; docformat = 'rst'

function calc_emiss_n_ii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, $
                      n_ii_rc_br=n_ii_rc_br, n_ii_rc_data=n_ii_rc_data
;+
;     This function calculates the emissivity 
;     for the given wavelength of N II recombination line 
;     by using the recombination coefficients from  
;     Escalante & Victor 1990ApJS...73..513E.
;
; :Returns:
;    type=double. This function returns the line emissivity.
;
; :Keywords:
;     temperature   :     in, required, type=float
;                         electron temperature
;     density       :     in, required, type=float
;                         electron density
;     wavelength    :     in, required, type=float
;                         Line Wavelength in Angstrom
;     n_ii_rc_br    :     in, required, type=array/object
;                         N II branching ratios (Br)
;     n_ii_rc_data  :     in, required, type=array/object
;                         N II recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> atom='h'
;     IDL> ion='ii' ; H I
;     IDL> h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> h_i_aeff_data=h_i_rc_data[0].Aeff
;     IDL> atom='n'
;     IDL> ion='iii' ; N II
;     IDL> n_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> n_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> wavelength=4442.02
;     IDL> emiss_n_ii=calc_emiss_n_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, $
;     IDL>                               n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data, $
;     IDL>                               h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'Emissivity:', emiss_n_ii
;        Emissivity:   3.0397397e-26
;
; :Categories:
;   Abundance Analysis, Recombination Lines, Emissivity
;
; :Dirs:
;  ./
;      Main routines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.3.0
;
; :History:
;     Based on Effective recombination coefficients for N II lines from
;     Escalante & Victor 1990ApJS...73..513E.
;     
;     Adopted from MIDAS Rnii script written by X.W.Liu.
;     
;     Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
;                       Ercolano et al. 2005MNRAS.362.1038E.
;     
;     10/05/2013, A. Danehkar, Translated to IDL code.
;     
;     25/04/2017, A. Danehkar, Integration with AtomNeb.
;     
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_n_ii_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_n_ii_rl().
;-

  ;  niiRLstructure ={Wave:double(0.0), $ ;REAL*8
  ;              Int:double(0.0),  $  
  ;              Obs:double(0.0), $ 
  ;              abundance:double(0.0), $             
  ;              g1:long(0), $ ;INTEGER 
  ;              g2:long(0), $ ;INTEGER
  ;              Mult1:'', $ ;CHARACTER*7
  ;              Term1:'', $ ;CHARACTER*9
  ;              Term2:'' $ ;CHARACTER*9
  ;              } 
  
  h_Planck = double(6.62606957e-27) ; erg s
  c_Speed = double(2.99792458e10) ; cm/s 
  
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(n_ii_rc_data) eq 0 then begin 
    print,'N II recombination coefficients (n_ii_rc_data) are not set'
    return, 0
  endif
  if keyword_set(n_ii_rc_br) eq 0 then begin 
    print,'N II branching ratios (n_ii_rc_br) are not set'
    return, 0
  endif
  if keyword_set(wavelength) eq 0 then begin 
    print,'Wavelength is not given'
    return, 0
  endif
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif
       
  Wave=double(0.0)
  RL_br=double(0.0)
  g1=double(0.0)
  g2=double(0.0)
  temp4 = temperature/10000.0
  loc1=where(abs(n_ii_rc_br.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(n_ii_rc_br[loc1].Wavelength)
    loc1=where(n_ii_rc_br.Wavelength eq  Wavelength_min)
  endif
  Wave=n_ii_rc_br[loc1].Wavelength
  RL_br=n_ii_rc_br[loc1].Br
  
  case 1 of
     ;---------------------------------------
     (loc1 ge 0) and (loc1 le 5): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p E1 3P* - 3D  : 03  row 
        ;i = 0002    ;case A
        i = 0003     ;case B
     end
     ;---------------------------------------
     (loc1 ge 6) and (loc1 le 8): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3S     : 04 row
        ;i = 0004    ;Case: A
        i = 0005     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 9) and (loc1 le 14): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3P     : 05 row
        ;i = 0006    ;Case: A
        i = 0007     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 15): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1P     : 08 row
        i = 0008     ;Case: A
        ;i = 0009    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 16): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1D     : 12 row
        i = 0010     ;Case: A
        ;i = 0011    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 17): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1S     : 13 row
        i = 0012     ;Case: A
        ;i = 0013    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 18): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1D*     : 15 row
        i = 0014     ;Case: A
        ;i = 0015    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 19): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1P*     : 17 row
        i = 0016     ;Case: A
        ;i = 0017    ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 20) and (loc1 le 25): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3F*     : 19 row
        ;i = 0018    ;Case: A
        i = 0019     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 26) and (loc1 le 32): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3D*     : 20 row
        ;i = 0020    ;Case: A
        i = 0021     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 33) and (loc1 le 38): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3P*     : 21 row
        ;i = 0022    ;Case: A
        i = 0023     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 39) and (loc1 le 44): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3D - 3P*     : 22 row
        ;i = 0024    ;Case: A
        i = 0025     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 45) and (loc1 le 47): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3S - 3P*     : 24 row
        ;i = 0026    ;Case: A
        i = 0027     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 48) and (loc1 le 50): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3S - 3P*     : 26 row
        ;i = 0028    ;Case: A
        i = 0029     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 51) and (loc1 le 56): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3D*     : 28 row
        ;i = 0030    ;Case: A
        i = 0031     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 57) and (loc1 le 62): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3P*     : 29 row
        ;i = 0032    ;Case: A
        i = 0033     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 63) and (loc1 le 68): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3P - 3P*     : 30 row
        ;i = 0034    ;Case: A
        i = 0035     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 69): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1D - 1F*     : 31 row
        i = 0036     ;Case: A
        ;i = 0037    ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 70) and (loc1 le 75): begin
        ; atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*).4p 3F* - 3D     : 36 row
        ;i = 0038     ;Case: A
        i = 0039     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 76) and (loc1 le 81): begin
        ; atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 3F* - 3G     : 39 row
        ;i = 0040    ;Case: A
        i = 0041     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 82): begin
        ; atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 1F* - 1G     : 58 row
        i = 0042     ;Case: A
        ;i = 0043    ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 83) and (loc1 le 88): begin
        ; atomic transitions: 3d 3D* - 4f 3F 4242     : 48 row
        ;i = 0044    ;Case: A
        i = 0045     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 89) and (loc1 le 94): begin
        ; atomic transitions: 3d 3P* - 4f 3D 4435     : 55 row
        ;i = 0046   ;Case: A
        i = 0047     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 95): begin
        ; atomic transitions: 3d 1D* - 4f 1F 4176     : 43 (RMT 42) row
        i = 0048     ;Case: A
        ;i = 0049     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 96): begin
        ; atomic transitions: 3d 1P* - 4f 1D 4677     : 61 (rmt 62) row
        i = 0049     ;Case: A
        ;i = 0051    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 97): begin
        ; Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
        ; atomic transitions: 3d 3F* - 4f 1G 4026     : 39b row
        ; Case: A 
        a0 = 0.108
        b0 = -0.754
        c0 = 2.587
        d0 = 0.719
        z0 = 2.
        br0 = 0.350   
        aeff = 1.e-13 * z0 * a0  * (temp4/z0^2)^(b0)
        aeff = aeff / (1. + c0 * (temp4/z0^2)^(d0)) * br0
     end
     ;---------------------------------------
     (loc1 eq 98): begin
        ; Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
        ; atomic transitions: 3d 1F* - 4f 3G 4552     : 58a row
        ; Case: A 
        a0 = 0.326
        b0 = -0.754
        c0 = 2.587
        d0 = 0.719
        z0 = 2.
        br0 = 0.074   
        aeff = 1.e-13 * z0 * a0  * (temp4/z0^2)^(b0)
        aeff = aeff / (1. + c0 * (temp4/z0^2)^(d0)) * br0
     end
  else: print, 'wavelength has an illegal value.'
  endcase
  if (loc1 ge 0) and (loc1 le 96) then begin
    aeff = 10.^(n_ii_rc_data[i-1].a + n_ii_rc_data[i-1].b * alog10(temp4) + n_ii_rc_data[i-1].c * alog10(temp4)^2) 
  endif
  
  ;niiRLs_Int = 100.0 * emissivity / hbeta_ems * abund 
  ;abund=line_flux/niiRLs_Int
  emissivity=(double(aeff*RL_br)/double(Wave))*double(h_Planck*c_Speed*1.e8)
  
  return,emissivity
end
