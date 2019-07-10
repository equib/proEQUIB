; docformat = 'rst'

function calc_emiss_o_ii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, $
                      o_ii_rc_br=o_ii_rc_br, o_ii_rc_data=o_ii_rc_data
;+
;     This function calculates the emissivity 
;     for the given wavelength of O II recombination line 
;     by using the recombination coefficients from  
;     Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.
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
;     o_ii_rc_br    :     in, required, type=array/object
;                         O II branching ratios (Br)
;     o_ii_rc_data  :     in, required, type=array/object
;                         O II recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> 
;     IDL> atom='o'
;     IDL> ion='iii' ; O II
;     IDL> o_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> o_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> wavelength=4613.68
;     IDL> emiss_o_ii=calc_emiss_o_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, $
;     IDL>                               o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data, $
;     IDL>                               h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'Emissivity:', emiss_o_ii
;        Emissivity:   5.9047319e-27
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
;     Based on recombination coefficients for O II lines from
;     Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.
;     
;     Adopted from MIDAS script Roii.prg written by X.W.Liu.
;     
;     Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
;                       Ercolano et al. 2005MNRAS.362.1038E.
;     
;     10/05/2013, A. Danehkar, Translated to IDL code.
;     
;     25/04/2017, A. Danehkar, Integration with AtomNeb.
;     
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_o_ii_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_o_ii_rl().
;-

  ;  oiiRLstructure ={Wave:double(0.0), $ ;REAL*8
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
  if keyword_set(o_ii_rc_data) eq 0 then begin 
    print,'O II recombination coefficients (o_ii_rc_data) are not set'
    return, 0
  endif
  if keyword_set(o_ii_rc_br) eq 0 then begin 
    print,'O II branching ratios (o_ii_rc_br) are not set'
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
  Br_A=double(0.0)
  Br_B=double(0.0)
  Br_C=double(0.0)
  g1=double(0.0)
  g2=double(0.0)
  temp4 = temperature/10000.0
  loc1=where(abs(o_ii_rc_br.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(o_ii_rc_br[loc1].Wavelength)
    loc1=where(o_ii_rc_br.Wavelength eq  Wavelength_min)
  endif
  
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then loc1=min(loc1)
  
  Wave=o_ii_rc_br[loc1].Wavelength
  Br_A=o_ii_rc_br[loc1].Br_A
  Br_B=o_ii_rc_br[loc1].Br_B
  Br_C=o_ii_rc_br[loc1].Br_C
  g1=o_ii_rc_br[loc1].g1
  g2=o_ii_rc_br[loc1].g2
  densi=double(density)
  log10ne=alog10(densi)
  
  case 1 of
     ;---------------------------------------
     (loc1 ge 0) and (loc1 le 182): begin
        ; 4f-3d transitions
        a = o_ii_rc_data[0].a4 ; 0.232
        b = o_ii_rc_data[0].b ;-0.92009
        c = o_ii_rc_data[0].c ; 0.15526
        d = o_ii_rc_data[0].d ; 0.03442
        ; an = [0.236, 0.232, 0.228, 0.222]
        an = [o_ii_rc_data[0].a2, o_ii_rc_data[0].a4, o_ii_rc_data[0].a5, o_ii_rc_data[0].a6]
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3] 
                        
        temp4 = temperature/10000.0
        aeff = 1.0e-14 * a * temp4^(b)
        aeff = aeff*(1. + c*(1. - temp4) + d * (1. - temp4)^2)
      
        ; data for 1000 < T < 5000 K, Ne = 100/cm3
        if (temp4 le 0.5) then begin
          a = 0.236
          b = -1.07552
          c = -0.04843
          aeff = 1.0e-14 * a * temp4^(b + c * alog(temp4))
        endif
        Br = Br_B
     end
     ;---------------------------------------
     (loc1 ge 183) and (loc1 le 218): begin
        ; 3d-3p ^4F transitions (Case A=B=C for a,b,c,d; Br diff. slightly, adopt Case B)
        a = o_ii_rc_data[1].a4 ; 0.876
        b = o_ii_rc_data[1].b ; -0.73465
        c = o_ii_rc_data[1].c ; 0.13689
        d = o_ii_rc_data[1].d ; 0.06220
        ; an = [0.876, 0.876, 0.877, 0.880] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a)
        an = [o_ii_rc_data[1].a2, o_ii_rc_data[1].a4, o_ii_rc_data[1].a5, o_ii_rc_data[1].a6]
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3] 
        
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 0.878
          b = -0.86175
          c = -0.02470
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
        endif
        Br = Br_B
     end
     ;---------------------------------------
     (loc1 ge 219) and (loc1 le 309): begin
        ; 3d-3p ^4D, ^4P transitions 
        a = o_ii_rc_data[3].a4 ; 0.745
        b = o_ii_rc_data[3].b ; -0.74621
        c = o_ii_rc_data[3].c ; 0.15710
        d = o_ii_rc_data[3].d ; 0.07059
        ;an = [0.727,0.726,0.725,0.726] ; Case: A
        ; an = [0.747, 0.745, 0.744, 0.745] ; Case: B
        ;an = [0.769,0.767,0.766,0.766] ; Case: C
        an = [o_ii_rc_data[3].a2, o_ii_rc_data[3].a4, o_ii_rc_data[3].a5, o_ii_rc_data[3].a6]
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3] 
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 0.747
          b = -0.89382
          c = -0.02906
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
        endif
        Br = Br_B
     end
     ;---------------------------------------
     (loc1 ge 310) and (loc1 le 327): begin
        ; 3d-3p ^2F transitions 
        a = o_ii_rc_data[5].a4 ; 0.745
        b = o_ii_rc_data[5].b ; -0.74621
        c = o_ii_rc_data[5].c ; 0.15710
        d = o_ii_rc_data[5].d ; 0.07059
        ;an = [0.727, 0.726, 0.725, 0.726] ; Case: A
        ;an = [0.747,0.745,0.744,0.745] ; Case: B
        ;an = [0.769,0.767,0.766,0.766] ; Case: C
        an = [o_ii_rc_data[5].a2, o_ii_rc_data[5].a4, o_ii_rc_data[5].a5, o_ii_rc_data[5].a6]
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 0.747
          b = -0.89382
          c = -0.02906
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
        endif
        Br = Br_A
     end
     ;---------------------------------------
     (loc1 ge 328) and (loc1 le 357): begin
        ; 3d-3p ^2D transitions 
        a = o_ii_rc_data[11].a4 ; 0.601
        b = o_ii_rc_data[11].b ; -0.79533
        c = o_ii_rc_data[11].c ; 0.15314
        d = o_ii_rc_data[11].d ; 0.05322
        ;an = [0.603, 0.601, 0.600, 0.599] ; Case: A
        ;an = [0.620,0.618,0.616,0.615] ; Case: C
        an = [o_ii_rc_data[11].a2, o_ii_rc_data[11].a4, o_ii_rc_data[11].a5, o_ii_rc_data[11].a6]
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 0.603
          b = -0.94025
          c = -0.03467
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
        endif
        Br = Br_A
     end
     ;---------------------------------------
     (loc1 ge 358) and (loc1 le 387): begin
        ; 3d-3p ^2P transitions 
        a = o_ii_rc_data[13].a4 ; 0.524
        b = o_ii_rc_data[13].b ; -0.78448
        c = o_ii_rc_data[13].c ; 0.13681
        d = o_ii_rc_data[13].d ; 0.05608
        ;an = [0.526, 0.524, 0.523, 0.524] ; Case: A
        ;an = [0.538,0.536,0.535,0.536] ; Case: C
        an = [o_ii_rc_data[13].a2, o_ii_rc_data[13].a4, o_ii_rc_data[13].a5, o_ii_rc_data[13].a6]
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 0.526
          b = -0.91758
          c = -0.03120
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
        endif
        Br = Br_A
     end
     ;---------------------------------------
     (loc1 ge 388) and (loc1 le 395): begin
        ; 3p-3s ^4D - ^4P transitions 
        ;an = [34.7,34.9,35.1,35.0] ;a  Case: A
        ;a =  36.2
        ;b =  -0.749
        ;c =  0.023
        ;d =  0.074
        an = [36.0, 36.2, 36.4, 36.3] ;a  Case: B
        a =  36.2
        b =  -0.736
        c =  0.033
        d =  0.077
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 36.288
          b = -0.75421
          c = 0.02883
          d = 0.01213
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
        endif
        Br = Br_B
     end
     ;---------------------------------------
     (loc1 ge 396) and (loc1 le 402): begin
        ; 3p-3s ^4P - ^4P transitions 
        ;an = [10.4,10.4,10.5,10.4] ; Case: A
        ;a =  10.4
        ;b =  -0.721
        ;c =  0.073
        ;d =  0.072
        an = [14.6, 14.6, 14.7, 14.6] ; Case: B
        a =  14.6
        b =  -0.732
        c =  0.081
        d =  0.066
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 14.656
          b = -0.80449
          c = 0.00018
          d = 0.00517
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
        endif
        Br = Br_B
     end
     ;---------------------------------------
     (loc1 ge 403) and (loc1 le 405): begin
        ; 3p-3s ^4S - ^4P transitions 
        ;an = [0.90,0.90,0.90,1.00] ; Case: A
        ;a =  0.90
        ;b =  -0.485
        ;c =  -0.047
        ;d =  0.140
        an = [4.80, 4.90, 4.90, 4.90] ; Case: B
        a =  4.90
        b =  -0.730
        c =  -0.003
        d =  0.057
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 4.8340
          b = -0.71947
          c = 0.02544
          d = 0.00936
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
        endif
        Br = Br_B
     end
     ;---------------------------------------
     (loc1 ge 406) and (loc1 le 408): begin
        ; 3p-3s ^2D - ^2P transitions 
        an = [2.40, 2.40, 2.50, 2.60] ; Case: A
        a =  2.40
        b =  -0.550
        c =  -0.051
        d =  0.178
        ;an = [14.5,14.6,14.5,14.3] ; Case: C
        ;a =  14.6
        ;b =  -0.736
        ;c =  0.068
        ;d =  0.066
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 2.3616
          b = -0.46263
          c = 0.14697
          d = 0.03856
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
        endif
        Br = Br_A
     end
     ;---------------------------------------
     (loc1 ge 409) and (loc1 le 412): begin
        ; 3p-3s ^2P - ^2P transitions 
        an = [1.10, 1.20, 1.20, 1.20] ; Case: A
        a =  1.20
        b =  -0.523
        c =  -0.044
        d =  0.173
        ;an = [1.30,1.40,1.40,1.40] ; Case: C
        ;a =  1.40
        ;b =  -0.565
        ;c =  -0.042
        ;d =  0.158
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 1.1198
          b = -0.44147
          c = 0.13837
          d = 0.03191
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
        endif
        Br = Br_A
     end
     ;---------------------------------------
     (loc1 ge 413) and (loc1 le 414): begin
        ; 3p-3s ^2S - ^2P transitions 
        an = [0.40, 0.40, 0.40, 0.40] ; Case: A
        a =  0.40
        b =  -0.461
        c =  -0.083
        d =  0.287
        ;an = [0.50,0.50,0.50,0.60] ; Case: C
        ;a =  0.50
        ;b =  -0.547
        ;c =  -0.074
        ;d =  0.244
        if (log10ne le 2) then a = an[0] $
        else if (log10ne gt 2 and log10ne le 4) then a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.) $
          else if (log10ne gt 4 and log10ne le 5) then a = an[1] + (an[2] - an[1]) * (log10ne - 2.) $
               else if (log10ne gt 5) and (log10ne le 6) then a = an[2] + (an[3] - an[2]) * (log10ne - 2.) $
                    else a = an[3]
        aeff = 1.e-14*a*temp4^(b) 
        aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
        ; for 1000 < T < 5000 K, Ne = 100/cm3
        if ( temp4 le 0.5) then begin
          a = 0.3922
          b = -0.35043
          c = 0.26366
          d = 0.06666
          aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
        endif
        Br = Br_A
     end
  else: print, 'wavelength has an illegal value.'
  endcase
  ; Ems1 = aeff * (h_Planck*c_Speed*1.e8) /Wave * g2 * Br_A;[loc1]
  ; oiiRLs_Int = 100. * Ems1 / hbeta_ems * abund 
  ;abund=line_flux/oiiRLs_Int
  
  emissivity=(double(aeff*g2*Br)/double(Wave))*double(h_Planck*c_Speed*1.e8)
  
  return,emissivity
end
