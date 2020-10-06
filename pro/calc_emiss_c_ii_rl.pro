; docformat = 'rst'

function calc_emiss_c_ii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, $
                      c_ii_rc_data=c_ii_rc_data
;+
;     This function calculates the emissivity 
;     for the given wavelength of C II recombination line 
;     by using the recombination coefficients from
;     from Davey et al. (2000) 2000A&AS..142...85D.
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
;     c_ii_rc_data  :     in, required, type=array/object
;                         C II recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> 
;     IDL> atom='c'
;     IDL> ion='iii' ; C II
;     IDL> c_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> wavelength=6151.43
;     IDL> emiss_c_ii=calc_emiss_c_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, $
;     IDL>                               c_ii_rc_data=c_ii_rc_data)
;     IDL> print, 'Emissivity:', emiss_c_ii
;        Emissivity:   5.4719511e-26
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
;     Based on recombination coefficients for C II lines from
;     Davey et al. 2000A&AS..142...85D.
;     
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.
;     
;     02/2003, Yong Zhang, added to MOCASSIN.
;     
;     10/05/2013, A. Danehkar, Translated to IDL code.
;     
;     15/04/2017, A. Danehkar, Integration with AtomNeb.
;
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_c_ii_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_c_ii_rl().
;-

  ; ciiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(c_ii_rc_data) eq 0 then begin 
    print,'C II recombination coefficients (c_ii_rc_data) are not set'
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
  
  lamb=double(0.0)
  a=double(0.0)
  b=double(0.0)
  c=double(0.0)
  d=double(0.0)
  f=double(0.0)
  aeff=double(0.0)
  br=double(1.0)
  temp4 = temperature/10000.0
  loc1=where(abs(c_ii_rc_data.Wavelength-wavelength) le 1.5)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(c_ii_rc_data[loc1].Wavelength)
    loc1=where(c_ii_rc_data.Wavelength eq  Wavelength_min)
  endif
  lamb=c_ii_rc_data[loc1].Wavelength
  a=c_ii_rc_data[loc1].a
  b=c_ii_rc_data[loc1].b
  c=c_ii_rc_data[loc1].c
  d=c_ii_rc_data[loc1].d
  f=c_ii_rc_data[loc1].f
  aeff = 1.0e-14 * (a*(temp4^f)) 
  aeff=aeff*(1. + (b*(1.-temp4)) + (c*((1.-temp4)^2) ) + (d * ((1.-temp4)^3) ) ) 
  ;ciiRLs_Int = 100.0*(aeff/hbeta_aeff)*br*(4861.33/lamb)*abund 
  ;abund=line_flux/ciiiRLs_Int
  emissivity=(double(aeff*br)/double(lamb))*double(h_Planck*c_Speed*1.e8)

  return,emissivity
end
