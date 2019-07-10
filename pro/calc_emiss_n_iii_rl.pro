; docformat = 'rst'

function calc_emiss_n_iii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, $
                      n_iii_rc_data=n_iii_rc_data
;+
;     This function calculates the emissivity 
;     for the given wavelength of N III recombination line 
;     by using the recombination coefficients from
;     Pequignot et al. 1991A&A...251..680P.
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
;     n_iii_rc_data  :     in, required, type=array/object
;                         N III recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_PPB91_file='/media/linux/proEQUIB/AtomNeb-idl/atomic-data-rc/rc_PPB91.fits'
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> 
;     IDL> atom='n'
;     IDL> ion='iv' ; N III
;     IDL> n_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> wavelength=4640.64
;     IDL> emiss_n_iii=calc_abund_n_iii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, $
;     IDL>                                 n_iii_rc_data=n_iii_rc_data)
;     IDL> print, 'Emissivity:', emiss_n_iii
;        Emissivity:   4.7908644e-24
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
;     Based on  effective radiative recombination coefficients for N III lines from
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.
;     
;     10/05/2013, A. Danehkar, IDL code written.
;     
;     20/04/2017, A. Danehkar, Integration with AtomNeb.
;     
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_n_iii_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_n_iii_rl().
;-

  ; niiiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

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
  if keyword_set(n_iii_rc_data) eq 0 then begin 
    print,'N III recombination coefficients (n_iii_rc_data) are not set'
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
  br=double(0.0)
  aeff=double(0.0)
  Ion=''
  
  z = 3.0 ; ion level c^3+
  ; equation (1) in 1991A&A...251..680P
  temp4 = 1.0e-4 * temperature /z^2
  loc1=where(abs(n_iii_rc_data.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(n_iii_rc_data[loc1].Wavelength)
    loc1=where(n_iii_rc_data.Wavelength eq  Wavelength_min)
    temp2=size(loc1,/DIMENSIONS)
    if temp2[0] ne 1 then begin
      loc1=where(n_iii_rc_data.Wavelength eq  Wavelength_min and n_iii_rc_data.Case1 eq 'B')
    endif
  endif
  lamb=n_iii_rc_data[loc1].Wavelength
  a=n_iii_rc_data[loc1].a
  b=n_iii_rc_data[loc1].b
  c=n_iii_rc_data[loc1].c
  d=n_iii_rc_data[loc1].d
  br=n_iii_rc_data[loc1].br
  ; equation (1) in 1991A&A...251..680P
  aeff = 1.0e-13 * z * br
  aeff=aeff*(a*(temp4^b))/(1.+c*(temp4^d)) 
  ;niiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund 
  ;abund=line_flux/niiiRLs_Int
  emissivity=(double(aeff)/double(lamb))*double(h_Planck*c_Speed*1.e8)
  
  return,emissivity
end
