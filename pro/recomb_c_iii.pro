function recomb_c_iii, ciii_rc_data, h_i_aeff_data, temperature, density, wavelength, iobs
;+
; NAME:
;     recomb_c_iii
; PURPOSE:
;     return the recombination coefficients of C III
;     lines from Pequignot et al. 1991A&A...251..680P
; EXPLANATION:
;
; CALLING SEQUENCE:
;     xiiiRLs=recomb_c_iii(temperature, density, Abund)
;
; INPUTS:
;     temperature  - electron temperature in K
;     density  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of C III
;          ciiiRLstructure
;          { Wave:0.0, 
;            Int:0.0, Obs:0.0, 
;            abundance:0.0, Ion:''}
; REVISION HISTORY:
;     Total and effective radiative recombination coefficients
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P
;     IDL code by A. Danehkar, 10/05/2013
;     Integration with AtomNeb, A. Danehkar, 06/04/2017
;- 
  common share1, Atomic_Data_Path
  
  ; ciiiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  TEh2=double(temperature)
  NEh2=double(density)
  abund=1.0

  nlines = 4
  hbeta_aeff= (10.0^gamma4861(h_i_aeff_data,TEh2,NEh2))*double(4861.33/(h_Planck*c_Speed*1.e8)) 
  
  ;ciiiRLs=REPLICATE(ciiiRLstructure, nlines)

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
  loc1=where(abs(ciii_rc_data.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(ciii_rc_data[loc1].Wavelength)
    loc1=where(ciii_rc_data.Wavelength eq  Wavelength_min)
  endif
  lamb=ciii_rc_data[loc1].Wavelength
  a=ciii_rc_data[loc1].a
  b=ciii_rc_data[loc1].b
  c=ciii_rc_data[loc1].c
  d=ciii_rc_data[loc1].d
  br=ciii_rc_data[loc1].br
  ; equation (1) in 1991A&A...251..680P
  aeff = 1.0e-13 * z * br
  aeff=aeff*(a*(temp4^b))/(1.+c*(temp4^d)) 
  ciiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund 
  
  abund=iobs/ciiiRLs_Int
  return,abund
end
