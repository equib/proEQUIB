function recomb_ne_ii, neii_rc_data, h_i_aeff_data, temperature, density, wavelength, iobs
;+
; NAME:
;     recomb_ne_ii
; PURPOSE:
;     return the recombination coefficients of Ne II line
;     from Kisielius et al. (1998) & Storey (unpublished)
; EXPLANATION:
;
; CALLING SEQUENCE:
;     neiiRLs=recomb_ne_ii(neii_rc_data, h_i_aeff_data, 
;     temperature, density, wavelength, iobs)
;
; INPUTS:
;     temperature  - electron temperature in K
;     density  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of Ne II
;          neiiRLstructure
;          { Wave:0.0, 
;            Int:0.0, Obs:0.0, 
;            Abundance:0.0}
; REVISION HISTORY:
;     from Kisielius et al. 1998A&AS..133..257K
;     & Storey (unpublished)
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
;     scripts added by Yong Zhang to MOCASSIN, 2003/02
;     Converted to IDL code by A. Danehkar, 10/05/2013
;     Integration with AtomNeb, A. Danehkar, 10/04/2017
;- 
  common share1, Atomic_Data_Path
  
  ; neiiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  TEh2=double(temperature)
  NEh2=double(density)
  abund=1.0

  nlines = 38
  hbeta_aeff= (10.0^gamma4861(h_i_aeff_data,TEh2,NEh2))*double(4861.33/(h_Planck*c_Speed*1.e8)) 
  
  ; neiiRLs=REPLICATE(neiiRLstructure, nlines)

  lamb=double(0.0)
  a=double(0.0)
  b=double(0.0)
  c=double(0.0)
  d=double(0.0)
  f=double(0.0)
  br=double(0.0)
  aeff=double(0.0)
  Ion=''
  
  z = 3.0 ; ion level c^3+
  ; equation (1) in 1991A&A...251..680P
  temp4 = temperature/10000.0
  loc1=where(abs(neii_rc_data.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(neii_rc_data[loc1].Wavelength)
    loc1=where(neii_rc_data.Wavelength eq  Wavelength_min)
  endif
  lamb=neii_rc_data[loc1].Wavelength
  a=neii_rc_data[loc1].a
  b=neii_rc_data[loc1].b
  c=neii_rc_data[loc1].c
  d=neii_rc_data[loc1].d
  f=neii_rc_data[loc1].f
  br=neii_rc_data[loc1].br
  ; equation (1) in 1991A&A...251..680P
  aeff=1.0e-14 * (a*(temp4^f)) * br
  aeff=aeff*(1.+b*(1.0-temp4)+c*(1.0-temp4)^2+d*(1.0-temp4)^3)
  neiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund
  
  abund=iobs/neiiRLs_Int
  return,abund
end
