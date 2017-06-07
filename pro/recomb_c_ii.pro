function recomb_c_ii, cii_rc_data, h_i_aeff_data, temp, dens, wavelength, iobs
;+
; NAME:
;     recomb_c_ii
; PURPOSE:
;     return the recombination coefficients of C II
;     from Davey et al. (2000) 2000A&AS..142...85D
; EXPLANATION:
;
; CALLING SEQUENCE:
;     ciiRLs=recomb_c_ii(tempi, densi, Abund)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of C II
;          ciiRLstructure
;          { Wave:0.0, 
;            Int:0.0, Obs:0.0, 
;            Abundance:0.0}
; REVISION HISTORY:
;     Recombination coefficients for C Ii lines from
;     Davey et al. 2000A&AS..142...85D
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
;     Originally added by Yong Zhang to MOCASSIN, 2003/02
;     Converted to IDL code by A. Danehkar, 10/05/2013
;     Integration with AtomNeb, A. Danehkar, 15/04/2017
;- 

  ; ciiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  TEh2=double(temp)
  NEh2=double(dens)
  abund=1.0
  nlines = 57 
  hbeta_aeff= (10.0^gamma4861(h_i_aeff_data,TEh2,NEh2))*double(4861.33/(h_Planck*c_Speed*1.e8)) 
  
  ; ciiRLs=REPLICATE(ciiRLstructure, nlines)
  
  lamb=double(0.0)
  a=double(0.0)
  b=double(0.0)
  c=double(0.0)
  d=double(0.0)
  f=double(0.0)
  aeff=double(0.0)
  br=double(1.0)
  temp4 = temp/10000.0
  loc1=where(abs(cii_rc_data.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(cii_rc_data[loc1].Wavelength)
    loc1=where(cii_rc_data.Wavelength eq  Wavelength_min)
  endif
  lamb=cii_rc_data[loc1].Wavelength
  a=cii_rc_data[loc1].a
  b=cii_rc_data[loc1].b
  c=cii_rc_data[loc1].c
  d=cii_rc_data[loc1].d
  f=cii_rc_data[loc1].f
  aeff = 1.0e-14 * (a*(temp4^f)) 
  aeff=aeff*(1. + (b*(1.-temp4)) + (c*((1.-temp4)^2) ) + (d * ((1.-temp4)^3) ) ) 
  ciiRLs_Int = 100.0*(aeff/hbeta_aeff)*br*(4861.33/lamb)*abund 
  
  abund=iobs/ciiRLs_Int
  return,abund
end
