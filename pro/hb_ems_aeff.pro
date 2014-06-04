function hb_ems_aeff, temp, dens
;+
; NAME:
;     hb_ems_aeff
; PURPOSE:
;     determine the value of Aeff and emissivity of H_beta
;     Table 4.4 in D. E. Osterbrock & G. J. Ferland, 
;     Astrophysics of Gaseius Nebulae and
;     Active Galactic Nuclei, 2nd Ed., 2006
; EXPLANATION:
;
; CALLING SEQUENCE:
;     hbeta=hb_ems_aeff(temp, dens)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
; OUTPUTS:
;     {aeff:double(0.0), ems:double(0.0)}
;     hbeta.aeff - effective recombination coefficient of H_beta
;     hbeta.em   - emissivity of H_beta
;
; REVISION HISTORY:
;     IDL code by A. Danehkar, 10/05/2013
;- 
  hbeta ={ems:double(0.0), aeff:double(0.0)}
  hbeta.aeff = hb_eff_rec_coef(temp, dens)
  logems = alog10(hbeta.aeff/double(4861.33/1.98648E-08))
  hbeta.ems = 10.0^logems
  return, hbeta
end
