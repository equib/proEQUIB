pro get_aeff_hb, temp, dens, aeff_hb, em_hb
;+
; NAME:
;     get_aeff_hb
; PURPOSE:
;     determine the value of Aeff and emissivity of H_beta
;     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41
; EXPLANATION:
;
; CALLING SEQUENCE:
;     get_aeff_hb, temp, dens, aeff_hb, em_hb
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     aeff_hb - effective recombination coefficient of H_beta
;     em_hb   - emissivity of H_beta
;
; REVISION HISTORY:
;     IDL code by A. Danehkar, 10/05/2013
;- 
  aeff_hb = hb_eff_rec_coef(temp, dens)
  logem = alog10(aeff_hb) - 11.38871 ; = log10(hc/lambda in cgs)
  em_hb = 10.0^logem
end
