function gamma4861, temp, dens
;+
; NAME:
;     gamm4861
; PURPOSE:
;     determine the value of Log10 (gamm(H Beta))
;     = Log10( 4*Pai*j(HBeta)/NpNe) at temperature Te and density Ne
;     Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41
; EXPLANATION:
;
; CALLING SEQUENCE:
;     gamm4861_theory = gamma4861(temp,dens)
;
; INPUTS:
;     temp -     electron temperature in K
;     dens -     electron density in cm-3
; RETURN:  Log10 (gamm(H Beta))
;
; REVISION HISTORY:
;     from Table 4.4 in D. E. Osterbrock & 
;     G. J. Ferland, Astrophysics of Gaseius Nebulae 
;     and Active Galactic Nuclei, 2nd Ed., 2006
;     IDL code by A. Danehkar, 31/08/2012
;- 
  hb_ems=alog10(hb_emissivity(temp, dens))
  return, hb_ems
end
