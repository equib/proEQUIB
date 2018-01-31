function he_ii_emissivity, he_ii_aeff_data, h_i_aeff_data, temperature, density, iobs
;+
; NAME:
;     he_i_emissivity_porter
; PURPOSE:
;     read helium emissivities from Porter et al.
;     2012MNRAS.425L..28P
; EXPLANATION:
;
; CALLING SEQUENCE:
;     heidata=he_i_emissivity_porter(temperature, density, linenum)
;
; INPUTS:
;     temperature -     electron temperature in K
;     density -     electron density in cm-3
;     linenum -  line number
; RETURN:  emissivity of He I line
;
; REVISION HISTORY:
;     Improved He I emissivities in the case B
;     from Porter et al. 2012MNRAS.425L..28P
;     IDL code by A. Danehkar, 15/12/2013
;     Integration with AtomNeb, A. Danehkar, 02/04/2017
;- 
  
  TEh2=double(temperature)
  NEh2=double(density)
  emiss_log=gamma4686(he_ii_aeff_data,TEh2, NEh2)
  
  ; wavl=he_i_aeff_data_Wavelength[line1]
  emissivity= 100.*10.^(emiss_log-gamma4861(h_i_aeff_data,TEh2,NEh2))
  
  abund=iobs/emissivity
  return,abund
end
