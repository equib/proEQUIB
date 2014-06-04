function he_i_emissivity_porter, temp, dens, linenum
;+
; NAME:
;     he_i_emissivity_porter
; PURPOSE:
;     read helium emissivities from Porter et al.
;     2012MNRAS.425L..28P
; EXPLANATION:
;
; CALLING SEQUENCE:
;     heidata=he_i_emissivity_porter(temp, dens, linenum)
;
; INPUTS:
;     temp -     electron temperature in K
;     dens -     electron density in cm-3
;     linenum -  line number
; RETURN:  emissivity of He I line
;
; REVISION HISTORY:
;     Improved He I emissivities in the case B
;     from Porter et al. 2012MNRAS.425L..28P
;     IDL code by A. Danehkar, 15/12/2013
;- 
  
  TEh2=double(temp)
  NEh2=double(dens)
  line1=long(linenum-1)
  emissivity=double(0.0)
  
  hei_ems=read_porter()
  
  ; restrict to the density & temperature ranges to 2012MNRAS.425L..28P  
  if (NEh2 lt 1.e1) then NEh2=1.e1
  if (NEh2 gt 1.e14) then NEh2=1.e14
  if (TEh2 lt 5000) then TEh2=5000.
  if (TEh2 gt 25000) then TEh2=25000.
  
  ; get logarithmic density
  dens_log=alog10(NEh2)
  
  dens_grid=double(indgen(14) + 1)
  temp_grid=double(1000*(indgen(21) + 5))

  hei_ems1=hei_ems[*,*,line1]
  ; Bilinearly interpolate density & temperature
  emiss_log = equib_interp2d(hei_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)

  emissivity= 100.*10.^(emiss_log-gamma4861(TEh2,NEh2))
  
  return,emissivity
end
