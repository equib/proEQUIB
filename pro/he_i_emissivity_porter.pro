function he_i_emissivity_porter, he_i_aeff_data, h_i_aeff_data, temperature, density, linenum, iobs
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
;     Integration with AtomNeb, A. Danehkar, 20/03/2017
;- 
  
  TEh2=double(temperature)
  NEh2=double(density)
  line1=long(linenum-1)
  emissivity=double(0.0)
  
  ; hei_ems=read_porter()
  hei_ems=dblarr(21,14) ;(21,14,44)
  temp1=dblarr(46)
  
  nlines = 294 
  
  for i=0, nlines-1 do begin 
    temp1=he_i_aeff_data[*,i]
    tpos=nint((temp1[0]/1000)-5)
    npos=nint(alog10(temp1[1])-1)
    hei_ems[tpos,npos]=temp1[line1+2];temp[2:45]
  endfor

  
  ; restrict to the density & temperature ranges to 2012MNRAS.425L..28P  
  if (NEh2 lt 1.e1) then NEh2=1.e1
  if (NEh2 gt 1.e14) then NEh2=1.e14
  if (TEh2 lt 5000) then TEh2=5000.
  if (TEh2 gt 25000) then TEh2=25000.
  
  ; get logarithmic density
  dens_log=alog10(NEh2)
  
  dens_grid=double(indgen(14) + 1)
  temp_grid=double(1000*(indgen(21) + 5))

  hei_ems1=hei_ems[*,*]
  ; Bilinearly interpolate density & temperature
  emiss_log = equib_interp2d(hei_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)
    
  ; wavl=he_i_aeff_data_Wavelength[line1]
  emissivity= 100.*10.^(emiss_log-gamma4861(h_i_aeff_data,TEh2,NEh2))
  
  abund=iobs/emissivity
  return,abund
end
