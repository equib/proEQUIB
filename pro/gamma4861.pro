function gamma4861, h_i_aeff_data, temperature, density
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
;     gamm4861_theory = gamma4861(temperature,density)
;
; INPUTS:
;     temperature -     electron temperature in K
;     density -     electron density in cm-3
; RETURN:  Log10 (gamm(H Beta))
;
; REVISION HISTORY:
;     from Table 4.4 in D. E. Osterbrock & 
;     G. J. Ferland, Astrophysics of Gaseius Nebulae 
;     and Active Galactic Nuclei, 2nd Ed., 2006
;     IDL code by A. Danehkar, 31/08/2012
;     Integration with AtomNeb, A. Danehkar, 11/03/2017
;- 
  ;h_a_col= find_aeff_sh95_column(3, 2)
  linenum= find_aeff_sh95_column(4, 2, 25)

  TEh2=double(temperature)
  NEh2=double(density)
  line1=long(linenum-1)
  emissivity=double(0.0)
  
  hi_ems=dblarr(10,13)
  temp1=dblarr(302)
  temp_grid=[500.,1000.,3000.,5000.,7500.,10000.,12500.,15000.,20000.,30000.]
  
  nlines = 130 
  
  for i=0, nlines-1 do begin 
    temp1=h_i_aeff_data[*,i]
    tpos=nint((where(temp_grid eq temp1[1])))
    npos=nint(alog10(temp1[0])-2)
    hi_ems[tpos,npos]=temp1[line1];temp[2:45]
  endfor
  
  ; restrict to the density & temperature ranges to 2012MNRAS.425L..28P  
  if (NEh2 lt 1.e2) then NEh2=1.e2
  if (NEh2 gt 1.e14) then NEh2=1.e14
  if (TEh2 lt 500.) then TEh2=500.
  if (TEh2 gt 30000.) then TEh2=30000.
  
  ; get logarithmic density
  dens_log=alog10(NEh2)
  
  dens_grid=double(indgen(13) + 2)

  hi_ems1=hi_ems[*,*]
  ; Bilinearly interpolate density & temperature
  emiss_log = _interp2d(hi_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)

  ;logems = alog10(hr_tmp/double(4861.33/1.98648E-08))
  ;hb_ems = 10.0^logems
  ;hb_ems=alog10(hb_emissivity(temp, density))

  ;hb_ems=alog10(emiss_log/double(4861.33/1.98648E-08))
  hb_ems = alog10(emiss_log)
  
  return, hb_ems
end
