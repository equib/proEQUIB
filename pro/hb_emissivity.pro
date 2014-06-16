function hb_emissivity, temp, dens
;+
; NAME:
;     hb_emissivity
; PURPOSE:
;     return Hydrogen Beta Balmer emissivity for
;     given electron temperature and density
;     Table 4.4 in D. E. Osterbrock & G. J. Ferland, 
;     Astrophysics of Gaseius Nebulae and
;     Active Galactic Nuclei, 2nd Ed., 2006
; EXPLANATION:
;
; CALLING SEQUENCE:
;     h_i = hb_emissivity(temp, dens)
;
; INPUTS:
;     temp -     electron temperature in K
;     dens -     electron density in cm-3
; RETURN:  Hydrogen Beta Balmer emissivity (Case B)
;
; REVISION HISTORY:
;     from Table 4.4 in D. E. Osterbrock & 
;     G. J. Ferland, Astrophysics of Gaseius Nebulae 
;     and Active Galactic Nuclei, 2nd Ed., 2006
;     IDL code by A. Danehkar, 31/08/2012
;-  
  dens_grid=[1.0e2, 1.0e4, 1.0e6]  
  temp_grid=[5000.0, 10000.0, 20000.0]
  hr_grid = [[2.20, 2.22, 2.29], [1.23, 1.24, 1.25], [0.658, 0.659, 0.661]] ; 4*pi*j_Hb/(ne*np)
  
  hr_grid = 1.0e-25 * hr_grid ; erg cm3 s-1
  ; Linearly interpolate extinction law in 1/lam
  if temp lt temp_grid[0] or temp gt temp_grid[2] then begin
    print, 'ouside temperature range!'
    return, 0
  endif
  if dens lt dens_grid[0] or dens gt dens_grid[2] then begin
    print, 'ouside density range!'
    return, 0
  endif
  
  ; Linearly interpolate density
  hr_dns0=lin_interp(hr_grid[*,0], dens_grid, dens)
  hr_dns1=lin_interp(hr_grid[*,1], dens_grid, dens)
  hr_dns2=lin_interp(hr_grid[*,2], dens_grid, dens)
  
  ; Linearly interpolate temperature
  hr_dns_grid=[hr_dns0, hr_dns1, hr_dns2]  
  hr_tmp=lin_interp(hr_dns_grid, temp_grid, temp)
  
  logems = alog10(hr_tmp/double(4861.33/1.98648E-08))
  hb_ems = 10.0^logems
  return, hb_ems
end
