function hb_eff_rec_coef, temp, dens
;+
; NAME:
;     hb_eff_rec_coef
; PURPOSE:
;     return effective recombination coefficient 
;     of Hydrogen Beta Balmer for
;     given electron temperature and density
;     Table 4.4 in D. E. Osterbrock & G. J. Ferland, 
;     Astrophysics of Gaseius Nebulae and
;     Active Galactic Nuclei, 2nd Ed., 2006
; EXPLANATION:
;
; CALLING SEQUENCE:
;     hb_eff = hb_rec_coef(temp, dens)
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
  hr_grid = [[5.37, 5.43, 5.59], [3.02, 3.03, 3.07], [1.61, 1.61, 1.62]] ; alpha_hb_eff
  
  hr_grid = 1.0e-14 * hr_grid ; cm3 s-1
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
  
  return, hr_tmp
end
