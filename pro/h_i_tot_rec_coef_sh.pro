function h_i_tot_rec_coef_sh, temp, dens, case_name=case_name
;+
; NAME:
;     h_i_tot_rec_coef_sh
; PURPOSE:
;     return total hydrogen recombination coefficient 
;     given electron temperature and density
;     Storey & Hummer, MNRAS, 272, 41S, 1995
; EXPLANATION:
;
; CALLING SEQUENCE:
;     hb_eff = h_i_tot_rec_coef_sh(temp, dens)
;
; INPUTS:
;     temp -     electron temperature in K
;     dens -     electron density in cm-3
; RETURN:  Hydrogen Beta Balmer emissivity (Case B)
;
; REVISION HISTORY:
;     Storey & Hummer, MNRAS, 272, 41S, 1995
;     1995MNRAS.272...41S
;     IDL code by A. Danehkar, 31/08/2012
;-  
  case_num=0
  if keyword_set(case_name) then begin
    case case_name of
       'A': case_num = 1
       'B': case_num = 0
       else: begin
        print, 'assume Case B'
        exit
       endelse
    endcase
  endif else begin
    case_num=0
  endelse
  
  if case_num eq 1 then begin
    ; Case A
    dens_grid=[1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9, 1.e10 ]  
    temp_grid=[5.e2, 1.e3, 3.e3, 5.e3, 7.5e3, 1.e4, 1.25e4, 1.5e4, 2.e4, 3.e4]
    ; alpha_hb_eff ; cm3 s-1
    hr_grid = [[3.251e-12, 3.351e-12, 3.536e-12, 3.880e-12, 4.552e-12, 5.943e-12, 9.129e-12, 1.769e-11, 4.694e-11], $
               [2.038e-12, 2.069e-12, 2.125e-12, 2.229e-12, 2.424e-12, 2.802e-12, 3.575e-12, 5.308e-12, 9.817e-12], $
               [9.690e-13, 9.735e-13, 9.819e-13, 9.973e-13, 1.026e-12, 1.079e-12, 1.179e-12, 1.374e-12, 1.778e-12], $
               [6.809e-13, 6.827e-13, 6.861e-13, 6.923e-13, 7.038e-13, 7.252e-13, 7.651e-13, 8.406e-13, 9.890e-13], $
               [5.120e-13, 5.128e-13, 5.145e-13, 5.174e-13, 5.230e-13, 5.333e-13, 5.524e-13, 5.883e-13, 6.569e-13], $
               [4.169e-13, 4.174e-13, 4.183e-13, 4.201e-13, 4.234e-13, 4.294e-13, 4.408e-13, 4.619e-13, 5.018e-13], $
               [3.547e-13, 3.550e-13, 3.557e-13, 3.568e-13, 3.590e-13, 3.630e-13, 3.705e-13, 3.844e-13, 4.107e-13], $
               [3.104e-13, 3.106e-13, 3.111e-13, 3.119e-13, 3.134e-13, 3.163e-13, 3.216e-13, 3.315e-13, 3.501e-13], $
               [2.507e-13, 2.509e-13, 2.511e-13, 2.516e-13, 2.525e-13, 2.541e-13, 2.572e-13, 2.630e-13, 2.737e-13], $
               [1.843e-13, 1.844e-13, 1.845e-13, 1.847e-13, 1.851e-13, 1.858e-13, 1.872e-13, 1.898e-13, 1.947e-13]]
  
    ; Linearly interpolate extinction law in 1/lam
    if temp lt temp_grid[0] or temp gt temp_grid[9] then begin
      print, 'ouside temperature range!'
      return, 0
    endif
    if dens lt dens_grid[0] or dens gt dens_grid[8] then begin
      print, 'ouside density range!'
      return, 0
    endif
    
    ; Linearly interpolate density
    hr_dns0=lin_interp(hr_grid[*,0], dens_grid, dens)
    hr_dns1=lin_interp(hr_grid[*,1], dens_grid, dens)
    hr_dns2=lin_interp(hr_grid[*,2], dens_grid, dens)
    hr_dns3=lin_interp(hr_grid[*,3], dens_grid, dens)
    hr_dns4=lin_interp(hr_grid[*,4], dens_grid, dens)
    hr_dns5=lin_interp(hr_grid[*,5], dens_grid, dens)
    hr_dns6=lin_interp(hr_grid[*,6], dens_grid, dens)
    hr_dns7=lin_interp(hr_grid[*,7], dens_grid, dens)
    hr_dns8=lin_interp(hr_grid[*,8], dens_grid, dens)
    hr_dns9=lin_interp(hr_grid[*,9], dens_grid, dens)
    
      
    ; Linearly interpolate temperature
    hr_dns_grid=[hr_dns0, hr_dns1, hr_dns2, hr_dns3, hr_dns4, hr_dns5, hr_dns6, hr_dns7, hr_dns8, hr_dns9]  
    hr_tmp=lin_interp(hr_dns_grid, temp_grid, temp)
  endif else begin
    ; Case B
    dens_grid=[1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14 ]  
    temp_grid=[5.e2, 1.e3, 3.e3, 5.e3, 7.5e3, 1.e4, 1.25e4, 1.5e4, 2.e4, 3.e4]
    ; alpha_hb_eff ; cm3 s-1
    hr_grid =  [[2.493E-12, 2.573E-12, 2.720E-12, 2.998E-12, 3.542E-12, 4.681E-12, 7.330E-12, 1.462E-11, 4.044E-11, 1.700E-10, 1.119E-09, 9.959E-09, 9.756E-08], $
                [1.512E-12, 1.535E-12, 1.579E-12, 1.658E-12, 1.810E-12, 2.106E-12, 2.717E-12, 4.111E-12, 7.823E-12, 2.037E-11, 7.981E-11, 4.937E-10, 4.266E-09], $
                [6.708E-13, 6.740E-13, 6.798E-13, 6.907E-13, 7.109E-13, 7.486E-13, 8.204E-13, 9.615E-13, 1.257E-12, 1.941E-12, 3.811E-12, 1.040E-11, 4.325E-11], $
                [4.522E-13, 4.534E-13, 4.556E-13, 4.597E-13, 4.674E-13, 4.816E-13, 5.083E-13, 5.592E-13, 6.601E-13, 8.726E-13, 1.371E-12, 2.755E-12, 7.792E-12], $
                [3.273E-13, 3.278E-13, 3.288E-13, 3.306E-13, 3.341E-13, 3.404E-13, 3.524E-13, 3.749E-13, 4.181E-13, 5.047E-13, 6.908E-13, 1.138E-12, 2.453E-12], $
                [2.585E-13, 2.588E-13, 2.594E-13, 2.604E-13, 2.623E-13, 2.658E-13, 2.724E-13, 2.848E-13, 3.083E-13, 3.541E-13, 4.477E-13, 6.551E-13, 1.200E-12], $
                [2.144E-13, 2.147E-13, 2.149E-13, 2.156E-13, 2.167E-13, 2.190E-13, 2.230E-13, 2.307E-13, 2.452E-13, 2.728E-13, 3.276E-13, 4.426E-13, 7.303E-13], $
                [1.836E-13, 1.837E-13, 1.839E-13, 1.843E-13, 1.851E-13, 1.866E-13, 1.893E-13, 1.944E-13, 2.040E-13, 2.221E-13, 2.571E-13, 3.278E-13, 5.036E-13], $
                [1.428E-13, 1.429E-13, 1.430E-13, 1.432E-13, 1.436E-13, 1.444E-13, 1.458E-13, 1.484E-13, 1.532E-13, 1.621E-13, 1.788E-13, 2.106E-13, 2.959E-13], $
                [9.911E-14, 9.913E-14, 9.917E-14, 9.924E-14, 9.937E-14, 9.962E-14, 1.001E-13, 1.009E-13, 1.025E-13, 1.054E-13, 1.103E-13, 1.190E-13, 1.528E-13]]
  
    ; Linearly interpolate extinction law in 1/lam
    if temp lt temp_grid[0] or temp gt temp_grid[9] then begin
      print, 'ouside temperature range!'
      return, 0
    endif
    if dens lt dens_grid[0] or dens gt dens_grid[12] then begin
      print, 'ouside density range!'
      return, 0
    endif
    
    ; Linearly interpolate density
    hr_dns0=lin_interp(hr_grid[*,0], dens_grid, dens)
    hr_dns1=lin_interp(hr_grid[*,1], dens_grid, dens)
    hr_dns2=lin_interp(hr_grid[*,2], dens_grid, dens)
    hr_dns3=lin_interp(hr_grid[*,3], dens_grid, dens)
    hr_dns4=lin_interp(hr_grid[*,4], dens_grid, dens)
    hr_dns5=lin_interp(hr_grid[*,5], dens_grid, dens)
    hr_dns6=lin_interp(hr_grid[*,6], dens_grid, dens)
    hr_dns7=lin_interp(hr_grid[*,7], dens_grid, dens)
    hr_dns8=lin_interp(hr_grid[*,8], dens_grid, dens)
    hr_dns9=lin_interp(hr_grid[*,9], dens_grid, dens)
      
    ; Linearly interpolate temperature
    hr_dns_grid=[hr_dns0, hr_dns1, hr_dns2, hr_dns3, hr_dns4, hr_dns5, hr_dns6, hr_dns7, hr_dns8, hr_dns9]  
    hr_tmp=lin_interp(hr_dns_grid, temp_grid, temp)
  endelse
  return, hr_tmp
end
