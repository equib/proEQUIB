function he_i_emissivity_smits, temp, dens, line
;+
; NAME:
;     he_i_emissivity_smits
; PURPOSE:
;     return helium emissivity for
;     given electron temperature and density
;     Smits D. P., 1996, MNRAS, 278, 683
;     1996MNRAS.278..683S
; EXPLANATION:
;
; CALLING SEQUENCE:
;     helium_emissivity = he_i_emissivity_smits(temp, dens, line)
;
; INPUTS:
;     temp -     electron temperature in K
;     dens -     electron density in cm-3
;     line -     Balmer line number
;        1 : He I 3889
;        2 : He I 4026
;        3 : He I 4387
;        4 : He I 4471
;        5 : He I 4922
;        6 : He I 5876
;        7 : He I 6678
;        8 : He I 7065
;        9 : He I 7281
;        10 : He I 10830

; RETURN:  Hydrogen Case B Balmer line ratio
;
; REVISION HISTORY:
;     Tables 1 and 2 in Smits 1996, MNRAS, 278, 683
;     1996MNRAS.278..683S
;     IDL code by A. Danehkar, 31/08/2012
;-  
  dens_grid=[1.0e2, 1.0e4, 1.0e6]
  temp_grid=[312.5, 625.0, 1250.0, 2500.0, 5000.0, 10000.0, 20000.0]  

  ; Table 1  in Smits 1996MNRAS.278..683S
  emiss4471_grid = [[9.364e-25, 5.854e-25, 3.578e-25, 2.104e-25, 1.174e-25, 6.146e-26, 3.001e-26], $
                  [1.049e-24, 6.195e-25, 3.674e-25, 2.129e-25, 1.179e-25, 6.155e-26, 3.001e-26], $
                  [1.512e-24, 7.501e-25, 4.048e-25, 2.228e-25, 1.202e-25, 6.196e-26, 3.010e-26]]

  ; Table 2  in Smits 1996MNRAS.278..683S
  ; Relative fluxes for He I 3889
  emiss3889_rel_grid = [[1.262, 1.340, 1.452, 1.617, 1.860, 2.215, 2.722], $
                   [1.261, 1.342, 1.455, 1.621, 1.865, 2.219, 2.727], $
                   [1.307, 1.382, 1.486, 1.646, 1.886, 2.238, 2.739]]

  ; Relative fluxes for He I 3889
  emiss4026_rel_grid = [[0.396, 0.406, 0.419, 0.434, 0.450, 0.466, 0.479], $
                    [0.401, 0.408, 0.420, 0.435, 0.451, 0.466, 0.479], $
                    [0.413, 0.418, 0.427, 0.439, 0.453, 0.467, 0.480]]
  
  ; Relative fluxes for He I 4387
  emiss4387_rel_grid = [[0.107, 0.110, 0.113, 0.116, 0.120, 0.123, 0.125], $
                  [0.109, 0.110, 0.113, 0.117, 0.120, 0.123, 0.125], $
                  [0.112, 0.113, 0.115, 0.118, 0.121, 0.124, 0.125]]
  
  ; Relative fluxes for He I 4922
  emiss4922_rel_grid = [[0.273, 0.273, 0.272, 0.271, 0.269, 0.266, 0.263], $
                  [0.274, 0.273, 0.272, 0.270, 0.269, 0.266, 0.263], $
                  [0.274, 0.273, 0.272, 0.271, 0.269, 0.266, 0.263]]
  
  ; Relative fluxes for He I 5876
  emiss5876_rel_grid = [[4.382, 4.035, 3.655, 3.295, 2.984, 2.752, 2.542], $
                  [4.047, 3.846, 3.548, 3.235, 2.952, 2.715, 2.534], $
                  [3.656, 3.532, 3.364, 3.131, 2.894, 2.686, 2.519]]
  
  ; Relative fluxes for He I 6678
  emiss6678_rel_grid = [[1.274, 1.171, 1.058, 0.950, 0.856, 0.777, 0.714], $
                  [1.181, 1.116, 1.026, 0.932, 0.847, 0.772, 0.712], $
                  [1.068, 1.024, 0.974, 0.902, 0.830, 0.764, 0.708]]
  
  ; Relative fluxes for He I 7065
  emiss7065_rel_grid = [[0.215, 0.230, 0.254, 0.293, 0.356, 0.461, 0.639], $
                  [0.216, 0.230, 0.254, 0.294, 0.356, 0.461, 0.639], $
                  [0.220, 0.234, 0.257, 0.295, 0.357, 0.461, 0.637]]
  
  ; Relative fluxes for He I 7281
  emiss7281_rel_grid = [[0.028, 0.031, 0.035, 0.042, 0.054, 0.075, 0.112], $
                  [0.028, 0.031, 0.035, 0.042, 0.054, 0.075, 0.112], $
                  [0.027, 0.030, 0.035, 0.042, 0.054, 0.075, 0.112]]
  
  ; Relative fluxes for He I 10830
  emiss10830_rel_grid = [[4.223, 4.072, 3.980, 3.990, 4.313, 5.584, 8.557], $
                  [4.033, 3.969, 3.928, 4.702, 13.10, 38.08, 90.87], $
                  [3.884, 3.853, 4.034, 7.482, 20.41, 51.69, 117.1]]
  
  case line of
     1: emiss_grid = emiss3889_rel_grid * emiss4471_grid ; He I 3889
     2: emiss_grid = emiss4026_rel_grid * emiss4471_grid ; He I 4026
     3: emiss_grid = emiss4387_rel_grid * emiss4471_grid ; He I 4387
     4: emiss_grid = emiss4471_grid ; He I 4471
     5: emiss_grid = emiss4922_rel_grid * emiss4471_grid ; He I 4922
     6: emiss_grid = emiss5876_rel_grid * emiss4471_grid ; He I 5876
     7: emiss_grid = emiss6678_rel_grid * emiss4471_grid ; He I 6678
     8: emiss_grid = emiss7065_rel_grid * emiss4471_grid ; He I 7065
     9: emiss_grid = emiss7281_rel_grid * emiss4471_grid ; He I 7281
     10: emiss_grid = emiss10830_rel_grid * emiss4471_grid ; He I 10830
     else: begin
      print, 'line cannnot find'
      exit
     endelse
  endcase
  
  ; Linearly interpolate extinction law in 1/lam
  if temp lt temp_grid[0] or temp gt temp_grid[6] then begin
    print, 'ouside temperature range!'
    return, 0
  endif
  if dens lt dens_grid[0] or dens gt dens_grid[2] then begin
    print, 'ouside density range!'
    return, 0
  endif
  
  ; Linearly interpolate density
  hr_tmp0=lin_interp(emiss_grid[*,0], temp_grid, temp)
  hr_tmp1=lin_interp(emiss_grid[*,1], temp_grid, temp)
  hr_tmp2=lin_interp(emiss_grid[*,2], temp_grid, temp)
  
  ; Linearly interpolate temperature
  hr_tmp_grid=[hr_tmp0, hr_tmp1, hr_tmp2]  
  hr_tmp=lin_interp(hr_tmp_grid, dens_grid, dens)
  
  return, hr_tmp
end
