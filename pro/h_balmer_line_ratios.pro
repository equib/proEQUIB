function h_balmer_line_ratios, temp, dens, line
;+
; NAME:
;     h_balmer_line_ratios
; PURPOSE:
;     return Hydrogen Case B Balmer line ratio for
;     given electron temperature and density
;     Table 4.4 in D. E. Osterbrock & G. J. Ferland, 
;     Astrophysics of Gaseius Nebulae and
;     Active Galactic Nuclei, 2nd Ed., 2006
; EXPLANATION:
;
; CALLING SEQUENCE:
;     H_balmer_theory = h_balmer_line_ratios(temp, dens, line)
;
; INPUTS:
;     temp -     electron temperature in K
;     dens -     electron density in cm-3
;     line -     Balmer line number
;        1 : Ha/Hb
;        2 : Hg/Hb
;        3 : Hd/Hb
;        4 : H10/Hb
;        5 : H15/Hb
;        6 : H20/Hb
; RETURN:  Hydrogen Case B Balmer line ratio
;
; REVISION HISTORY:
;     from Table 4.4 in D. E. Osterbrock & 
;     G. J. Ferland, Astrophysics of Gaseius Nebulae 
;     and Active Galactic Nuclei, 2nd Ed., 2006
;     IDL code by A. Danehkar, 31/08/2012
;-  
  dens_grid=[1.0e2, 1.0e4, 1.0e6]  
  temp_grid=[5000.0, 10000.0, 20000.0]
  case line of
     ; D. E. Osterbrock & G. J. Ferland, 
     ; Astrophysics of Gaseius Nebulae and Active Galactic Nuclei, 2nd Ed., 2006
     ; Table 4.4, Hydrogen Case B Balmer line ratio
     1: hr_grid = [[3.041, 3.001, 2.918], [2.863, 2.847, 2.806], [2.747, 2.739, 2.725]] ; Ha/Hb
     2: hr_grid = [[0.458, 0.460, 0.465], [0.468, 0.469, 0.471], [0.475, 0.476, 0.476]] ; Hg/Hb
     3: hr_grid = [[0.251, 0.253, 0.258], [0.259, 0.260, 0.262], [0.264, 0.264, 0.266]] ; Hd/Hb
     4: hr_grid = [[0.0515, 0.0520, 0.0616], [0.0530, 0.0533, 0.0591], [0.0540, 0.0541, 0.0575]] ; H10/Hb
     5: hr_grid = [[0.01534, 0.01628, 0.02602], [0.01561, 0.01620, 0.02147], [0.01576, 0.01612, 0.01834]] ; H15/Hb
     6: hr_grid = [[0.00657, 0.00819, 0.01394], [0.00662, 0.00755, 0.01058], [0.00664, 0.00717, 0.00832]] ; H20/Hb
     else: begin
      print, 'line cannnot find'
      exit
     endelse
  endcase
  
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
