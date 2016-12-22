function recomb_c_iii, temp, dens, abund
;+
; NAME:
;     recomb_c_iii
; PURPOSE:
;     return the recombination coefficients of C III and N III
;     lines from Pequignot et al. 1991A&A...251..680P
; EXPLANATION:
;
; CALLING SEQUENCE:
;     c_iii_rc=recomb_c_iii(tempi, densi, Abund)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of C III and N III
;          c_iii_rc_structure
;          { Wave:0.0, 
;            a:0.0, b:0.0, c:0.0, d:0.0, Br:0.0, aeff:0.0, 
;            Int:0.0, Obs:0.0, 
;            abundance:0.0, Ion:''}
; REVISION HISTORY:
;     Total and effective radiative recombination coefficients
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P
;     IDL code by A. Danehkar, 10/05/2013
;- 
  common share1, Atomic_Data_Path
  
  c_iii_rc_structure ={Wave:double(0.0), $
                  a:double(0.0), b:double(0.0), c:double(0.0), d:double(0.0), $
                  Br:double(0.0), aeff:double(0.0), $
                  Int:double(0.0), Obs:double(0.0), abundance:double(0.0), Ion:''}

  nlines = 6 
  get_aeff_hb, temp, dens, aeff_hb, em_hb
  
  c_iii_rc=REPLICATE(c_iii_rc_structure, nlines)
  
  orl_Wave=double(0.0)
  orl_a=double(0.0)
  orl_b=double(0.0)
  orl_c=double(0.0)
  orl_d=double(0.0)
  orl_Br=double(0.0)
  orl_Ion=''
    
  ; read C III atomic data
  ion1='R_c_iii'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, NLINES-1 do begin 
    readf,lun1,orl_Ion,orl_Wave, orl_a, orl_b, orl_c, orl_d, orl_Br, Format='(A3, F8.2, F6.3, F7.3, F6.3, F6.3, F6.3)'
    c_iii_rc[i].Ion = orl_Ion
    c_iii_rc[i].Wave = orl_Wave
    c_iii_rc[i].a = orl_a
    c_iii_rc[i].b = orl_b
    c_iii_rc[i].c = orl_c
    c_iii_rc[i].d = orl_d
    c_iii_rc[i].Br = orl_Br
  endfor
  free_lun, lun1
  
  z = 3.0 ; ion level c^3+
  ; equation (1) in 1991A&A...251..680P
  temp4 = 1.0e-4 * temperature /z^2
  for i = 0, 3 do begin 
    ; equation (1) in 1991A&A...251..680P
    c_iii_rc[i].aeff = c_iii_rc[i].Br * 1.0e-13 * z * $
                      (c_iii_rc[i].a*(temp4^c_iii_rc[i].b)) / $
                      (1.0 + (c_iii_rc[i].c * (temp4^c_iii_rc[i].d))) 
    c_iii_rc[i].Int = 100.0 * (c_iii_rc[i].aeff/aeff_hb) * (4861.33/c_iii_rc[i].Wave) * abund 
  endfor
  
  return,c_iii_rc
end
