function recomb_n_iii, temp, dens, abund
;+
; NAME:
;     recomb_n_iii
; PURPOSE:
;     return the recombination coefficients of C III and N III
;     lines from Pequignot et al. 1991A&A...251..680P
; EXPLANATION:
;
; CALLING SEQUENCE:
;     niii_rc=recomb_n_iii(tempi, densi, Abund)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of C III and N III
;          niii_rc_structure
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
  
  niii_rc_structure ={Wave:double(0.0), $
                  a:double(0.0), b:double(0.0), c:double(0.0), d:double(0.0), $
                  Br:double(0.0), aeff:double(0.0), $
                  Int:double(0.0), Obs:double(0.0), abundance:double(0.0), Ion:''}

  nlines = 6 
  getaeffhb, temp, dens, aeff_hb, em_hb
  
  niii_rc=REPLICATE(niii_rc_structure, nlines)
  
  orl_Wave=double(0.0)
  orl_a=double(0.0)
  orl_b=double(0.0)
  orl_c=double(0.0)
  orl_d=double(0.0)
  orl_Br=double(0.0)
  orl_Ion=''
    
  ; read C III & N III data from file
  ion1='Rxiii'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, NLINES-1 do begin 
    readf,lun1,orl_Ion,orl_Wave, orl_a, orl_b, orl_c, orl_d, orl_Br, Format='(A3, F8.2, F6.3, F7.3, F6.3, F6.3, F6.3)'
    niii_rc[i].Ion = orl_Ion
    niii_rc[i].Wave = orl_Wave
    niii_rc[i].a = orl_a
    niii_rc[i].b = orl_b
    niii_rc[i].c = orl_c
    niii_rc[i].d = orl_d
    niii_rc[i].Br = orl_Br
  endfor
  free_lun, lun1
  
  z = 3.0 ; ion level c^3+
  ; equation (1) in 1991A&A...251..680P
  temp4 = 1.0e-4 * temperature /z^2
  for i = 0, 1 do begin 
    ; equation (1) in 1991A&A...251..680P
    niii_rc[i].aeff = niii_rc[i].Br * 1.0e-13 * z * $
                      (niii_rc[i].a*(temp4^niii_rc[i].b)) / $
                      (1.0 + (niii_rc[i].c * (temp4^niii_rc[i].d)))
    niii_rc[i].Int = 100.0 * (niii_rc[i].aeff/aeff_hb) * (4861.33/niii_rc[i].Wave) * abund 
  endfor
  
  return,niii_rc
end
