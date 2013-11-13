function recomb_c_ii, temp, dens, abund
;+
; NAME:
;     recomb_c_ii
; PURPOSE:
;     return the recombination coefficients of C II
;     from Davey et al. (2000) 2000A&AS..142...85D
; EXPLANATION:
;
; CALLING SEQUENCE:
;     c_ii_orl=recomb_c_ii(tempi, densi, Abund)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of C II
;          c_ii_orl_structure
;          { Wave:0.0, 
;            a:0.0, b:0.0, c:0.0, d:0.0, f:0.0, aeff:0.0, 
;            Int:0.0, Obs:0.0, 
;            Abundance:0.0}
; REVISION HISTORY:
;     Recombination coefficients for C Ii lines from
;     Davey et al. 2000A&AS..142...85D
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
;     Originally added by Yong Zhang to MOCASSIN, 2003/02
;     Converted to IDL code by A. Danehkar, 10/05/2013
;- 
  common share1, Atomic_Data_Path
  
  c_ii_orl_structure ={Wave:double(0.0), $
                  a:double(0.0), b:double(0.0), c:double(0.0), d:double(0.0), $
                  f:double(0.0), aeff:double(0.0), $
                  Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  nlines = 57 
  get_aeff_hb, temp, dens, aeff_hb, em_hb
  
  c_ii_orl=REPLICATE(c_ii_orl_structure, nlines)
  
  orl_Wave=double(0.0)
  orl_a=double(0.0)
  orl_b=double(0.0)
  orl_c=double(0.0)
  orl_d=double(0.0)
  orl_f=double(0.0)
  
  ; read CII data from file
  ion1='R_c_ii'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, NLINES-1 do begin 
    readf,lun1,orl_Wave, orl_a, orl_b, orl_c, orl_d, orl_f
    c_ii_orl[i].Wave = orl_Wave
    c_ii_orl[i].a = orl_a
    c_ii_orl[i].b = orl_b
    c_ii_orl[i].c = orl_c
    c_ii_orl[i].d = orl_d
    c_ii_orl[i].f = orl_f
  endfor
  free_lun, lun1
  
  temp4 = temp/10000.0
  for i = 0, NLINES-1 do begin 
    c_ii_orl[i].aeff = 1.0e-14 * (c_ii_orl[i].a*(temp4^c_ii_orl[i].f)) * $
                    (1.0 + (c_ii_orl[i].b*(1.0-temp4)) + (c_ii_orl[i].c * $
                    ((1.0-temp4)^2) ) + (c_ii_orl[i].d * ((1.0-temp4)^3) ) ) 
    c_ii_orl[i].Int = 100.0 * (c_ii_orl[i].aeff/aeff_hb) * (4861.33/c_ii_orl[i].Wave) * abund 
  endfor
  
  return,c_ii_orl
end
