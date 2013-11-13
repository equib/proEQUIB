function recomb_ne_ii, temp, dens, abund
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
;     scripts added by Yong Zhang to MOCASSIN, 2003/02
;     Converted to IDL code by A. Danehkar, 10/05/2013
;- 
  common share1, Atomic_Data_Path
  
  ne_ii_orl_structure ={Wave:double(0.0), $
                  a:double(0.0), b:double(0.0), c:double(0.0), d:double(0.0), $
                  f:double(0.0), Br:double(0.0), aeff:double(0.0), $
                  Int:double(0.0), Obs:double(0.0), abundance:double(0.0)}

  nlines = 38 
  get_aeff_hb, temp, dens, aeff_hb, em_hb
  
  ne_ii_orl=REPLICATE(ne_ii_orl_structure, nlines)
  
  orl_Wave=double(0.0)
  orl_a=double(0.0)
  orl_b=double(0.0)
  orl_c=double(0.0)
  orl_d=double(0.0)
  orl_f=double(0.0)
  orl_Br=double(0.0)
  
  ; read Ne II data from file
  ion1='R_ne_ii'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, NLINES-1 do begin 
    readf,lun1,orl_Wave, orl_a, orl_b, orl_c, orl_d, orl_f, orl_Br
    ne_ii_orl[i].Wave = orl_Wave
    ne_ii_orl[i].a = orl_a
    ne_ii_orl[i].b = orl_b
    ne_ii_orl[i].c = orl_c
    ne_ii_orl[i].d = orl_d
    ne_ii_orl[i].f = orl_f
    ne_ii_orl[i].Br = orl_Br
  endfor
  free_lun, lun1
  
  temp4 = temp/10000.0
  for i = 0, NLINES-1 do begin 
    ne_ii_orl[i].aeff = ne_ii_orl[i].Br * 1.0e-14 * (ne_ii_orl[i].a*(temp4^ne_ii_orl[i].f)) * $
                    (1.0 + (ne_ii_orl[i].b*(1.0-temp4)) + (ne_ii_orl[i].c * ((1.0-temp4)^2) ) $
                     + (ne_ii_orl[i].d * ((1.0-temp4)^3) ) )
    ne_ii_orl[i].Int = 100.0 * (ne_ii_orl[i].aeff/aeff_hb) * (4861.33/ne_ii_orl[i].Wave) * abund 
  endfor
  
  return,ne_ii_orl
end
