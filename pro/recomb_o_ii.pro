function recomb_o_ii, temp, dens, abund
;+
; NAME:
;     recomb_o_ii
; PURPOSE:
;     return the recombination coefficients of O II lines
;     from Storey 1994A&A...282..999S 
;     and Liu et al. 1995MNRAS.272..369L
; EXPLANATION:
;
; CALLING SEQUENCE:
;     oiiRLs=recomb_o_ii(tempi, densi, Abund)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of O II
;    oiiRLstructure ={g1:0, g2:0, ION:0, bMult:0, 
;            Wave:0.0, E1:0.0, E2:0.0, Em:0.0, Int:0.0,
;            Br_A:0.0, Br_B:0.0, Br_C:0.0,
;            gf1:0.0, gf2:0.0, Obs:0.0, abundance:0.0,
;            Mult1:'', Term1:'', Term2:''} 
; REVISION HISTORY:
;     Recombination coefficients for O II lines at nebular 
;     temperatures and densities
;     Storey 1994A&A...282..999S 
;     The rich O II recombination spectrum of the planetary 
;     nebula NGC 7009: new observations and atomic data
;     Liu et al. 1995MNRAS.272..369L
;     Adopted from MIDAS script Roii.prg written by X.W.Liu
;     included in Fortran Program NEAT, 2012MNRAS.422.3516W
;     Revised based MOCASSIN, Ercolano et al. 2005
;     Converted to IDL code by A. Danehkar, 10/05/2013
;- 
  common share1, Atomic_Data_Path
 
  oiiRLstructure ={g1:long(0), $ ;INTEGER 
            g2:long(0), $ ;INTEGER 
            ION:long(0), $ ;INTEGER  
            bMult:0, $
            Wave:double(0.0), $ ;REAL*8
            E1:double(0.0), $ ;REAL*8 
            E2:double(0.0),  $ ;REAL*8    
            Em:double(0.0),  $    
            Int:double(0.0),  $  
            Br_A:double(0.0), $ ;REAL*8 
            Br_B:double(0.0), $ ;REAL*8 
            Br_C:double(0.0), $ ;REAL*8 
            gf1:double(0.0), $ ;REAL*8            
            gf2:double(0.0), $ ;REAL*8
            Obs:double(0.0), $ 
            abundance:double(0.0), $
            Mult1:'', $ ;CHARACTER*7
            Term1:'', $ ;CHARACTER*9
            Term2:'' $ ;CHARACTER*9
            } 

  nlines = 415 
  get_aeff_hb, temp, dens, aeff_hb, em_hb
  
  oiiRLs=REPLICATE(oiiRLstructure, nlines)
  
  orl_ION=long(0)
  orl_Wave=double(0.0)
  orl_Rem=''
  orl_gf1='' ; STR2FLOUT
  orl_q_gf1=''
  orl_gf2='' ; STR2FLOUT
  orl_q_gf2=''
  orl_Mult=''
  orl_E1='' ; STR2FLOUT
  orl_g1=long(0)
  orl_n_g1=''
  orl_Term1=''
  orl_E2='' ; STR2FLOUT
  orl_n_E2=''
  orl_g2=long(0)
  orl_Term2=''
  orl_Br_A=double(0.0)
  orl_Br_B=double(0.0)
  orl_Br_C=double(0.0)
          
  ; read OII atomic data
  ion1='R_o_ii'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, NLINES-1 do begin 
    readf,lun1,$
          orl_ION, $ ;INTEGER  
          orl_Wave, $ ;REAL*8
          ;oiiRLs(i)%Hyb, & CHARACTER*1
          ;&oiiRLs(i)%Rem1, CHARACTER*1
          ;oiiRLs(i)%Rem2, CHARACTER*1
          ;oiiRLs(i)%Rem3, CHARACTER*1
          ;oiiRLs(i)%Rem4, CHARACTER*1
          orl_Rem, $ ;CHARACTER*1
          orl_gf1, $ ;REAL*8   STR2FLOUT
          orl_q_gf1, $ ;CHARACTER*3
          orl_gf2, $ ;REAL*8   STR2FLOUT
          orl_q_gf2,$ ;CHARACTER*3
          orl_Mult, $ ;CHARACTER*7
          orl_E1, $ ;REAL*8    STR2FLOUT
          ;orl_n_E1, $ ;CHARACTER*1
          ;oiiRLs(i)%n_E1GA, CHARACTER*1
          orl_g1, $ ;INTEGER 
          orl_n_g1, $; CHARACTER*1
          orl_Term1, $ ;CHARACTER*9
          orl_E2,  $ ;REAL*8   STR2FLOUT
          orl_n_E2, $ ;CHARACTER*1
          ;oiiRLs(i)%n_E2GA, CHARACTER*1
          orl_g2, $ ;INTEGER 
          ;oiiRLs(i)%n_g2, CHARACTER*1
          orl_Term2, $ ;CHARACTER*9
          orl_Br_A, $ ;REAL*8 
          orl_Br_B, $ ;REAL*8 
          orl_Br_C,  $ ;REAL*8 
         ;Format='(I5, F11, A5, F9, A3, F9, I3, A10, F15, I3, A4, A9, F13, A2, I7, A9, F7, F8, F7)   
          Format='(I5, F11, A5, A9, A3, A9, A3, A10, A15, I3, A4, A9, A13, A2, I7, A9, F7, F8, F7)   
          
    orl_gf1=strtrim(orl_gf1,1)
    orl_gf2=strtrim(orl_gf2,1)
    orl_E1=strtrim(orl_E1,1)
    orl_E2=strtrim(orl_E2,1)
    oiiRLs[i].ION=orl_ION
    oiiRLs[i].Wave=orl_Wave
    ;oiiRLs[i].Rem=orl_Rem
    if orl_gf1 eq '*' then begin
       oiiRLs[i].gf1=double(0.0)
    endif else begin
       oiiRLs[i].gf1=double(orl_gf1)  ; STR2FLOUT
    endelse
    ;oiiRLs[i].q_gf1=orl_q_gf1
    if orl_gf2 eq '*' then begin
       oiiRLs[i].gf2=double(0.0)
    endif else begin
       oiiRLs[i].gf2=double(orl_gf2) ; STR2FLOUT
    endelse
    ;oiiRLs[i].q_gf2=orl_q_gf2
    oiiRLs[i].Mult1=orl_Mult
    if orl_E1 eq '*' then begin
       oiiRLs[i].E1=double(0.0)
    endif else begin
       oiiRLs[i].E1=double(orl_E1) ; STR2FLOUT
    endelse
    ;oiiRLs[i].n_E1=orl_n_E1
    oiiRLs[i].g1=orl_g1
    oiiRLs[i].Term1=orl_Term1
    if orl_E2 eq '*' then begin
       oiiRLs[i].E2=double(0.0)
    endif else begin
       oiiRLs[i].E2=double(orl_E2)  ; STR2FLOUT
    endelse
    ;oiiRLs[i].n_E2=orl_n_E2
    oiiRLs[i].g2=orl_g2
    oiiRLs[i].Term2=orl_Term2
    oiiRLs[i].Br_A=orl_Br_A
    oiiRLs[i].Br_B=orl_Br_B
    oiiRLs[i].Br_C=orl_Br_C  
  endfor
  free_lun, lun1
  
  temp4 = temp/10000.0
  
  ; 4f-3d transitions

  a = 0.232
  b = -0.92009
  c = 0.15526
  d = 0.03442
  an = [0.0, 0.236, 0.232, 0.228, 0.222]

  densi=double(dens)
  log10ne=alog10(densi)
  ;temp
  ;dens
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  
  temp4 = temp/10000.0
  aeff = 1.0e-14 * a * temp4^(b)
  aeff = aeff*(1. + c*(1. - temp4) + d * (1. - temp4)^2)

  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3

  if (temp4 le 0.5) then begin
    a = 0.236
    b = -1.07552
    c = -0.04843
    aeff = 1.0e-14 * a * temp4^(b + c * alog(temp4))
  endif

  for i = 1-1,183-1 do begin 
    oiiRLs[i].Em = aeff * 1.98648E-08 /oiiRLs[i].Wave * $
                   oiiRLs[i].g2 * oiiRLs[i].Br_B
    oiiRLs[i].Int = 100. * oiiRLs[i].Em / Em_hb * abund 
  endfor

  ; 3d-3p ^4F transitions (Case A=B=C for a,b,c,d; Br diff. slightly, adopt Case B)
  a =  0.876
  b =  -0.73465
  c =  0.13689
  d =  0.06220
  an = [0.0, 0.876, 0.876, 0.877, 0.880] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a)
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 0.878
    b = -0.86175
    c = -0.02470
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
  endif
  ;
  for i = 184-1,219-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_B
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3d-3p ^4D, ^4P transitions 
  a =  0.745
  b =  -0.74621
  c =  0.15710
  d =  0.07059
  ;      an = [0.727,0.726,0.725,0.726] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
  an = [0.0, 0.747, 0.745, 0.744, 0.745] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case B
  ;      an = [0.769,0.767,0.766,0.766] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 0.747
    b = -0.89382
    c = -0.02906
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
  endif
;
  for i = 220-1,310-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_B
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3d-3p ^2F transitions 
  a =  0.745
  b =  -0.74621
  c =  0.15710
  d =  0.07059
  an = [0.0, 0.727, 0.726, 0.725, 0.726] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
  ;      an = [0.747,0.745,0.744,0.745] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case B
  ;      an = [0.769,0.767,0.766,0.766] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 0.747
    b = -0.89382
    c = -0.02906
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
  endif
  ;
  for i = 311-1,328-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_A
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3d-3p ^2D transitions 
  a =  0.601
  b =  -0.79533
  c =  0.15314
  d =  0.05322
  an = [0.0, 0.603, 0.601, 0.600, 0.599] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
  ;      an = [0.620,0.618,0.616,0.615] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 0.603
    b = -0.94025
    c = -0.03467
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
  endif
  ;
  for i = 329-1,358-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_A
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3d-3p ^2P transitions 
  a =  0.524
  b =  -0.78448
  c =  0.13681
  d =  0.05608
  an = [0.0, 0.526, 0.524, 0.523, 0.524] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case A
  ;      an = [0.538,0.536,0.535,0.536] ;a for logNe = 2,4,5,6 (LSBC95, Tab.5a) Case C
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 0.526
    b = -0.91758
    c = -0.03120
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4))
  endif
  ;
  for i = 359-1,388-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_A
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3p-3s ^4D - ^4P transitions 
  ;      an = [34.7,34.9,35.1,35.0] ;a for logNe = 2,4,5,6 Case A
  ;      a =  36.2
  ;      b =  -0.749
  ;      c =  0.023
  ;      d =  0.074
  an = [0.0, 36.0, 36.2, 36.4, 36.3] ;a for logNe = 2,4,5,6 Case B
  a =  36.2
  b =  -0.736
  c =  0.033
  d =  0.077
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 36.288
    b = -0.75421
    c = 0.02883
    d = 0.01213
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
  endif
  ;
  for i = 389-1,396-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_B
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3p-3s ^4P - ^4P transitions 
  ;      an = [10.4,10.4,10.5,10.4] ;a for logNe = 2,4,5,6 Case A
  ;      a =  10.4
  ;      b =  -0.721
  ;      c =  0.073
  ;      d =  0.072
  an = [0.0, 14.6, 14.6, 14.7, 14.6] ;a for logNe = 2,4,5,6 Case B
  a =  14.6
  b =  -0.732
  c =  0.081
  d =  0.066
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 14.656
    b = -0.80449
    c = 0.00018
    d = 0.00517
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
  endif
  ;
  for i = 397-1,403-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_B
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3p-3s ^4S - ^4P transitions 
  ;      an = [0.90,0.90,0.90,1.00] ;a for logNe = 2,4,5,6 Case A
  ;      a =  0.90
  ;      b =  -0.485
  ;      c =  -0.047
  ;      d =  0.140
  an = [0.0, 4.80, 4.90, 4.90, 4.90] ;a for logNe = 2,4,5,6 Case B
  a =  4.90
  b =  -0.730
  c =  -0.003
  d =  0.057
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 4.8340
    b = -0.71947
    c = 0.02544
    d = 0.00936
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
  endif
;
  for i = 404-1,406-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_B
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3p-3s ^2D - ^2P transitions 
  an = [0.0, 2.40, 2.40, 2.50, 2.60] ;a for logNe = 2,4,5,6 Case A
  a =  2.40
  b =  -0.550
  c =  -0.051
  d =  0.178
  ;      an = [14.5,14.6,14.5,14.3] ;a for logNe = 2,4,5,6 Case C
  ;      a =  14.6
  ;      b =  -0.736
  ;      c =  0.068
  ;      d =  0.066
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 2.3616
    b = -0.46263
    c = 0.14697
    d = 0.03856
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
  endif
  ;
  for i = 407-1,409-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_A
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor 
  ;
  ; 3p-3s ^2P - ^2P transitions 
  an = [0.0, 1.10, 1.20, 1.20, 1.20] ;a for logNe = 2,4,5,6 Case A
  a =  1.20
  b =  -0.523
  c =  -0.044
  d =  0.173
  ;      an = [1.30,1.40,1.40,1.40] ;a for logNe = 2,4,5,6 Case C
  ;      a =  1.40
  ;      b =  -0.565
  ;      c =  -0.042
  ;      d =  0.158
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 1.1198
    b = -0.44147
    c = 0.13837
    d = 0.03191
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
  endif
  ;
  for i = 410-1,413-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_A
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor
  ;
  ; 3p-3s ^2S - ^2P transitions 
  an = [0.0, 0.40, 0.40, 0.40, 0.40] ;a for logNe = 2,4,5,6 Case A
  a =  0.40
  b =  -0.461
  c =  -0.083
  d =  0.287
  ;      an = [0.50,0.50,0.50,0.60] ;a for logNe = 2,4,5,6 Case C
  ;      a =  0.50
  ;      b =  -0.547
  ;      c =  -0.074
  ;      d =  0.244
  if (log10ne le 2) then begin
    a = an[1]
  endif else begin 
    if (log10ne gt 2 and log10ne le 4) then begin
      a = an[1] + (an[2] - an[1]) / 2. * (log10ne - 2.)
    endif else begin
      if (log10ne gt 4 and log10ne le 5) then begin
        a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
      endif else begin 
        if (log10ne gt 6) and (log10ne le 6) then begin
          a = an[3] + (an[4] - an[3]) * (log10ne - 2.)
        endif else begin
          a = an[4]
        endelse
      endelse
    endelse
  endelse
  ;
  aeff = 1.e-14*a*temp4^(b) 
  aeff = aeff*(1. + c*(1. - temp4) + d*(1. - temp4)^2)
  ;
  ; New eff. recomb. data for 1000 < T < 5000 K, Ne = 100/cm3
  if ( temp4 le 0.5) then begin
    a = 0.3922
    b = -0.35043
    c = 0.26366
    d = 0.06666
    aeff = 1.e-14*a*temp4^(b + c*alog(temp4) + d* alog(temp4)^2)
  endif
  ;
  for i = 414-1,415-1 do begin 
    oiiRLs[i].Em = aeff*1.98648E-08 / oiiRLs[i].Wave* $
                   oiiRLs[i].g2*oiiRLs[i].Br_A
    oiiRLs[i].Int = 100.*oiiRLs[i].Em / Em_hb*abund
  endfor 
  
  return,oiiRLs
end
