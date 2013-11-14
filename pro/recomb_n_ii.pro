function recomb_n_ii, temp, dens, abund
;+
; NAME:
;     recomb_n_ii
; PURPOSE:
;     return the recombination coefficients of N II lines
;     from Escalante & Victor 1990ApJS...73..513E
; EXPLANATION:
;
; CALLING SEQUENCE:
;     niiRLs=recomb_n_ii(tempi, densi, Abund)
;
; INPUTS:
;     temp  - electron temperature in K
;     dens  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of N II
;    niiRLstructure ={g1:0, g2:0, ION:0, bMult:0, 
;            Wave:0.0, E1:0.0, E2:0.0, Em:0.0, Int:0.0,
;            Br_LS:0.0,
;            gf1:0.0, gf2:0.0, Obs:0.0, abundance:0.0,
;            Mult1:'', Term1:'', Term2:''}
; REVISION HISTORY:
;     Effective recombination coefficients of neutral carbon 
;     and singly ionized nitrogen from
;     Escalante & Victor 1990ApJS...73..513E
;     Adopted from MIDAS script Rnii.prg written by X.W.Liu
;     Revised based MOCASSIN, Ercolano et al. 2005
;     and Fortran Program NEAT, 2012MNRAS.422.3516W
;     Converted to IDL code by A. Danehkar, 10/05/2013
;- 
  common share1, Atomic_Data_Path
 
  niiRLstructure ={g1:long(0), $ ;INTEGER 
            g2:long(0), $ ;INTEGER 
            ION:long(0), $ ;INTEGER  
            bMult:0, $
            Wave:double(0.0), $ ;REAL*8
            E1:double(0.0), $ ;REAL*8 
            E2:double(0.0),  $ ;REAL*8    
            Em:double(0.0),  $    
            Int:double(0.0),  $  
            Br_LS:double(0.0), $ ;REAL*8 
            gf1:double(0.0), $ ;REAL*8            
            gf2:double(0.0), $ ;REAL*8
            Obs:double(0.0), $ 
            abundance:double(0.0), $
            Mult1:'', $ ;CHARACTER*7
            Term1:'', $ ;CHARACTER*9
            Term2:'' $ ;CHARACTER*9
            } 

  nlines = 99 
  get_aeff_hb, temp, dens, aeff_hb, em_hb
  
  niiRLs=REPLICATE(niiRLstructure, nlines)
  
  orl_ION=long(0)
  orl_Wave=double(0.0)
  orl_Rem=''
  orl_gf1='' ; STR2FLOUT
  orl_q_gf1=''
  orl_gf2='' ; STR2FLOUT
  orl_q_gf2=''
  orl_Mult=''
  orl_E1='' ; STR2FLOUT
  ;orl_n_E1=''
  orl_g1=long(0)
  orl_n_g1=''
  orl_Term1=''
  orl_E2='' ; STR2FLOUT
  orl_n_E2=''
  orl_g2=long(0)
  orl_Term2=''
  Br_LS=double(0.0)
          
  ; read NII atomic data
  ion1='R_n_ii'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, NLINES-1 do begin 
    readf,lun1,$
          orl_ION, $ ;INTEGER  
          orl_Wave, $ ;REAL*8
          orl_Rem, $ ;CHARACTER*1
          orl_gf1, $ ;REAL*8   STR2FLOUT
          orl_q_gf1, $ ;CHARACTER*3
          orl_gf2, $ ;REAL*8   STR2FLOUT
          orl_q_gf2,$ ;CHARACTER*3
          orl_Mult, $ ;CHARACTER*7
          orl_E1, $ ;REAL*8    STR2FLOUT
          orl_g1, $ ;INTEGER 
          orl_n_g1, $; CHARACTER*1
          orl_Term1, $ ;CHARACTER*9
          orl_E2,  $ ;REAL*8   STR2FLOUT
          orl_n_E2, $ ;CHARACTER*1
          orl_g2, $ ;INTEGER 
          orl_Term2, $ ;CHARACTER*9
          Br_LS, $ ;REAL*8 
          Format='(I5, F11, A5, A9, A3, A9, A3, A10, A15, I3, A4, A9, A13, A2, I7, A9, F7, F8, F6)  
    orl_gf1=strtrim(orl_gf1,1)
    orl_gf2=strtrim(orl_gf2,1)
    orl_E1=strtrim(orl_E1,1)
    orl_E2=strtrim(orl_E2,1)
    niiRLs[i].ION=orl_ION
    niiRLs[i].Wave=orl_Wave
    ;niiRLs[i].Rem=orl_Rem
    if orl_gf1 eq '*' then begin
       niiRLs[i].gf1=double(0.0)
    endif else begin
       niiRLs[i].gf1=double(orl_gf1)  ; STR2FLOUT
    endelse
    ;niiRLs[i].q_gf1=orl_q_gf1
    if orl_gf2 eq '*' then begin
       niiRLs[i].gf2=double(0.0)
    endif else begin
       niiRLs[i].gf2=double(orl_gf2) ; STR2FLOUT
    endelse
    ;niiRLs[i].q_gf2=orl_q_gf2
    niiRLs[i].Mult1=orl_Mult
    if orl_E1 eq '*' then begin
       niiRLs[i].E1=double(0.0)
    endif else begin
       niiRLs[i].E1=double(orl_E1) ; STR2FLOUT
    endelse
    ;niiRLs[i].n_E1=orl_n_E1
    niiRLs[i].g1=orl_g1
    niiRLs[i].Term1=orl_Term1
    if orl_E2 eq '*' then begin
       niiRLs[i].E2=double(0.0)
    endif else begin
       niiRLs[i].E2=double(orl_E2)  ; STR2FLOUT
    endelse
    ;oiiRLs[i].n_E2=orl_n_E2
    niiRLs[i].g2=orl_g2
    niiRLs[i].Term2=orl_Term2
    niiRLs[i].Br_LS=Br_LS
  endfor
  free_lun, lun1
  
  temp4 = temp/10000.0

  ; 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p E1 3P* - 3D  M03  transitions 
  ; i = 0002     ;case A
  i = 0003     ;case B
  a = -12.7289
  b = -0.689816
  c = 0.022005
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2) 
  for ii = 1-1,6-1  do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund 
  endfor 
  ;      
  ; 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3S     M04 transitions
  ; i = 0004     ;case A
  i = 0005     ;case B
  a = -13.8161
  b = -0.778606
  c = -0.028944
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 7-1,9-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3P     M05 transitions
  ; i = 0006     ;case A
  i = 0007     ;case B
  a = -13.0765
  b = -0.734594
  c = -0.0251909
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 10-1,15-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1P     M08 transitions
  i = 0008     ;case A
  ; i = 0009     ;case B
  a = -14.1211
  b = -0.608107
  c = 0.0362301
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 16-1,16-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1D     M12 transitions
  i = 0010     ;case A
  ; i = 0011     ;case B
  a = -13.7473
  b = -0.509595
  c = 0.0255685
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 17-1,17-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1S     M13 transitions
  i = 0012     ;case A
  ; i = 0013     ;case B
  a = -14.3753
  b = -0.515547
  c = 0.0100966
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 18-1,18-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1D*     M15 transitions
  i = 0014     ;case A
  ; i = 0015     ;case B
  a = -14.3932
  b = -0.887946
  c = -0.0525855
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 19-1,19-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1P*     M17 transitions
  i = 0016     ;case A
  ; i = 0017     ;case B
  a = -15.0052
  b = -0.89811
  c = -0.0581789
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 20-1,20-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3F*     M19 transitions
  ; i = 0018     ;case A
  i = 0019     ;case B
  a = -12.6183
  b = -0.840727
  c = -0.0229685
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 21-1,26-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3D*     M20 transitions
  ; i = 0020     ;case A
  i = 0021     ;case B
  a = -13.3184
  b = -0.884034
  c = -0.0512093
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 27-1,33-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3P*     M21 transitions
  ; i = 0022     ;case A
  i = 0023     ;case B
  a = -14.5113
  b = -0.87792
  c = -0.0552785
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 34-1,39-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3D - 3P*     M22 transitions
  ; i = 0024     ;case A
  i = 0025     ;case B
  a = -14.1305
  b = -0.487037
  c = 0.0354135
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 40-1,45-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3S - 3P*     M24 transitions
  ; i = 0026     ;case A
  i = 0027     ;case B
  a = -13.3527
  b = -0.878224
  c = -0.0557112
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 46-1,48-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3S - 3P*     M26 transitions
  ; i = 0028     ;case A
  i = 0029     ;case B
  a = -14.9628
  b = -0.486746
  c = 0.0358261
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 49-1,51-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3D*     M28 transitions
  ; i = 0030     ;case A
  i = 0031     ;case B
  a = -13.0871
  b = -0.883624
  c = -0.0506882
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 52-1,57-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3P*     M29 transitions
  ; i = 0032     ;case A
  i = 0033     ;case B
  a = -13.5581
  b = -0.878488
  c = -0.0557583
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 58-1,63-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3P - 3P*     M30 transitions
  ; i = 0034     ;case A
  i = 0035     ;case B
  a = -14.3521
  b = -0.487527
  c = 0.0355516
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 64-1,69-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1D - 1F*     M31 transitions
  i = 0036     ;case A
  ; i = 0037     ;case B
  a = -15.0026
  b = -0.923093
  c = -0.0588371
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 70-1,70-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 2s2.2p.(2P*).3d - 2s2.2p.(2P*).4p 3F* - 3D     M36 transitions
  ; i = 0038     ;case A
  i = 0039     ;case B
  a = -13.8636
  b = -0.569144
  c = 0.0068655
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 71-1,76-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;       
  ;       2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 3F* - 3G M39 transitions
  ; i = 0040     ;case A
  i = 0041     ;case B
  a = -13.035
  b = -1.12035
  c = -0.10642
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 77-1,82-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;   
  ; 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 1F* - 1G M58 transitions
  i = 0042     ;case A
  ; i = 0043     ;case B
  a = -13.5484
  b = -1.11909
  c = -0.105123
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 83-1,83-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 3d 3D* - 4f 3F 4242 M48 transitions
  ; i = 0044     ;case A
  i = 0045     ;case B
  a = -13.2548
  b = -1.12902
  c = -0.110368
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 84-1,89-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 3d 3P* - 4f 3D 4435 M55 transitions
  ; i = 0046     ;case A
  i = 0047     ;case B
  a = -13.5656
  b = -1.11989
  c = -0.105818
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 90-1,95-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 3d 1D* - 4f 1F 4176 M43 (RMT M42) transitions
  i = 0048     ;case A
  ; i = 0049     ;case B
  a = -13.7426
  b = -1.13351
  c = -0.111146
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 96-1,96-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 3d 1P* - 4f 1D 4677 M61 (RMT M62) transitions
  i = 0049     ;case A
  ; i = 0051     ;case B
  a = -13.7373
  b = -1.12695
  c = -0.108158
  ;
  aeff = 10.^(a + b * alog10(temp4) + c * alog10(temp4)^2)
  for ii = 97-1,97-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;
  ; 3d 3F* - 4f 1G 4026 M39b transitions
  ; case A (PPB):
  a = 0.108
  b = -0.754
  c = 2.587
  d = 0.719
  z = 2.
  Br_term = 0.350
  ;      
  aeff = 1.e-13 * z * a  * (temp4/z^2)^(b)
  aeff = aeff / (1. + c * (temp4/z^2)^(d)) * Br_term
  for ii = 98-1,98-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  ;      
  ; 3d 1F* - 4f 3G 4552 M58a transitions
  ; case A (PPB):
  a = 0.326
  b = -0.754
  c = 2.587
  d = 0.719
  z = 2.
  Br_term = 0.074
  ;      
  aeff = 1.e-13 * z * a  * (temp4/z^2)^(b)
  aeff = aeff / (1. + c * (temp4/z^2)^(d)) * Br_term
  for ii = 99-1,99-1 do begin 
    niiRLs[ii].Em = aeff * 1.98648E-08 / niiRLs[ii].Wave * niiRLs[ii].Br_LS
    niiRLs[ii].Int = 100.0 * niiRLs[ii].Em / Em_Hb * abund
  endfor
  
  return,niiRLs
end
