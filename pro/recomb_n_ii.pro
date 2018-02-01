function recomb_n_ii, nii_rc_data_br, nii_rc_data, h_i_aeff_data, temperature, density, wavelength, iobs 
;+
; NAME:
;     recomb_n_ii
; PURPOSE:
;     return the recombination coefficients of N II lines
;     from Escalante & Victor 1990ApJS...73..513E
; EXPLANATION:
;
; CALLING SEQUENCE:
;     niiRLs=recomb_n_ii(nii_rc_data_br, nii_rc_data, 
;     h_i_aeff_data, temperature, density, wavelength, iobs )
;
; INPUTS:
;     temperature  - electron temperature in K
;     density  - electron density in cm-3
;     abund - abundance coefficient
; RETURN:  recombination coefficients of O II
;    niiRLstructure ={g1:0, g2:0, ION:0, bMult:0, 
;            Wave:0.0, E1:0.0, E2:0.0, Em:0.0, Int:0.0,
;            Br_A:0.0, Br_B:0.0, Br_C:0.0,
;            gf1:0.0, gf2:0.0, Obs:0.0, abundance:0.0,
;            Mult1:'', Term1:'', Term2:''} 
; REVISION HISTORY:
;     Effective recombination coefficients of neutral carbon 
;     and singly ionized nitrogen from
;     Escalante & Victor 1990ApJS...73..513E
;     Revised based on scripts by Yong Zhang added to MOCASSIN, 2003/02
;                       Ercolano et al. 2005MNRAS.362.1038E
;     Converted to IDL code by A. Danehkar, 10/05/2013
;     Integration with AtomNeb, A. Danehkar, 25/04/2017
;- 

  ;  niiRLstructure ={Wave:double(0.0), $ ;REAL*8
  ;              Int:double(0.0),  $  
  ;              Obs:double(0.0), $ 
  ;              abundance:double(0.0), $             
  ;              g1:long(0), $ ;INTEGER 
  ;              g2:long(0), $ ;INTEGER
  ;              Mult1:'', $ ;CHARACTER*7
  ;              Term1:'', $ ;CHARACTER*9
  ;              Term2:'' $ ;CHARACTER*9
  ;              } 
  
  h_Planck = double(6.62606957e-27) ; erg s
  c_Speed = double(2.99792458e10) ; cm/s 
  
  TEh2=double(temperature)
  NEh2=double(density)
  abund=1.0
  nlines = 99
  hbeta_ems= (10.0^gamma4861(h_i_aeff_data,TEh2,NEh2))
  ;hbeta_aeff= (10.0^atomneb_gamma4861(h_i_aeff_data,TEh2,NEh2))*double(4861.33/(h_Planck*c_Speed*1.e8)) 

  ; niiRLs=REPLICATE(niiRLstructure, nlines)
       
  Wave=double(0.0)
  RL_br=double(0.0)
  g1=double(0.0)
  g2=double(0.0)
  temp4 = temperature/10000.0
  loc1=where(abs(nii_rc_data_br.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(nii_rc_data_br[loc1].Wavelength)
    loc1=where(nii_rc_data_br.Wavelength eq  Wavelength_min)
  endif
  Wave=nii_rc_data_br[loc1].Wavelength
  RL_br=nii_rc_data_br[loc1].Br
  
  case 1 of
     ;---------------------------------------
     (loc1 ge 0) and (loc1 le 5): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p E1 3P* - 3D  : 03  row 
        ;i = 0002    ;case A
        i = 0003     ;case B
     end
     ;---------------------------------------
     (loc1 ge 6) and (loc1 le 8): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3S     : 04 row
        ;i = 0004    ;Case: A
        i = 0005     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 9) and (loc1 le 14): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3P     : 05 row
        ;i = 0006    ;Case: A
        i = 0007     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 15): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1P     : 08 row
        i = 0008     ;Case: A
        ;i = 0009    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 16): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1D     : 12 row
        i = 0010     ;Case: A
        ;i = 0011    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 17): begin
        ; atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1S     : 13 row
        i = 0012     ;Case: A
        ;i = 0013    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 18): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1D*     : 15 row
        i = 0014     ;Case: A
        ;i = 0015    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 19): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1P*     : 17 row
        i = 0016     ;Case: A
        ;i = 0017    ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 20) and (loc1 le 25): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3F*     : 19 row
        ;i = 0018    ;Case: A
        i = 0019     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 26) and (loc1 le 32): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3D*     : 20 row
        ;i = 0020    ;Case: A
        i = 0021     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 33) and (loc1 le 38): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3P*     : 21 row
        ;i = 0022    ;Case: A
        i = 0023     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 39) and (loc1 le 44): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3D - 3P*     : 22 row
        ;i = 0024    ;Case: A
        i = 0025     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 45) and (loc1 le 47): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3S - 3P*     : 24 row
        ;i = 0026    ;Case: A
        i = 0027     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 48) and (loc1 le 50): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3S - 3P*     : 26 row
        ;i = 0028    ;Case: A
        i = 0029     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 51) and (loc1 le 56): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3D*     : 28 row
        ;i = 0030    ;Case: A
        i = 0031     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 57) and (loc1 le 62): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3P*     : 29 row
        ;i = 0032    ;Case: A
        i = 0033     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 63) and (loc1 le 68): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3P - 3P*     : 30 row
        ;i = 0034    ;Case: A
        i = 0035     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 69): begin
        ; atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1D - 1F*     : 31 row
        i = 0036     ;Case: A
        ;i = 0037    ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 70) and (loc1 le 75): begin
        ; atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*).4p 3F* - 3D     : 36 row
        ;i = 0038     ;Case: A
        i = 0039     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 76) and (loc1 le 81): begin
        ; atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 3F* - 3G     : 39 row
        ;i = 0040    ;Case: A
        i = 0041     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 82): begin
        ; atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 1F* - 1G     : 58 row
        i = 0042     ;Case: A
        ;i = 0043    ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 83) and (loc1 le 88): begin
        ; atomic transitions: 3d 3D* - 4f 3F 4242     : 48 row
        ;i = 0044    ;Case: A
        i = 0045     ;Case: B
     end
     ;---------------------------------------
     (loc1 ge 89) and (loc1 le 94): begin
        ; atomic transitions: 3d 3P* - 4f 3D 4435     : 55 row
        ;i = 0046   ;Case: A
        i = 0047     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 95): begin
        ; atomic transitions: 3d 1D* - 4f 1F 4176     : 43 (RMT 42) row
        i = 0048     ;Case: A
        ;i = 0049     ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 96): begin
        ; atomic transitions: 3d 1P* - 4f 1D 4677     : 61 (rmt 62) row
        i = 0049     ;Case: A
        ;i = 0051    ;Case: B
     end
     ;---------------------------------------
     (loc1 eq 97): begin
        ; Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
        ; atomic transitions: 3d 3F* - 4f 1G 4026     : 39b row
        ; Case: A 
        a0 = 0.108
        b0 = -0.754
        c0 = 2.587
        d0 = 0.719
        z0 = 2.
        br0 = 0.350   
        aeff = 1.e-13 * z0 * a0  * (temp4/z0^2)^(b0)
        aeff = aeff / (1. + c0 * (temp4/z0^2)^(d0)) * br0
     end
     ;---------------------------------------
     (loc1 eq 98): begin
        ; Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
        ; atomic transitions: 3d 1F* - 4f 3G 4552     : 58a row
        ; Case: A 
        a0 = 0.326
        b0 = -0.754
        c0 = 2.587
        d0 = 0.719
        z0 = 2.
        br0 = 0.074   
        aeff = 1.e-13 * z0 * a0  * (temp4/z0^2)^(b0)
        aeff = aeff / (1. + c0 * (temp4/z0^2)^(d0)) * br0
     end
  else: print, 'wavelength has an illegal value.'
  endcase
  if (loc1 ge 0) and (loc1 le 96) then begin
    aeff = 10.^(nii_rc_data[i-1].a + nii_rc_data[i-1].b * alog10(temp4) + nii_rc_data[i-1].c * alog10(temp4)^2) 
  endif
  
  Ems1 = aeff * (h_Planck*c_Speed*1.e8) / Wave * RL_br
  niiRLs_Int = 100.0 * Ems1 / hbeta_ems * abund 
  
  abund=iobs/niiRLs_Int
  return,abund
end
