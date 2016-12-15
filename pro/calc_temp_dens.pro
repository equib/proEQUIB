function calc_temp_dens, ion, levu, levl, inratio, diagtype, fixedq
;+
; NAME:
;     calc_temp_dens
; PURPOSE:
;     determine electron density or temperature from given 
;     flux intensity ratio for specified ion with upper level(s)
;     lower level(s) by solving atomic level populations and 
;     line emissivities in statistical equilibrium 
;     for a fixed electron density or temperature.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     path='proEQUIB/atomic-data/'
;     set_atomic_data_path, path
;
;     ion='sii'
;     levu='1,2,1,3/'
;     levl='1,5/'
;     diagtype='T'
;     dens = double(2550)
;     niiTratio=double(10.753)
;     temp=calc_temp_dens(ion, levu, levl, niiTratio, diagtype, dens) 
;     print, temp
;
;     ion='sii'
;     levu='1,2/'
;     levl='1,3/'
;     diagtype='D'
;     temp=double(7000.0)
;     siiNratio=double(1.506)
;     dens=calc_temp_dens(ion, levu, levl, siiNratio, diagtype, temp)
;     print, dens
;
; INPUTS:
;     ion -       ion name e.g. 'sii', 'nii'
;     levu -      upper level(s) e.g '1,2/', '1,2,1,3/'
;     levl -      lower level(s) e.g '1,2/', '1,2,1,3/'
;     inratio -   flux intensity ratio
;     diagtype -  diagnostics type 
;                 'd' or 'D' for electron density
;                 't' or 'T' for electron temperature
;     fixedq -    fixed quantity 
;                 electron density when diagtype ='t' or 'T'
;                 electron temperature when diagtype ='d' or 'D'
; RETURN:  density or temperature
;                 electron density when diagtype ='d' or 'D'
;                 electron temperature when diagtype ='t' or 'T'
; REVISION HISTORY:
;     Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
;     Replaced str2int with strnumber, A. Danehkar, 20/10/2016
;     Replaced CFY, SPLMAT, and CFD with
;          IDL function INTERPOL( /SPLINE), A. Danehkar, 20/10/2016
;     Replaced LUSLV with IDL LAPACK function 
;                       LA_LINEAR_EQUATION, A. Danehkar, 20/10/2016
;     Replaced LA_LINEAR_EQUATION (not work in GDL)
;           with IDL function LUDC & LUSOL, A. Danehkar, 15/11/2016
;     Replaced INTERPOL (not accurate) with 
;                    SPL_INIT & SPL_INTERP, A. Danehkar, 19/11/2016
;     Made a new function calc_populations() for solving atomic 
;       level populations and separated it from
;       calc_abundance(), calc_temp_dens(), A. Danehkar, 20/11/2016
;     Integration with AtomNeb, now uses atomic data input elj_data,
;                      omij_data, aij_data, A. Danehkar, 10/03/2017
;     Cleaning the function, and remove unused varibales
;                           calc_temp_dens(), A. Danehkar, 12/06/2017
; 
; FORTRAN EQUIB HISTORY (F77/F90):
; 1981-05-03 I.D.Howarth  Version 1
; 1981-05-05 I.D.Howarth  Minibug fixed!
; 1981-05-07 I.D.Howarth  Now takes collision rates or strengths
; 1981-08-03 S.Adams      Interpolates collision strengths
; 1981-08-07 S.Adams      Input method changed
; 1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
;                         filenames given to SA's data files.
; 1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
; 1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
;                         modified such that matrix sizes (i.e. maximum
;                         of Te and maximum no of levels) can now be cha
;                         by modifying the parameters NDIM1, NDIM2 and N
;                         in the Main program. EASY!
;                         Now takes collision rates as well.
;                         All variables are declared explicitly
;                         Generate two extra files (ionpop.lis and ionra
;                         of plain stream format for plotting
; 1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
;                         Fixed readin bug for IBIG=2 case.
;                         Now reads reformatted upsilons (easier to see
;                         and the 0 0 0 data end is excluded for these c
;                         The A values have a different format for IBIG=
; 2006       B.Ercolano   Converted to F90
; 2009-04    R.Wesson     Misc updates and improvements, inputs from cmd line, 
;                         written purely to do diagnostics.
;- 
  common share1, Atomic_Data_Path
  
  NDIM1= long(35) 
  NDIM2= long(150)
  ; NDIM1T3 should be at least 3*NDIM1
  NDIM1T3= long(105)
  ; Maximum no. of Ne increments
  MAXND= long(100)
  
  GX= long(0)
  G=lonarr(NDIM2+1)
  ID=lonarr(2+1)
  JD=lonarr(2+1)
  ITRANA=lonarr(2+1,NDIM2+1)
  ITRANB=lonarr(2+1,NDIM2+1)
  ITRANC=lonarr(2+1,NDIM2+1)
  iteration= long(0)
     
  N=dblarr(NDIM2+1)
  TNIJ=dblarr(NDIM2+1,NDIM2+1)
  FINTIJ=dblarr(NDIM2+1,NDIM2+1)
  WAVA=dblarr(NDIM2+1)
  WAVB=dblarr(NDIM2+1)
  WAVC=dblarr(NDIM2+1)
  CS=dblarr(NDIM2+1,NDIM2+1)
  QEFF=dblarr(NDIM2+1,NDIM2+1)   
  QQ=dblarr(NDIM1+1)   
  QOM=dblarr(NDIM1+1,NDIM2+1,NDIM2+1)   
  A=dblarr(NDIM2+1,NDIM2+1)   
  E=dblarr(NDIM2+1)   
  T=dblarr(NDIM1+1)   
  ROOTT=dblarr(NDIM1+1) 
  X=dblarr(NDIM2+1,NDIM2+1)          
  Y=dblarr(NDIM2+1)   
  X2=dblarr(NDIM2+1,NDIM2+1)  
  XKEEP=dblarr(NDIM2+1,NDIM2+1)  
  Y2=dblarr(NDIM2+1)  
  YKEEP=dblarr(NDIM2+1)  
  HMH=dblarr(NDIM1+1,NDIM1+1)          
  D=dblarr(NDIM1+1)  
  check_value=dblarr(3+1)
     
  LABEL1=STRARR(NDIM2+1)
  
  I= long(0) 
  I1= long(0) 
  I2= long(0) 
  J= long(0) 
  K= long(0) 
  L= long(0) 
  KK= long(0) 
  LL= long(0) 
  JT= long(0) 
  JJD= long(0) 
  IONL= long(0) 
  NLINES= long(0) 
  NLEV= long(0) 
  NTEMP= long(0) 
  IBIG= long(0) 
  IRATS= long(0) 
  NTRA= long(0) 
  ITEMP= long(0) 
  IN= long(0) 
  NLEV1= long(0) 
  KP1= long(0) 
  INT= long(0) 
  IND= long(0) 
  IOPT= long(0) 
  IT= long(0) 
  IM1= long(0) 
  JM1= long(0) 
  IP1= long(0) 
  IAPR= long(0) 
  IBPR= long(0) 
  ICPR= long(0) 
  IKT= long(0) 
  IA= long(0) 
  IB= long(0) 
  IC= long(0) 
  IA1= long(0) 
  IA2= long(0) 
  IB1= long(0) 
  IB2= long(0) 
  IC1= long(0) 
  IC2= long(0) 
     
  TEMPI=double(0) 
  TINC=double(0)
  DENSI=double(0) 
  DINC=double(0)
  DENS=double(0)
  DLOGD=double(0)
  TEMP=double(0)
  TLOGT=double(0)
  TEMP2=double(0)
  DD=double(0)
  DELTEK=double(0)
  EXPE=double(0)
  VALUE=double(0)
  SUMN=double(0)
  TTT=double(0)
  TTP=double(0)
  AHB=double(0)
  EJI=double(0)
  WAV=double(0)
  RLINT=double(0)
  FINT=double(0)
  SUMA=double(0)
  SUMB=double(0)
  SUMC=double(0)
  QX=double(0)
  AX=double(0)
  EX=double(0)
  FRAT=double(0)
  DEE=double(0)
  LTEXT = '';
  
  result1=double(0)
     
  I= long(0)
  J= long(0)
  K= long(0)
  IP1= long(0)      
  ;A=dblarr(NR,NR)
  ;FACT=double(0)
  
  itrana[*,*]=0
  itranb[*,*]=0
  itranc[*,*]=0
  
  levu_str=strsplit(levu, ',', ESCAPE='/', /EXTRACT)
  levl_str=strsplit(levl, ',', ESCAPE='/', /EXTRACT)
  
  levu_num=size(levu_str, /N_ELEMENTS)
  levl_num=size(levl_str, /N_ELEMENTS)
  
  levu_i=0
  for i=1, 150 do begin 
    ITRANA[1,i]=equib_str2int(levu_str[levu_i])
    ITRANA[2,i]=equib_str2int(levu_str[levu_i+1])
    levu_i = levu_i + 2
    if levu_i ge levu_num then break
  endfor

  levl_i=0
  for i=1, 150 do begin 
    ITRANB[1,i]=equib_str2int(levl_str[levl_i])
    ITRANB[2,i]=equib_str2int(levl_str[levl_i+1])
    levl_i = levl_i + 2
    if levl_i ge levl_num then break;
  endfor
  
  ;READ(levu,*) ((ITRANA(LL,KK),LL=1,2),KK=1,150)
  ;READ(levl,*) ((ITRANB(LL,KK),LL=1,2),KK=1,150)
      
  ion1=strtrim(ion,1)
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  readf,lun1,NLINES
  for i=1, NLINES do begin 
    readf,lun1,LTEXT
  endfor
  ; Read no. of levels (max=NDIM2) NLEV,
  readf,lun1,NLEV, NTEMP
  ; no. of Te (max=NDIM1) NTEMP and the
  for i=1, NLEV do begin 
    ; input format (cf Readme)
    readf,lun1,LTEXT
    LABEL1[I]=LTEXT
  endfor
  ; be
  ibig=0 
  ; Read in Te's where coll. strengths are tabulated
  for i=1, NTEMP do begin 
    ddtemp=double(0.0)
    readf,lun1,ddtemp
    T[i] = ddtemp
    T[i] = alog10(T[i])
    ROOTT[i] = sqrt(T[i])
  endfor 
  ; If IRATS=0, what tabulated are collision strengths
  readf,lun1,IRATS
  ; Else Coll. rates = tabulated values * 10 ** IRATS
 
  if (IBIG eq 0) then begin
    QX = 1.0
    while (QX ne 0.D0) do begin 
      lontemp1=long(0)
      lontemp2=long(0)
      ddtemp=double(0)
      readf,lun1,lontemp1, lontemp2, ddtemp
      ID[2]=lontemp1
      JD[2]=lontemp2
      QX = ddtemp
      if QX eq 0 then break
      if (ID[2] eq 0) then begin
        ID[2] = ID[1]
        K = K + 1
      endif else begin
        ID[1]= ID[2]
        K = 1
      endelse
      if (JD[2] eq 0) then begin
        JD[2] = JD[1]
      endif else begin
        JD[1] = JD[2]
      endelse
      I = ID[2]
      J = JD[2]
      QOM[K,I,J] = QX
    endwhile
  endif
  
  if(IBIG eq 1) or (IBIG eq 2) then begin
    readf,lun1,NTRA
    for IN = 1, NTRA do begin
      readf,lun1,I,J,QOM[1:NTEMP,I,J]
      ;READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
    endfor
  endif
  ; Read transition probabilities
  NLEV1 = NLEV - 1
  if (IBIG eq 1) then begin
    readf,lun1,I,J,A[J,I];,L=K+1,NLEV),K=1,NLEV1
    ;READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
  endif else begin
    for K = 1, NLEV1 do begin
      KP1 = K + 1
      for L = KP1, NLEV do begin
        readf,lun1, I, J, AX
        A[J,I] = AX
      endfor
    endfor
  endelse
  ; Read statistical weights, energy levels (cm-1)
  for J = 1, NLEV do begin
    readf,lun1,  I, GX, EX
    G[I] = GX
    E[I] = EX
  endfor
      
  free_lun, lun1
  
  ITRANC[*,*] = 0
  
  ; start of iterations
  for iteration = 1, 9 do begin
    if (diagnostic_type eq 't') or (diagnostic_type eq 'T') then begin
      if (iteration eq 1) then begin
        TEMPI=5000.0
      endif else begin 
        TEMPI= check_value[1]
      endelse
      INT=4
      TINC=(15000.0)/((INT-1)^(iteration))
      densi=fixedq
      dinc=0
      ind=1
      ; ALLOCATE(RESULTS(3,INT))
      RESULTS=dblarr(3+1,INT+1)
    endif else begin
      if (iteration eq 1) then begin
        densi=0.0
      endif else begin
        densi=check_value[2]
      endelse
      IND=4
      DINC=(100000.0)/((IND-1)^(iteration))

      TempI=fixedq
      TINC=0
      INT=1
        
      ;allocate(results(3,IND))
      RESULTS=dblarr(3+1,IND+1)
    endelse
    if (densi le 0) then densi=1
    if (tempi lt 5000) then tempi=5000 ; add
    ; Start of Te iteration
    for JT = 1, INT do begin
      TEMP=TEMPI+(JT-1)*TINC 
      ; Start of Ne iteration=
      for JJD = 1, IND  do begin
        DENS=DENSI+(JJD-1)*DINC
        ; IF(DENSI.LT.30.D0) THEN
        ; DENS=10.D0**DENS
        ; ENDIF
        if (TEMP le 0.D0) or (DENS le 0.D0) then begin
            print,'Temp = ', TEMP, ', Dens = ', DENS
            return, 0
        endif
        DLOGD = alog10(DENS)
        TLOGT = alog10(TEMP)
        TEMP2= sqrt(TEMP)
        ; Form matrices
        X[*,*]=double(0.0)
        CS[*,*]=double(0.0)
        QEFF[*,*]=double(0.0)
        TNIJ[*,*]=double(0.0)
        Y[*]=double(0.0)

        IOPT=0
        if (NTEMP eq 1) then begin
          print, 'Coll. strengths available for 1 Te only - assuming const'
        endif else begin
          if (NTEMP eq 2) then begin
            print, 'Coll. strengths available for 2 Te only - linear interp'
          endif else begin
            equib_SPLMAT, T, NTEMP, IOPT, NDIM1, NDIM1T3, HMH
            equib_CFD, TLOGT, T, NTEMP, NDIM1, HMH, D
          endelse
        endelse
        for I = 2, NLEV do begin
          for J = I, NLEV do begin
            ;Negative!
            DELTEK = (E[I-1]-E[J])*1.4388463D0
            EXPE = exp(DELTEK/TEMP)
            for IT = 1, NTEMP do begin
            
              if (IRATS eq 0.D+00) then begin
                QQ[IT] = QOM[IT,I-1,J]
              endif else begin
                ;Take out the exp. depend.
                QQ[IT] = QOM[IT,I-1,J] / EXPE
                ; before interpolation
              endelse
              
            endfor
            
            if (NTEMP eq 1) then begin
              DD = QQ[1]
            endif else begin 
            
              if (NTEMP eq 2) then begin
                DD = QQ[1] +  (QQ[2] - QQ[1])/(T[2] - T[1]) * (TLOGT - T[1])
              endif else begin
                equib_CFY, TLOGT, DD, T, QQ, NTEMP, NDIM1, HMH, D
              endelse
            endelse
            if (IRATS eq 0.D+00) then begin
              CS[I-1,J] = DD
            endif else begin
              CS[I-1,J] = DD * EXPE
            endelse
              
            if (IRATS eq 0.D+00) then begin
              QEFF[I-1,J] = 8.63D-06*CS[I-1,J] * EXPE / (G[I-1]*TEMP2)
              QEFF[J,I-1] = 8.63D-06 * CS[I-1,J] / (G[J]*TEMP2)
            endif else begin
              QEFF[I-1,J] = CS[I-1,J] * 10.^IRATS
              ; Be careful
              QEFF[J,I-1] = G[I-1] * QEFF[I-1,J] / (EXPE * G[J])
              ; G integer!
            endelse
          endfor
        endfor
        for I = 2, NLEV do begin
          for J = 1, NLEV do begin
            if (J ne I) then begin
              X[I,J] = X[I,J] + DENS * QEFF[J,I]
              X[I,I] = X[I,I] - DENS * QEFF[I,J]
              if (J gt I) then begin
                X[I,J] = X[I,J] + A[J,I]
              endif else begin 
                X[I,I] = X[I,I] - A[I,J]
              endelse
            endif
          endfor
        endfor
        for I = 2, NLEV do begin
          IM1 = I - 1
          VALUE = 0.D0 - X[I,1]
          Y[IM1] = VALUE
          Y2[IM1] = VALUE
          YKEEP[IM1] = VALUE
          for J = 2, NLEV do begin
            JM1 = J - 1
            VALUE = X[I,J]
            X[IM1,JM1] = VALUE
            X2[IM1,JM1] = VALUE
            XKEEP[IM1,JM1] = VALUE
          endfor
        endfor
        ; Solve matrices for populations
        equib_LUSLV, X, Y, NLEV1, NDIM2
        for I = NLEV, 2, -1 do begin
          N[I] = Y[I-1]
        endfor
        SUMN = 1.D0
        for I = 2, NLEV do begin
          SUMN = SUMN + N[I]
        endfor
        for I = 2, NLEV do begin
          N[I] = N[I] / SUMN
        endfor
        N[1] = 1.D0 / SUMN
        ; Output data 
        TTT=TEMP*1.0D-4
        TTP=TTT^(-0.87D0)
        ; Eff. recombination coef. of Hb
        AHB=3.036D-14*TTP
        for I = 1, NLEV1 do begin
          IP1 = I + 1
          for J = IP1, NLEV do begin
            if (A[J,I] ne 0.D0) then begin
              EJI = E[J] - E[I]
              WAV = 1.D8 / EJI
              RLINT = A[J,I] * EJI
              RLINT = RLINT *N[J]
              TNIJ[I,J]=RLINT
              FINT=N[J]*A[J,I]*4861.D0/(DENS*AHB*WAV)
              FINTIJ[I,J]=FINT 
            endif
          endfor
        endfor
        ; Search ITRANA, ITRANB & ITRANC for transitions & sum up
        SUMA=0.D0
        SUMB=0.D0
        SUMC=0.D0
        IAPR=0
        IBPR=0
        ICPR=0
        for IKT = 1, NDIM2 do begin
          IA1=ITRANA[1,IKT]
          IA2=ITRANA[2,IKT]
          if(IA1 ne 0) and (IA2 ne 0) then begin
            SUMA=SUMA+TNIJ[IA1,IA2]
            IAPR=IAPR+1
          endif
          IB1=ITRANB[1,IKT]
          IB2=ITRANB[2,IKT]
          if(IB1 ne 0) and (IB2 ne 0) then begin
            IBPR=IBPR+1
            SUMB=SUMB+TNIJ[IB1,IB2]
          endif

          IC1=ITRANC[1,IKT]
          IC2=ITRANC[2,IKT]
          IF(IC1 ne 0) and (IC2 ne 0) then begin
            ICPR=ICPR+1
            SUMC=SUMC+FINTIJ[IC1,IC2]
          endif
        endfor
        FRAT=SUMA/SUMB
        ; SUMC = 1./SUMC
        ; TDRAT(1,JJD)=DENS  !are these lines necessary,
        ; TDRAT(2,JJD)=FRAT  !TDRAT is now never used again?
        ; write(6,*),jd,suma,sumb,sumc,dens,frat
        ; WRITE(7,1017) TEMP, DENS, SUMC
        ; WRITE(8,1017) TEMP, DENS, FRAT, FRAT-inratio
        if (diagnostic_type eq 't') or (diagnostic_type eq 'T') then begin
          RESULTS[1, JT] = TEMP
          RESULTS[2, JT] = DENS
          RESultS[3, JT] = FRAT-inratio
        endif else begin
          RESULTS[1, JJD] = TEMP
          RESULTS[2, JJD] = DENS
          RESultS[3, JJD] = FRAT-inratio
        endelse ;End of the Ne iteration
      endfor

      for IA = 1, IAPR do begin
        I1=ITRANA[1,IA]
        I2=ITRANA[2,IA]
        DEE=E[I2]-E[I1]
        WAVA[IA]=1.D8/DEE
      endfor
      for IB = 1, IBPR do begin
        I1=ITRANB[1,IB]
        I2=ITRANB[2,IB]
        DEE=E[I2]-E[I1]
        WAVB[IB]=1.D8/DEE
      endfor
      for IC = 1, ICPR do begin
        I1=ITRANC[1,IC]
        I2=ITRANC[2,IC]
        DEE=E[I2]-E[I1]
        WAVC[IC]=1.D8/DEE
      endfor
    ; End of the temperature iteration
    endfor
    if (diagnostic_type eq 'D') or (diagnostic_type eq 'd') then begin
      INT = ind
    endif

    ; iteration and detect the sign change.
    for I=2,INT do begin
      check=0
      if (equib_sign(results[3,I],results[3,1]) ne results[3,I]) then begin 
         ;if this condition, the values have a different sign
        check_value[*] = results[*,I-1] ; the value before the sign change returned
        check=1
        break
      endif
    endfor
    
    if(check eq 0) and (iteration lt 9) then begin ;check if no change of sign,
                             ;and checks if it should be upper or lower limit
      if(abs(results[3,1])) lt (abs(results[3,INT])) then begin
          check_value[*]=results[*,1]
      endif else begin 
                if(abs(results[3,INT]) lt abs(results[3,1])) then begin
                check_value[*]=results[*,INT-1]
            endif else begin
                print,'check_value failed'
                return, 0
           endelse
      endelse
    endif else begin 
      if (check eq 0) and (iteration eq 9) then begin ; check if there is any change of sign,
                             ; this kicks in then, and checks if it should be upper or lower limit
      if(abs(results[3,1]) lt abs(results[3,INT])) then begin
         check_value[*]=results[*,1]
      endif else begin 
                if (abs(results[3,INT]) lt abs(results[3,1])) then begin
                check_value[*]=results[*,INT]
            endif else begin
                print,'check_value failed'
                return, 0
            endelse
          endelse
      endif
    endelse
        
    ;iteration = iteration + 1
    ; DEALLOCATE(RESULTS) ; thanks Bruce!
  endfor
  ; end of iterations

  if (diagnostic_type eq 'D') or (diagnostic_type eq 'd') then begin
    result1 = check_value[2]
    ; print*,check_value(2)
  endif else begin
    result1 = check_value[1]
  endelse
  return, result1
end
