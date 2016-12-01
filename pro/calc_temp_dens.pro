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
  
  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  GX= long(0)
  ID=lonarr(2+1)
  JD=lonarr(2+1)
  iteration= long(0)
  
  I= long(0) 
  I1= long(0) 
  I2= long(0) 
  J= long(0) 
  K= long(0) 
  L= long(0) 
  KK= long(0)
  JT= long(0) 
  JJD= long(0)
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
  IP1= long(0) 
  IKT= long(0) 
  IA= long(0) 
  IB= long(0) 
  IA1= long(0) 
  IA2= long(0) 
  IB1= long(0) 
  IB2= long(0)
     
  TEMPI=double(0) 
  TINC=double(0)
  DENSI=double(0) 
  DINC=double(0)
  DENS=double(0)
  TEMP=double(0)
  EJI=double(0)
  WAV=double(0)
  SUMA=double(0)
  SUMB=double(0)
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
   
  ion1=strtrim(ion,1)
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  readf,lun1,NLINES
  for i=1, NLINES do begin 
    readf,lun1,LTEXT
  endfor
  ; Read no. of levels =  NLEV,
  readf,lun1,NLEV, NTEMP
  
  Glj=lonarr(NLEV+1)
  ITRANA=lonarr(2+1,NLEV+1)
  ITRANB=lonarr(2+1,NLEV+1)
  
  Nlj=dblarr(NLEV+1)
  WAVA=dblarr(NLEV+1)
  WAVB=dblarr(NLEV+1)
  Omij=dblarr(NTEMP+1,NLEV+1,NLEV+1)   
  Aij=dblarr(NLEV+1,NLEV+1)   
  Elj=dblarr(NLEV+1)   
  Telist=dblarr(NTEMP+1)
  check_value=dblarr(3+1)
     
  LABEL1=STRARR(NLEV+1)
  
  itrana[*,*]=0
  itranb[*,*]=0
  
  levu_str=strsplit(levu, ',', ESCAPE='/', /EXTRACT)
  levl_str=strsplit(levl, ',', ESCAPE='/', /EXTRACT)
  
  temp=size(levu_str, /N_ELEMENTS)
  levu_num=long(temp[0]/2)
  temp=size(levl_str, /N_ELEMENTS)
  levl_num=long(temp[0]/2)
  
  levu_i=0
  for i=1, levu_num do begin 
    res=equib_strnumber(levu_str[levu_i], val)
    if res eq 1 then ITRANA[1,i]=long(val)
    res=equib_strnumber(levu_str[levu_i+1], val)
    if res eq 1 then ITRANA[2,i]=long(val)
    levu_i = levu_i + 2
    ;if levu_i ge 2*levu_num then break
  endfor

  levl_i=0
  for i=1, levl_num do begin 
    res=equib_strnumber(levl_str[levl_i], val)
    if res eq 1 then ITRANB[1,i]=long(val)
    res=equib_strnumber(levl_str[levl_i+1], val)
    if res eq 1 then ITRANB[2,i]=long(val)
    levl_i = levl_i + 2
    ;if levl_i ge 2*levl_num then break;
  endfor
  
  ; no. of Te = NTEMP and the
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
    Telist[i] = ddtemp
    Telist[i] = alog10(Telist[i])
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
      Omij[K,I,J] = QX
    endwhile
  endif
  if(IBIG eq 1) or (IBIG eq 2) then begin
    readf,lun1,NTRA
    for IN = 1, NTRA do begin
      readf,lun1,I,J,Omij[1:NTEMP,I,J]
      ;READ(1,*) I,J,(Omij(ITEMP,I,J),ITEMP=1,NTEMP)
    endfor
  endif
  ; Read transition probabilities
  NLEV1 = NLEV - 1
  if (IBIG eq 1) then begin
    readf,lun1,I,J,Aij[J,I];,L=K+1,NLEV),K=1,NLEV1
    ;READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
  endif else begin
    for K = 1, NLEV1 do begin
      KP1 = K + 1
      for L = KP1, NLEV do begin
        readf,lun1, I, J, AX
        Aij[J,I] = AX
      endfor
    endfor
  endelse
  ; Read statistical weights, energy levels (cm-1)
  for J = 1, NLEV do begin
    readf,lun1,  I, GX, EX
    Glj[I] = GX
    Elj[I] = EX
  endfor
  free_lun, lun1
  ; start of iterations
  ; ****************************
  for iteration = 1, 9 do begin
    if (diagtype eq 't') or (diagtype eq 'T') then begin
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
        if (TEMP le 0.D0) or (DENS le 0.D0) then begin
            print,'Temp = ', TEMP, ', Dens = ', DENS
            return, 0
        endif

        Nlj=calc_populations(TEMP, DENS, Telist, Omij, Aij, Elj, Glj, NLEV, NTEMP, IRATS)
        
        ; Search ITRANA, ITRANB for transitions & sum up  
        SUMA=double(0.0)
        SUMB=double(0.0)
        for IKT=1, levu_num do begin 
          I=ITRANA[1,IKT]
          J=ITRANA[2,IKT]
          emissivity_line=double(0.0)
          if (Aij[J,I] ne 0.D0) then begin
            EJI = Elj[J] - Elj[I]
            WAV = 1.D8 / EJI
            emissivity_line=Nlj[J]*Aij[J,I]*h_Planck*c_Speed*1.e8/WAV
            SUMA=SUMA+emissivity_line
          endif
        endfor
        for IKT=1, levl_num do begin 
          I=ITRANB[1,IKT]
          J=ITRANB[2,IKT]
          emissivity_line=double(0.0)
          if (Aij[J,I] ne 0.D0) then begin
            EJI = Elj[J] - Elj[I]
            WAV = 1.D8 / EJI
            emissivity_line=Nlj[J]*Aij[J,I]*h_Planck*c_Speed*1.e8/WAV
            SUMB=SUMB+emissivity_line
          endif
        endfor
        FRAT=SUMA/SUMB
        if (diagtype eq 't') or (diagtype eq 'T') then begin
          RESULTS[1, JT] = TEMP
          RESULTS[2, JT] = DENS
          RESultS[3, JT] = FRAT-inratio
        endif else begin
          RESULTS[1, JJD] = TEMP
          RESULTS[2, JJD] = DENS
          RESultS[3, JJD] = FRAT-inratio
        endelse ;End of the Ne iteration
      endfor
      for IA = 1, levu_num do begin
        I1=ITRANA[1,IA]
        I2=ITRANA[2,IA]
        DEE=Elj[I2]-Elj[I1]
        WAVA[IA]=1.D8/DEE
      endfor
      for IB = 1, levl_num do begin
        I1=ITRANB[1,IB]
        I2=ITRANB[2,IB]
        DEE=Elj[I2]-Elj[I1]
        WAVB[IB]=1.D8/DEE
      endfor
    ; End of the temperature iteration
    endfor
    ; iteration and detect the sign change.
    if (diagtype eq 'D') or (diagtype eq 'd') then begin
      INT = ind
    endif
    ; iteration and detect the sign change.
    for I=2,INT do begin
      test=0
      if (equib_sign(results[3,I],results[3,1]) ne results[3,I]) then begin 
        ;if this condition, the values have a different sign
        check_value[*] = results[*,I-1] ; the value before the sign change returned
        test=1
        break
      endif
    endfor
    if(test eq 0) and (iteration lt 9) then begin ; check if there is any change of sign,
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
      if (test eq 0) and (iteration eq 9) then begin ;check if no change of sign,
                             ;and checks if it should be upper or lower limit
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
  endfor
  ; end of iterations
  ;****************************
  if (diagtype eq 'D') or (diagtype eq 'd') then begin
    result1 = check_value[2]
  endif else begin
    result1 = check_value[1]
  endelse
  return, result1
end
