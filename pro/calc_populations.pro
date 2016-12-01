function calc_populations, TEMP, DENS, Telist, Omij, Aij, Elj, Glj, NLEV, NTEMP, IRATS
;+
; NAME:
;     calc_populations
; PURPOSE:
;     solve atomic level populations in statistical equilibrium 
;     for given electron temperature and density.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     Nlj=calc_populations(TEMP, DENS, Telist, Omij, Aij, Elj, Glj, NLEV, NTEMP, IRATS)
;
; INPUTS:
;     TEMP -     electron temperature
;     DENS -     electron density
;     Telist -   temperature intervals (array)
;     Omij - Collision Strengths (Omega_ij)
;     Aij - Transition Probabilities (A_ij)
;     Elj - Energy Levels (E_j)
;     Glj - Ground Levels (G_j)
;     NLEV -Number of levels
;     NTEMP - Number of temperature intervals
;     IRATS - Else Coll. rates = tabulated values * 10 ** IRATS
; RETURN:  N_j (array): atomic level populations
;
; REVISION HISTORY:
;     Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
;     Replaced str2int with strnumber, A. Danehkar, 20/10/2016
;     Replaced CFY, SPLMAT, and CFD with
;          IDL function INTERPOL( /SPLINE), A. Danehkar, 20/10/2016
;     Replaced LUSLV with IDL LAPACK function 
;                       LA_LINEAR_EQUATION, A. Danehkar, 20/10/2016
;     Replaced LA_LINEAR_EQUATION (not work in GDL) with IDL function 
;                             LUDC & LUSOL, A. Danehkar, 15/11/2016
;     Replaced INTERPOL (not accurate) with 
;                    SPL_INIT & SPL_INTERP, A. Danehkar, 19/11/2016
;     Make a new function calc_populations() separated from 
;       calc_abundance(), calc_temp_dens(), A. Danehkar, 20/11/2016
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
; 2009-05    R.Wesson     Misc updates and improvements, inputs from cmd line, 
;                         written only for calculating ionic abundances.
;- 
  
  DD=double(0)
  I= long(0)
  J= long(0) 
  K= long(0)
  IM1= long(0) 
  JM1= long(0) 
  DELTEK=double(0)
  EXPE=double(0)
  SUMN=double(0)
  VALUE=double(0)
  
  CS=dblarr(NLEV+1,NLEV+1)
  QQ=dblarr(NTEMP+1)   
  QEFF=dblarr(NLEV+1,NLEV+1) 
  X=dblarr(NLEV+1,NLEV+1)  
  Y=dblarr(NLEV+1)
  Nlj=dblarr(NLEV+1)
  
  X[*,*]=double(0.0)
  CS[*,*]=double(0.0)
  QEFF[*,*]=double(0.0)
  Y[*]=double(0.0)
    
  NLEV1 = NLEV - 1
  
  TLOGT = alog10(TEMP)
  TEMP2= sqrt(TEMP)
  
  ;IOPT=0
  if (NTEMP eq 1) then begin
    print, 'Coll. strengths available for 1 Te only - assuming const'
  endif else begin
    if (NTEMP eq 2) then begin
      print, 'Coll. strengths available for 2 Te only - linear interp'
    endif
  endelse
  
  for I = 2, NLEV do begin
    for J = I, NLEV do begin
      ;Negative!
      DELTEK = (Elj[I-1]-Elj[J])*1.4388463D0
      EXPE = exp(DELTEK/TEMP)
      for IT = 1, NTEMP do begin
        if (IRATS eq 0.D+00) then begin
          QQ[IT] = Omij[IT,I-1,J]
        endif else begin
          ;Take out the exp. depend.
          QQ[IT] = Omij[IT,I-1,J] / EXPE
          ; before interpolation
        endelse
      endfor
      if (NTEMP eq 1) then begin
        DD = QQ[1]
      endif else begin
        if (NTEMP eq 2) then begin
          DD = QQ[1] +  (QQ[2] - QQ[1])/(Telist[2] - Telist[1]) * (TLOGT - Telist[1])
        endif else begin
          ;DD=interpol(QQ[1:NTEMP], T[1:NTEMP], TLOGT,/SPLINE)
          deriv1 = spl_init(Telist[1:NTEMP], QQ[1:NTEMP])
          DD=spl_interp(Telist[1:NTEMP], QQ[1:NTEMP], deriv1, TLOGT)
        endelse
      endelse
      if (IRATS eq 0.D+00) then begin
        CS[I-1,J] = DD
      endif else begin
        CS[I-1,J] = DD * EXPE
      endelse
      if (IRATS eq 0.D+00) then begin
        QEFF[I-1,J] = 8.63D-06*CS[I-1,J] * EXPE / (Glj[I-1]*TEMP2)
        QEFF[J,I-1] = 8.63D-06 * CS[I-1,J] / (Glj[J]*TEMP2)
      endif else begin
        QEFF[I-1,J] = CS[I-1,J] * 10.^IRATS
        ; Be careful
        QEFF[J,I-1] = Glj[I-1] * QEFF[I-1,J] / (EXPE * Glj[J])
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
          X[I,J] = X[I,J] + Aij[J,I]
        endif else begin 
          X[I,I] = X[I,I] - Aij[I,J]
        endelse
      endif
    endfor
  endfor
  for I = 2, NLEV do begin
    IM1 = I - 1
    VALUE = 0.D0 - X[I,1]
    Y[IM1] = VALUE
    for J = 2, NLEV do begin
      JM1 = J - 1
      VALUE = X[I,J]
      X[IM1,JM1] = VALUE
    endfor
  endfor
  ; Solve matrices for populations
  ; YY=la_linear_equation(transpose(X[1:NLEV1,1:NLEV1]), Y[1:NLEV1]); not work in GDL
  XX=transpose(X[1:NLEV1,1:NLEV1])
  ludc, XX, INDEX  ; supported by GDL
  YY = lusol(XX, INDEX, Y[1:NLEV1]) ; supported by GDL
  Y[1:NLEV1]=YY[0:NLEV1-1]
  for I = NLEV, 2, -1 do begin
    Nlj[I] = Y[I-1]
  endfor
  SUMN = 1.D0
  for I = 2, NLEV do begin
    SUMN = SUMN + Nlj[I]
  endfor
  for I = 2, NLEV do begin
    Nlj[I] = Nlj[I] / SUMN
  endfor
  Nlj[1] = 1.D0 / SUMN
  return, Nlj
end
