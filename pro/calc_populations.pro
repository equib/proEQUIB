; docformat = 'rst'

function calc_populations, temperature=temperature, density=density, $
                           temp_list=temp_list, $ 
                           Omij=Omij, Aij=Aij, Elj=Elj, $
                           Glj=Glj, level_num=level_num, $
                           temp_num=temp_num, irats=irats
;+
;     This function solves atomic level populations in statistical equilibrium 
;     for given electron temperature and density.
;
; :Returns:
;    type=array/object. This function returns the atomic level populations.
;
; :Keywords:
;     temperature :   in, required, type=float
;                     electron temperature
;     density     :   in, required, type=float
;                     electron density
;     temp_list   :   in, required, type=array
;                     temperature intervals (array)
;     Omij        :   in, required, type=array/object
;                     Collision Strengths (Omega_ij)
;     Aij         :   in, required, type=array/object
;                     Transition Probabilities (A_ij)
;     Elj         :   in, required, type=array
;                     Energy Levels (E_j)
;     Glj         :   in, required, type=array
;                     Ground Levels (G_j)
;     level_num   :   in, required, type=int
;                     Number of levels
;     temp_num    :   in, required, type=int
;                     Number of temperature intervals
;     irats       :     in, required, type=int
;                     Else Coll. rates = tabulated values * 10 ** irats
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.0.6
;
; :History:
;     15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.
;
;     20/10/2016, A. Danehkar, Replaced str2int with strnumber.
;
;     20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
;          IDL function INTERPOL( /SPLINE).
;
;     20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
;                       LA_LINEAR_EQUATION.
;
;     15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
;           with IDL function LUDC & LUSOL.
;
;     19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
;                    SPL_INIT & SPL_INTERP.
;
;     20/11/2016, A. Danehkar, Made a new function calc_populations()
;       for solving atomic level populations and separated it from
;       calc_abundance(), calc_density() and calc_temperature().
;
;     10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
;                      input elj_data, omij_data, aij_data.
;
;     12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
;                        from calc_populations().
;
; FORTRAN HISTORY:
;
;     03/05/1981, I.D.Howarth,  Version 1.
;
;     05/05/1981, I.D.Howarth,  Minibug fixed!
;
;     07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
;
;     03/08/1981, S.Adams,      Interpolates collision strengths.
;
;     07/08/1981, S.Adams,      Input method changed.
;
;     19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
;                               filenames given to SA's data files.
;
;     08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
;
;     02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
;                         modified such that matrix sizes (i.e. maximum
;                         of Te and maximum no of levels) can now be cha
;                         by modifying the parameters NDIM1, NDIM2 and N
;                         in the Main program. EASY!
;                         Now takes collision rates as well.
;                         All variables are declared explicitly
;                         Generate two extra files (ionpop.lis and ionra
;                         of plain stream format for plotting.
;
;     06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
;                         Fixed readin bug for IBIG=2 case.
;                         Now reads reformatted upsilons (easier to see
;                         and the 0 0 0 data end is excluded for these c
;                         The A values have a different format for IBIG=.
;
;     2006, B.Ercolano,   Converted to F90.
;-

;+
; NAME:
;     calc_populations
; PURPOSE:
;     This function solves atomic level populations in statistical equilibrium 
;     for given electron temperature and density.
;
; CALLING SEQUENCE:
;     Result = calc_populations(TEMPERATURE=temperature, DENSITY=density, $
;                        TEMP_LIST=temp_list, $ 
;                        OMIJ=Omij, AIJ=Aij, ELJ=Elj, GLJ=Glj, $
;                        LEVEL_NUM=level_num, TEMP_NUM=temp_num, 
;                        IRATS=irats)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE : in, required, type=float, electron temperature
;     DENSITY     : in, required, type=float, electron density
;     TEMP_LIST   : in, required, type=array, temperature intervals (array)
;     OMIJ        : in, required, type=array/object, Collision Strengths (Omega_ij)
;     AIJ         : in, required, type=array/object, Transition Probabilities (A_ij)
;     ELJ         : in, required, type=array, Energy Levels (E_j)
;     GLJ         : in, required, type=array, Ground Levels (G_j)
;     LEVEL_NUM   : in, required, type=int, Number of levels
;     TEMP_NUM    : in, required, type=int, Number of temperature intervals
;     IRATS       : in, required, type=int, Else Coll. rates = tabulated values * 10 ** irats
;     
; OUTPUTS:  This function returns a array/object as the atomic level populations (N_j)
;
; PROCEDURE: This function is called by calc_emissivity, calc_temperature and calc_density.
;
; MODIFICATION HISTORY:
;     15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.
;     20/10/2016, A. Danehkar, Replaced str2int with strnumber.
;     20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
;          IDL function INTERPOL( /SPLINE).
;     20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
;                       LA_LINEAR_EQUATION.
;     15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
;           with IDL function LUDC & LUSOL.
;     19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
;                    SPL_INIT & SPL_INTERP.
;     20/11/2016, A. Danehkar, Made a new function calc_populations()
;       for solving atomic level populations and separated it from
;       calc_abundance(), calc_density() and calc_temperature().
;     10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
;                      input elj_data, omij_data, aij_data.
;     12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
;                        from calc_populations().
; 
; FORTRAN HISTORY:
;     03/05/1981, I.D.Howarth,  Version 1.
;     05/05/1981, I.D.Howarth,  Minibug fixed!
;     07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
;     03/08/1981, S.Adams,      Interpolates collision strengths.
;     07/08/1981, S.Adams,      Input method changed.
;     19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
;                               filenames given to SA's data files.
;     08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
;     02/1996, X.W.Liu,    Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
;                          modified such that matrix sizes (i.e. maximum
;                          of Te and maximum no of levels) can now be cha
;                          by modifying the parameters NDIM1, NDIM2 and N
;                          in the Main program. EASY!
;                          Now takes collision rates as well.
;                          All variables are declared explicitly
;                          Generate two extra files (ionpop.lis and ionra
;                          of plain stream format for plotting.
;     06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
;                          Fixed readin bug for IBIG=2 case.
;                          Now reads reformatted upsilons (easier to see
;                          and the 0 0 0 data end is excluded for these c
;                          The A values have a different format for IBIG=.
;     2006, B.Ercolano,    Converted to F90.
;- 
  
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(temp_list) eq 0 then begin 
    print,'temp_list is not set'
    return, 0
  endif
  if keyword_set(Omij) eq 0 then begin 
    print,'Omij is not set'
    return, 0
  endif
  if keyword_set(Aij) eq 0 then begin 
    print,'Aij is not set'
    return, 0
  endif 
  if keyword_set(Elj) eq 0 then begin 
    print,'Elj is not set'
    return, 0
  endif  
  if keyword_set(Glj) eq 0 then begin 
    print,'Glj is not set'
    return, 0
  endif  
  if keyword_set(level_num) eq 0 then begin 
    print,'level_num is not set'
    return, 0
  endif
  if keyword_set(temp_num) eq 0 then begin 
    print,'temp_num is not set'
    return, 0
  endif 
  DD=double(0)
  I= long(0)
  J= long(0) 
  K= long(0)
  IM1= long(0) 
  JM1= long(0) 
  DELTEK=double(0)
  EXPE=double(0)
  pop_sum=double(0)
  VALUE=double(0)
  
  CS=dblarr(level_num+1,level_num+1)
  QQ=dblarr(temp_num+1)   
  QEFF=dblarr(level_num+1,level_num+1) 
  X=dblarr(level_num+1,level_num+1)  
  Y=dblarr(level_num+1)
  Nlj=dblarr(level_num+1)
  
  X[*,*]=double(0.0)
  CS[*,*]=double(0.0)
  QEFF[*,*]=double(0.0)
  Y[*]=double(0.0)
    
  level_num1 = level_num - 1
  
  TLOGT = alog10(temperature)
  TEMP2= sqrt(temperature)
  
  ;IOPT=0
  if (temp_num eq 1) then begin
    print, 'Coll. strengths available for 1 Te only - assuming const'
  endif else begin
    if (temp_num eq 2) then begin
      print, 'Coll. strengths available for 2 Te only - linear interp'
    endif
  endelse
  
  for I = 2, level_num do begin
    for J = I, level_num do begin
      ;Negative!
      DELTEK = (Elj[I-2]-Elj[J-1])*1.4388463D0
      EXPE = exp(DELTEK/temperature)
      for IT = 1, temp_num do begin
        if (irats eq 0.D+00) then begin
          QQ[IT] = Omij[IT-1,I-2,J-1]
        endif else begin
          ;Take out the exp. depend.
          QQ[IT] = Omij[IT-1,I-2,J-1] / EXPE
          ; before interpolation
        endelse
      endfor
      if (temp_num eq 1) then begin
        DD = QQ[1]
      endif else begin
        if (temp_num eq 2) then begin
          DD = QQ[1] +  (QQ[2] - QQ[1])/(temp_list[2-1] - temp_list[1-1]) * (TLOGT - temp_list[1-1])
        endif else begin
          ;DD=interpol(QQ[1:temp_num], T[1:temp_num], TLOGT,/SPLINE)
          deriv1 = spl_init(temp_list[0:temp_num-1], QQ[1:temp_num])
          DD=spl_interp(temp_list[0:temp_num-1], QQ[1:temp_num], deriv1, TLOGT)
        endelse
      endelse
      if (irats eq 0.D+00) then begin
        CS[I-1,J] = DD
      endif else begin
        CS[I-1,J] = DD * EXPE
      endelse
      if (irats eq 0.D+00) then begin
        QEFF[I-1,J] = 8.63D-06*CS[I-1,J] * EXPE / (Glj[I-2]*TEMP2)
        QEFF[J,I-1] = 8.63D-06 * CS[I-1,J] / (Glj[J-1]*TEMP2)
      endif else begin
        QEFF[I-1,J] = CS[I-1,J] * 10.^irats
        ; Be careful
        QEFF[J,I-1] = Glj[I-2] * QEFF[I-1,J] / (EXPE * Glj[J-1])
        ; G integer!
      endelse
    endfor
  endfor
  for I = 2, level_num do begin
    for J = 1, level_num do begin
      if (J ne I) then begin
        X[I,J] = X[I,J] + density * QEFF[J,I]
        X[I,I] = X[I,I] - density * QEFF[I,J]
        if (J gt I) then begin
          X[I,J] = X[I,J] + Aij[J-1,I-1]
        endif else begin 
          X[I,I] = X[I,I] - Aij[I-1,J-1]
        endelse
      endif
    endfor
  endfor
  for I = 2, level_num do begin
    IM1 = I - 1
    VALUE = 0.D0 - X[I,1]
    Y[IM1] = VALUE
    for J = 2, level_num do begin
      JM1 = J - 1
      VALUE = X[I,J]
      X[IM1,JM1] = VALUE
    endfor
  endfor
  ; Solve matrices for populations
  ; YY=la_linear_equation(transpose(X[1:level_num1,1:level_num1]), Y[1:level_num1]); not work in GDL
  XX=transpose(X[1:level_num1,1:level_num1])
  ludc, XX, INDEX  ; supported by GDL
  YY = lusol(XX, INDEX, Y[1:level_num1]) ; supported by GDL
  Y[1:level_num1]=YY[0:level_num1-1]
  for I = level_num, 2, -1 do begin
    Nlj[I] = Y[I-1]
  endfor
  pop_sum = 1.D0
  for I = 2, level_num do begin
    pop_sum = pop_sum + Nlj[I]
  endfor
  for I = 2, level_num do begin
    Nlj[I] = Nlj[I] / pop_sum
  endfor
  Nlj[1] = 1.D0 / pop_sum
  return, Nlj
end
