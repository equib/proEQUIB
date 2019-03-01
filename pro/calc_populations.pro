; docformat = 'rst'

function calc_populations, temperature=temperature, density=density, $
                           elj_data=elj_data, omij_data=omij_data, $
                           aij_data=aij_data, $
                           coeff_omij=coeff_omij, level_num=level_num, irats=irats
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
;     elj_data    :   in, required, type=array/object
;                            energy levels (Ej) data
;     omij_data   :   in, required, type=array/object
;                            collision strengths (omega_ij) data
;     aij_data    :   in, required, type=array/object
;                            transition probabilities (Aij) data
;     coeff_omij  :   in, type=array/object
;                     Collision Strengths (Omega_ij)
;     level_num   :   in, type=int
;                     Number of levels
;     irats       :   in, type=int
;                     Else Coll. rates = tabulated values * 10 ** irats
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_dir = ['atomic-data', 'chianti70']
;     IDL> Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> atom='s'
;     IDL> ion='ii'
;     IDL> s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej)
;     IDL> s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
;     IDL> s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)\
;     IDL> density = double(1000)
;     IDL> temperature=double(10000.0);
;     IDL> Nlj=calc_populations(temperature=temperature, density=density, $
;     IDL>                      elj_data=s_ii_elj, omij_data=s_ii_omij, $
;     IDL>                      aij_data=s_ii_aij)
;     IDL> print, 'Atomic Level Populations:', Nlj
;        Atomic Level Populations:    0.96992865    0.0070037131     0.023061848   2.6592895e-06   3.1276110e-06
;
; :Categories:
;   Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
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
;     27/06/2019, A. Danehkar, simplify the calc_populations() routine 
;                        for external usage.
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
;                           ELJ_DATA=elj_data, OMIJ_DATA=omij_data, $
;                           AIJ_DATA=aij_data, $
;                           COEFF_OMIJ=coeff_omij, LEVEL_NUM=level_num, IRATS=irats)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE : in, required, type=float, electron temperature
;     DENSITY     : in, required, type=float, electron density
;     ELJ_DATA    : in, required, type=array/object, energy levels (Ej) data
;     OMIJ_DATA   : in, required, type=array/object, collision strengths (omega_ij) data
;     AIJ_DATA    : in, required, type=array/object, transition probabilities (Aij) data
;     COEFF_OMIJ  : in, type=array/object, Collision Strengths (Omega_ij))
;     LEVEL_NUM   : in, type=int, Number of levels
;     IRATS       : in, type=int, Else Coll. rates = tabulated values * 10 ** irats
;     
; OUTPUTS:  This function returns a array/object as the atomic level populations (N_j)
;
; PROCEDURE: This function is called by calc_emissivity, calc_temperature and calc_density.
; 
; EXAMPLE:
;     base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     data_dir = ['atomic-data', 'chianti70']
;     Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
;     Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
;     Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )
;     atom='s'
;     ion='ii'
;     s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej)
;     s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
;     s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)\
;     density = double(1000)
;     temperature=double(10000.0);
;     Nlj=calc_populations(temperature=temperature, density=density, $
;                          elj_data=s_ii_elj, omij_data=s_ii_omij, $
;                          aij_data=s_ii_aij)
;      print, 'Atomic Level Populations:', Nlj
;     > Atomic Level Populations:   0.96992865    0.0070037131     0.023061848   2.6592895e-06   3.1276110e-06
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
;     27/06/2019, A. Danehkar, simplify the calc_populations() routine 
;                        for external usage.
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
  if keyword_set(elj_data) eq 0 then begin 
    print,'elj_data is not set'
    return, 0
  endif 
  if keyword_set(omij_data) eq 0 then begin 
    print,'omij_data is not set'
    return, 0
  endif 
  if keyword_set(aij_data) eq 0 then begin 
    print,'aij_data is not set'
    return, 0
  endif
  if keyword_set(level_num) eq 0 then begin 
    temp=size(elj_data,/DIMENSIONS)
    level_num=temp[0]
  endif
  temp=size(omij_data[0].strength,/DIMENSIONS)
  T_num=temp[0] ; Number of temperature intervals
  if keyword_set(coeff_omij) eq 0 then begin
    temp=size(omij_data,/DIMENSIONS)
    omij_num=temp[0]
    Omij=dblarr(T_num,level_num,level_num)   
    for k = 1, omij_num-1 do begin
      I = omij_data[k].level1
      J = omij_data[k].level2
      if I le level_num and J le level_num then begin
        Omij[0:T_num-1,I-1,J-1] = omij_data[k].strength
      endif
    endfor
  endif else begin
    Omij=coeff_omij
  endelse
  if keyword_set(irats) eq 0 then begin
    irats=0
  endif
  
  T_lin_list = omij_data[0].strength
  T_log_list = alog10(T_lin_list) ; temperature intervals (array)
  
  Aij =aij_data.AIJ ; Transition Probabilities (A_ij)
  Elj =elj_data.Ej ; Energy Levels (E_j)
  Glj =long(elj_data.J_v*2.+1.) ; Ground Levels (G_j)
  
  Q_vector=dblarr(T_num)   
  Q_eff=dblarr(level_num,level_num) 
  X=dblarr(level_num,level_num)  
  Y=dblarr(level_num)
  Nlj=dblarr(level_num)
  
  X[*,*]=double(0.0)
  Q_eff[*,*]=double(0.0)
  Y[*]=double(0.0)
  Nlj[*]=double(0.0)
  
  T_log = alog10(temperature)
  
  if (T_num eq 1) then begin
    print, 'Coll. strengths available for 1 Te only - assuming const'
  endif else begin
    if (T_num eq 2) then begin
      print, 'Coll. strengths available for 2 Te only - linear interp'
    endif
  endelse
  
  for I = 2, level_num do begin
    for J = I, level_num do begin
      d_E = (Elj[I-2]-Elj[J-1])*1.4388463D0 ;It is negative!
      exp_dE_T = exp(d_E/temperature)
      if (irats eq 0.D+00) then begin
        Q_vector[*] = Omij[*,I-2,J-1]
      endif else begin
        Q_vector[*] = Omij[*,I-2,J-1] / exp_dE_T ;Take out the exp. before interpolation
      endelse
      if (T_num eq 1) then begin
        Q_interp = Q_vector[0]
      endif else begin
        if (T_num eq 2) then begin
          Q_interp = Q_vector[0] +  (Q_vector[1] - Q_vector[0])/(T_log_list[1] - T_log_list[0]) * (T_log - T_log_list[0])
        endif else begin
          ;Q_interp=interpol(Q_vector[1:T_num], T[1:T_num], T_log,/SPLINE)
          Q_init = spl_init(T_log_list[0:T_num-1], Q_vector[0:T_num-1])
          Q_interp=spl_interp(T_log_list[0:T_num-1], Q_vector[0:T_num-1], Q_init, T_log)
        endelse
      endelse
      if (irats eq 0.D+00) then begin
        Q_eff[I-2,J-1] = 8.63D-06*Q_interp * exp_dE_T / (double(Glj[I-2])*sqrt(temperature))
        Q_eff[J-1,I-2] = 8.63D-06*Q_interp / (double(Glj[J-1])*sqrt(temperature))
      endif else begin
        Q_eff[I-2,J-1] = Q_interp * exp_dE_T * 10.^irats
        Q_eff[J-1,I-2] = double(Glj[I-2]) * Q_eff[I-2,J-1] / (exp_dE_T * double(Glj[J-1])) ; Be careful G integer!
      endelse
    endfor
  endfor
  for I = 2, level_num do begin
    for J = 1, level_num do begin
      if (J ne I) then begin
        X[I-1,J-1] = X[I-1,J-1] + density * Q_eff[J-1,I-1]
        X[I-1,I-1] = X[I-1,I-1] - density * Q_eff[I-1,J-1]
        if (J gt I) then begin
          X[I-1,J-1] = X[I-1,J-1] + Aij[J-1,I-1]
        endif else begin 
          X[I-1,I-1] = X[I-1,I-1] - Aij[I-1,J-1]
        endelse
      endif
    endfor
  endfor
  Y[0:level_num-2] = 0.D0 - X[1:level_num-1,0]
  X[0:level_num-2,0:level_num-2] = X[1:level_num-1,1:level_num-1]
  ; Solve matrices for populations
  ; YY=la_linear_equation(transpose(X[0:level_num-2,0:level_num-2]), Y[0:level_num-2]); this function does not work in GDL!
  XX=transpose(X[0:level_num - 2,0:level_num - 2])
  ludc, XX, INDEX  ; supported by GDL
  YY = lusol(XX, INDEX, Y[0:level_num - 2]) ; supported by GDL
  Nlj[1:level_num-1] = YY[0:level_num-2]
  Nlj[0] = 1.D0
  pop_sum = double(total(Nlj[0:level_num-1]))
  Nlj[0:level_num-1] = Nlj[0:level_num-1] / pop_sum
  return, Nlj
end
