; docformat = 'rst'

function calc_populations, temperature=temperature, density=density, $
                           elj_data=elj_data, omij_data=omij_data, $
                           aij_data=aij_data, $
                           eff_Omij=eff_Omij, $
                           level_num=level_num, irats=irats
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
;     eff_Omij    :   in, type=array/object
;                     effective collision strengths (Omij_T) at given temperature
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
;        Atomic Level Populations:    0.96992832    0.0070036315     0.023062261   2.6593671e-06   3.1277019e-06
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
;   0.3.0
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
;     27/02/2019, A. Danehkar, Simplify the calc_populations() routine 
;                        for external usage.
;                            
;     04/03/2019, A. Danehkar, Use the get_omij_temp() routine.
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
;     EFF_OMIJ    : in, type=array/object, effective collision strengths (Omij_T) at given temperature
;     LEVEL_NUM   : in, type=int, Number of levels
;     IRATS       : in, type=int, Else Coll. rates = tabulated values * 10 ** irats
;     
; OUTPUTS:  This function returns a array/object as the atomic level populations (N_j).
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
;     > Atomic Level Populations:   0.96992832    0.0070036315     0.023062261   2.6593671e-06   3.1277019e-06
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
;     27/02/2019, A. Danehkar, simplify the calc_populations() routine 
;                        for external usage.    
;     04/03/2019, A. Danehkar, use the get_omij_temp() routine.
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
  
  h_Planck = 4.13566766225D-15 ; eV.s ;6.62606957e-27 ; erg.s
  c_Speed = 2.99792458D10 ; cm/s
  k_B = 8.617330350D-5 ; eV/K ; 1.3806485279e-16 ; erg/K
  
  pi=3.1415926535897931D0
  h_bar_Planck = 1.054571800D-27 ;erg.s
  me = 9.10938356D-28 ; gram
  k_B_erg = 1.38064852D-16 ; erg/K
  Beta1 = h_bar_Planck^2*sqrt(2*pi/(k_B_erg*me))/me ; 8.629D-06 
  
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
  if keyword_set(irats) eq 0 then begin
    irats=0
  endif
  
  T_lin_list = omij_data[0].strength
  T_log_list = alog10(T_lin_list) ; temperature intervals (array)
  
  Aij =aij_data.Aij ; Transition Probabilities (A_ij)
  Ej =elj_data.Ej ; Energy Levels (E_j) in cm-1
  Gj =long(elj_data.J_v*2.+1.) ; Ground Levels (G_j)
  
  Qj=dblarr(T_num)   
  Qij=dblarr(level_num,level_num) 
  Equilib_Eqs=dblarr(level_num,level_num)  
  Nlj=dblarr(level_num)
  
  Qij[*,*]=double(0.0)
  Equilib_Eqs[*,*]=double(0.0)
  Nlj[*]=double(0.0)
  
  T_log = alog10(temperature)
  
  if (T_num eq 1) then begin
    print, 'Coll. strengths available for 1 Te only - assuming const'
  endif else begin
    if (T_num eq 2) then begin
      print, 'Coll. strengths available for 2 Te only - linear interp'
    endif
  endelse
  ; Derive the interpolated effective collision strengths (Omij_T) from collision strengths data (Omij)
  ; Obtain collisional de-excitation and excitation rates (Qij) from the effective collision strengths (Omij_T)
  if keyword_set(eff_Omij) eq 0 then begin
    Omij_T=get_omij_temp(temperature=temperature, omij_data=omij_data, level_num=level_num, irats=irats)
  endif else begin
    Omij_T=eff_Omij 
  endelse
  for I = 2, level_num do begin
    for J = I, level_num do begin
      d_E = double(Ej[J-1]-Ej[I-2])*h_Planck*c_Speed ; delta Energy in eV; convert from cm-1 to eV
      ; Calculate the Boltzmann factor
      exp_dE_kT = exp(-d_E/(k_B*temperature)) ; Maxwell-Boltzmann distribution      
      ; Obtain collisional de-excitation and excitation rates from the effective collision strengths Omij_T
      if (irats eq 0) then begin
        Qij[I-2,J-1] = Beta1*Omij_T[I-2,J-1]*exp_dE_kT / (double(Gj[I-2])*sqrt(temperature)) ; collisional excitation rates
        Qij[J-1,I-2] = Beta1*Omij_T[I-2,J-1] / (double(Gj[J-1])*sqrt(temperature)) ; collisional de-excitation rates
      endif else begin
        Qij[I-2,J-1] = Omij_T[I-2,J-1]*exp_dE_kT*10.^irats ; collisional excitation rates
        Qij[J-1,I-2] = double(Gj[I-2])*Qij[I-2,J-1] / (exp_dE_kT*double(Gj[J-1])) ; collisional de-excitation rates
      endelse
    endfor
  endfor
  ; Calculate the critical densities
  ; N_crit_i = Sum_{j} (Aij) / Sum_{j} (Qij)
  ;A_i_sum = TOTAL(Aij, 2) ; Sum each of the columns in Aij
  ;Q_i_sum = TOTAL(Qij, 2) ; Sum each of the columns in Qij
  ;N_crit=A_i_tot/Q_i_tot ; critical densities
  for I = 2, level_num do begin
    for J = 1, level_num do begin
      if (J ne I) then begin
        ; the equations for the equilibrium level populations:
        ; collisional de-excitation eqs -  collisional excitation eqs = 0
        ; Sum_{j ne i} (Ne * Nj * Qji) + Sum_{j > i} (Nj Aji)
        ;    - (Sum_{j ne i} (Ne * Ni * Qij) + Sum_{j < i} (Ni Aij)) = 0
        Equilib_Eqs[I-1,J-1] = Equilib_Eqs[I-1,J-1] + density * Qij[J-1,I-1] ; collisional de-excitation
        Equilib_Eqs[I-1,I-1] = Equilib_Eqs[I-1,I-1] - density * Qij[I-1,J-1] ; collisional excitation
        if (J gt I) then begin
          Equilib_Eqs[I-1,J-1] = Equilib_Eqs[I-1,J-1] + Aij[J-1,I-1] ; collisional de-excitation
        endif else begin
          Equilib_Eqs[I-1,I-1] = Equilib_Eqs[I-1,I-1] - Aij[I-1,J-1] ; collisional excitation
        endelse
      endif
    endfor
  endfor
  B0 = double(-Equilib_Eqs[1:level_num-1,0])
  Equilib_Eqs[0:level_num-2,0:level_num-2] = Equilib_Eqs[1:level_num-1,1:level_num-1]
  A0=transpose(Equilib_Eqs[0:level_num - 2,0:level_num - 2])
  ; Solve the equations for the equilibrium level populations
  ; A.X = B
  ; A: Matrix for the equilibrium level populations equations (i,j), Equilib_Eqs[*,*]
  ; B: Vector for the equilibrium level populations equations (i,j=0), Equilib_Eqs[*,0]
  ; X: Ionic population density (j), Nj
  ; X=la_linear_equation(A0, B0[0:level_num-2]); this function does not work in GDL!
  ludc, A0, Index0  ; Decompose A0; supported by GDL
  ; obtain the ionic population densities (j)
  X0 = lusol(A0, Index0, B0[0:level_num - 2]) ; Compute the solution X; supported by GDL
  ; Calculate the atomic level populations (Nlj) from the ionic population densities (Nj)
  Nlj[1:level_num-1] = X0[0:level_num-2]
  Nlj[0] = 1.D0
  N_tot = double(total(Nlj[0:level_num-1]))
  Nlj[0:level_num-1] = Nlj[0:level_num-1] / N_tot
  return, Nlj
end
