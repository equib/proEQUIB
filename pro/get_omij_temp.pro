; docformat = 'rst'

function get_omij_temp, temperature=temperature, $
                        omij_data=omij_data, elj_data=elj_data, $
                        level_num=level_num, irats=irats
;+
;     This function derives the effective collision strengths (Omij_T) from 
;     the collision strengths (omega_ij) data for the given temperature.
;
; :Returns:
;    type=array/object. This function returns the effective collision strengths (Omij_T).
;
; :Keywords:
;     temperature :   in, required, type=float
;                     electron temperature
;     omij_data   :   in, required, type=array/object
;                     collision strengths (omega_ij) data
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
;     IDL> temperature=double(10000.0);
;     IDL> Omij_T=get_omij_temp(temperature=temperature, omij_data=s_ii_omij)
;     IDL> print, 'Effective Collision Strengths: '
;     IDL> print, Omij_T
;        Effective Collision Strengths: 
;        0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        2.7800000       0.0000000       0.0000000       0.0000000       0.0000000
;        4.1600000       7.4600000       0.0000000       0.0000000       0.0000000
;        1.1700000       1.8000000       2.2000000       0.0000000       0.0000000
;        2.3500000       3.0000000       4.9900000       2.7100000       0.0000000
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
;   0.0.1
;
; :History:
;     04/03/2019, A. Danehkar, create the get_omij_temp() routine.
;
;-

  h_Planck = 4.13566766225D-15 ; eV.s ;6.62606957e-27 ; erg.s
  c_Speed = 2.99792458D10 ; cm/s
  k_B = 8.617330350D-5 ; eV/K ; 1.3806485279e-16 ; erg/K

  if keyword_set(temperature) eq 0 then begin
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(omij_data) eq 0 then begin
    print,'omij_data is not set'
    return, 0
  endif
  if keyword_set(level_num) eq 0 then begin
    if keyword_set(elj_data) eq 0 then begin
      level_num=max([max(omij_data[*].level1),max(omij_data[*].level2)])
    endif else begin
      temp=size(elj_data,/DIMENSIONS)
      level_num=temp[0]
    endelse
  endif
  if keyword_set(irats) eq 0 then begin
    irats=0
  endif else begin
    if keyword_set(elj_data) eq 0 then begin
        print,'elj_data is not set. It is required for irats'
        return, 0
    endif
    Ej =elj_data.Ej ; Energy Levels (E_j) in cm-1
  endelse
  T_log = double(alog10(temperature))
  temp=size(omij_data[0].strength,/DIMENSIONS)
  T_num=temp[0] ; Number of temperature intervals
  T_lin_list = double(omij_data[0].strength)
  T_log_list = alog10(T_lin_list) ; temperature intervals (array)
  temp=size(omij_data,/DIMENSIONS)
  omij_num=temp[0]
  Omij_T=dblarr(level_num,level_num)   
  for k = 1, omij_num-1 do begin
    I = omij_data[k].level1
    J = omij_data[k].level2
    if I le level_num and J le level_num then begin
      Qj=double(omij_data[k].strength)
      if (irats ne 0) then begin
        d_E = double(Ej[J-1]-Ej[I-2])*h_Planck*c_Speed ; delta Energy in eV; convert from cm-1 to eV
        ; Calculate the Boltzmann factor
        exp_dE_kT = exp(-d_E/(k_B*temperature)) ; Maxwell-Boltzmann distribution      
        Qj = Qj / exp_dE_kT ;Take out the exp. before interpolation
      endif
      if (T_num eq 1) then begin
        Omij_T[I-1,J-1]=Qj
      endif else begin
        if (T_num eq 2) then begin
          Omij_T[I-1,J-1] = Qj[0] +  (Qj[1] - Qj[0])/(T_log_list[1] - T_log_list[0]) * (T_log - T_log_list[0])
        endif else begin
          ;Qj_T=interpol(Qj, T_log_list, T_log, /SPLINE)
          ; Calculate interpolating cubic spline
          Qj_2 = spl_init(T_log_list, Qj)
          ; Calculate the interpolated Omij_T values corresponding to T_log
          ; Obtain the effective collision strengths Omij_T
          Qj_T=spl_interp(T_log_list, Qj, Qj_2, T_log, /Double)
          Omij_T[I-1,J-1] = Qj_T
        endelse
      endelse
    endif
  endfor
  return,Omij_T
end
