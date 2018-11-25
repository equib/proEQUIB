function calc_emissivity, temperature=temperature, density=density, $
                         atomic_levels=atomic_levels, $
                         elj_data=elj_data, omij_data=omij_data, $
                         aij_data=aij_data
;+
; NAME:
;     calc_emissivity
; PURPOSE:
;     calculate line emissivities for specified ion with level(s) by 
;     solving atomic level populations and in statistical equilibrium 
;     for given electron density and temperature.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     path='proEQUIB/atomic-data/'
;     set_atomic_data_path, path
;     Atom_Elj_file='idllib/AtomNeb/atomic-data/chianti70/AtomElj.fits'
;     Atom_Omij_file='idllib/AtomNeb/atomic-data/chianti70/AtomOmij.fits'
;     Atom_Aij_file='idllib/AtomNeb/atomic-data/chianti70/AtomAij.fits'
;     atom='o'
;     ion='iii'
;     o_iii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej) 
;     o_iii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
;     o_iii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)
;     temperature=double(10000.0)
;     density=double(5000.0)
;     atomic_levels='3,4/'
;     emiss5007=double(0.0) 
;     emiss5007=calc_emissivity(temperature=temperature, density=density, $
;                         atomic_levels=atomic_levels, $
;                         elj_data=o_iii_elj, omij_data=o_iii_omij, $
;                         aij_data=o_iii_aij
;     print, emiss5007
;
; INPUTS:
;     temperature   -     electron temperature
;     density       -     electron density
;     atomic_levels -     level(s) e.g '1,2/', '1,2,1,3/'
;     elj_data      -     energy levels (Ej) data
;     omij_data     -     collision strengths (omega_ij) data
;     aij_data      -     transition probabilities (Aij) data
;     
; RETURN:  ionic abundance
;
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
;       calc_abundance(), calc_density() and calc_temperature(), 
;                                            A. Danehkar, 20/11/2016
;     Made a new function calc_emissivity() for calculating 
;                      line emissivities and separated it 
;                      from calc_abundance(), A. Danehkar, 21/11/2016
;     Integration with AtomNeb, now uses atomic data input elj_data,
;                      omij_data, aij_data, A. Danehkar, 10/03/2017
;     Cleaning the function, and remove unused varibales
;                        calc_emissivity(), A. Danehkar, 12/06/2017
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
;- 
  common share1, Atomic_Data_Path

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s                    
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(elj_data) eq 0 then begin 
    print,'Energy Levels data (elj_data) are not set'
    return, 0
  endif
  if keyword_set(omij_data) eq 0 then begin 
    print,'Collision Strengths (omij_data) are not set'
    return, 0
  endif
  if keyword_set(aij_data) eq 0 then begin 
    print,'Transition Probabilities (aij_data) are not set'
    return, 0
  endif
  if keyword_set(atomic_levels) eq 0 then begin 
    print,'Atomic levels (atomic_levels) are not given'
    return, 0
  endif  
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif
  
  level_num= long(0) 
  temp_num= long(0) 
  IRATS= long(0) 
  ITEMP= long(0) 
  IKT= long(0) 
  
  EJI=double(0)
  WAV=double(0)
  
  temp=size(elj_data,/DIMENSIONS)
  level_num=temp[0]
  temp=size(omij_data[0].strength,/DIMENSIONS)
  temp_num=temp[0]
  temp=size(omij_data,/DIMENSIONS)
  omij_num=temp[0]
  
  Glj=lonarr(level_num)
  
  Nlj=dblarr(level_num+1)
  Omij=dblarr(temp_num,level_num,level_num)   
  Aij=dblarr(level_num+1,level_num+1)   
  Elj=dblarr(level_num)   
  Temp_List=dblarr(temp_num)
  
  LABEL1=STRARR(level_num+1) 
  
  Glj[*]=0
  levels_str=strsplit(atomic_levels, ',', ESCAPE='/', /EXTRACT)
  temp=size(levels_str, /N_ELEMENTS)
  levels_num=long(temp[0]/2)
  ITRANC=lonarr(2+1,levels_num+1)
  ITRANC[*,*]=0
  levels_i=0
  for i=1, levels_num do begin 
    res=_strnumber(levels_str[levels_i], val)
    if res eq 1 then ITRANC[1,i]=long(val)
    res=_strnumber(levels_str[levels_i+1], val)
    if res eq 1 then ITRANC[2,i]=long(val)
    levels_i = levels_i + 2
    ;if levels_i ge 2*levels_num then break
  endfor
  irats=0
  Temp_List = omij_data[0].strength
  Temp_List = alog10(Temp_List)
  for k = 1, omij_num-1 do begin
    I = omij_data[k].level1
    J = omij_data[k].level2
    if I le level_num and J le level_num then begin
      Omij[0:temp_num-1,I-1,J-1] = omij_data[k].strength
    endif
  endfor
  Aij =aij_data.AIJ
  Elj =elj_data.Ej
  Glj =long(elj_data.J_v*2.+1.)
  
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif

  Nlj=calc_populations(temperature=temperature, density=density, temp_list=temp_list, $
                       Omij=Omij, Aij=Aij, Elj=Elj, Glj=Glj, $
                       level_num=level_num, temp_num=temp_num, irats=irats)
  
  emissivity_all=double(0.0)
  for IKT=1, levels_num do begin 
    I=ITRANC[1,IKT]
    J=ITRANC[2,IKT]
    emissivity_line=double(0.0)
    if (Aij[J-1,I-1] ne 0.D0) then begin
      EJI = Elj[J-1] - Elj[I-1]
      WAV = 1.D8 / EJI
      emissivity_line=Nlj[J]*Aij[J-1,I-1]*h_Planck*c_Speed*1.e8/(WAV*density)
      emissivity_all=emissivity_all+emissivity_line
    endif
  endfor
  return, emissivity_all
end
