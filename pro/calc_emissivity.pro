; docformat = 'rst'

function calc_emissivity, temperature=temperature, density=density, $
                         atomic_levels=atomic_levels, $
                         elj_data=elj_data, omij_data=omij_data, $
                         aij_data=aij_data
;+
;     This function calculates line emissivities for specified ion with level(s) by 
;     solving atomic level populations and in statistical equilibrium 
;     for given electron density and temperature.
;
; :Returns:
;    type=double. This function returns the line emissivity.
;
; :Keywords:
;     temperature   :     in, required, type=float
;                         electron temperature
;     density       :     in, required, type=float
;                         electron density
;     atomic_levels :     In, required, type=string
;                         level(s) e.g '1,2/', '1,2,1,3/'
;     elj_data      :     in, required, type=array/object
;                         energy levels (Ej) data
;     omij_data     :     in, required, type=array/object
;                         collision strengths (omega_ij) data
;     aij_data      :     in, required, type=array/object
;                         transition probabilities (Aij) data
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_dir = ['atomic-data', 'chianti70']
;     IDL> Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> atom='o'
;     IDL> ion='iii'
;     IDL> o_iii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej)
;     IDL> o_iii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
;     IDL> o_iii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> atomic_levels='3,4/'
;     IDL> emiss5007=double(0.0)
;     IDL> emiss5007=calc_emissivity(temperature=temperature, density=density, $
;     IDL>                           atomic_levels=atomic_levels, $
;     IDL>                           elj_data=o_iii_elj, omij_data=o_iii_omij, $
;     IDL>                           aij_data=o_iii_aij
;     IDL> print, 'Emissivity(O III 5007):', emiss5007
;        Emissivity(O III 5007):   3.6041012e-21
;
; :Categories:
;   Abundance Analysis, Collisionally Excited Lines, Emissivity
;
; :Dirs:
;  ./
;      Main routines
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
;     21/11/2016, A. Danehkar, Made a new function calc_emissivity()
;                      for calculating line emissivities and separated it
;                      from calc_abundance().
;
;     10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
;                      input elj_data, omij_data, aij_data.
;
;     12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
;                        from calc_emissivity().
;                        
;     27/06/2019, A. Danehkar, Use the simplified calc_populations() routine.
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
;     calc_emissivity
; PURPOSE:
;     This function calculates line emissivities for specified ion with level(s) by 
;     solving atomic level populations and in statistical equilibrium 
;     for given electron density and temperature.
;
; CALLING SEQUENCE:
;     Result = calc_emissivity(TEMPERATURE=temperature, DENSITY=density, $
;                              ATOMIC_LEVELS=atomic_levels, $
;                              ELJ_DATA=elj_data, OMIJ_DATA=omij_data, $
;                              AIJ_DATA=aij_data)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE   :     in, required, type=float, electron temperature
;     DENSITY       :     in, required, type=float, electron density
;     ATOMIC_LEVELS :     In, required, type=string, level(s) e.g '1,2/', '1,2,1,3/'
;     ELJ_DATA      :     in, required, type=array/object, energy levels (Ej) data
;     OMIJ_DATA     :     in, required, type=array/object, collision strengths (omega_ij) data
;     AIJ_DATA      :     in, required, type=array/object, transition probabilities (Aij) data
;     
; OUTPUTS:  This function returns a double as the line emissivity.
;
; PROCEDURE: This function is called by calc_abundance. 
;            This function calls calc_populations.
;
; EXAMPLE:
;     base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     data_dir = ['atomic-data', 'chianti70']
;     Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
;     Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
;     Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )
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
;                               atomic_levels=atomic_levels, $
;                               elj_data=o_iii_elj, omij_data=o_iii_omij, $
;                               aij_data=o_iii_aij
;     print, 'Emissivity(O III 5007):', emiss5007
;     > Emissivity(O III 5007):   3.6041012e-21
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
;     21/11/2016, A. Danehkar, Made a new function calc_emissivity()
;                      for calculating line emissivities and separated it
;                      from calc_abundance().
;     10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
;                      input elj_data, omij_data, aij_data.
;     12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
;                        from calc_emissivity().
;     27/06/2019, A. Danehkar, use the simplified calc_populations() routine.
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
  
  IRATS= long(0) 
  ITEMP= long(0) 
  IKT= long(0) 
  
  EJI=double(0)
  WAV=double(0)
  
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
  Aij =aij_data.AIJ
  Elj =elj_data.Ej
  
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif
  
  Nlj=calc_populations(temperature=temperature, density=density, $
                       elj_data=elj_data, omij_data=omij_data, $
                       aij_data=aij_data, $
                       irats=irats)
  
  emissivity_all=double(0.0)
  for IKT=1, levels_num do begin 
    I=ITRANC[1,IKT]
    J=ITRANC[2,IKT]
    emissivity_line=double(0.0)
    if (Aij[J-1,I-1] ne 0.D0) then begin
      EJI = Elj[J-1] - Elj[I-1]
      WAV = 1.D8 / EJI
      emissivity_line=Nlj[J-1]*Aij[J-1,I-1]*h_Planck*c_Speed*1.e8/(WAV*density)
      emissivity_all=emissivity_all+emissivity_line
    endif
  endfor
  return, emissivity_all
end
