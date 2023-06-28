; docformat = 'rst'

;+
;     "Unit for Collisionally Excited Lines": 
;     This obejct library can be used
;     to determine electron temperature, electron density, 
;     ionic abundance from the observed flux of collisionally 
;     excited lines (CEL) for specified ion with level(s) by 
;     solving atomic level populations and line emissivities 
;     in statistical equilibrium for given electron density 
;     and temperature.
;
; :Examples:
;    For example::
;
;     IDL> s2=obj_new('collision')
;     IDL> s2->set,['s','ii']
;     IDL> 
;     IDL> upper_levels='1,2,1,3/'
;     IDL> lower_levels='1,5/'
;     IDL> density = double(2550)
;     IDL> line_flux_ratio=double(10.753)
;     IDL> temperature=s2->calc_temperature(line_flux_ratio=line_flux_ratio, density=density, $
;     IDL>   upper_levels=upper_levels, lower_levels=lower_levels)
;     IDL> print, "Electron Temperature:", temperature
;        Electron Temperature:       7920.2865
;        
;     IDL> upper_levels='1,2/'
;     IDL> lower_levels='1,3/'
;     IDL> diagtype='D'
;     IDL> temperature=double(7000.0);
;     IDL> line_flux_ratio=double(1.506);
;     IDL> density=s2->calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
;     IDL>   upper_levels=upper_levels, lower_levels=lower_levels)
;     IDL> print, "Electron Density:", density
;        Electron Density:       2312.6164
;        
;     IDL> density = double(1000)
;     IDL> temperature=double(10000.0);
;     IDL> Nlj=s2->calc_populations(temperature=temperature, density=density)
;     IDL> print, 'Atomic Level Populations:', Nlj
;        Atomic Level Populations:      0.96992796    0.0070037404     0.023062517   2.6594158e-06   3.1277593e-06
;        
;     IDL> temperature=double(10000.0)
;     IDL> N_crit=s2->calc_crit_density(temperature=temperature)
;     IDL> print, 'Critical Densities:', N_crit
;        Critical Densities:       0.0000000       5007.8396       1732.8414       1072685.0       2220758.1
;        
;     IDL> temperature=double(10000.0)
;     IDL> Omij_T=s2->get_omij_temp(temperature=temperature)
;     IDL> print, 'Effective Collision Strengths: '
;     IDL> print, Omij_T
;        Effective Collision Strengths:
;        0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        2.7800000       0.0000000       0.0000000       0.0000000       0.0000000
;        4.1600000       7.4600000       0.0000000       0.0000000       0.0000000
;        1.1700000       1.8000000       2.2000000       0.0000000       0.0000000
;        2.3500000       3.0000000       4.9900000       2.7100000       0.0000000
;        
;     IDL> s2->print_ionic, temperature=temperature, density=density
;        Temperature =   10000.0 K
;        Density =    1000.0 cm-3
;        
;        Level    Populations   Critical Densities
;        Level 1:   9.699E-01   0.000E+00
;        Level 2:   7.004E-03   5.008E+03
;        Level 3:   2.306E-02   1.733E+03
;        Level 4:   2.659E-06   1.073E+06
;        Level 5:   3.128E-06   2.221E+06
;        
;        1.231E-03
;        6732.69A
;        (2-->1)
;        2.544E-20
;        
;        3.338E-04   3.452E-07
;        6718.31A     314.47um
;        (3-->1)     (3-->2)
;        2.276E-20   5.029E-26
;        
;        1.076E-01   1.812E-01   7.506E-02
;        4077.51A       1.03um       1.04um
;        (4-->1)     (4-->2)     (4-->3)
;        1.394E-21   9.258E-22   3.823E-22
;        
;        2.670E-01   1.644E-01   1.938E-01   0.000E+00
;        4069.76A       1.03um       1.03um     214.14um
;        (5-->1)     (5-->2)     (5-->3)     (5-->4)
;        4.076E-21   9.927E-22   1.166E-21   0.000E+00
;        
;        H-beta emissivity: 1.237E-25 N(H+) Ne  [erg/s]
;        
;     IDL> o3=obj_new('collision')
;     IDL> o3->set,['o','iii']
;     IDL> 
;     IDL> levels5007='3,4/'
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> iobs5007=double(1200.0)
;     IDL> Abb5007=double(0.0)
;     IDL> 
;     IDL> emis=o3->calc_emissivity(temperature=temperature, density=density, $
;     IDL>   atomic_levels=levels5007)
;     IDL> print, 'Emissivity(O III 5007):', emis
;        Emissivity(O III 5007):   3.6041012e-21
;        
;     IDL> Abb5007=o3->calc_abundance(temperature=temperature, density=density, $
;     IDL>   line_flux=iobs5007, atomic_levels=levels5007)
;     IDL> print, 'N(O^2+)/N(H+):', Abb5007
;        N(O^2+)/N(H+):   0.00041256231
;        
;     IDL> Nlj=o3->calc_populations(temperature=temperature, density=density)
;     IDL> print, 'Atomic Level Populations:', Nlj
;        Atomic Level Populations:      0.15564960      0.42689831      0.41723001   0.00022205964   1.5224587e-08
;        
;     IDL> N_crit=o3->calc_crit_density(temperature=temperature)
;     IDL> print, 'Critical Densities:', N_crit
;        Critical Densities:       0.0000000       490.78115       3419.4864       685276.77       25472367.
;        
;     IDL> temperature=double(10000.0)
;     IDL> Omij_T=o3->get_omij_temp(temperature=temperature, level_num=8)
;     IDL> print, 'Effective Collision Strengths: '
;     IDL> print, Omij_T
;        Effective Collision Strengths:
;        0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        0.54300000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        0.27000000       1.2900000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        0.25300000      0.76000000       1.2700000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        0.032300000     0.097200000      0.16200000      0.57800000       0.0000000       0.0000000       0.0000000       0.0000000
;        0.13300000      0.39600000      0.66000000   1.9400000e-05       0.0000000       0.0000000       0.0000000       0.0000000
;        0.098800000       1.6300000      0.89000000      0.72700000    0.0029900000       1.4400000       0.0000000       0.0000000
;        0.66000000      0.62900000      0.28100000      0.29400000     0.024200000      0.46200000       1.0600000       0.0000000
;        
;     IDL> o3->print_ionic, temperature=temperature, density=density
;        Temperature =   10000.0 K
;        Density =    5000.0 cm-3
;        
;        Level    Populations   Critical Densities
;        Level 1:   1.556E-01   0.000E+00
;        Level 2:   4.269E-01   4.908E+02
;        Level 3:   4.172E-01   3.419E+03
;        Level 4:   2.221E-04   6.853E+05
;        Level 5:   1.522E-08   2.547E+07
;        
;        2.597E-05
;        88.34um
;        (2-->1)
;        4.986E-23
;        
;        0.000E+00   9.632E-05
;        32.66um      51.81um
;        (3-->1)     (3-->2)
;        0.000E+00   3.081E-22
;        
;        2.322E-06   6.791E-03   2.046E-02
;        4932.60A    4960.29A    5008.24A
;        (4-->1)     (4-->2)     (4-->3)
;        4.153E-25   1.208E-21   3.604E-21
;        
;        0.000E+00   2.255E-01   6.998E-04   1.685E+00
;        2315.58A    2321.67A    2332.12A    4364.45A
;        (5-->1)     (5-->2)     (5-->3)     (5-->4)
;        0.000E+00   5.875E-24   1.815E-26   2.335E-23
;        
;        H-beta emissivity: 1.239E-25 N(H+) Ne  [erg/s]
;        
;
; :Categories:
;   Abundance Analysis, Collisionally Excited Lines
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
;   0.2.0
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
;                        from calc_abundance(), calc_density(), and calc_temperature().
;     
;     27/02/2019, A. Danehkar, Fix a bug in the atomic level assumption, and 
;                        use the simplified calc_populations() routine.
;            
;     04/03/2019, A. Danehkar, Use the get_omij_temp() routine.
;     
;     24/05/2019, A. Danehkar, Add the optional density range to calc_density(), and
;                         the optional temperature range to calc_temperature().
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-
function collision::init
   self.data_dir = 'externals/atomneb/atomic-data/chianti70'
   self.data_rc_dir = 'externals/atomneb/atomic-data-rc'
   self.base_dir = file_dirname(file_dirname((routine_info('collision__define', /source)).path))
   self.Atom_Elj_file = filepath('AtomElj.fits', root_dir=self.base_dir, subdir=self.data_dir )
   self.Atom_Omij_file = filepath('AtomOmij.fits', root_dir=self.base_dir, subdir=self.data_dir )
   self.Atom_Aij_file = filepath('AtomAij.fits', root_dir=self.base_dir, subdir=self.data_dir )
   self.Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   self.level = 5
   return,1
end

function collision::free
  ptr_free, self.data_elj
  data_elj=ptr_new(/ALLOCATE_HEAP)
  ptr_free, self.data_omij
  data_omij=ptr_new(/ALLOCATE_HEAP)
  ptr_free, self.data_aij
  data_aij=ptr_new(/ALLOCATE_HEAP)
  ptr_free, self.hi_rc_data
  hi_rc_data=ptr_new(/ALLOCATE_HEAP)
  return,1
end

pro collision::set, atom_ion, level=level
  if n_elements(atom_ion) lt 2 then begin
    print, 'Error: atom and ionic level are not given'
    return
  endif
  if atom_ion[0] ne '' then self.atom=atom_ion[0] else print, 'Error: atom is not given'
  if atom_ion[1] ne '' then self.ion=atom_ion[1] else print, 'Error: ionic level is not given'
  ;if n_elements(atom_ion) gt 2 then if atom_ion[2] ne 0 then self.level = atom_ion[2]
  if keyword_set(level) then begin
    self.level = level
  endif
  
  data_elj=atomneb_read_elj(self.Atom_Elj_file, self.atom, self.ion, level_num=self.level) ; read Energy Levels (Ej) 
  data_omij=atomneb_read_omij(self.Atom_Omij_file, self.atom, self.ion) ; read Collision Strengths (Omegaij)
  data_aij=atomneb_read_aij(self.Atom_Aij_file, self.atom, self.ion) ; read Transition Probabilities (Aij)\ 
  hi_rc_data=atomneb_read_aeff_sh95(self.Atom_RC_SH95_file, 'h', 'ii') ; H I Rec

  self.data_elj=ptr_new(data_elj)
  self.data_omij=ptr_new(data_omij)
  self.data_aij=ptr_new(data_aij)
  self.hi_rc_data=ptr_new(hi_rc_data[0].Aeff)
  
  ;print, *(self.data_elj)
  ;print, *(self.data_omij)
  ;print, *(self.data_aij)
  ;print, *(self.hi_rc_data)
  return
end

function collision::calc_temperature, line_flux_ratio=line_flux_ratio, density=density, $
                                upper_levels=upper_levels, lower_levels=lower_levels, $
                                low_temperature=low_temperature, high_temperature=high_temperature, num_temperature=num_temperature, $
                                min_density=min_density
;+
;     This function determines electron temperature from given 
;     flux intensity ratio for specified ion with upper level(s)
;     lower level(s) by solving atomic level populations and 
;     line emissivities in statistical equilibrium 
;     for given electron density.
;
; :Returns:
;    type=double. This function returns the electron temperature.
;
; :Keywords:
;     line_flux_ratio  :     in, required, type=float
;                            flux intensity ratio
;     density          :     in, required, type=float
;                            electron density
;     upper_levels     :     in, required, type=string,
;                            upper atomic level(s) e.g '1,2/', '1,2,1,3/'
;     lower_levels     :     in, required, type=string
;                            lower atomic level(s) e.g '1,2/', '1,2,1,3/'
;     low_temperature  :     in, optional, type=float
;                            lower temperature range
;     high_temperature  :     in, optional, type=float
;                            upper temperature range
;     num_temperature  :     in, optional, type=integer
;                            number of the iteration step
;     min_density      :     in, optional, type=float
;                            lower density range
;
; :Examples:
;    For example::
;
;     IDL> s2=obj_new('collision')
;     IDL> s2->set,['s','ii']
;     IDL> 
;     IDL> upper_levels='1,2,1,3/'
;     IDL> lower_levels='1,5/'
;     IDL> density = double(2550)
;     IDL> line_flux_ratio=double(10.753)
;     IDL> temperature=s2->calc_temperature(line_flux_ratio=line_flux_ratio, density=density, $
;     IDL>   upper_levels=upper_levels, lower_levels=lower_levels)
;     IDL> print, "Electron Temperature:", temperature
;        Electron Temperature:       7920.2865
;     
;
; :Categories:
;   Plasma Diagnostics, Collisionally Excited Lines
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
;   0.2.0
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
;                        from calc_temperature().
;                        
;     27/02/2019, A. Danehkar, fix a bug in the atomic level assumption, and 
;                        use the simplified calc_populations() routine.                    
;          
;     04/03/2019, A. Danehkar, use the get_omij_temp() routine.
;     
;     24/05/2019, A. Danehkar, add the optional temperature range.
;     
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-
   value=calc_temperature(line_flux_ratio=line_flux_ratio, density=density, $
                          upper_levels=upper_levels, lower_levels=lower_levels, $
                          elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                          aij_data=*(self.data_aij), $
                          low_temperature=low_temperature, high_temperature=high_temperature, num_temperature=num_temperature, $
                          min_density=min_density)
   return, value                        
end

function collision::calc_density, line_flux_ratio=line_flux_ratio, temperature=temperature, $
                            upper_levels=upper_levels, lower_levels=lower_levels, $
                            low_density=low_density, high_density=high_density, num_density=num_density, $
                            min_temperature=min_temperature
;+
;     This function determines electron density from given
;     flux intensity ratio for specified ion with upper level(s)
;     lower level(s) by solving atomic level populations and
;     line emissivities in statistical equilibrium
;     for given electron temperature.
;
; :Returns:
;    type=double. This function returns the electron density.
;
; :Keywords:
;     line_flux_ratio  :     in, required, type=float
;                            flux intensity ratio
;     temperature      :     in, required, type=float
;                            electron temperature
;     upper_levels     :     in, required, type=string
;                            upper atomic level(s) e.g '1,2/', '1,2,1,3/'
;     lower_levels     :     in, required, type=string
;                            lower atomic level(s) e.g '1,2/', '1,2,1,3/'
;                            transition probabilities (Aij) data
;     low_density      :     in, optional, type=float
;                            lower density range
;     high_density      :     in, optional, type=float
;                            upper density range
;     num_density      :     in, optional, type=integer
;                            number of the iteration step
;     min_temperature  :     in, optional, type=float
;                            minimum temperature
;
; :Examples:
;    For example::
;
;     IDL> s2=obj_new('collision')
;     IDL> s2->set,['s','ii']
;     IDL> 
;     IDL> upper_levels='1,2/'
;     IDL> lower_levels='1,3/'
;     IDL> diagtype='D'
;     IDL> temperature=double(7000.0);
;     IDL> line_flux_ratio=double(1.506);
;     IDL> density=s2->calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
;     IDL>   upper_levels=upper_levels, lower_levels=lower_levels)
;     IDL> print, "Electron Density:", density
;        Electron Density:       2312.6164
;     
;
; :Categories:
;   Plasma Diagnostics, Collisionally Excited Lines
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
;   0.2.0
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
;                        from calc_density().
;
;     27/02/2019, A. Danehkar, fix a bug in the atomic level assumption, and
;                        use the simplified calc_populations() routine.
;
;     04/03/2019, A. Danehkar, use the get_omij_temp() routine.
;
;     24/05/2019, A. Danehkar, add the optional density range.
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-
    value=calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
                       upper_levels=upper_levels, lower_levels=lower_levels, $
                       elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                       aij_data=*(self.data_aij), $
                       low_density=low_density, high_density=high_density, num_density=num_density, $
                       min_temperature=min_temperature)
    return, value
end

function collision::calc_populations, temperature=temperature, density=density, $
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
;     IDL> s2=obj_new('collision')
;     IDL> s2->set,['s','ii']
;     IDL> 
;     IDL> density = double(1000)
;     IDL> temperature=double(10000.0);
;     IDL> Nlj=s2->calc_populations(temperature=temperature, density=density)
;     IDL> print, 'Atomic Level Populations:', Nlj
;        Atomic Level Populations:      0.96992796    0.0070037404     0.023062517   2.6594158e-06   3.1277593e-06
;     
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
;     27/02/2019, A. Danehkar, Simplify the calc_populations() routine 
;                        for external usage.
;                            
;     04/03/2019, A. Danehkar, Use the get_omij_temp() routine.
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-

  value=calc_populations(temperature=temperature, density=density, $
                         elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                         aij_data=*(self.data_aij), $
                         eff_Omij=eff_Omij, $
                         level_num=level_num, irats=irats)
  return, value
end

function collision::calc_crit_density, temperature=temperature, $
                                 level_num=level_num, irats=irats
;+
;     This function calculates critical densities in statistical equilibrium
;     for given electron temperature.
;
; :Returns:
;    type=array/object. This function returns the critical densities.
;
; :Keywords:
;     temperature :   in, required, type=float
;                     electron temperature
;     level_num   :   in, type=int
;                     Number of levels
;     irats       :   in, type=int
;                     Else Coll. rates = tabulated values * 10 ** irats
;
; :Examples:
;    For example::
;
;     IDL> s2=obj_new('collision')
;     IDL> s2->set,['s','ii']
;     IDL> 
;     IDL> temperature=double(10000.0)
;     IDL> N_crit=s2->calc_crit_density(temperature=temperature)
;     IDL> print, 'Critical Densities:', N_crit
;        Critical Densities:       0.0000000       5007.8396       1732.8414       1072685.0       2220758.1
;
; :Categories:
;   Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
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
;   0.2.0
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
;     27/02/2019, A. Danehkar, simplify the calc_populations() routine
;                        for external usage.
;
;     01/03/2019, A. Danehkar, create the calc_crit_density() routine
;                        from the calc_populations() routine.
;
;     04/03/2019, A. Danehkar, use the get_omij_temp() routine.
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-
  value=calc_crit_density(temperature=temperature, $
                         elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                         aij_data=*(self.data_aij), $
                         level_num=level_num, irats=irats)
  return, value
end

function collision::calc_emissivity, temperature=temperature, density=density, $
                               atomic_levels=atomic_levels
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
;
; :Examples:
;    For example::
;
;     IDL> o3=obj_new('collision')
;     IDL> o3->set,['o','iii']
;     IDL> 
;     IDL> levels5007='3,4/'
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> iobs5007=double(1200.0)
;     IDL> Abb5007=double(0.0)
;     IDL> 
;     IDL> emis=o3->calc_emissivity(temperature=temperature, density=density, $
;     IDL>   atomic_levels=levels5007)
;     IDL> print, 'Emissivity(O III 5007):', emis
;        Emissivity(O III 5007):   3.6041012e-21
;     
;
; :Categories:
;   Abundance Analysis, Collisionally Excited Lines
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
;   0.2.0
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
;     27/06/2019, A. Danehkar, use the simplified calc_populations() routine.
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-
  value=calc_emissivity(temperature=temperature, density=density, $
                        atomic_levels=atomic_levels, $
                        elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                        aij_data=*(self.data_aij))
  return, value
end

function collision::calc_abundance, temperature=temperature, density=density, $
                              line_flux=line_flux, atomic_levels=atomic_levels
;+
;     This function determines the ionic abundance from the observed
;     flux intensity for specified ion with level(s)
;     by solving atomic level populations and
;     line emissivities in statistical equilibrium
;     for given electron density and temperature.
;
; :Returns:
;    type=double. This function returns the ionic abundanc.
;
; :Keywords:
;     temperature   :     in, required, type=float
;                         electron temperature
;     density       :     in, required, type=float
;                         electron density
;     line_flux     :     in, required, type=float
;                         line flux intensity
;     atomic_levels :     in, required, type=string
;                         level(s) e.g '1,2/', '1,2,1,3/'
;
; :Examples:
;    For example::
;
;     IDL> o3=obj_new('collision')
;     IDL> o3->set,['o','iii']
;     IDL> 
;     IDL> levels5007='3,4/'
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> iobs5007=double(1200.0)
;     IDL> Abb5007=double(0.0)
;     IDL> 
;     IDL> Abb5007=o3->calc_abundance(temperature=temperature, density=density, $
;     IDL>   line_flux=iobs5007, atomic_levels=levels5007)
;     IDL> print, 'N(O^2+)/N(H+):', Abb5007
;        N(O^2+)/N(H+):   0.00041256231
;     
;
; :Categories:
;   Abundance Analysis, Collisionally Excited Lines
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
;   0.2.0
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
;                        from calc_abundance().
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
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
;
;     2009, R.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-

  value=calc_abundance(temperature=temperature, density=density, $
                       line_flux=line_flux, atomic_levels=atomic_levels, $
                       elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                       aij_data=*(self.data_aij), h_i_aeff_data=*(self.hi_rc_data))
  return, value
end
  

pro collision::print_ionic, temperature=temperature, density=density, $
                           printEmissivity=printEmissivity, $
                           printPopulations=printPopulations, $
                           printCritDensity=printCritDensity
;+
;    This function prints the atom's transitions information,
;    atomic level populations, critical densities, and emissivities
;    for given temperature and density.
;
; :Keywords:
;     temperature   :   in, required, type=float
;                       electron temperature
;     density       :   in, required, type=float
;                       electron density
;     printEmissivity  :   in, type=boolean
;                          Set for printing Emissivities
;     printPopulations :   in, type=boolean
;                          Set for printing Populations
;     printCritDensity  :  in, type=boolean
;                          Set for printing Critical Densities
;
; :Examples:
;    For example::
;
;     IDL> o3=obj_new('collision')
;     IDL> o3->set,['o','iii']
;     IDL>
;     IDL> levels5007='3,4/'
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> iobs5007=double(1200.0)
;     IDL> Abb5007=double(0.0)
;     IDL> 
;     IDL> o3->print_ionic, temperature=temperature, density=density
;        Temperature =   10000.0 K
;        Density =    5000.0 cm-3
;
;        Level    Populations   Critical Densities
;        Level 1:   1.556E-01   0.000E+00
;        Level 2:   4.269E-01   4.908E+02
;        Level 3:   4.172E-01   3.419E+03
;        Level 4:   2.221E-04   6.853E+05
;        Level 5:   1.522E-08   2.547E+07
;
;        2.597E-05
;        88.34um
;        (2-->1)
;        4.986E-23
;
;        0.000E+00   9.632E-05
;        32.66um      51.81um
;        (3-->1)     (3-->2)
;        0.000E+00   3.081E-22
;
;        2.322E-06   6.791E-03   2.046E-02
;        4932.60A    4960.29A    5008.24A
;        (4-->1)     (4-->2)     (4-->3)
;        4.153E-25   1.208E-21   3.604E-21
;
;        0.000E+00   2.255E-01   6.998E-04   1.685E+00
;        2315.58A    2321.67A    2332.12A    4364.45A
;        (5-->1)     (5-->2)     (5-->3)     (5-->4)
;        0.000E+00   5.875E-24   1.815E-26   2.335E-23
;
;        H-beta emissivity: 1.239E-25 N(H+) Ne  [erg/s]
;
;
; :Categories:
;   Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
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
;   0.2.0
;
; :History:
;     04/03/2019, A. Danehkar, create the print_ionic() routine.
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
;-
    print_ionic, temperature=temperature, density=density, $
                  elj_data=*(self.data_elj), omij_data=*(self.data_omij), $
                  aij_data=*(self.data_aij), h_i_aeff_data=*(self.hi_rc_data), $
                  printEmissivity=printEmissivity, $
                  printPopulations=printPopulations, $
                  printCritDensity=printCritDensity
end

function collision::get_omij_temp, temperature=temperature, $
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
;     level_num   :   in, type=int
;                     Number of levels
;     irats       :   in, type=int
;                     Else Coll. rates = tabulated values * 10 ** irats
;
; :Examples:
;    For example::
;
;     IDL> s2=obj_new('collision')
;     IDL> s2->set,['s','ii']
;     IDL>
;     IDL> temperature=double(10000.0)
;     IDL> Omij_T=s2->get_omij_temp(temperature=temperature)
;     IDL> print, 'Effective Collision Strengths: '
;     IDL> print, Omij_T
;        Effective Collision Strengths:
;        0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
;        2.7800000       0.0000000       0.0000000       0.0000000       0.0000000
;        4.1600000       7.4600000       0.0000000       0.0000000       0.0000000
;        1.1700000       1.8000000       2.2000000       0.0000000       0.0000000
;        2.3500000       3.0000000       4.9900000       2.7100000       0.0000000
;     
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
;     04/03/2019, A. Danehkar, create the get_omij_temp() routine.
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
;-

  value=get_omij_temp(temperature=temperature, $
                      omij_data=*(self.data_omij), elj_data=*(self.data_elj), $
                      level_num=level_num, irats=irats)
  return, value
end
                            
;-------------
pro collision::set_level_num, level_num
  if level_num ne '' then self.level_num=level_num else print, 'Error: level_num is not given'
  return
end

function collision::get_level_num
  if self.level_num ne '' then level_num=self.level_num else print, 'Error: level_num is not given'
  return, level_num
end
;-------------
pro collision::set_data_dir, data_dir
  if data_dir ne '' then self.data_dir=data_dir else print, 'Error: data_dir is not given'
  return
end

function collision::get_data_dir
  if self.data_dir ne '' then data_dir=self.data_dir else print, 'Error: data_dir is not given'
  return, data_dir
end
;-------------
pro collision::set_data_rc_dir, data_rc_dir
  if data_dir ne '' then self.data_dir=data_dir else print, 'Error: data_rc_dir is not given'
  return
end

function collision::get_data_rc_dir
  if self.data_rc_dir ne '' then data_rc_dir=self.data_rc_dir else print, 'Error: data_rc_dir is not given'
  return, data_rc_dir
end
;-------------
pro collision::set_Atom_Elj_file, Atom_Elj_file
  if Atom_Elj_file ne '' then self.Atom_Elj_file=Atom_Elj_file else print, 'Error: Atom_Elj_file is not given'
  return
end

function collision::get_Atom_Elj_file
  if self.Atom_Elj_file ne '' then Atom_Elj_file=self.Atom_Elj_file else print, 'Error: Atom_Elj_file is not given'
  return, Atom_Elj_file
end
;-------------
pro collision::set_Atom_Omij_file, Atom_Omij_file
  if Atom_Omij_file ne '' then self.Atom_Omij_file=Atom_Omij_file else print, 'Error: Atom_Omij_file is not given'
  return
end

function collision::get_Atom_Omij_file
  if self.Atom_Omij_file ne '' then Atom_Omij_file=self.Atom_Omij_file else print, 'Error: Atom_Omij_file is not given'
  return, Atom_Omij_file
end
;-------------
pro collision::set_Atom_Aij_file, Atom_Aij_file
  if Atom_Aij_file ne '' then self.Atom_Aij_file=Atom_Aij_file else print, 'Error: Atom_Aij_file is not given'
  return
end

function collision::get_Atom_Aij_file
  if self.Atom_Aij_file ne '' then Atom_Aij_file=self.Atom_Aij_file else print, 'Error: Atom_Aij_file is not given'
  return, Atom_Aij_file
end
;-------------
pro collision::set_Atom_RC_SH95_file, Atom_RC_SH95_file
  if Atom_RC_SH95_file ne '' then self.Atom_RC_SH95_file=Atom_RC_SH95_file else print, 'Error: Atom_RC_SH95_file is not given'
  return
end

function collision::get_Atom_RC_SH95_file
  if self.Atom_RC_SH95_file ne '' then Atom_RC_SH95_file=self.Atom_RC_SH95_file else print, 'Error: Atom_RC_SH95_file is not given'
  return, Atom_RC_SH95_file
end
;------------------------------------------------------------------

pro collision__define
  void={collision, level:0L, data_dir:'', data_rc_dir:'', $
          Atom_Elj_file:'', Atom_Omij_file:'', Atom_Aij_file:'', Atom_RC_SH95_file:'', $
          data_elj:ptr_new(/ALLOCATE_HEAP), data_omij:ptr_new(/ALLOCATE_HEAP), $
          data_aij:ptr_new(/ALLOCATE_HEAP), hi_rc_data:ptr_new(/ALLOCATE_HEAP), $
          inherits ion_unit}
  return 
end
