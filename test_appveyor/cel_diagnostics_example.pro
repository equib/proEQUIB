; Example: calc_temperature() and calc_density()
;     determine electron density or temperature from given 
;     flux intensity ratio for a fixed electron density 
;     or temperature using calc_temperature function
;     calc_density function from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 

; Locate datasets
;base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
base_dir='/home/appveyor/projects/proequib/externals/'
;data_dir = ['atomic-data', 'chianti70']
data_dir = ['atomneb', 'atomic-data', 'chianti70']
Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )

atom='s'
ion='ii'
s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej) 
s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)

upper_levels='1,2,1,3/'   
lower_levels='1,5/'
density = double(2550)
line_flux_ratio=double(10.753)
temperature=calc_temperature(line_flux_ratio=line_flux_ratio, density=density, $
                          upper_levels=upper_levels, lower_levels=lower_levels, $
                          elj_data=s_ii_elj, omij_data=s_ii_omij, $
                          aij_data=s_ii_aij)
print, "Electron Temperature:", temperature

upper_levels='1,2/'   
lower_levels='1,3/'
diagtype='D'
temperature=double(7000.0);
line_flux_ratio=double(1.506);
density=calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
                     upper_levels=upper_levels, lower_levels=lower_levels, $
                     elj_data=s_ii_elj, omij_data=s_ii_omij, $
                     aij_data=s_ii_aij)
print, "Electron Density:", density

density = double(1000)
temperature=double(10000.0);
Nlj=calc_populations(temperature=temperature, density=density, $
                     elj_data=s_ii_elj, omij_data=s_ii_omij, $
                     aij_data=s_ii_aij)
print, 'Atomic Level Populations:', Nlj

temperature=double(10000.0)
N_crit=calc_crit_density(temperature=temperature, $
                         elj_data=s_ii_elj, omij_data=s_ii_omij, $
                         aij_data=s_ii_aij)
print, 'Critical Densities:', N_crit

temperature=double(10000.0)
Omij_T=get_omij_temp(temperature=temperature, omij_data=s_ii_omij)
print, 'Effective Collision Strengths: '
print, Omij_T

; --- End $MAIN$ program. ---------------
exit

