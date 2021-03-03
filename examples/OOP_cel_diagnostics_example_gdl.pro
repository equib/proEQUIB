; Example: calc_temperature() and calc_density()
;     determine electron density or temperature from given 
;     flux intensity ratio for a fixed electron density 
;     or temperature using calc_temperature function
;     calc_density function from collision object in 
;     proEQUIB
;
; Example of object-oriented programming (OOP) for
;     collision object
;
; --- Begin $MAIN$ program. ---------------
; 
; 

s2=obj_new('collision')
s2->set,['s','ii']
print, s2->get()
print, s2->get_base_dir()
print, s2->get_data_dir()
print, s2->get_data_rc_dir()
print, s2->get_Atom_Elj_file()
print, s2->get_Atom_Omij_file()
print, s2->get_Atom_Aij_file()
print, s2->get_Atom_RC_SH95_file()

upper_levels='1,2,1,3/'
lower_levels='1,5/'
density = double(2550)
line_flux_ratio=double(10.753)
temperature=s2->calc_temperature(line_flux_ratio=line_flux_ratio, density=density, $
                                upper_levels=upper_levels, lower_levels=lower_levels)
print, "Electron Temperature:", temperature

upper_levels='1,2/'
lower_levels='1,3/'
temperature=double(7000.0);
line_flux_ratio=double(1.506);
density=s2->calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
                      upper_levels=upper_levels, lower_levels=lower_levels)
print, "Electron Density:", density

density = double(1000)
temperature=double(10000.0);
Nlj=s2->calc_populations(temperature=temperature, density=density)
print, 'Atomic Level Populations:', Nlj

temperature=double(10000.0)
N_crit=s2->calc_crit_density(temperature=temperature)
print, 'Critical Densities:', N_crit

temperature=double(10000.0)
Omij_T=s2->get_omij_temp(temperature=temperature)
print, 'Effective Collision Strengths: '
print, Omij_T

s2->print_ionic, temperature=temperature, density=density
  
;exit
exit
