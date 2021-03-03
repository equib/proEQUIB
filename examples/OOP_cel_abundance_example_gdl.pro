; Example: calc_abundance()
;     determine ionic abundance from observed 
;     flux intensity for gievn electron density 
;     and temperature using  calc_abundance function
;     from collision object in proEQUIB.
;
; Example of object-oriented programming (OOP) for
;     collision object
; 
; --- Begin $MAIN$ program. ---------------
; 
; 

o3=obj_new('collision')
o3->set,['o','iii']
print, o3->get()
print, o3->get_base_dir()
print, o3->get_data_dir()
print, o3->get_data_rc_dir()
print, o3->get_Atom_Elj_file()
print, o3->get_Atom_Omij_file()
print, o3->get_Atom_Aij_file()
print, o3->get_Atom_RC_SH95_file()


levels5007='3,4/'
temperature=double(10000.0)
density=double(5000.0)
iobs5007=double(1200.0)
Abb5007=double(0.0)

emis=o3->calc_emissivity(temperature=temperature, density=density, $
                    atomic_levels=levels5007)
print, 'Emissivity(O III 5007):', emis

Abb5007=o3->calc_abundance(temperature=temperature, density=density, $
                      line_flux=iobs5007, atomic_levels=levels5007)
print, 'N(O^2+)/N(H+):', Abb5007

Nlj=o3->calc_populations(temperature=temperature, density=density)
print, 'Atomic Level Populations:', Nlj

N_crit=o3->calc_crit_density(temperature=temperature)
print, 'Critical Densities:', N_crit

temperature=double(10000.0)
Omij_T=o3->get_omij_temp(temperature=temperature, level_num=8)
print, 'Effective Collision Strengths: '
print, Omij_T

o3->print_ionic, temperature=temperature, density=density

;exit
exit
