; Example: calc_temperature() and calc_density()
;     determine electron density or temperature from given 
;     flux intensity ratio for a fixed electron density 
;     or temperature using calc_temperature function
;     calc_density function from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 
Atom_Elj_file='/AtomNeb/atomic-data/chianti70/AtomElj.fits'
Atom_Omij_file='/AtomNeb/atomic-data/chianti70/AtomOmij.fits'
Atom_Aij_file='/AtomNeb/atomic-data/chianti70/AtomAij.fits'

atom='s'
ion='ii'
s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej) 
s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)\ 

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

end 
