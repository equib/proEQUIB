; Example: calc_abundance()
;     determine ionic abundance from observed 
;     flux intensity for gievn electron density 
;     and temperature using  calc_abundance function
;     from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 

; Locate datasets
base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
data_dir = ['externals','atomneb','atomic-data','chianti70']
data_rc_dir = ['externals','atomneb','atomic-data-rc']
Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )
Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )

atom='h'
ion='ii' ; H I Rec
hi_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)

atom='o'
ion='iii' ; [O III]
o_iii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej) 
o_iii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
o_iii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)

levels5007='3,4/'
temperature=double(10000.0)
density=double(5000.0)
iobs5007=double(1200.0)
Abb5007=double(0.0) 

emis=calc_emissivity(temperature=temperature, density=density, $
                         atomic_levels=levels5007, $
                         elj_data=o_iii_elj, omij_data=o_iii_omij, $
                         aij_data=o_iii_aij)
print, 'Emissivity(O III 5007):', emis
                      
Abb5007=calc_abundance(temperature=temperature, density=density, $
                       line_flux=iobs5007, atomic_levels=levels5007,$
                       elj_data=o_iii_elj, omij_data=o_iii_omij, $
                       aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data[0].Aeff)
print, 'N(O^2+)/N(H+):', Abb5007

Nlj=calc_populations(temperature=temperature, density=density, $
                     elj_data=o_iii_elj, omij_data=o_iii_omij, $
                     aij_data=o_iii_aij)
print, 'Atomic Level Populations:', Nlj


N_crit=calc_crit_density(temperature=temperature, $
                         elj_data=o_iii_elj, omij_data=o_iii_omij, $
                         aij_data=o_iii_aij)
print, 'Critical Densities:', N_crit

temperature=double(10000.0)
Omij_T=get_omij_temp(temperature=temperature, omij_data=o_iii_omij, level_num=8)
print, 'Effective Collision Strengths: '
print, Omij_T

print_ionic, temperature=temperature, density=density, $, $
             elj_data=o_iii_elj, omij_data=o_iii_omij, $
             aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data[0].Aeff
;     
; --- End $MAIN$ program. ---------------
;exit
end
