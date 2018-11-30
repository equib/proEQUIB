; Example: he_i_emissivity_porter(), he_ii_emissivity()
;          recomb_c_ii(), recomb_c_iii()
;          recomb_n_ii(), recomb_n_iii()
;          recomb_o_ii(), recomb_ne_ii()
;     determine ionic abundance from observed 
;     flux intensity for gievn electron density 
;     and temperature using he_i_emissivity_porter, 
;     he_ii_emissivity, recomb_c_ii, recomb_c_iii
;     recomb_o_ii, and recomb_ne_ii from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 

; Locate datasets
;base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
base_dir='/home/appveyor/projects/proequib/externals/'
;data_dir = ['atomic-data-rc']
data_rc_dir = ['atomneb', 'atomic-data-rc']
Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
Atom_RC_PPB91_file= filepath('rc_PPB91.fits', root_dir=base_dir, subdir=data_rc_dir )
Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )

atom='h'
ion='ii' ; H I
h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)

atom='he'
ion='ii' ; He I
he_i_rc_data=atomneb_read_aeff_he_i_pfsd12(Atom_RC_He_I_file, atom, ion)

atom='he'
ion='iii' ; He II
he_ii_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)

atom='c'
ion='iii' ; C II
c_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)

atom='c'
ion='iv' ; C III
c_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)

atom='n'
ion='iii' ; N II
n_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
n_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)

atom='n'
ion='iv' ; N III
n_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)

atom='o'
ion='iii' ; O II
o_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
o_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)

atom='ne'
ion='iii' ; Ne II
ne_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)

h_i_aeff_data=h_i_rc_data[0].Aeff
he_i_aeff_data=he_i_rc_data[0].Aeff
he_ii_aeff_data=he_ii_rc_data[0].Aeff

temperature=double(10000.0)
density=double(5000.0)

; 4120.84: linenum=7
; 4387.93: linenum=8
; 4437.55: linenum=9
; 4471.50: linenum=10
; 4921.93: linenum=12
; 5015.68: linenum=13
; 5047.74: linenum=14
; 5875.66: linenum=15
; 6678.16: linenum=16
; 7065.25: linenum=17
; 7281.35: linenum=18
he_i_4471_flux= 2.104
linenum=10; 4471.50
Abund_he_i=calc_abund_he_i_rl(temperature=temperature, density=density, $
                              linenum=linenum, line_flux=he_i_4471_flux, $
                              he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
print, 'N(He^+)/N(H^+):', Abund_he_i

he_ii_4686_flux = 135.833
Abund_he_ii=calc_abund_he_ii_rl(temperature=temperature, density=density, $
                                line_flux=he_ii_4686_flux, $
                                he_ii_aeff_data=he_ii_aeff_data, h_i_aeff_data=h_i_aeff_data)
print, 'N(He^2+)/N(H^+):', Abund_he_ii

c_ii_6151_flux = 0.028
wavelength=6151.43
Abund_c_ii=calc_abund_c_ii_rl(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=c_ii_6151_flux, $
                              c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
print, 'N(C^2+)/N(H+):', Abund_c_ii

c_iii_4647_flux = 0.107
wavelength=4647.42
Abund_c_iii=calc_abund_c_iii_rl(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=c_iii_4647_flux, $
                                c_iii_rc_data=c_iii_rc_data, h_i_aeff_data=h_i_aeff_data) 
print, 'N(C^3+)/N(H+):', Abund_c_iii

n_ii_4442_flux = 0.017
wavelength=4442.02
Abund_n_ii=calc_abund_n_ii_rl(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=n_ii_4442_flux, $
                              n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data, $
                              h_i_aeff_data=h_i_aeff_data)
print, 'N(N^2+)/N(H+):', Abund_n_ii

n_iii_4641_flux = 0.245
wavelength=4640.64
Abund_n_iii=calc_abund_n_iii_rl(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=n_iii_4641_flux, $
                                n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
print, 'N(N^3+)/N(H+):', Abund_n_iii

o_ii_4614_flux = 0.009
wavelength=4613.68
Abund_o_ii=calc_abund_o_ii_rl(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=o_ii_4614_flux, $
                              o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data, $ 
                              h_i_aeff_data=h_i_aeff_data)                      
print, 'N(O^2+)/N(H+):', Abund_o_ii

ne_ii_3777_flux = 0.056
wavelength=3777.14
Abund_ne_ii=calc_abund_ne_ii_rl(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=ne_ii_3777_flux, $
                                ne_ii_rc_data=ne_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
print, 'N(Ne^2+)/N(H+):', Abund_ne_ii

end 
