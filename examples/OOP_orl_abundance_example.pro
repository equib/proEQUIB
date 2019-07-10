; Example: calc_abund_he_i_rl(), calc_abund_he_ii_rl()
;          calc_abund_c_ii_rl(), calc_abund_c_iii_rl()
;          calc_abund_n_ii_rl(), calc_abund_n_iii_rl()
;          calc_abund_o_ii_rl(), calc_abund_ne_ii_rl()
;     determine ionic abundance from observed 
;     flux intensity for gievn electron density 
;     and temperature using calc_abund_he_i_rl, 
;     calc_abund_he_ii_rl, calc_abund_c_ii_rl, calc_abund_c_iii_rl
;     calc_abund_n_ii_rl, calc_abund_n_iii_rl
;     calc_abund_o_ii_rl, and calc_abund_ne_ii_rl from 
;     recombination object in proEQUIB
;
; Example of object-oriented programming (OOP) for
;     recombination object
; 
; --- Begin $MAIN$ program. ---------------
; 
; 

he1=obj_new('recombination')
he1->set,['he','ii'] ; He I

he2=obj_new('recombination')
he2->set,['he','iii'] ; He II

c2=obj_new('recombination')
c2->set,['c','iii'] ; C II

c3=obj_new('recombination')
c3->set,['c','iv'] ; C III

n2=obj_new('recombination')
n2->set,['n','iii'] ; N II

n3=obj_new('recombination')
n3->set,['n','iv'] ; N III

o2=obj_new('recombination')
o2->set,['o','iii'] ; O II

ne2=obj_new('recombination')
ne2->set,['ne','iii'] ; Ne II

print, he1->get()
print, he1->get_base_dir()
print, he1->get_data_rc_dir()
print, he1->get_Atom_RC_He_I_file()
print, he1->get_Atom_RC_PPB91_file()
print, he1->get_Atom_RC_SH95_file()

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
linenum=10; 4471.50
emiss_he_i=he1->calc_emissivity(temperature=temperature, density=density, linenum=linenum)
print, 'He I Emissivity:', emiss_he_i
he_i_4471_flux= 2.104
Abund_he_i=he1->calc_abundance(temperature=temperature, density=density, $
                              linenum=linenum, line_flux=he_i_4471_flux)
print, 'N(He^+)/N(H^+):', Abund_he_i

emiss_he_ii=he2->calc_emissivity(temperature=temperature, density=density)
print, 'He II Emissivity:', emiss_he_ii
he_ii_4686_flux = 135.833
Abund_he_ii=he2->calc_abundance(temperature=temperature, density=density, $
                                line_flux=he_ii_4686_flux)
print, 'N(He^2+)/N(H^+):', Abund_he_ii


wavelength=6151.43
emiss_c_ii=c2->calc_emissivity(temperature=temperature, density=density, $
                                  wavelength=wavelength)
print, 'C II Emissivity:', emiss_c_ii
c_ii_6151_flux = 0.028
Abund_c_ii=c2->calc_abundance(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=c_ii_6151_flux)
print, 'N(C^2+)/N(H+):', Abund_c_ii


wavelength=4647.42
emiss_c_iii=c3->calc_emissivity(temperature=temperature, density=density, $
                                    wavelength=wavelength)
print, 'C III Emissivity:', emiss_c_iii
c_iii_4647_flux = 0.107
Abund_c_iii=c3->calc_abundance(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=c_iii_4647_flux) 
print, 'N(C^3+)/N(H+):', Abund_c_iii

wavelength=4442.02
emiss_n_ii=n2->calc_emissivity(temperature=temperature, density=density, $
                                 wavelength=wavelength)
print, 'N II Emissivity:', emiss_n_ii
n_ii_4442_flux = 0.017
Abund_n_ii=n2->calc_abundance(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=n_ii_4442_flux)
print, 'N(N^2+)/N(H+):', Abund_n_ii

wavelength=4640.64
emiss_n_iii=n3->calc_emissivity(temperature=temperature, density=density, $
                                wavelength=wavelength)
print, 'N III Emissivity:', emiss_n_iii
n_iii_4641_flux = 0.245
Abund_n_iii=n3->calc_abundance(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=n_iii_4641_flux)
print, 'N(N^3+)/N(H+):', Abund_n_iii

wavelength=4613.68
emiss_o_ii=o2->calc_emissivity(temperature=temperature, density=density, $
                                  wavelength=wavelength)
print, 'O II Emissivity:', emiss_o_ii
o_ii_4614_flux = 0.009
Abund_o_ii=o2->calc_abundance(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=o_ii_4614_flux)                      
print, 'N(O^2+)/N(H+):', Abund_o_ii

wavelength=3777.14
emiss_ne_ii=ne2->calc_emissivity(temperature=temperature, density=density, $
                                wavelength=wavelength)
print, 'Ne II Emissivity:', emiss_ne_ii
ne_ii_3777_flux = 0.056
Abund_ne_ii=ne2->calc_abundance(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=ne_ii_3777_flux)
print, 'N(Ne^2+)/N(H+):', Abund_ne_ii

; --- End $MAIN$ program. ---------------
;exit
end

