; --- Begin $MAIN$ program. ---------------
; 
; 

mgunit, ['redlaw_ut', 'redlaw_gal_ut', 'redlaw_gal2_ut', 'redlaw_ccm_ut', $
        'redlaw_jbk_ut', 'redlaw_fm_ut', 'redlaw_smc_ut', 'redlaw_lmc_ut', $
        'deredden_relflux_ut', 'deredden_flux_ut', $
        'calc_abundance_ut', 'calc_emissivity_ut', $
        'calc_temperature_ut', 'calc_density_ut', $
        'calc_abund_he_i_rl_ut', 'calc_abund_he_ii_rl_ut', $
        'calc_abund_c_ii_rl_ut', 'calc_abund_c_iii_rl_ut', $
        'calc_abund_n_ii_rl_ut', 'calc_abund_n_iii_rl_ut', $
        'calc_abund_o_ii_rl_ut', 'calc_abund_ne_ii_rl_ut'], $
        filename='test-results.log'

mgunit, ['redlaw_ut', 'redlaw_gal_ut', 'redlaw_gal2_ut', 'redlaw_ccm_ut', $
        'redlaw_jbk_ut', 'redlaw_fm_ut', 'redlaw_smc_ut', 'redlaw_lmc_ut', $
        'deredden_relflux_ut', 'deredden_flux_ut', $
        'calc_abundance_ut', 'calc_emissivity_ut', $
        'calc_temperature_ut', 'calc_density_ut', $
        'calc_abund_he_i_rl_ut', 'calc_abund_he_ii_rl_ut', $
        'calc_abund_c_ii_rl_ut', 'calc_abund_c_iii_rl_ut', $
        'calc_abund_n_ii_rl_ut', 'calc_abund_n_iii_rl_ut', $
        'calc_abund_o_ii_rl_ut', 'calc_abund_ne_ii_rl_ut'], $
        filename='test-results.html', /html

; --- End $MAIN$ program. ---------------
exit
