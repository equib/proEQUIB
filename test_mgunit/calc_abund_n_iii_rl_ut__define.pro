function calc_abund_n_iii_rl_ut::test_basic
  compile_opt strictarr
  
  base_dir = file_dirname(file_dirname((routine_info('calc_abund_n_iii_rl_ut__define', /source)).path))
  data_rc_dir = ['externals', 'atomneb', 'atomic-data-rc']
  Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
  Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
  Atom_RC_PPB91_file= filepath('rc_PPB91.fits', root_dir=base_dir, subdir=data_rc_dir )
  Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )

  atom='h'
  ion='ii' ; H I
  h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
  
  atom='n'
  ion='iv' ; N III
  n_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)

  h_i_aeff_data=h_i_rc_data[0].Aeff
  
  temperature=double(10000.0)
  density=double(5000.0)
  
  n_iii_4641_flux = 0.245
  wavelength=4640.64
  Abund_n_iii=calc_abund_n_iii_rl(temperature=temperature, density=density, $
                                  wavelength=wavelength, line_flux=n_iii_4641_flux, $
                                  n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
  result= long(Abund_n_iii*1e7)
  assert, result eq 633, 'incorrect result: %d', result

  return, 1
end

pro calc_abund_n_iii_rl_ut__define
  compile_opt strictarr

  define = { calc_abund_n_iii_rl_ut, inherits proEquibUTTestCase}
end

