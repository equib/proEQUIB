function calc_abund_n_ii_rl_ut::test_basic
  compile_opt strictarr
  
  base_dir = file_dirname(file_dirname((routine_info('calc_abund_n_ii_rl_ut__define', /source)).path))
  data_rc_dir = ['externals', 'atomneb', 'atomic-data-rc']
  Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
  Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
  Atom_RC_PPB91_file= filepath('rc_PPB91.fits', root_dir=base_dir, subdir=data_rc_dir )
  Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )

  atom='h'
  ion='ii' ; H I
  h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
  
  atom='n'
  ion='iii' ; N II
  n_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
  n_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)

  h_i_aeff_data=h_i_rc_data[0].Aeff
  
  temperature=double(10000.0)
  density=double(5000.0)
  
  n_ii_4442_flux = 0.017
  wavelength=4442.02
  Abund_n_ii=calc_abund_n_ii_rl(temperature=temperature, density=density, $
                                wavelength=wavelength, line_flux=n_ii_4442_flux, $
                                n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data, $
                                h_i_aeff_data=h_i_aeff_data)
  result= long(Abund_n_ii*1e6)
  assert, result eq 692, 'incorrect result: %d', result

  return, 1
end

pro calc_abund_n_ii_rl_ut__define
  compile_opt strictarr

  define = { calc_abund_n_ii_rl_ut, inherits proEquibUTTestCase}
end

