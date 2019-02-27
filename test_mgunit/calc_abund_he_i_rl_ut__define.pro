function calc_abund_he_i_rl_ut::test_basic
  compile_opt strictarr
  
  base_dir = file_dirname(file_dirname((routine_info('calc_abundance_ut__define', /source)).path))
  data_rc_dir = ['externals', 'atomneb', 'atomic-data-rc']
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

  h_i_aeff_data=h_i_rc_data[0].Aeff
  he_i_aeff_data=he_i_rc_data[0].Aeff
  
  temperature=double(10000.0)
  density=double(5000.0)
  
  he_i_4471_flux= 2.104
  linenum=10; 4471.50
  Abund_he_i=calc_abund_he_i_rl(temperature=temperature, density=density, $
                                linenum=linenum, line_flux=he_i_4471_flux, $
                                he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
  result= long(Abund_he_i*1e4)
  assert, result eq 408, 'incorrect result: %d', result

  return, 1
end

pro calc_abund_he_i_rl_ut__define
  compile_opt strictarr

  define = { calc_abund_he_i_rl_ut, inherits proEquibUTTestCase}
end

