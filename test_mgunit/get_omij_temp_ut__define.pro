function get_omij_temp_ut::test_basic
  compile_opt strictarr
  
  base_dir = file_dirname(file_dirname((routine_info('get_omij_temp_ut__define', /source)).path))
  data_dir = ['externals', 'atomneb', 'atomic-data', 'chianti70']
  Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
  Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
  Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )

  atom='s'
  ion='ii'
  s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej)
  s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
  s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)
  density = double(1000)
  temperature=double(10000.0);
  Omij_T=get_omij_temp(temperature=temperature, omij_data=s_ii_omij)
  result= long(Omij_T[1,2]*1e2)
  assert, result eq 746, 'incorrect result: %d', result
  
  return, 1
end

pro get_omij_temp_ut__define
  compile_opt strictarr

  define = { get_omij_temp_ut, inherits proEquibUTTestCase}
end

