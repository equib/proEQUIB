function calc_populations_ut::test_basic
  compile_opt strictarr
  
  base_dir = file_dirname(file_dirname((routine_info('calc_temperature_ut__define', /source)).path))
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
  Nlj=calc_populations(temperature=temperature, density=density, $
                       elj_data=s_ii_elj, omij_data=s_ii_omij, $
                       aij_data=s_ii_aij)
  result= long(Nlj[0]*1e3)
  assert, result eq 969, 'incorrect result: %d', result
  
  return, 1
end

pro calc_populations_ut__define
  compile_opt strictarr

  define = { calc_populations_ut, inherits proEquibUTTestCase}
end

