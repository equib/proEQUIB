; docformat = 'rst'

function calc_abundance_ut::test_basic
  compile_opt strictarr
  
  base_dir = file_dirname(file_dirname((routine_info('calc_abundance_ut__define', /source)).path))
  data_dir = ['externals', 'atomneb', 'atomic-data', 'chianti70']
  data_rc_dir = ['externals', 'atomneb', 'atomic-data-rc']
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
  
  Abb5007=calc_abundance(temperature=temperature, density=density, $
                         line_flux=iobs5007, atomic_levels=levels5007,$
                         elj_data=o_iii_elj, omij_data=o_iii_omij, $
                         aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data[0].Aeff)
  result= long(Abb5007*1e7)
  assert, result eq 4125, 'incorrect result: %d', result
  
  return, 1
end

pro calc_abundance_ut__define
  compile_opt strictarr
  
  define = { calc_abundance_ut, inherits MGutLibTestCase }
end
