function redlaw_lmc_ut::test_basic
  compile_opt strictarr
  
  wavelength=6563.0
  m_ext=1.0
  flux=1.0
  R_V=3.1

  fl=redlaw_lmc(wavelength)
  
  result= long(fl*1e4)
  assert, result eq -3087, 'incorrect result: %d', result

  return, 1
end

pro redlaw_lmc_ut__define
  compile_opt strictarr

  define = { redlaw_lmc_ut, inherits MGutLibTestCase }
end

