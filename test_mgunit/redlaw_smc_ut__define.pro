function redlaw_smc_ut::test_basic
  compile_opt strictarr
  
  wavelength=6563.0
  m_ext=1.0
  flux=1.0
  R_V=3.1

  fl=redlaw_smc(wavelength)
  
  result= long(fl*1e4)
  assert, result eq -2265, 'incorrect result: %d', result

  return, 1
end

pro redlaw_smc_ut__define
  compile_opt strictarr

  define = { redlaw_smc_ut, inherits proEquibUTTestCase}
end

