function redlaw_gal_ut::test_basic
  compile_opt strictarr
  
  wavelength=6563.0
  m_ext=1.0
  flux=1.0
  R_V=3.1

  fl=redlaw_gal(wavelength, rv=R_V)
  
  result= long(fl*1e4)
  assert, result eq -3201, 'incorrect result: %d', result

  return, 1
end

pro redlaw_gal_ut__define
  compile_opt strictarr

  define = { redlaw_gal_ut, inherits proEquibUTTestCase}
end

