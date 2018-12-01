function redlaw_fm_ut::test_basic
  compile_opt strictarr
  
  wavelength=6563.0
  m_ext=1.0
  flux=1.0
  R_V=3.1
  fmlaw='AVGLMC'
  
  fl=redlaw_fm(wavelength, fmlaw=fmlaw, rv=R_V)
  
  result= long(fl*1e4)
  assert, result eq -3505, 'incorrect result: %d', result

  return, 1
end

pro redlaw_fm_ut__define
  compile_opt strictarr

  define = { redlaw_fm_ut, inherits MGutLibTestCase }
end

