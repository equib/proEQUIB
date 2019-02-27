function deredden_relflux_ut::test_basic
  compile_opt strictarr
  
  wavelength=6563.0
  m_ext=1.0
  flux=1.0
  R_V=3.1
  ext_law='GAL'

  flux_deredden=deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
  
  result= long(flux_deredden*1e5)
  assert, result eq 47847, 'incorrect result: %d', result
  
  return, 1
end

pro deredden_relflux_ut__define
  compile_opt strictarr

  define = { deredden_relflux_ut, inherits proEquibUTTestCase}
end

