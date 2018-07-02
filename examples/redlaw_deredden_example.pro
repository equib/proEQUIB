; Example: redlaw() and deredden_flux ()
;     determine reddening law function for
;     given wavelength and deredden flux
;     intensity (relative to Hb=100) based on
;     the reddening law
;
; --- Begin $MAIN$ program. ---------------
;
;
wavelength=6563.0
m_ext=1.0
flux=1.0
R_V=3.1

fl=redlaw_gal(wavelength, rv=R_V)
print, 'fl(6563)', fl

fl=redlaw_gal2(wavelength)
print, 'fl(6563)', fl

fl=redlaw_ccm(wavelength, rv=R_V)
print, 'fl(6563)', fl

fl=redlaw_jbk(wavelength)
print, 'fl(6563)', fl

;fmlaw='LMC2'
fmlaw='AVGLMC'
fl=redlaw_fm(wavelength, fmlaw=fmlaw, rv=R_V)
print, 'fl(6563)', fl

fl=redlaw_smc(wavelength)
print, 'fl(6563)', fl

fl=redlaw_lmc(wavelength)
print, 'fl(6563)', fl

fl=redlaw(wavelength, rv=R_V)
print, 'fl(6563)', fl

ext_law='GAL'
R_V=3.1
; deredden flux intensity relative to Hb=100
flux_deredden=deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
print, 'dereddened flux(6563)', flux_deredden

; deredden absolute flux intensity
flux_deredden=deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
print, 'dereddened flux(6563)', flux_deredden

end 