function deredden_relflux, wavelength, relflux, m_ext, ext_law=ext_law, rv=rv, fmlaw=fmlaw
;+
; NAME:
;     deredden_relflux
; PURPOSE:
;     determine deredden flux intensity relative to Hb=100, 
;     depending on the law used
; EXPLANATION:
;
; CALLING SEQUENCE:
;     relflux_deredden = deredden_relflux(wavelengths, relflux, m_ext)
;
; INPUTS:
;     wavelength -     Wavelength in Angstrom
;     relflux -      flux intensity relative to Hb=100
;     m_ext -      logarithmic extinction 
;     ext_law -     extinction law
;        ext_law='GAL' ; Howarth Galactic
;        ext_law='GAL2' ; Savage and Mathis
;        ext_law='CCM' ; CCM galactic
;        ext_law='JBK' ; Whitford, Seaton, Kaler
;        ext_law='FM' ; Fitxpatrick
;        ext_law='SMC' ; Prevot SMC
;        ext_law='LMC' ; Howarth LMC
; RETURN:  deredden relative intensity
;
; REVISION HISTORY:
;     IDL code by A. Danehkar, 31/08/2012
;-  
  if keyword_set(ext_law) then begin
    extlaw=ext_law
  endif else begin
    extlaw='GAL'
  endelse
  if keyword_set(rv) then begin
    R_V=rv
  endif else begin
    R_V=3.1
  endelse
  if keyword_set(fmlaw) then begin
    fm_law = fmlaw
  endif else begin
    fm_law = 'STANDARD'
  endelse
  fl = redlaw(wavelength,ext_law=extlaw,rv=R_V, fmlaw=fm_law)
  int_dered = relflux * 10.0^(m_ext*fl)
  return, int_dered
end
