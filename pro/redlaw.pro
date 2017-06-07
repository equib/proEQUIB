function redlaw, wavelength, ext_law=ext_law, rv=rv, fmlaw=fmlaw
;+
; NAME:
;     redlaw
; PURPOSE:
;     determine reddening law function for the line at the wavelength of lambda, 
;     depending on the law used
; EXPLANATION:
;
; CALLING SEQUENCE:
;     fl = redlaw(wavelength, ext_law)
;
; INPUTS:
;     wavelength -     Wavelength in Angstrom
;     ext_law -     extinction law
;        ext_law='GAL' ; Howarth Galactic
;        ext_law='GAL2' ; Savage and Mathis
;        ext_law='CCM' ; CCM galactic
;        ext_law='JBK' ; Whitford, Seaton, Kaler
;        ext_law='FM' ; Fitxpatrick
;        ext_law='SMC' ; Prevot SMC
;        ext_law='LMC' ; Howarth LMC
; RETURN:  f(lambda)
;
; REVISION HISTORY:
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;     Converted to IDL code by A. Danehkar, 31/08/2012
;-  
  if keyword_set(ext_law) then begin
    extlaw=ext_law
  endif else begin
    extlaw='GAL'
  endelse
  case extlaw of
     'GAL': fl = redlaw_gal(wavelength, rv=rv)
     'GAL2': fl = redlaw_gal2(wavelength)
     'CCM': fl = redlaw_ccm(wavelength,rv=rv)
     'JBK': fl = redlaw_jbk(wavelegth)
     'FM': fl = redlaw_fm(wavelength,fmlaw=fmlaw,rv=rv)
     'SMC': fl = redlaw_smc(wavelength)
     'LMC': fl = redlaw_lmc(wavelength)
     else: print, 'ext_law cannnot find'
  endcase
  return, fl
end
