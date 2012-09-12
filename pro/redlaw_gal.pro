function redlaw_gal, wave, rv=rv
;+
; NAME:
;     redlaw_smc
; PURPOSE:
;    reddening law function for Galactic Seaton1979+Howarth1983+CCM1983
; 
; EXPLANATION:
;
; CALLING SEQUENCE:
;     fl = redlaw_gal(wave)
;
; INPUTS:
;     wave[] -  wavelength of emission line, Angstroms
; RETURN: extl[] -  extinction evaluation array
;
; REVISION HISTORY:
;     Based on the UV Formulae from Seaton 1979, MNRAS, 187, 73
;     1979MNRAS.187P..73S, the opt/NIR from Howarth 1983, MNRAS, 203, 301
;     the FIR from Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
;     1989ApJ...345..245C
;     Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x, pyneb.extinction
;     Converted to IDL code by A. Danehkar, 31/08/2012
;-
; Tabulated inverse wavelengths in microns:

  xtable = [ 0.,  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, $
             1.6, 1.7, 1.8, 1.9, 2.0, 2.1, $
             2.2, 2.3, 2.4, 2.5, 2.6, 2.7 ]
  etable= [ 0.,   1.36, 1.64, 1.84, 2.04, 2.24, 2.44, $
           2.66, 2.88, 3.14, 3.36, 3.56, 3.77, $
           3.96, 4.15, 4.26, 4.40, 4.52, 4.64 ]
               
  temp=  size(wave,/DIMENSIONS)
  if temp[0] eq 0 then begin
    npts=1
    extl=double(0.0)
  endif else begin
    npts = temp[0]
    extl = dblarr(npts)
  endelse
  if keyword_set(rv) then begin
    R_V=rv
  endif else begin
    R_V=3.1
  endelse
  for pix = 0, npts-1 do begin
    ; Convert wavelength in angstroms to 1/microns
    x = 10000.D+0  / wave[pix]
    
    if (x le 1.1) then begin
      ; Infrared: extend optical results linearly to 0 at 1/lam = 0
      ; extl[pix] = etable[2] * x^2
      ; Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
      extl[pix] = R_V * (0.574 * x^1.61) - 0.527 * x^1.61
    endif else begin
      if (x lt 1.83) then begin
        ; Howarth 1983, Galactic
        extl[pix] = (R_V - 3.1)+(( (1.86*x^2.0) - (0.48*x^3.0) - (0.1*x)))   
      endif else begin
        if (x lt 2.75) then begin
          ; Optical region interpolates in Seaton's table 3
          ;extl[pix] = lin_interp(etable, xtable,  x)
          ; Howarth 1983, MNRAS, 203, 301
          extl[pix]  = R_V + 2.56*(x-1.83) - 0.993*(x-1.83)^2.0 ;
        endif else begin 
          if (x lt 3.65) then begin
            ; Ultraviolet uses analytic formulae from Seaton's table 2
            extl[pix] = (R_V - 3.2) + 1.56 + 1.048 * x + 1.01 /((x - 4.6)^2 + 0.280)
          endif else begin 
            if (x lt 7.14) then begin
              ; More ultraviolet
              extl[pix] = (R_V - 3.2) + 2.29 + 0.848 * x + 1.01 / ((x - 4.6)^2 + 0.280)
              ; And more ultraviolet still
            endif else begin
              x = min([x, 50.0])
              extl[pix] = (R_V - 3.2) + 16.17 -3.20 *x+ 0.2975 * x^2
            endelse
          endelse
        endelse
      endelse
    endelse
  endfor
  return, (extl / 3.63) - 1.0 
end
