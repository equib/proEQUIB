function redlaw_ccm, wave, rv=rv
;+
; NAME:
;     redlaw_ccm
; PURPOSE:
;    reddening law function of Cardelli, Clayton & Mathis
; 
; EXPLANATION:
;
; CALLING SEQUENCE:
;     fl = redlaw_ccm(wave, rv)
;
; INPUTS:
;     wave[] -  wavelength of emission line, Angstroms
;     rv - The ratio of extinction to reddening defined as R_V = A_V/E(B-V)
; RETURN: extl[] -  extinction evaluation array
;
;     Converted to IDL code by A. Danehkar, 31/08/2012
;-

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
    if (wave[pix] lt 1000.0) then print, "redlaw_ccm: Invalid wavelength"
    
    ; Convert input wavelength to inverse microns
    x = 10000.D+0 / wave[pix]
    
    ; For wavelengths longward of the L passband, linearly interpolate
    ; to 0. at 1/lambda = 0.  (a=0.08, b=-0.0734 @ x=0.29)
    if (x lt 0.29D+0) then begin
      a =  0.2759 * x
      b = -0.2531 * x
    endif else begin 
      if (x lt 1.1D+0) then begin
        y = x ^ 1.61
        a =  0.574 * y
        b = -0.527 * y
      endif else begin 
        if (x lt 3.3D+0) then begin
          y = x - 1.82
          a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 + $
              y * (0.72085 + y * (0.01979 + y * (-0.77530 + y * 0.32999))))))
          b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * $
              (-5.38434 + y * (-0.62251 + y * (5.30260 - y * 2.09002))))))
        endif else begin 
          if (x lt 5.9D+0) then begin
            y = (x - 4.67) ^ 2
            a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
            b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)
          endif else begin 
            if (x lt 8.0D+0) then begin 
              y = (x - 4.67) ^ 2
              a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
              b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)
              
              y = x - 5.9
              a = a - 0.04473 * y^2 - 0.009779 * y^3
              b = b + 0.2130 * y^2 + 0.1207 * y^3
              ; Truncate range of ISEF to that for 1000 Ang.
            endif else begin 
              if (x le 10.0D+0) then begin
                x = min([x, 10.0D+0])
                y = x - 8.D+0
                a = -1.072 - 0.628 * y + 0.137 * y^2 - 0.070 * y^3
                b = 13.670 + 4.257 * y - 0.420 * y^2 + 0.374 * y^3
              endif
            endelse
          endelse
        endelse
      endelse
    endelse  
    ; Compute A(lambda)/A(V)
    y = a*R_V + b 
    extl[pix] = y
  endfor
  return, (extl / ((1.015452*R_V) + 0.461000)) - 1.0
end
