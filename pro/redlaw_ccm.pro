; docformat = 'rst'

function redlaw_ccm, wavelength, rv=rv
;+
;    This function determines the reddening law function of Cardelli, Clayton & Mathis.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value for the given wavelength.
;
; :Params:
;     wavelength :  in, required, type=float/array
;                   Wavelength in Angstrom
;
; :Keywords:
;    RV       :  in, optional, type=float, default=3.1
;                the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> R_V=3.1
;     IDL> fl=redlaw_ccm(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.29756615
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.3.0
;
; :History:
;     Based on Formulae by Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256.
;     1989ApJ...345..245C
;     
;     Originally from IRAF STSDAS SYNPHOT redlaw.x
;     
;     18/05/1993, R. A. Shaw, Initial IRAF implementation, based upon CCM module
;         in onedspec.deredden.
;     
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

;+
; NAME:
;     redlaw_ccm
; 
; PURPOSE:
;    This function determines the reddening law function of Cardelli, Clayton & Mathis.
;
; CALLING SEQUENCE:
;     fl = redlaw_ccm(Wavelength, RV=rv)
;
; INPUTS:
;     Wavelength[] -  in, required, type=float/array, 
;               wavelength in Angstroms
; 
; KEYWORD PARAMETERS:
;    RV       :  in, optional, type=float, default=3.1, 
;                the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).
; 
; OUTPUTS: This function returns a double/array  as the reddening law function 
;                   value(s) f(lambda) for the given wavelength(s) lambda.
;
; PROCEDURE: This function is callsed by redlaw.
;
; EXAMPLE:
;     wavelength=6563.0
;     R_V=3.1
;     fl=redlaw_ccm(wavelength, rv=R_V)
;     print, 'fl(6563)', fl
;     > fl(6563)     -0.29756615
;
; MODIFICATION HISTORY:
;     Based on Formulae by Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256.
;     1989ApJ...345..245C
;     Originally from IRAF STSDAS SYNPHOT redlaw.x
;     18/05/1993, R. A. Shaw, Initial IRAF implementation, based upon CCM module
;         in onedspec.deredden.
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

  temp=  size(wavelength,/DIMENSIONS)
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
    if (wavelength[pix] lt 1000.0) then print, "redlaw_ccm: Invalid wavelength"
    
    ; Convert input wavelength to inverse microns
    x = 10000.D+0 / wavelength[pix]
    
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
