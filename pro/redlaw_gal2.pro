; docformat = 'rst'

function redlaw_gal2, wavelength
;+
;    This function determines the reddening law function of the line at the given wavelength
;    for Galactic Savage & Mathis 1979.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).
;
; :Params:
;     wavelength :  in, required, type=float
;                   Wavelength in Angstrom
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> fl=redlaw_gal2(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.30925984
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;     Subroutines
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
;     Based on Savage & Mathis 1979, ARA&A, vol. 17, 73-111
;     
;     Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x
;     
;     20/09/1994, R. A. Shaw, Initial IRAF implementation.
;     
;     04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;     
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

  ; Tabulated inverse wavelengths in microns:
  xtable=[ 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82, $
	   2.27,  2.50,  2.91,  3.65,  4.00,  4.17,  4.35,  4.57,  4.76,  5.00, $
	   5.26,  5.56,  5.88,  6.25,  6.71,  7.18,  8.00,  8.50,  9.00,  9.50, $
	  10.00 ]
  
  ; Tabulated extinction function, A(lambda)/E(B-V):
  etable=[0.00,  0.16,  0.38,  0.87,  1.50,  2.32,  3.10, $
	  4.10,  4.40,  4.90,  6.20,  7.29,  8.00,  8.87,  9.67,  9.33,  8.62, $
	  8.00,  7.75,  7.87,  8.12,  8.15,  8.49,  9.65, 10.55, 11.55, 12.90, $
	  14.40 ]


  temp=  size(wavelength,/DIMENSIONS)
  if temp[0] eq 0 then begin
    npts=1
    extl=double(0.0)
  endif else begin
    npts = temp[0]
    extl = dblarr(npts)
  endelse
  for pix = 0, npts-1 do begin
    if (wavelength[pix] lt 1000.0) then print, "redlaw_smc: Invalid wavelength"	
    ; Convert wavelength in angstroms to 1/microns
    x = 10000.D+0 / wavelength[pix]
    x = min([x, 10.0])

    ; Linearly interpolate extinction law in 1/lam
    val = lin_interp(etable, xtable,  x)
    ;deriv1 = spl_init(xtab, extab)
    ;val=spl_interp(xtab, extab, deriv1, x)

    extl[pix] = val
  endfor
  return, (extl / 3.63) - 1.0 
end
