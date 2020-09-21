; docformat = 'rst'

function redlaw_lmc, wavelength
;+
;    This function determines the reddening law function of the line at the given wavelength
;    for the Large Magellanic Cloud.
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
;     IDL> fl=redlaw_lmc(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.30871187
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
;     Based on Formulae by Howarth 1983, MNRAS, 203, 301
;     1983MNRAS.203..301H
;     
;     Originally from IRAF STSDAS SYNPHOT ebmvlfunc.x, redlaw.x
;     
;     18/10/1994, R. A. Shaw, Initial IRAF implementation.
;     
;     14/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;     
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

  ; Tabulated inverse wavelengths in microns:
  xtab=[ 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82 ]

  ; Tabulated extinction function, A(lambda)/E(B-V), from Savage & Mathis:
  extab=[ 0.00,  0.16,  0.38,  0.87,  1.50,  2.32,  3.10 ]
 
  temp=  size(wavelength,/DIMENSIONS)
  if temp[0] eq 0 then begin
    npts=1
    extl=double(0.0)
  endif else begin
    npts = temp[0]
    extl = dblarr(npts)
  endelse
  for pix = 0, npts-1 do begin
    if (wavelength[pix] lt 1000.0) then print, "redlaw_lmc: Invalid wavelength"
    
    ; Convert input wavelength to inverse microns
    x = 10000.D+0 / wavelength[pix]

    ; Infrared - optical
    if ( x le 1.82) then begin
        ; linear interpolation of Savage & Mathis 1979
        ;  val = lin_interp(extab, xtab,  x)
        ;  extl[pix] = val ;+ 3.1
        ; Infrared - extend optical results linearly to 0 at 1/lam = 0
        val = ((1.86 - 0.48 * x) * x - 0.1) * x
    endif else begin 
      ; The following polynomial evaluations assume R = 3.1
      ; Renormalize extinction function to A(lambda)/A(V)
      if ( x le 2.75 ) then begin
        ;  Violet
        val = 3.1 + (2.04 + 0.094 * (x - 1.82)) * (x - 1.82)
      endif else begin
        ;  Ultraviolet out to lambda = 1000 A
        x    = min([x, 10.0])
        val  = 3.1 - 0.236 + 0.462 * x + 0.105 * x^2 + 0.454 / ((x - 4.557)^2 + 0.293)
      endelse
    endelse
    extl[pix] = val
  endfor
  return,  (extl / 3.57) - 1.0 
end
