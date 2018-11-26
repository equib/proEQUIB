; docformat = 'rst'

function redlaw_smc, wavelength
;+
;    This function determines the reddening law function of the line at the given wavelength 
;    for Small Magellanic Cloud.
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
;     IDL> fl=redlaw_smc(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.22659261
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
;   0.0.1
;
; :History:
;     Based on Prevot et al. (1984), A&A, 132, 389-392
;     1984A%26A...132..389P
;     
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;     
;     20/09/1994, R. A. Shaw, Initial IRAF implementation.
;     
;     04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;     
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

;+
; NAME:
;     redlaw_smc
; PURPOSE:
;    This function determines the reddening law function of the line at the given wavelength 
;    for Small Magellanic Cloud.
;
; CALLING SEQUENCE:
;     Result = redlaw_smc(Wavelength)
;
; INPUTS:
;     Wavelength[] -  in, required, type=float/array, 
;               wavelength in Angstroms
; 
; OUTPUTS: This function returns a double/array  as the reddening law function 
;                   value(s) f(lambda) for the given wavelength(s) lambda.
;
; PROCEDURE: This function is callsed by redlaw.
;
; EXAMPLE:
;     wavelength=6563.0
;     fl=redlaw_smc(wavelength)
;     print, 'fl(6563)', fl
;     > fl(6563)     -0.22659261
; 
; MODIFICATION HISTORY:
;     Based on Prevot et al. (1984), A&A, 132, 389-392
;     1984A%26A...132..389P
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;     20/09/1994, R. A. Shaw, Initial IRAF implementation.
;     04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

  ; Tabulated inverse wavelengths in microns:
  xtab=[ 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82, $
         2.35,  2.70,  3.22,  3.34,  3.46,  3.60,  3.75,  3.92,  4.09,  4.28, $
         4.50,  4.73,  5.00,  5.24,  5.38,  5.52,  5.70,  5.88,  6.07,  6.27, $
         6.48,  6.72,  6.98,  7.23,  7.52,  7.84]
  
  ; Tabulated extinction function, E(lambda-V)/E(B-V):
  extab=[-3.10, -2.94, -2.72, -2.23, -1.60, -0.78,  0.00, $
          1.00,  1.67,  2.29,  2.65,  3.00,  3.15,  3.49,  3.91,  4.24,  4.53, $
          5.30,  5.85,  6.38,  6.76,  6.90,  7.17,  7.71,  8.01,  8.49,  9.06,  $
          9.28,  9.84, 10.80, 11.51, 12.52, 13.54]

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
    val = lin_interp(extab, xtab,  x)
    ;deriv1 = spl_init(xtab, extab)
    ;val=spl_interp(xtab, extab, deriv1, x)
          
    extl[pix] = val + 3.1
  endfor
  return, (extl / 3.242) - 1.0 
end
