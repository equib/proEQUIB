function redlaw_fm, wave, fmlaw=fmlaw, rv=rv
;+
; NAME:
;     redlaw_fm
; PURPOSE:
;    reddening law function of Fitzpatrick & Massa
; 
; EXPLANATION:
;
; CALLING SEQUENCE:
;     fl = redlaw_fm(wave)
;
; INPUTS:
;     wave[] -  wavelength of emission line, Angstroms
; RETURN: extl[] -  extinction evaluation array
;
; REVISION HISTORY:
;     Based on Formulae by Fitzpatrick 1999, PASP, 11, 63
;     1999PASP..111...63F, Fitzpatrick & Massa 1990, 
;     ApJS, 72, 163, 1990ApJS...72..163F
;     Adopted from NASA IDL Library & PyAstronomy
;     Revised in IDL code by A. Danehkar, 30/12/2016
;-
; Tabulated inverse wavelengths in microns:

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
  x0    =  4.596  
  gamma1 =  0.99 
  c3    =  3.23  
  c4   =  0.41    
  c2    = -0.824 + 4.717/R_V
  c1    =  2.030 - 3.007*c2
  if keyword_set(fmlaw) then begin
    if fmlaw eq 'LMC2' then begin
       x0    =  4.626
       gamma1 =  1.05 
       c4   =  0.42   
       c3    =  1.92  
       c2    = 1.31
       c1    =  -2.16
    endif else begin
      if fmlaw eq 'AVGLMC' then begin
         x0 = 4.596  
         gamma1 = 0.91
         c4   =  0.64  
         c3    =  2.73  
         c2    = 1.11
         c1    =  -1.28
      endif else begin
         x0    =  4.596  
         gamma1 =  0.99 
         c3    =  3.23  
         c4   =  0.41    
         c2    = -0.824 + 4.717/R_V
         c1    =  2.030 - 3.007*c2
      endelse
    endelse
  endif
  for pix = 0, npts-1 do begin
    ; Convert input wavelength to inverse microns
    x = 10000.D+0 / wave[pix]
    curve = x*0.
    
    ; Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
    ; R-dependent coefficients
    xcutuv = 10000.0/2700.0
    xspluv = 10000.0/[2700.0,2600.0]
     
    iuv = where(x ge xcutuv)
    temp = size(iuv,/DIMENSIONS)
    if temp[0] eq 0 then N_UV=0 else N_UV = temp[0]
    iopir = where(x le xcutuv)
    temp = size(iopir,/DIMENSIONS)
    if temp[0] eq 0 then Nopir=0 else Nopir = temp[0]
    
    if (N_UV gt 0) then xuv = [xspluv,x[iuv]] else  xuv = xspluv
  
    yuv = c1  + c2*xuv
    yuv = yuv + c3*xuv^2/((xuv^2-x0^2)^2 +(xuv*gamma1)^2)
    yuv = yuv + c4*(0.5392*((xuv>5.9)-5.9)^2+0.05644*((xuv>5.9)-5.9)^3)
    yuv = yuv + R_V
    yspluv  = yuv[0:1]                  ; save spline points
  
    if (N_UV gt 0) then curve[iuv] = yuv[2:*]      ; remove spline points
   
    ; Compute optical portion of A(lambda)/E(B-V) curve
    ; using cubic spline anchored in UV, optical, and IR
    xsplopir = [0,10000.0/[26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]]
    ysplir   = [0.0,0.26469,0.82925]*R_V/3.1 
    ysplop   = [poly(R_V, [-4.22809e-01, 1.00270, 2.13572e-04] ), $
               poly(R_V, [-5.13540e-02, 1.00216, -7.35778e-05] ), $
               poly(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05] ), $
               poly(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, $ 
                       -4.45636e-05] ) ]
    
    ysplopir = [ysplir,ysplop]
  
    if (Nopir GT 0) then begin
      curve[iopir] = cspline([xsplopir,xspluv],[ysplopir,yspluv],x[iopir])
    endif
    extl[pix] = curve 
  endfor
  return, (extl / 3.63) - 1.0 
end
