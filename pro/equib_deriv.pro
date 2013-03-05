pro equib_deriv, XY, D, X, N, NDIM
;+
; NAME:
;     equib_deriv
; PURPOSE:
;     Calculate the first derivative of the lagrangian interpolator
;     of a function F, tabulated at the N points XY(I), I=1 to N.
;     The derivative is given as the coefficients of F(I), I=1 to N,
;     in the array D(I), I=1 to N.   
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_deriv, XY, D, X, N, NDIM
;
; INPUTS:
;     XY -     XX parameter
;     D -      D parameter
;     X -      X parameter
;     N -      N parameter
;     NDIM -   NDIM parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;- 
  ;N= long(0)
  ;NDIM= long(0)
  I= long(0)
  J= long(0)
  K= long(0)  
  ;XY=dblarr(NDIM+1)
  ;D=dblarr(NDIM+1)
  ;X=double(0)
  P1=double(0)
  P2=double(0)
  S=double(0)  

  for I=1,N do begin
    P1=1.
    S=0.
    for J=1,N do begin
      if(J ne I) then begin
        P1=P1*(XY[I]-XY[J])
        P2=1.
        for K=1,N do begin
          if(K ne I) and (K ne J) then P2=P2*(X-XY[K])
        endfor
        S=S+P2
      endif
    endfor
    D[I]=S/P1
  endfor
end
