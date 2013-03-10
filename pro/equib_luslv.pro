pro equib_luslv, A, B, N, M
;+
; NAME:
;     equib_luslv
; PURPOSE:
;     Solving linear equations
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_luslv, A, B, N, M
;
; INPUTS:
;     A -     A parameter
;     B -     B parameter
;     N -     N parameter
;     M -     M parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;-  
  ;M= long(0)
  ;N= long(0)
  ;A=dblarr(M+1,M+1)
  ;B=dblarr(M+1)  
  equib_LURED, A, N, M
  equib_RESLV, A, B, N, M
end
