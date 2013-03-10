pro equib_lured, A, N, NR
;+
; NAME:
;     equib_lured
; PURPOSE:
;     
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_lured, A, N, NR
;
; INPUTS:
;     A -     A parameter
;     N -     N parameter
;     NR -     NR parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;-
  ; N= long(0)
  ; NR= long(0)
  NM1= long(0)
  I= long(0)
  J= long(0)
  K= long(0)
  IP1= long(0)      
  ;A=dblarr(NR+1,NR+1)
  FACT=double(0)
  if (N eq 1) then return
  NM1=N-1
  for I=1,NM1 do begin
    IP1=I+1
    for K=IP1,N do begin
      FACT=A[K,I]/A[I,I]
      for J=IP1,N do begin
        A[K,J]=A[K,J]-A[I,J]*FACT
      endfor
    endfor
  endfor
end
