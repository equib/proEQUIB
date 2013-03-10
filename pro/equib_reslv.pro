pro equib_reslv, A, B, N, NR
;+
; NAME:
;     equib_reslv
; PURPOSE:
;     Resolve A with B
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_reslv, A, B, N, NR
;
; INPUTS:
;     A -     A parameter
;     B -     B parameter
;     N -     N parameter
;     NR -    NR parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;-  
  ;N= long(0)
  ;NR= long(0)
  NM1= long(0)
  I= long(0)
  J= long(0)
  K= long(0)
  L= long(0)
  IP1= long(0)        
  ;A=dblarr(NR+1,NR+1)
  ;B=dblarr(NR+1)  
  IF(N eq 1) then begin
    B[N]=B[N]/A[N,N]
    return
  endif
  NM1=N-1
  for I=1,NM1 do begin
    IP1=I+1
    for J=IP1,N do begin
      B[J]=B[J]-B[I]*A[J,I]/A[I,I]
    endfor
  endfor
  B[N]=B[N]/A[N,N]
  for I=1,NM1 do begin
    K=N-I
    L=K+1
    for J=L,N do begin
      B[K]=B[K]-B[J]*A[K,J]
    endfor
    B[K]=B[K]/A[K,K]
  endfor
end
