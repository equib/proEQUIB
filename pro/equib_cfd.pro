pro equib_cfd, X, XX, NPT, NDIM, HMH, D
;+
; NAME:
;     equib_cfd
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_cfd, X, XX, NPT, NDIM, HMH, D
;
; INPUTS:
;     X -     X parameter
;     XX -    XX parameter
;     NPT -   NPT parameter
;     NDIM -  NDIM parameter
;     HMH -   HMH parameter
;     D -     D parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;-
  ;NPT= long(0)
  ;NDIM= long(0)
  NPTM= long(0)
  I= long(0)
  J= long(0)
  ;X= double(0) 
  ;XX=dblarr(NDIM+1)
  ;HMH=dblarr(NDIM+1,NDIM+1)
  ;D=dblarr(NDIM+1)
  X1= double(0)
  X2= double(0)
  A1= double(0) 
  A2= double(0)
  HI= double(0)
  if (X lt XX[1]) then begin
    ;print, XX[1]
    return
  endif
  if (X gt XX[NPT]) then begin
    ;print, XX[NPT]
    return
  endif
  NPTM=NPT-1
  for i=1, NPTM do begin 
    if(X lt XX[i+1]) then begin
      X1=XX[i+1]-X
      X2=X-XX[i]
      HI=XX[i+1]-XX[i]
      A1=X1*(X1*X1/(6*HI)-HI/6)
      A2=X2*(X2*X2/(6*HI)-HI/6)
      for j=1, NPT do begin 
        D[j]=A1*HMH[i,j]+A2*HMH[i+1,j]
      endfor
      D[i]=D[i]+X1/HI
      D[i+1]=D[i+1]+X2/HI
      return
    endif 
  endfor
end
