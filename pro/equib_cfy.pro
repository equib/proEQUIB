pro equib_cfy, X, Y, XX, YY, NPT, NDIM, HMH, D
;+
; NAME:
;     equib_cfy
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_cfy, X, Y, XX, YY, NPT, NDIM, HMH, D
;
; INPUTS:
;     X -     XX parameter
;     Y -     GH parameter
;     XX -    Y parameter
;     YY -    NPT parameter
;     NPT -   IOPT parameter
;     NDIM -  NDIM parameter
;     HMH -   NDIMT3 parameter
;     D -     HMH parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;- 
  ;NPT= long(0)
  ;NDIM= long(0)
  J= long(0)
  ;;XX=dblarr(NDIM+1)
  ;YY=dblarr(NDIM+1)
  ;HMH=dblarr(NDIM+1,NDIM+1)
  ;D=dblarr(NDIM+1)
  ;X= double(0) 
  ;Y= double(0)       
  TT= double(0) 
  if (X lt XX[1]) then begin
    Y=YY[1] 
  endif
  if(X gt XX[NPT]) then begin
    Y=YY[NPT] 
  endif
  TT=0.0
  for j=1, NPT do begin 
    TT=TT+D[J]*YY[J]
  endfor
  Y=TT
end
