pro equib_ghgen, GH, XX, NPT, IOPT, NDIM, NDIMT3
;+
; NAME:
;     equib_ghgen
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_ghgen, GH, XX, NPT, IOPT, NDIM, NDIMT3
;
; INPUTS:
;     GH -     GH parameter
;     XX -     XX parameter
;     NPT -    NPT parameter
;     IOPT -   IOPT parameter
;     NDIM -   NDIM parameter
;     NDIMT3 - NDIMT3 parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;- 
  ;NPT= long(0)
  ;IOPT= long(0)
  ;NDIM= long(0)
  ;NDIMT3= long(0)
  INDX= long(0)
  NPTM= long(0)
  I= long(0)
  J= long(0)
  IP= long(0)
  JP= long(0)
  IK= long(0)
  ;XX=dblarr(NDIM+1)
  ;GH=dblarr(NDIMT3+1)
  INDX=0
  NPTM=NPT-1
  for i=2, NPTM do begin
    IP=I-1
    for j=1, 3 do begin 
      JP=IP+J-2
      if(JP ge 1) and (JP le NPTM-1) then begin
        INDX=INDX+1
        if(J eq 2) then begin
          GH[INDX]=2*(XX[I+1]-XX[I-1])
        endif else begin
          IK=I+(J-1)/2
          GH[INDX]=XX[IK]-XX[IK-1]
        endelse
      endif   
    endfor
  endfor
  if(IOPT ge 1) then begin
    GH[1]=GH[1]-(XX[2]-XX[1])/2.
    GH[INDX]=GH[INDX]-(XX[NPT]-XX[NPT-1])/2.
  endif
end
