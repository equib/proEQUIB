pro equib_elu, GH, N, NDIM
;+
; NAME:
;     equib_elu
; PURPOSE:
; 
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_elu, GH, N, NDIM
;
; INPUTS:
;     GH -     GH parameter
;     N -      N parameter
;     NDIM -   NDIM parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;- 
  ;N= long(0)
  ;NDIM= long(0)
  INDX= long(0)
  I= long(0)
  J= long(0)
  JP= long(0)   
  ;GH=dblarr(NDIM+1)
  INDX=0
  for i=1, N do begin
    for j=1, 3 do begin 
      JP=I+J-2
      if(JP ge 1) and (JP le N) then begin
        INDX=INDX+1
        if(I gt 1) then begin
          if(J eq 1) then begin
            GH[INDX]=GH[INDX]/GH[INDX-2]
          endif
          if(J eq 2) then begin
            GH[INDX]=GH[INDX]-GH[INDX-1]*GH[INDX-2]
          endif
        endif
      endif
     endfor
  endfor
end
