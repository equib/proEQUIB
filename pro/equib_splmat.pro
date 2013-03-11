pro equib_splmat, XX, NPT, IOPT, NDIM, NDIMT3, HMH
;+
; NAME:
;     equib_splmat
; PURPOSE:
;     
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_splmat, XX, NPT, IOPT, NDIM, NDIMT3, HMH
;
; INPUTS:
;     XX -     XX parameter
;     NPT -    NPT parameter
;     IOPT -   IOPT parameter
;     NDIM -   NDIM parameter
;     NDIMT3 - NDIMT3 parameter
;     HMH -    HMH parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;-  
  ;NDIM= long(0)
  ;NDIMT3= long(0)
  ;NPT= long(0)
  ;IOPT= long(0)
  ;NPM= long(0)
  NELEM= long(0)
  ;XX=dblarr(NDIM)      
  GH=dblarr(NDIMT3+1)  
  Y=dblarr(NDIM+1)  
  ; HMH=dblarr(NDIM+1,NDIM+1)      
  NPM=NPT-2
  equib_GHGEN, GH, XX, NPT, IOPT, NDIM, NDIMT3
  NELEM=3*NPM-2
  equib_ELU, GH,NPM,NDIMT3
  equib_HGEN, XX, GH, Y, NPT, IOPT, NDIM, NDIMT3, HMH
end
