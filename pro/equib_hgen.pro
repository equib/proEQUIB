pro equib_hgen, XX, GH, Y, NPT, IOPT, NDIM, NDIMT3, HMH
;+
; NAME:
;     equib_hgen
; PURPOSE:
;     Cubic spline interpolation
;     The equation for the second derivatives at internal points
;     is of the form G*YPP=B, where G has been evaluated and LU
;     decomposed.
;     this routine writes B=HMH*Y and then solves YPP=G**(-1)*HMH*Y,
;     =HMH*Y.
;     Three options are provided for boundary conditions-
;     IOPT = 0  YPP=0 at end points
;     IOPT = 1  YP=0  at end points
;     IOPT = 2  YP at end points from lagarnge interpolant of a set of
;     internal points.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     equib_hgen, XX, GH, Y, NPT, IOPT, NDIM, NDIMT3, HMH
;
; INPUTS:
;     XX -     XX parameter
;     GH -     GH parameter
;     Y -      Y parameter
;     NPT -    NPT parameter
;     IOPT -   IOPT parameter
;     NDIM -   NDIM parameter
;     NDIMT3 - NDIMT3 parameter
;     HMH -    HMH parameter
; REVISION HISTORY:
;     Converted from FORTRAN EQUIB to IDL, 15/09/2013
;- 
  ;NPT= long(0)
  ;IOPT= long(0)
  ;NDIM= long(0)
  ;NDIMT3= long(0)
  NDIM3= long(0)
  NIP= long(0)
  I= long(0)
  J= long(0)
  K= long(0)
  NPM= long(0)
  INDX= long(0)
  ;XX=dblarr(NDIM+1)
  ;GH=dblarr(NDIMT3+1)
  ;Y=dblarr(NDIM+1)
  ;HMH=dblarr(NDIM+1,NDIM+1)
  XY=dblarr(5+1)     
  D=dblarr(5+1)     
  C=dblarr(2+1,5+1)     
  A0=double(0)
  AN1=double(0)
  H1=double(0)
  H2=double(0)
  ; Case of derivative boundary condition, with
  if (IOPT eq 2) then begin
    ; derivatives from NIP-point Lagrange at
    NDIM3=5
    ; internal points
    NIP=3
    for j=1, 2 do begin
      for i=1, NIP do begin 
        K=(NPT-NIP)*(J-1)
        XY[I]=XX[K+I]  
      endfor
      K=1+(NPT-1)*(J-1)
      equib_DERIV, XY, D, XX[K], NIP, NDIM3
      for i=1, NIP do begin 
        C[J,I]=D[I]  
      endfor
    endfor
  endif
  ; Set up matrix equation G*YPP=HMH*Y
  A0=XX[2]-XX[1]
  AN1=XX[NPT]-XX[NPT-1]
  NPM=NPT-2
  for I=1,NPM do begin
    H1=6./(XX[I+1]-XX[I])
    H2=6./(XX[I+2]-XX[I+1])
    for J=1,NPT do begin
      HMH[I,J]=0.
      if(J eq I) then HMH[I,J]=H1
      if(J eq I+2) then HMH[I,J]=H2
      if(J eq I+1) then HMH[I,J]=-H1-H2
    endfor
  endfor
  ;Correct matrix for case of
  if(IOPT eq 1) or (IOPT eq 2) then begin
    ; derivative boundary conditions
    HMH[1,1]=HMH[1,1]+3/A0
    HMH[1,2]=HMH[1,2]-3/A0
    HMH[NPM,NPT-1]=HMH[NPM,NPT-1]-3/AN1
    HMH[NPM,NPT]=HMH[NPM,NPT]+3/AN1
  endif
  if(IOPT eq 2) then begin
    for J=1,NIP do begin
      HMH[1,J]=HMH[1,J]+3*C[1,J]
      K=NPT+J-NIP
      HMH[NPM,K]=HMH[NPM,K]-3*C[2,J]
    endfor
  endif
  ;for I=1,NPM do begin
  ;endfor
  ; Solve matrix equation with results in the form
  for I=1,NPT do begin
    ; YPP=HMH*Y. matrix g has been LU decomposed
    Y[1]=HMH[1,I]
    INDX=0
    for J=2,NPM do begin
       INDX=INDX+3
       Y[J]=HMH[J,I]-GH[INDX]*Y[J-1]
    endfor
    INDX=INDX+1
    Y[NPM]=Y[NPM]/GH[INDX]
    for J=2,NPM do begin
      K=NPM-J+1
      INDX=INDX-3
      Y[K]=(Y[K]-GH[INDX+1]*Y[K+1])/GH[INDX]
    endfor
    for J=1,NPM do begin
      HMH[J+1,I]=Y[J]
    endfor
    ;Insert values for second derivative at end
    HMH[1,I]=0.
    ; points: first and last rows of the matrix
    HMH[NPT,I]=0.
  endfor
  ; Case of derivative boundary conditions
  if(IOPT gt 0) then begin
    for J=1,NPT do begin
      HMH[1,J]=-0.5*HMH[2,J]
      HMH[NPT,J]=-0.5*HMH[NPT-1,J]
    endfor
    HMH[1,1]=HMH[1,1]-3/(A0*A0)
    HMH[1,2]=HMH[1,2]+3/(A0*A0)
    HMH[NPT,NPT-1]=HMH[NPT,NPT-1]+3/(AN1*AN1)
    HMH[NPT,NPT]=HMH[NPT,NPT]-3/(AN1*AN1)
  endif
  IF(IOPT eq 2) then begin
    for J=1,NIP do begin
      HMH[1,J]=HMH[1,J]-3*C[1,J]/A0
      K=NPT+J-NIP
      HMH[NPT,K]=HMH[NPT,K]+3*C[2,J]/AN1
    endfor
  endif
  ;for I=1,NPT do begin
  ;endfor
end
