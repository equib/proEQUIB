function _sign, A, B
;+
; NAME:
;     _sign
; PURPOSE:
;     
; EXPLANATION:
;
; CALLING SEQUENCE:
;     ret= _sign(A, B)
;
; INPUTS:
;     A -     A parameter
;     B -     B parameter
; RETURN:  value
;-  
  if B lt 0 then return, -abs(A) else return, abs(A)
end
