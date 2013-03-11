function equib_sign, A, B
;+
; NAME:
;     equib_sign
; PURPOSE:
;     
; EXPLANATION:
;
; CALLING SEQUENCE:
;     ret= equib_sign(A, B)
;
; INPUTS:
;     A -     A parameter
;     B -     B parameter
; RETURN:  value
;-  
  if B lt 0 then return, -abs(A) else return, abs(A)
end
