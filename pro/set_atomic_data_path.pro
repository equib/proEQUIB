pro set_atomic_data_path, path
;+
; NAME:
;     set_atomic_data_path
; PURPOSE:
;     Set Atomic Data Path
; EXPLANATION:
;
; CALLING SEQUENCE:
;     path='proEQUIB/atomic-data/'
;     set_atomic_data_path, path
;
; INPUTS:
;     path -     Atomic Data Path
; REVISION HISTORY:
;     A. Danehkar, 5/09/2013
;-  
  common share1, Atomic_Data_Path
  Atomic_Data_Path = path
end
