function read_porter
;+
; NAME:
;     read_porter
; PURPOSE:
;     read in tables of helium emissivities from Porter et al.
;     http://cdsads.u-strasbg.fr/abs/2012MNRAS.425L..28P
; EXPLANATION:
;
; CALLING SEQUENCE:
;     heidata=read_porter()
;
; RETURN:  recombination coefficients of He I
;
; REVISION HISTORY:
;     IDL code by A. Danehkar, 15/12/2013
;- 
  common share1, Atomic_Data_Path
  
  hei_ems=dblarr(21,14,44)
  temp=dblarr(46)
  
  nlines = 294 
  
  ion1='R_he_i'
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  for i=0, nlines-1 do begin 
    readf,lun1, temp
    tpos=nint((temp[0]/1000)-5)
    npos=nint(temp[1]-1)
    hei_ems[tpos,npos,*]=temp[2:45]
  endfor
  free_lun, lun1 
  
  return, hei_ems
end
