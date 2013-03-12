function equib_str2int, str
;+
; NAME :
;       equib_str2int
; PURPOSE:
;       Convert string of number to integer.
; CALLING SEQUENCE:
;       xx = equib_str2int(str)
; INPUTS:
;       str   = string of number
; MODIFICATION HISTORY
;       Written April 10, 1992 by Bachtiar Anwar
;-
ref_num = indgen(10)

ref_str = string(format='(i0)',ref_num)

numb = lonarr(20)
st_len=strlen(str)
;
for i = 0, st_len - 1 do begin
  for j = 0, 9 do begin
  if (strmid(str,i,1) EQ ref_str(j)) then $
      numb(i) = ref_num(j)*10^(st_len - (i+1))
  endfor
endfor
n = total(numb)
num=intarr(2)
num(0) = n

return, num(0)
end

