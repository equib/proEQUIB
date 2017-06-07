function find_aeff_sh95_column, lo_lev, hi_lev, lev_num
  ;lev_num=25
  count=2
  ;lo_lev=3
  ;hi_lev=2
  for k=lev_num, 1, -1 do begin 
    for l=1, k-1 do begin 
        count=count+1
        if (k eq lo_lev) and (l eq hi_lev) then begin
          return, count
        endif
    endfor
  endfor
  return, 0
end

