function lin_interp, vv, xx, xout
;     linear interpolation/extrapolaton

  ; Make a copy so we don't overwrite the input arguments.
  v = vv
  x = xx
  m = n_elements(v)  ;# of input pnts

  type = size(v, /TYPE)
  
  s = value_locate(x, xout) > 0L < (m-2) ;Subscript intervals.
  ;Linear, not regular
  case (type) of
    1: diff = v[s+1] - fix(v[s])
    12: diff = v[s+1] - long(v[s])
    13: diff = v[s+1] - long64(v[s])
    15: diff = long64(v[s+1]) - long64(v[s])
  else: diff = v[s+1] - v[s]
  endcase
  p = (xout-x[s])*diff/(x[s+1] - x[s]) + v[s]
  return, p
end
