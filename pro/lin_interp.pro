; docformat = 'rst'

function lin_interp, vv, xx, xout
;+
;     This function perfoms a linear interpolation/extrapolaton.
;
; :Private:
;   
; :Returns:
;    type=double. This function returns the interpolated/extrapolated value.
;           
; :Params:
;     vv  :  in, required, type=float
;                   VV array to interpolate
;
;     xx  :  in, required, type=float
;                   X array that correspond to x(0), x(1), ...
;
;     xout :  in, required, type=float
;                   X values at which vv should be interpolated
;-
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

