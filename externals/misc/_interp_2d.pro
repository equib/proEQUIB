function _interp_2d,tab_in,x,y,x_tab,y_tab,xlog=xlog,ylog=ylog
;+
;   function interp_2d,tab_in,x,y,x_tab,y_tab,[xlog=xlog],[ylog=ylog]
;   return an interpolated 2 Dim less than tab_in (no extrapolation)
;   tab_in: fltarr([*,*,...],n_x,n_y)
;   x,y: values where we want to interpolate
;   x_tab,y_tab: fltarr(n_x),fltarr(n_y) : x and y vectors, must be
;      sorted, no necessarily in increasing order.
;   x_log and y_log: allows log interpolation
;   C. Morisset (LAS, Marseille, 2000)
;-   
   
   on_error,2
   
   
   size_tab_in =  size(tab_in)
   n_dim_tab =  size_tab_in[0]
   n_x =  n_elements(x_tab)
   n_y =  n_elements(y_tab)
   if n_x ne size_tab_in[n_dim_tab-1] or n_y ne size_tab_in[n_dim_tab] $
    then message,' Dimension of x_tab or y_tab incompatible with tab_in'
   tab_tmp = reform(tab_in,size_tab_in[n_dim_tab+2]/n_x/n_y,n_x,n_y)
   
; No extrapolation:
   if x ge max(x_tab) or x le min(x_tab) or y ge max(y_tab) $
    or y le min(y_tab) then message,' No extrapollation'
   
   i_x =  max(where(x_tab lt x))
   i_y =  max(where(y_tab lt y))
   
   if x_tab[0] gt x_tab[1] then x_incr =  -1 else x_incr = 1
   if y_tab[0] gt y_tab[1] then y_incr =  -1 else y_incr = 1
   
   if x_incr eq 1 then i_x =  max(where(x_tab lt x)) else $
    i_x =  max(where(x_tab gt x))
   if y_incr eq 1 then i_y =  max(where(y_tab lt y)) else $
    i_y =  max(where(y_tab gt y))
   
   if keyword_set(xlog) then $
    f_x =  1. - (alog10(x)-alog10(x_tab[i_x]))/ $
      (alog10(x_tab[i_x+1])-alog10(x_tab[i_x])) else $
    f_x =  1. - (x-x_tab[i_x])/(x_tab[i_x+1]-x_tab[i_x])
   if keyword_set(ylog) then $
    f_y =  1. - (alog10(y)-alog10(y_tab[i_y]))/ $
      (alog10(y_tab[i_y+1])-alog10(y_tab[i_y])) else $
    f_y =  1. - (y-y_tab[i_y])/(y_tab[i_y+1]-y_tab[i_y])
      
   tab_out =  tab_tmp[*,i_x,i_y]*f_x*f_y + $
    tab_tmp[*,i_x+1,i_y]*(1.-f_x)*f_y + $
    tab_tmp[*,i_x,i_y+1]*f_x*(1.-f_y) + $
    tab_tmp[*,i_x+1,i_y+1]*(1.-f_x)*(1.-f_y)
   
   if n_dim_tab gt 2 then $
   tab_out =  reform(tab_out,size_tab_in[1:n_dim_tab-2],/OVERWRITE)
   return,tab_out
end
