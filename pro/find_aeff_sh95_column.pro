; docformat = 'rst'

function find_aeff_sh95_column, lo_lev, hi_lev, lev_num
;+
;     This function locates and returns the data location 
;     of the given low energy level, high energy level,  
;     and the level number
;     within the database of H I emissivities given by
;     from Storey & Hummer, 1995MNRAS.272...41S.
;
; :Private:
;   
; :Returns:
;    type=double. This function returns the data location .
;           
; :Params:
;     lo_lev  :  in, required, type=float
;                   low energy level
;
;     hi_lev  :  in, required, type=float
;                   high energy level
;
;     lev_num :  in, required, type=float
;                   level number
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.0.1
;
; :History:
;     Based on H I emissivities
;     from Storey & Hummer, 1995MNRAS.272...41S.
;     
;     25/08/2012, A. Danehkar, IDL code written.
;     
;     11/03/2017, A. Danehkar, Integration with AtomNeb.
;-
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

