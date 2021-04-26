; docformat = 'rst'

function calc_emiss_h_alpha, temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data
;+
;     This function calculates the emissivity for H_alpha 6563A
;     Emis(HAlpha)= 4pi j(HAlpha 6563 A)/Np Ne) for the given temperature and density 
;     by using the hydrogen emissivities from 
;     Storey & Hummer, 1995MNRAS.272...41S.
;
; :Private:
;  
; :Returns:
;    type=double. This function returns the H beta emissivity 4pi j(HBeta 4861)/Np Ne).
;
; :Keywords:
;     temperature     :   in, required, type=float
;                         electron temperature
;     density         :   in, required, type=float
;                         electron density
;     h_i_aeff_data   :   in, required, type=array/object
;                         H I recombination coefficients
;
; :Categories:
;   Abundance Analysis, Recombination Lines, Emissivity
;
; :Dirs:
;  ./
;      Main routines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.3.0
;
; :History:
;     Based on H I emissivities
;     from Storey & Hummer, 1995MNRAS.272...41S.
;     
;     25/08/2012, A. Danehkar, IDL code written.
;     
;     11/03/2017, A. Danehkar, Integration with AtomNeb.
;
;     10/07/2019, A. Danehkar, Change from logarithmic to linear
;-

  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(h_i_aeff_data) eq 0 then begin 
    print,'H I recombination coefficients (h_i_aeff_data) are not set'
    return, 0
  endif

  ;h_a_col= find_aeff_sh95_column(3, 2)
  linenum= find_aeff_sh95_column(3, 2, 25)

  TEh2=double(temperature)
  NEh2=double(density)
  line1=long(linenum-1)
  emissivity=double(0.0)
  
  h_i_ems=dblarr(10,13)
  temp1=dblarr(302)
  temp_grid=[500.,1000.,3000.,5000.,7500.,10000.,12500.,15000.,20000.,30000.]
  
  nlines = 130 
  
  for i=0, nlines-1 do begin 
    temp1=h_i_aeff_data[*,i]
    tpos=nint((where(temp_grid eq temp1[1])))
    npos=nint(alog10(temp1[0])-2)
    h_i_ems[tpos,npos]=temp1[line1];temp[2:45]
  endfor
  
  ; restrict to the density & temperature ranges to 1995MNRAS.272...41S
  if (NEh2 lt 1.1e2) then NEh2=1.1e2
  if (NEh2 gt 1.e14) then NEh2=1.e14
  if (TEh2 lt 550.) then TEh2=550.
  if (TEh2 gt 30000.) then TEh2=30000.
  
  ; get logarithmic density
  dens_log=alog10(NEh2)
  
  dens_grid=double(indgen(13) + 2)
  ha_emissivity=double(0.0)
  ; Bilinearly interpolate density & temperature
  ; emiss_log =_interp2d(hi_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic);, /trigrid) not work on GDL
  ha_emissivity=_interp_2d(h_i_ems[*,*], TEh2, dens_log, temp_grid, dens_grid)

  ;ems_log = alog10(hr_tmp/double(4861.33/1.98648E-08))
  ;h_beta_emissivity = 10.0^ems_log
  ;h_beta_emissivity_log=alog10(h_beta_emissivity(temp, density))
  
  return, ha_emissivity
end
