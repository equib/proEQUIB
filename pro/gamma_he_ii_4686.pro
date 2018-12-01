; docformat = 'rst'

function gamma_he_ii_4686, temperature=temperature, density=density, he_ii_aeff_data=he_ii_aeff_data
;+
;     This function determines the value of gamma(He II 4686 A) = 
;     log10(4pi j(He II 4686 A)/Np Ne) for the given temperature and density 
;     by using the helium emissivities from 
;     Storey & Hummer, 1995MNRAS.272...41S.
;
; :Private:
;  
; :Returns:
;    type=double. This function returns the value of gamma(He II 4686) = log10(4pi j(He II 4686)/Np Ne).
;
; :Keywords:
;     temperature     :   in, required, type=float
;                         electron temperature
;     density         :   in, required, type=float
;                         electron density
;     he_ii_aeff_data :   in, required, type=array/object
;                         He II recombination coefficients
;
; :Categories:
;   Abundance Analysis, Recombination Lines
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
;   0.0.2
;
; :History:
;     Based on He II emissivities
;     from Storey & Hummer, 1995MNRAS.272...41S.
;     
;     31/08/2012, A. Danehkar, IDL code written.
;     
;     02/03/2017, A. Danehkar, Integration with AtomNeb.
;-

;+
; NAME:
;     gamma_he_ii_4686
;
; PURPOSE:
;     This function determines the value of gamma(He II 4686 A) = 
;     log10(4pi j(He II 4686 A)/Np Ne) for the given temperature and density 
;     by using the helium emissivities from 
;     Storey & Hummer, 1995MNRAS.272...41S.
;
; CALLING SEQUENCE:
;     Result = gamma_he_ii_4686(temperature=temperature, density=density, he_ii_aeff_data=he_ii_aeff_data)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE   :     in, required, type=float, electron temperature
;     DENSITY       :     in, required, type=float, electron density
;     HE_II_AEFF_DATA:     in, required, type=array/object, He II recombination coefficients
;
; OUTPUTS:  This function returns a double as the value of gamma(He II 4686 A) = log10(4pi j(He II 4686 A)/Np Ne).
;
; PROCEDURE: This function is called by calc_abund_he_ii_rl.
;
; MODIFICATION HISTORY:
;     Based on He II emissivities
;     from Storey & Hummer, 1995MNRAS.272...41S.
;     31/08/2012, A. Danehkar, IDL code written.
;     02/03/2017, A. Danehkar, Integration with AtomNeb.
;- 

  ;h_a_col= find_aeff_sh95_column(3, 2)
  linenum= find_aeff_sh95_column(4, 3, 25)

  TEh2=double(temperature)
  NEh2=double(density)
  line1=long(linenum-1)
  emissivity=double(0.0)
  
  heii_ems=dblarr(12,13)
  temp1=dblarr(302)
  temp_grid=[500.,1000.,3000.,5000.,7500.,10000.,12500.,15000.,20000.,30000.,50000.,100000.]
  
  nlines = 156 
  
  for i=0, nlines-1 do begin 
    temp1=he_ii_aeff_data[*,i]
    tpos=nint((where(temp_grid eq temp1[1])))
    npos=nint(alog10(temp1[0])-2)
    heii_ems[tpos,npos]=temp1[line1];temp[2:45]
  endfor
  
  ; restrict to the density & temperature ranges to 2012MNRAS.425L..28P  
  if (NEh2 lt 1.e2) then NEh2=1.e2
  if (NEh2 gt 1.e14) then NEh2=1.e14
  if (TEh2 lt 500.) then TEh2=500.
  if (TEh2 gt 100000.) then TEh2=100000.
  
  ; get logarithmic density
  dens_log=alog10(NEh2)
  
  dens_grid=double(indgen(13) + 2)

  heii_ems1=heii_ems[*,*]
  ; Bilinearly interpolate density & temperature
  ; emiss_log =_interp2d(heii_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic);, /trigrid) not work on GDL
  emiss_log=_interp_2d(heii_ems1, TEh2, dens_log, temp_grid, dens_grid)

  ;logems = alog10(hr_tmp/double(4861.33/1.98648E-08))
  ;hb_ems = 10.0^logems
  ;hb_ems=alog10(hb_emissivity(temperature, density))

  ;hb_ems=alog10(emiss_log/double(4861.33/1.98648E-08))
  heii_ems = alog10(emiss_log)
  
  return, heii_ems
end
