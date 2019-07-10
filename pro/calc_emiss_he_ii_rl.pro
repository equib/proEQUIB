; docformat = 'rst'

function calc_emiss_he_ii_rl, temperature=temperature, density=density, $
                      he_ii_aeff_data=he_ii_aeff_data
;+
;     This functioncalculates the emissivity
;     for the He II recombination line 4686 A 
;     by using the helium emissivities from 
;     Storey & Hummer, 1995MNRAS.272...41S.
;
; :Returns:
;    type=double. This function returns the line emissivity.
;
; :Keywords:
;     temperature     :   in, required, type=float
;                         electron temperature
;     density         :   in, required, type=float
;                         electron density
;     he_ii_aeff_data :   in, required, type=array/object
;                         He II recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> 
;     IDL> atom='he'
;     IDL> ion='iii' ; He II
;     IDL> he_ii_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> he_ii_aeff_data=he_ii_rc_data[0].Aeff
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> he_ii_4686_flux = 135.833
;     IDL> emiss_he_ii=calc_emiss_he_ii_rl(temperature=temperature, density=density, $
;     IDL>                                 he_ii_aeff_data=he_ii_aeff_data)
;     IDL> print, 'Emissivity:', emiss_he_ii
;        Emissivity:   1.4989134e-24
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
;     Based on He II emissivities
;     from Storey & Hummer, 1995MNRAS.272...41S.
;     
;     15/12/2013, A. Danehkar, IDL code written.
;     
;     02/04/2017, A. Danehkar, Integration with AtomNeb.
;     
;     10/07/2019, A. Danehkar, Change from logarithmic to linear
;     
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_he_ii_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_he_ii_rl().
;-
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(he_ii_aeff_data) eq 0 then begin 
    print,'He II recombination coefficients (he_ii_aeff_data) are not set'
    return, 0
  endif
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif
  
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
  ; emissivity =_interp2d(heii_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic);, /trigrid) not work on GDL
  emissivity=_interp_2d(heii_ems1, TEh2, dens_log, temp_grid, dens_grid)
  ; emissivity_log = alog10(emissivity)
  
  return,emissivity
end
