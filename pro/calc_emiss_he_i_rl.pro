; docformat = 'rst'

function calc_emiss_he_i_rl, temperature=temperature, density=density, $
                      linenum=linenum, $
                      he_i_aeff_data=he_i_aeff_data
;+
;     This function calculates the emissivity 
;     for the given wavelength of He I recombination line 
;     by using the recombination coefficients from Porter et al. 
;     2012MNRAS.425L..28P.
;
; :Returns:
;    type=double. This function returns the line emissivity.
;
; :Keywords:
;     temperature    :    in, required, type=float
;                         electron temperature
;     density        :    in, required, type=float
;                         electron density
;     linenum        :    in, required, type=int
;                         Line Number for Wavelength
;                         
;                         Wavelength=4120.84:linenum=7,  
;                         
;                         Wavelength=4387.93: linenum=8, 
;                         
;                         Wavelength=4437.55: linenum=9, 
;                         
;                         Wavelength=4471.50: linenum=10, 
;                         
;                         Wavelength=4921.93: linenum=12, 
;                         
;                         Wavelength=5015.68: linenum=13, 
;                         
;                         Wavelength=5047.74: linenum=14, 
;                         
;                         Wavelength=5875.66: linenum=15, 
;                         
;                         Wavelength=6678.16: linenum=16, 
;                         
;                         Wavelength=7065.25: linenum=17, 
;                         
;                         Wavelength=7281.35: linenum=18. 
;                         
;     line_flux      :    in, required, type=float
;                         line flux intensity
;     he_i_aeff_data :    in, required, type=array/object
;                         He I recombination coefficients
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
;     IDL> ion='ii' ; He I
;     IDL> he_i_rc_data=atomneb_read_aeff_he_i_pfsd12(Atom_RC_He_I_file, atom, ion)
;     IDL> he_i_aeff_data=he_i_rc_data[0].Aeff
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> linenum=10; 4471.50
;     IDL> emiss_he_i=calc_emiss_he_i_rl(temperature=temperature, density=density, $
;                                       linenum=linenum, $
;                                       he_i_aeff_data=he_i_aeff_data)
;     IDL> print, 'Emissivity:', emiss_he_i
;        Emissivity:   6.3822830e-26
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
;     Based on improved He I emissivities in the case B
;     from Porter et al. 2012MNRAS.425L..28P
;     
;     15/12/2013, A. Danehkar, IDL code written.
;     
;     20/03/2017, A. Danehkar, Integration with AtomNeb.
;     
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_he_i_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_he_i_rl().
;-
  
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(he_i_aeff_data) eq 0 then begin 
    print,'He I recombination coefficients (he_i_aeff_data) are not set'
    return, 0
  endif
  if keyword_set(linenum) eq 0 then begin 
    print,'Line Number for Wavelength is not given'
    return, 0
  endif
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif
  
  TEh2=double(temperature)
  NEh2=double(density)
  line1=long(linenum-1)
  emissivity=double(0.0)
  
  ; hei_ems=read_porter()
  hei_ems=dblarr(21,14) ;(21,14,44)
  temp1=dblarr(46)
  
  nlines = 294 
  
  for i=0, nlines-1 do begin 
    temp1=he_i_aeff_data[*,i]
    tpos=nint((temp1[0]/1000)-5)
    npos=nint(alog10(temp1[1])-1)
    hei_ems[tpos,npos]=temp1[line1+2];temp[2:45]
  endfor

  
  ; restrict to the density & temperature ranges to 2012MNRAS.425L..28P  
  if (NEh2 lt 1.e1) then NEh2=1.e1
  if (NEh2 gt 1.e14) then NEh2=1.e14
  if (TEh2 lt 5000) then TEh2=5000.
  if (TEh2 gt 25000) then TEh2=25000.
  
  ; get logarithmic density
  dens_log=alog10(NEh2)
  
  dens_grid=double(indgen(14) + 1)
  temp_grid=double(1000*(indgen(21) + 5))

  hei_ems1=hei_ems[*,*]
  ; Bilinearly interpolate density & temperature
  ; emiss_log =_interp2d(hei_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic);, /trigrid) not work on GDL
  emiss_log=_interp_2d(hei_ems1, TEh2, dens_log, temp_grid, dens_grid)
    
  ; wavl=he_i_aeff_data_Wavelength[line1]
  emissivity= 10.d^(emiss_log)
  
  return,emissivity
end
