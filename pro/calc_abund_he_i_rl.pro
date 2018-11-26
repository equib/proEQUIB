; docformat = 'rst'

function calc_abund_he_i_rl, temperature=temperature, density=density, $
                      linenum=linenum, line_flux=line_flux, $
                      he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of He I recombination line 
;     by using the recombination coefficients from Porter et al. 
;     2012MNRAS.425L..28P.
;
; :Returns:
;    type=double. This function returns the ionic abundanc.
;
; :Keywords:
;     temperature    :    in, required, type=float
;                         electron temperature
;     density        :    in, required, type=float
;                         electron density
;     linenum        :    in, required, type=int
;                         Line Number for Wavelength
;                         Wavelength=4120.84:linenum=7: 
;                         Wavelength=4387.93: linenum=8
;                         Wavelength=4437.55: linenum=9
;                         Wavelength=4471.50: linenum=10
;                         Wavelength=4921.93: linenum=12
;                         Wavelength=5015.68: linenum=13
;                         Wavelength=5047.74: linenum=14
;                         Wavelength=5875.66: linenum=15
;                         Wavelength=6678.16: linenum=16
;                         Wavelength=7065.25: linenum=17
;                         Wavelength=7281.35: linenum=18
;     line_flux      :    in, required, type=float
;                         line flux intensity
;     he_i_aeff_data :    in, required, type=array/object
;                         He I recombination coefficients
;     h_i_aeff_data  :    in, required, type=array/object
;                         H I recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> atom='h'
;     IDL> ion='ii' ; H I
;     IDL> h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> h_i_aeff_data=h_i_rc_data[0].Aeff
;     IDL> atom='he'
;     IDL> ion='ii' ; He I
;     IDL> he_i_rc_data=atomneb_read_aeff_he_i_pfsd12(Atom_RC_He_I_file, atom, ion)
;     IDL> he_i_aeff_data=he_i_rc_data[0].Aeff
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> he_i_4471_flux= 2.104
;     IDL> linenum=10; 4471.50
;     IDL> Abund_he_i=calc_abund_he_i_rl(temperature=temperature, density=density, $
;                                       linenum=linenum, line_flux=he_i_4471_flux, $
;                                       he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(He^+)/N(H^+):', Abund_he_i
;        N(He^+)/N(H^+):     0.040722892
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
;   0.0.3
;
; :History:
;     Based on improved He I emissivities in the case B
;     from Porter et al. 2012MNRAS.425L..28P
;     
;     15/12/2013, A. Danehkar, IDL code written.
;     
;     20/03/2017, A. Danehkar, Integration with AtomNeb.
;-

;+
; NAME:
;     calc_abund_he_i_rl
; 
; PURPOSE:
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of He I recombination line 
;     by using the recombination coefficients from Porter et al. 
;     2012MNRAS.425L..28P.
;
; CALLING SEQUENCE:
;     Result = calc_abund_he_i_rl(TEMPERATURE=temperature, DENSITY=density, $
;                                LINENUM=linenum, LINE_FLUX=line_flux, $
;                                HE_I_AEFF_DATA=he_i_aeff_data, H_I_AEFF_DATA=h_i_aeff_data)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE   :     in, required, type=float, electron temperature
;     DENSITY       :     in, required, type=float, electron density
;     LINENUM       :     in, required, type=int
;                         Line Number for Wavelength
;                         Wavelength=4120.84:linenum=7: 
;                         Wavelength=4387.93: linenum=8
;                         Wavelength=4437.55: linenum=9
;                         Wavelength=4471.50: linenum=10
;                         Wavelength=4921.93: linenum=12
;                         Wavelength=5015.68: linenum=13
;                         Wavelength=5047.74: linenum=14
;                         Wavelength=5875.66: linenum=15
;                         Wavelength=6678.16: linenum=16
;                         Wavelength=7065.25: linenum=17
;                         Wavelength=7281.35: linenum=18
;     LINE_FLUX     :     in, required, type=float, line flux intensity
;     HE_I_AEFF_DATA:     in, required, type=array/object, He I recombination coefficients
;     H_I_AEFF_DATA :     in, required, type=array/object, H I recombination coefficients
;
; OUTPUTS:  This function returns a double as the ionic abundance.
;
; PROCEDURE: This function calls gamma_hb_4861.
;
; EXAMPLE:
;     base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     data_rc_dir = ['atomic-data-rc']
;     Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
;     Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     atom='h'
;     ion='ii' ; H I
;     h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     h_i_aeff_data=h_i_rc_data[0].Aeff
;     atom='he'
;     ion='ii' ; He I
;     he_i_rc_data=atomneb_read_aeff_he_i_pfsd12(Atom_RC_He_I_file, atom, ion)
;     he_i_aeff_data=he_i_rc_data[0].Aeff
;     temperature=double(10000.0)
;     density=double(5000.0)
;     he_i_4471_flux= 2.104
;     linenum=10; 4471.50
;     Abund_he_i=calc_abund_he_i_rl(temperature=temperature, density=density, $
;                                   linenum=linenum, line_flux=he_i_4471_flux, $
;                                   he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
;     print, 'N(He^+)/N(H^+):', Abund_he_i
;     > N(He^+)/N(H^+):     0.040722892
; 
; MODIFICATION HISTORY:
;     Based on improved He I emissivities in the case B
;     from Porter et al. 2012MNRAS.425L..28P
;     15/12/2013, A. Danehkar, IDL code written.
;     20/03/2017, A. Danehkar, Integration with AtomNeb.
;- 
  
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
  emiss_log = _interp2d(hei_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)
    
  ; wavl=he_i_aeff_data_Wavelength[line1]
  emissivity= 100.*10.^(emiss_log-gamma_hb_4861(temperature=TEh2,density=NEh2,h_i_aeff_data=h_i_aeff_data))
  
  abund=line_flux/emissivity
  return,abund
end
