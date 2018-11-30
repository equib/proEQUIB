; docformat = 'rst'

function calc_abund_n_iii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, line_flux=line_flux, $
                      n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data
;
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of N III recombination line 
;     by using the recombination coefficients from
;     Pequignot et al. 1991A&A...251..680P.
;
; :Returns:
;    type=double. This function returns the ionic abundanc.
;
; :Keywords:
;     temperature   :     in, required, type=float
;                         electron temperature
;     density       :     in, required, type=float
;                         electron density
;     wavelength    :     in, required, type=float
;                         Line Wavelength in Angstrom
;     line_flux     :     in, required, type=float
;                         line flux intensity
;     n_iii_rc_data  :     in, required, type=array/object
;                         N III recombination coefficients
;     h_i_aeff_data :     in, required, type=array/object
;                         H I recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_PPB91_file='/media/linux/proEQUIB/AtomNeb-idl/atomic-data-rc/rc_PPB91.fits'
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> atom='h'
;     IDL> ion='ii' ; H I
;     IDL> h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> h_i_aeff_data=h_i_rc_data[0].Aeff
;     IDL> atom='n'
;     IDL> ion='iv' ; N III
;     IDL> n_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> n_iii_4641_flux = 0.245
;     IDL> wavelength=4640.64
;     IDL> Abund_n_iii=calc_abund_n_iii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, line_flux=n_iii_4641_flux, $
;     IDL>                                 n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(N^3+)/N(H+):', Abund_n_iii
;        N(N^3+)/N(H+):    6.3366175e-05
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
;     Based on  effective radiative recombination coefficients for N III lines from
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.
;     
;     10/05/2013, A. Danehkar, IDL code written.
;     
;     20/04/2017, A. Danehkar, Integration with AtomNeb.
;-

;+
; NAME:
;     calc_abund_n_iii_rl
; PURPOSE:
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of N III recombination line 
;     by using the recombination coefficients from
;     Pequignot et al. 1991A&A...251..680P.
;
; CALLING SEQUENCE:
;     Result = calc_abund_n_iii_rl(TEMPERATURE=temperature, DENSITY=density, $
;                                  WAVELENGTH=wavelength, LINE_FLUX=line_flux, $
;                                  N_III_RC_DATA=n_iii_rc_data, H_I_AEFF_DATA=h_i_aeff_data)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE   :     in, required, type=float, electron temperature
;     DENSITY       :     in, required, type=float, electron density
;     WAVELENGTH    :     in, required, type=float
;                         Line Wavelength in Angstrom
;     LINE_FLUX     :     in, required, type=float, line flux intensity
;     N_III_RC_DATA  :     in, required, type=array/object, N II recombination coefficients
;     H_I_AEFF_DATA :     in, required, type=array/object, H I recombination coefficients
;
; OUTPUTS:  This function returns a double as the ionic abundance.
;
; PROCEDURE: This function calls gamma_hb_4861.
;
; EXAMPLE:
;     base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     data_rc_dir = ['atomic-data-rc']
;     Atom_RC_PPB91_file='/media/linux/proEQUIB/AtomNeb-idl/atomic-data-rc/rc_PPB91.fits'
;     Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     atom='h'
;     ion='ii' ; H I
;     h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     h_i_aeff_data=h_i_rc_data[0].Aeff
;     atom='n'
;     ion='iv' ; N III
;     n_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)
;     temperature=double(10000.0)
;     density=double(5000.0)
;     n_iii_4641_flux = 0.245
;     wavelength=4640.64
;     Abund_n_iii=calc_abund_n_iii_rl(temperature=temperature, density=density, $
;                                     wavelength=wavelength, line_flux=n_iii_4641_flux, $
;                                     n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
;     print, 'N(N^3+)/N(H+):', Abund_n_iii
;     > N(N^3+)/N(H+):    6.3366175e-05
; 
; MODIFICATION HISTORY:
;     Based on  effective radiative recombination coefficients for N III lines from
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.
;     10/05/2013, A. Danehkar, IDL code written.
;     20/04/2017, A. Danehkar, Integration with AtomNeb.
;- 
  
  ; niiiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  TEh2=double(temperature)
  NEh2=double(density)
  abund=1.0

  nlines = 2
  hbeta_aeff= (10.0^gamma_hb_4861(temperature=TEh2,density=NEh2,h_i_aeff_data=h_i_aeff_data))*double(4861.33/(h_Planck*c_Speed*1.e8)) 
  
  ; niiiRLs=REPLICATE(niiiRLstructure, nlines)

  lamb=double(0.0)
  a=double(0.0)
  b=double(0.0)
  c=double(0.0)
  d=double(0.0)
  br=double(0.0)
  aeff=double(0.0)
  Ion=''
  
  z = 3.0 ; ion level c^3+
  ; equation (1) in 1991A&A...251..680P
  temp4 = 1.0e-4 * temperature /z^2
  loc1=where(abs(n_iii_rc_data.Wavelength-wavelength) le 0.01)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(n_iii_rc_data[loc1].Wavelength)
    loc1=where(n_iii_rc_data.Wavelength eq  Wavelength_min)
    temp2=size(loc1,/DIMENSIONS)
    if temp2[0] ne 1 then begin
      loc1=where(n_iii_rc_data.Wavelength eq  Wavelength_min and n_iii_rc_data.Case1 eq 'B')
    endif
  endif
  lamb=n_iii_rc_data[loc1].Wavelength
  a=n_iii_rc_data[loc1].a
  b=n_iii_rc_data[loc1].b
  c=n_iii_rc_data[loc1].c
  d=n_iii_rc_data[loc1].d
  br=n_iii_rc_data[loc1].br
  ; equation (1) in 1991A&A...251..680P
  aeff = 1.0e-13 * z * br
  aeff=aeff*(a*(temp4^b))/(1.+c*(temp4^d)) 
  niiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund 
  
  abund=line_flux/niiiRLs_Int
  return,abund
end
