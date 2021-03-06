; docformat = 'rst'

function calc_abund_he_ii_rl, temperature=temperature, density=density, $
                      line_flux=line_flux, $
                      he_ii_aeff_data=he_ii_aeff_data, h_i_aeff_data=h_i_aeff_data
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for the He II recombination line 4686 A 
;     by using the helium emissivities from 
;     Storey & Hummer, 1995MNRAS.272...41S.
;
; :Returns:
;    type=double. This function returns the ionic abundanc.
;
; :Keywords:
;     temperature     :   in, required, type=float
;                         electron temperature
;     density         :   in, required, type=float
;                         electron density
;     line_flux       :   in, required, type=float
;                         line flux intensity
;     he_ii_aeff_data :   in, required, type=array/object
;                         He II recombination coefficients
;     h_i_aeff_data   :   in, required, type=array/object
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
;     IDL> ion='iii' ; He II
;     IDL> he_ii_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> he_ii_aeff_data=he_ii_rc_data[0].Aeff
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> he_ii_4686_flux = 135.833
;     IDL> Abund_he_ii=calc_abund_he_ii_rl(temperature=temperature, density=density, $
;     IDL>                                 line_flux=he_ii_4686_flux, $
;     IDL>                                 he_ii_aeff_data=he_ii_aeff_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(He^2+)/N(H^+):', Abund_he_ii
;        N(He^2+)/N(H^+):      0.11228817
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
  if keyword_set(h_i_aeff_data) eq 0 then begin 
    print,'H I recombination coefficients (h_i_aeff_data) are not set'
    return, 0
  endif
  if keyword_set(line_flux) eq 0 then begin 
    print,'Line flux intensity (line_flux) is not given'
    return, 0
  endif  
  if (temperature le 0.D0) or (density le 0.D0) then begin
      print,'temperature = ', temperature, ', density = ', density
      return, 0
  endif
  
  emissivity_Hbeta=calc_emiss_h_beta(temperature=temperature,density=density,h_i_aeff_data=h_i_aeff_data)
  
  emissivity=calc_emiss_he_ii_rl(temperature=temperature, density=density, $
                      he_ii_aeff_data=he_ii_aeff_data)
  abund = (emissivity_Hbeta/emissivity)*double(line_flux/100.0)
  
  return,abund
end
