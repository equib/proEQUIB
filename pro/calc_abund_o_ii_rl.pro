; docformat = 'rst'

function calc_abund_o_ii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, line_flux=line_flux, $
                      o_ii_rc_br=o_ii_rc_br, o_ii_rc_data=o_ii_rc_data, $
                      h_i_aeff_data=h_i_aeff_data
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of O II recombination line 
;     by using the recombination coefficients from  
;     Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.
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
;     o_ii_rc_br    :     in, required, type=array/object
;                         O II branching ratios (Br)
;     o_ii_rc_data  :     in, required, type=array/object
;                         O II recombination coefficients
;     h_i_aeff_data :     in, required, type=array/object
;                         H I recombination coefficients
;
; :Examples:
;    For example::
;
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> atom='h'
;     IDL> ion='ii' ; H I
;     IDL> h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> h_i_aeff_data=h_i_rc_data[0].Aeff
;     IDL> atom='o'
;     IDL> ion='iii' ; O II
;     IDL> o_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> o_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> o_ii_4614_flux = 0.009
;     IDL> wavelength=4613.68
;     IDL> Abund_o_ii=calc_abund_o_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, line_flux=o_ii_4614_flux, $
;     IDL>                               o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data, $
;     IDL>                               h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(O^2+)/N(H+):', Abund_o_ii
;        N(O^2+)/N(H+):    0.0018886330
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
;     Based on recombination coefficients for O II lines from
;     Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.
;     
;     Adopted from MIDAS script Roii.prg written by X.W.Liu.
;     
;     Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
;                       Ercolano et al. 2005MNRAS.362.1038E.
;     
;     10/05/2013, A. Danehkar, Translated to IDL code.
;     
;     25/04/2017, A. Danehkar, Integration with AtomNeb.
;     
;     10/07/2019, A. Danehkar, Made a new function calc_emiss_o_ii_rl()
;                      for calculating line emissivities and separated it
;                      from calc_abund_o_ii_rl().
;-

  ;  oiiRLstructure ={Wave:double(0.0), $ ;REAL*8
  ;              Int:double(0.0),  $  
  ;              Obs:double(0.0), $ 
  ;              abundance:double(0.0), $             
  ;              g1:long(0), $ ;INTEGER 
  ;              g2:long(0), $ ;INTEGER
  ;              Mult1:'', $ ;CHARACTER*7
  ;              Term1:'', $ ;CHARACTER*9
  ;              Term2:'' $ ;CHARACTER*9
  ;              }
  
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(o_ii_rc_data) eq 0 then begin 
    print,'O II recombination coefficients (o_ii_rc_data) are not set'
    return, 0
  endif
  if keyword_set(o_ii_rc_br) eq 0 then begin 
    print,'O II branching ratios (o_ii_rc_br) are not set'
    return, 0
  endif
  if keyword_set(h_i_aeff_data) eq 0 then begin 
    print,'H I recombination coefficients (h_i_aeff_data) are not set'
    return, 0
  endif
  if keyword_set(wavelength) eq 0 then begin 
    print,'Wavelength is not given'
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
  emissivity=calc_emiss_o_ii_rl(temperature=temperature, density=density, wavelength=wavelength, $
                                o_ii_rc_br=o_ii_rc_br, o_ii_rc_data=o_ii_rc_data)
  abund = (emissivity_Hbeta/emissivity)*double(line_flux/100.0)
  
  return,abund
end
