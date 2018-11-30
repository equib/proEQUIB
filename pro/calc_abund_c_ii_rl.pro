; docformat = 'rst'

function calc_abund_c_ii_rl, temperature=temperature, density=density, $
                      wavelength=wavelength, line_flux=line_flux, $
                      c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of C II recombination line 
;     by using the recombination coefficients from
;     from Davey et al. (2000) 2000A&AS..142...85D.
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
;     c_ii_rc_data  :     in, required, type=array/object
;                         C II recombination coefficients
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
;     IDL> atom='c'
;     IDL> ion='iii' ; C II
;     IDL> c_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> c_ii_6151_flux = 0.028
;     IDL> wavelength=6151.43
;     IDL> Abund_c_ii=calc_abund_c_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, line_flux=c_ii_6151_flux, $
;     IDL>                               c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(C^2+)/N(H+):', Abund_c_ii
;        N(C^2+)/N(H+):    0.00063404650
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
;     Based on recombination coefficients for C II lines from
;     Davey et al. 2000A&AS..142...85D.
;     
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.
;     
;     02/2003, Yong Zhang, added to MOCASSIN.
;     
;     10/05/2013, A. Danehkar, Translated to IDL code.
;     
;     15/04/2017, A. Danehkar, Integration with AtomNeb.
;-

;+
; NAME:
;     calc_abund_c_ii_rl
;
; PURPOSE:
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of C II recombination line 
;     by using the recombination coefficients from
;     from Davey et al. (2000) 2000A&AS..142...85D.
;
; CALLING SEQUENCE:
;     Result = calc_abund_c_ii_rl(TEMPERATURE=temperature, DENSITY=density, $
;                                 WAVELENGTH=wavelength, LINE_FLUX=line_flux, $
;                                 C_II_RC_DATA=c_ii_rc_data, H_I_AEFF_DATA=h_i_aeff_data)
;
; KEYWORD PARAMETERS:
;     TEMPERATURE   :     in, required, type=float, electron temperature
;     DENSITY       :     in, required, type=float, electron density
;     WAVELENGTH    :     in, required, type=float
;                         Line Wavelength in Angstrom
;     LINE_FLUX     :     in, required, type=float, line flux intensity
;     C_II_RC_DATA  :     in, required, type=array/object, C II recombination coefficients
;     H_I_AEFF_DATA :     in, required, type=array/object, H I recombination coefficients
;
; OUTPUTS:  This function returns a double as the ionic abundance.
;
; PROCEDURE: This function calls gamma_hb_4861.
;
; EXAMPLE:
;     base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     data_rc_dir = ['atomic-data-rc']
;     Atom_RC_All_file= filepath('rc_collection.fits', root_dir=base_dir, subdir=data_rc_dir )
;     Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     atom='h'
;     ion='ii' ; H I
;     h_i_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     h_i_aeff_data=h_i_rc_data[0].Aeff
;     atom='c'
;     ion='iii' ; C II
;     c_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     temperature=double(10000.0)
;     density=double(5000.0)
;     c_ii_6151_flux = 0.028
;     wavelength=6151.43
;     Abund_c_ii=calc_abund_c_ii_rl(temperature=temperature, density=density, $
;                                   wavelength=wavelength, line_flux=c_ii_6151_flux, $
;                                   c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
;     print, 'N(C^2+)/N(H+):', Abund_c_ii
;     > N(C^2+)/N(H+):    0.00063404650
; 
; MODIFICATION HISTORY:
;     Based on recombination coefficients for C II lines from
;     Davey et al. 2000A&AS..142...85D.
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.
;     02/2003, Yong Zhang, added to MOCASSIN.
;     10/05/2013, A. Danehkar, Translated to IDL code.
;     15/04/2017, A. Danehkar, Integration with AtomNeb.
;- 

  ; ciiRLstructure ={Wave:double(0.0), Int:double(0.0), Obs:double(0.0), Abundance:double(0.0)}

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  TEh2=double(temperature)
  NEh2=double(density)
  abund=1.0
  nlines = 57 
  hbeta_aeff= (10.0^gamma_hb_4861(temperature=TEh2,density=NEh2,h_i_aeff_data=h_i_aeff_data))*double(4861.33/(h_Planck*c_Speed*1.e8)) 
  
  ; ciiRLs=REPLICATE(ciiRLstructure, nlines)
  
  lamb=double(0.0)
  a=double(0.0)
  b=double(0.0)
  c=double(0.0)
  d=double(0.0)
  f=double(0.0)
  aeff=double(0.0)
  br=double(1.0)
  temp4 = temperature/10000.0
  loc1=where(abs(c_ii_rc_data.Wavelength-wavelength) le 1.5)
  temp2=size(loc1,/DIMENSIONS)
  if temp2[0] ne 1 then begin
    Wavelength_min=min(c_ii_rc_data[loc1].Wavelength)
    loc1=where(c_ii_rc_data.Wavelength eq  Wavelength_min)
  endif
  lamb=c_ii_rc_data[loc1].Wavelength
  a=c_ii_rc_data[loc1].a
  b=c_ii_rc_data[loc1].b
  c=c_ii_rc_data[loc1].c
  d=c_ii_rc_data[loc1].d
  f=c_ii_rc_data[loc1].f
  aeff = 1.0e-14 * (a*(temp4^f)) 
  aeff=aeff*(1. + (b*(1.-temp4)) + (c*((1.-temp4)^2) ) + (d * ((1.-temp4)^3) ) ) 
  ciiRLs_Int = 100.0*(aeff/hbeta_aeff)*br*(4861.33/lamb)*abund 
  
  abund=line_flux/ciiRLs_Int
  return,abund
end
