; docformat = 'rst'

function calc_abundance, temperature=temperature, density=density, $
                         line_flux=line_flux, atomic_levels=atomic_levels, $
                         elj_data=elj_data, omij_data=omij_data, $
                         aij_data=aij_data, h_i_aeff_data=h_i_aeff_data              
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for specified ion with level(s)
;     by solving atomic level populations and 
;     line emissivities in statistical equilibrium 
;     for given electron density and temperature.
;   
; :Returns:
;    type=double. This function returns the ionic abundanc.
;
; :Keywords:
;     temperature   :     in, required, type=float 
;                         electron temperature
;     density       :     in, required, type=float 
;                         electron density
;     line_flux     :     in, required, type=float 
;                         line flux intensity
;     atomic_levels :     in, required, type=string 
;                         level(s) e.g '1,2/', '1,2,1,3/'
;     elj_data      :     in, required, type=array/object 
;                         energy levels (Ej) data
;     omij_data     :     in, required, type=array/object 
;                         collision strengths (omega_ij) data
;     aij_data      :     in, required, type=array/object 
;                         transition probabilities (Aij) data
;     h_i_aeff_data :     in, required, type=array/object 
;                         H I recombination coefficients
;
; :Examples:
;    For example::
;    
;     IDL> base_dir = file_dirname(file_dirname((routine_info('$MAIN$', /source)).path))
;     IDL> data_dir = ['atomic-data', 'chianti70']
;     IDL> Atom_Elj_file = filepath('AtomElj.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> Atom_Omij_file = filepath('AtomOmij.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> Atom_Aij_file = filepath('AtomAij.fits', root_dir=base_dir, subdir=data_dir )
;     IDL> data_rc_dir = ['atomic-data-rc']
;     IDL> Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=base_dir, subdir=data_rc_dir )
;     IDL> atom='o'
;     IDL> ion='iii'
;     IDL> o_iii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej)
;     IDL> o_iii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
;     IDL> o_iii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)
;     IDL> atom='h'
;     IDL> ion='ii' ; H I
;     IDL> hi_rc_data=atomneb_read_aeff_sh95(Atom_RC_SH95_file, atom, ion)
;     IDL> h_i_aeff_data=hi_rc_data[0].Aeff
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> atomic_levels='3,4/'
;     IDL> iobs5007=double(1200.0)
;     IDL> Abb5007=double(0.0)
;     IDL> Abb5007=calc_abundance(temperature=temperature, density=density, $
;     IDL>                        line_flux=iobs5007, atomic_levels=atomic_levels,$
;     IDL>                        elj_data=o_iii_elj, omij_data=o_iii_omij, $
;     IDL>                        aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data[0].Aeff)
;     IDL> print, 'N(O^2+)/N(H+):', Abb5007
;        N(O^2+)/N(H+):   0.00041256231                   
;
; :Categories:
;   Abundance Analysis, Collisionally Excited Lines
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
;     15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code. 
;     
;     20/10/2016, A. Danehkar, Replaced str2int with strnumber. 
;     
;     20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
;          IDL function INTERPOL( /SPLINE). 
;        
;     20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function 
;                       LA_LINEAR_EQUATION. 
;               
;     15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
;           with IDL function LUDC & LUSOL. 
;     
;     19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with 
;                    SPL_INIT & SPL_INTERP. 
;     
;     20/11/2016, A. Danehkar, Made a new function calc_populations() 
;       for solving atomic level populations and separated it from
;       calc_abundance(), calc_density() and calc_temperature(). 
;                         
;     21/11/2016, A. Danehkar, Made a new function calc_emissivity() 
;                      for calculating line emissivities and separated it 
;                      from calc_abundance(). 
;     
;     10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
;                      input elj_data, omij_data, aij_data. 
;     
;     12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
;                        from calc_abundance(). 
;        
; FORTRAN HISTORY:
; 
;     03/05/1981, I.D.Howarth,  Version 1.
; 
;     05/05/1981, I.D.Howarth,  Minibug fixed!
; 
;     07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
; 
;     03/08/1981, S.Adams,      Interpolates collision strengths.
; 
;     07/08/1981, S.Adams,      Input method changed.
; 
;     19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
;                               filenames given to SA's data files.
; 
;     08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
; 
;     02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
;                         modified such that matrix sizes (i.e. maximum
;                         of Te and maximum no of levels) can now be cha
;                         by modifying the parameters NDIM1, NDIM2 and N
;                         in the Main program. EASY!
;                         Now takes collision rates as well.
;                         All variables are declared explicitly
;                         Generate two extra files (ionpop.lis and ionra
;                         of plain stream format for plotting.
; 
;     06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
;                         Fixed readin bug for IBIG=2 case.
;                         Now reads reformatted upsilons (easier to see
;                         and the 0 0 0 data end is excluded for these c
;                         The A values have a different format for IBIG=.
; 
;     2006, B.Ercolano,   Converted to F90.
;
;     2009, B.Wesson,     Misc updates and improvements. 
;                         Converted to F90. Version written only for    
;                         calculating ionic abundances. Takes arguments
;                         from the command line.
;-

  AHB=double(0)
  
  h_Planck = 6.62606957e-27 ; erg.s
  c_Speed = 2.99792458e10 ; cm/s 
  
  if keyword_set(temperature) eq 0 then begin 
    print,'Temperature is not set'
    return, 0
  endif
  if keyword_set(density) eq 0 then begin 
    print,'Density is not set'
    return, 0
  endif
  if keyword_set(elj_data) eq 0 then begin 
    print,'Energy Levels data (elj_data) are not set'
    return, 0
  endif
  if keyword_set(omij_data) eq 0 then begin 
    print,'Collision Strengths (omij_data) are not set'
    return, 0
  endif
  if keyword_set(aij_data) eq 0 then begin 
    print,'Transition Probabilities (aij_data) are not set'
    return, 0
  endif
  if keyword_set(h_i_aeff_data) eq 0 then begin 
    print,'H I recombination coefficients (h_i_aeff_data) are not set'
    return, 0
  endif
  if keyword_set(atomic_levels) eq 0 then begin 
    print,'Atomic levels (atomic_levels) are not given'
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
  
;  T4=temperature*1.0D-4
;  AHB=3.036D-14*T4^(-0.87D0) ; Brocklehurt 1971; Aller (1984), Physics of Thermal Gaseous Nebulae, p. 76
;  WAVHB=4861.33D ;4861.D0
;  emissivity_Hbeta=AHB*h_Planck*c_Speed*1.e8/WAVHB ; N(H+) * N(e-) (erg/s) 
  ; emissivity_Hbeta=1.387D-25*T4^(-0.983D0)* 10.D0^(-0.0424D0/T4) ;  Brocklehurst (1971); Aller (1984)
  emissivity_Hbeta=calc_emiss_h_beta(temperature=temperature,density=density,h_i_aeff_data=h_i_aeff_data)
  
  emissivity_all=double(0.0)
  emissivity_all=calc_emissivity(temperature=temperature, density=density, $
                         atomic_levels=atomic_levels, $
                         elj_data=elj_data, omij_data=omij_data, $
                         aij_data=aij_data)
  
  if emissivity_all eq 0 then begin
    print,"cannot calculate emissivity" 
    return, 0
  endif
  abund = (emissivity_Hbeta/emissivity_all)*(line_flux/100.0)
  return, abund
end
