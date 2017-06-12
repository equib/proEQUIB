function calc_abundance, temperature=temperature, density=density, $
                         line_flux=line_flux, atomic_levels=atomic_levels, $
                         elj_data=elj_data, omij_data=omij_data, $
                         aij_data=aij_data, h_i_aeff_data=h_i_aeff_data
                        
;+
; NAME:
;     calc_abundance
; PURPOSE:
;     determine ionic abundance from observed 
;     flux intensity for specified ion with level(s)
;     by solving atomic level populations and 
;     line emissivities in statistical equilibrium 
;     for given electron density and temperature.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     path='proEQUIB/atomic-data/'
;     set_atomic_data_path, path
;
;     ion='oiii'
;     tempi=double(10000.0)
;     densi=double(5000.0)
;     levels5007='3,4/'
;     iobs5007=double(1200.0)
;     Abb5007=double(0.0) 
;     Abb5007=calc_abundance(ion, levels5007, tempi, densi, iobs5007)
;     print, Abb5007
;
; INPUTS:
;     temperature   -     electron temperature
;     density       -     electron density
;     line_flux     -     line flux intensity
;     atomic_levels -     level(s) e.g '1,2/', '1,2,1,3/'
;     elj_data      -     energy levels (Ej) data
;     omij_data     -     collision strengths (omega_ij) data
;     aij_data      -     transition probabilities (Aij) data
;     h_i_aeff_data -     H I recombination coefficients
;     
; RETURN:  ionic abundance
;
; REVISION HISTORY:
;     Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
;     Replaced str2int with strnumber,      A. Danehkar, 20/10/2016
;     Replaced CFY, SPLMAT, and CFD with
;          IDL function INTERPOL( /SPLINE), A. Danehkar, 20/10/2016
;     Replaced LUSLV with IDL LAPACK function 
;                       LA_LINEAR_EQUATION, A. Danehkar, 20/10/2016
;     Replaced LA_LINEAR_EQUATION (not work in GDL)
;           with IDL function LUDC & LUSOL, A. Danehkar, 15/11/2016
;     Replaced INTERPOL (not accurate) with 
;                    SPL_INIT & SPL_INTERP, A. Danehkar, 19/11/2016
;     Made a new function getpopulations() for solving atomic 
;       level populations and separated it from
;       calc_abundance() and do_diagnostic(), A. Danehkar, 20/11/2016
;     Made a new function calc_emissivity() for calculating 
;                      line emissivities and separated it 
;                      from calc_abundance(), A. Danehkar, 21/11/2016
;     Integration with AtomNeb, now uses atomic data input elj_data,
;                      omij_data, aij_data, A. Danehkar, 10/03/2017
;     Cleaning the function, and remove unused varibales
;                        calc_emissivity(), A. Danehkar, 12/06/2017
; 
; FORTRAN EQUIB HISTORY (F77/F90):
; 1981-05-03 I.D.Howarth  Version 1
; 1981-05-05 I.D.Howarth  Minibug fixed!
; 1981-05-07 I.D.Howarth  Now takes collision rates or strengths
; 1981-08-03 S.Adams      Interpolates collision strengths
; 1981-08-07 S.Adams      Input method changed
; 1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
;                         filenames given to SA's data files.
; 1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
; 1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
;                         modified such that matrix sizes (i.e. maximum
;                         of Te and maximum no of levels) can now be cha
;                         by modifying the parameters NDIM1, NDIM2 and N
;                         in the Main program. EASY!
;                         Now takes collision rates as well.
;                         All variables are declared explicitly
;                         Generate two extra files (ionpop.lis and ionra
;                         of plain stream format for plotting
; 1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
;                         Fixed readin bug for IBIG=2 case.
;                         Now reads reformatted upsilons (easier to see
;                         and the 0 0 0 data end is excluded for these c
;                         The A values have a different format for IBIG=
; 2006       B.Ercolano   Converted to F90
;- 

  AHB=double(0)
  
  h_Planck = 6.62606957e-27 ; erg s
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
  
  ; T4=TEMP*1.0D-4
  ; AHB=3.036D-14*T4^(-0.87D0) ; Brocklehurt 1971; Aller (1984), Physics of Thermal Gaseous Nebulae, p. 76
  ; WAVHB=4861.33D ;4861.D0
  ; emissivity_Hbeta=AHB*h_Planck*c_Speed*1.e8/WAVHB ; N(H+) * N(e-) (erg/s) 
  ; emissivity_Hbeta=1.387D-25*T4^(-0.983D0)* 10.D0^(-0.0424D0/T4) ;  Brocklehurst (1971); Aller (1984)
  emissivity_Hbeta=10.0^gamma4861(h_i_aeff_data, temperature, density)
  
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
