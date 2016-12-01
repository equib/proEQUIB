function calc_abundance, ion, levels, tempi, densi, iobs
;+
; NAME:
;     getabundance
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
;     Abb5007=getabundance(ion, levels5007, tempi, densi, iobs5007)
;     print, Abb5007
;
; INPUTS:
;     ion -       ion name e.g. 'oii', 'oiii'
;     levels -    level(s) e.g '1,2/', '1,2,1,3/'
;     tempi -     electron temperature
;     densi -     electron density
;     iobs -      observed flux intensity
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
;     Made a new function calc_populations() for solving atomic 
;       level populations and separated it from
;       calc_abundance(), calc_temp_dens(), A. Danehkar, 20/11/2016
;     Made a new function calc_emissivity() for calculating 
;                      line emissivities and separated it 
;                      from calc_abundance(), A. Danehkar, 21/11/2016
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
; 2009-05    R.Wesson     Misc updates and improvements, inputs from cmd line, 
;                         written only for calculating ionic abundances.
;- 

  AHB=double(0)
  
  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  TEMP=TEMPI
  ; Set Ne
  DENS=DENSI
  if (TEMP le 0.D0) or (DENS le 0.D0) then begin
      print,'Temp = ', TEMP, ', Dens = ', DENS
      return, 0
  endif
  
  ; Output data
  ; Eff. recombination coef. of Hb
  T4=TEMP*1.0D-4
  AHB=3.036D-14*T4^(-0.87D0) ; Brocklehurt 1971; Aller (1984), Physics of Thermal Gaseous Nebulae, p. 76
  WAVHB=4861.33D ;4861.D0
  emissivity_Hbeta=AHB*h_Planck*c_Speed*1.e8/WAVHB ; N(H+) * N(e-) (erg/s) 
  ; emissivity_Hbeta=1.387D-25*T4^(-0.983D0)* 10.D0^(-0.0424D0/T4) ;  Brocklehurst (1971); Aller (1984)
  
  emissivity_all=double(0.0)
  emissivity_all=calc_emissivity(ion, levels, tempi, densi)
  
  abund = (emissivity_Hbeta/emissivity_all)*(iobs/100.0)
  return, abund
end
