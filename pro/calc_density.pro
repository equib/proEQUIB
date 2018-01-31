function calc_density, line_flux_ratio=line_flux_ratio, temperature=temperature, $
                       upper_levels=upper_levels, lower_levels=lower_levels, $
                       elj_data=elj_data, omij_data=omij_data, $
                       aij_data=aij_data
;+
; NAME:
;     calc_density
; PURPOSE:
;     determine electron density from given 
;     flux intensity ratio for specified ion with upper level(s)
;     lower level(s) by solving atomic level populations and 
;     line emissivities in statistical equilibrium 
;     for given electron temperature.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     path='proEQUIB/atomic-data/'
;     set_atomic_data_path, path
;
;     Atom_Elj_file='idllib/AtomNeb/atomic-data/chianti70/AtomElj.fits'
;     Atom_Omij_file='idllib/AtomNeb/atomic-data/chianti70/AtomOmij.fits'
;     Atom_Aij_file='idllib/AtomNeb/atomic-data/chianti70/AtomAij.fits' 
;     atom='s'
;     ion='ii'
;     s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej) 
;     s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
;     s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)\
;     upper_levels='1,2/'   
;     lower_levels='1,3/'
;     temperature=double(7000.0);
;     line_flux_ratio=double(1.506);
;     density=calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
;                          upper_levels=upper_levels, lower_levels=lower_levels, $
;                          elj_data=s_ii_elj, omij_data=s_ii_omij, $
;                          aij_data=s_ii_aij)
;     print, "Electron Density:", density
;
; INPUTS:
;     line_flux_ratio  -     flux intensity ratio
;     temperature      -     electron temperature
;     upper_levels     -      upper atomic level(s) e.g '1,2/', '1,2,1,3/'
;     lower_levels     -      lower atomic level(s) e.g '1,2/', '1,2,1,3/'
;     elj_data         -     energy levels (Ej) data
;     omij_data        -     collision strengths (omega_ij) data
;     aij_data         -     transition probabilities (Aij) data
;     
; RETURN:  electron density
; 
; REVISION HISTORY:
;     Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
;     Replaced str2int with strnumber, A. Danehkar, 20/10/2016
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
;       calc_abundance(), calc_density() and calc_temperature(), 
;                                           A. Danehkar, 20/11/2016
;     Integration with AtomNeb, now uses atomic data input elj_data,
;                      omij_data, aij_data, A. Danehkar, 10/03/2017
;     Cleaning the function, and remove unused varibales
;                           calc_density(), A. Danehkar, 12/06/2017
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
  common share1, Atomic_Data_Path
  
  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
  
  iteration= long(0)
  
  I= long(0) 
  I1= long(0) 
  I2= long(0) 
  J= long(0) 
  K= long(0) 
  L= long(0) 
  JT= long(0) 
  JJD= long(0)
  level_num= long(0) 
  temp_num= long(0) 
  IRATS= long(0) 
  INT= long(0) 
  IND= long(0) 
  IT= long(0)
  IKT= long(0) 
  IA= long(0) 
  IB= long(0) 
     
  TEMPI=double(0) 
  TINC=double(0)
  DENSI=double(0) 
  DINC=double(0)
  density=double(0)
  EJI=double(0)
  WAV=double(0)
  emis_sum_a=double(0)
  emis_sum_b=double(0)
  QX=double(0)
  AX=double(0)
  EX=double(0)
  FRAT=double(0)
  DEE=double(0)
  LTEXT = '';
  
  result1=double(0)
     
  I= long(0)
  J= long(0)
  K= long(0)
  
  temp=size(elj_data,/DIMENSIONS)
  level_num=temp[0]
  temp=size(omij_data[0].strength,/DIMENSIONS)
  temp_num=temp[0]
  temp=size(omij_data,/DIMENSIONS)
  omij_num=temp[0]
  
  Glj=lonarr(level_num)

  Nlj=dblarr(level_num)
  WAVA=dblarr(level_num+1)
  WAVB=dblarr(level_num+1)
  Omij=dblarr(temp_num,level_num,level_num)   
  Aij=dblarr(level_num,level_num)   
  Elj=dblarr(level_num)   
  temp_list=dblarr(temp_num)
  check_value=dblarr(3+1)
     
  LABEL1=STRARR(level_num+1)
  
  upper_levels_str=strsplit(upper_levels, ',', ESCAPE='/', /EXTRACT)
  lower_levels_str=strsplit(lower_levels, ',', ESCAPE='/', /EXTRACT)
  
  temp=size(upper_levels_str, /N_ELEMENTS)
  upper_levels_num=long(temp[0]/2)
  temp=size(lower_levels_str, /N_ELEMENTS)
  lower_levels_num=long(temp[0]/2)
  
  ITRANA=lonarr(2,upper_levels_num)
  ITRANB=lonarr(2,lower_levels_num)
  
  ITRANA[*,*]=0
  ITRANB[*,*]=0
  
  upper_levels_i=0
  for i=0, upper_levels_num-1 do begin 
    res=equib_strnumber(upper_levels_str[upper_levels_i], val)
    if res eq 1 then ITRANA[0,i]=long(val)
    res=equib_strnumber(upper_levels_str[upper_levels_i+1], val)
    if res eq 1 then ITRANA[1,i]=long(val)
    upper_levels_i = upper_levels_i + 2
    ;if upper_levels_i ge 2*upper_levels_num then break
  endfor

  lower_levels_i=0
  for i=0, lower_levels_num-1 do begin 
    res=equib_strnumber(lower_levels_str[lower_levels_i], val)
    if res eq 1 then ITRANB[0,i]=long(val)
    res=equib_strnumber(lower_levels_str[lower_levels_i+1], val)
    if res eq 1 then ITRANB[1,i]=long(val)
    lower_levels_i = lower_levels_i + 2
    ;if lower_levels_i ge 2*lower_levels_num then break;
  endfor
  IRATS=0
  temp_list = omij_data[0].strength
  temp_list = alog10(temp_list)
  for k = 1, omij_num-1 do begin
    I = omij_data[k].level1
    J = omij_data[k].level2
    if I le level_num and J le level_num then begin
      Omij[0:temp_num-1,I-1,J-1] = omij_data[k].strength
    endif
  endfor
  level_max=max([max(ITRANA),max(ITRANB)])
  Aij =aij_data.AIJ
  Elj =elj_data.Ej
  Glj =long(elj_data.J_v*2.+1.)
  ; set density iterations
  ; start of iterations
  ; ****************************
  for iteration = 1, 9 do begin
    if (iteration eq 1) then begin
      densi=0.0
    endif else begin
      densi=check_value[2]
    endelse
;    IND=4
;    DINC=(100000.0)/((IND-1)^(iteration))
;    IND=8
;    DINC=(1000000.0)/((IND-1)^(iteration))
    IND=12
    DINC=(1000000.0)/((IND-1)^(iteration))
    TempI=temperature
    TINC=0
    INT=1
    
    RESULTS=dblarr(3+1,IND+1)
    if (densi le 0) then densi=1
    if (tempi lt 5000) then tempi=5000 ; add
    ; Start of temperature iteration
    for JT = 1, INT do begin
      temperature=TEMPI+(JT-1)*TINC 
      ; Start of density iteration=
      for JJD = 1, IND  do begin
        density=DENSI+(JJD-1)*DINC
        if (temperature le 0.D0) or (density le 0.D0) then begin
            print,'temperature = ', temperature, ', density = ', density
            return, 0
        endif
        if level_max gt level_num then begin
          print, "error outside level range"
          retunr, 0
        endif
        Nlj=calc_populations(temperature=temperature, density=density, $
                           temp_list=temp_list, $ 
                           Omij=Omij, Aij=Aij, Elj=Elj, $
                           Glj=Glj, level_num=level_max, $
                           temp_num=temp_num, irats=irats)
        
        ; Search ITRANA, ITRANB for transitions & sum up   
        emis_sum_a=double(0.0)
        emis_sum_b=double(0.0)
        for IKT=0, upper_levels_num-1 do begin 
          I=ITRANA[0,IKT]
          J=ITRANA[1,IKT]
          emissivity_line=double(0.0)
          if (Aij[J-1,I-1] ne 0.D0) then begin
            EJI = Elj[J-1] - Elj[I-1]
            WAV = 1.D8 / EJI
            emissivity_line=Nlj[J]*Aij[J-1,I-1]*h_Planck*c_Speed*1.e8/WAV
            emis_sum_a=emis_sum_a+emissivity_line
          endif
        endfor
        for IKT=0, lower_levels_num-1 do begin 
          I=ITRANB[0,IKT]
          J=ITRANB[1,IKT]
          emissivity_line=double(0.0)
          if (Aij[J-1,I-1] ne 0.D0) then begin
            EJI = Elj[J-1] - Elj[I-1]
            WAV = 1.D8 / EJI
            emissivity_line=Nlj[J]*Aij[J-1,I-1]*h_Planck*c_Speed*1.e8/WAV
            emis_sum_b=emis_sum_b+emissivity_line
          endif
        endfor
        FRAT=emis_sum_a/emis_sum_b
        RESULTS[1, JJD] = temperature
        RESULTS[2, JJD] = density
        RESultS[3, JJD] = FRAT-line_flux_ratio
      endfor
      for IA = 0, upper_levels_num-1 do begin
        I1=ITRANA[0,IA]
        I2=ITRANA[1,IA]
        DEE=Elj[I2-1]-Elj[I1-1]
        WAVA[IA]=1.D8/DEE
      endfor
      for IB = 0, lower_levels_num-1 do begin
        I1=ITRANB[0,IB]
        I2=ITRANB[1,IB]
        DEE=Elj[I2-1]-Elj[I1-1]
        WAVB[IB]=1.D8/DEE
      endfor
    ; End of the temperature iteration
    endfor
    INT = ind
    ; iteration and detect the sign change.
    for I=2,INT do begin
      check=0
      if (equib_sign(results[3,I],results[3,1]) ne results[3,I]) then begin 
        ;if this condition, the values have a different sign
        check_value[*] = results[*,I-1] ; the value before the sign change returned
        check=1
        break
      endif
    endfor
    if(check eq 0) and (iteration lt 9) then begin ; check if there is any change of sign,
                             ;and checks if it should be upper or lower limit
      if(abs(results[3,1])) lt (abs(results[3,INT])) then begin
          check_value[*]=results[*,1]
      endif else begin 
                if(abs(results[3,INT]) lt abs(results[3,1])) then begin
                check_value[*]=results[*,INT-1]
            endif else begin
                print,'check_value is wrong'
                return, 0
           endelse
      endelse
    endif else begin 
      if (check eq 0) and (iteration eq 9) then begin ;check if no change of sign,
                             ;and checks if it should be upper or lower limit
      if(abs(results[3,1]) lt abs(results[3,INT])) then begin
         check_value[*]=results[*,1]
      endif else begin 
                if (abs(results[3,INT]) lt abs(results[3,1])) then begin
                check_value[*]=results[*,INT]
            endif else begin
                print,'check_value is wrong'
                return, 0
            endelse
          endelse
      endif
    endelse
  endfor
  ; end of iterations
  ;****************************
  result1 = check_value[2]
  return, result1
end
