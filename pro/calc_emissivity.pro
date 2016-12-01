function calc_emissivity, ion, levels, tempi, densi
;+
; NAME:
;     calc_emissivity
; PURPOSE:
;     calculate line emissivities for specified ion with level(s) by 
;     solving atomic level populations and in statistical equilibrium 
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
;     emiss5007=double(0.0) 
;     emiss5007=calc_emissivity(ion, levels5007, tempi, densi)
;     print, emiss5007
;
; INPUTS:
;     ion -       ion name e.g. 'oii', 'oiii'
;     levels -    level(s) e.g '1,2/', '1,2,1,3/'
;     tempi -     electron temperature
;     densi -     electron density
; RETURN:  ionic abundance
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
  common share1, Atomic_Data_Path

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s 
    
  GX= long(0)
  ID=lonarr(2+1)
  JD=lonarr(2+1)
  
  I= long(0)
  J= long(0) 
  K= long(0) 
  L= long(0)
  NLINES= long(0) 
  NLEV= long(0) 
  NTEMP= long(0) 
  IBIG= long(0) 
  IRATS= long(0) 
  NTRA= long(0) 
  ITEMP= long(0) 
  IN= long(0) 
  KP1= long(0)
  IT= long(0) 
  IP1= long(0)
  IKT= long(0) 
  IC= long(0) 
  IC1= long(0) 
  IC2= long(0) 
  
  DENS=double(0)
  TEMP=double(0)
  
  EJI=double(0)
  WAV=double(0)
  SUMC=double(0)
  
  QX=double(0)
  AX=double(0)
  EX=double(0)
  LTEXT = '';

  abund=double(0)
     
  I= long(0)
  J= long(0)
  K= long(0)
  IP1= long(0)
  
  ion1=strtrim(ion,1)
  atomic_filename = Atomic_Data_Path+'/'+ion1+'.dat'
  openr, lun1, atomic_filename, /get_lun
  readf,lun1,NLINES
  for i=1, NLINES do begin 
    readf,lun1,LTEXT
  endfor
  ; Read no. of levels = NLEV,
  readf,lun1,NLEV, NTEMP
  
  Glj=lonarr(NLEV+1)
  
  Nlj=dblarr(NLEV+1)
  Omij=dblarr(NTEMP+1,NLEV+1,NLEV+1)   
  Aij=dblarr(NLEV+1,NLEV+1)   
  Elj=dblarr(NLEV+1)   
  Telist=dblarr(NTEMP+1)
  
  LABEL1=STRARR(NLEV+1) 
  
  Glj[*]=0
  levels_str=strsplit(levels, ',', ESCAPE='/', /EXTRACT)
  temp=size(levels_str, /N_ELEMENTS)
  levels_num=long(temp[0]/2)
  ITRANC=lonarr(2+1,levels_num+1)
  ITRANC[*,*]=0
  levels_i=0
  for i=1, levels_num do begin 
    res=equib_strnumber(levels_str[levels_i], val)
    if res eq 1 then ITRANC[1,i]=long(val)
    res=equib_strnumber(levels_str[levels_i+1], val)
    if res eq 1 then ITRANC[2,i]=long(val)
    levels_i = levels_i + 2
    ;if levels_i ge 2*levels_num then break
  endfor
  
  ; no. of Te = NTEMP and the
  for i=1, NLEV do begin 
    ; input format (cf Readme)
    readf,lun1,LTEXT
    LABEL1[I]=LTEXT
  endfor
  ; be
  ibig=0 
  ; Read in Te's where coll. strengths are tabulated
  for i=1, NTEMP do begin 
    ddtemp=double(0.0)
    readf,lun1,ddtemp
    Telist[i] = ddtemp
    Telist[i] = alog10(Telist[i])
  endfor 
  ; If IRATS=0, what tabulated are collision strengths
  readf,lun1,IRATS
  ; Else Coll. rates = tabulated values * 10 ** IRATS
 
  if (IBIG eq 0) then begin
    QX = 1.0
    while (QX ne 0.D0) do begin 
      lontemp1=long(0)
      lontemp2=long(0)
      ddtemp=double(0)
      readf,lun1,lontemp1, lontemp2, ddtemp
      ID[2]=lontemp1
      JD[2]=lontemp2
      QX = ddtemp
      if QX eq 0 then break
      if (ID[2] eq 0) then begin
        ID[2] = ID[1]
        K = K + 1
      endif else begin
        ID[1]= ID[2]
        K = 1
      endelse
      if (JD[2] eq 0) then begin
        JD[2] = JD[1]
      endif else begin
        JD[1] = JD[2]
      endelse
      I = ID[2]
      J = JD[2]
      Omij[K,I,J] = QX
    endwhile
  endif
  
  if(IBIG eq 1) or (IBIG eq 2) then begin
    readf,lun1,NTRA
    for IN = 1, NTRA do begin
      readf,lun1,I,J,Omij[1:NTEMP,I,J]
      ;READ(1,*) I,J,(Omij(ITEMP,I,J),ITEMP=1,NTEMP)
    endfor
  endif
  ; Read transition probabilities
  if (IBIG eq 1) then begin
    readf,lun1,I,J,Aij[J,I];,L=K+1,NLEV),K=1,NLEV1
    ;READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
  endif else begin
    for K = 1, NLEV - 1 do begin
      KP1 = K + 1
      for L = KP1, NLEV do begin
        readf,lun1, I, J, AX
        Aij[J,I] = AX
      endfor
    endfor
  endelse
  ; Read statistical weights, energy levels (cm-1)
  for J = 1, NLEV do begin
    readf,lun1,  I, GX, EX
    Glj[I] = GX
    Elj[I] = EX
  endfor
      
  free_lun, lun1
  
  ; Get levels for ratio 
  ; 150 large enough
  ; 
  ; Read in Te and Ne where the line 
  ; ratio is to be calculated
  ; Set Te
  TEMP=TEMPI
  ; Set Ne
  DENS=DENSI
  if (TEMP le 0.D0) or (DENS le 0.D0) then begin
      print,'Temp = ', TEMP, ', Dens = ', DENS
      return, 0
  endif

  Nlj=calc_populations(TEMP, DENS, Telist, Omij, Aij, Elj, Glj, NLEV, NTEMP, IRATS)
  
  emissivity_all=double(0.0)
  ; Search ITRANC for transitions & sum up
  for IKT=1, levels_num do begin 
    I=ITRANC[1,IKT]
    J=ITRANC[2,IKT]
    emissivity_line=double(0.0)
    if (Aij[J,I] ne 0.D0) then begin
      EJI = Elj[J] - Elj[I]
      WAV = 1.D8 / EJI
      emissivity_line=Nlj[J]*Aij[J,I]*h_Planck*c_Speed*1.e8/(WAV*DENS)
      emissivity_all=emissivity_all+emissivity_line
    endif
  endfor
  return, emissivity_all
end
