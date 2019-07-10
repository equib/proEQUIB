; docformat = 'rst'

pro print_ionic, temperature=temperature, density=density, $
                 elj_data=elj_data, omij_data=omij_data, $
                 aij_data=aij_data, h_i_aeff_data=h_i_aeff_data, $
                 printEmissivity=printEmissivity, $
                 printPopulations=printPopulations, $
                 printCritDensity=printCritDensity
;+
;    This function prints the atom's transitions information,
;    atomic level populations, critical densities, and emissivities 
;    for given temperature and density.
;
; :Keywords:
;     temperature   :   in, required, type=float
;                       electron temperature
;     density       :   in, required, type=float
;                       electron density
;     elj_data      :   in, required, type=array/object
;                       energy levels (Ej) data
;     omij_data     :   in, required, type=array/object
;                       collision strengths (omega_ij) data
;     aij_data      :   in, required, type=array/object
;                       transition probabilities (Aij) data
;     h_i_aeff_data :   in, type=array/object 
;                       H I recombination coefficients 
;     printEmissivity  :   in, type=boolean
;                          Set for printing Emissivities
;     printPopulations :   in, type=boolean
;                          Set for printing Populations
;     printCritDensity  :  in, type=boolean
;                          Set for printing Critical Densities
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
;     IDL> temperature=double(10000.0);
;     IDL> density = double(1000.)
;     IDL> print_ionic, temperature=temperature, density=density, $, $
;     IDL>              elj_data=o_iii_elj, omij_data=o_iii_omij, $
;     IDL>              aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data[0].Aeff
;        Temperature =   10000.0 K
;        Density =    1000.0 cm-3
;        
;        Level    Populations   Critical Densities 
;        Level 1:   3.063E-01   0.000E+00
;        Level 2:   4.896E-01   4.908E+02
;        Level 3:   2.041E-01   3.419E+03
;        Level 4:   4.427E-05   6.853E+05
;        Level 5:   2.985E-09   2.547E+07
;        
;         2.597E-05  
;             88.34um 
;            (2-->1) 
;         2.859E-22  
;        
;         0.000E+00   9.632E-05  
;             32.66um      51.81um 
;            (3-->1)     (3-->2) 
;         0.000E+00   7.536E-22  
;        
;         2.322E-06   6.791E-03   2.046E-02  
;           4932.60A    4960.29A    5008.24A 
;            (4-->1)     (4-->2)     (4-->3) 
;         4.140E-25   1.204E-21   3.593E-21  
;        
;         0.000E+00   2.255E-01   6.998E-04   1.685E+00  
;           2315.58A    2321.67A    2332.12A    4364.45A 
;            (5-->1)     (5-->2)     (5-->3)     (5-->4) 
;         0.000E+00   5.759E-24   1.779E-26   2.289E-23  
;        
;        H-beta emissivity: 1.237E-25 N(H+) Ne  [erg/s]
;        
; :Categories:
;   Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
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
;     04/03/2019, A. Danehkar, create the print_ionic() routine.
;
;-

  h_Planck = 6.62606957e-27 ; erg s
  c_Speed = 2.99792458e10 ; cm/s

  if keyword_set(temperature) eq 1 then begin
    print,temperature, FORMAT='("Temperature = ", F9.1," K")
  endif
  if keyword_set(density) eq 1 then begin
    print,density, FORMAT='("Density = ", F9.1, " cm-3")'
  endif
  if keyword_set(elj_data) eq 0 then begin
    print,'elj_data is not set'
    return
  endif
  if keyword_set(omij_data) eq 0 then begin
    print,'omij_data is not set'
    return
  endif
  if keyword_set(aij_data) eq 0 then begin
    print,'aij_data is not set'
    return
  endif
  if keyword_set(printEmissivity) eq 0 then begin
    printEmissivity=1 ; default
  endif
  if keyword_set(printPopulations) eq 0 then begin
    printPopulations=1 ; default
  endif
  if keyword_set(printCritDensity) eq 0 then begin
    printCritDensity=1 ; default
  endif
  if keyword_set(level_num) eq 0 then begin
    temp=size(elj_data,/DIMENSIONS)
    level_num=temp[0]
  endif
  if keyword_set(printPopulations) eq 1 then begin
    if keyword_set(temperature) eq 1 and keyword_set(density) eq 1 then begin
      Nlj=calc_populations(temperature=temperature, density=density, $
                           elj_data=elj_data, omij_data=omij_data, $
                           aij_data=aij_data, level_num=level_num)
    endif else begin
      print,'calc_populations needs temperature and density.'
    endelse
  endif
  if keyword_set(printCritDensity) eq 1 then begin
    if keyword_set(temperature) eq 1 then begin
      N_crit=calc_crit_density(temperature=temperature, $
                              elj_data=elj_data, omij_data=omij_data, $
                              aij_data=aij_data, level_num=level_num)
    endif else begin
      print,'calc_crit_density needs temperature.'
    endelse
  endif
  print, ''
  if keyword_set(temperature) eq 1 and (keyword_set(printCritDensity) eq 1 or keyword_set(printPopulations) eq 1) then begin
    if keyword_set(density) eq 1 then begin
      Print, 'Level    Populations   Critical Densities 
    endif else begin
      Print, 'Level    Critical Densities 
    endelse
    for I = 1, level_num do begin
      s=""
      Level_str=STRING(I, FORMAT='("Level ",I1,":")')
      s=s+Level_str+"   "
      if keyword_set(density) eq 1 then begin
        Nlj_str=STRING(Nlj[I-1], FORMAT='(E9.3)')
        s=s+Nlj_str+"   "
      endif
      N_crit_str=STRING(N_crit[I-1], FORMAT='(E9.3)')
      s=s+N_crit_str
      print, s
    endfor
    print, ''
  endif
  if keyword_set(printEmissivity) eq 1 then begin
    Aij =aij_data.AIJ
    Elj =elj_data.Ej
    for I = 2, level_num do begin
      Aij_str=''
      transition_str=''
      wavelength_str=''
      Emissivity_str=''
      for J = 1, I-1 do begin
        EJI =  Elj[I-1] - Elj[J-1]
        WAV = 1.D8 / EJI
        Aij_str=Aij_str+STRING(Aij[I-1,J-1],FORMAT='(E10.3,"  ")')
        if WAV lt 10000 then begin
          wavelength_str=wavelength_str+STRING(WAV,FORMAT='(F10.2,"A ")')
        endif else begin
          wavelength_str=wavelength_str+STRING(WAV*1.e-4,FORMAT='(F10.2,"um ")')
        endelse
        transition_str=transition_str+STRING(I, J,FORMAT='("    (",I1,"-->",I1,") ")')
        if keyword_set(temperature) eq 1 and keyword_set(density) eq 1 then begin
          emissivity_line=Nlj[I-1]*Aij[I-1,J-1]*h_Planck*c_Speed*1.e8/(WAV*density)
          Emissivity_str=Emissivity_str+STRING(emissivity_line,FORMAT='(E10.3,"  ")')
        endif else begin
          Emissivity_str=Emissivity_str+' '
        endelse
      endfor
      print, Aij_str
      print, wavelength_str
      print, transition_str
      print, Emissivity_str
      print, ''
    endfor
    if keyword_set(h_i_aeff_data) eq 1 then begin
       if keyword_set(temperature) eq 1 and keyword_set(density) eq 1 then begin
         emissivity_Hbeta=calc_emiss_h_beta(temperature=temperature,density=density,h_i_aeff_data=h_i_aeff_data)
         Emissivity_str=STRING(emissivity_Hbeta,FORMAT='(E10.3)')
         print, "H-beta emissivity:", Emissivity_str, " N(H+) Ne  [erg/s]"
       endif else begin
         print,'gamma_hb_4861 needs temperature and density.'
       endelse
    endif
  endif
end
