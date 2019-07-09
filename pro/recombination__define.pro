; docformat = 'rst'

;+
;     This obejct library can be used to determine 
;     the ionic abundance from the observed flux of
;     recombination lines (RL) by using the recombination 
;     coefficients.
;
; :Examples:
;    For example::
;
;     IDL> he1=obj_new('recombination')
;     IDL> he1->set,['he','ii'] ; He I
;     IDL> 
;     IDL> he2=obj_new('recombination')
;     IDL> he2->set,['he','iii'] ; He II
;     IDL> 
;     IDL> c2=obj_new('recombination')
;     IDL> c2->set,['c','iii'] ; C II
;     IDL> 
;     IDL> c3=obj_new('recombination')
;     IDL> c3->set,['c','iv'] ; C III
;     IDL> 
;     IDL> n2=obj_new('recombination')
;     IDL> n2->set,['n','iii'] ; N II
;     IDL> 
;     IDL> n3=obj_new('recombination')
;     IDL> n3->set,['n','iv'] ; N III
;     IDL> 
;     IDL> o2=obj_new('recombination')
;     IDL> o2->set,['o','iii'] ; O II
;     IDL> 
;     IDL> ne2=obj_new('recombination')
;     IDL> ne2->set,['ne','iii'] ; Ne II
;     IDL> 
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> 
;     IDL> ; 4120.84: linenum=7
;     IDL> ; 4387.93: linenum=8
;     IDL> ; 4437.55: linenum=9
;     IDL> ; 4471.50: linenum=10
;     IDL> ; 4921.93: linenum=12
;     IDL> ; 5015.68: linenum=13
;     IDL> ; 5047.74: linenum=14
;     IDL> ; 5875.66: linenum=15
;     IDL> ; 6678.16: linenum=16
;     IDL> ; 7065.25: linenum=17
;     IDL> ; 7281.35: linenum=18
;     IDL> he_i_4471_flux= 2.104
;     IDL> linenum=10; 4471.50
;     IDL> Abund_he_i=he1->calc_abund_he_i_rl(temperature=temperature, density=density, $
;     IDL>                               linenum=linenum, line_flux=he_i_4471_flux)
;     IDL> print, 'N(He^+)/N(H^+):', Abund_he_i
;        N(He^+)/N(H^+):     0.040848393
;        
;     IDL> he_ii_4686_flux = 135.833
;     IDL> Abund_he_ii=he2->calc_abund_he_ii_rl(temperature=temperature, density=density, $
;     IDL>                                 line_flux=he_ii_4686_flux)
;     IDL> print, 'N(He^2+)/N(H^+):', Abund_he_ii
;        N(He^2+)/N(H^+):      0.11228817
;        
;     IDL> c_ii_6151_flux = 0.028
;     IDL> wavelength=6151.43
;     IDL> Abund_c_ii=c2->calc_abund_c_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, line_flux=c_ii_6151_flux)
;     IDL> print, 'N(C^2+)/N(H+):', Abund_c_ii
;        N(C^2+)/N(H+):   0.00063404650
;        
;     IDL> c_iii_4647_flux = 0.107
;     IDL> wavelength=4647.42
;     IDL> Abund_c_iii=c3->calc_abund_c_iii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, line_flux=c_iii_4647_flux) 
;     IDL> print, 'N(C^3+)/N(H+):', Abund_c_iii
;        N(C^3+)/N(H+):   0.00017502840
;        
;     IDL> n_ii_4442_flux = 0.017
;     IDL> wavelength=4442.02
;     IDL> Abund_n_ii=n2->calc_abund_n_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, line_flux=n_ii_4442_flux)
;     IDL> print, 'N(N^2+)/N(H+):', Abund_n_ii
;        N(N^2+)/N(H+):   0.00069297541
;        
;     IDL> n_iii_4641_flux = 0.245
;     IDL> wavelength=4640.64
;     IDL> Abund_n_iii=n3->calc_abund_n_iii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, line_flux=n_iii_4641_flux)
;     IDL> print, 'N(N^3+)/N(H+):', Abund_n_iii
;        N(N^3+)/N(H+):   6.3366175e-05
;        
;     IDL> o_ii_4614_flux = 0.009
;     IDL> wavelength=4613.68
;     IDL> Abund_o_ii=o2->calc_abund_o_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, line_flux=o_ii_4614_flux)                      
;     IDL> print, 'N(O^2+)/N(H+):', Abund_o_ii
;        N(O^2+)/N(H+):    0.0018886330
;        
;     IDL> ne_ii_3777_flux = 0.056
;     IDL> wavelength=3777.14
;     IDL> Abund_ne_ii=ne2->calc_abund_ne_ii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, line_flux=ne_ii_3777_flux)
;     IDL> print, 'N(Ne^2+)/N(H+):', Abund_ne_ii
;        N(Ne^2+)/N(H+):   0.00043376850
;        
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
;   0.2.0
;
; :History:
;     calc_abund_n_ii_rl(), calc_abund_o_ii_rl(), calc_abund_ne_ii_rl()
;                   and calc_abund_c_ii_rl()  are mostly based on 
;                   scripts by Yong Zhang added to MOCASSIN, 02/2003
;                   Ercolano et al. 2005MNRAS.362.1038E.
;                   and MIDAS script written by X.W.Liu.
;
;     10/05/2013, A. Danehkar, Translated to IDL code.
;
;     25/04/2017, A. Danehkar, Integration with AtomNeb.
;     
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
;-
function recombination::init
   ;self.data_dir = 'externals/atomneb/atomic-data/chianti70'
   self.data_rc_dir = 'externals/atomneb/atomic-data-rc'
   self.base_dir = file_dirname(file_dirname((routine_info('recombination__define', /source)).path))
   self.Atom_RC_All_file= filepath('rc_collection.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   self.Atom_RC_He_I_file= filepath('rc_he_ii_PFSD12.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   self.Atom_RC_PPB91_file= filepath('rc_PPB91.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   self.Atom_RC_SH95_file= filepath('rc_SH95.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   self.Atom_RC_N_II_FSL13_file= filepath('rc_n_iii_FSL13.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   self.Atom_RC_O_II_SSB17_file= filepath('rc_o_iii_SSB17_orl_case_b.fits', root_dir=self.base_dir, subdir=self.data_rc_dir )
   return,1
end

function recombination::free
  ptr_free, self.data_elj
  data_elj=ptr_new(/ALLOCATE_HEAP)
  ptr_free, self.data_omij
  data_omij=ptr_new(/ALLOCATE_HEAP)
  ptr_free, self.data_aij
  data_aij=ptr_new(/ALLOCATE_HEAP)
  ptr_free, self.hi_rc_data
  hi_rc_data=ptr_new(/ALLOCATE_HEAP)
  return,1
end

pro recombination::set, atom_ion, new=new, wavelength_list=wavelength_list
    if n_elements(atom_ion) lt 2 then begin
      print, 'Error: atom and ionic level are not given'
      return
    endif
    if atom_ion[0] ne '' then self.atom=atom_ion[0] else print, 'Error: atom is not given'
    if atom_ion[1] ne '' then self.ion=atom_ion[1] else print, 'Error: ionic level is not given'
    if n_elements(atom_ion) gt 2 then if atom_ion[2] ne 0 then self.level = atom_ion[2]
    
    hi_rc_data=atomneb_read_aeff_sh95(self.Atom_RC_SH95_file, 'h', 'ii') ; H I
    ;h_i_aeff_data=hi_rc_data[0].Aeff
    self.hi_rc_data=ptr_new(hi_rc_data[0].Aeff)
    atom=self.atom
    ion=self.ion
    atom_ion=atom+'_'+ion
    if keyword_set(new) then begin
      atom_ion='_'+'new'
    endif
    case atom_ion of
     'he_ii':begin ; He I
                hei_rc_data=atomneb_read_aeff_he_i_pfsd12(self.Atom_RC_He_I_file, atom, ion)
                ;he_i_aeff_data=hei_rc_data[0].Aeff
                self.rc_data=ptr_new(hei_rc_data[0].Aeff)
             end
     'he_iii':begin ; He II
                heii_rc_data=atomneb_read_aeff_sh95(self.Atom_RC_SH95_file, atom, ion)
                ;he_ii_aeff_data=heii_rc_data[0].Aeff
                self.rc_data=ptr_new(heii_rc_data[0].Aeff)
             end
     'c_iii': begin ; C II
                cii_rc_data=atomneb_read_aeff_collection(self.Atom_RC_All_file, atom, ion)
                self.rc_data=ptr_new(cii_rc_data)
             end
     'c_iv': begin ; C III
                ciii_rc_data=atomneb_read_aeff_ppb91(self.Atom_RC_PPB91_file, atom, ion)
                self.rc_data=ptr_new(ciii_rc_data)
             end
     'n_iii': begin ; N II
                nii_rc_data_old=atomneb_read_aeff_collection(self.Atom_RC_All_file, atom, ion)
                nii_rc_data_old_br=atomneb_read_aeff_collection(self.Atom_RC_All_file, atom, ion, /br)
                self.rc_data=ptr_new(nii_rc_data_old)
                self.rc_data_br=ptr_new(nii_rc_data_old_br)
             end
     'n_iv': begin ; N III
                niii_rc_data=atomneb_read_aeff_ppb91(self.Atom_RC_PPB91_file, atom, ion)
                self.rc_data=ptr_new(niii_rc_data)
             end       
     'o_iii': begin ; O II
                oii_rc_data_old=atomneb_read_aeff_collection(self.Atom_RC_All_file, atom, ion)
                oii_rc_data_old_br=atomneb_read_aeff_collection(self.Atom_RC_All_file, atom, ion, /br)
                self.rc_data=ptr_new(oii_rc_data_old)
                self.rc_data_br=ptr_new(oii_rc_data_old_br)
             end 
     'ne_iii': begin ; Ne II
                neii_rc_data=atomneb_read_aeff_collection(self.Atom_RC_All_file, atom, ion)
                self.rc_data=ptr_new(neii_rc_data)
             end 
                                     
     'n_iii_new': begin  ; Ne II
                if keyword_set(wavelength_list) then begin
                  temp=size(wavelength_list,/DIMENSIONS)
                  n_line=temp[0]
                  if n_line ne -1 then begin
                    rc_element_template={Wavelength: float(0.0), Aeff:fltarr(7,4)}
                    nii_rc_data=replicate(rc_element_template, n_line)
                    for i=0,n_line-1 do  begin
                      wavelength=wavelength_list[i]
                      select_nii_aeff_data=atomneb_search_aeff_n_ii_fsl13(self.Atom_RC_N_II_FSL13_file, atom, ion, wavelength)
                      nii_rc_data[i].Wavelength=wavelength
                      nii_rc_data[i].Aeff=select_nii_aeff_data.Aeff
                    endfor
                   endif 
                   self.rc_data=ptr_new(nii_rc_data)
                endif else begin
                  print, 'Error: wavelength_list is not given for n_iii_new'
                endelse
             end  
     'o_iii_new': begin  ; O II
                if keyword_set(wavelength_list) then begin
                  temp=size(wavelength_list,/DIMENSIONS)
                  n_line=temp[0]
                  case1='B' ; O II
                  if n_line ne -1 then begin
                    rc_element_template={Wavelength: float(0.0), Aeff:fltarr(16,25)}
                    oii_rc_data=replicate(rc_element_template, n_line)
                    for i=0,n_line-1 do  begin
                      wavelength=orl_line[i].wavelength
                      select_oii_aeff_data=atomneb_search_aeff_o_ii_ssb17(Atom_RC_O_II_SSB17_file, atom, ion, case1, wavelength)
                      oii_rc_data[i].Wavelength=wavelength
                      oii_rc_data[i].Aeff=select_oii_aeff_data.Aeff
                    endfor
                  endif
                  self.rc_data=ptr_new(oii_rc_data)
                endif else begin
                  print, 'Error: wavelength_list is not given for o_iii_new'
                endelse
             end
    endcase
    ;print, *(self.rc_data)
    ;print, *(self.rc_data_br)
    ;print, *(self.hi_rc_data)
    return
end

function recombination::calc_abund_he_i_rl, temperature=temperature, density=density, $
                                            linenum=linenum, line_flux=line_flux
;+
;     This function determines the ionic abundance from the observed 
;     flux intensity for the given wavelength of He I recombination line 
;     by using the recombination coefficients from Porter et al. 
;     2012MNRAS.425L..28P.
;
; :Returns:
;    type=double. This function returns the ionic abundanc.
;
; :Keywords:
;     temperature    :    in, required, type=float
;                         electron temperature
;     density        :    in, required, type=float
;                         electron density
;     linenum        :    in, required, type=int
;                         Line Number for Wavelength
;                         
;                         Wavelength=4120.84:linenum=7,  
;                         
;                         Wavelength=4387.93: linenum=8, 
;                         
;                         Wavelength=4437.55: linenum=9, 
;                         
;                         Wavelength=4471.50: linenum=10, 
;                         
;                         Wavelength=4921.93: linenum=12, 
;                         
;                         Wavelength=5015.68: linenum=13, 
;                         
;                         Wavelength=5047.74: linenum=14, 
;                         
;                         Wavelength=5875.66: linenum=15, 
;                         
;                         Wavelength=6678.16: linenum=16, 
;                         
;                         Wavelength=7065.25: linenum=17, 
;                         
;                         Wavelength=7281.35: linenum=18. 
;                         
;     line_flux      :    in, required, type=float
;                         line flux intensity
;     he_i_aeff_data :    in, required, type=array/object
;                         He I recombination coefficients
;     h_i_aeff_data  :    in, required, type=array/object
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
;     IDL> ion='ii' ; He I
;     IDL> he_i_rc_data=atomneb_read_aeff_he_i_pfsd12(Atom_RC_He_I_file, atom, ion)
;     IDL> he_i_aeff_data=he_i_rc_data[0].Aeff
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> he_i_4471_flux= 2.104
;     IDL> linenum=10; 4471.50
;     IDL> Abund_he_i=calc_abund_he_i_rl(temperature=temperature, density=density, $
;                                       linenum=linenum, line_flux=he_i_4471_flux, $
;                                       he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(He^+)/N(H^+):', Abund_he_i
;        N(He^+)/N(H^+):     0.040848393
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
;   0.2.0
;
; :History:
;     Based on improved He I emissivities in the case B
;     from Porter et al. 2012MNRAS.425L..28P
;     
;     15/12/2013, A. Danehkar, IDL code written.
;     
;     20/03/2017, A. Danehkar, Integration with AtomNeb.
;-

    value=calc_abund_he_i_rl(temperature=temperature, density=density, $
                            linenum=linenum, line_flux=line_flux, $
                            he_i_aeff_data=*(self.rc_data), h_i_aeff_data=*(self.hi_rc_data))
    return, value                        
end

function recombination::calc_abund_he_ii_rl, temperature=temperature, density=density, $
                                             line_flux=line_flux
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
;   0.2.0
;
; :History:
;     Based on He II emissivities
;     from Storey & Hummer, 1995MNRAS.272...41S.
;
;     15/12/2013, A. Danehkar, IDL code written.
;
;     02/04/2017, A. Danehkar, Integration with AtomNeb.
;-
    value=calc_abund_he_ii_rl(temperature=temperature, density=density, $
                              line_flux=line_flux, $
                              he_ii_aeff_data=*(self.rc_data), h_i_aeff_data=*(self.hi_rc_data))
    return, value
end

function recombination::calc_abund_c_ii_rl, temperature=temperature, density=density, $
                                            wavelength=wavelength, line_flux=line_flux
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
;   0.2.0
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
    value=calc_abund_c_ii_rl(temperature=temperature, density=density, $
                             wavelength=wavelength, line_flux=line_flux, $
                             c_ii_rc_data=*(self.rc_data), h_i_aeff_data=*(self.hi_rc_data)) 
    return, value
end  
        
function recombination::calc_abund_c_iii_rl, temperature=temperature, density=density, $
                                             wavelength=wavelength, line_flux=line_flux
;+
;     This function determines the ionic abundance from the observed
;     flux intensity for the given wavelength of C III recombination line
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
;     c_iii_rc_data :     in, required, type=array/object
;                         C III recombination coefficients
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
;     IDL> atom='c'
;     IDL> ion='iv' ; C III
;     IDL> c_iii_rc_data=atomneb_read_aeff_ppb91(Atom_RC_PPB91_file, atom, ion)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> c_iii_4647_flux = 0.107
;     IDL> wavelength=4647.42
;     IDL> Abund_c_iii=calc_abund_c_iii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, line_flux=c_iii_4647_flux, $
;     IDL>                                 c_iii_rc_data=c_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(C^3+)/N(H+):', Abund_c_iii
;        N(C^3+)/N(H+):    0.00017502840
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
;   0.2.0
;
; :History:
;     Based on effective radiative recombination coefficients for C III lines from
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.
;
;     18/05/2013, A. Danehkar, Translated to IDL code.
;
;     06/04/2017, A. Danehkar, Integration with AtomNeb.
;-
    value=calc_abund_c_iii_rl(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=line_flux, $
                              c_iii_rc_data=*(self.rc_data), h_i_aeff_data=*(self.hi_rc_data)) 
    return, value
end

function recombination::calc_abund_n_ii_rl, temperature=temperature, density=density, $
                                            wavelength=wavelength, line_flux=line_flux
;+
;     This function determines the ionic abundance from the observed
;     flux intensity for the given wavelength of N II recombination line
;     by using the recombination coefficients from
;     Escalante & Victor 1990ApJS...73..513E.
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
;     n_ii_rc_br    :     in, required, type=array/object
;                         N II branching ratios (Br)
;     n_ii_rc_data  :     in, required, type=array/object
;                         N II recombination coefficients
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
;     IDL> atom='n'
;     IDL> ion='iii' ; N II
;     IDL> n_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> n_ii_rc_data_br=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion, /br)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> n_ii_4442_flux = 0.017
;     IDL> wavelength=4442.02
;     IDL> Abund_n_ii=calc_abund_n_ii_rl(temperature=temperature, density=density, $
;     IDL>                               wavelength=wavelength, line_flux=n_ii_4442_flux, $
;     IDL>                               n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data, $
;     IDL>                               h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(N^2+)/N(H+):', Abund_n_ii
;        N(N^2+)/N(H+):   0.00069297541
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
;   0.2.0
;
; :History:
;     Based on Effective recombination coefficients for N II lines from
;     Escalante & Victor 1990ApJS...73..513E.
;
;     Adopted from MIDAS Rnii script written by X.W.Liu.
;
;     Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
;                       Ercolano et al. 2005MNRAS.362.1038E.
;
;     10/05/2013, A. Danehkar, Translated to IDL code.
;
;     25/04/2017, A. Danehkar, Integration with AtomNeb.
;-
    value=calc_abund_n_ii_rl(temperature=temperature, density=density, $
                             wavelength=wavelength, line_flux=line_flux, $
                             n_ii_rc_br=*(self.rc_data_br), n_ii_rc_data=*(self.rc_data), $
                             h_i_aeff_data=*(self.hi_rc_data))
    return, value
end

function recombination::calc_abund_n_iii_rl, temperature=temperature, density=density, $
                                             wavelength=wavelength, line_flux=line_flux
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
;   0.2.0
;
; :History:
;     Based on  effective radiative recombination coefficients for N III lines from
;     Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.
;
;     10/05/2013, A. Danehkar, IDL code written.
;
;     20/04/2017, A. Danehkar, Integration with AtomNeb.
;-
    value=calc_abund_n_iii_rl(temperature=temperature, density=density, $
                             wavelength=wavelength, line_flux=line_flux, $
                             n_iii_rc_data=*(self.rc_data), h_i_aeff_data=*(self.hi_rc_data))
    return, value
end

function recombination::calc_abund_o_ii_rl, temperature=temperature, density=density, $
                                            wavelength=wavelength, line_flux=line_flux
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
;   0.2.0
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
;-
    value=calc_abund_o_ii_rl(temperature=temperature, density=density, $
                             wavelength=wavelength, line_flux=line_flux, $
                             o_ii_rc_br=*(self.rc_data_br), o_ii_rc_data=*(self.rc_data), $
                             h_i_aeff_data=*(self.hi_rc_data))
    return, value
end

function recombination::calc_abund_ne_ii_rl, temperature=temperature, density=density, $
                                             wavelength=wavelength, line_flux=line_flux
;+
;     This function determines the ionic abundance from the observed
;     flux intensity for the given wavelength of Ne II recombination line
;     by using the recombination coefficients from
;     Kisielius et al. (1998) & Storey (unpublished).
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
;     ne_ii_rc_data  :    in, required, type=array/object
;                         Ne II recombination coefficients
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
;     IDL> atom='ne'
;     IDL> ion='iii' ; Ne II
;     IDL> ne_ii_rc_data=atomneb_read_aeff_collection(Atom_RC_All_file, atom, ion)
;     IDL> temperature=double(10000.0)
;     IDL> density=double(5000.0)
;     IDL> ne_ii_3777_flux = 0.056
;     IDL> wavelength=3777.14
;     IDL> Abund_ne_ii=calc_abund_ne_ii_rl(temperature=temperature, density=density, $
;     IDL>                                 wavelength=wavelength, line_flux=ne_ii_3777_flux, $
;     IDL>                                 ne_ii_rc_data=ne_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
;     IDL> print, 'N(Ne^2+)/N(H+):', Abund_ne_ii
;        N(Ne^2+)/N(H+):    0.00043376850
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
;   0.2.0
;
; :History:
;     Based on effective radiative recombination coefficients for Ne II lines
;     from Kisielius et al. 1998A&AS..133..257K & Storey (unpublished).
;
;     Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.
;
;     02/2003, Yong Zhang, scripts added to MOCASSIN.
;
;     14/05/2013, A. Danehkar, Translated to IDL code.
;
;     10/04/2017, A. Danehkar, Integration with AtomNeb.
;-
    value=calc_abund_ne_ii_rl(temperature=temperature, density=density, $
                              wavelength=wavelength, line_flux=line_flux, $
                              ne_ii_rc_data=*(self.rc_data), h_i_aeff_data=*(self.hi_rc_data))
    return, value
end
;-------------
;pro recombination::set_data_dir, data_dir
;  if data_dir ne '' then self.data_dir=data_dir else print, 'Error: data_dir is not given'
;  return
;end

;function recombination::get_data_dir
;  if self.data_dir ne '' then data_dir=self.data_dir else print, 'Error: data_dir is not given'
;  return, data_dir
;end
;-------------
pro recombination::set_data_rc_dir, data_rc_dir
    if data_dir ne '' then self.data_dir=data_dir else print, 'Error: data_rc_dir is not given'
    return
end

function recombination::get_data_rc_dir
    if self.data_rc_dir ne '' then data_rc_dir=self.data_rc_dir else print, 'Error: data_rc_dir is not given'
    return, data_rc_dir
end
;-------------
pro recombination::set_Atom_RC_All_file, Atom_RC_All_file
    if Atom_RC_All_file ne '' then self.Atom_RC_All_file=Atom_RC_All_file else print, 'Error: Atom_RC_All_file is not given'
    return
end

function recombination::Atom_RC_All_file
    if self.Atom_RC_All_file ne '' then Atom_RC_All_file=self.Atom_RC_All_file else print, 'Error: Atom_RC_All_file is not given'
    return, Atom_RC_All_file
end
;-------------
pro recombination::set_Atom_Atom_RC_He_I_file, Atom_RC_He_I_file
    if Atom_RC_He_I_file ne '' then self.Atom_RC_He_I_file=Atom_RC_He_I_file else print, 'Error: Atom_RC_He_I_file is not given'
    return
end

function recombination::get_Atom_RC_He_I_file
    if self.Atom_RC_He_I_file ne '' then Atom_RC_He_I_file=self.Atom_RC_He_I_file else print, 'Error: Atom_RC_He_I_file is not given'
    return, Atom_RC_He_I_file
end
;-------------
pro recombination::set_Atom_Aij_file, Atom_RC_PPB91_file
    if Atom_RC_PPB91_file ne '' then self.Atom_RC_PPB91_file=Atom_RC_PPB91_file else print, 'Error: Atom_RC_PPB91_file is not given'
    return
end

function recombination::get_Atom_RC_PPB91_file
    if self.Atom_RC_PPB91_file ne '' then Atom_RC_PPB91_file=self.Atom_RC_PPB91_file else print, 'Error: Atom_RC_PPB91_file is not given'
    return, Atom_RC_PPB91_file
end
;-------------
pro recombination::set_Atom_RC_SH95_file, Atom_RC_SH95_file
    if Atom_RC_SH95_file ne '' then self.Atom_RC_SH95_file=Atom_RC_SH95_file else print, 'Error: Atom_RC_SH95_file is not given'
    return
end

function recombination::get_Atom_RC_SH95_file
    if self.Atom_RC_SH95_file ne '' then Atom_RC_SH95_file=self.Atom_RC_SH95_file else print, 'Error: Atom_RC_SH95_file is not given'
    return, Atom_RC_SH95_file
end
;-------------
pro recombination::set_Atom_RC_N_II_FSL13_file, Atom_RC_N_II_FSL13_file
    if Atom_RC_N_II_FSL13_file ne '' then self.Atom_RC_N_II_FSL13_file=Atom_RC_N_II_FSL13_file else print, 'Error: Atom_RC_N_II_FSL13_file is not given'
    return
end

function recombination::get_Atom_RC_N_II_FSL13_file
    if self.Atom_RC_N_II_FSL13_file ne '' then Atom_RC_N_II_FSL13_file=self.Atom_RC_N_II_FSL13_file else print, 'Error: Atom_RC_N_II_FSL13_file is not given'
    return, Atom_RC_N_II_FSL13_file
end
;-------------
pro recombination::set_Atom_RC_O_II_SSB17_file, Atom_RC_O_II_SSB17_file
    if Atom_RC_O_II_SSB17_file ne '' then self.Atom_RC_O_II_SSB17_file=Atom_RC_O_II_SSB17_file else print, 'Error: Atom_RC_O_II_SSB17_file is not given'
    return
end

function recombination::get_Atom_RC_O_II_SSB17_file
    if self.Atom_RC_O_II_SSB17_file ne '' then Atom_RC_O_II_SSB17_file=self.Atom_RC_O_II_SSB17_file else print, 'Error: Atom_RC_O_II_SSB17_file is not given'
    return, Atom_RC_O_II_SSB17_file
end
;------------------------------------------------------------------

pro recombination__define
    void={recombination, level:0L, data_rc_dir:'', $
            Atom_RC_All_file:'', Atom_RC_He_I_file:'', Atom_RC_PPB91_file:'', Atom_RC_SH95_file:'', $
            Atom_RC_N_II_FSL13_file:'', Atom_RC_O_II_SSB17_file:'', $
            rc_data:ptr_new(/ALLOCATE_HEAP), rc_data_br:ptr_new(/ALLOCATE_HEAP), $
            hi_rc_data:ptr_new(/ALLOCATE_HEAP), $
            inherits ion_unit}
    return 
end
