; Example: getdiagnostic()
;     determine electron density or temperature from given 
;     flux intensity ratio for a fixed electron density 
;     or temperature using  getdiagnostic function
;     from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 
Atom_Elj_file='/AtomNeb/atomic-data/chianti70/AtomElj.fits'
Atom_Omij_file='/AtomNeb/atomic-data/chianti70/AtomOmij.fits'
Atom_Aij_file='/AtomNeb/atomic-data/chianti70/AtomAij.fits'

atom='s'
ion='ii'
s_ii_elj=atomneb_read_elj(Atom_Elj_file, atom, ion, level_num=5) ; read Energy Levels (Ej) 
s_ii_omij=atomneb_read_omij(Atom_Omij_file, atom, ion) ; read Collision Strengths (Omegaij)
s_ii_aij=atomneb_read_aij(Atom_Aij_file, atom, ion) ; read Transition Probabilities (Aij)\ 

levu='1,2,1,3/'		
levl='1,5/'
diagtype='T'
dens = double(2550)
siiTratio=double(10.753)
temp=do_diagnostic(s_ii_elj, s_ii_omij, s_ii_aij, levu, levl, siiTratio, diagtype, dens) 
print, "Electron Temperature:", temp

levu='1,2/'		
levl='1,3/'
diagtype='D'
temp=double(7000.0);
siiNratio=double(1.506);
dens=do_diagnostic(s_ii_elj, s_ii_omij, s_ii_aij, levu, levl, siiNratio, diagtype, temp) 
print, "Electron Density:", dens

end 
