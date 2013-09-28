; Example: calc_temp_dens()
;     determine electron density or temperature from given 
;     flux intensity ratio for a fixed electron density 
;     or temperature using  calc_temp_dens function
;     from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 
path='proEQUIB/atomic-data/'
set_atomic_data_path, path

ion='s_ii'
levu='1,2,1,3/'
levl='1,5/'
diagtype='T'
dens = double(2550)
niiTratio=double(10.753)
temp=calc_temp_dens(ion, levu, levl, niiTratio, diagtype, dens) 
print, "Electron Temperature:",temp
      
ion='s_ii'
levu='1,2/'
levl='1,3/'
diagtype='D'
temp=double(7000.0);
siiNratio=double(1.506);;
dens=calc_temp_dens(ion, levu, levl, siiNratio, diagtype, temp) 
print, "Electron Density:",dens

end 
