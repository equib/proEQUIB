; Example: calc_abundance()
;     determine ionic abundance from observed 
;     flux intensity for gievn electron density 
;     and temperature using  calc_abundance function
;     from proEQUIB
; 
; --- Begin $MAIN$ program. ---------------
; 
; 
path='proEQUIB/atomic-data/'
set_atomic_data_path, path

ion='o_iii'
tempi=double(10000.0)
densi=double(5000.0)
levels5007='3,4/'
iobs5007=double(1200.0)
Abb5007=double(0.0) 
Abb5007=calc_abundance(ion, levels5007, tempi, densi, iobs5007)
print, Abb5007

end 
