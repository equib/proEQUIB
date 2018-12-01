;+
; Builds the proequib sav file.
;-

; clear any other compilations
.reset

; compile required code

@proequib_compile_all

; create the sav file
save, filename='proequib.sav', /routines, description='proEQUIB ' + proequib_version(/full)

exit
