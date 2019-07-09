function ion_unit::init
   self.base_dir = file_dirname(file_dirname((routine_info('ion_unit__define', /source)).path))
   return,1
end

;----------------------------------------------------------------

pro ion_unit::set_base_dir, base_dir
   if base_dir ne '' then self.base_dir=base_dir else print, 'Error: base_dir is not given'
   return
end

function ion_unit::get_base_dir
   if self.base_dir ne '' then base_dir=self.base_dir else print, 'Error: base_dir is not given'
   return, base_dir
end

;----------------------------------------------------------------

function ion_unit::get
   if self.atom ne '' and  self.ion ne '' then begin
      ion_name=[self.atom, self.ion]
   endif
   return,ion_name
end

;------------------------------------------------------------------

pro ion_unit__define
     void={ion_unit, atom:'',ion:'', base_dir:''}
   return
end
