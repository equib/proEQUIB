; docformat = 'rst'

;+
; Returns proEQUIB version. This file is automatically edited 
; by the builder to include the revision.
;
; :Returns:
;    string
;
; :Keywords:
;    full : in, optional, type=boolean
;       set to return Subversion revision as well
;-
function proequib_version, full=full
  compile_opt strictarr, hidden

  version = '0.3.0'
  revision = '-02b23a61'

  return, version + (keyword_set(full) ? (' ' + revision) : '')
end
