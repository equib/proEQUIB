orig_quiet = !quiet
!quiet = 1
@proequib__compile_all
!quiet = orig_quiet

profiler
profiler, /system

; You need to have IDLdoc (https://github.com/mgalloy/idldoc)
; create doc folder in docs/
idldoc, root='pro/', output='docs/doc/', $
  title='API Documentation for proEQUIB', $
  subtitle='IDL Library for Plasma Diagnostics and Abundance Analysis', index_level=1, /nosource, $
  overview='docs/overview', footer='docs/footer', /embed, $
  format_style='rst', markup_style='rst'

; create tex folder in docs/
idldoc, root='pro/', output='docs/tex/', $
  title='API Documentation for proEQUIB', $
  subtitle='IDL Library for Plasma Diagnostics and Abundance Analysis', index_level=1, /nosource, $
  overview='docs/overview', footer='docs/footer', /embed, $
  template_prefix='latex-', comment_style='latex', $
  format_style='rst', markup_style='rst'

exit
