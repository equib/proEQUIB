; docformat = 'rst'

function deredden_flux, wavelength, flux, m_ext, ext_law=ext_law, rv=rv, fmlaw=fmlaw
;+
;     This function dereddens absolute flux intensity based on the reddening law.
;
; :Returns:
;    type=double. This function returns the deredden flux intensity.
;
; :Params:
;     wavelength :  in, required, type=float/array
;                   Wavelength in Angstrom
;     flux       :  in, required, type=float,    
;                   absolute flux intensity
;     m_ext      :  in, required, type=float,    
;                   logarithmic extinction 
;
; :Keywords:
;    ext_law  :  in, optional, type=string, default='GAL'
;                the extinction law:
;
;                'GAL' for Howarth Galactic.
;
;                'GAL2' for Savage and Mathis.
;
;                'CCM' for CCM galactic.
;
;                'JBK' for Whitford, Seaton, Kaler.
;
;                'FM' for Fitxpatrick.
;
;                'SMC' for Prevot SMC.
;
;                'LMC' for Howarth LMC.
;
;    rv       :  in, optional, type=float, default=3.1
;                the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).
;
;    fmlaw    :  in, optional, type=string, default='GAL'
;                the fmlaw keyword is used only in the redlaw_fm function:
;
;                'GAL' for  the default fit parameters for the R-dependent
;                           Galactic extinction curve from Fitzpatrick & Massa
;                           (Fitzpatrick, 1999, PASP, 111, 63).
;
;                'LMC2' for the fit parameters are those determined for
;                              reddening the LMC2 field (inc. 30 Dor)
;                              from Misselt et al.  (1999, ApJ, 515, 128).
;
;                'AVGLMC' for  the fit parameters are those determined for
;                              reddening in the general Large Magellanic Cloud (LMC)
;                              field by Misselt et al.  (1999, ApJ, 515, 128).
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> ext_law='GAL'
;     IDL> R_V=3.1
;     IDL> m_ext=1.0
;     IDL> flux=1.0
;     IDL> flux_deredden=deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V) ; deredden absolute flux intensity
;     IDL> print, 'dereddened flux(6563):', flux_deredden
;        dereddened flux(6563):       4.7847785
;
; :Categories:
;   Interstellar Extinction
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
;     31/08/2012, A. Danehkar, IDL code.
;-
  if keyword_set(ext_law) then begin
    extlaw=ext_law
  endif else begin
    extlaw='GAL'
  endelse
  if keyword_set(rv) then begin
    R_V=rv
  endif else begin
    R_V=3.1
  endelse
  if keyword_set(fmlaw) then begin
    fm_law = fmlaw
  endif else begin
    fm_law = 'STANDARD'
  endelse
  fl = redlaw(wavelength,ext_law=extlaw,rv=R_V, fmlaw=fm_law)
  int_dered = flux * 10.0^(m_ext*(1+fl))
  return, int_dered
end
