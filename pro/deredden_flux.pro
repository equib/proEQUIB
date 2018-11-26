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
;   0.0.1
;
; :History:
;     31/08/2012, A. Danehkar, IDL code.
;-

;+
; NAME:
;     deredden_flux
; PURPOSE:
;     This function dereddens absolute flux intensity based on the reddening law.
;
; CALLING SEQUENCE:
;     flux_deredden = deredden_flux(Wavelength, Flux, M_ext, EXT_LAW=ext_law, RV=rv, FMLAW=fmlaw)
;
; INPUTS:
;     Wavelength  :  in, required, type=float, 
;                    wavelength in Angstroms
;     Flux        :  in, required, type=float,    
;                    absolute flux intensity
;     M_ext       :  in, required, type=float,    
;                    logarithmic extinction 
;
; KEYWORD PARAMETERS:
;     EXT_LAW :  in, optional, type=string, default='GAL'
;                extinction law:
;                ext_law='GAL' for Howarth Galactic
;                ext_law='GAL2' for Savage and Mathis
;                ext_law='CCM' for CCM galactic
;                ext_law='JBK' for Whitford, Seaton, Kaler
;                ext_law='FM' for Fitxpatrick
;                ext_law='SMC' for Prevot SMC
;                ext_law='LMC' for Howarth LMC
;    RV       :  in, optional, type=float, default=3.1,
;                the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).
;    FMLAW    :  in, optional, type=string
;                the fmlaw keyword is used only in the redlaw_fm function:
;                fmlaw=''   the default fit parameters for the R-dependent
;                           Galactic extinction curve from Fitzpatrick & Massa
;                           (Fitzpatrick, 1999, PASP, 111, 63).
;                fmlaw='LMC2'  the fit parameters are those determined for
;                              reddening the LMC2 field (inc. 30 Dor)
;                              from Misselt et al.  (1999, ApJ, 515, 128).
;                fmlaw='AVGLMC'   the fit parameters are those determined for
;                              reddening in the general Large Magellanic Cloud (LMC)
;                              field by Misselt et al.  (1999, ApJ, 515, 128).
; 
; OUTPUTS:  This function returns a double as the deredden flux intensity.
;
; PROCEDURE: This function calls redlaw.
;
; EXAMPLE:
;     wavelength=6563.0
;     ext_law='GAL'
;     R_V=3.1
;     m_ext=1.0
;     flux=1.0
;     flux_deredden=deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V) ; deredden absolute flux intensity
;     print, 'dereddened flux(6563):', flux_deredden
;     > dereddened flux(6563):       4.7847785
;     
; MODIFICATION HISTORY:
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
