; docformat = 'rst'

function redlaw, wavelength, ext_law=ext_law, rv=rv, fmlaw=fmlaw
;+
;     This function determines the reddening law function of the line at the given wavelength  
;     for the used extinction law.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value for the given wavelength.
;           
; :Params:
;     wavelength :  in, required, type=float/array
;                   Wavelength in Angstrom
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
;     IDL> R_V=3.1
;     IDL> fl=redlaw(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.32013816
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
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-

;+
; NAME:
;     redlaw
; PURPOSE:
;     This function determines the reddening law function of the line at the given wavelength  
;     for the used extinction law.
;
; CALLING SEQUENCE:
;     Result = redlaw(Wavelength, EXT_LAW=ext_law, RV=rv, FMLAW=fmlaw)
;
; INPUTS:
;     Wavelength[] -  in, required, type=float/array, 
;               wavelength in Angstroms
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
; OUTPUTS:  This function returns a double/array as the reddening law function 
;                   value(s) f(lambda) for the given wavelength(s) lambda.
;
; PROCEDURE: This function calls redlaw_gal or redlaw_gal2 or redlaw_ccm or redlaw_jbk
;            or redlaw_fm or redlaw_smc or redlaw_lmc.
;            This function is called by deredden_flux and deredden_relflux.
;
; EXAMPLE:
;     wavelength=6563.0
;     R_V=3.1
;     fl=redlaw(wavelength, rv=R_V)
;     print, 'fl(6563)', fl
;     > fl(6563)     -0.32013816
;     
; MODIFICATION HISTORY:
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-  
  if keyword_set(ext_law) then begin
    extlaw=ext_law
  endif else begin
    extlaw='GAL'
  endelse
  case extlaw of
     'GAL': fl = redlaw_gal(wavelength, rv=rv)
     'GAL2': fl = redlaw_gal2(wavelength)
     'CCM': fl = redlaw_ccm(wavelength,rv=rv)
     'JBK': fl = redlaw_jbk(wavelegth)
     'FM': fl = redlaw_fm(wavelength,rv=rv,fmlaw=fmlaw)
     'SMC': fl = redlaw_smc(wavelength)
     'LMC': fl = redlaw_lmc(wavelength)
     else: print, 'ext_law cannnot find'
  endcase
  return, fl
end
