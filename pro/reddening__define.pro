; docformat = 'rst'

;+
;     This obejct library can be used to determine the reddening 
;     law function of the line at the given wavelength
;     for the used extinction law.
;
; :Examples:
;    For example::
;
;     IDL> ext=obj_new('reddening')
;     IDL> wavelength=6563.0
;     IDL> m_ext=1.0
;     IDL> flux=1.0
;     IDL> R_V=3.1
;     IDL> 
;     IDL> fl=ext->redlaw_gal(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.32013816
;        
;     IDL> fl=ext->redlaw_gal2(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.30925984
;        
;     IDL> fl=ext->redlaw_ccm(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.29756615
;        
;     IDL> fl=ext->redlaw_jbk(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.33113684
;        
;     IDL> fmlaw='AVGLMC'
;     IDL> fl=ext->redlaw_fm(wavelength, fmlaw=fmlaw, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.35053032
;        
;     IDL> fl=ext->redlaw_smc(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.22659261
;        
;     IDL> fl=ext->redlaw_lmc(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.30871187
;        
;     IDL> fl=ext->redlaw(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.32013816
;        
;     IDL> ext_law='GAL'
;     IDL> R_V=3.1
;     IDL> flux_deredden=ext->deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
;     IDL> print, 'dereddened flux(6563)', flux_deredden
;        dereddened flux(6563)      0.47847785
;        
;     IDL> flux_deredden=ext->deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
;     IDL> print, 'dereddened flux(6563)', flux_deredden
;        dereddened flux(6563)       4.7847785
;        
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
;   0.2.0
;
; :History:
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;     
;
;     08/07/2019, A. Danehkar, Move to object-oriented programming (OOP).
;-
function reddening::init
   self.base_dir = file_dirname(file_dirname((routine_info('reddening__define', /source)).path))
   return,1
end

;----------------------------------------------------------------

function reddening::redlaw, wavelength, ext_law=ext_law, rv=rv, fmlaw=fmlaw
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
;   0.2.0
;
; :History:
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw(wavelength, ext_law=ext_law, rv=rv, fmlaw=fmlaw)
    return, value
end

;----------------------------------------------------------------

function reddening::redlaw_gal, wavelength, rv=rv
;+
;    This function determines the reddening law function of the line at the given wavelength
;    for Galactic Seaton1979+Howarth1983+CCM1983.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).
;
; :Params:
;    wavelength :  in, required, type=float
;                   Wavelength in Angstrom
;
; :Keywords:
;    rv       :  in, optional, type=float, default=3.1
;                the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> R_V=3.1
;     IDL> fl=redlaw_gal(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.32013816
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on the UV Formulae from Seaton 1979, MNRAS, 187, 73
;     1979MNRAS.187P..73S, the opt/NIR from Howarth 1983, MNRAS, 203, 301
;     the FIR from Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
;     1989ApJ...345..245C
;
;     Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x, pyneb.extinction
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw_gal(wavelength, rv=rv)
    return, value
end
  
function reddening::redlaw_gal2, wavelength
;+
;    This function determines the reddening law function of the line at the given wavelength
;    for Galactic Savage & Mathis 1979.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).
;
; :Params:
;     wavelength :  in, required, type=float
;                   Wavelength in Angstrom
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> fl=redlaw_gal2(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.30925984
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;     Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on Savage & Mathis 1979, ARA&A, vol. 17, 73-111
;
;     Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x
;
;     20/09/1994, R. A. Shaw, Initial IRAF implementation.
;
;     04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw_gal2(wavelength)
    return, value
end

function reddening::redlaw_ccm, wavelength, rv=rv
;+
;    This function determines the reddening law function of Cardelli, Clayton & Mathis.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value for the given wavelength.
;
; :Params:
;     wavelength :  in, required, type=float/array
;                   Wavelength in Angstrom
;
; :Keywords:
;    RV       :  in, optional, type=float, default=3.1
;                the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> R_V=3.1
;     IDL> fl=redlaw_ccm(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.29756615
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on Formulae by Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256.
;     1989ApJ...345..245C
;
;     Originally from IRAF STSDAS SYNPHOT redlaw.x
;
;     18/05/1993, R. A. Shaw, Initial IRAF implementation, based upon CCM module
;         in onedspec.deredden.
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw_ccm(wavelength, rv=rv)
    return, value
end

function reddening::redlaw_jbk, wavelength
;+
;    This function determines the reddening law function for Galactic Whitford1958 + Seaton1977 + Kaler1976.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).
;
; :Params:
;     wavelength :  in, required, type=float
;                   Wavelength in Angstrom
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> fl=redlaw_jbk(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.33113684
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on Whitford (1958), extended to the UV by Seaton (1977),
;     adapted by Kaler (1976).
;
;     Originally from IRAF STSDAS SYNPHOT redlaw.x
;
;     13/05/1993, R. A. Shaw, Initial IRAF implementation.
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw_jbk(wavelength)
    return, value
end

function reddening::redlaw_fm, wavelength, rv=rv, fmlaw=fmlaw
;+
;    This function determines the reddening law function by Fitzpatrick & Massa
;    for the line at the given wavelength.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value for the given wavelength.
;
; :Params:
;     wavelength :  in, required, type=float/array
;                   Wavelength in Angstrom
;
; :Keywords:
;    RV       :  in, optional, type=float, default=3.1
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
;     IDL> fl=redlaw_fm(wavelength, rv=R_V)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.35054942
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on Formulae by Fitzpatrick 1999, PASP, 11, 63
;     1999PASP..111...63F, Fitzpatrick & Massa 1990,
;     ApJS, 72, 163, 1990ApJS...72..163F
;
;     Adopted from NASA IDL Library & PyAstronomy.
;
;     30/12/2016, A. Danehkar, Revised in IDL code.
;-
    value=redlaw_fm(wavelength, rv=rv, fmlaw=fmlaw)
    return, value
end

function reddening::redlaw_smc, wavelength
;+
;    This function determines the reddening law function of the line at the given wavelength
;    for Small Magellanic Cloud.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).
;
; :Params:
;     wavelength :  in, required, type=float
;                   Wavelength in Angstrom
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> fl=redlaw_smc(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.22659261
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on Prevot et al. (1984), A&A, 132, 389-392
;     1984A%26A...132..389P
;
;     Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
;
;     20/09/1994, R. A. Shaw, Initial IRAF implementation.
;
;     04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw_smc(wavelength)
    return, value
end

function reddening::redlaw_lmc, wavelength
;+
;    This function determines the reddening law function of the line at the given wavelength
;    for the Large Magellanic Cloud.
;
; :Returns:
;    type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).
;
; :Params:
;     wavelength :  in, required, type=float
;                   Wavelength in Angstrom
;
; :Examples:
;    For example::
;
;     IDL> wavelength=6563.0
;     IDL> fl=redlaw_lmc(wavelength)
;     IDL> print, 'fl(6563)', fl
;        fl(6563)     -0.30871187
;
; :Categories:
;   Interstellar Extinction
;
; :Dirs:
;  ./
;      Subroutines
;
; :Author:
;   Ashkbiz Danehkar
;
; :Copyright:
;   This library is released under a GNU General Public License.
;
; :Version:
;   0.2.0
;
; :History:
;     Based on Formulae by Howarth 1983, MNRAS, 203, 301
;     1983MNRAS.203..301H
;
;     Originally from IRAF STSDAS SYNPHOT ebmvlfunc.x, redlaw.x
;
;     18/10/1994, R. A. Shaw, Initial IRAF implementation.
;
;     14/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.
;
;     31/08/2012, A. Danehkar, Converted to IDL code.
;-
    value=redlaw_lmc(wavelength)
    return, value
end

;------------------------------------------------------------------

function reddening::deredden_flux, wavelength, flux, m_ext, ext_law=ext_law, rv=rv, fmlaw=fmlaw
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
;   0.2.0
;
; :History:
;     31/08/2012, A. Danehkar, IDL code.
;-
    value=deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=rv, fmlaw=fmlaw)
    return, value
end

function reddening::deredden_relflux, wavelength, relflux, m_ext, ext_law=ext_law, rv=rv, fmlaw=fmlaw
;+
;     This function dereddens flux intensity relative to Hb=100,  based on the reddening law.
;
; :Returns:
;    type=double. This function returns the deredden flux intensity relative to Hb=100.
;
; :Params:
;     wavelength :  in, required, type=float/array
;                   Wavelength in Angstrom
;     relflux       :  in, required, type=float,
;                   flux intensity relative to Hb=100
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
;     IDL> flux_deredden=deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V) ; deredden absolute flux intensity
;     IDL> print, 'dereddened relative flux(6563):', flux_deredden
;        dereddened relative flux(6563):       0.47847785
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
;   0.2.0
;
; :History:
;     31/08/2012, A. Danehkar, IDL code.
;-
    value=deredden_relflux(wavelength, relflux, m_ext, ext_law=ext_law, rv=rv, fmlaw=fmlaw)
    return, value
end

;------------------------------------------------------------------

pro reddening__define
    void={reddening, base_dir:''}
    return
end
