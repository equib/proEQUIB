========
proEQUIB
========
    
.. image:: https://travis-ci.org/equib/proEQUIB.svg?branch=master
    :target: https://travis-ci.org/equib/proEQUIB
    :alt: Build Status

.. image:: https://ci.appveyor.com/api/projects/status/ab7ad315c6xejw3c?svg=true
    :target: https://ci.appveyor.com/project/danehkar/proequib
    :alt: Build Status

.. image:: http://mybinder.org/badge.svg
    :target: http://mybinder.org/repo/equib/proequib
    :alt: Binder

.. image:: https://img.shields.io/badge/license-GPL-blue.svg
    :target: https://github.com/equib/proEQUIB/blob/master/LICENSE
    :alt: GitHub license

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1890337.svg
    :target: https://doi.org/10.5281/zenodo.1890337
    :alt: Zenodo

.. image:: http://joss.theoj.org/papers/10.21105/joss.00899/status.svg
    :target: https://doi.org/10.21105/joss.00899
    :alt: JOSS


Description
===========

The **proEQUIB** library is a collection of `Interactive Data Language <http://www.harrisgeospatial.com/ProductsandSolutions/GeospatialProducts/IDL.aspx>`_ (IDL)/`GNU Data Language <http://gnudatalanguage.sourceforge.net/>`_ (GDL) programs developed to perform plasma diagnostics and abundance analysis using emission line fluxes measured in ionzed nebulae. It uses the `AtomNeb IDL library <https://github.com/atomneb/AtomNeb-idl>`_ to read collision strengths and transition probabilities for collisionally excited lines (CEL), and recombination coefficients for recombination lines (RL). This IDL package can be used to determine interstellar extinctions, electron temperatures, electron densities, and ionic abundances from the measured fluxes of emission lines. It mainly contains the follwing API functions written purely in IDL/GDL: 

* **API functions for collisionally excited lines (CEL)** have been developed based on the algorithm of the FORTRAN program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_ written in FORTRAN by `Howarth & Adams (1981) <http://adsabs.harvard.edu/abs/1981ucl..rept.....H>`_. The program EQUIB calculates atomic level populations and line emissivities in statistical equilibrium in multi-level atoms for different physical conditions of the stratification layers where the chemical elements are ionized. Using the IDL/GDL implementation of the program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_, electron temperatures, electron densities, and ionic abundances are determined from the measured fluxes of collisionally excited lines.

* **API functions for recombination lines (RL)** have been developed based on the algorithm of the recombination scripts by X. W. Liu and Y. Zhang included in the FORTRAN program `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_. These API functiosn are used to determine ionic abundances from recombination lines for some heavy element ions.
 
* **API functions for reddening and extinctions** have been developed according to the methods of the reddening law functions from `STSDAS IRAF Package <http://www.stsci.edu/institute/software_hardware/stsdas>`_, which are used to obtain interstellar extinctions and deredden measured fluxes based on different reddening laws.

Installation
============

Dependent IDL Packages
----------------------

* This package requires the following packages:

    - `The IDL Astronomy User's Library <https://idlastro.gsfc.nasa.gov/homepage.html>`_
    
    - `The AtomNeb IDL Library <https://github.com/atomneb/AtomNeb-idl>`_
    
    - `IDL MCMC Hammer library <https://github.com/mcfit/idl_emcee>`_ (currently not used!)
    
* To get this package with all the dependent packages, you can simply use ``git`` command as follows::

        git clone --recursive https://github.com/equib/proEQUIB.git


Installation in IDL
-------------------

* To install the **proEQUIB** library in the Interactive Data Language (IDL), you need to add the path of this package directory to your IDL path. For more information about the path management in IDL, read `the tips for customizing IDL program path <https://www.harrisgeospatial.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/16156/Quick-tips-for-customizing-your-IDL-program-search-path>`_ provided by Harris Geospatial Solutions or `the IDL library installation note <http://www.idlcoyote.com/code_tips/installcoyote.php>`_ by David Fanning in the Coyote IDL Library. 

* This package requires IDL version 7.1 or later. 


Installation in GDL
-------------------

*  You can install the GNU Data Language (GDL) if you do not have it on your machine:

    - Linux (Fedora)::

        sudo dnf install gdl
    
    - Linux (Ubuntu)::
    
        sudo apt-get install gnudatalanguage
    
    - OS X::
    
        brew install gnudatalanguage
    
    - Windows: You can use the `GNU Data Language for Win32 <https://sourceforge.net/projects/gnudatalanguage-win32/>`_ (Unofficial Version) or you can compile the `GitHub source <https://github.com/gnudatalanguage/gdl>`_ using Visual Studio 2015 as shown in `appveyor.yml <https://github.com/gnudatalanguage/gdl/blob/master/appveyor.yml>`_.

* To install the **proEQUIB** library in GDL, you need to add the path of this package directory to your ``.gdl_startup`` file in your home directory::

    !PATH=!PATH + ':/home/proEQUIB/pro/'
    !PATH=!PATH + ':/home/proEQUIB/externals/misc/'
    !PATH=!PATH + ':/home/proEQUIB/externals/astron/pro/'
    !PATH=!PATH + ':/home/proEQUIB/externals/atomneb/pro/'

  You may also need to set ``GDL_STARTUP`` if you have not done in ``.bashrc`` (bash)::

    export GDL_STARTUP=~/.gdl_startup

  or in ``.tcshrc`` (cshrc)::

    setenv GDL_STARTUP ~/.gdl_startup

* This package requires GDL version 0.9.8 or later.

How to Use
==========

The Documentation of the IDL functions provides in detail in the *API Documentation* (`equib.github.io/proEQUIB/doc <https://equib.github.io/proEQUIB/doc>`_). There are three main object units:

* **Collision Unit** which have the API functions for plasma diagnostics and abundance analysis of collisionally excited lines. Here are some examples of using *Collision* Unit:

    - *Temperature*::

        s2=obj_new('collision')
        s2->set,['s','ii']
        upper_levels='1,2,1,3/'
        lower_levels='1,5/'
        density = double(2550)
        line_flux_ratio=double(10.753)
        temperature=s2->calc_temperature(line_flux_ratio=line_flux_ratio, density=density, $
                                         upper_levels=upper_levels, lower_levels=lower_levels)
        print, "Electron Temperature:", temperature

      which gives::
    
        Electron Temperature:       7920.2865

    - *Density*::

        s2=obj_new('collision')
        s2->set,['s','ii']
        upper_levels='1,2/'
        lower_levels='1,3/'
        temperature=double(7000.0);
        line_flux_ratio=double(1.506);
        density=s2->calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, $
                                 upper_levels=upper_levels, lower_levels=lower_levels)
        print, "Electron Density:", density

      which gives::
      
        Electron Density:       2312.6395

    - *Ionic Abundance*::

        o3=obj_new('collision')
        o3->set,['o','iii']
        levels5007='3,4/'
        temperature=double(10000.0)
        density=double(5000.0)
        iobs5007=double(1200.0)
        Abb5007=o3->calc_abundance(temperature=temperature, density=density, $
                              line_flux=iobs5007, atomic_levels=levels5007)
        print, 'N(O^2+)/N(H+):', Abb5007

      which gives::
      
        N(O^2+)/N(H+):   0.00041256231 
        
    - *Emissivity*::
    
        o3=obj_new('collision')
        o3->set,['o','iii']
        levels5007='3,4/'
        temperature=double(10000.0)
        density=double(5000.0)
        iobs5007=double(1200.0)
        emis=o3->calc_emissivity(temperature=temperature, density=density, $
                            atomic_levels=levels5007)
        print, 'Emissivity(O III 5007):', emis

      which gives::
      
        Emissivity(O III 5007):   3.6041012e-21
        

    - *Atomic Level Population*::

        s2=obj_new('collision')
        s2->set,['s','ii']
        density = double(1000)
        temperature=double(10000.0);
        Nlj=s2->calc_populations(temperature=temperature, density=density)
        print, 'Populations:', Nlj

      which prints::
      
        Populations: 0.96992832 0.0070036315 0.023062261 2.6593671e-06 3.1277019e-06

    - *Critical Density*::
    
        s2=obj_new('collision')
        s2->set,['s','ii']
        temperature=double(10000.0)
        N_crit=s2->calc_crit_density(temperature=temperature)
        print, 'Critical Densities:', N_crit

      which gives::
      
        Critical Densities: 0.0000000 5007.8396 1732.8414 1072685.0 2220758.1

    - *All Ionic Level Information*::
    
        o3=obj_new('collision')
        o3->set,['o','iii']
        temperature=double(10000.0)
        density=double(5000.0)
        o3->print_ionic, temperature=temperature, density=density

      which gives::
      
        Temperature =   10000.0 K
        Density =    1000.0 cm-3
        
        Level    Populations   Critical Densities 
        Level 1:   3.063E-01   0.000E+00
        Level 2:   4.896E-01   4.908E+02
        Level 3:   2.041E-01   3.419E+03
        Level 4:   4.427E-05   6.853E+05
        Level 5:   2.985E-09   2.547E+07
          
         2.597E-05  
             88.34um 
             (2-->1) 
         2.859E-22  
        
         0.000E+00   9.632E-05  
             32.66um      51.81um 
             (3-->1)     (3-->2) 
         0.000E+00   7.536E-22  
        
         2.322E-06   6.791E-03   2.046E-02  
           4932.60A    4960.29A    5008.24A 
            (4-->1)     (4-->2)     (4-->3) 
         4.140E-25   1.204E-21   3.593E-21  
        
         0.000E+00   2.255E-01   6.998E-04   1.685E+00  
           2315.58A    2321.67A    2332.12A    4364.45A 
            (5-->1)     (5-->2)     (5-->3)     (5-->4) 
         0.000E+00   5.759E-24   1.779E-26   2.289E-23  
        
        H-beta emissivity: 1.237E-25 N(H+) Ne  [erg/s]


* **Recombination Unit** which have the API functions for plasma diagnostics and abundance analysis of recombination lines. Here are some examples of using *Recombination* Unit:

    - *He+ Ionic Abundance*::

        he1=obj_new('recombination')
        he1->set,['he','ii'] ; He I
        temperature=double(10000.0)
        density=double(5000.0)
        he_i_4471_flux= 2.104
        linenum=10; 4471.50
        Abund_he_i=he1->calc_abundance(temperature=temperature, density=density, $
                                      linenum=linenum, line_flux=he_i_4471_flux)
        print, 'N(He^+)/N(H^+):', Abund_he_i

      which gives::
      
        N(He^+)/N(H^+):     0.040848393

    - *He++ Ionic Abundance*::
    
        he2=obj_new('recombination')
        he2->set,['he','iii'] ; He II
        temperature=double(10000.0)
        density=double(5000.0)
        he_ii_4686_flux = 135.833
        Abund_he_ii=he2->calc_abundance(temperature=temperature, density=density, $
                                        line_flux=he_ii_4686_flux)
        print, 'N(He^2+)/N(H^+):', Abund_he_ii

      which gives::
      
        N(He^2+)/N(H^+):      0.11228817

    - *C++ Ionic Abundance*::
    
        c2=obj_new('recombination')
        c2->set,['c','iii'] ; C II
        temperature=double(10000.0)
        density=double(5000.0)
        wavelength=6151.43
        c_ii_6151_flux = 0.028
        Abund_c_ii=c2->calc_abundance(temperature=temperature, density=density, $
                                      wavelength=wavelength, line_flux=c_ii_6151_flux)
        print, 'N(C^2+)/N(H+):', Abund_c_ii

      which gives::
      
        N(C^2+)/N(H+):   0.00063404650 
      
    - *C3+ Ionic Abundance*::

        c3=obj_new('recombination')
        c3->set,['c','iv'] ; C III
        temperature=double(10000.0)
        density=double(5000.0)
        wavelength=4647.42
        c_iii_4647_flux = 0.107
        Abund_c_iii=c3->calc_abundance(temperature=temperature, density=density, $
                                        wavelength=wavelength, line_flux=c_iii_4647_flux) 
        print, 'N(C^3+)/N(H+):', Abund_c_iii

      which gives::
      
        N(C^3+)/N(H+):   0.00017502840

    - *N++ Ionic Abundance*::
    
        n2=obj_new('recombination')
        n2->set,['n','iii'] ; N II
        wavelength=4442.02
        n_ii_4442_flux = 0.017
        Abund_n_ii=n2->calc_abundance(temperature=temperature, density=density, $
                                      wavelength=wavelength, line_flux=n_ii_4442_flux)
        print, 'N(N^2+)/N(H+):', Abund_n_ii

      which gives::
      
        N(N^2+)/N(H+):   0.00069297541

    - *N3+ Ionic Abundance*::
    
        n3=obj_new('recombination')
        n3->set,['n','iv'] ; N III
        wavelength=4640.64
        n_iii_4641_flux = 0.245
        Abund_n_iii=n3->calc_abundance(temperature=temperature, density=density, $
                                        wavelength=wavelength, line_flux=n_iii_4641_flux)
        print, 'N(N^3+)/N(H+):', Abund_n_iii

      which gives::
      
        N(N^3+)/N(H+):   6.3366175e-05

    - *O++ Ionic Abundance*::

        o2=obj_new('recombination')
        o2->set,['o','iii'] ; O II
        wavelength=4613.68
        o_ii_4614_flux = 0.009
        Abund_o_ii=o2->calc_abundance(temperature=temperature, density=density, $
                                      wavelength=wavelength, line_flux=o_ii_4614_flux)                      
        print, 'N(O^2+)/N(H+):', Abund_o_ii
        
      which gives::
      
        N(O^2+)/N(H+):    0.0018886330

    - *Ne++ Ionic Abundance*::

        ne2=obj_new('recombination')
        ne2->set,['ne','iii'] ; Ne II
        wavelength=3777.14
        ne_ii_3777_flux = 0.056
        Abund_ne_ii=ne2->calc_abundance(temperature=temperature, density=density, $
                                        wavelength=wavelength, line_flux=ne_ii_3777_flux)
        print, 'N(Ne^2+)/N(H+):', Abund_ne_ii

      which gives::
      
        N(Ne^2+)/N(H+):   0.00043376850


    - *He I Emissivity*::

        he1=obj_new('recombination')
        he1->set,['he','ii'] ; He I
        temperature=double(10000.0)
        density=double(5000.0)
        linenum=10; 4471.50
        emiss_he_i=he1->calc_emissivity(temperature=temperature, density=density, $
                                        linenum=linenum)
        print, 'He I Emissivity:', emiss_he_i

      which gives::
      
        He I Emissivity:   6.3822830e-26

    - *He II Emissivity*::
    
        he2=obj_new('recombination')
        he2->set,['he','iii'] ; He II
        temperature=double(10000.0)
        density=double(5000.0)
        emiss_he_ii=he2->calc_emissivity(temperature=temperature, density=density)
        print, 'He II Emissivity:', emiss_he_ii

      which gives::
      
        He II Emissivity:   1.4989134e-24

    - *C II Emissivity*::
    
        c2=obj_new('recombination')
        c2->set,['c','iii'] ; C II
        temperature=double(10000.0)
        density=double(5000.0)
        wavelength=6151.43
        emiss_c_ii=c2->calc_emissivity(temperature=temperature, density=density, $
                                       wavelength=wavelength)
        print, 'C II Emissivity:', emiss_c_ii

      which gives::
      
        C II Emissivity:   5.4719511e-26
      
    - *C III Emissivity*::

        c3=obj_new('recombination')
        c3->set,['c','iv'] ; C III
        temperature=double(10000.0)
        density=double(5000.0)
        wavelength=4647.42
        emiss_c_iii=c3->calc_emissivity(temperature=temperature, density=density, $
                                        wavelength=wavelength)
        print, 'C III Emissivity:', emiss_c_iii

      which gives::
      
        C III Emissivity:   7.5749632e-25

    - *N II Emissivity*::
    
        n2=obj_new('recombination')
        n2->set,['n','iii'] ; N II
        wavelength=4442.02
        emiss_n_ii=n2->calc_emissivity(temperature=temperature, density=density, $
                                       wavelength=wavelength)
        print, 'N II Emissivity:', emiss_n_ii

      which gives::
      
        N II Emissivity:   3.0397397e-26

    - *N III Emissivity*::
    
        n3=obj_new('recombination')
        n3->set,['n','iv'] ; N III
        wavelength=4640.64
        emiss_n_iii=n3->calc_emissivity(temperature=temperature, density=density, $
                                        wavelength=wavelength)
        print, 'N III Emissivity:', emiss_n_iii

      which gives::
      
        N III Emissivity:   4.7908644e-24

    - *O II Emissivity*::

        o2=obj_new('recombination')
        o2->set,['o','iii'] ; O II
        wavelength=4613.68
        emiss_o_ii=o2->calc_emissivity(temperature=temperature, density=density, $
                                       wavelength=wavelength)
        print, 'O II Emissivity:', emiss_o_ii
        
      which gives::
      
        O II Emissivity:   5.9047319e-27

    - *Ne II Emissivity*::

        ne2=obj_new('recombination')
        ne2->set,['ne','iii'] ; Ne II
        wavelength=3777.14
        emiss_ne_ii=ne2->calc_emissivity(temperature=temperature, density=density, $
                                         wavelength=wavelength)
        print, 'Ne II Emissivity:', emiss_ne_ii

      which gives::
      
        Ne II Emissivity:   1.5996881e-25
        
* **Reddening Unit** which have the API functions for estimating logarithmic extinctions at H-beta and dereddening observed fluxes based on reddening laws and extinctions. Here are some examples of using *Reddening* Unit:

    - *Reddening Law Function*::

        ext=obj_new('reddening')
        wavelength=6563.0
        R_V=3.1
        fl=ext->redlaw(wavelength, rv=R_V, ext_law='GAL')
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.32013816

    - *Galactic Reddening Law Function based on Seaton (1979), Howarth (1983), & CCM (1983)*::

        ext=obj_new('reddening')
        wavelength=6563.0
        R_V=3.1
        fl=ext->redlaw_gal(wavelength, rv=R_V)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.32013816

    - *Galactic Reddening Law Function based on Savage & Mathis (1979)*::

        ext=obj_new('reddening')
        wavelength=6563.0
        fl=ext->redlaw_gal2(wavelength)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.30925984

    - *Reddening Law Function based on Cardelli, Clayton & Mathis (1989)*::
    
        ext=obj_new('reddening')
        wavelength=6563.0
        R_V=3.1
        fl=ext->redlaw_ccm(wavelength, rv=R_V)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.29756615

    - *Galactic Reddening Law Function based on Whitford (1958), Seaton (1977), & Kaler(1976)*::
    
        ext=obj_new('reddening')
        wavelength=6563.0
        fl=ext->redlaw_jbk(wavelength)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.33113684

    - *Reddening Law Function based on Fitzpatrick & Massa (1990), Fitzpatrick (1999), Misselt (1999)*::
    
        ext=obj_new('reddening')
        wavelength=6563.0
        R_V=3.1
        fmlaw='AVGLMC'
        fl=ext->redlaw_fm(wavelength, fmlaw=fmlaw, rv=R_V)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.35053032

    - *Reddening Law Function for the Small Magellanic Cloud*::
    
        ext=obj_new('reddening')
        wavelength=6563.0
        fl=ext->redlaw_smc(wavelength)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.22659261

    - *Reddening Law Function for the Large Magellanic Cloud*::
    
        ext=obj_new('reddening')
        wavelength=6563.0
        fl=ext->redlaw_lmc(wavelength)
        print, 'fl(6563):', fl

      which gives::
      
        fl(6563):     -0.30871187

    - *Dereddening Absolute Flux*::

        ext=obj_new('reddening')
        wavelength=6563.0
        m_ext=1.0
        flux=1.0
        ext_law='GAL'
        R_V=3.1
        flux_deredden=ext->deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
        print, 'dereddened flux(6563)', flux_deredden

      which gives::
      
        dereddened flux(6563)       4.7847785

    - *Dereddening Relative Flux*::

        ext=obj_new('reddening')
        wavelength=6563.0
        m_ext=1.0
        flux=1.0
        ext_law='GAL'
        R_V=3.1
        flux_deredden=ext->deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=R_V)
        print, 'dereddened flux(6563)', flux_deredden

      which gives::
      
        dereddened flux(6563)      0.47847785


Documentation
=============

For more information on how to use the API functions from the proEQUIB libray, please read the `API Documentation  <https://equib.github.io/proEQUIB/doc>`_ published on `equib.github.io/proEQUIB <https://equib.github.io/proEQUIB>`_.


References
==========
* Danehkar, A. (2020). pyEQUIB Python Package, an addendum to proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **5**, 2798. doi:`10.21105/joss.02798 <https://doi.org/10.21105/joss.02798>`_.

* Danehkar, A. (2018). proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **3**, 899. doi:`10.21105/joss.00899 <https://doi.org/10.21105/joss.00899>`_  ads:`2018JOSS....3..899D <https://ui.adsabs.harvard.edu/abs/2018JOSS....3..899D>`_.

* Danehkar, A. (2018). Bi-Abundance Ionisation Structure of the Wolf-Rayet Planetary Nebula PB 8, *PASA*, **35**, e005.  doi:`10.1017/pasa.2018.1 <https://doi.org/10.1017/pasa.2018.1>`_ ads:`2018PASA...35....5D <https://ui.adsabs.harvard.edu/abs/2018PASA...35....5D>`_.


Citation
========

Using **proEQUIB** in a scholarly publication? Please cite these papers:

.. code-block:: bibtex

   @article{Danehkar2020,
     author = {{Danehkar}, Ashkbiz},
     title = {pyEQUIB Python Package, an addendum to proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis},
     journal = {Journal of Open Source Software},
     volume = {5},
     number = {55},
     pages = {2798},
     year = {2020},
     doi = {10.21105/joss.02798}
   }

   @article{Danehkar2018,
     author = {{Danehkar}, Ashkbiz},
     title = {proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis},
     journal = {Journal of Open Source Software},
     volume = {3},
     number = {32},
     pages = {899},
     year = {2018},
     doi = {10.21105/joss.00899}
   }

Learn More
==========

==================  =============================================
**Documentation**   https://equib.github.io/proEQUIB/doc/
**Repository**      https://github.com/equib/proEQUIB
**Issues & Ideas**  https://github.com/equib/proEQUIB/issues
**DOI**             `10.21105/joss.00899 <https://doi.org/10.21105/joss.00899>`_
**Archive**         `10.5281/zenodo.3313731 <https://doi.org/10.5281/zenodo.3313731>`_
==================  =============================================

