=======
proEQUIB
=======
    
.. image:: https://travis-ci.org/equib/proEQUIB.svg?branch=master
    :target: https://travis-ci.org/equib/proEQUIB
    :alt: Build Status

.. image:: https://ci.appveyor.com/api/projects/status/ab7ad315c6xejw3c?svg=true
    :target: https://ci.appveyor.com/project/danehkar/proequib
    :alt: Build Status

.. image:: http://mybinder.org/badge.svg
    :target: http://mybinder.org/repo/equib/proequib
    :alt: Binder

.. image:: https://img.shields.io/aur/license/yaourt.svg
    :target: https://github.com/equib/proEQUIB/blob/master/LICENSE
    :alt: GitHub license

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1890337.svg
    :target: https://doi.org/10.5281/zenodo.1890337
    :alt: Zenodo

.. image:: http://joss.theoj.org/papers/10.21105/joss.00899/status.svg
    :target: https://doi.org/10.21105/joss.00899
    :alt: JOSS


Description
============

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

Documentation
=============

For more information on how to use the API functions from the proEQUIB libray, please read the `API Documentation  <https://equib.github.io/proEQUIB/doc>`_ published on `equib.github.io/proEQUIB <https://equib.github.io/proEQUIB>`_.


References
==========
* Danehkar, A. (2018). proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **3**, 899. doi:`10.21105/joss.00899 <https://doi.org/10.21105/joss.00899>`_

* Danehkar, A. (2018). Bi-Abundance Ionisation Structure of the Wolf-Rayet Planetary Nebula PB 8, *PASA*, **35**, e005.  doi:`10.1017/pasa.2018.1 <https://doi.org/10.1017/pasa.2018.1>`_ ads:`2018PASA...35....5D <http://adsabs.harvard.edu/abs/2018PASA...35....5D>`_.
