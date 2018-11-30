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

Description
============

The **proEQUIB** library is a collection of `Interactive Data Language <http://www.harrisgeospatial.com/ProductsandSolutions/GeospatialProducts/IDL.aspx>`_ (IDL)/`GNU Data Language <http://gnudatalanguage.sourceforge.net/>`_ (GDL) programs developed to calculate atomic level populations and line emissivities in statistical equilibrium in multi-level atoms for different physical conditions of the stratification layers where the chemical elements are ionized. This library includes the IDL/GDL implementation of the program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_, which was originally written in FORTRAN by `Howarth & Adams (1981) <http://adsabs.harvard.edu/abs/1981ucl..rept.....H>`_, and was recently converted to IDL/GDL. It also includes the IDL/GDL implementation of deredden functions from `STSDAS IRAF Package <http://www.stsci.edu/institute/software_hardware/stsdas>`_, and atomic data from the `AtomNeb database <https://github.com/atomneb/AtomNeb-idl>`_. It uses the IDL/GDL implementation of the `MIDAS <http://www.eso.org/~ohainaut/ccd/midas.html>`_ scripts by X. W. Liu and codes from `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_ for some heavy element recombination emissivities. 

History of codes can be found `here <https://physics.mq.edu.au/~ashkbiz/proequib/history/>`_.

Website: `physics.mq.edu.au/~ashkbiz/proequib <https://physics.mq.edu.au/~ashkbiz/proequib/>`_

Installation
============


To install the **proEQUIB** library, simply add the path of this package directory to your IDL path. 
This package requires IDL version 7.1 or later. For more information about the path management in IDL, read `the tips for customizing IDL program path <https://www.harrisgeospatial.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/16156/Quick-tips-for-customizing-your-IDL-program-search-path>`_
provided by Harris Geospatial Solutions or `the IDL library installation note <http://www.idlcoyote.com/code_tips/installcoyote.php>`_ by David Fanning in the Coyote IDL Library. 



References
==========
A Danehkar, `PASA, 35, e005, 2018 <http://adsabs.harvard.edu/abs/2018PASA...35....5D>`_

A Danehkar, Q A Parker & W Steffen, `AJ, 151, 38, 2016 <ttp://adsabs.harvard.edu/abs/2016AJ....151...38D>`_


