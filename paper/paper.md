---
title: 'proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis'
tags:
  - astrophysics
  - gaseous nebulae
  - plasma diagnostics
  - abundance analysis
  - IDL
  - GDL
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: "1, 2"
affiliations:
 - name: Research Centre in Astronomy, Astrophysics and Astrophotonics, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
date: 30 June 2018
bibliography: paper.bib
---

# Summary

The emission lines emitted from gaseous nebulae carry valuable
information about the physical conditions and chemical abundances of ionized gases 
in these objects, as well as the interstellar reddening. We determine the electron temperature, 
the electron density, and the ionic abundances from the dereddened fluxes of _collisionally excited lines_ (CEL) 
and _recombination lines_ (RL) identified in nebular spectra 
[see e.g. @Danehkar:2013; @Danehkar:2014; @Danehkar:2016; @Danehkar:2018].
 
``proEQUIB`` is a library including several application programming interface (API) functions developed 
in the Interactive Data Language (IDL), 
which can be used to determine temperatures, 
densities, and chemical abundances from emission lines of ionized nebulae. 
This IDL library can also be used with the GNU Data Language (GDL) [@Arabas:2010; @Coulais:2010], 
which is a free and open-source alternative IDL compiler. 
This IDL/GDL package employs the IDL library ``AtomNeb`` _Atomic Data for Ionized Nebulae_ [@Danehkar:2018b], 
which contains collision strengths and transition probabilities for collisional excitation calculations, 
and recombination coefficients for recombination calculations. 
This package includes several API functions to determine physical conditions and 
chemical abundances from CEL and RL, derive interstellar extinctions 
from Balmer lines, and deredden the observed fluxes:

- The API functions for _CEL analysis_ were developed in the IDL programming 
language based on the algorithm of the FORTRAN program EQUIB [@Howarth:1981; @Howarth:2016]. 
These API functions can be used to determine the electron temperature, the electron density, 
and the ionic abundances from the dereddened fluxes of _collisionally excited lines_ emitted from 
ionized gaseous nebulae.

- The API functions for _RL analysis_ were developed in IDL according to the algorithm 
of the recombination scripts by X. W. Liu and Y. Zhang included in the FORTRAN 
program MOCASSIN [@Ercolano:2003; @Ercolano:2005]. These API functions can be used to 
determine the ionic abundances from the dereddened fluxes of _recombination lines_ emitted from 
ionized nebulae.

- The API functions for _reddening analysis_ were developed based on the methods 
of the reddening functions in STSDAS IRAF Package [@Bushouse:1994; @Shaw:1994]. 
These API functions can be employed to obtain _interstellar extinctions_ for 
different reddening laws from the observed fluxes of Balmer lines detected 
in nebular spectra, and deredden the measured fluxes of emission lines. 

``proEQUIB`` has recently been used for plasma diagnostics and abundance analysis 
of some planetary nebulae [@Danehkar:2016; @Danehkar:2018]. This IDL/GDL package heavily 
relies on the IDL Astronomy User's library [@Landsman:1993; @Landsman:1995] 
and the IDL library AtomNeb [@Danehkar:2018b]. The API functions of this IDL library 
can easily be utilized to generate spatially-resolved maps of extinction, 
temperature, density, and chemical abundances from integral field spectroscopic observations 
[see e.g. @Danehkar:2013; @Danehkar:2014; @Danehkar:2014b].

# Acknowledgements

A.D. acknowledges the receipt of a Macquarie University Research Excellence Scholarship.

# References
