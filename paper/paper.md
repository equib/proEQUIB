---
title: proEQUIB - IDL Library for Plasma Diagnostics and Abundance Analysis
tags:
  - IDL
  - plasma diagnostics
  - abundance analysis
  - ionized nebulae
  - astrophysics
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: "1, 2"
affiliations:
 - name: Department of Physics and Astronomy, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
date: 29 June 2018
bibliography: paper.bib
---

# Summary

The emission line fluxes measured from the spectra of ionized nebulae provide 
insights into their physical conditions and elemental abundances, as well as 
interstellar extinctions. We can determine temperatures, densities, and chemical 
abundances from the measured fluxes of _collisionally excited lines_ (CEL) and 
_optical recombination lines_ (ORL) of ionized nebulae 
[see e.g. @Danehkar:2013; @Danehkar:2014; @Danehkar:2016; @Danehkar:2018].
 
``proEQUIB`` is a library including several functions developed 
in the Interactive Data Language (IDL), which can be used to derive 
electron temperatures, electron densities, and ionic abundances from the measured fluxes. 
It employs the IDL library ``AtomNeb`` _Atomic Data for Ionized Nebulae_ [@Danehkar:2018b], 
which contains collision strengths, transition probabilities, and recombination 
coefficients for collisional excitation and recombination calculations. 
It includes several IDL functions to determine physical conditions and 
chemical abundances from CEL and ORL, and estimate interstellar extinctions 
from Balmer lines.

- The CEL functions, which were completely developed in the IDL programming 
language based on the algorithm of the FORTRAN program EQUIB [@Howarth:1981; @Howarth:2016], 
allow to determine electron temperatures, electron densities, and ionic abundances 
from the measured fluxes of _collisionally excited lines_.

- The ORL functions, which were developed in pure-IDL according to the algorithm 
of the recombination scripts by X. W. Liu and  Y. Zhang included in the FORTRAN 
program MOCASSIN [@Ercolano:2003; @Ercolano:2005], allow to derive ionic 
abundances from the measured fluxes of _optical recombination lines_.

- The IDL functions based on the methods of the reddening law functions 
from STSDAS IRAF Package [@Bushouse:1994; @Shaw:1994] can be used to determine 
_interstellar extinctions_ for different reddening laws from the measured 
fluxes of Balmer lines.

``proEQUIB`` has recently been used for plasma diagnostics and abundance analysis 
of planetary nebulae [@Danehkar:2016; @Danehkar:2018]. ``proEQUIB`` heavily 
relies on the IDL Astronomy User's library [@Landsman:1993; @Landsman:1995] 
and the IDL library AtomNeb [@Danehkar:2018b]. The IDL functions of this package 
can easily be used to generate spatially-resolved maps of temperatures, density, 
and ionic abundances from integral field spectroscopic observations 
[see e.g. @Danehkar:2013; @Danehkar:2014; @Danehkar:2014b].

# Acknowledgements

A.D. acknowledges the receipt of a Macquarie University Research Excellence Scholarship.

# References
