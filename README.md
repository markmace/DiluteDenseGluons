# DILUTE-DENSE GLUONS #
# README #

### What is this repository for? ###

* A program for calculating gluon production and multi-particle correlations in ultra-relativistic nuclear collisions using the dilute-dense CGC framework
* v0.9

### How do I get set up? ###

* Compiles with C++ with MPI and OpenMP
* MPI is only used to run different seeds simulanteously on different cores. 
* Currently has no Open-MP optimization -- under construction with validation
* Dependencies: GSL, FFTW3 with multi-threading, LAPACK
* Has run on Mac OSX 10.8 and newer, NERSC

### How do I run? ###
* For the executable created, there are a large number of inputs which can be, and typically are, specified at run-time (see src/CommandlineParameters.cpp). If using the Glauber+IP-Sat module, there are additional input parameters in a separate input file (the name of which is specified also at run-time): an example is included with git repo. The lattice size and spacing are the only common parameters between the input file and the commandline inputs, both of which must be set.

### Who do I talk to? ###

* Author: Mark Mace -- https://sites.google.com/view/mark-mace 
* Source code: https://github.com/markfmace/DiluteDenseGluons
* Contact: mark.f.mace@jyu.fi OR markfmace@gmail.com
* Based on formalism from Kovchegov/Skokov: Phys.Rev. D97 (2018) no.9, 094021 (and McLerran/Skokov: Nucl.Phys. A959 (2017) 83-101)
* Glauber+IP-Sat portions of code modified from IP-Glasma code from B. Schenke, from Schenke, Tribedy, Venugopalan: Phys.Rev.Lett. 108 (2012) 252301, Phys.Rev. C86 (2012) 034908, Phys.Rev. C89 (2014) no.2, 024901 
