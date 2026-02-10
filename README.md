README file for Becker & Wolff X-ray spectral module.

Version 1.0.0

This is the README file for the Fortran and C-codes of the Becker & Wolff (2007) 
model for accretion onto magnetic neutron stars in accreting X-ray pulsar binary 
systems.

There are 2 Fortran codes in this release with their accumpanying .dat files:

bwphys.f (bwphys.dat)

bwsim.f (bwsim.dat)

There are also 2 C-codes in this release:

phys2sim.c

sim2phys.c

These two C-codes take the physical variable interface quantities to the 
similarity variable quantities (phys2sim) and from the similarity variable 
quantities to the physical variable formulation (sim2phys, see Becker & Wolff 
2007, ApJ, v.654, p.435 for details). The accompanying .dat files give the 
default starting parameters for the initial model and the ranges for the 
parameters.

Make sure that you are working on a system with the HEASoft tools installed 
from source! Otherwsie, you can not install the Becker & Wolff model codes 
included here.

Installation of the BW2007 codes is simple. Put the two Fortran codes bwphys.f 
and bwsim.f, along with bwphys.dat and bwsim.dat in a local model directory 
("local-model-directory-path"), start XSPEC as normal and enter at the 
command line

XSPEC> initpackage bwphys bwphys.dat /local-model-directory-path/

XSPEC> lmod bwphys /local-model-directory-path/

XSPEC> initpackage bwsim bwsim.dat /local-model-directory-path/

XSPEC> lmod bwsim /local-model-directory-path/

The file BWbuild.xcm is included in this model release and is useable as long 
as the user edits the "local-model-directory-path" to be what is true for their 
particular local model directory. Do not put the codes phys2sim.c and sim2phys.c 
in this local XSPEC model directory. Put these two C-codes in a local source 
directory that can be compiled separately!

After this the XSPEC BW2007 module will be ready to use for X-ray spectral 
modeling. Remember that you must have installed your HEASoft suite from source 
in order for this to work. "bwphys.f" will allow you to use the physical 
variables interface to set the model parameters. "bwsim.f" will allow you to 
use the similarity variables interface to set the model parameters.

You can go back and forth between the two variable interfaces using the codes 
"phys2sim" and "sim2phys". Because these later codes are written in C they must 
be installed in a separate directory and compiled with your local system C-compiler. 
These two C-codes do have help messages so that you can run them from the command 
line. Also, in order to compile these two C codes the libraries stdio.h, stdlib.h, 
string.h, unistd.h, getopt.h, and math.h, are required. These should be 
available on any computer with GCC installed.

Note that when the accuracy switch is set to low accuracy and fast convergence, 
fits in XSPEC move much faster than with the accuracy switch set to high accuracy.

When the accuracy switch is set to 0 (ASwitch=0) then up to 11 terms in the series 
are computed for all three Comptonization integrals, and, the convergence criteria 
for both EWfunc and EMfunc are stringent. When the accuracy switch is set to 1 
(ASwitch=1) a more limited set of series terms are computed for all three 
Comptonization integrals and the convergence criteria for both EWfunc and EMfunc 
are significantly cruder.

Good luck with your modeling of accreting X-ray pulsars!

February 5, 2026
