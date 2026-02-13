README file for Becker & Wolff X-ray spectral module.

Version 1.0.0

This is the README file for the Fortran and C-codes of the Becker & Wolff (2007) 
model for accretion onto magnetic neutron stars in accreting X-ray pulsar binary 
systems.

There are 2 Fortran codes in this release with their accumpanying .dat files in
the "localmodel" directory:

bwphys.f (bwphys.dat)

bwsim.f (bwsim.dat)

The accompanying .dat files give the default starting parameters for the initial 
model and the ranges for the parameters. The parameters can be separated into 
two different but equivalent sets of parameters that fully specify a spectral 
model. See the enclosed file BWparameters.pdf along with the publications 
Becker & Wolff 2007, ApJ, v.654, p.435, and Wolff et al. 2016, v.831, id.194, 
doi:10.3847/0004-637X/831/2/194, for explanations of the parameters.

There are also 2 C-codes in this release in the "source" directory along with
a Makefile that can be edited to compile the two .c codes on each particular 
user's system:

phys2sim.c

sim2phys.c

These two C-codes take the physical variable interface quantities to the 
similarity variable quantities (phys2sim) and from the similarity variable 
quantities to the physical variable formulation (sim2phys, see Becker & Wolff 
2007, ApJ, v.654, p.435 for details). 

Make sure that you are working on a system with the HEASoft tools installed 
from source! Otherwsie, you can not install as an XSPEC model the Becker & Wolff 
model codes included here.

Installation of the BW2007 codes is simple. The two Fortran codes bwphys.f 
and bwsim.f, along with bwphys.dat and bwsim.dat can be left in the local 
model directory ("localmodel-directory-path"). Edit the XSPEC script 
BWbuild.xcm with the actual "localmodel-directory-path" for your system. 
Then start XSPEC as normal and enter at the command line

XSPEC> @BWbuild.xcm

After this the XSPEC BW2007 module will be ready to use for X-ray spectral 
modeling. Remember that you must have installed your HEASoft suite from source 
in order for this to work. "bwphys.f" will allow you to use the physical 
variables interface to set the model parameters. "bwsim.f" will allow you to 
use the similarity variables interface to set the model parameters.

Do not put the codes phys2sim.c and sim2phys.c in this local XSPEC model 
directory. Leave these two C-codes in a local source directory so that they
can be compiled separately. A simple Makefile is included in this distribution
that works on a Macintosh with GNU make installed from the Macports project
web site. In order to compile these codes go into the "source" directory 
and type at the command line:

Unix%> make dirs

Unix%> make

If you use a different GNU compiler distribution such as the one from the
Homebrew web site the library and include paths in the Makefile will need
to be edited for your system. 

You can go back and forth between the two variable (similarity and physical) 
interfaces using the codes "phys2sim" and "sim2phys". These two C-codes have 
help messages so that you can run them from the command line. Also, in order 
to compile these two C codes the libraries stdio.h, stdlib.h, string.h, 
unistd.h, getopt.h, and math.h, are required. These should be available on 
any computer with GCC installed.

Note that when the accuracy switch is set to low accuracy and fast convergence, 
fits in XSPEC move much faster than with the accuracy switch set to high accuracy.

When the accuracy switch is set to 0 (ASwitch=0) then up to 11 terms in the series 
are computed for all three Comptonization integrals, and, the convergence criteria 
for both EWfunc and EMfunc are stringent. When the accuracy switch is set to 1 
(ASwitch=1) a more limited set of series terms are computed for all three 
Comptonization integrals and the convergence criteria for both EWfunc and EMfunc 
are significantly cruder.

Questions about this model can be addressed to the authors of this pacakge of codes:
Michael T. Wolff (mtwolff at cox dot net), and Peter A. Becker (pbecker at gmu dot edu).

If these model codes are used in a published paper then please cite the papers in the
Astrophysical Journal: Becker & Wolff 2007, ApJ, v.654, p.435, and Wolff et al. 2016,
v.831, id.194, doi:10.3847/0004-637X/831/2/194.

Good luck with your modeling of accreting X-ray pulsars!

February 11, 2026
