/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Program: phys2sim.c
%%  Purpose: Take the input values of the physical variable parameters of the 
%%           XSPEC implementation of the Becker & Wolff model and recast them 
%%           as the Greek letter similarity input variables. Either implementation
%%           can fully specify the BW spectral model flow parameters. 
%%
%%  Author:  Michael Wolff
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*
**	Standard C-library headers:
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>

#define STRLEN 1024
#define DEG2RAD (2.0*3.14159265358979323846/360.0)
#define RAD2DEG (360.0/(2.0*3.14159265358979323846))
#define PI 3.14159265358979323846

extern char *optarg;
extern int optopt,optind;

/*
**  MAIN code: Drive everything else!
*/
int main(int argc, char **argv)
{
	double	acrat,comptem,caprad,dist;
	double	NSmass,NSRadius,sigperp,sigpara,sigbar,Lcrit,w,Lambda;
	double	xi,protonmass,electronmass,speedlight,thomcs,alpha,delta;
	double	speedlight2;
	double	kboltz,b12,RAlfven,Lacc,Grav,LEdd,EddAccRat;
	double	vff,Lcoul,taustar,Nmax,Ecyc,TempEdd,alttauthmnd;
	double 	dummy,timeesc,zsonic,batzsonic,perpopdepth;
	double	umax,zmax,epsabs,rhoabs,chiabs,tauzsonic,tauthmmnd;
	double	thmdtem,thmdrho,thmdu,C1,taumax,batzmax,zthmmnd;
	double	TeffWMush,TeffDMush,thmdne,thmdthomopdth,tauatzmax;
	double	erg2kev,erg2ev,kev2erg,ev2erg,StefBoltz;
	double	firstvacenabs,bcrit,RMagnetopheric,plasmabeta;
	char	*convertablenumber;
	double	NSMagMom;
	double	thmdplasmaenergy,thmdplasmafreq,electroncharge;
	double 	ThmdCouplingGamma,Zabun,ThmdTemp6,Aabun;
	int		c,printhelp;

/*
**	Initialize some physical constants I will need:
*/
	protonmass=1.672623e-24;
	electronmass=9.1093897e-28;
	electroncharge=4.80325e-10;
	kboltz=1.380658e-16;
	Grav=6.67259e-8;
	thomcs=0.66524574e-24;
	speedlight=2.99792458e10;
	speedlight2=speedlight*speedlight;
	erg2kev=6.2414e8;
	erg2ev=6.2414e11;
	kev2erg=1.0/erg2kev;
	ev2erg=1.0/erg2ev;
	dummy=erg2kev*erg2ev*kev2erg*ev2erg;
	dummy=dummy*dummy;    /* Get rid of annoying compiler warning. */
	StefBoltz=5.6695e-5;
/*
**	Set defaults for some of the variables:
*/
	acrat=0.0;
	NSmass=1.4;
	NSRadius=10.0;
	comptem=0.0;
	caprad=0.0;
	sigperp=1.0;
	sigpara=0.0;
	sigbar=0.0;
	b12=0.0;
	xi=-1.0;
	delta=-1.0;
	w=1.0;   /*  Default value. */
	ThmdCouplingGamma=1.0;
	Zabun=1.0;
	ThmdTemp6=1.0;
	Aabun=1.0;
/*  
**	Default value. Parameter from Ghosh & Lamb treatment.
**		Lambda=1.00   Accretion from a stellar wind.
**		Lambda=0.50   Accretion from a Disk.
**		Lambda=0.10   Value used in Becker et al. (2012).
*/
	Lambda=0.5000;
	printhelp=0;
/*
**	Parse command line:
*/
	while ( ( c=getopt(argc,argv,"hva:b:c:d:m::p::r::s:t:z:") ) != -1 ) {
		switch(c) {
			case 'a': /* REQUIRED Accretion rate, units 1.0e17 grams/second */
				acrat=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( acrat==0.00 ) printhelp=1;
				break;
			case 'b': /* REQUIRED Magnetic Field Strength, units 1.0e12 Gauss; default=1.0 */
				b12=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
/* 				if ( b12==0.00 ) printhelp=1; */
				break;
			case 'c': /* REQUIRED Accretion polar cap radius, units meters. */
				caprad=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( caprad==0.00 ) printhelp=1;
				break;
			case 'd': /* REQUIRED Source distance, units kiloparsecs  */
				dist=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( dist==0.00 ) printhelp=1;
				break;
			case 'm': /* OPTIONAL Neutron star mass, units solar masses, default=1.4 */
				NSmass=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( NSmass==0.00 ) printhelp=1;
				break;
			case 'p': /* OPTIONAL Perpendicular scattering cross section, units Thomson value */
				sigperp=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( sigperp==0.00 ) printhelp=1;
				break;
			case 'r': /* OPTIONAL Neutron star radius, units kilometers, default=10.0 */
				NSRadius=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( NSRadius==0.00 ) printhelp=1;
				break;
			case 's': /* REQUIRED Parallel scattering cross section, units Thomson value */
				sigpara=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( sigpara==0.00 ) printhelp=1;
				break;
			case 't': /* REQUIRED Compton Temperature, units keV. */
				comptem=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( comptem==0.00 ) printhelp=1;
				break;
			case 'v': /* OPTIONAL Print version information and exit. */
				fprintf(stderr," PHYS2SIM Version: 1.47 Date: 2019/10/23\n");
				exit(1);
				break;
			case 'z': /* REQUIRED BW2007 Mean cross section (bar), units Thomson value.   */
				sigbar=strtod(optarg,&convertablenumber);
				if ( optarg==convertablenumber ) printhelp=1;
 				if ( sigbar==0.00 ) printhelp=1;
				break;
			case 'h': /* OPTIONAL Print help message and exit. */
				printhelp=1;
				break;
			case '?': /* Unrecognized option. Print help message and exit. */
				fprintf(stderr,"Unrecognized command line option: -%c\n",optopt);
				printhelp=1;
				break;
			default:  /* Print help message and exit is default switch action. */
				printhelp=1;
				break;
		}
	}
/*
** 	Certain variables must be defined. 
*/
	if ( acrat<=0.0 ) printhelp=1;
	if ( comptem<=0.0 ) printhelp=1;
	if ( caprad<=0.0 ) printhelp=1;
	if ( sigpara<=0.0 ) printhelp=1;
	if ( sigbar<=0.0 ) printhelp=1;
	if ( b12<0.0 ) printhelp=1;
/*
**	Print help message and die.
*/
	if ( printhelp==1 ) {
		fprintf(stderr," Usage: phys2sim [hva:b:c:d:m::p::r::s:t:z:] \n");
		fprintf(stderr,"\t-a\tAccretion Rate (/10^17).\n");
		fprintf(stderr,"\t-b\tMagnetic Field Strength (/10^12 Gauss; default=1.0)\n");
		fprintf(stderr,"\t-c\tAccretion polar cap radius (m).\n");
		fprintf(stderr,"\t-d\tSource distance (kpc).\n");
		fprintf(stderr,"\t-m\tNeutron star mass (solar masses; optional; default=1.4).\n");
		fprintf(stderr,"\t-p\tPerpendicular cross section (Thomson; optional; default=1.0).\n");
		fprintf(stderr,"\t-r\tNeutron star radius (km; optional; default=10.0).\n");
		fprintf(stderr,"\t-s\tParallel cross section (Thomson).\n");
		fprintf(stderr,"\t-t\tCompton Temperature (in keV).\n");
		fprintf(stderr,"\t-v\tPrint version information.\n");
		fprintf(stderr,"\t-z\tMean (bar) cross section (Thomson).\n");
		fprintf(stderr,"\t-h\tDisplay help message.\n");
		exit(1);
	}
/*
**	Determine what parameters have only their default values and print out this information:
*/
	fprintf(stdout,"\n");
	fprintf(stdout," PHYS2SIM Version: 1.48 Date: 2025/05/23\n");
	fprintf(stdout,"\n");
	if(NSmass==1.4) fprintf(stdout," Neutron Star Mass has Default value 1.4 solar masses\n");
	if(NSRadius==10.0) fprintf(stdout," Neutron Star Radius has Default value 1.0E6 cm.\n");
	if(sigperp==1.0) fprintf(stdout," Sigma Perpendicular has Default value 1.0.\n");
	if(b12==1.0) fprintf(stdout," Magnetic Field Strength has Default value 1.0E12 Gauss.\n");
/*
**	Re-scale inputs:
*/
	acrat=acrat*1.0e17;
	caprad=caprad*1.0e2;
	NSmass=NSmass*1.989e33;
	NSRadius=NSRadius*1.0e5;
	sigperp=sigperp*thomcs;
	sigpara=sigpara*thomcs;
	sigbar=sigbar*thomcs;
	comptem=comptem*11604832.2;
	NSMagMom=0.5*b12*1.0e12*NSRadius*NSRadius*NSRadius;
/*
**	Go from system of physical variables to system of similarity (Greek) variables:
*/
	xi=PI*caprad*protonmass*speedlight/(acrat*sqrt(sigperp*sigpara));
	alpha=xi*(32.0*sqrt(3.0)/(49.0*log(7.0/3.0)))*(Grav*NSmass/(NSRadius*speedlight*speedlight));
	delta=(alpha/3.0)*(sigpara/sigbar)*(electronmass*speedlight*speedlight/(kboltz*comptem));

	if (b12!=0.0) {
		Ecyc=11.58*b12;
	} else {
		Ecyc=0.0;
	}
	vff=sqrt(2.0*Grav*NSmass/NSRadius);
	Nmax=electronmass*vff*vff/(2.0*Ecyc*kev2erg);
	taustar=51.4*pow(NSmass/(1.4*1.989e33),2.0)*pow(NSRadius/1.0e6,-2.0)/log(2.0*Nmax);
	Lacc=Grav*NSmass*acrat/NSRadius;
	LEdd=4.0*PI*Grav*NSmass*protonmass*speedlight/thomcs;
	TempEdd=pow(LEdd/(4.0*PI*NSRadius*NSRadius*StefBoltz),0.25);
	EddAccRat=LEdd*NSRadius/(Grav*NSmass);
/*
**	Calculate some quantities I will need the magnetic field strength for.
**	First, test to make sure the magnetic field is non-zero.
*/
	if (b12!=0.0) {
		Lcrit=1.49e37*
			pow(Lambda/0.1,-7.0/5.0)*
			pow(w,-28.0/15.0)*
			pow(NSmass/(1.4*1.989e33),29.0/30.0)*
			pow(NSRadius/1.0e6,1.0/10.0)*
			pow(b12,16.0/15.0);
		RAlfven=2.73e7*(Lambda/0.1)*pow(NSmass/(1.4*1.989e33),1.0/7.0)*
			pow(NSRadius/1.0e6,10.0/7.0)*pow(b12,4.0/7.0)*
			pow(Lacc/1.0e37,-2.0/7.0);
		Lcoul=1.17e37*pow(Lambda/0.1,-7.0/12.0)*pow(taustar/20.0,7.0/12.0)*
			pow(NSmass/(1.4*1.989e33),11.0/8.0)*pow(NSRadius/1.0e6,-13.0/24.0)*pow(b12,-1.0/3.0);
	} else {
		Lcrit=0.0;
		RAlfven=0.0;
		Lcoul=0.0;
	}

	timeesc=acrat*thomcs/(3.14159265*protonmass*vff*speedlight/7.0);
	zsonic=caprad*sqrt(sigperp/sigpara)*0.2445938;
	batzsonic=b12/pow((NSRadius+zsonic)/NSRadius,3.0);
	perpopdepth=timeesc*speedlight/caprad;
/*
**	Calcuate quantities at the thermal mound:
*/
	thmdtem=2.318e3*(pow(acrat,0.4))*(pow(caprad,-0.666666667));
	thmdrho=4.05e-12*pow(thmdtem,1.75)/sqrt(caprad);
	thmdne=thmdrho/protonmass;
	thmdthomopdth=caprad*thmdne*sigperp;
	thmdu=2.62*acrat/(pow(caprad,1.5)*pow(thmdtem,1.75));
	alttauthmnd=2.64e28*acrat*NSRadius/(NSmass*pow(caprad,1.5)*pow(thmdtem,7.0/4.0)*xi);
	zthmmnd=5.44e15*acrat*NSRadius/(NSmass*caprad*pow(thmdtem,7.0/2.0)*xi*sigpara);

	C1=4.0*Grav*NSmass*caprad*xi/(alpha*speedlight*speedlight*NSRadius*NSRadius)*sqrt(sigperp/sigpara);
	zmax=NSRadius*(sqrt(1.0+C1)-1.0)/2.0;
	taumax=pow(sigpara/sigperp,0.25)*sqrt(2.0*zmax/(alpha*xi*caprad));
	batzmax=b12/pow((NSRadius+zmax)/NSRadius,3.0);

	umax=sqrt(2.0*Grav*NSmass/(NSRadius+zmax));
	tauatzmax=(sigperp*acrat/(PI*caprad*protonmass*speedlight))*
	 			pow(sigperp/sigpara,0.25)*sqrt(xi*caprad/(2.0*alpha*zmax));

	epsabs=(6.08e12*sqrt(thmdrho)/(pow(comptem,1.75)))*sqrt(sqrt(thmdu*speedlight/umax))*kboltz*comptem;
	chiabs=epsabs/(kboltz*comptem);
	rhoabs=chiabs*chiabs*pow(comptem,3.5)/3.696640e25;
/*
**	Do some plasma physics here:
*/
	thmdplasmafreq=sqrt(4.0*PI*thmdne*electroncharge*electroncharge/electronmass);
	thmdplasmaenergy=0.371*sqrt(thmdne/1.0e26);
	bcrit=4.412e13;
	firstvacenabs=3.0*sqrt(rhoabs/0.0167)*(0.1*bcrit/(1.0e12*b12));
	plasmabeta=(rhoabs*kboltz*(comptem*11604832.2)/protonmass)/(b12*b12*1.0e24/(8.0*PI));
/*
**	Get the Coulomb coupling parameter for the plasma here:
*/
	ThmdTemp6=thmdtem/1.0e6;
	ThmdCouplingGamma=22.75*(Zabun*Zabun/ThmdTemp6)*pow(thmdrho/(1.0e6*Aabun),1.0/3.0) ;
/*
**	Calculate some parameters from Ghosh & Lamb theory:
*/
	RMagnetopheric=2.7e8*pow(Lacc/1.0e37,-2.0/7.0)*pow(b12,4.0/7.0);
/*
**	Compute the vertical variation of the perpendicular optical thickness as
**	a function of altitude above the neutron star surface.
*/
	tauzsonic=(sigperp*acrat/(PI*caprad*protonmass*speedlight))*
			pow(sigperp/sigpara,0.250)*
			sqrt(xi*caprad/(2.0*alpha*zsonic));
	tauthmmnd=(sigperp*acrat/(PI*caprad*protonmass*speedlight))*
			pow(sigperp/sigpara,0.250)*
			sqrt(xi*caprad/(2.0*alpha*zthmmnd));
/*
**	Begin to work in the Mushtukov et al. (2015) framework
*/
	if (b12!=0.0) {
		TeffWMush=4.5*pow(Lambda,1.0/4.0)*
					pow(b12,1.0/7.0)*
					pow(Lacc/1.0e37,5.0/28.0)*
					pow(NSmass/1.989e33,1.0/28.0)*
					pow(NSRadius/1.0e6,-11.0/28.0);
		TeffDMush=6.6*pow(Lambda,7.0/32.0)*
					pow(NSmass/1.989e33,13.0/80.0)*
					pow(NSRadius/1.0e6,-19.0/40.0)*
					pow(b12,1.0/8.0)*
					pow(Lacc/1.0e37,3.0/20.0);
	} else {
		TeffWMush=0.0;
		TeffDMush=0.0;
	}

/*
**	Print out results:
*/
	fprintf(stdout,"\n");
	fprintf(stdout," AccRate=   %13.7f CapRad=    %13.7f B12=     %13.7f \n",acrat/1.0e17,caprad/100.0,b12);
	fprintf(stdout," Alpha=     %13.7f Xi=        %13.7f Delta=   %13.7f \n",alpha,xi,delta);
	fprintf(stdout," CompTemp=  %13.7f keV (= %13.7e degK)\n",comptem/11604832.2,comptem);
	fprintf(stdout," NSMass=    %10.3f SM NSRadius=  %10.3f km Distance=%9.3f kpc\n",NSmass/1.989e33,NSRadius/1.0e5,dist);
	fprintf(stdout," SigPerp= %15.12f SigPara= %15.12f SigBar=%15.12f  \n",
					sigperp/thomcs,sigpara/thomcs,sigbar/thomcs);
	fprintf(stdout,"\n");
	fprintf(stdout," Luminosity Estimates Scaling Parameters: Lambda=%7.3f w=%7.3f\n",Lambda,w);
	fprintf(stdout," 1+z = %12.7f \n", 1.0/sqrt(1.0-2.0*Grav*NSmass/(NSRadius*speedlight*speedlight)));
	fprintf(stdout," z   = %12.7f \n", (1.0/sqrt(1.0-2.0*Grav*NSmass/(NSRadius*speedlight*speedlight)))-1.0);
	fprintf(stdout," 1+z @ Z(max) = %12.7f \n", 1.0/sqrt(1.0-2.0*Grav*NSmass/((NSRadius+zmax)*speedlight*speedlight)));
	fprintf(stdout," z @ Z(max)   = %12.7f \n", (1.0/sqrt(1.0-2.0*Grav*NSmass/((NSRadius+zmax)*speedlight*speedlight)))-1.0);
	fprintf(stdout," Neutron Star Compactness: R/2M=%9.4f\n",(NSRadius/1.0e5)/(2.0*NSmass/(1.989e33)));
	fprintf(stdout," Neutron Star Compactness: M/R=%9.4f\n",(NSmass/1.989e33)/(NSRadius/1.0e5));
	fprintf(stdout," Neutron Star Compactness: GM/Rc^2=%9.4f\n",(Grav*NSmass)/(speedlight2*NSRadius));
	fprintf(stdout,"\n");
	fprintf(stdout," Accretion Rate=%12.5e g/s (%12.5e SM/year)\n",acrat,acrat*(365.25*86400.0)/1.989e33);
	fprintf(stdout," Accretion Flux=%12.5e g/s/cm2\n",acrat/(PI*caprad*caprad));
	fprintf(stdout," Accretion Luminosity=%12.5e ergs/s\n",Lacc);
	fprintf(stdout," Eddington Luminosity=%12.5e ergs/s\n",LEdd);
	fprintf(stdout," Critical Luminosity=%12.5e ergs/s\n",Lcrit);
	fprintf(stdout," Eddington Accretion Rate=%12.5e g/s (%12.5e SM/year)\n",EddAccRat,EddAccRat*(365.25*86400.0)/1.989e33);
	fprintf(stdout," Eddington Accretion Flux=%12.5e g/s/cm2\n",EddAccRat/(4.0*PI*NSRadius*NSRadius));
	fprintf(stdout," Eddington Temperature=%10.5f keV\n",TempEdd/11604832.2);
	fprintf(stdout," Thermal Mound Temperature=%10.5f keV\n",thmdtem/11604832.2);
	fprintf(stdout," Coulomb Luminosity=%12.5e ergs/s\n",Lcoul);
	fprintf(stdout," Coulomb Stopping Optical Depth=%10.4f \n",taustar);
	fprintf(stdout,"\n");
	fprintf(stdout," Alfven Radius= %12.5e cm (%8.3f NS Radii)\n",RAlfven,RAlfven/NSRadius);
	fprintf(stdout," Magnetopheric Radius= %12.5e cm (%8.3f NS Radii)\n",RMagnetopheric,RMagnetopheric/NSRadius);
	fprintf(stdout,"\n");
	fprintf(stdout," (Alt) Thermal Mound Parallel Thomson Optical Depth=%10.4f \n",alttauthmnd);
	fprintf(stdout," (Alt) Thermal Mound Perpendicular Thomson Optical Depth=%10.4f \n",thmdthomopdth);
	fprintf(stdout," (Pete) Thermal Mound Perpendicular Thomson Optical Depth=%10.4f \n",tauthmmnd);
	fprintf(stdout," (Pete) Sonic Point Perpendicular Thomson Optical Depth=%10.4f \n",tauzsonic);
	fprintf(stdout," Perpendicular Thomson Optical Depth (BW Eq. 17)=%10.4f \n",perpopdepth);
	fprintf(stdout," Perpendicular Thomson Optical Depth @ Z(max) (Appendix Eq. A7)=%10.4f \n",tauatzmax);
	fprintf(stdout,"\n");
	fprintf(stdout," Z(max) (BW Eq. 80)=%12.5f km\n",zmax/1.0e5);
	fprintf(stdout," Z(sonic) (BW Eq. 31)=%12.5f km\n",zsonic/1.0e5);
	fprintf(stdout," Tau(max) (BW Eq. 79)=%12.5f\n",taumax);
	fprintf(stdout," B12 at Z(max) =%12.5f\n",batzmax);
	fprintf(stdout," B12 at Z(sonic) =%12.5f\n",batzsonic);
	fprintf(stdout," Geometrical Bloom Fraction @Z(max) = %12.5f\n",pow((NSRadius+zmax)/NSRadius,3.0));
	fprintf(stdout," Geometrical Bloom Fraction @Z(sonic) = %12.5f\n",pow((NSRadius+zsonic)/NSRadius,3.0));
	fprintf(stdout," Escape Time through column walls=%12.5e s\n",timeesc);
	fprintf(stdout," First Vacuum Energy(abs)=%12.5e keV\n",firstvacenabs);
	fprintf(stdout," Plasma Beta=%14.7f\n",plasmabeta);
	fprintf(stdout," Thermal Mound Plasma Energy=%14.7f keV\n",thmdplasmaenergy);
	fprintf(stdout," Thermal Mound Plasma Frequency=%14.7e s^-1\n",thmdplasmafreq);
	fprintf(stdout,"\n");
	fprintf(stdout," Neutron Star Magnetic Moment (.5BR^3)=%12.5e \n",NSMagMom);
	fprintf(stdout," Polar Cap NS Area Fraction=%12.9f \n",PI*caprad*caprad/(4.0*PI*NSRadius*NSRadius));
	fprintf(stdout," Magnetic Field Pressure (B2/8PI)=%12.5e \n",b12*b12*1.0e24/(8.0*PI));
	fprintf(stdout," Bremsstrahlung Self-Absorption Density=%14.6f g/cm3 \n",rhoabs);
	fprintf(stdout," Thermal Mound (Radiation Energy Density / Thermal Energy Density) ~%11.5f \n",
					1.0/(0.028*rhoabs/protonmass/(thmdtem*thmdtem*thmdtem)));
	fprintf(stdout," Thermal Mound Plasma Coulomb Coupling Parameter ~%13.7f \n",ThmdCouplingGamma);
	fprintf(stdout," Free Fall Velocity at NS Surface=%12.5e cm/s\n",sqrt(2.0*Grav*NSmass/NSRadius));
	if ( RAlfven!=0.0 ) {
		fprintf(stdout," Orbit Velocity at Alfven Radius=%12.5e cm/s\n",sqrt(2.0*Grav*NSmass/RAlfven));
	}
	fprintf(stdout,"\n");
	fprintf(stdout," Mushtukov et al. Teff (Wind)=%12.5f keV \n",TeffWMush);
	fprintf(stdout," Mushtukov et al. Teff (Disk)=%12.5f keV \n",TeffDMush);
	fprintf(stdout,"\n");
/* 
**	All done: Git!
*/
  	exit(0);
}
