CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Subroutine bwsim(ear,ne,param,ifl,photar,photer)
C     Purpose:   Similarity Variable interface for Becker & Wolff (2007)
C                spectral model code for XSPEC.
C
C     Inputs to BW2007 routine for XSPEC using Similarity Variable values:
C     (1) Real*8 ear(0:ne)     - Energy array boundaries. Each pair of 
C                                boundaries comprises an energy bin. I will 
C                                ultimately return the value of the photon 
C                                flux in each energy bin (photons cm^-2 s^-1) 
C                                to XSPEC.
C     (2) Integer ne           - Size of flux array. This is the number of 
C                                input points in the ear array. Note that 
C                                the number of input energies in the variable 
C                                ear is ne+1.
C     (3) Real*8 param(nparam) - Parameter values. Right now, there are ten 
C                                parameters that specify a model. Also, each 
C                                parameter value needs to be close to 1.0 and 
C                                never deviate to much from 1.0 so scaling 
C                                parameters to some faducial or base value is
C                                preferred.
C     (4) Integer ifl          - Not sure what this is for.
C     (5) Real*8 photar(ne)    - Output photon flux array. I will ultimately 
C                                return the value of the photon flux in each 
C                                energy bin (photons cm^-2 s^-1) to XSPEC.
C     (6) Real*8 photer(ne)    - (Optional) Output flux error array. Right 
C                                now, I put nothing in this array. Do I need 
C                                to put 1's or something here even if I do not 
C                                use it?
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine bwsim(ear,ne,param,ifl,photar,photer)
      Implicit None
      Integer ne,ifl
      Integer maxpnt
      Parameter (maxpnt=8192)
      Double Precision ear(0:maxpnt),param(12),photar(1:maxpnt),
     1               photer(1:maxpnt)
      Integer jen
      Integer ESwitch,ASwitch
      Double Precision acrat,comptem,caprad,bfield,dist
      Double Precision nsmass,nsrad,sigperp,sigpara,sigbar,alpha
      Double Precision thcrsec,kpc2cm
      Double Precision epskev,PI,pmass,emass,c,c2,gravc,boltzc
      Double Precision CBrem,CCyc,CBb,TCSpec
      Double Precision TSpOut(0:maxpnt)
      Double Precision chiabsCB,sigperpCB,thmdtemCB,xiCB,wCB,alphaCB,
     1                 deltaCB,kappaCB,tauthCB,taumaxCB,
     2                 epscyckeVCB,NSmassCB,NSRadiusCB,AccLumCB
      Common /BW07CB2/ chiabsCB,sigperpCB,thmdtemCB,xiCB,wCB,alphaCB,
     1                 deltaCB,kappaCB,tauthCB,taumaxCB,
     2                 epscyckeVCB,NSmassCB,NSRadiusCB,AccLumCB
      Integer NintCBB1,NintCBB2
      Common /BW07CB5/ NintCBB1,NintCBB2
      Logical StartUp
      Data StartUp /.true./
C
C     At start up, only once, print out a message.
C
      if (StartUp) then
          NintCBB1=15
          NintCBB2=30
          write(6,1000) NintCBB1,NintCBB2
 1000 format(/,"BWSIM Version 1.0.0:",
     1       /,"Similarity Variable Interface.",
     2       /,"GLQ CBB Integration: NintCBB1=",i4," NintCBB2=",i4,
     3       /,"February 2, 2026")
          StartUp=.false.
      endif
C
C     Set some fundamental physical constants
C
      PI=3.14159265358979323846264338327950d0
      pmass=1.6726d-24
      emass=9.1094d-28
      c=2.99792458d10
      c2=c*c
      gravc=6.67384d-8
      boltzc=1.3807d-16
      thcrsec=6.6525d-25             ! Thomson Cross Section; Units: cm^-2
      kpc2cm=1000.0d0*3.086d18       ! kiloparsec in centimeters
C
C     Put test parameters into physical variables here. This is the section
C     that works with the combination of the similarity and physical input 
C     variables.  For each variable I need to convert it into the units that 
C     the model internal coding expects each parameter to take. 
C
      acrat=param(1)    *1.0d17       ! Mass accretion rate (Convert to g/s)
      comptem=param(2)  *11604832.2d0 ! Compton temperature (Convert to Kelvins)
      caprad=param(3)   *1.0d2        ! Accretion cap radius (Convert to cm)
      bfield=param(4)   *1.0d12       ! Magnetic field strength (Convert to gauss)
      dist=param(5)     *3.086d21     ! Distance (Convert to cm) 
      nsmass=param(6)   *1.989e33     ! Neutron star mass (Convert to grams)
      nsrad=param(7)    *1.0d5        ! Neutron star radius (Convert to cm)
      sigperp=param(8)  *6.6525d-25   ! Perpendicular Scattering cross section (Convert to cm^2)
      alpha=1.33499124d0*gravc*nsmass*param(9)/(nsrad*c2)
      sigpara=(PI*caprad*pmass*c/(acrat*param(9)*DSqrt(sigperp)))**2
      sigbar=alpha*sigpara*emass*c2/(3.0d0*param(10)*boltzc*comptem)
      ESwitch=NInt(param(11))         ! Switch to turn emission processes on/off.
      ASwitch=NInt(param(12))         ! Accuracy and speed switch (high-slow/low-fast)
C
C     ASwitch=0 is high accuracy, slow execution (default); 
C     ASwitch=1 is low accuracy, fast execution 
C
C     Get setup to loop over the energy grid.
C
      do jen=1,ne
C
C     Set up the energy at the center of each energy bin to input to 
C     spectral routine.
C
          epskev=0.5*(ear(jen)+ear(jen-1))
C
C     Actually compute the spectrum. The subroutine takes in the energy it is 
C     to compute the spectral point at in kev. This subroutine returns the 
C     spectral flux density in photons cm^-2 s^-1 keV^-1 at the input energy 
C     of each bin.
C
          Call BWMod(acrat,comptem,caprad,bfield,dist,
     1               sigperp,sigpara,sigbar,
     1               nsmass,nsrad,ESwitch,ASwitch,
     1               epskev,CBrem,CCyc,CBb,TCSpec)
C
C     Fill in the output spectral array for later use to construct the flux 
C     in each energy bin.
C
          TSpOut(jen)=TCSpec
      end do
C
C     Now I need to do the output photon flux array. Right now, all I do
C     is multiply the spectral flux density by the energy interval and return
C     that as a flux. In the future we may want to revisit this and make 
C     the flux determination more accurate (?)
C
      do jen=1,ne
        photar(jen)=TSpOut(jen)*(ear(jen)-ear(jen-1))
      end do
C
C     All done. Git!
C
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Subroutine f_bwsim(ear,ne,param,ifl,photar,photer)
C
C     Purpose: f_bwsim: Name wrapper for double precision usage in ISIS
C              of the Similarity Variable interface.
c
C     Parameters:
C              ear(0:ne)     --- Double Precision Input: Energy Bin Boundaries.
C              ne            --- Integer Input: Size of flux array.
C              param(nparam) --- Double Precision Input: Parameter values.
C              ifl           --- Integer Input: Dead parameter sometimes used
C                                for specific XSPEC implementations.
C              photar(ne)    --- Double Precision Output: Photon flux array.
C              photer(ne)    --- Double Precision Output: (Optional) flux 
C                                error array.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine f_bwsim(ear,ne,param,ifl,photar,photer)
C
C     FS: f_bwsim: Name wrapper for double precision usage in ISIS
C     of the Similarity Variable interface.
C
      Implicit None
      Integer ne,ifl
      Integer maxpnt
      Parameter (maxpnt=8192)
      Double Precision ear(0:maxpnt),param(12)
      Double Precision photar(1:maxpnt),photer(1:maxpnt)
      Call bwsim(ear,ne,param,ifl,photar,photer)
      Return
      End
