CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C              Becker & Wolff (2007) model XSPEC Physical Variables code.    C
C                                                                            C
C              This code can run spectral models based on the Becker &       C
C              Wolff (2007, ApJ., v.654, p.435, Hereafter BW2007) model      C
C              for accretion onto a magnetic neutron star.                   C
C                                                                            C
C              Date: June 26, 2024                                           C
C              Author:  Michael T. Wolff                                     C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine bwphys(ear,ne,param,ifl,photar,photer)
C
C      Inputs to BW2007 model routine for XSPEC using the Physical Variable values:
C      (1)  Real*8 ear(0:ne)      - Energy array boundaries. Each pair of boundaries 
C                                   comprises an energy bin. I will ultimately return 
C                                   the value of the photon flux in each energy bin 
C                                   (photons cm^-2 s^-1) to XSPEC.
C      (2)  Integer ne            - Size of flux array. This is the number of input
C                                   points in the ear array. Note that the number of 
C                                   input energies in the variable ear is ne+1.
C      (3)  Real*8 param(nparam)  - Parameter values. Right now, there are ten 
C                                   parameters that specify a model. Also, each 
C                                   parameter value needs to be close to 1.0 and 
C                                   never deviate to much from 1.0 so scaling 
C                                   parameters to some faducial or base value is
C                                   preferred.
C      (4)  Integer ifl           - Not sure what this is for.
C      (5)  Real*8 photar(ne)     - Output photon flux array. I will ultimately 
C                                   return the value of the photon flux in each 
C                                   energy bin (photons cm^-2 s^-1) to XSPEC.
C      (6)  Real*8 photer(ne)     - (Optional) Output flux error array. Right now, 
C                                   I put nothing in this array. Do I need to put 
C                                   1's or something here even if I do not use it?
C
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
 1000 format(/,"BWPHYS Version 1.0.0:",
     1       /,"Physical Variable Interface.",
     2       /,"GLQ CBB Integration: NintCBB1=",i4," NintCBB2=",i4,
     3       /,"February 2, 2026")
        StartUp=.false.
      endif
C
C     Set some fundamental physical constants
C
      PI=3.14159265358979323846264338327950D0
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
C     that works with the explicit physical input variables.  For each variable 
C     I need to do a conversion into the units that the model internal coding 
C     expects each parameter to take. 
C
      acrat=param(1)    *1.0d17       ! Mass accretion rate (Convert to g/s)
      comptem=param(2)  *11604832.2d0 ! Compton temperature (Convert to Kelvins)
      caprad=param(3)   *1.0d2        ! Accretion cap radius (Convert to cm)
      bfield=param(4)   *1.0d12       ! Magnetic field strength (Convert to gauss)
      dist=param(5)     *3.086d21     ! Distance (Convert to cm) 
      nsmass=param(6)   *1.989e33     ! Neutron star mass (Convert to grams)
      nsrad=param(7)    *1.0d5        ! Neutron star radius (Convert to cm)
      sigperp=param(8)  *6.6525d-25   ! Perpendicular Scattering cross section (Convert to cm^2)
      sigpara=param(9)  *6.6525d-25   ! Parallel Scattering cross section (Convert to cm^2)
      sigbar=param(10)  *6.6525d-25   ! Average Scattering cross section (Convert to cm^2)
      ESwitch=NInt(param(11))         ! Switch to turn emission processes on/off.
      ASwitch=NInt(param(12))         ! Accuracy and speed switch (high-slow/low-fast)
C
C      ASwitch=0 is high accuracy, slow execution (default); 
C      ASwitch=1 is low accuracy, fast execution 
C
C
C     Get setup to loop over the energy grid.
C
      do jen=1,ne
C
C      Set up the energy at the center of each energy bin to input to spectral routine.
C
        epskev=0.5*(ear(jen)+ear(jen-1))
C
C     Actually compute the spectrum. The subroutine takes in the energy it is to compute
C     the spectral point at in kev. This subroutine returns the spectral flux density
C     in photons cm^-2 s^-1 keV^-1 at the input energy of each bin.
C
        Call BWMod(acrat,comptem,caprad,bfield,dist,
     1             sigperp,sigpara,sigbar,
     1             nsmass,nsrad,ESwitch,ASwitch,
     1             epskev,CBrem,CCyc,CBb,TCSpec)
C
C     Fill in the output spectral array for later use to construct the flux in each energy bin.
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
C     Subroutine f_bwphys(ear,ne,param,ifl,photar,photer)
C
C     Purpose: f_bwphys: Name wrapper for double precision usage in ISIS
C     of the Physical Variable interface.
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
      Subroutine f_bwphys(ear,ne,param,ifl,photar,photer)
C
C     FS: f_bwphys: Name wrapper for double precision usage in ISIS
C     of the Physical Variable interface.
C
      Implicit None
      Integer ne,ifl
      Integer maxpnt
      Parameter (maxpnt=8192)
      Double Precision ear(0:maxpnt),param(12)
      Double Precision photar(1:maxpnt),photer(1:maxpnt)
      Call bwphys(ear,ne,param,ifl,photar,photer)
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Subroutine BWMod(acrat,comptem,caprad,bfield,dist,sigperp,sigpara,sigbar,
C     1                nsmass,nsrad,ESwitch,ASwitch,epskev,CBrem,CCyc,CBb,TCSpec)
C
C     Purpose: This subroutine is the main event. It will take in numerous
C              input physical parameters and return the X-ray flux for a 
C              Becker&Wolff (2007) model in each energy bin.
C
C     Parameters:
C              acrat    --- Double Precision Input: Mass Accretion Rate
C              comptem  --- Double Precision Input: Compton Temperature
C              caprad   --- Double Precision Input: Accretion Cap Radius
C              bfield   --- Double Precision Input: Magnetic Field Strength
C              dist     --- Double Precision Input: Distance
C              sigperp  --- Double Precision Input: Perpendicular Scattering 
C                           cross section
C              sigpara  --- Double Precision Input: Parallel Scattering 
C                           cross section
C              sigbar   --- Double Precision Input: Average Scattering 
C                           cross section
C              nsmass   --- Double Precision Input: Neutron Star Mass
C              nsrad    --- Double Precision Input: Neutron Star Radius
C              ESwitch  --- Integer Input: Switch to turn emission processes on/off.
C              ASwitch  --- Integer Input: Accuracy/Speed switch (high-slow/low-fast)
C              epskev   --- Double Precision Input: Energy for calculation.
C              CBrem    --- Double Precision Output: Bremsstrahlung spectral 
C                           component.
C              CCyc     --- Double Precision Output: Cyclotron spectral component.
C              CBb      --- Double Precision Output: Blackbody spectral component.
C              TCSpec   --- Double Precision Output: Total spectral flux.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine BWMod(acrat,comptem,caprad,bfield,dist,
     1                 sigperp,sigpara,sigbar,
     1                 nsmass,nsrad,ESwitch,ASwitch,
     1                 epskev,CBrem,CCyc,CBb,TCSpec)
      Implicit None
      Integer maxpnt
      Parameter (maxpnt=8192)
      double precision pmass,c,c2,gravc,boltzc,emass,plankc,plankc3
      double precision echrg,erg2kev,erg2ev,kev2erg,ev2erg
      double precision acrat,comptem,caprad,bfield,dist,chiabs,thcrsec
      double precision nsmass,nsrad,sigperp,sigpara,sigbar
      double precision epscyc,chicyc,thmdtem,xi,w,alpha,delta,kappa
      double precision thmdrho,thmdu,umax,tauth,c1,zmax,taumax,epskev
      double precision eps,chi,epsabs,chimin,chimax
      double precision BremCon1,lamlop,mun,brsum,WhitW,WhitM
      double precision Bfunc2,Xfunc,Afunc,FctrlExt
      double precision PI,kpc2cm
      Double Precision CycCon1,HCycFunc,cycsum,delcyc
      Double Precision BBCon1,bbsum,delbb,delbr
      Double Precision CBrem,CCyc,CBb,TCSpec
      Double Precision gn
      Double Precision gamfdum1,gamfdum2,gamfdum3
      Integer nloop
      Integer ESwitch,ASwitch
      Double Precision epsCB,TeCompCB,KapCB,MuCB,TeThermCB
      Common /BW07CB1/ epsCB,TeCompCB,KapCB,MuCB,TeThermCB
      Double Precision BBintLL,BBout1,BBintUL,BBout2
      Double Precision BBkernW,BBkernM
      Double Precision GammaExt
      Double Precision GLabsc(256),GLwait(256),BBout1GL,BBout2GL
      Integer gj
      Double Precision chiabsCB,sigperpCB,thmdtemCB,xiCB,wCB,alphaCB,
     1                 deltaCB,kappaCB,tauthCB,taumaxCB,
     2                 epscyckeVCB,NSmassCB,NSRadiusCB,AccLumCB
      Common /BW07CB2/ chiabsCB,sigperpCB,thmdtemCB,xiCB,wCB,alphaCB,
     1                 deltaCB,kappaCB,tauthCB,taumaxCB,
     2                 epscyckeVCB,NSmassCB,NSRadiusCB,AccLumCB
      Double Precision BBWDum,BBMDum,BBIDum
      Logical ldobrm,ldocyc,ldobb
      Integer NintCBB1,NintCBB2
      Common /BW07CB5/ NintCBB1,NintCBB2
      External BBkernW,BBkernM
      Common /BW07CB6/ EWconv,EMconv
      Common /BW07CB7/ Bloop,Cloop,BBloop
      Integer Bloop,Cloop,BBloop
      Double Precision EWconv,EMconv
      Double Precision BBtrap
      External BBtrap
C
C     Module to evaluate the Becker Wolff (2007) model for the spectrum of
C     an accreting X-ray pulsar.
C
C     Set some fundamental physical constants
C
      PI=3.14159265358979323846264338327950D0
      pmass=1.6726d-24
      c=2.99792458d10
      c2=c*c
      gravc=6.67384d-8
      boltzc=1.3807d-16
      emass=9.1094d-28
      plankc=6.6260755d-27
      plankc3=plankc*plankc*plankc
      echrg=4.8032068d-10
      erg2kev=6.2414d8
      erg2ev=6.2414d11 
      kev2erg=1.0/erg2kev
      ev2erg=1.0/erg2ev
      thcrsec=6.6525d-25             ! Thomson Cross Section; Units: cm^-2
      kpc2cm=1000.0d0*3.086d18       ! kiloparsec in centimeters
C
C      ASwitch=0 is high accuracy, slow execution; 
C      ASwitch=1 is low accuracy, fast execution 
C
      if (ASwitch.eq.0) then
        Bloop=10        ! ASwitch=0: High accuracy + slow convergence
        Cloop=10        ! ASwitch=0: High accuracy + slow convergence
        BBloop=10       ! ASwitch=0: High accuracy + slow convergence
        EWconv=1.0D-6   ! ASwitch=0: High accuracy + slow convergence
        EMconv=1.0D-6   ! ASwitch=0: High accuracy + slow convergence
      else
        Bloop=4         ! ASwitch=1: Low accuracy + fast convergence
        Cloop=4         ! ASwitch=1: Low accuracy + fast convergence
        BBloop=4        ! ASwitch=1: Low accuracy + fast convergence
        EWconv=1.0D-3   ! ASwitch=1: Low accuracy + fast convergence
        EMconv=1.0D-3   ! ASwitch=1: Low accuracy + fast convergence
      endif
C
C     Setup default physics switches here:
C
      ldobrm=.true.                   ! Turn Bremsstrahlung Comptonization on.
      ldocyc=.true.                   ! Turn Cyclotron Comptonization on.
      ldobb=.false.                   ! Turn BB Comptonization off.
C
C     Which components of the model do I compute? Parse the ESwitch parameter.
C
      if (ESwitch.eq.0) then
        ldobrm=.true.                   ! Turn Bremsstrahlung Comptonization on.
        ldocyc=.true.                   ! Turn Cyclotron Comptonization on.
        ldobb=.true.                    ! Turn BB Comptonization on.
        if (bfield.le.0.0d0) then
          ldocyc=.false.
        endif
      endif
      if (ESwitch.eq.1) then
        ldobrm=.true.                   ! Turn Bremsstrahlung Comptonization on.
        ldocyc=.true.                   ! Turn Cyclotron Comptonization on.
        ldobb=.false.                   ! Turn BB Comptonization off.
        if (bfield.le.0.0d0) then
          ldocyc=.false.
        endif
      endif
      if (ESwitch.eq.2) then
        ldobrm=.true.                    ! Turn Bremsstrahlung Comptonization on.
        ldocyc=.false.                   ! Turn Cyclotron Comptonization off.
        ldobb=.true.                     ! Turn BB Comptonization on.
      endif
      if (ESwitch.eq.3) then
        ldobrm=.false.                   ! Turn Bremsstrahlung Comptonization off.
        ldocyc=.true.                    ! Turn Cyclotron Comptonization on.
        ldobb=.true.                     ! Turn BB Comptonization on.
        if (bfield.le.0.0d0) then
          ldocyc=.false.
        endif
      endif
      if (ESwitch.eq.-1) then
        ldobrm=.true.                     ! Turn Bremsstrahlung Comptonization on.
        ldocyc=.false.                    ! Turn Cyclotron Comptonization off.
        ldobb=.false.                     ! Turn BB Comptonization off.
      endif
      if (ESwitch.eq.-2) then
        ldobrm=.false.                    ! Turn Bremsstrahlung Comptonization off.
        ldocyc=.true.                     ! Turn Cyclotron Comptonization on.
        ldobb=.false.                     ! Turn BB Comptonization off.
        if (bfield.le.0.0d0) then
          ldocyc=.false.
        endif
      endif
      if (ESwitch.eq.-3) then
        ldobrm=.false.                   ! Turn Bremsstrahlung Comptonization off.
        ldocyc=.false.                   ! Turn Cyclotron Comptonization off.
        ldobb=.true.                     ! Turn BB Comptonization on.
      endif
C
C     Now start the real calculation of the spectral point. first,
C     I need a number of important parameters that depend on the input
C     parameters of the model. Note: rather than making chiabs an
C     input parameter I calculate it based on the other input parameters. 
C
      epscyc=echrg*bfield*plankc/(2.0*PI*emass*c)
      chicyc=epscyc/(boltzc*comptem)
      thmdtem=2.318d3*(acrat**0.4)*(caprad**(-0.666666667))
      thmdrho=4.05d-12*(thmdtem**1.75)/DSqrt(caprad)
      thmdu=2.62d0*acrat/((caprad**1.5)*(thmdtem**1.75))
      xi=PI*caprad*pmass*c/(acrat*dsqrt(sigpara*sigperp))
      w=dsqrt(9.0d0+12.0d0*xi*xi)
      alpha=32.0d0*dsqrt(3.0d0)*gravc*nsmass*xi/
     1      (49.0d0*dlog(7.0d0/3.0d0)*nsrad*c2)
      delta=alpha*sigpara*emass*c2/(3.0d0*sigbar*boltzc*comptem)
      kappa=0.5d0*(delta+4.0d0)
      tauth=2.64d28*acrat*nsrad/(nsmass*(caprad**1.5)*
     1      (thmdtem**1.75)*xi)
      c1=4.0*gravc*nsmass*caprad*xi/(alpha*c2*nsrad*
     1      nsrad)*dsqrt(sigperp/sigpara)
      zmax=nsrad*(dsqrt(1.0d0+c1)-1.0d0)/2.0d0
      taumax=dsqrt(dsqrt(sigpara/sigperp))*dsqrt(2.0d0*zmax/
     1       (alpha*xi*caprad))
      umax=DSqrt(2.0d0*gravc*nsmass/(nsrad+zmax))
      epsabs=(6.08d12*DSqrt(thmdrho)/(comptem**1.75))*
     1       DSqrt(DSqrt(thmdu*c/umax))*boltzc*comptem
      chiabs=epsabs/(boltzc*comptem)
C
C     Transfer numbers to common block for black body integration module /BW07CB1/:
C
      TeCompCB=comptem
      KapCB=kappa
      TeThermCB=thmdtem
C 
C     Set up energy variables.
C
      eps=epskev*kev2erg
      epsCB=eps
      chi= eps/(boltzc*comptem)
      chimin=Min(chi,chicyc)
      chimax=Max(chi,chicyc)
C
C     Set the integration energy limits for Comptonized black body
C     integrals. 
C
      BBintLL=0.01d0*kev2erg
      BBintUL=200.0d0*kev2erg
      BBintUL=min(200.0d0*kev2erg,20.0d0*epskev*kev2erg)
C
C     Compute the Bremsstrahlung constant.
C
      brsum=1.0d-100
      if (ldobrm) then
        BremCon1=2.80D-12*acrat*xi*xi*DSqrt(alpha*alpha*alpha*w)*
     1           (eps**(kappa-2.d0))*
     2           DExp(-eps/(2.d0*boltzc*comptem))
     3           /(sigbar*((boltzc*comptem)**(kappa+.5d0)))
C
C     Loop through the first 11 terms of the expansion to output the 
C     Bremsstrahlung spectral value at this energy.
C
        do 100 nloop=0,Bloop
          lamlop=(4.0d0*dble(nloop)*w+w+3.0d0)/2.0d0
          mun=dsqrt((3.d0-delta)*(3.d0-delta)+4.d0*delta*lamlop)/
     1           2.d0
          gamfdum1=GammaExt(mun-kappa+0.5d0)
          gamfdum2=GammaExt(1.0d0+2.0d0*mun)
          gamfdum3=GammaExt(dble(nloop)+0.5d0)
          delbr=BremCon1*gamfdum1*((FctrlExt(nloop))/
     1           (gamfdum2*gamfdum3))*
     2           (Xfunc(nloop,w,alpha)*
     3           Afunc(nloop,w,alpha,tauth,taumax)*
     4           Bfunc2(kappa,mun,chiabs,chi))
          brsum=brsum+delbr
 100    continue
      endif
C
C     Compute the Cyclotron emisison front constant:
C
      cycsum=1.0d-100
      if (ldocyc) then
        CycCon1=3.43D-16*acrat*HCycFunc(chicyc)*xi*xi*
     1          DSqrt(alpha*alpha*alpha*w)*
     2          ((eps)**(kappa-2.0d0))/
     3          (sigbar*((epscyc)**(kappa+1.5d0))*
     4          DExp((epscyc+eps)/(2.0d0*boltzc*comptem)))
C
C     Loop to fill in the Comptonized cyclotron emisison:
C 
        do 110 nloop=0,Cloop
          lamlop=(4.0d0*dble(nloop)*w+w+3.0d0)/2.0d0
          mun=dsqrt((3.d0-delta)*(3.d0-delta)+4.d0*delta*lamlop)/2.d0
          gamfdum1=GammaExt(mun-kappa+0.5d0)
          gamfdum2=GammaExt(2.0d0*mun+1.0d0)
          gamfdum3=GammaExt(nloop+0.5d0)
          delcyc=CycCon1*(gamfdum1*FctrlExt(nloop)/(gamfdum2*
     1           gamfdum3))*Afunc(nloop,w,alpha,tauth,taumax)* 
     2           Xfunc(nloop,w,alpha)*
     3           WhitM(kappa,mun,chimin)*WhitW(kappa,mun,chimax)
          cycsum=cycsum+delcyc
          if ( nloop.gt.5 .and. 
     1          (DAbs(delcyc).lt.1.0d-4*DAbs(cycsum).or.
     1          DAbs(cycsum).lt.1.0d-4*DAbs(brsum)) ) goto 111
 110      continue
 111    continue
      endif
C
C     Compute black body integration front constant. Also, initialize some
C     variables so that I get rid of some annoying compilers warnings.
C     Contains Pete's typo correction.
C
      bbsum=1.0d-100
      if ( ldobb ) then
C
C     Execute second option for BB integration using trapezoid
C     integration scheme.
C
        bbsum=BBtrap(thmdtem,comptem,alpha,xi,delta,tauth,
     1               eps,caprad,BBloop)
      end if
C
C     Begin final processing of spectrum: Apply the distance to each
C     component.
C
      TCSpec=(brsum+cycsum+bbsum)*kev2erg/(4.0d0*PI*dist*dist)
      CBrem=brsum*kev2erg/(4.0d0*PI*dist*dist)
      CCyc=cycsum*kev2erg/(4.0d0*PI*dist*dist)
      CBb=bbsum*kev2erg/(4.0d0*PI*dist*dist)
C
C     Transfer some derived parameters to common block transfer storage /BW07CB2/:
C
      chiabsCB=chiabs
      sigperpCB=sigperp
      thmdtemCB=thmdtem
      xiCB=xi
      wCB=w
      alphaCB=alpha
      deltaCB=delta
      kappaCB=kappa
      tauthCB=tauth
      taumaxCB=taumax
      epscyckeVCB=epscyc*erg2kev
      NSmassCB=nsmass
      NSRadiusCB=nsrad
      AccLumCB=acrat*nsmass*gravc/nsrad
C
C     All done! Return!
C     
      Return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function FctrlExt(iin)
C
C     Purpose: This routine will calculate and return in Double Precision the 
C              value of n!. This routine utilizes the fundamental relation 
C              between the factorial function and the Gamma function if the 
C              value of the input integer is very large.
C
C     Parameters:
C              iin     --- Input: Input integer. Must be positive!
C
C     Returns: Double Precision value of FctrlExt(iin)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function FctrlExt(iin)
      Implicit None
      Integer iin,i
      Logical First,Pading   ! Pading is 32 bits to pad the common block
      Double Precision Result(0:32),lngamz
      Common /FacCom/ First,Pading,Result
      Data First,Result(0) /.True.,1.0d0/
C
C     Setup the Result array for later use. Set the logical First-call
C     variable to be a flag for subsequent calls to this function.
C
      If (First) Then
          Result(0)=1.0d0
          First=.False.
          Do i=1,32
              Result(i)=Result(i-1)*Dble(i)
          EndDo
      EndIf
C
C     If a negative integer is input then just return an error and a 
C     very large number.
C
      If (iin.lt.0) Then
          FctrlExt=1.0d90
          Return
      EndIf
C
C     Now set the answer and return it as a double precision real variable.
C
      If (IIN.le.32) Then
          FctrlExt=Result(iin)
          Return
      Else
          FctrlExt=Dexp(lngamz(dble(iin)+1.0d0))
          Return
      EndIf
      Return
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function BBkernM(eps0)
C
C     Purpose: This function evaluates the blackbody integration kernel at a
C              particular energy. See Equation 73 of BW2007.
C
C     Parameters:
C              eps0     --- Double Precision Input: Energy.
C 
C     Return : Double Precision value of BBkernM.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function BBkernM(eps0)
      Implicit None
      Double Precision eps0
      Double Precision epsCB,TeCompCB,KapCB,MuCB,TeThermCB
      Double Precision WhitM
      Common /BW07CB1/  epsCB,TeCompCB,KapCB,MuCB,TeThermCB
      Double Precision epdum0,HG
      WhitM=1.0d0
      epdum0=eps0/(1.3807d-16*TeCompCB)
      Call CHGM(0.5d0+MuCB-KapCB,1.0d0+2.0d0*MuCB,epdum0,HG)
      BBkernM=HG*(epdum0**(0.5d0+MuCB))*(eps0**(2.0d0-KapCB))/
     1        (DExp(eps0/(1.3807d-16*TeThermCB))-1.0d0)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function BBkernW(eps0)
C
C     Purpose: This function evaluates the blackbody integration kernel at a
C              particular energy. See Equation 73 of BW2007.
C
C     Parameters:
C              eps0     --- Double Precision Input: Energy.
C 
C     Return : Double Precision value of BBkernW.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function BBkernW(eps0)
      Implicit None
      Double Precision eps0
      Double Precision epsCB,TeCompCB,KapCB,MuCB,TeThermCB
      Double Precision WhitW
      Common /BW07CB1/  epsCB,TeCompCB,KapCB,MuCB,TeThermCB
      Integer MD
      Double Precision epdum0,HU
      WhitW=1.0d0
      epdum0=eps0/(1.3807d-16*TeCompCB)
      Call CHGU(0.5d0+MuCB-KapCB,1.0d0+2.0d0*MuCB,epdum0,HU,MD)
      BBkernW=HU*(epdum0**(0.5d0+MuCB))*(eps0**(2.0d0-KapCB))/
     1        (DExp(eps0/(1.3807d-16*TeThermCB))-1.0d0)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function Bfunc2(kappa,mun,chiabs,chi) 
C
C     Purpose: This function evaluates Equation 129 of BW2007.
C              This function utilizes the expansion version of Equation 129 
C              of BW2007. The input variables kappa and mun must be evaluated 
C              using Equation 47 of BW 2007. 
C
C     Parameters:
C              kappa     --- Double Precision Input Parameter.
C              mun       --- Double Precision Input Parameter.
C              chiabs    --- Double Precision Input: Normalized bremsstrahlung
C                            cutoff energy.
C              chi       --- Double Precision Input: Normalized energy.
C              
C     Return : Double Precision value of Bfunc2.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function Bfunc2(kappa,mun,chiabs,chi) 
      Implicit None
      double precision kappa,mun,chiabs,chi
      double precision IMfunc,IWfunc,WhitM,WhitW
      If( chi.le.chiabs ) then
        Bfunc2=WhitM(kappa,mun,chi)*IWfunc(kappa,mun,chiabs)
        return
      else
        Bfunc2=WhitW(kappa,mun,chi)*(IMfunc(kappa,mun,chi)-
     1         IMfunc(kappa,mun,chiabs))+
     2         WhitM(kappa,mun,chi)*IWfunc(kappa,mun,chi)
        return
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function IMfunc(kappa,mun,chi)
C
C     Purpose: This program computes part of the bremsstrahlung emission
C              integral expansion.
C
C     Parameters:
C              kappa    --- Double Precision Input Parameter
C              mun      --- Double Precision Input Parameter
C              chi      --- Double Precision Input Normalized Energy.
C
C     Returns:  Double Precision value of IMfunc (bremsstrahlung function)
C
C     Routines called: LMfunc, EMfunc.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function IMfunc(kappa,mun,chi)
      Implicit None
      double precision kappa,mun,chi
      double precision LMfunc,EMfunc
      if ( chi.lt.20.0d0 ) then
        IMfunc=LMfunc(kappa,mun,chi)
      else
        IMfunc=EMfunc(kappa,mun,chi)
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function IWfunc(kappa,mun,z)
C
C     Purpose: This program computes part of the bremsstrahlung emission
C              integral expansion.
C
C     Parameters:
C              kappa    --- Double Precision Input Parameter
C              mun      --- Double Precision Input Parameter
C              z        --- Double Precision Input Parameter
C
C     Returns:  Double Precision value of IWfunc (bremsstrahlung function)
C
C     Routines called: LWfunc, EWfunc.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function IWfunc(kappa,mun,z)
      Implicit None
      double precision kappa,mun,z
      double precision LWfunc,EWfunc
      if ( z.lt.15.0d0 ) then
        IWfunc=LWfunc(kappa,mun,z)
        return
      else
        IWfunc=EWfunc(kappa,mun,z)
        return
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function EWfunc(kappa,mun,z)
C
C     Purpose: This program computes part of the bremsstrahlung emission
C              integral expansion.
C
C     Parameters:
C              kappa    --- Double Precision Input Parameter
C              mun      --- Double Precision Input Parameter
C              z        --- Double Precision Input Parameter
C
C     Returns:  Double Precision value of EWfunc (bremsstrahlung function)
C
C     Routines called: CHGU, GammaExt.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function EWfunc(kappa,mun,z)
      Implicit None
      double precision kappa,mun,z
      double precision sum,front,gfout
      Double Precision HU1,fac1,gfout1
      double precision faclast
      Double Precision GammaExt
      integer MD,m
      Common /BW07CB6/ EWconv,EMconv
      Double Precision EWconv,EMconv
      gfout=GammaExt(mun-kappa+0.5d0)
      front=(z**(mun-kappa-0.5d0))*dexp(-z)*gfout 

C      write(6,1002) z,mun,kappa,mun-kappa+0.5d0
C 1002 format(1x,/,1x,"In EWfunc: z= ",f15.11,
C     2       " mun= ",f15.11," kappa= ",f15.11,
C     3       " mun-kappa+0.5= ",f15.11)
 
      sum=0.0d0
      faclast=0.0d0
      do m=0,50,1
        call CHGU(mun-kappa+0.5d0,2.0d0*mun-dble(m),z,
     1            HU1,MD)
        gfout1=GammaExt(mun-kappa-dble(m)+0.5d0)
        fac1=HU1/((z**m)*gfout1)
        sum=sum+fac1
        
C        write(6,1001) m,MD,HU1,fac1,
C     1                 dabs((faclast+fac1)/sum),sum
C 1001 format(1x,"EWfunc Loop:",2i5,5(1x,1pe15.8))

        if (dabs((faclast+fac1)/sum).lt.EWconv) then
          EWfunc=sum*front
          return
        endif
        faclast=fac1
      end do

      write(6,1000) z,dabs((faclast+fac1)/sum),mun,kappa
 1000 format(1x,"No convergence in EWfunc: z= ",f15.11,
     1 " abs((faclast+fac1)/sum)= ",1pe18.11,
     2 " mun= ",0pf15.11," kappa= ",f15.11)

      EWfunc=sum*front
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function LWfunc(kappa,mun,chi)
C
C     Purpose: This program computes part of the bremsstrahlung emission
C              integral expansion.
C
C     Parameters:
C              kappa    --- Double Precision Input Parameter
C              mun      --- Double Precision Input Parameter
C              chi      --- Double Precision Input Normalized Energy.
C
C     Returns:  Double Precision value of LWfunc (bremsstrahlung function)
C
C     Routines called: HG2F2, GammaExt
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function LWfunc(kappa,mun,chi)
      Implicit None
      double precision kappa,mun,chi
      double precision HGpFq2,HGpFq3
      double precision fac1,fac2,fac3
      double precision fac1ga1,fac1ga2,fac1ga3
      double precision fac2ga1,fac2ga2
      double precision fac3ga1,fac3ga2
      Double Precision GammaExt
      call HG2F2(mun+kappa+0.5d0,mun-kappa+0.5d0,
     1           1.0d0+2.0d0*mun,mun-kappa+1.5d0,-chi,
     2           HGpFq2)
      call HG2F2(-mun+kappa+0.5d0,-mun-kappa+0.5d0,
     1           1.0d0-2.0d0*mun,-mun-kappa+1.5d0,-chi,
     2           HGpFq3)
      fac1ga1=GammaExt(mun-kappa+0.5d0)
      fac1ga2=GammaExt(-mun-kappa+0.5d0)
      fac1ga3=GammaExt(1.0d0-2.0d0*kappa)
      fac1=fac1ga1*fac1ga2/fac1ga3
      fac2ga1=GammaExt(-2.0d0*mun)
      fac2ga2=GammaExt(-mun-kappa+0.5d0)
      fac2=(fac2ga1/fac2ga2)*(chi**(mun-kappa+0.5d0))*HGpFq2/
     1     (mun-kappa+0.5d0)
      fac3ga1=GammaExt(2.0d0*mun)
      fac3ga2=GammaExt(mun-kappa+0.5d0)
      fac3=(fac3ga1/fac3ga2)*(chi**(-mun-kappa+0.5d0))*HGpFq3/
     1     (mun+kappa-0.5d0)
      LWfunc=fac1-fac2+fac3
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function LMfunc(kappa,mun,chi)
C
C     Purpose: This program computes part of the bremsstrahlung emission
C              integral expansion.
C
C     Parameters:
C              kappa    --- Double Precision Input Parameter
C              mun      --- Double Precision Input Parameter
C              chi      --- Double Precision Input Normalized Energy.
C
C     Returns:  Double Precision value of LMfunc (bremsstrahlung function)
C
C     Routine called: HG2F2.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function LMfunc(kappa,mun,chi)
      Implicit None
      double precision kappa,mun,chi
      double precision HGpFq
      call HG2F2(mun+kappa+0.5d0,mun-kappa+0.5d0,
     1           1.0d0+2.0d0*mun,mun-kappa+1.5d0,-chi,
     2           HGpFq)
      LMfunc=HGpFq*(chi**(mun-kappa+0.5d0))/(mun-kappa+0.5d0)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function EMfunc(kappa,mun,chi)
C
C     Purpose: This program computes part of the bremsstrahlung emission
C              integral expansion.
C
C     Parameters:
C              kappa    --- Double Precision Input Parameter
C              mun      --- Double Precision Input Parameter
C              chi      --- Double Precision Input Normalized Energy.
C
C     Returns:  Double Precision value of EMfunc (bremsstrahlung function)
C
C     Routine called: GammaExt, FactrlExp
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function EMfunc(kappa,mun,chi)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      Implicit None
      double precision kappa,mun,chi
      double precision fac1,poch1,poch2,del,sum
      double precision FctrlExt,NLoInGamFunc
      Double Precision gamfdum1,gamfdum2,gamfdum3,gamfdum4,gamfdum5
      Double Precision GammaExt,pochrat,facrat
      integer mint
      Common /BW07CB6/ EWconv,EMconv
      Double Precision EWconv,EMconv
      Logical ieee_is_nan
      sum=0.0d0
      poch1=1.0d0
      poch2=1.0d0
      do mint=0,100
        gamfdum1=GammaExt(mun-kappa+0.5d0+dble(mint))

C        gamfdum1=GammaExt(mun-kappa+0.5d0+dble(mint))
C        gamfdum2=GammaExt(mun-kappa+0.5d0+dble(mint))
C        gamfdum2=gamfdum1
C        gamfdum3=GammaExt(mun-kappa+0.5d0)
C        gamfdum4=GammaExt(1.0d0+2.0d0*mun+dble(mint))
C        gamfdum5=GammaExt(1.0d0+2.0d0*mun)
C        fac1=NLoInGamFunc(mun-kappa+0.5d0+dble(mint),chi)*gamfdum1
C        poch1=gamfdum2/gamfdum3
C        poch2=gamfdum4/gamfdum5
C        del=fac1*poch1/(FctrlExt(mint)*poch2)

        fac1=NLoInGamFunc(mun-kappa+0.5d0+dble(mint),chi)*gamfdum1
      if (mint.ne.0) then
            poch1=poch1*(mun-kappa+0.5d0+dble(mint)-1.0d0)
            poch2=poch2*(1.0d0+2.0d0*mun+dble(mint)-1.0d0)
      endif

        if (fac1.gt.1.0d100 .or. mint.gt.70) then
          del=Dlog(poch1)+Dlog(fac1)-Dlog(poch2)-Dlog(FctrlExt(mint))
          del=Dexp(del)
        else
          pochrat=poch1/poch2
          facrat=fac1/FctrlExt(mint)
          del=facrat*pochrat
        endif

        if (ieee_is_nan(del)) then
            goto 299
        endif

        sum=sum+del
        if(Dabs(del/sum).lt.EMconv) then
          EMfunc=sum
          return
        endif
      end do

 299  continue
      write(6,1000) mint,kappa,mun,chi,del,sum
 1000 format(1x,/,1x,
     1 "Problem in EMfunc: mint,kappa,mun,chi,del,sum= ",
     1       i6,5(1x,1pe14.7))
      write(6,1001) fac1,poch1,FctrlExt(mint),poch2
 1001 format(1x,"fac1,poch1,FctrlExt(mint),poch2= ",
     1       4(1x,1pe14.7))
      write(6,1002) gamfdum1,pochrat,facrat
 1002 format(1x,"gamfdum1,pochrat,facrat= ",3(1x,1pe14.7))
      write(6,1003) mun-kappa+0.5d0+dble(mint),
     1   1.0d0+2.0d0*mun+dble(mint),Dabs(del/sum),EMconv
 1003 format(1x,"mun-kappa+1/2+mint,",
     1 "1+2mun+mint,Dabs(del/sum),EMconv= ",
     2 4(1x,1pe14.7))

      EMfunc=sum
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function HCycFunc(x)
C
C     Purpose: This is the Arons, Klein, and Lea (1987) special cyclotron function.
C
C     Parameters:
C              x        --- Double Precision Input Normalized Energy.
C
C     Returns:  Double Precision value of HCycFunc
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function HCycFunc(x)
      Implicit None
      Double Precision x
      if ( x.gt.7.5d0 ) then
          HCycFunc=0.410792d0
      else
          HCycFunc=0.15d0*dsqrt(x)
      end if
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function GammaExt(X)
C
C     Purpose: This is Gamma Function of the input parameter x. This subroutine
C              evaluates the expansion given in Abramowitz & Stegun Table 6.1.34.
C
C     Parameters: 
C              x        --- Double Precision Input Parameter.
C                           Can not equal 0, -1, -2, etc.
C
C     Returns: Double Precision value of GammaExt(X)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function GammaExt(X)
      Implicit None
      Double Precision CK(26)
      Double Precision X,GA,Z,R,GR,PI
      Integer M1,K,M,NX
C
C     Set up the coefficients from Abramowitz & Stegun Table 6.1.34:
C
      DATA CK/ 1.0D0,0.5772156649015329D0,
     1       -0.6558780715202538D0, -0.420026350340952D-1,
     2        0.1665386113822915D0,-0.421977345555443D-1,
     3       -0.96219715278770D-2, 0.72189432466630D-2,
     4       -0.11651675918591D-2, -0.2152416741149D-3,
     5        0.1280502823882D-3, -0.201348547807D-4,
     6       -0.12504934821D-5, 0.11330272320D-5,
     7       -0.2056338417D-6, 0.61160950D-8,
     8        0.50020075D-8, -0.11812746D-8,
     9        0.1043427D-9, 0.77823D-11,
     &       -0.36968D-11, 0.51D-12,
     1       -0.206D-13, -0.54D-14, 
     2        0.14D-14, 0.1D-15 /
      PI=3.14159265358979323846264338327950D0
C
C      Check to see if the input parameter is an integer:
C
      IF ( X.EQ.Int(X) ) THEN
C
C     We know that X has integer value. Now check to see if 
C     X is positive or negative. If X is positive then the
C     Gamma function of X is just X-1 factorial. If X is negative 
C     then the Gamma function of X is not defined and we set 
C     the returned value to a very large numbner and just return.
C
        IF ( X.GT.0.0D0 ) THEN
          GA=1.0D0
          NX=Int(X)
          M1=NX-1
          DO 10 K=2,M1
            GA=GA*Dble(K)
 10       CONTINUE
        ELSE
          GA=1.0D+300
        ENDIF
      ELSE
C
C      If I make it here then the input parameter is not an integer.
C      Check to see if the input parmater is greater or less than 1.0.
C
        IF ( DAbs(X).GT.1.0D0 ) THEN
          Z=DAbs(X)
          M=Int(Z)
          R=1.0D0
          DO 15 K=1,M
            R=R*(Z-K)
 15       CONTINUE
          Z=Z-M
        ELSE
          Z=X
        ENDIF
C
C     Set up series expansion for inverse of Gamma(z):
C
        GR=CK(26)
        DO 20 K=25,1,-1
          GR=GR*Z+CK(K)
 20     CONTINUE
        GA=1.0D0/(GR*Z)
        IF ( DAbs(X).GT.1.0D0 ) THEN
          GA=GA*R
C
C     Check: Is argument less than 0.0 and if it is then use reflection formula:
C
          IF ( X.LT.0.0D0 ) GA=-PI/(X*GA*DSIN(PI*X))
        ENDIF
      ENDIF
      GammaExt=GA    ! Transfer answer to output variable name.
C
C     All done: Git!
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function Afunc(n,w,alpha,tauth,taumax)
C
C     Purpose: This function does the vertical spatial integration of the Laguerre 
C     functions in equation 119 of BW2007. This function is really just equation 
C     119 from BW2007.
C
C     Parameters: 
C              n        --- Integer Input Parameter. 
C              w        --- Double Precision Input Parameter
C              alpha    --- Double Precision Input Parameter.
C              tauth    --- Double Precision Input Parameter
C              taumax   --- Double Precision Input Parameter.
C
C     Returns:  Double Precision value of Afunc 
C
C     Routine called: GammaExt, FactrlExp, Exp1Int.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function Afunc(n,w,alpha,tauth,taumax)
      implicit none
C
C     This function does the vertical spatial integration of the Laguerre 
C     functions in equation 119 of BW2007. This function is really just equation 
C     119 from BW2007.
C
      integer m,n
      double precision xm,xn,x,a
      double precision w,alpha,tauth,taumax
      double precision gammcf1,gammcf2
      double precision gln1,gln2,gmpre,dd,sum
      double precision FctrlExt,Exp1Int
      double precision NUpInGamFunc
      Double Precision GammaExt
      Double Precision gamfdum1
      xn=dble(n)
      dd=xn+0.5d0
      gamfdum1=GammaExt(dd)
      gmpre=0.5d0*gamfdum1
      sum=0.0d0
      do m=0,n
        xm=dble(m)
        a=dble(m)
        x=alpha*tauth*tauth*(w-3.0d0)/4.0d0
        if (m.ne.0) then
          gln1=Log(GammaExt(a))

          gammcf1=NUpInGamFunc(a,x)*dexp(gln1)

        else
          gammcf1=Exp1Int(x)
        end if
        x=alpha*taumax*taumax*(w-3.0d0)/4.0d0
        if (m.ne.0) then
          gln2=Log(GammaExt(a))

          gammcf2=NUpInGamFunc(a,x)*dexp(gln2)

        else
           gammcf2=Exp1Int(x)
        end if 
        gamfdum1=GammaExt(xm+0.5d0)
        sum=sum+(gmpre*((2.0d0*w/(3.0d0-w))**m)/
     1      (gamfdum1*FctrlExt(m)*FctrlExt(n-m)))*
     2      (gammcf1-gammcf2)
      end do
      Afunc=sum
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function Xfunc(n,w,alpha)
C
C     Purpose: This function does the vertical spatial integration of the Laguerre 
C     functions in equation 119 of BW2007. This function is really just equation 
C     75 from BW2007.
C
C     Parameters: 
C              n        --- Integer Input Parameter. 
C              w        --- Double Precision Input Parameter
C              alpha    --- Double Precision Input Parameter.
C
C     Returns:  Double Precision value of Xfunc 
C
C     Routine called: GammaExt, FactrlExp.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function Xfunc(n,w,alpha)
      Implicit None
      integer n
C
C     Calculate equation 75 of BW2007.
C
      double precision xn,w,alpha,FctrlExt
      Double Precision gamfdum1
      double Precision GammaExt
      xn=dble(n)
      gamfdum1=GammaExt(xn+0.5d0)
      xfunc=2.0d0*gamfdum1*((3.0d0-w)**(xn-1.0d0))*
     1      (3.0d0-w-4.0d0*xn*w)/(FctrlExt(n)*(alpha**1.5d0)*
     2      ((3.0d0+w)**(xn+1.5d0)))
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision FUNCTION gn(alpha,w,n,tau)
C
C     Purpose: This function computes the Spatial Eigenfunctions gn(tau). See
C              Equation 44 of BW2007.
C
C     Parameters: 
C              alpha  --- Double Precision Input Parameter.
C              w      --- Double Precision Input Parameter.
C              n      --- Integer Input Index
C              tau    --- Double Precision Input Argument.
C
C     Returns:  Double Precision value of gn(tau)
C
C     Routine called: GeneralLaguerre for Laguerre function
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision FUNCTION gn(alpha,w,n,tau)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION alpha,w,tau
      DOUBLE PRECISION CX(0:N)
      CALL GeneralLaguerre(n,-0.5d0,0.5d0*alpha*w*tau**2,cx)
      gn=DEXP(-0.25d0*alpha*(3.d0+w)*tau**2)*cx(n)
      END           ! End of gn
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Subroutine GeneralLaguerre(n,alpha,x,cx)
C
C     Purpose: This subroutine evaluates generalized Laguerre functions 
C              L(N,ALPHA)(X).
C
C     Parameters: 
C              n      --- Integer Input Parameter. N, the highest order 
C                         function to compute.
C              alpha  --- Double Precision Input Parameter. A parameter 
C                         which is part of the definition of the generalized 
C                         Laguerre functions.  ALPHA *MUST* be greater than -1.
C              x      --- Double Precision Input Parameter. X, the point at 
C                         which the functions are to be evaluated. 
C              cx     --- Double Precision Output Array, CX(0:N), the generalized 
C                         Laguerre functions of degrees 0 through N evaluated 
C                         at the point X.
C
C     Remark: Based on Abramowitz and Stegun, Handbook of Mathematical Functions,
C              US Department of Commerce, 1964.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GeneralLaguerre(n,alpha,x,cx)
      implicit none
      integer n
      double precision alpha
      double precision cx(0:n)
      integer i
      double precision x
      if ( alpha <= -1.0D+00 ) then
        write(10,1000)
 1000   format(1X)
        write(10,1001)
 1001   format(1x,"GeneralLaguerre Fatal error.")
        write(10,1002) ALPHA
 1002   format(1x,"Alpha= f12.5 must be > -1")
        stop
      end if
      if ( n < 0 ) then
        return
      end if
      cx(0) = 1.0D+00
      if ( n==0 ) then
        return
      end if
      cx(1) = 1.0D+00+alpha-x
      do 100 i = 2,n
         cx(i)=(((2.d0*i-1.d0)+alpha-x)*cx(i-1) 
     1           +(1.d0-i-alpha)*cx(i-2))/i
100   continue
      return
      end           ! End of GeneralLaguerre
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGM(A,B,X,HG)
C
C     Purpose: Compute confluent hypergeometric function
C              M(a,b,x) (Note: M(a,b,x) <=> 1F1(a;b;x) )
C
C     Parameters: 
C              a  --- Parameter
C              b  --- Parameter ( b <> 0,-1,-2,... )
C              x  --- Argument
C     Output:  HG --- M(a,b,x)
C
C     Routine called: GammaExt for computing Gamma(x)
C
C     Modified by Pete Becker 04/15/2024
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGM(A,B,X,HG)
      IMPLICIT None
      Double Precision A,B,A0,A1,HG,HG1,HG2,PI,X,XG,TA,TB,X0,TBA
      Double Precision SUM1,SUM2,R,RG,R1,R2,Y0,Y1,EPS
      Integer N,NL,LA,I,J,K,M,NMAX
      Double Precision GammaExt
      EPS=1.0d-8
C
C     Initialize variables
C
      PI=3.14159265358979323846264338327950D0
      A0=A
      A1=A
      X0=X
      HG=0.0D0
C
C     Check for special cases
C
      IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
          HG=1.0D+300      ! Whittaker-M is singular!
          RETURN
      ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
          HG=1.0D0         ! A=0 OR A,B!=0 and  X .EQ. 0 case. 
          RETURN
      ELSE IF (A.EQ.-1.0D0) THEN
          HG=1.0D0-X/B     ! A=-1, B=0, X=0 Case.
          RETURN
      ELSE IF (A.EQ.B) THEN
          HG=DEXP(X)       ! A=B, X!=0 Case. Exponential.
          RETURN
      ELSE IF (A-B.EQ.1.0D0) THEN
          HG=(1.0D0+X/B)*DEXP(X) ! A-B=1, X!= 0 Case.
          RETURN
      ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
          HG=(DEXP(X)-1.0D0)/X     ! A=1 and B=2 and X!= 0 Case.
          RETURN
      ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
          M=INT(-A)                ! A=Int < 0, B!=0, X!=0 Case.
          R=1.0D0
          HG=1.0D0
          DO 10 K=1,M
              R=R*(A+Dble(K)-1.0D0)/Dble(K)/(B+Dble(K)-1.0D0)*X
              HG=HG+R
 10      CONTINUE
         RETURN
      ENDIF
      IF (HG.NE.0.0D0) RETURN
C
C     Compute HG using the fundamental power series
C
      TA=1.d0
      HG=1.d0
      NMAX=2000
      Do 100 I=1,NMAX
         TA = TA*(A+DBLE(I)-1.0D0)*X/(DBLE(I)*(B+DBLE(I)-1.0D0))
         HG = HG + TA
         If (DABS(TA/HG).LT.EPS) return
 100  Continue
      WRITE(10,1000)
 1000 format(1x,"Failure to converge in CHGM.")
      STOP
      RETURN
      END           ! End of CHGM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGMOLD(A,B,X,HG)
C
C     Purpose: Compute confluent hypergeometric function
C              M(a,b,x) (Note: M(a,b,x) <=> 1F1(a;b;x) )
C
C     Parameters: 
C              a  --- Parameter
C              b  --- Parameter ( b <> 0,-1,-2,... )
C              x  --- Argument
C     Output:  HG --- M(a,b,x)
C
C     Routine called: GammaExt for computing Gamma(x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGMOLD(A,B,X,HG)
      IMPLICIT None
      Double Precision A,B,A0,A1,HG,HG1,HG2,PI,X,XG,TA,TB,X0,TBA
      Double Precision SUM1,SUM2,R,RG,R1,R2,Y0,Y1
      Integer N,NL,LA,I,J,K,M
      Double Precision GammaExt
      PI=3.14159265358979323846264338327950D0
      A0=A
      A1=A
      X0=X
      HG=0.0D0
      IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
          HG=1.0D+300
      ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
          HG=1.0D0
      ELSE IF (A.EQ.-1.0D0) THEN
          HG=1.0D0-X/B
      ELSE IF (A.EQ.B) THEN
          HG=DEXP(X)
      ELSE IF (A-B.EQ.1.0D0) THEN
          HG=(1.0D0+X/B)*DEXP(X)
      ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
          HG=(DEXP(X)-1.0D0)/X
      ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
          M=INT(-A)
          R=1.0D0
          HG=1.0D0
          DO 10 K=1,M
              R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
              HG=HG+R
 10      CONTINUE
C234567
      ENDIF
      IF (HG.NE.0.0D0) RETURN
      IF (X.LT.0.0D0) THEN
          A=B-A
          A0=A
          X=DABS(X)
      ENDIF
      IF (A.LT.2.0D0) NL=0
      IF (A.GE.2.0D0) THEN
          NL=1
          LA=INT(A)
          A=A-LA-1.0D0
      ENDIF
      DO 30 N=0,NL
          IF (A0.GE.2.0D0) A=A+1.0D0
          IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                  RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                  HG=HG+RG
                  IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
 15           CONTINUE
          ELSE
              TA=GammaExt(A)
              TB=GammaExt(B)
              XG=B-A
              TBA=GammaExt(XG)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8
                  R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                  R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                  SUM1=SUM1+R1
                  SUM2=SUM2+R2
 20           CONTINUE
              HG1=TB/TBA*X**(-A)*DCOS(PI*A)*SUM1
              HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2
              HG=HG1+HG2
          ENDIF
 25       Continue
          IF (N.EQ.0) Y0=HG
          IF (N.EQ.1) Y1=HG
 30   CONTINUE
      IF (A0.GE.2.0D0) THEN
          DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
              A=A+1.0D0
 35       CONTINUE
      ENDIF
      IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
      A=A1
      X=X0
      RETURN
      END           ! End of CHGMOLD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     DOUBLE PRECISION FUNCTION WhitM(kappa,mu,z)
C
C     Purpose: This program computes the Whittaker M function
C
C     Parameters: 
C              kappa  --- Parameter
C              mu     --- Parameter 
C              z      --- Argument
C
C     Returns:  Double Precision value of WhitM(a,b,z)
C
C     Routine called: CHGM for computing 1F1(a;b;z)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION WhitM(kappa,mu,z)
      IMPLICIT None
      DOUBLE PRECISION KAPPA,MU,Z,HG
      CALL CHGM(0.5d0+mu-kappa,1.d0+2.d0*mu,z,HG)
      WhitM=DEXP(-0.5d0*z)*(z**(mu+0.5d0))*HG
      END           ! End of WhitM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     DOUBLE PRECISION FUNCTION WhitW(kappa,mu,z)
C
C     Purpose: This program computes the Whittaker W
C              function
C
C     Parameters: 
C              kappa  --- Parameter
C              mu     --- Parameter 
C              z      --- Argument
C
C     Returns:  Double Precision value of WhitW(a,b,z)
C
C     Routine called: CHGU for computing U(a;b;z)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION WhitW(kappa,mu,z)
      IMPLICIT None
      DOUBLE PRECISION KAPPA,MU,Z,HU
      Integer MD                    ! Method code output by CHGU
      CALL CHGU(0.5d0+mu-kappa,1.d0+2.d0*mu,z,HU,MD)
      WhitW=DEXP(-0.5d0*z)*(z**(mu+0.5d0))*HU
      END           ! End of WhitW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGU(A,B,X,HU,MD)
C
C     Purpose: Compute the confluent hypergeometric function
C              U(a,b,x)
C
C     Parameters: 
C              a  --- Parameter
C              b  --- Parameter
C              x  --- Argument  ( x > 0 )
C              HU --- Output U(a,b,x)
C              MD --- Output Method code
C
C     Routines called:
C          (1) CHGUS for small x ( MD=1 )
C          (2) CHGUL for large x ( MD=2 )
C          (3) CHGUBI for integer b ( MD=3 )
C          (4) CHGUIT for numerical integration ( MD=4 )
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGU(A,B,X,HU,MD)
      IMPLICIT None
      LOGICAL IL1,IL2,IL3,BL1,BL2,BL3,BN
      Double Precision A,B,AA,X,HU,A00,B00,HU1
      Integer ID,ID1,MD
      AA=A-B+1.0D0
      IL1=A.EQ.INT(A).AND.A.LE.0.0D0
      IL2=AA.EQ.INT(AA).AND.AA.LE.0.0D0
      IL3=(DABS(A*(A-B+1.0D0))/X).LE.2.0D0
      BL1=X.LE.5.0D0.OR.(X.LE.10.0D0.AND.A.LE.2.0D0)
      BL2=(X.GT.5.0D0.AND.X.LE.12.5D0).AND.
     1    (A.GE.1.0D0.AND.B.GE.A+4.0D0)
      BL3=X.GT.12.5D0.AND.A.GE.5.0D0.AND.B.GE.A+5.0D0
      BN=B.EQ.INT(B).AND.B.NE.0.0D0
      ID1=-100
      IF (B.NE.INT(B)) THEN
          CALL CHGUS(A,B,X,HU,ID1)
          MD=1
          IF (ID1.GE.6) RETURN
          HU1=HU
      ENDIF
      IF (IL1.OR.IL2.OR.IL3) THEN
          CALL CHGUL(A,B,X,HU,ID)
          MD=2
          IF (ID.GE.6) RETURN
          IF (ID1.GT.ID) THEN
              MD=1
              ID=ID1
              HU=HU1
          ENDIF
      ENDIF
      IF (A.GE.0.0D0) THEN
          IF (BN.AND.(BL1.OR.BL2.OR.BL3)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
          ELSE
              CALL CHGUIT(A,B,X,HU,ID)
              MD=4
          ENDIF
      ELSE
          IF (B.LE.A) THEN
              A00=A
              B00=B
              A=A-B+1.0D0
              B=2.0D0-B
              CALL CHGUIT(A,B,X,HU,ID)
              HU=X**(1.0D0-B00)*HU
              A=A00
              B=B00
              MD=4
          ELSE IF (BN.AND.(.NOT.IL1)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
          ENDIF
      ENDIF
      IF (ID.LT.6) WRITE(10,1000)
 1000 FORMAT(1x,"No accurate result obtained in CHGU")
      RETURN
      END           ! End of CHGU
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGUS(A,B,X,HU,ID)
C
C     Purpose: Compute confluent hypergeometric function
C              U(a,b,x) for small argument x
C
C     Parameters: 
C              a  --- Parameter
C              b  --- Parameter ( b <> 0,-1,-2,...)
C              x  --- Argument
C              HU --- Output U(a,b,x)
C              ID --- Output Estimated number of significant digits
C
C     Routine called: GammaExt for computing gamma function
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGUS(A,B,X,HU,ID)
      IMPLICIT None
      Double Precision A,B,D1,D2,GA,GB,GAB,XG2,GB2,R1,R2,HMAX,HMIN
      double Precision PI,HU0,HUA,H0,HU,X,XG1
      Double Precision GammaExt
      Integer ID,J
      ID=-100
      PI=3.14159265358979323846264338327950D0
      GA=GammaExt(A)
      GB=GammaExt(B)
      XG1=1.0D0+A-B
      GAB=GammaExt(XG1)
      XG2=2.0D0-B
      GB2=GammaExt(XG2)
      HU0=PI/DSIN(PI*B)
      R1=HU0/(GAB*GB)
      R2=HU0*X**(1.0D0-B)/(GA*GB2)
      HU=R1-R2
      HMAX=0.0D0
      HMIN=1.0D+300
      DO 10 J=1,150
          R1=R1*(A+Dble(J)-1.0D0)/(Dble(J)*(B+Dble(J)-1.0D0))*X
          R2=R2*(A-B+Dble(J))/(Dble(J)*(1.0D0-B+Dble(J)))*X
          HU=HU+R1-R2
          HUA=DABS(HU)
          IF (HUA.GT.HMAX) HMAX=HUA
          IF (HUA.LT.HMIN) HMIN=HUA
          IF (DABS(HU-H0).LT.DABS(HU)*1.0D-15) GO TO 15
          H0=HU
 10   CONTINUE
 15   CONTINUE
      D1=DLOG10(HMAX)
      IF (HMIN.NE.0.0D0) D2=DLOG10(HMIN)
      ID=15-DABS(D1-D2)
      RETURN
      END           ! End of CHGUS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGUL(A,B,X,HU,ID)
C
C     Purpose: Compute the confluent hypergeometric function
C              U(a,b,x) for large argument x
C
C     Parameters: 
C              a  --- Parameter
C              b  --- Parameter
C              x  --- Argument
C              HU --- Output U(a,b,x)
C              ID --- Output Estimated number of significant digits
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGUL(A,B,X,HU,ID)
      IMPLICIT None
      Double Precision A,B,AA,HU,R,RA,R0,X
      Integer ID,K,NM
      LOGICAL IL1,IL2
      ID=-100
      AA=A-B+1.0D0
      IL1=A.EQ.INT(A).AND.A.LE.0.0D0
      IL2=AA.EQ.INT(AA).AND.AA.LE.0.0D0
      IF (IL1) NM=DABS(A)
      IF (IL2) NM=DABS(AA)
      IF (IL1.OR.IL2) THEN
          HU=1.0D0
          R=1.0D0
          DO 10 K=1,NM
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              HU=HU+R
 10       CONTINUE
          HU=X**(-A)*HU
          ID=10
      ELSE
          HU=1.0D0
          R=1.0D0
          DO 15 K=1,25
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              RA=DABS(R)
              IF (K.GT.5.AND.RA.GE.R0.OR.RA.LT.1.0D-15) GO TO 20
              R0=RA
              HU=HU+R
 15       CONTINUE
 20       Continue
          ID=DABS(DLOG10(RA))
          HU=X**(-A)*HU
      ENDIF
      RETURN
      END           ! End of CHGUL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGUBI(A,B,X,HU,ID)
C
C     Purpose: Compute confluent hypergeometric function
C              U(a,b,x) with integer b ( b = 1, 2,... )
C
C     Parameters: 
C              a  --- Parameter
C              b  --- Parameter
C              x  --- Argument
C              HU --- Output U(a,b,x)
C              ID --- Output Estimated number of significant digits
C
C     Routines called:
C              (1) GammaExt for computing gamma function Gamma(x)
C              (2) PSI2 for computing psi function
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGUBI(A,B,X,HU,ID)
      IMPLICIT None
      Double Precision GammaExt
      Integer ID,ID1,ID2,J,K,M,N
      Double Precision A,A0,A1,A2,B,DA1,DA2,DB1,DB2,EL,GA,GA1
      Double Precision H0,HM1,HM2,HM3,HU,HU1,HU2,HMIN,HMAX,HW,PS
      Double Precision R,RN,RN1,S0,S1,S2,SA,SB,UA,UB,X
      ID=-100
      EL=0.5772156649015329D0
      N=ABS(B-1)
      RN1=1.0D0
      RN=1.0D0
      DO 10 J=1,N
          RN=RN*J
          IF (J.EQ.N-1) RN1=RN
 10   CONTINUE
      CALL PSI2(A,PS)
      GA=GammaExt(A)
      IF (B.GT.0.0D0) THEN
          A0=A
          A1=A-N
          A2=A1
          GA1=GammaExt(A1)
          UA=(-1)**(N-1)/(RN*GA1)
          UB=RN1/GA*X**(-N)
      ELSE
          A0=A+N
          A1=A0
          A2=A
          GA1=GammaExt(A1)
          UA=(-1)**(N-1)/(RN*GA)*X**N
          UB=RN1/GA1
      ENDIF
      HM1=1.0D0
      R=1.0D0
      HMAX=0.0D0
      HMIN=1.0D+300
      DO 15 K=1,150
          R=R*(A0+K-1.0D0)*X/((N+K)*K)
          HM1=HM1+R
          HU1=DABS(HM1)
          IF (HU1.GT.HMAX) HMAX=HU1
          IF (HU1.LT.HMIN) HMIN=HU1
          IF (DABS(HM1-H0).LT.DABS(HM1)*1.0D-15) GO TO 20
          H0=HM1
 15   CONTINUE
 20   CONTINUE
      DA1=DLOG10(HMAX)
      IF (HMIN.NE.0.0D0) DA2=DLOG10(HMIN)
      ID=15-DABS(DA1-DA2)
      HM1=HM1*DLOG(X)
      S0=0.0D0
      DO 25 M=1,N
          IF (B.GE.0.0D0) S0=S0-1.0D0/M
          IF (B.LT.0.0D0) S0=S0+(1.0D0-A)/(M*(A+M-1.0D0))
 25   CONTINUE
      HM2=PS+2.0D0*EL+S0
      R=1.0D0
      HMAX=0.0D0
      HMIN=1.0D+300
      DO 50 K=1,150
          S1=0.0D0
          S2=0.0D0
          IF (B.GT.0.0) THEN
              DO 30 M=1,K
                  S1=S1-(M+2.0D0*A-2.0D0)/(M*(M+A-1.0D0))
 30           CONTINUE
              DO 35 M=1,N
                  S2=S2+1.0D0/(K+M)
 35           CONTINUE
          ELSE
              DO 40 M=1,K+N
                  S1=S1+(1.0D0-A)/(M*(M+A-1.0D0))
 40           CONTINUE
              DO 45 M=1,K
                  S2=S2+1.0D0/M
 45           CONTINUE
          ENDIF
          HW=2.0D0*EL+PS+S1-S2
          R=R*(A0+K-1.0D0)*X/((N+K)*K)
          HM2=HM2+R*HW
          HU2=DABS(HM2)
          IF (HU2.GT.HMAX) HMAX=HU2
          IF (HU2.LT.HMIN) HMIN=HU2
          IF (DABS((HM2-H0)/HM2).LT.1.0D-15) GO TO 55
          H0=HM2
 50   CONTINUE
 55   CONTINUE
      DB1=DLOG10(HMAX)
      IF (HMIN.NE.0.0D0) DB2=DLOG10(HMIN)
      ID1=15-DABS(DB1-DB2)
      IF (ID1.LT.ID) ID=ID1
      HM3=1.0D0
      IF (N.EQ.0) HM3=0.0D0
      R=1.0D0
      DO 60 K=1,N-1
          R=R*(A2+K-1.0D0)/((K-N)*K)*X
          HM3=HM3+R
 60   CONTINUE
      SA=UA*(HM1+HM2)
      SB=UB*HM3
      HU=SA+SB
      IF (SA.NE.0.0D0) ID1=INT(DLOG10(DABS(SA)))
      IF (HU.NE.0.0D0) ID2=INT(DLOG10(DABS(HU)))
      IF (SA*SB.LT.0.0D0) ID=ID-ABS(ID1-ID2)
      RETURN
      END           ! End of CHGUBI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE CHGUIT(A,B,X,HU,ID)
C
C     Purpose: Compute hypergeometric function U(a,b,x) by
C              using Gaussian-Legendre integration (n=60)
C
C     Parameters: 
C              a  --- Parameter ( a > 0 )
C              b  --- Parameter
C              x  --- Argument ( x > 0 )
C              HU --- Output U(a,b,z)
C              ID --- Output Estimated number of significant digits
C
C     Routine called: GammaExt for computing Gamma(x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHGUIT(A,B,X,HU,ID)
      IMPLICIT None
      Double Precision GammaExt
      Double Precision A,A1,B,B1,C,D,F1,F2,G,GA,HU,HU0,HU1,HU2
      Double Precision T(30),W(30)
      Double Precision S,T1,T2,T3,T4,X
      Integer ID,J,K,M
      DATA T/ .259597723012478D-01, .778093339495366D-01,
     &        .129449135396945D+00, .180739964873425D+00,
     &        .231543551376029D+00, .281722937423262D+00,
     &        .331142848268448D+00, .379670056576798D+00,
     &        .427173741583078D+00, .473525841761707D+00,
     &        .518601400058570D+00, .562278900753945D+00,
     &        .604440597048510D+00, .644972828489477D+00,
     &        .683766327381356D+00, .720716513355730D+00,
     &        .755723775306586D+00, .788693739932264D+00,
     &        .819537526162146D+00, .848171984785930D+00,
     &        .874519922646898D+00, .898510310810046D+00,
     &        .920078476177628D+00, .939166276116423D+00,
     &        .955722255839996D+00, .969701788765053D+00,
     &        .981067201752598D+00, .989787895222222D+00,
     &        .995840525118838D+00, .999210123227436D+00/
      DATA W/ .519078776312206D-01, .517679431749102D-01,
     &        .514884515009810D-01, .510701560698557D-01,
     &        .505141845325094D-01, .498220356905502D-01,
     &        .489955754557568D-01, .480370318199712D-01,
     &        .469489888489122D-01, .457343797161145D-01,
     &        .443964787957872D-01, .429388928359356D-01,
     &        .413655512355848D-01, .396806954523808D-01,
     &        .378888675692434D-01, .359948980510845D-01,
     &        .340038927249464D-01, .319212190192963D-01,
     &        .297524915007890D-01, .275035567499248D-01,
     &        .251804776215213D-01, .227895169439978D-01,
     &        .203371207294572D-01, .178299010142074D-01,
     &        .152746185967848D-01, .126781664768159D-01,
     &        .100475571822880D-01, .738993116334531D-02,
     &        .471272992695363D-02, .202681196887362D-02/
      ID=7
      A1=A-1.0D0
      B1=B-A-1.0D0
      C=12.0D0/X
      DO 20 M=10,100,5
          HU1=0.0D0
          G=0.5D0*C/M
          D=G
          DO 15 J=1,M
              S=0.0D0
              DO 10 K=1,30
                  T1=D+G*T(K)
                  T2=D-G*T(K)
                  F1=DEXP(-X*T1)*T1**A1*(1.0D0+T1)**B1
                  F2=DEXP(-X*T2)*T2**A1*(1.0D0+T2)**B1
                  S=S+W(K)*(F1+F2)
 10           CONTINUE
              HU1=HU1+S*G
              D=D+2.0D0*G
 15       CONTINUE
          IF (DABS(1.0D0-HU0/HU1).LT.1.0D-7) GO TO 25
          HU0=HU1
 20   CONTINUE
 25   CONTINUE
      GA=GammaExt(A)
      HU1=HU1/GA
      DO 40 M=2,10,2
          HU2=0.0D0
          G=0.5D0/M
          D=G
          DO 35 J=1,M
              S=0.0D0
              DO 30 K=1,30
                  T1=D+G*T(K)
                  T2=D-G*T(K)
                  T3=C/(1.0D0-T1)
                  T4=C/(1.0D0-T2)
                  F1=T3*T3/C*DEXP(-X*T3)*T3**A1*(1.0D0+T3)**B1
                  F2=T4*T4/C*DEXP(-X*T4)*T4**A1*(1.0D0+T4)**B1
                  S=S+W(K)*(F1+F2)
 30           CONTINUE
              HU2=HU2+S*G
              D=D+2.0D0*G
 35       CONTINUE
          IF (DABS(1.0D0-HU0/HU2).LT.1.0D-7) GO TO 45
          HU0=HU2
 40   CONTINUE
 45   CONTINUE
      GA=GammaExt(A)
      HU2=HU2/GA
      HU=HU1+HU2
      RETURN
      END           ! End of CHGUIT
C234567
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE PSI2(X,PS)
C
C     Purpose: Compute Psi function
C
C     Parameters: 
C              x  --- Argument of psi(x)
C              PS --- Output value of psi(x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PSI2(X,PS)
      IMPLICIT None
      Double Precision A1,A2,A3,A4,A5,A6,A7,A8
      Double Precision EL,PI,PS,S,X,XA,X2
      Integer K,N
      XA=DABS(X)
      PI=3.14159265358979323846264338327950D0
      EL=.5772156649015329D0
      S=0.0D0
      IF (X.EQ.INT(X).AND.X.LE.0.0D0) THEN
          PS=1.0D+300
          RETURN
      ELSE IF (XA.EQ.INT(XA)) THEN
          N=Int(XA)
          DO 10 K=1 ,N-1
              S=S+1.0D0/Dble(K)
 10       CONTINUE
          PS=-EL+S
      ELSE IF (XA+.5D0.EQ.INT(XA+.5D0)) THEN
          N=Int(XA-.5D0)
          DO 20 K=1,N
              S=S+1.0D0/(2.0D0*Dble(K)-1.0D0)
 20       CONTINUE
          PS=-EL+2.0D0*S-1.386294361119891D0
      ELSE
          IF (XA.LT.10.0D0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
                  S=S+1.0D0/(XA+Dble(K))
 30           CONTINUE
              XA=XA+Dble(N)
          ENDIF
          X2=1.0D0/(XA*XA)
          A1=-.8333333333333D-01
          A2=.83333333333333333D-02
          A3=-.39682539682539683D-02
          A4=.41666666666666667D-02
          A5=-.75757575757575758D-02
          A6=.21092796092796093D-01
          A7=-.83333333333333333D-01
          A8=.4432598039215686D0
          PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
          PS=PS-S
      ENDIF
      IF (X.LT.0.0D0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
      RETURN
      END           ! End of PSI2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE HG2F2(a1,a2,b1,b2,z,HGpFq)
C
C     Purpose: Compute the Hypergeometric function 2F2
C
C     Parameters: 
C              a1,a2  --- Double Precision Input Numerator Parameters
C              b1,b2  --- Double Precision Input Denominator Parameters 
C                         ( b <> 0,-1,-2,... )
C              z      --- Double Precision Input Argument
C              HGpFq  --- Ouptut Double Precision Result for 2F2(a1,a2;b1,b2;z)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HG2F2(a1,a2,b1,b2,z,HGpFq)
      Implicit None
      Double Precision a1,a2,b1,b2,z,HGpFq
      Double Precision eps,T
      Integer i,nmax
      eps=1.0d-10
      T=1.d0
      HGpFq=1.d0
      nmax=2000
      Do 100 i=1,nmax
         T = T*(a1+dble(i)-1.0d0)*(a2+dble(i)-1.0d0)*z
     1       /(dble(i)*(b1+dble(i)-1.0d0)*(b2+dble(i)-1.0d0))
         HGpFq = HGpFq + T
         If ( DABS(T/HGpFq).lt.eps ) return
 100  Continue
      WRITE(10,1000)
 1000 FORMAT(1X,"Failure to converge in HG2F2.")
      STOP
      END           ! End of HG2F2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function lngamz(xx)
C
C     Purpose: Computes the value of log[Gamma(xx)] (natural logarightm of the 
C     Gamma function) for real xx>0. This function uses the numerical algorithm 
C     worked out by Lanczos (1964, J. SIAM Numer. Anal. Ser. B, v.1, p.86).
C     See Lanczos (1964) equation (32) and page 95. This routine is also 
C     DoD-sourced to Sommer (1981).
C
C     Parameters:
C              xx      --- Double Precision Input: Argument.
C
C     Returns: Double Precision value of lngamz(xx)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function lngamz(xx)
      implicit none
      double precision x,xx,a5z,gam
      gam=5.0d0
      x=xx
      a5z=1.0000000001855653d0
     1    +76.18009172943101d0/(x+1.0d0)
     2    -86.5053203280994d0/(x+2.0d0)
     3    +24.01409823080286d0/(x+3.0d0)
     4    -1.231739543443382d0/(x+4.0d0)
     5    +0.1208615717587236d-2/(x+5.0d0)
     6    -0.5380068614613265d-5/(x+6.0d0)
C      print 1000,x,gam,a5z
C 1000 format(1x,3(1pd22.14,1x))
      lngamz=(x+0.5d0)*Dlog(x+gam+0.5d0)-(x+gam+0.5d0)
     1       +Dlog(2.506628274631d0*a5z/x)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function NUpInGamFunc(AIN,XIN)
C     Purpose: Compute and return the Normalized Upper Incomplete Gamma 
C              Function Q(a,x). For the Incomplete Gamma Function, there are 
C              really two functions: A Lower (lower case gamma) and Upper 
C              (upper case gamma) Incomplete Gamma Function. The sum of the 
C              lower and upper incomplete gamma functions is the full gamma 
C              function. This allows one to solve for either the upper or 
C              lower incomplete gamma function by simply first obtaining 
C              the full gamma function and then one or the other incomplete 
C              gamma function.
C
C     Parameters:
C              ain      --- Double Precision Input: Parameter.
C              xin      --- Double Precision Input: Argument.
C
C     Returns: Double Precision Normalized upper gamma function
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function NUpInGamFunc(AIN,XIN)
      Implicit None
      Double Precision X,XIN,TINYNUM,Dj,Djm1,Cj,Cjm1
      Double Precision Deltaj,Fj,Fjm1,Toler
      Double Precision a,b,c,d,f,g,lngam,lngamz
      Double Precision Euler,term1,term2,Delta,FctrlExt
      Double Precision term3,gam,ain,tmp,atmp
      Integer i,j
      Data Euler /0.577215664901532860D0/
      Data Toler /1.0D-9/
      Data TINYNUM /1.0D-34/
      X=XIN
      A=AIN
      lngam=lngamz(a)
      gam=Dexp(lngam)
      term1=Dexp(-x)
      term2=x**a
      if (x.lt.a+1) then
C
C     If x<a+1 then I use the series expression for the absolute
C     lower incomplete gamma function (lower case gamma) and then 
C     take its complement difference with 1.0. 
C
      atmp=a
      delta=1.0d0/atmp
      term3=delta
      do 100 i=1,100
          atmp=atmp+1.0d0
          delta=delta*x/atmp
          term3=term3+delta
          if (Dabs(delta).lt.Dabs(term3)*Toler) then
              NUpInGamFunc=1.0-term3*Dexp(-x+a*Dlog(x)-lngam)
              Return
          Endif
 100  Continue
      NUpInGamFunc=1.0-term3*Dexp(-x+a*Dlog(x)-lngam)
      Return
      else
C
C     If x>a+1 then  I use the continued fraction calculation for the 
C     absolute upper incomplete gamma function. 
C
      G=ain
      F=0.0d0+TINYNUM
      C=F
      D=0.0d0
      Djm1=D
      Cjm1=C
      Fjm1=F
      do 200 j=1,100
C
C     First, set up the values of a and b:
C
          if (mod(j,2).eq.0) then
              a=Dble(j/2)-g
              b=1.0d0
          else
              a=Dble(j/2)
              if (j.eq.1) a=1.0d0
              if (j.eq.3) a=1.0d0
              b=x
          endif
C
C     Now, calcuate the C and D parameters:
C
          Dj=b+a*Djm1
          If(Dj.eq.0.0d0) Dj=TINYNUM
          Cj=b+a/Cjm1
          If(Cj.eq.0.0d0) Cj=TINYNUM
          Dj=1.0d0/Dj
          Deltaj=Cj*Dj
          Fj=Fjm1*Deltaj
          If (Dabs(Deltaj-1.0d0).lt.Toler) then
              NUpInGamFunc=Dexp(-x+g*Dlog(x)-lngam)*Fj
              return
          Endif
          Cjm1=Cj
          Djm1=Dj
          Fjm1=Fj
 200  Continue
      NUpInGamFunc=Dexp(-x+g*Dlog(x)-lngam)*Fj
      return
      Endif
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function NLoInGamFunc(AIN,XIN)
C     Purpose: Compute and return the Normalized Lower Incomplete Gamma 
C              Function P(a,x). For the Incomplete Gamma Function, there are 
C              really two functions: A Lower (lower case gamma) and Upper 
C              (upper case gamma) Incomplete Gamma Function. The sum of the 
C              lower and upper incomplete gamma functions is the full gamma 
C              function. This allows one to solve for either the upper or 
C              lower incomplete gamma function by simply first obtaining 
C              the full gamma function and then one or the other incomplete 
C              gamma function.
C
C     Parameters:
C              ain      --- Double Precision Input: Parameter.
C              xin      --- Double Precision Input: Argument.
C
C     Returns: Double Precision Normalized lower gamma function
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function NLoInGamFunc(AIN,XIN)
      Implicit None
      Double Precision X,XIN,TINYNUM,Dj,Djm1,Cj,Cjm1
      Double Precision Deltaj,Fj,Fjm1,Toler
      Double Precision a,b,c,d,f,g,lngam,lngamz
      Double Precision Euler,term1,term2,Delta,FctrlExt
      Double Precision term3,gam,ain,tmp,atmp
      Integer i,j
      Data Euler /0.577215664901532860D0/
      Data Toler /1.0D-9/
      Data TINYNUM /1.0D-34/
      X=XIN
      A=AIN
      lngam=lngamz(a)
      gam=Dexp(lngam)
      term1=Dexp(-x)
      term2=x**a
      if (x.lt.a+1) then
C
C     If x<a+1 then I use the series expression for the absolute
C     lower incomplete gamma function (lower case gamma). 
C
      atmp=a
      delta=1.0d0/atmp
      term3=delta
      do 100 i=1,100
          atmp=atmp+1.0d0
          delta=delta*x/atmp
          term3=term3+delta
          if (Dabs(delta).lt.Dabs(term3)*Toler) then
              NLoInGamFunc=term3*Dexp(-x+a*Dlog(x)-lngam)
              Return
          Endif
 100  Continue
      NLoInGamFunc=term3*Dexp(-x+a*Dlog(x)-lngam)
      Return
      else
C
C     If x>a+1 then  I use the continued fraction calculation for the 
C     absolute upper incomplete gamma function taking the difference with 1.0.
C
      G=ain
      F=0.0d0+TINYNUM
      C=F
      D=0.0d0
      Djm1=D
      Cjm1=C
      Fjm1=F
      do 200 j=1,100
C
C     First, set up the values of a and b:
C
          if (mod(j,2).eq.0) then
              a=Dble(j/2)-g
              b=1.0d0
          else
              a=Dble(j/2)
              if (j.eq.1) a=1.0d0
              if (j.eq.3) a=1.0d0
              b=x
          endif
C
C     Now, calcuate the C and D parameters:
C
          Dj=b+a*Djm1
          If(Dj.eq.0.0d0) Dj=TINYNUM
          Cj=b+a/Cjm1
          If(Cj.eq.0.0d0) Cj=TINYNUM
          Dj=1.0d0/Dj
          Deltaj=Cj*Dj
          Fj=Fjm1*Deltaj
          If (Dabs(Deltaj-1.0d0).lt.Toler) then
              NLoInGamFunc=1.0d0-Dexp(-x+g*Dlog(x)-lngam)*Fj
              return
          Endif
          Cjm1=Cj
          Djm1=Dj
          Fjm1=Fj
 200  Continue
      NLoInGamFunc=1.0d0-Dexp(-x+g*Dlog(x)-lngam)*Fj
      return
      Endif
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Double Precision Function Exp1Int(xin)
C
C     Purpose: This function evaluates the 1st Exponential Integral based on
C              either a series or continued fraction expression. 
C
C     Parameters:
C              xin    --- Double Precision Input: Parameter.
C
C     Returns: Double Precision value for the 1st Exponential Integral for the
C              input value.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Double Precision Function Exp1Int(XIN)
      Implicit None
      Double Precision X,XIN,TINYNUM,Dj,Djm1,Cj,Cjm1
      Double Precision Deltaj,Fj,Fjm1,Toler
      Double Precision a,b,c,d,f
      Double Precision Euler,term1,term2,Delta,FctrlExt
      Integer i,j
      Data Euler /0.577215664901532860D0/
      Data Toler /1.0D-9/
      Data TINYNUM /1.0D-35/
      X=XIN
      If (X.LT.0.0D0) Then
          Write(6,2000) 
 2000     Format(1x," Error input to Exp1Int: Xin .lt. 0 ")
      Endif
C
C     Check to see if I use series expansion or continued fractions.
C
      If (X.lt.1.0d0) then   ! Use series expansion
          term1=-Dlog(x)
          term2=0.0d0
          do 200 i=1,50
              Delta=((-x)**i)/(Dble(i)*FctrlExt(i))
              term2=term2+Delta
              if (Dabs(Delta).lt.Toler) then
                  Exp1Int=-term2+term1-Euler
                  Return
              Endif
 200      Continue
          Exp1Int=term2+term1-Euler
          Return
      Endif
C
C     Crank up the continued fraction calculation.
C
      F=0.0d0+TINYNUM
      C=F
      D=0.0d0
      Djm1=D
      Cjm1=C
      Fjm1=F
      do 100 j=1,600
C
C     First, set up the values of a and b:
C
          if (j.eq.1) then
              a=1.0d0
              b=x
          elseif (j.eq.2) then
              a=1.0d0
              b=1.0d0
          elseif (j.eq.3) then
              a=1.0d0
              b=x
          elseif (j.gt.3) then
              if ( mod(j,2).eq.0 ) then
                  b=1.0d0
                  a=Dble(j/2)
              else
                  a=Dble((j-1)/2)
                  b=x
              endif
          endif
C
C     Now, calcuate the C and D parameters:
C
          Dj=b+a*Djm1
          If(Dj.le.0.0d0) Dj=TINYNUM
          Cj=b+a/Cjm1
          If(Cj.le.0.0d0) Cj=TINYNUM
          Dj=1.0d0/Dj
          Deltaj=Cj*Dj
          Fj=Fjm1*Deltaj
          If (Dabs(Deltaj-1.0d0).lt.Toler) then
              Exp1Int=Dexp(-x)*Fj
              return
          Endif
          Cjm1=Cj
          Djm1=Dj
          Fjm1=Fj
 100  Continue
      Exp1Int=Dexp(-x)*Fj
      return
      end

C
C     This function computes the unity-normalized column-integrated
C     Green's function from the BW07 paper.
C
      Double Precision Function BWGreenCol(Te,alpha,xi,
     1                          delta,tau0,eps0,eps,nmax)
      Implicit None
      Double precision Te,alpha,xi,delta,tau0,eps0,eps
      Double Precision w,kappa,mu,lambda,WhitM,WhitW,WhitProd
      Double Precision boltzc,GammaExt,FctrlExt,XFunc,gn
      Double Precision Ndot,NormVal,chiMin,chiMax,term,sum
      Double Precision gamfdum1,gamfdum2,gamfdum3
      Integer nmax,nloop
      boltzc=1.3807d-16
      kappa=0.5d0*(delta+4.0d0)
      w=dsqrt(9.0d0+12.0d0*xi*xi)
         Ndot=1.0d0
         chiMin=Dmin1(eps,eps0)/(boltzc*Te)
         chiMax=Dmax1(eps,eps0)/(boltzc*Te)
         NormVal=3*Ndot*delta*xi*xi*(boltzc*Te)*
     1        DSqrt(2.0d0*alpha*alpha*alpha*w)*
     2        DExp((3.0d0*alpha*tau0*tau0/2.0d0)+
     3        ((eps0-eps)/(2.0d0*boltzc*Te)))*(eps**(kappa-2.0d0))*
     4        (eps0**(-kappa))
      sum=1.0d-100
C
C     Loop to compute column-integrated Green's function, written by
C     PAB on 2/19/2025.
C
      do 120 nloop=0,nmax
        lambda=(4.0d0*dble(nloop)*w+w+3.0d0)/2.0d0
        mu=dsqrt((3.0d0-delta)*(3.0d0-delta)+4.0d0*delta*lambda)/2.0d0
        gamfdum1=GammaExt(mu-kappa+0.5d0)
        gamfdum2=GammaExt(2.0d0*mu+1.0d0)
        gamfdum3=GammaExt(nloop+0.5d0)
        WhitProd=WhitM(kappa,mu,chiMin)*WhitW(kappa,mu,chiMax)
        term=gamfdum1*FctrlExt(nloop)*
     1       XFunc(nloop,w,alpha)/(gamfdum2*
     2       gamfdum3)*gn(alpha,w,nloop,tau0)*WhitProd
        sum=sum+term
        if ( nloop.gt.3 .and. DAbs(term).lt.1.0d-3*DAbs(sum) )
     1       goto 121

 120    continue
 121  continue
      BWGreenCol=sum*NormVal
      return
      end
C
C     This function computes the integrand for the BB integration
C     in the BW07 paper, written by PAB on 2/19/2025.
C
      Double Precision Function BBintegrand(Tth,Te,alpha,xi,
     1                          delta,tauTH,eps0,eps,caprad,nmax)
      Implicit None
      Double precision Tth,Te,alpha,xi,delta,tauTH,eps0,eps,caprad
      Double Precision PI,boltzc,planckc,planckc3,clight
      Double Precision PlanckSource,BWGreenCol
      Integer nmax
      PI=3.14159265358979323846264338327950d0
      boltzc=1.3807d-16
      planckc=6.6260755d-27
      planckc3=planckc*planckc*planckc
      clight=2.99792458d10
      PlanckSource=2*PI*PI*caprad*caprad/(planckc3*clight*clight)/
     1        (DExp(eps0/(boltzc*Tth))-1.0d0)
      BBintegrand=eps0*eps0*PlanckSource*BWGreenCol(Te,alpha,xi,
     1                          delta,tauTH,eps0,eps,nmax)
      return
      end
C
C     This function computes the integral for the BB integration
C     in the BW07 paper using a trapezoid algorithm. The output is
C     the column-integrated photon number spectrum, written by PAB
C     on 2/19/2025.
C
      Double Precision Function BBtrap(Tth,Te,alpha,xi,
     1                          delta,tauTH,eps,caprad,nmax)
      Implicit None
      Double precision Tth,Te,alpha,xi,delta,tauTH,eps,caprad
      Double Precision PI,boltzc,planckc,planckc3,clight
      Double Precision eps0Min,eps0Max,DeltaEps0
      Double Precision BBintegrand,term,sum
      Double Precision eps0Left,eps0Right,funLeft,funRight
      Integer numVals,nmax,nloop
      PI=3.14159265358979323846264338327950d0
      boltzc=1.3807d-16
      planckc=6.6260755d-27
      planckc3=planckc*planckc*planckc
      clight=2.99792458d10
      numVals=50
      eps0Min=0.000001*eps
      eps0Max=Dmin1(7.0d0*eps,25.0d0*boltzc*Te)
      DeltaEps0=(eps0Max-eps0Min)/(numVals-1)
C
C     Loop to compute integral, written by PAB on 2/19/2025.
C
      sum=0.0d0
      eps0Left=eps0Min
      funLeft=BBintegrand(Tth,Te,alpha,xi,
     1                    delta,tauTH,eps0Left,eps,caprad,nmax)
      do 120 nloop=1,numVals-1
        eps0Right=eps0Min+nloop*DeltaEps0
        funRight=BBintegrand(Tth,Te,alpha,xi,
     1                       delta,tauTH,eps0Right,eps,caprad,nmax)
        term=0.5d0*(funLeft+funRight)*DeltaEps0
        sum=sum+term
        eps0Left=eps0Right
        funLeft=funRight
 120  continue
      BBtrap=sum
      return
      end


