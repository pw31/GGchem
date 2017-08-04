*********************************************************************
      SUBROUTINE NUCLEATION(species,T,V0,mass,n1in,SSin,Jstar,Nstar)
*********************************************************************
*****                                                           *****
*****  computes nucleation rate according to                    *****
*****  classical nucleation theory (Gail et al 1984)            *****
*****                                                           *****
*****  INPUT:  species = name of nucleating species             *****
*****          T  = Gastemperatur [K]                           *****
*****          V0 = monomer volume [cm-3]                       *****
*****          mass = mass of monomer [amu]                     *****
*****          n1 = monomer particle density [cm-3]             *****
*****          SS = supersaturation ratio                       ***** 
*****                                                           *****
*****  OUTPUT: Jstar = nucleation rate [cm^-3 s^-1]             *****
*****          Nstar = critical cluster size                    *****
*****                                                           *****
*********************************************************************
      use DUST_DATA,only: qp
      implicit none
      character(len=*),intent(in) :: species
      real*8,intent(in) :: T,V0,mass
      real(kind=qp),intent(in) :: n1in,SSin
      real*8,intent(out) :: Jstar,Nstar
      real*8 :: pi,bk,amu,a0,f0,sigma,Nf,alpha,slog
      real*8 :: thetun,thetaN,x0,x1,x2,x3,dgdn,nst,n1,SS
      real*8 :: zeldof,vth,beta,fNst
      data pi/3.14159265358979D+0/, bk/1.38066D-16/, amu/1.66055D-24/

      n1 = n1in
      SS = SSin
      if (trim(species)=='C') then
        Nf    = 5.0                        ! fit for sigma=sigma(N)
        sigma = 1400.                      ! erg/cm2 (Gail+1984)
      else if (trim(species)=='W') then
        Nf    = 10.0 
        sigma = 3340.                      ! erg/cm2 (R.Tran+2016)
      else
        print*,"*** surface tension not known in NUCLEATION.f"
        stop
      endif  

*     -----------------------------------------
*     ***  monomer radius and surface area  *** 
*     -----------------------------------------
      a0 = (3.d0*V0/(4.d0*pi))**(1.d0/3.d0)
      f0 = 4.d0*pi*a0**2
      !print*,T,V0,a0,f0,mass/amu,n1,SS

*     -------------------------------
*     ***  supersaturation ratio  ***
*     -------------------------------
      if (SS.le.1.d0) then
        Jstar = 0.d+0
        Nstar = 9.d+99
        goto 500
      end if  
      slog = LOG(SS)

*     -------------------------------------------------------------
*     ***  size of critical cluster according to droplet model  ***
*     -------------------------------------------------------------
      thetun = f0*sigma/bk
      x0     = 2.d0*thetun/(3.d0*T*slog)
      x1     = Nf**(1.d0/3.d0)
      x2     = x1/x0
      x3     = 0.5d0*(1.d0+DSQRT(1.d0+2.d0*x2)) - x2
      Nstar  = 1.d0 + (x0*x3)**3 
      if (Nstar<=1.d0) Nstar=1.000000001d0
*
*     --------------------------------------------
*     ***  number density of critical cluster  ***
*     --------------------------------------------
      x0     = x0*slog
      x2     = (Nstar-1.d0)**(1.d0/3.d0) + x1
      x3     = 1.d0/(Nstar-1.d0)**(2.d0/3.d0)
      dgdn   = x0*x3* ( 1.d0/x2**2 + x1/x2**3 ) / 3.d0
      zeldof = SQRT(dgdn/(2.d0*pi))
      thetaN = thetun/(1.d0+(Nf/(Nstar-1.d0))**(1.d0/3.d0))
      x1     = (Nstar-1.d0)*slog - (thetaN/T)
     &         *(Nstar-1.d0)**(2.d0/3.d0)
      nst    = n1*EXP(x1)
*
*     -------------------------
*     ***  growth velocity  ***
*     -------------------------
      alpha = 1.d0
      vth  = SQRT(bk*T/(2.d0*pi*mass))
      beta = alpha*vth*n1
      fNst = f0*Nstar**(2.d0/3.d0)
*
*     -------------------------
*     ***  nucleation rate  ***
*     -------------------------
      Jstar = nst*beta*fNst*zeldof

      !print*,T,SS,Nstar,n1,nst/n1,Jstar,vth
      !if (Nstar<20) stop

 500  continue
      RETURN
      end
