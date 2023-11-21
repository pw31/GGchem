************************************************************************
      module DISK
************************************************************************
      use OPACITY,ONLY: NLAMmax
      integer,parameter :: NZmax=1000
      real :: Mstar,Rstar,Mdot
      real,dimension(NZmax) :: zz,Td,nH,JJ,DD,kRoss,tauRoss,Hvis,dustgas
      real,dimension(NLAMmax,NZmax) :: kex
      integer :: Nx,Nz,izStart
      end
      
************************************************************************
      SUBROUTINE INIT_DISK
************************************************************************
      use DISK,ONLY: Nx,Nz,izStart,Mstar,Rstar,Mdot,
     >               zz,nH,Td,kRoss,tauRoss,Hvis,dustgas
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: grav=6.67259850000E-08
      real,parameter :: AU=1.495978700E+13
      real,parameter :: Msun=1.988922500E+33
      real,parameter :: Rsun=6.959900000E+10
      real,parameter :: yr=3.155760000E+7 
      integer :: i,ix,iz,ixSearch,Nsp,Nheat,Ncool,Nlam
      character(len=99999) :: line
      real :: rr,arr(1000),z1,z2,z3,f1,f2,f3,dz,Hcol,Hint,const,cc

      !----------------------------------------------
      ! ***  read a vertical column from ProDiMo  ***
      !----------------------------------------------
      print*
      print*,"reading ProDiMo.out ..."
      ixSearch = 35
      open(unit=12,file="ProDiMo.out",status='old')
      do i=1,21
        read(12,'(A99999)') line
        if (index(line,'Mstar')) read(line(20:),*) Mstar
        if (index(line,'Rstar')) read(line(20:),*) Rstar
      enddo
      Mstar = Mstar*Msun
      Rstar = Rstar*Rsun
      read(12,'(A99999)') line
      read(line(25:),*) Nx,Nz,Nsp,Nheat,Ncool,Nlam
      !print*,Nx,Nz,Nsp,Nheat,Ncool
      read(12,*)
      read(12,*)
      print'(A4,A4,2A11,A9,99A11)','ix','iz','r[AU]','z[AU]',
     >     'Td[K]','n<H>[cm-3]','kRoss','tauRoss','d/g'
      do iz=Nz,1,-1
        do ix=1,Nx
          read(12,'(A99999)') line
          read(line,*,end=100) arr
 100      if (ix==ixSearch) then
            rr     = arr(3)*AU
            zz(iz) = arr(4)*AU
            Td(iz) = arr(11)
            nH(iz) = arr(22+Nheat+Ncool)
            kRoss(iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+2)
            tauRoss(iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+3)
            dustgas(iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+5)
            if (tauRoss(iz)>10.0) then
              izStart=MAX(iz,izStart)
            endif
            print'(I4,I4,2(0pF11.6),0pF9.2,99(1pE11.4))',ix,iz,
     >           rr/AU,zz(iz)/AU,Td(iz),nH(iz),kRoss(iz),tauRoss(iz),
     >           dustgas(iz)     
          endif
        enddo
      enddo

      !---------------------------------------------
      ! ***  compute viscous heating per volume  ***
      !---------------------------------------------
      Mdot = 1.E-8*Msun/yr
      Hcol = grav*Mstar*Mdot/rr**3 * 3.0/(8.0*pi)   ! [erg/cm2/s]
     >     * (1.0-SQRT(Rstar/rr))
      Hint = 0.0
      do iz=1,Nz-1
        z1 = zz(iz)                                 ! [cm]
        z3 = zz(iz+1)
        z2 = 0.5*(z1+z3)
        f1 = nH(iz)*dustgas(iz)                     ! [cm-3]
        f3 = nH(iz+1)*dustgas(iz+1)
        !------ fit rho(z) = rho(z1)*exp(-(z-z1)^2/H^2) --------
        cc = LOG(f3/f1)/(z3-z1)**2                  ! -1/H^2
        f2 = f1*EXP(cc*(z2-z1)**2)   
        Hint = Hint + (f1+4*f2+f3)*(z3-z1)/6.0      ! [cm-2]
      enddo  
      const = Hcol/Hint                             ! [erg/s]
      do iz=1,Nz
        Hvis(iz) = const*nH(iz)*dustgas(iz)         ! [erg/cm3/s]
      enddo
      end

************************************************************************
      SUBROUTINE COMPUTE_DISK
************************************************************************
      use PARAMETERS,ONLY: verbose
      use DISK,ONLY: Nz,izStart,zz,nH,Td,kRoss,tauRoss,Hvis,dustgas
      use DUST_DATA,ONLY: NELEM,NDUST,eps0
      use OPACITY,ONLY: NLAMmax
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: hplanck=6.62607554E-27
      real,parameter :: bk=1.3806581200000E-16
      real,parameter :: cl=2.9979245800000E+10 
      real,parameter :: cPl1=2.0*hplanck*cl**2
      real,parameter :: cPl2=hplanck*cl/bk 
      real,parameter :: sig_SB=cPl1/cPl2**4*pi**5/15.0
      integer,parameter :: qp=selected_real_kind(33,4931)
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real,dimension(NLAMmax) :: kabs,ksca,kext
      real,dimension(Nz) :: Diff
      real :: nHtot,T,Told,dg,qual1,qual2,kapROSS
      real :: z1,z2,J0,J1,J2,D1,D2,dz1,dz2,term1,term2,f0,f1
      integer :: iz,verb,it
      character(len=1) :: char1

      !-------------------------------------------------------
      ! ***  check validity of 1d diffusion approximation  ***
      !-------------------------------------------------------
      do iz=1,Nz
        Diff(iz) = 4.0*pi/3.0/kRoss(iz)
      enddo
      do iz=izStart-1,2,-1
        z1  = 0.5*(zz(iz-1)+zz(iz))     ! zz(iz-1/2)
        z2  = 0.5*(zz(iz+1)+zz(iz))     ! zz(iz+1/2)
        J0  = sig_SB/pi*Td(iz)**4
        J1  = sig_SB/pi*Td(iz-1)**4
        J2  = sig_SB/pi*Td(iz+1)**4
        D1  = SQRT(Diff(iz)*Diff(iz-1))
        D2  = SQRT(Diff(iz)*Diff(iz+1))
        dz1 = zz(iz)-zz(iz-1)
        dz2 = zz(iz+1)-zz(iz)
        term1 = D1*(J0-J1)/dz1 + D2*(J0-J2)/dz2
        term2 = Hvis(iz)*(z2-z1)
        print'(I4,0pF8.1,9(1pE12.4))',iz,Td(iz),kRoss(iz),term1,term2
      enddo

      print*
      call INIT_OPAC
      print*
      verbose = -2
      verb    = -2
      do iz=izStart,izStart-2,-1
        nHtot = nH(iz)
        T = Td(iz)
        eps = eps0
        call EQUIL_COND(nHtot,T,eps,Sat,eldust,verbose)
        call CALC_OPAC(nHtot,eldust,kabs,ksca,kext,dg,-1)
        print'(I4,0pF9.1,9(1pE12.4))',iz,T,dustgas(iz),dg,
     >                         kRoss(iz),KapROSS(T,kext)
        kRoss(iz) = KapROSS(T,kext)
        Diff(iz) = 4.0*pi/3.0/kRoss(iz)
      enddo

      do iz=izStart-1,2,-1
        z1  = 0.5*(zz(iz-1)+zz(iz))     ! zz(iz-1/2)
        z2  = 0.5*(zz(iz+1)+zz(iz))     ! zz(iz+1/2)
        J0  = sig_SB/pi*Td(iz)**4
        J2  = sig_SB/pi*Td(iz+1)**4
        D2  = SQRT(Diff(iz)*Diff(iz+1))
        dz1 = zz(iz)-zz(iz-1)
        dz2 = zz(iz+1)-zz(iz)
        f0  = Hvis(iz)*(z2-z1) - D2*(J0-J2)/dz2
        nHtot = nH(iz-1)
        T = Td(iz)                      ! initial guess
        print*
        print'(I10,0pF11.4,2(1pE14.6))',iz-1,
     >       Td(iz-1),kRoss(iz-1),dustgas(iz-1)
        do it=1,10
          eps = eps0
          call EQUIL_COND(nHtot,T,eps,Sat,eldust,verbose)
          call CALC_OPAC(nHtot,eldust,kabs,ksca,kext,dg,-2)
          kRoss(iz-1) = KapROSS(T,kext)
          Diff(iz-1) = 4.0*pi/3.0/kRoss(iz-1)
          D1 = SQRT(Diff(iz)*Diff(iz-1))
          J1 = sig_SB/pi*T**4
          f1 = D1*(J0-J1)/dz1
          J1 = J0 - f0*dz1/D1 
          Told  = T
          T     = (J1*pi/sig_SB)**0.25
          qual1 = T/Told-1.0
          qual2 = (f0-f1)/(Hvis(iz)*(z2-z1))
          print'(I3,0pF11.4,2(1pE14.6),3(1pE9.1))',it,Told,
     >         kRoss(iz-1),dg,f0-f1,qual1,qual2
          if (ABS(qual2)<1.E-5) exit
        enddo
        Td(iz-1) = T
        call CALC_OPAC(nHtot,eldust,kabs,ksca,kext,dg,-1)
      enddo
      end
      
!-----------------------------------------------------------------------
      real function KapROSS(T,kext)
!-----------------------------------------------------------------------
      use OPACITY,ONLY: NLAMmax,NLAM,lam
      implicit none
      real,parameter :: mic=1.E-4
      real,parameter :: cl=2.9979245800000E+10
      real,intent(in) :: T,kext(NLAMmax)
      real,external :: DPLANCK
      real :: sum1,sum2,nu1,nu2,dnu,f1a,f1b,f2a,f2b
      integer :: i
      sum1 = 0.0
      sum2 = 0.0
      nu1 = cl/(lam(1)*mic)
      f1a = DPLANCK(T,nu1)
      f2a = f1a/kext(1)
      do i=2,NLAM
        nu2 = cl/(lam(i)*mic)
        f1b = DPLANCK(T,nu2)
        f2b = f1b/kext(i)
        dnu = nu1-nu2
        sum1 = sum1 + 0.5*(f1a+f1b)*dnu
        sum2 = sum2 + 0.5*(f2a+f2b)*dnu
        nu1 = nu2
        f1a = f1b
        f2a = f2b
        !print'(0pF9.2,9(1pE11.4))',lam(i),kext(i),f1a*dnu,f2a*dnu
      enddo
      KapROSS = sum1/sum2
      end
      
!-----------------------------------------------------------------------
      real function BPLANCK(T,nu)
!-----------------------------------------------------------------------
      implicit none
      real,parameter :: hplanck=6.62607554E-27
      real,parameter :: cl=2.9979245800000E+10
      real,parameter :: bk=1.3806581200000E-16 
      real,intent(in) :: T,nu
      real :: x
      real,save :: Bconst,one12,one720
      logical,save :: firstCall=.true.

      if (firstCall) then
        Bconst = 2.0*hplanck/cl**2
        one12  = 1.0/12.0
        one720 = 1.0/720.0
        firstCall=.false.
      endif  
      ! BPLANCK = 2*h/c^2 * v^3 / (exp(h*v/(k*T))-1)
      !         = 2*h/c^2 * v^3 / (exp(x)-1)
      x = hplanck*nu/(bk*T)
      if (x.gt.50.0) then
        BPLANCK = EXP(-x)
      else if (x.lt.0.1) then 
        BPLANCK = 1.0/x -0.5 +one12*x -one720*x**3
      else
        BPLANCK = 1.0 / ( EXP(x) - 1.0 )
      endif
      BPLANCK = Bconst * nu**3 * BPLANCK
      end

!-----------------------------------------------------------------------
      real function DPLANCK(T,nu)
!-----------------------------------------------------------------------
      implicit none
      real,parameter :: hplanck=6.62607554E-27
      real,parameter :: cl=2.9979245800000E+10
      real,parameter :: bk=1.3806581200000E-16 
      real,intent(in) :: T,NU
      real :: x,ex
      real,save :: Dconst,one12,one240
      logical,save :: firstCall=.true.

      if (firstCall) then
        Dconst = 2.0*bk/cl**2 
        one12  = 1.0/12.0
        one240 = 1.0/240.0
        firstCall=.false.
      endif  
      ! DPLANCK = 2*h/c^2 * v^3 / [EXP(h*v/(k*T))-1]^2
      !                          * EXP(h*v/(k*T)) * h*v/k / T^2
      !         = 2*h^2/(c^2*k*T^2) * v^4 * exp(x)/[exp(x)-1]^2
      !         = 2*h^2/(c^2*k*T^2) * v^4 * x^2/(hv/kT)^2 * exp(x)/[exp(x)-1]^2
      !         = 2*k/c^2 * v^2 * x^2*exp(x)/[exp(x)-1]^2
      x = hplanck*nu/(bk*T)
      if (x.gt.50.0) then
        DPLANCK = x**2 * EXP(-x)
      else if (x.lt.0.1) then
        DPLANCK = 1.0 -one12*x**2 +one240*x**4
      else
        ex = EXP(x) 
        DPLANCK = x**2 * ex/(ex-1.0)**2
      endif
      DPLANCK = Dconst * nu**2 * DPLANCK
      end
