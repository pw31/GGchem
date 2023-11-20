***********************************************************************
      module DISK
***********************************************************************
      use OPACITY,ONLY: NLAMmax
      integer,parameter :: NZmax=1000
      integer :: Nx,Nz,izStart
      real,dimension(NZmax) :: zz,Td,nH,JJ,DD,kRoss,tauRoss,Hvis
      real,dimension(NLAMmax,NZmax) :: kex
      end
      
***********************************************************************
      SUBROUTINE INIT_DISK
***********************************************************************
      use DISK,ONLY: Nx,Nz,izStart,zz,nH,Td,kRoss,tauRoss,Hvis
      implicit none
      real,parameter :: AU=1.495978700E+13
      integer :: i,ix,iz,ixSearch,Nsp,Nheat,Ncool,Nlam
      character(len=9999) :: line
      real :: rr,arr(1000)

      print*
      print*,"reading ProDiMo.out ..."
      ixSearch = 20
      izStart  = 0
      open(unit=12,file="ProDiMo.out",status='old')
      do i=1,21
        read(12,*)
      enddo
      read(12,'(A9999)') line
      read(line(25:),*) Nx,Nz,Nsp,Nheat,Ncool,Nlam
      !print*,Nx,Nz,Nsp,Nheat,Ncool
      read(12,*)
      read(12,*)
      print'(A4,A4,2A11,A9,99A11)','ix','iz','r[AU]','z[AU]',
     >     'Td[K]','n<H>[cm-3]','kRoss','tauRoss'
      do iz=Nz,1,-1
        do ix=1,Nx
          read(12,'(A9999)') line
          read(line,*,end=100) arr
 100      if (ix==ixSearch) then
            rr     = arr(3)*AU
            zz(iz) = arr(4)*AU
            Td(iz) = arr(11)
            nH(iz) = arr(22+Nheat+Ncool)
            kRoss(iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+2)
            tauRoss(iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+3)
            if (tauRoss(iz)>10.0) then
              izStart=MAX(iz,izStart)
            endif
            print'(I4,I4,2(0pF11.6),0pF9.2,99(1pE11.4))',ix,iz,
     >           rr/AU,zz(iz)/AU,Td(iz),nH(iz),kRoss(iz),tauRoss(iz)
          endif
        enddo
      enddo

      end

***********************************************************************
      SUBROUTINE COMPUTE_DISK
***********************************************************************
      use DISK,ONLY: izStart,zz,nH,Td,kRoss,tauRoss,Hvis
      use DUST_DATA,ONLY: NELEM,NDUST,eps0
      use OPACITY,ONLY: NLAMmax
      implicit none
      integer,parameter :: qp=selected_real_kind(33,4931)
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      integer :: iz,verbose
      real*8,dimension(NLAMmax) :: kabs,ksca,kext
      real :: nHtot,T
      
      call INIT_OPAC
      verbose = 0
      do iz=1,izStart
        nHtot = nH(iz)
        T = Td(iz)
        eps = eps0
        call EQUIL_COND(nHtot,T,eps,Sat,eldust,verbose)
        call CALC_OPAC(nHtot,eldust,kabs,ksca,kext)
      enddo
      end
      
!----------------------------------------------------------------------
      real function KapROSS(T,kext)
!----------------------------------------------------------------------
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
      
!----------------------------------------------------------------------
      real function BPLANCK(T,nu)
!----------------------------------------------------------------------
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

!----------------------------------------------------------------------
      real function DPLANCK(T,nu)
!----------------------------------------------------------------------
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
