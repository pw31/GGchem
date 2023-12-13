************************************************************************
      module DISK
************************************************************************
      integer,parameter :: NXmax=500, NZmax=500
      integer :: Nx,Nz
      real :: Mstar,Rstar,Mdot
      real,dimension(NXmax,NZmax) :: rr,zz,Td,nH,JJ,DD,kRoss,tauRossV,
     >                               tauRossR,Hvis,dustgas
      end
      
************************************************************************
      SUBROUTINE INIT_DISK
************************************************************************
      use DISK,ONLY: NXmax,Nx,Nz,Mstar,Rstar,Mdot,rr,zz,
     >               nH,Td,kRoss,tauRossV,tauRossR,Hvis,dustgas
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: grav=6.67259850000E-08
      real,parameter :: AU=1.495978700E+13
      real,parameter :: Msun=1.988922500E+33
      real,parameter :: Rsun=6.959900000E+10
      real,parameter :: yr=3.155760000E+7 
      integer :: i,ix,iz,Nsp,Nheat,Ncool,Nlam
      real :: r1,r2,z1,z2,z3,f1,f2,f3,Hcol,Hint,const,cc
      real :: arr(1000),tauR(NXmax),Hread(NXmax)
      character(len=99999) :: line,infile
      logical :: ex
      
      !----------------------------
      ! ***  read Parameter.in  ***
      !----------------------------
      inquire(file="Parameter.in",exist=ex)
      if (.not.ex) stop "*** need Parameter.in"
      print*,"reading Parameter.in ..."
      open(unit=12,file="Parameter.in",status='old')
      Mdot = 0.0
      infile = ''
      do i=1,99999
        read(12,'(A99999)',end=200) line
        if (index(line,'! Mdot')>0) read(line,*) Mdot
        if (index(line,'! fixed_surface_density')>0) then
          read(line,*) ex
          if (ex) then
            read(12,*) infile
          endif
        endif
      enddo
 200  continue
      if (infile=='') then
        print*,"found Mdot=",Mdot
        Mdot = Mdot*Msun/yr
      else
        print*,"viscous heating rates from "//trim(infile)
      endif

      !---------------------------
      ! ***  read ProDiMo.out  ***
      !---------------------------
      inquire(file="ProDiMo.out",exist=ex)
      if (.not.ex) stop "*** need ProDiMo.out"
      print*
      print*,"reading ProDiMo.out ..."
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
      !print'(A4,A4,2A11,A9,99A11)','ix','iz','r[AU]','z[AU]',
     >!     'Td[K]','n<H>[cm-3]','kRoss','tauRoss','d/g'
      do iz=Nz,1,-1
        tauR(:) = 0.0
        do ix=1,Nx
          read(12,'(A99999)') line
          read(line,*,end=100) arr
 100      rr(ix,iz) = arr(3)*AU
          zz(ix,iz) = arr(4)*AU
          Td(ix,iz) = arr(11)
          nH(ix,iz) = arr(22+Nheat+Ncool)
          kRoss(ix,iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+2)
          tauRossV(ix,iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+3)
          dustgas(ix,iz) = arr(22+Nheat+Ncool+Nsp+1+Nlam+5)
          if (ix>1) then
            r1 = SQRT(rr(ix-1,iz)**2+zz(ix-1,iz)**2)
            r2 = SQRT(rr(ix,iz)**2+zz(ix,iz)**2)
            tauRossR(ix,iz) = tauRossR(ix-1,iz)
     >           + 0.5*(kRoss(ix,iz)+kRoss(ix-1,iz))*(r2-r1)
          endif  
          !print'(I4,I4,2(0pF11.6),0pF9.2,99(1pE11.4))',ix,iz,
     >    !     rr(ix,iz)/AU,zz(ix,iz)/AU,Td(ix,iz),nH(ix,iz),
     >    !     kRoss(ix,iz),tauRossV(ix,iz),dustgas(ix,iz)     
        enddo
      enddo
      close(12)
      print*,"... done reading ProDiMo.out"
      
      !---------------------------------------------
      ! ***  compute viscous heating per volume  ***
      !---------------------------------------------
      Hread = 0.0
      if (infile.ne.'') then
        open(unit=12,file=infile,status='old')
        do
 400      read(12,'(A99999)',end=500) line
          read(line,*,err=400,end=400) arr(1:5)
          do ix=1,Nx
            if (ABS(arr(1)*AU/rr(ix,1)-1.0)<1.E-6) then
              Hread(ix) = arr(3)
              print*,ix,nH(ix,1),Hread(ix)
            endif
          enddo
        enddo
 500    close(12)
        do ix=Nx,1,-1
          if (Hread(ix)==0.0) then
            Hread(ix) = Hread(ix+1)*nH(ix,1)/nH(ix+1,1)
            print*,ix,nH(ix,1),Hread(ix)
          endif
        enddo
      endif
        
      do ix=1,Nx
        Hcol = Hread(ix)/2.0
        if (Hcol==0.0) then
          Hcol = grav*Mstar*Mdot/rr(ix,1)**3 * 3.0/(8.0*pi) ! [erg/cm2/s]
     >         * (1.0-SQRT(Rstar/rr(ix,1)))
        endif
        Hint = 0.0
        do iz=1,Nz-1
          z1 = zz(ix,iz)                                  ! [cm]
          z3 = zz(ix,iz+1)
          z2 = 0.5*(z1+z3)
          f1 = nH(ix,iz)*dustgas(ix,iz)                   ! [cm-3]
          f3 = nH(ix,iz+1)*dustgas(ix,iz+1)
          !------ fit rho(z) = rho(z1)*exp(-(z-z1)^2/H^2) --------
          cc = LOG(f3/f1)/(z3-z1)**2                      ! -1/H^2
          f2 = f1*EXP(cc*(z2-z1)**2)   
          Hint = Hint + (f1+4*f2+f3)*(z3-z1)/6.0          ! [cm-2]
        enddo  
        const = Hcol/Hint                                 ! [erg/s]
        do iz=1,Nz
          Hvis(ix,iz) = const*nH(ix,iz)*dustgas(ix,iz)    ! [erg/cm3/s]
        enddo
      enddo
        
      end

************************************************************************
      SUBROUTINE COMPUTE_DISK
************************************************************************
      use PARAMETERS,ONLY: verbose
      use DISK,ONLY: NXmax,Nx,Nz,rr,zz,nH,Td,kRoss,tauRossV,tauRossR,
     >               Hvis,dustgas,Mstar,Rstar,Mdot
      use DUST_DATA,ONLY: NELEM,NDUST,eps0,dust_nam,muH,
     >                    dust_mass,dust_Vol
      use OPACITY,ONLY: NLAMmax,NLAM,lam,xmin,xmax,kmin,kmax,nmin,nmax
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: grav=6.67259850000E-08
      real,parameter :: yr=3.155760000E+7 
      real,parameter :: hplanck=6.62607554E-27
      real,parameter :: bk=1.3806581200000E-16
      real,parameter :: cl=2.9979245800000E+10 
      real,parameter :: cPl1=2.0*hplanck*cl**2
      real,parameter :: cPl2=hplanck*cl/bk 
      real,parameter :: sig_SB=cPl1/cPl2**4*pi**5/15.0
      real,parameter :: AU=1.495978700E+13
      integer,parameter :: qp=selected_real_kind(33,4931)
      real,parameter :: Tminval=55.0,Tmaxval=5000.0
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real,dimension(NLAMmax) :: kabs,ksca,kext
      real,dimension(Nz) :: Diff,Fup
      real :: nHtot,rho,rhod,Vd,T,Told,dg,qual1,qual2,kRoss_GGchem
      real :: z1,z2,J0,J1,J2,D1,D2,dz1,dz2,term1,term2,f0,f1,qual,qold
      real :: Tmin,Tmax,fac,Tcrit,Twork,kRdust,kRgas
      real,external :: kapROSS,GAS_OPACITY
      integer :: iz0(NXmax),iz1(NXmax),i,ix,iz,it,iz_lastdust,totit
      character(len=1) :: char1,bem
      logical :: first

      iz0 = 0
      iz1 = 0
      do ix=1,Nx
        do iz=Nz,1,-1
          if (nH(ix,iz)>1.E+7) iz0(ix)=MAX(iz0(ix),iz)
          if (tauRossV(ix,iz)>10.0.and.tauRossR(ix,iz)>10.0)
     >                         iz1(ix)=MAX(iz1(ix),iz)
        enddo
        if (iz1(ix)>0) then
          do iz=iz1(ix),iz1(ix)+3
            !print*,iz,Td(ix,iz)/Td(ix,iz+1)
            if (Td(ix,iz)<1.25*Td(ix,iz+1)) exit
          enddo
          !print*,iz1(ix),iz
          iz1(ix)=iz
        endif
        !-------------------------------------------------------
        ! ***  check validity of 1d diffusion approximation  ***
        !-------------------------------------------------------
        do iz=1,Nz
          Diff(iz) = 4.0*pi/3.0/kRoss(ix,iz)
        enddo
        qual = 0.0
        do iz=iz1(ix)-1,2,-1
          z1  = 0.5*(zz(ix,iz-1)+zz(ix,iz))     ! zz(iz-1/2)
          z2  = 0.5*(zz(ix,iz+1)+zz(ix,iz))     ! zz(iz+1/2)
          J0  = sig_SB/pi*Td(ix,iz)**4
          J1  = sig_SB/pi*Td(ix,iz-1)**4
          J2  = sig_SB/pi*Td(ix,iz+1)**4
          D1  = SQRT(Diff(iz)*Diff(iz-1))
          D2  = SQRT(Diff(iz)*Diff(iz+1))
          dz1 = zz(ix,iz)-zz(ix,iz-1)
          dz2 = zz(ix,iz+1)-zz(ix,iz)
          term1 = D1*(J0-J1)/dz1 + D2*(J0-J2)/dz2
          term2 = Hvis(ix,iz)*(z2-z1)
          qual = qual+ABS(term2/term1-1.0)**0.3
          !print'(2(I4),0pF8.1,9(1pE12.4))',ix,iz,Td(ix,iz),
     >    !     kRoss(ix,iz),term1,term2
        enddo
        Tmin = MINVAL(Td(ix,1:Nz))
        if (iz1(ix)>2) then
          qual  = (qual/(iz1(ix)-2))**(1.0/0.3)
          Tcrit = Td(ix,1)/Tmin
        else
          qual  = 1.E+99
          Tcrit = 0.0
        endif
        first = (SUM(iz1(1:ix-1))==0)
        if (Td(ix,1)<30.0) then
          !--- too cold for GGchem ---
          iz0(ix) = 0
          iz1(ix) = 0
        else if (Tcrit<1.1.or.(first.and.qual>0.2*Tcrit)) then
          !--- vertical diffusion not valid ---
          iz1(ix) = 0
          print'(I4,0pF8.3,2I4," qual=",1pE9.2," Tcrit=",1pE9.2,
     >         " Tmin=",0pF7.1," Tmid=",0pF7.1," *")',
     >         ix,rr(ix,1)/AU,iz0(ix),iz1(ix),
     >         qual,Tcrit,Tmin,Td(ix,1)
        else
          print'(I4,0pF8.3,2I4," qual=",1pE9.2," Tcrit=",1pE9.2,
     >         " Tmin=",0pF7.1," Tstart=",0pF7.1," Tmax=",0pF7.1)',
     >         ix,rr(ix,1)/AU,iz0(ix),iz1(ix),qual,Tcrit,Tmin,
     >         Td(ix,iz1(ix)),MAXVAL(Td(ix,1:iz1(ix)))
        endif
      enddo
      print*,SUM(iz0(1:Nx))," phase.eq. & opacity calc."
      print*,SUM(iz1(1:Nx))," temp.iter., phase.eq. & opacity calc."      
      call INIT_OPAC      
      print*
      
      open(unit=9,file='disk.out',status='replace')
      write(9,*) Nx,Nz,NLAM,NDUST
      write(9,'(9999(1pE12.3))') lam(1:NLAM)
      write(9,'(9999(1pE12.3))') dust_mass(1:NDUST)
      write(9,'(9999(1pE12.3))') dust_Vol(1:NDUST)
      write(9,'(2(A4),A8,A10,9999(A16))')
     >     'ix','iz','z/r','T[K]','n<H>[cm-3]','rho[g/cm3]',
     >     'rhod[g/cm3]','Vdust[cm3/cm3]','kapR[cm2/g]',
     >     ('kabs[cm2/g]',i=1,NLAM),
     >     (trim(dust_nam(i)),i=1,NDUST)
      verbose = -2              ! avoid output from GGchem
      totit = 0
      do ix=1,Nx  !Nx
        print*
        print*,"==> new radius",ix,rr(ix,1)/AU,iz0(ix),iz1(ix)
        !Hcol = grav*Mstar*Mdot/rr(ix,1)**3 * 3.0/(8.0*pi) ! [erg/cm2/s]
     >  !     * (1.0-SQRT(Rstar/rr(ix,1)))
        iz_lastdust = 0
        do iz=iz0(ix),1,-1
          nHtot = nH(ix,iz)
          rho = nHtot*muH
          eps = eps0
          T = Td(ix,iz)
          Twork = MAX(MIN(T,Tmaxval),Tminval)
          call EQUIL_COND(nHtot,Twork,eps,Sat,eldust,verbose)
          call CALC_OPAC(nHtot,eldust,kabs,ksca,kext,dg,-1)
          rhod = 0.0
          Vd   = 0.0
          do i=1,NDUST
            rhod = rhod + nHtot*eldust(i)*dust_mass(i)
            Vd   = Vd   + nHtot*eldust(i)*dust_Vol(i)
          enddo
          kRdust = KapROSS(T,kext)
          kRgas  = GAS_OPACITY(T,nHtot)*rho
          kRoss_GGchem = kRdust+kRgas
          print'(2I4,0pF9.1,9(1pE12.4))',ix,iz,T,dustgas(ix,iz),dg,
     >                                kRoss(ix,iz)/rho,kRdust,kRgas
          kRoss(ix,iz) = kRoss_GGchem
          Diff(iz) = 4.0*pi/3.0/kRoss_GGchem
          write(9,'(2(I4),0pF8.5,0pF10.3,9999(1pE16.6E3))') ix,iz,
     >       zz(ix,iz)/rr(ix,iz),T,nHtot,rho,rhod,Vd,
     >       kRoss(ix,iz)/rho,kabs(1:NLAM)/rho,nHtot*eldust(1:NDUST)
          if (kRoss_GGchem>0.0) iz_lastdust=iz
          if (iz<=iz1(ix)-1.and.iz_lastdust>0) exit
        enddo
        if (kRoss_GGchem==0.0.and.iz1(ix)>0) then
          print*," !! no dust at start of integration"
          if (iz_lastdust==0) then
            print*," !!! no dust at all in this column!"
            iz1(ix) = 0
          else
            iz1(ix)=iz_lastdust+2
          endif
        endif
        Fup(1) = 0.0                            ! flux [erg/cm2/s] at iz
        do iz=2,iz1(ix)
          dz1 = zz(ix,iz)-zz(ix,iz-1)
          Fup(iz) = Fup(iz-1) + 0.5*(Hvis(ix,iz)+Hvis(ix,iz-1))*dz1
          !print*,iz,Fup(iz)/Hcol
        enddo
        do iz=iz1(ix)-1,2,-1
          z1  = 0.5*(zz(ix,iz-1)+zz(ix,iz))     ! zz(iz-1/2)
          z2  = 0.5*(zz(ix,iz+1)+zz(ix,iz))     ! zz(iz+1/2)
          J0  = sig_SB/pi*Td(ix,iz)**4
          J2  = sig_SB/pi*Td(ix,iz+1)**4
          D2  = SQRT(Diff(iz)*Diff(iz+1))
          dz1 = zz(ix,iz)-zz(ix,iz-1)
          dz2 = zz(ix,iz+1)-zz(ix,iz)
          f0  = Hvis(ix,iz)*(z2-z1) - D2*(J0-J2)/dz2
          nHtot = nH(ix,iz-1)
          rho = nHtot*muH
          T = Td(ix,iz)                         ! initial guess
          print*
          print'(I10,I4," T=",0pF11.4," kRos=",1pE9.3," d/g=",1pE9.3,
     >         " Hvis=",1pE9.3)',ix,iz-1,Td(ix,iz-1),
     >         kRoss(ix,iz-1)/rho,dustgas(ix,iz-1),Hvis(ix,iz)
          Tmin = 0.0
          Tmax = 9.e+99
          qual = 9.e+99
          do it=1,99
            eps = eps0
            Twork = MAX(MIN(T,Tmaxval),Tminval)
            call EQUIL_COND(nHtot,Twork,eps,Sat,eldust,verbose)
            call CALC_OPAC(nHtot,eldust,kabs,ksca,kext,dg,-2)
            kRdust = KapROSS(T,kext)
            kRgas  = GAS_OPACITY(T,nHtot)*rho
            kRoss(ix,iz-1) = kRdust + kRgas
            Diff(iz-1) = 4.0*pi/3.0/kRoss(ix,iz-1)
            D1 = SQRT(Diff(iz)*Diff(iz-1))
            J1 = J0 + 0.5*(Fup(iz-1)+Fup(iz))*dz1/D1   ! Fup = -D1*(J0-J1)/dz
            Told = T
            T    = (J1*pi/sig_SB)**0.25
            qold = qual
            qual = LOG(T/Told)
            print'(I3," T=",0pF10.4," kapR=",2(1pE12.6)," d/g=",1pE12.6,
     >             " q=",1pE9.2)',it,Told,kRoss(ix,iz-1)/rho,dg,qual
            if (T>Told) Tmin=MAX(Tmin,Told)            ! Tsolution>Told
            if (T<Told) Tmax=MIN(Tmax,Told)            ! Tsolution<Told
            if (it>1.and.qold*qual<0.0) then
              T = Told + (T-Told)*qold/(qold-qual)
              !print'(" qold,qual,Told,T=",2(1pE10.2),2(0pF13.6))',
     >        !       qold,qual,Told,T
            endif
            if (Tmin>0.0.and.Tmax<1.E+99) then
              fac = (T-Tmin)/(Tmax-Tmin)
              bem = " "
              if (it>3) then
                if (fac<0.3) T=0.7*Tmin+0.3*Tmax
                if (fac>0.7) T=0.3*Tmin+0.7*Tmax
                if (fac<0.3.or.fac>0.7) bem="*"
              endif  
              !print'(" Tmin,Tmax,fac=",2(0pF13.6),1pE10.2,
     >        !       " ==>",0pF13.6,A2)',Tmin,Tmax,fac,T,bem
            endif
            if (ABS(qual)<1.E-8) exit
          enddo
          if (it>99) print*,"*** WARNING T-iteration not converged"
          if (it>99) stop
          totit = totit + it
          Td(ix,iz-1) = T
          call CALC_OPAC(nHtot,eldust,kabs,ksca,kext,dg,-1)
          rhod = 0.0
          Vd   = 0.0
          do i=1,NDUST
            rhod = rhod + nHtot*eldust(i)*dust_mass(i)
            Vd   = Vd   + nHtot*eldust(i)*dust_Vol(i)
          enddo
          write(9,'(2(I4),0pF8.5,0pF10.3,9999(1pE16.6E3))') ix,iz-1,
     >         zz(ix,iz-1)/rr(ix,iz-1),T,nHtot,rho,rhod,Vd,
     >         kRoss(ix,iz-1)/rho,kabs(1:NLAM)/rho,nHtot*eldust(1:NDUST)
          read(*,'(A1)') char1
        enddo
      enddo
      close(9)
      print*
      print*,"--------------  D O N E  ---------------"
      print*,"nmin,nmax=",nmin,nmax
      print*,"kmin,kmax=",kmin,kmax
      print*,"xmin,xmax=",xmin,xmax
      print*,totit
      
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
