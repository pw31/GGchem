************************************************************************
      module OPACITY
************************************************************************
      use DUST_DATA,ONLY: NDUST
      integer,parameter :: NLAMmax=1000,NSIZEmax=500,NOPmax=200
      integer :: NLAM,NSIZE,NLIST
      real,dimension(NLAMmax) :: lam                ! wavelength[mic]
      real,dimension(NSIZEmax) :: aa,ff,aweight     ! size dist.func.
      real,allocatable,dimension(:,:) :: nn,kk      ! optical constants
      integer :: opind(NOPmax),duind(NOPmax)
      character(len=100) :: opfile(NOPmax)
      logical :: conducting(NOPmax)
      real :: xmin=1.0,xmax=0.0,kmin=1.0,kmax=0.0,nmin=1.0,nmax=0.0
      end 

************************************************************************
      subroutine INIT_OPAC
************************************************************************
      use DUST_DATA,ONLY: NDUST,dust_nam
      use OPACITY,ONLY: NLAM,lam,nn,kk,NLIST,opind,duind,
     >                  opfile,conducting
      implicit none
      real,parameter :: mic=1.E-4
      character(len=200) :: filename,line
      character(len=20) :: opname
      real :: lmin,lmax,kmin,kmax,nmin,nmax
      real,allocatable :: nread(:),kread(:)
      integer :: i,j,k
      logical :: ex

      print*
      print*,"INIT_OPAC ..."
      !----- wavelength grid -----
      NLAM = 300
      allocate(nn(NLAM,NDUST),kk(NLAM,NDUST),nread(NLAM),kread(NLAM))
      lmin = 0.3
      lmax = 300.0
      do i=1,NLAM
        lam(i) = EXP(LOG(lmin)+(i-1.0)/(NLAM-1.0)*LOG(lmax/lmin))
        !print*,i,lam(i)/mic
      enddo
      !----- default: set to vaccum -----
      nn(1:NLAM,:) = 1.0
      kk(1:NLAM,:) = 0.0
      opind(:) = 0
      duind(:) = 0

      !----- identify opacity species and read optical constants -----
      filename = "data/OpticalData/master.list"
      open(unit=12,file=filename,status='old')
      NLIST = 0
      nmin = 1.0
      nmax = 1.0
      kmin = 1.0
      kmax = 1.0
      do i=1,9999
        read(12,'(A200)',end=100) line
        opname = trim(line(1:24))
        read(line(25:26),*) conducting(i)
        opfile(i) = trim(line(27:))
        !print*,trim(opname),trim(opfile(i)),conducting(i)
        do j=1,NDUST
          if (dust_nam(j)==opname) exit
        enddo
        if (j>NDUST) then
          print*,"*** opacity species not found:",opname
          stop
        else
          print'(I3,I4,A40,A30," conducting=",L1)',i,j,
     >         " opacity species found: "//trim(opname),
     >         trim(opfile(i)),conducting(i)
          filename = "data/OpticalData/"//trim(opfile(i))
          inquire(file=filename,exist=ex)
          if (.not.ex) then
            print*,"*** file does not exist: "//filename
            stop
          endif
          call FETCH_OPTICALDATA(filename,conducting(i),nread,kread)
          nn(1:NLAM,i) = nread(1:NLAM) 
          kk(1:NLAM,i) = kread(1:NLAM)
          do k=1,NLAM
            nmax = MAX(nmax,nread(k))
            nmin = MIN(nmin,nread(k))
            kmax = MAX(kmax,kread(k))
            kmin = MIN(kmin,kread(k))
          enddo
          opind(j) = i
          duind(i) = j
          NLIST = i
        endif
      enddo
 100  close(12)
      !--------- add vaccum for porosity -----------
      NLIST = NLIST+1
      nn(1:NLAM,NLIST) = 1.0
      kk(1:NLAM,NLIST) = 0.0
      print*,"nmin,nmax=",nmin,nmax
      print*,"kmin,kmax=",kmin,kmax
      end

************************************************************************
      subroutine FETCH_OPTICALDATA(fileName,conducting,n_mono,k_mono)
************************************************************************
      use OPACITY,ONLY: NLAM,lam
      implicit none
      character(len=*),intent(IN) :: fileName
      logical,intent(IN) :: conducting
      real,dimension(10000) :: lread,nread,kread
      real,dimension(NLAM),intent(OUT) :: n_mono,k_mono
      real :: fac
      integer :: i,j,j0,j1,Ndat
      character(len=9999) :: line
      
      !-------------------------------
      !***  read the opacity file  ***
      !-------------------------------
      open(13,file=fileName,status='old')
      Ndat=1
 100  continue
      read(13,'(A9999)') line
      read(line,*,err=100,end=100) lread(1),nread(1),kread(1)
      !print*,lread(1),nread(1),kread(1)
      do
        read(13,*,end=200) lread(Ndat+1),nread(Ndat+1),kread(Ndat+1)
        if (lread(Ndat+1).le.lread(Ndat)) then
          write(*,*) "*** lambda not monoton increasing."
          write(*,*) Ndat,lread(Ndat),lread(Ndat+1)
          stop
        endif
        Ndat = Ndat+1
      enddo
 200  close(13)
      write(*,'(i5," datapoints between",0pF10.5," and ",
     >          0pF11.2," mic.")') Ndat,lread(1),lread(Ndat)
      
      !-----------------------------------------
      !***  interpolation and extrapolation  ***
      !-----------------------------------------
      j=1
      do i=1,NLAM
        do
          if (lread(j+1).gt.lam(i)) exit
          if (j.eq.Ndat) exit
          j=j+1
        enddo
        !write(*,*) lread(j),lam(i),lread(j+1)
        if (lam(i).lt.lread(1)) then
          !--- UV-extrapolation ---
          !write(*,*) "UV-extrapolation"
          !wp=(1.0-nread(1))/lread(1)**2
          !gamma=kread(1)/lread(1)**3
          !n_mono(i,kind) = 1.0-wp*lam(i)**2
          !k_mono(i,kind) = gamma*lam(i)**3
          n_mono(i) = nread(1)
          k_mono(i) = kread(1)
        else if (lam(i).gt.lread(Ndat)) then
          !--- cm-extrapolation ---
          !write(*,*) "cm-extrapolation"
          n_mono(i) = nread(Ndat)
          k_mono(i) = kread(Ndat)*lread(Ndat)/lam(i)
          if (conducting) then
            !--- log/log extraploation for conducting materials ---
            j0=Ndat
            do j1=Ndat,1,-1                         ! data can be noisy, so
              if (lread(j1).lt.0.7*lread(j0)) exit  ! it's safer to use larger
            enddo                                   ! region to get the slope
            !print *,j1,j0,lread(j1),lread(j0)
            fac = LOG(lam(i)/lread(j0))/LOG(lread(j1)/lread(j0))
            n_mono(i) = EXP(LOG(nread(j0))
     >                 +fac*LOG(nread(j1)/nread(j0)))
            k_mono(i) = EXP(LOG(kread(j0))
     >                 +fac*LOG(kread(j1)/kread(j0)))
          endif
        else
          !--- interpolation ---
          fac = LOG(lam(i)/lread(j))/LOG(lread(j+1)/lread(j))
          n_mono(i) = EXP(LOG(nread(j))
     >               +fac*LOG(nread(j+1)/nread(j)))
          if (kread(j).eq.0.0.or.kread(j+1).eq.0.0) then
            k_mono(i) = 0.0
          else
            k_mono(i) = EXP(LOG(kread(j))
     >                 +fac*LOG(kread(j+1)/kread(j)))
          endif
        endif
      enddo
      end
