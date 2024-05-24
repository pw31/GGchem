      module FASTMIEDIM
      use DATATYPE,ONLY: r2     ! from MIEX
      integer,parameter :: Nx=200,Nn=60,Nk=200
      real,parameter :: nmin=0.8,nmax=20.0
      real,parameter :: kmin=1.E-4,kmax=10.0      ! extrap to lower k works
      real,parameter :: xmin=1.E-3,xmax=10000.0   ! extrap to lower x works
      real*4,dimension(0:Nx,0:Nn,0:Nk) :: lQex,lQsc,lQab
      end
      
************************************************************************
      subroutine FASTMIETAB
************************************************************************
      use FASTMIEDIM,ONLY: Nx,Nn,Nk,nmin,nmax,kmin,kmax,xmin,xmax,
     >                     lQsc,lQab,lQex
      use DATATYPE,ONLY: r2            ! from MIEX
      use MIE_ROUTINES,ONLY: SHEXQNN2  ! from MIEX
      implicit none
      integer :: i,j,k
      real :: nnn,kkk,lam,aa,Qext,Qsca,Qabs,g
      !------------ variables for data exchange with MIEX ---------
      integer,parameter :: nang=2
      complex(kind=r2) :: ri
      real(kind=r2)    :: xxx,Qext2,Qsca2,Qabs2,Qbk2,Qpr2,albedo,g2
      integer          :: ier
      complex(kind=r2),dimension(nang) :: SA1,SA2

      open(unit=12,file='fastmie.dat',form="unformatted",
     >     status='replace')
      write(12) Nx,Nn,Nk
      write(12) xmin,xmax
      write(12) nmin,nmax
      write(12) kmin,kmax
      do i=0,Nx
        xxx = EXP(LOG(xmin)+i/REAL(Nx)*LOG(xmax/xmin))
!$omp parallel
!$omp& default(none)
!$omp& shared(i,xxx,lQsc,lQab,lQex)
!$omp& private(j,k,nnn,kkk,lam,aa,Qext,Qsca,Qabs,g)       
!$omp& private(ri,Qext2,Qsca2,Qabs2,Qbk2,Qpr2,albedo,g2,ier,SA1,SA2)       
!$omp do schedule(dynamic,1)
        do j=0,Nn
          nnn = EXP(LOG(nmin)+j/REAL(Nn)*LOG(nmax/nmin))
          do k=0,Nk
            kkk = EXP(LOG(kmin)+k/REAL(Nk)*LOG(kmax/kmin))
            ri = DCMPLX(nnn,kkk)
            call SHEXQNN2(ri,xxx,Qext2,Qsca2,Qabs2,Qbk2,Qpr2,
     >                    albedo,g2,ier,SA1,SA2,.false.,nang)
            Qext = Qext2
            Qabs = Qabs2
            Qsca = Qsca2
            !lam = 1.0
            !aa = xxx/6.2831853
            !call Q_MIE(nnn,kkk,lam,aa,Qext,Qsca,Qabs,g)            
            print'(4I4,1pE10.2)',i,j,k,ier,Qabs
            if (ier.ne.0) then
              print*,"*** error in MIEX, call Q_MIE ..."
              lam = 1.0
              aa = xxx/6.2831853
              call Q_MIE(nnn,kkk,lam,aa,Qext,Qsca,Qabs,g)
            endif
            lQex(i,j,k) = LOG(Qext)
            lQab(i,j,k) = LOG(Qabs)
            lQsc(i,j,k) = LOG(Qsca)
          enddo  
        enddo
!$omp end do      
!$omp end parallel        
      enddo
      write(12) lQsc,lQab,lQex
      close(12)
      end

************************************************************************
      subroutine FASTMIE(xin,nin,kin,Qext,Qsca,Qabs,debug)
************************************************************************
      use FASTMIEDIM,ONLY: Nx,Nn,Nk,nmin,nmax,kmin,kmax,xmin,xmax,
     >                     lQsc,lQab,lQex
      implicit none
      real,intent(in) :: xin,nin,kin
      real,intent(out) :: Qext,Qabs,Qsca
      logical,intent(in) :: debug
      integer :: i,j,k,idum(3),Nx_read,Nn_read,Nk_read
      real :: a1,a2,b1,b2,c1,c2
      real :: xval,nval,kval,xxx,kkk,nnn
      logical :: ex,match
      real :: n1_read,n2_read,k1_read,k2_read,x1_read,x2_read
      logical,save :: firstCall=.true.

      if (firstCall) then
        inquire(file='fastmie.dat',exist=ex)
 100    if (.not.ex) call FASTMIETAB
        open(unit=12,file='fastmie.dat',form="unformatted",status='old')
        read(12) Nx_read,Nn_read,Nk_read
        read(12) x1_read,x2_read
        read(12) n1_read,n2_read
        read(12) k1_read,k2_read
        match = (Nx_read==Nx).and.(Nn_read==Nn).and.(Nk_read==Nk).and.
     >          (x1_read==xmin).and.(x2_read==xmax).and. 
     >          (n1_read==nmin).and.(n2_read==nmax).and. 
     >          (k1_read==kmin).and.(k2_read==kmax)
        if (.not.match) then
          close(12)
          ex = .false.
          goto 100
        endif
        print*,"reading fastmie.dat ..."
        read(12) lQsc,lQab,lQex
        close(12)
        firstCall = .false.
        print*,"... done reading."
      endif

      xval = xin
      xval = MIN(xmax,xin)    ! don't extrapolate to larger x
      nval = nin
      kval = kin
      !if (kin<kmin) then
      !  kval = xval*kin/xin
      !endif
      xxx = LOG(xval/xmin)/LOG(xmax/xmin)*Nx
      nnn = LOG(nval/nmin)/LOG(nmax/nmin)*Nn
      kkk = LOG(kval/kmin)/LOG(kmax/kmin)*Nk
      i = INT(xxx)
      j = INT(nnn)
      k = INT(kkk)
      i = MAX(0,MIN(Nx-1,i))
      j = MAX(0,MIN(Nn-1,j))
      k = MAX(0,MIN(Nk-1,k))
      a1 = xxx-i
      a2 = 1.0-a1
      b1 = nnn-j
      b2 = 1.0-b1
      c1 = kkk-k
      c2 = 1.0-c1
      if (debug) then
        print'(1pE16.8,2I4,1pE16.8)',xxx,i,Nx,a1
        print'(3(1pE12.4))',EXP(LOG(xmin)+i/REAL(Nx)*LOG(xmax/xmin)),
     >        xval,EXP(LOG(xmin)+(i+1)/REAL(Nx)*LOG(xmax/xmin))
        print'(1pE16.8,2I4,1pE16.8)',nnn,j,Nn,b1
        print'(3(1pE12.4))',EXP(LOG(nmin)+j/REAL(Nn)*LOG(nmax/nmin)),
     >        nval,EXP(LOG(nmin)+(j+1)/REAL(Nn)*LOG(nmax/nmin))
        print'(1pE16.8,2I4,1pE16.8)',kkk,k,Nk,c1
        print'(3(1pE12.4))',EXP(LOG(kmin)+k/REAL(Nk)*LOG(kmax/kmin)),
     >        kval,EXP(LOG(kmin)+(k+1)/REAL(Nk)*LOG(kmax/kmin))
        print*,EXP(lQab(i,j,k)),EXP(lQsc(i,j,k))
      endif
      Qsca = EXP( lQsc(i  ,j,  k  )*a2*b2*c2
     >           +lQsc(i  ,j,  k+1)*a2*b2*c1
     >           +lQsc(i  ,j+1,k  )*a2*b1*c2
     >           +lQsc(i  ,j+1,k+1)*a2*b1*c1
     >           +lQsc(i+1,j  ,k  )*a1*b2*c2
     >           +lQsc(i+1,j  ,k+1)*a1*b2*c1
     >           +lQsc(i+1,j+1,k  )*a1*b1*c2
     >           +lQsc(i+1,j+1,k+1)*a1*b1*c1 )
      Qabs = EXP( lQab(i  ,j,  k  )*a2*b2*c2
     >           +lQab(i  ,j,  k+1)*a2*b2*c1
     >           +lQab(i  ,j+1,k  )*a2*b1*c2
     >           +lQab(i  ,j+1,k+1)*a2*b1*c1
     >           +lQab(i+1,j  ,k  )*a1*b2*c2
     >           +lQab(i+1,j  ,k+1)*a1*b2*c1
     >           +lQab(i+1,j+1,k  )*a1*b1*c2
     >           +lQab(i+1,j+1,k+1)*a1*b1*c1 )
      Qext = EXP( lQex(i  ,j,  k  )*a2*b2*c2
     >           +lQex(i  ,j,  k+1)*a2*b2*c1
     >           +lQex(i  ,j+1,k  )*a2*b1*c2
     >           +lQex(i  ,j+1,k+1)*a2*b1*c1
     >           +lQex(i+1,j  ,k  )*a1*b2*c2
     >           +lQex(i+1,j  ,k+1)*a1*b2*c1
     >           +lQex(i+1,j+1,k  )*a1*b1*c2
     >           +lQex(i+1,j+1,k+1)*a1*b1*c1 )
      if (kin<kmin.and.Qsca>Qabs) Qsca=Qext-Qabs
      if (kin<kmin.and.Qabs>Qsca) Qabs=Qext-Qsca
      end
