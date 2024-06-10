***********************************************************************
      PROGRAM CREATE_BURCAT
***********************************************************************
      use CHEMISTRY,ONLY: NMOLE,fit,cmol,catm,m_kind,m_anz,Natom,a,
     >                    th1,th2,th3,th4,TT1,TT2,TT3
      implicit none
      character(len=200) :: line,frmt
      character(len=2) :: atnam(80),atm
      real :: aatm(80,14)
      real,parameter :: bar=1.d+6,lnbar=LOG(bar)
      real,parameter :: Rgas=8.3144598d+0  ! J/K/mol
      integer,parameter :: Ndat=300
      integer :: it,i,j,k,at,Nat,Nit,Nel,jT1000,NX,Nstart
      integer :: Tfit,atnum(5000,10)
      logical :: has_atoms(5000)
      real,parameter :: Tmin=100.0, Tmax=6000.0
      real,dimension(2000) :: Tdat,Gdat,XX,YY
      real,external :: gk
      real :: lnk,dGRT,dGRTatom,dG,stoich(5000,10)
      real :: afit(7),eps0,ldum,delta,dev,qual
      
      !--- initialise GGchem molecular data ---
      call READ_PARAMETER
      call INIT
      call INIT_CHEMISTRY

      !--- read BURCAT atom fit coeffs wrt their standard state ---
      print*
      Nat = 1
      open(1,file='../Burcat_atoms.dat',status='old')
      do
        do 
          read(1,'(A200)',end=1000) line
          if (line(80:80)=="1") exit
        enddo
        atnam(Nat) = line(1:2)
        read(1,'(5(E15.8))') aatm(Nat,1:5)
        read(1,'(5(E15.8))') aatm(Nat,6:10)
        read(1,'(5(E15.8))') aatm(Nat,11:14)
        print*,"found BURCAT atom "//atnam(Nat)
        !print*,aatm(Natom,1:14)
        Nat = Nat+1
      enddo  
 1000 continue
      close(1)
      print*
      
      !--- create T-array for fiting ---
      do j=1,Ndat-1
        Tdat(j) = EXP(LOG(Tmin)+LOG(Tmax/Tmin)*(j-1.0)/(Ndat-2.0))
        if (Tdat(j)<1000.0) jT1000=j 
      enddo
      Tdat(jT1000+2:Ndat) = Tdat(jT1000+1:Ndat-1)
      Tdat(jT1000+1) = 1000.0
      
      !--- identify atoms in NASA data ---
      has_atoms = .true.
      do i=1,NMOLE
        !print'(I4,1x,A10,I2,I2,1x,99(1pE13.5))',
     >  !       i,cmol(i),Natom(i),fit(i),a(i,0:4)
        do k=1,m_kind(0,i)
          atm = catm(m_kind(k,i))
          do at=1,Nat
            if (atnam(at)==atm) exit
          enddo
          if (at>Nat) then
            print*,trim(cmol(i))//": "//atm//" BURCAT atom not found."
            has_atoms(i) = .false.
          else
            atnum(i,k) = at
            stoich(i,k) = m_anz(k,i)
            !print'(2(A3),I3)',atm,atnam(at),INT(stoich(k))
          endif
        enddo
      enddo
      
      !--- main loop ---
      open(unit=1,file="BURCAT_generated.dat")
      do i=1,NMOLE
        if (.not.has_atoms(i)) cycle
        print'(I4,1x,A10,I2,I2,1x,99(1pE13.5))',
     >         i,cmol(i),Natom(i),fit(i),a(i,0:4)
        do j=1,Ndat
          !--- compute dG/RT from GGchem's equil.const. ---
          TT1 = Tdat(j)
          TT2 = TT1*TT1
          TT3 = TT2*TT1
          th1 = 5040.d0/TT1
          th2 = th1*th1
          th3 = th2*th1
          th4 = th3*th1
          lnk = gk(i)                      ! mol. equilibrium constant
          dGRT = (1-Natom(i))*lnbar - lnk  ! dG/RT with respect to atoms
          !--- subtract NASA-atoms to get dGRT wrt BURCAT standard states ---
          do k=1,m_kind(0,i)
            dGRTatom = NASAatom(atnum(i,k))
            dGRT = dGRT - stoich(i,k)*dGRTatom
            !print*,atnam(atnum(i,k)),stoich(i,k),-dGRTatom*Rgas*TT1/1000
          enddo
          Gdat(j) = -dGRT
          dG = dGRT*Rgas*TT1/1000          ! dG [kJ/mol]
          !print*,TT1,dG
        enddo
        !if (trim(cmol(i))=='S2') then
        !afit(1:7) = (/ 3.83249656E+00, 8.88970881E-04,-2.59080844E-07,
     >  !               3.63847115E-11,-1.72606371E-15, 1.42836134E+04,
     >  !               5.33000845E+00/)
        afit(1:7) = (/ 2E+0, 1E-4, 1E-7, 1E-11, 1E-14, 1E+4, 1E+0 /)
        do Tfit=1,2
          if (Tfit==1) then   ! fit T>1000
            Nstart = jT1000-1
            NX = Ndat-Nstart+1
            XX(1:NX) = Tdat(Nstart:Ndat)
            YY(1:NX) = Gdat(Nstart:Ndat)
          else
            NX = jT1000+1
            XX(1:NX) = Tdat(1:NX)
            YY(1:NX) = Gdat(1:NX)
          endif
          call PARAFIND(NASAfit,7,NX,afit,XX,YY)
          delta = 1.E-6
          eps0  = 0.2d0
          do it=1,50
            Nit = 0
            call PARAM_LS(NASAfit,7,Nit,NX,dev,delta,eps0,afit,XX,YY)
            !print'("coefficients:",I4,7(1pE15.7),I7,1pE11.3)',
     >      !   it,afit(1:7),Nit,dev
          enddo
          if (Tfit==1) then   ! fit T>1000
            a(i,0:6) = afit(1:7)
          else
            a(i,7:13) = afit(1:7)
          endif
        enddo
        qual = 0.0
        do j=1,Ndat
          TT1 = Tdat(j)
          if (TT1>1000.0) then
            afit(1:7) = a(i,0:6)
          else
            afit(1:7) = a(i,7:13)
          endif
          dGRT = NASAfit(7,afit,TT1)
          qual = qual + (dGRT-Gdat(j))**2
          !print'(I4,0pF10.3,2(1pE12.4),1pE10.2)',
     >    !     j,Tdat(j),Gdat(j),dGRT,dGRT-Gdat(j)
        enddo
        qual = SQRT(qual/Ndat)   ! fit quality dG/RT
        Nel = m_kind(0,i)
        print*,i,cmol(i),qual
        if (qual>0.01) print*,"*** WARNING, fit is not great"
        write(frmt,'("(A18,I2,",I1,"(A3),",I1,"(I3),",I2,
     >               "x,0pF10.6)")') 
     >        Nel,Nel,6*(4-Nel)+1
        write(1,frmt) cmol(i),Nel,catm(m_kind(1:Nel,i)),
     >                int(m_anz(1:Nel,i)),qual
        write(1,'(14(1pE15.7))') a(i,0:13)
        !endif
      enddo

      contains

***************************************************************
      REAL FUNCTION NASAfit(Npara,afit,T)
***************************************************************
      implicit none
      integer,intent(in) :: Npara
      real*8,intent(in)  :: afit(Npara),T
      NASAfit = afit(1)*(log(T)-1.d0) 
     &        + afit(2)*T   / 2.d0
     &        + afit(3)*T**2/ 6.d0    
     &        + afit(4)*T**3/12.d0
     &        + afit(5)*T**4/20.d0    
     &        - afit(6)/T  
     &        + afit(7)
      end function NASAfit
      
***************************************************************
      real function NASAatom(i)
***************************************************************
      implicit none
      integer,intent(in) :: i
      real :: fit(7)
      if (TT1>1.d3) then
        fit = aatm(i,1:7)
      else
        fit = aatm(i,8:14)
      endif
      NASAatom = NASAfit(7,fit,TT1)
      end function NASAatom

      end

      
***************************************************************
      SUBROUTINE S200(func,Npara,l2,Ndat,d,A,X,Y) 
***************************************************************
      implicit none
      integer,parameter :: SIZE=2000
      integer,intent(in) :: Npara,Ndat
      real,intent(in) :: A(Npara),X(SIZE),Y(SIZE) 
      real,intent(out) :: d,l2
      real :: xx,yy,weight
      real,external :: func
      integer :: j
      !print*,Ndat
      !print*,X(1:Ndat)
      !print*,Y(1:Ndat)
      !print*,Npara,A
      l2 = 0.d0
      do j = 1,Ndat
         xx = X(j)
         yy = func(Npara,A,xx)
         weight = 1.0
         if (ABS(xx-1000.0)<0.01) weight=0.2*Ndat
         l2 = l2 + weight*(Y(j)-yy)**2
      enddo
      d = dsqrt(l2/(Ndat-Npara))
      end

!***************************************************************
!* Parametric least squares curve fit subroutine. This program *
!* least squares fits a function to a set of data values by    *
!* successively reducing the variance. Convergence depends on  *
!* the initial values and is not assured.                      *
!* n pairs of data values, X(i), Y(i), are given. There are l  *
!* parameters, A(j), to be optimized across.                   *
!* Required are initial values for the A(l) and e. Another     *
!* important parameter which affects stability is e1, which is *
!* initially converted to e1(l) for the first intervals.       *
!* The parameters are multiplied by (1 - e1(i)) on each pass.  *
!***************************************************************
      SUBROUTINE PARAM_LS(func,l,m,n,d,e,ee1,A,X,Y)  
      implicit none
      integer,parameter   :: SIZE=2000
      integer,intent(in)  :: l,n
      integer,intent(out) :: m
      real*8,intent(in)   :: e,ee1,X(SIZE),Y(SIZE)
      real*8,intent(inout)  :: d,A(l)
      real*8,external :: func
      real*8  :: a0,l1,l2,m0,m1
      real*8  :: E1(l)
      integer :: i
      do i = 1,l
         E1(i) = ee1
      end do	
      !Set up test residual
      l1 = 1.d+99
      !Make sweep through all parameters
 50   do i = 1,l
         a0 = A(i)
         !Get value of residual
         A(i) = a0
 100     call S200(func,l,l2,n,d,A,X,Y)
         !Store result in m0
         m0 = l2
         !Repeat for m1
         A(i) = a0 * (1.d0 - E1(i))
         call S200(func,l,l2,n,d,A,X,Y)
         m1 = l2
         !Change interval size if called for
         !If variance was increased, halve E1(i) 
         if (m1 > m0)  then
            E1(i) = -E1(i) / 2.d0
         end if
         !If variance was reduced, increase step size by increasing E1(i)
         if (m1 < m0)  then
            E1(i) = 1.2d0 * E1(i)
         end if 
         !If variance was increased, try to reduce it
         if (m1 > m0) then  
            A(i) = a0
         end if 
         if (m1 > m0) then
            goto 100
         end if
      enddo 
      !End of a complete pass
      !Test for convergence 
      m = m + 1
      if (l2.eq.0.d0) then
         return
      end if 
      if (dabs((l1-l2)/l2)<e)  then
         return
      end if 
      !If this point is reached, another pass is called for
      l1 = l2
      goto 50
      end

      SUBROUTINE PARAFIND(func,Npara,Ndat,A,X,Y)  
      implicit none
      integer,parameter   :: SIZE=2000
      integer,intent(in)  :: Npara,Ndat
      real,intent(in)     :: X(SIZE),Y(SIZE)
      real,intent(inout)  :: A(Npara)
      real,external :: func
      real :: dum,qfit,q0,qold,Asave,Abest(Npara)
      real :: dir(Npara),norm,rvar,eps,dqdp
      integer :: i,it,steps,ichange,irange
      integer :: range(7)=(/6,1,2,3,4,5,7/)

      call S200(func,Npara,dum,Ndat,q0,A,X,Y)
      !print'("preconditioning",I4,8(1pE10.3))',0,A(1:7),q0
      rvar = 1.E-7
      do it=1,100
        do irange=1,7
          i = range(irange)
          call S200(func,Npara,dum,Ndat,q0,A,X,Y)
          Asave = A(i)
          A(i) = Asave*(1.0+rvar)
          call S200(func,Npara,dum,Ndat,qfit,A,X,Y)
          A(i) = Asave
          dqdp = (qfit-q0)/(rvar*Asave)
          dir(:) = 0.0
          dir(i) = -q0/dqdp
          if (dir(i)/A(i)>+1.0) dir=+A(i)
          if (dir(i)/A(i)<-1.0) dir=-A(i)
          !print*,i,A(i),dir(i)
          eps = 0.2
          ichange = 0
          qold = q0
          do steps=1,100
            A(i) = A(i) + eps*dir(i)
            call S200(func,Npara,dum,Ndat,qfit,A,X,Y)
            if (qfit<qold) then
              Abest(i) = A(i)
            else
              eps=-0.6*eps
              ichange=ichange+1
            endif
            !print'(3(I4),4(1pE14.6))',it,i,steps,A(i),eps*dir(i),qfit
            if (ABS(eps)<1.E-5) exit
            qold=qfit
          enddo
          A(i) = Abest(i)
        enddo
        !print'("preconditioning",I4,8(1pE10.3))',it,A(1:7),qfit
      enddo
      return
      eps  = 0.01
      do it=1,20
        call S200(func,Npara,dum,Ndat,q0,A,X,Y)
        norm = 0.0
        do i=1,Npara
          Asave = A(i)
          A(i) = Asave*(1.0+rvar)
          call S200(func,Npara,dum,Ndat,qfit,A,X,Y)
          A(i) = Asave
          dqdp = (qfit-q0)/(rvar*Asave)
          dir(i) = -q0/dqdp
          print'(I4,3(1pE14.6))',i,A(i),qfit-q0,dir(i)
        enddo
        eps = ABS(eps)*2.0
        qold = q0
        print*,it,0,0.0,q0
        ichange = 0
        do steps=1,999
          A = A + eps*dir
          call S200(func,Npara,dum,Ndat,qfit,A,X,Y)
          if (qfit<qold) then
            Abest = A
          else
            eps=-0.6*eps
            ichange=ichange+1
          endif
          print*,it,steps,eps,qfit
          if (qfit<q0.and.qfit<qold.and.ichange==2) exit
          if (ichange>5) exit
          qold=qfit
        enddo
        A = Abest
        print'("coefficients after  :",7(1pE15.7),I7,1pE11.3)',
     >         A(1:7),it,qfit
        if (steps>10) eps=2*eps
        if (ABS(eps)<1.E-7.and.it>3) exit
      enddo
      end



      
