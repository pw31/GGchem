**********************************************************************
      MODULE DATABASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUSTmax,eps0,dust_nam
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: DMAX = 10**5
      integer :: NDAT=0
      TYPE ENTRY
        real*8 :: ln
        real*8 :: lT
        real(kind=qp) :: eps(NELEM)
        real(kind=qp) :: ddust(NDUSTmax)
      END TYPE ENTRY
      TYPE(ENTRY) :: dbase(DMAX)
      end MODULE DATABASE

**********************************************************************
      SUBROUTINE SAVE_DBASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST,eps0,dust_nam
      use DATABASE,ONLY: NDAT,dbase
      implicit none
      integer :: i
      character(len=80) :: filename="database.dat"
      open(unit=11,file=filename,form='unformatted',status='replace')
      write(11) NELEM,NDUST,NDAT
      write(11) eps0
      write(11) dust_nam
      do i=1,NDAT
        write(11) dbase(i)%ln 
        write(11) dbase(i)%lT
        write(11) dbase(i)%eps
        write(11) dbase(i)%ddust(1:NDUST)
      enddo 
      close(11)
      end

**********************************************************************
      SUBROUTINE LOAD_DBASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST,eps0,dust_nam
      use DATABASE,ONLY: qp,NDAT,dbase
      implicit none
      integer :: i,NELEM_read,NDUST_read
      logical :: ex
      real(kind=qp) :: eps0_read(NELEM)
      character(len=15) :: dust_nam_read(NDUST)
      character(len=80) :: filename="database.dat"

      NDAT = 0
      inquire(file=filename,exist=ex)
      if (.not.ex) goto 100
      open(unit=11,file=filename,form="unformatted",status="old")
      read(11) NELEM_read,NDUST_read,NDAT
      if (NELEM_read.ne.NELEM) goto 100
      if (NDUST_read.ne.NDUST) goto 100
      read(11) eps0_read
      do i=1,NELEM
        if (eps0(i).ne.eps0_read(i)) goto 100
      enddo
      read(11) dust_nam_read
      do i=1,NDUST
        if (dust_nam(i).ne.dust_nam_read(i)) goto 100
      enddo
      do i=1,NDAT
        read(11) dbase(i)%ln 
        read(11) dbase(i)%lT
        read(11) dbase(i)%eps
        read(11) dbase(i)%ddust(1:NDUST)
        !print*,i,EXP(dbase(i)%ln),EXP(dbase(i)%lT)
      enddo 
      close(11)
      print*,"... having read ",NDAT," datasets." 
      return
 100  close(11)
      NDAT = 0
      print*,"... no / unsuitable database."
      end

**********************************************************************
      SUBROUTINE PUT_DATA(nH,T,eps,ddust,qbest,ibest,active)
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST,eps0,dust_nam
      use DATABASE,ONLY: qp,NDAT,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T,qbest
      integer,intent(in) :: ibest
      real(kind=qp),intent(in) :: eps(NELEM),ddust(NDUST)
      logical,intent(in) :: active(0:NDUST)
      integer,save :: Nlast=0
      integer :: i,j
      
      if (qbest<1.d-8) then
        return 
      else if (qbest<1.d-3) then
        i = ibest
        write(*,'(" ... replacing database entry (",I6,
     >          ") nH,T=",2(1pE15.7))') i,nH,T
      else  
        NDAT = NDAT+1
        i = NDAT
        write(*,'(" ... adding database entry (",I6,
     >          ") nH,T=",2(1pE15.7))') i,nH,T
        if (NDAT>DMAX) then
          print*,"*** NDAT>DMAX in PUT_DATA",NDAT,DMAX
          stop
        endif  
      endif  
      dbase(i)%ln = LOG(nH)
      dbase(i)%lT = LOG(T)
      dbase(i)%eps = eps
      do j=1,NDUST
        dbase(i)%ddust(j) = ddust(j)
        if (.not.active(j)) dbase(i)%ddust(j)=0.Q0
      enddo  
      if (NDAT>Nlast+10) then
        call SAVE_DBASE
        print*,"... saved ",NDAT," datasets."
        Nlast = NDAT
      endif  
      end


**********************************************************************
      subroutine GET_DATA(nH,T,eps,ddust,qbest,ibest,active)
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST
      use DATABASE,ONLY: qp,NDAT,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T
      real*8,intent(out) :: qbest
      integer,intent(out) :: ibest
      real(kind=qp),intent(inout) :: eps(NELEM),ddust(NDUST)
      real*8 :: ln,lT,lnread,lTread,qual,pot
      logical,intent(out) :: active(0:NDUST)
      integer :: i
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        call LOAD_DBASE 
        firstCall = .false.
      endif  

      ln = LOG(nH)
      lT = LOG(T) 
      qbest  = 9.d+99
      ibest  = 0
      pot    = -0.03
      do i=NDAT,1,-1
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.1*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) exit
        endif  
      enddo
      active = .false.
      if (ibest>0) then
        eps    = dbase(ibest)%eps
        ddust  = dbase(ibest)%ddust(1:NDUST)
        do i=1,NDUST
          if (ddust(i)>0.Q0) active(i)=.true.
        enddo

        write(*,'(" ... found best dataset (",I6,
     >          ")  nH,T,qual=",3(1pE13.5))')
     >     ibest,EXP(dbase(ibest)%ln),EXP(dbase(ibest)%lT),qbest

      endif  
      end
