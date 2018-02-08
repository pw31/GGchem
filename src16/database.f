**********************************************************************
      MODULE DATABASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUSTmax,dust_nam
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: DMAX = 2*10**5
      integer :: NDAT=0,NLAST=0,NMODI=0,NPICK1=0,NPICK2=0
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
      use dust_data,ONLY: NELEM,NDUST,dust_nam
      use DATABASE,ONLY: NDAT,NLAST,dbase
      implicit none
      integer :: i
      character(len=80) :: filename="database.dat"
      if (NLAST==0) then
        open(unit=11,file=filename,form='unformatted',status='replace')
        write(11) NELEM,NDUST
        write(11) dust_nam
        do i=1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      else if (NDAT>NLAST) then 
        open(unit=11,file=filename,form='unformatted',position='append')
        do i=NLAST+1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      endif  
      NLAST = NDAT
      end

**********************************************************************
      SUBROUTINE LOAD_DBASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST,dust_nam
      use DATABASE,ONLY: qp,NDAT,NLAST,dbase
      implicit none
      integer :: i,NELEM_read,NDUST_read
      logical :: ex
      character(len=20) :: dust_nam_read(NDUST)
      character(len=80) :: filename="database.dat"

      NDAT = 0
      NLAST = 0
      inquire(file=filename,exist=ex)
      if (.not.ex) goto 200
      open(unit=11,file=filename,form="unformatted",status="old")
      read(11) NELEM_read,NDUST_read
      if (NELEM_read.ne.NELEM) goto 200
      if (NDUST_read.ne.NDUST) goto 200
      read(11) dust_nam_read
      do i=1,NDUST
        if (dust_nam(i).ne.dust_nam_read(i)) goto 200
      enddo
      do i=1,999999
        read(11,end=100) dbase(i)%ln 
        read(11) dbase(i)%lT
        read(11) dbase(i)%eps
        read(11) dbase(i)%ddust(1:NDUST)
        NDAT = NDAT+1
        !print*,i,EXP(dbase(i)%ln),EXP(dbase(i)%lT)
      enddo 
 100  close(11)
      print*,"... having read ",NDAT," datasets." 
      NLAST = NDAT
      return
 200  close(11)
      print*,"... no / unsuitable database."
      end

**********************************************************************
      SUBROUTINE PUT_DATA(nH,T,eps,ddust,qbest,ibest,active)
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST
      use DATABASE,ONLY: qp,NDAT,NLAST,NMODI,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T,qbest
      integer,intent(in) :: ibest
      real(kind=qp),intent(in) :: eps(NELEM),ddust(NDUST)
      logical,intent(in) :: active(0:NDUST)
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
      NMODI = i
      if (NDAT>NLAST+10) then
        call SAVE_DBASE
        print*,"... saved ",NDAT," datasets."
      endif  
      end


**********************************************************************
      subroutine GET_DATA(nH,T,eps,ddust,qbest,ibest,active)
**********************************************************************
      use dust_data,ONLY: NEPS,NELEM,NDUST,eps0,elnam,elnr,
     >                    dust_nel,dust_nu,dust_el,dust_nam
      use DATABASE,ONLY: qp,NDAT,NMODI,NPICK1,NPICK2,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T
      real*8,intent(out) :: qbest
      integer,intent(out) :: ibest
      real(kind=qp),intent(inout) :: eps(NELEM),ddust(NDUST)
      logical,intent(out) :: active(0:NDUST)
      real*8 :: ln,lT,lnread,lTread,qual,pot,rsort(NEPS)
      real(kind=qp) :: check(NELEM),error,errmax,corr,emain,del
      real(kind=qp) :: stoich(NEPS,NEPS),xx(NEPS),rest(NEPS),tmp
      integer :: i,j,k,it,el,elworst,b,bb,Nbuf,iloop
      integer :: isort(NEPS),jmain(NELEM),ibuf(NELEM)
      character(len=1) :: char
      character(len=80) :: frmt
      logical :: found,used(NDUST)
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        call LOAD_DBASE 
        firstCall = .false.
      endif  

      write(*,'("looking for nH,T=",2(1pE13.5)," ...")') nH,T
      ln = LOG(nH)
      lT = LOG(T) 
      qbest  = 9.d+99
      ibest  = 0
      pot    = -0.03
      !--- try last entry modified first ---
      if (NMODI>0) then
        i=NMODI
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        qbest = qual
        ibest = i
        if (qbest<1.d-3) goto 100
      endif  
      !--- try around entry picked last time ---  
      do i=MAX(1,NPICK1-1),MIN(NDAT,NPICK1+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      do i=MAX(1,NPICK2-1),MIN(NDAT,NPICK2+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      write(*,*) "entering full search ..."
      !--- check them all ---  
      do i=NDAT,1,-1
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
 100  active = .false.
      if (ibest>0) then
        eps    = dbase(ibest)%eps
        ddust  = dbase(ibest)%ddust(1:NDUST)
        do i=1,NDUST
          if (ddust(i)>0.Q0) active(i)=.true.
        enddo
        NPICK2 = NPICK1
        NPICK1 = ibest
        write(*,'(" ... found best dataset (",I6,
     >          ")  nH,T,qual=",3(1pE13.5))')
     >     ibest,EXP(dbase(ibest)%ln),EXP(dbase(ibest)%lT),qbest

        !----------------------------------------------------
        ! ***  adapt eps and ddust to reach current eps0  ***
        !----------------------------------------------------
        !--- 0. a few direct corrections ---
        do it=1,10
          check = eps
          do i=1,NDUST
            do j=1,dust_nel(i)
              el = dust_el(i,j)
              check(el) = check(el) + ddust(i)*dust_nu(i,j)    
            enddo
          enddo
          errmax = -1.Q0
          do el=1,NELEM
            error = ABS(1.Q0-check(el)/eps0(el))
            if (error.gt.errmax) then
              errmax = error
              elworst = el
              corr = eps0(el)/check(el)
            endif   
          enddo
          !print*,elnam(elworst),errmax,corr
          if (errmax<1.Q-25) return          ! perfect fit - nothing to do
          if (errmax<0.01) exit
          el = elworst
          eps(el) = eps(el)*corr
          do i=1,NDUST
            if (ddust(i)==0.Q0) cycle 
            do j=1,dust_nel(i)
              if (el==dust_el(i,j)) then
                ddust(i) = ddust(i)*corr
                exit
              endif
            enddo
          enddo  
        enddo  
        !--- 1. sort elements ---
        rsort = 9.d+99
        isort = 0
        do i=1,NEPS
          el = elnr(i) 
          do j=NEPS+1,2,-1
            if (eps0(el)>rsort(j-1)) exit
          enddo  
          isort(j+1:NEPS) = isort(j:NEPS-1)
          rsort(j+1:NEPS) = rsort(j:NEPS-1)
          isort(j) = el
          rsort(j) = eps0(el)
        enddo
        iloop = 1
 200    continue
        !--- 2. identify main reservoirs ---
        used = .false.
        Nbuf = 0
        do i=1,NEPS
          el = isort(i)
          check(el) = eps(el)
          emain = eps(el)
          jmain(el) = 0
          do j=1,NDUST
            if (ddust(j)==0.Q0) cycle
            do k=1,dust_nel(j)
              if (dust_el(j,k).ne.el) cycle
              del = ddust(j)*dust_nu(j,k)    
              check(el) = check(el) + del
              if (.not.used(j).and.del>emain) then
                emain = del
                jmain(el) = j
              endif
            enddo
          enddo  
          if (jmain(el)>0) then
            j = jmain(el) 
            used(j)=.true. 
            Nbuf = Nbuf+1
            ibuf(Nbuf) = el
            !print*,elnam(el)//" "//trim(dust_nam(j)),REAL(ddust(j))
          endif  
        enddo  
        !--- 3. setup linear equation system ---
        stoich = 0.Q0
        do b=1,Nbuf
          el = ibuf(b) 
          do bb=1,Nbuf
            j = jmain(ibuf(bb))
            do k=1,dust_nel(j)
              if (dust_el(j,k)==el) then
                 stoich(b,bb) = dust_nu(j,k) 
              endif
            enddo
          enddo
        enddo  
        check = eps0 - eps
        do j=1,NDUST
          if (ddust(j)==0.Q0) cycle
          if (used(j)) cycle 
          do k=1,dust_nel(j)
            el = dust_el(j,k) 
            check(el) = check(el) - ddust(j)*dust_nu(j,k)
          enddo
        enddo  
        do b=1,Nbuf
          el = ibuf(b) 
          j = jmain(el)
          rest(b) = check(el)
          write(frmt,'("(A2,2x,",I2,"(I2),A16,2(1pE13.6))")') Nbuf
          write(*,frmt) elnam(el),INT(stoich(b,1:Nbuf)),
     &           trim(dust_nam(j)),REAL(ddust(j)),REAL(rest(b))
        enddo  
        call GAUSS16( NEPS, Nbuf, stoich, xx, rest)
        do b=1,Nbuf
          el = ibuf(b)
          j = jmain(el)
          ddust(j) = xx(b)
          print'(A3,A16,1pE13.6)',elnam(el),trim(dust_nam(j)),ddust(j)
          if (xx(b)<0.Q0) then
            print*,"*** negative dust abundance in database.f"
            ddust(j) = 0.Q0
            active(j) = .false.
            qbest = 9.d+99
            return
          endif  
        enddo  
        !--- 4. correct gas element abundances ---
        check = 0.Q0
        do i=1,NDUST
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        do i=1,NEPS
          el = isort(i)
          if (jmain(el)==0) then
            tmp = eps0(el)-check(el)
            print*,elnam(el)//" gas",REAL(tmp)
            if (tmp<0.Q0) then
              print*,"*** negative element abundance in database.f" 
              qbest = 9.d+99
              return
              !if (jmain(el)==0) then
              !  isort(i) = isort(iloop) 
              !  isort(iloop) = el 
              !  iloop = iloop+1
              !  if (iloop>NEPS) iloop=1
              !else   
              !  j = jmain(el)
              !  ddust(j) = 0.Q0
              !  active(j) = .false.
              !endif  
              !goto 200
            endif  
            eps(el) = tmp
          endif
        enddo
        check = eps
        do i=1,NDUST
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        errmax = -1.Q0
        do el=1,NELEM
          error = ABS(1.Q0-check(el)/eps0(el))
          if (error.gt.errmax) then
            errmax = error
            elworst = el
          endif   
        enddo  
        !print*,"check ",elnam(elworst),errmax
        if (errmax>1.Q-10) then
          print*,"*** element conservation violation in database.f"
          stop
        endif  
      endif
  
      end
