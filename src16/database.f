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
        real*8 :: eprod
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
          write(11) dbase(i)%eprod
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      else if (NDAT>NLAST) then 
        open(unit=11,file=filename,form='unformatted',position='append')
        do i=NLAST+1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%eprod
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
        read(11) dbase(i)%eprod
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
      use dust_data,ONLY: NELEM,NDUST,NEPS,dust_nam,eps0,elnr
      use DATABASE,ONLY: qp,NDAT,NLAST,NMODI,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T,qbest
      integer,intent(in) :: ibest
      real*8 :: prod
      real(kind=qp),intent(in) :: eps(NELEM),ddust(NDUST)
      logical,intent(in) :: active(0:NDUST)
      integer :: i,j,e,el
      
      !if (qbest<1.d-8) then
      !  return 
      !else 
      if (qbest<1.d-3) then
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
      prod = 1.0
      do e=1,NEPS
        el = elnr(e) 
        prod = prod * eps0(el)
      enddo  
      dbase(i)%ln    = LOG(nH)
      dbase(i)%lT    = LOG(T)
      dbase(i)%eprod = LOG(prod)
      dbase(i)%eps   = eps
      do j=1,NDUST
        !if (ddust(j)>0.Q0) print*,active(j),dust_nam(j),real(ddust(j))
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
      subroutine GET_DATA(nH,T,eps,ddust,qbest,ibest,active,method)
**********************************************************************
      use dust_data,ONLY: NEPS,NELEM,NDUST,eps0,elnam,elnr,
     >                    dust_nel,dust_nu,dust_el,dust_nam
      use DATABASE,ONLY: qp,NDAT,NMODI,NPICK1,NPICK2,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T
      real*8,intent(out) :: qbest
      integer,intent(out) :: ibest,method
      real(kind=qp),intent(inout) :: eps(NELEM),ddust(NDUST)
      logical,intent(out) :: active(0:NDUST)
      real*8 :: prod,lp,ln,lT,lpread,lnread,lTread
      real*8 :: qual,pot,rsort(NEPS)
      real(kind=qp) :: check(NELEM),error,errmax,emain,del,sjk,sik
      real(kind=qp) :: stoich(NEPS,NEPS),xx(NEPS),rest(NEPS),tmp,val
      real(kind=qp) :: ecopy(NELEM),dcopy(NDUST),deps0(NELEM),echange
      integer :: i,j,k,it,el,el2,elworst,b,bb,Nbuf,iloop,jj,ii,Nact,e
      integer :: isort(NEPS),jmain(NELEM),ibuf(NELEM)
      character(len=1) :: char
      character(len=80) :: frmt
      character(len=999) :: condensates
      logical :: IS_NAN,found,contain1,contain2,corrected
      logical :: condensed(NELEM),used(NDUST)
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        call LOAD_DBASE 
        firstCall = .false.
      endif  

      prod = 1.0
      do e=1,NEPS
        el = elnr(e) 
        prod = prod * eps0(el)
      enddo  
      print'("looking for nH,T,eprod=",3(1pE13.5)," ...")',nH,T,prod
      ln = LOG(nH)
      lT = LOG(T) 
      lp = LOG(prod)
      qbest  = 9.d+99
      ibest  = 0
      pot    = -0.03
      !--- try last entry modified first ---
      if (NMODI>0) then
        i=NMODI
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 10.0*ABS(lpread-lp)
        qbest = qual
        ibest = i
        if (qbest<1.d-3) goto 100
      endif  
      !--- try around entry picked last time ---  
      do i=MAX(1,NPICK1-1),MIN(NDAT,NPICK1+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 10.0*ABS(lpread-lp)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      do i=MAX(1,NPICK2-1),MIN(NDAT,NPICK2+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 10.0*ABS(lpread-lp)
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
        lpread = dbase(i)%eprod
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 10.0*ABS(lpread-lp)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
 100  active = .false.
      condensates = ""
      if (ibest>0) then
        eps    = dbase(ibest)%eps
        ddust  = dbase(ibest)%ddust(1:NDUST)
        Nact = 0
        do i=1,NDUST
          if (ddust(i)>0.Q0) then
            Nact = Nact+1
            active(i)=.true.
            condensates = trim(condensates)//" "//trim(dust_nam(i))
          endif  
        enddo
        NPICK2 = NPICK1
        NPICK1 = ibest
        write(*,'(" ... found best dataset (",I6,
     >          ")  nH,T,eprod,qual=",4(1pE13.5))')
     >     ibest,EXP(dbase(ibest)%ln),EXP(dbase(ibest)%lT),
     >     EXP(dbase(ibest)%eprod),qbest
        write(*,*) "active condensates: "//trim(condensates)
        ecopy = eps
        dcopy = ddust

        condensed(:) = .false.
        check = eps
        do i=1,NDUST
          do j=1,dust_nel(i)
            if (ddust(i)==0.Q0) cycle
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
            condensed(el) = .true.
          enddo
        enddo
        deps0 = eps0-check
        errmax = -1.Q0
        do e=1,NEPS
          el = elnr(e) 
          error = ABS(1.Q0-check(el)/eps0(el))
          if (error.gt.errmax) then
            errmax = error
            elworst = el
          endif   
        enddo
        print*,"need fitting? "//elnam(elworst),SNGL(errmax)
        method = 0
        if (errmax<1.Q-25) return ! perfect fit - nothing to do

        !--------------------------------------------------------------
        ! ***  adapt eps and ddust to reach current eps0: method 1  ***
        !--------------------------------------------------------------
        print*,"adjust with method 1 ..."
        method = 1

        !--- 1. sort elements from rare to abundant ---
        rsort = 9.d+99
        isort = 0
        do i=1,NEPS
          el = elnr(i)
          val = eps(el)
          if (.not.condensed(el)) val=val+1.E+10  ! non-condensed to bottom
          do j=NEPS+1,2,-1
            if (val>rsort(j-1)) exit
          enddo  
          isort(j+1:NEPS) = isort(j:NEPS-1)
          rsort(j+1:NEPS) = rsort(j:NEPS-1)
          isort(j) = el
          rsort(j) = val
        enddo
        jmain = 0                                 ! list of condensates
        used  = .false.
        Nbuf  = 0
        do j=1,NDUST                              
          if (ddust(j)==0.Q0) cycle
          Nbuf = Nbuf+1          
          el = isort(Nbuf)
          ibuf(Nbuf) = el
          jmain(el) = j
          used(j)=.true.
          !print*,elnam(el),dust_nam(j),REAL(eps(el))
        enddo  

        !--------------------------------------------------------------
        !*** 2. setup linear equation system                        ***
        !*** the idea is to find modified ddust, such that eps is   ***
        !*** as in the database:  ((stoich))*(ddust) = (eps0-eps)   ***
        !*** ibuf(1...Nbuf) : sorted list of elements               ***
        !*** jmain(el) : dust index of main reservoir               ***
        !-------------------------------------------------------------- 
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
     &           trim(dust_nam(j)),ddust(j),rest(b)
        enddo  
        call GAUSS16( NEPS, Nbuf, stoich, xx, rest)
        do b=1,Nbuf
          el = ibuf(b)
          j = jmain(el)
          print'(A3,A16,1pE13.6," ->",1pE13.6)',
     &         elnam(el),trim(dust_nam(j)),ddust(j),xx(b)
          ddust(j) = xx(b)
          if (xx(b)<0.Q0.or.IS_NAN(DBLE(xx(b)))) then
            print*,"*** negative dust abund. in database.f method 1"
            goto 300
            !qbest = 9.d+99
            !return
          endif  
        enddo 
 
        !--- 3. correct remaining gas element abundances ---
        check = 0.Q0
        do i=1,NDUST
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        echange = 0.Q0
        do i=1,NEPS
          el = isort(i)
          if (jmain(el)==0) then
            tmp = eps0(el)-check(el)
            print'(A3,A16,1pE13.6," ->",1pE13.6)',
     &           elnam(el),"gas",eps(el),tmp
            if (tmp<0.Q0.or.IS_NAN(DBLE(tmp))) then
              print*,"*** negative el. abund. in database.f method 1" 
              goto 300
              !qbest = 9.d+99
              !return
            endif  
            eps(el) = tmp
            echange = MAX(echange,ABS(LOG10(eps(el)/ecopy(el))))
          endif
        enddo
        print*,'echange =',REAL(echange)
        if (echange<2.Q0) goto 500

 300    continue    
        !--------------------------------------------------------------
        ! ***  adapt eps and ddust to reach current eps0: method 2  ***
        !--------------------------------------------------------------
        print*,"adjust with method 2 ..."
        method = 2
        !--------------------------------------------------------------
        ! ***  the idea is to minimise the change of eps            ***
        ! ***      Sum_k (deps_k/eps_k)**2 -> Min                   ***
        ! ***  where deps_k = eps_k - eps_k(table)                  ***
        ! ***  with eps_k        = eps0_k        - Sum_j sjk xj     ***
        ! ***  and  eps_k(table) = eps0_k(table) - Sum_j sjk xj(table)
        ! ***       ---------------------------------------------------
        ! ***       deps_k       = deps0_k       - Sum_j sjk xx     ***
        ! ***       ---------------------------------------------------
        ! ***  the differentiate with respect to each xx_i = 0      ***
        ! -------------------------------------------------------------
        ddust = dcopy
        eps   = ecopy
        Nact = 0
        stoich = 0.Q0
        rest = 0.Q0
        i = 0
        do ii=1,NDUST
          if (ddust(ii)==0.Q0) cycle
          Nact = Nact+1
          i = i+1
          do k=1,NEPS
            el = elnr(k)
            do e=1,dust_nel(ii)
              if (el==dust_el(ii,e)) exit
            enddo  
            sik = dust_nu(ii,e)
            !print*,dust_nam(ii),elnam(el),sik,deps0(el)
            rest(i) = rest(i) + sik*deps0(el)/eps(el)**2
            j = 0
            do jj=1,NDUST
              if (ddust(jj)==0.Q0) cycle
              j = j+1
              do e=1,dust_nel(jj)
                if (el==dust_el(jj,e)) exit
              enddo  
              sjk = dust_nu(jj,e)
              stoich(i,j) = stoich(i,j) + sik*sjk/eps(el)**2 
            enddo  
          enddo  
        enddo
        !print*,NEPS,Nact,rest
        call GAUSS16( NEPS, Nact, stoich, xx, rest)
        eps = eps + deps0
        j = 0
        do jj=1,NDUST
          if (ddust(jj)==0.Q0) cycle
          j = j+1
          print'(A12,2(1pE12.4))',dust_nam(jj),dcopy(jj),xx(j)
          ddust(jj) = ddust(jj)+xx(j)
          do e=1,dust_nel(jj)
            el = dust_el(jj,e)
            eps(el) = eps(el)-dust_nu(jj,e)*xx(j)
          enddo
          if (ddust(jj)<0.Q0.or.IS_NAN(DBLE(ddust(jj)))) then
            method = 0
            print*,"*** negative dust abund. in database.f method 2"
            qbest = 9.d+99
            return
          endif
        enddo 
        do e=1,NEPS
          el = elnr(e)
          print'(A3,2(1pE12.4))',elnam(el),ecopy(el),eps(el)-ecopy(el)
          if (eps(el)<0.Q0.or.IS_NAN(DBLE(eps(el)))) then
            method = 0
            print*,"*** negative el. abund. in database.f method 2"
            qbest = 9.d+99
            return
          endif
        enddo

 500    continue
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
        if (errmax>1.Q-8) then
          print*,"*** element conservation violation in database.f"
          print*,elnam(elworst),errmax
          stop
        endif  
      endif

      end
