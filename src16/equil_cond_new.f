      module CONVERSION
      use DUST_DATA,ONLY: NELEM,NDUSTmax
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: Nind,Ndep,Iindex(NELEM),Dindex(NDUSTmax+NELEM)
      logical :: is_dust(NDUSTmax+NELEM)
      real(kind=qp) :: conv(NDUSTmax+NELEM,NELEM)
      end

!-------------------------------------------------------------------------
      SUBROUTINE EQUIL_COND(nHtot,T,eps,Sat,ddust,verbose)
!-------------------------------------------------------------------------
! ***  computes the gas element abundances eps, the saturation ratios  ***
! ***  Sat and the dust number densities ddust after equilibrium       ***
! ***  condensation. The routine determines the state in which each    ***
! ***  solid is either undersaturated or saturated S<=1 (all solids).  ***
! ***  The book-keeping is special: Although in principle redundant,   ***
! ***  eps(NELEM) and ddust(NDUST) are changed by dx(NELEM) separately ***
! ***  to avoid numerical problems close to complete condensation      ***
!-------------------------------------------------------------------------
      use PARAMETERS,ONLY: Tfast,useDatabase
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nam,dust_nel,dust_nu,dust_el,
     >                    eps0,elnam,elcode,NEPS,elnr
      use CHEMISTRY,ONLY: NewFastLevel,NELM,elnum,iel=>el
      use CONVERSION,ONLY: Nind,Ndep,Iindex,Dindex,is_dust,conv
      use EXCHANGE,ONLY: Fe,Mg,Si,Al,Ca,Ti,C,O,S,Na,Cl,H,Li,Mn,W,Ni,Cr,
     >                   Fluor=>F,Kalium=>K,Zr,V,
     >                   itransform,ieqcond,ieqconditer
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nHtot                ! H nuclei density [cm-3]
      real*8,intent(in) :: T                    ! temperature [K]
      real(kind=qp),intent(out) :: eps(NELEM)   ! gas element abundances
      real(kind=qp),intent(out) :: Sat(NDUST)   ! saturation ratio
      real(kind=qp),intent(out) :: ddust(NDUST) ! condensed units per <H> [-]
      integer,intent(inout) :: verbose
      real(kind=qp),dimension(NELEM) :: eps00,epsread,check,FF,Fsav,dx
      real(kind=qp),dimension(NELEM) :: eps_save,xstep,last_dx,null
      real(kind=qp),dimension(NELEM) :: scale,bvec
      real(kind=qp),dimension(NDUST) :: ddustread,dscale,pot,dust_save
      real(kind=qp),dimension(NDUST) :: Sat0,Sat1,Sat2,xvec,slin
      real(kind=qp),dimension(NELEM,NELEM) :: DF,DFsav
      real(kind=qp),dimension(NELEM,NDUST) :: AA
      real(kind=qp) :: worst,xx,xmin,Smax,Smin,qual,qold,SQUAL,del,abun
      real(kind=qp) :: turnon,maxon,minoff,fac,fac2,amount,Nt,dscalemax
      real(kind=qp) :: NRstep,dlim,target,dterm
      real(kind=qp) :: small=1.Q-30
      real(kind=qp) :: esum,emax,Qfinish,Sfinish
      integer,parameter :: itmax=5000
      integer,dimension(NELEM) :: eind
      integer,dimension(NDUST) :: dind,dlin,switchedON,switchedOFF
      integer :: it,i,j,k,el,elim,Nact,Nact_read,dk
      integer :: ii,jj,lastit,laston=0
      integer :: imaxon,iminoff,info,ipvt(NELEM)
      integer :: Nsolve,iread,ioff,method
      integer :: ifail,imax,irow,erow,Eact,Nlin,iback
      integer :: dust_to_act(NDUST),act_to_dust(NELEM)
      integer :: no_action
      logical,dimension(NELEM) :: e_act
      logical,dimension(0:NDUST) :: active,act_read,act_old
      logical,dimension(NDUST) :: itried
      logical :: action,changed,limited,ok
      logical :: limdust
      character(len=1) :: char1
      character(len=2) :: rem
      character(len=6) :: dum6
      character(len=500) :: txt,text
      real*8 :: time0,time1,qread

      open(unit=1,file='Last_abund.in',status='replace')
      do i=1,NEPS
        el = elnr(i)
        write(1,'(A2,2x,0pF25.20)') elnam(el),12.0+LOG10(eps0(el))
      enddo
      write(1,*) nHtot,T
      close(1)

      if (verbose>=0) then
        write(*,*)
        write(*,'("EQUIL_COND started")')
      endif  
      call CPU_TIME(time0)

      !------------------------
      ! ***  initial state  ***
      !------------------------
      ddust  = 0.Q0                 ! initial state dust-free
      eps    = eps0                 ! initial gas abundances 
      active = .false.              ! no solid condensing

      !--------------------------------------------
      ! ***  load initial state from database?  ***
      !--------------------------------------------
      Nact = 0
      Nact_read = 0
      act_read = .false.
      if (useDatabase) then
        call GET_DATA(nHtot,T,epsread,ddustread,qread,iread,
     >                act_read,method)
        if (qread.lt.1.0) then
          eps    = epsread
          ddust  = ddustread
          active = act_read
          text = "active solids:"
          Nact_read = 0
          do i=1,NDUST
            if (.not.act_read(i)) cycle
            Nact_read = Nact_read + 1
            text = trim(text)//" "//trim(dust_nam(i))
          enddo
          Nact = Nact_read
          !verbose = 0
          !if (qread>1.Q-3.and.Nact>0) verbose=2
          !if (qread>1.Q-3.and.iread==207) verbose=2
          !if (method==2) verbose=2
          if (verbose>0) then
            write(*,'(" ... using database entry (",I6,
     >            ") qual=",1pE15.7)') iread,qread
            write(*,*) trim(text)
          endif  
        endif  
      endif
  
      !----------------------------------------------------
      ! ***  recompute eps00 from initial state,        ***
      ! ***  because these can very slightly drift      ***
      ! ***  eps00 = total element abundances gas+dust  ***
      !----------------------------------------------------
      check = eps
      do i=1,NDUST
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          check(el) = check(el) + ddust(i)*dust_nu(i,j)    
        enddo
      enddo
      worst = 0.Q0
      do i=1,NEPS
        el = elnr(i)
        worst = MAX(worst,ABS(1.Q0-check(el)/eps0(el)))
      enddo
      !eps00 = check
      eps00 = eps0
      if (verbose>0) then
        write(*,*) "element conservation error 1:",worst
        write(*,*) "initial gas fractions ..."
        do i=1,NELEM
          if (elcode(i)==0) cycle
          print'(3x,A2,2(1pE15.6))',elnam(i),eps(i),eps(i)/eps0(i)
        enddo
      endif  
      if (worst>1.Q-8) stop "*** worst>1.Q-8 in equil_cond"

      !----------------------------------------------------------
      ! ***  compute maximum possible dust abundances dscale  ***
      ! ***  if all elements turn into selected dust species  ***
      !----------------------------------------------------------
      do i=1,NDUST
        xmin = 9.Q+99 
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          xmin = min(xmin,eps00(el)/REAL(dust_nu(i,j),kind=qp))    
        enddo
        dscale(i) = xmin                         ! max dust abundances
      enddo   
      dscalemax = MAXVAL(dscale(1:NDUST))

      null(:) = 0.Q0             
      call SUPER(nHtot,T,null,eps,Sat0,.false.) ! from scratch
      qual = SQUAL(Sat0,active)
      if (verbose>=0) then
        print'("it =",I4," qual =",1pE13.4E4)',0,qual
      endif  
      if (T>Tfast) then
        Qfinish = 1.Q-12
        Sfinish = 1.Q-09
      else
        Qfinish = 1.Q-20
        Sfinish = 1.Q-15
      endif
      switchedON(:) = 0
      switchedOFF(:) = 0
      if (verbose>-1) write(97,*)
      act_old = active
      lastit = -99
      iminoff = 0
      limited = .false.
      ifail = 0
      no_action = 0
      
      do it=1,itmax

        !---------------------------------------
        ! ***  selection of solids to solve  ***
        !---------------------------------------
        changed = .false.
        Smax = maxval(Sat0)
        ioff = 0
        if ((qread<0.5).and.(it<=3).and.(Nact_read>0)
     >                             .and.(qual>0.Q0)) then
          active = act_read
          Nact = Nact_read
          pot = 0.0
        else if (it>lastit+3) then
          maxon   = 0.Q0 
          minoff  = 0.Q0 
          imaxon  = 0
          if (qual==0.Q0) limited=.false.
          do i=1,NDUST
            xmin = 9.Q+99 
            esum = 0.Q0
            emax = 0.Q0
            do j=1,dust_nel(i)
              el = dust_el(i,j)
              esum = esum + dust_nu(i,j)
              if (eps(el)/REAL(dust_nu(i,j),kind=qp).lt.xmin) then    
                xmin = eps(el)/REAL(dust_nu(i,j),kind=qp)
                emax = REAL(dust_nu(i,j),kind=qp)
              endif
            enddo
            esum = esum**1.65            ! simple condensates first
            abun = dscalemax/dscale(i)   ! potentially abundant condensates first
            abun = MIN(abun,1.Q+10)
            pot(i)  = 1.Q0/(esum+0.05*abun+1.0*switchedOFF(i))
            Sat1(i) = Sat0(i)**pot(i)
            if (verbose>1.and.(.not.active(i)).and.Sat0(i)>1.Q0) then
              print'(A20," simplicity=",1pE10.3," abun=",1pE10.3,
     >                   " pot=",3(1pE10.3))',
     >           dust_nam(i),1.Q0/esum,1.Q0/abun,pot(i),Sat0(i),Sat1(i)
            endif  
          enddo 
          Smax = 0.Q0
          imax = 0
          do i=1,NDUST
            if (Sat1(i)>Smax) then
              Smax = Sat1(i)
              imax = i
            endif  
            if (Sat1(i)>1.Q0.and.(.not.active(i))) then
              turnon = Sat1(i)-1.Q0 
              if (turnon>maxon.and.(.not.limited)) then
                maxon  = turnon
                imaxon = i
              endif  
            endif  
          enddo
          if (qual>maxon.and.no_action<50) then  ! keep on iterating without switching on yet
            maxon = 0.0
            imaxon = 0
          endif  
          if (verbose>0) print'("limited=",L1,
     >                   "  Smax=",1pE10.3,2x,A18)',
     >                   limited,Smax,dust_nam(imax)
          if (verbose>0.and.maxon>0.Q0) print'("  maxon =",
     >                   1pE10.2,2x,A18)',maxon,dust_nam(imaxon)

          if (maxon>0.0*MAX(Smax-1.Q0,0.Q0)) then
            if (imaxon.ne.iminoff) then
              txt = dust_nam(imaxon)
              i = index(txt,'[l]')
              if (i>0) then
                txt = trim(txt(1:i-1))//'[s]'
                do i=1,NDUST
                  if (active(i).and.trim(dust_nam(i)).eq.trim(txt)) then
                    !print*,dust_nam(imaxon),Sat1(imaxon)
                    !print*,dust_nam(i),Sat1(i)
                    if (Sat1(imaxon)<Sat1(i)) imaxon=0
                  endif
                enddo  
              endif  
            endif
            if (imaxon>0) then
              active(imaxon) = .true.
              if (verbose>=0) then
                print*,"switch on ",trim(dust_nam(imaxon))
              endif  
            endif  
          endif  
          Nact = 0
          do i=1,NDUST
            active(i) = active(i).or.(ddust(i)>0.Q0)
            if (active(i).neqv.act_old(i)) changed=.true.
            if (active(i)) Nact=Nact+1
          enddo
          
          if (imaxon>0.and.Nact>1) then
            !----------------------------------------
            ! ***  eliminate linear combinations  ***
            !----------------------------------------
            eps_save = eps
            dust_save = ddust
            ioff = 0
            e_act(:) = .false.
            do i=1,NDUST
              if (.not.active(i)) cycle 
              do j=1,dust_nel(i)
                el = dust_el(i,j)
                e_act(el) = .true.
              enddo
            enddo  
            eind(:) = 0
            erow = 0
            do el=1,NELEM
              if (.not.e_act(el)) cycle
              erow = erow+1
              eind(el) = erow
            enddo  
            Eact = erow
            AA(:,:) = 0.Q0
            bvec(:) = 0.Q0
            irow = 0
            do i=1,NDUST
              if (.not.active(i)) cycle 
              if (i.ne.imaxon) then
                irow = irow+1
                dind(irow) = i
                do j=1,dust_nel(i)
                  el = dust_el(i,j)
                  erow = eind(el)
                  AA(erow,irow) = dust_nu(i,j)  
                enddo
              else
                dind(Nact) = i
                do j=1,dust_nel(i)
                  el = dust_el(i,j)
                  erow = eind(el)
                  bvec(erow) = dust_nu(i,j)  
                enddo
              endif  
            enddo  
            if (verbose>0) then
              print*,"searching for linear combination ..." 
              print'(2x,99(A10))',(trim(dust_nam(dind(i))),i=1,Nact) 
              do el=1,NELEM
                if (.not.e_act(el)) cycle
                erow = eind(el)
                print'(A2,99(F10.1))',trim(elnam(el)),AA(erow,1:Nact-1),
     >                                bvec(erow)
              enddo
            endif  
            call GAUSS_NM(NELEM,NDUST,Eact,Nact-1,AA,xvec,bvec,info)
            if (verbose>0) print'(" GAUSS_NM info =",I2)',info
            if (info.eq.0) then
              Nlin = 1
              dlin(1) = imaxon
              slin(imaxon) = -1.Q0
              do i=1,Nact-1
                if (ABS(xvec(i))<1.Q-25) cycle
                Nlin = Nlin+1
                dlin(Nlin) = dind(i)
                slin(dind(i)) = xvec(i)
              enddo  
              txt = trim(dust_nam(dlin(1)))//" <-> "
              do i=2,Nlin
                dk = dlin(i) 
                write(dum6,'(F6.3)') slin(dk)
                txt = trim(txt)//dum6//" "//trim(dust_nam(dk))
              enddo
              if (verbose>=0) then
                print*,"linear combination found: "//trim(txt)
              endif  
              itried(:) = .false.
              do
                Smin = 9.Q+99
                ioff = 0
                do i=1,Nlin
                  dk = dlin(i) 
                  if (dk==imaxon) cycle
                  if (Sat0(dk)<Smin.and.(.not.itried(dk))) then
                    Smin = Sat0(dk) 
                    ioff = dk
                  endif  
                enddo
                if (ioff==0) stop "*** ioff=0 should not occur"
                amount = ddust(ioff)
                ok = .true.
                do i=1,Nlin
                  dk = dlin(i)
                  if (dk==ioff) cycle
                  if (ddust(dk)-slin(dk)/slin(ioff)*amount<0.Q0) then
                    ok=.false.
                  endif
                enddo
                if (ok) exit
                itried(ioff) = .true.
              enddo
              changed = .true.
              active(ioff) = .false.  
              Nt = REAL(Nlin-1,kind=qp)
              amount = ddust(ioff)/Nt
              do i=1,Nlin
                dk = dlin(i)
                if (dk==ioff) cycle
                call TRANSFORM(ioff,dk,amount,-slin(dk)/slin(ioff)*Nt,
     >                         ddust,eps,dscale,ok)
              enddo  
              ddust(ioff) = 0.Q0
              eps = eps_save
              if (verbose>=0) then
                print*,"switch off ",dust_nam(ioff)
              endif  
            endif  
          endif
        endif  
        if (ioff>0) itransform=itransform+1

        action = .false.
        do i=1,NDUST
          if (active(i).and.(.not.act_old(i))) then
            action = .true.
            laston = i
            switchedON(i) = switchedON(i)+1
            if (verbose>-1) then
              write(97,*) it,"on  ",dust_nam(i),switchedON(i)
            endif  
          else if (act_old(i).and.(.not.active(i))) then
            action = .true.
            switchedOFF(i) = switchedOFF(i)+1
            if (verbose>-1) then
              write(97,*) it,"off ",dust_nam(i),switchedOFF(i)
            endif  
          endif
        enddo
        no_action = no_action+1
        if (action) no_action=0
        if (action) last_dx=0.Q0
  
        if (changed) then
          Nact = 0 
          do i=1,NDUST
            if (active(i)) Nact=Nact+1
          enddo
          call SUPER(nHtot,T,null,eps,Sat0,NewFastLevel<1)
          qual = SQUAL(Sat0,active)
          if (verbose>=0) then
            print'("it =",I4," qual =",1pE13.4E4)',it,qual
          endif  
          lastit = it
        endif
        if (verbose>0) then
          do i=1,NDUST
            rem = "  "
            if (active(i)) rem=" *"
            if (verbose>=0.and.(active(i).or.Sat0(i)>0.1)) then
              write(*,'(3x,A18,2(1pE11.3)1pE19.10,I3,1pE11.3,A2)') 
     >          dust_nam(i),ddust(i),ddust(i)/dscale(i),Sat0(i),
     >          switchedOFF(i),pot(i),rem
            endif  
          enddo
          do i=1,NDUST
            if (active(i).and.(.not.act_old(i))) then
              print*,"... switching on "//trim(dust_nam(i)) 
            else if (.not.active(i).and.act_old(i)) then
              print*,"... switching off "//trim(dust_nam(i)) 
            endif
          enddo   
        endif
        act_old = active
        if (Nact==0.and.qual<1.Q-30) exit   ! no solid supersaturated 

        !-----------------------------------------------
        ! ***  fill in r.h.s. vector FF and          ***
        ! ***  determine current quality of solution ***
        !-----------------------------------------------
        qual = SQUAL(Sat0,active)
        ii = 0
        do i=1,NDUST
          if (.not.active(i)) cycle
          ii = ii+1
          act_to_dust(ii) = i                ! index of active condensate 
          dust_to_act(i) = ii                ! to achieve Sat(i) = 1
          !FF(ii) = Sat0(i)-1.Q0
          FF(ii) = LOG(Sat0(i))              ! the function to be nullified
          elim = 0
          xmin = 9.Q+99
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            xx = eps(el)/REAL(dust_nu(i,j),kind=qp)
            if (xx<xmin) then
              xmin = xx
              elim = el
            endif  
          enddo
          !print*,dust_nam(i),"limited by",elnam(elim),xmin
          scale(ii) = xmin
        enddo 
        Nsolve = ii

        !------------------------------------------
        ! ***  compute numerical derivative DF  ***
        !------------------------------------------
        if (T>Tfast) then
          dlim = 1.Q-12
          target = 1.Q-5
        else
          dlim = 1.Q-30
          target = 1.Q-10
        endif  
        do jj=1,Nsolve
          j = act_to_dust(jj)                ! index of condensate as variable
          del = 1.Q-10*scale(jj) 
          !del = MAX(del,1.Q-12*last_dx(jj))
          !del = MAX(del,1.Q-25*scale(jj))
          xstep = 0.Q0                       ! evaporate a bit of that condensate
          do k=1,dust_nel(j)
            el = dust_el(j,k)
            xstep(el) = del*dust_nu(j,k)
          enddo
          call SUPER(nHtot,T,xstep,eps,Sat2,NewFastLevel<1)
          do ii=1,Nsolve
            i = act_to_dust(ii) 
            !dterm = Sat0(i)-Sat2(i)  !FF(ii)-Sat2(i)+1.Q0/Sat2(i)
            dterm = LOG(Sat0(i)/Sat2(i))
            DF(ii,jj) = dterm/del*scale(jj)
          enddo  
          print'("JAC:",A16,99(1pE10.2))',trim(dust_nam(j)),
     >                  scale(jj),del,last_dx(jj),
     >                  DF(1:Nsolve,jj)
        enddo
  
        !--------------------------------
        ! ***  Newton-Raphson step dx ***
        !--------------------------------
        Fsav  = FF
        DFsav = DF
        call GAUSS16( NELEM, Nsolve, DF, dx, FF )       ! the NR-step in condensates

        !call QGEFA( DF, NELEM, Nsolve, ipvt, info )
        !call QGESL( DF, NELEM, Nsolve, ipvt, FF, 0 )
        !dx = FF                                        
        !print*,"QGESL info=",info
        !if (info.ne.0) then
        !  FF = Fsav
        !  DF = DFsav
        !  call GAUSS16( NELEM, Nsolve, DF, dx, FF )
        !endif  

        xstep = 0.Q0                                    ! the corresponding NR-step in elements
        do ii=1,Nsolve
          dx(ii) = dx(ii)*scale(ii)
          i = act_to_dust(ii) 
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            xstep(el) = xstep(el) + dx(ii)*dust_nu(i,j)
          enddo
        enddo  

        !-----------------------------------
        ! ***  limit NR step physically  ***
        !-----------------------------------
        fac = 1.Q0
        iminoff = 0
        limdust = .false.
        do ii=1,Nsolve
          i = act_to_dust(ii)
          del = -dx(ii)
          !if (ddust(i)+del<0.Q0) 
     >    !  print'("dk,laston,it,lastit=",2(A20),2(i4))',
     >    !       trim(dust_nam(i)),trim(dust_nam(laston)),it,lastit

          if (del<0.Q0.and.ddust(i)>0.1*dscale(i)) then ! try small ddust before switching off
            fac2 = (-ddust(i)+0.05*dscale(i))/del       ! ddust+fac*del = 0.05*dscale
            if (fac2<1.Q0.and.verbose>0) print*,"*** limiting dust 1 "
     >                                   //dust_nam(i),REAL(fac2)
            if (fac2<fac) then
              fac = fac2
              iminoff = 0
              !limdust = .true.
            endif  
          else if (ddust(i)+del<0.Q0.and.
     >             i==laston.and.it<lastit+5) then      ! just switched on: keep on trying
            fac2 = -0.9*ddust(i)/del                    ! ddust+fac*del = ddust/10
            if (fac2<1.Q0.and.verbose>0) print*,"*** limiting dust 2 "
     >                                  //dust_nam(i),REAL(fac2)
            if (fac2<fac) then
              fac = fac2
              iminoff = 0
              limdust = .true.
            endif  
          else if (ddust(i)+del<0.Q0) then              ! preparation to switch off next
            fac2 = (-ddust(i)-small*dscale(i))/del      ! ddust+fac*del = -small*dscale
            if (fac2<1.Q0.and.verbose>0) print*,"*** limiting dust 3 "
     >                                   //dust_nam(i),REAL(fac2)
            if (fac2<fac) then
              fac = fac2 
              iminoff = i
              limdust = .true.
            endif
          endif  
        enddo
        do i=1,NELM
          if (i==iel) cycle
          el = elnum(i)
          del = xstep(el)
          if (eps(el)+del<0.02*eps(el)) then
            fac2 = (-0.98*eps(el))/del                  ! eps+fac*del = 0.02*eps
            if (verbose>0) print'(" *** limiting element1 ",A2,
     >        " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
            if (fac2<fac) then
              fac = fac2 
              iminoff = 0
              limdust = .false.
            endif
          endif
        enddo  
        dx    = dx*fac
        xstep = xstep*fac

        do ii=1,Nsolve
          i = act_to_dust(ii) 
          print'(A16,2(1pE18.10))',dust_nam(i),ddust(i),-dx(ii)
        enddo  
        do i=1,NELM
          if (i==iel) cycle
          el = elnum(i)
          print'(A3,2(1pE18.10))',elnam(el),eps(el),xstep(el)
        enddo

        limited = (fac<1.Q0)
        if (iminoff>0) then
          if (verbose>=0) print*,"switch off ",dust_nam(iminoff) 
          active(iminoff) = .false.
          lastit = -99
        endif

        !------------------------------------
        ! ***  apply dx to ddust and eps  ***
        !------------------------------------
        eps_save = eps
        dust_save = ddust
        NRstep = 1.Q0
        qold = qual
        do iback=1,7
          eps = eps_save
          ddust = dust_save
          do ii=1,Nsolve
            i = act_to_dust(ii)
            ddust(i) = ddust(i) - NRstep*dx(ii)
          enddo
          do i=1,NELM
            if (i==iel) cycle
            el = elnum(i)
            eps(el) = eps(el) + NRstep*xstep(el)
          enddo  
          call SUPER(nHtot,T,null,eps,Sat0,NewFastLevel<1)
          qual = SQUAL(Sat0,active)
          if (verbose>0) print'("--> pullback",i3,0pF6.3," Q=",1pE11.3,
     >                          " ->",1pE11.3)',iback,NRstep,qold,qual
          if (qual<qold*1.5) exit
          if (qual<1.0) exit
          if (limdust) exit
          if (iback==1) then
            NRstep = 0.9*NRstep
          else  
            NRstep = 0.5*NRstep
          endif  
        enddo
        last_dx = NRstep*ABS(dx)

        !-------------------------------------
        ! ***  check element conservation  ***
        !-------------------------------------
        !check = eps
        !do i=1,NDUST
        !  do j=1,dust_nel(i)
        !    el = dust_el(i,j)
        !    check(el) = check(el) + ddust(i)*dust_nu(i,j)    
        !  enddo
        !enddo
        !worst = 0.d0
        !do i=1,NEPS
        !  el = elnr(i)
        !  worst = MAX(worst,ABS(1.Q0-check(el)/eps00(el)))
        !enddo
        !if (verbose>1.or.worst>1.Q-8) write(*,*) 
     >  !   "element conservation error 2:",worst
        !if (worst>1.Q-8) stop

        Smax = maxval(Sat0)
        print'("it =",I4," qual =",1pE13.4E4)',it,qual
        if ((Smax<1.Q0+Sfinish).and.(qual<Qfinish)) exit
        if (verbose>0) read(*,'(a1)') char1

      enddo  
      Sat = Sat0

      call CPU_TIME(time1)
      if (it.lt.itmax) then
        if (verbose>=0) then
          write(*,'("EQUIL_COND converged after ",I3," iter, time =",
     >          0pF7.3," CPU sec.")') it,time1-time0
        endif  
      else
        write(*,'("*** EQUIL_COND failed after ",I3," iter,  time =",
     >            0pF9.4," CPU sec.")') it,time1-time0 
        stop
      endif

      !-------------------------
      ! ***  check solution  ***
      !-------------------------
      do i=1,NDUST
        if (ddust(i)>0.Q0.and.Sat(i)<0.9999) then
          print*,"*** error: ddust>0 but S<1"
          print*,dust_nam(i),REAL(ddust(i)),REAL(Sat(i))
          stop
        endif
        if (Sat(i)>1.00001) then
          print*,"*** error: S>1"
          print*,dust_nam(i),REAL(ddust(i)),REAL(Sat(i))
          stop
        endif
        if (ddust(i)<-10*small*dscale(i)) then
          print*,"*** error: ddust<0"
          print*,dust_nam(i),REAL(ddust(i)),REAL(Sat(i))
          stop
        endif  
      enddo  

      !-------------------------------------
      ! ***  check element conservation  ***
      !-------------------------------------
      check = eps
      do i=1,NDUST
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          check(el) = check(el) + ddust(i)*dust_nu(i,j)
        enddo
      enddo
      do i=1,NEPS
        el = elnr(i)
        if (ABS(1.Q0-check(el)/eps00(el))>1.Q-8) then
          print*,"*** element conservation error"
          print*,elnam(el),check(el),eps00(el)
          stop
        endif  
      enddo

      !----------------------------------
      ! ***  save result to database  ***
      !----------------------------------
      !if (qual<1.Q-10.and.useDatabase) then
      if (useDatabase) then
        call PUT_DATA(nHtot,T,eps,ddust,qread,iread,active)
      endif  
      ieqcond = ieqcond + 1
      ieqconditer = ieqconditer + it

      end
            

!-------------------------------------------------------------------------
      subroutine SUPER(nHtot,T,xx,eps,Sat,merk)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,
     >                    dust_nam,elnam
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in)  :: nHtot,T
      real(kind=qp),intent(in) :: xx(NELEM),eps(NELEM)
      real(kind=qp),intent(out) :: Sat(NDUST)
      logical,intent(in) :: merk
      real(kind=qp) :: eps1(NELEM)

      eps1 = eps + xx
      !----------------------------------------------
      ! ***  compute chemistry & supersaturation  ***
      !----------------------------------------------
      call GGCHEM(nHtot,T,eps1,merk,0)
      call SUPERSAT(T,nat,nmol,Sat)
      
      end
      

!-------------------------------------------------------------------------
      function SQUAL(Sat,active)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NDUST
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),intent(in) :: Sat(NDUST)
      real(kind=qp) :: SQUAL,qual
      logical,intent(in) :: active(0:NDUST)
      integer :: i

      qual = 0.d0
      do i=1,NDUST
        if (active(i)) then
         !qual = qual + (1.Q0-Sat(i))**2
          qual = qual + (Sat(i)-1.Q0/Sat(i))**2
         !qual = qual + LOG(Sat(i))**2
        else if (Sat(i).gt.1.Q0) then
         !qual = qual + MIN(Sat(i)-1.Q0,1.Q0)
        endif  
      enddo
      SQUAL = qual
      end

!-------------------------------------------------------------------------
      subroutine VAPORIZE(i,ddust,eps)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp) :: del
      integer :: j,el
      
      del = ddust(i)
      print*," ==>  vaporize "//trim(dust_nam(i)),REAL(del)
      ddust(i) = 0.Q0
      do j=1,dust_nel(i)
        el = dust_el(i,j)
        eps(el) = eps(el) + del*dust_nu(i,j)    
      enddo
      end

!-------------------------------------------------------------------------
      subroutine TRANSFORM(i1,i2,del,fac,ddust,eps,dscale,ok)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam,
     >                    eps0
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i1,i2
      real(kind=qp),parameter :: dsmall=1.Q-30
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp),intent(in) :: del,fac,dscale(NDUST)
      logical,intent(inout) :: ok
      integer :: j,el
      
      print*," ==>  transform "//trim(dust_nam(i1))//" -> "
     &       //trim(dust_nam(i2)),REAL(fac*del/dscale(i1))
      ddust(i1) = ddust(i1)-del
      ddust(i2) = ddust(i2)+fac*del
      do j=1,dust_nel(i1)
        el = dust_el(i1,j)
        eps(el) = eps(el) + del*dust_nu(i1,j)    
      enddo
      do j=1,dust_nel(i2)
        el = dust_el(i2,j)
        eps(el) = eps(el) - fac*del*dust_nu(i2,j)    
      enddo

      if (ddust(i1)<-dsmall.or.ddust(i2)<-dsmall) ok=.false.
      !if (ddust(i1)<-dsmall) then
      !  call VAPORIZE(i1,ddust,eps)
      !  active(i1) = .false.
      !  vap = .true.
      !endif  
      !if (ddust(i2)<-dsmall) then
      !  call VAPORIZE(i2,ddust,eps)
      !  active(i2) = .false.
      !  vap = .true.
      !endif
  
      !-------------------------------------
      ! ***  check element conservation  ***
      !-------------------------------------
      !check = eps
      !do i=1,NDUST
      !  do j=1,dust_nel(i)
      !    el = dust_el(i,j)
      !    check(el) = check(el) + ddust(i)*dust_nu(i,j)    
      !  enddo
      !enddo
      !worst = 0.d0
      !do i=1,NEPS
      !  el = elnr(i)
      !  worst = MAX(worst,ABS(1.Q0-check(el)/eps00(el)))
      !enddo
      !write(*,*) "element conservation error 3:",worst
      !if (worst>1.Q-8) stop

      end

!-------------------------------------------------------------------------
      subroutine GET_KNOWNS(eq,mat,emat,e_act,d_resolved,e_resolved,
     >                      elem,Nslot,dustkind,stoich,vec,verbose)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nam,dust_nel,dust_nu,dust_el,
     >                    elnam 
      use CONVERSION,ONLY: Nind,Ndep,Iindex,Dindex,is_dust
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: eq
      integer,intent(in),dimension(NELEM) :: elem,Nslot
      integer,intent(in),dimension(NELEM,NDUST) :: dustkind,stoich
      logical,intent(in),dimension(NDUST) :: d_resolved
      logical,intent(in),dimension(NELEM) :: e_resolved,e_act
      integer,intent(in) :: verbose
      real(kind=qp),intent(in) :: mat(NDUST,NELEM)
      real(kind=qp),intent(in) :: emat(NELEM,NELEM)
      real(kind=qp),intent(out):: vec(NELEM)
      integer :: i,j,el,el2,slots,sl

      vec(:) = 0.Q0
      slots = Nslot(eq)
      do sl=1,slots
        i = dustkind(eq,sl)
        if (.not.d_resolved(i)) cycle
        do j=1,Nind
          el2 = Iindex(j)
          if (mat(i,el2).eq.0.Q0) cycle
          if (verbose>1) print*,"out1 ",trim(dust_nam(i))//" "
     >           //elnam(el2),REAL(mat(i,el2)),stoich(eq,sl)
          vec(el2) = vec(el2)+ mat(i,el2)*stoich(eq,sl)
        enddo
      enddo  
      el = elem(eq)
      if (e_act(el)) then
        vec(el) = vec(el) + 1.Q0 
      else if (e_resolved(el)) then
        do j=1,Nind
          el2 = Iindex(j)
          if (emat(el,el2).eq.0.Q0) cycle
          vec(el) = vec(el) + emat(el,el2)
          if (verbose>1) print*,"out2 ",trim(elnam(el))//" "
     >           //elnam(el2),REAL(emat(el,el2))
        enddo  
      endif
      end
