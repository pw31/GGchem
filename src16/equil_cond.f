      module CONVERSION
      use DUST_DATA,ONLY: NELEM,NDUST
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: Nind,Ndep,Iindex(NELEM),Dindex(NDUST+NELEM)
      logical :: is_dust(NDUST+NELEM)
      real(kind=qp) :: conv(NDUST+NELEM,NELEM)
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
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nam,dust_nel,dust_nu,dust_el,
     >                    eps0,elnam 
      use CONVERSION,ONLY: Nind,Ndep,Iindex,Dindex,is_dust,conv
      use EXCHANGE,ONLY: Fe,Mg,Si,Al,Ca,Ti,O,S,Na,Kalium=>K,Cl,H
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nHtot                ! H nuclei density [cm-3]
      real*8,intent(in) :: T                    ! temperature [K]
      real(kind=qp),intent(out) :: eps(NELEM)   ! gas element abundances
      real(kind=qp),intent(out) :: Sat(NDUST)   ! saturation ratio
      real(kind=qp),intent(out) :: ddust(NDUST) ! density of solid units [cm-3]
      integer,intent(inout) :: verbose
      real(kind=qp),dimension(NELEM) :: eps00,epsread,check,FF,Fsav,dx
      real(kind=qp),dimension(NELEM) :: eps_save,vec,xstep,Iabund,work
      real(kind=qp),dimension(NELEM) :: scale
      real(kind=qp),dimension(NDUST) :: ddustread,dscale,pot
      real(kind=qp),dimension(NDUST) :: Sat0,Sat1,Sat2
      real(kind=qp),dimension(NELEM,NELEM) :: DF,DFsav,emat,vecs
      real(kind=qp),dimension(NDUST,NELEM) :: mat
      real(kind=qp) :: worst,xmin,Smax,qual,SQUAL,del
      real(kind=qp) :: turnon,turnoff,maxon,minoff,fac,fac2,amount
      real(kind=qp) :: deps1,deps2,deps
      real(kind=qp) :: det(2),converge(500,NELEM),crit,cbest
      real(kind=qp) :: small=1.Q-30
      integer,parameter :: itmax=200
      integer,dimension(NELEM) :: elem,Nslot
      integer,dimension(NELEM,NDUST) :: dustkind,stoich
      integer :: it,i,j,k,el,el2,Nact,Nact_read,Neq,slots,sl,dk,eq
      integer :: itry,knowns,unknowns,unknown,ii,jj,lastit
      integer :: imaxon,iminoff,info,ipvt(NELEM),e_num(NELEM)
      integer :: Nunsolved,unsolved(NELEM),Nvar1,Nvar2,var(NELEM)
      integer :: Nsolve,ebest,dbest,nonzero,itrivial,iread
      integer :: act_to_elem(NELEM),act_to_dust(NELEM)
      integer :: Nzero,Ntrivial,etrivial(NELEM),dtrivial(NELEM)
      logical,dimension(NELEM) :: e_resolved,e_act,e_taken,is_esolved
      logical,dimension(NDUST) :: active,act_read,act_old,d_resolved
      logical,dimension(NDUST) :: is_dsolved
      logical :: action,changed,solved,limited
      character(len=1) :: char1,txt0
      character(len=2) :: rem
      character(len=500) :: txt,txt1,txt2,text,filename
      logical,save :: firstCall=.true.
      integer,save :: iAl2O3,iFe,iFeS,iNa2SiO3,iMgSiO3,iMg2SiO4,iTi4O7
      integer,save :: iCaSiO3,iCaMgSi2O6,iNaAlSi3O8,iMgAl2O4,iCaTiO3
      integer,save :: iSiO,iSiO2,iTiO2,iMgTi2O5,iSiC,iCaS,iFe2SiO4,iFeO
      integer,save :: iNaCl,iKCl,iKAlSi3O8,iH2O,iH2O_l,iFe_l,iMgAl2O4_l
      integer,save :: it_tot=0, sit_tot=0, fail_tot=0
      real*8 :: time0,time1,qread

      if (firstCall) then
        do i=1,NDUST
          if (dust_nam(i).eq.'Al2O3[s]')     iAl2O3=i 
          if (dust_nam(i).eq.'SiO[s]')       iSiO=i 
          if (dust_nam(i).eq.'SiO2[s]')      iSiO2=i 
          if (dust_nam(i).eq.'Fe[s]')        iFe=i 
          if (dust_nam(i).eq.'Fe[l]')        iFe_l=i 
          if (dust_nam(i).eq.'FeS[s]')       iFeS=i 
          if (dust_nam(i).eq.'Na2SiO3[s]')   iNa2SiO3=i 
          if (dust_nam(i).eq.'MgSiO3[s]')    iMgSiO3=i 
          if (dust_nam(i).eq.'Mg2SiO4[s]')   iMg2SiO4=i 
          if (dust_nam(i).eq.'CaSiO3[s]')    iCaSiO3=i
          if (dust_nam(i).eq.'CaTiO3[s]')    iCaTiO3=i
          if (dust_nam(i).eq.'CaMgSi2O6[s]') iCaMgSi2O6=i
          if (dust_nam(i).eq.'NaAlSi3O8[s]') iNaAlSi3O8=i
          if (dust_nam(i).eq.'MgAl2O4[s]')   iMgAl2O4=i
          if (dust_nam(i).eq.'MgAl2O4[l]')   iMgAl2O4_l=i
          if (dust_nam(i).eq.'Ti4O7[s]')     iTi4O7=i
          if (dust_nam(i).eq.'TiO2[s]')      iTiO2=i
          if (dust_nam(i).eq.'MgTi2O5[s]')   iMgTi2O5=i
          if (dust_nam(i).eq.'SiC[s]')       iSiC=i
          if (dust_nam(i).eq.'CaS[s]')       iCaS=i
          if (dust_nam(i).eq.'Fe2SiO4[s]')   iFe2SiO4=i
          if (dust_nam(i).eq.'FeO[s]')       iFeO=i
          if (dust_nam(i).eq.'NaCl[s]')      iNaCl=i
          if (dust_nam(i).eq.'KCl[s]')       iKCl=i
          if (dust_nam(i).eq.'KAlSi3O8[s]')  iKAlSi3O8=i
          if (dust_nam(i).eq.'H2O[s]')       iH2O=i
          if (dust_nam(i).eq.'H2O[l]')       iH2O_l=i
        enddo
        firstCall = .false. 
      endif

      write(*,*)
      write(*,'("EQUIL_COND started")') 
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
      call GET_DATA(nHtot,T,epsread,ddustread,qread,iread,act_read)
      if (qread.lt.0.5) then
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
        verbose = 0
        !if (qread>1.Q-3.and.Nact>0) verbose=2
        !if (iread==2775) verbose=2
        if (verbose>0) then
          write(*,'(" ... using database entry (",I6,
     >          ") qual=",1pE15.7)') iread,qread
          write(*,*) trim(text)
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
      do i=1,NELEM
        worst = MAX(worst,ABS(1.Q0-check(i)/eps0(i)))
      enddo
      eps00 = check
      if (verbose>0) then
        write(*,*) "element conservation error:",worst
        write(*,*) "initial gas fractions ..."
        do i=1,NELEM
          if (eps00(i).gt.1.Q-20) then
            print'(3x,A2,2(1pE15.6))',elnam(i),eps(i),eps(i)/eps00(i)
          endif  
        enddo
      endif  
      if (worst>1.Q-8) stop

      !----------------------------------------------------------
      ! ***  compute maximum possible dust abundances dscale  ***
      ! ***  if all elements turn into selected dust species  ***
      !----------------------------------------------------------
      do i=1,NDUST
        xmin = 9.Q+99 
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          xmin = min(xmin,eps00(el)/dust_nu(i,j))    
        enddo
        dscale(i) = xmin                        ! max dust abundances
      enddo   

      call GGCHEM(nHtot,T,eps,.false.,0)        ! one call from scratch
      xstep(:) = 0.Q0             
      call SUPER(nHtot,T,xstep,eps,Sat0)
      qual = SQUAL(Sat0,active)
      print'("it=",I3," qual=",1pE12.4)',0,qual
      act_old = active
      lastit = -99
      iminoff = 0
      limited = .false.

      do it=1,itmax
        
        !---------------------------------------
        ! ***  selection of solids to solve  ***
        !---------------------------------------
        changed = .false.
        Smax = maxval(Sat0)
        if ((qread<0.5).and.(it<=3).and.(Nact_read>0)) then
          active = act_read
          Nact = Nact_read
        else if (it>lastit+3) then
          maxon   = 0.Q0 
          minoff  = 0.Q0 
          imaxon  = 0
          do i=1,NDUST
            xmin = 9.Q+99 
            pot(i) = 0.Q0
            do j=1,dust_nel(i)
              el = dust_el(i,j)
              !if (el.ne.O) pot(i)=pot(i)+dust_nu(i,j)
              if (eps(el)/dust_nu(i,j).lt.xmin) then    
                xmin = eps(el)/dust_nu(i,j)
                pot(i) = 1.Q0/DBLE(dust_nu(i,j))
              endif
            enddo
            Sat1(i)=Sat0(i)**(1.Q0/pot(i))
            !print'(A15,0pF8.3)',dust_nam(i),pot(i)
          enddo 
          do i=1,NDUST
            if (Sat1(i)>1.Q0.and.(.not.active(i))) then
              turnon = Sat1(i)-1.Q0 
              if (i==iMg2SiO4) turnon=turnon*10.0
              if (turnon>maxon.and..not.limited) then
                maxon  = turnon
                imaxon = i
              endif  
            endif  
          enddo  
          if (verbose>0.and.maxon>0.Q0) then
            print'("  maxon =",1pE10.2,2x,A12)',maxon,dust_nam(imaxon)
          endif  
          Smax = maxval(Sat1)
          if (maxon>0.1*MAX(Smax-1.Q0,0.Q0)) then
            if (imaxon.ne.iminoff) then 
              active(imaxon) = .true.
            endif  
          endif  
          Nact = 0
          do i=1,NDUST
            active(i) = active(i).or.(ddust(i)>0.Q0)
            if (active(i).neqv.act_old(i)) changed=.true.
            if (active(i)) Nact=Nact+1
          enddo

          !----------------------------------------------------
          ! ***  some explicit help with sudden conversions ***
          !----------------------------------------------------
          if (active(iKCl).and.active(iKalSi3O8).and.
     >        active(iNaCl).and.active(iMgAl2O4).and.
     >        active(iMgSiO3).and.active(iMg2SiO4).and.
     >        active(iNa2SiO3)) then
            changed = .true.
            eps_save = eps
            !--- decide ---
            if (Sat0(iKCl)>Sat0(iKAlSi3O8)) then
              active(iKAlSi3O8) = .false.
              amount = ddust(iKalSi3O8)/6.Q0
              call TRANSFORM(iKAlSi3O8,iKCl,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKAlSi3O8,iNaCl,amount,-1.Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKAlSi3O8,iNa2SiO3,amount,0.5Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKAlSi3O8,iMgAl2O4,amount,0.5Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKAlSi3O8,iMgSiO3,amount,5.5Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKAlSi3O8,iMg2SiO4,amount,-3.Q0*6.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iKCl) = .false.
              amount = ddust(iKCl)/6.Q0
              call TRANSFORM(iKCl,iKAlSi3O8,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKCl,iNaCl,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKCl,iNa2SiO3,amount,-0.5Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKCl,iMgAl2O4,amount,-0.5Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKCl,iMgSiO3,amount,-5.5Q0*6.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iKCl,iMg2SiO4,amount,3.Q0*6.Q0,
     >                       ddust,eps,dscale)
            endif
            eps(K)  = eps_save(K)
            eps(Cl) = eps_save(Cl)
            eps(Al) = eps_save(Al)
            eps(Na) = eps_save(Na)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif               
          if (active(iAl2O3).and.active(iMgAl2O4).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            eps_save = eps
            !--- decide ---
            if (Sat0(iAl2O3)>Sat0(iMgAl2O4)) then
              active(iMgAl2O4) = .false.
              amount = ddust(iMgAl2O4)/3.Q0
              call TRANSFORM(iMgAl2O4,iAl2O3,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgAl2O4,iMgSiO3,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgAl2O4,iMg2SiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iAl2O3) = .false.
              amount = ddust(iAl2O3)/3.Q0
              call TRANSFORM(iAl2O3,iMgAl2O4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iAl2O3,iMgSiO3,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iAl2O3,iMg2SiO4,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale)
            endif  
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif  
          !if (active(iCaMgSi2O6).and.active(iCaSiO3).and.
     >    !    .not.e_act(Mg).and..not.e_act(Si)) then
          !  changed = .true.
          !  !--- decide ---
          !  if (Sat0(iCaSiO3).gt.Sat0(iCaMgSi2O6)) then
          !    active(iCaMgSi2O6) = .false.  
          !    amount = ddust(iCaMgSi2O6)
          !    call TRANSFORM(iCaMgSi2O6,iCaSiO3,amount,1.Q0,
     >    !                   ddust,eps,dscale)
          !  else  
          !    active(iCaSiO3) = .false.  
          !    amount = ddust(iCaSiO3)
          !    call TRANSFORM(iCaSiO3,iCaMgSi2O6,amount,1.Q0,
     >    !                   ddust,eps,dscale)
          !  endif  
          !endif
          if (active(iFeO).and.active(iFe2SiO4).and.active(iSiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFe2SiO4).gt.Sat0(iFeO)) then
              active(iFeO) = .false.  
              amount = ddust(iFeO)/2.Q0
              call TRANSFORM(iFeO,iFe2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iFeO,iSiO2,amount,-0.5Q0*2.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iFe2SiO4) = .false.  
              amount = ddust(iFe2SiO4)/2.Q0
              call TRANSFORM(iFe2SiO4,iFeO,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iFe2SiO4,iSiO2,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif   
          if (active(iSiO).and.active(iMgSiO3).and.active(iMg2SiO4)
     >        .and.(.not.e_act(Mg).or..not.e_act(Si))) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgSiO3).gt.Sat0(iSiO)) then
              active(iSiO) = .false.  
              amount = ddust(iSiO)/2.Q0
              call TRANSFORM(iSiO,iMgSiO3,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iSiO,iMg2SiO4,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iMgSiO3) = .false.  
              amount = ddust(iMgSiO3)/2.Q0
              call TRANSFORM(iMgSiO3,iSiO,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgSiO3,iMg2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif
          if (active(iSiO2).and.active(iMgSiO3).and.
     >        active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgSiO3).gt.Sat0(iSiO2)) then
              active(iSiO2) = .false.  
              amount = ddust(iSiO2)/2.Q0
              call TRANSFORM(iSiO2,iMgSiO3,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iSiO2,iMg2SiO4,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iMgSiO3) = .false.  
              amount = ddust(iMgSiO3)/2.Q0
              call TRANSFORM(iMgSiO3,iSiO2,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgSiO3,iMg2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif
          if (active(iTi4O7).and.active(iCaTiO3).and.
     >        active(iCaMgSi2O6).and..not.e_act(Mg)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaTiO3).gt.Sat0(iTi4O7)) then
              active(iTi4O7) = .false.  
              amount = ddust(iTi4O7)/2.Q0
              call TRANSFORM(iTi4O7,iCaTiO3,amount,4.Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iTi4O7,iCaMgSi2O6,amount,-4.Q0*4.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iCaTiO3) = .false.  
              amount = ddust(iCaTiO3)/2.Q0
              call TRANSFORM(iCaTiO3,iTi4O7,amount,0.25Q0*2.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iCaTiO3,iCaMgSi2O6,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif  
          if (active(iTi4O7).and.active(iTiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTi4O7).gt.Sat0(iTiO2)) then
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)
              call TRANSFORM(iTiO2,iTi4O7,amount,0.25Q0,
     >                       ddust,eps,dscale)
            else  
              active(iTi4O7) = .false.  
              amount = ddust(iTi4O7)
              call TRANSFORM(iTi4O7,iTiO2,amount,4.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif   
          if (active(iMgTi2O5).and.active(iTiO2).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgTi2O5).gt.Sat0(iTiO2)) then
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)/3.Q0
              call TRANSFORM(iTiO2,iMgTi2O5,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iTiO2,iMgSiO3,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iTiO2,iMg2SiO4,amount,-0.5Q0*3.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iMgTi2O5) = .false.  
              amount = ddust(iMgTi2O5)/3.Q0
              call TRANSFORM(iMgTi2O5,iTiO2,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgTi2O5,iMgSiO3,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgTi2O5,iMg2SiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif   
          if (active(iMgTi2O5).and.active(iTiO2).and.
     >        active(iSiO2).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgTi2O5).gt.Sat0(iTiO2)) then
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)/3.Q0
              call TRANSFORM(iTiO2,iMgTi2O5,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iTiO2,iSiO2,amount,0.25Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iTiO2,iMg2SiO4,amount,-0.25Q0*3.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iMgTi2O5) = .false.  
              amount = ddust(iMgTi2O5)/3.Q0
              call TRANSFORM(iMgTi2O5,iTiO2,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgTi2O5,iSiO2,amount,-0.5Q0*3.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iMgTi2O5,iMg2SiO4,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale)
            endif  
          endif   
          if (active(iNa2SiO3).and.active(iNaAlSi3O8).and.
     >        active(iMgSiO3).and.active(iMg2SiO4).and.
     >        active(iMgAl2O4)) then
            changed = .true.
            eps_save = eps
            !--- decide ---
            if (Sat0(iNa2SiO3).gt.Sat0(iNaAlSi3O8)) then
              active(iNaAlSi3O8) = .false.  
              amount = ddust(iNaAlSi3O8)/4.Q0
              ddust(iNa2SiO3) = MAX(ddust(iNa2SiO3),0.Q0)
              !print*,eps(Na),eps(Al),eps(Mg),eps(Si)
              call TRANSFORM(iNaAlSi3O8,iNa2SiO3,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iNaAlSi3O8,iMgAl2O4,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iNaAlSi3O8,iMgSiO3,amount,5.5Q0*4.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iNaAlSi3O8,iMg2SiO4,amount,-3.Q0*4.Q0,
     >                       ddust,eps,dscale)
              !print*,eps(Na),eps(Al),eps(Mg),eps(Si)
            else  
              active(iNa2SiO3) = .false.  
              amount = ddust(iNa2SiO3)/4.Q0
              !print*,eps(Na),eps(Al),eps(Mg),eps(Si)
              call TRANSFORM(iNa2SiO3,iNaAlSi3O8,amount,2.0Q0*4.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iNa2SiO3,iMgAl2O4,amount,-1.0Q0*4.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iNa2SiO3,iMgSiO3,amount,-11.0Q0*4.Q0,
     >                       ddust,eps,dscale)
              call TRANSFORM(iNa2SiO3,iMg2SiO4,amount,6.Q0*4.Q0,
     >                       ddust,eps,dscale)
              !print*,eps(Na),eps(Al),eps(Mg),eps(Si)
            endif  
            eps(Na) = eps_save(Na)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iH2O).and.active(iH2O_l)) then
            changed = .true.
            eps_save = eps
            !--- decide ---
            if (Sat0(iH2O).gt.Sat0(iH2O_l)) then
              active(iH2O_l) = .false.  
              amount = ddust(iH2O_l)
              call TRANSFORM(iH2O_l,iH2O,amount,1.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iH2O) = .false.  
              amount = ddust(iH2O)
              call TRANSFORM(iH2O,iH2O_l,amount,1.Q0,
     >                       ddust,eps,dscale)
            endif  
            eps(O) = eps_save(O)
            eps(H) = eps_save(H)
          endif   
          if (active(iFe).and.active(iFe_l)) then
            changed = .true.
            eps_save = eps
            !--- decide ---
            if (Sat0(iFe).gt.Sat0(iFe_l)) then
              active(iFe_l) = .false.  
              amount = ddust(iFe_l)
              call TRANSFORM(iFe_l,iFe,amount,1.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iFe) = .false.  
              amount = ddust(iFe)
              call TRANSFORM(iFe,iFe_l,amount,1.Q0,
     >                       ddust,eps,dscale)
            endif  
            eps(Fe) = eps_save(Fe)
          endif   
          if (active(iMgAl2O4).and.active(iMgAl2O4_l)) then
            changed = .true.
            eps_save = eps
            !--- decide ---
            if (Sat0(iMgAl2O4).gt.Sat0(iMgAl2O4_l)) then
              active(iMgAl2O4_l) = .false.  
              amount = ddust(iMgAl2O4_l)
              call TRANSFORM(iMgAl2O4_l,iMgAl2O4,amount,1.Q0,
     >                       ddust,eps,dscale)
            else  
              active(iMgAl2O4) = .false.  
              amount = ddust(iMgAl2O4)
              call TRANSFORM(iMgAl2O4,iMgAl2O4_l,amount,1.Q0,
     >                       ddust,eps,dscale)
            endif  
            eps(Mg) = eps_save(Mg)
            eps(Al) = eps_save(Al)
            eps(O) = eps_save(O)
          endif   
        endif  

        if (changed) then
          Nact = 0 
          do i=1,NDUST
            if (active(i)) Nact=Nact+1
          enddo
          xstep(:)= 0.Q0             
          call SUPER(nHtot,T,xstep,eps,Sat0)
          qual = SQUAL(Sat0,active)
          print'("it=",I3," qual=",1pE12.4)',it,qual
          lastit = it
        endif
        if (verbose>0) then
          do i=1,NDUST
            rem = "  "
            if (active(i)) rem=" *"
            if (active(i).or.Sat0(i)>0.01) then
              write(*,'(3x,A12,2(1pE19.10),A2)') 
     >          dust_nam(i),ddust(i)/dscale(i),Sat0(i),rem
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
        if (it>1.and.(.not.changed)) goto 100

        !--------------------------------------
        ! ***  select independent elements  ***
        !--------------------------------------
        Nind = 1
        Iabund(1) = 9.E+99
        e_act(:) = .false.
        do i=1,NDUST
          if (.not.active(i)) cycle
          if (dust_nel(i)>1) cycle          ! include pure metals first
          el = dust_el(i,1)
          if (e_act(el)) cycle
          e_act(el) = .true.
          Iindex(2:Nind+1) = Iindex(1:Nind)
          Iabund(2:Nind+1) = Iabund(1:Nind)
          Iindex(1) = el
          Iabund(1) = 0.Q0
          !print*,1,elnam(Iindex(1:Nind))
          Nind = Nind+1
        enddo  
        do i=1,NDUST
          if (.not.active(i)) cycle         ! include elements in compounds
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            if (e_act(el)) cycle
            e_act(el) = .true.
            do k=1,Nind                  
              if (eps(el)<Iabund(k)) exit   ! sort by element abundance 
            enddo
            Iindex(k+1:Nind+1) = Iindex(k:Nind)
            Iabund(k+1:Nind+1) = Iabund(k:Nind)
            Iindex(k) = el
            Iabund(k) = eps(el) 
            !print*,2,elnam(Iindex(1:Nind))
            Nind = Nind+1
          enddo  
        enddo  
        e_act(:) = .false.
        e_num(:) = 0
        do i=1,Nact
          el = Iindex(i)
          e_act(el) = .true.
        enddo 
        do i=1,NDUST 
          if (.not.active(i)) cycle
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            e_num(el) = e_num(el)+1
          enddo  
        enddo  
        !do i=1,NDUST 
        !  if (active(i)) print*,dust_nam(i)
        !enddo  
        !print*,elnam(Iindex(1:Nind-1)) 
        !print*,e_num(Iindex(1:Nind-1)) 
        !print*,e_act(Iindex(1:Nind-1))
        Nind = Nact                         ! truncate at number of condensates
        !-------------------------------------
        ! ***  some explict special cases  ***
        !-------------------------------------
        if (active(iAl2O3).and.active(iMgAl2O4).and..not.e_act(Mg)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Mg"
          e_act(Iindex(Nact)) = .false.
          Iindex(Nact) = Mg
          e_act(Mg) = .true.
        endif
        if (active(iCaS).and.e_act(Ca).and.e_act(S).and.
     >      (e_num(Ca)==1).and.(e_num(S)==1)) then
          do i=1,Nind
            if (Iindex(i)==S) exit
          enddo  
          if ((e_num(Si)>0).and.(.not.e_act(Si))) then
            Iindex(i) = Si 
          else if ((e_num(Mg)>0).and.(.not.e_act(Mg))) then
            Iindex(i) = Mg
          else if ((e_num(O)>0).and.(.not.e_act(O))) then
            Iindex(i) = O
          else
            stop "CaS" 
          endif
        endif 
        if (active(iAl2O3).and..not.e_act(Al)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Al"
          e_act(Iindex(Nact)) = .false.
          Iindex(Nact) = Al
          e_act(Al) = .true.
        endif  
        if (active(iMgAl2O4).and..not.e_act(Mg).and..not.e_act(Al)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Al"
          e_act(Iindex(Nact)) = .false.
          Iindex(Nact) = Al
          e_act(Al) = .true.
        endif  
        if (active(iCaMgSi2O6).and.e_act(Ca).and.e_act(Si).and.
     >      e_num(Ca)==1.and.e_num(Si)==1) then
          e_act(Iindex(Nact)) = .false.
          e_act(Iindex(Nact+1)) = .true.
          Iindex(Nact) = Iindex(Nact+1)
        endif   
        if (verbose>1) print*,"solving for ... ",
     >                      (elnam(Iindex(i))//" ",i=1,Nind)
        if (verbose>1) print'(99(1pE11.3))',(Iabund(i),i=1,Nind)

        !------------------------------------------------
        ! ***  determine dependent dust and elements  ***
        !------------------------------------------------
        d_resolved(1:NDUST) = .false.
        e_resolved(1:NELEM) = .false.
        e_act(:) = .false.
        e_taken(:) = .false.
        do i=1,Nind
          el = Iindex(i)
          e_act(el) = .true.
          e_resolved(el) = .true.
        enddo  
        Ndep = 0
        do i=1,NDUST
          if (.not.active(i)) cycle
          Ndep = Ndep+1
          Dindex(Ndep) = i
          is_dust(Ndep) =.true.
        enddo
        do i=1,NDUST
          if (.not.active(i)) cycle
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            if (e_resolved(el)) cycle
            if (e_taken(el)) cycle
            Ndep = Ndep+1
            Dindex(Ndep) = el
            is_dust(Ndep) = .false.
            e_taken(el) = .true.
          enddo  
        enddo

        !----------------------------
        ! ***  conversion matrix  ***
        !----------------------------
        Neq = 1
        do el=1,NELEM
          slots = 0
          do i=1,NDUST
            if (.not.active(i)) cycle
            do j=1,dust_nel(i)
              if (el.ne.dust_el(i,j)) cycle
              elem(Neq) = el
              slots = slots+1
              dustkind(Neq,slots) = i
              stoich(Neq,slots) = dust_nu(i,j)
            enddo
          enddo
          if (slots>0) then
            Nslot(Neq) = slots
            Neq=Neq+1
          endif  
        enddo  
        Neq = Neq-1
        do itry=1,99
          action = .false. 
          solved = .true.
          Nunsolved = 0
          do eq=1,Neq
            el = elem(eq)
            slots = Nslot(eq)
            knowns = 0
            unknown = 0
            if (e_resolved(el)) knowns=1
            do sl=1,slots
              dk = dustkind(eq,sl) 
              if (d_resolved(dk)) then
                knowns = knowns + 1
              else
                unknown = sl
              endif  
            enddo 
            unknowns = slots+1-knowns
            text = ""
            do sl=1,slots
              dk = dustkind(eq,sl)  
              write(txt1,'(I2)') stoich(eq,sl)
              txt0 = " "
              if (d_resolved(dk)) txt0="*"
              text = trim(text)//" "//trim(txt1)//" "
     >             //trim(dust_nam(dk))//txt0
            enddo
            write(txt1,'(I2,": ")') eq 
            write(txt2,'(" (",I2,")")') unknowns 
            txt0 = " "
            if (e_resolved(el)) txt0="*"
            if (verbose>1) print*,trim(txt1)//trim(text)//" + "
     >                    //trim(elnam(el))//txt0//trim(txt2)
            if (unknowns>1) solved=.false.
            if (unknowns==1) then
              action = .true. 
              call GET_KNOWNS(eq,mat,emat,e_act,d_resolved,e_resolved,
     >                        elem,Nslot,dustkind,stoich,vec,verbose)
              if (unknown>0) then 
                dk = dustkind(eq,unknown)
                d_resolved(dk) = .true.
                mat(dk,:) = -vec(:)/REAL(stoich(eq,unknown))
                do j=1,Nind
                  el2 = Iindex(j)
                  if (mat(dk,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in1 ",trim(dust_nam(dk))
     >                      //" "//elnam(el2),REAL(mat(dk,el2))
                enddo  
              else
                e_resolved(el) = .true.
                emat(el,:) = -vec(:)
                do j=1,Nind
                  el2 = Iindex(j)
                  if (emat(el,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in2 ",trim(elnam(el))//" "
     >                      //elnam(el2),REAL(emat(el,el2))
                enddo  
              endif
              exit
            else if (unknowns>1) then 
              Nunsolved = Nunsolved+1 
              unsolved(Nunsolved) = eq
            endif  
          enddo
          if (.not.action.and.Nunsolved==0) exit
          if (.not.action.and.Nunsolved>1) then
            !--------------------------------------------
            ! ***  solve N equations with N unknowns?  *** 
            !--------------------------------------------
            Nvar1 = 0
            do i=1,NDUST
              if (.not.active(i)) cycle
              if (d_resolved(i)) cycle
              Nvar1 = Nvar1+1
              var(Nvar1) = i
            enddo  
            Nvar2 = 0            
            do el=1,NELEM
              if (.not.e_taken(el)) cycle
              if (e_resolved(el)) cycle
              Nvar2 = Nvar2+1
              var(Nvar1+Nvar2) = el
            enddo  
            !write(*,*) Nunsolved,Nvar1,Nvar2
            !write(*,*) unsolved(1:Nunsolved)
            !write(*,*) dust_nam(var(1:Nvar1)), 
     >      !           elnam(var(Nvar1+1:Nvar1+Nvar2))
            if (verbose>1) print'("solving",I2," equations with",I2,
     >                 " unknowns ",99(A12))',Nunsolved,Nvar1+Nvar2,
     >                 dust_nam(var(1:Nvar1)), 
     >                 elnam(var(Nvar1+1:Nvar1+Nvar2))
            if (Nunsolved/=Nvar1+Nvar2) then
              print*,"... is impossible"
              stop
            endif  
            if (Nunsolved==Nvar1+Nvar2) then
              DF = 0.Q0 
              FF = 0.Q0
              do i=1,Nunsolved
                eq = unsolved(i)
                call GET_KNOWNS(eq,mat,emat,e_act,d_resolved,e_resolved,
     >                     elem,Nslot,dustkind,stoich,vecs(i,:),verbose)
                slots = Nslot(eq)
                do j=1,Nvar1
                  dk = var(j)
                  do sl=1,slots
                    if (dk==dustkind(eq,sl)) then  
                      DF(i,j) = stoich(eq,sl)
                    endif  
                  enddo
                enddo
                el = elem(eq)
                do j=Nvar1+1,Nvar1+Nvar2
                  el2 = var(j)
                  if (el2==el) then
                    DF(i,j) = 1.Q0
                  endif  
                enddo
              enddo
              !do i=1,Nunsolved
              !  print'(99(1pE12.3))',(DF(i,j),j=1,Nunsolved)
              !enddo
              !--- compute inverse matrix ---
              DFsav = DF
              call QGEFA ( DF, NELEM, Nunsolved, ipvt, info )
              call QGEDI ( DF, NELEM, Nunsolved, ipvt, det, work, 1 )
              !do i=1,Nunsolved
              !  print'(99(1pE12.3))',(DF(i,j),j=1,Nunsolved)
              !enddo
              do i=1,Nvar1
                dk = var(i)
                vec(:) = 0.Q0
                do j=1,Nunsolved
                  vec = vec + DF(i,j)*vecs(j,:)
                enddo
                mat(dk,:) = -vec(:)
                d_resolved(dk) = .true.
                do j=1,Nind
                  el2 = Iindex(j)
                  if (mat(dk,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in3 ",trim(dust_nam(dk))//" "
     >                                //elnam(el2),REAL(mat(dk,el2))
                enddo  
              enddo
              do i=Nvar1+1,Nvar1+Nvar2
                el = var(i)
                vec(:) = 0.Q0
                do j=1,Nunsolved
                  vec = vec + DF(i,j)*vecs(j,:)
                enddo
                e_resolved(el) = .true.
                emat(el,:) = -vec(:)
                do j=1,Nind
                  el2 = Iindex(j)
                  if (emat(el,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in4 ",trim(elnam(el))//" "
     >                            //elnam(el2),REAL(emat(el,el2))
                enddo  
              enddo  
            endif  
          endif    
          if (itry==100) stop "*** itry==100"
        enddo
        if (.not.solved) then
          write(*,*) "*** couldn't resolve the conversion matrix."
          stop
        endif   
        do i=1,Nind
          el = Iindex(i) 
          do j=1,Ndep 
            if (is_dust(j)) then         
              dk  = Dindex(j)
              conv(j,i) = mat(dk,el)
            else
              el2 = Dindex(j) 
              conv(j,i) = emat(el2,el)
            endif
          enddo
        enddo  
        if (verbose>0) then
          print'(A24,99(A7))',"solving for ...",elnam(Iindex(1:Nind))
          do j=1,Ndep 
            if (is_dust(j)) then
              dk  = Dindex(j)
              txt = dust_nam(dk)
            else
              el  = Dindex(j) 
              txt = elnam(el)
            endif   
            print'(" conv.mat ",A14,99(0pF7.3))',trim(txt),
     >             (conv(j,i),i=1,Nind)
          enddo  
        endif  

 100    continue
        !-----------------------------------------------
        ! ***  stop iteration of parts of solution?  ***
        !-----------------------------------------------
        if (changed.or.it==1) then
          converge(:,:) = 9.Q+99
          is_esolved(:) = .false.
          is_dsolved(:) = .false.
        else if (it>3) then
          Ntrivial = 0
          dtrivial = 0
          etrivial = 0
          do ii=1,Nsolve
            Nzero = 0
            nonzero = 0
            do jj=1,Nsolve
              if (ABS(DFsav(jj,ii))<1.Q-3) then
                Nzero = Nzero+1
              else
                nonzero = jj
              endif
            enddo
            if (Nzero==Nsolve-1) then
              Ntrivial = Ntrivial+1
              dtrivial(Ntrivial) = nonzero
              etrivial(Ntrivial) = ii
            endif  
          enddo  
          ebest = 0
          dbest = 0
          cbest = 9.Q+99
          do itrivial=1,Ntrivial
            ii = etrivial(itrivial) 
            jj = dtrivial(itrivial)
            i  = act_to_elem(ii)
            j  = act_to_dust(jj) 
            el = Iindex(i)
            dk = Dindex(j)
            crit = MAXVAL(ABS(converge(it-3:it-1,i)))  
            !print'(A2,A13,3(1pE12.3))',elnam(el),trim(dust_nam(dk)),
     >      !     converge(it-1,i),crit
            if (crit<1.Q-15.and.crit<cbest) then
              cbest = crit
              ebest = i
              dbest = j
            endif
          enddo
          if (ebest.ne.0) then
            el = Iindex(ebest) 
            dk = Dindex(dbest)
            is_esolved(ebest) = .true.
            is_dsolved(dbest) = .true.
            !print*,elnam(el)//" (->"//trim(dust_nam(dk))//
     >      !       ") has converged."
          endif  
        endif  
        !-----------------------------------------------
        ! ***  fill in r.h.s. vector FF and          ***
        ! ***  determine current quality of solution ***
        !-----------------------------------------------
        qual = SQUAL(Sat0,active)
        ii = 0
        do i=1,Nind
          if (is_dsolved(i)) cycle
          ii = ii+1
          act_to_dust(ii) = i
          dk = Dindex(i)
          !FF(ii) = 1.Q0-Sat0(dk)
          FF(ii) = LOG(Sat0(dk))
        enddo 
        Nsolve = ii

        !------------------------------------------
        ! ***  compute numerical derivative DF  ***
        !------------------------------------------
        jj = 0
        do j=1,Nind
          if (is_esolved(j)) cycle
          jj = jj+1
          act_to_elem(jj) = j
          el = Iindex(j) 
          deps1 = +1.Q-3*eps(el)            ! limited by el abundance
          deps2 = -1.Q-3*eps(el)            ! limited by el abundance
          do i=1,Ndep
            if (conv(i,j)==0.Q0) cycle
            if (is_dust(i)) cycle
            el2 = Dindex(i)                 ! limited by dep. element?
            del = 1.Q-2*eps(el2)/conv(i,j)
            if (del>0.Q0) deps2=MAX(deps2,-del)
            if (del<0.Q0) deps1=MIN(deps1,-del)
            !if (verbose>1) print*,elnam(el)//" "//elnam(el2),
     >      !               REAL((/conv(i,j),deps1,deps2/))
          enddo
          deps = deps2
          if (ABS(deps1)>ABS(deps2)) deps=deps1
          xstep(:) = 0.Q0
          xstep(j) = deps
          scale(j) = eps(el)
          call SUPER(nHtot,T,xstep,eps,Sat2)
          do ii=1,Nsolve
            i  = act_to_dust(ii) 
            dk = Dindex(i)
            !DF(ii,jj) = (Sat2(dk)-Sat0(dk))/deps*scale(j)
            DF(ii,jj) = LOG(Sat0(dk)/Sat2(dk))/deps*scale(j)
          enddo  
        enddo            
        if (verbose>1) then
          print'(12x,99(A11))',elnam(Iindex(act_to_elem(1:Nsolve)))
          do ii=1,Nsolve
            i  = act_to_dust(ii) 
            dk = Dindex(i)
            print'(A12,99(1pE11.3))',dust_nam(dk),DF(ii,1:Nsolve),FF(ii)
          enddo  
        endif

        !--------------------------------
        ! ***  Newton-Raphson step dx ***
        !--------------------------------
        Fsav  = FF
        DFsav = DF
        call GAUSS16( NELEM, Nsolve, DF, dx, FF)
        !--- re-scale ---
        if (it>1) converge(it,:) = converge(it-1,:)
        do ii=1,Nsolve
          i = act_to_elem(ii) 
          el = Iindex(i) 
          dx(ii) = dx(ii)*scale(i)
          converge(it,i) = dx(ii)/eps(el)
        enddo  

        !-----------------------------------
        ! ***  limit NR step physically  ***
        !-----------------------------------
        fac = 1.Q0
        iminoff = 0
        do ii=1,Nsolve
          i  = act_to_elem(ii) 
          el = Iindex(i)
          if (eps(el)+dx(ii)<0.3*eps(el)) then
            fac2 = (-0.7*eps(el))/dx(ii)        ! eps+fac*dx = 0.3*eps
            if (verbose>0) print*,"*** limiting element "
     >                            //elnam(el),REAL(fac2)
            if (fac2<fac) then
              fac = fac2 
            endif
          endif  
        enddo
        do j=Ndep,1,-1
          del = 0.Q0 
          do ii=1,Nsolve
            i = act_to_elem(ii) 
            del = del + conv(j,i)*dx(ii)
          enddo 
          if (is_dust(j)) then
            dk = Dindex(j)
            if (ddust(dk)+del<0.Q0) then
              fac2 = (-ddust(dk)-small*dscale(dk))/del
              if (verbose>0) print*,"*** limiting dust "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2 
                iminoff = dk
              endif
            endif  
          else  
            el = Dindex(j)
            if (eps(el)+del<0.1*eps(el)) then
              fac2 = (-0.9*eps(el))/del        ! eps+fac*dx = 0.1*eps
              if (verbose>0) print*,"*** limiting element "
     >                              //elnam(el),REAL(fac2)
              if (fac2<fac) then
                fac = fac2 
              endif  
            endif
          endif  
        enddo  
        dx = dx*fac
        if (iminoff>0) then
          print*,"switch off ",dust_nam(iminoff) 
          active(iminoff) = .false.
          lastit = -99
        endif
        limited = (fac<1.Q0)

        !------------------------------------
        ! ***  apply dx to ddust and eps  ***
        !------------------------------------
        do ii=1,Nsolve
          i = act_to_elem(ii) 
          el = Iindex(i)
          eps(el) = eps(el) + dx(ii)     ! direct effect
        enddo  
        do j=1,Ndep
          del = 0.Q0 
          do ii=1,Nsolve
            i = act_to_elem(ii) 
            del = del + conv(j,i)*dx(ii) ! effect of indep.element i
          enddo 
          if (is_dust(j)) then
            dk = Dindex(j)
            ddust(dk) = ddust(dk) + del
          else  
            el = Dindex(j)
            eps(el) = eps(el) + del
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
        worst = 0.d0
        do i=1,NELEM
          worst = MAX(worst,ABS(1.Q0-check(i)/eps00(i)))
        enddo
        if (verbose>1) write(*,*) "element conservation error:",worst
        if (worst>1.Q-8) stop

        xstep(:) = 0.Q0
        call SUPER(nHtot,T,xstep,eps,Sat0)
        qual = SQUAL(Sat0,active)
        print'("it=",I3," qual=",1pE12.4)',it,qual
        if (qual<1.Q-20) exit
        if (verbose>0) read(*,'(a1)') char1

      enddo  
      Sat = Sat0

      call CPU_TIME(time1)
      if (it.lt.itmax) then
        write(*,'("EQUIL_COND converged after ",I3," iter, time =",
     >            0pF7.3," CPU sec.")') it,time1-time0 
      else
        write(*,'("*** EQUIL_COND failed after ",I3," iter,  time =",
     >            0pF9.4," CPU sec.")') it,time1-time0 
        fail_tot = fail_tot+1
        stop
      endif   

      !----------------------------------
      ! ***  save result to database  ***
      !----------------------------------
      if (qual<1.Q-10) then
        !if (ddust(iTiO2)>0.Q0.and.ddust(iMgTi2O5)>0.Q0) stop
        call PUT_DATA(nHtot,T,eps,ddust,qread,iread,active)
      endif  
      it_tot = it_tot + it

      end
            

!-------------------------------------------------------------------------
      subroutine SUPER(nHtot,T,xx,eps,Sat)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,
     >                    dust_nam,elnam
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      use CONVERSION,ONLY: Ndep,Nind,Dindex,Iindex,is_dust,conv
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in)  :: nHtot,T
      real(kind=qp),intent(in) :: xx(NELEM),eps(NELEM)
      real(kind=qp),intent(out) :: Sat(NDUST)
      real(kind=qp) :: eps1(NELEM),dx
      integer :: i,j,el

      !-------------------------------------------
      ! ***  compute remaining gas abundances  ***
      !-------------------------------------------
      eps1 = eps
      do i=1,Nind
        el = Iindex(i)
        eps1(el) = eps1(el) + xx(i)  ! direct effect
      enddo  
      do j=1,Ndep
        if (is_dust(j)) cycle
        dx = 0.Q0 
        do i=1,Nind
          dx = dx + conv(j,i)*xx(i)  ! effect of indep.element i on el j
        enddo 
        el = Dindex(j)
        eps1(el) = eps1(el) + dx
      enddo

      do el=1,NELEM
        if (eps1(el).le.0.Q0) then
          write(*,*) "*** negative el.abund. SUPER",elnam(el),eps1(el)
          stop
        endif  
      enddo
      
      !----------------------------------------------
      ! ***  compute chemistry & supersaturation  ***
      !----------------------------------------------
      call GGCHEM(nHtot,T,eps1,.true.,0)
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
      logical,intent(in) :: active(NDUST)
      integer :: i

      qual = 0.d0
      do i=1,NDUST
        if (active(i).or.(Sat(i).gt.1.Q0)) then
          qual = qual + (1.Q0-Sat(i))**2
         !qual = qual + (Sat(i)-1.Q0/Sat(i))**2
         !qual = qual + LOG(Sat(i))**2
        endif  
      enddo
      SQUAL = qual
      end

!-------------------------------------------------------------------------
      subroutine VAPORIZE(i,ddust,eps,dscale)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp),intent(in) :: dscale(NDUST)
      real(kind=qp) :: del
      integer :: j,el
      
      del = ddust(i)
      print*," ==>  vaporize "//trim(dust_nam(i)),REAL(del/dscale(i))
      ddust(i) = 0.Q0
      do j=1,dust_nel(i)
        el = dust_el(i,j)
        eps(el) = eps(el) + del*dust_nu(i,j)    
      enddo
      end

!-------------------------------------------------------------------------
      subroutine TRANSFORM(i1,i2,del,fac,ddust,eps,dscale)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam,
     >                    eps0
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i1,i2
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp),intent(in) :: del,fac,dscale(NDUST)
      real(kind=qp) :: check(NELEM),worst
      integer :: i,j,el
      
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
      worst = 0.d0
      do i=1,NELEM
        worst = MAX(worst,ABS(1.Q0-check(i)/eps0(i)))
      enddo
      write(*,*) "element conservation error:",worst
      if (worst>1.Q-8) stop

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
