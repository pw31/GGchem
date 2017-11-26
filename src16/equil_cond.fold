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
      use PARAMETERS,ONLY: Tfast,CzuOvari
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nam,dust_nel,dust_nu,dust_el,
     >                    eps0,elnam,elcode
      use CONVERSION,ONLY: Nind,Ndep,Iindex,Dindex,is_dust,conv
      use EXCHANGE,ONLY: Fe,Mg,Si,Al,Ca,Ti,C,O,S,Na,Cl,H,Li,Mn,W,Ni,Cr,
     >                   Kalium=>K,Zr,V,itransform,ieqcond
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
      real(kind=qp),dimension(NELEM) :: scale,save_eps
      real(kind=qp),dimension(NDUST) :: ddustread,dscale,pot,dust_save
      real(kind=qp),dimension(NDUST) :: Sat0,Sat1,Sat2,save_ddust
      real(kind=qp),dimension(NELEM,NELEM) :: DF,DFsav,emat,vecs
      real(kind=qp),dimension(NDUST,NELEM) :: mat
      real(kind=qp) :: worst,xmin,Smax,qual,SQUAL,del,save_qual
      real(kind=qp) :: turnon,turnoff,maxon,minoff,fac,fac2,amount
      real(kind=qp) :: deps1,deps2,deps
      real(kind=qp) :: det(2),converge(500,NELEM),crit,cbest
      real(kind=qp) :: small=1.Q-30
      integer,parameter :: itmax=200
      integer,dimension(NELEM) :: elem,Nslot
      integer,dimension(NELEM,NDUST) :: dustkind,stoich
      integer :: it,i,j,el,el2,Nact,Nact_read,Neq,slots,sl,dk,eq
      integer :: itry,knowns,unknowns,unknown,ii,jj,lastit,laston
      integer :: imaxon,iminoff,info,ipvt(NELEM),e_num(NELEM)
      integer :: Nunsolved,unsolved(NELEM),Nvar1,Nvar2,var(NELEM)
      integer :: Nsolve,ebest,dbest,nonzero,itrivial,iread,ioff
      integer :: ifail,Nact_save,Nall,imax,swap
      integer :: act_to_elem(NELEM),act_to_dust(NELEM)
      integer :: Nzero,Ntrivial,etrivial(NELEM),dtrivial(NELEM)
      logical,dimension(NELEM) :: e_resolved,e_act,e_taken,is_esolved
      logical,dimension(0:NDUST) :: active,act_read,act_old
      logical,dimension(0:NDUST) :: save_active,active_save
      logical,dimension(NDUST) :: is_dsolved,d_resolved
      logical :: action,changed,solved,limited,ok,conserved,exch
      character(len=1) :: char1,txt0
      character(len=2) :: rem
      character(len=500) :: txt,txt1,txt2,text,filename
      logical,save :: firstCall=.true.
      integer,save :: iAl2O3=0,iFe=0,iFeS=0,iNa2SiO3=0,iMgSiO3=0
      integer,save :: iMg2SiO4=0,iTi4O7=0,iCaSiO3=0,iCaMgSi2O6=0
      integer,save :: iNaAlSi3O8=0,iMgAl2O4=0,iCaTiO3=0,iSiO=0,iSiO2=0
      integer,save :: iTiO2=0,iMgTi2O5=0,iSiC=0,iCaS=0,iFe2SiO4=0,iFeO=0
      integer,save :: iNaCl=0,iKCl=0,iKAlSi3O8=0,iFe_l=0,iH2O=0,iH2O_l=0
      integer,save :: iFeS_l=0,iNaCl_l=0,iTiO2_l=0,iSiO2_l=0
      integer,save :: iNa2SiO3_l=0,iMgAl2O4_l=0,iMg2SiO4_l=0,iMgSiO3_l=0
      integer,save :: iAl2O3_l=0,iCaAl2Si2O8=0,iC=0,iTiC=0,iFe2O3=0
      integer,save :: iMgO=0,iNa=0,iS=0,iMgS=0,iLiCl=0,iSiS2=0,iFeS2=0
      integer,save :: iH2SO4_l=0,iNa2S=0,iAlCl3=0,iNH3=0,iCaO=0,iNa_l=0
      integer,save :: iKCl_l=0,iCaCl2_l=0,iLiCl_l=0,iTi4O7_l=0,iFeO_l=0
      integer,save :: iCH4=0,iCaCl2=0,iLiH=0,iLiH_l=0,iMgTi2O5_l=0
      integer,save :: iMgO_l=0,iAlCl3_l=0,iMgTiO3=0,iMgTiO3_l=0,iCaO_l=0
      integer,save :: iS_l=0,iK2SiO3=0,iK2SiO3_l=0,iTiC_l=0,iTi=0
      integer,save :: iTi_l=0,iTiO=0,iTiO_l=0,iSiS2_l=0,iLiOH=0
      integer,save :: iLiOH_l=0,iMnS=0,iW=0,iW_l=0,iZrO2=0,iZrSiO4=0
      integer,save :: iVO=0,iV2O3=0,iCr=0
      integer,save :: iNi=0,iNi_l,iNi3S2=0,iFe3O4=0,iKMg3AlSi3O12H2=0
      integer,save :: iKFe3AlSi3O12H2=0,iMg3Si2O9H4=0,iFe3Si2O9H4=0
      integer,save :: iMgCr2O4=0,iCr2O3=0,iMn3Al2Si3O12=0,iMn2SiO4=0
      integer,save :: iKAlSi2O6=0,iCa3Al2Si3O12=0,iFeAl2SiO7H2=0
      integer,save :: iNaMg3AlSi3O12H2=0,iNaAlSiO4=0,iCa2MgSi2O7=0
      integer,save :: iCa2Al2SiO7=0
      integer,save :: iCaTiSiO5=0,iNaAlSi2O6=0,iKAlSiO4=0,iMg3Si4O12H2=0
      integer,save :: iCa3MgSi2O8=0,iCaMgSiO4=0
      integer,save :: it_tot=0, sit_tot=0, fail_tot=0
      real*8 :: time0,time1,qread

      if (firstCall) then
        do i=1,NDUST
          if (dust_nam(i).eq.'Al2O3[s]')      iAl2O3=i 
          if (dust_nam(i).eq.'Al2O3[l]')      iAl2O3_l=i 
          if (dust_nam(i).eq.'Fe2O3[s]')      iFe2O3=i 
          if (dust_nam(i).eq.'SiO[s]')        iSiO=i 
          if (dust_nam(i).eq.'SiO2[s]')       iSiO2=i
          if (dust_nam(i).eq.'SiO2[l]')       iSiO2_l=i 
          if (dust_nam(i).eq.'SiS2[s]')       iSiS2=i
          if (dust_nam(i).eq.'SiS2[l]')       iSiS2_l=i
          if (dust_nam(i).eq.'Fe[s]')         iFe=i 
          if (dust_nam(i).eq.'Fe[l]')         iFe_l=i 
          if (dust_nam(i).eq.'FeS[s]')        iFeS=i 
          if (dust_nam(i).eq.'FeS[l]')        iFeS_l=i 
          if (dust_nam(i).eq.'FeS2[s]')       iFeS2=i 
          if (dust_nam(i).eq.'Na2SiO3[s]')    iNa2SiO3=i
          if (dust_nam(i).eq.'Na2SiO3[l]')    iNa2SiO3_l=i 
          if (dust_nam(i).eq.'MgSiO3[s]')     iMgSiO3=i 
          if (dust_nam(i).eq.'MgSiO3[l]')     iMgSiO3_l=i 
          if (dust_nam(i).eq.'Mg2SiO4[s]')    iMg2SiO4=i 
          if (dust_nam(i).eq.'Mg2SiO4[l]')    iMg2SiO4_l=i 
          if (dust_nam(i).eq.'CaSiO3[s]')     iCaSiO3=i 
          if (dust_nam(i).eq.'CaTiO3[s]')     iCaTiO3=i 
          if (dust_nam(i).eq.'CaMgSi2O6[s]')  iCaMgSi2O6=i
          if (dust_nam(i).eq.'CaAl2Si2O8[s]') iCaAl2Si2O8=i
          if (dust_nam(i).eq.'Ca2Al2SiO7[s]') iCa2Al2SiO7=i
          if (dust_nam(i).eq.'NaAlSi3O8[s]')  iNaAlSi3O8=i
          if (dust_nam(i).eq.'MgAl2O4[s]')    iMgAl2O4=i
          if (dust_nam(i).eq.'MgAl2O4[l]')    iMgAl2O4_l=i
          if (dust_nam(i).eq.'Ti4O7[s]')      iTi4O7=i
          if (dust_nam(i).eq.'Ti4O7[l]')      iTi4O7_l=i
          if (dust_nam(i).eq.'TiO2[s]')       iTiO2=i
          if (dust_nam(i).eq.'TiO2[l]')       iTiO2_l=i
          if (dust_nam(i).eq.'MgTi2O5[s]')    iMgTi2O5=i
          if (dust_nam(i).eq.'MgTi2O5[l]')    iMgTi2O5_l=i
          if (dust_nam(i).eq.'MgTiO3[s]')     iMgTiO3=i
          if (dust_nam(i).eq.'MgTiO3[l]')     iMgTiO3_l=i
          if (dust_nam(i).eq.'SiC[s]')        iSiC=i
          if (dust_nam(i).eq.'CaS[s]')        iCaS=i
          if (dust_nam(i).eq.'Fe2SiO4[s]')    iFe2SiO4=i
          if (dust_nam(i).eq.'Fe3O4[s]')      iFe3O4=i
          if (dust_nam(i).eq.'FeO[s]')        iFeO=i
          if (dust_nam(i).eq.'FeO[l]')        iFeO_l=i
          if (dust_nam(i).eq.'NaCl[s]')       iNaCl=i
          if (dust_nam(i).eq.'NaCl[l]')       iNaCl_l=i
          if (dust_nam(i).eq.'LiCl[s]')       iLiCl=i
          if (dust_nam(i).eq.'LiCl[l]')       iLiCl_l=i
          if (dust_nam(i).eq.'KCl[s]')        iKCl=i
          if (dust_nam(i).eq.'KCl[l]')        iKCl_l=i
          if (dust_nam(i).eq.'KAlSi3O8[s]')   iKAlSi3O8=i
          if (dust_nam(i).eq.'KAlSi2O6[s]')   iKAlSi2O6=i
          if (dust_nam(i).eq.'K2SiO3[s]')     iK2SiO3=i
          if (dust_nam(i).eq.'K2SiO3[l]')     iK2SiO3_l=i
          if (dust_nam(i).eq.'H2O[s]')        iH2O=i
          if (dust_nam(i).eq.'H2O[l]')        iH2O_l=i
          if (dust_nam(i).eq.'H2SO4[l]')      iH2SO4_l=i
          if (dust_nam(i).eq.'C[s]')          iC=i
          if (dust_nam(i).eq.'S[s]')          iS=i
          if (dust_nam(i).eq.'S[l]')          iS_l=i
          if (dust_nam(i).eq.'Ti[s]')         iTi=i
          if (dust_nam(i).eq.'Ti[l]')         iTi_l=i
          if (dust_nam(i).eq.'TiC[s]')        iTiC=i
          if (dust_nam(i).eq.'TiC[l]')        iTiC_l=i
          if (dust_nam(i).eq.'MgO[s]')        iMgO=i
          if (dust_nam(i).eq.'MgO[l]')        iMgO_l=i
          if (dust_nam(i).eq.'CaO[s]')        iCaO=i
          if (dust_nam(i).eq.'CaO[l]')        iCaO_l=i
          if (dust_nam(i).eq.'MgS[s]')        iMgS=i
          if (dust_nam(i).eq.'Na[s]')         iNa=i 
          if (dust_nam(i).eq.'Na[l]')         iNa_l=i 
          if (dust_nam(i).eq.'Na2S[s]')       iNa2S=i
          if (dust_nam(i).eq.'AlCl3[s]')      iAlCl3=i 
          if (dust_nam(i).eq.'AlCl3[l]')      iAlCl3_l=i 
          if (dust_nam(i).eq.'NH3[s]')        iNH3=i
          if (dust_nam(i).eq.'CaCl2[s]')      iCaCl2=i
          if (dust_nam(i).eq.'CaCl2[l]')      iCaCl2_l=i
          if (dust_nam(i).eq.'CH4[s]')        iCH4=i
          if (dust_nam(i).eq.'LiH[s]')        iLiH=i
          if (dust_nam(i).eq.'LiH[l]')        iLiH_l=i
          if (dust_nam(i).eq.'TiO[s]')        iTiO=i
          if (dust_nam(i).eq.'TiO[l]')        iTiO_l=i
          if (dust_nam(i).eq.'LiOH[s]')       iLiOH=i
          if (dust_nam(i).eq.'LiOH[l]')       iLiOH_l=i
          if (dust_nam(i).eq.'MnS[s]')        iMnS=i 
          if (dust_nam(i).eq.'W[s]')          iW=i 
          if (dust_nam(i).eq.'W[l]')          iW_l=i 
          if (dust_nam(i).eq.'Ni[s]')         iNi=i 
          if (dust_nam(i).eq.'Ni[l]')         iNi_l=i 
          if (dust_nam(i).eq.'ZrO2[s]')       iZrO2=i 
          if (dust_nam(i).eq.'ZrSiO4[s]')     iZrSiO4=i 
          if (dust_nam(i).eq.'KAlSiO4[s]')    iKAlSiO4=i
          if (dust_nam(i).eq.'Ni3S2[s]')      iNi3S2=i
          if (dust_nam(i).eq.'Mg3Si2O9H4[s]') iMg3Si2O9H4=i
          if (dust_nam(i).eq.'Mg3Si4O12H2[s]') iMg3Si4O12H2=i
          if (dust_nam(i).eq.'Fe3Si2O9H4[s]') iFe3Si2O9H4=i
          if (dust_nam(i).eq.'MgCr2O4[s]')    iMgCr2O4=i
          if (dust_nam(i).eq.'Cr[s]')         iCr=i
          if (dust_nam(i).eq.'Cr2O3[s]')      iCr2O3=i
          if (dust_nam(i).eq.'Mn2SiO4[s]')    iMn2SiO4=i
          if (dust_nam(i).eq.'NaAlSi2O6[s]')  iNaAlSi2O6=i
          if (dust_nam(i).eq.'NaAlSiO4[s]')   iNaAlSiO4=i
          if (dust_nam(i).eq.'CaMgSiO4[s]')   iCaMgSiO4=i
          if (dust_nam(i).eq.'CaTiSiO5[s]')   iCaTiSiO5=i
          if (dust_nam(i).eq.'Ca2MgSi2O7[s]') iCa2MgSi2O7=i
          if (dust_nam(i).eq.'Ca3MgSi2O8[s]') iCa3MgSi2O8=i
          if (dust_nam(i).eq.'Ca3Al2Si3O12[s]') iCa3Al2Si3O12=i
          if (dust_nam(i).eq.'Mn3Al2Si3O12[s]') iMn3Al2Si3O12=i
          if (dust_nam(i).eq.'FeAl2SiO7H2[s]')  iFeAl2SiO7H2=i
          if (dust_nam(i).eq.'KMg3AlSi3O12H2[s]') iKMg3AlSi3O12H2=i
          if (dust_nam(i).eq.'KFe3AlSi3O12H2[s]') iKFe3AlSi3O12H2=i
          if (dust_nam(i).eq.'NaMg3AlSi3O12H2[s]') iNaMg3AlSi3O12H2=i
          if (dust_nam(i).eq.'VO[s]')         iVO=i
          if (dust_nam(i).eq.'V2O3[s]')       iV2O3=i
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
      Nact = 0
      if (qread.lt.0.5.and.(.not.CzuOvari)) then
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
        !if (qread>1.Q-3.and.iread==50655) verbose=2
        if (verbose>0) then
          write(*,'(" ... using database entry (",I6,
     >          ") qual=",1pE15.7)') iread,qread
          write(*,*) trim(text)
        endif  
      endif
      !verbose=2
  
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
        write(*,*) "element conservation error 1:",worst
        write(*,*) "initial gas fractions ..."
        do i=1,NELEM
          if (elcode(i)==0) cycle
          print'(3x,A2,2(1pE15.6))',elnam(i),eps(i),eps(i)/eps00(i)
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
        dscale(i) = xmin                        ! max dust abundances
      enddo   

      call GGCHEM(nHtot,T,eps,.false.,0)  ! one call from scratch
      xstep(:) = 0.Q0             
      call SUPER(nHtot,T,xstep,eps,Sat0)
      qual = SQUAL(Sat0,active)
      print'("it=",I3," qual=",1pE12.4)',0,qual
      act_old = active
      lastit = -99
      iminoff = 0
      limited = .false.
      ifail = 0

      do it=1,itmax
        
        !-------------------------
        ! ***  fall-back copy  ***
        !-------------------------
        !save_active = active
        !save_eps = eps
        !save_ddust = ddust

        !---------------------------------------
        ! ***  selection of solids to solve  ***
        !---------------------------------------
        changed = .false.
        Smax = maxval(Sat0)
        ioff = 0
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
              if (eps(el)/REAL(dust_nu(i,j),kind=qp).lt.xmin) then    
                xmin = eps(el)/REAL(dust_nu(i,j),kind=qp)
                pot(i) = 1.Q0/REAL(dust_nu(i,j),kind=qp)
              endif
            enddo
            Sat1(i)=Sat0(i)**(1.Q0/pot(i))
            !print'(A20,0pF8.3)',dust_nam(i),pot(i)
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
              !if (i==iMg2SiO4) turnon=turnon*10.0
              if (turnon>maxon.and.(.not.limited)) then
                !if (.not.(Nact>0.and.qual>1.E+10)) then 
                  maxon  = turnon
                  imaxon = i
                !endif  
              endif  
            endif  
          enddo  
          if (verbose>0) print'("limited=",L1,
     >                   "  Smax=",1pE10.3,2x,A18)',
     >                   limited,Smax,dust_nam(imax)
          if (verbose>0.and.maxon>0.Q0) print'("  maxon =",
     >                   1pE10.2,2x,A18)',maxon,dust_nam(imaxon)
          active_save = active 
          Nact_save = Nact
          if (maxon>0.0*MAX(Smax-1.Q0,0.Q0)) then
            if (imaxon.ne.iminoff) then 
              active(imaxon) = .true.
              print*,"switch on ",trim(dust_nam(imaxon))
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
          eps_save = eps
          dust_save = ddust
          ioff = 0
          ok = .true.
          if (active(iVO).and.active(iV2O3)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iVO)>Sat0(iV2O3)) then
              ioff = iV2O3
              active(iV2O3) = .false.  
              amount = ddust(iV2O3)/1.Q0
              call TRANSFORM(iV2O3,iVO,amount,2.0*1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else
              ioff = iVO
              active(iVO) = .false.  
              amount = ddust(iVO)/1.Q0
              call TRANSFORM(iVO,iV2O3,amount,0.5*1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            eps(V) = eps_save(V)
          endif  
          if (active(iMg3Si2O9H4).and.active(iMg3Si4O12H2).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMg3Si2O9H4)>Sat0(iMg3Si4O12H2)) then
              ioff = iMg3Si4O12H2
              active(iMg3Si4O12H2) = .false.  
              amount = ddust(iMg3Si4O12H2)/3.Q0
              call TRANSFORM(iMg3Si4O12H2,iMg3Si2O9H4,amount,
     >                      0.5Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMg3Si4O12H2,iMgSiO3,amount,
     >                      4.5Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMg3Si4O12H2,iMg2SiO4,amount,
     >                     -1.5Q0*3.Q0,ddust,eps,dscale,active,ok)
            else
              ioff = i
              active(iMg3Si2O9H4) = .false.  
              amount = ddust(iMg3Si2O9H4)/3.Q0
              call TRANSFORM(iMg3Si2O9H4,iMg3Si4O12H2,amount,
     >                       2.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMg3Si2O9H4,iMgSiO3,amount,
     >                      -9.Q0*2.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMg3Si2O9H4,iMg2SiO4,amount,
     >                       3.Q0*2.Q0,ddust,eps,dscale,active,ok)
            endif   
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif  
          if (active(iNaAlSiO4).and.active(iNaAlSi2O6).and.
     >        active(iMg2SiO4).and.active(iMgSiO3)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNaAlSiO4)>Sat0(iNaAlSi2O6)) then
              ioff = iNaAlSi2O6
              active(iNaAlSi2O6) = .false.  
              amount = ddust(iNaAlSi2O6)/3.Q0
              call TRANSFORM(iNaAlSi2O6,iNaAlSiO4,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSi2O6,iMgSiO3,amount,
     >                       2.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSi2O6,iMg2SiO4,amount,
     >                      -1.Q0*3.Q0,ddust,eps,dscale,active,ok)
            else
              ioff = iNaAlSiO4
              active(iNaAlSiO4) = .false.  
              amount = ddust(iNaAlSiO4)/3.Q0
              call TRANSFORM(iNaAlSiO4,iNaAlSi2O6,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSiO4,iMgSiO3,amount,
     >                      -2.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSiO4,iMg2SiO4,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Na),eps(Al),eps(Mg),eps(Si)
            !print*,eps_save(Na),eps_save(Al),eps_save(Mg),
     >      !       eps_save(Si)
            eps(Na) = eps_save(Na)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif  
          if (active(iNaMg3AlSi3O12H2).and.active(iNaAlSiO4).and.
     >        active(iMg2SiO4).and.active(iFe2SiO4).and.
     >        active(iFe)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNaMg3AlSi3O12H2)>Sat0(iNaAlSiO4)) then
              ioff = iNaAlSiO4
              active(iNaAlSiO4) = .false.  
              amount = ddust(iNaAlSiO4)/4.Q0
              call TRANSFORM(iNaAlSiO4,iNaMg3AlSi3O12H2,amount,
     >                       1.Q0*4.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSiO4,iMg2SiO4,amount,
     >                     -1.5Q0*4.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSiO4,iFe,amount,
     >                       1.Q0*4.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSiO4,iFe2SiO4,amount,
     >                     -0.5Q0*4.Q0,ddust,eps,dscale,active,ok)
            else
              ioff = iNaMg3AlSi3O12H2
              active(iNaMg3AlSi3O12H2) = .false.  
              amount = ddust(iNaMg3AlSi3O12H2)/4.Q0
              call TRANSFORM(iNaMg3AlSi3O12H2,iNaAlSiO4,amount,
     >                       1.Q0*4.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaMg3AlSi3O12H2,iMg2SiO4,amount,
     >                      1.5Q0*4.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaMg3AlSi3O12H2,iFe,amount,
     >                      -1.Q0*4.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaMg3AlSi3O12H2,iFe2SiO4,amount,
     >                      0.5Q0*4.Q0,ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Na),eps(Al),eps(Mg),eps(Si),eps(Fe)
            !print*,eps_save(Na),eps_save(Al),eps_save(Mg),
     >      !       eps_save(Si),eps_save(Fe)
            eps(Na) = eps_save(Na)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Fe) = eps_save(Fe)
          endif  
          if (active(iCa3MgSi2O8).and.active(iCaMgSiO4).and.
     >        active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaMgSiO4)>Sat0(iCa3MgSi2O8)) then
              ioff = iCa3MgSi2O8
              active(iCa3MgSi2O8) = .false.  
              amount = ddust(iCa3MgSi2O8)/2.Q0
              call TRANSFORM(iCa3MgSi2O8,iCaMgSiO4,amount,
     >                    3.Q0*2.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iCa3MgSi2O8,iMg2SiO4,amount,
     >                   -1.Q0*2.Q0,ddust,eps,dscale,active,ok)
            else
              ioff = iCaMgSiO4
              active(iCaMgSiO4) = .false.  
              amount = ddust(iCaMgSiO4)/2.Q0
              call TRANSFORM(iCaMgSiO4,iCa3MgSi2O8,amount,
     >                    1.Q0/3.Q0*2.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaMgSiO4,iMg2SiO4,amount,
     >                    1.Q0/3.Q0*2.Q0,ddust,eps,dscale,active,ok)
            endif
            eps(Ca) = eps_save(Ca)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif     
          if (active(iCa3Al2Si3O12).and.active(iCaMgSi2O6).and.
     >        active(iFeAl2SiO7H2).and.active(iFe2SiO4).and.
     >        active(iFe3O4).and.active(iMg3Si2O9H4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCa3Al2Si3O12)>Sat0(iCaMgSi2O6)) then
              if (ddust(iCaMgSi2O6)/3<ddust(iFeAl2SiO7H2)) then
                ioff = iCaMgSi2O6
                active(iCaMgSi2O6) = .false.  
                amount = ddust(iCaMgSi2O6)/5.Q0
                call TRANSFORM(iCaMgSi2O6,iCa3Al2Si3O12,amount,
     >                    1.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iCaMgSi2O6,iFeAl2SiO7H2,amount,
     >                   -1.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iCaMgSi2O6,iFe2SiO4,amount,
     >                    2.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iCaMgSi2O6,iFe3O4,amount,
     >                   -1.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iCaMgSi2O6,iMg3Si2O9H4,amount,
     >                    1.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
              else
                ioff = iFeAl2SiO7H2
                active(iFeAl2SiO7H2) = .false.  
                amount = ddust(iFeAl2SiO7H2)/5.Q0
                call TRANSFORM(iFeAl2SiO7H2,iCa3Al2Si3O12,amount,
     >                         1.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iFeAl2SiO7H2,iCaMgSi2O6,amount,
     >                        -3.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iFeAl2SiO7H2,iFe2SiO4,amount,
     >                         2.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iFeAl2SiO7H2,iFe3O4,amount,
     >                        -1.Q0*5.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iFeAl2SiO7H2,iMg3Si2O9H4,amount,
     >                         1.Q0*5.Q0,ddust,eps,dscale,active,ok)
              endif  
            else
              ioff = iCa3Al2Si3O12
              active(iCa3Al2Si3O12) = .false.  
              amount = ddust(iCa3Al2Si3O12)/5.Q0
              call TRANSFORM(iCa3Al2Si3O12,iCaMgSi2O6,amount,
     >                       3.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iCa3Al2Si3O12,iFeAl2SiO7H2,amount,
     >                       1.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iCa3Al2Si3O12,iFe2SiO4,amount,
     >                      -2.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iCa3Al2Si3O12,iFe3O4,amount,
     >                       1.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iCa3Al2Si3O12,iMg3Si2O9H4,amount,
     >                      -1.Q0*5.Q0,ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Ca),eps(Al),eps(Mg),eps(Si),eps(Fe)
            !print*,eps_save(Ca),eps_save(Al),eps_save(Mg),
     >      !       eps_save(Si),eps_save(Fe)
            eps(Ca) = eps_save(Ca)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Fe) = eps_save(Fe)
          endif  
          if (active(iKAlSi3O8).and.active(iKAlSi2O6).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iKAlSi3O8)>Sat0(iKAlSi2O6)) then
              ioff = iKAlSi2O6
              active(iKAlSi2O6) = .false.  
              amount = ddust(iKAlSi2O6)/3.Q0
              call TRANSFORM(iKAlSi2O6,iKAlSi3O8,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi2O6,iMgSiO3,amount,
     >                      -2.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi2O6,iMg2SiO4,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
            else  
              ioff = iKAlSi3O8
              active(iKAlSi3O8) = .false.  
              amount = ddust(iKAlSi3O8)/3.Q0
              call TRANSFORM(iKAlSi3O8,iKAlSi2O6,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi3O8,iMgSiO3,amount,
     >                      -2.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi3O8,iMg2SiO4,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Kalium),eps(Mg),eps(Al),eps(Si)
            !print*,eps_save(Kalium),eps_save(Mg),eps_save(Al),
     >      !       eps_save(Si)
            eps(Kalium) = eps_save(Kalium)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif  
          if (active(iMn3Al2Si3O12).and.active(iMn2SiO4).and.
     >        active(iCaAl2Si2O8).and.active(iCaMgSi2O6).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMn3Al2Si3O12)>Sat0(iMn2SiO4)) then
              ioff = iMn2SiO4
              active(iMn2SiO4) = .false.  
              amount = ddust(iMn2SiO4)/5.Q0
              call TRANSFORM(iMn2SiO4,iMn3Al2Si3O12,amount,
     >                  2.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn2SiO4,iCaAl2Si2O8,amount,
     >                 -2.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn2SiO4,iCaMgSi2O6,amount,
     >                  2.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn2SiO4,iMgSiO3,amount,
     >                 -4.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn2SiO4,iMg2SiO4,amount,
     >                  1.Q0/3.Q0*5.Q0,ddust,eps,dscale,active,ok)
            else
              ioff = iMn3Al2Si3O12
              active(iMn3Al2Si3O12) = .false.  
              amount = ddust(iMn3Al2Si3O12)/5.Q0
              call TRANSFORM(iMn3Al2Si3O12,iMn2SiO4,amount,
     >                  3.Q0/2.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn3Al2Si3O12,iCaAl2Si2O8,amount,
     >                       1.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn3Al2Si3O12,iCaMgSi2O6,amount,
     >                      -1.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn3Al2Si3O12,iMgSiO3,amount,
     >                       2.Q0*5.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iMn3Al2Si3O12,iMg2SiO4,amount,
     >                 -1.Q0/2.Q0*5.Q0,ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Mn),eps(Al),eps(Mg),eps(Si),eps(Ca)
            !print*,eps_save(Mn),eps_save(Al),eps_save(Mg),
     >      !       eps_save(Si),eps_save(Ca)
            eps(Mn) = eps_save(Mn)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Ca) = eps_save(Ca)
            eps(O)  = eps_save(O)
          endif
          if (active(iKMg3AlSi3O12H2).and.active(iKFe3AlSi3O12H2).and.
     >        active(iMg3Si2O9H4).and.active(iFe3Si2O9H4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iKMg3AlSi3O12H2)>Sat0(iKFe3AlSi3O12H2)) then
              ioff = iKFe3AlSi3O12H2
              active(iKFe3AlSi3O12H2) = .false.  
              amount = ddust(iKFe3AlSi3O12H2)/3.Q0
              call TRANSFORM(iKFe3AlSi3O12H2,iKMg3AlSi3O12H2,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKFe3AlSi3O12H2,iMg3Si2O9H4,amount,
     >                      -1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKFe3AlSi3O12H2,iFe3Si2O9H4,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
            else
              ioff = iKMg3AlSi3O12H2
              active(iKMg3AlSi3O12H2) = .false.  
              amount = ddust(iKMg3AlSi3O12H2)/3.Q0
              call TRANSFORM(iKMg3AlSi3O12H2,iKFe3AlSi3O12H2,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iMg3Si2O9H4,amount,
     >                       1.Q0*3.Q0,ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iFe3Si2O9H4,amount,
     >                      -1.Q0*3.Q0,ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Kalium),eps(Mg),eps(Al),eps(Si),eps(Fe)
            !print*,eps_save(Kalium),eps_save(Mg),eps_save(Al),
     >      !       eps_save(Si),eps_save(Fe)
            eps(Kalium) = eps_save(Kalium)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Fe) = eps_save(Fe)
            eps(O)  = eps_save(O)
            eps(H)  = eps_save(H)
          endif  
          if (active(iCr).and.active(iCr2O3)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCr)>Sat0(iCr2O3)) then
              ioff = iCr2O3
              active(iCr2O3) = .false.  
              amount = ddust(iCr2O3)/1.Q0
              call TRANSFORM(iCr2O3,iCr,amount,2.Q0*1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iCr
              active(iCr) = .false.  
              amount = ddust(iCr)/1.Q0
              call TRANSFORM(iCr,iCr2O3,amount,0.5Q0*1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif 
            !print*,eps(Cr),eps(Mg),eps(Si),eps(O)
            !print*,eps_save(Cr),eps_save(Mg),eps_save(Si),eps_save(O)
            eps(Cr) = eps_save(Cr)
          endif    
          if (active(iMgCr2O4).and.active(iCr2O3).and.
     >        active(iMg2SiO4).and.active(iMgSiO3)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgCr2O4)>Sat0(iCr2O3)) then
              ioff = iCr2O3
              active(iCr2O3) = .false.  
              amount = ddust(iCr2O3)/3.Q0
              call TRANSFORM(iCr2O3,iMgCr2O4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCr2O3,iMgSiO3,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCr2O3,iMg2SiO4,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgCr2O4
              active(iMgCr2O4) = .false.  
              amount = ddust(iMgCr2O4)/3.Q0
              call TRANSFORM(iMgCr2O4,iCr2O3,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgCr2O4,iMgSiO3,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgCr2O4,iMg2SiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif 
            !print*,eps(Cr),eps(Mg),eps(Si),eps(O)
            !print*,eps_save(Cr),eps_save(Mg),eps_save(Si),eps_save(O)
            eps(Cr) = eps_save(Cr)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif    
          if (active(iKMg3AlSi3O12H2).and.active(iKAlSiO4).and.
     >        active(iMg2SiO4).and.active(iFe2SiO4).and.
     >        active(iFe)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iKMg3AlSi3O12H2)>Sat0(iKAlSiO4)) then
              ioff = iKAlSiO4
              active(iKAlSiO4) = .false.  
              amount = ddust(iKAlSiO4)/4.Q0
              call TRANSFORM(iKAlSiO4,iKMg3AlSi3O12H2,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSiO4,iMg2SiO4,amount,-1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSiO4,iFe2SiO4,amount,-0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSiO4,iFe,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iKMg3AlSi3O12H2
              active(iKMg3AlSi3O12H2) = .false.  
              amount = ddust(iKMg3AlSi3O12H2)/4.Q0
              call TRANSFORM(iKMg3AlSi3O12H2,iKAlSiO4,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iMg2SiO4,amount,1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iFe2SiO4,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iFe,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Al),eps(Mg),eps(Si),eps(Fe),eps(Kalium)
            !print*,eps_save(Al),eps_save(Mg),eps_save(Si),
     >      !       eps_save(Fe),eps_save(Kalium)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Fe) = eps_save(Fe)
            eps(Kalium) = eps_save(Kalium)
            eps(O)  = eps_save(O)
          endif  
          if (active(iKMg3AlSi3O12H2).and.active(iKAlSi3O8).and.
     >        active(iMg2SiO4).and.active(iFe2SiO4).and.
     >        active(iFe3O4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iKMg3AlSi3O12H2)>Sat0(iKAlSi3O8)) then
              ioff = iKAlSi3O8
              active(iKAlSi3O8) = .false.  
              amount = ddust(iKAlSi3O8)/4.Q0
              call TRANSFORM(iKAlSi3O8,iKMg3AlSi3O12H2,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi3O8,iMg2SiO4,amount,-1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi3O8,iFe3O4,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKAlSi3O8,iFe2SiO4,amount,1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iKMg3AlSi3O12H2
              active(iKMg3AlSi3O12H2) = .false.  
              amount = ddust(iKMg3AlSi3O12H2)/4.Q0
              call TRANSFORM(iKMg3AlSi3O12H2,iKAlSi3O8,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iMg2SiO4,amount,1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKMg3AlSi3O12H2,iFe3O4,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
             call TRANSFORM(iKMg3AlSi3O12H2,iFe2SiO4,amount,-1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Al),eps(Mg),eps(Si),eps(Fe),eps(Kalium)
            !print*,eps_save(Al),eps_save(Mg),eps_save(Si),
     >      !       eps_save(Fe),eps_save(Kalium)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Fe) = eps_save(Fe)
            eps(Kalium) = eps_save(Kalium)
            eps(O)  = eps_save(O)
          endif  
          if (active(iFeS).and.active(iFe).and.
     >        active(iNi3S2).and.active(iNi)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNi)>Sat0(iNi3S2)) then
              ioff = iNi3S2
              active(iNi3S2) = .false.  
              amount = ddust(iNi3S2)/3.Q0
              call TRANSFORM(iNi3S2,iNi,amount,3.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNi3S2,iFeS,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNi3S2,iFe,amount,-2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iNi
              active(iNi) = .false.  
              amount = ddust(iNi)/3.Q0
              call TRANSFORM(iNi,iNi3S2,amount,1.Q0/3.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNi,iFeS,amount,-2.Q0/3.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNi,iFe,amount,2.Q0/3.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Fe),eps(S),eps(Ni)
            !print*,eps_save(Fe),eps_save(S),eps_save(Ni)
            eps(Fe) = eps_save(Fe)
            eps(Ni) = eps_save(Ni)
            eps(S)  = eps_save(S)
          endif   
          if (active(iAl2O3).and.active(iCaAl2Si2O8).and.
     >        active(iCaMgSi2O6).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaAl2Si2O8)>Sat0(iAl2O3)) then
              ioff = iAl2O3
              active(iAl2O3) = .false.
              amount = ddust(iAl2O3)/4.Q0
              call TRANSFORM(iAl2O3,iCaAl2Si2O8,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iCaMgSi2O6,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iMgSiO3,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iMg2SiO4,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iCaAl2Si2O8
              active(iCaAl2Si2O8) = .false.
              amount = ddust(iCaAl2Si2O8)/4.Q0
              call TRANSFORM(iCaAl2Si2O8,iAl2O3,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaAl2Si2O8,iCaMgSi2O6,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaAl2Si2O8,iMgSiO3,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaAl2Si2O8,iMg2SiO4,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Ca),eps(Mg),eps(Si),eps(Al)
            !print*,eps_save(Ca),eps_save(Mg),eps_save(Si),eps_save(Al)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Ca) = eps_save(Ca)
          endif  
          if (active(iMgAl2O4).and.active(iCaAl2Si2O8).and.
     >        active(iCaMgSi2O6).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaAl2Si2O8)>Sat0(iMgAl2O4)) then
              ioff = iMgAl2O4
              active(iMgAl2O4) = .false.
              amount = ddust(iMgAl2O4)/4.Q0
              call TRANSFORM(iMgAl2O4,iCaAl2Si2O8,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iCaMgSi2O6,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iMgSiO3,amount,-2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iMg2SiO4,amount,2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iCaAl2Si2O8
              active(iCaAl2Si2O8) = .false.
              amount = ddust(iCaAl2Si2O8)/4.Q0
              call TRANSFORM(iCaAl2Si2O8,iMgAl2O4,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaAl2Si2O8,iCaMgSi2O6,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaAl2Si2O8,iMgSiO3,amount,2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaAl2Si2O8,iMg2SiO4,amount,-2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Ca),eps(Mg),eps(Si),eps(Al)
            !print*,eps_save(Ca),eps_save(Mg),eps_save(Si),eps_save(Al)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Ca) = eps_save(Ca)
          endif  
          if (active(iZrO2).and.active(iZrSiO4).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iZrO2)>Sat0(iZrSiO4)) then
              ioff = iZrSiO4
              active(iZrSiO4) = .false.
              amount = ddust(iZrSiO4)/3.Q0
              call TRANSFORM(iZrSiO4,iZrO2,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iZrSiO4,iMgSiO3,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iZrSiO4,iMg2SiO4,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iZrO2
              active(iZrO2) = .false.
              amount = ddust(iZrO2)/3.Q0
              call TRANSFORM(iZrO2,iZrSiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iZrO2,iMgSiO3,amount,-2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iZrO2,iMg2SiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Zr),eps(Mg),eps(Si)
            !print*,eps_save(Zr),eps_save(Mg),eps_save(Si)
            eps(Zr) = eps_save(Zr)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif  
          if (active(iZrO2).and.active(iZrSiO4).and.
     >        active(iSiO)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iZrO2)>Sat0(iZrSiO4)) then
              ioff = iZrSiO4
              active(iZrSiO4) = .false.
              amount = ddust(iZrSiO4)/2.Q0
              call TRANSFORM(iZrSiO4,iZrO2,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iZrSiO4,iSiO,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iZrO2
              active(iZrO2) = .false.
              amount = ddust(iZrO2)/2.Q0
              call TRANSFORM(iZrO2,iZrSiO4,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iZrO2,iSiO,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Zr),eps(O),eps(Si)
            !print*,eps_save(Zr),eps_save(O),eps_save(Si)
            eps(Zr) = eps_save(Zr)
            eps(Si) = eps_save(Si)
          endif  
          if (active(iKCl).and.active(iKAlSi3O8).and.
     >        active(iNaCl).and.active(iMgAl2O4).and.
     >        active(iMgSiO3).and.active(iMg2SiO4).and.
     >        active(iNa2SiO3)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iKCl)>Sat0(iKAlSi3O8)) then
              if (ddust(iKAlSi3O8)<ddust(iNaCl)) then 
                ioff = iKAlSi3O8 
                active(iKAlSi3O8) = .false.
                amount = ddust(iKAlSi3O8)/6.Q0
                call TRANSFORM(iKAlSi3O8,iKCl,amount,1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iNaCl,amount,-1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iNa2SiO3,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iMgAl2O4,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iMgSiO3,amount,5.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iMg2SiO4,amount,-3.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
              else
                ioff = iNaCl 
                active(iNaCl) = .false.
                amount = ddust(iNaCl)/6.Q0
                call TRANSFORM(iNaCl,iKCl,amount,1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iKAlSi3O8,amount,-1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iNa2SiO3,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iMgAl2O4,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iMgSiO3,amount,5.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iMg2SiO4,amount,-3.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
              endif   
            else  
              ioff = iKCl 
              active(iKCl) = .false.
              amount = ddust(iKCl)/6.Q0
              call TRANSFORM(iKCl,iKAlSi3O8,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iNaCl,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iNa2SiO3,amount,-0.5Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iMgAl2O4,amount,-0.5Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iMgSiO3,amount,-5.5Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iMg2SiO4,amount,3.Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            eps(Kalium) = eps_save(Kalium)
            eps(Cl) = eps_save(Cl)
            eps(Al) = eps_save(Al)
            eps(Na) = eps_save(Na)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif               
          if (active(iKCl).and.active(iKAlSi3O8).and.
     >        active(iNaCl).and.active(iAl2O3).and.
     >        active(iMgSiO3).and.active(iMg2SiO4).and.
     >        active(iNa2SiO3)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iKCl)>Sat0(iKAlSi3O8)) then
              if (ddust(iKAlSi3O8)<ddust(iNaCl)) then 
                ioff = iKAlSi3O8 
                active(iKAlSi3O8) = .false.
                amount = ddust(iKAlSi3O8)/6.Q0
                call TRANSFORM(iKAlSi3O8,iKCl,amount,1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iNaCl,amount,-1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iNa2SiO3,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iAl2O3,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iMgSiO3,amount,5.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iKAlSi3O8,iMg2SiO4,amount,-2.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
              else
                ioff = iNaCl 
                active(iNaCl) = .false.
                amount = ddust(iNaCl)/6.Q0
                call TRANSFORM(iNaCl,iKCl,amount,1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iKAlSi3O8,amount,-1.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iNa2SiO3,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iAl2O3,amount,0.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iMgSiO3,amount,5.Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaCl,iMg2SiO4,amount,-2.5Q0*6.Q0,
     >                         ddust,eps,dscale,active,ok)
              endif   
            else  
              ioff = iKCl 
              active(iKCl) = .false.
              amount = ddust(iKCl)/6.Q0
              call TRANSFORM(iKCl,iKAlSi3O8,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iNaCl,amount,1.Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iNa2SiO3,amount,-0.5Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iAl2O3,amount,-0.5Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iMgSiO3,amount,-5.0Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iKCl,iMg2SiO4,amount,2.5Q0*6.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            eps(Kalium) = eps_save(Kalium)
            eps(Cl) = eps_save(Cl)
            eps(Al) = eps_save(Al)
            eps(Na) = eps_save(Na)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif               
          if (active(iAl2O3).and.active(iMgAl2O4).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iAl2O3)>Sat0(iMgAl2O4)) then
              ioff = iMgAl2O4
              active(iMgAl2O4) = .false.
              amount = ddust(iMgAl2O4)/3.Q0
              call TRANSFORM(iMgAl2O4,iAl2O3,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iMgSiO3,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iMg2SiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              if (ddust(iAl2O3)<ddust(iMg2SiO4)) then 
                ioff = iAl2O3
                active(iAl2O3) = .false.
                amount = ddust(iAl2O3)/3.Q0
                call TRANSFORM(iAl2O3,iMgAl2O4,amount,1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iAl2O3,iMgSiO3,amount,1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iAl2O3,iMg2SiO4,amount,-1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
              else  
                ioff = iMg2SiO4
                active(iMg2SiO4) = .false.
                amount = ddust(iMg2SiO4)/3.Q0
                call TRANSFORM(iMg2SiO4,iMgAl2O4,amount,1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iMg2SiO4,iAl2O3,amount,-1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iMg2SiO4,iMgSiO3,amount,1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
              endif  
            endif  
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif  
          if (active(iAl2O3).and.active(iMgAl2O4).and.
     >        active(iCaMgSi2O6).and.active(iCaSiO3).and.
     >        active(iSiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iAl2O3)>Sat0(iMgAl2O4)) then
              ioff = iMgAl2O4
              active(iMgAl2O4) = .false.
              amount = ddust(iMgAl2O4)/4.Q0
              call TRANSFORM(iMgAl2O4,iAl2O3,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iCaMgSi2O6,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iCaSiO3,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iSiO2,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iAl2O3
              active(iAl2O3) = .false.
              amount = ddust(iAl2O3)/4.Q0
              call TRANSFORM(iAl2O3,iMgAl2O4,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iCaMgSi2O6,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iCaSiO3,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iSiO2,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Al),eps(Mg),eps(Si),eps(Ca),eps(O)
            !print*,eps_save(Al),eps_save(Mg),eps_save(Si),
     >      !       eps_save(Ca),eps_save(O)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(Ca) = eps_save(Ca)
            eps(O)  = eps_save(O)
          endif  
          if (active(iFeO).and.active(iFe2SiO4).and.active(iSiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFe2SiO4)>Sat0(iFeO)) then
              ioff = iFeO
              active(iFeO) = .false.  
              amount = ddust(iFeO)/2.Q0
              call TRANSFORM(iFeO,iFe2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFeO,iSiO2,amount,-0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFe2SiO4
              active(iFe2SiO4) = .false.  
              amount = ddust(iFe2SiO4)/2.Q0
              call TRANSFORM(iFe2SiO4,iFeO,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iSiO2,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Fe) = eps_save(Fe)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iFeO_l).and.active(iFe2SiO4).and.
     >        active(iSiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFe2SiO4)>Sat0(iFeO_l)) then
              ioff = iFeO_l
              active(iFeO_l) = .false.  
              amount = ddust(iFeO_l)/2.Q0
              call TRANSFORM(iFeO_l,iFe2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFeO_l,iSiO2,amount,-0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFe2SiO4
              active(iFe2SiO4) = .false.  
              amount = ddust(iFe2SiO4)/2.Q0
              call TRANSFORM(iFe2SiO4,iFeO_l,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iSiO2,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            eps(Fe) = eps_save(Fe)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iFeO_l).and.active(iFe2SiO4).and.
     >        active(iSiO2_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFe2SiO4)>Sat0(iFeO_l)) then
              ioff = iFeO_l
              active(iFeO_l) = .false.  
              amount = ddust(iFeO_l)/2.Q0
              call TRANSFORM(iFeO_l,iFe2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFeO_l,iSiO2_l,amount,-0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFe2SiO4
              active(iFe2SiO4) = .false.  
              amount = ddust(iFe2SiO4)/2.Q0
              call TRANSFORM(iFe2SiO4,iFeO_l,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iSiO2_l,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Fe) = eps_save(Fe)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iSiO).and.active(iMgSiO3).and.
     >        active(iMg2SiO4).and.active(iC).and..false.) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgSiO3)>Sat0(iSiO)) then
              ioff = iSiO
              active(iSiO) = .false.  
              amount = ddust(iSiO)/2.Q0
              call TRANSFORM(iSiO,iMgSiO3,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iSiO,iMg2SiO4,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              active(iMgSiO3) = .false.  
              amount = ddust(iMgSiO3)/2.Q0
              call TRANSFORM(iMgSiO3,iSiO,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgSiO3,iMg2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif
          if (active(iSiO2).and.active(iMgSiO3).and.
     >        active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgSiO3)>Sat0(iSiO2)) then
              ioff = iSiO2
              active(iSiO2) = .false.  
              amount = ddust(iSiO2)/2.Q0
              call TRANSFORM(iSiO2,iMgSiO3,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iSiO2,iMg2SiO4,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgSiO3
              active(iMgSiO3) = .false.  
              amount = ddust(iMgSiO3)/2.Q0
              call TRANSFORM(iMgSiO3,iSiO2,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgSiO3,iMg2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif
          if (active(iSiO2_l).and.active(iMgSiO3_l).and.
     >        active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgSiO3_l)>Sat0(iSiO2_l)) then
              ioff = iSiO2_l
              active(iSiO2_l) = .false.  
              amount = ddust(iSiO2_l)/2.Q0
              call TRANSFORM(iSiO2_l,iMgSiO3_l,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iSiO2_l,iMg2SiO4,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgSiO3_l
              active(iMgSiO3_l) = .false.  
              amount = ddust(iMgSiO3_l)/2.Q0
              call TRANSFORM(iMgSiO3_l,iSiO2_l,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgSiO3_l,iMg2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Mg),eps(Si)
            !print*,eps_save(Mg),eps_save(Si)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif
          if (active(iSiO2_l).and.active(iMgSiO3).and.
     >        active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgSiO3)>Sat0(iSiO2_l)) then
              ioff = iSiO2_l
              active(iSiO2_l) = .false.  
              amount = ddust(iSiO2_l)/2.Q0
              call TRANSFORM(iSiO2_l,iMgSiO3,amount,2.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iSiO2_l,iMg2SiO4,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgSiO3
              active(iMgSiO3) = .false.  
              amount = ddust(iMgSiO3)/2.Q0
              call TRANSFORM(iMgSiO3,iSiO2_l,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgSiO3,iMg2SiO4,amount,0.5Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Mg),eps(Si)
            !print*,eps_save(Mg),eps_save(Si)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif
          if (active(iTi4O7).and.active(iCaTiO3).and.
     >        active(iCaMgSi2O6).and.(.not.active(iMg2SiO4))) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaTiO3)>Sat0(iTi4O7)) then
              ioff = iTi4O7
              active(iTi4O7) = .false.  
              amount = ddust(iTi4O7)/2.Q0
              call TRANSFORM(iTi4O7,iCaTiO3,amount,4.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTi4O7,iCaMgSi2O6,amount,-4.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iCaTiO3
              active(iCaTiO3) = .false.  
              amount = ddust(iCaTiO3)/2.Q0
              call TRANSFORM(iCaTiO3,iTi4O7,amount,0.25Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaTiO3,iCaMgSi2O6,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ca) = eps_save(Ca)
            eps(Ti) = eps_save(Ti)
          endif  
          if (active(iCaMgSi2O6).and.active(iCaSiO3).and.
     >        active(iMgSiO3_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaMgSi2O6)>Sat0(iMgSiO3_l)) then
              if (ddust(iMgSiO3_l)<ddust(iCaSiO3)) then 
                ioff = iMgSiO3_l
                active(iMgSiO3_l) = .false.  
                amount = ddust(iMgSiO3_l)/2.Q0
                call TRANSFORM(iMgSiO3_l,iCaMgSi2O6,amount,1.Q0*2.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iMgSiO3_l,iCaSiO3,amount,-1.Q0*2.Q0,
     >                         ddust,eps,dscale,active,ok)
              else  
                ioff = iCaSiO3
                active(iCaSiO3) = .false.  
                amount = ddust(iCaSiO3)/2.Q0
                call TRANSFORM(iCaSiO3,iCaMgSi2O6,amount,1.Q0*2.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iCaSiO3,iMgSiO3_l,amount,-1.Q0*2.Q0,
     >                         ddust,eps,dscale,active,ok)
              endif  
            else  
              ioff = iCaMgSi2O6
              active(iCaMgSi2O6) = .false.  
              amount = ddust(iCaMgSi2O6)/2.Q0
              call TRANSFORM(iCaMgSi2O6,iMgSiO3_l,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaMgSi2O6,iCaSiO3,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ca) = eps_save(Ca)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif  
          if (active(iTi4O7).and.active(iTiC).and.active(iC)
     >        .and..false.) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTi4O7)>Sat0(iTiC)) then
              ioff = iTiC
              active(iTiC) = .false.  
              amount = ddust(iTiC)/2.Q0
              call TRANSFORM(iTiC,iTi4O7,amount,0.25Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiC,iC,amount,1Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iTi4O7
              active(iTi4O7) = .false.  
              amount = ddust(iTi4O7)/2.Q0
              call TRANSFORM(iTi4O7,iTiC,amount,4.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTi4O7,iC,amount,-4.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ti) = eps_save(Ti)
            eps(C)  = eps_save(C)
          endif   
          if (active(iTi4O7).and.active(iTiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTi4O7)>Sat0(iTiO2)) then
              ioff = iTiO2
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)
              call TRANSFORM(iTiO2,iTi4O7,amount,0.25Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iTi4O7
              active(iTi4O7) = .false.  
              amount = ddust(iTi4O7)
              call TRANSFORM(iTi4O7,iTiO2,amount,4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ti) = eps_save(Ti)
          endif   
          if (active(iTi4O7).and.active(iTiO)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTi4O7)>Sat0(iTiO)) then
              ioff = iTiO
              active(iTiO) = .false.  
              amount = ddust(iTiO)
              call TRANSFORM(iTiO,iTi4O7,amount,0.25Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iTi4O7
              active(iTi4O7) = .false.  
              amount = ddust(iTi4O7)
              call TRANSFORM(iTi4O7,iTiO,amount,4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ti) = eps_save(Ti)
          endif   
          if (active(iTi4O7_l).and.active(iTiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTi4O7_l)>Sat0(iTiO2)) then
              ioff = iTiO2
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)
              call TRANSFORM(iTiO2,iTi4O7_l,amount,0.25Q0,
     >                       ddust,eps,dscale,active,ok)
            else
              ioff = iTi4O7_l
              active(iTi4O7_l) = .false.  
              amount = ddust(iTi4O7_l)
              call TRANSFORM(iTi4O7_l,iTiO2,amount,4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ti) = eps_save(Ti)
          endif   
          if (active(iMgTi2O5).and.active(iTiO2).and.
     >        active(iMgSiO3).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgTi2O5)>Sat0(iTiO2)) then
              ioff = iTiO2
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)/3.Q0
              call TRANSFORM(iTiO2,iMgTi2O5,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iMgSiO3,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iMg2SiO4,amount,-0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgTi2O5
              active(iMgTi2O5) = .false.  
              amount = ddust(iMgTi2O5)/3.Q0
              call TRANSFORM(iMgTi2O5,iTiO2,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iMgSiO3,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iMg2SiO4,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            eps(Ti) = eps_save(Ti)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iMgTi2O5).and.active(iTiO2).and.
     >        active(iSiO2).and.active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgTi2O5)>Sat0(iTiO2)) then
              ioff = iTiO2
              active(iTiO2) = .false.  
              amount = ddust(iTiO2)/3.Q0
              call TRANSFORM(iTiO2,iMgTi2O5,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iSiO2,amount,0.25Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iMg2SiO4,amount,-0.25Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgTi2O5
              active(iMgTi2O5) = .false.  
              amount = ddust(iMgTi2O5)/3.Q0
              call TRANSFORM(iMgTi2O5,iTiO2,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iSiO2,amount,-0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iMg2SiO4,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ti) = eps_save(Ti)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iNa2SiO3).and.active(iNaAlSi3O8).and.
     >        active(iMgSiO3).and.active(iMg2SiO4).and.
     >        active(iMgAl2O4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNa2SiO3)>Sat0(iNaAlSi3O8)) then
              if (3.Q0*ddust(iNaAlSi3O8)<ddust(iMg2SiO4)) then 
                ioff = iNaAlSi3O8
                active(iNaAlSi3O8) = .false.  
                amount = ddust(iNaAlSi3O8)/4.Q0
                call TRANSFORM(iNaAlSi3O8,iNa2SiO3,amount,0.5Q0*4.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaAlSi3O8,iMgAl2O4,amount,0.5Q0*4.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaAlSi3O8,iMgSiO3,amount,5.5Q0*4.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iNaAlSi3O8,iMg2SiO4,amount,-3.Q0*4.Q0,
     >                         ddust,eps,dscale,active,ok)
              else
                ioff = iMg2SiO4
                active(iMg2SiO4) = .false.  
                amount = ddust(iMg2SiO4)/4.Q0
                call TRANSFORM(iMg2SiO4,iNaAlSi3O8,amount,-1.Q0/3.Q0
     >                         *4.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iMg2SiO4,iNa2SiO3,amount,1.Q0/6.Q0
     >                         *4.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iMg2SiO4,iMgAl2O4,amount,1.Q0/6.Q0
     >                         *4.Q0,ddust,eps,dscale,active,ok)
                call TRANSFORM(iMg2SiO4,iMgSiO3,amount,11.Q0/6.Q0
     >                         *4.Q0,ddust,eps,dscale,active,ok)
              endif   
            else  
              ioff = iNa2SiO3
              active(iNa2SiO3) = .false.  
              amount = ddust(iNa2SiO3)/4.Q0
              call TRANSFORM(iNa2SiO3,iNaAlSi3O8,amount,2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3,iMgAl2O4,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3,iMgSiO3,amount,-11.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3,iMg2SiO4,amount,6.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Na) = eps_save(Na)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif   
          if (active(iNa2SiO3_l).and.active(iNaAlSi3O8).and.
     >        active(iAl2O3_l).and.active(iSiO2)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNa2SiO3_l)>Sat0(iNaAlSi3O8)) then
              ioff = iNaAlSi3O8
              active(iNaAlSi3O8) = .false.  
              amount = ddust(iNaAlSi3O8)/3.Q0
              call TRANSFORM(iNaAlSi3O8,iNa2SiO3_l,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSi3O8,iAl2O3_l,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSi3O8,iSiO2,amount,2.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iNa2SiO3_l
              active(iNa2SiO3_l) = .false.  
              amount = ddust(iNa2SiO3_l)/3.Q0
              call TRANSFORM(iNa2SiO3_l,iNaAlSi3O8,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3_l,iAl2O3_l,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3_l,iSiO2,amount,-5.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Na),eps(Al),eps(Si)
            !print*,eps_save(Na),eps_save(Al),eps_save(Si)
            eps(Al) = eps_save(Al)
            eps(Na) = eps_save(Na)
            eps(Si) = eps_save(Si)
          endif  
          if (active(iNa2SiO3_l).and.active(iNaAlSi3O8).and.
     >        active(iAl2O3_l).and.active(iSiO2_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNa2SiO3_l)>Sat0(iNaAlSi3O8)) then
              ioff = iNaAlSi3O8
              active(iNaAlSi3O8) = .false.  
              amount = ddust(iNaAlSi3O8)/3.Q0
              call TRANSFORM(iNaAlSi3O8,iNa2SiO3_l,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSi3O8,iAl2O3_l,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNaAlSi3O8,iSiO2_l,amount,2.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iNa2SiO3_l
              active(iNa2SiO3_l) = .false.  
              amount = ddust(iNa2SiO3_l)/3.Q0
              call TRANSFORM(iNa2SiO3_l,iNaAlSi3O8,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3_l,iAl2O3_l,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iNa2SiO3_l,iSiO2_l,amount,-5.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif
            !print*,eps(Na),eps(Al),eps(Si),eps(O)
            !print*,eps_save(Na),eps_save(Al),eps_save(Si),eps_save(O)
            eps(Al) = eps_save(Al)
            eps(Na) = eps_save(Na)
            eps(Si) = eps_save(Si)
          endif  
          if (active(iH2O).and.active(iH2O_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iH2O)>Sat0(iH2O_l)) then
              ioff = iH2O_l
              active(iH2O_l) = .false.  
              amount = ddust(iH2O_l)
              call TRANSFORM(iH2O_l,iH2O,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iH2O
              active(iH2O) = .false.  
              amount = ddust(iH2O)
              call TRANSFORM(iH2O,iH2O_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(O) = eps_save(O)
            eps(H) = eps_save(H)
          endif   
          if (active(iTiO).and.active(iTiO_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTiO)>Sat0(iTiO_l)) then
              ioff = iTiO_l
              active(iTiO_l) = .false.  
              amount = ddust(iTiO_l)
              call TRANSFORM(iTiO_l,iTiO,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iTiO
              active(iTiO) = .false.  
              amount = ddust(iTiO)
              call TRANSFORM(iTiO,iTiO_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ti) = eps_save(Ti)
            eps(O)  = eps_save(O)
          endif   
          if (active(iFe).and.active(iFe_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFe)>Sat0(iFe_l)) then
              ioff = iFe_l
              active(iFe_l) = .false.  
              amount = ddust(iFe_l)
              call TRANSFORM(iFe_l,iFe,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFe
              active(iFe) = .false.  
              amount = ddust(iFe)
              call TRANSFORM(iFe,iFe_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Fe) = eps_save(Fe)
          endif   
          if (active(iNi).and.active(iNi_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNi)>Sat0(iNi_l)) then
              ioff = iNi_l
              active(iNi_l) = .false.  
              amount = ddust(iNi_l)
              call TRANSFORM(iNi_l,iNi,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iNi
              active(iNi) = .false.  
              amount = ddust(iNi)
              call TRANSFORM(iNi,iNi_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Ni) = eps_save(Ni)
          endif   
          if (active(iFeO).and.active(iFeO_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFeO)>Sat0(iFeO_l)) then
              ioff = iFeO_l
              active(iFeO_l) = .false.  
              amount = ddust(iFeO_l)
              call TRANSFORM(iFeO_l,iFeO,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFeO
              active(iFeO) = .false.  
              amount = ddust(iFeO)
              call TRANSFORM(iFeO,iFeO_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Fe) = eps_save(Fe)
            eps(O)  = eps_save(O)
          endif   
          if (active(iAl2O3).and.active(iAl2O3_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iAl2O3)>Sat0(iAl2O3_l)) then
              ioff = iAl2O3_l
              active(iAl2O3_l) = .false.  
              amount = ddust(iAl2O3_l)
              call TRANSFORM(iAl2O3_l,iAl2O3,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iAl2O3
              active(iAl2O3) = .false.  
              amount = ddust(iAl2O3)
              call TRANSFORM(iAl2O3,iAl2O3_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Al) = eps_save(Al)
            eps(O) = eps_save(O)
          endif   
          if (active(iMgAl2O4).and.active(iMgAl2O4_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgAl2O4)>Sat0(iMgAl2O4_l)) then
              ioff = iMgAl2O4_l
              active(iMgAl2O4_l) = .false.  
              amount = ddust(iMgAl2O4_l)
              call TRANSFORM(iMgAl2O4_l,iMgAl2O4,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgAl2O4
              active(iMgAl2O4) = .false.  
              amount = ddust(iMgAl2O4)
              call TRANSFORM(iMgAl2O4,iMgAl2O4_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Mg) = eps_save(Mg)
            eps(Al) = eps_save(Al)
            eps(O)  = eps_save(O)
          endif   
          if (active(iSiO2).and.active(iSiO2_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iSiO2)>Sat0(iSiO2_l)) then
              ioff = iSiO2_l
              active(iSiO2_l) = .false.  
              amount = ddust(iSiO2_l)
              call TRANSFORM(iSiO2_l,iSiO2,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iSiO2
              active(iSiO2) = .false.  
              amount = ddust(iSiO2)
              call TRANSFORM(iSiO2,iSiO2_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif   
          if (active(iFeS).and.active(iFeS_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFeS)>Sat0(iFeS_l)) then
              ioff = iFeS_l
              active(iFeS_l) = .false.  
              amount = ddust(iFeS_l)
              call TRANSFORM(iFeS_l,iFeS,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFeS
              active(iFeS) = .false.  
              amount = ddust(iFeS)
              call TRANSFORM(iFeS,iFeS_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Fe) = eps_save(Fe)
            eps(S)  = eps_save(S)
          endif   
          if (active(iLiCl).and.active(iLiCl_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iLiCl)>Sat0(iLiCl_l)) then
              ioff = iLiCl_l
              active(iLiCl_l) = .false.  
              amount = ddust(iLiCl_l)
              call TRANSFORM(iLiCl_l,iLiCl,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iLiCl
              active(iLiCl) = .false.  
              amount = ddust(iLiCl)
              call TRANSFORM(iLiCl,iLiCl_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Li) = eps_save(Li)
            eps(Cl) = eps_save(Cl)
          endif   
          if (active(iNaCl).and.active(iNaCl_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iNaCl)>Sat0(iNaCl_l)) then
              ioff = iNaCl_l
              active(iNaCl_l) = .false.  
              amount = ddust(iNaCl_l)
              call TRANSFORM(iNaCl_l,iNaCl,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iNaCl
              active(iNaCl) = .false.  
              amount = ddust(iNaCl)
              call TRANSFORM(iNaCl,iNaCl_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Li) = eps_save(Li)
            eps(Cl) = eps_save(Cl)
          endif   
          if (active(iW).and.active(iW_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iW)>Sat0(iW_l)) then
              ioff = iW_l
              active(iW_l) = .false.  
              amount = ddust(iW_l)
              call TRANSFORM(iW_l,iW,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iW
              active(iW) = .false.  
              amount = ddust(iW)
              call TRANSFORM(iW,iW_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(W) = eps_save(W)
          endif   
          if (active(iMg2SiO4).and.active(iMg2SiO4_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMg2SiO4)>Sat0(iMg2SiO4_l)) then
              ioff = iMg2SiO4_l
              active(iMg2SiO4_l) = .false.  
              amount = ddust(iMg2SiO4_l)
              call TRANSFORM(iMg2SiO4_l,iMg2SiO4,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMg2SiO4
              active(iMg2SiO4) = .false.  
              amount = ddust(iMg2SiO4)
              call TRANSFORM(iMg2SiO4,iMg2SiO4_l,amount,1.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            eps(Li) = eps_save(Li)
            eps(Cl) = eps_save(Cl)
          endif   
          if (active(iMg2SiO4).and.active(iMgSiO3).and.
     >        active(iFe).and.active(iFe2SiO4).and.
     >        .not.active(iFeS)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFe)>Sat0(iFe2SiO4)) then
              ioff = iFe2SiO4
              active(iFe2SiO4) = .false.
              amount = ddust(iFe2SiO4)/3.Q0
              call TRANSFORM(iFe2SiO4,iFe,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iMgSiO3,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iMg2SiO4,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iFe
              active(iFe) = .false.
              amount = ddust(iFe)/3.Q0
              call TRANSFORM(iFe,iFe2SiO4,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe,iMgSiO3,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe,iMg2SiO4,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Fe),eps(Mg),eps(Si)
            !print*,eps_save(Fe),eps_save(Mg),eps_save(Si)
            eps(Fe) = eps_save(Fe)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif 
          if (active(iMg2SiO4).and.active(iMgSiO3).and.
     >        active(iFeO).and.active(iFe2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iFeO)>Sat0(iFe2SiO4)) then
              ioff = iFe2SiO4
              active(iFe2SiO4) = .false.
              amount = ddust(iFe2SiO4)/3.Q0
              call TRANSFORM(iFe2SiO4,iFeO,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iMgSiO3,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iFe2SiO4,iMg2SiO4,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              if (ddust(iFeO)<ddust(iMgSiO3)) then 
                ioff = iFeO
                active(iFeO) = .false.
                amount = ddust(iFeO)/3.Q0
                call TRANSFORM(iFeO,iFe2SiO4,amount,0.5Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iFeO,iMgSiO3,amount,-1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iFeO,iMg2SiO4,amount,0.5Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
              else
                ioff = iMgSiO3
                active(iMgSiO3) = .false.
                amount = ddust(iMgSiO3)/3.Q0
                call TRANSFORM(iMgSiO3,iFe2SiO4,amount,0.5Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iMgSiO3,iFeO,amount,-1.Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
                call TRANSFORM(iMgSiO3,iMg2SiO4,amount,0.5Q0*3.Q0,
     >                         ddust,eps,dscale,active,ok)
              endif   
            endif  
            !print*,eps(Fe),eps(Mg),eps(Si),eps(O)
            !print*,eps_save(Fe),eps_save(Mg),eps_save(Si),eps_save(O)
            eps(Fe) = eps_save(Fe)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
            eps(O)  = eps_save(O)
          endif 
          if (active(iMgTi2O5).and.active(iTiO2).and.
     >        active(iMgSiO3_l).and.active(iSiO2_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTiO2)>Sat0(iMgTi2O5)) then
              ioff = iMgTi2O5
              active(iMgTi2O5) = .false.
              amount = ddust(iMgTi2O5)/3.Q0
              call TRANSFORM(iMgTi2O5,iTiO2,amount,2.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iMgSiO3_l,amount,1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iSiO2_l,amount,-1.Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iTiO2
              active(iTiO2) = .false.
              amount = ddust(iTiO2)/3.Q0
              call TRANSFORM(iTiO2,iMgTi2O5,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iMgSiO3_l,amount,-0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iSiO2_l,amount,0.5Q0*3.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Ti),eps(Mg),eps(Si)
            !print*,eps_save(Ti),eps_save(Mg),eps_save(Si)
            eps(Ti) = eps_save(Ti)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif 
          if (active(iMgTi2O5).and.active(iTiO2).and.
     >        active(iCaMgSi2O6).and.active(iCaSiO3).and.
     >        active(iSiO2_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iTiO2)>Sat0(iMgTi2O5)) then
              ioff = iMgTi2O5
              active(iMgTi2O5) = .false.
              amount = ddust(iMgTi2O5)/4.Q0
              call TRANSFORM(iMgTi2O5,iTiO2,amount,2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iCaMgSi2O6,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iCaSiO3,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iSiO2_l,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iTiO2
              active(iTiO2) = .false.
              amount = ddust(iTiO2)/4.Q0
              call TRANSFORM(iTiO2,iMgTi2O5,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iCaMgSi2O6,amount,-0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iCaSiO3,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iTiO2,iSiO2_l,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Ti),eps(Mg),eps(Si),eps(Ca)
            !print*,eps_save(Ti),eps_save(Mg),eps_save(Si),eps_save(Ca)
            eps(Ca) = eps_save(Ca)
            eps(Ti) = eps_save(Ti)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif 
          if (active(iMgTi2O5).and.active(iCaTiO3).and.
     >        active(iCaSiO3).and.active(iMgSiO3_l).and.
     >        active(iMg2SiO4)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iCaTiO3)>Sat0(iMgTi2O5)) then
              ioff = iMgTi2O5
              active(iMgTi2O5) = .false.
              amount = ddust(iMgTi2O5)/4.Q0
              call TRANSFORM(iMgTi2O5,iCaTiO3,amount,2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iCaSiO3,amount,-2.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iMgSiO3_l,amount,3.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgTi2O5,iMg2SiO4,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iCaTiO3
              active(iCaTiO3) = .false.
              amount = ddust(iCaTiO3)/4.Q0
              call TRANSFORM(iCaTiO3,iMgTi2O5,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaTiO3,iCaSiO3,amount,1.0Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaTiO3,iMgSiO3_l,amount,-1.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iCaTiO3,iMg2SiO4,amount,0.5Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Ti),eps(Mg),eps(Si),eps(Ca)
            !print*,eps_save(Ti),eps_save(Mg),eps_save(Si),eps_save(Ca)
            eps(Ca) = eps_save(Ca)
            eps(Ti) = eps_save(Ti)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif 
          if (active(iAl2O3).and.active(iMgAl2O4).and.
     >        active(iCaMgSi2O6).and.active(iCaSiO3).and.
     >        active(iSiO2_l)) then
            changed = .true.
            !--- decide ---
            if (Sat0(iMgAl2O4)>Sat0(iAl2O3)) then
              ioff = iAl2O3
              active(iAl2O3) = .false.
              amount = ddust(iAl2O3)/4.Q0
              call TRANSFORM(iAl2O3,iMgAl2O4,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iCaMgSi2O6,amount,-1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iCaSiO3,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iAl2O3,iSiO2_l,amount,1.Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            else  
              ioff = iMgAl2O4
              active(iMgAl2O4) = .false.
              amount = ddust(iMgAl2O4)/4.Q0
              call TRANSFORM(iMgAl2O4,iAl2O3,amount,1.0Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iCaMgSi2O6,amount,1.0Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iCaSiO3,amount,-1.0Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iMgAl2O4,iSiO2_l,amount,-1.0Q0*4.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif  
            !print*,eps(Al),eps(Mg),eps(Si),eps(Ca)
            !print*,eps_save(Al),eps_save(Mg),eps_save(Si),eps_save(Ca)
            eps(Ca) = eps_save(Ca)
            eps(Al) = eps_save(Al)
            eps(Mg) = eps_save(Mg)
            eps(Si) = eps_save(Si)
          endif 
          if (active(iSiO).and.active(iSiC).and.active(iC)
     >        .and..false.) then
            changed = .true.
            !--- decide ---
            if (Sat0(iSiC)>Sat0(iSiO)) then
              ioff = iSiO
              active(iSiO) = .false.
              amount = ddust(iSiO)/2.Q0
              call TRANSFORM(iSiO,iSiC,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iSiO,iC,amount,-1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            else
              ioff = iSiC
              active(iSiC) = .false.
              amount = ddust(iSiC)/2.Q0
              call TRANSFORM(iSiC,iSiO,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
              call TRANSFORM(iSiC,iC,amount,1.Q0*2.Q0,
     >                       ddust,eps,dscale,active,ok)
            endif   
            !print*,eps(C),eps(O),eps(Si)
            !print*,eps_save(C),eps_save(O),eps_save(Si)
            eps(Si) = eps_save(Si)
            eps(C)  = eps_save(C)  
          endif  
          if (.not.ok) then
            print*,"TRANSFORM resulted in negative dust abundance"
            ifail = ifail+1
            if (ifail>=3) stop "*** too many fails."
            print*,"try to carry on without ..."
            call SLEEP(1)
            eps = eps_save 
            ddust = dust_save 
            !active(ioff) = .true.
            active = active_save 
            Nact = Nact_save
            changed = .false.
            lastit = it
          endif  
        endif  
        if (ioff>0) itransform=itransform+1

        laston = 0
        do i=1,NDUST
          if (active(i).and.(.not.act_old(i))) then
            laston = i
          endif
        enddo
  
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
            if (active(i).or.Sat0(i)>0.1) then
              write(*,'(3x,A18,2(1pE11.3)1pE19.10,A2)') 
     >          dust_nam(i),ddust(i),ddust(i)/dscale(i),Sat0(i),rem
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
            do ii=1,Nind                  
              if (eps(el)<Iabund(ii)) exit  ! sort by element abundance 
            enddo
            Iindex(ii+1:Nind+1) = Iindex(ii:Nind)
            Iabund(ii+1:Nind+1) = Iabund(ii:Nind)
            Iindex(ii) = el
            Iabund(ii) = eps(el) 
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
        if (verbose>1) print'(99(A3))',
     >                    (trim(elnam(Iindex(j))),j=1,Nind-1)
        if (verbose>1) print'(99(I3))',e_num(Iindex(1:Nind-1)) 
        if (Nind-1<Nact) stop "*** Nind<Nact in equil_cond."
        Nall = Nind-1
        Nind = Nact                         ! truncate at number of condensates

        !-------------------------------------
        ! ***  some explict special cases  ***
        !-------------------------------------
        if (active(iCaS).and.active(iMnS).and.(Nall>Nact)
     >      .and.e_act(Ca).and.e_act(S).and.e_act(Mn)) then
          exch = (e_num(Ca)+e_num(S)+e_num(Mn)==4)
     >      .or.((e_num(Ca)+e_num(S)+e_num(Mn)==5).and.(e_num(Al)==1))
     >      .or.((e_num(Ca)+e_num(S)+e_num(Mn)==5).and.(e_num(Al)==2)
     >           .and.(e_num(Na)==1).and.active(iNaAlSi3O8))
          if (exch) then       
            print*,"... exchanging S for "//elnam(Iindex(Nact+1))
            do i=1,Nind
              if (Iindex(i)==S) exit
            enddo  
            swap = Iindex(Nact+1) 
            Iindex(Nact+1) = S
            Iindex(i) = swap
            e_act(S) = .false.
            e_act(swap) = .true.
          endif  
        endif
        if (active(iCaS).and.e_act(Ca).and.e_act(S).and.Nall>Nact) then
          exch = ((e_num(Ca)==1).and.(e_num(S)==1)) 
          exch = ((e_num(Ca)==2).and.(e_num(S)==1).and.
     >            (e_num(Al)==2).and.(e_num(Si)==1).and. 
     >            active(iCa2Al2SiO7).and.active(iAl2O3).and.
     >            e_act(Si))
          if (exch) then
            do i=1,Nind
              if (Iindex(i)==Ca) exit
            enddo  
            do j=1,Nind
              if (Iindex(j)==S) exit
            enddo  
            if (i>j) then
              e_act(Iindex(i)) = .false.
              j=i
            else
              e_act(Iindex(j)) = .false.
              i=j
            endif
            print*,"... exchanging "//elnam(Iindex(i))//
     >                       " for "//elnam(Iindex(Nact+1))
            swap = Iindex(i)
            Iindex(i) = Iindex(Nact+1) 
            Iindex(Nact+1) = swap
            e_act(Iindex(i)) = .true.
          endif  
        endif
        if (active(iCa2Al2SiO7).and.
     >      e_act(Ca).and.e_act(Al).and.(.not.e_act(Si)).and.
     >      (e_num(Ca)==1).and.(e_num(Al)==1)) then
          print*,"... exchanging Al for Si"
          do i=1,Nind
            if (Iindex(i)==Al) exit
          enddo  
          do j=Nind+1,Nall
            if (Iindex(j)==Si) exit
          enddo  
          Iindex(i) = Si
          Iindex(j) = Al
          e_act(Al) = .false.
          e_act(Si) = .true.
        endif   
        if (active(iTi4O7).and.active(iTiC).and.active(iCa2Al2SiO7).and.
     >      (e_num(Ti)==2).and.(e_num(Ca)==1).and.(e_num(Si)==1).and.
     >      e_act(Si).and.(.not.e_act(O)).and.(.not.e_act(C))) then
          print*,"... exchanging Si for O"
          do i=1,Nind
            if (Iindex(i)==Si) exit
          enddo  
          do j=Nind+1,Nall
            if (Iindex(j)==O) exit
          enddo  
          Iindex(i) = O
          Iindex(j) = Si
          e_act(Si) = .false.
          e_act(O) = .true.
        endif   
        if (active(iCaMgSi2O6).and.active(iCaAl2Si2O8).and.
     >      (.not.active(iCaTiSiO5)).and.e_act(Si).and.
     >      (.not.e_act(Mg))) then
          print*,"... exchanging Si for Mg"
          do i=1,Nind
            if (Iindex(i)==Si) exit
          enddo  
          do j=Nind+1,Nall
            if (Iindex(j)==Mg) exit
          enddo  
          Iindex(i) = Mg
          Iindex(j) = Si
          e_act(Si) = .false.
          e_act(Mg) = .true.
        endif   
        if (active(iAl2O3).and.active(iMgAl2O4).and..not.e_act(Mg)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Mg"
          print*,(trim(elnam(Iindex(j))),j=1,Nall)
          print*,Nact
          do i=Nact+1,Nall
            if (Iindex(i)==Mg) exit
          enddo  
          Iindex(i) = Iindex(Nact)
          e_act(Iindex(Nact)) = .false.
          Iindex(Nact) = Mg
          e_act(Mg) = .true.          
          print*,(trim(elnam(Iindex(j))),j=1,Nall)
        endif
        if (active(iSiO).and..not.e_act(Si)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Si"
          do i=Nact+1,Nall
            if (Iindex(i)==Si) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Si
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Si) = .true.
        endif   
        if (active(iAl2O3).and..not.e_act(Al)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Al"
          do i=Nact+1,Nall
            if (Iindex(i)==Si) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Al
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Al) = .true.
        endif  
        if (active(iMgAl2O4).and..not.e_act(Mg).and..not.e_act(Al)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Al"
          do i=Nact+1,Nall
            if (Iindex(i)==Al) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Al
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Al) = .true.
        endif  
        if (active(iCaMgSi2O6).and.e_act(Ca).and.e_act(Si).and.
     >      e_num(Ca)==1.and.e_num(Si)==1) then
          do i=1,Nind
            if (Iindex(i)==Ca) exit
          enddo  
          do j=1,Nind
            if (Iindex(j)==Si) exit
          enddo  
          if (i>j) then
            e_act(Iindex(i)) = .false.
            j=i
          else
            e_act(Iindex(j)) = .false.
            i=j
          endif
          print*,"... exchanging "//elnam(Iindex(i))//
     >                     " for "//elnam(Iindex(Nact+1))
          swap = Iindex(i)
          Iindex(i) = Iindex(Nact+1) 
          Iindex(Nact+1) = swap
          e_act(Iindex(i)) = .true.
        endif   
        if (active(iCaTiO3).and.e_act(Ca).and.e_act(Ti).and.
     >      e_num(Ca)==1.and.e_num(Ti)==1) then
          do i=1,Nind
            if (Iindex(i)==Ca) exit
          enddo  
          do j=1,Nind
            if (Iindex(j)==Ti) exit
          enddo  
          if (i>j) then
            e_act(Iindex(i)) = .false.
            j=i
          else
            e_act(Iindex(j)) = .false.
            i=j
          endif
          print*,"... exchanging "//elnam(Iindex(i))//
     >                     " for "//elnam(Iindex(Nact+1))
          swap = Iindex(i)
          Iindex(i) = Iindex(Nact+1) 
          Iindex(Nact+1) = swap
          e_act(Iindex(i)) = .true.
        endif   
        if (active(iFeS).and..not.e_act(S).and.e_num(Fe)>1) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for S"
          do i=Nact+1,Nall
            if (Iindex(i)==S) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = S
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(S) = .true.
        endif   
        if (active(iFeS).and.e_act(Fe).and.e_act(S).and.
     >      e_num(Fe)==1.and.e_num(S)==1) then
          do i=1,Nind
            if (Iindex(i)==Fe) exit
          enddo  
          do j=1,Nind
            if (Iindex(j)==S) exit
          enddo  
          if (i>j) then
            e_act(Iindex(i)) = .false.
            j=i
          else
            e_act(Iindex(j)) = .false.
            i=j
          endif
          print*,"... exchanging "//elnam(Iindex(i))//
     >                     " for "//elnam(Iindex(Nact+1))
          swap = Iindex(i)
          Iindex(i) = Iindex(Nact+1) 
          Iindex(Nact+1) = swap 
          e_act(Iindex(i)) = .true.
        endif
        if (active(iFeS_l).and.e_act(Fe).and.e_act(S).and.
     >      e_num(Fe)==1.and.e_num(S)==1) then
          do i=1,Nind
            if (Iindex(i)==Fe) exit
          enddo  
          do j=1,Nind
            if (Iindex(j)==S) exit
          enddo  
          if (i>j) then
            e_act(Iindex(i)) = .false.
            j=i
          else
            e_act(Iindex(j)) = .false.
            i=j
          endif
          print*,"... exchanging "//elnam(Iindex(i))//
     >                     " for "//elnam(Iindex(Nact+1))
          swap = Iindex(i)
          Iindex(i) = Iindex(Nact+1) 
          Iindex(Nact+1) = swap
          e_act(Iindex(i)) = .true.
        endif   
        if (active(iMnS).and.e_act(Mn).and.e_act(S).and.
     >      e_num(Mn)==1.and.e_num(S)==1) then
          do i=1,Nind
            if (Iindex(i)==Mn) exit
          enddo  
          do j=1,Nind
            if (Iindex(j)==S) exit
          enddo  
          if (i>j) then
            e_act(Iindex(i)) = .false.
            j=i
          else
            e_act(Iindex(j)) = .false.
            i=j
          endif
          print*,"... exchanging "//elnam(Iindex(i))//
     >                     " for "//elnam(Iindex(Nact+1))
          swap = Iindex(i)
          Iindex(i) = Iindex(Nact+1) 
          Iindex(Nact+1) = swap
          e_act(Iindex(i)) = .true.
        endif
        if (active(iLiCl).and.e_act(Li).and.e_act(Cl).and.
     >      e_num(Li)==1.and.e_num(Cl)==1) then
          do i=1,Nind
            if (Iindex(i)==Li) exit
          enddo  
          do j=1,Nind
            if (Iindex(j)==Cl) exit
          enddo  
          if (i>j) then
            e_act(Iindex(i)) = .false.
            j=i
          else
            e_act(Iindex(j)) = .false.
            i=j
          endif
          print*,"... exchanging "//elnam(Iindex(i))//
     >                     " for "//elnam(Iindex(Nact+1))
          swap = Iindex(i)
          Iindex(i) = Iindex(Nact+1) 
          Iindex(Nact+1) = swap
          e_act(Iindex(i)) = .true.
        endif   
        if (active(iH2O_l).and..not.e_act(O).and..not.e_act(H)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for O"
          do i=Nact+1,Nall
            if (Iindex(i)==O) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = O
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(O) = .true.
        endif   
        if (active(iH2O).and..not.e_act(O).and..not.e_act(H)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for O"
          do i=Nact+1,Nall
            if (Iindex(i)==O) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = O
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(O) = .true.
        endif
        if (active(iFeO_l).and..not.e_act(Fe).and..not.e_act(O)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Fe"
          do i=Nact+1,Nall
            if (Iindex(i)==Fe) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Fe
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Fe) = .true.
        endif   
        if (active(iMgO).and..not.e_act(Mg).and..not.e_act(O)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Mg"
          do i=Nact+1,Nall
            if (Iindex(i)==Mg) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Mg
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Mg) = .true.
        endif   
        if (active(iNaAlSi3O8).and..not.e_act(Na).and.e_num(Na)==1) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Na"
          do i=Nact+1,Nall
            if (Iindex(i)==Na) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Na
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Na) = .true.
        endif   
        if (active(iCaSiO3).and.active(iCa2MgSi2O7).and.
     >      active(iCaAl2Si2O8).and..not.e_act(Mg)) then
          print*,"... exchanging "//elnam(Iindex(Nact))//" for Mg"
          do i=Nact+1,Nall
            if (Iindex(i)==Mg) exit
          enddo  
          swap = Iindex(Nact)
          Iindex(Nact) = Mg
          Iindex(i) = swap
          e_act(swap) = .false.
          e_act(Mg) = .true.
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
                mat(dk,:) = -vec(:)/REAL(stoich(eq,unknown),kind=qp)
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
     >                 " unknowns ",99(A18))',Nunsolved,Nvar1+Nvar2,
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
              if (info.ne.0) then
                print*,"*** singular matrix in QGEFA: info=",info
                stop
              endif   
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
          deps1 = +1.Q-6*eps(el)            ! limited by el abundance
          deps2 = -1.Q-6*eps(el)            ! limited by el abundance
          if (T<Tfast) then
            deps1 = +1.Q-12*eps(el)         ! quadrupole precision chemistry calls
            deps2 = -1.Q-12*eps(el) 
          endif  
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
            print'(A18,99(1pE11.3))',dust_nam(dk),DF(ii,1:Nsolve),FF(ii)
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
        limited = (fac<1.Q0)
        !if (iminoff>0.and.(iminoff.ne.laston)) then
        if (iminoff>0) then
          print*,"switch off ",dust_nam(iminoff) 
          active(iminoff) = .false.
          lastit = -99
          !if (iminoff.eq.laston) then
          !  print*,"=> fall back"
          !  active = save_active
          !  eps = save_eps
          !  ddust = save_ddust
          !endif  
        endif

        !------------------------------------
        ! ***  apply dx to ddust and eps  ***
        !------------------------------------
        do ii=1,Nsolve
          i = act_to_elem(ii) 
          el = Iindex(i)
          eps(el) = eps(el) + dx(ii)          ! direct effect
        enddo  
        do j=1,Ndep
          del = 0.Q0 
          do ii=1,Nsolve
            i = act_to_elem(ii) 
            del = del + conv(j,i)*dx(ii)      ! effect of indep.element i
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
        if (verbose>1.or.worst>1.Q-8) write(*,*) 
     >     "element conservation error 2:",worst
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
      if (qual<1.Q-10.and.(.not.CzuOvari)) then
        call PUT_DATA(nHtot,T,eps,ddust,qread,iread,active)
      endif  
      it_tot  = it_tot + it
      ieqcond = ieqcond+1

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
      logical,intent(in) :: active(0:NDUST)
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
      subroutine TRANSFORM(i1,i2,del,fac,ddust,eps,dscale,active,ok)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam,
     >                    eps0
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i1,i2
      real(kind=qp),parameter :: dsmall=1.Q-30
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp),intent(in) :: del,fac,dscale(NDUST)
      real(kind=qp) :: check(NELEM),worst
      logical,intent(inout) :: active(0:NDUST)
      logical,intent(inout) :: ok
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
      write(*,*) "element conservation error 3:",worst
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
