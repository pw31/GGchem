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
      real(kind=qp),dimension(NELEM) :: eps_save,vec,xstep,Iabund,work
      real(kind=qp),dimension(NELEM) :: scale,bvec
      real(kind=qp),dimension(NDUST) :: ddustread,dscale,pot,dust_save
      real(kind=qp),dimension(NDUST) :: Sat0,Sat1,Sat2,xvec,slin
      real(kind=qp),dimension(NELEM,NELEM) :: DF,DFsav,emat,vecs
      real(kind=qp),dimension(NDUST,NELEM) :: mat
      real(kind=qp),dimension(NELEM,NDUST) :: AA
      real(kind=qp) :: worst,xmin,Smax,Smin,qual,qold,SQUAL,del
      real(kind=qp) :: turnon,maxon,minoff,fac,fac2,amount,Nt
      real(kind=qp) :: deps1,deps2,deps,esum,emax,NRstep
      real(kind=qp) :: det(2),converge(5000,NELEM),crit,cbest
      real(kind=qp) :: deplete1,deplete2,small=1.Q-30
      real(kind=qp) :: dterm,dtmax,dlim,target,Qfinish,Sfinish
      integer,parameter :: itmax=5000
      integer,dimension(NELEM) :: elem,Nslot,eind
      integer,dimension(NDUST) :: dind,dlin,switchedON,switchedOFF
      integer,dimension(NELEM,NDUST) :: dustkind,stoich
      integer :: it,i,j,m,el,el2,Nact,Nact_read,Neq,slots,sl,dk,eq
      integer :: itry,knowns,unknowns,unknown,ii,jj,lastit,laston=0
      integer :: dtry,dtry_last=0,imaxon,iminoff,info,ipvt(NELEM)
      integer :: e_num(NELEM),e_num_save(NELEM)
      integer :: Nunsolved,unsolved(NELEM),Nvar1,Nvar2,var(NELEM)
      integer :: Nsolve,ebest,dbest,nonzero,itrivial,iread,ioff,method
      integer :: ifail,Nall,imax,swap,irow,erow,Eact,Nlin,iback,e1,e2
      integer :: act_to_elem(NELEM),act_to_dust(NELEM)
      integer :: Nzero,Ntrivial,etrivial(NELEM),dtrivial(NELEM)
      logical,dimension(NELEM) :: e_resolved,e_act,e_taken,is_esolved
      logical,dimension(NELEM) :: e_eliminated,eblocked
      logical,dimension(0:NDUST) :: active,act_read,act_old
      logical,dimension(NDUST) :: is_dsolved,d_resolved,d_eliminated
      logical,dimension(NDUST) :: itried
      logical :: action,changed,solved,limited,ok,found,all_two
      logical :: limdust,dtry_break
      character(len=1) :: char1,txt0
      character(len=2) :: rem,tnum
      character(len=6) :: dum6
      character(len=500) :: txt,txt1,txt2,text
      logical,save :: firstCall=.true.
      integer,save :: iAl2O3=0,iFe=0,iFeS=0,iNa2SiO3=0,iMgSiO3=0
      integer,save :: iMg2SiO4=0,iTi4O7=0,iCaSiO3=0,iCaMgSi2O6=0
      integer,save :: iNaAlSi3O8=0,iMgAl2O4=0,iCaTiO3=0,iSiO=0,iSiO2=0
      integer,save :: iTiO2=0,iMgTi2O5=0,iSiC=0,iCaS=0,iFe2SiO4=0,iFeO=0
      integer,save :: iNaCl=0,iKCl=0,iKAlSi3O8=0,iFe_l=0,iH2O=0,iH2O_l=0
      integer,save :: iFeS_l=0,iNaCl_l=0,iTiO2_l=0,iSiO2_l=0,iCaSO4=0
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
      integer,save :: iVO=0,iV2O3=0,iCr=0,iH2SO4=0,iS2=0
      integer,save :: iNi=0,iNi_l,iNi3S2=0,iFe3O4=0,iKMg3AlSi3O12H2=0
      integer,save :: iKFe3AlSi3O12H2=0,iMg3Si2O9H4=0,iFe3Si2O9H4=0
      integer,save :: iMgCr2O4=0,iCr2O3=0,iMn3Al2Si3O12=0,iMn2SiO4=0
      integer,save :: iKAlSi2O6=0,iCa3Al2Si3O12=0,iFeAl2SiO7H2=0
      integer,save :: iNaMg3AlSi3O12H2=0,iNaAlSiO4=0,iCa2MgSi2O7=0
      integer,save :: iCa2Al2SiO7=0
      integer,save :: iCaTiSiO5=0,iNaAlSi2O6=0,iKAlSiO4=0,iMg3Si4O12H2=0
      integer,save :: iCa3MgSi2O8=0,iCaMgSiO4=0
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
          if (dust_nam(i).eq.'CaSO4[s]')      iCaSO4=i
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
          if (dust_nam(i).eq.'H2SO4[s]')      iH2SO4=i
          if (dust_nam(i).eq.'S2[s]')         iS2=i
        enddo
        firstCall = .false. 
      endif

      !open(unit=1,file='Last_abund.in',status='replace')
      !do i=1,NEPS
      !  el = elnr(i)
      !  write(1,'(A2,2x,0pF25.20)') elnam(el),12.0+LOG10(eps0(el))
      !enddo
      !write(1,*) nHtot,T
      !close(1)

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

      xstep(:) = 0.Q0             
      call SUPER(nHtot,T,xstep,eps,Sat0,.false.) ! from scratch
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
      if (verbose>0) write(97,*)
      act_old = active
      lastit = -99
      iminoff = 0
      limited = .false.
      ifail = 0
      
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
            esum = esum**1.75   ! simple condensates first
            pot(i)  = 1.0/(0.0*emax+1.0*esum+1.0*switchedOFF(i))
            Sat1(i) = Sat0(i)**pot(i)
            !print'(A20,3(0pF8.3),1pE10.3)',
     >      !     dust_nam(i),emax,esum,pot(i),Sat1(i)
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
          if (qual>maxon) then  ! keep on iterating without switching on yet
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

        do i=1,NDUST
          if (active(i).and.(.not.act_old(i))) then
            laston = i
            switchedON(i) = switchedON(i)+1
            if (verbose>0) then
              write(97,*) it,"on  ",dust_nam(i),switchedON(i)
            endif  
          endif
        enddo
  
        if (changed) then
          Nact = 0 
          do i=1,NDUST
            if (active(i)) Nact=Nact+1
          enddo
          xstep(:)= 0.Q0             
          call SUPER(nHtot,T,xstep,eps,Sat0,NewFastLevel<1)
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
        if (Nind-1<Nact) stop "*** Nind<Nact in equil_cond."
        Nall = Nind-1
        Nind = Nact                         ! truncate at number of condensates
        if (verbose>1) print'(99(A3))',(elnam(Iindex(j)),j=1,Nall)
        if (verbose>1) print'(99(I3))',e_num(Iindex(1:Nall)) 

        !-----------------------------------------------
        ! ***  check and correct choice of elements  ***
        !-----------------------------------------------
        eblocked(:) = .false.
        do dk=1,NDUST
          if (.not.active(dk)) cycle
          e1 = 0
          e2 = 0
          do j=1,dust_nel(dk)
            el=dust_el(dk,j)
            if (e_num(el)==1) then
              e2 = e1
              e1 = el
            endif
          enddo
          if (e1>0.and.e2>0) then
            eblocked(e1) = .true.
            eblocked(e2) = .true.
            if (verbose>1) then
              print*,"*** blocked ",dust_nam(dk),elnam(e1),elnam(e2)
            endif  
          endif  
        enddo  
        e_num_save(:) = e_num(:)
        dtry = 0
 200    continue
        if (Nact<Nall) then
          e_num(:) = e_num_save
          e_eliminated(:) = .false.
          d_eliminated(:) = .false.
          do 
            if (verbose==2) then
              txt  = ''
              txt1 = ''
              txt2 = ''
              do dk=1,NDUST
                if (active(dk).and.(.not.d_eliminated(dk))) then
                  txt = trim(txt)//" "//trim(dust_nam(dk))
                endif  
              enddo
              do i=1,Nall
                el = Iindex(i)
                if (.not.e_eliminated(el)) then
                  txt1 = trim(txt1)//" "//elnam(el)
                  write(tnum,'(I2)') e_num(el) 
                  txt2 = trim(txt2)//" "//tnum
                endif  
              enddo
              print*,trim(txt)
              print*," "//trim(txt1)
              print*,trim(txt2)
            endif  
            found = .false. 
            do i=1,Nact
              el = Iindex(i)
              if (e_eliminated(el)) cycle
              if (e_num(el)==1) then
                do dk=1,NDUST
                  if (.not.active(dk)) cycle
                  if (d_eliminated(dk)) cycle
                  do j=1,dust_nel(dk)
                    if (el.eq.dust_el(dk,j)) then
                      if (verbose==2) then
                        print*,elnam(el)//" "//trim(dust_nam(dk)) 
                      endif  
                      found = .true.
                      exit
                    endif   
                  enddo  
                  if (found) exit
                enddo   
              endif   
              if (found) exit
            enddo  
            if (.not.found) exit
            e_eliminated(el) = .true.
            d_eliminated(dk) = .true.
            do j=1,dust_nel(dk)
              el = dust_el(dk,j)
              e_num(el) = e_num(el)-1
            enddo   
          enddo   
          !--- is there still a selected element with e_num=0?  ---
          found = .false.
          do i=1,Nact
            el = Iindex(i)
            if (e_eliminated(el)) cycle
            if (e_num(el)==0) then
              found = .true. 
              exit
            endif   
          enddo
          if (found) then
            found = .false. 
            do j=Nact+1,Nall 
              el = Iindex(j)
              if (e_num(el)>0) then
                found=.true.
                exit
              endif
            enddo  
            if (found) then
              if (verbose>=0) then
                print*,"... exchanging1 "//elnam(Iindex(i))//
     >               " for "//elnam(Iindex(j))
              endif  
              swap = Iindex(i)   
              Iindex(i) = Iindex(j)
              Iindex(j) = swap
              e_act(Iindex(i)) = .true.
              e_act(Iindex(j)) = .false.
            endif
            if (.not.found) then
              print*,"*** no alternative element selection found."
              stop 
            endif   
            goto 200 
          endif 
          !--- is there an unselected element with e_num=1?
          found = .false.
          !print*,"Nall,Nact=",Nall,Nact
          !if (Nall>Nact+1) then
          if (Nall>Nact) then
            do i=Nact+1,Nall
              el = Iindex(i)
              if (e_num(el)==1) then
                found = .true. 
                exit
              endif   
            enddo          
          endif
          !if (found) then
          !  all_two = (Nact>1)
          !  do m=1,Nact
          !    el2 = Iindex(m) 
          !    if (e_eliminated(el2)) cycle
          !    if (e_num(el2).ne.2) all_two=.false. 
          !  enddo
          !  if (all_two) found=.false.
          !endif  
          if (found) then
            found = .false. 
            do j=Nact,1,-1
              el = Iindex(j)
              if (e_num(el)>0) then 
                found = .true.
                exit
              endif 
            enddo
            deplete1 = eps(Iindex(i))**2/eps0(Iindex(i))
            deplete2 = eps(Iindex(j))**2/eps0(Iindex(j))
            if (verbose>1) then
              print*,"deplete",elnam(Iindex(i)),elnam(Iindex(j))
              print*,deplete1,deplete2
            endif  
            if (dtry==0.and.deplete1>1.Q+2*deplete2) found=.false.
            if (found) then
              if (verbose>=0) then
                print*,"... exchanging2 "//elnam(Iindex(j))//
     >               " for "//elnam(Iindex(i))
              endif  
              swap = Iindex(i)   
              Iindex(i) = Iindex(j)
              Iindex(j) = swap
              e_act(Iindex(i)) = .true.
              e_act(Iindex(j)) = .false.
              goto 200 
            endif  
          endif
          !--- special cases ---
          found = (Nall==Nact+1.and.Iindex(Nall)==H.and.Iindex(Nact)==O)
     >       .and.active(iFe3O4).and.active(iFe2SiO4)
     >       .and.active(iFeAl2SiO7H2).and.active(iMg3Si2O9H4)
     >       .and.active(iCaMgSi2O6).and.active(iCa3Al2Si3O12)     
          if (found) then
            ! there is a linear-combination disregarding hydrogen
            i = Nall
            j = Nact
            if (verbose>=0) then
              print*,"... exchanging3 "//elnam(Iindex(j))//
     >             " for "//elnam(Iindex(i))
            endif  
            swap = Iindex(i)   
            Iindex(i) = Iindex(j)
            Iindex(j) = swap
            e_act(Iindex(i)) = .true.
            e_act(Iindex(j)) = .false.
            goto 200 
          endif
          !print*,Nall,Nact,elnam(Iindex(Nall)),e_num(Iindex(Nall))
          found = (Nall>=Nact+1.and.Iindex(Nall)==O
     >         .and.e_num(Iindex(Nall))==1)
     >         .and.active(iCaS).and.active(iCaSO4)
          if (found) then
            found = .false.
            do j=Nact,1,-1
              !print*,j,elnam(Iindex(j)),S,e_num(Iindex(j))
              if (Iindex(j)==S.and.e_num(Iindex(j))==2) then
                found = .true.
                exit
              endif 
            enddo
          endif  
          if (found) then
            ! can't have O as independent variable in this case
            i = Nall
            if (verbose>=0) then
              print*,"... exchanging4 "//elnam(Iindex(j))//
     >           " for "//elnam(Iindex(i))
            endif  
            swap = Iindex(i)   
            Iindex(i) = Iindex(j)
            Iindex(j) = swap
            e_act(Iindex(i)) = .true.
            e_act(Iindex(j)) = .false.
            stop
            goto 200 
          endif  
        endif   

 210    continue
        dtry_break = .false.
        if (verbose>1) then
          print*,"dependent elements ... ",
     >         (elnam(Iindex(i))//" ",i=Nind+1,Nall)
          print*,"solving for ... ",
     >         (elnam(Iindex(i))//" ",i=1,Nind)
          print'(99(1pE11.3))',(Iabund(i),i=1,Nind)
        endif

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
              dtry = dtry+1
              dtry_break = .true.
              if (dtry<5) exit
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
                dtry = dtry+1
                dtry_break = .true.
                if (dtry<5) exit
                do i=1,Nunsolved
                  print'(99(1pE12.3))',(DFsav(i,j),j=1,Nunsolved)
                enddo
                stop
              endif   
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
            if (dtry>0.and.dtry_break) exit
          endif    
          if (itry==100) stop "*** itry==100"
        enddo
        if (.not.solved) then
          if (dtry==1) goto 200   ! may work by relaxing the depletion-criterium
          i = Nind
          dtry_last = dtry_last+1
          do j=Nind+dtry_last,Nall
            if (.not.eblocked(Iindex(j))) exit
          enddo
          !print*,i,j,Nind,Nall
          !print*,elnam(Iindex(i))
          !print*,elnam(Iindex(j))
          if (j<=Nall) then
            if (verbose>1) then
              print*,"... dtry_last=",dtry_last
              print*,"... exchanging5 "//elnam(Iindex(i))//
     >               " for "//elnam(Iindex(j))
            endif  
            swap = Iindex(i)   
            Iindex(i) = Iindex(j)
            Iindex(j) = swap
            e_act(Iindex(i)) = .true.
            e_act(Iindex(j)) = .false.
            goto 210
          endif
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
          deps1 = +1.Q-7*eps(el)            ! limited by el abundance
          deps2 = -1.Q-7*eps(el)            ! limited by el abundance
          if (T<Tfast) then
            deps1 = +1.Q-15*eps(el)         ! quadrupole precision chemistry calls
            deps2 = -1.Q-15*eps(el) 
          endif  
          do i=1,Ndep
            if (conv(i,j)==0.Q0) cycle
            if (is_dust(i)) cycle
            el2 = Dindex(i)                 ! limited by dep. element?
            del = 1.Q-7*eps(el2)/conv(i,j)
            if (T<Tfast) del=1.Q-15*eps(el2)/conv(i,j)
            if (del>0.Q0) deps2=MAX(deps2,-del)
            if (del<0.Q0) deps1=MIN(deps1,-del)
            !if (verbose>1) print*,elnam(el)//" "//elnam(el2),
     >      !               REAL((/conv(i,j),deps1,deps2/))
          enddo
          deps = deps2
          if (ABS(deps1)>ABS(deps2)) deps=deps1
          scale(j) = eps(el)
          if (T>Tfast) then
            dlim = 1.Q-12
            target = 1.Q-5
          else
            dlim = 1.Q-30
            target = 1.Q-10
          endif  
          do
            xstep(:) = 0.Q0
            xstep(j) = deps
            call SUPER(nHtot,T,xstep,eps,Sat2,NewFastLevel<1)
            dtmax = 0.Q0
            do ii=1,Nsolve
              i  = act_to_dust(ii) 
              dk = Dindex(i)
              !DF(ii,jj) = (Sat2(dk)-Sat0(dk))/deps*scale(j)
              dterm = LOG(Sat0(dk)/Sat2(dk))
              dtmax = MAX(dtmax,ABS(dterm))
              DF(ii,jj) = dterm/deps*scale(j)
            enddo  
            !print*,"JAC:",elnam(el),REAL(eps(el)),REAL(deps),REAL(dtmax)
            if (dtmax<target.or.ABS(deps)<dlim*eps(el)) exit
            !--- decrease deps to get more precice DF-entry ---
            deps = deps * 2.Q-1
          enddo  
        enddo            
        !if (verbose>1) then
        !  print'(12x,99(A11))',elnam(Iindex(act_to_elem(1:Nsolve)))
        !  do ii=1,Nsolve
        !    i  = act_to_dust(ii) 
        !    dk = Dindex(i)
        !    print'(A18,99(1pE11.3))',dust_nam(dk),DF(ii,1:Nsolve),FF(ii)
        !  enddo  
        !endif

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
        limdust = .false.
        do ii=1,Nsolve
          i  = act_to_elem(ii) 
          el = Iindex(i)
          if (eps(el)+dx(ii)<0.05*eps(el)) then
            fac2 = (-0.95*eps(el))/dx(ii)               ! eps+fac*dx = 0.05*eps
            if (verbose>0) print'(" *** limiting element1 ",A2,
     >        " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
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
            if (ddust(dk)+del<0.Q0) 
     >        print'("dk,laston,it,lastit=",2(A20),2(i4))',
     >             trim(dust_nam(dk)),trim(dust_nam(laston)),it,lastit
            if (del<0.Q0.and.ddust(dk)>0.1*dscale(dk)) then
              fac2 = (-ddust(dk)+0.05*dscale(dk))/del   ! ddust+fac*del = 0.05*dscale
              if (fac2<1.0.and.verbose>0) print*,"*** limiting dust 1 "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2
                iminoff = 0
                limdust = .true.
              endif  
            else if (ddust(dk)+del<0.Q0.and.
     >               dk==laston.and.it<lastit+5) then
              fac2 = -0.9*ddust(dk)/del                 ! ddust+fac*del = ddust/10
              if (fac2<1.0.and.verbose>0) print*,"*** limiting dust 2 "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2
                iminoff = 0
                limdust = .true.
              endif  
            else if (ddust(dk)+del<0.Q0) then
              fac2 = (-ddust(dk)-small*dscale(dk))/del  ! ddust+fac*del = -small*dscale
              if (fac2<1.0.and.verbose>0) print*,"*** limiting dust 3 "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2 
                iminoff = dk
                limdust = .true.
              endif
            endif  
          else  
            el = Dindex(j)
            if (eps(el)+del<0.05*eps(el)) then
              fac2 = (-0.95*eps(el))/del                ! eps+fac*dx = 0.05*eps
              if (verbose>0) print'(" *** limiting element2 ",A2,
     >        " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
              if (fac2<fac) then
                fac = fac2 
                iminoff = 0
                limdust = .false.
              endif  
            endif
          endif  
        enddo  
        dx = dx*fac
        limited = (fac<1.Q0)
        !if (iminoff>0.and.(iminoff.ne.laston)) then
        if (iminoff>0) then
          if (verbose>=0) print*,"switch off ",dust_nam(iminoff) 
          active(iminoff) = .false.
          lastit = -99
          switchedOFF(iminoff) = switchedOFF(iminoff)+1
          if (verbose>0) then
            write(97,*) it,"off ",dust_nam(iminoff),switchedOFF(iminoff)
          endif  
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
        eps_save = eps
        dust_save = ddust
        NRstep = 1.Q0
        qold = qual
        do iback=1,8
          eps = eps_save
          ddust = dust_save
          do ii=1,Nsolve
            i = act_to_elem(ii) 
            el = Iindex(i)
            eps(el) = eps(el) + NRstep*dx(ii) ! direct effect
          enddo  
          do j=1,Ndep
            del = 0.Q0 
            do ii=1,Nsolve
              i = act_to_elem(ii) 
              del = del + conv(j,i)*dx(ii) ! effect of indep.element i
            enddo 
            if (is_dust(j)) then
              dk = Dindex(j)
              ddust(dk) = ddust(dk) + NRstep*del
            else  
              el = Dindex(j)
              eps(el) = eps(el) + NRstep*del
            endif  
          enddo
          xstep(:) = 0.Q0
          call SUPER(nHtot,T,xstep,eps,Sat0,NewFastLevel<1)
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
        !del = 0.Q0
        !do i=1,NELM
        !  if (i==iel) cycle
        !  el = elnum(i)
        !  del = MAX(del,ABS(LOG10(work(el)/eps(el))))
        !  print'(A3,1pE10.3," -> ",1pE10.3,0pF5.2)',
     >  !       elnam(el),work(el),eps(el),LOG10(work(el)/eps(el))
        !enddo
        !print*,"epsilon max-log-change =",REAL(del)

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

        !xstep(:) = 0.Q0
        !call SUPER(nHtot,T,xstep,eps,Sat0,NewFastLevel<1)
        !qual = SQUAL(Sat0,active)
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
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      use CONVERSION,ONLY: Ndep,Nind,Dindex,Iindex,is_dust,conv
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in)  :: nHtot,T
      real(kind=qp),intent(in) :: xx(NELEM),eps(NELEM)
      real(kind=qp),intent(out) :: Sat(NDUST)
      logical,intent(in) :: merk
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

      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        if (eps1(el).le.0.Q0) then
          write(*,*) "*** negative el.abund. SUPER",elnam(el),eps1(el)
          stop
        endif  
      enddo
      
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
