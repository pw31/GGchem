***********************************************************************
      PROGRAM EQ_CHEMISTRY
***********************************************************************
      use PARAMETERS,ONLY: model_dim,model_struc,model_eqcond,
     >                     useDatabase
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform,preEst,preUse,preIter
      use DATABASE,ONLY: NLAST
      implicit none

      call READ_PARAMETER
      call INIT
      call INIT_CHEMISTRY
      call INIT_DUSTCHEM
      
      if (model_dim==0) then
        call DEMO_CHEMISTRY
      else if (model_dim==1) then  
        if (model_struc==0) then 
          call DEMO_SWEEP
        else  
          call DEMO_STRUCTURE
        endif 
      else if (model_dim==2) then  
        call DEMO_PHASEDIAGRAM
      else
        print*,'*** model_dim=',model_dim,' ???'
        stop
      endif   
      
      print*
      print'("            smchem calls = ",I8)',chemcall
      print'("         iterations/call = ",0pF8.2)',
     >                     REAL(chemiter)/REAL(chemcall)
      print'("     pre-iterations/call = ",0pF12.3)',
     >                      REAL(preIter)/REAL(chemcall)
      print'("usage of saved estimates = ",0pF12.3,"%")',
     >                       REAL(preUse)/REAL(preEst)*100.0
      if (model_eqcond) then
        print'("   eq condensation calls = ",I8)',ieqcond
        print'("      eq iterations/call = ",0pF8.2)',
     >                   REAL(ieqconditer)/REAL(ieqcond)
        print'("         transform calls = ",I8)',itransform
        NLAST=0         ! also save replaced database entries
        if (useDatabase) call SAVE_DBASE
      endif

      end


***********************************************************************
      SUBROUTINE DEMO_CHEMISTRY
***********************************************************************
      use PARAMETERS,ONLY: nHmax,Tmax,pmax,model_pconst,model_eqcond,
     >                     verbose 
      use CHEMISTRY,ONLY: NMOLE,NELM,m_kind,m_anz,elnum,cmol,el
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,
     >                    muH,mass,mel,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol,C,N,O,Si
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: nmax,threshold,deps
      real*8,parameter :: pi=3.14159265358979D+0
      real*8  :: Tg,nHges,p,mu,muold,pgas,fold,ff,dfdmu,dmu,mugas,Vol
      real*8  :: rhog,rhoc,rhod,d2g,mcond,mgas,Vcon,Vcond,Vgas,ngas
      real*8  :: nkey,nkeyt,tchem,tchemtot,AoverV,mic,atyp,alpha,vth
      real*8  :: yr,Myr,molm,molmass,stoich
      integer :: i,imol,iraus,e,aIraus,aIIraus,j,verb,dk,it,stindex
      integer :: k,keyel,imax,dustst
      logical :: included,haeufig,raus(NMOLE)
      logical :: rausI(NELEM),rausII(NELEM)
      character(len=10) :: sp
      character(len=20) :: limcond


      !deps = eps0(Si)*0.1*(1.0-108./200.)
      !eps0(Si) = eps0(Si) + deps
      !eps0(O)  = eps0(O)  + deps*2
      
      Tg    = Tmax
      nHges = nHmax
      p     = pmax
      eps   = eps0
      mu    = muH
      eldust = 0.Q0
      verb = verbose
      if (model_eqcond) verb=0

      do it=1,999
        if (model_pconst) nHges = p*mu/(bk*Tg)/muH
        if (model_eqcond) then
          call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
        endif
        call GGCHEM(nHges,Tg,eps,.false.,verb)
        ngas = nel
        do j=1,NELEM
          ngas = ngas + nat(j)
        enddo
        do j=1,NMOLE
          ngas = ngas + nmol(j)
        enddo
        pgas  = ngas*bk*Tg
        ff    = p-pgas
        if (it==1) then
          muold = mu
          mu = nHges/pgas*(bk*Tg)*muH
          dmu = mu-muold
          if (.not.model_pconst) exit
        else
          dfdmu = (ff-fold)/(mu-muold)
          dmu   = -ff/dfdmu
          !write(98,'(I3,99(1pE14.7))') it,muold,mu,fold,ff,dfdmu,dmu/mu
          muold = mu
          if ((dmu>0.0).or.ABS(dmu/mu)<0.7) then
            mu = muold+dmu
          else
            mu = nHges/pgas*(bk*Tg)*muH
          endif  
        endif
        fold = ff
        print '("p-it=",i3,"  mu=",2(1pE20.12))',it,mu/amu,dmu/mu
        if (ABS(dmu/mu)<1.E-10) exit
      enddo  

      print*
      write(*,*) '----- total nuclei dens. and fractions in gas -----'
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
        write(*,'(" n<",A2,">=",1pE10.4,2x,1pE10.4)')
     >      elnam(i),nHges*eps(i),eps(i)/eps0(i)
      enddo  

      ngas = nel      ! cm-3
      rhog = nel*mel  ! g cm-3
      do j=1,NELEM
        ngas = ngas + nat(j)
        rhog = rhog + nat(j)*mass(j)
      enddo
      do j=1,NMOLE
        ngas = ngas + nmol(j)
        rhog = rhog + nmol(j)*mmol(j)
      enddo
      rhod = 0.d0     ! g cm-3  
      Vcon = 0.d0     ! cm3
      do i=1,NDUST
        if (eldust(i)<=0.Q0) cycle 
        rhod = rhod + nHges*eldust(i)*dust_mass(i)
        Vcon = Vcon + nHges*eldust(i)*dust_Vol(i)
      enddo  
      mugas = rhog/ngas             ! mean molecular weight [g]
      d2g   = rhod/rhog             ! dust/gas mass ratio [-]
      rhoc  = rhod/Vcon             ! condensed matter density [g cm-3]
      mcond = 1.0                   ! mass of condensates [= 1 g]
      Vcond = mcond/rhoc            ! volume of condensates [cm3]
      mgas  = mcond/d2g             ! mass of gas [g]
      Vgas  = mgas/mugas*bk*Tg/pgas ! volume of gas [cm3]
      
      print*
      write(*,*) '----- bulk properties -----'
      print'("for Tg[K]=",0pF8.2," and n<H>[cm-3]=",1pE10.3)',Tg,nHges
      print'(16x,3(A12))',"condensates","gas","check"
      print'("         mass[g]",2(1pE12.4))',mcond,mgas
      print'("  density[g/cm3]",2(1pE12.4))',rhoc,rhog
      print'("tot.dens.[g/cm3]",12x,2(1pE12.4))',nHges*muH,
     >                                   (mcond+mgas)/Vgas
      print'("     volume[cm3]",2(1pE12.4))',Vcond,Vgas
      print'("   pressure[bar]",12x,1pE12.4)',pgas/bar
      print'("  el.press.[bar]",12x,1pE12.4)',nel*bk*Tg/bar
      print'(" mol.weight[amu]",12x,1pE12.4)', mugas/amu
      print'("         mu[amu]",12x,2(1pE12.4))', mu/amu,
     >                (mcond+mgas)/Vgas/pgas*(bk*Tg)/amu
      print'("        muH[amu]",12x,2(1pE12.4))',muH/amu,
     >                     (mcond+mgas)/(nHges*Vgas)/amu
      
      print*
      write(*,*) '----- condensates [cm3] [mfrac] [Vfrac] -----'
      raus = .false.
      do 
        iraus = 0
        nmax  = 0.Q0
        do i=1,NDUST
          if (raus(i).or.eldust(i)<=0.Q0) cycle 
          if (eldust(i)>nmax) then
            iraus = i
            nmax = eldust(i)
          endif
        enddo
        if (iraus==0) exit
        raus(iraus) = .true.
        write(*,1020) ' n'//trim(dust_nam(iraus))//'=',
     >                eldust(iraus)*nHges,
     >                eldust(iraus)*dust_mass(iraus)*nHges/rhod,
     >                eldust(iraus)*dust_Vol(iraus)*nHges/Vcon
      enddo
  
      print*
      write(*,*) '----- atoms and ions [cm3] -----'
      write(*,1000) ' nel=',nel
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
        write(*,1010) ' n'//trim(elnam(i))//'I=',nat(i),
     >               '  n'//trim(elnam(i))//'II=',nion(i)
      enddo
  
      write(*,*) '----- most abundant species [cm3] [molfrac] -----'
      raus   = .false.
      rausI  = .false.
      rausII = .false.
      do
        iraus   = 0
        aIraus  = 0
        aIIraus = 0
        nmax  = 0.Q0
        do i=1,NMOLE
          if ((nmol(i).gt.nmax).and.(.not.raus(i))) then
            iraus = i
            nmax  = nmol(i)
          endif
        enddo
        do e=1,NELM
          if (e==el) cycle
          i = elnum(e) 
          if ((nat(i).gt.nmax).and.(.not.rausI(i))) then
            iraus = 0
            aIraus = i
            aIIraus = 0
            nmax  = nat(i)
          endif
          if ((nion(i).gt.nmax).and.(.not.rausII(i))) then
            iraus = 0
            aIraus = 0
            aIIraus = i
            nmax  = nion(i)
          endif
        enddo
        haeufig = (nmax.gt.ngas*1.Q-9)
        if (.not.haeufig) exit
        if (iraus>0) then
          raus(iraus) = .true.
          write(*,4010) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
        else if (aIraus>0) then 
          rausI(aIraus) = .true.
          write(*,4010) elnam(aIraus)//"I       ",
     >                  nat(aIraus),nat(aIraus)/ngas
        else if (aIIraus>0) then 
          rausII(aIIraus) = .true.
          write(*,4010) elnam(aIIraus)//"II     ",
     >                  nat(aIIraus),nion(aIIraus)/ngas
        endif  
      enddo
      iraus = stindex(cmol,NMOLE,'O2')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'NO')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'S3')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'S4')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'S5')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'H2SO4')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'SF6')
      if (.not.raus(iraus))
     >   write(*,4020) cmol(iraus),nmol(iraus),nmol(iraus)/ngas
          
      print*
      write(*,*) '-----  where are the elements?  -----'
      do e=1,NELM
        i = elnum(e)
        if (e==el) then
          write(*,'("    Element ",A2,1pE15.3)') 'el',0.Q0
          write(*,'(1x,A18,1pE10.3)') "nel",nel
          threshold = 1.Q-3*nel
        else   
          write(*,'("    Element ",A2,1pE15.3)') elnam(i),eps0(i)*nHges 
          threshold = eps(i)*nHges*1.D-2
          if (nat(i).gt.eps(i)*nHges*1.D-2) then
            write(*,'(1x,A18,1pE10.3)') "n"//trim(elnam(i)), nat(i) 
          endif  
        endif  

        raus = .false.
        do 
          iraus = 0
          nmax  = 0.Q0
          do dk=1,NDUST
            if (eldust(dk)<=0.Q0) cycle 
            included = .false. 
            do j=1,dust_nel(dk)
              if (i==dust_el(dk,j)) then
                included = .true.
              endif
            enddo  
            if (included) then
              if ((eldust(dk).gt.nmax).and.(.not.raus(dk))) then
                iraus = dk
                raus(dk) = .true.
                nmax = eldust(dk)
              endif  
            endif
          enddo  
          if (nmax==0.Q0) exit
          write(*,'(1x,A18,1pE10.3)') 
     >          "n"//trim(dust_nam(iraus)),eldust(iraus)*nHges 
        enddo  

        raus = .false.
        do 
          iraus = 0
          nmax  = 0.Q0
          do imol=1,NMOLE
            sp = cmol(imol) 
            if ((nmol(imol).gt.nmax).and.(.not.raus(imol))) then
              included = .false. 
              do j=1,m_kind(0,imol)
                if (e==m_kind(j,imol)) included=.true.
              enddo  
              if (included) then
                iraus = imol
                nmax = nmol(imol)
              endif  
            endif
          enddo  
          haeufig = (nmax.gt.threshold)
          if (.not.haeufig) exit
          write(*,'(1x,A18,1pE10.3)') "n"//trim(cmol(iraus)),nmol(iraus)
          raus(iraus) = .true.
        enddo
      enddo  

*     ------------------------------
      call SUPERSAT(Tg,nat,nmol,Sat)
*     ------------------------------
      print*
      write(*,*) '----- supersaturation ratios -----'
      do i=1,NDUST
        if (Sat(i)<1.Q-2) cycle 
        write(*,5000) dust_nam(i),Sat(i) 
      enddo

*     -----------------------------------------------------
*     ***  Calculation of the condenstation timescales  ***
*     -----------------------------------------------------
      print*
      print'("----- condensation timescales -----")'
      yr  = 365.25*24.0*3600.0
      Myr = 1.E+6*yr
      mic = 1.E-4                ! one micrometer [cm]
      atyp = 1.0*mic             ! typical particle size
      AoverV = 3.0*atyp          ! A/V = 4pi a^2/(4pi/3 a^3)
      alpha = 1.0                ! sticking probability
      tchemtot = 1.E-99          ! chem.timescale of slowest condensate
      do i=1,NDUST
        !--- loop over present condensates ---
        if (eldust(i)<=0.0) cycle
        write(*,*) 'n_'//trim(dust_nam(i))//' = ',eldust(i)*nHges
        !--- find least abundant element in the gas phase   ---
        !--- print stoichiometry and gas element abundances ---
        !--- nkey = min (eps(e )*nHges/s_e )  ---
        nkey = 1.E+99
        do k=1,dust_nel(i)
          e = dust_el(i,k)
          !write(*,'("    Element ",A2,1pE11.4,I2)')
     >    !        elnam(e), eps(e)*nHges, dust_nu(i,k)
          nkeyt = eps(e)*nHges/dust_nu(i,k)
          if (nkeyt<nkey) then
            nkey = nkeyt
            keyel = e
            dustst = dust_nu(i,k)
          endif
        enddo
        !--- find atom/molecule which contains most ---
        !--- of the key element                     ---
        sp = elnam(keyel)
        nmax = nat(keyel)/dustst
        molmass = mass(keyel)
        do imol=1,NMOLE
          included = .false. 
          molm = 0.0
          do j=1,m_kind(0,imol)
            if (m_kind(j,imol)==el) cycle      ! avoid free electrons
            e = elnum(m_kind(j,imol))
            if (e==keyel) then
              included = .true.
              stoich = m_anz(j,imol)
            endif  
            molm = molm + m_anz(j,imol)*mass(e)
          enddo  
          if (included.and.nmol(imol)*stoich/dustst>nmax) then
            sp = cmol(imol)
            nmax = nmol(imol)*stoich/dustst
            molmass = molm
          endif  
        enddo
        vth = SQRT(8.0*bk*Tg/(pi*molmass))     ! [cm/s]
        Vol = eldust(i)*nHges*dust_Vol(i)      ! [cm3/cm3]
        tchem = 1.d0/(vth*alpha*AoverV*Vol)    ! [s]
        write(*,'(" limiting element = ",A2)') elnam(keyel)
        write(*,'(" mostly present as ",A10)') sp
        write(*,'("        nmax,nkey = ",2(1pE11.4)," cm-3")')
     >        nmax,nkey
        write(*,'("             mass = ",2(1pE11.4)," amu")')
     >        molmass/amu,mass(keyel)/amu
        write(*,'("        timescale = ",1pE11.4," yr")')
     >        tchem/yr
        if (tchem>tchemtot) then
          tchemtot = tchem
          limcond  = dust_nam(i)
        endif
      enddo  
      write(*,'("Limiting condensate ",A22,"  timescale/yr = ",
     >          1pE11.3)') limcond, tchemtot/yr
      
 1000 format(a6,1pE9.3)
 1010 format(a6,1pE9.3,a8,1pE9.3)
 1020 format(a22,1pE15.9,2(0pF10.5))
 1030 format(a22,1pE12.4)
 4000 format(a7,1pE10.4,a5,1pE10.4)     
 4010 format(' n',a8,1pE12.4,0pF13.9)
 4020 format(' n',a8,1pE12.4,1pE13.3)
 5000 format(1x,a20,' S=',1pE9.3)
      RETURN
      end      
