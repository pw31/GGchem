***********************************************************************
      SUBROUTINE ADAPT_CONDENSATES
***********************************************************************
      use PARAMETERS,ONLY: nHmax,Tmax,pmax,model_pconst,model_eqcond,
     >                     verbose,Nseq,Tseq,adapt_file
      use CHEMISTRY,ONLY: NMOLE,NELM,m_kind,m_anz,elnum,cmol,el
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,
     >                    muH,mass,mel,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol,H,C,N,O,Si
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: SPMAX=100
      real(kind=qp) :: eps(NELEM),eps00(NELEM)
      real(kind=qp) :: Sat(NDUST),eldust(NDUST)
      real*8 :: Tg,nHges,p,mu,muold,pgas,ngas,ntot,nn,ctot
      real*8 :: ff,dfdmu,dmu,fold,conc,q,qual,deps,fac,emax
      real*8,dimension(SPMAX) :: CCadapt,delta,left,right,estep
      integer,dimension(SPMAX) :: SPadapt
      integer :: verb,iseq,it,ad,i,j,e,Nadapt,worst,oworst,iter
      logical,dimension(SPMAX) :: is_atom
      logical :: bracketed
      character(len=20) :: mol
      character(len=1) :: char

      !--------------------------------------------------
      ! ***  read the gas concentrations to adapt to  ***
      !--------------------------------------------------
      print*,"read concentrations to adapt to from "//trim(adapt_file)
      open(unit=1,file=adapt_file,status='old')
      Nadapt = 0
      ctot = 0.0
      do 
        read(1,*,end=1000) conc,mol
        ctot = ctot + conc/100.0
        Nadapt = Nadapt + 1
        SPadapt(Nadapt) = 0
        do j=1,NELM
          e = elnum(j)
          if (elnam(e).eq.mol) then
            SPadapt(Nadapt) = e
            CCadapt(Nadapt) = conc/100.0
            is_atom(Nadapt) = .true.
          endif
        enddo
        do j=1,NMOLE
          if (cmol(j).eq.mol) then
            SPadapt(Nadapt) = j
            CCadapt(Nadapt) = conc/100.0
            is_atom(Nadapt) = .false.
          endif
        enddo
        if (SPadapt(Nadapt)==0) then
          print*,"*** molecule not found: "//mol
          stop
        endif
        i = SPadapt(Nadapt)
        if (is_atom(Nadapt)) then
          estep(Nadapt) = eps0(SPadapt(Nadapt))
        else
          estep(Nadapt) = 9.E+99
          do j=1,m_kind(0,i)
            if (m_kind(j,i)==el) cycle ! avoid free electrons
            e = elnum(m_kind(j,i))
            estep(Nadapt) = MIN(estep(Nadapt),eps0(e)/m_anz(j,i))  
            !print'(A2,I3,1pE13.5)',elnam(e),m_anz(j,i),eps0(e)
          enddo
        endif
        estep(Nadapt) = MIN(0.1,0.1*estep(Nadapt))
        if (is_atom(Nadapt)) then
          print'("fit at  ",A2,8x,"conc=",1pE11.4,"  estep=",1pE9.3)',
     >       elnam(i),CCadapt(Nadapt),estep(Nadapt)
        else
          print'("fit mol ",A10,"conc=",1pE11.4,"  estep=",1pE9.3)',
     >       cmol(i),CCadapt(Nadapt),estep(Nadapt)
        endif
      enddo
 1000 CCadapt = CCadapt/ctot
      print*

      eps00 = eps0
      mu    = muH
      nHges = nHmax
      p     = pmax
      Tseq(Nseq) = Tmax
      delta = 0.0
      left  = 9.e+99
      right = 9.e+99
      worst = 0
      verb  = 0

      do iter=1,999
        !-------------------------------------------
        ! ***  add amounts of molecules to eps0  ***
        !-------------------------------------------
        eps0 = eps00
        do ad=1,Nadapt
          i = SPadapt(ad)
          if (is_atom(ad)) then
            print'(A2,8x,3(1pE15.7),"  estep=",1pE9.3)',
     >         elnam(i),left(ad),delta(ad),right(ad),estep(ad)
            eps0(i) = eps0(i) + delta(ad)
          else
            print'(A10,3(1pE15.7),"  estep=",1pE9.3)',
     >         cmol(i),left(ad),delta(ad),right(ad),estep(ad)
            do j=1,m_kind(0,i)
              if (m_kind(j,i)==el) cycle
              e = elnum(m_kind(j,i))
              eps0(e) = eps0(e) + delta(ad)*m_anz(j,i)
            enddo
          endif
        enddo
        open(unit=2,file='abund_test.in',status='replace')
        do i=1,NELM
          e = elnum(i)
          if (i==el) cycle
          write(2,*) elnam(e),LOG10(eps0(e)/eps0(1))+12.Q0
        enddo
        close(2)
        !read(*,'(A1)') char

        !---------------------------------------------
        ! ***  solve the chemical and phase equil. ***
        !---------------------------------------------
        if (iter==1.or..not.is_atom(worst)) call ERASE_DBASE
        do iseq=1,Nseq          ! possibility to set a sequence of T-points hot->cold
          Tg = Tseq(iseq)       ! (not used when default Nseq=1)
          do it=1,999
            if (model_pconst) nHges = p*mu/(bk*Tg)/muH
            if (model_eqcond) then
              call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verb)
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
        enddo

        print*
        print*,"ITERATION = ",iter
        print*,"comparison ..."
        ntot = 0.0
        do ad=1,Nadapt
          i = SPadapt(ad)
          if (is_atom(ad)) then
            ntot = ntot+nat(i)
          else
            ntot = ntot+nmol(i)
          endif
        enddo
        oworst = worst
        worst = 0
        qual  = 0.0
        do ad=1,Nadapt
          i = SPadapt(ad)
          if (is_atom(ad)) then
            q = (nat(i)/ntot-CCadapt(ad))/(nat(i)/ntot+CCadapt(ad))
            print'(A2,8x,2(1pE14.6),0pF8.4)',
     >         elnam(i),nat(i)/ntot,CCadapt(ad),q
          else
            q = (nmol(i)/ntot-CCadapt(ad))/(nmol(i)/ntot+CCadapt(ad))
            print'(A10,2(1pE14.6),0pF8.4)',
     >         cmol(i),nmol(i)/ntot,CCadapt(ad),q
          endif
          if (iter.gt.1.and.ad.ne.oworst) then
            right(ad) = 9.e+99
            left(ad) = 9.e+99
          endif
          if (q>0.0) then
            right(ad) = delta(ad)
          else
            left(ad)  = delta(ad)
          endif
          fac = 1.0
          if (ad.ne.oworst) fac=0.25    ! stay with oworst for a bit longer
          if (ABS(q)*fac>qual) then
            qual = ABS(q)*fac
            worst = ad
          endif
        enddo
        bracketed = (left(oworst)<1.e+99.and.right(oworst)<1.e+99)
        if (iter>1.and..not.bracketed) worst=oworst
        i = SPadapt(worst)
        if (is_atom(worst)) then
          print'(" it=",I4,"  worst=",A10,"  qual=",0pF8.5)',
     >         iter,elnam(i),qual
          nn = nat(i)
        else
          print'(" it=",I4,"  worst=",A10,"  qual=",0pF8.5)',
     >         iter,cmol(i),qual
          nn = nmol(i)
        endif  
        if (qual<1.E-3) exit

        !--- avoid too large estep leading to negative eps0 ---
        if (is_atom(worst)) then
          emax = eps0(i)
        else
          emax = 9.E+99
          do j=1,m_kind(0,i)
            if (m_kind(j,i)==el) cycle 
            e = elnum(m_kind(j,i))
            emax = MIN(emax,eps0(e)/m_anz(j,i))  
          enddo
        endif
        emax = emax/2

        !--- iterate ---
        if (left(worst)<1.e+99.and.right(worst)<1.e+99) then
          delta(worst) = (right(worst)+left(worst))/2
          estep(worst) = (right(worst)-left(worst))/2
        else if (nn/ntot<CCadapt(worst)) then
          delta(worst) = delta(worst) + MIN(emax,estep(worst))
          estep(worst) = MIN(emax,estep(worst)*2)
        else
          delta(worst) = delta(worst) - MIN(emax,estep(worst))
          estep(worst) = MIN(emax,estep(worst)*2)
        endif
 
      enddo  

      call DEMO_CHEMISTRY

      end      
