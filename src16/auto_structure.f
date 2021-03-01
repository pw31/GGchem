***********************************************************************
      SUBROUTINE AUTO_STRUCTURE
***********************************************************************
      use PARAMETERS,ONLY: Rpl,Mpl,Tmax,pmin,pmax,gamma,verbose,
     >                     model_eqcond,remove_condensates,Npoints,
     >                     useDataBase,stop_after_init
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,grav,
     >                    muH,mass,mel,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol,H,C,N,O,Si
      use DATABASE,ONLY: NLAST
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      character(len=20) :: name,short_name(NDUST)
      real(kind=qp) :: eps(NELEM),eps00(NELEM),eadd
      real(kind=qp) :: Sat(NDUST),eldust(NDUST),out(NDUST)
      real(kind=qp) :: fac,e_reservoir(NELEM),d_reservoir(NDUST)
      real :: rho,Tg,p,kT,nHges,nges,pgas,ff,fold,mu,xx,xold,dx,dfdx
      real :: muold,dfdmu,dmu
      real :: p1,p2,T1,T2,g1,g2,mu1,mu2,rho1,rho2,dg,dgold,limit
      real :: zz,dz,Hp,kappa,rhog,rhod,dustV,dustM,km=1.D+5
      integer :: i,j,k,e,ip,jj,it,NOUT,jbest
      logical :: out_at(NELEM),out_mol(NMOLE),atom,ex
      logical :: outAllHistory=.false.
      
      !-----------------------------------------------------------
      ! ***  compute base point with equilibrium condensation  ***
      !-----------------------------------------------------------
      !eadd = 0.2*eps0(Si)
      !eps0(H) = eps0(H)+1.3Q0*eadd
      !eps0(O) = eps0(O)+eadd
      Tg = Tmax
      p  = pmax
      mu = muH
      inquire(file='LastMu.dat',exist=ex)
      if (ex) then
        open(12,file='LastMu.dat')
        read(12,*) mu
        close(12)
        mu = mu*muH
      endif  
      do it=1,99
        nHges = p*mu/(bk*Tg)/muH
        call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
        call GGCHEM(nHges,Tg,eps,.false.,0)
        kT   = bk*Tg
        nges = nel
        rho  = nel*mel
        do j=1,NELEM
          nges = nges + nat(j)
          rho  = rho  + nat(j)*mass(j)
        enddo
        do j=1,NMOLE
          nges = nges + nmol(j)
          rho  = rho  + nmol(j)*mmol(j)
        enddo
        pgas = nges*kT
        ff   = p-pgas
        if (it==1) then
          muold = mu
          mu = nHges/pgas*(bk*Tg)*muH
          mu = MAX(0.3*muold,MIN(3.0*muold,mu))
          dmu = mu-muold
        else
          dfdmu = (ff-fold)/(mu-muold)
          dmu   = -ff/dfdmu
          muold = mu
          if ((dmu>0.0).or.ABS(dmu/mu)<0.7) then
            mu = muold+dmu
          else
            mu = nHges/pgas*(bk*Tg)*muH
          endif  
          mu = MAX(0.3*muold,MIN(3.0*muold,mu))
        endif
        fold = ff
        print '("p-it=",i3,"  mu=",2(1pE20.12))',it,mu/amu,dmu/mu
        if (ABS(dmu/mu)<1.E-10) exit
      enddo
      open(12,file='LastMu.dat',status='replace')
      write(12,*) mu/muH
      close(12)
      
      dustV = 0.0
      dustM = 0.0
      do jj=1,NDUST
        dustV = dustV + eldust(jj)*dust_Vol(jj)
        dustM = dustM + eldust(jj)*dust_mass(jj)
      enddo
      print*
      print'("----- the planet surface consists of ... -----")'
      do jj=1,NDUST
        if (eldust(jj)==0.Q0) cycle
        print'(A18,"  Vfrac=",0pF12.8," %  Mfrac=",0pF12.8," %")',
     >         dust_nam(jj),eldust(jj)*dust_Vol(jj)/dustV*100,
     >         eldust(jj)*dust_mass(jj)/dustM*100
      enddo

      !------------------------------------------------
      ! ***  reset total <- gas element abundances  ***
      !------------------------------------------------
      nHges  = nHges*eps(H)   ! in the gas phase over the crust
      fac    = 1.Q0/eps(H)
      eps    = fac*eps        ! nomalise to have eps(H)=1
      eps0   = eps
      muH    = rho/nHges      ! in the gas phase over the crust
      mu     = rho/nges       ! in the gas phase over the crust
      eldust = 0.Q0           ! no condensates in that gas
      e_reservoir = 0.Q0
      d_reservoir = 0.Q0
      eps00 = eps0
      print*
      print'("----- gas properties ... -----")'
      do e=1,NELM
        if (e==el) cycle
        j = elnum(e) 
        print'("eps(",A2,")=",1pE11.3)',elnam(j),eps(j)
      enddo
      print*
      print'("----- most abundant molecules [mol-fractions] -----")'
      out_at=.false.
      out_mol=.false.
      do it=1,12
        limit = 0.0
        do j=1,NELEM
          if (nat(j)>limit.and.(.not.out_at(j))) then
            limit = nat(j)
            jbest = j
            atom = .true.
          endif
        enddo  
        do j=1,NMOLE
          if (nmol(j)>limit.and.(.not.out_mol(j))) then
            limit = nmol(j)
            jbest = j
            atom = .false.
          endif
        enddo
        if (limit<1.E-9*nges) exit
        if (atom) then
          print 3000,elnam(jbest),nat(jbest)/nges
          out_at(jbest) = .true.
        else  
          print 3000,cmol(jbest),nmol(jbest)/nges
          out_mol(jbest) = .true.
        endif
      enddo
      print'("mu[amu] =",0pF8.4)',mu/amu
      call SUPERSAT(Tg,nat,nmol,Sat)
      print*
      print'("----- supersaturation ratios -----")'
      do i=1,NDUST
        !if (Sat(i)<1.Q-2) cycle
        if (index(dust_Nam(i),"H2O")<=0) cycle
        write(*,'(A18,1pE9.2)') dust_nam(i),Sat(i) 
      enddo
      NLAST=0      
      if (useDatabase) call SAVE_DBASE
      if (stop_after_init) return
      
      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
      enddo
      NOUT = NELM
      if (charge) NOUT=NOUT-1
      open(unit=70,file='Static_Conc.dat',status='replace')
      write(70,1000) 'H',eps( H), 'C',eps( C),
     &               'N',eps( N), 'O',eps( O)
      write(70,*) NOUT,NMOLE,NDUST,Npoints
      write(70,2000) 'Tg','nHges','pgas','el',
     &               (trim(elnam(elnum(j))),j=1,el-1),
     &               (trim(elnam(elnum(j))),j=el+1,NELM),
     &               (trim(cmol(i)),i=1,NMOLE),
     &               ('S'//trim(short_name(i)),i=1,NDUST),
     &               ('n'//trim(short_name(i)),i=1,NDUST),
     &               ('eps'//trim(elnam(elnum(j))),j=1,el-1),
     &               ('eps'//trim(elnam(elnum(j))),j=el+1,NELM),
     &               'dust/gas','dustVol/H'
      open(unit=71,file='Structure.dat',status='replace')
      write(71,2000) 'z[km]','rho[g/cm3]','pgas[dyn/cm2]','T[K]',
     &               'n<H>[cm-3]','mu[amu]','g[m/s2]','Hp[km]'

      !-----------------------------------------------------------
      ! ***  construct a polytopic atmosphere above the crust  ***
      ! ***  in hydrostatic equilibrium using dp/dz = - rho g  ***
      ! ***  and  T = const p**kappa.                          ***
      !-----------------------------------------------------------
      kappa = (gamma-1.0)/gamma      ! polytrope index T ~ p^kappa 
      zz    = 0.0                    ! height [cm] 
      T1    = Tg                     ! temperature [K]
      p1    = pgas                   ! pressure [dyn/cm2]
      rho1  = rho                    ! gas mass density [g/cm3]
      mu1   = mu                     ! mean molecular weight [g]
      g1    = grav*Mpl/(Rpl+zz)**2   ! gravity [cm/s2]
      !useDataBase = .false.
      !verbose = 1

      do ip=1,Npoints
        Hp  = bk*T1/(mu1*g1)         ! pressure scale height

        !--- compute properties of condensates ---
        rhod  = 0.0
        dustV = 0.0
        do jj=1,NDUST
          rhod  = rhod  + nHges*eldust(jj)*dust_mass(jj)
          dustV = dustV + eldust(jj)*dust_Vol(jj)
          out(jj) = LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj))))
          if (ABS(Sat(jj)-1.Q0)<1.E-10) out(jj)=0.Q0
        enddo  
        call SUPERSAT(T1,nat,nmol,Sat)

        !--- output stuff ---
        write(70,2010) T1,nHges,p1,
     &       LOG10(MAX(1.Q-300, nel)),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=1,el-1),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=el+1,NELM),
     &      (LOG10(MAX(1.Q-300, nmol(jj))),jj=1,NMOLE),
     &      (out(jj),jj=1,NDUST),
     &      (LOG10(MAX(1.Q-300, eldust(jj))),jj=1,NDUST),
     &      (LOG10(eps(elnum(jj))),jj=1,el-1),
     &      (LOG10(eps(elnum(jj))),jj=el+1,NELM),
     &       LOG10(MAX(1.Q-300, rhod/rho1)),
     &       LOG10(MAX(1.Q-300, dustV))
        write(71,2011) zz/km,rho1,p1,T1,nHges,mu1/amu,g1/100.0,Hp/km
        
        !-----------------------------------------------------------------
        ! ***  solve chemistry, hydrostatic equilibrium and T~p^kappa  ***
        ! ***  dT/dz = dT/dp dp/dz = -kappa T/p rho g = -kappa mu/k g. ***
        ! ***  eps0(:) and hence muH = rho/n<H> do not change.         ***
        ! ***  During the iteration, muH is with respect to total      ***
        ! ***  rho = rho_gas+rho_dust, so muH = rho/n<H>, but mu is    ***
        ! ***  with respect to rho_gas, so mu = rho_gas/p kT           *** 
        ! ***  => rho = rho_gas*(1+dg) = p*mu/kT*(1+dg) = p/kT * xx    ***
        ! ***  The iteration is done in xx                             *** 
        !-----------------------------------------------------------------
        dz  = LOG(pmax/pmin)*Hp/Npoints                 ! step in height
        g2  = grav*Mpl/(Rpl+zz+dz)**2                   ! gravity there
        mu2 = mu1                                       ! initial guess
        dg  = 0.0                                       ! dust/gas ratio
        xx  = (1.0+dg)*mu2
        ff  = p2
        do it=1,99
          T2    = T1 - kappa/bk*0.5*(mu2*g2+mu1*g1)*dz
          p2    = p1*(T2/T1)**(1.0/kappa)               ! target pressure
          kT    = bk*T2
          rho2  = p2/kT*xx                              ! total mass density
          nHges = rho2/muH                              ! n<H>
          if (model_eqcond) then
            call EQUIL_COND(nHges,T2,eps,Sat,eldust,verbose)
          else
            eps(:) = eps0(:)
          endif  
          call GGCHEM(nHges,T2,eps,.false.,0)
          nges = nel
          rhog = nel*mel
          do j=1,NELEM
            nges = nges + nat(j)
            rhog = rhog + nat(j)*mass(j)
          enddo
          do j=1,NMOLE
            nges = nges + nmol(j)
            rhog = rhog + nmol(j)*mmol(j)
          enddo
          pgas  = nges*kT
          ff    = p2-pgas                               ! remaining error
          dgold = dg
          dg    = rho2/rhog-1.0
          if (it==1) then
            muold = mu2
            xold  = xx
            mu2   = rhog*kT/pgas
            xx    = (1.0+dg)*mu2
            dx    = xx-xold
            dmu   = mu2-muold
            xold  = xx
            fold  = ff
            muold = mu2
          else
            mu2   = rhog*kT/pgas
            xx    = (1.0+dg)*mu2
            dfdx  = (ff-fold)/(xx-xold)
            dx    = -ff/dfdx
            dfdmu = (ff-fold)/(mu2-muold)
            dmu   = -ff/dfdmu
            xold  = xx
            fold  = ff
            muold = mu2
            if (ABS(dx/xold)<0.5.and.ABS(dmu/mu2)<0.5) then
              xx   = xx+dx
              !mu2 = xx/(1.0+dg)
              mu2  = mu2+dmu 
            else
              xx   = xx *EXP(dx/xx)
              mu2  = mu2*EXP(dmu/mu2)
            endif
          endif
          !write(98,'(I3,99(2(1pE14.7),2x))')
     >    !        it,xx/amu,dx/xx,mu2/amu,dmu/mu2,dg-dgold
          print '("p-it=",i3,"  mu=",2(1pE20.12))',it,mu2/amu,dx/xx
          if (ABS(dmu/mu2)<1.E-7.and.ABS(dx/xx)<1.E-7) exit
        enddo  
        if (it.ge.99) then
          stop "mu-iteration in auto_structure.f not converging."
        endif  
        !print*,p2,pgas
        !print*,bk/kappa*(T2-T1),-0.5*(mu1*g1+mu2*g2)*dz
        !print*,p1**(1-gamma)*T1**gamma,p2**(1-gamma)*T2**gamma
        !print*,LOG(p2/p1)/dz,-0.5/bk*(mu1*g1/T1+mu2*g2/T2)
        
        !--- make step ---
        zz   = zz+dz
        p1   = p2
        T1   = T2
        rho1 = rhog
        mu1  = mu2
        g1   = g2

        !--- remove all condensates and put them into the reservoir? ---
        if (remove_condensates) then
          fac = 1.Q+0
          do j=1,NDUST
            d_reservoir(j) = d_reservoir(j) + fac*eldust(j)
            do jj=1,dust_nel(j)
              k = dust_el(j,jj)
              e_reservoir(k) = e_reservoir(k) 
     &                       + fac*dust_nu(j,jj)*eldust(j)
            enddo
          enddo  
          if (verbose>0) then
            print*,"fractions of elements still in the gas phase ..."
            do j=1,NELM
              if (j==el) cycle 
              k = elnum(j)
              print'(A3,2(1pE18.10))',elnam(k),eps(k)/eps00(k),
     &                      (eps(k)+e_reservoir(k))/eps00(k)
            enddo
          endif  
          eps0(:) = eps(:) + (1.Q0-fac)*e_reservoir(:) ! remove elements
          !--- choice of output ---
          if (outAllHistory) then                      ! output will contain:
            eldust(:) = d_reservoir(:)                 ! all condensates ever
          else 
            eldust(:) = eldust(:)                      ! only local condensates
          endif
          muH = 0.d0
          do j=1,NELM
            if (j==el) cycle
            k = elnum(j)
            muH = muH + mass(k)*eps0(k)                ! re-compute muH
          enddo
          nHges = rho1/muH
        endif  

      enddo
      close(70)
      close(71)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,2(1pE20.6),9999(0pF20.7))
 2011 format(9999(1x,1pE19.10))
 3000 format(A8,0pF13.9)
      end
