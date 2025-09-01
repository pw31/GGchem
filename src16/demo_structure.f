***********************************************************************
      SUBROUTINE DEMO_STRUCTURE
***********************************************************************
      use PARAMETERS,ONLY: Mpl,Rpl,Tmin,Tmax,pmin,pmax,nHmin,nHmax,
     >                     model_eqcond,model_pconst,Npoints,
     >                     model_struc,struc_file,remove_condensates,
     >                     model_eqcond,Nseq,Tseq
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,
     >                    muH,mass,mel,amu,grav,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol,H,C,N,O,W
      use STRUCTURE,ONLY: Npmax,Tgas,press,pelec,dens,nHtot,estruc,zz
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),eps00(NELEM)
      real(kind=qp) :: Sat(NDUST),eldust(NDUST),out(NDUST)
      real(kind=qp) :: fac,e_reservoir(NELEM),d_reservoir(NDUST)
      real :: dat(1000),ddust
      real :: tau,p,pe,rho,nHges,nges,kT,pgas,mu,gg,Hp,mugas
      real :: muold,dmu,ff,fold,dfdmu,km=1.D+5,AU=1.4959787d+13 
      real :: Jstar,Nstar,rhog,dustV,rhod,L3,bmix,emono,rdum,zdum
      real :: Tg,Td,cT,rcut,Kzz
      integer :: i,j,k,l,e,jj,dk,NOUT,Nfirst,Nlast,Ninc,iW,idum
      integer :: it,n1,n2,n3,n4,n5,Ndat,dind(1000),ek,eind(1000)
      integer :: Nx,Nz,ix,iz,Nfound,e_source(100),e_target(100)
      integer :: stindex,cut,iseq,verbose=0
      character(len=20000) :: header
      character(len=200) :: line,filename
      character(len=20) :: name,short_name(NDUST),dname,ename
      character(len=1) :: char
      logical :: hasW,efound,take
      logical :: outAllHistory=.false.

      !-----------------------------
      ! ***  read the structure  ***
      !-----------------------------
      if ((index(struc_file,'~')==1).or.
     &    (index(struc_file,'/')==1)) then 
        filename = trim(struc_file)
      else
        filename = 'structures/'//trim(struc_file)
      endif  
      print*,"reading "//trim(filename)//" ..."

      !--------------------------------------------------------
      if (model_struc==1) then
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        Npoints = 256 
        do i=1,5+Npoints+2 
          read(3,'(A200)') line
        enddo  
        do i=1,Npoints 
          read(3,*) iz,tau,Tg,p,pe,rho
          print*,iz
          Tgas(i)  = Tg
          press(i) = p
          pelec(i) = pe
          dens(i)  = rho
          nHtot(i) = rho/muH
          estruc(i,:) = eps0(:)
        enddo
        close(3)
        Nfirst = Npoints
        Nlast  = 2
        Ninc   = -1  ! botton to top

      !--------------------------------------------------------
      else if (model_struc==2) then    ! tp_w18.dat
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        Npoints = 48
        do i=1,2
          read(3,'(A200)') line
        enddo  
        do i=1,Npoints 
          read(3,*) Tg,p
          Tgas(i)  = Tg
          press(i) = p*bar
          estruc(i,:) = eps0(:)
        enddo
        close(3)
        model_pconst = .true.
        Nfirst = Npoints
        Nlast  = 2
        Ninc   = -1            ! bottom to top

      !--------------------------------------------------------
      else if (model_struc==3) then    ! weather_*.dat
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        read(3,'(A200)') line
        read(3,*) n1,n2,n3,n4,n5,Npoints  ! NELM,NMOLE,NDUST,NEPS,NNUC,Npoints
        Ndat = 6 + 2*n1 + n2 + 2*n3
        read(3,'(A20000)') header
        do j=Ndat-n1-n3+1,Ndat-n1
          dname = adjustl(header(j*20-19:j*20))
          dname = trim(dname(2:))//"[s]"
          dind(j) = 0
          do dk=1,NDUST
            if (trim(dust_nam(dk))==trim(dname)) then
              print*,dk,trim(dust_nam(dk))//" = "//trim(dname) 
              dind(j) = dk
            endif  
          enddo
          if (dind(j)==0) then
            print*,"*** dust kind "//trim(dname)//" not found."
            stop
          endif  
        enddo  
        Nfound = 0
        do k=1,NELM-1
          e = elnum(k)
          efound = .false.
          do j=6+n1+n2+2*n3+1,6+n1+n2+2*n3+n1
            ename = adjustl(header(j*20-19:j*20))
            !print*,e,"eps"//trim(elnam(e))//" = "//trim(ename)
            if (trim(ename)=="eps"//trim(elnam(e))) then
              Nfound = Nfound+1
              efound = .true.
              e_source(Nfound) = j
              e_target(Nfound) = e
              exit
            endif  
          enddo
          if (efound) then
            print*,"found "//trim(elnam(e))
          else
            print*,"not found "//trim(elnam(e))
          endif
        enddo
        do i=1,Npoints
          print*,i,Npoints
          read(3,*) dat(1:Ndat) 
          Tgas(i)  = dat(1)
          nHtot(i) = dat(2)
          press(i) = dat(3)
          dens(i)  = dat(5)
          pelec(i) = 10.d0**dat(7)*bk*Tgas(i)
          estruc(i,:) = eps0(:)          
          do k=1,Nfound
            j = e_source(k)
            e = e_target(k)
            estruc(i,e) = 10.Q0**dat(j)
          enddo   
          if (model_eqcond) then
            do j=Ndat-n1-n3+1,Ndat-n1
              ddust = 10.Q0**dat(j)
              dk = dind(j)
              do k=1,dust_nel(dk)
                e = dust_el(dk,k)
                estruc(i,e) = estruc(i,e) + ddust*dust_nu(dk,k)    
              enddo
            enddo  
          endif  
          !do k=1,NELM-1
          !  e = elnum(k) 
          !  print'(I3,A3,2(1pE18.10))',i,elnam(e),eps0(e),estruc(i,e) 
          !enddo  
        enddo
        close(3)
        Nfirst = 1
        Nlast  = Npoints
        Ninc   = 1             ! bottom to top

      !--------------------------------------------------------
      else if (model_struc==4) then
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        do i=1,1
          read(3,'(A200)') line
        enddo  
        do i=1,9999
          read(3,*,end=444) idum,p,Tg
          Tgas(i)  = Tg
          press(i) = p*bar
          estruc(i,:) = eps0(:)
        enddo
 444    continue
        close(3)
        Npoints = i-1
        model_pconst = .true.
        Nfirst = Npoints
        Nlast  = 1
        Ninc   = -1            ! bottom to top

      !--------------------------------------------------------
      else if (model_struc==5) then   ! Jup-T-p.dat
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        do i=1,1
          read(3,'(A200)') line
        enddo  
        do i=1,9999
          read(3,*,end=555) p,Tg
          Tgas(i)  = Tg
          press(i) = p*bar
          estruc(i,:) = eps0(:)
        enddo
 555    continue
        close(3)
        Npoints = i-1
        model_pconst = .true.
        Nfirst = 1
        Nlast  = Npoints
        Ninc   = 1             ! bottom to top

      !--------------------------------------------------------
      else if (model_struc==6) then   ! StaticWeather output
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        read(3,'(A200)') line
        read(3,*) n1,n2
        read(3,'(A200)') line
        read(3,'(A200)') line
        do i=1,999999
          read(3,*,end=666) dat(1:5)
          zz(i)    = dat(1)
          Tgas(i)  = dat(2)
          nHtot(i) = dat(3)
          dens(i)  = dat(4)
          press(i) = dat(5)
          muH = dens(i)/nHtot(i)
          !print*,i,Tgas(i),nHtot(i)
        enddo
 666    Npoints=i-1
        close(3)
        filename = 'structures/out3_dust.dat'
        open(3,file=filename,status='old')
        read(3,'(A200)') line
        read(3,*) n1,n2
        read(3,'(A200)') line
        read(3,'(A20000)') header
        do k=1,n1
          j = 8+k + 2*n2 + n1
          ename = adjustl(header(j*14-13:j*14))
          ename = trim(ename(5:))
          do ek=1,NELM
            e = elnum(ek)
            if (trim(ename)==trim(elnam(e))) then
              print*,e,trim(ename)//" found"
              eind(j) = e
            endif
          enddo
        enddo
        do k=1,n2
          j = 8+k
          dname = adjustl(header(j*14-13:j*14))
          dname = trim(dname(3:))
          dind(j) = 0
          do dk=1,NDUST
            if (trim(dust_nam(dk))==trim(dname)) then
              print*,dk,trim(dust_nam(dk))//" = "//trim(dname) 
              dind(j) = dk
            endif  
          enddo
        enddo
        Ndat = 8 + 2*n1 + 2*n2
        do i=1,Npoints
          read(3,*) dat(1:Ndat) 
          rho = dens(i)                       ! gas density g/cm3
          L3 = dat(6)                         ! dust volume cm3/g 
          estruc(i,:) = eps0(:) 
          do k=1,n1
            j = 8+k + 2*n2 + n1
            e = eind(j) 
            estruc(i,e) = dat(j)
            !print*,elnam(e),dat(j)
          enddo   
          !print'(99(A12))',elnam(1:NELEM)
          !print'(99(1pE12.3))',estruc(i,:)/eps0(:)
          do k=1,n2
            j = 8+k
            dk = dind(j) 
            bmix = dat(j)                     ! volume ratio
            emono = bmix*rho*L3/dust_vol(dk)/nHtot(i)
            do l=1,dust_nel(dk)
              e = dust_el(dk,l)
              estruc(i,e) = estruc(i,e) + emono*dust_nu(dk,l)    
            enddo
          enddo  
          print'(99(A12))',(elnam(elnum(j)),j=1,el-1),
     &                     (elnam(elnum(j)),j=el+1,NELM)
          print'(99(1pE12.3))',(eps0(elnum(j)),j=1,el-1),
     &                         (eps0(elnum(j)),j=el+1,NELM)
        enddo  
        close(3)
        Nfirst = Npoints
        Nlast  = 1
        Ninc   = -1            ! top to bottom 

      !--------------------------------------------------------
      else if (model_struc==7) then     ! VenusHighFit.dat
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        do i=1,1
          read(3,'(A200)') line
        enddo
        read(3,*) Npoints
        do i=1,Npoints
          read(3,*) idum,zdum,p,Tg
          zz(i) = zdum*km
          Tgas(i)  = Tg
          press(i) = p*bar
          estruc(i,:) = eps0(:)
          if (p*bar<pmin) exit
        enddo
        close(3)
        model_pconst = .true.
        Nfirst = 1
        Nlast  = MIN(Npoints,i)
        Ninc   = 1             ! bottom to top

      !--------------------------------------------------------
      else if (model_struc==8) then     ! ProDiMo.dat
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        do i=1,999
          read(3,'(A200)') line
          print*,trim(line)
          if (index(line,"Nx,Nz,Nsp")>0) exit  
        enddo
        i = index(line,"=")
        read(line(i+1:),*) Nx,Nz
        read(3,'(A200)') line
        read(3,'(A200)') line
        print*
        print'("Would you like a radial (1) or vertical (2) cut? ",$)'
        read*,cut
        if (cut==1) then
          !--------- midplane cut --------
          Npoints = 0
          do iz=Nz,1,-1
            take = .true.
            do ix=1,Nx
              read(3,*) i,j,rdum,zdum,dat(1:5),Tg,Td,cT,rho,pgas
              take = take.and.(Td>Tmin)
              if (iz==1.and.take) then
                !print*,i,rdum,Td
                Npoints = Npoints+1
                zz(i) = rdum*AU
                Tgas(i) = Td
                dens(i) = rho
                press(i) = pgas
                nHtot(i) = rho/muH
                estruc(i,:) = eps0(:)
              endif
            enddo
          enddo
          Nfirst = 1
          Nlast  = Npoints
          Ninc   = 1            ! high-to-low T
        else
          !--------- vertical cut --------
          print'("vertical cut at which radius[AU]? ",$)'
          read*,rcut
          do iz=Nz,1,-1
            do ix=1,Nx
              read(3,*) i,j,rdum,zdum,dat(1:5),Tg,Td,cT,rho,pgas
              if (rdum>rcut) cycle
              !print*,i,j,rdum,Td
              zz(Nz+1-iz) = zdum*AU
              Tgas(Nz+1-iz) = Td
              dens(Nz+1-iz) = rho
              press(Nz+1-iz) = pgas
              nHtot(Nz+1-iz) = rho/muH
              estruc(Nz+1-iz,:) = eps0(:)
            enddo
          enddo
          Nfirst  = Nz
          Npoints = 0
          do i=1,Nz
            if (nHtot(i)<1.E+4) cycle 
            print*,i,Tgas(i),nHtot(i)
            Nfirst = MIN(Nfirst,i)
          enddo
          Npoints = Nlast-Nfirst+1
          Nlast  = Nz
          Ninc   = 1            ! high-to-low T
        endif
        
      !--------------------------------------------------------
      else if (model_struc==9) then   ! ARCiS output
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        do i=1,1
          read(3,'(A200)') line
        enddo  
        do i=1,99999
          read(3,*,end=556) rdum,zdum,rho,dat(1),Tg,pgas
          zz(i) = zdum
          Tgas(i)  = Tg
          press(i) = pgas*bar
          dens(i) = rho
          nHtot(i) = rho/muH
          estruc(i,:) = eps0(:)
        enddo
 556    continue
        close(3)
        Npoints = i-1
        Nfirst = 1
        Nlast  = Npoints
        Ninc   = 1             ! bottom to top

      !--------------------------------------------------------
      else if (model_struc==10) then     ! Paul's pT-*.dat
      !--------------------------------------------------------
        open(3,file=filename,status='old')
        do i=1,2
          read(3,'(A200)') line
        enddo
        Npoints = 0
        do i=1,99999
          read(3,*,end=557) zdum,Tg,Kzz,rho,pgas
          if (pgas*bar<pmin) cycle
          Npoints = Npoints+1
          zz(Npoints) = zdum*km
          Tgas(Npoints)  = Tg
          press(Npoints) = pgas*bar
          dens(Npoints) = rho
          estruc(Npoints,:) = eps0(:)
        enddo
 557    close(3)
        Nfirst = Npoints
        Nlast  = 1
        Ninc   = -1             ! data given top to bottom
        
      else
        print*,"*** unknown file format =",model_struc
        stop
      endif  

      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
      enddo
      iW = stindex(dust_nam,NDUST,'W[s]')
      hasW = (iW>0)
      eps  = eps0
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
     &               'dust/gas','dustVol/H','Jstar(W)','Nstar(W)'
      open(unit=71,file='Structure.dat',status='replace')
      write(71,2000) 'z[km]','rho[g/cm3]','pgas[dyn/cm2]','T[K]',
     &               'n<H>[cm-3]','mu[amu]','g[m/s2]','Hp[km]'

      !-------------------------------------
      ! ***  run chemistry on structure  ***
      !-------------------------------------
      e_reservoir = 0.Q0
      d_reservoir = 0.Q0
      eps00 = eps0
      eldust = 0.Q0
      mu = muH
      do i=Nfirst,Nlast,Ninc
        Tg      = Tgas(i)
        p       = press(i) 
        nHges   = nHtot(i)
        eps0(:) = estruc(i,:)

        print*
        print'("new point",I4,"  n<H>=",1pE12.5," cm-3  T=",0pF8.3,
     &         " K  p=",1pE12.5," bar")',i,nHges,Tg,p/bar
        print'(99(A12))',(elnam(elnum(j)),j=1,el-1),
     &                   (elnam(elnum(j)),j=el+1,NELM)
        print'(99(1pE12.3))',(eps0(elnum(j)),j=1,el-1),
     &                       (eps0(elnum(j)),j=el+1,NELM)

        if (i==Nfirst.and.Nseq>1) then  ! on base point, use Tseq() to get solution
          iseq = 1
          Tg = Tseq(1)
        endif
        !--- run chemistry (+phase equilibrium)    ---
        !--- iterate to achieve requested pressure ---
 100    continue
        do it=1,99
          if (model_pconst) nHges = p*mu/(bk*Tg)/muH
          if (model_eqcond) then
            call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
          else
            eps(:) = eps0(:)
          endif  
          call GGCHEM(nHges,Tg,eps,.false.,0)
          kT = bk*Tg
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
          pgas = nges*kT
          ff   = p-pgas
          if (it==1) then
            muold = mu
            mu = nHges/pgas*(bk*Tg)*muH
            dmu = mu-muold
            if (.not.model_pconst) exit
          else
            dfdmu = (ff-fold)/(mu-muold)
            dmu   = -ff/dfdmu
            !write(98,'(I3,99(1pE14.7))')
     >      !     it,muold,mu,fold,ff,dfdmu,dmu/mu
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
        if (i==Nfirst.and.Nseq>1) then  ! on base point, use Tseq() to get solution
          if (iseq<Nseq) then
            iseq = iseq+1
            Tg = Tseq(iseq)
            if (iseq==Nseq) Tg=Tgas(i)
            goto 100
          endif
        endif
          
        !--- remove all condensates and put them in reservoir? ---
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
          do j=1,NELM
            if (j==el) cycle 
            k = elnum(j)
            print'(A3,2(1pE18.10))',elnam(k),eps(k)/eps00(k),
     &                      (eps(k)+e_reservoir(k))/eps00(k)
          enddo
          eps0(:) = eps(:) + (1.Q0-fac)*e_reservoir(:)
          estruc(i+Ninc,:) =  eps0(:)                  ! for next layer
          if (outAllHistory) then                      ! output will contain:
            eldust(:) = d_reservoir(:)                 ! all condensates ever
          else 
            eldust(:) = eldust(:)                      ! only local condensates
          endif
        endif  

        !--- compute supersat ratios and nucleation rates ---
        call SUPERSAT(Tg,nat,nmol,Sat)
        if (hasW) then
          call NUCLEATION('W',Tg,dust_vol(iW),nat(W),
     &                    Sat(iW),Jstar,Nstar)
        else
          Jstar = 0
          Nstar = 9.e+99
        endif  

        !--- compute dust/gas density ratio ---
        rhod  = 0.0
        dustV = 0.0
        do jj=1,NDUST
          rhod  = rhod  + nHges*eldust(jj)*dust_mass(jj)
          dustV = dustV + eldust(jj)*dust_Vol(jj)
          out(jj) = LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj))))
          if (ABS(Sat(jj)-1.Q0)<1.E-10) out(jj)=0.Q0
        enddo  

        print'(i4," Tg[K] =",0pF8.2,"  n<H>[cm-3] =",1pE10.3)',
     >        i,Tg,nHges

        write(*,1010) ' Tg=',Tg,' n<H>=',nHges,
     &                ' p=',pgas/bar,' mu=',mu/amu,
     &                ' dust/gas=',rhod/rhog
        print*
        write(70,2010) Tg,nHges,pgas,
     &       LOG10(MAX(1.Q-300, nel)),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=1,el-1),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=el+1,NELM),
     &      (LOG10(MAX(1.Q-300, nmol(jj))),jj=1,NMOLE),
     &      (out(jj),jj=1,NDUST),
     &      (LOG10(MAX(1.Q-300, eldust(jj))),jj=1,NDUST),
     &      (LOG10(eps(elnum(jj))),jj=1,el-1),
     &      (LOG10(eps(elnum(jj))),jj=el+1,NELM),
     &       LOG10(MAX(1.Q-300, rhod/rhog)),
     &       LOG10(MAX(1.Q-300, dustV)),
     &       LOG10(MAX(1.Q-300, Jstar)), 
     &       MIN(999999.99999,Nstar)
        gg = grav*Mpl/(Rpl+zz(i))**2
        mugas = rhog/nges
        Hp = bk*Tg/(mugas*gg)
        write(71,2011) zz(i)/km,rhog,pgas,Tg,nHges,mugas/amu,gg,Hp

        if (verbose>0) read(*,'(a1)') char

      enddo  

      close(70)
      close(71)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,2(1pE20.6),9999(0pF20.7))
 2011 format(9999(1x,1pE19.10))
      end  

