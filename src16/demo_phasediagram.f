***********************************************************************
      SUBROUTINE DEMO_PHASEDIAGRAM
***********************************************************************
      use PARAMETERS,ONLY: Tmin,Tmax,pmin,pmax,nHmin,nHmax,Tfast,
     >                     model_eqcond,model_pconst,Npoints
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,muH,mass,
     >                    amu,dust_nam,dust_mass,dust_Vol
      use EXCHANGE,ONLY: nel,nat,nion,nmol,H,C,N,O,W
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real :: p,Tg,nHges,nges,kT,pgas,mu,muold,fac
      real :: ff,fold,dmu,dfdmu,xx,yy,epsH,epsO,epsC,epsN
      real :: rhog,rhod,Jstar,Nstar,xmin,xmax,ymin,ymax
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST),out(NDUST)
      integer :: i,it,ix,iy,j,jj,NOUT,stindex,iC,iW,iH2O,iH2Ol
      character(len=20) :: name,short_name(NDUST)
      integer :: verbose=0
      logical :: hasW,has_basecond

      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      iC = 0
      iH2O = 0
      iH2Ol = 0
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
        if (name=='C[s]') iC=i
        if (name=='H2O[s]') iH2O=i
        if (name=='H2O[l]') iH2Ol=i
      enddo
      has_basecond = (iC>0.and.iH2O>0.and.iH2Ol>0)
      iW = stindex(dust_nam,NDUST,'W[s]')
      hasW = (iW>0)
      eps = eps0
      NOUT = NELM
      if (charge) NOUT=NOUT-1
      open(unit=70,file='Static_Conc_2D.dat',status='replace')
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
     &               'dust/gas','Jstar(W)','Nstar(W)'

      !-------------------------------------
      ! ***  run chemistry on structure  ***
      !-------------------------------------
      mu = muH
      do ix=1,Npoints
        do iy=1,Npoints 
          fac = REAL(ix-1)/REAL(Npoints-1) 
          if (model_pconst) then
            p = EXP(LOG(pmax)+fac*LOG(pmin/pmax))
          else  
            nHges = EXP(LOG(nHmax)+fac*LOG(nHmin/nHmax))
          endif  
          Tg = EXP(LOG(Tmax)+LOG(Tmin/Tmax)*REAL(iy-1)/REAL(Npoints-1))
          if (Tmax==Tmin.and.pmax==pmin) then
            !--- interprete as variation of epsC and epsO ---
            p  = pmax
            Tg = Tmax
            xmin = 1.E-6
            xmax = 0.4
            ymin =-1.0+1.E-6
            ymax = 1.0-1.E-6
            xx = xmin + (xmax-xmin)*REAL(ix-1)/REAL(Npoints-1)
            yy = ymin + (ymax-ymin)*REAL(iy-1)/REAL(Npoints-1)
            epsC = xx                              ! x = C
            epsO = 0.5*((1.0-xx)*(yy+1.0))         ! y = (O-H)/(O+H)
            epsH = 0.5*((xx-1.0)*(yy-1.0))         ! 1 = H+O+C
            epsN = 1.E-3
            eps0(O) = epsO/epsH
            eps0(C) = epsC/epsH
            eps0(N) = epsN/epsH
            !print*,epsH,epsC,epsO
            !print*,yy,(epsO-epsH)/(epsO+epsH)
            !print*,xx,epsC/(epsH+epsO+epsC)
            !print*,1.0,(epsH+epsO+epsC)
            !print*,xx,eps0(C)/(eps0(H)+eps0(O)+eps0(C))
            !print*,yy,(eps0(O)-eps0(H))/(eps0(O)+eps0(H))
            !print*,eps0(H),eps0(C),eps0(N),eps0(O)
            !stop
            write(60,'(2(i4),99(1pE13.6))') ix,iy,
     &            eps0(H),eps0(C),eps0(N),eps0(O)
            muH = 0.0
            do i=1,NELM
              j = elnum(i)
              if (j==el) cycle
              muH = muH + eps0(j)*mass(j)
            enddo  
          endif  
          mu = muH     ! only useful to make it reproducible
          eldust = 0.0
          eps = eps0
          !--- iterate to achieve requested pressure ---
          do it=1,99
            if (model_pconst) nHges = p*mu/(bk*Tg)/muH
            if (model_eqcond) then
              call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
            endif  
            call GGCHEM(nHges,Tg,eps,.false.,0)
            kT = bk*Tg
            nges = nel
            do j=1,NELEM
              nges = nges + nat(j)
            enddo
            do j=1,NMOLE
              nges = nges + nmol(j)
            enddo
            pgas  = nges*bk*Tg
            ff    = p-pgas
            if (it==1) then
              muold = mu
              mu = nHges/pgas*(bk*Tg)*muH
              dmu = mu-muold
              if (.not.model_pconst) exit
            else
              dfdmu = (ff-fold)/(mu-muold)
              dmu   = -ff/dfdmu
              !write(98,'(I3,99(1pE14.7))')
     >        !      it,muold,mu,fold,ff,dfdmu,dmu/mu
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
          
          !--- compute supersat ratios and nucleation rates ---
          call SUPERSAT(Tg,nat,nmol,Sat)
          if (hasW) then
            call NUCLEATION('W',Tg,dust_vol(iW),nat(W),
     &                      Sat(iW),Jstar,Nstar)
          else
            Jstar = 0
            Nstar = 9.e+99
          endif
          !--- compute dust/gas density ratio ---
          rhog = nHges*muH
          rhod = 0.0
          do jj=1,NDUST
            rhod = rhod + nHges*eldust(jj)*dust_mass(jj)
            out(jj) = LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj))))
            if (ABS(Sat(jj)-1.Q0)<1.E-10) out(jj)=0.Q0
          enddo  

          print'(i4,i4," Tg[K] =",0pF8.2,"  n<H>[cm-3] =",1pE10.3)',
     >          ix,iy,Tg,nHges
          write(*,1010) ' Tg=',Tg,' n<H>=',nHges,
     &                  ' p=',pgas/bar,' mu=',mu/amu,
     &                  ' dust/gas=',rhod/rhog
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
     &       LOG10(MAX(1.Q-300, Jstar)), 
     &       MIN(999999.99999,Nstar)

          if (has_basecond) then
            if (model_eqcond) then
              write(60,'(2(i4),99(1pE13.6))') ix,iy,
     &               eps(H),eps(C),eps(N),eps(O),
     &               eldust(iC),eldust(iH2O),eldust(iH2Ol)
            else  
              write(60,'(2(i4),99(1pE13.6))') ix,iy,
     &               eps(H),eps(C),eps(N),eps(O),
     &               Sat(iC),Sat(iH2O),Sat(iH2Ol)
            endif  
          endif
          
        enddo
      enddo  
      close(70)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(a4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,2(1pE20.6),9999(0pF20.7))
      end  
