***********************************************************************
      SUBROUTINE CALC_OPAC(nHges,eldust,T)
***********************************************************************
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,
     >                    muH,mass,mel,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      use OPACITY,ONLY: NLAM,lam,NLIST,opind,duind,nn,kk,
     >                  NSIZE,aa,ff,aweight,kabs,ksca,kext
      use DATATYPE,ONLY: r2            ! from MIEX
      use MIE_ROUTINES,ONLY: SHEXQNN2  ! from MIEX
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: mic=1.E-4
      integer,parameter :: qp=selected_real_kind(33,4931)
      real,intent(in) :: nHges,T
      real(kind=qp),intent(in) :: eldust(NDUST)
      real :: rhod1,rhod2,Vcon1,Vcon2,rhogr1,rhogr2,rho
      real :: neff,keff,Vs(NDUST),porosity
      integer :: i,j
      !------------ variables for data exchange with MIEX ---------
      complex(kind=r2) :: ri
      real(kind=r2)    :: xx,Qext,Qsca,Qabs,Qbk,Qpr,albedo,g,mm1
      integer          :: ier,nang
      complex(kind=r2),dimension(2) :: SA1,SA2
      
      rhod1 = 0.d0               ! dust mass density   [g/cm3]  
      Vcon1 = 0.d0               ! dust volume density [cm3/cm3]
      rhod2 = 0.d0               
      Vcon2 = 0.d0               ! same, but only opacity species
      do i=1,NDUST
        if (eldust(i)<=0.Q0) cycle
        j = duind(i)
        rhod1 = rhod1 + nHges*eldust(i)*dust_mass(i)
        Vcon1 = Vcon1 + nHges*eldust(i)*dust_Vol(i)
        if (j==0) cycle
        rhod2 = rhod2 + nHges*eldust(i)*dust_mass(i)
        Vcon2 = Vcon2 + nHges*eldust(i)*dust_Vol(i)
      enddo
      porosity = 0.25            ! add 25% porosity
      Vcon1 = Vcon1*(1.0+porosity)
      Vcon2 = Vcon2*(1.0+porosity)
      rhogr1 = rhod1/Vcon1       
      rhogr2 = rhod2/Vcon2
      Vs(1:NLIST) = 0.0
      do i=1,NDUST
        if (eldust(i)<=0.Q0) cycle
        j = duind(i)
        if (j==0) then
          print'(" ***",A16,1pE11.3,"  not an opacity species ***")',
     >         trim(dust_nam(i)),nHges*eldust(i)*dust_Vol(i)/Vcon1
        else
          Vs(j) = nHges*eldust(i)*dust_Vol(i)/Vcon2
          if (eldust(i)<=0.Q0) cycle
          print'(A20,1pE11.3)',trim(dust_nam(i)),(1.0-porosity)*Vs(j)
        endif
      enddo
      Vs = (1.0-porosity)*Vs
      Vs(NLIST) = porosity
      print'(A20,1pE11.3)',"vacuum",porosity
      rho = nHges*muH
      print*,"dust material density [g/cm3]",rhogr1,rhogr2
      print*,"dust/gas mass ratio   [g/cm3]",rhod1/rho,rhod2/rho
      print*,"dust volume density [cm3/cm3]",Vcon1,Vcon2

      !-------------------------------
      ! ***  dust size dist.model  ***
      !-------------------------------
      call SIZE_DIST(rho,rhogr2,Vcon2)

      !---------------------------------------------
      ! ***  effective medium and Mie opacities  ***
      !---------------------------------------------
      kabs(1:NLAM) = 0.0
      ksca(1:NLAM) = 0.0
      do j=1,NLAM
        call effMedium(j,Vs,neff,keff)
        print*,lam(j),neff,keff
        nang = 3
        ri = DCMPLX(neff,keff)
        do i=1,NSIZE
          xx = 2.0*pi*aa(i)/(lam(j)*mic)     ! Mie size parameter
          call SHEXQNN2(ri,xx,Qext,Qsca,Qabs,Qbk,Qpr,
     >                  albedo,g,ier,SA1,SA2,.false.,nang)
          kabs(j) = kabs(j) + ff(i)*pi*aa(i)**2 * Qabs * aweight(i)
          ksca(j) = ksca(j) + ff(i)*pi*aa(i)**2 * Qsca * aweight(i)
        enddo  
      enddo
      kext(1:NLAM) = kabs(1:NLAM)+ksca(1:NLAM)

      do j=1,NLAM
        print'(0pF8.2,3(1pE12.4))',lam(j),
     >       kabs(j)/rhod2,ksca(j)/rhod2,kext(j)/rhod2
      enddo
      
      end

***********************************************************************
      SUBROUTINE SIZE_DIST(rho,rhogr,Vcon)
***********************************************************************
      use OPACITY,ONLY: NSIZE,aa,ff,aweight
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: mic=1.E-4
      real,intent(in) :: rho    ! gas mass density [g/cm3]
      real,intent(in) :: rhogr  ! dust material density [g/cm3]
      real,intent(in) :: Vcon   ! dust volume/H-nucleus [cm3]
      real :: rhoref,a1ref,a2ref,pp,VV,mm,Vref,mref,ndref,dg,scale
      real :: ndtest,Vtest,mtest,da
      integer :: i

      !------------------------------------------------
      ! ***  set up reference dust size dist.model  ***
      !------------------------------------------------
      NSIZE  = 100
      rhoref = 2.0              ! g/cm3
      a1ref  = 0.05             ! mic
      a2ref  = 1000.0           ! mic
      pp     = -3.5
      do i=1,NSIZE
        aa(i) = EXP(LOG(a1ref)+(i-1.0)/(NSIZE-1.0)*LOG(a2ref/a1ref))*mic
        ff(i) = aa(i)**pp
      enddo
      aweight = 0.0
      do i=2,NSIZE
        da = 0.5*(aa(i)-aa(i-1))
        aweight(i)   = aweight(i)   + da
        aweight(i-1) = aweight(i-1) + da
      enddo
      ndref = 0.0
      Vref  = 0.0
      mref  = 0.0
      do i=1,NSIZE
        VV = 4.0*pi/3.0*aa(i)**3
        mm = VV*rhoref
        ndref = ndref + ff(i)*aweight(i)
        Vref  = Vref  + ff(i)*VV*aweight(i)
        mref  = mref  + ff(i)*mm*aweight(i)
      enddo
      dg    = mref/rho             ! dust/gas mass ratio ...
      scale = 0.004/dg             ! ... which should be 0.004
      ff    = scale*ff             ! [cm-4]
      ndref = scale*ndref          ! [cm-3]
      Vref  = scale*Vref           ! [cm3/H-nucleus]

      !------------------------------------------------------
      ! ***  adjust dust sizes and volume while nd=const  ***
      ! ***  we assume that each grain looses or gains    ***
      ! ***  the same volume fraction.                    ***               
      !------------------------------------------------------
      scale   = (Vcon/Vref)**(1.0/3.0)
      aa      = scale*aa
      ff      = ff/scale
      aweight = scale*aweight
      !--- done, this is just a check ---
      ndtest  = 0.0
      Vtest   = 0.0
      mtest   = 0.0
      do i=1,NSIZE
        VV = 4.0*pi/3.0*aa(i)**3
        mm = VV*rhogr
        ndtest = ndtest + ff(i)*aweight(i)
        Vtest  = Vtest  + ff(i)*VV*aweight(i)
        mtest  = mtest  + ff(i)*mm*aweight(i)
      enddo
      dg = mtest/rho
      print*,"scale=",scale
      print*,"   nd=",ndref,ndtest
      print*,"  d/g=",dg
      print*,"Vdust=",Vref,Vtest,Vcon

      end
      
      
