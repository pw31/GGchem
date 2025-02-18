***********************************************************************
      SUBROUTINE CALC_OPAC(nHges,eldust,kabs,ksca,kext,dustgas,verb)
***********************************************************************
***                                                                 ***
***   calculates the dust opacities kabs,ksca,kext [1/cm] using     ***
***                                                                 ***
***   (1) optical constants for various solid materials (nn,kk)     ***
***       see README in data/OpticalData/README                     ***
***       note that many materials have no opacity data             ***
***   (2) effective mixing theory (neff,keff)                       ***
***   (3) a model for the size distibution function f(a)            ***
***   (4) Mie theory (routine MIEX.f90)                             ***
***                                                                 ***
***   input is n<H> [cm-3] and solid unit concentrations eldust(:)  ***
***                                                                 ***
***********************************************************************
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_mass,dust_Vol,muH
      use OPACITY,ONLY: NLAMmax,NLAM,lam,NLIST,opind,nn,kk,
     >                  NSIZE,aa,ff,aweight,
     >                  nmin,nmax,kmin,kmax,xmin,xmax
      use DATATYPE,ONLY: r2            ! from MIEX
      use MIE_ROUTINES,ONLY: SHEXQNN2  ! from MIEX
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: mic=1.E-4
      integer,parameter :: qp=selected_real_kind(33,4931)
      real,intent(in) :: nHges
      real(kind=qp),intent(in) :: eldust(NDUST)
      real,intent(out),dimension(NLAMmax) :: kabs,ksca,kext
      real,intent(out) :: dustgas
      integer,intent(in) :: verb
      real :: rhod1,rhod2,Vcon1,Vcon2,rhogr1,rhogr2,rho
      real :: neff,keff,porosity
      real,dimension(NDUST) :: qsave,qsort,Vs
      integer :: i,j,k
      logical,parameter :: use_fastmie=.true.
      !------------ variables for data exchange with MIEX ---------
      complex(kind=r2) :: ri
      real(kind=r2)    :: xx,Qext,Qsca,Qabs,Qbk,Qpr,albedo,g,mm1
      integer          :: ier
      complex(kind=r2),dimension(2) :: SA1,SA2
      
      kabs(1:NLAM) = 0.0         ! [1/cm]
      ksca(1:NLAM) = 0.0         ! [1/cm]
      kext(1:NLAM) = 0.0         ! [1/cm]
      rhod1 = 0.d0               ! dust mass density   [g/cm3]  
      Vcon1 = 0.d0               ! dust volume density [cm3/cm3]
      rhod2 = 0.d0               
      Vcon2 = 0.d0               ! same, but only opacity species
      dustgas = 0.d0
      do i=1,NDUST
        if (eldust(i)<=0.Q0) cycle
        j = opind(i)
        rhod1 = rhod1 + nHges*eldust(i)*dust_mass(i)
        Vcon1 = Vcon1 + nHges*eldust(i)*dust_Vol(i)
        if (j==0) cycle
        rhod2 = rhod2 + nHges*eldust(i)*dust_mass(i)
        Vcon2 = Vcon2 + nHges*eldust(i)*dust_Vol(i)
      enddo
      if (rhod1==0.d0) return    ! nothing to do
      porosity = 0.25            ! add 25% porosity
      Vcon1 = Vcon1/(1.0-porosity)
      Vcon2 = Vcon2/(1.0-porosity)
      rhogr1 = rhod1/Vcon1       
      rhogr2 = rhod2/Vcon2
      do i=1,NDUST
        qsort(i) = nHges*eldust(i)*dust_Vol(i)/Vcon2
      enddo
      qsave = qsort
      call sort(NDUST,qsort)
      Vs(1:NLIST) = 0.0
      do k=NDUST,1,-1
        if (qsort(k)<=0.0) cycle
        do i=1,NDUST
          if (qsort(k)==qsave(i)) then
            j = opind(i)
            if (j==0) then
              if (nHges*eldust(i)*dust_Vol(i)/Vcon1>3.E-3) then
                if (verb>=-1) print'(" ***",A16,1pE11.3,
     >              "  not an opacity species ***")',trim(dust_nam(i)),
     >              nHges*eldust(i)*dust_Vol(i)/Vcon1
              endif
            else 
              Vs(j) = qsave(i)
              if (verb>=-1) print'(A20,1pE11.3)',trim(dust_nam(i)),Vs(j)
            endif
          endif
        enddo
      enddo
      Vs(NLIST) = porosity
      rho = nHges*muH
      dustgas = rhod1/rho
      if (verb>=-1) then
        print'(A20,1pE11.3)',"vacuum",porosity
      endif  
      if (verb>=0) then
        print*,"dust material density [g/cm3]",rhogr1,rhogr2
        print*,"dust/gas mass ratio   [g/cm3]",rhod1/rho,rhod2/rho
        print*,"dust volume density [cm3/cm3]",Vcon1,Vcon2
      endif
      if (Vcon1==0.0) then
        kabs(1:NLAM) = 0.0
        ksca(1:NLAM) = 0.0
        kext(1:NLAM) = 0.0
        return
      endif
      
      !--------------------------------------
      ! ***  dust size distibution model  ***
      !--------------------------------------
      call SIZE_DIST(nHges,rhogr1,Vcon1,verb)

      !---------------------------------------------
      ! ***  effective medium and Mie opacities  ***
      !---------------------------------------------
      kext(1:NLAM) = 0.0                    ! [1/cm]
      kabs(1:NLAM) = 0.0                    ! [1/cm]
      ksca(1:NLAM) = 0.0                    ! [1/cm]
      do j=1,NLAM
        call effMedium(j,Vs,neff,keff)
        !print*,lam(j),neff,keff
        nmax = MAX(nmax,neff)
        nmin = MIN(nmin,neff)
        kmax = MAX(kmax,keff)
        kmin = MIN(kmin,keff)
        ri = DCMPLX(neff,keff)
        do i=1,NSIZE
          xx = 2.0*pi*aa(i)/(lam(j)*mic)    ! Mie size parameter
          xmin = MIN(xmin,xx)
          xmax = MAX(xmax,xx)
          if (use_fastmie) then
            call FASTMIE(xx,neff,keff,Qext,Qsca,Qabs,.false.)
          else
            call SHEXQNN2(ri,xx,Qext,Qsca,Qabs,Qbk,Qpr,
     >           albedo,g,ier,SA1,SA2,.false.,2)
          endif
          kext(j) = kext(j) + ff(i)*pi*aa(i)**2 * Qext * aweight(i)
          kabs(j) = kabs(j) + ff(i)*pi*aa(i)**2 * Qabs * aweight(i)
          ksca(j) = ksca(j) + ff(i)*pi*aa(i)**2 * Qsca * aweight(i)
        enddo  
      enddo
      end

***********************************************************************
      SUBROUTINE SIZE_DIST(nHges,rhogr,Vcon,verb)
***********************************************************************
      use DUST_DATA,ONLY: muH
      use OPACITY,ONLY: NSIZEmax,NSIZE,aa,ff,aweight
      implicit none
      real,parameter :: pi=ACOS(-1.0)
      real,parameter :: mic=1.E-4
      real,intent(in) :: nHges  ! H-nuclei particle density [cm-3]
      real,intent(in) :: rhogr  ! dust material density [g/cm3]
      real,intent(in) :: Vcon   ! dust volume/H-nucleus [cm3]
      integer,intent(in) :: verb
      real :: VV,mm,dg,scale,ndtest,Vtest,mtest,da,pp
      integer :: i
      logical,save :: firstCall=.true.
      real,save :: rhoref,a1ref,a2ref,Vref,mref,ndref,nHref,dgref
      real,dimension(NSIZEmax),save :: aref,fref,awref

      !------------------------------------------------
      ! ***  set up reference dust size dist.model  ***
      !------------------------------------------------
      if (firstCall) then
        NSIZE  = 100
        rhoref = 2.0              ! g/cm3
        a1ref  = 0.05             ! mic
        a2ref  = 3000.0           ! mic
        dgref  = 0.004            ! fully condensed solor dust/gas ratio
        pp     = -3.5
        do i=1,NSIZE
          aref(i) = EXP(LOG(a1ref)+(i-1.0)/(NSIZE-1.0)*LOG(a2ref/a1ref))
     >              *mic   
          fref(i) = aref(i)**pp
        enddo
        awref = 0.0
        do i=2,NSIZE
          da = 0.5*(aref(i)-aref(i-1))
          awref(i)   = awref(i)   + da
          awref(i-1) = awref(i-1) + da
        enddo
        ndref = 0.0
        Vref  = 0.0
        mref  = 0.0
        do i=1,NSIZE
          VV = 4.0*pi/3.0*aref(i)**3
          mm = VV*rhoref
          ndref = ndref + fref(i)*awref(i)
          Vref  = Vref  + fref(i)*VV*awref(i)
          mref  = mref  + fref(i)*mm*awref(i)
        enddo
        nHref = 1.0                  ! reference n<H>
        dg    = mref/(nHref*muH)     ! dust/gas mass ratio ...
        scale = dgref/dg             ! ... which should be 0.004
        fref  = scale*fref           ! [cm-4]
        ndref = scale*ndref          ! [cm-3]
        Vref  = scale*Vref           ! [cm3/H-nucleus]
        firstCall = .false.
      endif
        
      !------------------------------------------------------
      ! ***  adjust dust sizes and volume while nd=const  ***
      ! ***  we assume that each grain looses or gains    ***
      ! ***  the same volume fraction.                    ***               
      !------------------------------------------------------
      scale   = (Vcon/(Vref*nHges/nHref))**(1.0/3.0)
      aa      = scale*aref
      ff      = fref*nHges/nHref/scale
      aweight = scale*awref
      !--- done, this is just a check ---
      !ndtest  = 0.0
      !Vtest   = 0.0
      !mtest   = 0.0
      !do i=1,NSIZE
      !  VV = 4.0*pi/3.0*aa(i)**3
      !  mm = VV*rhogr
      !  ndtest = ndtest + ff(i)*aweight(i)
      !  Vtest  = Vtest  + ff(i)*VV*aweight(i)
      !  mtest  = mtest  + ff(i)*mm*aweight(i)
      !enddo
      !dg = mtest/(nHges*muH)
      !if (verb>=0) then
      !  print*,"scale=",scale
      !  print*,"   nd=",ndref*nHges/nHref,ndtest
      !  print*,"  d/g=",dgref,dg
      !  print*,"Vdust=",Vref*nHges/nHref,Vtest,Vcon
      !endif
      
      end
      
      
!---------------------------------------------
      SUBROUTINE sort(n,arr)
!---------------------------------------------
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
