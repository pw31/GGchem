***********************************************************************
      PROGRAM EQ_CHEMISTRY
***********************************************************************
      use PARAMETERS,ONLY: model_dim,model_struc,model_eqcond
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform
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
      print'("         smchem calls = ",I8)',chemcall
      print'("      iterations/call = ",0pF8.2)',
     >                     REAL(chemiter)/REAL(chemcall)
      if (model_eqcond) then
        print'("eq condensation calls = ",I8)',ieqcond
        print'("   eq iterations/call = ",0pF8.2)',
     >                   REAL(ieqconditer)/REAL(ieqcond)
        print'("      transform calls = ",I8)',itransform
        NLAST=0         ! also save replaced database entries
        call SAVE_DBASE
      endif

      end


***********************************************************************
      SUBROUTINE DEMO_CHEMISTRY
***********************************************************************
      use PARAMETERS,ONLY: nHmax,Tmax,pmax,model_pconst,model_eqcond
      use CHEMISTRY,ONLY: NMOLE,NELM,m_kind,elnum,cmol,el
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,muH,
     >                    dust_nam,dust_mass,dust_Vol,dust_nel,dust_el
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: nges,kT,nmax,threshold
      real*8  :: Tg,nHges,p,mu,muold,pgas
      integer :: i,imol,iraus,e,j,verbose,dk
      logical :: included,haeufig,raus(NMOLE)
      character(len=10) :: sp

      Tg    = Tmax
      nHges = nHmax
      p     = pmax
      eps   = eps0
      mu    = muH
      eldust  = 0.Q0
      verbose = 2
      if (model_eqcond) verbose=0

      do 
        if (model_pconst) nHges = p*mu/(bk*Tg)/muH
        if (model_eqcond) then
          call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
        endif  
        call GGCHEM(nHges,Tg,eps,.false.,0)
        nges = nel
        do j=1,NELEM
          nges = nges + nat(j)
        enddo
        do j=1,NMOLE
          nges = nges + nmol(j)
        enddo
        pgas = nges*bk*Tg
        muold = mu
        mu = nHges/pgas*(bk*Tg)*muH
        if (.not.model_pconst) exit
        print '("mu=",2(1pE12.5))',muold/amu,mu/amu
        if (ABS(mu/muold-1.0)<1.E-5) exit
      enddo  

      write(*,*)
      write(*,'("Tg=",0pF8.2,"  rho=",1pE10.3,"  n<H>=",1pE10.3)') 
     >        Tg,nHges*muH,nHges
      write(*,*) '----- total particle densities -----'
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
        write(*,'(" n<",A2,">=",1pE10.4,2x,1pE10.4)')
     >      elnam(i),nHges*eps(i),eps(i)/eps0(i)
      enddo  

      write(*,*) '----- condensates -----'
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
     >                eldust(iraus)*nHges
      enddo
  
      write(*,*) '----- atoms and ions -----'
      write(*,1000) ' nel=',nel
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
        write(*,1010) ' n'//trim(elnam(i))//'I=',nat(i),
     >               '  n'//trim(elnam(i))//'II=',nion(i)
      enddo
  
      write(*,*) '----- some abundant molecules -----'
      raus = .false.
      do
        iraus = 0
        nmax  = 0.Q0
        do i=1,NMOLE
          if ((nmol(i).gt.nmax).and.(.not.raus(i))) then
            iraus = i
            nmax  = nmol(i)
          endif
        enddo
        haeufig = (nmax.gt.nHges*1.Q-5)
        if (.not.haeufig) exit
        raus(iraus) = .true.
        write(*,4010) cmol(iraus), nmol(iraus)
      enddo
  
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

      write(*,*) '----- gas and electron pressure -----'
      write(*,'("pgas[bar]=",1pE10.3,"  pe[bar]=",1pE10.3)') 
     >      nges*kT/bar,nel*kT/bar
     
*     ------------------------------
      call SUPERSAT(Tg,nat,nmol,Sat)
*     ------------------------------
      write(*,*)
      write(*,*) '----- supersaturation ratios -----'
      do i=1,NDUST
        if (Sat(i)<1.Q-2) cycle 
        write(*,5000) dust_nam(i),Sat(i) 
      enddo  

 1000 format(a6,1pE9.3)
 1010 format(a6,1pE9.3,a8,1pE9.3)
 1020 format(a20,1pE9.3)
 4000 format(a7,1pE10.4,a5,1pE10.4)     
 4010 format(' n',a8,1pE12.4)
 5000 format(1x,a20,' S=',1pE9.3)
      RETURN
      end      
