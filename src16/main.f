***********************************************************************
      PROGRAM EQ_CHEMISTRY
***********************************************************************
      use EXCHANGE,ONLY: chemcall,chemiter
      chemcall = 0
      chemiter = 0

      call INIT
      call INIT_CHEMISTRY
      call INIT_DUSTCHEM
      !call DEMO_CHEMISTRY
      call DEMO_SWEEP
      !call DEMO_STRUCTURE
      !call DEMO_PHASEDIAGRAM
      
      print*
      print'("   smchem calls = ",I7)',chemcall
      print'("iterations/call = ",0pF7.2)',
     >     REAL(chemiter)/REAL(chemcall)

      call SAVE_DBASE
      end


***********************************************************************
      SUBROUTINE DEMO_CHEMISTRY
***********************************************************************
      use CHEMISTRY,ONLY: NMOLE,NELM,m_kind,elnum,cmol,el
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,
     >                    dust_nam,dust_mass,dust_Vol
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: nges,kT,nmax,threshold
      real*8  :: Tg,nHges,Vtot,mtot
      integer :: i,imol,iraus,e,j,verbose
      logical :: included,haeufig,raus(NMOLE)
      character(len=2) :: search 
      character(len=10) :: sp
      character(len=2),external :: upper

   10 continue
      write(*,*)
      write(*,*)
      write(*,*) 'Eingabe von Tgas [K], nHtot [cm^-3]'
      read (*,*) Tg, nHges
      eps = eps0
      verbose = 0

*     ------------------------------------------------
      call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
      call GGCHEM(nHges,Tg,eps,.false.,verbose)
*     ------------------------------------------------

      write(*,*) '----- total particle densities -----'
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
        write(*,'(" n<",A2,">=",1pE10.4,2x,1pE10.4)')
     >      elnam(i),nHges*eps(i),eps(i)/eps0(i)
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
 200  continue
      iraus = 0
      nmax  = 0.Q0
      do i=1,NMOLE
        if ((nmol(i).gt.nmax).and.(.not.raus(i))) then
          iraus = i
          nmax  = nmol(i)
        endif
      enddo
      haeufig = (nmax.gt.nHges*1.Q-8)
      if (haeufig) then
        write(*,4010) cmol(iraus), nmol(iraus)
        raus(iraus) = .true.
        goto 200
      endif     

      write(*,*) '-----  where are the elements?  -----'
      do e=1,NELM
        i = elnum(e)
        if (e==el) then
          write(*,'("Element ",A2,1pE12.4)') 'el',0.Q0
          write(*,'(1x,A10,1pE11.4)') "nel       ",nel
          threshold = 1.Q-3*nel
        else   
          write(*,'("Element ",A2,1pE12.4)') elnam(i),eps(i)*nHges 
          threshold = eps(i)*nHges*1.D-3
          if (nat(i).gt.eps(i)*nHges*1.D-3) then
            write(*,'(1x,A10,1pE11.4)') 
     >       "n"//trim(elnam(i))//"       ", nat(i) 
          endif  
        endif  
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
          write(*,4010) cmol(iraus), nmol(iraus)
          raus(iraus) = .true.
        enddo
      enddo  

      write(*,*) '----- gas and electron pressure -----'
      kT = bk*Tg
      nges = nel
      do i=1,NELEM
        nges = nges + nat(i)
      enddo
      do i=1,NMOLE
        nges = nges + nmol(i)
      enddo
      write(*,1010) 'pe=',nel*kT,' pges=',nges*kT
     
*     ------------------------------
      call SUPERSAT(Tg,nat,nmol,Sat)
*     ------------------------------
      write(*,*)
      write(*,*) '----- supersaturation ratios -----'
      do i=1,NDUST
        write(*,5000) dust_nam(i),Sat(i) 
      enddo  

      if (.true.) goto 10

 1000 format(a6,1pE9.3)
 1010 format(a6,1pE9.3,a8,1pE9.3)
 4000 format(a7,1pE10.4,a5,1pE10.4)     
 4010 format(' n',a8,1pE12.4)
 5000 format(1x,a10,' S=',1pE9.3)
      RETURN
      end      
