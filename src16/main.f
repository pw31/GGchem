      PROGRAM EQ_CHEMISTRY
        use EXCHANGE,ONLY: chemcall,chemiter
        chemcall = 0
        chemiter = 0

        call INIT
        call INIT_DUSTCHEM
        !call DEMO_CHEMISTRY
        call DEMO_SWEEP
        !call DEMO_STRUCTURE

        print*
        print'("   smchem calls = ",I7)',chemcall
        print'("iterations/call = ",0pF7.2)',
     >              REAL(chemiter)/REAL(chemcall)
        call SAVE_DBASE

      end


***********************************************************************
      SUBROUTINE DEMO_CHEMISTRY
***********************************************************************
      use DUST_DATA,ONLY: NELEM,NMOLE,NDUST,elnam,cmol,eps0,bk,
     >                    dust_nam,dust_mass,dust_Vol
      use EXCHANGE,ONLY: nel,nat,nion,nmol,
     >                   H,He,Li,C,N,O,Fl,Ne,Na,Mg,Al,Si,S,Cl,K,Ca,Ti,
     >                   Cr,Mn,Fe,Ni
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: nges,kT,nmax
      real*8  :: Tg,nHges,Vtot,mtot
      integer :: i,iraus,el,j,verbose
      logical :: haeufig,raus(NMOLE)
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
      !call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
      call GGCHEM(nHges,Tg,eps,.false.,2)
*     ------------------------------------------------

      write(*,*) '----- total particle densities -----'
      do i=1,NELEM
        if (eps0(i).gt.1.E-20) then 
          write(*,'(" n<",A2,">=",1pE10.4,2x,1pE10.4)')
     >      elnam(i),nHges*eps(i),eps(i)/eps0(i)
        endif
      enddo  

      write(*,*) '----- atoms and ions -----'
      write(*,1000) ' nel=',nel
      write(*,1010) ' nHI=',nat(H) ,'  nHII=',nion(H)
      write(*,1010) 'nHeI=',nat(He),' nHeII=',nion(He)
      write(*,1010) ' nCI=',nat(C) ,'  nCII=',nion(C)
      write(*,1010) ' nNI=',nat(N) ,'  nNII=',nion(N)
      write(*,1010) ' nOI=',nat(O) ,'  nOII=',nion(O)
      write(*,1010) 'nSiI=',nat(Si),' nSiII=',nion(Si)
      write(*,1010) 'nMgI=',nat(Mg),' nMgII=',nion(Mg)
      write(*,1010) 'nAlI=',nat(Al),' nAlII=',nion(Al)
      write(*,1010) 'nFeI=',nat(Fe),' nFeII=',nion(Fe)
      write(*,1010) ' nSI=',nat(S) ,'  nSII=',nion(S)
      write(*,1010) 'nNaI=',nat(Na),' nNaII=',nion(Na)
      write(*,1010) ' nKI=',nat(K) ,'  nKII=',nion(K)
      write(*,1010) 'nTiI=',nat(Ti),' nTiII=',nion(Ti)
      write(*,1010) 'nCaI=',nat(Ca),' nCaII=',nion(Ca)
      write(*,1010) 'nLiI=',nat(Li),' nLiII=',nion(Li)
      write(*,1010) 'nClI=',nat(Cl),' nClII=',nion(Cl)

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
      do el=1,NELEM
        if (eps0(el).lt.1.E-20) cycle 
        search = upper(elnam(el))        
        write(*,'("Element ",A2,1pE12.4)') elnam(el),eps(el)*nHges 
        if (nat(el).gt.eps(el)*nHges*1.D-3) then
          write(*,'(1x,A10,1pE11.4)') 
     >     "n"//trim(elnam(el))//"       ", nat(el) 
        endif  
        raus = .false.
        do 
          iraus = 0
          nmax  = 0.Q0
          do i=1,NMOLE
            sp = cmol(i) 
            if ((nmol(i).gt.nmax).and.(.not.raus(i))) then
              j = index(sp,trim(search))
              if ((j.gt.0).and.(search.eq.'C')) then
                if (sp(j+1:j+1).eq.'L') j=0   ! CL
                if (sp(j+1:j+1).eq.'A') j=0   ! CA
              endif  
              if ((j.gt.0).and.(search.eq.'N')) then
                if (sp(j+1:j+1).eq.'A') j=0   ! NA
              endif  
              if ((j.gt.0).and.(search.eq.'S')) then
                if (sp(j+1:j+1).eq.'I') j=0   ! SI
              endif  
              if (j.gt.0) then
                iraus = i
                nmax = nmol(i)
              endif  
            endif
          enddo  
          haeufig = (nmax.gt.eps(el)*nHges*1.D-3)
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

 1000 format(a5,1pE9.3)
 1010 format(a5,1pE9.3,a7,1pE9.3)
 4000 format(a7,1pE10.4,a5,1pE10.4)     
 4010 format(' n',a8,1pE12.4)
 5000 format(1x,a10,' S=',1pE9.3)
      RETURN
      end      
