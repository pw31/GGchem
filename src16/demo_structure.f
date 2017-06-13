***********************************************************************
      SUBROUTINE DEMO_STRUCTURE
***********************************************************************
      use STRUCTURE,ONLY: Npmax,Npoints,Tgas,press,pelec,dens,nHtot
      use DUST_DATA,ONLY: NELEM,NMOLE,NDUST,elnam,cmol,eps0,bk,bar,amu,
     >                    dust_nam,dust_mass,dust_Vol,muH
      use EXCHANGE,ONLY: nel,nat,nion,nmol,
     >                   H,He,Li,C,N,O,F,Ne,Na,Mg,Al,Si,S,Cl,K,Ca,Ti,
     >                   Cr,Mn,Fe,Ni
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real :: tau,p,pe,Tg,rho,nHges,nges,kT,pges,mu
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      integer :: i,j,jj,iz
      character(len=200) :: line
      character(len=20) :: name,short_name(NDUST)
      character(len=1) :: char
      integer :: verbose=0

      !-----------------------------
      ! ***  read the structure  ***
      !-----------------------------
      open(3,file='lte_1200_3.0-0.0_cutforpegas',status='old')
      Npoints = 256 
      do i=1,5+Npoints+3 
        read(3,'(A200)') line
      enddo  
      do i=1,Npoints 
        read(3,*) iz,tau,Tg,p,pe,rho
        Tgas(i)  = Tg
        press(i) = p
        pelec(i) = pe
        dens(i)  = rho
        nHtot(i) = rho/muH
      enddo
      close(3)

      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[")
        short_name(i) = name(1:j-1)
      enddo
      eps = eps0
      eps = eps0
      open(unit=70,file='Static_Conc.dat',status='replace')
      write(70,1000) 'H',eps( H), 'C',eps( C),
     &               'N',eps( N), 'O',eps( O)
      write(70,*) NELEM,NMOLE,NDUST,Npoints
      write(70,2000) 'Tg','nHges','pges','el',
     &               'H','He','Li','C','N','O','F','Ne','Na','Mg','Al',
     &               'Si','S','Cl','K','Ca','Ti','Cr','Mn','Fe','Ni',
     &               (trim(cmol(i)),i=1,NMOLE),
     &               ('S'//trim(short_name(i)),i=1,NDUST),
     &               ('n'//trim(short_name(i)),i=1,NDUST)

      !-------------------------------------
      ! ***  run chemistry on structure  ***
      !-------------------------------------
      do i=Npoints,2,-1     ! inside-out, avoid outermost point 
        nHges = nHtot(i)
        Tg = Tgas(i)
        print*
        print '(i4," Tg[K] =",0pF8.2," n<H>[cm-3] =",1pE10.3)',
     >        i,Tg,nHges

        eldust = 0.0
        call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
        call GGCHEM(nHges,Tg,eps,.true.,0)
        call SUPERSAT(Tg,nat,nmol,Sat)

        kT = bk*Tg
        nges = nel
        do j=1,NELEM
          nges = nges + nat(j)
        enddo
        do j=1,NMOLE
          nges = nges + nmol(j)
        enddo
        pges = nges*kT
        mu = nHges/pges*(bk*Tg)*muH
        write(*,1010) ' Tg=',Tg,' pe=',nel*kT,' n<H>=',nHges,
     &                ' p=',pges/bar,' mu=',mu/amu

        write(70,2010) Tg,nHges,pges,
     &                LOG10(MAX(1.Q-300, nel)),
     &                LOG10(MAX(1.Q-300, nat( H))),
     &                LOG10(MAX(1.Q-300, nat(He))),
     &                LOG10(MAX(1.Q-300, nat(Li))),
     &                LOG10(MAX(1.Q-300, nat( C))),
     &                LOG10(MAX(1.Q-300, nat( N))),
     &                LOG10(MAX(1.Q-300, nat( O))),
     &                LOG10(MAX(1.Q-300, nat( F))),
     &                LOG10(MAX(1.Q-300, nat(Ne))),
     &                LOG10(MAX(1.Q-300, nat(Na))),
     &                LOG10(MAX(1.Q-300, nat(Mg))),
     &                LOG10(MAX(1.Q-300, nat(Al))),
     &                LOG10(MAX(1.Q-300, nat(Si))),
     &                LOG10(MAX(1.Q-300, nat( S))),
     &                LOG10(MAX(1.Q-300, nat(Cl))),
     &                LOG10(MAX(1.Q-300, nat( K))),
     &                LOG10(MAX(1.Q-300, nat(Ca))),
     &                LOG10(MAX(1.Q-300, nat(Ti))),
     &                LOG10(MAX(1.Q-300, nat(Cr))),
     &                LOG10(MAX(1.Q-300, nat(Mn))),
     &                LOG10(MAX(1.Q-300, nat(Fe))),
     &                LOG10(MAX(1.Q-300, nat(Ni))),
     &               (LOG10(MAX(1.Q-300, nmol(jj))),jj=1,NMOLE),
     &               (LOG10(Sat(jj)),jj=1,NDUST),
     &               (LOG10(MAX(1.Q-300, eldust(jj))),jj=1,NDUST)

        if (verbose>0) read(*,'(a1)') char

      enddo  

      close(70)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(a5,0pF8.2,99(a6,1pE10.3))
 2000 format(999(1x,a12))
 2010 format(0pF13.4,2(1pE13.4),999(0pF13.5))
      end  

