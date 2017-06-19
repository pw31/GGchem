***********************************************************************
      SUBROUTINE DEMO_PHASEDIAGRAM
***********************************************************************
      use DUST_DATA,ONLY: NELEM,NMOLE,NDUST,elnam,cmol,eps0,bk,bar,muH,
     >                    amu,dust_nam,dust_mass,dust_Vol
      use EXCHANGE,ONLY: nel,nat,nion,nmol,
     >                   H,He,Li,C,N,O,F,Ne,Na,Mg,Al,Si,S,Cl,K,Ca,Ti,
     >                   Cr,Mn,Fe,Ni
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: Npoints=100
      real,dimension(Npoints) :: nHtot,Tgas
      real :: T1,T2,p1,p2,p,pe,Tg,rho,nHges,nges,kT,pges,mu,muold
      real :: nTEA,pTEA
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      integer :: i,ii,j,jj,l,iz,stindex
      character(len=5000) :: species,NISTspecies,elnames
      character(len=200) :: line
      character(len=20) :: frmt,name,short_name(NDUST),test1,test2
      character(len=4) :: sup
      character(len=2) :: test3
      character(len=1) :: char
      integer :: verbose=0
      logical :: isOK

      !---------------------------------
      ! ***  setup sweep parameters  ***
      !---------------------------------
      print*,"start and end temperature [K] (decreasing)"
      !read*,T1,T2
      T1 = 4000
      T2 = 100
      print*,"start and end p [bar] (any)"
      !read*,p1,p2
      p1 = 1E+09*bar
      p2 = 0.1*bar
      mu = 2.3*amu

      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
      enddo
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
      do i=1,Npoints
        p  = EXP(LOG(p1)+LOG(p2/p1)*REAL(i-1)/REAL(Npoints-1)) 
        do ii=1,Npoints
          Tg = EXP(LOG(T1)+LOG(T2/T1)*REAL(ii-1)/REAL(Npoints-1))        
          do 
            nHges = p*mu/(bk*Tg)/muH
            eldust = 0.0
            call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
            call GGCHEM(nHges,Tg,eps,.false.,0)
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
            muold = mu
            mu = nHges/pges*(bk*Tg)*muH
            print '("mu=",2(1pE12.5))',muold/amu,mu/amu
            if (ABS(mu/muold-1.0)<1.E-10) exit
          enddo  

          print*
          print '(i4," Tg[K] =",0pF8.2,"  n<H>[cm-3] =",1pE10.3)',
     >          ii,Tg,nHges

          write(*,1010) ' Tg=',Tg,' pe=',nel*kT,' n<H>=',nHges,
     &                  ' p=',pges/bar,' mu=',mu/amu
          write(70,2010) Tg,nHges,pges,
     &                  LOG10(MAX(1.Q-300, nel)),
     &                  LOG10(MAX(1.Q-300, nat( H))),
     &                  LOG10(MAX(1.Q-300, nat(He))),
     &                  LOG10(MAX(1.Q-300, nat(Li))),
     &                  LOG10(MAX(1.Q-300, nat( C))),
     &                  LOG10(MAX(1.Q-300, nat( N))),
     &                  LOG10(MAX(1.Q-300, nat( O))),
     &                  LOG10(MAX(1.Q-300, nat( F))),
     &                  LOG10(MAX(1.Q-300, nat(Ne))),
     &                  LOG10(MAX(1.Q-300, nat(Na))),
     &                  LOG10(MAX(1.Q-300, nat(Mg))),
     &                  LOG10(MAX(1.Q-300, nat(Al))),
     &                  LOG10(MAX(1.Q-300, nat(Si))),
     &                  LOG10(MAX(1.Q-300, nat( S))),
     &                  LOG10(MAX(1.Q-300, nat(Cl))),
     &                  LOG10(MAX(1.Q-300, nat( K))),
     &                  LOG10(MAX(1.Q-300, nat(Ca))),
     &                  LOG10(MAX(1.Q-300, nat(Ti))),
     &                  LOG10(MAX(1.Q-300, nat(Cr))),
     &                  LOG10(MAX(1.Q-300, nat(Mn))),
     &                  LOG10(MAX(1.Q-300, nat(Fe))),
     &                  LOG10(MAX(1.Q-300, nat(Ni))),
     &                 (LOG10(MAX(1.Q-300, nmol(jj))),jj=1,NMOLE),
     &        (LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj)))),jj=1,NDUST),
     &                 (LOG10(MAX(1.Q-300, eldust(jj))),jj=1,NDUST)

        enddo
      enddo  
      close(70)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(a5,0pF8.2,99(a6,1pE10.3))
 2000 format(999(1x,a12))
 2010 format(0pF13.4,2(1pE13.4),999(0pF13.5))
      end  
