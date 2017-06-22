**********************************************************************
      SUBROUTINE INIT
**********************************************************************
*****                                                            *****
*****   Initialisiert Elementhaeufigkeiten                       *****
*****   - Anders + Grevesse (1989):                              *****
*****     Geochimica et Cosmochemica Acta Vol 53, pp 197--214    *****
*****     ("Photosphere")                                        *****
*****   - wie in MARCS-Code                                      *****
*****   - wie in Tsuji-Chemie                                    *****
*****                                                            *****
**********************************************************************
      use DUST_DATA,ONLY: NELEM,eps=>eps0,mass,muH,elnam,amu
      use EXCHANGE,ONLY: H,He,Li,C,N,O,F,Ne,Na,Mg,Al,Si,S,Cl,K,Ca,Ti,
     >                   Cr,Mn,Fe,Ni
      implicit none
      integer :: i,nr,pick
      real*8 :: m,abund(NELEM,4)
      character(len=2) :: el
      character(len=20) :: elname
      character(len=20) :: source(4)
      character(len=200) :: line

      write(*,*) 
      write(*,*) "elemental abundances ..."
      write(*,*) "========================"
      do i=1,NELEM
        elnam(i) = '  '
        eps(i)  = -99.d0
        mass(i) = 0.d0
      enddo
      elnam(1)  = 'H '
      elnam(2)  = 'He'
      elnam(3)  = 'Li'
      elnam(6)  = 'C '
      elnam(7)  = 'N '
      elnam(8)  = 'O '
      elnam(9)  = 'F '
      elnam(10) = 'Ne'
      elnam(11) = 'Na'
      elnam(12) = 'Mg'
      elnam(13) = 'Al'
      elnam(14) = 'Si'
      elnam(16) = 'S '
      elnam(17) = 'Cl'
      elnam(19) = 'K '
      elnam(20) = 'Ca'
      elnam(22) = 'Ti'
      elnam(24) = 'Cr'
      elnam(25) = 'Mn'
      elnam(26) = 'Fe'
      elnam(28) = 'Ni'
       
*     ---------------------------------------
*     Grevesse + Noels (1996, "photosphere"):
*     ---------------------------------------
      eps(H)  = 12.00 D0
      eps(He) = 10.99 D0
      eps(Li) =  1.16 D0
      eps(C)  =  8.55 D0
      eps(N)  =  7.97 D0
      eps(O)  =  8.87 D0
      eps(Ne) =  8.08 D0
      eps(Na) =  6.33 D0
      eps(Mg) =  7.58 D0
      eps(Al) =  6.47 D0
      eps(Si) =  7.55 D0
      eps(S)  =  7.33 D0
      eps(Cl) =  5.50 D0
      eps(K)  =  5.12 D0
      eps(Ca) =  6.36 D0
      eps(Ti) =  5.02 D0
      eps(Cr) =  5.67 D0
      eps(Mn) =  5.39 D0
      eps(Fe) =  7.50 D0
      eps(Ni) =  6.25 D0

*     ---------------------------------
*     Grevesse, Asplund, Sauval (2007):
*     ---------------------------------
      eps(H)  = 12.00 D0
      eps(He) = 10.93 D0
      eps(Li) =  1.10 D0  ! Lodders, Palme Gail 2009
      eps(C)  =  8.39 D0  
      eps(N)  =  7.78 D0  
      eps(O)  =  8.66 D0  
      eps(F)  =  4.56 D0
      eps(Ne) =  7.84 D0
      eps(Na) =  6.17 D0
      eps(Mg) =  7.53 D0
      eps(Al) =  6.37 D0
      eps(Si) =  7.51 D0
      eps(S)  =  7.14 D0
      eps(Cl) =  5.50 D0
      eps(K)  =  5.08 D0
      eps(Ca) =  6.31 D0
      eps(Ti) =  4.90 D0
      eps(Cr) =  5.64 D0
      eps(Mn) =  5.39 D0
      eps(Fe) =  7.45 D0
      eps(Ni) =  6.23 D0

      !eps(C)  = eps(O) + LOG10(2.0)   ! try C/O=2

      !eps(Si) = eps(H)+20.0           ! pure SiO2 modelling ...
      !eps(O)  = eps(Si)+LOG10(2.0)    ! Si:O = 1:2

      !eps(:) = eps(:)-40.0            ! pure H2O modelling ...
      !eps(H) = 12.0
      !eps(O) = eps(H)-LOG10(2.0)      ! H:O = 2:1

      !do i=1,NELEM
      !  if (elnam(i).ne.'  ') then
      !    if ((elnam(i).ne.'H' ).and.(elnam(i).ne.'He').and.
     >!        (elnam(i).ne.'C' ).and.(elnam(i).ne.'O' ).and.
     >!        (elnam(i).ne.'N' ).and.(elnam(i).ne.'Si').and.
     >!        (elnam(i).ne.'S' ).and.
     >!        (elnam(i).ne.'Na').and.(elnam(i).ne.'Ca').and.
     >!        (elnam(i).ne.'Cl').and.(elnam(i).ne.'Ti')) then
      !      eps(i) = eps(i)-30.0
      !    endif
      !  endif
      !enddo

*     --------------------
*     ***  Atommassen  ***
*     --------------------
      mass(H)  = 1.008  D0 * amu
      mass(He) = 4.0026 D0 * amu
      mass(Li) = 6.94   D0 * amu
      mass(C)  = 12.011 D0 * amu
      mass(N)  = 14.007 D0 * amu
      mass(O)  = 15.999 D0 * amu
      mass(F)  = 19.000 D0 * amu      
      mass(Ne) = 20.18  D0 * amu            
      mass(Na) = 22.990 D0 * amu
      mass(Mg) = 24.312 D0 * amu
      mass(Al) = 26.98  D0 * amu
      mass(Si) = 28.086 D0 * amu
      mass(S)  = 32.064 D0 * amu
      mass(Cl) = 35.45  D0 * amu
      mass(K)  = 39.10  D0 * amu  
      mass(Ca) = 40.08  D0 * amu
      mass(Ti) = 47.90  D0 * amu
      mass(Cr) = 52.00  D0 * amu
      mass(Mn) = 54.94  D0 * amu
      mass(Fe) = 55.847 D0 * amu
      mass(Ni) = 58.71  D0 * amu

      muH = 0.d0
      do i=1,NELEM
        eps(i) = 10.d0 ** (eps(i)-12.d0)
        if (elnam(i).ne.'  ') then
          write(*,'(1x,a2,1x,0pF8.3,1pE12.4)') 
     &          elnam(i),12.d0+LOG10(eps(i)),eps(i)
        endif  
        muH = muH + mass(i)*eps(i)
      enddo
      write(*,*) 'rho = n<H> *',muH/amu,' amu'

*     ------------------------------------
*     ***  read abundances from file?  ***      
*     ------------------------------------
      if (.true.) then
        source = (/'EarthCrust','Ocean','Solar','Meteorites'/)
        pick = 1
        write(*,*)
        write(*,*) "replacing from file Abundances.dat ("
     &             //trim(source(pick))//") ..."
        open(1,file='Abundances.dat',status='old')
        do i=1,5
          read(1,'(A200)') line
        enddo  
        do i=1,NELEM
          read(1,*) nr,elname,el,m,abund(nr,1:4)
          !print*,nr,el,real(eps(nr)),abund(nr,pick)/abund(1,pick)
          mass(nr)  = m*amu
          elnam(nr) = el
          eps(nr) = MAX(1.e-99,abund(nr,pick)/abund(1,pick))
        enddo  
        close(1)
        muH = 0.d0
        do i=1,NELEM
          write(*,'(1x,a2,1x,0pF8.3,1pE12.4)') 
     &          elnam(i),12.d0+LOG10(eps(i)),eps(i)
          muH = muH + mass(i)*eps(i)
        enddo
        write(*,*) 'rho = n<H> *',muH/amu,' amu'
      endif  

      RETURN
      end
