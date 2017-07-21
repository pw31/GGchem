**********************************************************************
      SUBROUTINE INIT_DUSTCHEM
**********************************************************************
      use CHEMISTRY,ONLY: NMOLE,NELM,catm
      use DUST_DATA,ONLY: NEPS,NELEM,NDUST,eps0,
     &                    dust_nam,dust_rho,dust_vol,dust_mass,
     &                    dust_nel,dust_nu,dust_el,
     &                    elnr,elcode,elnam,mass
      implicit none
      integer :: i,imax,j,k,el
      real*8 :: dmass
      character(len=200):: zeile
      character(len=2)  :: name
      logical :: is_atom,found,allfound

      write(*,*) 
      write(*,*) "reading DustChem.dat ..."
      write(*,*) "========================"

      open(12, file='DustChem.dat', status='old')
 
      write(*,*) '--- dust species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) imax
      NDUST = 1
      do i=1,imax
        read(12,1000) zeile
        read(12,*) dust_nam(NDUST)
        read(12,*) dust_rho(NDUST)
        read(12,*) dust_nel(NDUST)
        dmass = 0.d0
        allfound = .true.
        do j=1,dust_nel(NDUST)
          read(12,1030) dust_nu(NDUST,j),name
          found = .false. 
          do k=1,NELEM
            if (elnam(k).eq.name) then
              dust_el(NDUST,j) = k
              dmass = dmass + dust_nu(NDUST,j)*mass(k)
              found = .true.
            endif
          enddo
          if (.not.found) stop 'Element in dust species not found'
          found = .false.
          do k=1,NELM
            if (catm(k).eq.name) then
              found = .true.
              exit
            endif
          enddo
          if (.not.found) allfound=.false.
        enddo
        if (allfound) then
          dust_mass(NDUST) = dmass
          dust_vol(NDUST) = dmass/dust_rho(NDUST)
          write(*,1060) dust_nam(NDUST),dust_rho(NDUST),dust_vol(NDUST), 
     &    (dust_nu(NDUST,j),elnam(dust_el(NDUST,j)),j=1,dust_nel(NDUST))
          NDUST = NDUST+1
        endif  
      enddo
      NDUST=NDUST-1

      write(*,*) '--- involved elements ---'
      NEPS=0
      elcode(:)=0
      do i=1,NDUST
        do j=1,dust_nel(i)
          name = elnam(dust_el(i,j)) 
          do k=1,NELEM
            if (elnam(k).eq.name) then
              el = k
              exit
            endif
          enddo
          found = .false.           
          do k=1,NEPS
            if (el==elnr(k)) found=.true.
          enddo
          if (.not.found) then
            NEPS = NEPS+1 
            elnr(NEPS) = el
            elcode(el) = NEPS
            write(*,*) elcode(elnr(NEPS)),' ',name,el
          endif
        enddo
      enddo
      write(*,*)
 
      RETURN 
 1000 format(a200)
 1010 format(a2)
 1020 format(2(l1,1x),i1,1x,a10)
 1030 format(i2,1x,a2)
 1040 format(i2,1x,a10)
 1050 format(1x,a10,i4,' mass=',0pf7.3," amu")
 1060 format(1x,a14," rhod=",0pf6.3," V0=",1pe11.3,2x,99(i1,"x",A2,1x))
 1070 format(1x,a10,99(i1,"x",i2,1x))
 2011 format(1(I2,1x,a8),22x,'->',I2,1x,a10,99(I2,1x,a8))
 2021 format(2(I2,1x,a8),11x,'->',I2,1x,a10,99(I2,1x,a8))
 2031 format(3(I2,1x,a8)    ,'->',I2,1x,a10,99(I2,1x,a8))
      end 
