**********************************************************************
      SUBROUTINE INIT_DUSTCHEM
**********************************************************************
      use dust_data,ONLY: NEPS,NELEM,NDUST,NMOLE,eps0,
     &                    dust_nam,dust_rho,dust_vol,dust_mass,
     &                    dust_nel,dust_nu,dust_el,
     &                    elnr,elcode,elnam,mass
      implicit none
      integer :: i,j,k
      real*8 :: dmass
      character(len=200):: zeile
      character(len=2)  :: name2
      logical :: is_atom,found

      call GGCHEM(1.d+15,2000.d0,eps0,.false.,0)   ! damit cmol vorliegt
 
      write(*,*) 
      write(*,*) "reading DustChem57.dat ..."
      write(*,*) "=========================="

      open(12, file='DustChem57.dat', status='old')
 
      write(*,*) '--- involved elements ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NEPS.ne.i) stop 'Parameter NEPS inconsistent.'
      do i=1,NEPS
        read(12,1010) name2
        found = .false. 
        do j=1,NELEM
          if (name2.eq.elnam(j)) then
            found = .true. 
            elnr(i) = j
            elcode(j) = i
          endif
        enddo
        if (.not.found) stop 'Element not found.'
        write(*,*) elcode(elnr(i)),' ',name2,elnr(i)
      enddo
 
      write(*,*) '--- dust species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) i
      if (NDUST.ne.i) stop 'Parameter NDUST incorrect.'
      do i=1,NDUST
        read(12,*) dust_nam(i)
        read(12,*) dust_rho(i)
        read(12,*) dust_nel(i)
        dmass = 0.d0
        do j=1,dust_nel(i)
          read(12,1030) dust_nu(i,j),name2
          found = .false. 
          do k=1,NELEM
            if (elnam(k).eq.name2) then
              dust_el(i,j) = k
              dmass = dmass + dust_nu(i,j)*mass(k)
              found = .true.
            endif
          enddo
          if (.not.found) stop 'Element in dust species not found'
        enddo
        dust_mass(i) = dmass
        dust_vol(i) = dmass/dust_rho(i)
        write(*,1060) dust_nam(i), dust_rho(i), dust_vol(i), 
     &      (dust_nu(i,j),elcode(dust_el(i,j)),j=1,dust_nel(i))
      enddo

      RETURN 
 1000 format(a200)
 1010 format(a2)
 1020 format(2(l1,1x),i1,1x,a10)
 1030 format(i2,1x,a2)
 1040 format(i2,1x,a10)
 1050 format(1x,a10,i4,' mass=',0pf7.3," amu")
 1060 format(1x,a10," rhod=",0pf6.3," V0=",1pe11.3,2x,99(i1,"x",i2,1x))
 1070 format(1x,a10,99(i1,"x",i2,1x))
 2011 format(1(I2,1x,a8),22x,'->',I2,1x,a10,99(I2,1x,a8))
 2021 format(2(I2,1x,a8),11x,'->',I2,1x,a10,99(I2,1x,a8))
 2031 format(3(I2,1x,a8)    ,'->',I2,1x,a10,99(I2,1x,a8))
      end 
