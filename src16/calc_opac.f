***********************************************************************
      SUBROUTINE CALC_OPAC(nHges,eldust,T)
***********************************************************************
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,
     >                    muH,mass,mel,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      implicit none
      integer,parameter :: qp=selected_real_kind(33,4931)
      real,intent(in) :: nHges,T
      real(kind=qp),intent(in) :: eldust(NDUST)
      real :: rhod,Vcon,Vrel(NDUST)
      integer :: i
      
      rhod = 0.d0               ! dust mass density [g/cm3]  
      Vcon = 0.d0               ! dust volume [cm3/cm3]
      do i=1,NDUST
        if (eldust(i)<=0.Q0) cycle 
        rhod = rhod + nHges*eldust(i)*dust_mass(i)
        Vcon = Vcon + nHges*eldust(i)*dust_Vol(i)
      enddo
      do i=1,NDUST
        Vrel(i) = nHges*eldust(i)*dust_Vol(i)/Vcon
      enddo

      end
