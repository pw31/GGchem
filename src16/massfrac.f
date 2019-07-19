************************************************************************
      subroutine eps2mf(eps,mf)
************************************************************************
      use PARAMETERS,ONLY: elements
      use DUST_DATA,ONLY: NELEM,elnam,mass,amu
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),intent(in) :: eps(NELEM)   ! with respect to H nuclei
      real(kind=qp),intent(out) :: mf(NELEM)   ! mass fractions
      real(kind=qp) :: msum
      integer :: e,el

      msum = 0.Q0
      do el=1,NELEM
        if (index(elements," "//trim(elnam(el))//" ")<=0) cycle
        msum = msum + eps(el)*mass(el)
      enddo  
      mf = 0.Q0
      do el=1,NELEM
        if (index(elements," "//trim(elnam(el))//" ")<=0) cycle
        mf(el) = eps(el)*mass(el)/msum
        !print'(A3,2(1pE12.4))',elnam(el),eps(el),mf(el)
      enddo
      end

************************************************************************
      subroutine mf2eps(mf,eps)
************************************************************************
      use PARAMETERS,ONLY: elements
      use DUST_DATA,ONLY: NELEM,elnam,mass
      use EXCHANGE,ONLY: H
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),intent(in)  :: mf(NELEM)  ! mass fractions
      real(kind=qp),intent(out) :: eps(NELEM) ! with respect to H nuclei
      real(kind=qp) :: epsH
      integer :: el

      eps = 0.Q0
      do el=1,NELEM
        if (index(elements," "//trim(elnam(el))//" ")<=0) cycle
        eps(el) = mf(el)/mass(el) 
      enddo  
      epsH = eps(H)
      eps  = eps/epsH
      end
