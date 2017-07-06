***********************************************************************
      SUBROUTINE GGCHEM(nHges,Tg,eps,merk,verbose)
***********************************************************************
*****                                                             *****
*****  Ruft lediglich SMCHEM auf (mit kompatiber Datenstruktur)   *****
*****                                                             *****
***********************************************************************
      use PARAMETERS,ONLY: Tfast
      use CHEMISTRY,ONLY: NMOLE,NELEM,NELM,elnum,elion,el,charge,catm
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nHges,Tg
      real*8 :: epsi8(NELM),anmono8(NELM),nmol8(NMOLE)
      real(kind=qp),intent(in) :: eps(NELEM)
      real(kind=qp) :: epsi(NELM),anmono(NELM)
      logical,intent(in) :: merk
      integer,intent(in) :: verbose
      integer :: i,verb
      logical :: mk

      if (charge) epsi(el)=0.Q0
      do i=1,NELM
        if (i==el) cycle
        epsi(i) = eps(elnum(i))
        !print*,i,catm(i),elnum(i),epsi(i)
      enddo  
      verb = verbose
      mk = merk

      if (Tg<Tfast) then
        call SMCHEM16(nHges, Tg, epsi, anmono, nmol, mk, verb)
      else
        epsi8 = epsi
        call SMCHEM8 (nHges, Tg,epsi8,anmono8,nmol8, mk, verb)
        anmono = anmono8
        nmol   = nmol8
      endif  

      if (charge) nel=anmono(el)
      do i=1,NELM
        if (i==el) cycle                   ! slot for electrons
        nat(elnum(i))  = anmono(i)         ! neutral atoms
        if (charge) then
          nion(elnum(i)) = nmol(elion(i))  ! positive ions
        endif  
      enddo  
             
      RETURN
      end 




