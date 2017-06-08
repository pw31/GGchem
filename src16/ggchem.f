***********************************************************************
      SUBROUTINE GGCHEM (nHges,Tg,eps,merk,verbose)
***********************************************************************
*****                                                             *****
*****  Ruft lediglich SMCHEM auf (mit kompatiber Datenstruktur)   *****
*****                                                             *****
***********************************************************************
      use DUST_DATA,ONLY: NELEM,NELM,NMOLE,cmol 
      use EXCHANGE,ONLY: nel,nat,nion,nmol,
     &                   HII,CII,NII,OII,NaII,MgII,AlII,KII,TiII,SII,
     &                   SiII,FeII,CaII,LiII,ClII,HeII,
     &                   H,He,Li,C,N,O,Fl,Ne,Na,Mg,Al,Si,S,Cl,K,Ca,Ti,
     &                   Cr,Mn,Fe,Ni
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nHges,Tg
      real*8 :: epsi8(NELM),anmono8(NELM),nmol8(NMOLE)
      real(kind=qp),intent(in) :: eps(NELEM)
      real(kind=qp) :: epsi(NELM),anmono(NELM)
      logical,intent(in) :: merk
      integer,intent(in) :: verbose
      integer :: verb
      logical :: mk

      epsi( 1) = eps(He)
      epsi( 2) = 0.Q0
      epsi( 3) = eps( H)
      epsi( 4) = eps( C)
      epsi( 5) = eps( N)
      epsi( 6) = eps( O)
      epsi( 7) = eps(Si)
      epsi( 8) = eps(Mg)
      epsi( 9) = eps(Al)
      epsi(10) = eps(Fe)
      epsi(11) = eps( S)
      epsi(12) = eps(Na)
      epsi(13) = eps( K)
      epsi(14) = eps(Ti)
      epsi(15) = eps(Ca)
      epsi(16) = eps(Li)
      epsi(17) = eps(Cl)
      verb = verbose
      mk = merk

      if (Tg<1000.d0) then
        call SMCHEM16(nHges, Tg, epsi, anmono, nmol, mk, verb)
      else
        epsi8 = epsi
        call SMCHEM8 (nHges, Tg,epsi8,anmono8,nmol8, mk, verb)
        anmono = anmono8
        nmol   = nmol8
      endif  

      nat(He)  = anmono( 1)
      nel      = anmono( 2)
      nat( H)  = anmono( 3)
      nat( C)  = anmono( 4)
      nat( N)  = anmono( 5)
      nat( O)  = anmono( 6)
      nat(Si)  = anmono( 7)
      nat(Mg)  = anmono( 8)
      nat(Al)  = anmono( 9)
      nat(Fe)  = anmono(10)
      nat( S)  = anmono(11)
      nat(Na)  = anmono(12)
      nat( K)  = anmono(13)
      nat(Ti)  = anmono(14)
      nat(Ca)  = anmono(15)
      nat(Li)  = anmono(16)
      nat(Cl)  = anmono(17)

      nion( H) = nmol(HII)
      nion(He) = nmol(HeII)
      nion( C) = nmol(CII)
      nion( N) = nmol(NII)
      nion( O) = nmol(OII)
      nion(Si) = nmol(SiII)
      nion(Mg) = nmol(MgII)
      nion(Al) = nmol(AlII)
      nion(Fe) = nmol(FeII)
      nion(S ) = nmol(SII)
      nion(Na) = nmol(NaII)
      nion( K) = nmol(KII)
      nion(Ti) = nmol(TiII)
      nion(Ca) = nmol(CaII)
      nion(Li) = nmol(LiII)
      nion(Cl) = nmol(ClII)
             
      RETURN
      end 




