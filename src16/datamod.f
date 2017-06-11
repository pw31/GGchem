      module DUST_DATA
        integer,parameter :: qp = selected_real_kind ( 33, 4931 )
        integer,parameter :: NELEM=37        ! number of elements
        integer,parameter :: NELM=17         ! no of elements treated in smchem
        integer,parameter :: NMOLE=205       ! number of molecules
        integer,parameter :: NDUST=42        ! number of dust species
        integer,parameter :: NEPS=14         ! number of affected elements

        character(len=2)  :: elnam(NELEM)    ! names of elements
        character(len=10) :: cmol(NMOLE)     ! names of molecules
        character(len=15) :: dust_nam(NDUST) ! names of dust species
        integer :: elnr(NEPS),elcode(NELEM)  ! element cross-indices
        real(kind=qp) :: eps0(NELEM)         ! element abundances
        real*8  :: mass(NELEM)               ! element masses
        real*8  :: dust_rho(NDUST)           ! dust material densities
        real*8  :: dust_mass(NDUST)          ! dust monomer volume
        real*8  :: dust_vol(NDUST)           ! dust monomer volume
        integer :: dust_nel(NDUST)           ! no of elements in dust
        integer :: dust_el(NDUST,5)          ! indices of elements
        integer :: dust_nu(NDUST,5)          ! stoichiometric coeffs

        real(kind=qp) :: bk=1.380662Q-16     ! Boltzman constant
        real(kind=qp) :: bar=1.Q+6           ! 1 bar in dyn/cm2
        real(kind=qp) :: amu=1.66055Q-24     ! atomar mass unit
        real(kind=qp) :: atm=1.013Q+6        ! standard atmosphere pressure
        real(kind=qp) :: rgas=8.31434Q+0     ! gas constant 
        real(kind=qp) :: muH                 ! rho/n<H>
      end

      module STRUCTURE
        integer,parameter :: Npmax=10000 
        integer :: Npoints
        real*8,dimension(Npmax) :: Tgas,press,pelec,dens,nHtot
      end

      module EXCHANGE
        use DUST_DATA,ONLY:NELEM,NMOLE 
        integer,parameter :: qp = selected_real_kind ( 33, 4931 )
        real(kind=qp) :: nel,nat(NELEM),nion(NELEM),nmol(NMOLE)
        integer :: HII,HeII,CII,NII,OII,NaII,MgII,LiII,ClII
        integer :: AlII,KII,TiII,SII,SiII,FeII,CaII
        integer,parameter :: H=1,He=2,Li=3,C=6,N=7,O=8,Fl=9,Ne=10,Na=11
        integer,parameter :: Mg=12,Al=13,Si=14,S=16,Cl=17,K=19,Ca=20
        integer,parameter :: Ti=22,Cr=24,Mn=25,Fe=26,Ni=28
        integer*8 :: chemcall,chemiter
      end
