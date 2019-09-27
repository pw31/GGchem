************************************************************************
      module PARAMETERS
************************************************************************
      character(len=200) :: elements,abund_file,struc_file
      integer :: abund_pick,model_dim,Npoints,model_struc,verbose
      logical :: model_eqcond,model_pconst,pick_mfrac,initchem_info
      logical :: useDataBase,remove_condensates,phyllosilicates
      real*8  :: Tfast,Tmin,Tmax,pmin,pmax,nHmin,nHmax
      end

************************************************************************
      module DUST_DATA
************************************************************************
      character(len=200) :: DustChem_file
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: NELEM=41        ! number of elements (up to Zr + W)
      integer,parameter :: NDUSTmax=500    ! max number of condensed species
      integer :: NDUST                     ! number of condensed species
      integer :: NEPS                      ! number of affected elements
      
      character(len=2)  :: elnam(NELEM)       ! names of elements
      character(len=20) :: dust_nam(NDUSTmax) ! names of dust species
      integer :: elnr(NELEM),elcode(NELEM)    ! element cross-indices
      real(kind=qp) :: eps0(NELEM)            ! element abundances
      real*8  :: mass(NELEM)                  ! element masses
      real*8  :: dust_rho(NDUSTmax)           ! dust material densities
      real*8  :: dust_mass(NDUSTmax)          ! dust monomer volume
      real*8  :: dust_vol(NDUSTmax)           ! dust monomer volume
      real*8  :: Tmelt(NDUSTmax)              ! melting points
      real*8  :: Tcorr(NDUSTmax)
      logical :: is_liquid(NDUSTmax)
      integer :: dust_nel(NDUSTmax)           ! no of elements in dust
      integer :: dust_el(NDUSTmax,8)          ! indices of elements
      integer :: dust_nu(NDUSTmax,8)          ! stoichiometric coeffs
      
      integer :: fit(NDUSTmax)                ! fit-formular identifier
      real*8  :: cfit(NDUSTmax,0:4)           ! pvap fit coefficients
      
      real(kind=qp) :: bk=1.380662Q-16        ! Boltzman constant
      real(kind=qp) :: bar=1.Q+6              ! 1 bar in dyn/cm2
      real(kind=qp) :: amu=1.66055Q-24        ! atomar mass unit
      real(kind=qp) :: atm=1.013Q+6           ! standard atmosphere pressure
      real(kind=qp) :: rgas=8.3144598Q+0      ! gas constant 
      real(kind=qp) :: mel=9.109389754Q-28    ! electron mass
      real(kind=qp) :: muH                    ! rho/n<H>
      end

************************************************************************
      module CHEMISTRY
************************************************************************
      use DUST_DATA,ONLY: NELEM
      character(len=200) :: dispol_file(4)
      logical :: NewFullIt
      integer :: NewBackIt,NewFastLevel,NewPreMethod
      real*8  :: NewBackFac
      integer :: NMOLdim         ! max number of molecules
      integer :: NMOLE           ! number of molecules found
      integer :: NELM            ! number of elements found
      integer :: el=0,H=0,He=0,Li=0,Be=0,B=0,C=0,N=0,O=0,F=0,Ne=0
      integer :: Na=0,Mg=0,Al=0,Si=0,P=0,S=0,Cl=0,Ar=0,K=0,Ca=0
      integer :: Sc=0,Ti=0,V=0,Cr=0,Mn=0,Fe=0,Co=0,Ni=0,Cu=0,Zn=0
      integer :: Ga=0,Ge=0,As=0,Se=0,Br=0,Kr=0,Rb=0,Sr=0
      integer :: Y=0,Zr=0,W=0
      logical :: charge
      character(len=2) :: catm(NELEM)           ! names of elements
      character(len=20),allocatable :: cmol(:)  ! names of molecules
      integer :: elnum(NELEM)                   ! indices of found elements
      integer :: elion(NELEM)                   ! indices of ions
      integer,allocatable :: fit(:)             ! fit-formular identifier
      integer,allocatable :: natom(:)           ! no of atoms in molecule    
      integer,allocatable :: source(:)          ! no of source file
      integer,allocatable :: m_kind(:,:)        ! index of elements
      integer,allocatable :: m_anz(:,:)         ! stoichiometric coeffs
      real*8,allocatable  :: a(:,:)             ! kp fit-coeffs
      real*8,allocatable  :: error(:)           ! kp fit errors
      real*8 :: b_nasa(NELEM,0:13)              ! kp fit-coeffs Added by Yui Kawashima
      integer :: i_nasa,c_nasa(NELEM)           ! Added by Yui Kawashima
      real*8 :: th1,th2,th3,th4,TT1,TT2,TT3     
      end

************************************************************************
      module STRUCTURE
************************************************************************
      use DUST_DATA,ONLY: NELEM
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: Npmax=10000 
      real*8,dimension(Npmax) :: Tgas,press,pelec,dens,nHtot
      real(kind=qp) :: estruc(Npmax,NELEM)
      end

************************************************************************
      module EXCHANGE
************************************************************************
      use CHEMISTRY,ONLY: NMOLE
      use DUST_DATA,ONLY: NELEM
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: nel,nat(NELEM),nion(NELEM)
      real(kind=qp),allocatable :: nmol(:),mmol(:)
      integer :: HII,HeII,CII,NII,OII,NaII,MgII,LiII,ClII
      integer :: AlII,KII,TiII,SII,SiII,FeII,CaII
      integer,parameter :: H=1,He=2,Li=3,Be=4,B=5,C=6,N=7,O=8,F=9
      integer,parameter :: Ne=10,Na=11,Mg=12,Al=13,Si=14,P=15,S=16
      integer,parameter :: Cl=17,Ar=18,K=19,Ca=20,Sc=21,Ti=22
      integer,parameter :: V=23,Cr=24,Mn=25,Fe=26,Co=27,Ni=28
      integer,parameter :: Cu=29,Zn=30,Ga=31,Ge=32,As=33,Se=34
      integer,parameter :: Br=35,Kr=36,Rb=37,Sr=38,Y=39,Zr=40,W=41
      integer*8 :: chemcall=0,chemiter=0,itransform=0,ieqcond=0
      integer*8 :: preIter=0,preEst=0,preUse=0
      integer*8 :: ieqconditer=0
      end
