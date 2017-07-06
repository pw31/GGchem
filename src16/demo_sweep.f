***********************************************************************
      SUBROUTINE DEMO_SWEEP
***********************************************************************
      use PARAMETERS,ONLY: Tmin,Tmax,pmin,pmax,nHmin,nHmax,
     >                     model_eqcond,model_pconst,Npoints
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,el,charge
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,muH,
     >                    amu,dust_nam,dust_mass,dust_Vol
      use EXCHANGE,ONLY: nel,nat,nion,nmol,
     >                   H,He,Li,C,N,O,F,Ne,Na,Mg,Al,Si,S,Cl,K,Ca,Ti,
     >                   Cr,Mn,Fe,Ni
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real :: p,pe,Tg,rho,nHges,nges,kT,pges
      real :: nTEA,pTEA,mu,muold,fac
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      integer :: i,j,jj,l,iel,NOUT
      character(len=5000) :: species,NISTspecies,elnames
      character(len=200) :: line
      character(len=20) :: frmt,name,short_name(NDUST),test1,test2
      character(len=4) :: sup
      character(len=2) :: test3
      character(len=1) :: char
      integer :: verbose=0
      logical :: isOK

      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
      enddo
      eps  = eps0
      NOUT = NELM
      if (charge) NOUT=NOUT-1
      open(unit=70,file='Static_Conc.dat',status='replace')
      write(70,1000) 'H',eps( H), 'C',eps( C),
     &               'N',eps( N), 'O',eps( O)
      write(70,*) NOUT,NMOLE,NDUST,Npoints
      write(70,2000) 'Tg','nHges','pges','el',
     &               (trim(elnam(elnum(j))),j=1,el-1),
     &               (trim(elnam(elnum(j))),j=el+1,NELM),
     &               (trim(cmol(i)),i=1,NMOLE),
     &               ('S'//trim(short_name(i)),i=1,NDUST),
     &               ('n'//trim(short_name(i)),i=1,NDUST)

      !--- TEA automated choice from dispol_large.dat ---
      species = "H_g He_ref C_g N_g O_g Si_g S_g Na_g "
     &        //"Ca_g Cl_g Ti_g K_g Al_g Mg_g Fe_g Li_g "
      elnames = " H He C N O Si S Na Ca Cl Ti K Al Mg Fe Li " 
      do i=1,NMOLE
        name = trim(cmol(i))
        isOK = .true.
        if (index(name,"+")>0) isOK=.false.
        if (index(name,"-")>0) isOK=.false.
        !print*,trim(name),isOK
        jj   = len(trim(name))
        j    = 1
        do while(j<=jj) 
          test1 = name(j:j)
          test2 = name(j:j+1)
          l = iachar(name(j:j))
          if (l<iachar("0").or.l>iachar("9")) then
            test3 = ''
            if (trim(test2)=='LI') test3='Li'
            if (trim(test2)=='MG') test3='Mg'
            if (trim(test2)=='AL') test3='Al'
            if (trim(test2)=='SI') test3='Si'
            if (trim(test2)=='CL') test3='Cl'
            if (trim(test2)=='CA') test3='Ca'
            if (trim(test2)=='TI') test3='Ti'
            if (trim(test2)=='FE') test3='Fe'
            if (trim(test2)=='HE') test3='He'
            if (test3.ne.'') then
              name(j:j+1) = test3(1:2) 
              test1 = test3
              j = j+1
            endif
            if (index(elnames," "//trim(test1)//" ")<=0) isOK=.false.
            !print*,trim(name),j,trim(test1),isOK
          endif  
          j = j+1
        enddo  
        if (trim(name)=="C3H") isOK=.false.  
        if (isOK) then
          sup = "_g"
          if (trim(name)=="HO2")   name="HOO" 
          if (trim(name)=="C2N2")  name="(CN)2" 
          if (trim(name)=="CHNO")  name="HNCO" 
          if (trim(name)=="HNO3")  name="HONO2"
          if (trim(name)=="N2H2")  name="HNNH"
          if (trim(name)=="CNO")   name="NC101"
          if (trim(name)=="H2SO4") name="O2S(OH)2"
          if (trim(name)=="SN")    name="NS"
          if (trim(name)=="S2O")   name="SSO"
          if (trim(name)=="H2") sup="_ref"
          if (trim(name)=="O2") sup="_ref"
          if (trim(name)=="N2") sup="_ref"
          species=trim(species)//" "//trim(name)//trim(sup)
        endif  
      enddo

      !--- TEA explicit choice ---  
      species = "H_g He_ref H2_ref C_g N_g O_g Si_g S_g "
     &  //"Na_g Ca_g Cl_g Ti_g K_g Al_g Mg_g Fe_g Li_g "
     &  //"H2O_g O2_ref CO_g CO2_g CH4_g C2H2_g N2_ref NH3_g " 
     &  //"OH_g CH_g CN_g HCO_g HCN_g NH_g NO_g NH2_g C2_g CH2_g CH3_g " 
     &  //"C2H_g HNO_g (CN)2_g C2H4_g C3_g C2N_g C2O_g "
     &  //"SiO_g SiH4_g SiS_g SiO2_g SiH_g SiN_g SiC_g Si2C_g SiC2_g "
     &  //"Si2N_g Si2_g Si3_g Si(CH3)4_g SiH3Cl_g SiH2Cl2_g "            
     &  //"H2S_g HS_g S2_g COS_g SO_g SO2_g SSO_g NS_g CS_g CS2_g "
     &  //"O2S(OH)2_g "
     &  //"Na2_g NaO_g NaOH_g (NaOH)2_g NaH_g NaCN_g Na2SO4_g "
     &  //"HCl_g NaCl_g (NaCl)2_g CaCl_g CaCl2_g SiCl_g "
     &  //"CaOH_g Ca(OH)2_g CaS_g CaO_g Ca2_g "                          ! CaH missing
     &  //"TiO_g TiO2_g TiCl_g TiCl2_g TiCl3_g TiCl4_g OTiCl_g "         
     &  //"TiOCl2_g "                                                    ! TiC,TiC2,TiH,TiS missing
     &  //"KCl_g KOH_g KH_g (KOH)2_g KCN_g (KCl)2_g K2SO4_g KO_g K2_g "
     &  //"AlOH_g OAlOH_g Al2O_g (AlO)2_g AlCl_g AlH_g AlS_g "
     &  //"AlCl2_g AlCl3_g (AlCl3)2_g AlO_g AlO2_g OAlH_g AlC_g "
     &  //"OAlCl_g AlN_g "
     &  //"Mg(OH)2_g MgOH_g MgH_g MgCl2_g MgS_g MgCl_g MgO_g MgN_g "
     &  //"(MgCl2)2_g Fe(OH)2_g FeCl_g FeCl2_g FeO_g FeS_g "             ! FeH missing
     &  //"LiCl_g LiOH_g LiH_g LiO_g (LiOH)2_g (LiCl)2_g (LiCl)3_g "  
     &  //"LiN_g "    

      !--- TEA complete ---
      NISTspecies = trim(species)
     &  //"Li2_g Li2O_g Li2SO4_g LiN_g (LiO)2_g LiOCl_g LiONa_g LiON_g " ! new Li
     &  //"(FeCl2)2_g (FeCl3)2_g FeCl3_g Fe(CO)5_g "                     ! new Fe
     &  //"Mg2_g "                                                       ! new Mg
     &  //"Al2_g "                                                       ! new Al
     &  //"(NaCN)2_g "                                                   ! new Na
     &  //"(KCN)2_g "                                                    ! new K
     &  //"SiCH3Cl3_g SiCl2_g SiCl3_g SiCl4_g "                          ! new Si
     &  //"ClSSCl_g S2Cl_g S3_g S4_g S5_g S6_g S7_g S8_g SCl2_g SCl_g "  ! new S
     &  //"SO2Cl2_g SO3_g "
     &  //"C2Cl2_g C2Cl4_g C2Cl6_g C2HCl_g CCl2_g CCl3_g CCl4_g CCl_g "  ! new Cl
     &  //"CH2Cl2_g CH3Cl_g CHCl3_g CHCl_g Cl2O_g ClCN_g ClO2_g ClO_g "
     &  //"COCl2_g COCl_g HOCl_g NO2Cl_g ONCl_g "
     &  //"C2H4O_g C3O2_g H2CO_g HNCO_g HONO2_g HOO_g HOOH_g N2O3_g "    ! new O
     &  //"N2O4_g N2O5_g N2O_g NO2_g NO3_g "                             
     &  //"C4N2_g CNN_g HNNH_g N2H4_g N3_g NC101_g NCN_g "               ! new N
     &  //"C4_g C5_g "                                                   ! new C

      !species = trim(NISTspecies)     ! want full NIST species list?

      write(frmt,'("(A",I4.4,")")') len(trim(species))
      open(unit=71,file='ggchem.atm',status='replace')
      write(71,'(A8)') "#SPECIES"
      write(71,frmt) trim(species) 
      write(71,*)
      write(71,'(A8)') "#TEADATA"
      write(71,'(A10,A8,99(A18))') 
     &      "#Pressure ","Temp","H","He","C","N","O","Si","S",
     &      "Na","Ca","Cl","Ti","K","Al","Mg","Fe","Li"

      !----------------------------------------
      ! ***  run chemistry on linear track  ***
      !----------------------------------------
      mu = muH
      do i=1,Npoints
        fac = REAL(i-1)/REAL(Npoints-1) 
        Tg  = EXP(LOG(Tmax)+fac*LOG(Tmin/Tmax))
        if (model_pconst) then
          p = EXP(LOG(pmax)+fac*LOG(pmin/pmax))
        else  
          nHges = EXP(LOG(nHmax)+fac*LOG(nHmin/nHmax))
        endif  
        eldust = 0.0
        do 
          if (model_pconst) nHges = p*mu/(bk*Tg)/muH
          if (model_eqcond) then
            call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
          endif  
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
          if (.not.model_pconst) exit
          print '("mu=",2(1pE12.5))',muold/amu,mu/amu
          if (ABS(mu/muold-1.0)<1.E-5) exit
        enddo  

        print'(i4," Tg[K] =",0pF8.2,"  n<H>[cm-3] =",1pE10.3)',
     >        i,Tg,nHges

        write(*,1010) ' Tg=',Tg,' pe=',nel*kT,' n<H>=',nHges,
     &                ' p=',pges/bar,' mu=',mu/amu
        print*
        write(70,2010) Tg,nHges,pges,
     &       LOG10(MAX(1.Q-300, nel)),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=1,el-1),
     &      (LOG10(MAX(1.Q-300, nat(elnum(jj)))),jj=el+1,NELM),
     &      (LOG10(MAX(1.Q-300, nmol(jj))),jj=1,NMOLE),
     &    (LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj)))),jj=1,NDUST),
     &             (LOG10(MAX(1.Q-300, eldust(jj))),jj=1,NDUST)

        nTEA = 0.0                           ! no electrons
        do j=1,NELEM
          nTEA = nTEA + nat(j)               ! atoms
        enddo
        do j=1,NMOLE
          if (index(cmol(j),"+")>0) cycle    ! no ions
          if (index(cmol(j),"-")>0) cycle    ! no cations
          nTEA = nTEA + nmol(j)              ! molecules
        enddo
        pTEA = nTEA*bk*Tg
        write(71,'(1pE10.4,0pF8.2,99(1pE18.11))') 
     &     pTEA/bar,Tg,eps(H),eps(He),eps(C),eps(N),eps(O),
     &     eps(Si),eps(S),eps(Na),eps(Ca),eps(Cl),eps(Ti),
     &     eps(K),eps(Al),eps(Mg),eps(Fe),eps(Li)

        if (verbose>0) read(*,'(a1)') char

      enddo  

      close(70)
      close(71)

      write(*,*)
      write(*,frmt) trim(species)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(a5,0pF8.2,99(a6,1pE10.3))
 2000 format(999(1x,a12))
 2010 format(0pF13.4,2(1pE13.4),999(0pF13.5))
      end  
