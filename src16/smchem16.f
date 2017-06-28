************************************************************************
      SUBROUTINE smchem16 (anHges,Tg,eps,anmono,anmol,merk,verbose)
************************************************************************
*                                                                      *
*     small chemistry                                                  *
*     ---------------                                                  *
*     Diese Routine berechnet eine GG-Chemie                           *
*                                                                      *
*   - Wird der Parameter "nachit"=FALSE gesetzt, so berechnet          *
*     sie kein exaktes Dissoziationsgleichgewicht, sondern nutzt       *
*     das Wissen aus dem Irsee-paper. Nur die wichtigen Dissozia-      *
*     tionsgleichgewichte werden dann geloest und so die Dichten       *
*     der haeufigsten Molekuele bestimmt. Dabei loest die Routine      *
*     zunaechst die Gleichgewichte der haeufigsten Elemente und        *
*     geht dann weiter zu weniger haeufigen. Die bereits geloesten     *
*     GG werden als von selteneren Elementen unbeeinflusst angenommen. *
*     Hierin liegt natuerlich auch die Schwaeche der Routine: Sie geht *
*     von bestimmten chemischen Haeufigkeiten aus. Weichen diese       *
*     stark von den Vorstellungen der Routine ab, so werden die        *
*     Ergebnisse falsch (z.B. wenn  [C]/[O] ~ 1 ist!!!)                *
*                                                                      *
*     Welche Molekuele in diesem Fall beruecksichtigt werden, steht    *
*     unten im Text.                                                   *
*                                                                      *
*     Die Abweichungen der wichtigen Molekuele liegen in der Regel     *
*     unter einem Faktor 2, oft unter 10%.                             *
*     Die Namen der beteiligten Atome und Molekuele werden ueber den   *
*     commonblock chemnam an die Aussenwelt uebergeben. Mit Hilfe      *
*     der Funktion stindex (siehe unten) kann man sich leicht nach     *
*     dem ersten Aufruf der Routine smchem die Indices der einzelnen   *
*     Molekuele im array anmol verschaffen; z.B.                       *
*     nnco = stindex(cmol,dim,'CO')                                    *
*                                                                      *
*     Die Routine ist extrem schnell und beansprucht wenig Platz.      *
*                                                                      *
*                                                                      *
*   - Wird der Parameter "nachit"=TRUE gesetzt, so berechnet die       *
*     Routine die GG-Chemie korrekt nach Newton-Raphson, wobei         *
*     obige Konzentrationen als Startwerte verwendet werden.           * 
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     e i n g a b e :                                                  *
*     anHges : totale Dichte der Wasserstoffatome (~ nh+2*nh2)         *
*     tg     : Gastemperatur                                           *
*     eps    : Vektor mit linearen chemischen Haeufigkeiten ([H]=1)    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     a u s g a b e :                                                  *
*     anmono : vektor mit den dichten der monomere  falls relevant     *
*     anmol  : vektor mit den dichten der molekuele falls relevant     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     b e n o e t i g t e  u n t e r p r o g r a m m e :               *
*     stindex                                                          *
*                                                                      *
*     b e n o e t i g t e  d a t e i e n :                             *
*     dispol.dat                                                       *
*                                                                      *
************************************************************************
*                                                                      *
*     Die Nummerierung der Elemente ist wie folgt:                     *
*     __________________________________________________________________ 
*     Element: He e- H C N O Si Mg Al Fe S  Na K  Ti Ca Li Cl Fl
*     Nummer:  1  2  3 4 5 6  7  8  9 10 11 12 13 14 15 16 17 18
*     __________________________________________________________________
*                                                                      *
************************************************************************
*     (c)       Carsten Dominik                     Do 11. Maerz 1993  *
*     Revision und nachit-Erweiterung Peter Woitke  Do  7. Maerz 1996  *
*     GG-Konstanten fuer TiC von Andreas Gauger gewonnen               *
*     auf der Grundlage spekrtoskopischer Messungen          11.09.97  *
************************************************************************
      use DUST_DATA,ONLY: NewChemIt,NewBackIt
      use CHEMISTRY,ONLY: nml=>NMOLE,nel=>NELM,cmol,catm,
     >                    m_kind,m_anz,a,natom,charge,elion,
     >                    He,el,H,C,N,O,Si,Mg,Al,Fe,S,Na,K,Ti,Ca,Li,Cl,
     >                    th1,th2,th3,th4,fit,TT1,TT2,TT3
      use EXCHANGE,ONLY: HII,CII,NII,OII,NaII,MgII,AlII,KII,TiII,SII,
     >                   SiII,FeII,CaII,LiII,ClII,HeII,chemcall,chemiter
      implicit none
*-----------------------------------------------------------------------
*  Dimensionierung fuer die Molekuel- und Atom Felder. Hier ist auf
*  Konsistenz mit dem aufrufenden Programm zu achten.
*-----------------------------------------------------------------------
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in)  :: anHges,Tg
      real(kind=qp),intent(in)  :: eps(nel)
      real(kind=qp),intent(out) :: anmono(nel),anmol(nml)
      integer,intent(inout) :: verbose
      real(kind=qp),parameter :: bk=1.380662d-16
*-----------------------------------------------------------------------
*  Die Variable "alle" entscheidet, ob die nicht unmittelbar 
*  beruecksichtigten Molekuele dennoch inkonsistent mitgerechnet werden. 
      logical alle
      data alle/.true./
*-----------------------------------------------------------------------
*  Die Variable "nachit" entscheidet, ob die GG-Chemie zur Ermittlung
*  der korrekten Loesung nach Newton-Raphson nachiteriert wird.
*  (braucht laenger, manchmal Konvergenzprobleme)
      logical nachit
      data nachit/.true./      
*-----------------------------------------------------------------------
*  Bei merk=.true. merkt sich die Routine die letzte konvergiert Loesung
*  und geht beim naechsten Mal von diesen Startwerten aus.
      logical,intent(INOUT) :: merk
*-----------------------------------------------------------------------
*  Die Variable "ngestst" entscheidet, ob die Elementerhaltung ueber-
*  prueft werden soll.
      logical ngestst
      data ngestst/.false./
*-----------------------------------------------------------------------
*  Die Variable tdispol bestimmt, welches die niedrigste Temperatur 
*  ist, die in die Dissoziationspolynome eingesetzt werden darf.
      real(kind=qp),parameter :: tdispol=100.Q0
*-----------------------------------------------------------------------
      integer stindex,at,info,ipvt(nel),Nconv,switch,iredo,ido
      integer Nact,all_to_act(nel),act_to_all(nel),switchoff(nel)
      integer e,i,j,j1,ii,jj,kk,l,it,m1,m2,piter,iraus,itry
      integer Nseq,imin,imax,enew,e2,eseq(nel)
      integer,parameter :: itmax=200
      real(kind=qp),parameter :: finish=1.Q-25
      real(kind=qp) :: ppp,qqq
      real(kind=qp) :: g(0:nml),limit,condnum1,work(nel)
      real(kind=qp) :: kT,kT1,nelek,ng,Sa,Nenner,fak,lth,arg,term
      real(kind=qp) :: f0,f1,f2,f3,ff1,ff2,ff3,f,fs
      real(kind=qp) :: VIETA,VIETA2
      real(kind=qp) :: pH,pHe,pC,pN,pO,pSi,pMg,pAl,pFe,pS,pNa,pK,pTi,pCa
      real(kind=qp) :: pLi,pCl,pel,worst
      real(kind=qp) :: pHges,pHeges,pCges,pNges,pOges,pSiges,pMgges,
     &                 pAlges,pFeges,pSges,pNages,pKges,pTiges,pCages,
     &                 pLiges,pClges,pHalt,pCalt,pOalt,pNalt,
     &                 pNaalt,pCaalt,pClalt,pKalt,pTialt,pSialt,pSalt
      real(kind=qp) :: aa,bb,cc,dd,ee,gg,hh,a2,a3,delta,pat,dpat
      real(kind=qp) :: nges(nel),pmono1(nel),coeff(-1:12),atmax,atfrac
      real(kind=qp) :: DF(nel,nel),dp(nel),FF(nel),pmol,q0,qq,crit
      real(kind=qp) :: DF0(nel,nel),FF0(nel),scale(nel),conv(0:500,nel)
      real(kind=qp) :: converge(0:500),delp,nold,soll,haben,abw,sum
      real(kind=qp) :: dpp(4),func(4),dfunc(4,4),sca(4)
      real(kind=qp) :: pbefore1(4),pbefore2(4),pbefore3(2),pbefore4(2)
      real(kind=qp) :: pbefore(nel)
      real(kind=qp) :: emax,pges,pwork
      logical :: from_merk,eact(nel),redo(nel),done(nel),affect,known
      character(len=100) :: txt,line
      character(len=1) :: char
      integer,save :: ilauf=0 
      integer,save :: Al2O,AlH,AlO2H,AlOH,C2,C3,C2H,C3H,C2H2,CH4,CN,CO
      integer,save :: CO2,CS,NO,FeS,H2,H2O,H2S,HCN,HS,MgH,MgO,MgOH,MgS
      integer,save :: N2,NH3,O2,OH,SO,SO2,Si2C,SiC,SiC2,SiH,SiH4,SiN
      integer,save :: SiO,SiO2,SiS,SiC4H12,FeO,FeO2H2,TiO,TiO2,TiS,TiC
      integer,save :: TiC2,TiCl,TiCl2,TiCl4,MGO2H2,NAOH,NA2O2H2,CaOH
      integer,save :: CaO2H2,KOH,K2O2H2,CKN,C2K2N2,FeH,CaH,TiH,LiH
      integer,save :: LiO,LiOH,LiCl,Li2Cl2,Li3Cl3,Li2O2H2,KCl,CaCl2
      integer,save :: CaCl,NaCl,HCl,NaH,H2SO4,MgCl,MgCl2,FeCl,FeCl2
      integer,save :: AlCl,AlCl2,AlCl3,Na2CL2,K2CL2,TiOCl2
      real(kind=qp),allocatable,save :: amerk(:),ansave(:)
      real(kind=qp),allocatable,save :: badness(:),pcorr(:,:) 
      real(kind=qp),save :: pcorr1(4),pcorr2(4),pcorr3(2),pcorr4(2)
*-----------------------------------------------------------------------
*  Die Formelfunktion zur Loesung quadratische Gleichungen mit Vieta
      VIETA(ppp,qqq)  = qqq/(-ppp/2.Q0-SQRT(ppp**2/4.Q0-qqq))
      VIETA2(ppp,qqq) = -ppp/2.Q0+SQRT(ppp**2/4.Q0-qqq)
*-----------------------------------------------------------------------      

      ilauf = ilauf + 1
      if ( ilauf .eq. 1 ) then
        TiC = stindex(cmol,nml,'TIC    ')
        if (.not.NewChemIt) then
          Al2O   = stindex(cmol,nml,'AL2O   ')
          AlH    = stindex(cmol,nml,'ALH    ')
          AlO2H  = stindex(cmol,nml,'ALO2H  ')
          AlOH   = stindex(cmol,nml,'ALOH   ')
          C2     = stindex(cmol,nml,'C2     ')
          C3     = stindex(cmol,nml,'C3     ')
          C2H    = stindex(cmol,nml,'C2H    ')
          C3H    = stindex(cmol,nml,'C3H    ')
          C2H2   = stindex(cmol,nml,'C2H2   ')
          CH4    = stindex(cmol,nml,'CH4    ')
          CN     = stindex(cmol,nml,'CN     ')
          CO     = stindex(cmol,nml,'CO     ')
          CO2    = stindex(cmol,nml,'CO2    ')
          CS     = stindex(cmol,nml,'CS     ')
          FeO    = stindex(cmol,nml,'FEO    ')
          FeO2H2 = stindex(cmol,nml,'FE(OH)2')
          FeS    = stindex(cmol,nml,'FES    ')
          H2     = stindex(cmol,nml,'H2     ')
          H2O    = stindex(cmol,nml,'H2O    ')
          H2S    = stindex(cmol,nml,'H2S    ')
          HCN    = stindex(cmol,nml,'HCN    ')
          HS     = stindex(cmol,nml,'HS     ')
          MgH    = stindex(cmol,nml,'MGH    ')
          MgO    = stindex(cmol,nml,'MGO    ')
          MgOH   = stindex(cmol,nml,'MGOH   ')
          MgO2H2 = stindex(cmol,nml,'MG(OH)2')
          MgS    = stindex(cmol,nml,'MGS    ')
          N2     = stindex(cmol,nml,'N2     ')
          NO     = stindex(cmol,nml,'NO     ')
          NH3    = stindex(cmol,nml,'NH3    ')
          O2     = stindex(cmol,nml,'O2     ')
          OH     = stindex(cmol,nml,'OH     ')
          SO     = stindex(cmol,nml,'SO     ')
          SO2    = stindex(cmol,nml,'SO2    ')
          Si2C   = stindex(cmol,nml,'SI2C   ')
          SiC    = stindex(cmol,nml,'SIC    ')
          SiC2   = stindex(cmol,nml,'SIC2   ')
          SiH    = stindex(cmol,nml,'SIH    ')
          SiH4   = stindex(cmol,nml,'SIH4   ')
          SiN    = stindex(cmol,nml,'SIN    ')
          SiO    = stindex(cmol,nml,'SIO    ')
          SiO2   = stindex(cmol,nml,'SIO2   ')
          SiS    = stindex(cmol,nml,'SIS    ')
          TiO    = stindex(cmol,nml,'TIO    ')
          TiO2   = stindex(cmol,nml,'TIO2   ')
          TiS    = stindex(cmol,nml,'TIS    ')
          TiC2   = stindex(cmol,nml,'TIC2   ')
          TiCl   = stindex(cmol,nml,'TICL   ')
          TiCl2  = stindex(cmol,nml,'TICL2  ')
          TiCl4  = stindex(cmol,nml,'TICL4  ')
          NaH    = stindex(cmol,nml,'NAH    ')
          NaOH   = stindex(cmol,nml,'NAOH   ')
          Na2O2H2= stindex(cmol,nml,'NA2O2H2')
          CaOH   = stindex(cmol,nml,'CAOH   ')
          CaO2H2 = stindex(cmol,nml,'CA(OH)2')
          KOH    = stindex(cmol,nml,'KOH    ')
          K2O2H2 = stindex(cmol,nml,'K2O2H2 ')
          C2K2N2 = stindex(cmol,nml,'C2K2N2 ')
          CKN    = stindex(cmol,nml,'CKN    ')
          HeII   = stindex(cmol,nml,'HE+    ')
          HII    = stindex(cmol,nml,'H+     ')
          CII    = stindex(cmol,nml,'C+     ')
          NII    = stindex(cmol,nml,'N+     ')
          OII    = stindex(cmol,nml,'O+     ')
          NaII   = stindex(cmol,nml,'NA+    ')
          MgII   = stindex(cmol,nml,'MG+    ')
          AlII   = stindex(cmol,nml,'AL+    ')
          KII    = stindex(cmol,nml,'K+     ')
          TiII   = stindex(cmol,nml,'TI+    ')
          SII    = stindex(cmol,nml,'S+     ')
          SiII   = stindex(cmol,nml,'SI+    ')
          FeII   = stindex(cmol,nml,'FE+    ')
          CaII   = stindex(cmol,nml,'CA+    ')
          LiII   = stindex(cmol,nml,'LI+    ')
          ClII   = stindex(cmol,nml,'CL+    ')
          FeH    = stindex(cmol,nml,'FEH    ')
          CaH    = stindex(cmol,nml,'CAH    ')
          TiH    = stindex(cmol,nml,'TIH    ')
          LiH    = stindex(cmol,nml,'LIH    ')
          LiO    = stindex(cmol,nml,'LIO    ')
          LiOH   = stindex(cmol,nml,'LIOH   ')
          LiCl   = stindex(cmol,nml,'LICL   ')
          Li2Cl2 = stindex(cmol,nml,'LI2CL2 ')
          Li3Cl3 = stindex(cmol,nml,'LI3CL3 ')
          Li2O2H2= stindex(cmol,nml,'LI2O2H2')
          KCl    = stindex(cmol,nml,'KCL    ')
          CaCl2  = stindex(cmol,nml,'CACL2  ')
          CaCl   = stindex(cmol,nml,'CACL   ')
          NaCl   = stindex(cmol,nml,'NACL   ')
          HCl    = stindex(cmol,nml,'HCL    ')
          H2SO4  = stindex(cmol,nml,'H2SO4  ')
          MgCl   = stindex(cmol,nml,'MGCL   ')
          MgCl2  = stindex(cmol,nml,'MGCL2  ')
          FeCl   = stindex(cmol,nml,'FECL   ')
          FeCl2  = stindex(cmol,nml,'FECL2  ')
          AlCl   = stindex(cmol,nml,'ALCL   ')
          AlCl2  = stindex(cmol,nml,'ALCL2  ')
          AlCl3  = stindex(cmol,nml,'ALCL3  ')
          Na2Cl2 = stindex(cmol,nml,'NA2CL2 ')
          K2Cl2  = stindex(cmol,nml,'K2CL2  ')
          TiOCL2 = stindex(cmol,nml,'TIOCL2 ')
          SiC4H12= stindex(cmol,nml,'SI(CH3)4')
        endif  
        allocate(badness(nel),pcorr(nel,nel),amerk(nel),ansave(nel))
        badness = 1.d0
        pcorr1  = 1.d0
        pcorr2  = 1.d0
        pcorr3  = 1.d0
        pcorr4  = 1.d0
        pcorr   = 1.d0
      endif
*-----------------------------------------------------------------------
*     ! zu niedrige Temperaturen abfangen und
*     ! Variable fuer die Dissoziationskonstanten berechnen
*     =====================================================
      TT1 = MAX(tdispol,Tg)
      TT2 = TT1*TT1
      TT3 = TT2*TT1
      th1 = 5040.Q0/TT1
      th2 = th1*th1
      th3 = th2*th1
      th4 = th3*th1
      kT  = bk*TT1
      kT1 = 1.Q0/kT
*      
*-----------------------------------------------------------------------
*     ! Vektoren initialisieren
*     =========================
      do j = 1 , nel 
        anmono(j) = 0.Q0
      enddo
      do j = 1 , nml 
        anmol(j)  = 0.Q0
        g(j)      = 0.Q0
      enddo
*
* --------------------------------------------------------------------------
*    TiC Gleichgewichtskonstante von Andreas Gauger ist anders
*        definiert als die Gleichgewichtskonstanten von Gail
*  Gauger: 
*  log Kp = 12.75293-5.4485*th1-1.56672*log(th1)+1.56041*(log(th1))**2
*           - 0.93275(log(th1))**3
*         = log ( p(A)p(B)/p(AB) )
*  Gail:
*   ln Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         =  ln ( p(AB)/p(A)p(B) )
*  Tsuji:
*  log Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         = log ( p(A)p(B)/p(AB) )
*
*  Umrechnung der Gauger-TiC-GG-Konstante in das Gail'sche System
*  -log(10)*log Kp(Gauger) = -2.30256*log Kp(Gauger) = ln Kp(Gail)

        lth = LOG10(th1)
        arg = 12.75293 - 5.44850*th1    - 1.56672*lth
     &                 + 1.56041*lth**2 - 0.93275*lth**3
        g(TiC) = EXP(MIN(1.Q+4,-2.30256*arg))

*  Umrechnen der Tsuji-CaH-GG-Konstante in das Gail'sche System
*  -log(10)*log Kp(Tsuji) = -2.30256*log Kp(Tsuji) = ln Kp(Gail)

      !   arg = 1.13401E+01 
      !&       -3.01442E+00*th1 
      !&       +4.23487E-01*th2
      !&       -6.14674E-02*th3
      !&       +3.16392E-03*th4
      !  g(CaH) = EXP(MIN(1.Q+4,-2.30256*arg))
c
      !   arg = 0.0
      !   arg =  1.17824E+01
      !&        -3.44431E+00*th1  ! this term differs between Tsuji 1073 paper and molecBP data file!  
      !&        +3.27412E-01*th2 
      !&        -4.84559E-02*th3  
      !&        +2.53356E-03*th4
      !   g(LiH) =  MIN(1.Q+300, EXP(-2.30256*arg))
c
      !   arg = 0.0
      !   arg =  1.23622E+01
      !&        -4.54966E+00*th1
      !&        +3.40687E-01*th2
      !&        -5.00589E-02*th3
      !&        +2.60132E-03*th4
      !   g(LiO) =  MIN(1.Q+300, EXP(-2.30256*arg))
c
      !   arg = 0.0
      !   arg =  1.24491E+01
      !&        -6.82101E+00*th1
      !&        +2.85775E-01*th2
      !&        -4.16715E-02*th3
      !&        +2.15772E-03*th4
      !   g(LiF) =  MIN(1.Q+300, EXP(-2.30256*arg))
c fit to Janaf tables made 27/02/2017
      !   arg  = 0.0
      !   arg  = 1.20145E+01
      !&        -5.66908E+00*th1
      !&        +2.55583E-01*th2
      !&        -3.69367E-02*th3
      !&        +1.90539E-03*th4
      !   g(LiCl) =  MIN(1.Q+300, EXP(-2.30256*arg))
c  fit to Janaf tables made 27/02/2017  
      !   arg  = 0.0
      !   arg  = 2.51219E+01
      !&        -1.06248E+01*th1
      !&        +5.03575E-01*th2
      !&        -7.21409E-02*th3
      !&        +3.71830E-03*th4
      !   g(LiOH) =  MIN(1.Q+300, EXP(-2.30256*arg))

*  Umrechnen der Burrow-FeH-deltaG polynome in das Gail'sche System
*  lnKp = - DeltaG_Burrows/RcalT - lnPstd
*  Achtung: [DeltaG_Burrows]= cal/mol ; [Rcal] = 1.987 cal/(mol K)
*  Pstd = 1atm
      !RcalT = Rcal*TT1
      !  arg = (-3.87740E+04 
      !&        +2.47290E+01*TT1
      !&        -4.64016E-04*TT2
      !&        +6.79410E-08*TT3)/RcalT
      !  g(FeH) = EXP(MIN(1.Q+4,-arg))/atm

c Sharp & Huebner 1990
      !  arg = (-3.04561E+05/TT1 
      !&        -4.91489E+04 
      !&        +2.32539E+01*TT1 
      !&        -1.56846E-04*TT2 
      !&        +4.53896E-08*TT3)/RcalT
      !  g(TiH) = EXP(MIN(1.Q+4,-arg))/atm
c fit to Janaf tables made 27/02/2017
      !   arg = ( 3.07543E+05/TT1 
      !&         -1.00087E+05 
      !&         +2.36342E+01*TT1 
      !&         +1.24759E-05*TT2 
      !&         +1.23522E-08*TT3)/RcalT
      !   g(NaCl) = MIN(1.Q+300, EXP(-arg)/atm)
c fit to Janaf tables made 27/02/2017
      !   arg = ( 4.57163e+05/TT1 
      !&         -1.04386E+05 
      !&         +2.38671E+01*TT1 
      !&         -3.71773E-04*TT2 
      !&         +6.22590E-08*TT3)/RcalT
      !   g(KCl) = MIN(1.Q+300, EXP(-arg)/atm)
c fit to Janaf tables made 27/02/2017
      !   arg = ( 2.89596e+05/TT1 
      !&         -9.88707E+04 
      !&         +2.13234E+01*TT1 
      !&         -1.27956E-04*TT2 
      !&         +4.49210E-08*TT3)/RcalT
      !   g(CaCl) = MIN(1.Q+300, EXP(-arg)/atm)
c fit to Janaf tables made 27/02/2017
      !   arg = ( 2.75428e+05/TT1 
      !&         -2.15484E+05 
      !&         +4.91645E+01*TT1 
      !&         -4.41153E-04*TT2 
      !&         +7.17853E-08*TT3)/RcalT
      !   g(CaCl2) = MIN(1.Q+300, EXP(-arg)/atm)
c fit to Janaf tables made 27/02/2017
      !   arg = ( 4.30684e+05/TT1 
      !&         -1.06291E+05 
      !&         +2.61097e+01*TT1 
      !&         +4.66915E-04*TT2 
      !&         -2.90088E-08*TT3)/RcalT
      !   g(HCl) = MIN(1.Q+300, EXP(-arg)/atm)

*---------------------------------------------------------------------------
      if ((ilauf.gt.10).and.merk) then
c       write(*,*) 'benutze Konzentrationen von vorher'
        do i=1,nel
          anmono(i) = amerk(i) * anhges
        enddo
        from_merk = .true.
        goto 200
      endif

*---------------------------------------------------------------------------
*     ! Elektronendichte
*     ================== 
 100  continue
      from_merk = .false.
      if (charge) then
        nelek = 0.Q0 
        do i=1,nel
          if (i==el) cycle 
          ng = anHges * eps(i)
          Sa = gk(elion(i))*kT1
          nelek = nelek + ng/(0.5Q0 + SQRT(0.25Q0 + ng/Sa))
        enddo
        anmono(el) = nelek
        pel = nelek*kT
      endif  

      if (.not.NewChemIt) then
*
*-----------------------------------------------------------------------
*     ! He: atomar
*     ============
      anmono(He) = anHges * eps(He) 
*
*     ! Gleichgewicht zwischen H, H+, H2
*     ==================================
      g(HII)    = gk( HII )
      g(H2)     = gk( H2 )
      pHges     = eps(H) * anHges * kT
      ppp       = ( 1.Q0 + g(HII)/pel )  / (2.Q0*g(H2))
      qqq       = -pHges / (2.Q0*g(H2))
      pH        = vieta(ppp,qqq)
*
*     ! Gleichgewicht  H, H+, H2, O, O+, O2, H2O
*     ==========================================
      g(OII) = gk( OII )
      g(O2)  = gk( O2  )
      g(H2O) = gk( H2O )
      pOges  = eps(O)*anHges*kT
      if (Tg<500.Q0) then
        pOges = max(pOges*0.5,pOges-(eps(Fe)+eps(Mg))*anHges*kT)
      endif  
      aa  = 2.Q0*g(H2)
      bb  = g(H2O)
      cc  = 1.Q0 + g(HII)/pel
      dd  = 1.Q0 + g(OII)/pel
      ee  = 2.Q0*g(O2)
      ppp = (dd + pH**2*bb) / ee
      qqq = -pOges / ee
      pO  = VIETA(ppp,qqq)
      ppp = cc / (aa + 2.Q0*pO*bb)
      qqq = -pHges / (aa + 2.Q0*pO*bb)
      pH  = VIETA(ppp,qqq)
      pbefore4(1:2) = (/pO,pH/)
      pO = pO * pcorr4(1)
      pH = pH * pcorr4(2)
      piter = 0
      sca(1) = pHges
      sca(2) = pOges
      do
        func(1)    = 2.Q0*pH**2*pO*bb + pH**2*aa + pH*cc - pHges
        dfunc(1,1) = 4.Q0*pH*pO*bb + 2.Q0*pH*aa + cc
        dfunc(1,2) = 2.Q0*pH**2*bb
        func(2)    = pH**2*pO*bb + pO**2*ee + pO*dd - pOges
        dfunc(2,1) = 2.Q0*pH*pO*bb
        dfunc(2,2) = pH**2*bb + 2.Q0*pO*ee + dd 
        dfunc(1:2,1) = dfunc(1:2,1)*sca(1)
        dfunc(1:2,2) = dfunc(1:2,2)*sca(2)
        call GAUSS16(4,2,dfunc,dpp,func)
        pHalt = pH
        pOalt = pO
        fak = 100.Q0
        pH  = MAX(MIN(pH-dpp(1)*sca(1),pH*fak),pH/fak)
        pO  = MAX(MIN(pO-dpp(2)*sca(2),pO*fak),pO/fak)
        delta = MAX(ABS(pHalt/pH-1.Q0),ABS(pOalt/pO-1.Q0))
        piter = piter + 1
        if (verbose>1) write(*,'(a11,i3,3(1pE11.4))') 
     &     'pH/pO-iter:',piter,pHalt,pOalt,delta
        if ((piter>99).or.(delta<1.Q-3)) exit
      enddo
      if (delta>1.Q-3) print*,"*** no conv. pH/pO",delta
      pcorr4(1:2) = (/pO,pH/) / pbefore4(1:2)  
*              
*     ! Gleichgewicht  O, O+, C, C+, CO, CO2, H2O, C2H2, CH4
*     ======================================================
      g(OII) = gk( OII )
      g(CII) = gk( CII )
      g(CH4) = gk( CH4 )
      g(C2H2)= gk( C2H2)
      g(H2O) = gk( H2O )
      g(CO)  = gk( CO  )
      g(CO2) = gk( CO2 )
      pCges  = eps(C) * anHges * kT
      aa     = 1.Q0 + g(OII)/pel + pH**2*g(H2O)
      bb     = 1.Q0 + g(CII)/pel + pH**4*g(CH4)
      cc     = 2.Q0*pH**2*g(C2H2)
      ppp    = (bb+pO*g(CO)+pO**2*g(CO2))/cc
      qqq    = -pCges/cc
      pC     = VIETA(ppp,qqq)
      sca(1) = pOges
      sca(2) = pCges
      pbefore3(1:2) = (/pC,pO/)
      pC = pC * pcorr3(1)
      pO = pO * pcorr3(2)
      piter = 0
      do 
        func(1)    = 2.Q0*pO**2*pC*g(CO2) + pO*(pC*g(CO)+aa) - pOges 
        dfunc(1,1) = 4.Q0*pO*pC*g(CO2)    + (pC*g(CO)+aa)
        dfunc(1,2) = 2.Q0*pO**2*g(CO2)    + pO*g(CO)
        func(2)    = pC**2*cc + pC*(pO**2*g(CO2)+pO*g(CO)+bb) - pCges
        dfunc(2,1) = pC*(2.Q0*pO*g(CO2)+g(CO))
        dfunc(2,2) = 2.Q0*pC*cc +  (pO**2*g(CO2)+pO*g(CO)+bb) 
        dfunc(1:2,1) = dfunc(1:2,1)*sca(1)
        dfunc(1:2,2) = dfunc(1:2,2)*sca(2)
        call GAUSS16(4,2,dfunc,dpp,func)
        pCalt = pC
        pOalt = pO
        fak = 100.Q0
        pO  = MAX(MIN(pO-dpp(1)*sca(1),pO*fak),pO/fak)
        pC  = MAX(MIN(pC-dpp(2)*sca(2),pC*fak),pC/fak)        
        delta = MAX(ABS(pCalt/pC-1.Q0),ABS(pOalt/pO-1.Q0))
        piter = piter+1
        if (verbose>1) write(*,'(a11,i3,3(1pE11.4))') 
     &       'pC/pO-iter:',piter,pCalt,pOalt,delta
        if ((piter>99).or.(delta<1.Q-3)) exit
      enddo  
      if (delta>1.Q-3) print*,"*** no conv. pC/pO",delta
      pcorr3(1:2) = (/pC,pO/) / pbefore3(1:2)  
      !write(98,'(0pF8.2,I4,4(1pE11.3))') Tg,piter,pcorr3

*     ! Gleichgewicht N, N+, N2, CN, HCN, NH3, NO
*     ===========================================
      g(NII)     = gk(NII)
      g(N2)      = gk(N2)
      g(CN)      = gk(CN)
      g(NO)      = gk(NO)
      g(HCN)     = gk(HCN)
      g(NH3)     = gk(NH3)
      pNges      = eps(N) * anHges * kT
      ppp        = (1.Q0 + g(CN)*pC + g(NO)*pO + g(HCN)*pC*pH 
     &              + g(NH3)*pH**3 + g(NII)/pel)/(2.Q0*g(N2))
      qqq        = -pNges / (2.Q0*g(N2))
      pN         = vieta(ppp,qqq)

*     ! pH/pC/pO/pN-iteration:  H,H+,O,O+,H2,H2O,O2,CO,CO2,
*     !                         C,C+,CH4,C2H2,C2H,C2,C3,C3H
*     !                         N,N+,N2,CN,HCN,NH3,NO
*     =====================================================
      g(C2)  = gk( C2  )
      g(C2H) = gk( C2H )
      g(C3)  = gk( C3  )
      g(C3H) = gk( C3H )
      aa = 1.Q0 + g(CII)/pel
      bb = 1.Q0 + g(OII)/pel
      cc = 1.Q0 + g(HII)/pel
      dd = 1.Q0 + g(NII)/pel
      sca(1) = pCges
      sca(2) = pOges
      sca(3) = pHges
      sca(4) = pNges
      pbefore1(1:4) = (/pC,pO,pH,pN/)
      pC = pC * pcorr1(1)
      pO = pO * pcorr1(2)
      pH = pH * pcorr1(3)
      pN = pN * pcorr1(4)
      piter = 0
      do
        !------------------------------------------------------------------------
        func(1)    = 3.Q0*pC**3*( g(C3) + pH*g(C3H) ) - pCges          ! carbon
     &             + 2.Q0*pC**2*( pH**2*g(C2H2) + pH*g(C2H) + g(C2) )
     &             + pC*( aa + pO*g(CO) + pO**2*g(CO2) + pH**4*g(CH4) 
     &                    + pN*g(CN) + pH*pN*g(HCN) )
        dfunc(1,1) = 9.Q0*pC**2*( g(C3) + pH*g(C3H) )
     &             + 4.Q0*pC*( pH**2*g(C2H2) + pH*g(C2H) + g(C2) )
     &             + ( aa + pO*g(CO) + pO**2*g(CO2) + pH**4*g(CH4) 
     &                    + pN*g(CN) + pH*pN*g(HCN) )
        dfunc(1,2) = pC*( g(CO) + 2.Q0*pO*g(CO2) )
        dfunc(1,3) = 3.Q0*pC**3*g(C3H)
     &             + 2.Q0*pC**2*( 2.Q0*pH*g(C2H2) + g(C2H) )
     &             + pC*( 4.Q0*pH**3*g(CH4) + pN*g(HCN) )
        dfunc(1,4) = pC*( g(CN) + pH*g(HCN) )
        !------------------------------------------------------------------------
        func(2)    = 2.Q0*pO**2*( g(O2) + pC*g(CO2) ) - pOges          ! oxygen
     &             + pO*( bb + pH**2*g(H2O) + pC*g(CO) + pN*g(NO) )               
        dfunc(2,1) = 2.Q0*pO**2*g(CO2)
     &             + pO*g(CO) 
        dfunc(2,2) = 4.Q0*pO*( g(O2) + pC*g(CO2) )
     &             + ( bb + pH**2*g(H2O) + pC*g(CO) + pN*g(NO) )
        dfunc(2,3) = pO*2.Q0*pH*g(H2O)
        dfunc(2,4) = pO*g(NO)
        !------------------------------------------------------------------------
        func(3)    = 4.Q0*pH**4*pC*g(CH4) - pHges                      ! hydrogen
     &             + 3.Q0*pH**3*pN*g(NH3)
     &             + 2.Q0*pH**2*( g(H2) + pO*g(H2O) + pC**2*g(C2H2) )
     &             + pH*( cc + pC**2*g(C2H) + pC**3*g(C3H) 
     &                       + pC*pN*g(HCN) )
        dfunc(3,1) = 4.Q0*pH**4*g(CH4)
     &             + 2.Q0*pH**2*( 2.Q0*pC*g(C2H2) )
     &             + pH*(2.Q0*pC*g(C2H) + 3.Q0*pC**2*g(C3H) + pN*g(HCN))
        dfunc(3,2) = 2.Q0*pH**2*g(H2O)
        dfunc(3,3) = 16.Q0*pH**3*pC*g(CH4)
     &             + 9.Q0*pH**2*pN*g(NH3)
     &             + 4.Q0*pH*( g(H2) + pO*g(H2O) + pC**2*g(C2H2) )
     &             + ( cc + pC**2*g(C2H) + pC**3*g(C3H) + pC*pN*g(HCN) )
        dfunc(3,4) = 3.Q0*pH**3*g(NH3)
     &             + pH*( pC*g(HCN) )
        !------------------------------------------------------------------------
        func(4)    = 2.Q0*pN**2*g(N2) - pNges                          ! nitrogen
     &             + pN*( dd + pC*g(CN) + pH*pC*g(HCN) + pH**3*g(NH3) 
     &                  + pO*g(NO) )
        dfunc(4,1) = pN*( g(CN) + pH*g(HCN) )
        dfunc(4,2) = pN*g(NO)
        dfunc(4,3) = pN*( pC*g(HCN) + 3.Q0*pH**2*g(NH3) )
        dfunc(4,4) = 4.Q0*pN*g(N2) 
     &             + ( dd + pC*g(CN) + pH*pC*g(HCN) + pH**3*g(NH3)
     &               + pO*g(NO) )
        dfunc(1:4,1) = dfunc(1:4,1)*sca(1)
        dfunc(1:4,2) = dfunc(1:4,2)*sca(2)
        dfunc(1:4,3) = dfunc(1:4,3)*sca(3)
        dfunc(1:4,4) = dfunc(1:4,4)*sca(4)

        call GAUSS16(4,4,dfunc,dpp,func)
        pCalt = pC
        pOalt = pO
        pHalt = pH
        pNalt = pN
        fak = 10.Q0
        pC  = MAX(MIN(pC-dpp(1)*sca(1),pC*fak),pC/fak)
        pO  = MAX(MIN(pO-dpp(2)*sca(2),pO*fak),pO/fak)
        pH  = MAX(MIN(pH-dpp(3)*sca(3),pH*fak),pH/fak)
        pN  = MAX(MIN(pN-dpp(4)*sca(4),pN*fak),pN/fak)
        delta = MAX(ABS(pCalt/pC-1.Q0),ABS(pOalt/pO-1.Q0),
     &              ABS(pHalt/pH-1.Q0),ABS(pNalt/pN-1.Q0))          
        piter = piter + 1
        if (verbose>1.or.piter>80) write(*,'(a15,i3,5(1pE11.4))') 
     &       'pC/pO/pH/pN-it:',piter,pCalt,pOalt,pHalt,pNalt,delta
        if ((piter>99).or.(delta<1.Q-3)) exit
      enddo  
      if (delta>1.Q-3) print*,"*** no conv. pC/pO/pH/pN",delta 
      pcorr1(1:4) = (/pC,pO,pH,pN/) / pbefore1(1:4)  
      !write(98,'(0pF8.2,I4,4(1pE11.3))') Tg,piter,pcorr1
      anmono(H)   = pH / kT
      anmono(C)   = pC / kT
      anmono(N)   = pN / kT
      anmono(O)   = pO / kT
      anmol(H2)   = g(H2)   * pH**2 / kT         
      anmol(H2O)  = g(H2O)  * pH**2 * pO / kT
      anmol(CO)   = g(CO)   * pC * pO / kT
      anmol(CO2)  = g(CO2)  * pC * pO**2 / kT
      anmol(O2)   = g(O2)   * pO**2 / kT
      anmol(C2)   = g(C2)   * pC**2 / kT
      anmol(C3)   = g(C3)   * pC**3 / kT
      anmol(C2H)  = g(C2H)  * pC**2 * pH / kT
      anmol(C3H)  = g(C3H)  * pC**3 * pH / kT
      anmol(C2H2) = g(C2H2) * pC**2 * pH**2 / kT
      anmol(CH4)  = g(CH4)  * pC * pH**4 / kT
      anmol(N2)   = g(N2)   * pN**2 / kT
      anmol(CN)   = g(CN)   * pC * pN / kT
      anmol(HCN)  = g(HCN)  * pH * pC * pN / kT
      anmol(NH3)  = g(NH3)  * pN * pH**3 / kT
*
*     ! Gleichgewicht Si,S,SiS,SiC,SiO,SiO2,Si2C,SiC2,SiH,SiH4,SiN,CS,HS,H2S,SO,SO2,H2SO4
*     ===================================================================================
      g(SiII)  = gk( SiII )
      g(SII)   = gk( SII )
      g(SiS)   = gk( SiS )
      g(SiC)   = gk( SiC )
      g(SiO)   = gk( SiO )
      g(SiO2)  = gk( SiO2 )
      g(Si2C)  = gk( Si2C )
      g(SiC2)  = gk( SiC2 )
      g(SiH)   = gk( SiH )
      g(SiH4)  = gk( SiH4 )
      g(SiN)   = gk( SiN )
      g(CS)    = gk( CS )
      g(HS)    = gk( HS )
      g(H2S)   = gk( H2S )
      g(SO)    = gk( SO )
      g(SO2)   = gk( SO2 )
      g(H2SO4) = gk( H2SO4 )
      g(SiC4H12) = gk( SiC4H12 )
      pSiges   = eps(Si) * anHges * kT
      pSges    = eps( S) * anHges * kT
      aa       = 1.Q0 + pC*g(SiC) + pC**2*g(SiC2) + pH*g(SiH) 
     &         + pH**4*g(SiH4) + pN*g(SiN) + pO*g(SiO) + pO**2*g(SiO2)
     &         + pC**4*pH**12*g(SiC4H12) + g(SiII)/pel
      bb       = 1.Q0 + pC*g(CS) + pH*g(HS) + pH**2*g(H2S) 
     &         + pO*g(SO) + pO**2*g(SO2) + pH**2*pO**4*g(H2SO4)
     &         + g(SII)/pel
      ppp      = aa/g(SiS) + (pSiges-pSges)/bb
      qqq      = -pSges * aa / bb / g(SiS)
      pS       = vieta(ppp,qqq)
      pSi      = pSiges / ( aa + pS * g(SiS) )
*
*     ! Nachiteration wegen Si2C
*     ==========================                
      sca(1) = pSiges
      sca(2) = pSges
      piter = 0
      do 
        func(1)    = 2.Q0*pSi**2*pC*g(Si2C) + pSi*(pS*g(SiS)+aa)   ! Si
     &             - pSiges
        dfunc(1,1) = 4.Q0*pSi*pC*g(Si2C) + (pS*g(SiS)+aa)
        dfunc(1,2) = pSi*g(SiS) 
        func(2)    = pS*(pSi*g(SiS)+bb) - pSges                    ! S
        dfunc(2,1) = pS*g(SiS)
        dfunc(2,2) = (pSi*g(SiS)+bb)
        dfunc(1:2,1) = dfunc(1:2,1)*sca(1)
        dfunc(1:2,2) = dfunc(1:2,2)*sca(2)
        call GAUSS16(4,2,dfunc,dpp,func)
        pSialt = pSi
        pSalt  = pS
        fak = 100.Q0
        pSi = MAX(MIN(pSi-dpp(1)*sca(1),pSi*fak),pSi/fak)
        pS  = MAX(MIN(pS -dpp(2)*sca(2),pS *fak),pS /fak)        
        delta = MAX(ABS(pSialt/pSi-1.Q0),ABS(pSalt/pS-1.Q0))
        piter = piter+1
        if (verbose>1) write(*,'(a11,i3,3(1pE11.4))') 
     &       'pSi/pS-iter:',piter,pSialt,pSalt,delta
        if ((piter>99).or.(delta<1.Q-3)) exit
      enddo  
      if (delta>1.Q-3) print*,"*** no conv. pSi/pS",delta
      anmono(Si)  = pSi / kT
      anmono(S)   = pS / kT
      anmol(Si2C) = g(Si2C) * pSi**2 * pC    / kT
      anmol(SiC2) = g(SiC2) * pSi    * pC**2 / kT
      anmol(SiH)  = g(SiH)  * pSi    * pH    / kT
      anmol(CS)   = g(CS)   * pS     * pC    / kT
      anmol(HS)   = g(HS)   * pS     * pH    / kT
      anmol(H2S)  = g(H2S)  * pS     * pH**2 / kT
      anmol(SiS)  = g(SiS)  * pSi    * pS    / kT
      anmol(SiC)  = g(SiC)  * pSi    * pC    / kT
      anmol(SiH4) = g(SiH4) * pSi    * pH**4 / kT
      anmol(SiO2) = g(SiO2) * pSi    * pO**2 / kT
      anmol(SiN)  = g(SiN)  * pSi    * pN    / kT
      anmol(SO)   = g(SO)   * pS     * pO / kT
      anmol(SO2)  = g(SO2)  * pS     * pO**2 / kT
      anmol(H2SO4)= g(H2SO4)* pS     * pH**2 * pO**4 / kT
*
*     ! Gleichgewicht Na, Na+, NaOH, (NaOH)2, NaH
*                     Ca, Ca+, CaH, CaOH, Ca(OH)2
*                     Cl, Cl+, HCl, NaCl, CaCl, CaCl2, Na2Cl2
*                     Ti, Ti+, TiC2, TiS, TiO, TiO2, TiCl, TiCl2, TiCl4, TiOCl2
*     =========================================================================
      g(NaII)    = gk( NaII )
      g(NaH)     = gk( NaH )
      g(NaOH)    = gk( NaOH )
      g(Na2O2H2) = gk( Na2O2H2 )
      g(CaH)     = gk( CaH )
      g(CaII)    = gk( CaII )
      g(CaOH)    = gk( CaOH )
      g(CaO2H2)  = gk( CaO2H2 )
      g(ClII)    = gk( ClII )
      g(NaCl)    = gk( NaCl )
      g(Na2Cl2)  = gk( Na2Cl2 )
      g(HCl)     = gk( HCl )
      g(CaCl)    = gk( CaCl )
      g(CaCl2)   = gk( CaCl2 )
c     g(TiC)   : siehe oben!
      g(TiH)     = gk( TiH )
      g(TiII)    = gk( TiII )        
      g(TiC2)    = gk( TiC2 )
      g(TiS)     = gk( TiS )
      g(TiO)     = gk( TiO )
      g(TiO2)    = gk( TiO2 )
      g(TiCl)    = gk( TiCl )
      g(TiCl2)   = gk( TiCl2 )
      g(TiCl4)   = gk( TiCl4 )
      g(TiOCL2)  = gk( TiOCL2 )
      pNages = eps(Na) * anHges * kT
      pCages = eps(Ca) * anHges * kT
      pClges = eps(Cl) * anHges * kT
      pTiges = eps(Ti) * anHges * kT    
      aa  = 1.Q0 + g(ClII)/pel + pH*g(HCl)
      bb  = 1.Q0 + g(NaII)/pel + pO*pH*g(NaOH) + pH*g(NaH) 
      ppp = bb/g(NaCl) + (pNages-pClges)/aa
      qqq = -pClges * bb / aa / g(NaCl)
      if (ppp>0.Q0) then
        pCl = vieta(ppp,qqq)
      else  
        pCl = vieta2(ppp,qqq)
      endif  
      aa  = 2.Q0*pO**2*pH**2*g(Na2O2H2)
      ppp = (bb+pCl*g(NaCl))/aa
      qqq = -pNages/aa
      pNa = vieta(ppp,qqq)
      pCa = pCages / ( 1.Q0 + g(CaII)/pel + pO*pH*g(CaOH) 
     &                 + pH*g(CaH) + pO**2*pH**2*g(CaO2H2) 
     &                 + pCl*g(CaCl) + pCl**2*g(CaCl2) ) 
      pTi = pTiges / ( 1.Q0 + g(TiII)/pel + pC**2*g(TiC2) + pS*g(TiS) 
     &               + pO*g(TiO) + pO**2*g(TiO2) + pCl*g(TiCl) 
     &               + pCl**2*g(TiCl2) + pCl**4*g(TiCl4)
     &               + pO*pCl**2*g(TiOCl2) ) 
      aa  = 2.Q0*(pTi*pO*g(TiOCl2)+pNa**2*g(Na2Cl2)+pCa*g(CaCl2))
      ppp = (1.Q0+g(ClII)/pel+pNa*g(NaCl)+pH*g(HCl))/aa
      qqq = -pClges/aa
      pCl = vieta(ppp,qqq)      
      aa  = 1.Q0 + g(NaII)/pel + pO*pH*g(NaOH) + pH*g(NaH)
      bb  = 1.Q0 + g(CaII)/pel + pO*pH*g(CaOH) + pH*g(CaH) 
     &                         + pO**2*pH**2*g(CaO2H2) 
      cc  = 1.Q0 + g(ClII)/pel + pH*g(HCl)
      dd  = 1.Q0 + g(TiII)/pel + pC**2*g(TiC2) + pS*g(TiS) + pO*g(TiO)
     &                         + pO**2*g(TiO2)
      piter = 0
      pbefore2(1:4) = (/pNa,pCa,pCl,pTi/)
      pNa = pNa * pcorr2(1)
      pCa = pCa * pcorr2(2)
      pCl = pCl * pcorr2(3)
      pTi = pTi * pcorr2(4)
      sca(1) = pNages
      sca(2) = pCages
      sca(3) = pClges
      sca(4) = pTiges
      do
        !-------------------------------------------------------------------------
        func(1)    = 2.Q0*pNa**2*( pCl**2*g(Na2Cl2)                       ! sodium
     &                            +pO**2*pH**2*g(Na2O2H2) ) - pNages
     &             + pNa*( aa + pCl*g(NaCl) )
        dfunc(1,1) = 4.Q0*pNa*(pCl**2*g(Na2Cl2)+pO**2*pH**2*g(Na2O2H2))
     &             + ( aa + pCl*g(NaCl) )
        dfunc(1,2) = 0.Q0
        dfunc(1,3) = 4.Q0*pNa**2*pCl*g(Na2Cl2) + pNa*g(NaCl)
        dfunc(1,4) = 0.Q0
        !-------------------------------------------------------------------------
        func(2)    = pCa*( bb + pCl*g(CaCl) + pCl**2*g(CaCl2) ) - pCages ! calcium
        dfunc(2,1) = 0.Q0
        dfunc(2,2) = ( bb + pCl*g(CaCl) + pCl**2*g(CaCl2) )
        dfunc(2,3) = pCa*( g(CaCl) + 2.Q0*pCl*g(CaCl2) )
        dfunc(2,4) = 0.Q0
        !-------------------------------------------------------------------------
        func(3)    = 4.Q0*pCl**4*( pTi*g(TiCl4) ) - pClges               ! clorine
     &             + 2.Q0*pCl**2*( pTi*g(TiCl2) + pCa*g(CaCl2) 
     &                            +pNa**2*g(Na2Cl2) + pTi*pO*g(TiOCl2) )
     &             + pCl*( cc + pNa*g(NaCl) + pCa*g(CaCl) + pTi*g(TiCl)) 
        dfunc(3,1) = 4.Q0*pCl**2*pNa*g(Na2Cl2) + pCl*g(NaCl)
        dfunc(3,2) = 2.Q0*pCl**2*g(CaCl2) + pCl*g(CaCl)
        dfunc(3,3) = 16.Q0*pCl**3*( pTi*g(TiCl4) )
     &             + 4.Q0*pCl*( pTi*g(TiCl2) + pCa*g(CaCl2)
     &                         +pNa**2*g(Na2Cl2) + pTi*pO*g(TiOCl2) )
     &             + ( cc + pNa*g(NaCl) + pCa*g(CaCl) + pTi*g(TiCl) ) 
        dfunc(3,4) = 4.Q0*pCl**4*g(TiCl4)
     &             + 2.Q0*pCl**2*( g(TiCl2) + pO*g(TiOCl2) )
     &             + pCl*g(TiCl)  
        !-------------------------------------------------------------------------
        func(4)    = pTi*( dd + pCl*g(TiCl) + pCl**2*g(TiCl2)              ! titan
     &                   + pCl**4*g(TiCl4) + pO*pCl**2*g(TiOCl2) ) 
     &             - pTiges
        dfunc(4,1) = 0.Q0
        dfunc(4,2) = 0.Q0
        dfunc(4,3) = pTi*( g(TiCl) + 2.Q0*pCl*g(TiCl2) 
     &                   + 4.Q0*pCl**3*g(TiCl4) + pO*2.Q0*pCl*g(TiOCl2))
        dfunc(4,4) = ( dd + pCl*g(TiCl) + pCl**2*g(TiCl2) 
     &               + pCl**4*g(TiCl4) + pO*pCl**2*g(TiOCl2) )
        dfunc(1:4,1) = dfunc(1:4,1)*sca(1)
        dfunc(1:4,2) = dfunc(1:4,2)*sca(2)
        dfunc(1:4,3) = dfunc(1:4,3)*sca(3)
        dfunc(1:4,4) = dfunc(1:4,4)*sca(4)

        call GAUSS16(4,4,dfunc,dpp,func)
        pNaalt = pNa
        pCaalt = pCa
        pClalt = pCl
        pTialt = pTi
        fak = 5.0
        pNa = MAX(MIN(pNa-dpp(1)*sca(1),pNa*fak),pNa/fak)
        pCa = MAX(MIN(pCa-dpp(2)*sca(2),pCa*fak),pCa/fak)
        pCl = MAX(MIN(pCl-dpp(3)*sca(3),pCl*fak),pCl/fak)
        pTi = MAX(MIN(pTi-dpp(4)*sca(4),pTi*fak),pTi/fak)
        delta = MAX(ABS(pNaalt/pNa-1.Q0),ABS(pCaalt/pCa-1.Q0),
     &              ABS(pClalt/pCl-1.Q0),ABS(pTialt/pTi-1.Q0))
        piter = piter + 1
        if (verbose>1.or.piter>80) write(*,'(a19,i3,5(1pE11.4))') 
     &       'pNa/pCa/pCl/pTi-it:',piter,pNaalt,pCaalt,
     &                             pClalt,pTialt,delta
        if ((piter>99).or.(delta<1.Q-3)) exit
      enddo  
      if (delta>1.Q-3) print*,"*** no conv. pNa/pCa/pCl/pTi",delta
      pcorr2(1:4) = (/pNa,pCa,pCl,pTi/) / pbefore2(1:4)  
      !write(98,'(0pF8.2,I4,4(1pE11.3))') Tg,piter,pcorr2
      anmono(Na)     = pNa / kT
      anmono(Ca)     = pCa / kT
      anmono(Cl)     = pCl / kT 
      anmono(Ti)     = pTi / kT
      anmol(NaH)     = g(NaH)     * pNa * pH / kT
      anmol(NaOH)    = g(NaOH)    * pNa * pO * pH / kT
      anmol(Na2O2H2) = g(Na2O2H2) * pNa**2 * pO**2 * pH**2 / kT
      anmol(CaH)     = g(CaH)     * pCa * pH / kT
      anmol(CaOH)    = g(CaOH)    * pCa * pO * pH / kT
      anmol(CaO2H2)  = g(CaO2H2)  * pCa * pO**2 * pH**2 / kT
      anmol(NaCl)    = g(NaCl)    * pNa * pCl / kT
      anmol(HCl)     = g(HCl)     * pCl * pH / kT
      anmol(CaCl)    = g(CaCl)    * pCl * pCa / kT
      anmol(CaCl2)   = g(CaCl2)   * pCl**2 * pCa / kT
      anmol(TiH)     = g(TiH)     * pTi * pH / kT
      anmol(TiC)     = g(TiC)     * pTi * pC / kT
      anmol(TiC2)    = g(TiC2)    * pTi * pC**2 / kT
      anmol(TiS)     = g(TiS)     * pTi * pS / kT
      anmol(TiO)     = g(TiO)     * pTi * pO / kT 
      anmol(TiO2)    = g(TiO2)    * pTi * pO**2 / kT 
      anmol(TiCl)    = g(TiCl)    * pTi * pCl / kT
      anmol(TiCl2)   = g(TiCl2)   * pTi * pCl**2 / kT
      anmol(TiCl4)   = g(TiCl4)   * pTi * pCl**4 / kT
*        
*     ! Gleichgewicht Mg, MgH, MgO, MgOH, Mg(OH)2, MgS, Mg+, MgCl, MgCl2
*     ==================================================================
      g(MgII)   = gk( MgII )
      g(MgO)    = gk( MgO )
      g(MgOH)   = gk( MgOH )
      g(MgO2H2) = gk( MgO2H2 )
      g(MgH)    = gk( MgH )
      g(MgS)    = gk( MgS )
      g(MgCl)   = gk( MgCl )
      g(MgCl2)  = gk( MgCl2 )
      pMgges    = eps(Mg) * anHges * kT
      pMg       = pMgges / ( 1.Q0 + g(MgH)*pH + g(MgO)*pO + g(MgS)*pS 
     &          + g(MgOH)*pH*pO + g(MgO2H2)*pO**2*pH**2 
     &          + g(MgCl)*pCl + g(MgCl2)*pCl**2 + g(MgII)/pel )
      anmono(Mg)    = pMg / kT
      anmol(MgH)    = g(MgH)    * pMg * pH            / kT
      anmol(MgO)    = g(MgO)    * pMg * pO            / kT
      anmol(MgOH)   = g(MgOH)   * pMg * pO    * pH    / kT
      anmol(MgO2H2) = g(MgO2H2) * pMg * pO**2 * pH**2 / kT
      anmol(MgS)    = g(MgS)    * pMg * pS            / kT
      anmol(MgCl)   = g(MgCl)   * pMg * pCl           / kT
      anmol(MgCl2)  = g(MgCl2)  * pMg * pCl**2        / kT
*
*     ! Gleichgewicht  Fe, Fe+, FeO, FeS, FeH, Fe(OH)2, FeCl, FeCl2
*     =============================================================
      g(FeH)        = gk( FeH    )
      g(FeII)       = gk( FeII   )
      g(FeO)        = gk( FeO    )
      g(FeS)        = gk( FeS    )
      g(FeO2H2)     = gk( FeO2H2 )
      g(FeCl)       = gk( FeCl   )
      g(FeCl2)      = gk( FeCl2  )
      pFeges        = eps(Fe) * anHges * kT
      pFe           = pFeges / ( 1.Q0 + pO*g(FeO) + g(FeII)/pel 
     &                + pH**2*pO**2*g(FeO2H2) + pH*g(FeH) + pS*g(FeS)
     &                + pCl*g(FeCl) + pCl**2*g(FeCl2) )
      anmono(Fe)    = pFe / kT
      anmol(FeO)    = g(FeO)    * pFe * pO / kT 
      anmol(FeS)    = g(FeS)    * pFe * pS / kT 
      anmol(FeH)    = g(FeH)    * pFe * pH / kT 
      anmol(FeO2H2) = g(FeO2H2) * pFe * pO**2 * pH**2 / kT 
      anmol(FeCl)   = g(FeCl)   * pFe * pCl / kT 
      anmol(FeCl2)  = g(FeCl2)  * pFe * pCl**2 / kT 
*
*     ! Gleichgewicht Al, Al+, AlOH, AlO2H, Al2O, AlH, AlCl, AlCl2, AlCl3
*     ===================================================================
      g(AlII)  = gk( AlII  )
      g(AlH)   = gk( AlH   )
      g(AlOH)  = gk( AlOH  )
      g(AlO2H) = gk( AlO2H )
      g(Al2O)  = gk( Al2O  )
      g(AlCl)  = gk( AlCl  )
      g(AlCl2) = gk( AlCl2 )
      g(AlCl3) = gk( AlCl3 )
      pAlges   = eps(Al) * anHges * KT
      ppp      = 1.Q0 + g(AlII)/pel + pO*pH*g(AlOH) + pO**2*pH*g(AlO2H) 
     &         + pCl*g(AlCl) + pCl**2*g(AlCl2) + pCl**3*g(AlCl3)  
     &         + pH*g(AlH)
      ppp      = ppp / pO / g(Al2O) / 2.Q0
      qqq      = -pAlges / pO / g(Al2O) / 2.Q0
      pAl      = vieta(ppp,qqq)
      anmono(Al)   = pAl / kT
      anmol(AlH)   = g(AlH)   * pAl * pH / kT
      anmol(AlOH)  = g(AlOH)  * pAl * pO * pH / kT
      anmol(AlO2H) = g(AlO2H) * pAl * pO**2 * pH / kT
      anmol(Al2O)  = g(Al2O)  * pAl**2 * pO / kT
      anmol(AlCl)  = g(AlCl)  * pAl * pCl / kT
      anmol(AlCl2) = g(AlCl2) * pAl * pCl**2 / kT
      anmol(AlCl3) = g(AlCl3) * pAl * pCl**3 / kT
*
*     ! Gleichgewicht  K, K+, KCl, KOH, (KOH)2, CKN, C2K2N2, K2Cl2
*     ============================================================
      g(KII)    = gk( KII )
      g(KCl)    = gk( KCl )
      g(KOH)    = gk( KOH )
      g(K2O2H2) = gk( K2O2H2 )
      g(CKN)    = gk( CKN )
      g(C2K2N2) = gk( C2K2N2 )
      g(K2CL2)  = gk( K2Cl2 )
      pKges     = eps(K) * anHges * kT
      aa  = 2.Q0*( pO**2*pH**2*g(K2O2H2) + pC**2*pN**2*g(C2K2N2) 
     &            +pCl**2*g(K2Cl2) )
      ppp = ( 1.Q0 + g(KII)/pel + pCl*g(KCl) + pO*pH*g(KOH) 
     &        + pC*pN*g(CKN) ) /aa 
      qqq = -pKges/aa
      pK  = VIETA(ppp,qqq)
      anmono(K)     = pK / kT
      anmol(KCl)    = g(KCl)    * pK * pCl / kT
      anmol(KOH)    = g(KOH)    * pK * pO * pH / kT
      anmol(K2O2H2) = g(K2O2H2) * (pK*pO*pH)**2 / kT
      anmol(CKN)    = g(CKN)    * pK * pC * pN / kT
      anmol(C2K2N2) = g(C2K2N2) * (pK*pC*pN)**2 / kT
      anmol(KCl)    = g(KCl)    * pK * pCl / kT
*
*     ! Gleichgewicht  Li, LiO, LiH, LiOH, LiCl, Li2Cl2, Li2O2H2
*     ==========================================================
      g(LiII)    = gk( LiII )
      g(LiH)     = gk( LiH  )
      g(LiO)     = gk( LiO  )
      g(LiCl)    = gk( LiCl )
      g(Li2Cl2)  = gk( Li2Cl2 )
      g(LiOH)    = gk( LiOH )
      g(Li2O2H2) = gk( Li2O2H2 )
      pLiges  = eps(Li) * anHges * kT
      aa  = 2.Q0*(pO**2*pH**2*g(Li2O2H2)+pCl**2*g(Li2Cl2))
      ppp = ( 1.Q0 + g(LiII)/pel + pH*g(LiH) + pCl*g(LiCl)
     &        + pO*g(LiO)  + pO*pH*g(LiOH) )/aa
      qqq = -pLiges/aa
      pLi = VIETA(ppp,qqq)
      anmono(Li)   = pLi / kT
      anmol(LiH)   = g(LiH)   * pLi * pH / kT
      anmol(LiO)   = g(LiO)   * pLi * pO / kT
      anmol(LiOH)  = g(LiOH)  * pLi * pO * pH / kT
      anmol(LiCl)  = g(LiCl)  * pLi * pCl / kT

      else
*
*     ! estimate atomic pressures: new method
*     =======================================
      do i=1,nml
        if (i.ne.TiC) g(i)=gk(i)       ! compute all equil.constants
      enddo  
      Nseq = nel
      done(:) = .false.                ! all elements to be estimated here
      if (charge) done(el)=.true.      ! ... except for the electrons      
      if (charge) Nseq=nel-1
      eseq(:) = 0                      ! hirachical sequence of elements
      do ido=1,Nseq
        !---------------------------------------------------------
        ! search for the most abundant element not yet considered 
        !---------------------------------------------------------
        emax = 0.Q0 
        enew = 0
        do e=1,nel
          if (done(e)) cycle   
          if (eps(e)<emax) cycle
          emax = eps(e)
          enew = e                     
        enddo  
        if (verbose>1) print*,'estimate p'//trim(catm(enew))//' ...'
        done(enew) = .true.
        eseq(ido) = enew               ! add to hirachical sequence 
        pges = eps(enew)*anHges*kT
        pwork = pges
        !-------------------------------------------
        ! store coeff for Sum_l coeff(l) p^l = pges 
        !-------------------------------------------
        coeff(:) = 0.Q0          
        do i=1,nml
          affect = .false. 
          known  = .true. 
          pmol = g(i)
          do j=1,m_kind(0,i)
            e = m_kind(j,i) 
            if (.not.done(e)) then
              known = .false.
              exit
            endif  
            pat = anmono(e)*kT
            if (e==enew) then
              l = m_anz(j,i)   
              affect = .true.
            else if (m_anz(j,i).gt.0) then
              do kk=1,m_anz(j,i)
                pmol = pmol*pat
              enddo
            else
              do kk=1,-m_anz(j,i)
                pmol = pmol/pat
              enddo
            endif  
          enddo  
          if (.not.affect) cycle  
          if (.not.known) cycle
          coeff(l) = coeff(l) + l*pmol
          !------------------------------------
          ! for initial guess, consider this 
          ! molecule to have all of element e2 
          !------------------------------------
          pwork = MIN(pwork,(pges/(l*pmol))**(1.d0/REAL(l)))
          !if (verbose>1) print'(A10,1pE10.3)',cmol(i),pwork
        enddo  
        !----------------------------------------------
        ! solve 1d equation above with Newton's method 
        !----------------------------------------------
        do piter=1,99                  
          f  = pwork-pges
          fs = 1.Q0
          do l=1,12
            if (coeff(l)==0.d0) cycle
            f  = f  + coeff(l)*pwork**l
            fs = fs + coeff(l)*l*pwork**(l-1)
          enddo
          delta = f/fs
          pwork = pwork-delta
          if (verbose>1) print'(A2,I3,1pE25.15,1pE10.2)',
     >                   catm(enew),piter,pwork,delta/pwork
          if (ABS(delta)<1.Q-4*ABS(pwork)) exit 
        enddo  
        if (piter>=99) then
          write(*,*) "*** no convergence in 1D pre-it "//catm(enew)
          write(*,*) coeff
          stop
        endif  
        anmono(enew) = pwork*kT1
        !-----------------------------------------------------------
        ! take into account feedback on elements considered before,
        ! unless they are much more abundant, with Newton-Raphson 
        !-----------------------------------------------------------
        eact(:) = .false.
        Nact = 0
        do iredo=MAX(1,ido-NewBackIt),ido
          e = eseq(iredo)
          if (eps(e)<100*eps(enew)) then
            eact(e) = .true. 
            Nact = Nact+1
            all_to_act(e) = Nact
            act_to_all(Nact) = e
            pbefore(e) = anmono(e)
            anmono(e) = anmono(e)*pcorr(ido,e) 
          endif
        enddo
        if (verbose>1) print*,catm(eseq(1:ido))
        if (verbose>1) print*,eact(eseq(1:ido))
        if (verbose>1) print'("corr",99(1pE11.2))',
     >                 pcorr(ido,act_to_all(1:Nact))
        do it=1,99
          do ii=1,Nact
            i = act_to_all(ii) 
            FF(ii) = anHges*eps(i)*kT - anmono(i)*kT
            scale(i) = anmono(i)  
            DF(ii,:) = 0.Q0
            DF(ii,ii) = -scale(i)
            pmono1(i) = scale(i) / (anmono(i)*kT)
          enddo	
          do i=1,nml
            affect = .false. 
            known  = .true.
            do j=1,m_kind(0,i)
              if (.not.done(m_kind(j,i))) then
                known = .false.
                exit
              endif
              if (eact(m_kind(j,i))) affect=.true.
            enddo  
            if (known.and.affect) then
              pmol = g(i)
              do j=1,m_kind(0,i)
                pat = anmono(m_kind(j,i))*kT
                if (m_anz(j,i).gt.0) then
                  do kk=1,m_anz(j,i)
                    pmol = pmol*pat
                  enddo
                else
                  do kk=1,-m_anz(j,i)
                    pmol = pmol/pat
                  enddo
                endif
              enddo
              do j=1,m_kind(0,i)
                m1 = m_kind(j,i)
                if (.not.eact(m1)) cycle
                m1 = all_to_act(m1)
                term   = m_anz(j,i) * pmol
                FF(m1) = FF(m1) - term
                do l=1,m_kind(0,i)
                  m2 = m_kind(l,i)
                  if (.not.eact(m2)) cycle
                  jj = all_to_act(m2)
                  DF(m1,jj) = DF(m1,jj) - m_anz(l,i)*term*pmono1(m2)
                enddo	    
              enddo
              !if (it==1) print'(A12,1pE12.4)',cmol(i),pmol
            endif  
          enddo
          call GAUSS16(nel,Nact,DF,dp,FF)
          do ii=1,Nact
            i = act_to_all(ii) 
            dp(ii) = dp(ii)*scale(i)
          enddo
          fak = 1.Q0+4.Q0*EXP(-(MAX(0,it-20))/13.Q0)
          delta = 0.Q0
          do ii=1,Nact
            i = act_to_all(ii)
            delp = -dp(ii)/(anmono(i)*kT)
            delta = MAX(delta,ABS(delp))
            delp = -dp(ii)*kT1
            nold = anmono(i)
            anmono(i) = MAX(nold/fak,MIN(nold*fak,nold+delp))
          enddo
          if (verbose>1) print'(I3,99(1pE12.4))',
     >                   it,anmono(act_to_all(1:Nact))*kT,delta
          if (delta<1.Q-4) exit
        enddo  
        if (delta>1.Q-4) then
          write(*,*) "*** no convergence in NR pre-it "
          print*,"Tg=",Tg
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
          stop
        endif  
        pcorr(ido,:) = 1.Q0
        do ii=1,Nact
          i = act_to_all(ii)
          pcorr(ido,i) = anmono(i)/pbefore(i)    ! save after/before for next run
        enddo 
        if (verbose>1) print'("corr",99(1pE11.2))',
     >                 pcorr(ido,act_to_all(1:Nact))
        if (verbose>1) read(*,'(A1)') char
      enddo  
      endif
*
*     ! redo electron density
*     =======================
      if (charge) then
        coeff(:) = 0.Q0
        do i=1,nml
          if (i.ne.TiC) g(i)=gk(i)
          pmol = g(i)
          l=0
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_kind(j,i)==el) then
              l = m_anz(j,i)   
            else if (m_anz(j,i).gt.0) then
              do kk=1,m_anz(j,i)
                pmol = pmol*pat
              enddo
            else
              do kk=1,-m_anz(j,i)
                pmol = pmol/pat
              enddo
            endif  
          enddo
          if (l.ne.0) coeff(l)=coeff(l)+pmol
        enddo
        pel = SQRT(coeff(-1)/(1.Q0+coeff(+1)))     ! 0 = pel - a/pel + b*pel
        anmono(el) = pel/kT
      endif  

*     ! use memory of deviations between predicted atom pressures 
*     ! and converged atom pressures to improve the predictions
*     ============================================================
      ansave = anmono
      anmono = anmono*badness
*     
*-----------------------------------------------------------------------
 200  continue

      if ( alle ) then
        ! alle Molekuele mitrechnen
*       ===========================
        do i=1,nml
          if (i.ne.TiC) g(i)=gk(i)
          if (verbose>1.and.g(i)>1.Q+300) then
            print'("huge kp",A12,1pE12.3E4,I2)',cmol(i),g(i),fit(i)
          else if (g(i)>exp(9999.Q0)) then
            print'("*** limited kp",A12,1pE12.3E4,I2)',
     >                                          cmol(i),g(i),fit(i)
          endif 
          pmol = g(i)
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_anz(j,i).gt.0) then
              do kk=1,m_anz(j,i)
                pmol = pmol*pat
              enddo
            else
              do kk=1,-m_anz(j,i)
                pmol = pmol/pat
              enddo
            endif
          enddo
  	  anmol(i) = pmol*kT1
	enddo
        if (verbose>1) then
          imin = MINLOC(g(1:nml),1) 
          imax = MAXLOC(g(1:nml),1) 
          print'("min kp: ",A12,1pE12.3E4)',cmol(imin),g(imin)
          print'("max kp: ",A12,1pE12.3E4)',cmol(imax),g(imax)
        endif  
      endif  

*-----------------------------------------------------------------------
      if (nachit) then
*       ! Jacobi matrix and rhs vector for Newton-Raphson
*       =================================================
        it = 0
        eact(:) = .true.
        conv(:,:) = 9.Q+99
        switchoff(:) = 0
 300    continue
        Nact = 0
        ii = 0
        do i=1,nel
          if (.not.eact(i)) cycle
          Nact = Nact+1
          ii = ii+1
          all_to_act(i) = ii
          act_to_all(ii) = i
          FF(ii) = anHges*eps(i)*kT - anmono(i)*kT
          scale(i)  = anmono(i)  
          DF(ii,:)  = 0.Q0
          DF(ii,ii) = -scale(i)
          pmono1(i) = scale(i) / (anmono(i)*kT)
        enddo	
        do i=1,nml
          pmol = g(i)
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_anz(j,i).gt.0) then
              do kk=1,m_anz(j,i)
                pmol = pmol*pat
              enddo
            else
              do kk=1,-m_anz(j,i)
                pmol = pmol/pat
              enddo
            endif
          enddo
          anmol(i) = pmol*kT1
          do j=1,m_kind(0,i)
            m1 = m_kind(j,i)
            if (.not.eact(m1)) cycle
            m1 = all_to_act(m1)
            term   = m_anz(j,i) * pmol
            FF(m1) = FF(m1) - term
            do l=1,m_kind(0,i)
              m2 = m_kind(l,i)
              if (.not.eact(m2)) cycle
              jj = all_to_act(m2)
              DF(m1,jj) = DF(m1,jj) - m_anz(l,i)*term*pmono1(m2)
            enddo	    
          enddo
        enddo

*       ! compute the Newton-Naphson step	  
*       =================================
        FF0 = FF
        DF0 = DF
        call GAUSS16(nel,Nact,DF,dp,FF)
        !--- re-scale ---
        do ii=1,Nact
          i = act_to_all(ii) 
          dp(ii) = dp(ii)*scale(i)
        enddo  

*       ! limit NR-step and check convergence
*       =====================================
        limit = 1.Q0                                   ! limit step, keep direction
        converge(it) = 0.Q0
        Nconv = 0
        txt = ""
        do i=1,nel
          if (.not.eact(i)) then
            Nconv = Nconv+1 
            txt = trim(txt)//" "//catm(i)
          else 
            ii = all_to_act(i) 
            delp = -dp(ii)/(anmono(i)*kT)              ! relative change dx/x
            conv(it,i) = delp
            converge(it) = MAX(converge(it),ABS(delp))
            if (ABS(delp)<finish) then
              Nconv = Nconv+1 
              txt = trim(txt)//" "//catm(i)
            endif  
            !if (delp.gt.0.Q0) then
            !  limit = min(limit,(fak-1.Q0)/delp)       ! such that xnew=xold*fac 
            !else if (delp.lt.0.Q0) then
            !  limit = min(limit,(1.Q0/fak-1.Q0)/delp)  ! such that xnew=xold/fac
            !endif
          endif  
        enddo
        if (verbose>1.and.it==0) then
          write(*,*) 
          print'(7x,A14,A14,A14)',"natom","dnatom","badness" 
          do ii=1,Nact
            i = act_to_all(ii) 
            print'(A7,2(1pE14.6),0pF14.9)',catm(i),anmono(i),
     >                      -dp(ii)/(anmono(i)*kT),badness(i)
          enddo
        endif
        
*       ! apply limited NR step
*       =======================
        fak = 5.Q0
        !fak = 1.Q0+4.Q0*EXP(-(MAX(0,it-20))/13.Q0)
        do ii=1,nact
          i = act_to_all(ii)
          delp = -dp(ii)*kT1
          nold = anmono(i)
          anmono(i) = MAX(nold/fak,MIN(nold*fak,nold+delp))
        enddo
        if (it>itmax-10) then
          verbose=2
          do ii=1,nact
            i = act_to_all(ii) 
            print'(A3,2(1pE12.3))',catm(i),
     >           anmono(i),-dp(ii)/(anmono(i)*kT) 
          enddo  
        endif  
        crit = MAXVAL(converge(MAX(0,it-1):it))
        !if (verbose>1) print'(99(1pE36.28E3))',
     >  !                  anmono(act_to_all(1:nact))
        if (verbose>1) print'(i3,i3,2(1pE9.1)," converged(",i2,"):",
     >                    A50)',it,Nact,converge(it),limit,Nconv,txt
        if (it==itmax) then 
          write(*,*) '*** keine Konvergenz in SMCHEM!'
          write(*,*) 'it, converge, ind =',it,converge(it),limit
          write(*,*) '  n<H>, T =',anhges,Tg
          !if (from_merk) then
            chemiter  = chemiter + it
            from_merk = .false.
            verbose = 2             
            goto 100                   ! try again from scratch before giving up
          !endif  
          do ii=1,nact
            i = act_to_all(ii)
            write(*,*) catm(i),eps(i),-dp(ii)/(anmono(i)*kT) 
          enddo  
          stop                         ! give up.
        endif
        if (it>=5) then
          j = 0 
          do ii=1,Nact
            i = act_to_all(ii)
            if (MAXVAL(ABS(conv(it-5:it,i)))<finish) then
              switchoff(i) = it
              eact(i) = .false.
              j = j+1
              if (verbose>1) then
                print*,"switching off "//catm(i)//" ..."
              endif  
            endif
          enddo
          Nact = Nact-j
        endif  
        it = it + 1
        if (verbose.gt.1) read(*,'(a1)') char
        if (crit>finish.and.Nact>0) goto 300       ! continue iterating
*
*       ! redo rare elements
*       ====================
        redo(:) = .false.
        do iredo=1,nel
          atmax = 0.Q0 
          e = 0
          do i=1,nel
            atfrac = anmono(i)/anHges
            if (redo(i)) cycle   
            if (atfrac>1.Q-100) cycle   
            if (atfrac<atmax) cycle
            atmax = atfrac
            e = i
          enddo  
          if (e==0) exit
          redo(e) = .true.
          coeff(:) = 0.Q0
          do i=1,nml
            pmol = g(i)
            l=0
            do j=1,m_kind(0,i)
              pat = anmono(m_kind(j,i))*kT
              if (m_kind(j,i)==e) then
                l = m_anz(j,i)   
              else if (m_anz(j,i).gt.0) then
                do kk=1,m_anz(j,i)
                  pmol = pmol*pat
                enddo
              else
                do kk=1,-m_anz(j,i)
                  pmol = pmol/pat
                enddo
              endif  
            enddo
            if (l.ne.0) coeff(l)=coeff(l)+pmol
          enddo
          pat = anmono(e)*kT
          if (verbose>1) print'(2x,A25,A10)',
     >                   "redo rare element patom","dp/patom"
          do piter=1,99
            f  = pat-eps(e)*anHges*kT
            fs = 1.Q0
            do l=-1,12
              if (coeff(l)==0.Q0) cycle
              f  = f  + coeff(l)*l*pat**l
              fs = fs + coeff(l)*l**2*pat**(l-1)
            enddo
            delta = f/fs
            pat = pat-delta
            if (verbose>1) print'(A2,1pE25.15,1pE10.2)',
     >                            catm(e),pat,delta/pat
            if (ABS(delta/pat)<finish) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** no convergence in post-it "//catm(e)
            write(*,*) coeff
          endif  
          anmono(e) = pat/kT  
        enddo
*
*       ! how bad was the initial guess?
*       ================================
        if (.not.from_merk) then
          if (verbose>1) print'(7x,3(A14))',
     >                   "natom","conv.","init.guess"
          do i=1,nel
            badness(i) = anmono(i)/ansave(i)
            switch = switchoff(i)
            if (switch==0) switch=it-1
            if (verbose>1) then
              print'(A7,2(1pE14.6),0pF14.9)',catm(i),anmono(i),
     >                               conv(switch,i),badness(i)
            endif  
          enddo
          if (ilauf==1) write(99,'(A9,A10,A4,99(A10))') 
     >          'Tg','n<H>','it',catm(1:nel)
          write(99,'(0pF9.3,1pE10.3,I4,99(1pE10.3))') 
     >          Tg,anHges,it,badness
        endif  
*
*       ! final anmol determination
*       ===========================
        amerk = anmono/anHges
        do i=1,nml
          pmol = g(i)
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_anz(j,i).gt.0) then
              do kk=1,m_anz(j,i)
                pmol = pmol*pat
              enddo
            else
              do kk=1,-m_anz(j,i)
                pmol = pmol/pat
              enddo
            endif
          enddo
          anmol(i) = pmol*kT1
        enddo
        if (charge) pel=anmono(el)*kT

      endif      ! nachit

      if (ngestst) then
*       ! Test auf Elementerhaltung
*       ===========================
        do i=1,nel
          nges(i) = anmono(i)
        enddo
        do i=1,nml
          do j=1,m_kind(0,i)
            j1 = m_kind(j,i)
            nges(j1) = nges(j1) + m_anz(j,i)*anmol(i)
          enddo
        enddo
        do e=3,nel
          soll  = anHges * eps(e)
          haben = nges(e)
          abw   = ABS(soll-haben)/MAX(soll,haben)
          if (abw>1.Q-15) then
            if (verbose>1) then
              print'("*** element conservation error ",A2)',catm(e)
              print'(A12,1pE14.7)',catm(e),anmono(e)/(eps(e)*anHges)
            endif  
            sum = anmono(e)/(eps(e)*anHges)
            do i=1,nml
              do j=1,m_kind(0,i)
                j1 = m_kind(j,i)
                cc = m_anz(j,i)*anmol(i)/(eps(e)*anHges)
                if (j1==e.and.cc>1.Q-7) then
                  if (verbose>1) print'(A12,1pE14.7)',cmol(i),cc
                  sum = sum+cc
                endif  
              enddo
            enddo
            if (verbose>1) then
              print'(3(1pE14.7))',soll/anHges,haben/anHges,sum
            endif  
            nachit = .true.
            from_merk = .false.
            ansave = anmono
            goto 200
          endif
        enddo
      endif
      
      if (verbose.gt.0) print '("  ==> smchem used it=",I3,
     &                          " conv=",1pE9.2)',it,crit

      if (verbose.gt.1) read(*,'(a1)') char

      chemcall = chemcall + 1
      chemiter = chemiter + it


      CONTAINS       ! internal functions - not visible to other units 
************************************************************************
      FUNCTION gk(i)
************************************************************************
*****  kp [cgs] for different fit formula                          *****
*****  fit=1  :  Gail's polynom                                    *****
*****  fit=2  :  Tsuji's polynom                                   *****
************************************************************************
      use CHEMISTRY,ONLY: a,th1,th2,th3,th4,TT1,TT2,TT3,fit,natom,cmol
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),parameter :: bar=1.Q+6, atm=1.013Q+6, Rcal=1.987Q+0
      real(kind=qp),parameter :: Rgas=8.3144598Q+0
      real(kind=qp),parameter :: ln10=LOG(10.Q0)
      real(kind=qp),parameter :: lnatm=LOG(atm), lnbar=LOG(bar)
      integer,intent(in) :: i    ! index of molecule
      real(kind=qp) :: gk,dG,lnk ! return kp in [cgs]
      if (i.eq.0) then
        gk = 1.Q-300             ! tiny kp for unassigned molecules
        return
      endif
      if (fit(i).eq.1) then
        !---------------------
        ! ***  Gail's fit  *** 
        !---------------------
        lnk = a(i,0) + a(i,1)*th1 + a(i,2)*th2 
     &               + a(i,3)*th3 + a(i,4)*th4 
      else if (fit(i).eq.2) then
        !---------------------------
        ! ***  Tsuji (1973) fit  *** 
        !---------------------------
        lnk = ln10*( - a(i,0) - a(i,1)*th1 - a(i,2)*th2
     &                        - a(i,3)*th3 - a(i,4)*th4 ) 
      else if (fit(i).eq.3) then  
        !---------------------------------
        ! ***  Sharp & Huebner (1990)  ***
        !---------------------------------
        dG  = a(i,0)/TT1 + a(i,1) + a(i,2)*TT1 + a(i,3)*TT2 + a(i,4)*TT3
        lnk = -dG/(Rcal*TT1) + (1-Natom(i))*lnatm

      else if (fit(i).eq.4) then
        !-----------------------------------
        ! ***  Stock (2008) & Kietzmann  ***
        !-----------------------------------
        dG  = a(i,0)/TT1+a(i,1)*LOG(TT1)+a(i,2)+a(i,3)*TT1+a(i,4)*TT2
        lnk = dG + (1-Natom(i))*lnbar

      else if (fit(i).eq.5) then
        !--------------------
        ! ***  dG(T)-fit  ***
        !--------------------
        dG  = a(i,0)/TT1 + a(i,1) + a(i,2)*TT1 + a(i,3)*TT2 + a(i,4)*TT3
        lnk = -dG/(Rgas*TT1) + (1-Natom(i))*lnbar
         
      else if (fit(i).eq.6) then
        !-------------------------------
        ! ***  Barklem & Collet fit  ***
        !-------------------------------
        lnk = a(i,0)/TT3 + a(i,1)/TT2 + a(i,2)/TT1 + a(i,3)/TT1**0.05d0
     &      + a(i,4)*LOG(TT1) + a(i,5) + a(i,6)*TT1 + a(i,7)*TT2
         
      else
        print*,cmol(i),"i,fit=",i,fit(i)
        stop "???"
      endif  
      gk = EXP(MIN(1.Q+4,lnk))
      end FUNCTION gk

      end SUBROUTINE smchem16
