************************************************************************
      SUBROUTINE smchem8 (anHges,Tg,eps,anmono,anmol,merk,verbose)
************************************************************************
*                                                                      *
*     small chemistry                                                  *
*     ---------------                                                  *
*     Diese Routine berechnet eine GG-Chemie                           *
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
************************************************************************
      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,
     >                    NewFastLevel,NewPreMethod,
     >                    nml=>NMOLE,nel=>NELM,cmol,catm,
     >                    m_kind,m_anz,charge,elion,el,
     >                    Natom,Ncmax,STOImax,H,C,N,O,
     >                    th1,th2,th3,th4,TT1,TT2,TT3
      use EXCHANGE,ONLY: chemcall,chemiter,preUse,preIter,preEst,
     >                   DUALcorr,HCOcorr
      implicit none
*-----------------------------------------------------------------------
*  Dimensionierung fuer die Molekuel- und Atom Felder. Hier ist auf
*  Konsistenz mit dem aufrufenden Programm zu achten.
*-----------------------------------------------------------------------
      real*8,intent(in)  :: anHges,Tg
      real*8,intent(in)  :: eps(nel)
      real*8,intent(out) :: anmono(nel),anmol(nml)
      integer,intent(inout) :: verbose
      real*8,parameter :: bk=1.380662d-16
*-----------------------------------------------------------------------
*  Bei merk=.true. merkt sich die Routine die letzte konvergiert Loesung
*  und geht beim naechsten Mal von diesen Startwerten aus.
      logical,intent(INOUT) :: merk
*-----------------------------------------------------------------------
*  Die Variable "ngestst" entscheidet, ob die Elementerhaltung ueber-
*  prueft werden soll.
      logical,parameter :: ngestst=.false.
*-----------------------------------------------------------------------
      integer :: stindex,Nconv,switch,ido,iredo,nu,tmp
      integer :: Nact,all_to_act(nel),act_to_all(nel),switchoff(nel)
      integer :: e,i,j,j1,ii,jj,l,it,m1,m2,piter,ifatal,ipull,pullmax
      integer :: lmin,lmax,Nseq,iloop,imethod,enew,i1,i2,i3,s1,s2,s3,s4
      integer :: eseq(nel),imaj(nel),imaj2(nel)
      integer :: NpreLoop,NpreIt,Ntaken,Nestim,DUALco,HCOco
      integer :: ind,indx(nel),info,ipvt(nel)
      integer,parameter :: itmax=200
      real*8 :: finish,qual,qual0,qual1,qual2,qual3,qmost,qq
      real*8 :: g(0:nml),limit
      real*8 :: pH2O,pCO2,pCH4,lpH2O,lpCO2,lpCH4,lpH,lpO,lpC
      real*8 :: kT,kT1,cc,nelek,ng,Sa,fak,lth,arg,term,f,fs
      real*8 :: pel,delta,pat,atfrac,atmax,lnp,lnp3
      real*8 :: nges(nel),coeff(-3:Ncmax),lnc(-3:Ncmax)
      real*8 :: DF(nel,nel),dp(nel),FF(nel),pmol,crit,delp,nold
      real*8 :: DF0(nel,nel),FF0(nel),xscale(nel),fscale(nel)
      real*8 :: nsave(nel),null(nel),sortq(nel)
      real*8 :: conv(0:itmax,nel),converge(0:itmax)
      real*8 :: soll,haben,abw,sum,maxs,mcharge
      real*8 :: pbefore(nel),norm(nel),xx(nel)
      real*8 :: emax,pges,pwork,pwork2,pp,psc,ptest,aim
      real*8 :: work(nel*(nel+1))
      real*8 :: condnum1,work2(nel)
      logical :: from_merk,eact(nel),redo(nel),done(nel),affect,known
      logical :: relevant(nml)
      logical :: possible,ptake,HCOproblem,hasH,hasC,hasO,IS_NAN
      character(len=5000) :: mols
      character(len=100) :: txt
      character(len=1) :: char,bem
      integer,save :: TiC,H2O,CH4,CO2,ilauf=0
      real*8,allocatable,save :: amerk(:),ansave(:)
      real*8,allocatable,save :: badness(:),pcorr(:,:) 
      integer,allocatable,save :: pkey(:)
!$omp threadprivate(TiC,ilauf,amerk,ansave,badness,pcorr,pkey)
*-----------------------------------------------------------------------      

      ifatal = 0
      if (.not.allocated(amerk)) then
        allocate(badness(nel),pcorr(nel,nel),pkey(nel),
     >           amerk(nel),ansave(nel))
        TiC = stindex(cmol,nml,'TIC    ')
        H2O = stindex(cmol,nml,'H2O    ')
        CH4 = stindex(cmol,nml,'CH4    ')
        CO2 = stindex(cmol,nml,'CO2    ')
        badness = 1.d0
        pcorr   = 1.d0
        pkey    = 0
      endif

*-----------------------------------------------------------------------
*     ! zu niedrige Temperaturen abfangen und
*     ! Variable fuer die Dissoziationskonstanten berechnen
*     =====================================================
      TT1 = Tg
      TT2 = TT1*TT1
      TT3 = TT2*TT1
      th1 = 5040.d0/TT1
      th2 = th1*th1
      th3 = th2*th1
      th4 = th3*th1
      kT  = bk*TT1
      kT1 = 1.d0/kT
*      
*-----------------------------------------------------------------------
*     ! init vectors
*     ==============
      anmono = 0.d0
      anmol  = 0.d0

* --------------------------------------------------------------------------
*     ! compute equilibrium constants
*     ===============================
      do i=1,nml
        if (i.ne.TiC) g(i)=gk(i)       ! compute all equil.constants lnk
      enddo  

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
     &               + 1.56041*lth**2 - 0.93275*lth**3
      g(TiC) = -2.30256*arg

*---------------------------------------------------------------------------
      !print'("smchem8 called ilauf,merk",i5,l2)',ilauf,merk
      !print'("T,n<H>,eps=",0pF10.2,99(1pE11.3))',Tg,anHges,eps
      NpreLoop = 0
      NpreIt = 0
      Ntaken = 0
      Nestim = 0
      DUALco = 0
      HCOco  = 0
      if ((ilauf.gt.10).and.merk) then
        do i=1,nel
          anmono(i) = amerk(i) * anhges
        enddo
        from_merk = .true.
        goto 200
      endif

*---------------------------------------------------------------------------
*     ! estimate electron density
*     =========================== 
 100  continue
      if (verbose>1) print'("SMCHEM_8:")'
      from_merk = .false.
      if (charge) then
        nelek = 0.d0 
        do i=1,nel
          if (i==el) cycle 
          ng = anHges * eps(i)
          Sa = EXP(g(elion(i)))*kT1
          nelek = nelek + ng/(0.5d0 + SQRT(0.25d0 + ng/Sa))
        enddo
        anmono(el) = nelek
        pel = nelek*kT
        mcharge = pel
        !peest = pel
        !pel = pecorr*pel
        !anmono(el) = pecorr*anmono(el) 
        if (verbose>1) print'(" estimate pel=",1pE10.3)',pel
      endif  

*-----------------------------------------------------------------------
*     ! estimate atomic pressures: new method
*     =======================================
      sortq(:) = 1.Q0
 110  Nseq = nel
      done(:) = .false.                ! all elements to be estimated here
      !if (charge) done(el)=.true.     ! ... except for the electrons      
      !if (charge) Nseq=nel-1
      eseq(:) = 0                      ! hirachical sequence of elements
      ptake = .true.
      do ido=1,Nseq
        !---------------------------------------------------------
        ! search for the most abundant element not yet considered 
        !---------------------------------------------------------
        emax = 0.d0
        enew = 0
        do e=1,nel
          if (done(e)) cycle
          norm(e) = eps(e)*sortq(e)
          if (e==el) norm(e)=anmono(el)*sortq(el)/anHges
          if (norm(e)<emax.or.(ido==1.and.e==el)) cycle
          emax = norm(e)
          enew = e
        enddo
        if (enew==0) then
          print*,catm 
          print*,eps
          print*,ido
          print*,done
          stop "*** should not occur."
        endif  
        if (verbose>1) print*,'estimate p'//trim(catm(enew))//' ...'
        if (enew.ne.pkey(ido)) ptake=.false.
        done(enew) = .true.
        eseq(ido) = enew               ! add to hirachical sequence 
        pges = eps(enew)*anHges*kT
        if (verbose>0) mols = ''
        !--------------------------------------------------------------------
        ! make rough estimate of atom pressure, considering that one molecule 
        ! takes all the element.  Store coeff for Sum_l coeff(l) p^l = pges 
        !--------------------------------------------------------------------
        coeff(:) = 0.d0   
        lnc(:) = 0.d0
        pwork = pges
        pwork2 = 0.d0
        imaj(enew) = 0
        imaj2(enew) = 0
        lmin = +99
        lmax = -99
        do e=1,nel
          if (.not.done(e)) cycle
          xx(e) = LOG(anmono(e)*kT)
        enddo  
        do i=1,nml
          affect = .false.
          known  = .true.
          lnp = g(i)
          do j=1,m_kind(0,i)
            e  = m_kind(j,i)
            nu = m_anz(j,i)
            if (.not.done(e)) then
              known = .false.
              exit
            endif  
            if (e==enew) then
              l = nu
              affect = .true.
            else
              lnp = lnp + nu*xx(e)
            endif  
          enddo  
          if (.not.affect) cycle  
          if (.not.known) cycle
          if (verbose>0) mols = trim(mols)//" "//cmol(i)
          lmin = MIN(lmin,l)
          lmax = MAX(lmax,l)
          coeff(l) = coeff(l) + l*EXP(lnp)
          lnc(l) = MAX(lnc(l),LOG(REAL(l))+lnp)
          if (l>0) then
            sum = LOG(pges) - lnp - LOG(REAL(l))
            ptest = EXP(sum/l)
            if (ptest<pwork) then
              pwork2 = pwork
              pwork = ptest
              imaj2(enew) = imaj(enew)
              imaj(enew) = i
              !print*,cmol(i),pwork,lnc(l)
            else if (ptest<pwork2) then  
              pwork2 = ptest
              imaj2(enew) = i
            endif
          endif
        enddo
        if (verbose>1) print*,trim(mols)
        if (enew==el) then
          !print*,coeff(-1),coeff(1)
          pel = SQRT(-coeff(-1)/(1.d0+coeff(+1)))    ! 0 = pel - a/pel + b*pel
          if (verbose>1) print*,'pel=',anmono(el)*kT,pel
          anmono(el) = pel*kT1
          mcharge = MAX(pel,MAX(coeff(-1)/pel,coeff(1)*pel))
        else   
          !---------------------------------------
          ! scale to make coeff(:) fit into real*8  
          !---------------------------------------
          aim = LOG(1.d+10)
          psc = -LOG(pges)                           ! when atom dominates
          do l=lmin,lmax
            if (lnc(l)==0.d0) cycle
            ptest = (aim-lnc(l))/l                   ! lnc(l)+l*ln(psc) = aim
            psc = MAX(psc,-ptest)
          enddo  
          psc = -psc
          if (verbose>1) print*,"pwork=",pwork,"  psc=",psc
          !------------------------------------------------
          ! store coeff for Sum_l coeff(l) (p*psc)^l = pges 
          !------------------------------------------------
          coeff(:) = 0.d0
          do l=lmin,lmax
            if (lnc(l)==0.d0) cycle
            coeff(l) = EXP(lnc(l)+l*psc)
          enddo
          psc = EXP(psc)
          !if (verbose>1) print'(" coeff(",I2,":",I2,")=",99(1pE10.2))',
     >    !               lmin,lmax,coeff(lmin:lmax)
          !----------------------------------------------
          ! solve 1d equation above with Newton's method 
          !----------------------------------------------
          pp = pwork/psc
          do piter=1,99                  
            f  = pp*psc-pges
            fs = psc
            do l=lmin,lmax
              if (coeff(l)==0.d0) cycle
              f  = f  + coeff(l)*pp**l
              fs = fs + coeff(l)*l*pp**(l-1)
            enddo
            if (fs==0.d0) stop "*** fs=0 in smchem8 1d-pre-it."
            delta = f/fs
            pp = pp-delta
            if (verbose>1) print'(A2,I3,1pE25.15,1pE10.2)',
     >                     catm(enew),piter,pp*psc,delta/pp
            if (ABS(delta)<1.d-4*ABS(pp)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** smchem8 no conv in 1D pre-it "//catm(enew)
            write(*,*) anHges,Tg
            do i=1,nel
              write(*,'(A2,1x,0pF30.26)') catm(i),eps(i)
            enddo
            write(*,*) "psc=",real(psc),lmin,lmax
            write(*,*) "  lnc:",lnc(lmin:lmax)
            write(*,*) "coeff:",coeff(lmin:lmax)
            goto 1000
          endif  
          anmono(enew) = pp*psc*kT1
        endif  

        !-----------------------------------------------------------
        ! take into account feedback on elements considered before,
        ! unless they are much more abundant, with Newton-Raphson 
        !-----------------------------------------------------------
        nsave = anmono
 150    continue
        eact(:) = .false.
        Nact = 0
        !if (verbose>0) print*,"NewFastLevel,ptake=",NewFastLevel,ptake
        do iredo=MAX(1,ido-NewBackIt),ido
          e = eseq(iredo)
          if (norm(e)<NewBackFac*norm(enew)) then
            eact(e) = .true. 
            Nact = Nact+1
            all_to_act(e) = Nact
            act_to_all(Nact) = e
            pbefore(e) = anmono(e)
            if (NewFastLevel<3.and.ptake) then
              anmono(e) = anmono(e)*pcorr(enew,e)
              Ntaken = Ntaken+1
            else
              pcorr(enew,e) = 1.d0
            endif
            Nestim = Nestim+1
          endif
        enddo
        if (verbose>1) then
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
        endif
        if (Nact==1) then
          pkey(ido) = enew
          e = eseq(ido)
          xx(e) = LOG(anmono(e)*kT)
          cycle
        endif  
        do i=1,nml
          affect = .false. 
          known  = .true.
          do j=1,m_kind(0,i)
            e = m_kind(j,i) 
            if (.not.done(e)) then
              known = .false.
              exit
            endif
            if (eact(e)) affect=.true.
          enddo  
          relevant(i) = (known.and.affect)
        enddo

        !---------------------------------------------------------
        !***  check if new important molecule has huge impact  ***
        !---------------------------------------------------------
        i = imaj(enew)
        pmol = g(i)
        do j=1,m_kind(0,i)
          pmol = pmol + m_anz(j,i)*LOG(anmono(m_kind(j,i))*kT)
        enddo
        pmol = EXP(pmol)
        i1 = 0
        s1 = 0
        qmost = 0.Q0
        do j=1,m_kind(0,i)
          e = m_kind(j,i) 
          if (e==enew) then
            i1 = enew
            s1 = m_anz(j,i) 
            cycle
          endif
          qq = pmol*m_anz(j,i)/(eps(e)*anHges*kT)
          !print*,catm(e),cmol(i),REAL(qq)
          if (qq>qmost) then
            qmost = qq
            i2 = e
            s2 = m_anz(j,i)
          endif  
        enddo
        possible = (qmost>0.95)
        if (possible) then
          imaj2(i2) = imaj(i2)
          imaj(i2) = i
          pmol = eps(i2)*anHges*kT/s2
          if (eps(i1)*anHges*kT-s1*pmol<0.d0) then
            tmp= i1
            i1 = i2
            i2 = tmp
            tmp= s1
            s1 = s2
            s2 = tmp
          endif
          pmol = eps(i2)*anHges*kT/s2
          if (eps(i1)*anHges*kT-s1*pmol<0.d0) possible=.false.
          if (i2==el) possible=.false.
        endif  
        if (possible) then
          !print*,catm(i1),cmol(i),cmol(imaj2(i1)),
     >    !       REAL(eps(i1)*anHges*kT-s1*pmol)
          !print*,catm(i2),cmol(i),cmol(imaj2(i2)),
     >    !       REAL(eps(i2)*anHges*kT-s2*pmol)
          i3 = imaj2(i1)   ! second most important carrier of i1
          s3 = 0
          s4 = 0
          lnp3 = g(i3)
          do j=1,m_kind(0,i3)
            e = m_kind(j,i3) 
            if (e==i1) then
              s3 = m_anz(j,i3)
              cycle
            endif  
            if (e==i2) then
              s4 = m_anz(j,i3)
              cycle
            endif  
            lnp3 = lnp3 + m_anz(j,i3)*LOG(anmono(e)*kT)
          enddo
          lnp = g(i)
          do j=1,m_kind(0,i)
            e = m_kind(j,i)
            if (e==i1.or.e==i2) cycle
            lnp = lnp + m_anz(j,i)*LOG(anmono(e)*kT)
          enddo
          if (verbose>1) print*,"readjusting "//catm(i1)//" "//catm(i2)
     >                   //" with molecules "//trim(cmol(i))//" "
     >                   //trim(cmol(i3))//" ..."
          !--- solve eps(i1)*anHges*kT = s1*pmol + s3*p_second ---
          !--- with  eps(i2)*anHges*kT = s2*pmol               ---
          FF(1) = LOG((eps(i1)*anHges*kT-s1*pmol)/s3)-lnp3
          FF(2) = LOG((eps(i2)*anHges*kT)/s2)-lnp
          DF(1,1) = s3
          DF(1,2) = s4
          DF(2,1) = s1
          DF(2,2) = s2
          call GAUSS8(nel,2,DF,dp,FF)
          possible = (dp(1)<1000.0).and.(dp(2)<1000.0)
        endif
        if (possible) then
          anmono(i1) = EXP(dp(1))*kT1
          anmono(i2) = EXP(dp(2))*kT1
          !print'(A2,9(1pE16.8))',catm(i1),eps(i1)*anHges*kT,
     >    !                       s1*EXP(lnp+s1*dp(1)+s2*dp(2))
     >    !                      +s3*EXP(lnp3+s3*dp(1)+s4*dp(2))
          !print'(A2,9(1pE16.8))',catm(i2),eps(i2)*anHges*kT,
     >    !                       s2*EXP(lnp+s1*dp(1)+s2*dp(2))
          pbefore(i1) = anmono(i1)
          pbefore(i2) = anmono(i2)
          if (NewFastLevel<3.and.ptake) then
            anmono(i1) = anmono(i1)*pcorr(enew,i1)
            anmono(i2) = anmono(i2)*pcorr(enew,i2)
          endif  
          DUALco = DUALco+1
        endif  
        !--- special case CH4+CO2+H2O ---
        HCOproblem = (Nact>2).and.(enew==H.or.enew==C.or.enew==O)
     >               .and.(CH4>0).and.(H2O>0).and.(CO2>0)
        if (HCOproblem) then
          hasH = .false.
          hasC = .false.
          hasO = .false.
          do ii=1,Nact
            i = act_to_all(ii)
            if (i==H) hasH=.true.
            if (i==C) hasC=.true.
            if (i==O) hasO=.true.
          enddo  
          HCOproblem = hasH.and.hasO.and.hasC
        endif  
        if (HCOproblem) then
          pH2O = (2*eps(O)-4*eps(C)+eps(H))*anHges*kT/4
          pCO2 = (2*eps(O)+4*eps(C)-eps(H))*anHges*kT/8
          pCH4 = (eps(H)+4*eps(C)-2*eps(O))*anHges*kT/8
          HCOproblem = (pH2O>0).and.(pCO2>0).and.(pCH4>0)
        endif  
        if (HCOproblem) then
          lpH2O = LOG(pH2O)
          lpCO2 = LOG(pCO2)
          lpCH4 = LOG(pCH4)
          lpH = (lpCH4+lpH2O*2+g(CO2)-lpCO2-g(CH4)-g(H2O)*2)/8
          lpO = (lpCO2+lpH2O*2+g(CH4)-lpCH4-g(CO2)-g(H2O)*2)/4
          lpC = (lpCH4+lpCO2+g(H2O)*2-lpH2O*2-g(CH4)-g(CO2))/2
          if (verbose>1) print*,"applying H2O-CH4-CO2 correction ..."
          !print'(2(1pE16.8))',2*pH2O+4*pCH4,anHges*kT*eps(H)
          !print'(2(1pE16.8))',  pH2O+2*pCO2,anHges*kT*eps(O)
          !print'(2(1pE16.8))',  pCH4+  pCO2,anHges*kT*eps(C)
          !print'(2(1pE16.8))',pH2O,EXP(lpH*2+lpO+g(H2O))
          !print'(2(1pE16.8))',pCO2,EXP(lpC+lpO*2+g(CO2))
          !print'(2(1pE16.8))',pCH4,EXP(lpC+lpH*4+g(CH4))
          anmono(H) = EXP(lpH)*kT1
          anmono(O) = EXP(lpO)*kT1
          anmono(C) = EXP(lpC)*kT1
          pbefore(H) = anmono(H)
          pbefore(O) = anmono(O)
          pbefore(C) = anmono(C)
          if (NewFastLevel<3.and.ptake) then
            anmono(H) = anmono(H)*pcorr(enew,H)
            anmono(O) = anmono(O)*pcorr(enew,O)
            anmono(C) = anmono(C)*pcorr(enew,C)
          endif  
          HCOco = HCOco + 1
        endif  
          
        !--- solve with different methods ---
        if (verbose>1) then 
          print'("corr",99(1pE11.2))',pcorr(enew,act_to_all(1:Nact))
        endif  
        qual  = 9.d+99
        do iloop=1,4
          if (iloop==1) then
            imethod = NewPreMethod                    ! first choice
          elseif (iloop==2) then
            imethod = NewPreMethod                    ! second try with pcorr=1.0
          elseif (iloop==3) then
            imethod = MOD(NewPreMethod,3)+1
          elseif (iloop==4) then
            imethod = MOD(NewPreMethod+1,3)+1
          endif
          if (verbose>1) print*,"imethod=",iloop,imethod
          !write(97,*) catm(eseq(1:ido))
          !write(97,*) eact(eseq(1:ido))
          !write(97,*) imethod
          !write(97,*) real(pcorr(enew,act_to_all(1:Nact)))          
          if (imethod==1) then
            !-------- method 1: xx=log(patm)-variables --------
            null = anmono
            qual0 = qual 
            qual1 = qual 
            qual2 = qual 
            do it=1,299
              do ii=1,Nact
                i = act_to_all(ii)
                xx(i) = LOG(anmono(i)*kT)
                xscale(i) = anmono(i)*kT
                if (i==el) then
                  fscale(ii) = 1.d0/mcharge
                else  
                  fscale(ii) = 1.d0/(anHges*eps(i)*kT)
                endif  
                FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                DF(ii,:)  = 0.d0
                DF(ii,ii) = -fscale(ii)*anmono(i)*kT
              enddo  
              do i=1,nml
                if (.not.relevant(i)) cycle 
                pmol = g(i)
                do j=1,m_kind(0,i)
                  pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                enddo
                pmol = EXP(pmol)
                do j=1,m_kind(0,i)
                  m1 = m_kind(j,i)
                  if (.not.eact(m1)) cycle
                  ii = all_to_act(m1)
                  term = fscale(ii)*m_anz(j,i)*pmol
                  FF(ii) = FF(ii) - term
                  do l=1,m_kind(0,i)
                    m2 = m_kind(l,i)
                    if (.not.eact(m2)) cycle
                    jj = all_to_act(m2)
                    DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                  enddo	    
                enddo
              enddo
              call GAUSS8(nel,Nact,DF,dp,FF)
              bem = " "
              qual3 = qual2
              qual2 = qual1
              qual1 = qual0
              qual0 = qual
              qual  = 0.d0
              do ii=1,Nact
                qual = qual + ABS(dp(ii))
              enddo  
              maxs = 3.0
              if (it>30.and.(qual >qual0.or.qual0>qual1.or.
     >                       qual1>qual2.or.qual2>qual3)) then
                maxs = 3.0*exp(-MAX(0,it-30)/70.0)
                bem = "*"
              endif  
              do ii=1,Nact
                i = act_to_all(ii) 
                xx(i) = xx(i) - MAX(-maxs,MIN(maxs,dp(ii)))
                anmono(i) = exp(xx(i))*kT1
              enddo
              !write(97,'(I4,A2,99(1pE11.3E3))')
     >        !               it,bem,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,A2,99(1pE11.3E3))',
     >                       it,bem,anmono(act_to_all(1:Nact))*kT,qual
              NpreIt = NpreIt+1
              if (it>1.and.qual<1.d-10) exit
            enddo
            NpreLoop = NpreLoop+1

          else if (imethod==2) then
            !-------- method 2: lin-variables with pullback --------
            dp(:) = 0.d0
            qual  = 9.d+99
            fak   = 1.d0
            do it=1,199
              qual0 = qual 
              null  = anmono
              pullmax = 1
              if (it>100) pullmax=10
              do ipull=1,pullmax  ! pullback if quality gets worse
                !--- make a step ---
                do ii=1,Nact
                  i = act_to_all(ii)
                  anmono(i) = null(i)-fak*dp(ii)*kT1
                  xx(i) = LOG(anmono(i)*kT)
                enddo  
                !--- determine new FF and DF ---
                do ii=1,Nact
                  i = act_to_all(ii) 
                  xscale(i) = anmono(i)*kT
                  if (i==el) then
                    fscale(ii) = 1.d0/mcharge
                  else  
                    fscale(ii) = 1.d0/(anHges*eps(i)*kT)
                  endif  
                  FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                  DF(ii,:) = 0.d0
                  DF(ii,ii) = -fscale(ii)*anmono(i)*kT
                enddo
                do i=1,nml
                  if (.not.relevant(i)) cycle 
                  pmol = g(i)
                  do j=1,m_kind(0,i)
                    pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                  enddo
                  pmol = EXP(pmol)
                  do j=1,m_kind(0,i)
                    m1 = m_kind(j,i)
                    if (.not.eact(m1)) cycle
                    ii = all_to_act(m1)
                    term = fscale(ii)*m_anz(j,i)*pmol
                    FF(ii) = FF(ii) - term
                    do l=1,m_kind(0,i)
                      m2 = m_kind(l,i)
                      if (.not.eact(m2)) cycle
                      jj = all_to_act(m2)
                      DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                    enddo	    
                  enddo
                enddo
                !--- determine new quality ---
                qual = 0.d0
                do ii=1,Nact
                  i = act_to_all(ii)
                  qual = qual + FF(ii)**2
                enddo  
                if (qual<qual0) exit
                if (ipull==pullmax) exit
                !write(97,'("pullback",3(1pE11.3))') fak,qual0,qual
                if (verbose>1) print'("pullback",3(1pE11.3))',
     >                         fak,qual0,qual
                fak = 0.5*fak   ! reduce NR-step
              enddo  
              !write(97,'(I4,99(1pE11.3))')
     >        !             it,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,99(1pE11.3))',
     >                     it,anmono(act_to_all(1:Nact))*kT,qual
              NpreIt = NpreIt+1
              if (it>1.and.qual<1.d-10) exit
              !--- determine new NR-vector ---
              call GAUSS8(nel,Nact,DF,dp,FF)
              do ii=1,Nact
                i = act_to_all(ii) 
                dp(ii) = dp(ii)*xscale(i)
              enddo
              null = anmono
              !--- limit step physically, keep direction ---
              fak = 1.d0
              do ii=1,Nact
                i = act_to_all(ii)
                if (null(i)*kT-fak*dp(ii)>5.d0*null(i)*kT) then
                  fak=MIN(fak,-4.d0*null(i)*kT/dp(ii))
                endif
                if (null(i)*kT-fak*dp(ii)<0.2d0*null(i)*kT) then
                  fak=MIN(fak,0.8d0*null(i)*kT/dp(ii))
                endif
              enddo
            enddo  
            NpreLoop = NpreLoop+1

          else if (imethod==3) then
            !-------- method 3: xx=log(patm)-variables --------
            null = anmono
            do it=1,299
              do ii=1,Nact
                i = act_to_all(ii)
                xx(i) = LOG(anmono(i)*kT)
                xscale(i) = anmono(i)*kT
                if (i==el) then
                  fscale(ii) = 1.d0/mcharge
                else  
                  fscale(ii) = 1.d0/(anHges*eps(i)*kT)
                endif  
                FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                DF(ii,:)  = 0.d0
                DF(ii,ii) = -fscale(ii)*anmono(i)*kT
              enddo  
              do i=1,nml
                if (.not.relevant(i)) cycle 
                pmol = g(i)
                do j=1,m_kind(0,i)
                  pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                enddo
                pmol = EXP(pmol)
                do j=1,m_kind(0,i)
                  m1 = m_kind(j,i)
                  if (.not.eact(m1)) cycle
                  ii = all_to_act(m1)
                  term = fscale(ii)*m_anz(j,i)*pmol
                  FF(ii) = FF(ii) - term
                  do l=1,m_kind(0,i)
                    m2 = m_kind(l,i)
                    if (.not.eact(m2)) cycle
                    jj = all_to_act(m2)
                    DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                  enddo	    
                enddo
              enddo
              call GAUSS8(nel,Nact,DF,dp,FF)
              qual = 0.d0
              do ii=1,Nact
                qual = MAX(qual,ABS(dp(ii)))
              enddo  
              fak = MIN(1.0,3.0/qual)
              do ii=1,Nact
                i = act_to_all(ii) 
                xx(i) = xx(i) - fak*dp(ii)
                anmono(i) = exp(xx(i))*kT1
                if (IS_NAN(anmono(i))) then
                  qual=2.d+99
                  exit
                endif  
              enddo
              !write(97,'(I4,A2,99(1pE11.3E3))')
     >        !               it,bem,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,99(1pE11.3E3))',
     >                       it,anmono(act_to_all(1:Nact))*kT,fak,qual
              NpreIt = NpreIt+1
              if (qual<1.d-10.or.qual>1.d+99) exit
            enddo  
            NpreLoop = NpreLoop+1
          endif
          if (qual<1.d-5) exit
          do ii=1,Nact
            i = act_to_all(ii)
            anmono(i) = pbefore(i)    ! do not use pre-it corrections again
          enddo 
          if (verbose>1) read(*,'(A1)') char
        enddo  
        if (qual>1.d-4) then
          if (ptake) then
            anmono = nsave 
            ptake = .false.
            goto 150
          endif  
          write(*,*) "*** no convergence in NR pre-it "
          print*,"Tg=",Tg
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
          print*,pcorr(enew,act_to_all(1:Nact))
          goto 1000
        endif  
        !--- save ratio after/before for next run ---
        !print*,"ido,pkey(ido),enew=",ido,pkey(ido),enew
        pkey(ido) = enew
        pcorr(enew,:) = 1.d0
        do ii=1,Nact
          i = act_to_all(ii)
          pcorr(enew,i) = anmono(i)/pbefore(i)    
        enddo 
        if (verbose>1) print'("corr",99(1pE11.2))',
     >                 pcorr(enew,act_to_all(1:Nact))
        if (verbose>1) read(*,'(A1)') char
      enddo  
*
*     ! redo electron density
*     =======================
      if (charge) then
        coeff(:) = 0.d0
        relevant = .false.
        do i=1,nml
          do j=1,m_kind(0,i)
            if (m_kind(j,i)==el) relevant(i)=.true. 
          enddo  
          if (.not.relevant(i)) cycle
          pmol = g(i)
          l = 0
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_kind(j,i)==el) then
              l = m_anz(j,i)   
            else 
              pmol = pmol + m_anz(j,i)*LOG(pat)
            endif  
          enddo
          coeff(l) = coeff(l)+EXP(pmol)
        enddo
        pel = SQRT(coeff(-1)/(1.d0+coeff(+1)))     ! 0 = pel - a/pel + b*pel
        mcharge = MAX(pel,MAX(coeff(-1)/pel,coeff(1)*pel))
        if (verbose>1) print'(" redo pel =",1pE17.10," ->",1pE17.10)',
     >                        anmono(el)*kT,pel
        anmono(el) = pel*kT1
        !pecorr = pel/peest
      endif  

*     ! use memory of deviations between predicted atom pressures 
*     ! and converged atom pressures to improve the predictions
*     ============================================================
      ansave = anmono
      if (NewFastLevel<2.and.ptake) anmono = anmono*badness
*     
 200  continue
*-----------------------------------------------------------------------
      if (NewFullIt) then
*       ! Jacobi matrix and rhs vector for Newton-Raphson
*       =================================================
        it = 0
        eact(:) = .true.
        conv(:,:) = 9.d+99
        switchoff(:) = 0
        finish = 1.d-12
 300    continue
        if (it>30) finish=10.d0**(-12.0+3.0*(it-30.0)/(itmax-30.0))
        Nact = 0
        ii = 0
        do i=1,nel
          if (.not.eact(i)) cycle
          Nact = Nact+1
          ii = ii+1
          all_to_act(i) = ii
          act_to_all(ii) = i
          xx(i) = LOG(anmono(i)*kT)
          xscale(i) = anmono(i)*kT
          if (i==el) then
            fscale(ii) = 1.d0/mcharge
          else  
            fscale(ii) = 1.d0/(anHges*eps(i)*kT)
          endif  
          FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
          DF(ii,:)  = 0.d0
          DF(ii,ii) = -fscale(ii)*anmono(i)*kT
        enddo	
        do i=1,nml
          pmol = g(i)
          do j=1,m_kind(0,i)
            pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
          enddo
          pmol = EXP(pmol)
          anmol(i) = pmol*kT1
          do j=1,m_kind(0,i)
            m1 = m_kind(j,i)
            if (.not.eact(m1)) cycle
            ii = all_to_act(m1)
            term = m_anz(j,i)*pmol*fscale(ii)
            FF(ii) = FF(ii) - term
            do l=1,m_kind(0,i)
              m2 = m_kind(l,i)
              if (.not.eact(m2)) cycle
              jj = all_to_act(m2)
              DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
            enddo	    
          enddo
        enddo

*       ! compute the Newton-Naphson step	  
*       =================================
        FF0 = FF
        DF0 = DF
        !call QGEFA ( DF, nel, nel, ipvt, info )
        !call QGECO ( DF, nel, nel, ipvt, condnum1, work2 )
        !call QGESL ( DF, nel, nel, ipvt, FF, 0 )
        !dp  = FF
        !print'("condnum1 = ",1pE12.2E3)',condnum1
        call SGEIR(DF,nel,Nact,FF,1,ind,work,indx)
        dp = FF
        if (ind<0) then
          FF = FF0
          DF = DF0
          call GAUSS8(nel,Nact,DF,dp,FF)
        endif  
        !do i=1,Nact
        !  print'(99(1pE11.3E3))',DF0(i,1:Nact),dp(i),FF0(i)
        !enddo  
        !print*,ind

*       ! apply limited NR step and check convergence 	  
*       =============================================
        if (.true.) then        
          qual = 0.d0
          do ii=1,Nact
            qual = MAX(qual,ABS(dp(ii)))
          enddo  
          fak = MIN(1.0,3.0/qual)
          converge(it) = 0.d0
          Nconv = 0
          if (verbose>0) txt = ""
          do i=1,nel
            if (.not.eact(i)) then
              Nconv = Nconv+1 
              if (verbose>0) txt = trim(txt)//" "//catm(i)
            else 
              ii = all_to_act(i) 
              delp = -dp(ii)                             ! relative change dx/x
              conv(it,i) = delp
              converge(it) = MAX(converge(it),ABS(delp))
              if (ABS(delp)<finish) then
                Nconv = Nconv+1 
                if (verbose>0) txt = trim(txt)//" "//catm(i)
              endif  
            endif  
          enddo
          if (verbose>1.and.it<100) then
            print*,"ind,fak=",ind,fak
            print'(7x,A14,A14,A14)',"patom","dpatom","badness" 
            do ii=1,Nact
              i = act_to_all(ii) 
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >             -dp(ii),badness(i)
            enddo
          endif
          do ii=1,Nact
            i = act_to_all(ii) 
            xx(i) = xx(i) - fak*dp(ii)
            anmono(i) = exp(xx(i))*kT1
          enddo
          limit = fak

        else
          !--- re-scale ---
          do ii=1,Nact
            i = act_to_all(ii) 
            dp(ii) = dp(ii)*xscale(i)
          enddo  
          fak = 5.d0
          limit = 1.d0                                   ! limit step, keep direction
          converge(it) = 0.d0
          Nconv = 0
          if (verbose>0) txt = ""
          do i=1,nel
            if (.not.eact(i)) then
              Nconv = Nconv+1 
              if (verbose>0) txt = trim(txt)//" "//catm(i)
            else 
              ii = all_to_act(i) 
              delp = -dp(ii)/(anmono(i)*kT)              ! relative change dx/x
              conv(it,i) = delp
              converge(it) = MAX(converge(it),ABS(delp))
              if (ABS(delp)<finish) then
                Nconv = Nconv+1 
                if (verbose>0) txt = trim(txt)//" "//catm(i)
              endif  
              if (1.d0+delp>fak) then
                limit = min(limit,(fak-1.d0)/delp)       ! such that xnew=xold*fac 
              else if (1.d0+delp<1.d0/fak) then
                limit = min(limit,(1.d0/fak-1.d0)/delp)  ! such that xnew=xold/fac
              endif
            endif  
          enddo
          if (it<=10) then
            limit = 1.d0
          else
            dp = dp*limit
          endif  
          if (verbose>1.and.it==0) then
            write(*,*) 
            print'(7x,A14,A14,A14)',"patom","dpatom","badness" 
            do ii=1,Nact
              i = act_to_all(ii) 
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >             -dp(ii)/(anmono(i)*kT),badness(i)
            enddo
          endif

*         ! apply limited NR step
*         =======================
          !fak = 1.d0+4.d0*EXP(-(MAX(0,it-20))/13.d0)
          do ii=1,nact
            i = act_to_all(ii)
            delp = -dp(ii)*kT1
            nold = anmono(i)
            anmono(i) = MAX(nold/fak,MIN(nold*fak,nold+delp))
          enddo
          if (it>itmax-10) then
            verbose = 2
            do ii=1,Nact
              i = act_to_all(ii) 
              delp = -dp(ii)/(anmono(i)*kT)
              print'(A3,99(1pE12.4))',catm(i),anmono(i),delp
            enddo  
          endif  

        endif
        crit = MAXVAL(converge(MAX(0,it-1):it))
        if (verbose>1) print'(3(i3),2(1pE9.1)," converged(",i2,"):",
     >                 A50)',it,Nact,ind,converge(it),limit,Nconv,txt
        if (it==itmax) then 
          write(*,*) '*** keine Konvergenz in SMCHEM8!'
          write(*,*) 'it, converge, ind =',it,converge(it),limit
          write(*,*) '  n<H>, T =',anhges,Tg
          write(*,*) 'from_merk,NewFastLevel=',from_merk,NewFastLevel
          if (ifatal<=1) then
            chemiter  = chemiter + it
            from_merk = .false.
            badness   = 1.Q0
            pcorr     = 1.Q0
            ifatal    = ifatal+1
            if (ifatal==2) verbose=2
            goto 100        ! try again from scratch before giving up
          endif  
          goto 1000
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
          atmax = 0.d0 
          e = 0
          do i=1,nel
            xx(i) = LOG(anmono(i)*kT)
            if (redo(i)) cycle   
            atfrac = anmono(i)/anHges
            if (atfrac>1.d-30) cycle   
            if (atfrac<atmax) cycle
            atmax = atfrac
            e = i
          enddo  
          if (e==0) exit
          redo(e) = .true.
          psc = xx(e)
          coeff(:) = 0.d0
          lmin = +99
          lmax = -99
          do i=1,nml
            pmol = g(i)
            affect = .false.
            do j=1,m_kind(0,i)
              nu = m_anz(j,i)   
              if (m_kind(j,i)==e) then
                l = nu
                affect = .true.
                pmol = pmol + l*psc
              else
                pmol = pmol + nu*xx(m_kind(j,i))
              endif  
            enddo
            if (.not.affect) cycle
            lmin = MIN(lmin,l)
            lmax = MAX(lmax,l)
            coeff(l) = coeff(l)+EXP(pmol)
          enddo
          psc = EXP(psc)
          if (verbose>1) then
            print'("redo rare element patom dp/patom ",A2)',catm(e)
            print*,"psc=",psc
            print'(" coeff(",I2,":",I2,")=",99(1pE10.2))',
     >                   lmin,lmax,coeff(lmin:lmax)
          endif  
          pat = anmono(e)*kT/psc
          do piter=1,99
            f  = pat*psc-eps(e)*anHges*kT
            fs = psc
            do l=lmin,lmax
              if (coeff(l)==0.d0) cycle
              f  = f  + coeff(l)*l*pat**l
              fs = fs + coeff(l)*l**2*pat**(l-1)
            enddo
            delta = f/fs
            if (verbose>1) print'(1x,A2,1pE25.15,1pE10.2)',
     >                            catm(e),pat*psc,delta/pat
            pat = pat-delta
            if (ABS(delta)<finish*ABS(pat)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** no convergence in post-it "//catm(e)
            write(*,*) anmono(e),eps(e)
            write(*,*) psc,lmin,lmax
            write(*,*) anHges,Tg
            write(*,*) coeff(lmin:lmax)
            goto 1000
          endif  
          anmono(e) = pat*psc*kT1
        enddo
*
*       ! how bad was the initial guess?
*       ================================
        if (.not.from_merk) then
          if (verbose>1) print'(7x,3(A14))',
     >                   "patom","conv.","init.guess"
          do i=1,nel
            badness(i) = anmono(i)/ansave(i)
            switch = switchoff(i)
            if (switch==0) switch=it-1
            if (verbose>1) then
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >              conv(switch,i),badness(i)
            endif  
          enddo
!$omp critical(fort99)
          ilauf = ilauf+1
          !if (ilauf==1) write(99,'(A9,A10,A4,99(A10))') 
     >    !      'Tg','n<H>','it',catm(1:nel)
          !write(99,'(0pF9.3,1pE10.3,I4,99(1pE10.3))') 
     >    !      Tg,anHges,it,badness
!$omp end critical(fort99)
        endif  

      endif     ! NewFullIt

*     ! final anmol determination
*     ===========================
      do i=1,nel
        xx(i) = LOG(anmono(i)*kT)
      enddo  
      amerk = anmono/anHges
      anmol = 0.d0
      do i=1,nml
        pmol = g(i)
        do j=1,m_kind(0,i)
          pat = anmono(m_kind(j,i))*kT
          pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
        enddo
        anmol(i) = EXP(pmol)*kT1
      enddo
      if (charge) pel=anmono(el)*kT

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
        do e=1,nel
          if (e==el) cycle 
          soll  = anHges * eps(e)
          haben = nges(e)
          abw   = ABS(soll-haben)/MAX(soll,haben)
          if (abw>1.d-5) then
            if (verbose>1) then
              print'("*** element conservation error ",A2)',catm(e)
              print'(A12,1pE14.7)',catm(e),anmono(e)/(eps(e)*anHges)
            endif  
            sum = anmono(e)/(eps(e)*anHges)
            do i=1,nml
              do j=1,m_kind(0,i)
                j1 = m_kind(j,i)
                cc = m_anz(j,i)*anmol(i)/(eps(e)*anHges)
                if (j1==e.and.cc>1.d-7) then
                  if (verbose>1) print'(A12,1pE14.7)',cmol(i),cc
                  sum = sum+cc
                endif  
              enddo
            enddo
            if (verbose>1) then
              print'(3(1pE14.7))',soll/anHges,haben/anHges,sum
            endif  
            from_merk = .false.
            ansave = anmono
            verbose=2
            goto 200
          endif
        enddo
      endif
      
      if (verbose.gt.0) print '("  ==> smchem used it=",I3,
     &                          " conv=",1pE9.2)',it,crit
      if (verbose.gt.0) print'("number of pre-iterations",I4,
     &     " -- used saved initial guesses",0pF5.1,"%")',
     &     NpreIt,REAL(Ntaken)/REAL(Nestim+1)*100
      if (verbose.gt.1) read(*,'(a1)') char

!$omp critical(counters)
      chemcall = chemcall + 1
      chemiter = chemiter + it
      preIter  = preIter  + NpreIt
      preUse   = preUse   + Ntaken
      preEst   = preEst   + Nestim
      DUALcorr = DUALcorr + DUALco
      HCOcorr  = HCOcorr  + HCOco
!$omp end critical(counters)

      return

 1000 continue
      open(unit=12,file='fatal.case')
      do i=1,nel
        write(12,'(A2,1x,0pF22.16)') catm(i),12.d0+log10(eps(i))
      enddo  
      write(12,*) anhges,Tg
      close(12)
      stop "***  giving up."


      CONTAINS       ! internal functions - not visible to other units 
************************************************************************
      REAL*8 FUNCTION gk(i)
************************************************************************
*****  returns  ln(kp) [cgs] for different fit formula             *****
************************************************************************
      use CHEMISTRY,ONLY: a,th1,th2,th3,th4,TT1,TT2,TT3,fit,natom,cmol,
     >                    NELEM,elnum,b_nasa,c_nasa
      implicit none
      integer,intent(in) :: i          ! index of molecule
      real*8,parameter :: bar=1.d+6, atm=1.013d+6, Rcal=1.987d+0
      real*8,parameter :: Rgas=8.3144598d+0
      real*8,parameter :: ln10=LOG(10.d0)
      real*8,parameter :: lnatm=LOG(atm), lnbar=LOG(bar)
      real*8 :: lnk,dG
      real*8 :: h_rt,s_r               !Added by Yui Kawashima
      real*8 :: dG_rt_ref(NELEM),dG_rt !Added by Yui Kawashima
      integer:: k,j                    !Added by Yui Kawashima

      if (i.eq.0) then
        gk = -1000.d0                  ! tiny kp for unassigned molecules
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

      else if (fit(i).eq.7) then
        !-----------------------------------------------------
        ! ***  NASA polynomial fit added by Yui Kawashima  ***
        !-----------------------------------------------------         
        if (TT1>1.d3) then
          h_rt = a(i,0) + a(i,1)*TT1/2.d0 
     &         + a(i,2)*TT1**2/3.d0 + a(i,3)*TT1**3/4.d0
     &         + a(i,4)*TT1**4/5.d0 + a(i,5)/TT1
          s_r  = a(i,0)*log(TT1) + a(i,1)*TT1
     &         + a(i,2)*TT1**2/2.d0 + a(i,3)*TT1**3/3.d0
     &         + a(i,4)*TT1**4/4.d0 + a(i,6)
        else
          h_rt = a(i,7) + a(i,8)*TT1/2.d0
     &         + a(i,9) *TT1**2/3.d0 + a(i,10)*TT1**3/4.d0
     &         + a(i,11)*TT1**4/5.d0 + a(i,12)/TT1
          s_r  = a(i,7)*log(TT1) + a(i,8)*TT1
     &         + a(i,9) *TT1**2/2.d0 + a(i,10)*TT1**3/3.d0
     &         + a(i,11)*TT1**4/4.d0 + a(i,13)           
        endif
        dG_rt = h_rt - s_r
        do k=1,m_kind(0,i)
          j = elnum(m_kind(k,i))
          if (c_nasa(j)==0) then
            print*,"Provide the data in data/Burcat_ref-elements.dat"
     &            ," and edit nasa_polynomial.f for "
     &            ,trim(catm(m_kind(k,i)))
            stop
          else
            if (TT1>1.d3) then
              h_rt = b_nasa(j,0) + b_nasa(j,1)*TT1/2.d0 
     &             + b_nasa(j,2)*TT1**2/3.d0
     &             + b_nasa(j,3)*TT1**3/4.d0
     &             + b_nasa(j,4)*TT1**4/5.d0 + b_nasa(j,5)/TT1
              s_r  = b_nasa(j,0)*log(TT1) + b_nasa(j,1)*TT1
     &             + b_nasa(j,2)*TT1**2/2.d0
     &             + b_nasa(j,3)*TT1**3/3.d0
     &             + b_nasa(j,4)*TT1**4/4.d0 + b_nasa(j,6)
            else
              h_rt = b_nasa(j,7) + b_nasa(j,8)*TT1/2.d0
     &             + b_nasa(j,9) *TT1**2/3.d0
     &             + b_nasa(j,10)*TT1**3/4.d0
     &             + b_nasa(j,11)*TT1**4/5.d0
     &             + b_nasa(j,12)/TT1
              s_r  = b_nasa(j,7)*log(TT1) + b_nasa(j,8)*TT1
     &             + b_nasa(j,9) *TT1**2/2.d0
     &             + b_nasa(j,10)*TT1**3/3.d0
     &             + b_nasa(j,11)*TT1**4/4.d0 + b_nasa(j,13)           
            endif
            dG_rt_ref(j) = h_rt - s_r
            dG_rt = dG_rt - m_anz(k,i)*dG_rt_ref(j)
          endif
        enddo
        dG = -dG_rt
        lnk = dG + (1-Natom(i))*lnbar

      else if (fit(i).eq.8) then
        !-----------------------------------------------
        ! ***  NASA 7-polynomial fit for BURCAT data ***
        !-----------------------------------------------
        if (TT1>1.d3) then
          dG = a(i,0) *(log(TT1)-1.d0) 
     &       + a(i,1) *TT1   / 2.d0
     &       + a(i,2) *TT1**2/ 6.d0    
     &       + a(i,3) *TT1**3/12.d0
     &       + a(i,4) *TT1**4/20.d0    
     &       - a(i,5) /TT1  
     &       + a(i,6)
        else
          dG = a(i,7) *(log(TT1)-1.d0) 
     &       + a(i,8) *TT1   / 2.d0
     &       + a(i,9) *TT1**2/ 6.d0    
     &       + a(i,10)*TT1**3/12.d0
     &       + a(i,11)*TT1**4/20.d0    
     &       - a(i,12)/TT1  
     &       + a(i,13)
        endif
        lnk = dG + (1-Natom(i))*lnbar

      else
        print*,cmol(i),"i,fit=",i,fit(i)
        stop "???"
      endif  

      gk = lnk
      end FUNCTION gk

      end SUBROUTINE smchem8
