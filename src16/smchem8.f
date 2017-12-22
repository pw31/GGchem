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
     >                    NewFastLevel,nml=>NMOLE,nel=>NELM,cmol,catm,
     >                    m_kind,m_anz,charge,elion,el,
     >                    th1,th2,th3,th4,TT1,TT2,TT3
      use EXCHANGE,ONLY: chemcall,chemiter
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
*  Die Variable "alle" entscheidet, ob die nicht unmittelbar 
*  beruecksichtigten Molekuele dennoch inkonsistent mitgerechnet werden. 
      logical,parameter :: alle=.true.
*-----------------------------------------------------------------------
*  Bei merk=.true. merkt sich die Routine die letzte konvergiert Loesung
*  und geht beim naechsten Mal von diesen Startwerten aus.
      logical,intent(INOUT) :: merk
*-----------------------------------------------------------------------
*  Die Variable "ngestst" entscheidet, ob die Elementerhaltung ueber-
*  prueft werden soll.
      logical,parameter :: ngestst=.false.
*-----------------------------------------------------------------------
*  Die Variable tdispol bestimmt, welches die niedrigste Temperatur 
*  ist, die in die Dissoziationspolynome eingesetzt werden darf.
      real*8,parameter :: tdispol=300.d0
*-----------------------------------------------------------------------
      integer stindex,info,ipvt(nel),Nconv,switch,ido,iredo
      integer Nact,all_to_act(nel),act_to_all(nel),switchoff(nel)
      integer e,i,j,j1,ii,jj,kk,l,it,m1,m2,piter,ifatal,ipull,pullmax
      integer Nseq,imin,imax,enew,eseq(nel)
      integer,parameter :: itmax=200,Ncmax=16
      real*8,parameter :: finish=1.d-12
      real*8 :: qual0,qual
      real*8 :: g(0:nml),limit
      real*8 :: work(nel*(nel+1))
      integer:: ind,indx(nel)
      real*8 :: condnum1,work2(nel)
      real*8 :: kT,kT1,cc,nelek,ng,Sa,fak,lth,arg,term,f,fs
      real*8 :: pel,delta,pat,atfrac,atmax
      real*8 :: nges(nel),pmono1(nel),coeff(-1:Ncmax)
      real*8 :: DF(nel,nel),dp(nel),FF(nel),pmol,crit
      real*8 :: DF0(nel,nel),FF0(nel),scale(nel),conv(0:500,nel)
      real*8 :: converge(0:500),delp,nold,null(nel),nsave(nel)
      real*8 :: soll,haben,abw,sum
      real*8 :: pbefore(nel),norm(nel)
      real*8 :: emax,pges,pwork
      logical :: from_merk,eact(nel),redo(nel),done(nel),affect,known
      logical :: relevant(nml)
      logical :: ptake,isOK,IS_NAN
      character(len=5000) :: mols
      character(len=100) :: txt
      character(len=1) :: char
      integer,save :: TiC,ilauf=0
      real*8,allocatable,save :: amerk(:),ansave(:)
      real*8,allocatable,save :: badness(:),pcorr(:,:) 
      integer,allocatable,save :: pkey(:)
!$omp threadprivate(TiC,amerk,ansave,badness,pcorr,pkey)
*-----------------------------------------------------------------------      

      ifatal = 0
      if (.not.allocated(amerk)) then
        allocate(badness(nel),pcorr(nel,nel),pkey(nel),
     >           amerk(nel),ansave(nel))
        TiC = stindex(cmol,nml,'TIC    ')
        badness = 1.d0
        pcorr   = 1.d0
        pkey    = 0
      endif

*-----------------------------------------------------------------------
*     ! zu niedrige Temperaturen abfangen und
*     ! Variable fuer die Dissoziationskonstanten berechnen
*     =====================================================
      TT1 = MAX(tdispol,Tg)
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
        if (i.ne.TiC) g(i)=gk(i)       ! compute all equil.constants
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
      g(TiC) = EXP(MIN(700.d0,-2.30256*arg))

*---------------------------------------------------------------------------
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
      from_merk = .false.
      if (charge) then
        nelek = 0.Q0 
        do i=1,nel
          if (i==el) cycle 
          ng = anHges * eps(i)
          Sa = g(elion(i))*kT1
          nelek = nelek + ng/(0.5d0 + SQRT(0.25d0 + ng/Sa))
        enddo
        anmono(el) = nelek
        pel = nelek*kT
        !peest = pel
        !pel = pecorr*pel
        !anmono(el) = pecorr*anmono(el) 
        if (verbose>1) print'(" estimate pel=",1pE10.3)',pel
      endif  

*-----------------------------------------------------------------------
*     ! estimate atomic pressures: new method
*     =======================================
      Nseq = nel
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
          norm(e) = eps(e)
          if (e==el) norm(e)=anmono(el)/anHges
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
        pwork = pges
        !-------------------------------------------
        ! store coeff for Sum_l coeff(l) p^l = pges 
        !-------------------------------------------
        coeff(:) = 0.d0   
        if (verbose>0) mols = ''
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
          if (verbose>0) mols = trim(mols)//" "//cmol(i)
          coeff(l) = coeff(l) + l*pmol
          !------------------------------------
          ! for initial guess, consider this 
          ! molecule to have all of element e2 
          !------------------------------------
          if (pmol>0.d0.and.l>0) then
            pwork = MIN(pwork,(pges/(l*pmol))**(1.d0/REAL(l)))
            !if (verbose>1) print'(A10,1pE10.3)',cmol(i),pwork
          endif  
        enddo  
        if (verbose>1) print*,trim(mols)
        if (enew==el) then
          pel = SQRT(-coeff(-1)/(1.Q0+coeff(+1)))     ! 0 = pel - a/pel + b*pel
          anmono(el) = pel*kT1
        else   
          !----------------------------------------------
          ! solve 1d equation above with Newton's method 
          !----------------------------------------------
          do piter=1,99                  
            f  = pwork-pges
            fs = 1.d0
            do l=1,Ncmax
              if (coeff(l)==0.d0) cycle
              f  = f  + coeff(l)*pwork**l
              fs = fs + coeff(l)*l*pwork**(l-1)
            enddo
            if (fs==0.d0) stop "*** fs=0 in smchem8 1d-pre-it."
            delta = f/fs
            pwork = pwork-delta
            if (verbose>1) print'(A2,I3,1pE25.15,1pE10.2)',
     >                     catm(enew),piter,pwork,delta/pwork
            if (ABS(delta)<1.E-4*ABS(pwork)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** smchem8 no conv in 1D pre-it "//catm(enew)
            write(*,*) coeff
            goto 1000
          endif  
          anmono(enew) = pwork*kT1
        endif  

        !-----------------------------------------------------------
        ! take into account feedback on elements considered before,
        ! unless they are much more abundant, with Newton-Raphson 
        !-----------------------------------------------------------
        nsave = anmono
 150    continue
        eact(:) = .false.
        Nact = 0
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
            else
              pcorr(enew,e) = 1.d0 
            endif  
          endif
        enddo
        if (verbose>1) then
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
          print'("corr",99(1pE11.2))',pcorr(enew,act_to_all(1:Nact))
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
        qual  = 9.d+99
        dp(:) = 0.d0
        fak   = 1.d0
        null  = anmono
        do it=1,99
          qual0 = qual 
          pullmax = 1
          if (it>30) pullmax=10
          do ipull=1,pullmax    ! pullback if quality gets worse
            !--- make a step ---
            do ii=1,Nact
              i = act_to_all(ii)
              anmono(i) = null(i)-fak*dp(ii)*kT1
            enddo  
            !--- determine new FF and DF ---
            do ii=1,Nact
              i = act_to_all(ii) 
              FF(ii) = anHges*eps(i)*kT - anmono(i)*kT
              scale(i) = anmono(i)  
              DF(ii,:) = 0.d0
              DF(ii,ii) = -scale(i)
              pmono1(i) = scale(i) / (anmono(i)*kT)
            enddo
            do i=1,nml
              if (.not.relevant(i)) cycle
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
            enddo
            !--- determine new quality ---
            qual = 0.d0
            do ii=1,Nact
              i = act_to_all(ii)           
              qual = qual + (FF(ii)/(anHges*norm(i)*kT))**2
            enddo  
            if (qual<qual0) exit
            if (ipull==pullmax) exit
            if (verbose>1) print'("pullback",3(1pE11.3))',fak,qual0,qual
            fak = 0.5*fak   ! reduce NR-step
          enddo  
          if (verbose>1) print'(I4,99(1pE11.3))',
     >                   it,anmono(act_to_all(1:Nact))*kT,qual
          if (it>1.and.qual<1.d-4) exit
          !--- determine new NR-vector ---
          call GAUSS8(nel,Nact,DF,dp,FF)
          do ii=1,Nact
            i = act_to_all(ii) 
            dp(ii) = dp(ii)*scale(i)
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
        pkey(ido) = enew
        pcorr(enew,:) = 1.Q0
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
        do i=1,nml
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
        pel = SQRT(coeff(-1)/(1.d0+coeff(+1)))     ! 0 = pel - a/pel + b*pel
        anmono(el) = pel/kT
        !if (verbose>1) print'(" pecorr =",3(1pE10.3))',pecorr,pel/peest
        !pecorr = pel/peest
      endif  

*     ! use memory of deviations between predicted atom pressures 
*     ! and converged atom pressures to improve the predictions
*     ============================================================
      ansave = anmono
      if (NewFastLevel<2.and.ptake) anmono = anmono*badness
*     
*-----------------------------------------------------------------------
 200  continue

      if ( alle ) then
        ! alle Molekuele mitrechnen
*       ===========================
        do i=1,nml
          if (g(i)>1.d+300) then
            print'("huge kp",A12,0pF11.3,1pE12.3E4)',cmol(i),Tg,g(i)
            stop
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
      if (NewFullIt) then
*       ! Jacobi matrix and rhs vector for Newton-Raphson
*       =================================================
        it = 0
        eact(:) = .true.
        conv(:,:) = 9.d+99
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
          DF(ii,:)  = 0.d0
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
            ii = all_to_act(m1)
            term   = m_anz(j,i) * pmol
            FF(ii) = FF(ii) - term
            do l=1,m_kind(0,i)
              m2 = m_kind(l,i)
              if (.not.eact(m2)) cycle
              jj = all_to_act(m2)
              DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term*pmono1(m2)
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
        !--- re-scale ---
        isOK = .true.
        do ii=1,Nact
          i = act_to_all(ii) 
          dp(ii) = dp(ii)*scale(i)
          if (IS_NAN(dp(ii))) then
            print*,catm(i),eps(i),anmono(i)/anHges,dp(ii)/kT/anHges,
     >             scale(i),FF0(ii)
            isOK = .false.
          endif  
        enddo  
        if (.not.isOK) then 
          print*,Nact,ind 
          stop "*** NaN in smchem8"
        endif  

*       ! limit NR-step and check convergence
*       =====================================
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
          print'(7x,A14,A14,A14)',"natom","dnatom","badness" 
          do ii=1,Nact
            i = act_to_all(ii) 
            print'(A7,3(1pE14.6))',catm(i),anmono(i),
     >           -dp(ii)/(anmono(i)*kT),badness(i)
          enddo
        endif
        
*       ! apply limited NR step
*       =======================
        !fak = 1.Q0+4.Q0*EXP(-(MAX(0,it-20))/13.Q0)
        do ii=1,nact
          i = act_to_all(ii)
          delp = -dp(ii)*kT1
          nold = anmono(i)
          anmono(i) = MAX(nold/fak,MIN(nold*fak,nold+delp))
        enddo
        if (it>itmax-10) then
          verbose=2
          do ii=1,Nact
            i = act_to_all(ii) 
            print'(A3,2(1pE12.3))',catm(i),
     >           anmono(i),-dp(ii)/(anmono(i)*kT) 
          enddo  
        endif  
        crit = MAXVAL(converge(MAX(0,it-1):it))
        if (verbose>1) print'(i3,i3,2(1pE9.1)," converged(",i2,"):",
     >                    A50)',it,Nact,converge(it),limit,Nconv,txt
        if (it==itmax) then 
          write(*,*) '*** keine Konvergenz in SMCHEM8!'
          write(*,*) 'it, converge, ind =',it,converge(it),limit
          write(*,*) '  n<H>, T =',anhges,Tg
          if (ifatal==0) then
            chemiter  = chemiter + it
            from_merk = .false.
            ifatal  = 1
            verbose = 2             
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
            atfrac = anmono(i)/anHges
            if (redo(i)) cycle   
            if (atfrac>1.d-30) cycle   
            if (atfrac<atmax) cycle
            atmax = atfrac
            e = i
          enddo  
          if (e==0) exit
          redo(e) = .true.
          coeff(:) = 0.d0
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
            fs = 1.d0
            do l=-1,Ncmax
              if (coeff(l)==0.d0) cycle
              f  = f  + coeff(l)*l*pat**l
              fs = fs + coeff(l)*l**2*pat**(l-1)
            enddo
            delta = f/fs
            pat = pat-delta
            if (verbose>1) print'(A2,1pE25.15,1pE10.2)',
     >                            catm(e),pat,delta/pat
            if (ABS(delta)<finish*ABS(pat)) exit 
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
              print'(A7,3(1pE14.6))',catm(i),anmono(i),
     >              conv(switch,i),badness(i)
            endif  
          enddo
!$omp critical(fort99)
          ilauf = ilauf+1
          if (ilauf==1) write(99,'(A9,A10,A4,99(A10))') 
     >          'Tg','n<H>','it',catm(1:nel)
          write(99,'(0pF9.3,1pE10.3,I4,99(1pE10.3))') 
     >          Tg,anHges,it,badness
!$omp end critical(fort99)
        endif  

      endif     ! NewFullIt

*     ! final anmol determination
*     ===========================
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

      if (verbose.gt.1) read(*,'(a1)') char

      chemcall = chemcall + 1
      chemiter = chemiter + it
      return

 1000 continue
      open(unit=12,file='fatal.case')
      do i=1,nel
        write(12,'(A2,1x,0pF30.26)') catm(i),12+log10(eps(i))
      enddo  
      write(12,*) anhges,Tg
      close(12)
      stop "***  giving up."


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
      real*8,parameter :: bar=1.d+6, atm=1.013d+6, Rcal=1.987d+0
      real*8,parameter :: Rgas=8.3144598d+0
      real*8,parameter :: ln10=DLOG(10.d0)
      real*8,parameter :: lnatm=DLOG(atm), lnbar=DLOG(bar)
      integer,intent(in) :: i    ! index of molecule
      real*8 :: lnk,gk,dG            ! return kp in [cgs]
      if (i.eq.0) then
        gk = 1.d-300             ! tiny kp for unassigned molecules
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
        !-----------------
        ! ***  dG-fit  ***
        !-----------------
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
      gk = EXP(MIN(700.d0,lnk))
      end FUNCTION gk

      end SUBROUTINE smchem8
