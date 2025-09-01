*********************************************************************
      SUBROUTINE SUPERSAT(T,nat,nmol,Sat)
*********************************************************************
      use CHEMISTRY,ONLY: NMOLE,cmol
      use DUST_DATA,ONLY: NELEM,NDUST,bk,atm,rgas,bar,fit,cfit,
     &                    dust_nam,dust_nel,dust_el,dust_nu,dust_mass,
     &                    is_liquid,Tcorr,elnam,Nfit,Tfit,Bfit,fitTmax
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: T
      real(kind=qp),intent(in) :: nat(NELEM),nmol(NMOLE)
      real(kind=qp),intent(out):: Sat(NDUST)
      real(kind=qp),parameter :: cal=4.184Q+0       ! 1 cal in J
      real(kind=qp),parameter :: mmHg=1.3328Q+3     ! 1 mmHg in dyn/cm2
      real(kind=qp),parameter :: Joule=1.Q+7        ! 1 Joule in erg
      real(kind=qp),parameter :: mol=6.02214076Q+23 ! 1 mol
      real(kind=qp),parameter :: eV=1.60218Q-12     ! 1 eV in erg
      real(kind=qp) :: T1,T2,T3,TC,kT,RT,dG,lbruch,pst,psat,dGRT
      real(kind=qp) :: a(0:4),term,n1,natom,aa(0:6)
      real(kind=qp) :: TT1,TT2,TT3,ffff,dfdT
      !real(kind=qp) :: tiny16=TINY(T1),huge16=HUGE(T1)
      integer :: i,j,l,STINDEX,el,imol,imol1,imol2,ifit
      character(len=20) :: search,leer='                    '

      T1  = T
      T2  = T1**2
      T3  = T1**3
      kT  = bk*T1
      RT  = rgas*T1
      Sat = 0.Q0
      do i=1,NDUST
        a(:) = cfit(i,:)
        if (fit(i)==1) then
          !-------------------------------------
          !***  dG-fit Sharp & Huebner 1990  ***
          !-------------------------------------
          pst = atm
          dG = a(0)/T1 + a(1) + a(2)*T1 + a(3)*T2 + a(4)*T3
          dG = dG/(rgas/cal*T1)
          lbruch = 0.Q0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
          enddo
          Sat(i) = EXP(lbruch-dG)

        else if (fit(i)==2) then
          !-----------------------------------
          !***  dG polynom-fit NIST-Janaf  ***
          !-----------------------------------
          pst = bar
          dG = a(0)/T1 + a(1) + a(2)*T1 + a(3)*T2 + a(4)*T3
          dGRT = dG/RT
          lbruch = 0.Q0
          Natom = 0.0          
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
            Natom  = Natom + dust_nu(i,j)            
          enddo
          Sat(i) = EXP(lbruch-dGRT)
          !print*,dust_nam(i),T1,-dG*Joule/mol/eV/Natom,      ! eV/atom
     >    !                      -dG*Joule/mol/dust_mass(i)   ! erg/g          

        else if (fit(i)==3) then
          !-----------------------------------------
          !***  ln(pvap) polynom-fit NIST-Janaf  ***
          !-----------------------------------------
          psat = EXP(a(0)/T1 + a(1) + a(2)*T1 + a(3)*T2 + a(4)*T3)
          if (dust_nel(i)==1.and.dust_nu(i,1)==1) then
            n1 = nat(dust_el(i,1))
          else
            search = trim(dust_nam(i))
            call upper(search)
            l = index(search,'[')
            search = search(1:l-1)//leer(l:20)
            imol = STINDEX(cmol,NMOLE,search)
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            n1 = nmol(imol)
          endif
          Sat(i) = n1*kT/psat

        else if (fit(i)==4) then
          !----------------------------------------
          !***  ln(pvap) 3-para-fit NIST-Janaf  ***
          !----------------------------------------
          psat = EXP(a(0) + a(1)/(T1 + a(2)))
          if (dust_nel(i)==1.and.dust_nu(i,1)==1) then
            n1 = nat(dust_el(i,1))
          else
            search = trim(dust_nam(i))
            call upper(search)
            l = index(search,'[')
            search = search(1:l-1)//leer(l:20)
            imol = STINDEX(cmol,NMOLE,search)
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            n1 = nmol(imol)
          endif
          Sat(i) = n1*kT/psat

        else if (fit(i)==5) then
          !--------------------------
          !***  -dG/RT Stock-fit  ***
          !--------------------------
          pst = bar
          dGRT = a(0)/T1 + a(1)*LOG(T1) + a(2) + a(3)*T1 + a(4)*T2
          if (T1>fitTmax(i)) then
            TT1  = fitTmax(i)
            TT2  = TT1*TT1
            TT3  = TT1*TT2
            !-- we assume that dGRT is linear beyond Tmax --
            !ffff = a(0)/TT1 + a(1)*LOG(TT1) + a(2) + a(3)*TT1 + a(4)*TT2
            !dfdT =-a(0)/TT2 + a(1)/TT1 + a(3) + 2.0*a(4)*TT1
            !dGRT = ffff + dfdT*(T1-TT1)
            !-- we assume that dG = dGRT*T is linear beyond Tmax 
            ffff = a(0) + a(1)*LOG(TT1)*TT1 + a(2)*TT1 + a(3)*TT2
     >                                                 + a(4)*TT3
            dfdT = a(1) + a(1)*LOG(TT1) + a(2) + 2*a(3)*TT1
     >                                         + 3*a(4)*TT2
            dGRT = (ffff + dfdT*(T1-TT1))/T1
          endif
          lbruch = 0.Q0
          Natom = 0.0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
            Natom  = Natom + dust_nu(i,j)
          enddo
          Sat(i) = EXP(lbruch+dGRT)
          !dG = dGRT*RT                                      ! Joule/mol
          !print*,dust_nam(i),T1,dG*Joule/mol/eV/Natom,      ! eV/atom
     >    !                      dG*Joule/mol/dust_mass(i)   ! erg/g

        else if (fit(i)==6) then
          !---------------------------------------------------------------
          !***  Yaws' Chemical Properties Handbook (McGraw-Hill 1999)  ***
          !---------------------------------------------------------------
          psat = a(0) + a(1)/T1 + a(2)*LOG10(T1) + a(3)*T1 + a(4)*T2
          psat = 10.Q0**psat * mmHg
          if (dust_nel(i)==1.and.dust_nu(i,1)==1) then
            n1 = nat(dust_el(i,1))
          else
            search = trim(dust_nam(i))
            call upper(search)
            l = index(search,'[')
            search = search(1:l-1)//leer(l:20)
            imol = STINDEX(cmol,NMOLE,search)
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            n1 = nmol(imol)
          endif
          Sat(i) = n1*kT/psat

        else if (fit(i)==7) then
          !---------------------
          !***  BURCAT data  ***
          !---------------------
          do ifit=1,Nfit(i)
            if (T1<Tfit(i,ifit+1)) exit
            if (ifit==Nfit(i)) exit
          enddo
          if (T1>1000.0) then
            aa(0:6) = Bfit(i,ifit,1:7)
          else
            aa(0:6) = Bfit(i,ifit,8:14)
          endif  
          dGRT = aa(0) *(log(T1)-1.d0) 
     &         + aa(1) *T1   / 2.d0
     &         + aa(2) *T1**2/ 6.d0    
     &         + aa(3) *T1**3/12.d0
     &         + aa(4) *T1**4/20.d0    
     &         - aa(5) /T1  
     &         + aa(6)
          lbruch = 0.Q0
          Natom = 0.0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
            Natom  = Natom + dust_nu(i,j)
          enddo
          Sat(i) = EXP(lbruch+dGRT)
          !print*,dust_nam(i),"BURCAT:",T1,ifit,dGRT,Sat(i)
          if (T1>1.2*Tfit(i,Nfit(i)+1)) then
            Sat(i) = 1.Q-999   ! BURCAT-fits go crasy beyond fit-regime
          endif  
          
        else if (fit(i)==99) then
          !-----------------------
          !***  special cases  ***
          !-----------------------
          if (dust_nam(i).eq.'H2O[l]') then
            !--- Ackerman & Marley 2001 ---
            TC   = MIN(2000.0,T1)-273.15         ! T[degree Celsius]
            psat = 6112.1*exp((18.729*TC - TC**2/227.3)/(TC + 257.87))
            imol = STINDEX(cmol,NMOLE,"H2O")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'H2O[s]') then
            !--- Ackerman & Marley 2001 ---
            TC   = MIN(2000.0,T1)-273.15         ! T[degree Celsius]
            psat = 6111.5*exp((23.036*TC - TC**2/333.7)/(TC + 279.82))
            imol = STINDEX(cmol,NMOLE,"H2O")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'NH3[s]') then
            !--- CRC Handbook of Chemistry and Physics (Weast 1971) ---
            psat = exp(10.53 - 2161.Q0/T1 - 86596.Q0/T2)*bar
            imol = STINDEX(cmol,NMOLE,"NH3")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'CH4[s]') then
            !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
            psat = 10.0**(3.9895 - 443.028/(T1-0.49))*bar
            imol = STINDEX(cmol,NMOLE,"CH4")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'NH4SH[s]') then
            !--- G.Lee's fit to Walker & Lumsden (1897) ---
            psat  = 10.0**(7.8974 - 2409.4/T1)*bar
            imol1 = STINDEX(cmol,NMOLE,"NH3")
            imol2 = STINDEX(cmol,NMOLE,"H2S")
            if (imol1<=0.or.imol2<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = SQRT(nmol(imol1)*kT/(psat*0.5))
     &             * SQRT(nmol(imol2)*kT/(psat*0.5))

          else if (dust_nam(i).eq.'H2S[s]') then
            !--- Stull (1947) ---
            if (T1 < 30.0) then ! Limiter for very cold T
              T1 = 30.0
            end if
            if (T1 < 212.8) then
              psat = 10.0**(4.43681 - 829.439/(T1-25.412))*bar
            else
              psat = 10.0**(4.52887 - 958.587/(T1-0.539))*bar
            end if
            imol = STINDEX(cmol,NMOLE,"H2S")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'S2[s]') then
            !--- Zahnle et al. (2016) ---
            psat = exp(27.0 - 18500.0/T1)*bar
            !if (T1 < 413.0) then
            !  psat = exp(27.0 - 18500.0/T1)*bar
            !else
            !  psat = exp(16.1 - 14000.0/T1)*bar
            !end if
            !--- Lyons 2008 ---
            !write(50,'(F8.1,2(1pE13.4))')
     &      !         T1,psat,10.0**(7.024 - 6091.0/T1)*bar
            !psat = 10.0**(7.024 - 6091.0/T1)*bar
            imol = STINDEX(cmol,NMOLE,"S2")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'S2[l]') then
            !--- Zahnle et al. (2016) ---
            psat = exp(16.1 - 14000.0/T1)*bar
            imol = STINDEX(cmol,NMOLE,"S2")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'S8[s]') then
            !--- Zahnle et al. (2016) ---
            psat = exp(20.0 - 11800.0/T1)*bar
            !if (T1 < 413.0) then
            !  psat = exp(20.0 - 11800.0/T1)*bar
            !else
            !  psat = exp(9.6 - 7510.0/T1)*bar
            !end if
            !--- Lyons 2008 ---
            !psat = 10.0**(4.188 - 3269.0/T1)*bar
            !Zahnle et al.(1989)
            psat = 1.316E-3*exp(4.969 - 2201.0/T1)*atm
            imol = STINDEX(cmol,NMOLE,"S8")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'S8[l]') then
            !--- Zahnle et al. (2016) ---
            psat = exp(9.6 - 7510.0/T1)*bar
            imol = STINDEX(cmol,NMOLE,"S8")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else
            print*,"*** supersat.f fit=",fit(i)," ??? ",dust_nam(i)
            stop
          endif

        else
          print*,"*** supersat.f fit=",fit(i)," ??? ",dust_nam(i)
          stop
        endif

        if (Tcorr(i)>0.0) then
          if (is_liquid(i).and.T1<Tcorr(i)) then
            Sat(i) = Sat(i)/EXP(0.1*(Tcorr(i)-T1))
            Sat(i) = MIN(Sat(i),0.99*Sat(i-1))
          endif
        endif
        if (i>1) then
          if (Tcorr(i-1)>0.0) then
            if ((.not.is_liquid(i-1)).and.T1>Tcorr(i-1)) then
              Sat(i-1) = Sat(i-1)/EXP(0.1*(T1-Tcorr(i)))
              Sat(i-1) = MIN(Sat(i-1),0.99*Sat(i))
            endif
          endif
        endif

        !Sat(i) = MAX(Sat(i),tiny16) 
        !Sat(i) = MIN(Sat(i),huge16) 

      enddo

      RETURN
      end
