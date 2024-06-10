************************************************************************
      REAL*8 FUNCTION gk(i)
************************************************************************
*****  returns  ln(kp) [cgs] for different fit formula             *****
************************************************************************
      use CHEMISTRY,ONLY: a,th1,th2,th3,th4,TT1,TT2,TT3,fit,natom,cmol,
     >                    NELEM,elnum,b_nasa,c_nasa,m_kind,m_anz,catm
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
