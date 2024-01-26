**********************************************************************
      SUBROUTINE INIT_DUSTCHEM
**********************************************************************
      use PARAMETERS,ONLY: model_eqcond,phyllosilicates,use_SiO,
     &                     metal_sulphates
      use CHEMISTRY,ONLY: NMOLE,NELM,catm
      use DUST_DATA,ONLY: NDUSTmax,NEPS,NELEM,NDUST,eps0,amu,
     &                    dust_nam,dust_rho,dust_vol,dust_mass,
     &                    dust_nel,dust_nu,dust_el,fit,cfit,
     &                    Nfit,Tfit,Bfit,
     &                    elnr,elcode,elnam,mass,Tmelt,Tcorr,
     &                    DustChem_file
      use EXCHANGE,ONLY: H,Si,Al,Ca
      implicit none
      integer :: i,imax,j,k,el,j1,j2
      real*8 :: dmass,prec(NDUSTmax)
      character(len=10000) :: allcond
      character(len=200):: zeile,lastzeile
      character(len=100) :: trivial(NDUSTmax),tmp
      character(len=2)  :: name
      logical :: found,allfound,hasH,hasSi,hasAl,hasCa

      write(*,*) 
      write(*,*) "reading "//trim(DustChem_file)//" ..."
      write(*,*) "========================"
      trivial(:)=' '

      open(12, file='data/'//trim(DustChem_file), status='old')
 
      write(*,*) '--- dust species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) imax
      read(12,1000) zeile
      allcond = " "
      NDUST = 1
      do i=1,imax
        read(12,1000) zeile
        read(zeile,*) dust_nam(NDUST)
        j1 = index(zeile,' ')
        read(zeile(j1+1:),*) trivial(NDUST)
        if (index(zeile,'[l]')>0) then
          j2 = index(zeile,trim(trivial(NDUST)))
     &       + len(trim(trivial(NDUST)))
          read(zeile(j2+1:),*) Tmelt(NDUST)
          read(zeile(j1+1:j2),*) trivial(NDUST)
        endif
        print*,trim(dust_nam(NDUST)),trim(trivial(NDUST))
        read(12,*) dust_rho(NDUST)
        read(12,*) dust_nel(NDUST)
        dmass = 0.d0
        allfound = .true.
        hasH  = .false.
        hasSi = .false.
        hasCa = .false.
        hasAl = .false.
        do j=1,dust_nel(NDUST)
          read(12,1030) dust_nu(NDUST,j),name
          found = .false. 
          do k=1,NELEM
            if (elnam(k).eq.name) then
              dust_el(NDUST,j) = k
              dmass = dmass + dust_nu(NDUST,j)*mass(k)
              found = .true.
              if (k==H)  hasH =.true.
              if (k==Si) hasSi=.true.
              if (k==Al) hasAl=.true.
              if (k==Ca) hasCa=.true.
            endif
          enddo
          if (.not.found) then
            print*,trim(dust_nam(NDUST)),name
            print*,elnam(1:NELEM)
            stop 'Element in dust species not found'
          endif  
          found = .false.
          do k=1,NELM
            if (catm(k).eq.name) then
              found = .true.
              exit
            endif
          enddo
          if (.not.found) allfound=.false.
        enddo
        found = .false.
        do 
          lastzeile = zeile 
          read(12,1000) zeile
          if (trim(zeile)=='') exit
          if (zeile(1:1)=='#') cycle
          read(zeile,*) fit(NDUST)
          prec(NDUST) = 0.0
          if (fit(NDUST).ne.7) then
            read(zeile,*) fit(NDUST),cfit(NDUST,0:4)
            j1 = index(lastzeile,'+/-')
            j2 = index(lastzeile,':')
            if (j1>0) then
              tmp = lastzeile(j1+3:)
              if (j2>j1) tmp=lastzeile(j1+3:j2-1)            
              read(tmp,*) prec(NDUST)
            endif  
            !print*,trim(tmp),prec(NDUST)
          else
            read(zeile,*) fit(NDUST),Nfit(NDUST),
     >                    Tfit(NDUST,1:Nfit(NDUST)+1)
            do k=1,Nfit(NDUST)
              read(12,*) Bfit(NDUST,k,1:14)
            enddo  
          endif  
          found = .true.
        enddo
        if (.not.found) then
          print*,"*** syntax error in DustChem.dat, condensate=",
     &         dust_nam(NDUST)
          stop
        endif  
        j1 = index(allcond," "//trim(dust_nam(NDUST)))
        if (j1>0) then
          print*,"*** double condensate in DustChem.dat"
          print*,dust_nam(NDUST)
          stop
        endif  
        if ((.not.phyllosilicates).and.hasH
     &       .and.(hasSi.or.hasAl.or.hasCa)) allfound=.false.
        if ((.not.use_SiO).and.
     &       trim(dust_nam(NDUST))=='SiO[s]') allfound=.false.
        if ((.not.metal_sulphates).and.
     &       index(trivial(NDUST),'SULFATE')>0) allfound=.false.   
        if (allfound) then
          dust_mass(NDUST) = dmass
          dust_vol(NDUST) = dmass/dust_rho(NDUST)
          write(*,1060) NDUST,dust_nam(NDUST),dust_rho(NDUST),
     &                  dust_vol(NDUST), (dust_nu(NDUST,j),
     &                  elnam(dust_el(NDUST,j)),j=1,dust_nel(NDUST))
          allcond = " "//trim(allcond)//" "//trim(dust_nam(NDUST))
          NDUST = NDUST+1
        endif
      enddo
      NDUST=NDUST-1
      write(*,*) NDUST," condensed species"
      write(*,*)
      write(*,*) '--- involved elements ---'
      NEPS=0
      elcode(:)=0
      do i=1,NDUST
        do j=1,dust_nel(i)
          name = elnam(dust_el(i,j)) 
          do k=1,NELEM
            if (elnam(k).eq.name) then
              el = k
              exit
            endif
          enddo
          found = .false.           
          do k=1,NEPS
            if (el==elnr(k)) found=.true.
          enddo
          if (.not.found) then
            NEPS = NEPS+1 
            elnr(NEPS) = el
            elcode(el) = NEPS
            write(*,*) elcode(elnr(NEPS)),' ',name,el
          endif
        enddo
      enddo

      Tcorr(:) = -1.d0
      if (model_eqcond) call CHECK_MELTING
      write(*,*)
      
      !open(unit=1,file='condensates.tex')
      !do i=1,NDUST
      !  limit = ' '
      !  j = index(dust_nam(i),"[l]")
      !  if (Tcorr(i)>0.and.j>0) then
      !    write(limit,'("$>$",I4,"K")') int(Tcorr(i))
      !  else if (Tcorr(i)>0) then
      !    write(limit,'("$<$",I4,"K")') int(Tcorr(i))
      !  endif  
      !  if (prec(i)>0.0) then
      !    write(1,3000)
     &!      i,dust_nam(i),trivial(i),dust_rho(i),
     &!      fit(i),limit,cfit(i,0:4),prec(i)
      !  else  
      !    write(1,3001)
     &!      i,dust_nam(i),trivial(i),dust_rho(i),
     &!      fit(i),limit,cfit(i,0:4)
      !  endif  
      !enddo  
      !close(1)

      RETURN 
 1000 format(a200)
 1010 format(a2)
 1020 format(2(l1,1x),i1,1x,a10)
 1030 format(i2,1x,a2)
 1040 format(i2,1x,a10)
 1050 format(1x,a10,i4,' mass=',0pf7.3," amu")
 1060 format(I4,1x,a20," rhod=",0pf6.3," V0=",1pe11.3,2x,
     &       99(i2,"x",A2,1x))
 1070 format(1x,a10,99(i1,"x",i2,1x))
 2011 format(1(I2,1x,a8),22x,'->',I2,1x,a10,99(I2,1x,a8))
 2021 format(2(I2,1x,a8),11x,'->',I2,1x,a10,99(I2,1x,a8))
 2031 format(3(I2,1x,a8)    ,'->',I2,1x,a10,99(I2,1x,a8))
 3000 format(I3," & ",A20," & ",A25," & ",0pF5.2," & ",I2," & ",A8,
     &       " & ",5(1pE12.5," & "),"$\pm$",0pF4.2,"\\")
 3001 format(I3," & ",A20," & ",A25," & ",0pF5.2," & ",I2," & ",A8,
     &       " & ",5(1pE12.5," & "),9x,"\\")
      end 

***********************************************************************
      SUBROUTINE CHECK_MELTING
***********************************************************************
      use CHEMISTRY,ONLY: NMOLE,NELM,catm
      use DUST_DATA,ONLY: qp,NELEM,NDUST,dust_nam,Tmelt,Tcorr,is_liquid
      implicit none
      real*8 :: T
      real(kind=qp) :: nat(NELEM),nmol(NMOLE),Sat(NDUST)
      real(kind=qp) :: old,new,S(NDUST,10000)
      integer :: i,j,k,iT,Ncheck,il,is
      integer :: iliq(NDUST),isol(NDUST)
      character(len=15) :: search

      !--------------------------------------
      ! ***  identify solid/liquid pairs  ***
      !--------------------------------------
      is_liquid(:) = .false.
      Ncheck = 0
      do i=1,NDUST
        k = index(dust_nam(i),'[l]')
        if (k>0) then
          is_liquid(i) = .true. 
          Ncheck = Ncheck+1 
          iliq(Ncheck) = i
          isol(Ncheck) = 0
          search = dust_nam(i)
          search = search(1:k-1)//'[s]'
          do j=1,NDUST
            if (search==dust_nam(j)) then
              isol(Ncheck) = j
            endif
          enddo
          if (isol(Ncheck)==0) then
            print*,"*** liquid without solid "//trim(dust_nam(i))
            stop
          endif  
        endif
      enddo  
      if (Ncheck==0) return

      !-------------------------------
      ! ***  check melting points  ***
      !-------------------------------
      print*
      print*,'auto-correction for spurious liquid <-> solid '//
     &       'phase transitions ...'
      nat = 1.Q-100
      nmol = 1.Q-100
      do iT=100,10000
        T = DBLE(iT) 
        call SUPERSAT(T,nat,nmol,Sat)
        S(:,iT) = Sat(:)
      enddo  
      do i=1,Ncheck
        il = iliq(i)
        is = isol(i)
        do iT=101,10000
          T = DBLE(iT) 
          old = S(is,iT-1)/S(il,iT-1)
          new = S(is,iT)/S(il,iT)
          if (old>1.Q0.and.new<1.Q0) then
            !print'(A15,"-> ",A15,":",2(0pF8.1))',
     &      !     dust_nam(is),dust_nam(il),T,Tmelt(il)
          else if (old<1.Q0.and.new>1.Q0) then
            !print'(A15,"<- ",A15,":",0pF8.1,
     &      !     " false intersection point")',
     &      !     dust_nam(is),dust_nam(il),T
            if (T<Tmelt(il)) then
              Tcorr(il) = 0.5*(T+Tmelt(il))  
              print'(" ... correct ",A15," T <",0pF7.1)',
     &             dust_nam(il),Tcorr(il) 
            else  
              Tcorr(is) = 0.5*(T+Tmelt(il))  !correct solid
              print'(" ... correct ",A15," T >",0pF7.1)',
     &             dust_nam(is),Tcorr(is) 
            endif  
          endif  
        enddo   
      enddo
      do iT=100,10000
        T = DBLE(iT) 
        call SUPERSAT(T,nat,nmol,Sat)
        S(:,iT) = Sat(:)
      enddo  
      print'(26x,"melting point[K]  should be")'
      do i=1,Ncheck
        il = iliq(i)
        is = isol(i)
        do iT=101,10000
          T = DBLE(iT) 
          old = S(is,iT-1)/S(il,iT-1)
          new = S(is,iT)/S(il,iT)
          if (old>1.Q0.and.new<1.Q0) then
            print'(A15,"-> ",A15,":",2(0pF8.1))',
     &           dust_nam(is),dust_nam(il),T,Tmelt(il)
          else if (old<1.Q0.and.new>1.Q0) then
            print'(A15,"<- ",A15,":",0pF8.1,
     &           " false intersection point")',
     &           dust_nam(is),dust_nam(il),T
            stop
          endif  
        enddo   
      enddo
      stop
      end
