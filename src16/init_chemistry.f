************************************************************************
      subroutine INIT_CHEMISTRY
************************************************************************
      use PARAMETERS,ONLY: elements,initchem_info
      use CHEMISTRY,ONLY: NMOLdim,NMOLE,NELM,catm,cmol,el,
     &    dispol_file,source,fit,natom,a,error,i_nasa,
     &    m_kind,m_anz,elnum,elion,charge,
     &    el,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,
     &    Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,W
      use DUST_DATA,ONLY: mass,mel,amu
      use EXCHANGE,ONLY: nmol,mmol
      implicit none
      integer :: loop,i,ii,j,iel,e,smax,ret
      character(len=2) :: cel(40),elnam
      character(len=20) :: molname,upper,leer='                    '
      character(len=200) :: filename
      character(len=300) :: line
      logical :: found,charged
      real*8 :: fiterr

      cel(:) = '.'
      read(elements,*,end=100) cel
 100  NELM = 0
      charge = .false.
      do i=1,99
        if (cel(i)=='.') exit
        elnam = cel(i)
        found = .false.
        do e=1,NELM
          if (elnam.eq.catm(e)) then
            found=.true.
            exit
          endif  
        enddo
        NELM = NELM+1
        catm(NELM) = elnam
        if     (elnam=='el') then; el=NELM ; charge=.true.
        elseif (elnam=='H')  then;  H=NELM ; elnum(NELM)=1 
        elseif (elnam=='He') then; He=NELM ; elnum(NELM)=2 
        elseif (elnam=='Li') then; Li=NELM ; elnum(NELM)=3
        elseif (elnam=='Be') then; Be=NELM ; elnum(NELM)=4 
        elseif (elnam=='B')  then;  B=NELM ; elnum(NELM)=5 
        elseif (elnam=='C')  then;  C=NELM ; elnum(NELM)=6 
        elseif (elnam=='N')  then;  N=NELM ; elnum(NELM)=7 
        elseif (elnam=='O')  then;  O=NELM ; elnum(NELM)=8 
        elseif (elnam=='F')  then;  F=NELM ; elnum(NELM)=9 
        elseif (elnam=='Ne') then; Ne=NELM ; elnum(NELM)=10
        elseif (elnam=='Na') then; Na=NELM ; elnum(NELM)=11 
        elseif (elnam=='Mg') then; Mg=NELM ; elnum(NELM)=12
        elseif (elnam=='Al') then; Al=NELM ; elnum(NELM)=13
        elseif (elnam=='Si') then; Si=NELM ; elnum(NELM)=14
        elseif (elnam=='P')  then;  P=NELM ; elnum(NELM)=15 
        elseif (elnam=='S')  then;  S=NELM ; elnum(NELM)=16
        elseif (elnam=='Cl') then; Cl=NELM ; elnum(NELM)=17
        elseif (elnam=='Ar') then; Ar=NELM ; elnum(NELM)=18
        elseif (elnam=='K')  then;  K=NELM ; elnum(NELM)=19
        elseif (elnam=='Ca') then; Ca=NELM ; elnum(NELM)=20
        elseif (elnam=='Sc') then; Sc=NELM ; elnum(NELM)=21
        elseif (elnam=='Ti') then; Ti=NELM ; elnum(NELM)=22
        elseif (elnam=='V')  then;  V=NELM ; elnum(NELM)=23
        elseif (elnam=='Cr') then; Cr=NELM ; elnum(NELM)=24
        elseif (elnam=='Mn') then; Mn=NELM ; elnum(NELM)=25
        elseif (elnam=='Fe') then; Fe=NELM ; elnum(NELM)=26
        elseif (elnam=='Co') then; Co=NELM ; elnum(NELM)=27
        elseif (elnam=='Ni') then; Ni=NELM ; elnum(NELM)=28
        elseif (elnam=='Cu') then; Cu=NELM ; elnum(NELM)=29
        elseif (elnam=='Zn') then; Zn=NELM ; elnum(NELM)=30
        elseif (elnam=='Ga') then; Ga=NELM ; elnum(NELM)=31
        elseif (elnam=='Ge') then; Ge=NELM ; elnum(NELM)=32 
        elseif (elnam=='As') then; As=NELM ; elnum(NELM)=33 
        elseif (elnam=='Se') then; Se=NELM ; elnum(NELM)=34 
        elseif (elnam=='Br') then; Br=NELM ; elnum(NELM)=35 
        elseif (elnam=='Kr') then; Kr=NELM ; elnum(NELM)=36 
        elseif (elnam=='Rb') then; Rb=NELM ; elnum(NELM)=37 
        elseif (elnam=='Sr') then; Sr=NELM ; elnum(NELM)=38 
        elseif (elnam=='Y')  then;  Y=NELM ; elnum(NELM)=39 
        elseif (elnam=='Zr') then; Zr=NELM ; elnum(NELM)=40
        elseif (elnam=='W')  then;  W=NELM ; elnum(NELM)=41
        else
          stop "*** unknown element "
        endif
        if (initchem_info) then
          print*,'element '//elnam,elnum(NELM)
        endif
      enddo

      NMOLdim = 10000
      allocate(cmol(NMOLdim),fit(NMOLdim),natom(NMOLdim))
      allocate(error(NMOLdim),a(NMOLdim,0:13))
      allocate(source(NMOLdim),m_kind(0:6,NMOLdim),m_anz(6,NMOLdim))
      i=1
      i_nasa = 0
      do loop=1,4
        filename = trim(dispol_file(loop))
        if (filename=='') exit
        filename = 'data/'//trim(filename)
        if (initchem_info) write(*,*)
        write(*,*) 'reading kp-data from '
     &             //trim(filename)//" ..."
        open(unit=12, file=filename, status='old')
        read(12,*) NMOLdim
        do ii=1,NMOLdim
          read(12,'(A300)') line
          read(line,*) molname,iel,cel(1:iel),m_anz(1:iel,i)
          molname=trim(molname)
          fiterr = 0.0
          j = index(line,"+/-")
          if (j>0) read(line(j+3:),*) fiterr
          error(i) = fiterr
          read(12,'(A300)') line
          read(line,*) fit(i)
          if (fit(i)==6) then
            read(line,*) fit(i),(a(i,j),j=0,7)
          elseif(fit(i)==7) then
             i_nasa = 1
             read(line,*) fit(i),(a(i,j),j=0,13)
          else   
            read(line,*) fit(i),(a(i,j),j=0,4)
          endif  
          m_kind(0,i) = iel
          natom(i) = 0
          found = .true.
          smax  = 0
          do j=1,m_kind(0,i)
            natom(i) = natom(i)+m_anz(j,i)
            if (index(elements,cel(j))<=0) found=.false. 
            smax = MAX(smax,ABS(m_anz(j,i)))
          enddo  
          if (.not.found) cycle    ! molecule has unselected element 
          if (smax>16) cycle       ! stoichiometric coefficient > 16
          if (m_kind(0,i)==1.and.natom(i)==1) cycle  ! pure atom
          j = index(molname,"_")
          if (j>1) then
            cmol(i) = upper(molname(j+1:)//leer(1:j))
          else
            cmol(i) = upper(molname)
          endif
          charged = .false.
          do j=1,m_kind(0,i)
            elnam = cel(j)
            found = .false.
            do e=1,NELM
              if (elnam.eq.catm(e)) then
                found=.true.
                exit
              endif  
            enddo
            if (.not.found) stop "*** should not occur"
            m_kind(j,i) = e
            if (e==el) charged=.true. 
          enddo  
          if (fit(i)==6.and.charged) cycle ! old charged BarklemCollet
          source(i) = loop
          call CHECK_DOUBLE(cmol(i),m_kind(:,i),m_anz(:,i),i,loop,ret)
          if (ret>0) then
            source(ret) = loop
            cmol(ret) = cmol(i)
            fit(ret)  = fit(i)
            a(ret,:)  = a(i,:)
            error(ret)= error(i)
            if (initchem_info) then
              write(line,'(I4,A20,1x,99(I3,1x,A2,1x))')
     &           ret,trim(cmol(ret)),(m_anz(j,ret),cel(j),j=1,iel)
              print*,trim(line)//"    OVERWRITE"
            endif  
          else  
            if (initchem_info) then
              write(line,'(I4,A20,1x,99(I3,1x,A2,1x))')
     &            i,trim(cmol(i)),(m_anz(j,i),catm(m_kind(j,i)),j=1,iel)
              if (loop==1) then
                print*,trim(line)
              else
                print*,trim(line)//"    --> NEW" 
              endif
            endif  
            if (iel==2.and.
     >       ((m_kind(1,i)==el.and.m_anz(1,i)==-1.and.m_anz(2,i)==1).or.
     >        (m_kind(2,i)==el.and.m_anz(2,i)==-1.and.m_anz(1,i)==1))
     >        ) then
              e = m_kind(1,i)
              if (e==el) e=m_kind(2,i)
              elion(e) = i
            endif
            i = i+1
          endif  
        enddo
 200    close(12)
      enddo  
      NMOLE = i-1
      allocate(nmol(NMOLE),mmol(NMOLE))

      if(i_nasa==1) call NASA_POLYNOMIAL !Added by Yui Kawashima

      if (loop>1.and.initchem_info) then
        print* 
        do i=1,NMOLE
          print*,i,cmol(i),' ->  '//trim(dispol_file(source(i)))
        enddo
      endif  
  
      !open(unit=1,file='chemicals.tex')
      !write(1,*) NMOLE
      !do i=1,NMOLE
      !  if (error(i)>0.0) then 
      !    write(1,3000)
     &!      i,cmol(i),source(i),fit(i),a(i,0:4),error(i)
      !  else  
      !    write(1,3010)
     &!      i,cmol(i),source(i),fit(i),a(i,0:4)
      !  endif  
      !enddo  
      !close(1)
      !stop

      do i=1,NMOLE
        mmol(i) = 0.d0
        do j=1,m_kind(0,i)
          if (m_kind(j,i)==el) then
            mmol(i) = mmol(i) + m_anz(j,i)*mel
          else
            mmol(i) = mmol(i) + m_anz(j,i)*mass(elnum(m_kind(j,i)))
          endif
        enddo
        !print*,cmol(i),mmol(i)/amu
      enddo  

      print* 
      print*,NMOLE,' species'
      print*,NELM,' elements'
      print'(99(A4))',(trim(catm(j)),j=1,NELM)
      print'(99(I4))',elnum(1:NELM)
      !print'(99(I4))',H,He,C,N,O,Si,Mg,Fe,Na,Al,S,Ca,Ti,Cl,K,Li,el
      if (charge) then
        print'(1x,99(A4))',(trim(cmol(elion(j))),j=1,el-1),'  ',
     >                     (trim(cmol(elion(j))),j=el+1,NELM)
      endif  
      
 3000 format(I4," & ",A12," & (",I1,") & ",I1," & ",
     &       5(1pE12.5," & "),"$\pm$",0pF4.2,"\\")
 3010 format(I4," & ",A12," & (",I1,") & ",I1," & ",
     &       5(1pE12.5," & "),"\\")
      end

************************************************************************
      subroutine CHECK_DOUBLE(molname,kind,anz,N,loop,ret)
************************************************************************
      use PARAMETERS,ONLY: initchem_info
      use CHEMISTRY,ONLY: cmol,m_kind,m_anz,dispol_file,source
      implicit none
      character(len=20) :: molname
      integer,intent(IN) :: kind(0:6),anz(6),N,loop
      integer,intent(OUT) :: ret
      integer :: i,j,jj,el,ambi
      logical :: found,allfound,eqname,eqsource

      ret  = 0
      ambi = 0
      do i=1,N-1
        if (kind(0).ne.m_kind(0,i)) cycle   ! different no(elements)
        allfound=.true.
        do j=1,kind(0)
          el = kind(j)
          found = .false.
          do jj=1,m_kind(0,i)
            if (el.ne.m_kind(jj,i)) cycle
            found = .true.
            exit
          enddo
          if (.not.found) then
            allfound = .false.
            exit                            ! different elements
          else if (anz(j).ne.m_anz(jj,i)) then
            allfound = .false.
            exit                            ! different stoich.fac.
          endif
        enddo
        if (.not.allfound) cycle
        eqname = (trim(molname)==trim(cmol(i)))
        eqsource = (loop==source(i))
        if (eqname.and.eqsource) then
          print*,"*** double molecule in "//dispol_file(loop)
          print*,trim(molname)//", "//trim(cmol(i))
          stop
        else if ((.not.eqname).and.eqsource.and.loop==1) then
          if (initchem_info) then
            print*,trim(molname)//", "//trim(cmol(i))//
     &           " different isomere in first source is OK"
          endif
          return  
        else if (eqname.and.(.not.eqsource)) then  
          ret = i
          return
        else
          ambi = i 
        endif
      enddo
      if (ambi>0) then
        if (source(ambi)==loop) then 
          if (initchem_info) then
            print*,trim(molname)//", "//trim(cmol(ambi))//
     &           " different isomere in subsequent source is OK"
          endif  
          ret = 0
          return
        else  
          print*,"*** "//trim(molname)//", "//trim(cmol(ambi))//
     &         " ambiguous names in ..."
          print*,trim(dispol_file(loop))//
     &         ", "//trim(dispol_file(source(ambi)))
          print*,"please equalise in both data files."
          stop 
        endif  
      endif  
      end
