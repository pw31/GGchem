************************************************************************
      subroutine INIT_CHEMISTRY
************************************************************************
      use CHEMISTRY,ONLY: NMOLdim,NMOLE,NELM,catm,cmol,fit,natom,a,
     &    m_kind,m_anz,elnum,elion,charge,
     &    el,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,
     &    Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr
      use EXCHANGE,ONLY: nmol
      implicit none
      integer :: i,ii,j,iel,e
      character(len=2) :: cel(40),elnam
      character(len=200) :: line,elements
      logical :: found,allfound

      open(unit=12, file='dispol_new.dat', status='old')
      write(*,*)
      write(*,*) 'reading molecules and kp-data from dispol_new.dat'
      read(12,'(A200)') elements
      elements = ' '//trim(elements)//' '
      read(12,*) NMOLdim
      allocate(cmol(NMOLdim),fit(NMOLdim),natom(NMOLdim),a(NMOLdim,0:7))
      allocate(m_kind(0:6,NMOLdim),m_anz(6,NMOLdim))
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
        else
          stop "*** unknown element "
        endif   
        print*,'element '//elnam,elnum(NELM)
      enddo
      i=1
      do ii=1,NMOLdim
        read(12,'(A200)',end=200) line
        read(line,'(A10)') cmol(i)
        line = line(11:)
        read(line,*) iel,cel(1:iel),m_anz(1:iel,i)
        read(12,'(A200)') line
        read(line,*) fit(i)
        print*,trim(line),fit(i)
        if (fit(i)==6) then
          read(line,*) fit(i),(a(i,j),j=0,7)
        else   
          read(line,*) fit(i),(a(i,j),j=0,4)
        endif  
        m_kind(0,i) = iel
        natom(i) = 0
        found = .true.
        do j=1,m_kind(0,i)
          natom(i) = natom(i)+m_anz(j,i)
          if (index(elements,cel(j))<=0) found=.false. 
        enddo  
        if (.not.found) cycle
        print'(I4,A10,1x,99(I3,1x,A2,1x))',
     &        i,trim(cmol(i)),(m_anz(j,i),cel(j),j=1,iel)
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
        enddo  
        if (iel==2.and.
     >      ((m_kind(1,i)==el.and.m_anz(1,i)==-1.and.m_anz(2,i)==1).or.
     >       (m_kind(2,i)==el.and.m_anz(2,i)==-1.and.m_anz(1,i)==1))
     >     ) then
          e = m_kind(1,i)
          if (e==el) e=m_kind(2,i)
          elion(e) = i
        endif
        i = i + 1
      enddo
 200  close(12)
      NMOLE = i-1
      allocate(nmol(NMOLdim))
  
      print*,NMOLE,' species'
      print*,NELM,' elements'
      print'(99(A4))',(trim(catm(j)),j=1,NELM)
      print'(99(I4))',elnum(1:NELM)
      !print'(99(I4))',H,He,C,N,O,Si,Mg,Fe,Na,Al,S,Ca,Ti,Cl,K,Li,el
      if (charge) then
        print'(1x,99(A4))',(trim(cmol(elion(j))),j=1,el-1),'  ',
     >                     (trim(cmol(elion(j))),j=el+1,NELM)
      endif  

      end

