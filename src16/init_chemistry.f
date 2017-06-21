************************************************************************
      subroutine INIT_CHEMISTRY
************************************************************************
      use CHEMISTRY,ONLY: NMOLE,NELM,catm,cmol,fit,natom,a,
     &    m_kind,m_anz,elnum,elion,charge,
     &    el,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,
     &    Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr
      use EXCHANGE,ONLY: nmol
      implicit none
      integer :: i,j,iel,e
      character(len=2) :: cel(8),elnam
      character(len=100) :: line
      logical :: found

      open(unit=12, file='dispol_new.dat', status='old')
      write(*,*)
      write(*,*) 'reading molecules and kp-data from dispol_new.dat'
      read(12,*) NMOLE
      allocate(nmol(NMOLE))
      allocate(cmol(NMOLE),fit(NMOLE),natom(NMOLE),a(NMOLE,0:4))
      allocate(m_kind(0:8,NMOLE),m_anz(8,NMOLE))
      charge=.false.
      NELM = 0
      !catm(1:3)  = (/'He','el','H'/)
      !elnum(1:3) = (/  2 ,  0 , 1 /) 
      !He = 1
      !el = 2
      !H  = 3
      do i=1,NMOLE
        read(12,'(A100)') line
        read(line,'(A10)') cmol(i)
        line = line(11:)
        read(line,*) iel,cel(1:iel),m_anz(1:iel,i)
        read(12,*) fit(i),(a(i,j),j=0,4)
        m_kind(0,i) = iel
        natom(i) = 0
        do j=1,m_kind(0,i)
          natom(i) = natom(i)+m_anz(j,i)
        enddo  
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
          m_kind(j,i) = e
          if (.not.found) then
            print*,'new element '//elnam
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
          endif  
        enddo
        if (iel==2.and.
     >      ((m_kind(1,i)==el.and.m_anz(1,i)==-1.and.m_anz(2,i)==1).or.
     >       (m_kind(2,i)==el.and.m_anz(2,i)==-1.and.m_anz(1,i)==1))
     >     ) then
          e = m_kind(1,i)
          if (e==el) e=m_kind(2,i)
          elion(e) = i
        endif    
      enddo
      close(12)
  
      print*,NMOLE,' species'
      print*,NELM,' elements'
      print'(99(A3))',(trim(catm(j)),j=1,NELM)
      print'(99(I3))',elnum(1:NELM)
      !print'(99(I3))',H,C,O,N,el,Si,Al,Mg,S,K,Na,Fe,Ti,Cl,Ca,Li,He
      if (charge) then
        print'(1x,99(A3))',(trim(cmol(elion(j))),j=1,el-1),'  ',
     >                     (trim(cmol(elion(j))),j=el+1,NELM)
      endif  
      end

