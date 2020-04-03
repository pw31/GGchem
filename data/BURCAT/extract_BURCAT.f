      program EXTRACT
      implicit none
      real*8 :: aatm(80,14),amol(3000,14),amol9(3000,5,10),aa(14)
      real*8 :: Trange(3000,5,2),stoich(3000,5),st,sto(4)
      integer :: Natom,Nmol,Nvalid,i,j,k,ifit,Nfit(3000),Nat(3000)
      integer :: nomatch
      logical :: ivalid(3000)
      character(len=200) :: line,frmt,outmol,info,text(10),tmp,oldi
      character(len=18) :: spnam(3000),sp
      character(len=60) :: fullname(3000)
      character(len=2) :: atnam(80),elnam(3000,5),enam,eln(4)
      character(len=1) :: type(3000),ty,oldt

      !-----------------------------------------------------------------
      Natom = 1
      open(1,file='Burcat_atoms.dat',status='old')
      do
        do 
          read(1,'(A200)',end=1000) line
          if (line(80:80)=="1") exit
        enddo
        atnam(Natom) = line(1:2)
        read(1,'(5(E15.8))') aatm(Natom,1:5)
        read(1,'(5(E15.8))') aatm(Natom,6:10)
        read(1,'(5(E15.8))') aatm(Natom,11:14)
        print*,"found atom "//atnam(Natom)
        !print*,aatm(Natom,1:14)
        Natom = Natom+1
      enddo  
 1000 continue
      close(1)
      print*

      !-----------------------------------------------------------------
      Nmol = 1
      stoich = 0.0
      open(1,file='THERM.DAT',status='old')
      do
        do 
          read(1,'(A200)',end=2000) line
          if (line(80:80)=="1") exit
        enddo
        read(line,'(24x,4(A2,F3.0),A1)') 
     >       elnam(Nmol,1),stoich(Nmol,1),
     >       elnam(Nmol,2),stoich(Nmol,2),
     >       elnam(Nmol,3),stoich(Nmol,3),
     >       elnam(Nmol,4),stoich(Nmol,4),
     >       type(Nmol)                         ! stoichiometry
        i=index(line,' ')
        spnam(Nmol) = line(1:i)                 ! name of species
        if (type(Nmol)=="G") spnam(Nmol) = upper(line(1:i))        
        fullname(Nmol) = line(i+1:18)           ! further information
        do i=1,4
          if (elnam(Nmol,i)=="  ") exit
          enam = elnam(Nmol,i)
          elnam(Nmol,i) = enam(1:1)//lower(enam(2:2))
          if (elnam(Nmol,i)=='E ') elnam(Nmol,i)='el'
        enddo
        Nat(Nmol) = i-1                         ! number of elements
        !print'(I5,1x,A12," (",A1,") ",3(A2,1x,F4.1,2x))',
     >  !      Nmol,spnam(Nmol),type(Nmol),
     >  !      (elnam(Nmol,i),stoich(Nmol,i),i=1,Nat(Nmol))
        read(1,'(5(E15.8))') amol(Nmol,1:5)
        read(1,'(5(E15.8))') amol(Nmol,6:10)
        read(1,'(5(E15.8))') amol(Nmol,11:14)
        !print*,amol(Nmol,1:14)
        Nmol = Nmol + 1
      enddo  
 2000 continue
      Nmol = Nmol - 1
      close(1)

      !-----------------------------------------------------------------
      ivalid = .false.
      open(1,file='BURCAT.THR',status='old')
      do j=1,120
        read(1,'(A200)',end=2500) line
      enddo  
      text(:) = " " 
      info = " "
      oldi = " "
      Nvalid = 1
      do 
        do
          read(1,'(A200)',end=2500) line
          do i=10,2,-1
            text(i) = text(i-1)
          enddo  
          text(1) = line
          if (line(80:80)=="1") exit
          i=len(trim(text(2)))
          if (i>2.and.i<20) then 
            tmp = text(2)
            print*,trim(tmp)
            tmp = tmp(index(tmp,"-")+1:)
            !print*,trim(tmp)
            if (index(tmp,"-")>0) info=text(1)
          endif  
        enddo
        sp = line(1:index(line,' '))
        ty = line(45:45)
        if (ty.ne."G".and.ty.ne."L".and.ty.ne."S") then
          print*,"type not recognised."
          print*,trim(line)
          print*,line(43:47)
          stop
        endif  
        if (ty=="G") sp = upper(sp)
        i=index(info," by")
        if (i>0) info=info(1:i-1)
        i=index(info,"from")
        if (i>0) info=info(1:i-1)
        i=index(info,"From")
        if (i>0) info=info(1:i-1)
        i=index(info,"based")
        if (i>0) info=info(1:i-1)
        i=index(info,"STATW")
        if (i>0) info=info(1:i-1)
        i=index(info,"1.")
        if (i>0) info=info(1:i-1)
        i=index(info,"Ia=")
        if (i>0) info=info(1:i-1)
        i=index(info,"IA=")
        if (i>0) info=info(1:i-1)
        i=index(info,"IAIB=")
        if (i>0) info=info(1:i-1)
        i=index(info,"IAIBIC=")
        if (i>0) info=info(1:i-1)
        i=index(info,"Max Lst")
        if (i>0) info=info(1:i-1)
        i=index(info,"according")
        if (i>0) info=info(1:i-1)
        print*,"info = "//trim(info)
        if (info==oldi.and.ty==oldt.and.ty.eq.'G') then
          print*,"something went wrong."
          stop
        endif
        oldt=ty
        oldi=info
        print*,trim(sp)//" "//trim(spnam(Nvalid))//"  "//
     >         ty//" "//type(Nvalid)
        if (trim(sp)==trim(spnam(Nvalid)).and.ty==type(Nvalid)) then
          print*,Nvalid,"... match."
          read(1,'(5(E15.8))') aa(1:5)
          read(1,'(5(E15.8))') aa(6:10)
          read(1,'(5(E15.8))') aa(11:14)
          do i=1,14
            if (aa(i).ne.amol(Nvalid,i)) then
              print*,"data does not agree",i
              print*,aa(1:14)
              print*,amol(Nvalid,1:14)
              stop
            endif
          enddo  
          read(line,'(24x,4(A2,F3.0),A1)') 
     >       eln(1),sto(1),
     >       eln(2),sto(2),
     >       eln(3),sto(3),
     >       eln(4),sto(4)
          print*,eln(1:4)
          print*,elnam(Nvalid,1:4)
          do i=1,4
            enam = eln(i)
            eln(i) = enam(1:1)//lower(enam(2:2))
            if (eln(i)=='E ') eln(i)='el'
            if (eln(i).ne.elnam(Nvalid,i).or.
     >          sto(i).ne.stoich(Nvalid,i)) then
              print*,"stoichiometry does not match."
              print*,eln
              print*,sto
              stop
            endif  
          enddo
          fullname(Nvalid) = info
          ivalid(Nvalid) = .true.
          Nvalid = Nvalid+1
          nomatch = 0
        else
          print*,"... no match"
          nomatch = nomatch+1
          if (nomatch>0) then
            print*,"*** too many mismatches"
            stop
            goto 2500
          endif  
        endif  
      enddo
 2500 continue
      close(1)

      if (.false.) then
      !-----------------------------------------------------------------
      Nmol = 1
      open(1,file='NEWNASA.TXT',status='old')
      do i=1,39
        read(1,*) line
      enddo  
      do
        read(1,'(A200)',end=3000) line
        if (line(1:1)==' ') read(1,'(A200)') line
        if (line(1:1)==' ') read(1,'(A200)') line
        !print*,trim(line)
        i=index(line,' ')
        spnam(Nmol) = line(1:i)                 ! name of species
        read(1,'(I2,8x,5(A2,F6.2))') Nfit,    
     >       elnam(Nmol,1),stoich(Nmol,1),
     >       elnam(Nmol,2),stoich(Nmol,2),
     >       elnam(Nmol,3),stoich(Nmol,3),
     >       elnam(Nmol,4),stoich(Nmol,4),
     >       elnam(Nmol,5),stoich(NMol,5)       ! stoichiometry
        do i=1,5
          if (elnam(Nmol,i)=="  ") exit
          enam = elnam(Nmol,i)
          elnam(Nmol,i) = enam(1:1)//lower(enam(2:2))
        enddo
        Nat(Nmol) = i-1                         ! number of elements
        do ifit=1,Nfit(Nmol)
          read(1,*) Trange(Nmol,ifit,1:2)
          !print*,Trange(ifit,1:2)
          read(1,'(5(D16.9))') amol9(Nmol,ifit,1:5)
          read(1,'(5(D16.9))') amol9(Nmol,ifit,6:10)
        enddo  
        if (Nat(Nmol)==1.and.stoich(Nmol,1)==1) then
          print*,trim(line)
          print*,spnam(Nmol),Nat,Nfit
        endif  
        Nmol = Nmol + 1
      enddo 
 3000 continue
      close(1)
      endif

      !----------------------------
      ! ***  create output file ***
      !----------------------------
      Nvalid = 0
      outmol = ' CO H2 N2 H- H+ SH ALO ALO2- CNO- CH '
      do i=1,Nmol
        if (type(i).ne.'G') then
          ivalid(i)=.false.
        else if (Nat(i)==1.and.stoich(i,1)==1) then  
          ivalid(i)=.false.
        !else if (index(outmol," "//trim(spnam(i))//" ")<=0) cycle
        else if (ivalid(i)) then
          Nvalid = Nvalid+1
        endif  
      enddo  
      open(1,file='../dispol_BURCAT.dat')
      write(1,*) Nvalid
      Nvalid = 0
      do i=1,Nmol
        if (.not.ivalid(i)) cycle
        Nvalid = Nvalid+1
        sp = checkname(spnam(i))
        sp = check_iso(sp,i,spnam,stoich,elnam)
        spnam(i) = sp
        print'(I5,1x,A18," (",a1,") ",4(A2,1x,F4.1,2x))',
     >        Nvalid,sp,type(i),
     >        (elnam(i,j),stoich(i,j),j=1,Nat(i))
        aa(1:14) = amol(i,1:14)
        !print*,aa
        do j=1,Nat(i)
          enam = elnam(i,j)
          st   = stoich(i,j)
          do k=1,Natom
            if (enam==atnam(k)) exit
          enddo
          if (k>Natom) stop "element not found."
          !print*,atnam(k),aatm(k,1:14)
          aa(1:14) = aa(1:14) - st*aatm(k,1:14)
        enddo
        write(frmt,'("(A18,I2,",I1,"(A3),",I1,"(I3),",I2,"x,A60)")') 
     >        Nat(i),Nat(i),6*(4-Nat(i))+1
        write(1,frmt) sp,Nat(i),
     >                elnam(i,1:Nat(i)),int(stoich(i,1:Nat(i))),
     >                fullname(i)
        write(1,'(I2,14(1pE16.8))') 8,aa(1:14)
      enddo
      close(1)

      contains

!-------------------------------------------------------------------------
      function lower(str) result(string)
!-------------------------------------------------------------------------
      implicit none
      character(*),intent(In) :: str
      character(len(str)) :: string
      integer :: ic,i
      character(26),parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26),parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
      string = str
      do i = 1,LEN_TRIM(str)
        ic = index(cap, str(i:i))
        if (ic>0) string(i:i) = low(ic:ic)
      enddo
      end function lower

!-------------------------------------------------------------------------
      function upper(str) result(string)
!-------------------------------------------------------------------------
      implicit none
      character(*),intent(In) :: str
      character(len(str)) :: string
      integer :: ic,i
      character(26),parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26),parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
      string = str
      do i = 1,LEN_TRIM(str)
        ic = index(low, str(i:i))
        if (ic>0) string(i:i) = cap(ic:ic)
      enddo
      end function upper

!-------------------------------------------------------------------------
      function checkname(str) result(string)
!-------------------------------------------------------------------------
      implicit none
      character(*),intent(In) :: str
      character(len(str)) :: string
      integer :: i,ic,N
      character(3),parameter :: nono = '*=-'
      string = ''
      N = LEN_TRIM(str)
      !print*,"###"//trim(str)//"###  ",N
      do i=1,N
        ic = index(nono, str(i:i))
        if (i==N.and.str(i:i)=='-') ic=0
        if (i==N-1.and.str(i:i+1)=='--') ic=0
        if (i==N-2.and.str(i:i+2)=='---') ic=0
        if (ic==0) string=trim(string)//str(i:i)
      enddo
      ic = index(string,'(O)')
      if (ic>1) then
        string=string(1:ic-1)//"O"//string(ic+3:)
      endif  
      if (str.ne.string) then
        print*,"*** changed species name from "//trim(str)//
     >         " -> "//trim(string)
      endif  
      end function checkname

!-------------------------------------------------------------------------
      function check_iso(str,N,spnam,stoich,elnam) result(string)
!-------------------------------------------------------------------------
      implicit none
      character(*),intent(IN) :: str
      integer,intent(IN) :: N
      character(len=18),intent(IN) :: spnam(3000)
      real*8,intent(IN) :: stoich(3000,5)
      character(len(str)) :: string
      character(len=2) :: elnam(3000,5)
      character(len=5) :: add
      character(len=1) :: char
      logical :: double,is_double
      integer :: i,j
      string = str
      is_double=.false.
      add = ':1'
      do i=1,N-1
        if (str==spnam(i)) then
          double=.true.
          do j=1,4
            if ( elnam(i,j).ne. elnam(N,j)) double=.false.
            if (stoich(i,j).ne.stoich(N,j)) double=.false.
          enddo
          if (double) then
            !print*,i,N
            !print*,trim(str),trim(spnam(i))
            is_double=.true.
          endif  
        endif  
        if (trim(str)//':1'==spnam(i)) add=':2'
        if (trim(str)//':2'==spnam(i)) add=':3'
        if (trim(str)//':3'==spnam(i)) add=':4'
        if (trim(str)//':4'==spnam(i)) add=':5'
        if (trim(str)//':5'==spnam(i)) add=':6'
        if (trim(str)//':5'==spnam(i)) add=':6'
        if (trim(str)//':6'==spnam(i)) add=':7'
        if (trim(str)//':7'==spnam(i)) add=':8'
        if (trim(str)//':8'==spnam(i)) add=':9'
        if (trim(str)//':9'==spnam(i)) add=':a'
        if (trim(str)//':a'==spnam(i)) add=':b'
        if (trim(str)//':b'==spnam(i)) add=':c'
        if (trim(str)//':c'==spnam(i)) stop
      enddo
      if (is_double) then
        string=trim(str)//add
        print*,"*** changed species name from "//trim(str)//
     >         " -> "//trim(string)
        !read(*,'(A1)') char
      endif  
      end function check_iso

      end
