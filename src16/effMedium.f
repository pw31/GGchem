!----------------------------------------------------------------------
      subroutine effMedium(ilam,Vs,neff,keff)
!----------------------------------------------------------------------
      use DUST_DATA,ONLY: dust_nam
      use OPACITY,ONLY: NLIST,opind,nn,kk
      implicit none
      integer,intent(IN) :: ilam
      real,intent(IN) :: Vs(NLIST)
      real,intent(OUT) :: neff,keff
      complex(kind=8) :: mm(NLIST),ee(NLIST),mmean,emean
      complex(kind=8) :: m2e,e2m
      real :: Fold(2),FF(2),DF(2,2),corr(2),xx(2),xnew(2)
      real :: FF1(2),FF2(2),FF3(2),FF4(2),de1,de2,qual
      integer,parameter :: itmax=30
      integer :: i,it,start
      logical :: unphysical

      mmean = DCMPLX(0.d0,0.d0)
      emean = DCMPLX(0.d0,0.d0)
      do i=1,NLIST
        mm(i) = DCMPLX( nn(ilam,i), kk(ilam,i) )
        ee(i) = m2e(mm(i))
        mmean = mmean + Vs(i)*mm(i)
        emean = emean + Vs(i)*ee(i)
        !write(*,*) i,dust_nam(opind(i)),Vs(i)
        !write(*,*) "   (n,k)=",mm(i)
        !write(*,*) " (er,ei)=",ee(i)
        !write(*,*) "   Probe=",e2m(ee(i))
      enddo
      !write(*,*) "   <n,k>=",mmean
      !write(*,*) "         ",e2m(emean)
      !write(*,*) " <ei,er>=",emean
      !write(*,*) "         ",m2e(mmean)
      !write(*,*)
      
      !----- use Bruggeman method (see Bohren & Huffman) -----
      do start=1,2+NLIST
        if (start.eq.1) then
          mmean=mmean
        elseif (start.eq.2) then
          mmean=e2m(emean)
        else
          mmean=mm(start-2)
        endif
        !write(*,*) "trying initial value no.",start
        do it=1,itmax
          xx(1) = REAL(mmean,kind=8)
          xx(2) = DIMAG(mmean)
          call EFF_FUNC(mmean,FF,NLIST,Vs,mm) ! this complex function=!0
          qual = FF(1)**2+FF(2)**2
          de1 = xx(1)*1.d-5
          de2 = xx(2)*1.d-5
          call EFF_FUNC(mmean+DCMPLX(de1,0.d0),FF1,NLIST,Vs,mm)
          call EFF_FUNC(mmean-DCMPLX(de1,0.d0),FF2,NLIST,Vs,mm)
          call EFF_FUNC(mmean+DCMPLX(0.d0,de2),FF3,NLIST,Vs,mm)
          call EFF_FUNC(mmean-DCMPLX(0.d0,de2),FF4,NLIST,Vs,mm)
          DF(1,1) = (FF1(1)-FF2(1)) / (2.d0*de1)
          DF(1,2) = (FF3(1)-FF4(1)) / (2.d0*de2)
          DF(2,1) = (FF1(2)-FF2(2)) / (2.d0*de1)
          DF(2,2) = (FF3(2)-FF4(2)) / (2.d0*de2)
          Fold = FF
          call GAUSS8(2,2,DF,corr,FF)
          corr = -corr
          call EFF_PULLBACK(2,xx,corr,Fold,xnew,unphysical,NLIST,Vs,mm)
          !write(*,'(i3," mm=(",1pE11.4,",",1pE11.4,") qual=",1pE10.4)') 
     &    !   it,xx(1),xx(2),qual
          if (unphysical) exit
          mmean = DCMPLX(xnew(1),xnew(2))
          if (ABS(qual).lt.1.d-13) exit
        enddo           
        if (unphysical.or.(it.ge.itmax)) then
          unphysical=.true.
          if (start<2+NLIST) then
            !write(*,*) "no convergence in effMedium (",start,")"//
     &      !           " - trying next init.val."
          else
            write(*,*) "no convergence in effMedium."
          endif  
        else
          exit
        endif
      enddo  
         
      if (unphysical) then
        write(*,*) "  effMedium could not find the physical solution"
        stop
      else  
        !write(*,*) "  effMedium converged after ",it," iterations"
      endif  
         
      !write(*,*) " (er,ei)eff=",emean
      !write(*,*) "   (n,k)eff=",mmean

      neff = REAL(mmean,kind=8)
      keff = DIMAG(mmean)
      emean = m2e(mmean)   

      end

        
!========================================================================
      subroutine EFF_PULLBACK(N,xx,dx,Fold,xnew,unphys,NLIST,Vs,mm)
      implicit none
      integer,intent(in) :: N,NLIST
      real,intent(in) :: xx(N),Fold(N),Vs(NLIST)
      complex(kind=8),intent(in) :: mm(NLIST)
      real,intent(inout) :: dx(N)
      real,intent(out) :: xnew(N)
      logical,intent(out) :: unphys
      real(kind=8) :: Fnew(2),fac,qold,qnew
      complex(kind=8) :: mwork
      integer,parameter :: itmax=20 
      integer :: it

      qold=Fold(1)**2+Fold(2)**2
      fac=1.d0
      do it=1,itmax
        xnew = xx + fac*dx
        if ((xnew(1).gt.0.d0).and.(xnew(2).gt.0.d0)) then 
          mwork = DCMPLX(xnew(1),xnew(2))
          call EFF_FUNC(mwork,Fnew,NLIST,Vs,mm)
          qnew = Fnew(1)**2+Fnew(2)**2
          !write(*,*) it,qold,qnew
          unphys = .false.
          if (qnew<qold) exit
        else
          !write(*,*) it,"negative (n,k)",xnew
          unphys=.true.
        endif
        fac=fac*0.7
      enddo
      return
      end      

!========================================================================
      function e2m(e)
      implicit none
      complex(kind=8) :: e2m
      complex(kind=8), intent(in) :: e
      real(kind=8) :: ereal, eimag, n, k
      real(kind=8) :: sqrte2
      ereal = real(e,kind=8)
      eimag = dimag(e)
      sqrte2 = dsqrt(ereal**2 + eimag**2)
      n = dsqrt( 0.5d0 *( ereal + sqrte2))
      k = dsqrt( 0.5d0 *(-ereal + sqrte2))
      e2m = dcmplx(n, k)      
      return
      end 
!========================================================================
      function m2e(m)
      implicit none
      complex(kind=8) :: m2e
      complex(kind=8), intent(in) :: m
      real :: ereal, eimag, n, k
      n = real(m,kind=8)
      k = dimag(m)
      ereal = n**2 - k**2
      eimag = 2.d0 * n * k
      m2e = dcmplx(ereal, eimag)       
      return
      end 
!========================================================================
      subroutine EFF_FUNC(mmean,FF,NLIST,Vs,mm)
!     -----------------  Formel nach Bruggeman (1935)  ------------------
      implicit none
      complex(kind=8),intent(in) :: mmean
      real,intent(out) :: FF(2) 
      integer,intent(in) :: NLIST
      real,intent(in) :: Vs(NLIST)
      complex(kind=8),intent(in) :: mm(NLIST)
      complex(kind=8) :: Fcplx,mm2,mmi2
      integer :: i,kind
      real,parameter :: gamma=2.d0
      mm2 = mmean**2
      Fcplx = dcmplx(0.d0,0.d0)
      do i=1,NLIST
        mmi2 = mm(i)**2
        Fcplx = Fcplx + Vs(i)*(mmi2-mm2)/(mmi2+gamma*mm2)
      enddo
      FF(1) = REAL(Fcplx,kind=8)
      FF(2) = DIMAG(Fcplx)
      return
      end

