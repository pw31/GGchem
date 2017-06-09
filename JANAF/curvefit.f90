!*******************************************************
!*    Program to demonstrate the Parafit subroutine    *
!* --------------------------------------------------- *
!*  Reference: BASIC Scientific Subroutines, Vol. II   *
!*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
!*                                                     *
!*                  F90 version by J-P Moreau, Paris   *
!*                          (www.jpmoreau.fr)          *
!* --------------------------------------------------- *
!* SAMPLE RUN:                                         *
!*                                                     *
!* Parametric least squares fit                        *
!*                                                     *
!* The input data are:                                 *
!*                                                     *
!* X( 1) =  1.000000   Y( 1) =  0.033702               *
!* X( 2) =  2.000000   Y( 2) =  0.249029               *
!* X( 3) =  3.000000   Y( 3) =  0.944733               *
!* X( 4) =  4.000000   Y( 4) =  1.840089               *
!* X( 5) =  5.000000   Y( 5) =  1.840089               *
!* X( 6) =  6.000000   Y( 6) =  0.944733               *
!* X( 7) =  7.000000   Y( 7) =  0.249029               *
!* X( 8) =  8.000000   Y( 8) =  0.033702               *
!* X( 9) =  9.000000   Y( 9) =  0.002342               *
!* X(10) = 10.000000   Y(10) =  0.000084               *
!*                                                     *
!* The coefficients are:                               *
!*  2.000000000000000                                  *
!*  4.500000000000000                                  *
!*  3.000000000000000                                  *
!*                                                     *
!* The standard deviation of the fit is 0.000000       *
!*                                                     *
!* The number of iterations was 49                     *
!*******************************************************
PROGRAM Parafit

  parameter(LMAX=10,SIZE=25)

  real*8 :: A(LMAX),X(SIZE),Y(SIZE)

  real*8 :: dev,e,ee1
  integer :: Nit,Ndat,Npara,i

  Ndat = 10; Npara = 3; Nit = 0; dev = 0.d0

  print *, 'PARAMETRIC LEAST SQUARES FIT'
  print *, ' '
  print *, 'The input data are:'
  print *, ' '
  do i = 1,Ndat
    X(i) = dfloat(i)
    Y(i) = 2.d0 * dexp(-(X(i) - 4.5d0) * (X(i) - 4.5d0) / 3.d0);
    write(*,50) i,X(i),i,Y(i)
  end do
  
  e = 0.1d0; ee1 = 0.5d0; A(1) = 10.d0; A(2) = 10.d0; A(3) = 10.d0

  call Param_LS(Npara,Nit,Ndat,dev,e,ee1,A,X,Y)

  print *, ' '
  print *,'The coefficients are:'
  print *, A(1:Npara) 
  print *, ' '
  write (*,60)  dev
  write (*,70)  Nit
  print *, ' '

  stop

50 format(' X(',I2,') = ',F9.6,'   Y(',I2,') = ',F9.6)
60 format(' The standard deviation of the fit is ',F10.6)
70 format(' The number of iterations was ',I2)
 
End !main program


  !Function subroutine
  Subroutine S500(l,x,y,A)  
  real*8 x,y,A(l)
    y = A(1) * dexp(-(x - A(2)) * (x - A(2)) / A(3))
    return
  end

  !Residual generation subroutine
  Subroutine S200(l,l2,n,d,A,X,Y) 
  parameter(SIZE=25)
  real*8 d,l2,A(l),X(SIZE),Y(SIZE) 
  integer l,n,j
  real*8 xx,yy
    l2 = 0.d0
    do j = 1, n
      xx = X(j)
      ! Obtain function
      call S500(l,xx,yy,A)
      l2 = l2 + (Y(j) - yy) * (Y(j) - yy)
    end do
    d = dsqrt(l2 / dfloat(n - l))
    return
  end


!***************************************************************
!* Parametric least squares curve fit subroutine. This program *
!* least squares fits a function to a set of data values by    *
!* successively reducing the variance. Convergence depends on  *
!* the initial values and is not assured.                      *
!* n pairs of data values, X(i), Y(i), are given. There are l  *
!* parameters, A(j), to be optimized across.                   *
!* Required are initial values for the A(l) and e. Another     *
!* important parameter which affects stability is e1, which is *
!* initially converted to e1(l)for the first intervals.        *
!* The parameters are multiplied by (1 - e1(i)) on each pass.  *
!***************************************************************
Subroutine Param_LS(l,m,n,d,e,ee1,A,X,Y)  
  !Labels: 50,100
  parameter(SIZE=25)
  real*8 d,e,ee1,A(l),X(SIZE),Y(SIZE)
  real*8 E1(l)
  integer i,l,m,n
  real*8 a0,l1,l2,m0,m1
  do i = 1, l
    E1(i) = ee1
  end do	;
  !Set up test residual
  l1 = 1.d6
  !Make sweep through all parameters
50 do i = 1, l
    a0 = A(i)
    !Get value of residual
    A(i) = a0
100 call S200(l,l2,n,d,A,X,Y)
    !Store result in m0
    m0 = l2
    !Repeat for m1
    A(i) = a0 * (1.d0 - E1(i))
    call S200(l,l2,n,d,A,X,Y)
    m1 = l2
    !Change interval size if called for
    !If variance was increased, halve E1(i) 
    if (m1 > m0)  then
       E1(i) = -E1(i) / 2.d0
    end if
    !If variance was reduced, increase step size by increasing E1(i)
    if (m1 < m0)  then
       E1(i) = 1.2d0 * E1(i)
    end if 
    !If variance was increased, try to reduce it
    if (m1 > m0) then  
       A(i) = a0
    end if 
    if (m1 > m0)  then
       goto 100
    end if
  end do !i loop
  !End of a complete pass
  !Test for convergence 
  m = m + 1
  if (l2.eq.0.d0) then
    return
  end if 
  if (dabs((l1 - l2) / l2) < e)  then
    return
  end if 
  !If this point is reached, another pass is called for
  l1 = l2
  goto 50

end
