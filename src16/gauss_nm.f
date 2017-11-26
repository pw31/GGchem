**********************************************************************
      SUBROUTINE GAUSS_NM(Ndim,Mdim,N,M,A,x,b,info)
**********************************************************************
*****  tries to solve a system of N equations for M unknowns     *****
*****                   A x = b                                  *****
*****  info = 0  means that a solution was found                 *****
*****  info = 1  no solution in case N>M                         *****
*****  info = 2  NaNs produced (which normally occurs if N<M)    *****
**********************************************************************
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: Ndim,Mdim,N,M
      real(kind=qp),intent(inout) :: A(Ndim,Mdim),b(Ndim)
      real(kind=qp),intent(out) :: x(Mdim)
      integer,intent(out) :: info
      integer :: i,j,k,kmax,D
      real(kind=qp) :: c,Amax

      D = min(N-1,M)

      do i=1,D

        !do k=1,N
        !  print'(99(1pE10.3))',A(k,1:M),b(k) 
        !enddo
        !print*

        !-------------------------------------
        !***  MAX-exchange of i-th column  ***      
        !-------------------------------------
        kmax = i
        Amax = ABS(A(i,i))
        do k=i+1,N
          if (ABS(A(k,i))>Amax) then
            kmax = k
            Amax = ABS(A(k,i))
          endif
        enddo  
        if (kmax.ne.i) then
          do j=1,M
            c         = A(i,j)
            A(i,j)    = A(kmax,j)
            A(kmax,j) = c 
          enddo  
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
        !-------------------------------
        !***  make triangular shape  ***
        !-------------------------------
        do k=i+1,N
          if (A(i,i)==0.Q0) then
            info = 2
            return
          endif  
          c = A(k,i) / A(i,i)
          A(k,i) = 0.Q0
          do j=i+1,M
            A(k,j) = A(k,j) - c * A(i,j)
          enddo  
          b(k) = b(k) - c * b(i)
        enddo  
      enddo  

      !do i=1,N
      !  print'(99(1pE10.3))',A(i,1:M),b(i) 
      !enddo
      !print*

      !-------------------
      !***  resolve x  ***
      !-------------------
      info = 0
      x(:) = 0.Q0
      do i=M,1,-1
        if (A(i,i)==0.Q0) then
          info = 2
          return
        endif  
        c = 0.Q0
        do j=i+1,M
          c = c + A(i,j) * x(j)
        enddo  
        x(i) = (b(i) - c) / A(i,i)
        if (ABS(x(i))>1.Q+25) info=2
      enddo  

      !print*,x(1:M)
      !print*

      if (N>M) then
        !--- more equations than unknowns --- 
        do i=M+1,N
          c = 0.Q0         
          do j=1,M
            c = c + A(i,j)*x(j)
          enddo
          !print'(2(1pE18.11))',c,b(i)
          if (ABS(c-b(i))>1.Q-25) info=1 
        enddo 
      endif   

      !if (info==0) then
      !  print*,'solution found.'
      !else if (info==1) then
      !  print*,'more eqs than unknowns, no linear combi'
      !else if (info==2) then
      !  print*,'more unknowns than eqs, no linear combi'
      !endif  

      end
