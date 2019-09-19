************************************************************************
      subroutine NASA_POLYNOMIAL
************************************************************************
      use CHEMISTRY,ONLY: NELEM,b_nasa,c_nasa
      implicit none
      integer*4 :: j,k
      

      c_nasa(1:NELEM) = 0
      open(10, file = 'data/Burcat_ref-elements.dat')

      do k=1,4
         read(10,*)
      end do
      !H as a reference state for H
      c_nasa(1) = 1
      read(10, '(5E15.8)') (b_nasa(1,j),j=0,4)
      read(10, '(5E15.8)') (b_nasa(1,j),j=5,9)
      read(10, '(4E15.8)') (b_nasa(1,j),j=10,13)
      
      do k=1,7
         read(10,*)
      end do
      !C as a reference state for C
      c_nasa(6) = 1
      read(10, '(5E15.8)') (b_nasa(6,j),j=0,4)
      read(10, '(5E15.8)') (b_nasa(6,j),j=5,9)
      read(10, '(4E15.8)') (b_nasa(6,j),j=10,13)
          
      do k=1,5
         read(10,*)
      end do
      !O as a reference state for O
      c_nasa(8) = 1
      read(10, '(5E15.8)') (b_nasa(8,j),j=0,4)
      read(10, '(5E15.8)') (b_nasa(8,j),j=5,9)
      read(10, '(4E15.8)') (b_nasa(8,j),j=10,13)

      do k=1,6
         read(10,*)
      end do
      !N as a reference state for N
      c_nasa(7) = 1
      read(10, '(5E15.8)') (b_nasa(7,j),j=0,4)
      read(10, '(5E15.8)') (b_nasa(7,j),j=5,9)
      read(10, '(4E15.8)') (b_nasa(7,j),j=10,13)
      
      close(10)
      end
