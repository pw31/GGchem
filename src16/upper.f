      function upper(strIn) result(strOut)
      implicit none
      character(len=*),intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j,l
      logical :: change

      change = .true.
      l = len(strIn)
      do i=1,l
        if (i<l-4) then 
          if (strIn(i:i+4)=='trans') change=.false. 
        endif  
        if (i<l-2) then 
          if (strIn(i:i+2)=='cis'  ) change=.false. 
        endif  
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") .and. change) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      enddo

      end
