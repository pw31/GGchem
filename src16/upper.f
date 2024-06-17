      subroutine upper(str)
      implicit none
      character(*),intent(inout) :: str
      integer :: i,j,l
      logical :: change
      change = .true.
      l = len(str)
      do i=1,l
        if (i<l-4) then 
          if (str(i:i+4)=='trans') change=.false. 
        endif
        if (i<l-2) then 
          if (str(i:i+2)=='cis'  ) change=.false. 
        endif
        j = iachar(str(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") .and. change) then
          str(i:i) = achar(iachar(str(i:i))-32)
        else
          str(i:i) = str(i:i)
        end if
      enddo
      end subroutine upper

      !function upper(strIn) result(strOut)
      !implicit none
      !character(len=*),intent(in) :: strIn
      !character(len=len(strIn)) :: strOut
      !integer :: i,j,l
      !logical :: change
      !change = .true.
      !l = len(strIn)
      !do i=1,l
      !  if (i<l-4) then 
      !    if (strIn(i:i+4)=='trans') change=.false. 
      !  endif  
      !  if (i<l-2) then 
      !    if (strIn(i:i+2)=='cis'  ) change=.false. 
      !  endif  
      !  j = iachar(strIn(i:i))
      !  if (j>= iachar("a") .and. j<=iachar("z") .and. change) then
      !    strOut(i:i) = achar(iachar(strIn(i:i))-32)
      !  else
      !    strOut(i:i) = strIn(i:i)
      !  end if
      !enddo
      !end
