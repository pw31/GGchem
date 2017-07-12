      function upper(strIn) result(strOut)
      implicit none
      character(len=*),intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j
      logical :: change

      change = .true.
      do i = 1, len(strIn)
        if (strIn(i:i+4)=='trans') change=.false. 
        if (strIn(i:i+2)=='cis'  ) change=.false. 
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") .and. change) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      enddo

      end
