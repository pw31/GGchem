
************************************************************************
      integer function stindex(array,dim,string)
************************************************************************
      integer dim
      character*(*) array(dim),string
      integer i
      
      do i=1,dim
        if ( array(i) .eq. string ) then
          stindex = i
          return
        endif
      enddo
      write(6,*) 'not found: ' ,string
      stindex = 0
      end
