
************************************************************************
      integer function stindex(array,dim,string)
************************************************************************
*                                                                      *
*     Diese Funktion sucht in einem Character array nach einem String  *
*     und gibt den zugehoerigen index zurueck.                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     I N P U T                                                        *
*     array   : character array                                        *
*     dim     : dimensionierung des arrays                             *
*     string  : character sting, der gesucht werden soll               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     O U T P U T                                                      *
*     stindex : der index des Elementes in array, das mit string       *
*               uebereinstimmt.                                        *
*               Falls string nicht gefunden wird, wird 0 zurueckgegeben*
*                                                                      *
************************************************************************
*     (c)       Carsten Dominik                   Do 11. Maerz 1993    *
*     Revision: Carsten Dominik                   Do 11. Maerz 1993    *
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
