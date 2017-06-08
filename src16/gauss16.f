**********************************************************************
           SUBROUTINE GAUSS16 (Ndim,N,a,x,b)
**********************************************************************
*****                                                            *****
*****   Diese Routine loesst ein lineares Gleichungssystem       *****
*****   der Form    (( a )) * ( x ) = ( b )     nach x auf.      *****
*****   Der Algorithmus funktioniert, indem die Matrix a         *****
*****   auf Dreiecksform gebracht wird.                          *****
*****                                                            *****
*****   EINGABE:  n = Dimension der Vektoren, der Matrix         *****
*****             a = (N x N)-Matrix                             *****
*****             b = (N)-Vektor                                 *****
*****   AUSGABE:  x = (N)-Vektor                                 *****
*****                                                            *****
**********************************************************************
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: Ndim,N,i,j,k,kmax
      real(kind=qp) :: a(Ndim,Ndim),x(Ndim),b(Ndim),c,amax
*
      do 500 i=1,N-1
*       ------------------------------------------
*       ***  MAX-Zeilentausch der i-ten Zeile  ***      
*       ------------------------------------------
        kmax = i
        amax = ABS(a(i,i))
        do 200 k=i+1,N
          if ( ABS(a(k,i)) .gt. amax ) then
            kmax = k
            amax = ABS(a(k,i))
          endif
  200   continue
        if (kmax.ne.i) then
          do 210 j=1,N
            c         = a(i,j)
            a(i,j)    = a(kmax,j)
            a(kmax,j) = c 
  210     continue
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
*
*       ---------------------------------
*       ***  bringe auf Dreiecksform  ***
*       ---------------------------------
        do 310 k=i+1,N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.Q0
          do 300 j=i+1,N
            a(k,j) = a(k,j) - c * a(i,j)
  300     continue        
          b(k) = b(k) - c * b(i)
  310   continue
*
  500 continue
*
*     --------------------------
*     ***  loese nach x auf  ***
*     --------------------------
      do 610 i=N,1,-1
        c = 0.Q0
        if (i.lt.N) then
          do 600 j=i+1,N
            c = c + a(i,j) * x(j)
  600     continue
        end if
        x(i) = (b(i) - c) / a(i,i)
  610 continue
      RETURN
      end
