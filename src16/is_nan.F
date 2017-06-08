#ifndef IBM
!======================================================================
      LOGICAL FUNCTION IS_NAN(X)
!======================================================================
      IMPLICIT NONE
      REAL*8 X,Y
      INTEGER*4 I(2)
#ifdef IFORT
      EQUIVALENCE(Y,I)
      Y = X
      IS_NAN = ((I(2) .AND. Z'7ff00000') .EQ. Z'7ff00000') .AND.  
     &        (((I(2) .AND. Z'000fffff') .NE. Z'00000000') .OR.   
     &          (I(1) .NE. 0))
#else
      IS_NAN = (.not.(x.eq.x))
!     IS_NAN = ((I(1) .AND. Z'7ff00000') .EQ. Z'7ff00000') .AND. 
!    &        (((I(1) .AND. Z'000fffff') .NE. Z'00000000') .OR.   
!    &          (I(2) .NE. 0))
#endif
      RETURN
      END FUNCTION IS_NAN
#endif
