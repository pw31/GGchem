***********************************************************************
*       New Mie subroutine that approximates for big grains           *
*                                                                     *
***********************************************************************
      SUBROUTINE Q_MIE(E1,E2,LAM,RAD,QEX,QSC,QAB,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION,INTENT(IN) ::  LAM,RAD,E1,E2
      DOUBLE PRECISION,INTENT(OUT) :: QEX,QSC,QAB,G
      DOUBLE PRECISION :: T
C
C     MIE THEORY EFFICIENCY FACTORS FOR SPHERICAL PARTICLES OF
C     RADIUS 'RAD' AT WAVELENGTH 'LAM'.
C     E=E1 + I*E2 IS THE SQUARE OF THE COMPLEX REFRACTIVE INDEX.
C     THE REFRACTIVE INDEX IS GIVEN BY SUBROUTINE 'EPS'
C
C     MOVE AN out of the stack to avoid problem in openmp version
C     http://gcc.gnu.org/onlinedocs/gfortran/OpenMP.html
C 
      COMPLEX*16,ALLOCATABLE :: AN(:)
C
      COMPLEX*16 E,RM,Y,ZN,ZN1,ZN2,C,A,B,AO,RRAT,A1,ANM1,BNM1
      ALLOCATE(AN(700000))
      T=0.d0
      X=6.2831853*RAD/LAM
      E=DCMPLX(E1,-E2)
      E=E**2.
      IF(X.LT.0.001)THEN
C
C        USE SMALL PARTICLE FORMULAE.
C        CHANGED CRITERIION FROM X < 0.01 TO 0.001 BECAUSE SILICATE
C        SCATTERING WAS NOT CORRECT.
C        15-8-2001: Changed scattering formula from
C        QSC=(X**4/.375)*DBLE(C**2)
C        into the correct formula QSC=(X**4/.375)*DABS(C)**2
C        Michiel Min
C
         C=(E-1.)/(E+2.)
         QSC=(X**4/.375)*ABS(C)**2
         A=DIMAG(-4.*C)
         B=DIMAG(-C*(E*E+27.*E+38.)/(2.*E+3.)/3.75)
         QAB=X*(A+X*X*B)
         QEX=QAB+QSC
C
C        G THE ASYMMETRY PARAMETER IS ALWAYS NEGLIGIBLE FOR SMALL PARTICLES.
C
         G=0.0
         RETURN
      END IF
C
C     FULL MIE THEORY CALCULATION.
C     RM - COMPLEX REFRACTIVE INDEX
C
      RM=CDSQRT(E)
      EN1=DBLE(RM)
      EN2=DIMAG(RM)
      Y=X*RM
      ZN2=DCMPLX(DCOS(X),-DSIN(X))
      ZN1=DCMPLX(DSIN(X),DCOS(X))
      RIND=EN1**2+EN2**2     ! Rind = |rm|²
      NTIL=1.5*SQRT(RIND)*X+1
      NTOT=MAX0(20,NTIL)

c     Number of iterations changed to improve for small |m| (Michiel Min)
      if(real(ntil).lt.(1.5*x))ntil=1.5*x
c
      if (ntot.le.700000) then    ! go ahead with full Mie theory
c
      AN(NTOT)=DCMPLX(0,0)
      SUME=0.
      SUMS=0.
      SUMG1=0.
      SUMG2=0.
      PSG1=0.
      PSG2=0.
      NTOTA=NTOT
  100 P=DFLOAT(NTOTA)
      AN(NTOTA-1)=P/Y-(1./(P/Y+AN(NTOTA)))
      NTOTA=NTOTA-1
      IF(NTOTA.EQ.1) GOTO 101
      GOTO 100
  101 AO1=DSIN(EN1*X)*DCOS(EN1*X)
      EN2P=-EN2
c      IF(EN2P*X.GE.44.)WRITE(6,*)'EN2P,X,LAM,RAD,E1,E2',EN2P,X,LAM,
c     >RAD,E1,E2
      if(EN2P*X.GE.350.) then
         AO=dcmplx(0.0,1.0)
      else
        AO2=DSINH(EN2P*X)*DCOSH(EN2P*X)
        AO3=(DSIN(EN1*X))**2+(DSINH(EN2P*X))**2
        AO=DCMPLX(AO1,AO2)
        AO=AO/AO3
      endif
      A1=-1./Y+(1./(1./Y-AO))
      RRAT=A1/AN(1)
      f=2.0/(x*x)
      DO 4 N=1,NTOT
         AN(N)=AN(N)*RRAT
    4 CONTINUE 
      DO 2 N=1,NTOT
         P=DFLOAT(N)
         ZN=DFLOAT(2*N-1)*ZN1/X-ZN2
         C=AN(N)/RM+P/X
         A=C*DBLE(ZN)-DBLE(ZN1)
         A=A/(C*ZN-ZN1)
         C=RM*AN(N)+P/X
         B=C*DBLE(ZN)-DBLE(ZN1)
         B=B/(C*ZN-ZN1)
C
C        PP, PPG1, PPG2 ARE CONSTANTS CONTAINING THE N TERMS IN THE 
C        SUMMATIONS.
C
         PP=DFLOAT(2*N+1)
C
         PSS=PP*(A*dCONJG(A)+B*dCONJG(B))
         PSE=PP*DBLE(A+B)
         IF(N.GT.1)THEN
C
C           CALCULATE G USING FORMULA ON P.128 OF VAN DE HULST'S BOOK.
C           HAVE REPLACED N BY (N-1) IN THE FORMULA SO THAT WE CAN USE
C           PREVIOUS A(N) AND B(N) INSTEAD OF A(N+1) AND B(N+1)
C
            REN=DFLOAT(N)
            PPG1=(REN-1.)*(REN+1.)/REN
            PPG2=(2.*REN-1.)/((REN-1.)*REN)
            PSG1=PPG1*DBLE(ANM1*dCONJG(A)+BNM1*dCONJG(B))
            PSG2=PPG2*DBLE(ANM1*dCONJG(BNM1))
         END IF
         SUME=SUME+PSE
         SUMS=SUMS+PSS
         SUMG1=SUMG1+PSG1
         SUMG2=SUMG2+PSG2
         D1=ABS(PSE/SUME)
         D2=ABS(PSS/SUMS)
C        IF(D1.LT.1.E-7.AND.D2.LT.1.E-7) GO TO 5
         PT=ABS(PSS/PP)
         IF(PT.LE.1.E-20) GOTO 5
C
C        SAVE PREVIOUS A AND B FOR CALCULATION OF G THE ASYMMETRY PARAMETER
C
         ANM1=A
         BNM1=B
         ZN2=ZN1
         ZN1=ZN
    2 CONTINUE
    5 F=2.0/(X*X)
      QEX=F*SUME
      QSC=F*SUMS
      QAB=F*(SUME-SUMS)
      G=2.0*F*(SUMG1+SUMG2)/QSC
      RETURN
      else
c               Geometrical optics for big spheres
      call geopt(rm,ans)
      qex =2.0d0
      g=9.23d-01   !approx true for D&L silicate.......
      qsc=ans
      end if
      return
      END          
c******************************************************************************
      subroutine geopt(m,ans)
c      intgrates the reflection coefficient
c      trapezium rule integration from 0 to pi/2
      implicit real*8 (a-h,o-z)
      complex*16 m
      a=0.0d0
      b=1.570796327d0
      nstrip = 5000
      tot=0
      h=(b-a)/dfloat(nstrip)   !strip width
      tot=tot+0.5*ref(m,a)*h   !1st term
      do i=1,nstrip-1
       x=a+h*dfloat(i)
       tot=tot+ref(m,x)*h      !middle terms
      end do
      tot=tot+0.5*ref(m,b)*h   !last term
      ans=1.+2.*tot    !ans is Qsca
      return
      end
      
c******************************************************************************
                                                                               
      function ref(m,thetai)
c         Calculates Reflection coeffs
      implicit real*8 (a-h,o-z)
      complex*16 sinTHETAt,cosTHETAt ,m,rpll,rper          
      sinTHETAt=sin(THETAi)/m
      cosTHETAt=cdsqrt(1-(sinTHETAt*sinTHETAt))
c       r for E parallel to plane
      rpll = (cosTHETAt-m*cos(THETAi)) / (cosTHETAt+m*cos(THETAi))
c       r for E perp. to plane
      rper = (cos(THETAi)-m*cosTHETAt) / (cos(THETAi)+m*cosTHETAt)
C       R = ½(|rpll|²+|rper|²)
      R= (abs(rpll)*abs(rpll) + abs(rper)*abs(rper))/2.0
      ref=r*sin(THETAi)*cos(THETAi)
      return                                            
      end                 

C
C
      SUBROUTINE INTERP(X,Y,NPTS,NTERMS,XIN,YOUT)
      REAL*8 DELTAX,PROD,SUM,X(5000),Y(5000),XIN,YOUT
      REAL*8 DELTA(10),A(10)
      REAL*8 DENOM
**************************************************   
*     SEARCH FOR AN APPROPRIATE VALUE OF X(1)    *
**************************************************
   11 DO 19 I=1,NPTS
      IF (XIN-X(I)) 13,17,19
   13 I1=I-NTERMS/2
      IF(I1) 15,15,21
   15 I1=1
      GOTO 21
   17 YOUT=Y(I)
   18 GOTO 61
   19 CONTINUE
      I1=NPTS-NTERMS+1
   21 I2=I1+NTERMS-1
      IF (NPTS-I2) 23,31,31
   23 I2=NPTS
      I1=I2-NTERMS+1
   25 IF (I1) 26,26,31
   26 I1=1
   27 NTERMS=I2-I1+1
C
C  EVALUATE DEVIATIONS DELTA
C
   31 DENOM=X(I1+1)-X(I1)
      DELTAX=(XIN-X(I1))/DENOM
      DO 35 I=1,NTERMS
         IX=I1+I-1
   35 DELTA(I)=(X(IX)-X(I1)) / DENOM
**********************************************
*           ACCUMULATE COEFFICIENTS A        *
**********************************************
   40 A(1)=Y(I1)
   41 DO 50 K=2,NTERMS
         PROD=1.
         SUM=0.
         IMAX=K-1
         IXMAX=I1+IMAX
         DO 49 I=1,IMAX
            J=K-I
            PROD=PROD*(DELTA(K)-DELTA(J))
   49    SUM=SUM-A(J)/PROD
   50 A(K)=SUM+Y(IXMAX)/PROD
***********************************************
*         ACCUMULATE SUM OF EXPANSION         *
***********************************************
   51 SUM=A(1)
      DO 57 J=2,NTERMS
         PROD=1.
         IMAX=J-1
         DO 56 I=1,IMAX
   56    PROD=PROD*(DELTAX-DELTA(I))
   57 SUM=SUM+A(J)* PROD
   60 YOUT=SUM
   61 RETURN
      END                
