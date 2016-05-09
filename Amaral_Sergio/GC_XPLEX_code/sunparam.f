C $Id: sunparam.f,v 1.1.1.1 2009/06/09 21:51:53 daven Exp $
      SUBROUTINE SUNPARAM(X)

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

C===============================================
C the sequence is lai,suncos,cloud fraction
C===============================================
C  NN = number of variables (lai,suncos,cloud fraction)
      INTEGER NN
      PARAMETER(NN=3)
C  ND = scaling factor for each variable
      INTEGER ND(NN),I
      DATA ND /55,20,11/
C  X0 = maximum for each variable
      TYPE (XPLEX) X(NN),X0(NN),XLOW
      DATA X0 /xplex(11.d0,0d0),xplex(1.d0,0d0),xplex(1.d0,0d0)/

      DO I=1,NN
        X(I)=MIN(X(I),X0(I))
C XLOW = minimum for each variable
        IF (I.NE.3) THEN
          XLOW=X0(I)/XPLX(ND(I))
        ELSE
          XLOW= 0.d0
        END IF
        X(I)=MAX(X(I),XLOW)
        X(I)=X(I)/X0(I)
      END DO

      RETURN
      END
