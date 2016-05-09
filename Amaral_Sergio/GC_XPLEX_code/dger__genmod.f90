        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGER__genmod
          INTERFACE 
            SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              TYPE (XPLEX) :: ALPHA
              TYPE (XPLEX) :: X(*)
              INTEGER(KIND=4) :: INCX
              TYPE (XPLEX) :: Y(*)
              INTEGER(KIND=4) :: INCY
              TYPE (XPLEX) :: A(LDA,*)
            END SUBROUTINE DGER
          END INTERFACE 
        END MODULE DGER__genmod
