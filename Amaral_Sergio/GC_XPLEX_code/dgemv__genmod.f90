        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGEMV__genmod
          INTERFACE 
            SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              TYPE (XPLEX) :: ALPHA
              TYPE (XPLEX) :: A(LDA,*)
              TYPE (XPLEX) :: X(*)
              INTEGER(KIND=4) :: INCX
              TYPE (XPLEX) :: BETA
              TYPE (XPLEX) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DGEMV
          END INTERFACE 
        END MODULE DGEMV__genmod
