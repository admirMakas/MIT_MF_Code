        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DTRSM__genmod
          INTERFACE 
            SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB&
     &)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              TYPE (XPLEX) :: ALPHA
              TYPE (XPLEX) :: A(LDA,*)
              TYPE (XPLEX) :: B(LDB,*)
            END SUBROUTINE DTRSM
          END INTERFACE 
        END MODULE DTRSM__genmod
