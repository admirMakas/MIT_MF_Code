        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DTBSV__genmod
          INTERFACE 
            SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              TYPE (XPLEX) :: A(LDA,*)
              TYPE (XPLEX) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE DTBSV
          END INTERFACE 
        END MODULE DTBSV__genmod
