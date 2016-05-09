        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGEMM__genmod
          INTERFACE 
            SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,&
     &C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: TRANSB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              TYPE (XPLEX) :: ALPHA
              TYPE (XPLEX) :: A(LDA,*)
              TYPE (XPLEX) :: B(LDB,*)
              TYPE (XPLEX) :: BETA
              TYPE (XPLEX) :: C(LDC,*)
            END SUBROUTINE DGEMM
          END INTERFACE 
        END MODULE DGEMM__genmod
