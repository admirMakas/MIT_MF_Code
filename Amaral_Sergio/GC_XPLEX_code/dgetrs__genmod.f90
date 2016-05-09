        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGETRS__genmod
          INTERFACE 
            SUBROUTINE DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              TYPE (XPLEX) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              TYPE (XPLEX) :: B(LDB,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGETRS
          END INTERFACE 
        END MODULE DGETRS__genmod
