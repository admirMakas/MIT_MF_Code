        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGBTRS__genmod
          INTERFACE 
            SUBROUTINE DGBTRS(TRANS,N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO&
     &)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDAB
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KL
              INTEGER(KIND=4) :: KU
              INTEGER(KIND=4) :: NRHS
              TYPE (XPLEX) :: AB(LDAB,*)
              INTEGER(KIND=4) :: IPIV(*)
              TYPE (XPLEX) :: B(LDB,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGBTRS
          END INTERFACE 
        END MODULE DGBTRS__genmod
