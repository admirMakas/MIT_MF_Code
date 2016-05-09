        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGBTRF__genmod
          INTERFACE 
            SUBROUTINE DGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
              INTEGER(KIND=4) :: LDAB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KL
              INTEGER(KIND=4) :: KU
              TYPE (XPLEX) :: AB(LDAB,*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGBTRF
          END INTERFACE 
        END MODULE DGBTRF__genmod
