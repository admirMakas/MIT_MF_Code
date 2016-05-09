        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DGBTF2__genmod
          INTERFACE 
            SUBROUTINE DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
              INTEGER(KIND=4) :: LDAB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KL
              INTEGER(KIND=4) :: KU
              TYPE (XPLEX) :: AB(LDAB,*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGBTF2
          END INTERFACE 
        END MODULE DGBTF2__genmod
