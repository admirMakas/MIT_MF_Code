        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DLASWP__genmod
          INTERFACE 
            SUBROUTINE DLASWP(N,A,LDA,K1,K2,IPIV,INCX)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              TYPE (XPLEX) :: A(LDA,*)
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: K2
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE DLASWP
          END INTERFACE 
        END MODULE DLASWP__genmod
