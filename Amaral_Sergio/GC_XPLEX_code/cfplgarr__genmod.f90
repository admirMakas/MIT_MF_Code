        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE CFPLGARR__genmod
          INTERFACE 
            SUBROUTINE CFPLGARR(MAXMOM,L,M,X,CFPLG)
              INTEGER(KIND=4), INTENT(IN) :: MAXMOM
              INTEGER(KIND=4), INTENT(IN) :: L
              INTEGER(KIND=4), INTENT(IN) :: M
              TYPE (XPLEX), INTENT(IN) :: X
              TYPE (XPLEX), INTENT(OUT) :: CFPLG(0:MAXMOM)
            END SUBROUTINE CFPLGARR
          END INTERFACE 
        END MODULE CFPLGARR__genmod
