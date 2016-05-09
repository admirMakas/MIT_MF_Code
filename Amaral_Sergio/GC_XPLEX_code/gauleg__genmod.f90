        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE GAULEG__genmod
          INTERFACE 
            SUBROUTINE GAULEG(X1,X2,X,W,N)
              INTEGER(KIND=4), INTENT(IN) :: N
              TYPE (XPLEX), INTENT(IN) :: X1
              TYPE (XPLEX), INTENT(IN) :: X2
              TYPE (XPLEX), INTENT(OUT) :: X(N)
              TYPE (XPLEX), INTENT(OUT) :: W(N)
            END SUBROUTINE GAULEG
          END INTERFACE 
        END MODULE GAULEG__genmod
