        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE DSWAP__genmod
          INTERFACE 
            SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
              INTEGER(KIND=4) :: N
              TYPE (XPLEX) :: DX(*)
              INTEGER(KIND=4) :: INCX
              TYPE (XPLEX) :: DY(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DSWAP
          END INTERFACE 
        END MODULE DSWAP__genmod
