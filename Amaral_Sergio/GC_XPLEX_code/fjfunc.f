! $Id: fjfunc.f,v 1.1.1.1 2009/06/09 21:51:53 daven Exp $
      TYPE (XPLEX) FUNCTION FJFUNC( I, J, L, NREAC, BRCH, NAME )
!
!*****************************************************************************
!  Subroutine FJFUNC supplies J-values to SMVGEAR solver.
!  (ppm, 4/98, bmy, 9/99, 10/15/02)
!
!  Arguments as input:
!  ===========================================================================
!  (1-3) I, J, L : Latitude, Longitude, Altitude indices of CTM grid box
!  (4  ) NREAC   : SMVGEAR photo reaction number (read from "chem.dat")
!  (5  ) BRCH    : SMVGEAR branch index (computed from "chem.dat") 
!  (6  ) NAME    : SMVGEAR species name (read from "chem.dat")
!
!  NOTES:
!  (1  ) "cmn_fj.h" also includes "CMN_SIZE" and "define.h".
!  (2  ) J-values are stored in array "ZPJ" from "cmn_fj.h".
!  (3  ) Now references ERROR_STOP from "error_mod.f".  Updated comments,
!         and made some cosmetic changes. (bmy, 10/15/02)
!*****************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP,GEOS_CHEM_STOP

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"  

      ! Arguments
      INTEGER,           INTENT(IN) :: I, J, L, NREAC, BRCH
      CHARACTER (LEN=4), INTENT(IN) :: NAME

      ! Local variables
      INTEGER                       :: N

      !=================================================================
      ! FJFUNC begins here!
      !
      ! If your compiler has subscript-range checking (-C or 
      ! -check_bounds) then it is recommended to use this option to 
      ! test for the validity of (I,J,L), since repeated IF statements 
      ! are computationally expensive.
      !
      ! If your compiler does not have subscript-range checking, then 
      ! uncomment the following lines to do a manual test for the 
      ! validity of (I,J,L).  
      !=================================================================
      !IF ( I > IPAR .OR. J > JPAR .OR. L > JPNL ) THEN 
      !   STOP 'invalid grid-box # in call to fjfunc - check fjfunc.f'
      !ENDIF

      !=================================================================
      ! RINDEX converts the J-value index as read from "chem.dat" to 
      ! the J-value index as read from "ratj.d". (bmy, 10/5/98)
      !
      ! Make sure that we have taken the proper reaction! 
      !=================================================================
      N = RINDEX(NREAC)

      IF ( N > JPPJ ) THEN
         WRITE(6,*) 'RXN for ',name,', branch ',brch,' not found!'
         CALL ERROR_STOP( 'Check FJFUNC.F', 'fjfunc.f' )
      ENDIF

      !=================================================================
      ! Return the appropriate J-value as the value of the function 
      !=================================================================
      FJFUNC = ZPJ(L,N,I,J)
!      if (isnan(zpj(L,n,I,j))) then
!         print*,'zpj is nan in fjfunc end',zpj(l,n,i,j)
!         CALL GEOS_CHEM_STOP
!      endif
      ! Return to calling program
      END FUNCTION FJFUNC
