! $Id: tcorr.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      FUNCTION TCORR( TEMP )
!
!******************************************************************************
!  Function TCORR applies the temperature correction for isoprene emissions, 
!  according to Guenther et al.(92) (yhw, 11/15/93; bmy, 4/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TEMP (TYPE (XPLEX)) : Temperature [K]
!
!  References:
!  ============================================================================
!  Guenther et al, 1992, ... 
!
!  NOTES:
!  (1 ) Removed DATA statements, replaced w/ F90 syntax.  Updated comments
!        and made cosmetic changes (bmy, 4/4/03)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: TEMP

      ! Local variables
      TYPE (XPLEX), PARAMETER  :: R   = xplex(8.314d0,0d0)
      TYPE (XPLEX), PARAMETER  :: CT1 = xplex(95000.d0,0d0)
      TYPE (XPLEX), PARAMETER  :: CT2 = xplex(230000.d0,0d0)
      TYPE (XPLEX), PARAMETER  :: T1  = xplex(303.d0,0d0)
      TYPE (XPLEX), PARAMETER  :: T3  = xplex(314.d0,0d0)
      
      ! Function value
      TYPE (XPLEX)             :: TCORR
      
      !=================================================================
      ! TCORR begins here!
      !=================================================================
      TCORR =
     &     EXP( CT1/(R*T1*TEMP) * (TEMP-T1) ) /
     &     ( 1 + EXP( CT2/(R*T1*TEMP) * (TEMP-T3) ) )

      ! Return to calling program
      END FUNCTION TCORR
