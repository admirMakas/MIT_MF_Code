C $Id: emmonot.f,v 1.1.1.1 2009/06/09 21:51:51 daven Exp $
      FUNCTION EMMONOT( IJLOOP, TMMP, XNUMOL )
!
!******************************************************************************
!  Subroutine EMMONOT computes the BIOGENIC MONOTERPENE EMISSIONS for each 
!  grid box in units of [atoms C/box/step]. (yhw, bdf, bmy, 9/4/01, 11/26/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IJLOOP    (INTEGER ) : 1-D grid box index
!  (2 ) TMMP      (TYPE (XPLEX)  ) : Local air temperature (K)
!  (3 ) XNUMOL    (TYPE (XPLEX)  ) : Number of atoms C / kg C 
!
!  Important Common Block Variables:
!  ============================================================================
!  (1 ) XYLAI     (CMN_VEL ) : Leaf Area Index of land type for current MONTH
!  (2 ) IJREG     (CMN_VEL ) : Number of Olson land types per grid box
!  (3 ) IJLAND+1  (CMN_VEL ) : Olson land type index
!  (4 ) IJUSE     (CMN_VEL ) : Olson land type fraction per box (in mils)
!  (5 ) BASEMONOT (CMN_ISOP) : Baseline MONOTERPENE emissions [kg C/box/step]
!
!  NOTES:
!  (1 ) Now use F90 syntax.  Use "D" exponents to force TYPE (XPLEX).
!        Updated comments, and mad cosmetic changes (bmy, 9/4/01) 
!  (2 ) Removed obsolete, commented-out code from 8/01 (bmy, 11/26/01)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_VEL"   ! XYLAI, IJREG, IJLAND, IJUSE
#     include "CMN_MONOT" ! BASEMONOT

      ! Arguments
      INTEGER, INTENT(IN) :: IJLOOP
      TYPE (XPLEX),  INTENT(IN) :: TMMP, XNUMOL

      ! Local variables
      INTEGER             :: INVEG
      TYPE (XPLEX)              :: EMMONOT 
      TYPE (XPLEX),PARAMETER::TS=xplex(303d0,0d0),BETA=xplex(0.09d0,0d0)

      !=================================================================
      ! EMMONOT begins here!
      !=================================================================

      ! Initialize
      EMMONOT = 0d0

      ! Loop over all Olson land types in this grid box
      DO INVEG = 1, IJREG(IJLOOP)

         ! Compute monoterpene emissions at box IJLOOP in [kg C/box/step].  
         ! Monoterpenes are now scaled to leaf area index XYLAI.  Also 
         ! multiply by the fraction of grid box IJLOOP occupied
         ! by this Olson land type. (bdf, bmy, 8/2/01)
         EMMONOT = EMMONOT + 
     &             ( BASEMONOT(IJLOOP,INVEG) * XYLAI(IJLOOP,INVEG) * 
     &               XPLX( IJUSE(IJLOOP,INVEG) ) / 1000d0 )

      ENDDO

      !=================================================================
      ! Temperature correction from Guenther et al. (1995)
      ! BETA is an empirical coefficient given by Guenther. (.09 K-1)
      ! TS is leaf temperature at standard conditions, (303 K)
      ! foliar density is accounted for in monotemis.table.
      !=================================================================

      ! Temp-corrected MONOTERPENE emissions in [kg C/box/step]
      EMMONOT = EMMONOT * EXP( BETA * ( TMMP - TS ) )

      ! Convert MONOTERPENE emissions to [atoms C/box/step]
      EMMONOT = EMMONOT * XNUMOL

      ! Return to calling program
      END FUNCTION EMMONOT
