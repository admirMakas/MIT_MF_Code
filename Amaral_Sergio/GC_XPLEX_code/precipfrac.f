! $Id: precipfrac.f,v 1.1.1.1 2009/06/09 21:51:53 daven Exp $
      SUBROUTINE PRECIPFRAC( I, J, RATE, FRAC )
!
!*****************************************************************************
!  Subroutine PRECIPFRAC computes the fraction of a grid box that is 
!  actually precipitating, along with the precipitation rate.
!  (djj, hyl, bmy, 10/18/99, 2/11/03)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) I (INTEGER)  : Longitude index of grid box 
!  (2 ) J (INTEGER)  : Latitude  index of grid box
!
!  Arguments as output:
!  ===========================================================================
!  (1) RATE (TYPE (XPLEX)) : Rate of precipitation for grid box (I,J) [mm/day  ]
!  (2) FRAC (TYPE (XPLEX)) : Fraction of grid box undergoing precip   [unitless]
!
!  Inputs passed via "CMN_PRECIP"
!  ===========================================================================
!  (3 ) PREACC (TYPE (XPLEX)) : DAO total      precipitation at ground [mm/day]
!  (4 ) PRECON (TYPE (XPLEX)) : DAO convective precipitation at ground [mm/day]
!
!  References:
!  ===========================================================================
!  Liu, H. Y., D. J. Jacob, I. Bey, R. M. Yantosca, and D. M. Koch,
!  Three-dimensional simulation of $210Pb$ and $7Be$ in the Harvard-DAO
!  tropospheric chemistry model, Eos Trans. AGU, 80 (17), S32, 1999a.
!
!  NOTES:
!  (1 ) PRECIPFRAC is written in Fixed-Form Fortran 90.
!  (2 ) This version of PRECIPFRAC replaces Yuhang Wang's original version,
!        as used in the GEOS-CTM prior to 10/18/99.
!  (3 ) Be sure to force TYPE (XPLEX) with the "D" exponent.
!  (4 ) Now reference PREACC, PRECON from "dao_mod.f" instead of from
!        common block header file "CMN_PRECIP".  (bmy, 6/26/00)
!  (5 ) Removed obsolete code from 6/26/00 (bmy, 8/31/00)
!  (6 ) Replaced JMX with JGLOB.  Updated comments, cosmetic changes.
!        (bmy, 6/25/02)
!  (7 ) Now use function GET_YOFFSET from "grid_mod.f" (bmy, 2/11/03)
!*****************************************************************************
!
      ! Reference to F90 modules
      USE DAO_MOD,  ONLY : PREACC, PRECON
      USE GRID_MOD, ONLY : GET_YOFFSET

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"   ! JGLOB

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J
      TYPE (XPLEX),  INTENT(OUT) :: RATE, FRAC

      ! Local variables
      INTEGER              :: JREF
      TYPE (XPLEX)               :: FRAC_LS, FRAC_CONV

      !=================================================================
      ! PRECIPFRAC begins here! 
      !
      ! For the polar boxes there is no precipitation.  
      ! Set RATE = 0, FRAC = 0 and return.
      !=================================================================
      JREF = J + GET_YOFFSET()

      IF ( JREF == 1 .OR. JREF == JGLOB ) THEN
         FRAC = 0.0d0
         RATE = 0.0d0
         RETURN
      ENDIF

      !=================================================================
      ! Large scale precipitation at (I,J) = PREACC(I,J) - PRECON(I,J).
      !
      ! If there is large-scale precipitation at grid box (I,J), then 
      ! assume that it covers 7% of the area of grid box(I,J).  Store 
      ! this value in the variable FRAC_LS.
      !=================================================================
      IF ( ( PREACC(I,J) - PRECON(I,J) ) > 0.0d0 ) THEN
         FRAC_LS = 7.0d-2
      ELSE
         FRAC_LS = 0.0d0
      ENDIF

      !=================================================================
      ! Convective precipitation at (I,J) = PRECON(I,:J)
      !
      ! If there is convective precipitation at (I,J), then 
      !  assume that it covers 0.3% of the area of grid box (I,J). 
      !  Store this value in the variable FRAC_CONV.
      !=================================================================
      IF ( PRECON(I,J) > 0.0d0 ) THEN
         FRAC_CONV = 3.0d-3
      ELSE
         FRAC_CONV = 0.0d0
      ENDIF

      !=================================================================
      ! FRAC = total fraction of grid box (I,J) covered by precip
      !      = FRAC_LS + FRAC_CONV
      !
      ! The possible values of FRAC are: 0.0%, 0.3%, 7.0%, or 7.3%.  
      !=================================================================
      FRAC = FRAC_LS + FRAC_CONV

      !=================================================================
      ! RATE = total precipitation rate in mm/day, adjusted for the 
      !        fraction of the grid box that is precipitating.  
      !
      ! To get RATE, take total precip at (I,J) and divide it by FRAC.
      !=================================================================
      IF ( FRAC > 0.0d0 ) THEN
         RATE = PREACC(I,J) / FRAC
      ELSE   
         RATE = 0.0d0
      ENDIF

      ! Return to calling program
      END SUBROUTINE PRECIPFRAC


