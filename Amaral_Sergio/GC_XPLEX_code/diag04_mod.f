! $Id: diag04_mod.f,v 1.2 2010/05/07 20:39:47 daven Exp $
      MODULE DIAG04_MOD
!
!******************************************************************************
!  Module DIAG04_MOD contains arrays and routines for archiving the ND04
!  diagnostic -- CO2 emissions and fluxes (bmy, 7/26/05, 9/5/06) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD04         (TYPE (XPLEX)) : Array for 2-D CO2 emissions/uptake 
!  (2 ) AD04_plane   (TYPE (XPLEX)) : Array for 3-D CO2 emissions from aircraft 
!  (3 ) AD04_chem    (TYPE (XPLEX)) : Array for 3-D CO2 emissions from chemical oxidation 
!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG04           : Sets all module arrays to zero
!  (2 ) WRITE_DIAG04          : Writes data in module arrays to bpch file
!  (3 ) INIT_DIAG04           : Allocates all module arrays
!  (4 ) CLEANUP_DIAG04        : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag04_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f           : Module w/ NaN and other error check routines
!  (3 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f            : Module w/ horizontal grid information
!  (5 ) time_mod.f            : Module w/ routines to compute date & time
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (3 ) Modified for ship emissions (2-D), aircraft emissions (3-D) and 
!       chemical source for CO2 (3-D) (RayNassar, 2009-12-23)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag04_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND04, LD04
      INTEGER, PARAMETER   :: PD04 = 10

      ! Arrays
      TYPE (XPLEX),  ALLOCATABLE :: AD04(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD04_plane(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD04_chem(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG04
!
!******************************************************************************
!  Subroutine ZERO_DIAG04 zeroes the ND04 diagnostic array (bmy, 7/26/05)
!******************************************************************************
!
!      ! References to F90 modules

#     include "CMN_SIZE"  ! Size parameters

      !=================================================================
      ! ZERO_DIAG04 begins here!
      !=================================================================

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Zero 2-D array (for N=7 tracers) and 3-D plane and chem arrays
      AD04(:,:,:)       = 0d0
      AD04_plane(:,:,:) = 0d0
      AD04_chem(:,:,:)  = 0d0

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG04

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG04
!
!******************************************************************************
!  Subroutine WRITE_DIAG04 writes the ND04 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 7/26/05, 9/3/06)
!
!   # : Field     : Description                  : Units       : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) CO2-SRCE  : CO2 fossil fuel emissions    : molec/cm2/s : SCALE
!  (2 ) CO2-SRCE  : CO2 ocean emissions          : molec/cm2/s : SCALE
!  (3 ) CO2-SRCE  : CO2 balanced biosphere       : molec/cm2/s : SCALE
!  (4 ) CO2-SRCE  : CO2 biomass emissions        : molec/cm2/s : SCALE
!  (5 ) CO2-SRCE  : CO2 biofuel emissions        : molec/cm2/s : SCALE
!  (6 ) CO2-SRCE  : CO2 net terrestrial exchange : molec/cm2/s : SCALE
!  (7 ) CO2-SRCE  : CO2 ship emissions           : molec/cm2/s : SCALE
!  (8 ) CO2-SRCE  : CO2 aircraft emissions (3-D) : molec/cm2/s : SCALE
!  (9 ) CO2-SRCE  : CO2 chemical source (3-D)    : molec/cm2/s : SCALE
!  (10) CO2-SRCE  : CO2 chem source surf correct : molec/cm2/s : SCALE
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,  ONLY : IU_BPCH
      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,  ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! TINDEX

      ! Local variables
      INTEGER            :: CENTER180, HALFPOLAR, IFIRST, JFIRST 
      INTEGER            :: LFIRST,    LMAX,      M,      N       
      TYPE (XPLEX)             :: ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)             :: LONRES,    LATRES
      TYPE (XPLEX)             :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20)  :: MODELNAME 
      CHARACTER(LEN=40)  :: CATEGORY,  RESERVED,    UNIT

      !=================================================================
      ! WRITE_DIAG04 begins here!
      !=================================================================

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      SCALE     = DBLE( GET_CT_EMIS() ) + 1d-32

      !=================================================================
      ! Write data to the bpch file
	! Note: if any of the ARRAY or AD04* dimensions are wrong, the 
	! run will crash with "ERROR RUNNING GEOS-CHEM" at the end.
      !=================================================================

      ! Loop over ND04 diagnostic tracers
      DO M = 1, TMAX(4)

         ! Get quantities
         N            = TINDEX(4,M)
 
	   IF (N <= 7) THEN

         	CATEGORY     = 'CO2-SRCE'
         	UNIT         = 'molec/cm2/s'
         	!UNIT         = ''                     ! Let GAMAP pick the unit
         	LMAX = 1
		ARRAY(:,:,1) = AD04(:,:,N) / SCALE

	   ELSEIF (N == 8) THEN

         	CATEGORY     = 'CO2-SRCE'
         	UNIT         = 'molec/cm3/s'
         	LMAX = LD04
         	ARRAY(:,:,1:LMAX) = AD04_plane(:,:,1:LMAX) / SCALE

	   ELSEIF (N == 9) THEN

         	CATEGORY     = 'CO2-SRCE'
         	UNIT         = 'molec/cm3/s'
         	LMAX = LD04
		ARRAY(:,:,1:LMAX) = AD04_chem(:,:,1:LMAX) / SCALE
		
	   ELSEIF (N == 10) THEN

         	CATEGORY     = 'CO2-SRCE'
         	UNIT         = 'molec/cm2/s'
         	LMAX = 1
		ARRAY(:,:,1) = AD04(:,:,N) / SCALE

         ELSE 

            CYCLE

	   ENDIF

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LMAX,    IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )

      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG04

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG04
!
!******************************************************************************
!  Subroutine INIT_DIAG04 allocates all module arrays (bmy, 7/26/05)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
   
#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG04 begins here!
      !=================================================================

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Get number of levels for 3-D arrays
	LD04 = MIN( ND04, LLPAR )
	
      ! 2-D array ("CO2-SRCE")

!      ALLOCATE( AD04( IIPAR, JJPAR, PD04-2 ), STAT=AS )
      ALLOCATE( AD04( IIPAR, JJPAR, PD04 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD04' )

      ! 3-D arrays ("CO2-SRCE")

      ALLOCATE( AD04_plane( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD04_plane' )

      ALLOCATE( AD04_chem( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD04_chem' )

      ! Zero arrays
      CALL ZERO_DIAG04

      ! Return to calling program
      END SUBROUTINE INIT_DIAG04

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG04
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG04 deallocates all module arrays (bmy, 7/26/05)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG04 begins here!
      !=================================================================

      IF ( ALLOCATED( AD04       ) ) DEALLOCATE( AD04 ) 
      IF ( ALLOCATED( AD04_plane ) ) DEALLOCATE( AD04_plane )
      IF ( ALLOCATED( AD04_chem  ) ) DEALLOCATE( AD04_chem  ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG04

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG04_MOD
