! $Id: lai_mod.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      MODULE LAI_MOD
!
!******************************************************************************
!  Module LAI_MOD reads and stores AVHRR LAI for calculating MEGAN biogenic 
!  VOC emissions. (dsa, tmf, bmy, 10/20/05, 11/6/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) ISOLAI     (TYPE (XPLEX) ) : AVHRR LAI data for the current day
!  (2 ) MISOLAI    (TYPE (XPLEX) ) : AVHRR LAI data for the current month
!  (3 ) NMISOLAI   (TYPE (XPLEX) ) : AVHRR LAI data for the next month 
!  (4 ) PMISOLAI   (TYPE (XPLEX) ) : AVHRR LAI data for the previous month
!  (5 ) DAYS_BTW_M (INTEGER) : days btw the current & previous months for LAI
!
!  Module Routines:
!  ============================================================================
!  (1 ) READISOLAI           : Reads monthly AVHRR LAI data
!  (2 ) RDISOLAI             : Calls READISOLAI and interpolates to daily LAI
!  (8 ) INIT_LAI             : Allocate and initialize data array
!  (9 ) CLEANUP_LAI          : Deallocate data array
!
!  GEOS-CHEM modules referenced by megan_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f          : Module w/ routines for binary punch file I/O
!  (2 ) error_mod.f          : Module w/ I/O error and NaN check routines
!  (3 ) transfer_mod.f       : Module w/ routines to cast & resize arrays
!
!  References:
!  ============================================================================
!
!  NOTES:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (7/8/03).  Updated  
!        and modified for the standard code by May Fu (11/2004).
!  (2 ) MEGAN is currently locked to use AVHRR LAI data.  
!        The LAVHRRLAI logical switch controls whether the AVHRR LAI data 
!        is used for the GEIA inventory and dry deposition.
!  (3 ) Modifications for 0.5 x 0.667 nested grid.  Added routine 
!        READISOLAI_05x0666 to read finer-resolution data for GEOS-5 nested
!        grids. (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "lai_mod.f"
      !=================================================================

      ! PRIVATE module variables

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      TYPE (XPLEX), ALLOCATABLE :: ISOLAI(:,:)
      TYPE (XPLEX), ALLOCATABLE :: MISOLAI(:,:)
      TYPE (XPLEX), ALLOCATABLE :: NMISOLAI(:,:)
      TYPE (XPLEX), ALLOCATABLE :: PMISOLAI(:,:)
      INTEGER             :: DAYS_BTW_M

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READISOLAI( MM )
!
!******************************************************************************
!  Subroutine READISOLAI reads AVHRR LAI data from bpch file for the current 
!  month, the previous month, and the next month. (dsa, tmf, bmy, 10/18/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MM (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (7/8/03).  Updated  
!        and modified for the standard code by May Fu (11/2004).
!******************************************************************************
!
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2, GET_RES_EXT
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
    
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters
   
      ! Arguments
      INTEGER, INTENT(IN)    :: MM

      ! Local variables
      INTEGER                :: I, J, K, INDEX, MMM, PMM, IJLOOP
      TYPE (XPLEX)                 :: ARRAY(I1x1,J1x1,1)
      TYPE (XPLEX)                 :: TAU0
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READISOLAI begins here!
      !=================================================================     

      ! Zero arrays
      MISOLAI  = 0.d0
      NMISOLAI = 0.d0
      ARRAY    = 0.d0

      !------------------------------------
      ! Read current month's lai at (I,J) 
      !------------------------------------
      
      ! Filename
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'leaf_area_index_200412/avhrrlai.global.geos.1x1.2000'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - READISOLAI: Reading ', a )

      ! Get TAU0 value
      TAU0 = GET_TAU0( MM, 1, 2000 )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      CALL DO_REGRID_1x1( 'cm2/cm2', ARRAY, MISOLAI )

      !------------------------------------
      ! Read next month's lai at (I,J) 
      !------------------------------------

      ! MMM is next month
      MMM = MM + 1
      IF ( MMM == 13 ) MMM = 1

      ! TAU0 for 1st day of next month
      TAU0 = GET_TAU0( MMM, 1, 2000 ) 
      
      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      CALL DO_REGRID_1x1( 'cm2/cm2', ARRAY, NMISOLAI )

      !------------------------------------
      ! Read previous month's lai at (I,J) 
      !------------------------------------

      ! PMM is previous month
      PMM = MM - 1
      IF ( PMM == 0 ) PMM = 12

      ! TAU0 for 1st day of previous month
      TAU0 = GET_TAU0( PMM, 1, 2000 ) 
      
      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      CALL DO_REGRID_1x1( 'cm2/cm2', ARRAY, PMISOLAI )

      ! Return to calling program
      END SUBROUTINE READISOLAI

!------------------------------------------------------------------------------

      SUBROUTINE READISOLAI_05x0666( MM )
!
!******************************************************************************
!  Subroutine READISOLAI reads AVHRR LAI data from bpch file for the current 
!  month, the previous month, and the next month.  Specially constructed to
!  read hi-res data for the GEOS-5 0.5 x 0.666 nested grid simulations.
!  (yxw, bmy, dan, 11/6/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MM (INTEGER) : Current month number (1-12)
!
!  NOTES:
!******************************************************************************
!
      ! Modules
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2, GET_RES_EXT
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      USE DIRECTORY_MOD,  ONLY : DATA_DIR

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)     :: MM

      ! Local variables
      INTEGER                 :: I, J, K, INDEX, MMM, PMM, IJLOOP
      TYPE (XPLEX)                  :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                  :: TAU0
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READISOLAI begins here!
      !=================================================================     

      ! Zero arrays
      MISOLAI  = 0.d0
      NMISOLAI = 0.d0
      ARRAY    = 0.d0

      !------------------------------------
      ! Read current month's lai at (I,J) 
      !------------------------------------

      ! Filename
      FILENAME = TRIM( DATA_DIR ) //
     &      'leaf_area_index_200412/avhrrlai.global.geos.05x0666.2000'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - READISOLAI: Reading ', a )

      ! Get TAU0 value
      TAU0 = GET_TAU0( MM, 1, 2000 )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,
     &                 TAU0,      IGLOB,   JGLOB,
     &                 1,         ARRAY,  QUIET=.TRUE. )


      CALL TRANSFER_2D( ARRAY(:,:,1), MISOLAI )


      !------------------------------------
      ! Read next month's lai at (I,J) 
      !------------------------------------

      ! MMM is next month
      MMM = MM + 1
      IF ( MMM == 13 ) MMM = 1

      ! TAU0 for 1st day of next month
      TAU0 = GET_TAU0( MMM, 1, 2000 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,
     &                 TAU0,      IGLOB,   JGLOB,
     &                 1,         ARRAY,  QUIET=.TRUE. )

      CALL TRANSFER_2D( ARRAY(:,:,1), NMISOLAI )


      !------------------------------------
      ! Read previous month's lai at (I,J) 
      !------------------------------------

      ! PMM is previous month
      PMM = MM - 1
      IF ( PMM == 0 ) PMM = 12

      ! TAU0 for 1st day of previous month
      TAU0 = GET_TAU0( PMM, 1, 2000 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,
     &                 TAU0,      IGLOB,   JGLOB,
     &                 1,         ARRAY,  QUIET=.TRUE. )

      CALL TRANSFER_2D( ARRAY(:,:,1), PMISOLAI )


      ! Return to calling program
      END SUBROUTINE READISOLAI_05x0666

!------------------------------------------------------------------------------

      SUBROUTINE RDISOLAI( JDAY, MONTH )
!
!******************************************************************************
!  Subroutine RDISOLAI sets ISOLAI daily.  The stored monthly LAI are used for
!  the middle day in the month and LAIs are interpolated for other days.
!  (dsa, tmf, bmy, 10/20/05, 11/6/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) JDAY  (INTEGER) : Julian Day
!  (2 ) MONTH (INTEGER) : Calendar month JDAY is in.  
!
!  NOTES:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (7/8/03).  Updated  
!        and modified for the standard code by May Fu (11/2004).
!  (2 ) Now call READISOLAI_05x0666 to read hi-res LAI data if we are doing a 
!        GEOS-5 0.5 x 0.666 nested grid simulation. (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,   ONLY : ITS_A_LEAPYEAR

#     include "CMN_SIZE"   ! Size parameters
   
      ! Arguments
      INTEGER, INTENT(IN) :: JDAY, MONTH

      ! Local variables
      INTEGER             :: I, J, IMUL, ITD, IJLOOP, MM
      INTEGER, SAVE       :: LAST_MM = -1 
      TYPE (XPLEX)              :: FRACTION

      ! specify midmonth day for year 2000
      INTEGER, PARAMETER  :: STARTDAY(13) = 
     &                         (/  15,  45,  74, 105, 135, 166,
     &                            196, 227, 258, 288, 319, 349, 380/)

      !=================================================================
      ! RDISOLAI begins here!
      !=================================================================

      ! Find the month if we index by midmonth
      CALL FINDMON( JDAY, MONTH, MM, STARTDAY )

      ! Read new data if it's a new LAI month
      IF ( MM /= LAST_MM ) THEN

#if   defined( GRID05x0666 )
         CALL READISOLAI_05x0666( MM )   ! GEOS-5 nested grid simulation
#else 
         CALL READISOLAI( MM )           ! Global simulations
#endif

         ! Save for next month
         LAST_MM = MM
      ENDIF
      
      ! IMUL is days since midmonth
      ! ITD  is days between midmonths
      IF ( JDAY < STARTDAY(1) ) THEN
         IMUL = 365 + JDAY - STARTDAY(12) 
         ITD  = 31
      ELSE
         IMUL = JDAY           - STARTDAY(MM)
         ITD  = STARTDAY(MM+1) - STARTDAY(MM)
      END IF

      ! Archive the days between midmonths in the LAI data
      DAYS_BTW_M     = ITD

      ! Fraction of the LAI month that we are in
      FRACTION       = DBLE( IMUL ) / DBLE( ITD ) 
       
      ! Interpolate to daily LAI value
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ISOLAI(I,J) = MISOLAI(I,J) + 
     &                ( FRACTION * ( NMISOLAI(I,J) - MISOLAI(I,J) ) )
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE RDISOLAI

!------------------------------------------------------------------------------

      SUBROUTINE INIT_LAI
!
!******************************************************************************
!  Subroutine INIT_ISOLAI allocates and initializes arrays for AVHRR LAI.
!  (dsa, tmf, 7/8/03, 11/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY    : ALLOC_ERR

#     include "CMN_SIZE"     ! Size parameters

      ! Local Variables
      INTEGER               :: AS

      !=================================================================
      ! INIT_ISOLAI begins here!
      !=================================================================

      ALLOCATE( ISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ISOLAI' )
      ISOLAI = 0d0

      ALLOCATE( MISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MISOLAI' )
      MISOLAI = 0d0

      ALLOCATE( NMISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NMISOLAI' )
      NMISOLAI = 0d0

      ALLOCATE( PMISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PMISOLAI' )
      PMISOLAI = 0d0        

      ! Return to calling program
      END SUBROUTINE INIT_LAI

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_LAI
!
!******************************************************************************
!  Subroutine CLEANUP_ISOLAI deallocates all allocated arrays at the
!  end of a GEOS-CHEM model run. (dsa 7/8/03)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_ISOLAI begins here!
      !=================================================================
      IF ( ALLOCATED( ISOLAI   ) ) DEALLOCATE( ISOLAI   )
      IF ( ALLOCATED( MISOLAI  ) ) DEALLOCATE( MISOLAI  )
      IF ( ALLOCATED( NMISOLAI ) ) DEALLOCATE( NMISOLAI )
      IF ( ALLOCATED( PMISOLAI ) ) DEALLOCATE( PMISOLAI )

      ! Return to calling program
      END SUBROUTINE CLEANUP_LAI

!------------------------------------------------------------------------------

      ! End of module
      END MODULE LAI_MOD
