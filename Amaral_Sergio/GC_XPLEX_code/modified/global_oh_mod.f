! $Id: global_oh_mod.f,v 1.2 2012/03/01 22:00:27 daven Exp $
      MODULE GLOBAL_OH_MOD
!
!******************************************************************************
!  Module GLOBAL_OH_MOD contains variables and routines for reading the
!  global monthly mean OH concentration from disk. (bmy, 7/28/00, 10/3/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) OH (TYPE (XPLEX))       : stores global monthly mean OH field
!  
!  Module Routines:
!  ============================================================================
!  (1 ) GET_OH            : Wrapper for GET_GLOBAL_OH
!  (2 ) GET_GLOBAL_OH     : Reads global monthly mean OH from disk
!  (3 ) INIT_GLOBAL_OH    : Allocates & initializes the OH array
!  (4 ) CLEANUP_GLOBAL_OH : Deallocates the OH array
!
!  GEOS-CHEM modules referenced by global_nox_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f  : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f  : Module containing NaN and other error-check routines
!
!  NOTES:
!  (1 ) Updated comments (bmy, 9/4/01)
!  (2 ) Now use routines from "transfer_mod.f" to regrid OH to 30 levels
!        for reduced GEOS-3 grid.  Also size OH array properly. (bmy, 1/14/02)
!  (3 ) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (4 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (5 ) Now use updated OH fields (bmy, 10/2/02)
!  (6 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (7 ) Minor bug fixes in FORMAT statements (bmy, 3/23/03)
!  (8 ) Cosmetic changes to simplify output (bmy, 3/27/03)
!  (9 ) Bug fix: OH should be (IIPAR,JJPAR,LLPAR) (bmy, 5/4/04)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!     
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean OH field
      TYPE (XPLEX), ALLOCATABLE :: OH(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_OH( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_OH reads global OH from binary punch files stored
!  in the /data/ctm/GEOS_MEAN directory.  This OH data is needed as oxidant
!  for various chemistry mechanisms (HCN, Tagged CO, etc...)  
!  (bmy, 7/28/00, 10/3/05)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) GET_GLOBAL_OH assumes that we are reading global OH data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate OH files.
!  (2 ) Now use version of GET_TAU0 with 3 arguments.  Now call READ_BPCH2
!        with IGLOB,JGLOB,LGLOB.  Call TRANSFER_3D to cast from TYPE (XPLEX) to
!        TYPE (XPLEX) and to regrid to 30 levels for GEOS-3 (if necessary).
!        ARRAY should now be of size (IGLOB,JGLOB,LGLOB). (bmy, 1/11/02)
!  (3 ) Now point to new OH files in the v4-26 subdirectory.  Also eliminated
!        obsolete code from 11/01. (bmy, 2/27/02)
!  (4 ) Now point to OH files in the v4-33 subdirectory. (bmy, 10/2/02)
!  (5 ) Replace missing commas in the FORMAT statement (bmy, 3/23/03)
!  (6 ) Cosmetic changes to simplify output (bmy, 3/27/03)
!  (7 ) Add Mat's OH as an option.  Also read bpch file quietly (bmy, 5/4/04)
!  (8 ) Now use OH_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : OH_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE TIME_MOD,      ONLY : GET_YEAR
      USE TRACER_MOD,    ONLY : ITS_A_CH4_SIM

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH

      ! Local variables
      INTEGER              :: I, J, L, THISYEAR
      TYPE (XPLEX)               :: ARRAY(IGLOB,JGLOB,LGLOB)
      TYPE (XPLEX)               :: ARRAY2(IIPAR,JJPAR,LGLOB)
      TYPE (XPLEX)               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_OH begins here!
      !=================================================================

      ! Allocate OH array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_OH
         FIRST = .FALSE.
      ENDIF

      ! Filename
        ! dkh hack: hardwire to geos4.  The geos5 doesn't work becuase
        ! it assumed 72 vertical levels. 
      FILENAME = TRIM( OH_DIR ) // 'OH_3Dglobal.' // GET_NAME_EXT() // 
     &                              '.'           // GET_RES_EXT()


      ! (kjw, dkh, 02/12/12, adj32_023) 
      IF ( ITS_A_CH4_SIM() ) THEN 
#if defined( GRID05x0666 ) && defined( NESTED_NA )
               FILENAME = '/home/kjw/GEOS-Chem/files/OH/'     //
     &                    'OH_3Dglobal.geos5.05x0666_NA'
#endif
      ENDIF 

!      FILENAME = TRIM( OH_DIR ) // 'OH_3Dglobal.geos4' //
!     &                              '.'           // GET_RES_EXT()
!          print*, 'dkh hack: use GEOS_4 OH fields!!!'
!          print*, 'dkh hack: use GEOS_4 OH fields!!!'
!          print*, 'dkh hack: use GEOS_4 OH fields!!!'
!          print*, 'dkh hack: use GEOS_4 OH fields!!!'

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_OH: Reading OH from: ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      THISYEAR = GET_YEAR()

      ! (kjw, dkh, 02/12/12, adj32_023) 
      IF ( ITS_A_CH4_SIM() ) THEN 
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
#if   defined( GRID05x0666 ) && defined( NESTED_NA )    
         FILENAME = '/home/kjw/GEOS-Chem/files/OH/'        //
     &              'OH_3Dglobal.geos5.05x0666_NA'
      
         CALL READ_BPCH2( TRIM(FILENAME), 'CHEM-L=$', 1,
     &                    XTAU,      IIPAR,     JJPAR,      
     &                    LGLOB,     ARRAY2,    QUIET=.FALSE.)
      
         ! Assign data from ARRAY to the module variable OH
         CALL TRANSFER_3D( ARRAY2, OH )

#else 
         ! Read OH data from the binary punch file
         CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 1,    
     &                    XTAU,      IGLOB,     JGLOB,       
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )
      
         ! Assign data from ARRAY to the module variable OH
         CALL TRANSFER_3D( ARRAY, OH )

#endif

      ELSE 
#if defined( GEOS_5 ) 
         XTAU = GET_TAU0( THISMONTH, 1, THISYEAR )  !(zhe 11/28/10)
#else
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
#endif  

         ! Read OH data from the binary punch file
         CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 1,     
     &                    XTAU,      IGLOB,     JGLOB,      
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

         ! Assign data from ARRAY to the module variable OH
         CALL TRANSFER_3D( ARRAY, OH )

      ENDIF  

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_OH

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_OH
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_OH allocates and zeroes the OH array, which holds 
!  global monthly mean OH concentrations. (bmy, 7/28/00, 5/4/04)
!
!  NOTES:
!  (1 ) OH array now needs to be sized (IGLOB,JGLOB,LGLOB) (bmy, 1/14/02)
!  (2 ) Also eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (3 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (4 ) OH should be (IIPAR,JJPAR,LLPAR): avoid subscript errors (bmy, 5/4/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_OH begins here!
      !=================================================================

      ! Allocate OH array
      ALLOCATE( OH( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH' )

      ! Zero OH array
      OH = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_OH
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_OH
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_OH deallocates the OH array. (bmy, 7/28/00)
!
!  NOTES:
!******************************************************************************
!        
      !=================================================================
      ! CLEANUP_GLOBAL_OH begins here!
      !=================================================================
      IF ( ALLOCATED( OH ) ) DEALLOCATE( OH ) 
     
      ! Return to calling program
      END SUBROUTINE CLEANUP_GLOBAL_OH

!------------------------------------------------------------------------------

      END MODULE GLOBAL_OH_MOD
