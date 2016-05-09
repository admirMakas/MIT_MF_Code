      MODULE PARTNER_EXTRA_MOD

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "aircraft_nox_mod.f"
      !=================================================================
     

      ! PRIVATE module variables

      TYPE (XPLEX),  ALLOCATABLE :: SUM_STT_ADJ(:,:,:,:) 
      TYPE (XPLEX),  ALLOCATABLE :: SUM_STT(:,:,:,:) 

       ! reading in input.aircraft
       LOGICAL            :: VERBOSE  = .FALSE.
       INTEGER, PARAMETER :: FIRSTCOL = 26
       INTEGER, PARAMETER :: MAXDIM   = 255
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!=========================================================================

      SUBROUTINE CUSTOM_WEIGHT
! Initializing variables
      
      USE ADJ_ARRAYS_MOD,    ONLY : WEIGHT
      USE DAO_MOD,           ONLY : IS_LAND
      USE FILE_MOD,          ONLY : IU_GEOS, IOERROR
      USE POPULATION_MOD,    ONLY : WGT_MORT 
#     include "CMN_SIZE"
     
      INTEGER   :: WGT_IMAX, WGT_JMAX, WGT_IMIN, WGT_JMIN
      LOGICAL   :: WGT_LAND
      INTEGER   :: I, J
      INTEGER   :: N, IOS
      LOGICAL  :: EOF
      CHARACTER(LEN=255) :: FILE_PARTNER
      CHARACTER(LEN=255) :: TITLE
      CHARACTER(LEN=255) :: SUBSTRS(255)
      
      FILE_PARTNER = 'input.partner'
      WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'P A R T N E R   I N P U T'
      WRITE( 6, 200   ) TRIM( FILE_PARTNER )
 200  FORMAT( 'READ_AIRCRAFT_INPUT_FILE: Reading ', a )

      WGT_IMAX = 0
      WGT_JMAX = 0  
      WGT_IMIN = 0 
      WGT_JMIN = 0
      OPEN( IU_GEOS, FILE=TRIM( FILE_PARTNER ), IOSTAT=IOS )
      IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS, IU_GEOS, 'read_partner_file:1' )
      PRINT *, 'test0'
      TITLE = READ_ONE_LINE( EOF  )
      PRINT *, 'test1'
      TITLE = READ_ONE_LINE( EOF  )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) WGT_IMIN
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) WGT_IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) WGT_JMIN
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) WGT_JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) WGT_LAND
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) WGT_MORT
    
      TITLE = READ_ONE_LINE( EOF  )

      CLOSE(IU_GEOS)
 
      PRINT *, 'WGT_IMIN, WGT_IMAX, WGT_JMIN, WGT_JMAX, WGT_LAND', 
     &      WGT_IMIN, WGT_IMAX, WGT_JMIN, WGT_JMAX, WGT_LAND 
      PRINT *, 'SUM WEIGHT', SUM(WEIGHT) 

      !cost function in sub area
      IF ( WGT_IMIN .ne. 0 .AND. WGT_IMIN .ne. 1) 
     &             WEIGHT(1:WGT_IMIN-1,:,:,:)=0d0
      IF ( WGT_IMAX .ne. 0 .AND. WGT_IMAX .ne. IIPAR) 
     &             WEIGHT(WGT_IMAX+1:IIPAR,:,:,:)=0d0
 
      IF ( WGT_JMIN .ne. 0 .AND. WGT_JMIN .ne. 1 ) 
     &             WEIGHT(:,1:WGT_JMIN-1,:,:)=0d0
      IF ( WGT_JMAX .ne. 0 .AND. WGT_JMAX .ne. JJPAR ) 
     &             WEIGHT(:,WGT_JMAX+1:JJPAR,:,:)=0d0

      ! only land has non zero weight
      IF ( WGT_LAND ) THEN
         DO J = 1, JJPAR
            DO I = 1, IIPAR
                IF ( .not. IS_LAND(I,J) ) THEN
                    WEIGHT(I,J,:,:)=0d0
                ENDIF
            ENDDO
         ENDDO
      ENDIF
      PRINT *, 'SUM WEIGHT', SUM(WEIGHT) 
!      CALL READ_PARNER_MENU
  
      ! Return to calling program
      END SUBROUTINE CUSTOM_WEIGHT

!-----------------------------------------------------------------------------

      SUBROUTINE INIT_PARTNER
! Initializing variables
      USE TRACER_MOD,    ONLY : N_TRACERS, TRACER_NAME
      USE ERROR_MOD,    ONLY : ALLOC_ERR
#     include "CMN_SIZE"
#     include "define.h"
      
      ! Local variables
      INTEGER :: AS
      INTEGER :: L
     
      ALLOCATE( SUM_STT_ADJ( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUM_STT_ADJ' )
      SUM_STT_ADJ = 0d0
      ALLOCATE( SUM_STT( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUM_STT' )
      SUM_STT = 0d0
 

!      CALL READ_PARNER_MENU
  
      ! Return to calling program
      END SUBROUTINE INIT_PARTNER

!------------------------------------------------------------------------
      SUBROUTINE SAVE_FINAL_OUTPUTS

      USE AIRCRAFT_ADJ_MOD,      ONLY : EMS_AC_ADJ
      USE DIRECTORY_ADJ_MOD,     ONLY : OPTDATA_DIR, DIAGADJ_DIR
      USE PARTNER_SAVE_MOD,      ONLY : SAVE_EMS_AC_FILE,
     &                                  SAVE_NETCDF_FILE,
     &                                  SAVE_NETCDF_FILE_X,
     &                                  SAVE_EMS_AC_FILE_X
      CHARACTER(LEN=255) :: FILENAME

         FILENAME='gdt_ac_ems.nc' !jkoo
         FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )
         CALL SAVE_EMS_AC_FILE( FILENAME, EMS_AC_ADJ)
         
         FILENAME='gdt_ac_ems_x.nc' !bvconsta
         FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )
         CALL SAVE_EMS_AC_FILE_X( FILENAME, EMS_AC_ADJ)         

         FILENAME='gdt_con.nc' !jko
         FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )
         CALL SAVE_NETCDF_FILE( FILENAME, SUM_STT_ADJ )
     
         FILENAME='gdt_con_x.nc' !bvconsta
         FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )
         CALL SAVE_NETCDF_FILE_X( FILENAME, SUM_STT_ADJ )

       IF ( ALLOCATED( SUM_STT_ADJ ) ) DEALLOCATE( SUM_STT_ADJ  )
       IF ( ALLOCATED( SUM_STT     ) ) DEALLOCATE( SUM_STT      )

      END SUBROUTINE SAVE_FINAL_OUTPUTS

!------------------------------------------------------------------------------

      FUNCTION READ_ONE_LINE( EOF, LOCATION ) RESULT( LINE )
!
!******************************************************************************
!  Subroutine READ_ONE_LINE reads a line from the input file.  If the global 
!  variable VERBOSE is set, the line will be printed to stdout.  READ_ONE_LINE
!  can trap an unexpected EOF if LOCATION is passed.  Otherwise, it will pass
!  a logical flag back to the calling routine, where the error trapping will
!  be done. (bmy, 7/20/04)
! 
!  Arguments as Output:
!  ===========================================================================
!  (1 ) EOF      (CHARACTER) : Logical flag denoting EOF condition
!  (2 ) LOCATION (CHARACTER) : Name of calling routine; traps premature EOF
!
!  Function value:
!  ===========================================================================
!  (1 ) LINE     (CHARACTER) : A line of text as read from the file
!
!  NOTES:
!******************************************************************************
!      
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_GEOS, IOERROR

      ! Arguments
      LOGICAL,          INTENT(OUT)          :: EOF
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: LOCATION

      ! Local variables
      INTEGER                                :: IOS
      CHARACTER(LEN=255)                     :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file
      READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.

         ! Trap unexpected EOF -- stop w/ error msg if LOCATION is passed
         ! Otherwise, return EOF to the calling program
         IF ( PRESENT( LOCATION ) ) THEN
            MSG = 'READ_ONE_LINE: error at: ' // TRIM( LOCATION )
            WRITE( 6, '(a)' ) MSG
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_ONE_LINE (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            STOP
         ELSE
            RETURN
         ENDIF
      ENDIF

      ! IO Status > 0: true I/O error condition
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_one_line:1' )

      ! Print the line (if necessary)
      IF ( VERBOSE ) WRITE( 6, '(a)' ) TRIM( LINE )

      ! Return to calling program
      END FUNCTION READ_ONE_LINE

!------------------------------------------------------------------------------

      SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, N_EXP, LOCATION ) 
!
!******************************************************************************
!  Subroutine SPLIT_ONE_LINE reads a line from the input file (via routine 
!  READ_ONE_LINE), and separates it into substrings. (bmy, 7/20/04)
!
!  SPLIT_ONE_LINE also checks to see if the number of substrings found is 
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori, 
!  you can skip the error check.
! 
!  Arguments as Input:
!  ===========================================================================
!  (3 ) N_EXP     (INTEGER  ) : Number of substrings we expect to find
!                               (N_EXP < 0 will skip the error check!)
!  (4 ) LOCATION  (CHARACTER) : Name of routine that called SPLIT_ONE_LINE
!
!  Arguments as Output:
!  ===========================================================================
!  (1 ) SUBSTRS   (CHARACTER) : Array of substrings (separated by " ")
!  (2 ) N_SUBSTRS (INTEGER  ) : Number of substrings actually found
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY: STRSPLIT
      
      ! Arguments
      CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM)
      INTEGER,            INTENT(OUT) :: N_SUBSTRS
      INTEGER,            INTENT(IN)  :: N_EXP
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION 

      ! Local varaibles
      LOGICAL                         :: EOF
      CHARACTER(LEN=255)              :: LINE, MSG

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================      

      ! Create error msg
      MSG = 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )

      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = READ_ONE_LINE( EOF )

      ! STOP on End-of-File w/ error msg
      IF ( EOF ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, '(a)' ) 'End of file encountered!' 
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong 
      IF ( N_EXP /= N_SUBSTRS ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
      ENDIF
       
      ! Return to calling program
      END SUBROUTINE SPLIT_ONE_LINE

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_PARNER
!
!******************************************************************************
!  Subroutine CLEANUP_AIRCRAFT deallocates module variables. (srhb, 08/27/09)
!  Added additional allocated variables (jkoo,03/02/09)
!******************************************************************************
!

      ! Return to calling program
      END SUBROUTINE CLEANUP_PARNER

!-----------------------------------------------------------------------------
      ! End of module
      END MODULE PARTNER_EXTRA_MOD
