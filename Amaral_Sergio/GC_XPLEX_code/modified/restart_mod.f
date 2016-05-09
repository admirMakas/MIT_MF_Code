! $Id: restart_mod.f,v 1.2 2012/03/01 22:00:27 daven Exp $
      MODULE RESTART_MOD
!
!******************************************************************************
!  Module RESTART_MOD contains variables and routines which are used to read
!  and write GEOS-CHEM restart files, which contain tracer concentrations
!  in [v/v] mixing ratio. (bmy, 6/25/02, 12/16/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) INPUT_RESTART_FILE   : Full path name of the restart file to be read
!  (2 ) OUTPUT_RESTART_FILE  : Full path name (w/ tokens!) of output file
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_RESTART_FILE    : Writes restart file to disk 
!  (2 ) READ_RESTART_FILE    : Reads restart file from disk 
!  (3 ) CONVERT_TRACER_TO_VV : Converts from [ppbv], [ppmv], etc to [v/v]
!  (4 ) CHECK_DIMENSIONS     : Ensures that restart file contains global data
!  (5 ) COPY_STT             : Converts [v/v] to [kg] and stores in STT
!  (6 ) CHECK_DATA_BLOCKS    : Makes sure we have read in data for each tracer
!  (7 ) SET_RESTART          : Gets restart filenames from "input_mod.f"
!
!  GEOS-CHEM modules referenced by restart_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f          : Module w/ routines for binary punch file I/O
!  (2 ) error_mod.f          : Module w/ NaN and other error check routines
!  (3 ) file_mod.f           : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f           : Module w/ horizontal grid information
!  (5 ) logical_mod.f        : Module w/ GEOS-CHEM logical switches
!  (6 ) time_mod.f           : Module w/ routines for computing time & date
!  (7 ) tracer_mod.f         : Module w/ GEOS-CHEM tracer array STT etc.
!
!  NOTES:
!  (1 ) Moved routines "make_restart_file.f"" and "read_restart_file.f" into
!        this module.  Also now internal routines to "read_restart_file.f"
!        are now a part of this module.  Now reference "file_mod.f" to get
!        file unit numbers and error checking routines. (bmy, 6/25/02)
!  (2 ) Now reference AD from "dao_mod.f".  Now reference "error_mod.f".
!        Also added minor bug fix for ALPHA platform. (bmy, 10/15/02)
!  (3 ) Now references "grid_mod.f" and the new "time_mod.f" (bmy, 2/11/03)
!  (4 ) Added error-check and cosmetic changes (bmy, 4/29/03)
!  (5 ) Removed call to COPY_STT_FOR_OX, it's obsolete (bmy, 8/18/03)
!  (6 ) Add fancy output (bmy, 4/26/04)
!  (7 ) Added routine SET_RESTART.  Now reference "logical_mod.f" and
!        "tracer_mod.f" (bmy, 7/20/04)
!  (8 ) Removed obsolete routines TRUE_TRACER_INDEX and COPY_DATA_FOR_CO_OH
!        (bmy, 6/28/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now pass TAU via the arg list in MAKE_RESTART_FILE (bmy, 12/15/05)
!  (11) Add MAKE_CSPEC_FILE and READ_CSPEC_FILE routines to save and read
!        CSPEC_FULL restart files (dkh, 02/12/09)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "restart_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC  :: MAKE_RESTART_FILE
      PUBLIC  :: READ_RESTART_FILE
      PUBLIC  :: SET_RESTART
      PUBLIC  :: MAKE_CSPEC_FILE
      PUBLIC  :: READ_CSPEC_FILE
      ! adj_group (dkh, ks, mak, cs, 06/08/09) 
      PUBLIC  :: CHECK_DIMENSIONS



      !=================================================================
      ! MODULE VARIABLES
      !=================================================================    
      CHARACTER(LEN=255) :: INPUT_RESTART_FILE  
      CHARACTER(LEN=255) :: OUTPUT_RESTART_FILE 

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_RESTART_FILE( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine MAKE_RESTART_FILE creates GEOS-CHEM restart files of tracer 
!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!
!  NOTES:
!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
!        (bmy, 6/22/00)
!  (3 ) Now do not write more than NTRACE data blocks to disk.  
!        Also updated comments. (bmy, 7/17/00)
!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
!        restart file. (bmy, 6/24/02)
!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
!        Now references function GET_TAU from "time_mod.f".  Now added a call 
!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!  (8 ) Cosmetic changes (bmy, 4/29/03)
!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
!        remove hardwired output restart filename.   Now references LPRT
!        from "logical_mod.f". (bmy, 7/20/04)
!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
!        grids. (bmy, 6/28/05)
!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (12) Add TAU to the argument list (bmy, 12/16/05)
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE DAO_MOD,     ONLY : AD
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
      TYPE (XPLEX),  INTENT(IN)  :: TAU

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1
      
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      !=================================================================
      ! MAKE_RESTART_FILE begins here!
      !=================================================================

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Restart File: ' // 
     &           'Instantaneous Tracer Concentrations (v/v)'
      UNIT     = 'v/v'
      CATEGORY = 'IJ-AVG-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Call GET_HALFPOLAR to return the proper value
      ! for either GCAP or GEOS grids (bmy, 6/28/05)
      HALFPOLAR = GET_HALFPOLAR()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the restart file for output -- binary punch format
      !=================================================================

      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( OUTPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_RESTART_FILE: Writing ', a )

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each tracer to the restart file
      !=================================================================
      DO N = 1, N_TRACERS
         
         ! Convert from [kg] to [v/v] and store in the TRACER array
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TRACER(I,J,L) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Convert STT from [kg] to [v/v] mixing ratio 
         ! and store in temporary variable TRACER
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      TAU,       TAU,       RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
     &               J0+1,      1,         TRACER )
      ENDDO  

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_RESTART_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_RESTART_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_RESTART_FILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_RESTART_FILE initializes GEOS-CHEM tracer concentrations 
!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  NOTES:
!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
!        Also reorganize some print statements  (bmy, 10/25/99)
!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!  (5 ) Broke up sections of code into internal subroutines.  Also updated
!        comments & cleaned up a few things. (bmy, 7/17/00)
!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
!  (10) Added updates from amf for tagged Ox run.  Also updated comments
!        and made some cosmetic changes (bmy, 7/3/01)
!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
!        accordingly. (amf, bmy, 9/6/01)
!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
!  (13) Updated comments (bmy, 1/25/02)
!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
!  (17) Add fancy output string (bmy, 4/26/04)
!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
!        and "tracer_mod.f" (bmy, 7/20/04)
!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
!        (bmy, 6/24/05)
!  (20) Updated comments (bmy, 12/16/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE DAO_MOD,     ONLY : AD
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N
      INTEGER             :: NCOUNT(NNPAR)
      REAL*4              :: D_TRACER(IIPAR,JJPAR,LLPAR) 
      TYPE (XPLEX)        :: TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: D_LONRES, D_LATRES
      REAL*8              :: D_ZTAU0, D_ZTAU1
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_RESTART_FILE begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0d0
      D_TRACER(:,:,:) = 0d0

      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RESTART_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      ! Echo more output
      WRITE( 6, 110 )
 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
     &        /, '(in volume mixing ratio units: v/v)' )

      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO 
         READ( IU_RST, IOSTAT=IOS ) 
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LATRES%r = dble(D_LATRES)
         LONRES%r = dble(D_LONRES)
         !LATRES%i = dimag(D_LATRES)
         !LONRES%i = dimag(D_LONRES) 
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )

         READ( IU_RST, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r= dble(D_ZTAU0)
         ZTAU1%r = dble(D_ZTAU1)
         !ZTAU0%i= dimag(D_ZTAU0)
         !ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')

         READ( IU_RST, IOSTAT=IOS ) 
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         TRACER(:,:,:)%r = dble(D_TRACER(:,:,:))
         !TRACER(:,:,:)%i = dimag(D_TRACER(:,:,:)) 
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
  
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-AVG-$' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

            ! Convert TRACER from its native units to [v/v] mixing ratio
            CALL CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )

            ! Convert TRACER from [v/v] to [kg] and copy into STT array
            CALL COPY_STT( NTRACER, TRACER, NCOUNT )

         ENDIF
      ENDDO

      !=================================================================
      ! Examine data blocks, print totals, and return
      !=================================================================

      ! Check for missing or duplicate data blocks
      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )

      ! Close file
      CLOSE( IU_RST )      

      ! Print totals atmospheric mass for each tracer
      WRITE( 6, 120 )
 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 

      DO N = 1, N_TRACERS

         ! For tracers in kg C, be sure to use correct unit string
         IF ( INT( TRACER_MW_G(N) + 0.5 ) == 12 ) THEN
            UNIT = 'kg C'
         ELSE
            UNIT = 'kg  '
         ENDIF

         ! Print totals
         WRITE( 6, 130 ) N,                   TRACER_NAME(N), 
     &                   SUM( STT(:,:,:,N) ), ADJUSTL( UNIT )
 130     FORMAT( 'Tracer ', i3, ' (', a10, ') ', 2es12.5, 1x, a4)
      ENDDO

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_RESTART_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_RESTART_FILE

!------------------------------------------------------------------------------

      SUBROUTINE CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )
!
!******************************************************************************
!  Subroutine CONVERT_TRACER_TO_VV converts the TRACER array from its
!  natural units (e.g. ppbv, ppmv) as read from the restart file to v/v
!  mixing ratio. (bmy, 6/25/02, 6/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER)   : Tracer number
!  (2 ) TRACER  (TYPE (XPLEX) )   : Array containing tracer concentrations
!  (3 ) UNIT    (CHARACTER) : Unit of tracer as read in from restart file
!
!  NOTES:
!  (1 ) Added to "restart_mod.f".  Can now also convert from ppm or ppmv
!        to v/v mixing ratio. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (3 ) Remove obsolete reference to CMN (bmy, 6/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : TRANUC
      USE ERROR_MOD,   ONLY : GEOS_CHEM_STOP

#     include "CMN_SIZE"              ! Size parameters

      ! Arguments
      INTEGER,          INTENT(IN)    :: NTRACER
      TYPE (XPLEX),           INTENT(INOUT) :: TRACER(IIPAR,JJPAR,LLPAR) 
      CHARACTER(LEN=*), INTENT(IN)    :: UNIT

      !=================================================================
      ! CONVERT_TRACER_TO_VV begins here!
      !=================================================================

      ! Convert UNIT to uppercase
      CALL TRANUC( UNIT )
      
      ! Convert from the current unit to v/v
      SELECT CASE ( TRIM( UNIT ) )

         CASE ( '', 'V/V' )
            ! Do nothing, TRACER is already in v/v
            
         CASE ( 'PPM', 'PPMV', 'PPMC' ) 
            TRACER = TRACER * 1d-6

         CASE ( 'PPB', 'PPBV', 'PPBC' ) 
            TRACER = TRACER * 1d-9

         CASE ( 'PPT', 'PPTV', 'PPTC' )
            TRACER = TRACER * 1d-12

         CASE DEFAULT
            WRITE( 6, '(a)' ) 'Incompatible units in punch file!'
            WRITE( 6, '(a)' ) 'STOP in CONVERT_TRACER_TO_VV'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP

      END SELECT

      ! Print the min & max of each tracer as it is read from the file
      WRITE( 6, 110 ) NTRACER,  MINVAL( TRACER ), MAXVAL( TRACER )
 110  FORMAT( 'Tracer ', i3, ': Min = ', 2es12.5, '  Max = ',  2es12.5 )

      ! Return to READ_RESTART_FILE
      END SUBROUTINE CONVERT_TRACER_TO_VV

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DIMENSIONS( NI, NJ, NL ) 
!
!******************************************************************************
!  Subroutine CHECK_DIMENSIONS makes sure that the dimensions of the
!  restart file extend to cover the entire grid. (bmy, 6/25/02, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NI (INTEGER) : Number of longitudes read from restart file
!  (2 ) NJ (INTEGER) : Number of latitudes  read from restart file
!  (3 ) NL (INTEGER) : Numbef of levels     read from restart file
!
!  NOTES:
!  (1 ) Added to "restart_mod.f".  Now no longer allow initialization with 
!        less than a globally-sized data block. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NI, NJ, NL

#     include "CMN_SIZE"

      !=================================================================
      ! CHECK_DIMENSIONS begins here!
      !=================================================================

      ! Error check longitude dimension: NI must equal IIPAR
      IF ( NI /= IIPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
         WRITE( 6, '(a)' ) 'Wrong number of longitudes encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Error check latitude dimension: NJ must equal JJPAR
      IF ( NJ /= JJPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
         WRITE( 6, '(a)' ) 'Wrong number of latitudes encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
      
      ! Error check vertical dimension: NL must equal LLPAR
      IF ( NL /= LLPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
         WRITE( 6, '(a)' ) 'Wrong number of levels encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_DIMENSIONS

!------------------------------------------------------------------------------
     
      SUBROUTINE COPY_STT( NTRACER, TRACER, NCOUNT )
!
!******************************************************************************
!  Subroutine COPY_STT converts tracer concetrations from [v/v] to [kg] and 
!  then copies the results into the STT tracer array. (bmy, 6/25/02, 6/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Tracer number
!  (2 ) NCOUNT  (INTEGER) : Ctr array - # of data blocks read for each tracer
!  (3 ) TRACER  (TYPE (XPLEX) ) : Tracer concentrations from restart file [v/v]
!
!  NOTES:
!  (1 ) Added to "restart_mod.f".  Also added parallel loops. (bmy, 6/25/02)
!  (2 ) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!  (3 ) Now exit if N is out of range (bmy, 4/29/03)
!  (4 ) Now references N_TRACERS, STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  (5 ) Remove call to TRUE_TRACER_INDEX (bmy, 6/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : AD
      USE TRACER_MOD, ONLY : N_TRACERS, STT, TCVV
      USE ERROR_MOD, ONLY: GEOS_CHEM_STOP 
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      TYPE (XPLEX),  INTENT(IN)    :: TRACER(IIPAR,JJPAR,LLPAR)
      INTEGER, INTENT(INOUT) :: NCOUNT(NNPAR)

      ! Local variables
      INTEGER                :: I, J, L, N
      
      !=================================================================
      ! COPY_STT begins here!
      !=================================================================

      ! Tracer number
      N = NTRACER

      ! Exit if N is out of range
      IF ( N < 1 .or. N > N_TRACERS ) RETURN

      ! Convert from [v/v] to [kg] and store in STT
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT(I,J,L,N) = TRACER(I,J,L) * AD(I,J,L) / TCVV(N) 
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Increment the # of records found for tracer N
      NCOUNT(N) = NCOUNT(N) + 1

      END SUBROUTINE COPY_STT

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DATA_BLOCKS( NTRACE, NCOUNT )
!
!******************************************************************************
!  Subroutine CHECK_DATA_BLOCKS checks to see if we have multiple or 
!  missing data blocks for a given tracer. (bmy, 6/25/02, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACE (INTEGER) : Number of tracers
!  (2 ) NCOUNT (INTEGER) : Ctr array - # of data blocks found per tracer
!
!  NOTES:
!  (1 ) Added to "restart_mod.f".  Also now use F90 intrinsic REPEAT to
!        write a long line of "="'s to the screen. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!******************************************************************************
!      
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP


#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: NTRACE, NCOUNT(NNPAR)
  
      ! Local variables
      INTEGER             :: N

      !=================================================================
      ! CHECK_DATA_BLOCKS begins here! 
      !=================================================================

      ! Loop over all tracers
      DO N = 1, NTRACE

         ! Stop if a tracer has more than one data block 
         IF ( NCOUNT(N) > 1 ) THEN 
            WRITE( 6, 100 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
         
         ! Stop if a tracer has no data blocks 
         IF ( NCOUNT(N) == 0 ) THEN
            WRITE( 6, 110 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDDO

      ! FORMAT statements
 100  FORMAT( 'More than one record found for tracer : ', i4 )
 110  FORMAT( 'No records found for tracer : ',           i4 ) 
 120  FORMAT( 'STOP in CHECK_DATA_BLOCKS (restart_mod.f)'    )

      ! Return to calling program
      END SUBROUTINE CHECK_DATA_BLOCKS

!------------------------------------------------------------------------------

      SUBROUTINE SET_RESTART( INFILE, OUTFILE )
!
!******************************************************************************
!  Subroutine SET_RESTART initializes the variables INPUT_RESTART_FILE and
!  OUTPUT_RESTART_FILE with the values read from the "input.geos" file.
!  (bmy, 7/9/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INFILE  (CHAR*255) : Input restart file name from "input.geos"
!  (2 ) OUTFILE (CHAR*255) : Output restart file name from "input.geos"
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      CHARACTER(LEN=255) :: INFILE, OUTFILE
      
      !=================================================================
      ! SET_RESTART begins here
      !=================================================================
      INPUT_RESTART_FILE  = INFILE
      OUTPUT_RESTART_FILE = OUTFILE
     
      ! Return to calling program
      END SUBROUTINE SET_RESTART

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_CSPEC_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_CSPEC_FILE creates GEOS-CHEM checkpt files of species
!  concentrations. 
!  (dkh, 8/27/04)!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Passed via comode_mod
!  ============================================================================
!  (1 ) CSPEC	: Array of quantities to be checkpointed     
!                      
!
!  NOTES:
! (1) Based on MAKE_RESTART_FILE 
!
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,  ONLY : DEBUG_MSG,   ERROR_STOP
      USE FILE_MOD,   ONLY : IU_RST,      IOERROR
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,   ONLY : EXPAND_DATE, GET_TAU
      USE COMODE_MOD, ONLY : CSPEC,       JLOP
      ! write CSPEC_FULL hotp 2/25/09
      USE COMODE_MOD,  ONLY : CSPEC_FULL
      USE COMODE_MOD,  ONLY : IXSAVE, IYSAVE, IZSAVE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! TAU , NSRCX, LSOILNOX
#     include "comode.h"   ! IGAS

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
      CHARACTER(LEN=255)   :: FILENAME
      CHARACTER(LEN=255)   :: OUTPUT_CHECKPT_FILE


      ! Temporary storage arrays for checkpointed variables
      TYPE (XPLEX)               :: TMP(ILONG, ILAT, IPVERT)

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      ! make HALFPOLAR variable (hotp 2/25/09)
      !INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1

      INTEGER 		   :: MAX_nitr_max
      INTEGER 		   :: NSOFAR
      
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 


      !=================================================================
      ! MAKE_CSPEC_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_CHECKPT_FILE = 'restart.cspec.YYYYMMDDhh'

      ! Clear some arrays 
      ! use minimum value instead of zero hotp 2/25/09
      !TMP(:,:,:)   = 0e0
      TMP(:,:,:)   = 1d-30


      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
     &           'Instantaneous Species Concentrations (#/cm3)'
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      ! get value of HALFPOLAR hotp 2/25/09
      HALFPOLAR = GET_HALFPOLAR()

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CSPEC_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================

      ! Checkpt additional values for full chem simulation
  
      ! Write the final species concetrations after full chemistry
      UNIT = 'molec/cm3/box'
          
      ! replace NCS with 1 for safety hotp 2/25/09
      !DO N = 1 , NTSPEC(NCS)
      DO N = 1 , NTSPEC(1)
 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, IPVERT
         DO J = 1, ILAT
         DO I = 1, ILONG
        
            ! fix recommended by pls, save CSPEC_FULL (hotp 2/18/09)
            ! now set conc no smalle than 1d-30
            IF ( CSPEC_FULL(I,J,L,N) .LE. 1d-30 ) THEN
                 TMP(I,J,L) = 1d-30
            ELSE
                 TMP(I,J,L) = CSPEC_FULL(I,J,L,N)
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               ILONG,     ILAT,      IPVERT,    I0+1,
     &               J0+1,      1,         TMP(:,:,:) )
 
      ENDDO

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_CSPEC_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_CSPEC_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_CSPEC_FILE( YYYYMMDD, HHMMSS, IT_EXISTS ) 
!
!******************************************************************************
!  Subroutine READ_CSPEC_FILE initializes GEOS-CHEM species concentrations 
!  from a checkpoint file (binary punch file format) 
!  (dkh, 8/30/04)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Notes
!  (1 ) Based on READ_RESTART
!  
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD,   ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,    ONLY : IU_RST, IOERROR
      USE FILE_MOD,    ONLY : FILE_EXISTS
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE COMODE_MOD,  ONLY : CSPEC_FULL, JLOP


#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"	   ! LPRT, NSRCX, LSOILNOX
#     include "comode.h"   ! ITLOOP, IGAS

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, JLOOP, NN, NTL
      INTEGER             :: NCOUNT(NNPAR) 
      TYPE (XPLEX)		  :: TMP(ILONG,ILAT,IPVERT)
      LOGICAL             :: IT_EXISTS 

      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: INPUT_CHECKPT_FILE
     

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_CSPEC_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_CHECKPT_FILE = 'restart.cspec.YYYYMMDDhh'

      ! Initialize some variables
      !TMP(:,:,:) = 0e0
      ! use 1e-30 as min (hotp 2/25/09)
      TMP(:,:,:) = 1d-30

      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CSPEC_FILE: Reading ', a )
 
      ! Check to see if cspec restart file exists
      IT_EXISTS = FILE_EXISTS( FILENAME )
      IF ( .not. IT_EXISTS ) THEN 
         RETURN
      ENDIF 

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      ! force NCS to one (hotp 2/25/09)
      !DO N = 1, NTSPEC(NCS)
      DO N = 1, NTSPEC(1)

         ! Read the values of CSPEC
         READ( IU_RST, IOSTAT=IOS )
     &       MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_cspec_file:13' )

         READ( IU_RST, IOSTAT=IOS )
     &         CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &         NTL,      NN,       NL,   IFIRST, JFIRST, LFIRST,
     &         NSKIP

         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_cspec_file:14' )

         READ( IU_RST, IOSTAT=IOS )
     &       ( ( ( TMP(I,J,L), I= 1, NTL), J=1,NN ), L = 1, NL)

         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_cspec_file:16' )

         !==============================================================
         ! Assign data from the TMP array to CSPEC
         !==============================================================

         ! Only process checkpoint data 
         IF ( CATEGORY(1:8) == 'IJ-CHK-$' .and.
     &        NTL           == ILONG      .and. 
     &        NN            == ILAT       .and. 
     &        NL            == IPVERT            ) THEN

            CSPEC_FULL(:,:,:,N) = TMP(:,:,:)


         ELSE
            CALL ERROR_STOP(' Restart data is not correct ', 
     &                   ' reading CSPEC, restart_mod')

         ENDIF
 
      ENDDO ! N 

 555  CONTINUE

      ! Close file
      CLOSE( IU_RST )      


      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_CSPEC_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_CSPEC_FILE

!----------------------------------------------------------------------
      ! End of module
      END MODULE RESTART_MOD
