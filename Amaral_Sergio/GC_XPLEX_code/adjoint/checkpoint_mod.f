! $Id: checkpoint_mod.f,v 1.6 2012/03/01 22:00:26 daven Exp $
      MODULE CHECKPOINT_MOD
!
!******************************************************************************
!  Module CHECKPOINT_MOD contains variables and routines which are used to read
!  and write GEOS-CHEM restart files, which contain tracer concentrations
!  in [v/v] mixing ratio. (bmy, 6/25/02, 12/16/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) INPUT_CHECKPOINT_FILE   : Full path name of the restart file to be read
!  (2 ) OUTPUT_CHECKPOINT_FILE  : Full path name (w/ tokens!) of output file
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_CHECKPOINT_FILE    : Writes restart file to disk 
!  (2 ) READ_CHECKPOINT_FILE    : Reads restart file from disk 
!  (3 ) CONVERT_TRACER_TO_VV : Converts from [ppbv], [ppmv], etc to [v/v]
!  (4 ) CHECK_DIMENSIONS     : Ensures that restart file contains global data
!  (5 ) COPY_STT             : Converts [v/v] to [kg] and stores in STT
!  (6 ) CHECK_DATA_BLOCKS    : Makes sure we have read in data for each tracer
!  (7 ) SET_CHECKPOINT          : Gets restart filenames from "input_mod.f"
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
!  (7 ) Added routine SET_CHECKPOINT.  Now reference "logical_mod.f" and
!        "tracer_mod.f" (bmy, 7/20/04)
!  (8 ) Removed obsolete routines TRUE_TRACER_INDEX and COPY_DATA_FOR_CO_OH
!        (bmy, 6/28/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now pass TAU via the arg list in MAKE_CHECKPOINT_FILE (bmy, 12/15/05)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================    
      CHARACTER(LEN=255) :: INPUT_CHECKPOINT_FILE   
      CHARACTER(LEN=255) :: OUTPUT_CHECKPOINT_FILE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_CONVECTION_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
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
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
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
      CHARACTER*10         :: SUFFIX1
      CHARACTER*1          :: SUFFIX2(4)
      INTEGER              :: T,MULT,IT,LT
      TYPE (XPLEX),  PARAMETER   :: SMALLNUM = xplex(1d-12,0d0)
      !=================================================================
      ! MAKE_CHECKPOINT_FILE begins here!
      !=================================================================

      WRITE (SUFFIX1,'(I8)')YYYYMMDD 

      T = HHMMSS/100

      DO IT = 1, 4
         LT = T-(T/10)*10
         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
         T = T/10
      END DO

      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
     &     //TRIM('CONV_CHK.')//TRIM(SUFFIX1)//TRIM('.')
     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
     &     //TRIM(SUFFIX2(4))

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
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
      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )

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
      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
      

      ! Return to calling program
      END SUBROUTINE MAKE_CONVECTION_CHKFILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_CONVECTION_CHKFILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
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
      COMPLEX*16 :: D_TRACER(IIPAR,JJPAR,LLPAR) 
      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER*10         :: SUFFIX1
      CHARACTER*1          :: SUFFIX2(4)
      INTEGER              :: T,MULT,IT,LT

      !=================================================================
      ! READ_CHECKPOINT_FILE begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) =0e0
      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================

      WRITE (SUFFIX1,'(I8)')YYYYMMDD 

      T = HHMMSS/100

      DO IT = 1, 4
         LT = T-(T/10)*10
         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
         T = T/10
      END DO

      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
     &     //TRIM('CONV_CHK.')//TRIM(SUFFIX1)//TRIM('.')
     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
     &     //TRIM(SUFFIX2(4))

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )

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
         LATRES%i = dimag(D_LATRES)
         LONRES%i = dimag(D_LONRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )

         READ( IU_RST, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
         ZTAU1%r = dble(D_ZTAU1)
         ZTAU0%i = dimag(D_ZTAU0)
         ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')

         READ( IU_RST, IOSTAT=IOS ) 
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         TRACER%r = dble(D_TRACER)
         TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
  
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-AVG-$' ) THEN 

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

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')

      ! Return to calling program
      END SUBROUTINE READ_CONVECTION_CHKFILE

c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------

      SUBROUTINE READ_CHEMISTRY_CHKFILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_CHEMISTRY_CHKFILE  (ks, ???)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  NOTES:
!  (1 ) Based on READ_RESTART_FILE
!  (2 ) Updated for v8 adjoint (mak, dkh, 06/23/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,        ONLY : OPEN_BPCH2_FOR_READ
      USE DAO_MOD,          ONLY : AD
      USE DIRECTORY_ADJ_MOD,ONLY : ADJTMP_DIR 
      USE ERROR_MOD,        ONLY : DEBUG_MSG
      USE FILE_MOD,         ONLY : IU_RST,      IOERROR
      USE LOGICAL_MOD,      ONLY : LSPLIT,      LPRT
      USE TIME_MOD,         ONLY : EXPAND_DATE
      USE TRACER_MOD,       ONLY : N_TRACERS,   STT
      USE TRACER_MOD,       ONLY : TRACER_NAME, TRACER_MW_G

#     include "CMN_SIZE"         ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)       :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER                   :: I, IOS, J, L, N
      INTEGER                   :: NCOUNT(NNPAR) 
      TYPE (XPLEX)                    :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16                    :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                    :: SUMTC
      CHARACTER(LEN=255)        :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER                   :: NI,     NJ,     NL
      INTEGER                   :: IFIRST, JFIRST, LFIRST
      INTEGER                   :: NTRACER,   NSKIP
      INTEGER                   :: HALFPOLAR, CENTER180
      TYPE (XPLEX)                    :: LONRES,    LATRES
      TYPE (XPLEX)                    :: ZTAU0,     ZTAU1
      COMPLEX*16                    :: D_LONRES,    D_LATRES
      COMPLEX*16                    :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)         :: MODELNAME
      CHARACTER(LEN=40)         :: CATEGORY
      CHARACTER(LEN=40)         :: UNIT     
      CHARACTER(LEN=40)         :: RESERVED
      CHARACTER*10              :: SUFFIX1
      CHARACTER*1               :: SUFFIX2(4)
      INTEGER                   :: T,MULT,IT,LT

      !=================================================================
      ! READ_CHEMISTRY_CHKFILE begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) = 0e0
      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================

      ! Use EXPAND_DATE instead of this (dkh, 06/23/09) 
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('CHEM_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))

      INPUT_CHECKPOINT_FILE = TRIM('CHECK_CHK.YYYYMMDD.hhmm')

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      FILENAME =  TRIM( ADJTMP_DIR ) 
     &         // TRIM( FILENAME ) 


      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C H E C K P O I N T   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_CHEMISTRY_CHKFILE: Reading ', a )

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
         LATRES%i = dimag(D_LATRES)
         LONRES%i = dimag(D_LONRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_chemistry_chk:4')

         READ( IU_RST, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU1%r = dble(D_ZTAU1)
         ZTAU0%r = dble(D_ZTAU0)
         ZTAU1%i = dimag(D_ZTAU1)
         ZTAU0%i = dimag(D_ZTAU0)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'readchemistry_chk:5')

         READ( IU_RST, IOSTAT=IOS ) 
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         TRACER%r = dble(D_TRACER)
         TRACER%i = dimag(D_TRACER)
         !-------------------------------------------
         !     *****TESTING CHECKPOINTING*****
         !-------------------------------------------
         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
         
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'readchemistry_chk:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================

         CALL COPY_STT( NTRACER, TRACER, NCOUNT )

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

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      !### Debug
      IF ( LPRT ) 
     &   CALL DEBUG_MSG('### READ_CHEMISTRY_CHKFILE: read file')

      ! Return to calling program
      END SUBROUTINE READ_CHEMISTRY_CHKFILE

c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE_CSP1( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CSP,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$      USE COMODE_MOD,   ONLY : CSPEC
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$#     include "comode.h"
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      INTEGER              :: JLOOP,JJ, KK
c$$$      TYPE (XPLEX)              :: TRACER(ITLOOP,IGAS)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_CSP1.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      !PRINT*,'ITLOOP, IGAS = ',ITLOOP,IGAS
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, IGAS
c$$$      DO I = 1, ITLOOP
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         TRACER(I,J) = CSPEC(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      !-------------------------------------------
c$$$      !     *****TESTING CHECKPOINTING*****
c$$$      !-------------------------------------------
c$$$
c$$$      CALL BPCH2_CSP( IU_RST, ITLOOP, IGAS, TRACER )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE_CSP1
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_CHEMISTRY_CHKFILE_CSP1( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$      USE COMODE_MOD,   ONLY : CSPEC,  JLOP
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$#     include "comode.h"
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(ITLOOP,IGAS)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_CSP1.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS )
c$$$     &        NI,       NJ
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( TRACER(I,J), I=1,ITLOOP ), J=1,IGAS )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, IGAS
c$$$      DO I = 1, ITLOOP
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         CSPEC(I,J) = TRACER(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO         
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_CHEMISTRY_CHKFILE_CSP1
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE_CSP2( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CSP,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$      USE COMODE_MOD,   ONLY : CSPEC
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$#     include "comode.h"
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      INTEGER              :: JLOOP,JJ, KK
c$$$      TYPE (XPLEX)              :: TRACER(ITLOOP,IGAS)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_CSP2.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      !PRINT*,'ITLOOP, IGAS = ',ITLOOP,IGAS
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, IGAS
c$$$      DO I = 1, ITLOOP
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         TRACER(I,J) = CSPEC(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      !-------------------------------------------
c$$$      !     *****TESTING CHECKPOINTING*****
c$$$      !-------------------------------------------
c$$$
c$$$      CALL BPCH2_CSP( IU_RST, ITLOOP, IGAS, TRACER )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE_CSP2
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_CHEMISTRY_CHKFILE_CSP2( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$      USE COMODE_MOD,   ONLY : CSPEC,  JLOP
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$#     include "comode.h"
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(ITLOOP,IGAS)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_CSP2.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS )
c$$$     &        NI,       NJ
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( TRACER(I,J), I=1,ITLOOP ), J=1,IGAS )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, IGAS
c$$$      DO I = 1, ITLOOP
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         CSPEC(I,J) = TRACER(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO         
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_CHEMISTRY_CHKFILE_CSP2
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine CONVERT_TRACER_TO_VV converts the TRACER array from its
c$$$!  natural units (e.g. ppbv, ppmv) as read from the restart file to v/v
c$$$!  mixing ratio. (bmy, 6/25/02, 6/24/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) NTRACER (INTEGER)   : Tracer number
c$$$!  (2 ) TRACER  (TYPE (XPLEX) )   : Array containing tracer concentrations
c$$$!  (3 ) UNIT    (CHARACTER) : Unit of tracer as read in from restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Added to "restart_mod.f".  Can now also convert from ppm or ppmv
c$$$!        to v/v mixing ratio. (bmy, 6/25/02)
c$$$!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
c$$$!        allocated memory before stopping the run. (bmy, 10/15/02)
c$$$!  (3 ) Remove obsolete reference to CMN (bmy, 6/24/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE CHARPAK_MOD, ONLY : TRANUC
c$$$      USE ERROR_MOD,   ONLY : GEOS_CHEM_STOP
c$$$
c$$$#     include "CMN_SIZE"              ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER,          INTENT(IN)    :: NTRACER
c$$$      TYPE (XPLEX),           INTENT(INOUT) :: TRACER(IIPAR,JJPAR,LLPAR) 
c$$$      CHARACTER(LEN=*), INTENT(IN)    :: UNIT
c$$$
c$$$      !=================================================================
c$$$      ! CONVERT_TRACER_TO_VV begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Convert UNIT to uppercase
c$$$      CALL TRANUC( UNIT )
c$$$
c$$$      ! Convert from the current unit to v/v
c$$$      SELECT CASE ( TRIM( UNIT ) )
c$$$
c$$$         CASE ( '', 'V/V' )
c$$$            ! Do nothing, TRACER is already in v/v
c$$$            
c$$$         CASE ( 'PPM', 'PPMV', 'PPMC' ) 
c$$$            TRACER = TRACER * 1d-6
c$$$
c$$$         CASE ( 'PPB', 'PPBV', 'PPBC' ) 
c$$$            TRACER = TRACER * 1d-9
c$$$
c$$$         CASE ( 'PPT', 'PPTV', 'PPTC' )
c$$$            TRACER = TRACER * 1d-12
c$$$
c$$$         CASE DEFAULT
c$$$            WRITE( 6, '(a)' ) 'Incompatible units in punch file!'
c$$$            WRITE( 6, '(a)' ) 'STOP in CONVERT_TRACER_TO_VV'
c$$$            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$            CALL GEOS_CHEM_STOP
c$$$
c$$$      END SELECT
c$$$
c$$$      ! Print the min & max of each tracer as it is read from the file
c$$$      WRITE( 6, 110 ) NTRACER,  MINVAL( TRACER ), MAXVAL( TRACER )
c$$$ 110  FORMAT( 'Tracer ', i3, ': Min = ', es12.5, '  Max = ',  es12.5 )
c$$$
c$$$      ! Return to READ_CHECKPOINT_FILE
c$$$      END SUBROUTINE CONVERT_TRACER_TO_VV
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE CHECK_DIMENSIONS( NI, NJ, NL ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine CHECK_DIMENSIONS makes sure that the dimensions of the
c$$$!  restart file extend to cover the entire grid. (bmy, 6/25/02, 10/15/02)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) NI (INTEGER) : Number of longitudes read from restart file
c$$$!  (2 ) NJ (INTEGER) : Number of latitudes  read from restart file
c$$$!  (3 ) NL (INTEGER) : Numbef of levels     read from restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Added to "restart_mod.f".  Now no longer allow initialization with 
c$$$!        less than a globally-sized data block. (bmy, 6/25/02)
c$$$!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
c$$$!        allocated memory before stopping the run. (bmy, 10/15/02)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: NI, NJ, NL
c$$$
c$$$#     include "CMN_SIZE"
c$$$
c$$$      !=================================================================
c$$$      ! CHECK_DIMENSIONS begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Error check longitude dimension: NI must equal IIPAR
c$$$      IF ( NI /= IIPAR ) THEN
c$$$         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
c$$$         WRITE( 6, '(a)' ) 'Wrong number of longitudes encountered!'
c$$$         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
c$$$         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$         CALL GEOS_CHEM_STOP
c$$$      ENDIF
c$$$
c$$$      ! Error check latitude dimension: NJ must equal JJPAR
c$$$      IF ( NJ /= JJPAR ) THEN
c$$$         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
c$$$         WRITE( 6, '(a)' ) 'Wrong number of latitudes encountered!'
c$$$         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
c$$$         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$         CALL GEOS_CHEM_STOP
c$$$      ENDIF
c$$$      
c$$$      ! Error check vertical dimension: NL must equal LLPAR
c$$$      IF ( NL /= LLPAR ) THEN
c$$$         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
c$$$         WRITE( 6, '(a)' ) 'Wrong number of levels encountered!'
c$$$         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
c$$$         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$         CALL GEOS_CHEM_STOP
c$$$      ENDIF
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE CHECK_DIMENSIONS
c$$$
!------------------------------------------------------------------------------
     
      SUBROUTINE COPY_STT( NTRACER, TRACER, NCOUNT )
!
!******************************************************************************
!  Subroutine COPY_STT copies the results into the STT tracer array. 
!  (Kumaresh, 01/24/08)
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

      ! store Tracers into GEOS-CHEM tracer arry
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
     
      SUBROUTINE COPY_STT_TMP( NTRACER, TRACER, NCOUNT )
!
!******************************************************************************
!  Subroutine COPY_STT copies the results into the STT tracer array. 
!  (Kumaresh, 01/24/08)
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
      USE TRACER_MOD, ONLY : N_TRACERS, STT_TMP, TCVV
      
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

      ! store Tracers into GEOS-CHEM tracer arry
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT_TMP(I,J,L,N) = TRACER(I,J,L) 
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Increment the # of records found for tracer N
      NCOUNT(N) = NCOUNT(N) + 1

      END SUBROUTINE COPY_STT_TMP

c$$$!------------------------------------------------------------------------------
c$$$     
c$$$      SUBROUTINE COPY_STT_ADJ( NTRACER, TRACER, NCOUNT )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine COPY_STT_ADJ converts tracer concetrations copies the results into 
c$$$!  the STT_ADJ tracer array. (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) NTRACER (INTEGER) : Tracer number
c$$$!  (2 ) NCOUNT  (INTEGER) : Ctr array - # of data blocks read for each tracer
c$$$!  (3 ) TRACER  (TYPE (XPLEX) ) : Tracer concentrations from restart file [v/v]
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Added to "restart_mod.f".  Also added parallel loops. (bmy, 6/25/02)
c$$$!  (2 ) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (3 ) Now exit if N is out of range (bmy, 4/29/03)
c$$$!  (4 ) Now references N_TRACERS, STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (5 ) Remove call to TRUE_TRACER_INDEX (bmy, 6/24/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE DAO_MOD,    ONLY : AD
c$$$      USE TRACER_MOD, ONLY : N_TRACERS, STT_ADJ, TCVV
c$$$      
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)    :: NTRACER
c$$$      TYPE (XPLEX),  INTENT(IN)    :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      INTEGER, INTENT(INOUT) :: NCOUNT(NNPAR)
c$$$
c$$$      ! Local variables
c$$$      INTEGER                :: I, J, L, N
c$$$      
c$$$      !=================================================================
c$$$      ! COPY_STT_ADJ begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Tracer number
c$$$      N = NTRACER
c$$$
c$$$      ! Exit if N is out of range
c$$$      IF ( N < 1 .or. N > N_TRACERS ) RETURN
c$$$
c$$$      ! Store Tracer values in GEOS-CHEM tracers
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$      DO L = 1, LLPAR
c$$$      DO J = 1, JJPAR
c$$$      DO I = 1, IIPAR
c$$$         STT_ADJ(I,J,L,N) = TRACER(I,J,L) !* AD(I,J,L) / TCVV(N) 
c$$$      ENDDO
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      ! Increment the # of records found for tracer N
c$$$      NCOUNT(N) = NCOUNT(N) + 1
c$$$
c$$$      END SUBROUTINE COPY_STT_ADJ
c$$$
c$$$!------------------------------------------------------------------------------
c$$$     
c$$$      SUBROUTINE COPY_F( NTRACER, TRACER, NCOUNT )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine COPY_STT copies the results into the STT tracer array. 
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) NTRACER (INTEGER) : Tracer number
c$$$!  (2 ) NCOUNT  (INTEGER) : Ctr array - # of data blocks read for each tracer
c$$$!  (3 ) TRACER  (TYPE (XPLEX) ) : Tracer concentrations from restart file [v/v]
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Added to "restart_mod.f".  Also added parallel loops. (bmy, 6/25/02)
c$$$!  (2 ) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (3 ) Now exit if N is out of range (bmy, 4/29/03)
c$$$!  (4 ) Now references N_TRACERS, STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (5 ) Remove call to TRUE_TRACER_INDEX (bmy, 6/24/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE DAO_MOD,    ONLY : AD
c$$$      USE TRACER_MOD, ONLY : N_TRACERS, F, TCVV
c$$$      
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)    :: NTRACER
c$$$      TYPE (XPLEX),  INTENT(IN)    :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      INTEGER, INTENT(INOUT) :: NCOUNT(NNPAR)
c$$$
c$$$      ! Local variables
c$$$      INTEGER                :: I, J, L, N
c$$$      
c$$$      !=================================================================
c$$$      ! COPY_STT begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Tracer number
c$$$      N = NTRACER
c$$$
c$$$      ! Exit if N is out of range
c$$$      IF ( N < 1 .or. N > N_TRACERS ) RETURN
c$$$
c$$$      ! store Tracers into GEOS-CHEM tracer arry
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$      DO L = 1, LLPAR
c$$$      DO J = 1, JJPAR
c$$$      DO I = 1, IIPAR
c$$$         F(I,J,L,N) = TRACER(I,J,L) !* AD(I,J,L) / TCVV(N) 
c$$$      ENDDO
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      ! Increment the # of records found for tracer N
c$$$      NCOUNT(N) = NCOUNT(N) + 1
c$$$
c$$$      END SUBROUTINE COPY_F
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_RRATE_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CSP,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$      USE COMODE_MOD,  ONLY : R_KPP
c$$$      USE gckpp_Global, ONLY : NTT, IND
c$$$      USE gckpp_Parameters
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      INTEGER              :: JLOOP,JJ, KK
c$$$      TYPE (XPLEX)              :: TRACER(NTT,NREACT)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('RRATE_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      !PRINT*,'ITLOOP, IGAS = ',ITLOOP,IGAS
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, NREACT
c$$$      DO I = 1, NTT
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         TRACER(I,J) = R_KPP(I,IND(J))
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      !-------------------------------------------
c$$$      !     *****TESTING CHECKPOINTING*****
c$$$      !-------------------------------------------
c$$$
c$$$      CALL BPCH2_CSP( IU_RST, NTT, NREACT, TRACER )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_RRATE_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_RRATE_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$      USE COMODE_MOD,  ONLY : R_KPP
c$$$      USE gckpp_Global
c$$$      USE gckpp_Parameters
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(NTT,NREACT)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: II, JJ
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('RRATE_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS )
c$$$     &        II,       JJ
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( TRACER(I,J), I=1,NTT ), J=1,NREACT )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, NREACT
c$$$      DO I = 1, NTT
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         R_KPP(I,IND(J)) = TRACER(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO     
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_RRATE_CHKFILE
c$$$
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
      USE LOGICAL_MOD,   ONLY : LLINOZ

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
 120  FORMAT( 'STOP in CHECK_DATA_BLOCKS (checkpoint_mod.f)' )

      ! Return to calling program
      END SUBROUTINE CHECK_DATA_BLOCKS

c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_ADJOINT_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT_ADJ,       N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLOX  = 1d-6
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLNOX = 1d-8
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLCO  = 1d-9
c$$$      
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('ADJ.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_ADJ_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT_ADJ(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$c$$$            IF(N==1)THEN
c$$$c$$$               IF(STT_ADJ(I,J,L,N).ne.0d0) THEN
c$$$c$$$                  IF( STT_ADJ(I,J,L,N) < SMALLNOX )THEN
c$$$c$$$                     TRACER(I,J,L) = 0d0
c$$$c$$$                  ENDIF
c$$$c$$$               ENDIF
c$$$c$$$            ELSEIF(N==2)THEN
c$$$c$$$               IF(STT_ADJ(I,J,L,N).ne.0d0) THEN
c$$$c$$$                  IF( STT_ADJ(I,J,L,N) < SMALLOX )THEN
c$$$c$$$                     TRACER(I,J,L) = 0d0
c$$$c$$$                  ENDIF
c$$$c$$$               ENDIF
c$$$c$$$            ELSE
c$$$c$$$               IF(STT_ADJ(I,J,L,N).ne.0d0) THEN
c$$$c$$$                  IF( STT_ADJ(I,J,L,N) < SMALLCO )THEN
c$$$c$$$                     TRACER(I,J,L) = 0d0
c$$$c$$$                  ENDIF
c$$$c$$$               ENDIF
c$$$c$$$            ENDIF
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_ADJOINT_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P creates GEOS-CHEM restart files of tracers
c$$$!  in binary punch file format for perturbed chemistry concentrations. 
c$$$!  (Kumaresh, 01/24/08)
c$$$
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_P.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_CHEMISTRY_CHKFILE_P( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHEMISTRY_CHKFILE_P initializes GEOS-CHEM tracer concentrations 
c$$$!  from a binary punch file for perturbed chemistry concentrations. 
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_P.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      DO 
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS < 0 is end-of-file, so exit
c$$$         IF ( IOS < 0 ) EXIT
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &        NSKIP
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$         CALL COPY_STT( NTRACER, TRACER, NCOUNT )
c$$$
c$$$      ENDDO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Check for missing or duplicate data blocks
c$$$      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_CHEMISTRY_CHKFILE_P
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P1( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P1 creates GEOS-CHEM restart files of tracers
c$$$!  in binary punch file format for chemistry checkpoints of type1 information. 
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_P1.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P1
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_CHEMISTRY_CHKFILE_P1( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHEMISTRY_CHKFILE_P1 initializes GEOS-CHEM tracer concentrations 
c$$$!  from a binary punch file for chemistry checkpoints of type1 informations. 
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_P1.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      DO 
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS < 0 is end-of-file, so exit
c$$$         IF ( IOS < 0 ) EXIT
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &        NSKIP
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$         CALL COPY_STT( NTRACER, TRACER, NCOUNT )
c$$$
c$$$      ENDDO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Check for missing or duplicate data blocks
c$$$      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_CHEMISTRY_CHKFILE_P1
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_HSAVE_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_HSAVE_CHKFILE creates GEOS-CHEM restart files of KPP chemistry 
c$$$!  step size in binary punch file format. (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE COMODE_MOD,  ONLY : IXSAVE, IYSAVE, IZSAVE, HSAVE_KPP
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$      USE GCKPP_Global    
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS, JJLOOP
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('HSAVE_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      N = 1
c$$$         
c$$$      ! Store KPP Chemistry step size in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( JJLOOP )
c$$$      DO JJLOOP = 1,NTT
c$$$         I = IXSAVE(JJLOOP)
c$$$         J = IYSAVE(JJLOOP)
c$$$         L = IZSAVE(JJLOOP)
c$$$         TRACER(I,J,L) = HSAVE_KPP(I,J,L)
c$$$      END DO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_HSAVE_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_HSAVE_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_HSAVE_CHKFILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a binary punch file with KPP chemistry step size. (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$      USE GCKPP_Global  
c$$$      USE COMODE_MOD,  ONLY : IXSAVE, IYSAVE, IZSAVE, HSAVE_KPP
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR), JJLOOP 
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('HSAVE_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$      
c$$$      READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &     NI,       NJJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &     NSKIP
c$$$      
c$$$      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$      
c$$$      READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJJ ), L=1,NL )
c$$$         
c$$$      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( JJLOOP )
c$$$      DO JJLOOP = 1,NTT
c$$$         I = IXSAVE(JJLOOP)
c$$$         J = IYSAVE(JJLOOP)
c$$$         L = IZSAVE(JJLOOP)
c$$$         HSAVE_KPP(I,J,L) = TRACER(I,J,L) !* TCVV(N) / AD(I,J,L)
c$$$      END DO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_HSAVE_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_PART_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_HSAVE_CHKFILE creates GEOS-CHEM restart files of KPP chemistry 
c$$$!  step size in binary punch file format. (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE COMODE_MOD,  ONLY : IXSAVE, IYSAVE, IZSAVE, PART_CASE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$      USE GCKPP_Global    
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS, JJLOOP
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('PART_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      N = 1
c$$$         
c$$$      ! Store KPP Chemistry step size in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( JJLOOP )
c$$$      DO JJLOOP = 1,NTT
c$$$         I = IXSAVE(JJLOOP)
c$$$         J = IYSAVE(JJLOOP)
c$$$         L = IZSAVE(JJLOOP)
c$$$         TRACER(I,J,L) = PART_CASE(JJLOOP)
c$$$      END DO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_PART_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_PART_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_HSAVE_CHKFILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a binary punch file with KPP chemistry step size. (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$      USE GCKPP_Global  
c$$$      USE COMODE_MOD,  ONLY : IXSAVE, IYSAVE, IZSAVE, PART_CASE
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR), JJLOOP 
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10        :: SUFFIX1
c$$$      CHARACTER*1         :: SUFFIX2(4)
c$$$      INTEGER             :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('PART_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$      
c$$$      READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &     NI,       NJJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &     NSKIP
c$$$      
c$$$      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$      
c$$$      READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJJ ), L=1,NL )
c$$$         
c$$$      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( JJLOOP )
c$$$      DO JJLOOP = 1,NTT
c$$$         I = IXSAVE(JJLOOP)
c$$$         J = IYSAVE(JJLOOP)
c$$$         L = IZSAVE(JJLOOP)
c$$$         PART_CASE(JJLOOP) = TRACER(I,J,L) !* TCVV(N) / AD(I,J,L)
c$$$      END DO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_PART_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P2( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P3 creates GEOS-CHEM restart files of tracers 
c$$$!  in binary punch file format. Used to checkpoint tracers for type2 information
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_P2.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P2
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P3( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P3 creates GEOS-CHEM restart files of tracers 
c$$$!  in binary punch file format. Used to checkpoint tracers for type3 information
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CHEM_CHK_P3.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! store GEOS-CHEM tracer in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CHEMISTRY_CHKFILE_P3
c$$$
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_PRESSURE_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine MAKE_PRESSURE_CHKFILE make pressure checkpoint file. 
!  Originally from v7 adj (ks), updated (dkh, 03/07/10) 
!
!  Based on fwd model code MAKE_RESTART_FILE 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : BPCH2, GET_MODELNAME
      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE DAO_MOD,     ONLY : AD
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE,   GET_TAU
      USE TRACER_MOD,  ONLY : STT, N_TRACERS, TCVV, TMP_PRESS
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR


#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
      TYPE (XPLEX),  INTENT(IN)  :: TAU

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      INTEGER              :: JLOOP,JJ, KK
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,1)
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
      CHARACTER*10         :: SUFFIX1
      CHARACTER*1          :: SUFFIX2(4)
      INTEGER              :: T,MULT,IT,LT
      !=================================================================
      ! MAKE_PRESSURE_CHKFILE begins here!
      !=================================================================

! v7 adj kludge
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('PRESS_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))

      OUTPUT_CHECKPOINT_FILE = TRIM('press.chk.YYYYMMDD.hhmm')

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM CHECKPOINT File: '
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE    
      LATRES   = DJSIZE
      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()
      UNIT     = 'hPa'

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )


      !=================================================================
      ! Open the restart file for output -- binary punch format
      !=================================================================

      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      FILENAME =  TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_PRESSURE_CHKFILE: Writing ', a )

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each tracer to the restart file
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
            ! Compute tracer concentration [molec/cm3/box] by
            ! looping over all species belonging to this tracer
         TRACER(I,J,1) = TMP_PRESS(I,J)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1, 
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,         I0+1,
     &            J0+1,      1,         TRACER )

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT )CALL DEBUG_MSG('### MAKE_PRESSURE_CHKFILE: wrote file')
      

      ! Return to calling program
      END SUBROUTINE MAKE_PRESSURE_CHKFILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_PRESSURE_CHKFILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_PRESSURE_CHKFILE reads PRESS_CHK files.  
!   (ks, 2008; dkh, 03/07/10) 
!
!  Based on READ_RESTART_FILE from fwd model 
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  NOTES:
!  (1 ) Now remove files after they have been read (dkh, 05/02/10) 
!  (2 ) Now delete the press.chk.* files after reading (dkh, 05/02/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DAO_MOD,           ONLY : AD
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE LOGICAL_MOD,       ONLY : LSPLIT,      LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : N_TRACERS, STT, TMP_PRESS
      USE TRACER_MOD,        ONLY : TRACER_NAME, TRACER_MW_G
      USE COMODE_MOD,        ONLY : JLOP
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD


#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N
      INTEGER             :: NCOUNT(NNPAR) 
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,1)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,1)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: REMOVE_CHK_FILE_CMD


      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER*10        :: SUFFIX1
      CHARACTER*1         :: SUFFIX2(4)
      INTEGER             :: T,MULT,IT,LT

      !=================================================================
      ! READ_PRESSURE_CHKFILE begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:)=0e0
      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================

      ! v7 adjoint kludge
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('PRESS_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))

      INPUT_CHECKPOINT_FILE =  TRIM('press.chk.YYYYMMDD.hhmm')

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
!      ! Echo more output
!      WRITE( 6, 110 )
! 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
!     &        /, '(in volume mixing ratio units: v/v)' )

      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      READ( IU_RST, IOSTAT=IOS )
     &  MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
        LONRES%r = dble(D_LONRES)
        LATRES%r = dble(D_LATRES)
        LONRES%i = dimag(D_LONRES)
        LATRES%i = dimag(D_LATRES)
      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 )
     &   CALL IOERROR( IOS,IU_RST,'READ_PRESSURE_CHKFILE:4' )

      READ( IU_RST, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP
        ZTAU0%r = dble(D_ZTAU0)
        ZTAU1%r = dble(D_ZTAU1)
        ZTAU0%i = dimag(D_ZTAU0)
        ZTAU1%i = dimag(D_ZTAU1)
      IF ( IOS /= 0 )
     &   CALL IOERROR( IOS,IU_RST,'READ_PRESSURE_CHKFILE:5')

      READ( IU_RST, IOSTAT=IOS )
     &     ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
       TRACER%r = dble(D_TRACER)
       TRACER%i = dimag(D_TRACER)
      IF ( IOS /= 0 )
     &   CALL IOERROR( IOS,IU_RST,'READ_PRESSURE_CHKFILE:6')

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
            ! Compute tracer concentration [molec/cm3/box] by
            ! looping over all species belonging to this tracer
         TMP_PRESS(I,J) = TRACER(I,J,1)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO         

      !=================================================================
      ! Examine data blocks, print totals, and return
      !=================================================================

      ! Close file
      CLOSE( IU_RST )      


      ! Remove files if L_CHK_DEL = TRUE 
      IF ( LDEL_CHKPT ) THEN

        REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                         TRIM ( FILENAME )

        CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

        WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102    FORMAT( '     - READ_PRESSURE_CHKFILE: Executing: ',a )

      ENDIF 

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG('### READ_PRESSURE_CHKFILE: read file')

      ! Return to calling program
      END SUBROUTINE READ_PRESSURE_CHKFILE

! now obsolete (dkh, 03/07/10) 
!!------------------------------------------------------------------------------
!
!      SUBROUTINE MAKE_FPBL_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!!
!!******************************************************************************
!!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
!!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Date 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!!
!!  NOTES:
!!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
!!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
!!        (bmy, 6/22/00)
!!  (3 ) Now do not write more than NTRACE data blocks to disk.  
!!        Also updated comments. (bmy, 7/17/00)
!!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
!!        restart file. (bmy, 6/24/02)
!!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
!!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
!!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
!!        Now references function GET_TAU from "time_mod.f".  Now added a call 
!!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!!  (8 ) Cosmetic changes (bmy, 4/29/03)
!!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
!!        remove hardwired output restart filename.   Now references LPRT
!!        from "logical_mod.f". (bmy, 7/20/04)
!!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
!!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
!!        grids. (bmy, 6/28/05)
!!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (12) Add TAU to the argument list (bmy, 12/16/05)
!!******************************************************************************
!!     
!      ! References to F90 modules
!      USE BPCH2_MOD,   ONLY : BPCH2_CSP,         GET_MODELNAME
!      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
!      USE DAO_MOD,     ONLY : AD
!      USE ERROR_MOD,   ONLY : DEBUG_MSG
!      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
!      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
!      USE LOGICAL_MOD, ONLY : LPRT
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV, FP
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "comode.h"
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
!      TYPE (XPLEX),  INTENT(IN)  :: TAU
!
!      ! Local Variables      
!      INTEGER              :: I,    I0, IOS, J,  J0, L, N
!      INTEGER              :: YYYY, MM, DD,  HH, SS
!      INTEGER              :: JLOOP,JJ, KK
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR)
!      CHARACTER(LEN=255)   :: FILENAME
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER              :: HALFPOLAR
!      INTEGER, PARAMETER   :: CENTER180 = 1
!      
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!      CHARACTER*10         :: SUFFIX1
!      CHARACTER*1          :: SUFFIX2(4)
!      INTEGER              :: T,MULT,IT,LT
!      !=================================================================
!      ! MAKE_CHECKPOINT_FILE begins here!
!      !=================================================================
!
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('FPBL_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))
!
!      !PRINT*,'ITLOOP, IGAS = ',ITLOOP,IGAS
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
!     &           'Instantaneous Tracer Concentrations (v/v)'
!      !=================================================================
!      ! Open the restart file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output restart file name into a local variable
!      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
!
!      ! Open restart file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      !=================================================================
!      ! Write each tracer to the restart file
!      !=================================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!            ! Compute tracer concentration [molec/cm3/box] by
!            ! looping over all species belonging to this tracer
!         TRACER(I,J) = FP(I,J)
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      !-------------------------------------------
!      !     *****TESTING CHECKPOINTING*****
!      !-------------------------------------------
!
!      CALL BPCH2_CSP( IU_RST, IIPAR, JJPAR, TRACER )
!
!      ! Close file
!      CLOSE( IU_RST )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
!      
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_FPBL_CHKFILE
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_FPBL_CHKFILE( YYYYMMDD, HHMMSS ) 
!!
!!******************************************************************************
!!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
!!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
!!
!!  Arguments as input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Day 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!!
!!  NOTES:
!!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
!!        Also reorganize some print statements  (bmy, 10/25/99)
!!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
!!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
!!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!!  (5 ) Broke up sections of code into internal subroutines.  Also updated
!!        comments & cleaned up a few things. (bmy, 7/17/00)
!!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
!!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
!!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
!!  (10) Added updates from amf for tagged Ox run.  Also updated comments
!!        and made some cosmetic changes (bmy, 7/3/01)
!!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
!!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
!!        accordingly. (amf, bmy, 9/6/01)
!!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
!!  (13) Updated comments (bmy, 1/25/02)
!!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
!!  (17) Add fancy output string (bmy, 4/26/04)
!!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
!!        and "tracer_mod.f" (bmy, 7/20/04)
!!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
!!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
!!        (bmy, 6/24/05)
!!  (20) Updated comments (bmy, 12/16/05)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
!      USE DAO_MOD,     ONLY : AD
!      USE ERROR_MOD,   ONLY : DEBUG_MSG
!      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
!      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!      USE TRACER_MOD,  ONLY : N_TRACERS
!      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G, FP
!      USE COMODE_MOD,  ONLY : JLOP
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "comode.h"
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
!
!      ! Local Variables
!      INTEGER             :: I, IOS, J, L, N
!      INTEGER             :: NCOUNT(NNPAR) 
!      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR)
!      TYPE (XPLEX)              :: SUMTC
!      CHARACTER(LEN=255)  :: FILENAME
!
!      ! For binary punch file, version 2.0
!      INTEGER             :: NI,     NJ,     NL
!      INTEGER             :: IFIRST, JFIRST, LFIRST
!      INTEGER             :: NTRACER,   NSKIP
!      INTEGER             :: HALFPOLAR, CENTER180
!      TYPE (XPLEX)              :: LONRES,    LATRES
!      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
!      CHARACTER(LEN=20)   :: MODELNAME
!      CHARACTER(LEN=40)   :: CATEGORY
!      CHARACTER(LEN=40)   :: UNIT     
!      CHARACTER(LEN=40)   :: RESERVED
!      CHARACTER*10         :: SUFFIX1
!      CHARACTER*1          :: SUFFIX2(4)
!      INTEGER              :: T,MULT,IT,LT
!
!      !=================================================================
!      ! READ_CHECKPOINT_FILE begins here!
!      !=================================================================
!
!      ! Initialize some variables
!      NCOUNT(:)     = 0
!      TRACER(:,:) = 0e0
!
!      !=================================================================
!      ! Open restart file and read top-of-file header
!      !=================================================================
!
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('FPBL_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      ! Echo more output
!      WRITE( 6, 110 )
! 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
!     &        /, '(in volume mixing ratio units: v/v)' )
!
!      !=================================================================
!      ! Read concentrations -- store in the TRACER array
!      !=================================================================
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        NI,       NJ
!
!         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
!
!         READ( IU_RST, IOSTAT=IOS ) 
!     &        ( ( TRACER(I,J), I=1,IIPAR ), J=1,JJPAR )
!
!         !-------------------------------------------
!         !     *****TESTING CHECKPOINTING*****
!         !-------------------------------------------
!         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
!         
!         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!            ! Compute tracer concentration [molec/cm3/box] by
!            ! looping over all species belonging to this tracer
!         FP(I,J) = TRACER(I,J)
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO         
!
!      !=================================================================
!      ! Examine data blocks, print totals, and return
!      !=================================================================
!
!      ! Close file
!      CLOSE( IU_RST )      
!
!      ! Print totals atmospheric mass for each tracer
!      WRITE( 6, 120 )
! 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
!
!      ! Fancy output
!      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
!
!      ! Return to calling program
!      END SUBROUTINE READ_FPBL_CHKFILE
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE MAKE_IMIX_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!!
!!******************************************************************************
!!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
!!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Date 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!!
!!  NOTES:
!!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
!!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
!!        (bmy, 6/22/00)
!!  (3 ) Now do not write more than NTRACE data blocks to disk.  
!!        Also updated comments. (bmy, 7/17/00)
!!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
!!        restart file. (bmy, 6/24/02)
!!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
!!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
!!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
!!        Now references function GET_TAU from "time_mod.f".  Now added a call 
!!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!!  (8 ) Cosmetic changes (bmy, 4/29/03)
!!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
!!        remove hardwired output restart filename.   Now references LPRT
!!        from "logical_mod.f". (bmy, 7/20/04)
!!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
!!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
!!        grids. (bmy, 6/28/05)
!!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (12) Add TAU to the argument list (bmy, 12/16/05)
!!******************************************************************************
!!     
!      ! References to F90 modules
!      USE BPCH2_MOD,   ONLY : BPCH2_INT,         GET_MODELNAME
!      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
!      USE DAO_MOD,     ONLY : AD
!      USE ERROR_MOD,   ONLY : DEBUG_MSG
!      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
!      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
!      USE LOGICAL_MOD, ONLY : LPRT
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV, IM
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "comode.h"
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
!      TYPE (XPLEX),  INTENT(IN)  :: TAU
!
!      ! Local Variables      
!      INTEGER              :: I,    I0, IOS, J,  J0, L, N
!      INTEGER              :: YYYY, MM, DD,  HH, SS
!      INTEGER              :: JLOOP,JJ, KK
!      INTEGER              :: TRACER(IIPAR,JJPAR)
!      CHARACTER(LEN=255)   :: FILENAME
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER              :: HALFPOLAR
!      INTEGER, PARAMETER   :: CENTER180 = 1
!      
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!      CHARACTER*10         :: SUFFIX1
!      CHARACTER*1          :: SUFFIX2(4)
!      INTEGER              :: T,MULT,IT,LT
!      !=================================================================
!      ! MAKE_CHECKPOINT_FILE begins here!
!      !=================================================================
!
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('IMIX_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))
!
!      !PRINT*,'ITLOOP, IGAS = ',ITLOOP,IGAS
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
!     &           'Instantaneous Tracer Concentrations (v/v)'
!      !=================================================================
!      ! Open the restart file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output restart file name into a local variable
!      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
!
!      ! Open restart file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      !=================================================================
!      ! Write each tracer to the restart file
!      !=================================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!            ! Compute tracer concentration [molec/cm3/box] by
!            ! looping over all species belonging to this tracer
!         TRACER(I,J) = IM(I,J)
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      !-------------------------------------------
!      !     *****TESTING CHECKPOINTING*****
!      !-------------------------------------------
!
!      CALL BPCH2_INT( IU_RST, IIPAR, JJPAR, TRACER )
!
!      ! Close file
!      CLOSE( IU_RST )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
!      
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_IMIX_CHKFILE
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_IMIX_CHKFILE( YYYYMMDD, HHMMSS ) 
!!
!!******************************************************************************
!!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
!!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
!!
!!  Arguments as input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Day 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!!
!!  NOTES:
!!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
!!        Also reorganize some print statements  (bmy, 10/25/99)
!!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
!!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
!!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!!  (5 ) Broke up sections of code into internal subroutines.  Also updated
!!        comments & cleaned up a few things. (bmy, 7/17/00)
!!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
!!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
!!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
!!  (10) Added updates from amf for tagged Ox run.  Also updated comments
!!        and made some cosmetic changes (bmy, 7/3/01)
!!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
!!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
!!        accordingly. (amf, bmy, 9/6/01)
!!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
!!  (13) Updated comments (bmy, 1/25/02)
!!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
!!  (17) Add fancy output string (bmy, 4/26/04)
!!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
!!        and "tracer_mod.f" (bmy, 7/20/04)
!!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
!!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
!!        (bmy, 6/24/05)
!!  (20) Updated comments (bmy, 12/16/05)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
!      USE DAO_MOD,     ONLY : AD
!      USE ERROR_MOD,   ONLY : DEBUG_MSG
!      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
!      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!      USE TRACER_MOD,  ONLY : N_TRACERS
!      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G, IM
!      USE COMODE_MOD,  ONLY : JLOP
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "comode.h"
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
!
!      ! Local Variables
!      INTEGER             :: I, IOS, J, L, N
!      INTEGER             :: NCOUNT(NNPAR) 
!      INTEGER             :: TRACER(IIPAR,JJPAR)
!      TYPE (XPLEX)              :: SUMTC
!      CHARACTER(LEN=255)  :: FILENAME
!
!      ! For binary punch file, version 2.0
!      INTEGER             :: NI,     NJ,     NL
!      INTEGER             :: IFIRST, JFIRST, LFIRST
!      INTEGER             :: NTRACER,   NSKIP
!      INTEGER             :: HALFPOLAR, CENTER180
!      TYPE (XPLEX)              :: LONRES,    LATRES
!      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
!      CHARACTER(LEN=20)   :: MODELNAME
!      CHARACTER(LEN=40)   :: CATEGORY
!      CHARACTER(LEN=40)   :: UNIT     
!      CHARACTER(LEN=40)   :: RESERVED
!      CHARACTER*10         :: SUFFIX1
!      CHARACTER*1          :: SUFFIX2(4)
!      INTEGER              :: T,MULT,IT,LT
!
!      !=================================================================
!      ! READ_CHECKPOINT_FILE begins here!
!      !=================================================================
!
!      ! Initialize some variables
!      NCOUNT(:)     = 0
!      TRACER(:,:) = 0
!
!      !=================================================================
!      ! Open restart file and read top-of-file header
!      !=================================================================
!
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('IMIX_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      ! Echo more output
!      WRITE( 6, 110 )
! 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
!     &        /, '(in volume mixing ratio units: v/v)' )
!
!      !=================================================================
!      ! Read concentrations -- store in the TRACER array
!      !=================================================================
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        NI,       NJ
!
!         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
!
!         READ( IU_RST, IOSTAT=IOS ) 
!     &        ( ( TRACER(I,J), I=1,IIPAR ), J=1,JJPAR )
!
!         !-------------------------------------------
!         !     *****TESTING CHECKPOINTING*****
!         !-------------------------------------------
!         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
!         
!         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!            ! Compute tracer concentration [molec/cm3/box] by
!            ! looping over all species belonging to this tracer
!         IM(I,J) = TRACER(I,J)
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO         
!
!      !=================================================================
!      ! Examine data blocks, print totals, and return
!      !=================================================================
!
!      ! Close file
!      CLOSE( IU_RST )      
!
!      ! Print totals atmospheric mass for each tracer
!      WRITE( 6, 120 )
! 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
!
!      ! Fancy output
!      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
!
!      ! Return to calling program
!      END SUBROUTINE READ_IMIX_CHKFILE
!
!------------------------------------------------------------------------------
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_EMISRATE_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CSP,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$      USE COMODE_MOD,  ONLY : EMIS_RATE
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$#     include "comode.h"
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      INTEGER              :: JLOOP,JJ, KK
c$$$      INTEGER, PARAMETER   :: IND = 40
c$$$      TYPE (XPLEX)              :: TRACER(ITLOOP,IND)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('EMISRATE_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      !PRINT*,'ITLOOP, IGAS = ',ITLOOP,IGAS
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_EMISRATE_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, IND
c$$$      DO I = 1, ITLOOP
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         TRACER(I,J) = EMIS_RATE(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      !-------------------------------------------
c$$$      !     *****TESTING CHECKPOINTING*****
c$$$      !-------------------------------------------
c$$$
c$$$      CALL BPCH2_CSP( IU_RST, ITLOOP, IND, TRACER )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_EMISRATE_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_EMISRATE_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_EMISRATE_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$      USE COMODE_MOD,  ONLY : EMIS_RATE,  JLOP
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$#     include "comode.h"
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      INTEGER, PARAMETER  :: IND = 40
c$$$      TYPE (XPLEX)             :: TRACER(ITLOOP,IND)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('EMISRATE_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_EMISRATE_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS )
c$$$     &        NI,       NJ
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( TRACER(I,J), I=1,ITLOOP ), J=1,IND )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J )
c$$$      DO J = 1, IND
c$$$      DO I = 1, ITLOOP
c$$$            ! Compute tracer concentration [molec/cm3/box] by
c$$$            ! looping over all species belonging to this tracer
c$$$         EMIS_RATE(I,J) = TRACER(I,J)
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO         
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_EMISRATE_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_EMISRATE_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_F_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : F,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLNUM = 1d-12
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('F_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = F(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_F_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_F_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   F
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('F_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$      DO 
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS < 0 is end-of-file, so exit
c$$$         IF ( IOS < 0 ) EXIT
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &        NSKIP
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$  
c$$$         ! Only process concentration data (i.e. mixing ratio)
c$$$         IF ( CATEGORY(1:8) == 'IJ-AVG-$' ) THEN 
c$$$
c$$$            ! Convert TRACER from [v/v] to [kg] and copy into STT array
c$$$            CALL COPY_F( NTRACER, TRACER, NCOUNT )
c$$$
c$$$         ENDIF
c$$$      ENDDO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Check for missing or duplicate data blocks
c$$$      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_F_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_EMISDEP_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,       N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLOX  = 1d-6
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLNOX = 1d-8
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLCO  = 1d-9
c$$$      
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('EMISDEP.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_ADJ_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_EMISDEP_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_SRCEMIS_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,       N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLOX  = 1d-6
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLNOX = 1d-8
c$$$      TYPE (XPLEX),  PARAMETER   :: SMALLCO  = 1d-9
c$$$      
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('SRCEMIS.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_ADJ_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_SRCEMIS_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_OBS_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('OBS_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_OBS_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_OBS_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('OBS_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      DO 
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS < 0 is end-of-file, so exit
c$$$         IF ( IOS < 0 ) EXIT
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &        NSKIP
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$         CALL COPY_STT( NTRACER, TRACER, NCOUNT )
c$$$
c$$$      ENDDO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Check for missing or duplicate data blocks
c$$$      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_OBS_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_CURR_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CURR_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_CURR_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_CURR_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('CURR_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      DO 
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS < 0 is end-of-file, so exit
c$$$         IF ( IOS < 0 ) EXIT
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &        NSKIP
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$         CALL COPY_STT( NTRACER, TRACER, NCOUNT )
c$$$
c$$$      ENDDO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Check for missing or duplicate data blocks
c$$$      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_CURR_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_BG_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHECKPOINT_FILE creates GEOS-CHEM restart files of tracer 
c$$$!  mixing ratios (v/v), in binary punch file format. (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
c$$$!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
c$$$!        (bmy, 6/22/00)
c$$$!  (3 ) Now do not write more than NTRACE data blocks to disk.  
c$$$!        Also updated comments. (bmy, 7/17/00)
c$$$!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
c$$$!        restart file. (bmy, 6/24/02)
c$$$!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
c$$$!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
c$$$!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
c$$$!        Now references function GET_TAU from "time_mod.f".  Now added a call 
c$$$!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (8 ) Cosmetic changes (bmy, 4/29/03)
c$$$!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
c$$$!        remove hardwired output restart filename.   Now references LPRT
c$$$!        from "logical_mod.f". (bmy, 7/20/04)
c$$$!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
c$$$!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
c$$$!        grids. (bmy, 6/28/05)
c$$$!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
c$$$!  (12) Add TAU to the argument list (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2_CHK,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('BG_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Convert from [kg] to [v/v] and store in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) !* TCVV(N) / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         ! Convert STT from [kg] to [v/v] mixing ratio 
c$$$         ! and store in temporary variable TRACER
c$$$         CALL BPCH2_CHK( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_BG_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE READ_BG_CHKFILE( YYYYMMDD, HHMMSS ) 
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine READ_CHECKPOINT_FILE initializes GEOS-CHEM tracer concentrations 
c$$$!  from a restart file (binary punch file format) (bmy, 5/27/99, 12/16/05)
c$$$!
c$$$!  Arguments as input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Day 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
c$$$!
c$$$!  NOTES:
c$$$!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
c$$$!        Also reorganize some print statements  (bmy, 10/25/99)
c$$$!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
c$$$!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
c$$$!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
c$$$!        Y2K compliant string for all data sets. (bmy, 6/22/00)
c$$$!  (5 ) Broke up sections of code into internal subroutines.  Also updated
c$$$!        comments & cleaned up a few things. (bmy, 7/17/00)
c$$$!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
c$$$!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
c$$$!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
c$$$!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
c$$$!  (10) Added updates from amf for tagged Ox run.  Also updated comments
c$$$!        and made some cosmetic changes (bmy, 7/3/01)
c$$$!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
c$$$!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
c$$$!        accordingly. (amf, bmy, 9/6/01)
c$$$!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
c$$$!  (13) Updated comments (bmy, 1/25/02)
c$$$!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
c$$$!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
c$$$!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
c$$$!  (17) Add fancy output string (bmy, 4/26/04)
c$$$!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
c$$$!        and "tracer_mod.f" (bmy, 7/20/04)
c$$$!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
c$$$!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to A10.
c$$$!        (bmy, 6/24/05)
c$$$!  (20) Updated comments (bmy, 12/16/05)
c$$$!******************************************************************************
c$$$!
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,      IOERROR
c$$$      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : N_TRACERS,   STT,  TCVV
c$$$      USE TRACER_MOD,  ONLY : TRACER_NAME, TRACER_MW_G
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
c$$$
c$$$      ! Local Variables
c$$$      INTEGER             :: I, IOS, J, L, N
c$$$      INTEGER             :: NCOUNT(NNPAR) 
c$$$      TYPE (XPLEX)             :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      TYPE (XPLEX)              :: SUMTC
c$$$      CHARACTER(LEN=255)  :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      INTEGER             :: NI,     NJ,     NL
c$$$      INTEGER             :: IFIRST, JFIRST, LFIRST
c$$$      INTEGER             :: NTRACER,   NSKIP
c$$$      INTEGER             :: HALFPOLAR, CENTER180
c$$$      TYPE (XPLEX)              :: LONRES,    LATRES
c$$$      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
c$$$      CHARACTER(LEN=20)   :: MODELNAME
c$$$      CHARACTER(LEN=40)   :: CATEGORY
c$$$      CHARACTER(LEN=40)   :: UNIT     
c$$$      CHARACTER(LEN=40)   :: RESERVED
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$
c$$$      !=================================================================
c$$$      ! READ_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      ! Initialize some variables
c$$$      NCOUNT(:)     = 0
c$$$      TRACER(:,:,:) = 0e0
c$$$
c$$$      !=================================================================
c$$$      ! Open restart file and read top-of-file header
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
c$$$     &     //TRIM('BG_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Copy input file name to a local variable
c$$$      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      ! Echo some input to the screen
c$$$      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
c$$$      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( 'READ_CHECKPOINT_FILE: Reading ', a )
c$$$
c$$$      ! Open the binary punch file for input
c$$$      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
c$$$      
c$$$      ! Echo more output
c$$$      WRITE( 6, 110 )
c$$$ 110  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
c$$$     &        /, '(in volume mixing ratio units: v/v)' )
c$$$
c$$$      !=================================================================
c$$$      ! Read concentrations -- store in the TRACER array
c$$$      !=================================================================
c$$$
c$$$      DO 
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
c$$$
c$$$         ! IOS < 0 is end-of-file, so exit
c$$$         IF ( IOS < 0 ) EXIT
c$$$
c$$$         ! IOS > 0 is a real I/O error -- print error message
c$$$         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
c$$$     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
c$$$     &        NSKIP
c$$$
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')
c$$$
c$$$         READ( IU_RST, IOSTAT=IOS ) 
c$$$     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
c$$$
c$$$         !-------------------------------------------
c$$$         !     *****TESTING CHECKPOINTING*****
c$$$         !-------------------------------------------
c$$$         !PRINT*,'TRACER(2,2,2)=',TRACER(2,2,2)
c$$$         
c$$$         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')
c$$$
c$$$         !==============================================================
c$$$         ! Assign data from the TRACER array to the STT array.
c$$$         !==============================================================
c$$$
c$$$         CALL COPY_STT( NTRACER, TRACER, NCOUNT )
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            STT(I,J,L,NTRACER) = STT(I,J,L,NTRACER) !* AD(I,J,L) / 
c$$$!     &                           TCVV(NTRACER) 
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$
c$$$      ENDDO
c$$$
c$$$      !=================================================================
c$$$      ! Examine data blocks, print totals, and return
c$$$      !=================================================================
c$$$
c$$$      ! Check for missing or duplicate data blocks
c$$$      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )      
c$$$
c$$$      ! Print totals atmospheric mass for each tracer
c$$$      WRITE( 6, 120 )
c$$$ 120  FORMAT( /, 'Total atmospheric masses for each tracer: ' ) 
c$$$
c$$$      ! Fancy output
c$$$      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### READ_CHECKPOINT_FILE: read file')
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE READ_BG_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
! This still needs updating, plus I don't think it's even used
! (dkh, 06/23/09) 
!      SUBROUTINE MAKE_ORIG_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!!
!!******************************************************************************
!!  Subroutine MAKE_ORIG_CHKFILE ??? ks ???
!!  
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Date 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!!******************************************************************************
!!     
!      ! References to F90 modules
!      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
!      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
!      USE DAO_MOD,     ONLY : AD
!      USE ERROR_MOD,   ONLY : DEBUG_MSG
!      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
!      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
!      USE LOGICAL_MOD, ONLY : LPRT
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
!
!#     include "CMN_SIZE"   ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
!      TYPE (XPLEX),  INTENT(IN)  :: TAU
!
!      ! Local Variables      
!      INTEGER              :: I,    I0, IOS, J,  J0, L, N
!      INTEGER              :: YYYY, MM, DD,  HH, SS
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
!      CHARACTER(LEN=255)   :: FILENAME
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER              :: HALFPOLAR
!      INTEGER, PARAMETER   :: CENTER180 = 1
!      
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!      CHARACTER*10         :: SUFFIX1
!      CHARACTER*1          :: SUFFIX2(4)
!      INTEGER              :: T,MULT,IT,LT
!      !=================================================================
!      ! MAKE_CHECKPOINT_FILE begins here!
!      !=================================================================
!
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      OUTPUT_CHECKPOINT_FILE = TRIM('opt/')
!     &     //TRIM('ORIG_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
!     &           'Instantaneous Tracer Concentrations (v/v)'
!      UNIT     = 'v/v'
!      CATEGORY = 'IJ-AVG-$'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Call GET_HALFPOLAR to return the proper value
!      ! for either GCAP or GEOS grids (bmy, 6/28/05)
!      HALFPOLAR = GET_HALFPOLAR()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the restart file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output restart file name into a local variable
!      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
!
!      ! Open restart file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      !=================================================================
!      ! Write each tracer to the restart file
!      !=================================================================
!
!      DO N = 1, N_TRACERS
!         
!         ! Store GEOS-CHEM tracers in the TRACER array
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            TRACER(I,J,L) = STT(I,J,L,N) * TCVV(N) * 1d9  / AD(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!         
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  N,
!     &               UNIT,      TAU,       TAU,       RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         TRACER )
!      ENDDO  
!
!      ! Close file
!      CLOSE( IU_RST )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
!      
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_ORIG_CHKFILE
!
!!------------------------------------------------------------------------------

c$$$      SUBROUTINE MAKE_PERT_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P3 creates GEOS-CHEM restart files of tracers 
c$$$!  in binary punch file format. Used to checkpoint tracers for type2 information
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('opt/')
c$$$     &     //TRIM('PERT_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) * TCVV(N) * 1d9 / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_PERT_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_OPTZ_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P3 creates GEOS-CHEM restart files of tracers 
c$$$!  in binary punch file format. Used to checkpoint tracers for type2 information
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('opt/')
c$$$     &     //TRIM('OPTZ_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) * TCVV(N) * 1d9 / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_OPTZ_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_DIFFPERT_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P3 creates GEOS-CHEM restart files of tracers 
c$$$!  in binary punch file format. Used to checkpoint tracers for type2 information
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('opt/')
c$$$     &     //TRIM('DIFFPERT_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) * TCVV(N) * 1d9 / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_DIFFPERT_CHKFILE
c$$$
c$$$!------------------------------------------------------------------------------
c$$$
c$$$      SUBROUTINE MAKE_DIFFOPTZ_CHKFILE( YYYYMMDD, HHMMSS, TAU )
c$$$!
c$$$!******************************************************************************
c$$$!  Subroutine MAKE_CHEMISTRY_CHKFILE_P3 creates GEOS-CHEM restart files of tracers 
c$$$!  in binary punch file format. Used to checkpoint tracers for type2 information
c$$$!  (Kumaresh, 01/24/08)
c$$$!
c$$$!  Arguments as Input:
c$$$!  ============================================================================
c$$$!  (1 ) YYYYMMDD : Year-Month-Date 
c$$$!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
c$$$!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
c$$$!******************************************************************************
c$$$!     
c$$$      ! References to F90 modules
c$$$      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
c$$$      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
c$$$      USE DAO_MOD,     ONLY : AD
c$$$      USE ERROR_MOD,   ONLY : DEBUG_MSG
c$$$      USE FILE_MOD,    ONLY : IU_RST,        IOERROR
c$$$      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
c$$$      USE LOGICAL_MOD, ONLY : LPRT
c$$$      USE TIME_MOD,    ONLY : EXPAND_DATE
c$$$      USE TRACER_MOD,  ONLY : STT,           N_TRACERS,  TCVV
c$$$
c$$$#     include "CMN_SIZE"   ! Size parameters
c$$$
c$$$      ! Arguments
c$$$      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
c$$$      TYPE (XPLEX),  INTENT(IN)  :: TAU
c$$$
c$$$      ! Local Variables      
c$$$      INTEGER              :: I,    I0, IOS, J,  J0, L, N
c$$$      INTEGER              :: YYYY, MM, DD,  HH, SS
c$$$      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
c$$$      CHARACTER(LEN=255)   :: FILENAME
c$$$
c$$$      ! For binary punch file, version 2.0
c$$$      TYPE (XPLEX)               :: LONRES, LATRES
c$$$      INTEGER              :: HALFPOLAR
c$$$      INTEGER, PARAMETER   :: CENTER180 = 1
c$$$      
c$$$      CHARACTER(LEN=20)    :: MODELNAME
c$$$      CHARACTER(LEN=40)    :: CATEGORY
c$$$      CHARACTER(LEN=40)    :: UNIT     
c$$$      CHARACTER(LEN=40)    :: RESERVED = ''
c$$$      CHARACTER(LEN=80)    :: TITLE 
c$$$      CHARACTER*10         :: SUFFIX1
c$$$      CHARACTER*1          :: SUFFIX2(4)
c$$$      INTEGER              :: T,MULT,IT,LT
c$$$      !=================================================================
c$$$      ! MAKE_CHECKPOINT_FILE begins here!
c$$$      !=================================================================
c$$$
c$$$      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
c$$$
c$$$      T = HHMMSS/100
c$$$
c$$$      DO IT = 1, 4
c$$$         LT = T-(T/10)*10
c$$$         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
c$$$         T = T/10
c$$$      END DO
c$$$
c$$$      OUTPUT_CHECKPOINT_FILE = TRIM('opt/')
c$$$     &     //TRIM('DIFFOPTZ_CHK.')//TRIM(SUFFIX1)//TRIM('.')
c$$$     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
c$$$     &     //TRIM(SUFFIX2(4))
c$$$
c$$$      ! Define variables for BINARY PUNCH FILE OUTPUT
c$$$      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
c$$$     &           'Instantaneous Tracer Concentrations (v/v)'
c$$$      UNIT     = 'v/v'
c$$$      CATEGORY = 'IJ-AVG-$'
c$$$      LONRES   = DISIZE
c$$$      LATRES   = DJSIZE
c$$$
c$$$      ! Call GET_MODELNAME to return the proper model name for
c$$$      ! the given met data being used (bmy, 6/22/00)
c$$$      MODELNAME = GET_MODELNAME()
c$$$
c$$$      ! Call GET_HALFPOLAR to return the proper value
c$$$      ! for either GCAP or GEOS grids (bmy, 6/28/05)
c$$$      HALFPOLAR = GET_HALFPOLAR()
c$$$
c$$$      ! Get the nested-grid offsets
c$$$      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
c$$$      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
c$$$
c$$$      !=================================================================
c$$$      ! Open the restart file for output -- binary punch format
c$$$      !=================================================================
c$$$
c$$$      ! Copy the output restart file name into a local variable
c$$$      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
c$$$
c$$$      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
c$$$      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
c$$$
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$ 100  FORMAT( '     - MAKE_CHECKPOINT_FILE: Writing ', a )
c$$$
c$$$      ! Open restart file for output
c$$$      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
c$$$
c$$$      !=================================================================
c$$$      ! Write each tracer to the restart file
c$$$      !=================================================================
c$$$
c$$$      DO N = 1, N_TRACERS
c$$$         
c$$$         ! Store GEOS-CHEM tracers in the TRACER array
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, L )
c$$$         DO L = 1, LLPAR
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            TRACER(I,J,L) = STT(I,J,L,N) * TCVV(N) * 1d9 / AD(I,J,L)
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$!$OMP END PARALLEL DO
c$$$         
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
c$$$     &               HALFPOLAR, CENTER180, CATEGORY,  N,
c$$$     &               UNIT,      TAU,       TAU,       RESERVED,   
c$$$     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
c$$$     &               J0+1,      1,         TRACER )
c$$$      ENDDO  
c$$$
c$$$      ! Close file
c$$$      CLOSE( IU_RST )
c$$$
c$$$      !### Debug
c$$$      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_CHECKPOINT_FILE: wrote file')
c$$$      
c$$$
c$$$      ! Return to calling program
c$$$      END SUBROUTINE MAKE_DIFFOPTZ_CHKFILE

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_UPBDFLX_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine MAKE_UPBDFLX_CHKFILE saves STT values for LINOZE adjoint 
!  (ks, dkh, 05/02/10) 
!
!  Based on MAKE_RESTART_FILE (bmy, 5/27/99, 12/16/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!
!  NOTES:
!  ( 1) Add date tokens, clean up (dkh, 05/02/10) 
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : BPCH2,         GET_MODELNAME
      USE BPCH2_MOD,         ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE DAO_MOD,           ONLY : AD
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST,        IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : STT_TMP,        N_TRACERS,  TCVV
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 


#     include "CMN_SIZE"          ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS
      TYPE (XPLEX),  INTENT(IN)        :: TAU

      ! Local Variables      
      INTEGER                    :: I,    I0, IOS, J,  J0, L, N
      INTEGER                    :: YYYY, MM, DD,  HH, SS
      TYPE (XPLEX)                     :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)         :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)                     :: LONRES, LATRES
      INTEGER                    :: HALFPOLAR
      INTEGER, PARAMETER         :: CENTER180 = 1
      
      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT     
      CHARACTER(LEN=40)          :: RESERVED = ''
      CHARACTER(LEN=80)          :: TITLE 
! old code from ks 
!      CHARACTER*10         :: SUFFIX1
!      CHARACTER*1          :: SUFFIX2(4)
!      INTEGER              :: T,MULT,IT,LT
!      TYPE (XPLEX),  PARAMETER   :: SMALLNUM = 1d-12

      !=================================================================
      ! MAKE_UPBDFLX_CHKFILE begins here!
      !=================================================================

! old code from ks 
!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      OUTPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('UPBD_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))
     
      ! now use date tokens (dkh, 05/02/10) 
      OUTPUT_CHECKPOINT_FILE = 'upbd.chk.YYYYMMDD.hhmm'


      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' // 
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
      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      FILENAME = TRIM( ADJTMP_DIR ) // 
     &           TRIM( FILENAME ) 

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_UPBDFLX_CHKFILE: Writing ', a )

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each tracer to the restart file
      !=================================================================
      DO N = 1, 2
         
         ! Convert from [kg] to [v/v] and store in the TRACER array
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TRACER(I,J,L) = STT_TMP(I,J,L,N)
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
      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_UPBDFLX_CHKFILE: wrote file')
      

      ! Return to calling program
      END SUBROUTINE MAKE_UPBDFLX_CHKFILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_UPBDFLX_CHKFILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_UPBDFLX_CHKFILE reads in STT_TMP for LINOZE. 
!  (ks, dkh, 05/02/10) 
!
!  Based on READ_RESTART_FILE (bmy, 5/27/99, 12/16/05)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  NOTES:
!  ( 1) Now use date tokens to make filename (dkh, 05/02/10) 
!  ( 2) Now delete the upbd.chk.* files after reading (dkh, 05/02/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DAO_MOD,           ONLY : AD
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE LOGICAL_MOD,       ONLY : LSPLIT,      LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : N_TRACERS,   STT_TMP
      USE TRACER_MOD,        ONLY : TRACER_NAME, TRACER_MW_G
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD



#     include "CMN_SIZE"          ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER                    :: I, IOS, J, L, N
      INTEGER                    :: NCOUNT(NNPAR) 
      TYPE (XPLEX)                     :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16 :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: SUMTC
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=255)         :: REMOVE_CHK_FILE_CMD


      ! For binary punch file, version 2.0
      INTEGER                    :: NI,     NJ,     NL
      INTEGER                    :: IFIRST, JFIRST, LFIRST
      INTEGER                    :: NTRACER,   NSKIP
      INTEGER                    :: HALFPOLAR, CENTER180
      TYPE (XPLEX)                     :: LONRES,    LATRES
      TYPE (XPLEX)                     :: ZTAU0,     ZTAU1
      COMPLEX*16 :: D_LONRES, D_LATRES, D_ZTAU0,D_ZTAU1
      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT     
      CHARACTER(LEN=40)          :: RESERVED

      !=================================================================
      ! READ_UPBDFLX_CHKFILE begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:)=0e0
      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================

!      WRITE (SUFFIX1,'(I8)')YYYYMMDD 
!
!      T = HHMMSS/100
!
!      DO IT = 1, 4
!         LT = T-(T/10)*10
!         WRITE (SUFFIX2(4-IT+1),'(I1)')LT
!         T = T/10
!      END DO
!
!      INPUT_CHECKPOINT_FILE = TRIM('adjtmp/')
!     &     //TRIM('UPBD_CHK.')//TRIM(SUFFIX1)//TRIM('.')
!     &     //TRIM(SUFFIX2(1))//TRIM(SUFFIX2(2))//TRIM(SUFFIX2(3))
!     &     //TRIM(SUFFIX2(4))

      INPUT_CHECKPOINT_FILE = 'upbd.chk.YYYYMMDD.hhmm'

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_UPBDFLX_CHKFILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO 
         READ( IU_RST, IOSTAT=IOS ) 
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
           LONRES%r = dble(D_LONRES)
        LATRES%r = dble(D_LATRES)
        LONRES%i = dimag(D_LONRES)
        LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'READ_UPBDFLX:4' )

         READ( IU_RST, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
        ZTAU0%r = dble(D_ZTAU0)
        ZTAU1%r = dble(D_ZTAU1)
        ZTAU0%i = dimag(D_ZTAU0)
        ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'READ_UPBDFLX:5')

         READ( IU_RST, IOSTAT=IOS ) 
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
        TRACER%r = dble(D_TRACER)
        TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'READ_UPBDFLX:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
  
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-AVG-$' ) THEN 

            ! Convert TRACER from [v/v] to [kg] and copy into STT array
            CALL COPY_STT_TMP( NTRACER, TRACER, NCOUNT )

         ENDIF
      ENDDO

      !=================================================================
      ! Examine data blocks, print totals, and return
      !=================================================================

      ! Check for missing or duplicate data blocks
      CALL CHECK_DATA_BLOCKS( 2, NCOUNT )

      ! Close file
      CLOSE( IU_RST )      

      ! Remove files if L_CHK_DEL = TRUE 
      IF ( LDEL_CHKPT ) THEN

        REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                         TRIM ( FILENAME )

        CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

        WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102    FORMAT( '     - READ_UPBDFLX_CHKFILE: Executing: ',a )

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG('### READ_UPBDFLX_CHKFILE: read file')

      ! Return to calling program
      END SUBROUTINE READ_UPBDFLX_CHKFILE

!------------------------------------------------------------------------------
      
      SUBROUTINE MAKE_BEFSTRAT_CHKFILE( YYYYMMDD, HHMMSS, TAU )
!     
!******************************************************************************
!  Subroutine MAKE_BEFSTRAT_CHKFILE saves STT values for STRAT_CHEM adjoint 
!  (hml, 07/28/11, adj32_025)
!     
!  Based on MAKE_UPBDFLX_FILE (bmy, 5/27/99, 12/16/05) 
!     
!  Arguments as Input:           
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!     
!  NOTES:
!  ( 1) Add date tokens, clean up (dkh, 05/02/10) 
!******************************************************************************
!     
      ! References to F90 modules 
      USE BPCH2_MOD,         ONLY : BPCH2,         GET_MODELNAME
      USE BPCH2_MOD,         ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE DAO_MOD,           ONLY : AD
      USE ERROR_MOD,         ONLY : DEBUG_MSG      
      USE FILE_MOD,          ONLY : IU_RST,        IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : STT_STRAT_TMP, N_TRACERS,  TCVV
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE TROPOPAUSE_MOD,    ONLY : GET_MIN_TPAUSE_LEVEL
      USE TROPOPAUSE_MOD,    ONLY : ITS_IN_THE_STRAT
      
      ! for new strat chem (hml, 10/07/11)
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACERID_MOD,      ONLY : IDTOX



#     include "CMN_SIZE"          ! Size parameters
      
      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS
      TYPE (XPLEX),  INTENT(IN)        :: TAU
      
      ! Local Variables          
      INTEGER                    :: I,    I0, IOS, J,  J0, L, N
      INTEGER                    :: YYYY, MM, DD,  HH, SS
      INTEGER                    :: LMIN
      TYPE (XPLEX)                     :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)         :: FILENAME
      

      ! For binary punch file, version 2.0
      TYPE (XPLEX)                     :: LONRES, LATRES
      INTEGER                    :: HALFPOLAR 
      INTEGER, PARAMETER         :: CENTER180 = 1
      
      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT
      CHARACTER(LEN=40)          :: RESERVED = ''
      CHARACTER(LEN=80)          :: TITLE
      
      !=================================================================
      ! MAKE_BEFSTRAT_CHKFILE begins here!
      !=================================================================
      
      ! now use date tokens (hml, 07/31/11) 
      OUTPUT_CHECKPOINT_FILE = 'befstrat.chk.YYYYMMDD.hhmm'

      
      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM CHECKPOINT File: ' //
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
      FILENAME = TRIM( OUTPUT_CHECKPOINT_FILE )
      
      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
      
      FILENAME = TRIM( ADJTMP_DIR ) //
     &           TRIM( FILENAME )
      
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_BEFSTRAT_CHKFILE: Writing ', a )
      
      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
      
      !=================================================================
      ! Write each tracer to the restart file
      !=================================================================
      DO N = 1,N_TRACERS

         ! Now use GMI rate for Ox (hml)
!         IF ( ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGOX_SIM() ) .and.
!     &        ( N .eq. IDTOx ) ) CYCLE
         
         ! Get the minimum level extent of the tropopause
         LMIN = GET_MIN_TPAUSE_LEVEL()
         
         ! Convert from [kg] to [v/v] and store in the TRACER array
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = LMIN, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TRACER(I,J,L) = STT_STRAT_TMP(I,J,L,N)
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
      IF ( LPRT ) CALL DEBUG_MSG
     &                 ('### MAKE_BEFSTRAT_CHKFILE: wrote file')

      
      ! Return to calling program
      END SUBROUTINE MAKE_BEFSTRAT_CHKFILE

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_BEFSTRAT_CHKFILE( YYYYMMDD, HHMMSS )

!
!******************************************************************************
!  Subroutine READ_BEFSTRAT_CHKFILE reads in STT_STRAT_TMP for STRAT_CHEM_ADJ. 
!  (hml, 07/28/11, adj32_025) 
!
!  Based on READ_UPDBFLX_FILE (hml, 07/28/11)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  NOTES:
!  ( 1) Now use date tokens to make filename (dkh, 05/02/10) 
!  ( 2) Now delete the upbd.chk.* files after reading (dkh, 05/02/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DAO_MOD,           ONLY : AD
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE LOGICAL_MOD,       ONLY : LSPLIT,      LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : N_TRACERS,   STT_STRAT_TMP
      USE TRACER_MOD,        ONLY : TRACER_NAME, TRACER_MW_G
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD


#     include "CMN_SIZE"          ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER                    :: I, IOS, J, L, N
      INTEGER                    :: NCOUNT(NNPAR)
      TYPE (XPLEX)                     :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16:: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: SUMTC
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=255)         :: REMOVE_CHK_FILE_CMD


      ! For binary punch file, version 2.0
      INTEGER                    :: NI,     NJ,     NL
      INTEGER                    :: IFIRST, JFIRST, LFIRST
      INTEGER                    :: NTRACER,   NSKIP
      INTEGER                    :: HALFPOLAR, CENTER180
      TYPE (XPLEX)                     :: LONRES,    LATRES
      TYPE (XPLEX)                     :: ZTAU0,     ZTAU1
      COMPLEX*16:: D_LONRES,D_LATRES,D_ZTAU0,D_ZTAU1
      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT
      CHARACTER(LEN=40)          :: RESERVED
      !=================================================================
      ! READ_BEFSTRAT_CHKFILE begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) = 0e0
      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================

      INPUT_CHECKPOINT_FILE = 'befstrat.chk.YYYYMMDD.hhmm'

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPOINT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
      
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )
      
      ! Echo some input to the screen
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_BEFSTRAT_CHKFILE: Reading ', a )
      
      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r = dble(D_LONRES)
        LATRES%r = dble(D_LATRES)
        LONRES%i = dimag(D_LONRES)
        LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT     
      
         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'READ_BEFSTRAT:4' )
      
         READ( IU_RST, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
        ZTAU0%r = dble(D_ZTAU0)
        ZTAU1%r = dble(D_ZTAU1)
        ZTAU0%i = dimag(D_ZTAU0)
        ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'READ_BEFSTRAT:5')
      
         READ( IU_RST, IOSTAT=IOS ) 
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
        TRACER%r = dble(D_TRACER)
        TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'READ_BEFSTRAT:6')
      
         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-AVG-$' ) THEN

            ! Convert TRACER from [v/v] to [kg] and copy into STT array
            CALL COPY_STT_STRAT_TMP( NTRACER, TRACER, NCOUNT )

         ENDIF
      ENDDO

      !=================================================================
      ! Examine data blocks, print totals, and return
      !=================================================================

      ! Check for missing or duplicate data blocks
      CALL CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )

      ! Close file
      CLOSE( IU_RST )

      ! Remove files if L_CHK_DEL = TRUE 
      IF ( LDEL_CHKPT ) THEN

        REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                         TRIM ( FILENAME )

        CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

        WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102    FORMAT( '     - READ_BEFSTRAT_CHKFILE: Executing: ',a )

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG('### READ_BEFSTRAT_CHKFILE: read file')

      ! Return to calling program
      END SUBROUTINE READ_BEFSTRAT_CHKFILE

!------------------------------------------------------------------------------

      SUBROUTINE COPY_STT_STRAT_TMP( NTRACER, TRACER, NCOUNT )
!
!******************************************************************************
!  Subroutine COPY_STT copies the results into the STT tracer array. 
!  Based on code by Kumaresh, 01/24/08. 
!  (hml, dkh, 02/14/12, adj32_025) 
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
      USE DAO_MOD,      ONLY : AD
      USE TRACER_MOD,   ONLY : N_TRACERS, STT_STRAT_TMP, TCVV
         
      ! for new strat chem (hml, 10/07/11)
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACERID_MOD, ONLY : IDTOX
      
      
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

      ! Now use GMI rate for Ox (hml)
!      IF ( ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGOX_SIM() ) .and.
!     &     ( N .eq. IDTOx ) ) RETURN

      ! store Tracers into GEOS-CHEM tracer arry
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT_STRAT_TMP(I,J,L,N) = TRACER(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Increment the # of records found for tracer N
      NCOUNT(N) = NCOUNT(N) + 1

      END SUBROUTINE COPY_STT_STRAT_TMP

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CHECKPOINT_MOD
