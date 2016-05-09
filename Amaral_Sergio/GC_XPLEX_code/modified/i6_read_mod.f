! $Id: i6_read_mod.f,v 1.2 2009/06/17 07:39:04 daven Exp $
      MODULE I6_READ_MOD
!
!******************************************************************************
!  Module I6_READ_MOD contains subroutines that unzip, open, and read
!  GEOS-CHEM I-6 (instantaneous 6-hr) met fields from disk. 
!  (bmy, 6/23/03, 10/7/08)
! 
!  Module Routines:
!  =========================================================================
!  (1 ) UNZIP_I6_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) OPEN_I6_FIELDS  : Opens met field files residing in the temp dir
!  (3 ) GET_I6_FIELDS_1 : Wrapper for routine READ_I6
!  (4 ) GET_I6_FIELDS_2 : Wrapper for routine READ_I6
!  (5 ) GET_N_I6        : Returns # of A-3 fields for each DAO data set 
!  (6 ) READ_I6         : Reads A-3 fields from disk
!  (7 ) I6_CHECK        : Checks if we have found all of the fields
! 
!  GEOS-CHEM modules referenced by i6_read_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f     : Module containing NaN and other error check routines
!  (6 ) logical_mod.f   : Module containing GEOS-CHEM logical switches 
!  (7 ) file_mod.f      : Module containing file unit #'s and error checks
!  (8 ) time_mod.f      : Module containing routines for computing time & date
!  (9 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!  (10) tropopause_mod.f: Module containing routines for tropopause (only for 
!                         variable one here)
!  (11) unix_cmds_mod.f : Module containing Unix commands for unzipping etc.
!
!  NOTES:
!  (1 ) Adapted from "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Now use TIMESTAMP_STRING for formatted date/time output (bmy, 10/28/03)
!  (3 ) Now reads either zipped or unzipped files (bmy, 12/11/03)
!  (4 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (5 ) Now references "directory_mod.f", "unix_cmds_mod.f", and
!        "logical_mod.f" (bmy, 7/20/04)
!  (6 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  (7 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (8 ) Now account for GEOS-4 coastal boxes in LWI properly (bmy, 8/10/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now make LWI TYPE (XPLEX) for near-land formulation (ltm, bmy, 5/9/06)
!  (11) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (12) Now set negative SPHU to a very small positive # (bmy, 9/8/06)
!  (13) Now read TROPP files for GEOS-4, and check tropopause level 
!       in case of a variable tropopause (phs, bmy, bdf, 9/14/06)
!  (14) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "i6_read_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: GET_I6_FIELDS_1
      PUBLIC :: GET_I6_FIELDS_2
      PUBLIC :: OPEN_I6_FIELDS 
      PUBLIC :: UNZIP_I6_FIELDS 
      ! adj_group (dkh, ks, mak, cs  06/12/09) 
      PUBLIC :: OPEN_I6_FIELDS_ADJ
      PUBLIC :: GET_I6_FIELDS_1_ADJ

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Number of I6 fields in the file
      INTEGER :: N_I6_FIELDS

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_I6_FIELDS( OPTION, NYMD )
!
!******************************************************************************
!  Subroutine UNZIP_I6_FIELDS invokes a FORTRAN system call to uncompress
!  GEOS-CHEM I-6 met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (bmy, bdf, 6/15/98, 8/4/06)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of A-6 file to be unzipped (optional)
!
!  NOTES:
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR,   GEOS_3_DIR 
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE UNIX_CMDS_MOD, ONLY : BACKGROUND, REDIRECT,   REMOVE_CMD 
      USE UNIX_CMDS_MOD, ONLY : UNZIP_CMD,  WILD_CARD,  ZIP_SUFFIX

#     include "CMN_SIZE"

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD

      ! Local variables
      CHARACTER(LEN=255)            :: GEOS_DIR,   I6_STR
      CHARACTER(LEN=255)            :: I6_FILE_GZ, I6_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_MET_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN

#if   defined( GEOS_3 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         I6_STR   = 'YYYYMMDD.i6.' // GET_RES_EXT() 

#elif defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         I6_STR   = 'YYYYMMDD.i6.' // GET_RES_EXT() 

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         I6_STR   = 'YYYYMMDD.i6.' // GET_RES_EXT() 

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         I6_STR   = 'YYYYMMDD.i6.' // GET_RES_EXT() 

#endif
         
         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
         CALL EXPAND_DATE( I6_STR,   NYMD, 000000 )

         ! Location of zipped A-3 file in data dir
         I6_FILE_GZ  = TRIM( DATA_DIR   ) // TRIM( GEOS_DIR   ) // 
     &                 TRIM( I6_STR     ) // TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir
         I6_FILE     = TRIM( TEMP_DIR   ) // TRIM( I6_STR     )
         
         ! Remove A-3 files for this date from temp dir
         REMOVE_DATE = TRIM( REMOVE_CMD ) // ' '                // 
     &                 TRIM( TEMP_DIR   ) // TRIM( I6_STR     )

         !==============================================================
         ! Define the foreground and background UNZIP commands
         !==============================================================

         ! Foreground unzip
         UNZIP_FG = TRIM( UNZIP_CMD ) // ' ' // TRIM( I6_FILE_GZ ) // 
     &              TRIM( REDIRECT  ) // ' ' // TRIM( I6_FILE    )  

         ! Background unzip
         UNZIP_BG  = TRIM( UNZIP_FG ) // TRIM( BACKGROUND )
      ENDIF

      !=================================================================
      ! Define command to remove all I-6 files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( REMOVE_CMD ) // ' '    // TRIM( TEMP_DIR  ) // 
     &             TRIM( WILD_CARD )  // '.i6.' // TRIM( WILD_CARD ) 

      !=================================================================
      ! Perform an F90 system call to do the desired operation
      !=================================================================
      SELECT CASE ( TRIM( OPTION ) )
         
         ! Unzip A-3 fields in the Unix foreground
         CASE ( 'unzip foreground' )
            WRITE( 6, 100 ) TRIM( I6_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_FG ) )

         ! Unzip A-3 fields in the Unix background
         CASE ( 'unzip background' )
            WRITE( 6, 100 ) TRIM( I6_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_BG ) )

         ! Remove A-3 field for this date in temp dir
         CASE ( 'remove date' )
            WRITE( 6, 110 ) TRIM( I6_FILE )
            CALL SYSTEM( TRIM( REMOVE_DATE ) )
            
         ! Remove all A-3 fields in temp dir
         CASE ( 'remove all' )
            WRITE( 6, 120 ) TRIM( REMOVE_ALL )
            CALL SYSTEM( TRIM( REMOVE_ALL ) )

         ! Error -- bad option!
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid value for OPTION!', 
     &                       'UNZIP_I6_FIELDS (i6_read_mod.f)' )
            
      END SELECT

      ! FORMAT strings
 100  FORMAT( '     - Unzipping: ', a )
 110  FORMAT( '     - Removing: ', a )
 120  FORMAT( '     - About to execute command: ', a )

      ! Return to calling program
      END SUBROUTINE UNZIP_I6_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_I6_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine OPEN_I6_FIELDS opens the I-6 met fields file for date NYMD and 
!  time NHMS. (bmy, bdf, 6/15/98, 10/7/08)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER) : Current value of YYYYMMDD
!  (2 ) NHMS (INTEGER) : Current value of HHMMSS
!
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/13/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (4 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (5 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit IU_I6
!        refers to a valid file on disk (bmy, 3/23/05
!  (6 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (9 ) Updated for variable tropopause (phs, bmy, 9/14/06)
!  (10) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR,   GEOS_3_DIR 
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LUNZIP, LVARTROP
      USE FILE_MOD,      ONLY : IU_I6, IOERROR, FILE_EXISTS, IU_TP
      USE TIME_MOD,      ONLY : EXPAND_DATE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMD, NHMS

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      LOGICAL                :: IT_EXISTS
      INTEGER                :: IOS, IUNIT
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: I6_FILE, TP_FILE
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_I6_FIELDS begins here!
      !=================================================================

      ! Open the A-3 file 0 GMT of each day, or on the first call
      IF ( NHMS == 000000 .or. FIRST ) THEN

#if   defined( GEOS_3 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

#elif defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()
         TP_FILE  = 'YYYYMMDD.tropp.' //  GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )
         CALL EXPAND_DATE( I6_FILE,  NYMD, NHMS )

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( I6_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( I6_FILE )
         ENDIF

         ! Close previously opened A-3 file
         CLOSE( IU_I6 )

         ! Make sure the file unit is valid before we open it 
         IF ( .not. FILE_EXISTS( IU_I6 ) ) THEN 
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_I6_FIELDS (i6_read_mod.f)' )
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_I6,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_I6, 'open_i6_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a )
         
         ! Set the proper first-time-flag false
         FIRST = .FALSE.

#if   !defined( GEOS_3 )

         ! Skip past the ident string
         READ( IU_I6, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_I6, 'open_i6_fields:2' )
         ENDIF
         
         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_I6_FIELDS    

#endif

#if   defined( GEOS_4 ) 

         ! Test if variable tropopause is turned on
         IF ( LVARTROP ) THEN

            !===========================================================
            ! ALSO NEED TO OPEN THE TROPOPAUSE PRESSURE FILE
            !===========================================================

            ! Replace date tokens
            CALL EXPAND_DATE( TP_FILE, NYMD, NHMS )

            ! If unzipping, open GEOS-1 file in TEMP dir
            ! If not unzipping, open GEOS-1 file in DATA dir
            IF ( LUNZIP ) THEN
               PATH = TRIM( TEMP_DIR ) // TRIM( TP_FILE )
            ELSE
               PATH = TRIM( DATA_DIR ) // 
     &              TRIM( GEOS_DIR ) // TRIM( TP_FILE )
            ENDIF

            ! Close previously opened A-3 file
            CLOSE( IU_TP )

            ! Make sure the file unit is valid before we open it 
            IF ( .not. FILE_EXISTS( IU_TP ) ) THEN 
               CALL ERROR_STOP( 'Could not find TROPP file!', 
     &                          'OPEN_I6_FIELDS (i6_read_mod.f)' )
            ENDIF

            ! Open the file
            OPEN( UNIT   = IU_TP,         FILE   = TRIM( PATH ),
     &            STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &            FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
            IF ( IOS /= 0 ) THEN
               CALL IOERROR( IOS, IU_TP, 'open_i6_fields:3' )
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( PATH )
         
            ! Set the proper first-time-flag false
            FIRST = .FALSE.

            ! Skip past the ident string
            READ( IU_TP, IOSTAT=IOS ) IDENT

            IF ( IOS /= 0 ) THEN
               CALL IOERROR( IOS, IU_I6, 'open_i6_fields:4' )
            ENDIF
         ENDIF
#endif

      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_I6_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_I6_FIELDS_ADJ( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine OPEN_I6_FIELDS_ADJ opens the I-6 met fields file for date NYMD and 
!  time NHMS. (bmy, bdf, 6/15/98, 10/7/08)
!
!  Just like OPEN_I6_FIELDS but don't restrict to beginning of day or FIRST
!  (dkh, 9/14/04)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER) : Current value of YYYYMMDD
!  (2 ) NHMS (INTEGER) : Current value of HHMMSS
!
!  NOTES:
!  (1 ) Adapted from OPEN_I6_FIELDS
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR,   GEOS_3_DIR 
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LUNZIP, LVARTROP
      USE FILE_MOD,      ONLY : IU_I6, IOERROR, FILE_EXISTS, IU_TP
      USE TIME_MOD,      ONLY : EXPAND_DATE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMD, NHMS

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      LOGICAL                :: IT_EXISTS
      INTEGER                :: IOS, IUNIT
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: I6_FILE, TP_FILE
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_I6_FIELDS_ADJ begins here!
      !=================================================================

      ! Open the A-3 file 0 GMT of each day, or on the first call
      !IF ( NHMS == 000000 .or. FIRST ) THEN

#if   defined( GEOS_3 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

#elif defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()
         TP_FILE  = 'YYYYMMDD.tropp.' //  GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )
         CALL EXPAND_DATE( I6_FILE,  NYMD, NHMS )

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( I6_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( I6_FILE )
         ENDIF

         ! Close previously opened A-3 file
         CLOSE( IU_I6 )

         ! Make sure the file unit is valid before we open it 
         IF ( .not. FILE_EXISTS( IU_I6 ) ) THEN 
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_I6_FIELDS_ADJ (i6_read_mod.f)' )
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_I6,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_I6, 'open_i6_fields_adj:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a )
         
         ! Set the proper first-time-flag false
         FIRST = .FALSE.

#if   !defined( GEOS_3 )

         ! Skip past the ident string
         READ( IU_I6, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_I6, 'open_i6_fields_adj:2' )
         ENDIF
         
         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_I6_FIELDS    

#endif

#if   defined( GEOS_4 ) 

         ! Test if variable tropopause is turned on
         IF ( LVARTROP ) THEN

            !===========================================================
            ! ALSO NEED TO OPEN THE TROPOPAUSE PRESSURE FILE
            !===========================================================

            ! Replace date tokens
            CALL EXPAND_DATE( TP_FILE, NYMD, NHMS )

            ! If unzipping, open GEOS-1 file in TEMP dir
            ! If not unzipping, open GEOS-1 file in DATA dir
            IF ( LUNZIP ) THEN
               PATH = TRIM( TEMP_DIR ) // TRIM( TP_FILE )
            ELSE
               PATH = TRIM( DATA_DIR ) // 
     &              TRIM( GEOS_DIR ) // TRIM( TP_FILE )
            ENDIF

            ! Close previously opened A-3 file
            CLOSE( IU_TP )

            ! Make sure the file unit is valid before we open it 
            IF ( .not. FILE_EXISTS( IU_TP ) ) THEN 
               CALL ERROR_STOP( 'Could not find TROPP file!', 
     &                          'OPEN_I6_FIELDS_ADJ (i6_read_mod.f)' )
            ENDIF

            ! Open the file
            OPEN( UNIT   = IU_TP,         FILE   = TRIM( PATH ),
     &            STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &            FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
            IF ( IOS /= 0 ) THEN
               CALL IOERROR( IOS, IU_TP, 'open_i6_fields_adj:3' )
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( PATH )
         
            ! Set the proper first-time-flag false
            FIRST = .FALSE.

            ! Skip past the ident string
            READ( IU_TP, IOSTAT=IOS ) IDENT

            IF ( IOS /= 0 ) THEN
               CALL IOERROR( IOS, IU_I6, 'open_i6_fields_adj:4' )
            ENDIF
         ENDIF
#endif

      !ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_I6_FIELDS_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE GET_I6_FIELDS_1( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_I6_FIELDS_1 is a wrapper for routine READ_I6.  This routine
!  calls READ_I6 properly for reading I-6 fields from GEOS-3, GEOS-4, GEOS-5,
!  or GCAP met data sets at the START of a GEOS-CHEM run. 
!  (bmy, 6/23/03, 1/16/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of I-6 fields to be read from disk
!
!  NOTES:
!  (1 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Now also read TO3 and TTO3 for GEOS-5 (bmy, 1/16/07)
!******************************************************************************
!      
      ! References to F90 modules
      USE DAO_MOD, ONLY : ALBD1, LWI,  PS1,   SLP,    SPHU1, T
      USE DAO_MOD, ONLY : TO3,   TTO3, TMPU1, TROPP1, UWND1, VWND1

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

#if   defined( GEOS_3 )

      !=================================================================
      ! GEOS-3: read PS1,   ALBD1, LWI,   SLP,  TROPP, 
      !              UWND1, VWND1, TMPU1, SPHU1    
      !=================================================================
      CALL READ_I6( NYMD=NYMD,    NHMS=NHMS,   
     &              PS=PS1,       ALBD=ALBD1,    LWI=LWI,    
     &              SLP=SLP,      TROPP=TROPP1,  UWND=UWND1,  
     &              VWND=VWND1,   TMPU=TMPU1,    SPHU=SPHU1 ) 

      ! Initialize T with TMPU1
      T = TMPU1

#elif defined( GEOS_4 )

      !=================================================================
      ! GEOS-4: read LWI, PS1, SLP      
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, 
     &              LWI=LWI,   PS=PS1,  SLP=SLP, TROPP=TROPP1 )

#elif defined( GEOS_5 )

      !=================================================================
      ! GEOS-5: read LWI, PS1, SLP, TO3, TTO3
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, 
     &              LWI=LWI,   PS=PS1,  SLP=SLP, TO3=TO3, TTO3=TTO3 )

#elif defined( GCAP )

      !=================================================================
      ! GCAP: read PS1, SLP
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, PS=PS1, SLP=SLP )

#endif
      
      ! Return to calling program
      END SUBROUTINE GET_I6_FIELDS_1

!-----------------------------------------------------------------------------
      SUBROUTINE GET_I6_FIELDS_1_ADJ( NYMD, NHMS )

!
!******************************************************************************
!  Subroutine GET_I6_FIELDS_1_ADJ is a wrapper for routine READ_I6.  This routine
!  calls READ_I6 properly for reading I-6 fields from GEOS-3, GEOS-4, GEOS-5,
!  or GCAP met data sets at the START of a GEOS-CHEM run. 
!  (bmy, 6/23/03, 1/16/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of I-6 fields to be read from disk
!
!  NOTES:
!  (1 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Now also read TO3 and TTO3 for GEOS-5 (bmy, 1/16/07)
!******************************************************************************
!      
      ! References to F90 modules
      USE DAO_MOD, ONLY : ALBD1, LWI,  PS1,   SLP,    SPHU1, T
      USE DAO_MOD, ONLY : TO3,   TTO3, TMPU1, TROPP1, UWND1, VWND1
   
      USE DAO_MOD, ONLY : SLP_TMP, LWI_TMP, TO3_TMP, TTO3_TMP

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

#if   defined( GEOS_3 )

      !=================================================================
      ! GEOS-3: read PS1,   ALBD1, LWI,   SLP,  TROPP, 
      !              UWND1, VWND1, TMPU1, SPHU1    
      !=================================================================
      CALL READ_I6( NYMD=NYMD,    NHMS=NHMS,   
! adj_group: USE SLP_TMP, LWI_TMP
!     &              PS=PS1,       ALBD=ALBD1,    LWI=LWI,    
!     &              SLP=SLP,  TROPP=TROPP1,  UWND=UWND1,  
     &              PS=PS1,       ALBD=ALBD1,    LWI=LWI_TMP,    
     &              SLP=SLP_TMP,  TROPP=TROPP1,  UWND=UWND1,  
     &              VWND=VWND1,   TMPU=TMPU1,    SPHU=SPHU1 ) 

      ! Initialize T with TMPU1
      T = TMPU1

#elif defined( GEOS_4 )

      !=================================================================
      ! GEOS-4: read LWI, PS1, SLP      
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, 
! adj_group: USE SLP_TMP, LWI_TMP
!     &              LWI=LWI,   PS=PS1,  SLP=SLP, TROPP=TROPP1 )
     &              LWI=LWI_TMP,   PS=PS1,  SLP=SLP_TMP, TROPP=TROPP1 )

#elif defined( GEOS_5 )

      !=================================================================
      ! GEOS-5: read LWI, PS1, SLP, TO3, TTO3
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, 
! adj_group: USE LWI_TMP, SLP_TMP, TO3_TMP, TTO3_TMP
!     &              LWI=LWI,   PS=PS1,  SLP=SLP, TO3=TO3, TTO3=TTO3 )
     &              LWI=LWI_TMP,   PS=PS1,  SLP=SLP_TMP, TO3=TO3_TMP, 
     &              TTO3=TTO3_TMP )

#elif defined( GCAP )

      !=================================================================
      ! GCAP: read PS1, SLP
      !=================================================================
! adj_group: USE SLP_TMP
!      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, PS=PS1, SLP=SLP )
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, PS=PS1, SLP=SLP_TMP )

#endif
      
      ! Return to calling program
      END SUBROUTINE GET_I6_FIELDS_1_ADJ

!-----------------------------------------------------------------------------

      SUBROUTINE GET_I6_FIELDS_2( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_I6_FIELDS_2 is a wrapper for routine READ_I6.  This routine
!  calls READ_I6 properly for reading I-6 fields from GEOS-3, GEOS-4, GEOS_5,
!  or GCAP met data sets every 6 hours during a GEOS-CHEM run. 
!  (bmy, 6/23/03, 1/16/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-3 fields to be read from disk
!
!  NOTES:
!  (1 ) Now modified for GEOS-5 and GCAP met fields
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Now read TO3 and TTO3 for GEOS-5 (bmy, 1/16/07)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : ALBD2,  LWI,   PS2,   SLP,    SPHU2
      USE DAO_MOD, ONLY : T,      TO3,   TTO3,  TMPU1,  TMPU2
      USE DAO_MOD, ONLY : TROPP2, UWND2, VWND2

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

#if   defined( GEOS_3 )

      !=================================================================
      ! GEOS-3: read PS2,   ALBD2, LWI,   SLP, TROPP, 
      !              UWND2, VWND2, TMPU2, SPHU2     
      !=================================================================
      CALL READ_I6( NYMD=NYMD,    NHMS=NHMS,   
     &              PS=PS2,       ALBD=ALBD2,    LWI=LWI,    
     &              SLP=SLP,      TROPP=TROPP2,  UWND=UWND2,  
     &              VWND=VWND2,   TMPU=TMPU2,    SPHU=SPHU2 ) 

#elif defined( GEOS_4 )

      !=================================================================
      ! GEOS-4: read LWI, PS2, SLP
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, 
     &              LWI=LWI,   PS=PS2, SLP=SLP, TROPP=TROPP2 )

#elif defined( GEOS_5 )

      !=================================================================
      ! GEOS-5: read LWI, PS2, SLP, TO3, TTO3
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, 
     &              LWI=LWI,   PS=PS2,   SLP=SLP, TO3=TO3, TTO3=TTO3 )


#elif defined( GCAP )

      !=================================================================
      ! GCAP: read PS1, SLP
      !=================================================================
      CALL READ_I6( NYMD=NYMD, NHMS=NHMS, PS=PS2, SLP=SLP )

#endif

      ! Return to calling program
      END SUBROUTINE GET_I6_FIELDS_2

!------------------------------------------------------------------------------

      FUNCTION GET_N_I6( NYMD ) RESULT( N_I6 )
!
!******************************************************************************
!  Function GET_N_I6 returns the number of I-6 fields per met data set
!  (GEOS-3, GEOS-4, GEOS-5, or GCAP). (bmy, 6/23/03, 1/17/06) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD for which to read in I-6 fields
!
!  NOTES:
!  (1 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/25/05)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Increase # of I-6 fields to 5 for GEOS-5 (bmy, 1/17/06)
!******************************************************************************
!
#     include "CMN_SIZE" 

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD

      ! Function value
      INTEGER             :: N_I6

      !=================================================================
      ! GET_N_I6 begins here!
      !=================================================================

#if   defined( GEOS_3 ) 
      
      ! N_I6 is a function of year for GEOS-3
      SELECT CASE ( NYMD / 10000 ) 
         CASE ( 1997, 1998, 1999 )
            N_I6 = 8

         CASE DEFAULT
            N_I6 = 9
            
      END SELECT

#elif defined( GEOS_4 )

      ! GEOS-4 has 3 I-6 fields
      N_I6 = 3 

#elif defined( GEOS_5 )

      ! GEOS-5 has 5 I-6 fields
      N_I6 = 5

#elif defined( GCAP )

      ! GCAP has 2 I-6 fields
      N_I6 = 2 

#endif

      ! Return to calling program
      END FUNCTION GET_N_I6

!---------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the A-3 field just
!  read from disk matches the current time.  If so, then it's time to return
!  the A-3 field to the calling program. (bmy, 6/23/03, 8/4/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) XYMD (INTEGER) : YYYYMMDD timestamp for A-3 field in file
!  (2 ) XHMS (INTEGER) : HHMMSS   timestamp for A-3 field in file
!  (3 ) NYMD (INTEGER) : YYYYMMDD at which A-3 field is to be read
!  (4 ) NHMS (INTEGER) : HHMMSS   at which A-3 field is to be read
!
!  NOTES:
!  (1 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments 
      INTEGER, INTENT(IN) :: XYMD, XHMS, NYMD, NHMS
      
      ! Function value
      LOGICAL             :: ITS_TIME

      !=================================================================
      ! CHECK_TIME begins here!
      !=================================================================
      IF ( XYMD == NYMD .AND. XHMS == NHMS ) THEN
         ITS_TIME = .TRUE.
      ELSE
         ITS_TIME = .FALSE.
      ENDIF

      ! Return to calling program
      END FUNCTION CHECK_TIME

!-----------------------------------------------------------------------------

      SUBROUTINE READ_I6( NYMD, NHMS, 
     &                    ALBD, LWI,   PS,   SLP,  SPHU, TMPU, 
     &                    TO3,  TROPP, TTO3, UWND, VWND         )
!
!******************************************************************************
!  Subroutine READ_I6 reads GEOS-CHEM I-6 (inst. 6-hr) met fields from disk.
!  (bmy, 5/8/98, 10/7/08)
! 
!  Arguments as Input:
!  ===========================================================================
!  (1 ) NYMD   : YYYYMMDD
!  (2 ) NHMS   :  and HHMMSS of I-6 fields to be read from disk
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) ALBD  : (2-D) GMAO Surface albedo            [unitless]
!  (4 ) LWI   : (2-D) GMAO Land-water indices        [unitless]
!  (5 ) PS    : (2-D) GMAO Surface pressure          [hPa]
!  (6 ) SLP   : (2-D) GMAO Sea-level pressures       [hPa]
!  (7 ) SPHU  : (3-D) GMAO Specific humidity field   [g H2O/kg air]
!  (8 ) TMPU  : (3-D) GMAO Temperature field         [K]
!  (9 ) TO3   : (2-D) GMAO GEOS-5 column ozone       [DU]
!  (10) TROPP : (2-D) GMAO tropopause pressure pressures     [hPa]
!  (11) TTO3  : (2-D) GMAO GEOS-5 trop column ozone  [DU]
!  (12) UWND  : (3-D) GMAO U-wind (zonal wind)       [m/s]
!  (13) VWND  : (3-D) GMAO V-wind (meridional wind)  [m/s]
!
!  NOTES:
!  (1 ) Adapted from "READ_I6" of "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!  (3 ) Round up to account for GEOS-4 coastal boxes properly (bmy, 8/10/05)
!  (4 ) For near-land formulation: (a) make LWI a TYPE (XPLEX) and (b) do not round 
!        up LWI for GEOS-4 meteorology (ltm, bmy, 5/9/06)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Now set negative SPHU to a small positive number (1d-32) instead of 
!        zero, so as not to blow up logarithms (bmy, 9/8/06)
!  (7 ) Now read TROPP files for GEOS-4 (phs, bmy, bdf, 9/12/06)
!  (8 ) Now read TO3 and TTO3 for GEOS-5 (bmy, 1/16/07)
!  (9 ) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,       ONLY : AD66,            AD67
      USE FILE_MOD,       ONLY : IOERROR,         IU_I6, IU_TP
      USE LOGICAL_MOD,    ONLY : LVARTROP
      USE TIME_MOD,       ONLY : SET_CT_I6,       TIMESTAMP_STRING
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D,     TRANSFER_3D

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! NDxx flags

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: LWI  (IIPAR,JJPAR      )    
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: ALBD (IIPAR,JJPAR      )  
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: PS   (IIPAR,JJPAR      )     
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: SLP  (IIPAR,JJPAR      )      
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: SPHU (IIPAR,JJPAR,LLPAR) 
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: TMPU (IIPAR,JJPAR,LLPAR)  
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: TO3  (IIPAR,JJPAR      )  
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: TROPP(IIPAR,JJPAR      )  
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: TTO3 (IIPAR,JJPAR      )  
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: UWND (IIPAR,JJPAR,LLPAR)   
      TYPE (XPLEX),  INTENT(OUT), OPTIONAL :: VWND (IIPAR,JJPAR,LLPAR)    

      ! Local Variables
      INTEGER                        :: IOS, NFOUND, N_I6
      REAL*4                         :: D_Q2(IGLOB,JGLOB)
      REAL*4                         :: D_Q3(IGLOB,JGLOB,LGLOB)      
      TYPE (XPLEX)                         :: Q2(IGLOB,JGLOB)
      TYPE (XPLEX)                         :: Q3(IGLOB,JGLOB,LGLOB)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP
      CHARACTER(LEN=255)             :: FILENAME, I6_FILE !for var. tropopause
      INTEGER                        :: XYMD, XHMS

      !=================================================================
      ! READ_I6 begins here!
      !=================================================================

      ! Get the number of I-6 fields
#if   defined( GEOS_5 ) && defined( IN_CLOUD_OD )
      N_I6 = N_I6_FIELDS
#else
      N_I6 = GET_N_I6( NYMD )
#endif

      ! Zero the number of I-6 fields we have already found
      NFOUND = 0

      !=================================================================
      ! Read the I-6 fields from disk
      !=================================================================
      DO 

         ! I-6 field name
         READ( IU_I6, IOSTAT=IOS ) NAME

         ! IOS < 0: End-of-file, but make sure we have 
         ! found all I-6 fields before exiting loop!
         IF ( IOS < 0 ) THEN
            CALL I6_CHECK( NFOUND, N_I6 )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error, stop w/ error msg
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:1' )

         ! CASE statement for met fields
         SELECT CASE ( TRIM( NAME ) )

            !------------------------------------
            ! ALBD: GEOS-3 surface albedo
            !------------------------------------
            CASE ( 'ALBD', 'ALBEDO', 'ALBVISDF' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:2' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ALBD ) ) CALL TRANSFER_2D( Q2, ALBD )
                  NFOUND = NFOUND + 1
               ENDIF
         
            !------------------------------------
            ! LWI: land-water flags
            !------------------------------------
            CASE ( 'LWI', 'SURFTYPE' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:3' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( LWI ) ) CALL TRANSFER_2D( Q2, LWI )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! PS: surface pressure
            !------------------------------------
            CASE ( 'PS' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:4' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PS ) ) CALL TRANSFER_2D( Q2, PS )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! SLP: sea-level pressure 
            !------------------------------------
            CASE ( 'SLP' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:5' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SLP ) ) CALL TRANSFER_2D( Q2, SLP )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! SPHU: GEOS-3 specific humidity
            !------------------------------------
            CASE ( 'SPHU' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:6' )
               Q3(:,:,:) = (D_Q3(:,:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SPHU ) ) CALL TRANSFER_3D( Q3, SPHU )
                  !=====================================================
                  ! NOTE: Now set negative SPHU to a small positive # 
                  ! instead of zero, so as not to blow up logarithms
                  ! (bmy, 9/8/06)
                  WHERE ( SPHU%r < 0d0 ) 
                          SPHU%r = 1d-32
                          SPHU%i = 0d0
                  endwhere
                  !=====================================================
                  NFOUND = NFOUND + 1
               ENDIF

            !---------------------------------
            ! TMPU: GEOS-3 temperature
            !---------------------------------
            CASE ( 'TMPU' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:7' )
               Q3(:,:,:) = (D_Q3(:,:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TMPU ) ) CALL TRANSFER_3D( Q3, TMPU )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------            
            ! TO3: GEOS-5 column ozone
            !------------------------------------
            CASE ( 'TO3' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:8' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TO3 ) ) CALL TRANSFER_2D( Q2, TO3 )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! TROPP: GEOS-3 tropopause pressure
            !------------------------------------
            CASE ( 'TROPP' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:9' )
               Q2(:,:) = (D_Q2(:,:)) 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TROPP ) ) CALL TRANSFER_2D( Q2, TROPP )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------            
            ! TTO3: GEOS-5 trop column ozone
            !------------------------------------
            CASE ( 'TTO3' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:10' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TTO3 ) ) CALL TRANSFER_2D( Q2, TTO3 )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! UWND: GEOS-3 zonal wind field
            !------------------------------------
            CASE ( 'UWND' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:11' )
               Q3(:,:,:) = (D_Q3(:,:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( UWND ) ) CALL TRANSFER_3D( Q3, UWND )
                  NFOUND = NFOUND + 1
               ENDIF
            !UWND(:,:,:)%i = 1d-120*UWND(:,:,:)%r
            !------------------------------------            
            ! VWND: GEOS-3 meridional wind field
            !------------------------------------
            CASE ( 'VWND' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:12' )
               Q3(:,:,:) = (D_Q3(:,:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( VWND ) ) CALL TRANSFER_3D( Q3, VWND )
                  NFOUND = NFOUND + 1
               ENDIF
               
            !------------------------------------
            ! TKE: Just skip over this
            !------------------------------------
            CASE ( 'TKE' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:13' )
               Q3(:,:,:) = (D_Q3(:,:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------  
            ! RH: just skip over this
            !--------------------------------  
            CASE ( 'RH' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, D_Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:14' )
               Q3(:,:,:) = (D_Q3(:,:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! Field not found
            !--------------------------------
            CASE DEFAULT
               WRITE ( 6, '(a)' ) 'Searching for next I-6 field!'
               
         END SELECT

         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop and return to the calling program.  Otherwise, 
         ! go to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .AND. 
     &        NFOUND == N_I6 ) THEN
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, ' I-6 met fields for ', a )
            EXIT
         ENDIF
      ENDDO


      !======== CASE OF VARIABLE TROPOPAUSE ============================
      ! We need to read TROPP from offline files in case of GEOS-4.
      !=================================================================
      IF ( LVARTROP ) THEN

#if   defined( GEOS_4 )

         NFOUND = 0   ! need reset because dealing with new file

         DO 
            ! I-6 field name
            READ( IU_TP, IOSTAT=IOS ) NAME

            ! IOS < 0: End-of-file, but make sure we have 
            ! found all I-6 fields before exiting loop (only 4 here)!
            IF ( IOS < 0 ) THEN
               CALL I6_CHECK( NFOUND, 1 )
               EXIT
            ENDIF

            ! IOS > 0: True I/O error, stop w/ error msg
            IF ( IOS > 0 ) CALL IOERROR( IOS, IU_TP, 'read_i6:15' )

            ! CASE statement for met fields
            SELECT CASE ( TRIM( NAME ) )

            !------------------------------------
            ! TROPP: GEOS-4 tropopause pressure
            !------------------------------------
            CASE ( 'TROPP' ) 
               READ( IU_TP, IOSTAT=IOS ) XYMD, XHMS, D_Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_TP, 'read_i6:16' )
               Q2(:,:) = (D_Q2(:,:))
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TROPP ) ) CALL TRANSFER_2D( Q2, TROPP )
                  NFOUND = NFOUND + 1
               ENDIF
            
            CASE DEFAULT 
               ! Nothing

            END SELECT

            !==========================================================
            ! If we have found the one field for this time, then exit 
            ! the loop. Otherwise, go to the next iteration.
            !==========================================================
            IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .AND. 
     &           NFOUND == 1 ) THEN
               STAMP = TIMESTAMP_STRING( NYMD, NHMS )
               WRITE( 6, 230 ) NFOUND, STAMP
 230           FORMAT( '     - Found all ', i3, 
     &                 ' TROPP met fields for ', a )

               EXIT
            ENDIF
         ENDDO

#endif

!         ! Testing:
!         write( 6, 220 ) MINVAL( TROPP ), MAXVAL( TROPP ) 
! 220     FORMAT( ' CHECKING TROPP: MIN=', f10.3 , ' MAX=', f10.3 )


      ENDIF
      !======== END CASE OF VARIABLE TROPOPAUSE ================



      ! Increment KDI6FLDS -- this is the # of times READ_I6 is called. 
      CALL SET_CT_I6( INCREMENT=.TRUE. )

      !=================================================================
      ! ND66 diagnostic: I-6 fields (3-dimensional)
      !
      ! (1 ) UWND : Instantaneous U-winds             [m/s]
      ! (2 ) VWND : Instantaneous V-winds             [m/s]
      ! (3 ) TMPU : Instantaneous Temperature         [K]
      ! (4 ) SPHU : Instantaneous Specific humidity   [g H20/kg air] 
      !=================================================================
      IF ( ND66 > 0 ) THEN
         IF ( PRESENT( UWND ) ) THEN
            AD66(:,:,1:LD66,1) = AD66(:,:,1:LD66,1) + UWND(:,:,1:LD66)
         ENDIF

         IF ( PRESENT( VWND ) ) THEN
            AD66(:,:,1:LD66,2) = AD66(:,:,1:LD66,2) + VWND(:,:,1:LD66)
         ENDIF

         IF ( PRESENT( TMPU ) ) THEN
            AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + TMPU(:,:,1:LD66)
         ENDIF

         IF ( PRESENT( SPHU ) ) THEN
            AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + SPHU(:,:,1:LD66)
         ENDIF
      ENDIF

      !=================================================================
      ! ND67 diagnostic: I-6 fields (at surface)
      !
      ! (14) ALBD   : Surface Albedo                     [unitless]
      ! (17) TROPP  : Tropopause pressure                [hPa] 
      ! (18) SLP    : Sea Level pressure                 [hPa]
      !=================================================================
      IF ( ND67 > 0 ) THEN 
         IF( PRESENT( ALBD  ) ) AD67(:,:,14) = AD67(:,:,14) + ALBD
         IF( PRESENT( TROPP ) ) AD67(:,:,17) = AD67(:,:,17) + TROPP
         IF( PRESENT( SLP   ) ) AD67(:,:,18) = AD67(:,:,18) + SLP
      ENDIF 

      ! Return to calling program      
      END SUBROUTINE READ_I6

!------------------------------------------------------------------------------

      SUBROUTINE I6_CHECK( NFOUND, N_I6 )
!
!******************************************************************************
!  Subroutine I6_CHECK prints an error message if not all of the I-6 met 
!  fields are found.  The run is also terminated. (bmy, 10/27/00, 6/23/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NFOUND (INTEGER) : # of I-6 met fields read from disk
!  (2 ) N_I6   (INTEGER) : # of I-6 met fields expected to be read from disk
!
!  NOTES
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/23/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NFOUND, N_I6

      !=================================================================
      ! I6_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_I6 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough I-6 fields found!'      

         WRITE( 6, 120   ) N_I6, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in I6_CHECK (i6_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE I6_CHECK

!------------------------------------------------------------------------------

      ! End of module
      END MODULE I6_READ_MOD
