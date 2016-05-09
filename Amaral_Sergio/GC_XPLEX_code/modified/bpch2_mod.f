! $Id: bpch2_mod.f,v 1.3 2010/03/09 15:03:47 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bpch2_mod.f
!
! !DESCRIPTION: Module BPCH2\_MOD contains the routines used to read data 
!  from and write data to binary punch (BPCH) file format (v. 2.0).
!\\
!\\
! !INTERFACE: 
!
      MODULE BPCH2_MOD
! 
! !USES:
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: OPEN_BPCH2_FOR_READ 
      PUBLIC  :: OPEN_BPCH2_FOR_WRITE 
      PUBLIC  :: BPCH2_HDR           
      PUBLIC  :: BPCH2               
      PUBLIC  :: READ_BPCH2          
      PUBLIC  :: GET_MODELNAME       
      PUBLIC  :: GET_NAME_EXT        
      PUBLIC  :: GET_NAME_EXT_2D     
      PUBLIC  :: GET_RES_EXT         
      PUBLIC  :: GET_HALFPOLAR       
      PUBLIC  :: GET_TAU0
      ! adj_group
      PUBLIC  :: BPCH3               

      INTERFACE GET_TAU0
         MODULE PROCEDURE GET_TAU0_6A
      END INTERFACE

!      INTERFACE READ_BPCH2
!         MODULE PROCEDURE READ_BPCH2_C
!         MODULE PROCEDURE READ_BPCH2_C2
!         MODULE PROCEDURE READ_BPCH2_R
!         MODULE PROCEDURE READ_BPCH2_R2
!         MODULE PROCEDURE READ_BPCH2_SCALAR
!      END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_TAU0_6A
!      PRIVATE :: READ_BPCH2_R
!      PRIVATE :: READ_BPCH2_R2
!      PRIVATE :: READ_BPCH2_C
!      PRIVATE :: READ_BPCH2_C2
!      PRIVATE :: READ_BPCH2_SCALAR

! !REVISION HISTORY:
!  (1 ) Added routine GET_TAU0 (bmy, 7/20/00)
!  (2 ) Added years 1985-2001 for routine GET_TAU0 (bmy, 8/1/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF (bmy, 9/12/00)
!  (4 ) Removed obsolete code in "read_bpch2.f" (bmy, 12/18/00)
!  (5 ) Correct error for 1991 TAU values in GET_TAU0 (bnd, bmy, 1/4/01)
!  (6 ) BPCH2_MOD is now independent of any GEOS-CHEM size parameters.
!        (bmy, 4/18/01)
!  (7 ) Now have 2 versions of "GET_TAU0" overloaded by an interface.  The
!        original version takes 2 arguments (MONTH, YEAR).  The new version
!        takes 3 arguments (MONTH, DAY, YEAR). (bmy, 8/22/01)
!  (8 ) Updated comments (bmy, 9/4/01)
!  (9 ) Renamed GET_TAU0_3A to GET_TAU0_6A, and updated the GET_TAU0 
!        interface.  Also updated comments (bmy, 9/26/01)
!  (10) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (11) Minor bug fix in GET_TAU0_2A.  Also deleted obsolete code from 9/01.
!        (bmy, 11/15/01)
!  (12) Moved routines JULDAY, MINT, CALDATE to "julian_mod.f".  Now 
!        references routine JULDAY from "julday_mod.f".  Also added code
!        for GEOS-4/fvDAS model type. (bmy, 11/20/01)
!  (23) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (24) Added OPEN_BPCH2_FOR_READ and OPEN_BPCH2_FOR_WRITE.  Also now 
!        reference IU_FILE and IOERROR from "file_mod.f". (bmy, 7/30/02)
!  (25) Now references "error_mod.f".  Also obsoleted routine GET_TAU0_2A.
!        (bmy, 10/15/02)
!  (26) Made modification in READ_BPCH2 for 1x1 nested grids (bmy, 3/11/03)
!  (27) Modifications for GEOS-4, 30-layer grid (bmy, 11/3/03)
!  (28) Added cpp switches for GEOS-4 1x125 grid (bmy, 12/1/04)
!  (29) Modified for GCAP and GEOS-5 met fields.  Added function
!        GET_HALFPOLAR. (bmy, 6/28/05)
!  (30) Added GET_NAME_EXT_2D to get filename extension for files which do
!        not contain any vertical information (bmy, 8/16/05)
!  (31) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (32) Renamed GRID30LEV to GRIDREDUCED.  Also increase TEMPARRAY in
!        READ_BPCH2 for GEOS-5 vertical levels. (bmy, 2/16/07)
!  (33) Modifications for GEOS-5 nested grids (bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_bpch2_for_read
!
! !DESCRIPTION: Subroutine OPEN\_BPCH2\_FOR\_READ opens a binary punch file 
!  (version 2.0 format) for reading only.  Also reads FTI and TITLE strings. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_BPCH2_FOR_READ( IUNIT, FILENAME, TITLE )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER,           INTENT(IN)            :: IUNIT     ! LUN for file I/O
      CHARACTER(LEN=*),  INTENT(IN)            :: FILENAME  ! Name of file
!
! !OUTPUT PARAMETERS:
!
      CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: TITLE     ! File title string
!
! !REVISION HISTORY: 
!  (1 ) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                                  :: IOS
      CHARACTER(LEN=40)                        :: FTI
      CHARACTER(LEN=80)                        :: TMP_TITLE

      !=================================================================
      ! OPEN_BPCH2_FOR_READ begins here!
      !=================================================================

      ! Open file for input -- readonly
      OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) THEN
         WRITE(6,*)'Error opening filename=',trim(filename)
         CALL FLUSH(6)
         CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:1')
      ENDIF

      
      ! Read file type identifier
      READ( IUNIT, IOSTAT=IOS ) FTI

      ! Error check
      IF ( IOS /= 0 ) THEN
         WRITE(6,*)'Error reading FTI for filename=',trim(filename)
         CALL FLUSH(6)
         CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:2' )
      ENDIF
         
      ! Stop if this is not a binary punch file
      IF ( TRIM( FTI ) /= 'CTM bin 02' ) THEN
         WRITE(6,*)'Error filename=',trim(filename)
         CALL FLUSH(6)
         CALL ERROR_STOP( 'Invalid file format!', 
     &                    'OPEN_BPCH2_FOR_READ (bpch2_mod.f)')
      ENDIF

      
      ! Read top title
      READ( IUNIT, IOSTAT=IOS ) TMP_TITLE

      ! Error check
      IF ( IOS /= 0 ) THEN
         WRITE(6,*)'Error reading filename=',trim(filename)
         CALL FLUSH(6)
         CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:3' )
      ENDIF
   

      ! Copy value of TMP_TITLE to TITLE for return 
      IF ( PRESENT( TITLE ) ) TITLE = TMP_TITLE

      ! Return to calling program
      END SUBROUTINE OPEN_BPCH2_FOR_READ
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_bpch2_for_write
!
! !DESCRIPTION: Subroutine OPEN\_BPCH2\_FOR\_WRITE opens a binary punch file
!  (version 2.0) for writing.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_BPCH2_FOR_WRITE( IUNIT, FILENAME, TITLE )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER,           INTENT(IN)            :: IUNIT     ! LUN for file I/O
      CHARACTER(LEN=*),  INTENT(IN)            :: FILENAME  ! Name of file
!
! !OUTPUT PARAMETERS:
!
      CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: TITLE     ! File title string
!
! !REVISION HISTORY:
!  30 Jul 2002 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER                                 :: IOS
      CHARACTER(LEN=80)                       :: TMP_TITLE

      !=================================================================
      ! OPEN_BPCH2_FOR_WRITE begins here!
      !=================================================================

      ! If TITLE is not passed, create a default title string
      IF ( PRESENT( TITLE ) ) THEN
         TMP_TITLE = TITLE
      ELSE
         TMP_TITLE = 'GEOS-CHEM binary punch file v. 2.0'
      ENDIF

      ! Open file for output
      OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) ' '
         WRITE(6,*) "CANNOT WRITE : " // FILENAME
         CALL IOERROR( IOS, IUNIT,'open_bpch2_for_write:1')
      ENDIF
         

      ! Write the top-of-file title to disk
      CALL BPCH2_HDR( IUNIT, TMP_TITLE )

      ! Return to calling program
      END SUBROUTINE OPEN_BPCH2_FOR_WRITE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bpch2_hdr
!
! !DESCRIPTION: Subroutine BPCH2\_HDR writes a header at the top of the binary
!  punch file, version 2.0.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BPCH2_HDR ( IUNIT, TITLE )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER,           INTENT(IN) :: IUNIT   ! LUN for file I/O
      CHARACTER(LEN=80), INTENT(IN) :: TITLE   ! Top-of-file title string
!
! !REVISION HISTORY:
!  (1 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (2 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (3 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER                       :: IOS
      CHARACTER(LEN=40)             :: FTI = 'CTM bin 02'

      !=================================================================
      ! BPCH2_HDR begins here!
      !
      ! Write header information to binary punch file 
      ! Also be sure to trap I/O Error conditions
      !=================================================================
      WRITE ( IUNIT, IOSTAT=IOS ) FTI
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:1' )

      WRITE ( IUNIT, IOSTAT=IOS ) TITLE
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:2' )

      ! Return to calling program    
      END SUBROUTINE BPCH2_HDR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bpch2
!
! !DESCRIPTION: Subroutine BPCH2 writes binary punch file (version 2.0) to 
!  disk.  Information about the model grid is also stored with each data block.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BPCH2( IUNIT,     MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NTRACER,    
     &                  UNIT,      TAU0,      TAU1,     RESERVED,   
     &                  NI,        NJ,        NL,       IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER,           INTENT(IN) :: IUNIT            ! LUN for file I/O
      CHARACTER(LEN=20), INTENT(IN) :: MODELNAME        ! Met field type
      TYPE (XPLEX),            INTENT(IN) :: LONRES           ! Lon resolution [deg]
      TYPE (XPLEX),            INTENT(IN) :: LATRES           ! Lat resolution [deg]
      INTEGER,           INTENT(IN) :: HALFPOLAR        ! 1/2-size polar boxes?
      INTEGER,           INTENT(IN) :: CENTER180        ! 1st box center -180?
      CHARACTER(LEN=40), INTENT(IN) :: CATEGORY         ! Diag. category name
      INTEGER,           INTENT(IN) :: NTRACER          ! Tracer index #
      CHARACTER(LEN=40), INTENT(IN) :: UNIT             ! Unit string
      TYPE (XPLEX),            INTENT(IN) :: TAU0             ! TAU values @ start &
      TYPE (XPLEX),            INTENT(IN) :: TAU1             !  end of diag interval
      CHARACTER(LEN=40), INTENT(IN) :: RESERVED         ! Extra string
      INTEGER,           INTENT(IN) :: NI, NJ, NL       ! Dimensions of ARRAY
      INTEGER,           INTENT(IN) :: IFIRST           ! (I,J,L) indices of
      INTEGER,           INTENT(IN) :: JFIRST           !  the first grid box
      INTEGER,           INTENT(IN) :: LFIRST           !  in Fortran notation
      TYPE (XPLEX),            INTENT(IN) :: ARRAY(NI,NJ,NL)  ! Data array
!
! !REVISION HISTORY:
!  (1 ) Added indices to IOERROR calls (e.g. "bpch2:1", "bpch2:2", etc.) 
!        (bmy, 10/4/99)
!  (2 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (4 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER                       :: I, J, L, NSKIP, IOS
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER            :: BYTES_PER_NUMBER = 16
      INTEGER, PARAMETER            :: END_OF_RECORD    = 8

      !=================================================================
      ! BPCH2 begins here!!  
      !
      ! Compute the number of bytes to skip between the end of one 
      ! data block and the beginning of the next data header line
      !=================================================================
      NSKIP = ( BYTES_PER_NUMBER * ( NI * NJ * NL ) ) + END_OF_RECORD

      !=================================================================
      ! Write data block to binary punch file
      ! Check for I/O errors
      !=================================================================
      !print*,'LONRES',LONRES
      !print*,'LATRES',LATRES
      !print*,'TAU0',TAU0
      !print*,'TAU1',TAU1
      WRITE( IUNIT, IOSTAT=IOS ) 
     & MODELNAME, cmplx(LONRES%r,LONRES%i),cmplx(LATRES%r,LATRES%i), 
     & HALFPOLAR,
     & CENTER180

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:1' )

      WRITE( IUNIT, IOSTAT = IOS ) 
     &     CATEGORY, NTRACER,  UNIT, cmplx(TAU0%r,TAU0%i),
     & cmplx(TAU1%r,TAU1%i),   RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:2' )

      WRITE( IUNIT, IOSTAT=IOS ) 
     &     ( ( ( cmplx(ARRAY(I,J,L)%r,ARRAY(I,J,L)%i), I=1,NI ),
     & J=1,NJ ), L=1,NL )

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:3' )

      !=================================================================
      ! Return to calling program      
      !=================================================================
      END SUBROUTINE BPCH2
!EOC
! adj_group (dkh, 03/07/10) 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bpch3
!
! !DESCRIPTION: Subroutine BPCH3 writes binary punch file (version 2.0) to 
!  disk.  Information about the model grid is also stored with each data block.
!  Just like BPCH2, except use TYPE (XPLEX). Based on BPCH2. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BPCH3( IUNIT,     MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NTRACER,    
     &                  UNIT,      TAU0,      TAU1,     RESERVED,   
     &                  NI,        NJ,        NL,       IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER,           INTENT(IN) :: IUNIT            ! LUN for file I/O
      CHARACTER(LEN=20), INTENT(IN) :: MODELNAME        ! Met field type
      TYPE (XPLEX),            INTENT(IN) :: LONRES           ! Lon resolution [deg]
      TYPE (XPLEX),            INTENT(IN) :: LATRES           ! Lat resolution [deg]
      INTEGER,           INTENT(IN) :: HALFPOLAR        ! 1/2-size polar boxes?
      INTEGER,           INTENT(IN) :: CENTER180        ! 1st box center -180?
      CHARACTER(LEN=40), INTENT(IN) :: CATEGORY         ! Diag. category name
      INTEGER,           INTENT(IN) :: NTRACER          ! Tracer index #
      CHARACTER(LEN=40), INTENT(IN) :: UNIT             ! Unit string
      TYPE (XPLEX),            INTENT(IN) :: TAU0             ! TAU values @ start &
      TYPE (XPLEX),            INTENT(IN) :: TAU1             !  end of diag interval
      CHARACTER(LEN=40), INTENT(IN) :: RESERVED         ! Extra string
      INTEGER,           INTENT(IN) :: NI, NJ, NL       ! Dimensions of ARRAY
      INTEGER,           INTENT(IN) :: IFIRST           ! (I,J,L) indices of
      INTEGER,           INTENT(IN) :: JFIRST           !  the first grid box
      INTEGER,           INTENT(IN) :: LFIRST           !  in Fortran notation
      TYPE (XPLEX),            INTENT(IN) :: ARRAY(NI,NJ,NL)  ! Data array
!
! !REVISION HISTORY:
!  (1 ) See BPCH2 
!  (2 ) Updated to v8 (dkh, 03/07/10) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER                       :: I, J, L, NSKIP, IOS
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER            :: BYTES_PER_NUMBER = 16
      INTEGER, PARAMETER            :: END_OF_RECORD    = 8

      !=================================================================
      ! BPCH3 begins here!!  
      !
      ! Compute the number of bytes to skip between the end of one 
      ! data block and the beginning of the next data header line
      !=================================================================
      NSKIP = ( BYTES_PER_NUMBER * ( NI * NJ * NL ) ) + END_OF_RECORD

      !=================================================================
      ! Write data block to binary punch file
      ! Check for I/O errors
      !=================================================================
      WRITE( IUNIT, IOSTAT=IOS ) 
     &     MODELNAME, cmplx(LONRES%r,LONRES%i),
     &   cmplx(LATRES%r,LATRES%i), HALFPOLAR, CENTER180

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch3:1' )

      WRITE( IUNIT, IOSTAT = IOS ) 
     &     CATEGORY, NTRACER,  UNIT, cmplx(TAU0%r,TAU0%i),
     &    cmplx(TAU1%r,TAU1%i),   RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch3:2' )

      WRITE( IUNIT, IOSTAT=IOS ) 
     & (((cmplx(ARRAY(I,J,L)%r,ARRAY(I,J,L)%i),I=1,NI),J=1,NJ),L=1,NL)

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch3:3' )

      !=================================================================
      ! Return to calling program      
      !=================================================================
      END SUBROUTINE BPCH3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_bpch2
!
! !DESCRIPTION: Subroutine READ\_BPCH2 reads a binary punch file (v. 2.0) 
!  and extracts a data block that matches the given category, tracer, and 
!  tau value.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BPCH2( FILENAME, CATEGORY_IN, TRACER_IN, 
     &                       TAU0_IN,  IX,          JX,          
     &                       LX,       ARRAY,       QUIET ) 
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE, IOERROR

#     include "define.h" 
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME         ! Bpch file to read
      CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN      ! Diag. category name
      INTEGER,           INTENT(IN)  :: TRACER_IN        ! Tracer index #
      TYPE (XPLEX),        INTENT(IN)  :: TAU0_IN          ! TAU timestamp 
      INTEGER,           INTENT(IN)  :: IX, JX, LX       ! Dimensions of ARRAY
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET            ! Don't print output
!
! !OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),            INTENT(OUT) :: ARRAY(IX,JX,LX)  ! Data array from file
!
! !REVISION HISTORY:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy, 3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy, 4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!  (10) Now reference IU_FILE and IOERROR from "file_mod.f".  Now call 
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use IU_FILE
!        as the unit number instead of a locally-defined IUNIT. (bmy, 7/30/02)
!  (11) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (12) Now set IFIRST=1, JFIRST=1 for 1x1 nested grids.  Now needs to
!        reference "define.h".  Added OPTIONAL QUIET flag. (bmy, 3/14/03)
!  (13) Now separate off nested grid code in an #ifdef block using
!        NESTED_CH or NESTED_NA cpp switches (bmy, 12/1/04)
!  (14) Make TEMPARRAY big enough for GEOS-5 72 levels (and 73 edges) 
!        (bmy, 2/15/07)
!  (15) Make TEMPARRAY large enough for 0.5 x 0.666 arrays -- but only if we
!        are doing a 0.5 x 0.666 nested simulation. (yxw, dan, bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG
      
      ! Make TEMPARRAY big enough for a global grid.  For 0.5 x 0.666 nested
      ! grid simulations we need to define this as 540x361x73.  However, this
      ! may cause memory problems on some Linux boxes for people who want to
      ! run only the global simulations.  Therefore increase the size of 
      ! TEMPARRAY only if we are doing a 0.5 x 0.666 nested simulation.
      ! (yxw, bmy, dan, 11/6/08)
#if   defined( GRID05x0666 ) 
      REAL*4             :: TEMPARRAY(540,361,73)   
#else
      REAL*4             :: TEMPARRAY(360,181,73)   
#endif

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: D_LONRES, D_LATRES
      TYPE (XPLEX)             :: LONRES,    LATRES
      REAL*8             :: D_ZTAU0, D_ZTAU1
      TYPE (XPLEX)             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT     
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY(:,:,:)     = 0d0
      TEMPARRAY(:,:,:) = 0d0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS ) 
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES = (D_LONRES)
         LATRES = (D_LATRES) 
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0 = (D_ZTAU0)
         ZTAU1 = (D_ZTAU1) 
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )
         
         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and. 
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN 

#if   defined( GRID1x1 )   || defined( GRID05x0666 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     
 
         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1
                  
         ARRAY(I1:I2,J1:J2,L1:L2)=(TEMPARRAY(1:NI,1:NJ,1:NL))

         ! Flag to decide whether or not we will echo info (bmy, 3/14/03)
         IF ( .not. TMP_QUIET ) THEN 
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2, 
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2

      SUBROUTINE READ_BPCH2_SCALAR( FILENAME, CATEGORY_IN, TRACER_IN,
     &                       TAU0_IN,  IX,          JX,
     &                       LX,       ARRAY,       QUIET )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE, IOERROR

#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME         ! Bpch file to      Cread
      CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN      ! Diag.      Ccategory name
      INTEGER,           INTENT(IN)  :: TRACER_IN        ! Tracer index#
      TYPE (XPLEX),        INTENT(IN)  :: TAU0_IN          ! TAU timestamp 
      INTEGER,           INTENT(IN)  :: IX, JX, LX       ! Dimensions ofARRAY
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET            ! Don't printoutput
!
! !OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),            INTENT(OUT) :: ARRAY  ! Data      Carray from file
!
! !REVISION HISTORY:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy,
!  3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy,
!  9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy,
!  4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!  (10) Now reference IU_FILE and IOERROR from "file_mod.f".  Now call 
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use
!        IU_FILE
!        as the unit number instead of a locally-defined IUNIT. (bmy,
!        7/30/02)
!  (11) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (12) Now set IFIRST=1, JFIRST=1 for 1x1 nested grids.  Now needs to
!        reference "define.h".  Added OPTIONAL QUIET flag. (bmy,
!        3/14/03)
!  (13) Now separate off nested grid code in an #ifdef block using
!        NESTED_CH or NESTED_NA cpp switches (bmy, 12/1/04)
!  (14) Make TEMPARRAY big enough for GEOS-5 72 levels (and 73 edges) 
!        (bmy, 2/15/07)
!  (15) Make TEMPARRAY large enough for 0.5 x 0.666 arrays -- but only
!  if we
!        are doing a 0.5 x 0.666 nested simulation. (yxw, dan, bmy,
!        11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG

      ! Make TEMPARRAY big enough for a global grid.  For 0.5 x 0.666
      ! nested
      ! grid simulations we need to define this as 540x361x73.  However,
      ! this
      ! may cause memory problems on some Linux boxes for people who
      ! want to
      ! run only the global simulations.  Therefore increase the size of 
      ! TEMPARRAY only if we are doing a 0.5 x 0.666 nested simulation.
      ! (yxw, bmy, dan, 11/6/08)
#if   defined( GRID05x0666 ) 
      REAL*4             :: TEMPARRAY
#else
      REAL*4             :: TEMPARRAY
#endif

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: D_LONRES, D_LATRES
      TYPE (XPLEX)             :: LONRES,    LATRES
      REAL*8             :: D_ZTAU0, D_ZTAU1
      TYPE (XPLEX)             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY     = 0d0
      TEMPARRAY = 0d0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right
      ! format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LATRES = (D_LATRES)
         LONRES = (D_LONRES)
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0 = (D_ZTAU0)
         ZTAU1 = (D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS )
     &        TEMPARRAY

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )

         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and.
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN

#if   defined( GRID1x1 )   || defined( GRID05x0666 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     

         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1

         ARRAY = (TEMPARRAY)

         ! Flag to decide whether or not we will echo info (bmy,
         ! 3/14/03)
         IF ( .not. TMP_QUIET ) THEN
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2,
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2_SCALAR


      SUBROUTINE READ_BPCH2_C2( FILENAME, CATEGORY_IN, TRACER_IN,
     &                       TAU0_IN,  IX,          JX,
     &                       LX,       ARRAY,       QUIET )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE, IOERROR

#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME         ! Bpch file to      Cread
      CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN      ! Diag.      Ccategory name
      INTEGER,           INTENT(IN)  :: TRACER_IN        ! Tracer index#
      TYPE (XPLEX),        INTENT(IN)  :: TAU0_IN          ! TAU timestamp 
      INTEGER,           INTENT(IN)  :: IX, JX, LX       ! Dimensions ofARRAY
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET            ! Don't printoutput
!
! !OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),            INTENT(OUT) :: ARRAY(IX,JX)  ! Data      Carray from file
!
! !REVISION HISTORY:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy,
!  3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy,
!  9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy,
!  4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!  (10) Now reference IU_FILE and IOERROR from "file_mod.f".  Now call 
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use
!        IU_FILE
!        as the unit number instead of a locally-defined IUNIT. (bmy,
!        7/30/02)
!  (11) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (12) Now set IFIRST=1, JFIRST=1 for 1x1 nested grids.  Now needs to
!        reference "define.h".  Added OPTIONAL QUIET flag. (bmy,
!        3/14/03)
!  (13) Now separate off nested grid code in an #ifdef block using
!        NESTED_CH or NESTED_NA cpp switches (bmy, 12/1/04)
!  (14) Make TEMPARRAY big enough for GEOS-5 72 levels (and 73 edges) 
!        (bmy, 2/15/07)
!  (15) Make TEMPARRAY large enough for 0.5 x 0.666 arrays -- but only
!  if we
!        are doing a 0.5 x 0.666 nested simulation. (yxw, dan, bmy,
!        11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG

      ! Make TEMPARRAY big enough for a global grid.  For 0.5 x 0.666
      ! nested
      ! grid simulations we need to define this as 540x361x73.  However,
      ! this
      ! may cause memory problems on some Linux boxes for people who
      ! want to
      ! run only the global simulations.  Therefore increase the size of 
      ! TEMPARRAY only if we are doing a 0.5 x 0.666 nested simulation.
      ! (yxw, bmy, dan, 11/6/08)
#if   defined( GRID05x0666 ) 
      REAL*4             :: TEMPARRAY(540,361)
#else
      REAL*4             :: TEMPARRAY(360,181)
#endif

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: D_LONRES, D_LATRES
      TYPE (XPLEX)             :: LONRES,    LATRES
      REAL*8             :: D_ZTAU1, D_ZTAU0
      TYPE (XPLEX)             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY(:,:)     = 0d0
      TEMPARRAY(:,:) = 0d0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right
      ! format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LATRES = (D_LATRES)
         LONRES = (D_LONRES)
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0 = (D_ZTAU0)
         ZTAU1 = (D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS )
     &         ( ( TEMPARRAY(I,J), I=1,NI ), J=1,NJ ) 

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )

         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and.
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN

#if   defined( GRID1x1 )   || defined( GRID05x0666 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     

         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1

         ARRAY( I1:I2, J1:J2) = (TEMPARRAY( 1:NI, 1:NJ))

         ! Flag to decide whether or not we will echo info (bmy,
         ! 3/14/03)
         IF ( .not. TMP_QUIET ) THEN
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2,
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2_C2

      SUBROUTINE READ_BPCH2_R2( FILENAME, CATEGORY_IN, TRACER_IN,
     &                       TAU0_IN,  IX,          JX,
     &                       LX,       ARRAY,       QUIET )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE, IOERROR

#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME         ! Bpch file to      CCread
      CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN      ! Diag.      CCcategory name
      INTEGER,           INTENT(IN)  :: TRACER_IN        ! Tracer index#
      REAL*8,        INTENT(IN)  :: TAU0_IN          ! TAU timestamp 
      INTEGER,           INTENT(IN)  :: IX, JX, LX       ! DimensionsofARRAY
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET            ! Don'tprintoutput
!
! !OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),            INTENT(OUT) :: ARRAY(IX,JX)  ! Data      CCarray from file
!
! !REVISION HISTORY:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy,
!  3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy,
!  9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy,
!  4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!  (10) Now reference IU_FILE and IOERROR from "file_mod.f".  Now call 
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use
!        IU_FILE
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use
!        IU_FILE
!        as the unit number instead of a locally-defined IUNIT. (bmy,
!        7/30/02)
!  (11) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (12) Now set IFIRST=1, JFIRST=1 for 1x1 nested grids.  Now needs to
!        reference "define.h".  Added OPTIONAL QUIET flag. (bmy,
!        3/14/03)
!  (13) Now separate off nested grid code in an #ifdef block using
!        NESTED_CH or NESTED_NA cpp switches (bmy, 12/1/04)
!  (14) Make TEMPARRAY big enough for GEOS-5 72 levels (and 73 edges) 
!        (bmy, 2/15/07)
!  (15) Make TEMPARRAY large enough for 0.5 x 0.666 arrays -- but only
!  if we
!        are doing a 0.5 x 0.666 nested simulation. (yxw, dan, bmy,
!        11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG

      ! Make TEMPARRAY big enough for a global grid.  For 0.5 x 0.666
      ! nested
      ! grid simulations we need to define this as 540x361x73.  However,
      ! this
      ! may cause memory problems on some Linux boxes for people who
      ! want to
      ! run only the global simulations.  Therefore increase the size of 
      ! TEMPARRAY only if we are doing a 0.5 x 0.666 nested simulation.
      ! (yxw, bmy, dan, 11/6/08)
#if   defined( GRID05x0666 ) 
      REAL*4             :: TEMPARRAY(540,361)
#else
      REAL*4             :: TEMPARRAY(360,181)
#endif

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4            :: D_LONRES, D_LATRES
      TYPE (XPLEX)             :: LONRES,    LATRES
      REAL*8            :: D_ZTAU0, D_ZTAU1
      TYPE (XPLEX)             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY(:,:)     = 0d0
      TEMPARRAY(:,:) = 0d0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right
      ! format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LATRES = (D_LATRES) 
         LONRES = (D_LONRES)
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0 = (D_ZTAU0)
         ZTAU1 = (D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS )
     &         ( ( TEMPARRAY(I,J), I=1,NI ), J=1,NJ )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )

         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and.
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN

#if   defined( GRID1x1 )   || defined( GRID05x0666 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     

         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1

         ARRAY( I1:I2, J1:J2) = (TEMPARRAY( 1:NI, 1:NJ))

         ! Flag to decide whether or not we will echo info (bmy,
         ! 3/14/03)
         IF ( .not. TMP_QUIET ) THEN
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2,
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2_R2

      SUBROUTINE READ_BPCH2_R( FILENAME, CATEGORY_IN, TRACER_IN,
     &                       TAU0_IN,  IX,          JX,
     &                       LX,       ARRAY,       QUIET )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE, IOERROR

#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME         ! Bpch file to      Cread
      CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN      ! Diag.      Ccategory name
      INTEGER,           INTENT(IN)  :: TRACER_IN        ! Tracer index#
      REAL*8,            INTENT(IN)  :: TAU0_IN          ! TAU  Ctimestamp 
      INTEGER,           INTENT(IN)  :: IX, JX, LX       ! Dimensions of ARRAY
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET            ! Don't print output
!
! !OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),            INTENT(OUT) :: ARRAY(IX,JX,LX)  ! Data      Carray from file
!
! !REVISION HISTORY:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy,
!  3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy,
!  9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy,
!  4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!  (10) Now reference IU_FILE and IOERROR from "file_mod.f".  Now call 
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use
!        IU_FILE
!        as the unit number instead of a locally-defined IUNIT. (bmy,
!        7/30/02)
!  (11) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (12) Now set IFIRST=1, JFIRST=1 for 1x1 nested grids.  Now needs to
!        reference "define.h".  Added OPTIONAL QUIET flag. (bmy,
!        3/14/03)
!  (13) Now separate off nested grid code in an #ifdef block using
!        NESTED_CH or NESTED_NA cpp switches (bmy, 12/1/04)
!  (14) Make TEMPARRAY big enough for GEOS-5 72 levels (and 73 edges) 
!        (bmy, 2/15/07)
!  (15) Make TEMPARRAY large enough for 0.5 x 0.666 arrays -- but only
!  if we
!        are doing a 0.5 x 0.666 nested simulation. (yxw, dan, bmy,
!        11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG

      ! Make TEMPARRAY big enough for a global grid.  For 0.5 x 0.666
      ! nested
      ! grid simulations we need to define this as 540x361x73.  However,
      ! this
      ! may cause memory problems on some Linux boxes for people who
      ! want to
      ! run only the global simulations.  Therefore increase the size of 
      ! TEMPARRAY only if we are doing a 0.5 x 0.666 nested simulation.
      ! (yxw, bmy, dan, 11/6/08)
#if   defined( GRID05x0666 ) 
      REAL*4             :: TEMPARRAY(540,361,73)
#else
      REAL*4             :: TEMPARRAY(360,181,73)
#endif

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: D_LONRES, D_LATRES
      TYPE (XPLEX)             :: LONRES,    LATRES
      REAL*8             :: D_ZTAU0, D_ZTAU1
      TYPE (XPLEX)             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY(:,:,:)     = 0d0
      TEMPARRAY(:,:,:) = 0d0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right
      ! format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LATRES = (D_LATRES)
         LONRES = (D_LONRES)
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0 = (D_ZTAU0)
         ZTAU1 = (D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS )
     &        ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )

         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and.
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN

#if   defined( GRID1x1 )   || defined( GRID05x0666 )

#if   defined( NESTED_CH ) || defined( NESTED_NA )
         ! *** NOTE: now use NESTED_CH or NESTED_NA cpp switches ***
         ! *** to block off this section of code (bmy, 12/1/04)  ***
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#endif

#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     

         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1

         ARRAY(I1:I2,J1:J2,L1:L2)=(TEMPARRAY(1:NI,1:NJ,1:NL))

         ! Flag to decide whether or not we will echo info (bmy,
         ! 3/14/03)
         IF ( .not. TMP_QUIET ) THEN
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2,
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2_R



!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_modelname
!
! !DESCRIPTION: Function GET\_MODELNAME returns the proper value of MODELNAME 
!  for current GEOS or GCAP met field type.  MODELNAME is written to the 
!  binary punch file and is also used by the GAMAP package.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MODELNAME() RESULT( MODELNAME )
!
! !USES:
!
#     include "CMN_SIZE"
!
! !RETURN_VALUE:
!
      CHARACTER(LEN=20) :: MODELNAME   ! Model name for the current met field
!
! !REVISION HISTORY:
!  (1 ) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (2 ) Added modelname for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (3 ) Added "GEOS4_30L" for reduced GEOS-4 grid.  Also now use C-preprocessor
!        switch "GRID30LEV" instead of IF statements. (bmy, 11/3/03)
!  (4 ) Updated for GCAP and GEOS-5 met fields.  Rearranged coding for
!        simplicity. (swu, bmy, 5/24/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Rename GRID30LEV to GRIDREDUCED (bmy, 2/7/07)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_3 ) && defined( GRIDREDUCED )
      MODELNAME = 'GEOS3_30L'

#elif defined( GEOS_3 )
      MODELNAME = 'GEOS3'

#elif defined( GEOS_4 ) && defined( GRIDREDUCED )
      MODELNAME = 'GEOS4_30L'

#elif defined( GEOS_4 )
      MODELNAME = 'GEOS4'

#elif defined( GEOS_5 ) && defined( GRIDREDUCED )
      MODELNAME = 'GEOS5_47L'
      
#elif defined( GEOS_5 ) 
      MODELNAME = 'GEOS5'

#elif defined( GCAP )
      MODELNAME = 'GCAP'

#endif

      ! Return to calling program
      END FUNCTION GET_MODELNAME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_name_ext
!
! !DESCRIPTION: Function GET\_NAME\_EXT returns the proper filename extension 
!  the current GEOS-Chem met field type (e.g. "geos3", "geos4", "geos5", or 
!  "gcap").  
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NAME_EXT() RESULT( NAME_EXT )
!
! !USES:
!
#     include "define.h"
!
! !RETURN VALUE:
!
#if   defined( GEOS_3 )
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos3'

#elif defined( GEOS_4 )
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos4'

#elif defined( GEOS_5 )
      CHARACTER(LEN=5) :: NAME_EXT
      NAME_EXT = 'geos5'

#elif defined( GCAP )
      CHARACTER(LEN=4) :: NAME_EXT
      NAME_EXT = 'gcap'

#endif
!
! !REVISION HISTORY:
!  (1 ) Added name string for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (2 ) Remove obsolete "geos2" model name strning (bmy, 11/3/03)
!  (3 ) Modified for GCAP and GEOS-5 met fields (bmy, 5/24/05)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      END FUNCTION GET_NAME_EXT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_name_ext_2d
!
! !DESCRIPTION: Function GET\_NAME\_EXT\_2D returns the proper filename 
!  extension for CTM model name for files which do not contain any vertical 
!  information (i.e. "geos" or "gcap").
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NAME_EXT_2D() RESULT( NAME_EXT_2D )
!
! !RETURN VALUE:
!
      CHARACTER(LEN=4) :: NAME_EXT_2D
!
! !REVISION HISTORY:
!  (1 ) Added name string for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (2 ) Remove obsolete "geos2" model name strning (bmy, 11/3/03)
!  (3 ) Modified for GCAP and GEOS-5 met fields (bmy, 5/24/05)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
!
      ! Local variables
      CHARACTER(LEN=5) :: TEMP_NAME

      !=================================================================
      ! GET_NAME_EXT_2D begins here!
      !=================================================================

      ! Get the name extension
      TEMP_NAME   = GET_NAME_EXT()

      ! Take the 1st 4 characters ("geos" or "gcap") and return
      NAME_EXT_2D = TEMP_NAME(1:4)

      ! Return to calling program
      END FUNCTION GET_NAME_EXT_2D 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_res_ext
!
! !DESCRIPTION: Function GET\_RES\_EXT returns the proper filename extension 
!  for the GEOS-Chem horizontal grid resolution (e.g. "1x1", "2x25", "4x5").
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_RES_EXT() RESULT( RES_EXT )
!
! !USES:
!
#     include "define.h"
!
! !RETURN VALUE:
!
#if   defined( GRID4x5 )
      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '4x5'
     
#elif defined( GRID2x25 ) 
      CHARACTER(LEN=4) :: RES_EXT
      RES_EXT = '2x25'

#elif defined( GRID1x125 )
      CHARACTER(LEN=5) :: RES_EXT
      RES_EXT = '1x125'

#elif defined( GRID1x1 ) 
      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '1x1'

#elif defined( GRID05x0666 )
      CHARACTER(LEN=7) :: RES_EXT
      RES_EXT = '05x0666'

#endif
!
! !REVISION HISTORY:
!  (1 ) Added extension for 1 x 1.25 grid (bmy, 12/1/04)
!  (2 ) Added extension for 0.5 x 0.666 grid (yxw, dan, bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      END FUNCTION GET_RES_EXT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_halfpolar
!
! !DESCRIPTION: Function GET\_HALFPOLAR returns 1 if the current grid has 
!  half-sized polar boxes (e.g. GEOS), or zero otherwise (e.g. GCAP).  
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HALFPOLAR() RESULT( HALFPOLAR )
!
! !USES:
!
#     include "define.h"
!
! !RETURN VALUE:
!
      INTEGER :: HALFPOLAR  ! =1 if we have half-sized polar boxes, =0 if not
!
! !REVISION HISTORY:
!  28 Jun 2005 - S. Wu & R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca         - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
! 
#if   defined( GCAP ) 

      ! GCAP grid does not have half-sized polar boxes
      HALFPOLAR = 0

#else

      ! All GEOS grids have half-sized polar boxes
      HALFPOLAR = 1

#endif

      ! Return to calling program
      END FUNCTION GET_HALFPOLAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_tau0_6a
!
! !DESCRIPTION: Function GET\_TAU0\_6A returns the corresponding TAU0 value 
!  for the first day of a given MONTH of a given YEAR.  This is necessary to 
!  index monthly mean binary punch files, which are used as input to GEOS-Chem.
!\\
!\\
!  This function takes 3 mandatory arguments (MONTH, DAY, YEAR) and 3 
!  optional arguments (HOUR, MIN, SEC).  It is intended to replace the current 
!  2-argument version of GET\_TAU0.  The advantage being that GET\_TAU0\_6A 
!  can compute a TAU0 for any date and time in the GEOS-Chem epoch, rather 
!  than just the first day of each month.  Overload this w/ an interface so 
!  that the user can also choose the version of GET\_TAU0 w/ 2 arguments 
!  (MONTH, YEAR), which is the prior version.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TAU0_6A( MONTH, DAY, YEAR, 
     &                      HOUR,  MIN, SEC  ) RESULT( THIS_TAU0 )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE JULDAY_MOD, ONLY : JULDAY
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)           :: MONTH       
      INTEGER, INTENT(IN)           :: DAY         
      INTEGER, INTENT(IN)           :: YEAR        
      INTEGER, INTENT(IN), OPTIONAL :: HOUR        
      INTEGER, INTENT(IN), OPTIONAL :: MIN
      INTEGER, INTENT(IN), OPTIONAL :: SEC
!
! !RETURN VALUE:
!
      TYPE (XPLEX)                        :: THIS_TAU0   ! TAU0 timestamp
!
! !REMARKS:
!  TAU0 is hours elapsed since 00:00 GMT on 01 Jan 1985.
!
! !REVISION HISTORY:
!  (1 ) 1985 is the first year of the GEOS epoch.
!  (2 ) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!  (3 ) Correct error for 1991 TAU values.  Also added 2002 and 2003.
!        (bnd, bmy, 1/4/01)
!  (4 ) Updated comments  (bmy, 9/26/01)
!  (5 ) Now references JULDAY from "julday_mod.f" (bmy, 11/20/01)
!  (6 ) Now references ERROR_STOP from "error_mod.f"  (bmy, 10/15/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: TMP_HOUR, TMP_MIN, TMP_SEC
      TYPE (XPLEX)                        :: DAYS

      ! Return value
      
      !=================================================================
      ! GET_TAU0_6A begins here!
      !=================================================================

      ! Error checking 
      IF ( MONTH < 1 .or. MONTH > 12 ) THEN
         CALL ERROR_STOP ( 'Invalid MONTH selection!', 'GET_TAU0' )
      ENDIF

      ! Error checking 
      IF ( DAY < 1 .or. DAY > 31 ) THEN
         CALL ERROR_STOP ( 'Invalid DAY selection!', 'GET_TAU0' )
      ENDIF

      ! If HOUR isn't passed, default to 0
      IF ( PRESENT( HOUR ) ) THEN
         TMP_HOUR = HOUR
      ELSE
         TMP_HOUR = 0
      ENDIF 

      ! If MIN isn't passed, default to 0
      IF ( PRESENT( MIN ) ) THEN
         TMP_MIN = MIN
      ELSE
         TMP_MIN = 0 
      ENDIF 

      ! If SEC isn't passed, default to 0
      IF ( PRESENT( SEC ) ) THEN
         TMP_SEC = SEC
      ELSE
         TMP_SEC = 0 
      ENDIF 

      ! Number of days since midnight on 1/1/1985
      THIS_TAU0 = JULDAY( YEAR, MONTH, XPLEX( DAY,0d0 ) ) - 2446066.5d0

      ! Multiply by 24 to get hours since 1/1/1985
      ! Also add in the hours elapsed since midnight on this date
      THIS_TAU0 = ( THIS_TAU0 * 24d0 ) + ( TMP_HOUR         ) + 
     &            ( TMP_MIN   / 60d0 ) + ( TMP_SEC / 3600d0 )

      ! Return to calling program
      END FUNCTION GET_TAU0_6A
!EOC

      ! End of module
      END MODULE BPCH2_MOD
