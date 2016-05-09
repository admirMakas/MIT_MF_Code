! $Id: HdfVdModule.f90,v 1.2 2012/03/01 22:00:27 daven Exp $
MODULE HdfVdModule

  !===========================================================================
  ! Module "HdfVdModule" contains variables and methods that are used to
  ! read data fields stored in HDF-VDATA format. (bmy, 7/3/03, 12/12/05)
  !
  ! In order to use HdfVdModule, you must first install the HDF-4 library 
  ! on your system.  You may download the library source code from:
  !
  !     http://hdf.ncsa.uiuc.edu/hdf4.html
  !
  ! There is also a good online tutorial about the HDF-VD interface at:
  !
  !     http://hdf.ncsa.uiuc.edu/training/HDFtraining/tutorial/vd/vds.html
  !
  ! Module Variables: 
  ! --------------------------------------------------------------------------
  ! (1 ) fileId            : ID number for the HDF file 
  ! (2 ) saveFileName      : Shadow variable for filename 
  !
  ! Module Methods:
  ! --------------------------------------------------------------------------
  ! (1 ) vdOpen            : Opens the HDF file
  ! (2 ) vdClose           : Closes the HDF file
  ! (3 ) vdOpenField       : Opens access to a HDF-VDATA field w/in the file
  ! (4 ) vdCloseField      : Closes acess to a HDF-VDATA field w/in the file
  ! (5 ) vdPrintInfo       : Prints information about all HDF-VDATA fields
  ! (6 ) vdGetFieldDim     : Gets dimensions of a given HDF-VDATA field
  ! (7 ) vdGetDataR4       : Reads a 1-D TYPE (XPLEX) HDF-VDATA field from the file
  ! (8 ) vdGetDataR8       : Reads a 1-D TYPE (XPLEX) HDF-VDATA field from the file
  ! (9 ) vdShift           : Shifts a 1-D TYPE (XPLEX)  data field by 180 degrees
  ! (10) getTauFromDate    : Converts a date to a TAU value
  ! (11) calDate           : Converts Julian day to NYMD, NHMS
  ! (12) julDay            : Converts Year/month/day to Julian day
  ! (13) mint              : Function required by routine julDay
  !
  ! Module Interfaces: 
  ! --------------------------------------------------------------------------
  ! (1 ) vdGetData         : vdGetDataR4, vdGetDataR8
  !
  ! NOTES:
  ! (1 ) Based on HdfSdModule.f90 (bmy, 7/3/03)
  ! (2 ) Added function getTauFromDate (bmy, 12/12/05)
  !===========================================================================
  USE HdfIncludeModule
 
  USE MYTYPE
  USE COMPLEXIFY
  IMPLICIT NONE

  !=====================================================================
  ! MODULE PRIVATE DECLARATIONS
  !=====================================================================

  ! Make everything PRIVATE ...
  PRIVATE

  ! ... except these variables ...
  PUBLIC :: fileId
  PUBLIC :: saveFileName

  ! ... and these routines
  PUBLIC :: vdOpen       
  PUBLIC :: vdClose       
  PUBLIC :: vdCloseField   
  PUBLIC :: vdGetData
  PUBLIC :: vdGetFieldDim
  PUBLIC :: vdOpenField    
  PUBLIC :: vdPrintInfo  
  PUBLIC :: vdShift      
  PUBLIC :: getTauFromDate 

  !=====================================================================
  ! Private module variables -- visible only within HdfModule
  !=====================================================================
  INTEGER            :: fileId
  CHARACTER(LEN=255) :: saveFileName

  !=======================================================================
  ! Module interfaces: allow you to associate a name w/ several routines
  ! with different numbers of arguments or different argument types
  !=======================================================================
  INTERFACE vdGetData
     MODULE PROCEDURE vdGetDataInt
     MODULE PROCEDURE vdGetDataR4
     MODULE PROCEDURE vdGetDataR8
  END INTERFACE

  INTERFACE vdShift
     MODULE PROCEDURE vdShift1d
  END INTERFACE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE vdOpen( fileName )

    !=====================================================================
    ! Subroutine "vdOpen" opens an HDF file and initializes the 
    ! HDF-VDATA interface. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fileName : CHARACTER name of the HDF-EOS file to be opened
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) hOpen    : returns INTEGER value ( fileId )
    ! (2) vfStart  : returns INTEGER value ( vdId   )
    !
    ! NOTES: 
    !=====================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    CHARACTER(LEN=*)       :: fileName

    ! Local Variables
    INTEGER                :: status
    CHARACTER(LEN=255)     :: message

    ! External functions
    INTEGER,      EXTERNAL :: hOpen, vfStart

    !=====================================================================
    ! vdOpen begins here!
    !=====================================================================

    ! Save file name to a private shadow variable for error msgs
    saveFileName = TRIM( fileName )

    ! Open the HDF file
    fileId = hopen( TRIM( fileName ), DFACC_READ, 16 )
  
    ! Error check fileId
    IF ( fileId == FAIL ) THEN 
       message = 'ERROR: Could not open HDF file ' // TRIM( fileName )
       CALL ERROR_STOP( message, 'vdOpen' )
    ENDIF
    
    ! Start the VDATA interface for this file
    status = vfstart( fileId )

    ! Error check
    IF ( status == FAIL ) THEN 
       message = 'ERROR: Could not start HDF-VDATA interface for file ' // &
                  TRIM( fileName )
       CALL ERROR_STOP( message, 'vdOpen' )
    ENDIF

  END SUBROUTINE vdOpen

!------------------------------------------------------------------------------
  
  SUBROUTINE vdClose( fileName )

    !=====================================================================
    ! Subroutine "vdClose" terminates the HDF Scientific Dataset 
    ! (HDF-VD) interface and closes the HDF file. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fileName : CHARACTER name of the HDF-EOS file to be opened
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) hClose   : takes INTEGER value ( fileId )
    ! (2) vfEnd    : takes INTEGER value ( fileId )
    !
    ! NOTES:
    !=====================================================================
    
    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: fileName

    ! Local variables
    INTEGER                      :: status
    CHARACTER(LEN=255)           :: message

    ! External functions
    INTEGER, EXTERNAL            :: hClose, vfEnd

    !=====================================================================
    ! vdClose begins here!
    !=====================================================================

    ! Close VDATA interface to the file
    status = vfEnd( fileId )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not close HDF-VDATA interface for ' // &
                  TRIM( fileName )
       CALL ERROR_STOP( message, 'vdClose' )
    ENDIF
        
    ! Close the HDF file
    status = hClose( fileId )

    ! Error check status
    IF ( status == FAIL ) THEN 
       message = 'ERROR: Could not close the HDF file ' // TRIM( fileName )
       CALL ERROR_STOP( message, 'vdClose' )
    ENDIF
    
  END SUBROUTINE vdClose

!-----------------------------------------------------------------------------
  
  SUBROUTINE vdOpenField( name, vdId )

    !=====================================================================
    ! Subroutine "vdOpenField" initializes the HDF-VDATA interface 
    ! for a given VDATA field w/in the HDF file. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) name : CHARACTER name of the VDATA field to be initialized  
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (2) vdId : INTEGER VDATA ID# of the field
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) vsfAtch : returns INTEGER value ( vdId   )
    ! (2) vsfDtch : returns INTEGER value ( status )
    ! (3) vsfGid  : returns INTEGER value ( status )
    ! (4) vsfInq  : returns INTEGER value ( status )
    !
    ! NOTES: 
    !=====================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER,          INTENT(OUT) :: vdId

    ! Local Variables
    LOGICAL                       :: found 
    INTEGER                       :: intMode, nRec  
    INTEGER                       :: numType, status,  vDataRef
    CHARACTER(LEN=255)            :: list,    message, vdName

    ! External functions
    INTEGER, EXTERNAL             :: vsfAtch, vsfDtch, vsfGid, vsfInq

    !=====================================================================
    ! vdOpenField begins here!
    !=====================================================================

    ! Initialize
    found    = .FALSE.
    vDataRef = -1

    ! Loop thru file
    DO

       ! Look for HDF-VDATA field
       vDataRef = vsfGid( fileId, vDataRef )

       ! Exit if we are have come to EOF
       if ( vDataRef == FAIL ) EXIT

       ! Attach to this HDF-VDATA field 
       vdId = vsfAtch( fileId, vDataRef, 'r' )
       
       ! Get the name of this HDF-VDATA field
       status = vsfInq( vdId, nRec, intMode, list, numType, vdName )

       ! If the name of the field matches the name that we are 
       ! looking for, exit and return vdID to the calling routine
       IF ( TRIM( vdName ) == TRIM( name ) ) THEN
          found = .TRUE.
          EXIT
       ENDIF

       ! Otherwise, detach from this HDF-VDATA and loop again
       status = vsfDtch( vdId )
    ENDDO
   
    ! Error check if no files were found
    IF ( .not. found ) THEN 
       message = 'ERROR: Could not HDF-VDATA field ' // TRIM( name ) // &
                 ' in file ' // TRIM( saveFileName ) 
       CALL ERROR_STOP( message, 'vdOpenField' )
    ENDIF
    
  END SUBROUTINE vdOpenField

!-----------------------------------------------------------------------------

  SUBROUTINE vdCloseField( vdId )

    !=====================================================================
    ! Subroutine "vdCloseField" terminates the HDF-VDATA interface 
    ! for a given VDATA field w/in the HDF file. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) name : CHARACTER name of the VDATA field to be closed  
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) vsfDtch : returns INTEGER value 
    !
    ! NOTES: 
    !=====================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    INTEGER, INTENT(IN)    :: vdId

    ! Local Variables
    INTEGER                :: status
    CHARACTER(LEN=255)     :: message

    ! External functions
    INTEGER, EXTERNAL      :: vsfDtch

    !=====================================================================
    ! vdCloseField begins here!
    !=====================================================================

    ! Terminate VDATA interface
    status = vsfDtch( vdId )
   
    ! Error check 
    IF ( status == FAIL ) THEN 
       message = 'ERROR: Could not terminate HDF-VDATA interface!'
       CALL ERROR_STOP( message, 'vdCloseField' )
    ENDIF
  END SUBROUTINE vdCloseField

!-----------------------------------------------------------------------------

  SUBROUTINE vdPrintInfo

    !=====================================================================
    ! Subroutine "vdPrintInfo: obtains and prints information about 
    ! each HDF-VDATA field stored in the HDF file. (bmy, 7/3/03)
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) vsfAtch : returns INTEGER value ( vdId   )
    ! (2) vsfDtch : returns INTEGER value ( status )
    ! (3) vsfEx   : returns INTEGER value ( status )
    ! (4) vsfGid  : returns INTEGER value ( status )
    ! (5) vsfInq  : returns INTEGER value ( status )
    !=====================================================================

    ! Local variables
    INTEGER                  :: vdId,    vDataRef, numType
    INTEGER                  :: nRec,    intMode,  status
    CHARACTER(LEN=6)         :: numStr
    CHARACTER(LEN=255)       :: message, name,     list

    ! External functions
    INTEGER, EXTERNAL        :: vsfAtch, vsfDtch, vsfEx, vsfGid, vsfInq

    !=====================================================================
    ! vdPrintInfo begins here!
    !=====================================================================

    ! Start at beginning of file
    vDataRef = -1
  
    ! Loop thru file
    DO

       ! Look for VDATA Reference
       vDataRef = vsfGid( fileId, vDataRef )
       if ( VDataRef == FAIL ) EXIT

       ! Attach to this VDATA 
       vdId = vsfAtch( fileId, vDataRef, 'r' )

       ! If attach was successful, continue...
       IF ( status == SUCCEED ) THEN

          ! Get information about field #N from the HDF File
          status = vsfInq( vdId, nRec, intMode, list, numType, name )
       
          ! If status is successful, then print info
          IF ( status == SUCCEED ) THEN

             ! Pick number string
             SELECT CASE ( numType )
                CASE( 4 )
                   numStr = 'TYPE (XPLEX)'
                CASE( 8 )
                   numStr = 'TYPE (XPLEX)'
                CASE DEFAULT
                   numStr = 'N/A   '
             END SELECT

             ! Print info
             PRINT*, '--------------------------------------'
             PRINT*, 'HDF-VDATA # : ', vdId
             PRINT*, 'Name        : ', TRIM( name )
             PRINT*, '# records   : ', nRec
             PRINT*, 'Number Type : ', numStr
             PRINT*, 'Interlace   : ', intMode

          ENDIF
       ENDIF

       ! Detach from this VDATA and try again
       status = vsfDtch( vdId )
    ENDDO

  END SUBROUTINE vdPrintInfo

!-----------------------------------------------------------------------------

  SUBROUTINE vdGetFieldDim( vdId, vdDim )
    
    !===================================================================
    ! Subroutine vdGetFieldDim returns dimension information for
    ! a given HDF-VDATA field stored in the HDF file. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1) vdId  (INTEGER) : HDF-VD ID # for the given field
    !
    ! Arguments as Output:
    ! ------------------------------------------------------------------
    ! (2) vdDim (INTEGER) : Dimension (# of elements) of the VDATA
    !
    ! HDF-EOS library routines referenced:
    ! ------------------------------------------------------------------
    ! (1) vsQfNelt        : returns INTEGER ( status )
    !===================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)    :: vdId
    INTEGER, INTENT(OUT)   :: vdDim

    ! Local variables
    INTEGER                :: status
    CHARACTER(LEN=255)     :: message

    ! External functions
    INTEGER, EXTERNAL      :: vsQfNelt

    !===================================================================
    ! vdGetFieldSize begins here!
    !===================================================================

    ! Get information about field #N from the HDF File
    status = vsQfNelt( vdId, vdDim )
       
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not get the dimensions of HDF-VDATA field!'
       CALL ERROR_STOP( message, 'vdGetFieldDim' )
    ENDIF

  END SUBROUTINE vdGetFieldDim

!-----------------------------------------------------------------------------
  SUBROUTINE vdGetDataInt( vdId, nX, tData )

    !=====================================================================
    ! Subroutine vdGetDataInt reads a 1-D data array (INTEGER) from the
    ! HDF file.  The entire array will be returned. (zhe, 14/6/11)
    ! Added to standard code (zhej, dkh, 01/17/12, adj32_016) 
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) vdId  (INTEGER) : HDF-VD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (3 ) tData (INTEGER ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) vsfInq          : Returns INTEGER ( status            )
    ! (2 ) vsfRd           : Returns INTEGER ( # of records read )
    !=====================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER,  INTENT(IN)   :: vdId, nX
    INTEGER,  INTENT(OUT)  :: tdata(nX)

    ! Local variables
    INTEGER                :: intMode, nRec, numType, status
    CHARACTER(LEN=255)     :: message, name, list

    ! External functions
    INTEGER, EXTERNAL      :: vsfInq, vsfRd

    !===================================================================
    ! vdGetDataR8 begins here!
    !===================================================================
 
    ! Get information about the HDF-VDATA field
    status = vsfInq( vdId, nRec, intMode, list, numType, name )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not find HDF-VDATA field ' // &
                  TRIM( name ) // ' in file ' // TRIM( saveFileName )
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

    ! Also make sure the dimensions are compatible
    IF ( nX /= nRec ) THEN
       message = 'ERROR: nX does not match number of records in file!'
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

    ! Read the HDF-VDATA field from the file
    ! (status returns the # of records read)
    status = vsfRd( vdId, tData, nRec, intMode )
    
    ! Error check
    IF ( status <= 0 ) THEN
       message = 'ERROR: Did not read any records for HDF-VDATA field ' // &
                  TRIM( name )
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

  END SUBROUTINE vdGetDataInt

!-----------------------------------------------------------------------------

  SUBROUTINE vdGetDataR4( vdId, nX, tData )

    !=====================================================================
    ! Subroutine vdGetDataR4 reads a 1-D data array (TYPE (XPLEX)) from the
    ! HDF file.  The entire array will be returned. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) vdId  (INTEGER) : HDF-VD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (3 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) vsfInq          : Returns INTEGER ( status            )
    ! (2 ) vsfRd           : Returns INTEGER ( # of records read )
    !=====================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER,  INTENT(IN)   :: vdId, nX
    TYPE (XPLEX),   INTENT(OUT)  :: tdata(nX)

    ! Local variables
    INTEGER                :: intMode, nRec, numType, status
    CHARACTER(LEN=255)     :: message, name, list

    ! External functions
    INTEGER, EXTERNAL      :: vsfInq, vsfRd

    !===================================================================
    ! vdGetDataR8 begins here!
    !===================================================================
 
    ! Get information about the HDF-VDATA field
    status = vsfInq( vdId, nRec, intMode, list, numType, name )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not find HDF-VDATA field ' // &
                  TRIM( name ) // ' in file ' // TRIM( saveFileName )
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

    ! Also make sure the dimensions are compatible
    IF ( nX /= nRec ) THEN
       message = 'ERROR: nX does not match number of records in file!'
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

    ! Read the HDF-VDATA field from the file
    ! (status returns the # of records read)
    status = vsfRd( vdId, tData, nRec, intMode )
    
    ! Error check
    IF ( status <= 0 ) THEN
       message = 'ERROR: Did not read any records for HDF-VDATA field ' // &
                  TRIM( name )
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

  END SUBROUTINE vdGetDataR4

!-----------------------------------------------------------------------------

  SUBROUTINE vdGetDataR8( vdId, nX, tData )

    !=====================================================================
    ! Subroutine vdGetData reads a 1-D data array (TYPE (XPLEX)) from the
    ! HDF file.  The entire array will be returned. (bmy, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) vdId  (INTEGER) : HDF-VD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (3 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) vsfInq          : Returns INTEGER ( status            )
    ! (2 ) vsfRd           : Returns INTEGER ( # of records read )
    !=====================================================================

    ! References to F90 modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER,  INTENT(IN)  :: vdId, nX
    TYPE (XPLEX),   INTENT(OUT) :: tdata(nX)

    ! Local variables
    INTEGER               :: intMode, nRec, numType, status
    CHARACTER(LEN=255)    :: message, name, list

    ! External functions
    INTEGER, EXTERNAL     :: vsfInq, vsfRd

    !===================================================================
    ! vdGetDataR8 begins here!
    !===================================================================
 
    ! Get information about this field
    status = vsfInq( vdId, nRec, intMode, list, numType, name )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not find HDF-VDATA field ' // &
                  TRIM( name ) // ' in file ' // TRIM( saveFileName )
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

    ! Also make sure the dimensions are compatible
    IF ( nX /= nRec ) THEN
       message = 'ERROR: nX does not match number of records in file!'
       CALL ERROR_STOP( message, 'vdGetDataR4' )
    ENDIF

    ! Read the HDF-VDATA field from the file
    ! (status returns the # of records read)
    status = vsfRd( vdId, tData, nRec, intMode )
    
    ! Error check
    IF ( status <= 0 ) THEN
       message = 'ERROR: Did not read any records for HDF-VDATA field ' // &
                  TRIM( name )
       CALL ERROR_STOP( message, 'vdGetDataR8' )
    ENDIF

  END SUBROUTINE vdGetDataR8

!-----------------------------------------------------------------------------

  SUBROUTINE vdShift1d( nX, tData )

    !=====================================================================
    ! Subroutine vdShift1d shifts a 1-D data array by 180 degrees.  
    ! This is necessary since fvDAS data starts at 0 longitude, but 
    ! GEOS-CHEM needs the first box to be at -180 longitude. (bmy, 4/3/02)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) nX    (INTEGER) : Number of elements in the X-dimension
    ! (2) tData (TYPE (XPLEX) ) : Data array 
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (2) tData (TYPE (XPLEX) ) : Data array (shifted by 180 degrees)
    !=====================================================================

    ! Arguments
    INTEGER, INTENT(IN)    :: nX
    TYPE (XPLEX),  INTENT(INOUT) :: tData(nX)

    !===================================================================
    ! vdShift1d begins here!
    !===================================================================

    ! Shift the longitude dimension by nX/2 elements
    tData = CSHIFT( tdata, nX/2, 1 )

    END SUBROUTINE vdShift1d

!-----------------------------------------------------------------------------

  FUNCTION getTauFromDate( year, month, day ) RESULT( tau )

    !=====================================================================
    ! Function getTauFromDate returns the TAU value (hours since 0 GMT
    ! on Jan 1, 1985) at the beginning of the given date. (bmy, 12/12/05)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) year (INTEGER) : Current YYYY year  value
    ! (2 ) year (INTEGER) : Current MM   month value
    ! (3 ) year (INTEGER) : Current DD   day   value
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    INTEGER           :: year, month, day

    ! Local variables
    TYPE (XPLEX)            :: tau, jdToday

    ! Astronomical Julian Date at 0 GMT, 1 Jan 1985
    TYPE (XPLEX), PARAMETER :: JD85 = 2446066.5d0

    !=====================================================================
    ! getTauFromDate begins here!
    !=====================================================================

    ! Get today's astronomical Julian date
    jdToday = julDay( year, month, DCMPLX( day ) )
    
    ! Get Tau0 value
    tau     = ( jdToday - jd85 ) * 24d0 

  END FUNCTION getTauFromDate

!-----------------------------------------------------------------------------

  SUBROUTINE calDate( julDay, nymd, nhms )
    
    !=====================================================================
    ! Subroutine "calDate" converts an astronomical Julian day to 
    ! the NYMD (e.g. YYYYMMDD) and NHMS (i.e. HHMMSS) format.
    !
    ! Algorithm taken from "Practical Astronomy With Your Calculator",
    ! Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) julDay : TYPE (XPLEX)  : Astronomical julian day
    !
    ! Arguments as output:
    ! --------------------------------------------------------------------
    ! (1) nymd   : INTEGER : YYYYMMDD corresponding to JDAY
    ! (2) nhms   : INTEGER : HHMMSS   corresponding to JDAY
    !=====================================================================

    ! Arguments
    TYPE (XPLEX),  INTENT(IN)  :: julDay
    INTEGER, INTENT(OUT) :: nymd, nhms

    ! Local variables
    TYPE (XPLEX)               :: a, b, c, d, day, e, f 
    TYPE (XPLEX)               :: fDay, g, i, j, jd, m, y

    !=====================================================================
    ! "calDate begins here!
    ! See "Practical astronomy with your calculator", Peter Duffett-Smith
    ! 1992, for an explanation of the following algorithm.
    !=====================================================================
    jd = julDay + 0.5d0
    i  = INT( jd )
    f  = jd - INT( I )

    IF ( i > 2299160d0 ) THEN
       a = INT( ( I - 1867216.25d0 ) / 36524.25 )
       b = i + 1 + a - INT( a / 4 )
    ELSE
       b = i
    ENDIF

    c = b + 1524d0
    
    d = INT( ( c - 122.1d0 ) / 365.25d0 )

    e = INT( 365.25d0 * d )

    g = INT( ( c - e ) / 30.6001d0 )

    ! Day is the day number
    day  = c - e + f - INT( 30.6001d0 * g ) 

    ! fDay is the fractional day number
    fDay = day - int( day )
    
    ! M is the month number
    IF ( g < 13.5d0 ) THEN
       m = g - 1d0
    ELSE
       m = g - 13d0
    ENDIF

    ! Y is the year number
    IF ( m > 2.5d0 ) THEN
       y = d - 4716d0
    ELSE
       y = d - 4715d0
    ENDIF

    ! NYMD is YYYYMMDD
    nymd = ( INT( y ) * 10000 ) + ( INT( m ) * 100 ) + INT( day )
    
    ! NHMS is HHMMSS
    nhms = INT( fday * 24 ) * 10000 

  END SUBROUTINE calDate

!-----------------------------------------------------------------------------

  FUNCTION julDay( year, month, day ) RESULT( julianDay )

    !===================================================================
    ! Function JULDAY returns the astronomical Julian day.
    !
    ! Algorithm taken from "Practical Astronomy With Your Calculator",
    ! Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
    ! 
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1) YEAR  : (INTEGER) Current year
    ! (2) MONTH : (INTEGER) Current month
    ! (3) DAY   : (TYPE (XPLEX) ) Current day (can be fractional, e.g. 17.25)
    !
    ! NOTES:
    ! (2) JULDAY requires the external function MINT.F.
    ! 
    ! (3) JULDAY will compute the correct Julian day for any 
    !     BC or AD date.
    !      
    ! (4) For BC dates, subtract 1 from the year and append a minus 
    !     sign.  For example, 1 BC is 0, 2 BC is -1, etc.  This is 
    !     necessary for the algorithm.  
    !===================================================================

    ! Arguments
    INTEGER   :: year, month
    TYPE (XPLEX)    :: day,  julianDay

    ! Local variables
    INTEGER   :: year1, month1
    TYPE (XPLEX)    :: x1, a, b, c, d
    LOGICAL   :: isGregorian

    !===================================================================
    ! JULDAY begins here!
    !
    ! Follow algorithm from Peter Duffett-Smith (1992)
    !===================================================================

    ! Compute YEAR and MONTH1
    IF ( ( month == 1 ) .OR. ( month == 2 ) ) THEN
       year1  = year  - 1
       month1 = month + 12 
    ELSE
       year1  = year
       month1 = month
    ENDIF

    ! Compute the "A" term. 
    x1 = DCMPLX( year ) / 100.0d0
    a  = mint( x1 )

    ! The Gregorian calendar begins on 10 October 1582
    ! Any dates prior to this will be in the Julian calendar
    IF ( year > 1582 ) THEN
       isGregorian = .TRUE.
    ELSE
       IF ( ( year   == 1582 )  .AND. &
            ( month1 >= 10   )  .AND. &
            ( day    >= 15.0 ) ) THEN 
          isGregorian = .TRUE.
       ELSE
          isGregorian = .FALSE.
       ENDIF
    ENDIF
           
    ! Compute the "B" term according to Gregorian or Julian calendar
    IF ( isGregorian ) THEN
       b = 2.0d0 - a + mint( a / 4.0d0 )
    ELSE
       b = 0.0d0
    ENDIF

    ! Compute the "C" term for BC dates (YEAR1 <= 0 ) 
    ! or AD dates (YEAR1 > 0)
    IF ( year1 < 0 ) THEN
       x1 = ( 365.25d0 * year1 ) - 0.75d0
       c  = mint( x1 )
    ELSE
       x1 = 365.25d0 * year1
       c  = mint( x1 ) 
    ENDIF

    ! Compute the "D" term    
    x1 = 30.6001d0 * DCMPLX( month1 + 1 )
    d  = mint( x1 )

    
    ! Add the terms to get the Julian Day number 
    julianDay = b + c + d + day + 1720994.5d0

  END FUNCTION julDay

!-----------------------------------------------------------------------------

  FUNCTION mint( x ) RESULT ( value )

    !===================================================================
    ! Function MINT is defined as follows:
    ! 
    ! MINT = -INT( ABS( X ) ), X <  0
    ! MINT =  INT( ABS( X ) ), X >= 0
    ! 
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1) X : (TYPE (XPLEX)) Argument for the function MINT
    !
    ! NOTES:
    ! (1) MINT is primarily intended for use with routine JULDAY.
    !===================================================================

    ! Arguments
    TYPE (XPLEX), INTENT(IN) :: x
        
    ! Return value
    TYPE (XPLEX)             :: value

    !===================================================================
    ! MINT begins here!
    !===================================================================
    IF ( x < 0d0 ) THEN 
       value = -INT( ABS( x ) )        
    ELSE
       value =  INT( ABS( x ) )        
    ENDIF

  END FUNCTION MINT

!------------------------------------------------------------------------------

END MODULE HdfVdModule
