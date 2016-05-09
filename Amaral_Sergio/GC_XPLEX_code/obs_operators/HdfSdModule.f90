! $Id: HdfSdModule.f90,v 1.1 2009/06/18 19:53:07 daven Exp $
MODULE HdfSdModule

  !===========================================================================
  ! Module "HdfSdModule" contains variables and methods that are used to
  ! read data fields stored in HDF-SD format. (bmy, 4/26/02, 4/27/05)
  !
  ! In order to use HdfSdModule, you must first install the HDF-4 library 
  ! on your system.  You may download the library source code from:
  !
  !     http://hdf.ncsa.uiuc.edu/hdf4.html
  !
  ! There is also a good online tutorial about the HDF-SD interface at:
  !
  !     http://hdf.ncsa.uiuc.edu/training/HDFtraining/tutorial/sd/sds.html
  !
  ! Module Variables: 
  ! --------------------------------------------------------------------------
  ! (1 ) fileId            : ID number for the HDF file 
  ! (2 ) nDataSets         : # of data fields contained w/in a HDF file
  ! (3 ) nAttributes       : # of global atttributes contained w/in a HDF file
  ! (4 ) saveFileName      : Shadow variable used to store HDF file name
  !
  ! Module Methods:
  ! --------------------------------------------------------------------------
  ! (1 ) sdOpen            : Opens the HDF file
  ! (2 ) sdClose           : Closes the HDF file
  ! (3 ) sdName2Index      : Locates position of field w/in a HDF file by name
  ! (4 ) sdOpenField       : Opens a data field w/in a HDF file by index
  ! (5 ) sdOpenFieldByName : Opens a data field w/in a HDF file by name
  ! (6 ) sdCloseField      : Closes access to a data field w/in a HDF file
  ! (7 ) sdPrintInfo       : Prints information about fields w/in a HDF file
  ! (8 ) sdGetFieldDims    : Gets dimensions of a given field w/in a HDF file
  ! (9 ) sdGetData1dI4     : Reads  a 1-D INTEGER data field from the HDF file
  ! (10) sdGetData1d       : Reads  a 1-D TYPE (XPLEX)  data field from the HDF file
  ! (11) sdGetData2d       : Reads  a 2-D TYPE (XPLEX)  data field from the HDF file
  ! (12) sdGetData3d       : Reads  a 3-D TYPE (XPLEX)  data field from the HDF file
  ! (13) sdGetData4d       : Reads  a 4-D TYPE (XPLEX)  data field from the HDF file
  ! (14) sdGetData1d1Time  : Reads one time value of a 2-D TYPE (XPLEX) data field
  ! (15) sdGetData2d1Time  : Reads one time value of a 2-D TYPE (XPLEX) data field
  ! (16) sdGetData2d1TimeR8: Reads one time value of a 2-D TYPE (XPLEX) data fiedl
  ! (17) sdGetData3d1Time  : Reads one time value of a 3-D TYPE (XPLEX) data field
  ! (18) sdGetData3d1TimeR8: Reads one time value of a 3-D TYPE (XPLEX) data fiedl
  ! (19) sdShift1d         : Shifts a 1-D TYPE (XPLEX)  data field by 180 degrees
  ! (20) sdShift2d         : Shifts a 2-D TYPE (XPLEX)  data field by 180 degrees
  ! (21) sdShift3d         : Shifts a 3-D TYPE (XPLEX)  data field by 180 degrees
  ! (22) sdShift4d         : Shifts a 4-D TYPE (XPLEX)  data field by 180 degrees
  ! (23) sdGetMaxDims      : Returns the value of MAX_DIMS to outside routines
  !
  ! Module Interfaces: 
  ! --------------------------------------------------------------------------
  ! (1 ) sdGetData         : sdGetData1d, sdGetData1d_i4, sdGetData2d, 
  !                          sdGetData3d, sdGetData4d
  ! (2 ) sdGetData1        : sdGetData1d1Time, sdGetData2d1Time, 
  !                          sdGetData2d1TimeR8, sdGetData3d1Time,
  !                          sdGetData3d1TimeR8
  ! (3 ) sdShift           : sdShift1d, sdShift2d, sdShift3d, sdShift4d 
  !
  ! NOTES:
  ! (1 ) HdfSdModule is designed to only have one HDF file open at a time.  
  !       Once you have opened a file, you may attach/detach from as many
  !       individual data fields as you want.  (bmy, 4/9/02)
  ! (2 ) Added routines sdGetData1d1Time, sdGetData2d1Time, and 
  !       sdGetData3d1Time to read an array for only one time value
  !       from the HDF file, instead of reading the whole array (bmy, 4/25/02)
  ! (3 ) Declared internal routines and variables PRIVATE (bmy, 7/19/02)
  ! (4 ) Added routines sdGetData2d1TimeR8 and sdGetData3d1TimeR8 in order
  !       to read TYPE (XPLEX) data from the HDF file.  Also changed the name of
  !       sdGetData1d_i4 to sdGetData1dI4. (bmy, 7/19/02)
  ! (5 ) Minor updates.  Improved documentation and error/warning messages.
  !       (bmy, 7/3/03)
  ! (6 ) Modified for inclusion into GEOS-CHEM (bmy, 4/27/05)
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

  ! ... except these routines
  PUBLIC :: sdClose           
  PUBLIC :: sdCloseField      
  PUBLIC :: sdGetMaxDims 
  PUBLIC :: sdGetData     
  PUBLIC :: sdGetData1
  PUBLIC :: sdGetFieldDims    
  PUBLIC :: sdName2Index     
  PUBLIC :: sdOpen           
  PUBLIC :: sdOpenField       
  PUBLIC :: sdOpenFieldByName 
  PUBLIC :: sdPrintInfo       
  PUBLIC :: sdShift
  
  !=====================================================================
  ! MODULE VARIABLES
  !=====================================================================
  INTEGER                       :: fileId, nDataSets, nAttributes
  INTEGER,            PARAMETER :: MAX_DIMS = 4

  ! Shadow variable for file name
  CHARACTER(LEN=255), PRIVATE   :: saveFileName

  !=======================================================================
  ! Module interfaces: allow you to associate a name w/ several routines
  ! with different numbers of arguments or different argument types
  !=======================================================================
  INTERFACE sdGetData
     MODULE PROCEDURE sdGetData1dI4
     MODULE PROCEDURE sdGetData1d
     MODULE PROCEDURE sdGetData2d
     MODULE PROCEDURE sdGetData3d
     MODULE PROCEDURE sdGetData4d
  END INTERFACE 
  
  INTERFACE sdGetData1
     MODULE PROCEDURE sdGetData1d1Time
     MODULE PROCEDURE sdGetData2d1Time
     MODULE PROCEDURE sdGetData3d1Time
     MODULE PROCEDURE sdGetData2d1TimeR8
     MODULE PROCEDURE sdGetData3d1TimeR8
  END INTERFACE 

  INTERFACE sdShift
     MODULE PROCEDURE sdShift1d
     MODULE PROCEDURE sdShift2d
     MODULE PROCEDURE sdShift3d
     MODULE PROCEDURE sdShift4d 
  END INTERFACE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE sdOpen( fileName )

    !=====================================================================
    ! Subroutine "sdOpenFile" opens an HDF file and initializes the 
    ! scientific datasaet (HDF-SD) interface. (bmy, 4/3/02)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fileName : CHARACTER name of the HDF-EOS file to be opened
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) sfStart  : returns INTEGER value ( fileId )
    !
    ! NOTES: 
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    CHARACTER(LEN=*)       :: fileName

    ! Local Variables
    CHARACTER(LEN=255)     :: message

    ! External functions
    INTEGER, EXTERNAL :: sfStart

    !=====================================================================
    ! sdOpen begins here!
    !=====================================================================

    ! Save file name to a private shadow variable
    saveFileName = TRIM( fileName )

    ! Open the HDF file
    fileId = sfStart( TRIM( fileName ), DFACC_RDONLY )
  
    ! Error check fileId
    IF ( fileId == FAIL ) THEN 
       message = 'ERROR: Could not open the HDF file ' // TRIM( fileName )
       CALL ERROR_STOP( message, 'sdOpen' )
    ENDIF
    
  END SUBROUTINE sdOpen

!------------------------------------------------------------------------------
  
  SUBROUTINE sdClose( fileName )

    !=====================================================================
    ! Subroutine "sdClose" terminates the HDF Scientific Dataset 
    ! (HDF-SD) interface and closes the HDF file. (bmy, 4/3/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fileName : CHARACTER name of the HDF-EOS file to be opened
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) sfEnd    : takes INTEGER value ( fileId )
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP    

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: fileName

    ! Local variables
    INTEGER                 :: status
    CHARACTER(LEN=255)           :: message

    ! External functions
    INTEGER, EXTERNAL       :: sfEnd

    !=====================================================================
    ! sdClose begins here!
    !=====================================================================
        
    ! Close the HDF file
    status = sfEnd( fileId )

    ! Error check status
    IF ( status == FAIL ) THEN 
       message = 'ERROR: could not close the HDF file: ' // TRIM( fileName )
       CALL ERROR_STOP( message, 'sdClose' )
    ENDIF
    
  END SUBROUTINE sdClose

!-----------------------------------------------------------------------------

  SUBROUTINE sdName2Index( name, n ) 

    !===================================================================
    ! Subroutine sdName2Index finds out the number of a given data set
    ! within an HDF file given its name. (bmy, 4/3/02)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) name (CHARACTER)    : Name of the field to search for
    !
    ! Arguments as Output:
    ! ------------------------------------------------------------------
    ! (2 ) n    (INTEGER)      : Number of the field in the HDF file
    !
    ! External Functions:
    ! ------------------------------------------------------------------
    ! (1 ) sfN2Index (INTEGER) : Returns index based on the name
    !===================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER,          INTENT(OUT) :: n
    
    ! Local variables
    CHARACTER(LEN=255)            :: message

    ! External functions
    INTEGER,           EXTERNAL   :: sfN2Index

    !=====================================================================
    ! sdName2Index begins here!
    !=====================================================================

    ! Translate NAME to the HDF-EOS index number
    n = sfN2Index( fileId, TRIM( name ) )

    ! Make sure INDEX is valid
    IF ( n == FAIL ) THEN
       message = 'ERROR: '           // TRIM( name         ) // &
                 ' is not found in ' // TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdName2Index' ) 
    ENDIF

  END SUBROUTINE sdName2Index

!-----------------------------------------------------------------------------

  SUBROUTINE sdOpenField( n, sdId )

    !===================================================================
    ! Function sdOpenField initializes the HDF-SD interface for the Nth
    ! field in the HDF file.  The SD ID # is returned to the calling
    ! program. (bmy, 4/3/02, 7/3/03)    
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) n (INTEGER) : # of the scientific dataset in the HDF file
    !
    ! Arguments as Output:
    ! ------------------------------------------------------------------
    ! (2 ) sdId (INTEGER) : ID # for the corresponding SD
    !
    ! External Functions:
    ! ------------------------------------------------------------------
    ! (1 ) sfSelect (INTEGER) : Returns ID # for scientific dataset
    !===================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP    

    ! Arguments
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(OUT) :: sdId
    
    ! Local variables
    CHARACTER(LEN=255)   :: message

    ! External functions
    INTEGER              :: sfSelect

    !===================================================================
    ! sdOpenField begins here!
    !===================================================================

    ! Get the SD Index for the Nth dataset in the file
    sdId = sfSelect( fileId, n )
 
    ! Make sure data set ID is valid
    IF ( sdId == FAIL ) then 
       message = 'ERROR: Invalid ID # for HDF-SDATA field!'
       CALL ERROR_STOP( message, 'sdGetFieldId' )
    ENDIF
   
  END SUBROUTINE sdOpenField

!-----------------------------------------------------------------------------

  SUBROUTINE sdOpenFieldByName( name, sdId )

    !===================================================================
    ! Function sdOpenFieldbyName initializes the HDF-SD interface for
    ! given the field name.  The SD ID # is returned to the calling
    ! program. (bmy, 4/3/02, 7/3/03)    
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) name (CHARACTER) : name  of the SD field in the HDF file
    !
    ! Arguments as Output:
    ! ------------------------------------------------------------------
    ! (2 ) sdId (INTEGER) : ID # for the corresponding SD
    !
    ! External Functions:
    ! ------------------------------------------------------------------
    ! (1 ) sfSelect (INTEGER) : Returns ID # for scientific dataset
    !===================================================================
    
    ! Arguments
    !CHARACTER(LEN=255)   :: name
    CHARACTER(LEN=*)     :: name
    INTEGER, INTENT(OUT) :: sdId
    
    ! Local variables
    INTEGER              :: n
    CHARACTER(LEN=255)   :: message

    !===================================================================
    ! sdOpenFieldByName begins here!
    !===================================================================
    
    ! Convert name to index
    CALL sdName2Index( name, n )

    ! Open field w/ via the index
    CALL sdOpenField( n, sdId )
   
  END SUBROUTINE sdOpenFieldByName

!-----------------------------------------------------------------------------

  SUBROUTINE sdCloseField( sdId )

    !=====================================================================
    ! Subroutine sdCloseField terminates the HDF-SD interface for a 
    ! given field. (bmy, 4/3/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) N (INTEGER) : Number of the scientific dataset in the HDF file
    !
    ! External Functions:
    ! --------------------------------------------------------------------
    ! (1 ) sfSelect (INTEGER) : Returns ID # if successful or FAIL if not
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP
    
    ! Arguments
    INTEGER, INTENT(IN) :: sdId
    
    ! Local variables
    INTEGER             :: status
    CHARACTER(LEN=255)  :: message

    ! External functions
    INTEGER             :: sfEndAcc

    !=====================================================================
    ! sdCloseField begins here!
    !=====================================================================

    ! Terminate the ID # for this field
    status = sfEndAcc( sdId )

    ! Make sure data set ID is valid
    IF ( status == FAIL ) then 
       message = 'ERROR: Could not terminate HDF-SDATA interface!'
       CALL ERROR_STOP( message, 'sdCloseField' )
    ENDIF
   
  END SUBROUTINE sdCloseField

!-----------------------------------------------------------------------------

  SUBROUTINE sdPrintInfo

    !=====================================================================
    ! Subroutine sdPrintInfo obtains and prints information about each
    ! data field stored in the HDF file. (bmy, 4/3/02, 7/3/03)
    !
    ! HDF-EOS library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) sdFInfo : returns nAttributes, nDataSets 
    ! (2) sdGInfo : returns name, rank, dims, numType, nAttrs for fields
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Local variables
    INTEGER             :: sdId,    rank, dims(MAX_DIMS), numType
    INTEGER             :: n,       i,    nAttrs,         status
    CHARACTER(LEN=9)    :: numStr
    CHARACTER(LEN=255)  :: message, name

    ! External functions
    INTEGER, EXTERNAL   :: sfFInfo, sfGInfo
   
    !=====================================================================
    ! sdPrintInfo begins here!
    !
    ! Get global file information: # of data sets and attributes
    !=====================================================================
    status = sfFInfo( fileId, nDataSets, nAttributes )
    
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'Could not get info for ' // TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetInfo' )
    ENDIF
 
    !=====================================================================
    ! Get information about each field in the file
    !=====================================================================
    DO n = 0, nDataSets-1

       ! Initialize the HDF-SD interface for field # N
       CALL sdOpenField( n, sdId )

       ! Get information about field #N from the HDF File
       status = sfGInfo( sdId, name, rank, dims, numType, nAttrs )
       
       ! Print info if successful
       IF ( status == SUCCEED ) THEN

          ! Define string for number type
          SELECT CASE ( numType )
             CASE( 5 )
                numStr = 'TYPE (XPLEX)   '
             CASE( 6 )
                numStr = 'TYPE (XPLEX)   '
             CASE( 24 )
                numStr = 'INTEGER*4'
             CASE DEFAULT
                numStr = 'N/A      '
          END SELECT

          ! Print information
          PRINT*, '--------------------------------------'
          PRINT*, 'HDF-SDATA # : ', n
          PRINT*, 'Name        : ', TRIM( name )
          PRINT*, 'Rank        : ', rank
          PRINT*, 'Dimensions  : ', (dims(i), i=1,rank)
          PRINT*, 'Number Type : ', numStr
          PRINT*, 'Attributes  : ', nAttrs
       ENDIF

       ! Terminate the HDF-SD interface for field # N
       CALL sdCloseField( sdId )

    ENDDO
           
  END SUBROUTINE sdPrintInfo

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetFieldDims( sdId, nDims, dims )
    
    !===================================================================
    ! Subroutine sdGetFieldDims returns dimension information for
    ! the given field stored in the HDF file. (bmy, 4/3/02, 7/3/03)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1) sdId  (INTEGER) : HDF-SD ID # for the given field
    !
    ! Arguments as Output:
    ! ------------------------------------------------------------------
    ! (2) nDims (INTEGER) : Number of dimensions
    ! (3) dims  (INTEGER) : Array containing dimension information
    !
    ! HDF-EOS library routines referenced:
    ! ------------------------------------------------------------------
    ! (1) sdGInfo : returns name, rank, dims, numType, nAttrs for fields
    !===================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId
    INTEGER, INTENT(OUT) :: nDims, dims(MAX_DIMS)

    ! Local variables
    INTEGER              :: numType, nAttrs, status
    CHARACTER(LEN=255)   :: message, name

    ! External functions
    INTEGER, EXTERNAL    :: sfFInfo, sfGInfo

    !===================================================================
    ! sdGetFieldDims begins here!
    !===================================================================

    ! Zero out dimension array
    dims(:) = 0

    ! Get information about field #N from the HDF File
    status = sfGInfo( sdId, name, nDims, dims, numType, nAttrs )
       
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not locate HDF-SDATA field ' // &
                 TRIM( name ) // ' in ' // TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdPrintInfo' )
    ENDIF

  END SUBROUTINE sdGetFieldDims

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData1dI4( sdId, nX, tData )

    !=====================================================================
    ! Subroutine sdGetData1dI4 reads a 1-D data array (INTEGER) 
    ! from the HDF file.  The entire array will be returned. 
    ! (bmy, 7/19/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (3 ) tData (INTEGER) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, tData(nX)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message

    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData1d_i4 begins here!
    !===================================================================
 
    ! Read the data for the given field
    status = sfRData( sdId, (/0/), (/1/), (/nX/), tData )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData1d' )
    ENDIF

  END SUBROUTINE sdGetData1dI4

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData1d( sdId, nX, tData )

    !=====================================================================
    ! Subroutine sdGetData1d reads a 1-D data array from the HDF file.
    ! The entire array will be returned. (bmy, 4/3/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (3 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message

    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData1d begins here!
    !===================================================================
 
    ! Read the data for the given field
    status = sfRData( sdId, (/0/), (/1/), (/nX/), tData )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData1d' )
    ENDIF

  END SUBROUTINE sdGetData1d
  
!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData2d( sdId, nX, nY, tData )

    !=====================================================================
    ! Subroutine sdGetData2d reads a 2-D data array from the HDF file.
    ! The entire array will be returned. (bmy, 4/3/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (3 ) nY    (INTEGER) : Number of elements in the Y-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (4 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, nY
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message

    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData2d begins here!
    !===================================================================
 
    ! Read the data for the given field
    status = sfRData( sdId, (/0,0/), (/1,1/), (/nX,nY/), tData )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData2d' )
    ENDIF

  END SUBROUTINE sdGetData2d

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData3d( sdId, nX, nY, nZ, tData )

    !=====================================================================
    ! Subroutine sdGetData3d reads a 3-D data array from the HDF file.
    ! The entire array will be returned. (bmy, 4/3/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (3 ) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (4 ) nZ    (INTEGER) : Number of elements in the Z-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (5 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, nY, nZ
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY,nZ)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message

    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData3d begins here!
    !===================================================================
 
    ! Read the data for the given field
    status = sfRData( sdId, (/0,0,0/), (/1,1,1/), (/nX,nY,nZ/), tData )

    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData3d' )
    ENDIF

  END SUBROUTINE sdGetData3d

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData4d( sdId, nX, nY, nZ, nW, tData )

    !=====================================================================
    ! Subroutine sdGetData4d reads a 3-D data array from the HDF file.
    ! The entire array will be returned. (bmy, 4/9/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (3 ) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (4 ) nZ    (INTEGER) : Number of elements in the Z-dimension
    ! (5 ) nW    (INTEGER) : Number of elements in the W-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (6 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !
    ! NOTES:
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, nY, nZ, nW
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY,nZ,nW)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message
  
    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData4d begins here!
    !===================================================================

    ! Read the data for the given field
    status = sfRData( sdId, (/0,0,0,0/), (/1,1,1,1/), (/nX,nY,nZ,nW/), tData )
       
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData4d' )
    ENDIF

  END SUBROUTINE sdGetData4d

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData1d1Time( sdId, nTime, nX, tData )

    !=====================================================================
    ! Subroutine sdGetData2d1Time reads a 1-D data array for a single
    ! time.  We assume that the first dimension of the array in the
    ! HDF file is a spatial dimension, and the 2nd dimension is a time 
    ! dimension. (bmy, 4/26/02, 7/3/03)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nTime (INTEGER) : Time index (starting from 0) 
    ! (3 ) nX    (INTEGER) : Number of elements in the X-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (4 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !
    ! NOTES:
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, nTime
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message
  
    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData1d1Time begins here!
    !===================================================================

    ! Read the data for the given field
    status = sfRData( sdId, (/0,nTime/), (/1,1/), (/nX,1/), tData )
       
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData2d1Time' )
    ENDIF

  END SUBROUTINE sdGetData1d1Time

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData2d1Time( sdId, nTime, nX, nY, tData )

    !=====================================================================
    ! Subroutine sdGetData2d1Time reads a 2-D data array for a single
    ! time.  We assume that the first 2 dimensions of the array in the
    ! HDF file are spatial dimensions, and the 3rd dimension is a time 
    ! dimension. (bmy, 4/26/02, 7/3/03) 
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nTime (INTEGER) : Time index (starting from 0) 
    ! (3 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (4 ) nY    (INTEGER) : Number of elements in the Y-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (5 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !
    ! NOTES:
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, nY, nTime
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)    :: message
  
    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData2d1Time begins here!
    !===================================================================

    ! Read the data for the given field
    status = sfRData( sdId, (/0,0,nTime/), (/1,1,1/), (/nX,nY,1/), tData )
       
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData2d1Time' )
    ENDIF

  END SUBROUTINE sdGetData2d1Time

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData2d1TimeR8( sdId, nTime, nX, nY, tData )

    !=====================================================================
    ! Subroutine sdGetData2d1TimeR8 reads a 2-D data array (TYPE (XPLEX)) for  
    ! a single time.  We assume that the first 2 dimensions of the array 
    ! in the HDF file are spatial dimensions, and the 3rd dimension is 
    ! a time dimension. (bmy, 7/19/02, 7/3/03) 
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2 ) nTime (INTEGER) : Time index (starting from 0) 
    ! (3 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (4 ) nY    (INTEGER) : Number of elements in the Y-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (5 ) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1 ) sfRData         : Reads numeric data from the HDF file
    !
    ! NOTES:
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nX, nY, nTime
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message
  
    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData2d1Time begins here!
    !===================================================================

    ! Read the data for the given field
    status = sfRData( sdId, (/0,0,nTime/), (/1,1,1/), (/nX,nY,1/), tData )
       
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData2d1Time' )
    ENDIF

  END SUBROUTINE sdGetData2d1TimeR8

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData3d1Time( sdId, nTime, nX, nY, nZ, tData )

    !=====================================================================
    ! Subroutine sdGetData3d1Time reads a 3-D data array for a single
    ! time.  We assume that the first 3 dimensions of the array in the
    ! HDF file are spatial dimensions, and the 4th dimension is a time 
    ! dimension. (bmy, 4/26/02, 7/3/03) 
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2) nX    (INTEGER) : Number of elements in the X-dimension
    ! (3) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (4) nZ    (INTEGER) : Number of elements in the Z-dimension
    ! (5) nW    (INTEGER) : Number of elements in the W-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (6) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) sfRData         : Reads numeric data from the HDF file
    !
    ! NOTES:
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nTime, nX, nY, nZ
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY,nZ)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message
  
    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData3d1Time begins here!
    !===================================================================

    ! Read the data for the given field
    status = sfRData( sdId, (/0,0,0,nTime/), (/1,1,1,1/), &
                            (/nX,nY,nZ,1/),  tData )
     
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDAT field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData3d1Time' )
    ENDIF

  END SUBROUTINE sdGetData3d1Time

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetData3d1TimeR8( sdId, nTime, nX, nY, nZ, tData )

    !=====================================================================
    ! Subroutine sdGetData3d1TimeR8 reads a 3-D data array (TYPE (XPLEX)) for 
    ! a single time.  We assume that the first 3 dimensions of the array 
    ! in the HDF file are spatial dimensions, and the 4th dimension is 
    ! a time dimension. (bmy, 7/19/02, 7/3/03) 
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) sdId  (INTEGER) : HDF-SD # of the data field in the HDF file
    ! (2) nX    (INTEGER) : Number of elements in the X-dimension
    ! (3) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (4) nZ    (INTEGER) : Number of elements in the Z-dimension
    ! (5) nW    (INTEGER) : Number of elements in the W-dimension
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (6) tData (TYPE (XPLEX) ) : Data array
    ! 
    ! HDF library routines referenced:
    ! --------------------------------------------------------------------
    ! (1) sfRData         : Reads numeric data from the HDF file
    !
    ! NOTES:
    !=====================================================================

    ! References to F90 Modules
    USE ERROR_MOD, ONLY : ERROR_STOP

    ! Arguments
    INTEGER, INTENT(IN)  :: sdId, nTime, nX, nY, nZ
    TYPE (XPLEX),  INTENT(OUT) :: tData(nX,nY,nZ)

    ! Local variables
    INTEGER              :: status
    CHARACTER(LEN=255)   :: message
  
    ! External functions
    INTEGER, EXTERNAL    :: sfRData

    !===================================================================
    ! sdGetData3d1Time begins here!
    !===================================================================

    ! Read the data for the given field
    status = sfRData( sdId, (/0,0,0,nTime/), (/1,1,1,1/), &
                            (/nX,nY,nZ,1/),  tData )
     
    ! Error check
    IF ( status == FAIL ) THEN
       message = 'ERROR: Could not read HDF-SDATA field from ' // &
                  TRIM( saveFileName )
       CALL ERROR_STOP( message, 'sdGetData3d1Time' )
    ENDIF

  END SUBROUTINE sdGetData3d1TimeR8

!-----------------------------------------------------------------------------

  SUBROUTINE sdShift1d( nX, tData )

    !=====================================================================
    ! Subroutine sdShift1d shifts a 1-D data array by 180 degrees.  
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
    INTEGER,  INTENT(IN)    :: nX
    TYPE (XPLEX),   INTENT(INOUT) :: tData(nX)

    !===================================================================
    ! sdShift1d begins here!
    !===================================================================

    ! Shift the longitude dimension by nX/2 elements
    tData = CSHIFT( tdata, nX/2, 1 )

    END SUBROUTINE sdShift1d

!-----------------------------------------------------------------------------

  SUBROUTINE sdShift2d( nX, nY, tData )

    !=====================================================================
    ! Subroutine sdShift2d shifts a 2-D data array by 180 degrees.  
    ! This is necessary since fvDAS data starts at 0 longitude, but 
    ! GEOS-CHEM needs the first box to be at -180 longitude. (bmy, 4/3/02)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) nX    (INTEGER) : Number of elements in the X-dimension
    ! (2) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (3) tData (TYPE (XPLEX) ) : Data array 
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (3) tData (TYPE (XPLEX) ) : Data array (shifted by 180 degrees)
    !=====================================================================

    ! Arguments
    INTEGER,  INTENT(IN)    :: nX, nY
    TYPE (XPLEX),   INTENT(INOUT) :: tData(nX,nY)

    !===================================================================
    ! sdShift2d begins here!
    !===================================================================

    ! Shift the longitude dimension by nX/2 elements
    tData = CSHIFT( tdata, nX/2, 1 )

    END SUBROUTINE sdShift2d

!-----------------------------------------------------------------------------

  SUBROUTINE sdShift3d( nX, nY, nZ, tData )

    !=====================================================================
    ! Subroutine sdShift3d shifts a 3-D data array by 180 degrees.  
    ! This is necessary since fvDAS data starts at 0 longitude, but 
    ! GEOS-CHEM needs the first box to be at -180 longitude. (bmy, 4/3/02)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (2 ) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (3 ) nZ    (INTEGER) : Number of elements in the Z-dimension
    ! (4 ) tData (TYPE (XPLEX) ) : Data array 
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (5 ) tData (TYPE (XPLEX) ) : Data array (shifted by 180 degrees)
    !=====================================================================

    ! Arguments
    INTEGER,  INTENT(IN)    :: nX, nY, nZ
    TYPE (XPLEX),   INTENT(INOUT) :: tData(nX,nY,nZ)

    !===================================================================
    ! sdShift3d begins here!
    !===================================================================

    ! Shift the longitude dimension by nX/2 elements
    tData = CSHIFT( tdata, nX/2, 1 )

    END SUBROUTINE sdShift3d

!-----------------------------------------------------------------------------

  SUBROUTINE sdShift4d( nX, nY, nZ, nW, tData )

    !=====================================================================
    ! Subroutine sdShift4d shifts a 4-D data array by 180 degrees.  
    ! This is necessary since fvDAS data starts at 0 longitude, but 
    ! GEOS-CHEM needs the first box to be at -180 longitude. (bmy, 4/3/02)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1 ) nX    (INTEGER) : Number of elements in the X-dimension
    ! (2 ) nY    (INTEGER) : Number of elements in the Y-dimension
    ! (3 ) nZ    (INTEGER) : Number of elements in the Z-dimension
    ! (4 ) nW    (INTEGER) : Number of elements in the W-dimension
    ! (5 ) tData (TYPE (XPLEX) ) : Data array 
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------
    ! (5 ) tData (TYPE (XPLEX) ) : Data array (shifted by 180 degrees)
    !=====================================================================

    ! Arguments
    INTEGER,  INTENT(IN)    :: nX, nY, nZ, nW
    TYPE (XPLEX),   INTENT(INOUT) :: tData(nX,nY,nZ,nW)

    !===================================================================
    ! sdShift4d begins here!
    !===================================================================

    ! Shift the longitude dimension by nX/2 elements
    tData = CSHIFT( tData, nX/2, 1 )

    END SUBROUTINE sdShift4d

!-----------------------------------------------------------------------------

  SUBROUTINE sdGetMaxDims( maxDims )

    !===================================================================
    ! Subroutine sdGetMaxDims returns the value of MAX_DIMS to the 
    ! calling program.  This allows us to keep MAX_DIMS private.
    ! (bmy, 4/3/02)
    !
    ! Arguments as Output:
    ! ------------------------------------------------------------------
    ! (1 ) maxDims : Maximum # of dimensions for arrays in the HDF file
    !===================================================================

    ! Arguments
    INTEGER, INTENT(OUT) :: maxDims

    !===================================================================
    ! sdGetMaxDims begins here!
    !===================================================================
    maxDims = MAX_DIMS

  END SUBROUTINE sdGetMaxDims

!------------------------------------------------------------------------------

END MODULE HdfSdModule
