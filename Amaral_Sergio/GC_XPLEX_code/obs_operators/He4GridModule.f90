! $Id: He4GridModule.f90,v 1.1 2009/06/18 19:53:07 daven Exp $
MODULE He4GridModule

  !========================================================================
  ! Module "He4GridModule" contains variable declarations and routines
  ! that are used for reading data from HDF-EOS4 "Grid" structures.
  ! (bmy, 1/17/06, 4/10/08)
  !
  ! Module Variables:
  ! --------------------------------------------------------------------------
  ! fId          (INTEGER  ) : ID number for the HDF_EOS file 
  ! gId          (INTEGER  ) : ID number for the HDF_EOS grid 
  ! start3D      (INTEGER  ) : Starting index for 3-D arrays
  ! stride3D     (INTEGER  ) : Stride         for 3-D arrays
  ! edge3D       (INTEGER  ) : Ending   index for 3-D arrays
  ! start4D      (INTEGER  ) : Starting index for 4-D arrays
  ! stride4D     (INTEGER  ) : Stride         for 4-D arrays
  ! edge4D       (INTEGER  ) : Ending   index for 4-D arrays
  ! xStart       (INTEGER  ) : Starting index for the XDIM array 
  ! xStride      (INTEGER  ) : Stride         for the XDIM array 
  ! xEdge        (INTEGER  ) : Ending   index for the XDIM array
  ! yStart       (INTEGER  ) : Starting index for the YDIM array 
  ! yStride      (INTEGER  ) : Stride         for the YDIM array 
  ! yEdge        (INTEGER  ) : Ending   index for the YDIM array
  ! zStart       (INTEGER  ) : Starting index for the ZDIM array 
  ! zStride      (INTEGER  ) : Stride         for the ZDIM array 
  ! zEdge        (INTEGER  ) : Ending   index for the ZDIM array
  ! tStart       (INTEGER  ) : Starting index for the TIME, TAU arrays 
  ! tStride      (INTEGER  ) : Stride         for the TIME, TAU arrays 
  ! tEdge        (INTEGER  ) : Ending   index for the TIME, TAU arrays
  ! nFields      (INTEGER  ) : Number of fields in the current file
  ! fieldRank    (INTEGER  ) : Array of Argument for "gdinqfld" routine
  ! fieldType    (INTEGER  ) : Array of data types for each field in the file
  ! fieldName    (CHARACTER) : Array of names for each field in the file
  ! nAttrs       (INTEGER  ) : Number of attributes defined in the file
  ! attrName     (CHARACTER) : Array of attribute names
  ! attrValue    (TYPE (XPLEX)   ) : Array of attribute values
  ! time         (TYPE (XPLEX)   ) : Array of times (# of seconds since 1/1/1993)
  ! saveFileName (CHARACTER) : Shadow variable for the file name
  ! timefrom1993 (LOGICAL  ) : = T if TIME starts from 1/1/1993
  ! VERBOSE      (LOGICAL  ) : = T if we are printing info to the std output
  ! xDimSize     (INTEGER  ) : # of grid boxes in the X (longitude) dimension
  ! yDimSize     (INTEGER  ) : # of grid boxes in the Y (latitude ) dimension
  ! zDimSize     (INTEGER  ) : # of grid boxes in the Z (altitude ) dimension
  ! tDimSize     (INTEGER  ) : # of time intervals contained within this file
  ! xDim         (TYPE (XPLEX)   ) : Array of longitude centers
  ! yDim         (TYPE (XPLEX)   ) : Array of latitude  centers
  ! zdim         (TYPE (XPLEX)   ) : Array of alitude   centers
  ! nymd         (INTEGER  ) : Array of YYYYMMDD values -- date indices
  ! nhms         (INTEGER  ) : Array of HHMMSS   values -- hour indices
  !
  ! Module Methods:
  ! --------------------------------------------------------------------------
  ! (1 ) He4SetVerbose             : toggles information display on/off
  ! (2 ) He4GridOpen               : opens file & attaches to HDF-EOS4 grid
  ! (3 ) He4GridClose              : detaches from HDF-EOS4 grid & closes file
  ! (4 ) He4GridGetDimInfo         : gets dimensions of data fields in file
  ! (5 ) He4GridGetFldInfo         : gets info about each field
  ! (6 ) He4GridGetAttrInfo        : gets global attributes for HDF-EOS4 grid
  ! (7 ) He4GridReadAttrChar       : Reads CHARACTER attribute from grid
  ! (8 ) He4GridReadAttrI2         : Reads INTEGER*2 attribute from grid  
  ! (9 ) He4GridReadAttrI4         : Reads INTEGER*4 attribute from grid 
  ! (10) He4GridReadAttrR4         : Reads TYPE (XPLEX)    attribute from grid 
  ! (11) He4GridReadAttrR8         : Reads TYPE (XPLEX)    attribute from grid 
  ! (12) He4GridGetFillValue       : gets missing data "fill" value for fields
  ! (13) He4GridReadData3D         : reads a 3-D data block from the file
  ! (14) He4GridReadData4D         : reads a 4-D data block from the file
  ! (15) He4GridReadX              : gets the longitudes  (X) for the grid
  ! (16) He4GridReadY              : gets the latitudes   (Y) for the grid
  ! (17) He4GridReadZ              : gets the altitudes   (Z) for the grid    
  ! (18) He4GridReadT              : gets the time values (T) for the grid
  ! (19) He4GetNymdNhms            : converts T to NYMD, NHMS
  ! (20) He4CleanUpIndexFields     : deallocates index arrays 
  ! (21) makeCharArrayFromCharList : separates a string into a string array
  ! (22) calDate                   : converts Julian day to NYMD, NHMS
  ! (24) julDay                    : converts Year/month/day to Julian day
  ! (25) mint                      : function required by routine julDay
  !
  ! Module Interfaces:
  ! --------------------------------------------------------------------------
  ! (1 ) He4ReadGridAttr -- overloads these routines
  !      (a) He4GridReadAttrChar
  !      (b) He4GridReadAttrI2  
  !      (c) He4GridReadAttrI4
  !      (d) He4GridReadAttrR4
  !      (e) He4GridReadAttrR8
  !
  ! (2 ) He4GridReadData overloads the following routines
  !       (a) He4GridReadData3D
  !       (b) He4GridReadData4D
  !
  ! NOTES:
  ! (1 ) Updated for more consistency
  ! (2 ) Now declare "makeCharArrayFromCharList" public. (bmy, 11/8/06)
  ! (3 ) Updated comments. TYPEARRAY is now a global variable. (bmy, 8/14/07)
  ! (4 ) Added interface to read attribute data (bmy, 4/10/08)
  !===========================================================================

  ! References to F90 modules
  USE He4ErrorModule
  USE He4IncludeModule
 
  ! Force explicit data types
  USE MYTYPE
  USE COMPLEXIFY
  IMPLICIT NONE

  !---------------------------------------------------------------------
  ! PUBLIC / PRIVATE declarations
  !---------------------------------------------------------------------

  ! Make everything PRIVATE ...
  PRIVATE 

  ! ... and these routines
  PUBLIC :: xDimSize
  PUBLIC :: yDimSize
  PUBLIC :: zDimSize 
  PUBLIC :: tDimSize
  PUBLIC :: nymd
  PUBLIC :: nhms
  PUBLIC :: xDim
  PUBLIC :: yDim 
  PUBLIC :: zDim

  ! ... and these routines
  PUBLIC :: He4SetVerbose                
  PUBLIC :: He4GridOpen                  
  PUBLIC :: He4GridClose                 
  PUBLIC :: He4GridGetDimInfo      
  PUBLIC :: He4GridGetFldInfo          
  PUBLIC :: He4GridGetAttrInfo
  PUBLIC :: He4GridReadAttr
  PUBLIC :: He4GridGetFillValue          
  PUBLIC :: He4GridReadData          
  PUBLIC :: He4GridReadX                 
  PUBLIC :: He4GridReadY                 
  PUBLIC :: He4GridReadZ                     
  PUBLIC :: He4GridReadT                 
  PUBLIC :: He4GetNymdNhms               
  PUBLIC :: He4CleanUpIndexFields      
  PUBLIC :: makeCharArrayFromCharList

  !------------------------------------------------------------------------
  ! MODULE VARIABLES
  !------------------------------------------------------------------------

  ! Switch for printing output to the screen
  LOGICAL                     :: VERBOSE

  ! ID's
  INTEGER                     :: fId, gId

  ! Data extent
  INTEGER                     :: start3D(3), stride3D(3), edge3D(3)
  INTEGER                     :: start4D(4), stride4D(4), edge4D(4) 
  INTEGER                     :: xStart(1),  xStride(1),  xEdge(1)
  INTEGER                     :: yStart(1),  yStride(1),  yEdge(1)
  INTEGER                     :: zStart(1),  zStride(1),  zEdge(1)
  INTEGER                     :: tStart(1),  tStride(1),  tEdge(1)

  ! Fields
  INTEGER                     :: nFields
  INTEGER                     :: fieldRank(HE4_MAX_FLDS)
  INTEGER                     :: fieldType(HE4_MAX_FLDS)   
  CHARACTER(LEN=HE4_MAX_CHAR) :: fieldName(HE4_MAX_FLDS)

  ! Attributes
  INTEGER                     :: nAttrs
  TYPE (XPLEX)                      :: attrValue(HE4_MAX_ATRS)   
  CHARACTER(LEN=HE4_MAX_CHAR) :: attrName(HE4_MAX_ATRS)

  ! Variables for timing
  LOGICAL                     :: timeFrom1993 = .TRUE.
  TYPE (XPLEX),  ALLOCATABLE        :: time(:)

  ! Shadow variable for file name
  CHARACTER(LEN=HE4_MAX_CHAR) :: saveFileName

  ! Index arrays for grid
  INTEGER                     :: xDimSize, yDimSize, zDimSize, tDimSize
  INTEGER, ALLOCATABLE        :: nymd(:), nhms(:)
  TYPE (XPLEX),  ALLOCATABLE        :: xDim(:), yDim(:), zDim(:)

  ! Array for number types
  CHARACTER(LEN=10)           :: typeArray(10) =                         &
     (/ '         ', '         ', '         ', '         ', 'TYPE (XPLEX)   ', &
        'TYPE (XPLEX)   ', '         ', '         ', '         ', '         ' /)

  !------------------------------------------------------------------------
  ! MODULE INTERFACES
  !------------------------------------------------------------------------
  INTERFACE He4GridReadData
     MODULE PROCEDURE He4GridReadData3D
     MODULE PROCEDURE He4GridReadData4D
  END INTERFACE 
  
  INTERFACE He4GridReadAttr
     MODULE PROCEDURE He4GridReadAttrChar      
     MODULE PROCEDURE He4GridReadAttrI2           
     MODULE PROCEDURE He4GridReadAttrI4          
     MODULE PROCEDURE He4GridReadAttrR4          
     MODULE PROCEDURE He4GridReadAttrR8          
  END INTERFACE

  !------------------------------------------------------------------------
  ! MODULE ROUTINES
  !------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE He4SetVerbose( v )

    !======================================================================
    ! Subroutine setVerbose sets the value of module variable "verbose"
    ! which determines if information is echoed to the standard output
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) v (LOGICAL) : TRUE or FALSE value
    !
    ! NOTES:
    !======================================================================

    ! Arguments 
    LOGICAL, INTENT(IN) :: v

    ! Set the value of verbose
    VERBOSE = v

  END SUBROUTINE He4SetVerbose

!------------------------------------------------------------------------------

  SUBROUTINE He4GridOpen( fileName )

    !======================================================================
    ! Subroutine "gridOpen" opens the HDF_EOS file and attaches to the
    ! grid structure contained in the file. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fileName (CHARACTER) : Name of the HDF-EOS file to be opened
    !
    ! NOTES: 
    ! (1) The DAO uses the generic name "EOSGRID" for the grid 
    !     structure in all products.
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*)            :: fileName

    ! Local Variables
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                     :: gdAttach, gdOpen

    !--------------------------
    ! He4GridOpen begins here!
    !--------------------------

    ! Save file name to a private shadow variable
    saveFileName = TRIM( fileName )

    ! Call HDF library routine "gdopen" to open the HDF file
    fId = gdOpen( TRIM( fileName ), DFACC_RDONLY )
  
    ! Error check fId
    IF ( fId == FAIL ) THEN 
       msg = 'ERROR: Could not open file ' // TRIM( fileName )
       loc = 'He4GridOpen ("He4GridModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Call HDF library routine "gdAttach" to attach to the 
    ! grid structure contained in the HDF-EOS file.
    gId = gdAttach( fId, 'EOSGRID' )  

    ! Error check gId
    IF ( gId == FAIL ) THEN
       msg = 'ERROR: Could not attach to grid structure!'
       loc = 'He4GridOpen ("He4GridModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
    
  END SUBROUTINE He4GridOpen

!------------------------------------------------------------------------------
  
  SUBROUTINE He4GridClose( fileName )

    !=====================================================================
    ! Subroutine He4GridClose detaches from the currently opened grid
    ! and closes the HDF-EOS file that contains the grid. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fileName (CHARACTER) : Name of the HDF-EOS4 file to be closed
    !
    ! NOTES:
    !=====================================================================
    
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: fileName

    ! Local variables
    INTEGER                      :: status
    CHARACTER(LEN=HE4_MAX_CHAR)  :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                      :: gdDetach, gdClose

    !---------------------------
    ! He4GridClose begins here!
    !---------------------------

    ! Call HDF library routine "gdDetach" to detach from the 
    ! grid structure in the HDF-EOS file.
    status = gdDetach( gId )

    ! Error check status
    IF ( status == FAIL ) THEN 
       msg = 'ERROR detaching from grid structure!'
       loc = 'He4GridClose ("He4GridModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
    
    ! Call HDF library routine "gdClose" to close the HDF-EOS file.
    status = gdClose( fId )

    ! Error check status
    IF ( status == FAIL ) THEN 
       msg = 'ERROR closing the file: ' // TRIM( fileName )
       loc = 'He4GridClose ("He4GridModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
    
  END SUBROUTINE He4GridClose

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridGetDimInfo

    !=====================================================================
    ! Subroutine He4GridGetDimInfo obtains information about each 
    ! dimension of the grid contained in the HDF-EOS4 file. (bmy, 1/17/06)  
    !
    ! He4GridGetDimInfo also creates the various START, STRIDE, and EDGE 
    ! arrays needed to read data from the grid structure.
    !
    ! NOTES:
    !=====================================================================

    ! Local variables
    INTEGER                     :: dims(4),  status
    INTEGER                     :: xNumType, yNumType, zNumType, tNumType
    CHARACTER(LEN=HE4_MAX_CHAR) :: dimList
   
    ! HDF-EOS4 library routines
    INTEGER                     :: gdFldInfo

    !---------------------------------------------------------------------
    ! He4GridGetDimInfo begins here
    !
    ! Call HDF library routine "gdFldInfo" to get the grid dimensions.
    !
    ! xDimSize < 0 denotes missing X-dimension
    ! yDimSize < 0 denotes missing Y-dimension
    ! zDimSize < 0 denotes missing Z-dimension
    ! tDimSize < 0 denotes missing Time-dimension
    !---------------------------------------------------------------------
    status = gdFldInfo( gId, 'XDim',   dims, xDimSize, xNumType, dimList )
    status = gdFldInfo( gId, 'YDim',   dims, yDimSize, yNumType, dimList )
    status = gdFldInfo( gId, 'Height', dims, zDimSize, zNumType, dimList )
    status = gdFldInfo( gId, 'Time',   dims, tDimSize, tNumType, dimList )

    ! Create START, STRIDE, EDGE arrays for 3-D data fields (X,Y,Time)
    IF ( xDimSize > 0 .and. yDimSize > 0 .and. tDimSize > 0 ) THEN
       start3D  = 0
       stride3D = 1
       edge3D   = (/ xDimSize, yDimSize, tDimSize /)
    ENDIF

    ! Create START, STRIDE, EDGE arrays for 4-D data fields (X,Y,Z,Time)
    IF ( xDimSize > 0 .and. yDimSize > 0  .and. &
         zDimSize > 0 .and. tDimSize > 0 ) THEN
       start4D  = 0
       stride4D = 1
       edge4D   = (/ xDimSize, yDimSize, zDimSize, tDimSize /)
    endif

    ! Create START, STRIDE, EDGE arrays for index field xDim
    IF ( xDimSize > 0 ) THEN
       xStart  = 0
       xStride = 1
       xEdge   = xDimSize 
    ENDIF

    ! Create START, STRIDE, EDGE arrays for index field yDim
    IF ( yDimSize > 0 ) THEN
       yStart  = 0
       yStride = 1
       yEdge   = yDimSize
    ENDIF

    ! Create START, STRIDE, EDGE arrays for index field zDim
    IF ( zDimSize > 0 ) THEN
       zStart  = 0
       zStride = 1
       zEdge   = zDimSize 
    ENDIF

    ! Create START, STRIDE, EDGE arrays for index field Time
    IF ( tDimSize > 0 ) THEN
       tStart  = 0
       tStride = 1
       tEdge   = tDimSize 
    ENDIF

    ! Echo dimension information to the screen
    IF ( VERBOSE ) THEN
       WRITE( 6, '(a)' )
       WRITE( 6, '(a)' ) '      Index    Quantity   (Units)     Number of    Number'
       WRITE( 6, '(a)' ) '      Field                           Elements     Type'
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       
       IF ( xDimSize > 0 ) THEN
          WRITE( 6, '(6x, ''XDim     Longitude  (degrees) '', i8,7x,a10 )' ) &
               xDimSize, typeArray(xNumType) 
       ENDIF
    
       IF ( yDimSize > 0 ) THEN
          WRITE( 6, '(6x, ''YDim     Latitude   (degrees) '', i8,7x,a10 )' ) &
               yDimSize, typeArray(yNumType)
       ENDIF
       
       IF ( zDimSize > 0 ) THEN
          WRITE( 6, '(6x, ''Height   Altitude   (mb     ) '', i8,7x,a10 )' ) &
               zDimSize, typeArray(zNumType)
       ENDIF
       
       IF ( tDimSize > 0 ) THEN
          WRITE( 6, '(6x, ''Time     Time Index (seconds) '', i8,7x,a10 )' ) &
               tDimSize, typeArray(tNumType)
       ENDIF
    ENDIF
    
  END SUBROUTINE He4GridGetDimInfo

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridGetFldInfo
  
    !======================================================================
    ! Subroutine He4GridGetFldInfo obtains information about each of
    ! the fields stored in the HDF-EOS4 file.  Some information
    ! is echoed to the standard output.
    !
    ! NOTES:
    !======================================================================

    ! Local variables
    INTEGER                     :: i, dims(4), rank, numType, status
    TYPE (XPLEX)                      :: fillValue
    CHARACTER(LEN=HE4_MAX_CHAR) :: fieldList, dimList, msg, loc

    ! HDF-EOS4 library routines 
    INTEGER                     :: gdInqFlds, gdFldInfo

    !--------------------------------
    ! He4GridGetFldInfo begins here!
    !--------------------------------
    
    ! Call HDF library routine gdInqFlds to get information about
    ! each of the fields contained in the HDF-EOS file.
    nFields = gdInqFlds( gId, fieldList, fieldRank, fieldType )

    ! Call "makeCharArrayFromCharList" to create a character array
    ! using the comma-separated list of field names, FIELDLIST.
    CALL makeCharArrayFromCharList( fieldList, ',', fieldName )

    ! Write some header lines to the standard output
    IF ( VERBOSE ) THEN
       WRITE( 6, '(a)' )
       WRITE( 6, '(a)' )'      Field      Number      Fill          Dimensions'
       WRITE( 6, '(a)' )'      Name       Type        Value         of Field  '
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    ! Loop over each field
    DO i = 1, nFields

       ! Call HDF-EOS library routine "gdFldInfo" to read 
       ! information about each field's dimensions and type
       status = gdFldInfo( gId, TRIM( fieldName(i) ), &
                           dims, rank, numType, dimList )
 
       ! Get the "missing data" fill value for each field
       CALL He4GridGetFillValue( fieldName(i), fillValue )

       ! Echo information to the standard output       
       IF ( VERBOSE ) THEN
          WRITE( 6, '( i4,'') '', a10, 1x, a10, 1x, es9.2, 6x, a )' ) &
               i, fieldName(i), typeArray( fieldType(i) ), &
               fillValue, TRIM( dimList ) 
       ENDIF

    ENDDO
        
  END SUBROUTINE He4GridGetFldInfo

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridGetAttrInfo
  
    !======================================================================
    ! Subroutine He4GridGetAttrInfo obtains information about each of
    ! the global attributes for the HDF-EOS4 grid.  Some information
    ! is echoed to the standard output. (bmy, 1/17/06)
    !
    ! NOTES
    !======================================================================

    ! Local variables
    INTEGER                     :: i, status, strBufSize
    TYPE (XPLEX)                        :: attrValue
    CHARACTER(LEN=HE4_MAX_CHAR) :: attrList, message
    
    ! HDF-EOS4 library routines
    INTEGER                     :: gdInqAttrs, gdRdAttr

    !---------------------------------
    ! He4GridGetAttrInfo begins here!
    !---------------------------------

    ! Call HDF library routine gdInqAttrs to get information about
    ! each of the global attributes for the HDF-EOS grid
    nAttrs = gdInqAttrs( gId, attrList, strBufSize )

    ! Call "makeCharArrayFromCharList" to create a character array
    ! using the comma-separated list of attribute names, ATTRLIST.
    CALL makeCharArrayFromCharList( attrList, ',', attrName )

    ! Write some header lines to the standard output
    IF ( VERBOSE ) then
       WRITE( 6, '(a)' )
       WRITE( 6, '(a)' ) '      Attribute             Attribute'
       WRITE( 6, '(a)' ) '      Name                  Value'
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    
    ! Loop over all field names
    DO i = 1, nAttrs

       ! Read the value of each attribute
       status = gdRdAttr( gId, TRIM( attrName(i) ), attrValue )
       
       ! Echo information to the standard output       
       IF ( verbose ) then
          WRITE( 6, '( i4,'') '', a20, 1x, es13.6 )' ) &
               i, attrName(i), attrValue
       ENDIF
    ENDDO
    
  END SUBROUTINE He4GridGetAttrInfo

!------------------------------------------------------------------------------

 SUBROUTINE He4GridReadAttrChar( gId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4GridAttrChar returns a global attributes of type
    ! CHARACTER associated with the grid data structure. (bmy, 4/10/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) gId       (INTEGER  ) : HDF-EOS4 Grid ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (CHARACTER) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: gId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    CHARACTER(LEN=*), INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: GdRdAttr
    
    !-----------------------------------
    ! He4GridReadAttrChar begins here!
    !-----------------------------------

    ! Read attribute
    status = GdRdAttr( gId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4GridReadAttrChar

!------------------------------------------------------------------------------

  SUBROUTINE He4GridReadAttrI2( gId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4GridAttrI2 returns a global attributes of type
    ! INTEGER*2 associated with the grid data structure. (bmy, 4/10/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) gId       (INTEGER  ) : HDF-EOS4 grid ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (INTEGER*2) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: gId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    INTEGER*2,        INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: GdRdAttr
    
    !-----------------------------------
    ! He4GridReadAttrI2 begins here!
    !-----------------------------------

    ! Read attribute
    status = GdRdAttr( gId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4GridReadAttrI2

!------------------------------------------------------------------------------

  SUBROUTINE He4GridReadAttrI4( gId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4GridAttrI4 returns a global attributes of type
    ! INTEGER*4 associated with the grid data structure. (bmy, 4/10/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) gId       (INTEGER  ) : HDF-EOS4 grid ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (INTEGER*2) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: gId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    INTEGER,          INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: GdRdAttr
    
    !-----------------------------------
    ! He4GridReadAttrI4 begins here!
    !-----------------------------------

    ! Read attribute
    status = GdRdAttr( gId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4GridReadAttrI4

!------------------------------------------------------------------------------

  SUBROUTINE He4GridReadAttrR4( gId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4GridAttrR4 returns a global attributes of type
    ! TYPE (XPLEX) associated with the grid data structure. (bmy, 4/10/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) gId       (INTEGER  ) : HDF-EOS4 grid ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (TYPE (XPLEX)   ) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: gId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    TYPE (XPLEX),           INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: GdRdAttr
    
    !-----------------------------------
    ! He4SwathReadAttrR4 begins here!
    !-----------------------------------

    ! Read attribute
    status = GdRdAttr( gId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4GridReadAttrR4

!------------------------------------------------------------------------------

  SUBROUTINE He4GridReadAttrR8( gId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4GridAttrR8 returns a global attributes of type
    ! TYPE (XPLEX) associated with the grid data structure. (bmy, 4/10/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (TYPE (XPLEX)   ) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: gId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    TYPE (XPLEX),           INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: GdRdAttr
    
    !-----------------------------------
    ! He4GridReadAttrR8 begins here!
    !-----------------------------------

    ! Read attribute
    status = GdRdAttr( gId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4GridReadAttrR8

!------------------------------------------------------------------------------

  SUBROUTINE He4GridGetFillValue( fieldName, fillValue )
    
    !=====================================================================
    ! Subroutine He4GridGetFillValue reads the missing data "fill" value
    ! for a field contained in the HDF-EOS file. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fieldName (CHARACTER) : Name of the field to read in
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------    
    ! (2) fillValue (TYPE (XPLEX)   ) : Fill value for missing data
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: fieldName
    TYPE (XPLEX),           INTENT(OUT) :: fillValue

    ! Local variables
    INTEGER                       :: status

    ! HDF-EOS4 library routines
    INTEGER                       :: gdGetFill

    !----------------------------------
    ! He4GridGetFillValue begins here!
    !----------------------------------

    ! Call HDF library routine "gdrdfld" to read the data field
    status = gdGetFill( gId, TRIM( fieldName ), fillValue )

    ! Assign a large negative number to FILLVALUE if 
    IF ( status == FAIL ) fillValue = 0.0

  END SUBROUTINE He4GridGetFillValue

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridReadData3D( fldName, data3D )
    
    !=====================================================================
    ! Subroutine He4GridReadData3D" reads a 3-dimensional data field
    ! (X,Y,Time) from the HDF-EOS4 file. (bmy, 1/17/06) 
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fldName (CHARACTER) : Name of the field to read in
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------    
    ! (2) data3D  (TYPE (XPLEX)   ) : Data array (3 dimensions)
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: data3D(:,:,:)

    ! Local variables
    INTEGER                       :: status
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                       :: gdRdFld

    !--------------------------------
    ! He4GridReadData3D begins here!
    !--------------------------------

    ! Call HDF library routine "gdrdfld" to read the data field
    status = gdRdFld( gId, TRIM( fldName ), start3D, stride3D, edge3D, data3D )

    ! Error check
    IF ( status == FAIL ) THEN
       msg = 'ERROR reading data for field ' // TRIM( fldName )
       loc = 'He4GridReadData3D ("He4GridModule.f90")' 
       CALL He4ErrMsg( msg, loc )
    ENDIF

  END SUBROUTINE He4GridReadData3D

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridReadData4D( fldName, data4D )
    
    !=====================================================================
    ! Subroutine He4GridReadData4D reads a 4-dimensional data field
    ! (X,Y,Z,Time) from the HDF-EOS file. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) fldName (CHARACTER) : Name of the field to read in
    !
    ! Arguments as Output:
    ! --------------------------------------------------------------------    
    ! (1) data4D  (TYPE (XPLEX)   ) : Data array (4 dimensions)
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: data4D(:,:,:,:)

    ! Local variables
    INTEGER                       :: status
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                       :: gdRdFld

    !--------------------------------
    ! He4GridReadData4D begins here!
    !--------------------------------

    ! Call HDF library routine "gdrdfld" to read the data field
    status = gdRdFld( gId, TRIM( fldName ), start4D, stride4D, edge4D, data4D )

    ! Error check
    IF ( status == FAIL ) THEN
       msg = 'ERROR reading data for field ' // TRIM( fldName )
       loc = 'He4GridReadData4D ("He4GridReadModule.f90")' 
       CALL He4ErrMsg( msg, loc )
    ENDIF

  END SUBROUTINE He4GridReadData4D

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridReadX
    
    !=====================================================================
    ! Subroutine He4GridReadX reads XDIM, the index field for the
    ! X-dimension of the HDF-EOS4 grid structure. (bmy, 1/17/06)
    !
    ! NOTES:
    !=====================================================================

    ! Local variables    
    INTEGER                     :: status, as
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg,    loc

    ! HDF-EOS4 library routines
    INTEGER                     :: gdRdFld

    !---------------------------
    ! He4GridReadX begins here!
    !---------------------------

    ! Allocate the XDIM array of longitude centers (if necessary)
    IF ( xDimSize > 0 .and. .not. ALLOCATED( xDim ) ) THEN
       ALLOCATE( xDim( xEdge(1) ), stat=as )
       IF ( as /= 0 ) CALL He4AllocErr( 'xDim' )
    ELSE
       RETURN
    ENDIF

    ! Read XDIM, the vector of longitudes (in degrees), if present
    status = gdRdFld( gId, 'XDim', xStart, xStride, xEdge, xDim )
    
    ! Error check
    IF ( status == FAIL ) THEN
       msg = 'ERROR reading data for field xDim!' 
       loc = 'gridReadX ("He4GridReadModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
 
    ! Echo information to the standard output
    IF ( VERBOSE ) THEN
       WRITE( 6, '(a)'     )
       WRITE( 6, '(a)'     ) ' XDim: Longitude Centers (in degrees)'
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       WRITE( 6, '(9f8.2)' ) xDim
    ENDIF

  END SUBROUTINE He4GridReadX

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridReadY
    
    !=====================================================================
    ! Subroutine He4GridReadY reads YDIM, the index field for the
    ! Y-dimension of the HDF-EOS grid structure. (bmy, 1/17/06)
    ! 
    ! NOTES:
    !=====================================================================

    ! Local variables    
    INTEGER                     :: status, as
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg,    loc

    ! HDF-EOS4 library routines
    INTEGER                     :: gdRdFld

    !---------------------------
    ! He4GridReadY begins here!
    !---------------------------

    ! Allocate the YDIM array of longitude centers (if necessary)
    IF ( yDimSize > 0 .and. .not. ALLOCATED( yDim ) ) THEN
       ALLOCATE( yDim( yedge(1) ), stat=as )
       IF ( as /= 0 ) CALL He4AllocErr( 'yDim' ) 
    ELSE
       RETURN
    ENDIF

    ! Read YDIM, the vector of latitudes (in degrees), if present
    status = gdRdFld( gId, 'YDim', yStart, yStride, yEdge, yDim )
       
    ! Error check
    IF ( status == FAIL ) THEN
       msg = 'ERROR reading data for field yDim!' 
       loc = 'gridReadY ("He4GridModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
 
    ! Echo information to the standard output
    IF ( verbose ) THEN
       WRITE( 6, '(a)'     )
       WRITE( 6, '(a)'     ) ' YDim: Latitude Centers (in degrees)'
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       WRITE( 6, '(9f8.2)' ) yDim
    ENDIF

  END SUBROUTINE He4GridReadY

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridReadZ
    
    !=====================================================================
    ! Subroutine He4GridReadZ reads ZDIM, the index field for the
    ! Z-dimension of the HDF-EOS grid structure. (bmy, 1/17/06)
    !
    ! 
    ! NOTES:
    !=====================================================================

    ! Local variables    
    INTEGER                     :: status, as
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg,    loc

    ! HDF-EOS4 library routines
    INTEGER                     :: gdRdFld

    !---------------------------
    ! He4GridReadZ begins here!
    !---------------------------

    ! Allocate the ZDIM array of altitude centers (if necessary)
    IF ( zDimSize > 0 .and. .not. ALLOCATED( zDim ) ) THEN
       ALLOCATE( zDim( zedge(1) ), stat=as )
       IF ( as /= 0 ) CALL He4AllocErr( 'zDim' )
    ELSE
       RETURN
    ENDIF
    
    ! Read ZDIM, the vector of pressures, if present
    status = gdRdFld( gId, 'Height', zStart, zStride, zEdge, zDim )
       
    ! Error check
    IF ( status == FAIL ) THEN
       msg = 'ERROR reading data for field zDim!' 
       loc = 'gridReadZ ("He4GridReadModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
 
    ! Echo information to the standard output
    IF ( verbose ) THEN
       WRITE( 6, '(a)'     )
       WRITE( 6, '(a)'     ) ' ZDim: Altitude Indices'
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       WRITE( 6, '(9f8.2)' ) zDim
    ENDIF

  END SUBROUTINE He4GridReadZ

!-----------------------------------------------------------------------------

  SUBROUTINE He4GridReadT
    
    !=====================================================================
    ! Subroutine He4GridReadT reads TIME, the index field for the
    ! time-dimension of the HDF-EOS grid structure. (bmy, 1/17/06)
    !
    ! NOTES:
    !=====================================================================

    ! Local variables    
    INTEGER                     :: status, as, t
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg,    loc

    ! HDF-EOS4 library routines
    INTEGER                     :: gdRdFld

    !---------------------------
    ! He4gridReadT begins here!
    !---------------------------

    ! Allocate the time array, if it has not been previously allocated
    IF ( tDimSize > 0 .and. .not. ALLOCATED( time ) ) THEN
       ALLOCATE( time( tEdge(1) ), stat=as )
       IF ( as /= 0 ) CALL He4AllocErr( 'time' )
    ELSE
       RETURN
    ENDIF

    ! Read TIME, the vector of longitudes (in degrees), if present
    status = gdRdFld( gId, 'Time', tStart, tStride, tEdge, time )

    ! Error check TIME
    IF ( status == FAIL ) THEN
       msg = 'ERROR reading data for field Time!' 
       loc = 'He4GridReadT ("He4GridModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
    
  END SUBROUTINE He4GridReadT

!-----------------------------------------------------------------------------

  SUBROUTINE He4GetNymdNhms

    !=====================================================================
    ! Subroutine He4getNymdNhms converts the "Time" array into YYYYMMDD
    ! (NYMD) and HHMMSS (NHMS) values.  
    !
    ! NOTES: 
    ! (1) Some HDF-EOS files index time as seconds since 1993.  If this
    !     is the case, then set module variable TIMEFROM1993 = .TRUE.
    !
    ! (2) The HDF-EOS files for the GEOS-3/Terra assimilation index time
    !     from the starting date and time contained in the file name.
    !     To read these files, first set TIMEFROM1993 = .FALSE.
    !
    ! (3) Call routines JULDAY and CALDAT to compute the Year/Month/Day
    !     and Hour/Minute/Second.  These will account for time periods
    !     that straddle a month change.
    !
    ! (4) Now trim excess spaces from SAVEFILENAME before splitting 
    !     it up into segments.  Also error check NYMD and NHMS for
    !     negative values. (bmy, 9/21/00)
    !=====================================================================
    
    ! Local variables
    INTEGER                     :: t, as, year0, month0, day0, hour0
    TYPE (XPLEX)                      :: julianDay, julianDay0
    CHARACTER(LEN=HE4_MAX_CHAR) :: tmpStr, suffix(10)

    !-----------------------------
    ! He4GetNymdNhms begins here!
    !-----------------------------

    ! Get the time index array
    CALL He4GridReadT

    ! Allocate the nymd and nhms arrays
    IF ( tDimSize > 0 ) THEN

       ! NYMD array 
       IF ( .not. ALLOCATED( nymd ) ) THEN
          ALLOCATE( nymd( tEdge(1) ), stat=as )
          IF ( as /= 0 ) CALL He4AllocErr( 'nymd' )
       ENDIF

       ! NHMS array
       IF ( .not. ALLOCATED( nhms ) ) THEN
          ALLOCATE( nhms( tEdge(1) ), stat=as )
          IF ( as /= 0 ) CALL He4AllocErr( 'nhms' )
       ENDIF
    ELSE
       RETURN
    ENDIF

    !=====================================================================
    ! If TIME is measured from 1/1/1993, call "calDate" to 
    ! convert time to YYYYMMDD and HHMMSS values.
    !=====================================================================
    IF ( timeFrom1993 ) THEN 
       DO t = 1, tDimSize
          julianDay = 2448988.5d0 + ( time(t) / 86400d0 ) 
          CALL calDate( julianDay, nymd(t), nhms(t) )
       ENDDO

    !=====================================================================
    ! If TIME is NOT measured from 1/1/1993, then read YYYYMMDD from
    ! from the file name.  HDF-EOS containing GEOS-3/Terra data stick
    ! to the following naming convention:
    ! 
    ! (1) Assimilation files -- the suffix "tYYYYMMDD" indicates the
    !     starting date.  The starting time is always 0h GMT.
    !
    ! (2) Forecast files -- the suffix "bYYYYMMDDHH" indicates the
    !     starting date and GMT time.  
    ! 
    ! Therefore, extract the time/date info from the appropriate suffix.
    !=====================================================================
    ELSE
 
       ! Initialize
       suffix = ''

       ! Separate file name into individual segments
       ! Trim excess spaces from SAVEFILENAME (bmy, 9/21/00)
       CALL makeCharArrayfromCharList( TRIM( saveFileName ), '.', suffix )

       ! Initialize TMPSTR, for safety's sake
       tmpStr = ''

       ! Loop thru the file name segments from right to left
       DO t = 1, 10
          
          ! Save each segment into a temp string
          tmpStr = suffix(t)
             
          ! Skip null strings
          IF ( LEN_TRIM( tmpStr ) == 0 ) CYCLE
          
          ! Assimilation files have the date listed as "tYYYYMMDD"
          ! Extract starting year, month, day, and hour
          IF ( tmpStr(1:1) == 't' ) THEN
             IF ( LEN_TRIM( tmpStr ) == 9 ) THEN
                READ( tmpStr, '(1x,i4,i2,i2)' ) year0, month0, day0
                hour0 = 0
                EXIT
             ENDIF
          ENDIF

          ! Forecast files have the starting date/time listed as
          ! "bYYYYMMDDHH". Extract starting year, month, day, and hour
          IF ( tmpStr(1:1) == 'b' ) THEN
             IF ( LEN_TRIM( tmpStr ) == 11 ) THEN
                READ( tmpStr, '(1x,i4,i2,i2,i2)' ) year0, month0, day0, hour0
                EXIT
             ENDIF
          ENDIF
       ENDDO
       
       ! Compute starting Julian day
       julianDay0 = julDay( year0, month0, DCMPLX( day0 ) )

       ! Loop over all the elements of TIME
       DO t = 1, tDimSize 

          ! Compute the julian day corresponding to each element of TIME
          julianDay = julianDay0 + ( time(t)       / 1440.d0 ) + &
                                   ( DCMPLX( hour0 ) / 24d0    )

          ! Convert Julian day to NYMD, NHMS
          ! This will work for days that straddle the 1st of the month
          CALL calDate( julianDay, nymd(t), nhms(t) )
       ENDDO
    ENDIF

    ! Error check NYMD
    IF ( ANY( nymd < 0 ) ) THEN
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
       WRITE( 6, '(a)'   ) 'ERROR: NYMD is negative!'
       WRITE( 6, '(8i9)' ) nymd
       WRITE( 6, '(a)'   ) 'STOP in getNymdNhms (HdfModule)'
       STOP
    ENDIF

    ! Error check NHMS
    IF ( ANY( nhms < 0 ) ) THEN
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
       WRITE( 6, '(a)'   ) 'ERROR: NYMD is negative!'
       WRITE( 6, '(8i9)' ) nymd
       WRITE( 6, '(a)'   ) 'STOP in getNymdNhms (HdfModule)'
       STOP
    ENDIF

    ! Echo information values to the standard output
    IF ( VERBOSE ) THEN
       WRITE( 6, '(a)' )
       WRITE( 6, '(a)' ) ' Time      Time from           Date          Time'
    IF ( timeFrom1993 ) THEN
       WRITE( 6, '(a)' ) ' Index     1993 (s)       (YYYYMMDD)      (HHMMSS)'  
    ELSE
       WRITE( 6, '(a)' ) ' Index     start (s)      (YYYYMMDD)      (HHMMSS)'  
    ENDIF
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )

       DO t = 1, tDimSize
          WRITE( 6, '(i5,1x,f13.1,8x,i8.8,8x,i6.6)' ) &
               t, time(t), nymd(t), nhms(t)
       ENDDO
    ENDIF

  END SUBROUTINE He4GetNymdNhms

!-----------------------------------------------------------------------------

  SUBROUTINE He4CleanUpIndexFields

    !=====================================================================
    ! Subroutine He4CleanUpIndexFields deallocates the HDF-EOS index 
    ! fields xDim, YDim, ZDim, time, nymd, and nhms.  (bmy, 1/17/06)
    !=====================================================================
    IF ( ALLOCATED( xDim ) ) DEALLOCATE( xDim )
    IF ( ALLOCATED( yDim ) ) DEALLOCATE( yDim )
    IF ( ALLOCATED( zDim ) ) DEALLOCATE( zDim )
    IF ( ALLOCATED( time ) ) DEALLOCATE( time )   
    IF ( ALLOCATED( nymd ) ) DEALLOCATE( nymd )
    IF ( ALLOCATED( nhms ) ) DEALLOCATE( nhms )

  END SUBROUTINE He4CleanUpIndexFields

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

    !======================================================================
    ! Function JULDAY returns the astronomical Julian day. (bmy, 1/17/06)
    !
    ! Algorithm taken from "Practical Astronomy With Your Calculator",
    ! Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
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
    INTEGER, INTENT(IN) :: year, month
    TYPE (XPLEX),  INTENT(IN) :: day

    ! Local variables
    LOGICAL             :: isGregorian
    INTEGER             :: year1, month1
    TYPE (XPLEX)              :: x1, a, b, c, d, julianDay

    !======================================================================
    ! JULDAY begins here!
    !
    ! Follow algorithm from Peter Duffett-Smith (1992)
    !======================================================================

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

    !======================================================================
    ! Function MINT is defined as follows:
    ! 
    ! MINT = -INT( ABS( X ) ), X <  0
    ! MINT =  INT( ABS( X ) ), X >= 0
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) X : (TYPE (XPLEX)) Argument for the function MINT
    !
    ! NOTES:
    ! (1) MINT is primarily intended for use with routine JULDAY.
    !======================================================================

    ! Arguments
    TYPE (XPLEX), INTENT(IN) :: x
	
    ! Function value
    TYPE (XPLEX)             :: value

    !--------------------
    ! MINT begins here!
    !--------------------
    IF ( x < 0d0 ) THEN 
       value = -INT( ABS( x ) )        
    ELSE
       value =  INT( ABS( x ) )        
    ENDIF

  END FUNCTION MINT

!-----------------------------------------------------------------------------

  SUBROUTINE makeCharArrayFromCharList( list, separator, array )

    !======================================================================
    ! Subroutine makeCharArrayFromCharList takes a comma-separated word 
    ! list, and places each word into a separate element of a character 
    ! array. (bmy, 1/17/06, 11/8/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) list      (CHARACTER) : String with comma-separated words
    ! (2) separator (CHARACTER) : String for separator text
    !
    ! Arguments as output:
    ! ---------------------------------------------------------------------
    ! (3) array     (CHARACTER) : Array of substrings
    !
    ! NOTES:
    ! (1) Now set the output "array" argument to '' (bmy, 11/8/06)
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*),            INTENT(IN)  :: list
    CHARACTER(LEN=1           ), INTENT(IN)  :: separator
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: array(:)

    ! local variables
    INTEGER                                  :: P, N, ind(HE4_MAX_CHAR)
    CHARACTER(LEN=1)                         :: C

    !----------------------------------------
    ! makeCharArrayFromCharList begins here!
    !----------------------------------------

    ! Initialize
    N     = 1
    ind   = 0
    array = ''

    ! Find the positions of all the commas in LIST
    DO P = 1, LEN( list )

       ! Look at each character individually
       C = list(P:P)

       ! If a comma... 
       IF ( C == separator ) THEN 

          ! Increment comma
          N      = N + 1
          ind(N) = P 
       ENDIF
    ENDDO

    ! Add the position of the end of the string into IND
    ind(N+1) = LEN( list )

    ! Save text between the commas into ARRAY
    DO P = 1, N
       IF ( P == N ) THEN 
          array(P) = list( ind(P)+1:ind(P+1)   )
       ELSE
          array(P) = list( ind(P)+1:ind(P+1)-1 )
       ENDIF
    ENDDO

  END SUBROUTINE makeCharArrayFromCharList

!-----------------------------------------------------------------------------

END MODULE He4GridModule
