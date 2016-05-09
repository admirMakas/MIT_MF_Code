! $Id: He4SwathModule.f90,v 1.1 2009/06/18 19:53:07 daven Exp $
MODULE He4SwathModule

  !========================================================================
  ! Module He4SwathModule contains routines for reading data from swath
  ! data structures in HDF-EOS4 data files. (bmy, 1/17/06, 4/8/08)
  !
  ! Module Variables:
  ! -----------------------------------------------------------------------
  ! (1 ) VERBOSE       (LOGICAL  ) : Flag for toggling verbose output
  ! (2 ) dataTypeName  (CHARACTER) : Array w/ names of HDF-EOS4 data types
  ! (3 ) saveFileName  (CHARACTER) : Shadow variable for filename
  ! (4 ) saveSwathName (CHARACTER) : Shadow variable for swath name
  !
  ! Module Routines:
  ! -----------------------------------------------------------------------
  ! (1 ) He4VerboseOutput          : Toggles verbose output for file I/O
  ! (2 ) He4FileOpen               : Opens HDF4-EOS file; gets file ID #
  ! (3 ) He4FileClose              : Closes HDF4-EOS file
  ! (4 ) He4SwathAttach            : Attaches to swath; gets swath ID #
  ! (5 ) He4SwathDetach            : Detaches from swath
  ! (6 ) He4SwathDimInfo           : Gets dimension names, types, sizes
  ! (7 ) He4SwathGeoFldInfo        : Gets info about swath geoloc fields
  ! (8 ) He4SwathDataFldInfo       : Gets info about swath data fields
  ! (9 ) He4SwathFldInfo           : Gets info about an individual field
  ! (10) He4SwathFillValue         : Gets missing data fill values
  ! (11) He4SwathAttrs             : Gets attributes from HDF4-EOS swath
  ! (11) He4SwathReadAttrChar      : Reads CHARACTER attribute from swath
  ! (11) He4SwathReadAttrI2        : Reads INTEGER*2 attribute from swath   
  ! (11) He4SwathReadAttrI4        : Reads INTEGER*4 attribute from swath  
  ! (11) He4SwathReadAttrR4        : Reads TYPE (XPLEX)    attribute from swath  
  ! (11) He4SwathReadAttrR8        : Reads TYPE (XPLEX)    attribute from swath  
  ! (12) He4SwathReadData1dI2      : Reads 1-D INTEGER*2 data array 
  ! (12) He4SwathReadData1dI4      : Reads 1-D INTEGER*4 data array 
  ! (12) He4SwathReadData1dR4      : Reads 1-D TYPE (XPLEX)    data array 
  ! (13) He4SwathReadData1dR8      : Reads 1-D TYPE (XPLEX)    data array 
  ! (12) He4SwathReadData2dI2      : Reads 2-D INTEGER*2 data array 
  ! (12) He4SwathReadData2dI4      : Reads 2-D INTEGER*4 data array 
  ! (14) He4SwathReadData2dR4      : Reads 2-D TYPE (XPLEX)    data array 
  ! (15) He4SwathReadData2dR8      : Reads 2-D TYPE (XPLEX)    data array 
  ! (16) He4SwathReadData3dI2      : Reads 3-D INTEGER*2 data array 
  ! (17) He4SwathReadData3dI4      : Reads 3-D INTEGER*4 data array 
  ! (18) He4SwathReadData3dR4      : Reads 3-D TYPE (XPLEX)    data array 
  ! (19) He4SwathReadData3dR8      : Reads 3-D TYPE (XPLEX)    data array
  ! (20) He4SwathReadData4dI2      : Reads 3-D INTEGER*2 data array 
  ! (21) He4SwathReadData4dI4      : Reads 3-D INTEGER*4 data array 
  ! (22) He4SwathReadData4dR4      : Reads 3-D TYPE (XPLEX)    data array 
  ! (23) He4SwathReadData4dR8      : Reads 3-D TYPE (XPLEX)    data array
  ! (24) makeCharArrayFromCharList : Splits char list into char array
  ! (25) He4DataTypeName           : Returns data type name from type #
  !
  ! Module Interfaces:
  ! -----------------------------------------------------------------------
  ! (1 ) He4ReadSwathData -- overloads these routines
  !      (a) He4SwathReadData1dI2     
  !      (b) He4SwathReadData1dI4     
  !      (c) He4SwathReadData1dR4     
  !      (d) He4SwathReadData1dR8     
  !      (e) He4SwathReadData2dI2     
  !      (f) He4SwathReadData2dI4     
  !      (g) He4SwathReadData2dR4     
  !      (h) He4SwathReadData2dR8     
  !      (i) He4SwathReadData3dI2     
  !      (j) He4SwathReadData3dI4     
  !      (k) He4SwathReadData3dR4     
  !      (l) He4SwathReadData3dR8      
  !      (m) He4SwathReadData4dI2     
  !      (n) He4SwathReadData4dI4     
  !      (o) He4SwathReadData4dR4     
  !      (p) He4SwathReadData4dR8            
  !                               
  ! (2 ) He4ReadSwathAttr -- overloads these routines
  !      (a) He4SwathReadAttrChar
  !      (b) He4SwathReadAttrI2  
  !      (c) He4SwathReadAttrI4
  !      (d) He4SwathReadAttrR4
  !      (e) He4SwathReadAttrR8
  !
  ! Other Information:
  ! -----------------------------------------------------------------------
  ! (1 ) The data type HE4-INTEGER (represented by parameter HE4_INT in
  !       He4IncludeModule) is either INTEGER*4 or INTEGER*8 depending
  !       on platform.  However, most HDF4-EOS applications require
  !       INTEGER*4 dimension variables, etc.
  ! (2 ) You must select your machine type in the file "He4Define.h".
  !       This will automatically set parameter HE4_INT accordingly.
  ! (3 ) Data arrays which are passed to the HDF-EOS4 library function
  !       SwRdFld must be dimensioned with values of type HE4_INT.
  ! (4 ) Created interface for He4SwathReadAttr* functions.  This is 
  !       necessary to read attributes directly from the HDF4-EOS swath
  !       data structure. (bmy, 4/8/08)
  !
  ! References:
  ! -----------------------------------------------------------------------
  ! (1 ) http://hdf.ncsa.uiuc.edu/HDF4/
  !         -- HDF4 home page
  ! (2 ) http://newsroom.gsfc.nasa.gov/sdptoolkit/toolkit.html
  !         -- ECS toolkit home page (home of HDF-EOS4)
  !
  ! NOTES:
  ! (1 ) Added routines for 1dI2, IdI4, 2dI2, 2dI4, 3dI2, 3dI4 data types
  !       (bmy, 8/20/07)
  ! (2 ) Added routines for 4dI2, 4dI4, 4dR4, 4dR8 data types (bmy, 9/20/07)
  !========================================================================

  ! References to F90 modules
  USE He4ErrorModule
  USE He4IncludeModule

  ! Force explicit data types
  USE MYTYPE
  USE COMPLEXIFY
  IMPLICIT NONE

  !------------------------------------------------------------------------
  ! PRIVATE / PUBLIC DECLARATIONS
  !------------------------------------------------------------------------

  ! Make everything PRIVATE ...
  PRIVATE

  ! ... except these routines
  PUBLIC :: He4VerboseOutput
  PUBLIC :: He4FileOpen
  PUBLIC :: He4FileClose
  PUBLIC :: He4SwathAttach
  PUBLIC :: He4SwathDetach
  PUBLIC :: He4SwathDimInfo
  PUBLIC :: He4SwathGeoFldInfo
  PUBLIC :: He4SwathDataFldInfo
  PUBLIC :: He4SwathFldInfo
  PUBLIC :: He4SwathFillValue
  PUBLIC :: He4SwathAttrs
  PUBLIC :: He4SwathReadAttr
  PUBLIC :: He4SwathReadData

  !------------------------------------------------------------------------
  ! MODULE VARIABLES 
  !------------------------------------------------------------------------
  LOGICAL                     :: VERBOSE      = .FALSE.
  CHARACTER(LEN=HE4_MAX_CHAR) :: dataTypeName(57)
  CHARACTER(LEN=HE4_MAX_CHAR) :: saveFileName
  CHARACTER(LEN=HE4_MAX_CHAR) :: saveSwathName
  
  !------------------------------------------------------------------------
  ! MODULE INTERFACES
  !------------------------------------------------------------------------
  INTERFACE He4SwathReadAttr
     MODULE PROCEDURE He4SwathReadAttrChar      
     MODULE PROCEDURE He4SwathReadAttrI2           
     MODULE PROCEDURE He4SwathReadAttrI4          
     MODULE PROCEDURE He4SwathReadAttrR4          
     MODULE PROCEDURE He4SwathReadAttrR8          
  END INTERFACE

  INTERFACE He4SwathReadData
     MODULE PROCEDURE He4SwathReadData1dI2
     MODULE PROCEDURE He4SwathReadData1dI4
     MODULE PROCEDURE He4SwathReadData1dR4
     MODULE PROCEDURE He4SwathReadData1dR8
     MODULE PROCEDURE He4SwathReadData2dI2
     MODULE PROCEDURE He4SwathReadData2dI4
     MODULE PROCEDURE He4SwathReadData2dR4
     MODULE PROCEDURE He4SwathReadData2dR8
     MODULE PROCEDURE He4SwathReadData3dI2
     MODULE PROCEDURE He4SwathReadData3dI4
     MODULE PROCEDURE He4SwathReadData3dR4
     MODULE PROCEDURE He4SwathReadData3dR8
     MODULE PROCEDURE He4SwathReadData4dI2
     MODULE PROCEDURE He4SwathReadData4dI4
     MODULE PROCEDURE He4SwathReadData4dR4
     MODULE PROCEDURE He4SwathReadData4dR8
  END INTERFACE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE He4VerboseOutput( v )

    !======================================================================
    ! Subroutine He4VerboseOutput is used to trigger "extra" output from
    ! the routines in this module. (bmy, 1/17/06)
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

  END SUBROUTINE He4VerboseOutput

!------------------------------------------------------------------------------

  SUBROUTINE He4FileOpen( fileName, fId )

    !======================================================================
    ! Subroutine He4FileOpen opens an HDF-EOS4 file and returns the
    ! file Id number. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fileName (CHARACTER) : Name of HDF-EOS4 file to open
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) fId      (INTEGER  ) : HDF-EOS4 file ID number
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: fileName
    INTEGER,          INTENT(OUT) :: fId      

    ! Local variables
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwOpen

    !--------------------------
    ! He4FileOpen begins here!
    !--------------------------

    ! Store filename in a shadow variable
    saveFileName = fileName

    ! Open HDF-EOS4 file and get file ID #
    fId          = SwOpen( fileName, DFACC_RDONLY )

    ! Error check
    IF ( fId == FAILURE ) THEN
       msg = 'Error opening file ' // TRIM( saveFileName )
       loc = 'He4FileOpen ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
 
    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveFileName )
       WRITE( 6, 110 ) fId 
100    FORMAT( '===> HDF-EOS4 file name : "', a, '"' )
110    FORMAT( '===> HDF-EOS4 file ID   : ', i10     )
    ENDIF

  END SUBROUTINE He4FileOpen

!------------------------------------------------------------------------------

  SUBROUTINE He4FileClose( fId )

    !======================================================================
    ! Subroutine He4FileClose closes an HDF-EOS4 file. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fId (INTEGER) : HDF-EOS4 file ID number
    !
    ! NOTES:
    !======================================================================
    
    ! Arguments
    INTEGER, INTENT(IN)         :: fId      

    ! Local variables
    INTEGER                     :: status
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                     :: SwClose

    !---------------------------
    ! He4FileClose begins here!
    !---------------------------

    ! Get HDF-EOS4 file ID
    status = SwClose( fId )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error closing file ' // TRIM( savefileName )
       loc = 'He4FileClose ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveFileName )
100    FORMAT( '===> Closed file "', a, '"' )
    ENDIF

  END SUBROUTINE He4FileClose

!------------------------------------------------------------------------------

  SUBROUTINE He4SwathAttach( fId, swathName, sId )

    !======================================================================
    ! Subroutine He4SwathAttach attaches to an HDF-EOS4 swath data
    ! structure and returns the swath ID number. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fId       (INTEGER)   : HDF-EOS4 file ID (see He4FileOpen)
    ! (2) swathName (CHARACTER) : Name of HDF-EOS4 swath to attach to
    ! 
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) sId       (INTEGER)   : HDF-EOS4 swath ID number
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: fId
    CHARACTER(LEN=*), INTENT(IN)  :: swathName
    INTEGER,          INTENT(OUT) :: sId  

    ! Local variables
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwAttach

    !-----------------------------
    ! He4SwathAttach begins here!
    !-----------------------------

    ! Save swathname in a shadow variable
    saveSwathName = swathName

    ! Attach to swath
    sId           = SwAttach( fId, swathName )

    ! Error check
    IF ( sId == FAILURE ) THEN
       msg = 'Error attaching to swath ' // TRIM( saveSwathName )
       loc = 'He4SwathAttach ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF
  
    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveSwathName )
       WRITE( 6, 110 ) sId 
100    FORMAT( '===> HDF-EOS4 swath name: "', a, '"' )
110    FORMAT( '===> HDF-EOS4 swath ID  : ', i10     )  
    ENDIF

  END SUBROUTINE He4SwathAttach

!------------------------------------------------------------------------------

  SUBROUTINE He4SwathDetach( sId )

    !======================================================================
    ! Subroutine He4SwathDetach detaches from an HDF-EOS4 swath
    ! data structure. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId (INTEGER)   : HDF-EOS4 swath ID number (see He4SwathAttach)
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER, INTENT(IN)         :: sId  

    ! Local variables
    INTEGER                     :: status
    CHARACTER(LEN=HE4_MAX_CHAR) :: msg, loc

    ! HDF-EOS4 library routines
    INTEGER                     :: SwDetach

    !-----------------------------
    ! He4SwathDetach begins here!
    !-----------------------------

    ! Detach from swath
    status = SwDetach( sId )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error detaching from swath ' // TRIM( saveSwathName )
       loc = 'He4SwathDetach ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveSwathName )
100    FORMAT( '===> Detached from swath "', a, '"' )
    ENDIF

  END SUBROUTINE He4SwathDetach

!------------------------------------------------------------------------------

  SUBROUTINE He4SwathDimInfo( sId, nDims, dims, dimNames )

    !======================================================================
    ! Subroutine He4SwathDetach detaches from an HDF-EOS4 swath
    ! data structure. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId (INTEGER)   : HDF-EOS4 swath ID number (see He4SwathAttach)
    !
    ! NOTES:
    !======================================================================
    
    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId  
    INTEGER,                     INTENT(OUT) :: nDims  
    INTEGER,                     INTENT(OUT) :: dims(HE4_MAX_DIMS) 
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: dimNames(HE4_MAX_DIMS)
    
    ! Local variables
    INTEGER                                  :: N, C
    CHARACTER(LEN=HE4_MAX_CHAR)              :: msg, loc, dimList

    ! HDF-EOS4 library routines
    INTEGER                                  :: SwInqDims

    !------------------------------
    ! He4SwathDimInfo begins here!
    !------------------------------

    ! Initialize
    nDims       = 0
    dims(:)     = 0
    dimNames(:) = ''

    ! Get dimension info for this swath
    nDims       = SwInqDims( sId, dimList, dims )

    ! Make an array from the dimension list
    CALL makeCharArrayFromCharList( dimList, ',', dimNames )

! Comment out for now (bmy, 8/20/07)
!    ! NOTE: Sometimes every other element of DIMS is zero.  
!    ! I don't know why but we can just pack the array to be 
!    ! on the safe side. (bmy, 1/17/06)
!    C = 0
!    DO N = 1, 2*nDims+1 
!       IF ( dims(N) > 0  ) THEN
!          C       = C + 1
!          dims(C) = dims(N)
!          IF ( N > 1 ) dims(N) = 0
!       ENDIF
!    ENDDO
    
    ! Error check
    IF ( nDims <= 0 ) THEN
       msg = 'Error getting dim info from swath ' // TRIM( saveSwathName )
       loc = 'He4SwathDetach ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nDims

       DO N = 1, nDims 
          WRITE( 6, 110 ) TRIM( dimNames(N) ), dims(N)
       ENDDO

100    FORMAT( '===> There are ', i4, ' dimensions'      )
110    FORMAT( '===> ', a25,' is of size ', i10 )
    ENDIF

  END SUBROUTINE He4SwathDimInfo

!------------------------------------------------------------------------------

  SUBROUTINE He4SwathGeoFldInfo( sId, nGeo, geoRank, geoName, geoType )

    !======================================================================
    ! Subroutine He4SwathGeoFieldInfo obtains information about the 
    ! geolocation fields in the HDF-EOS4 swath data structure. 
    ! (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId     (INTEGER  ) : HDF-EOS4 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nGeo    (INTEGER  ) : Number of geolocation fields
    ! (3) geoRank (INTEGER  ) : Number of dimensions for each geoloc field
    ! (4) geoName (CHARACTER) : Name of each geolocation field
    ! (5) geoType (CHARACTER) : Data type of each geolocation field
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    INTEGER,                     INTENT(OUT) :: nGeo
    INTEGER,                     INTENT(OUT) :: geoRank(HE4_MAX_FLDS)
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: geoName(HE4_MAX_FLDS)
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: geoType(HE4_MAX_FLDS)

    ! Local variables
    INTEGER                                  :: N
    INTEGER                                  :: typeNum(HE4_MAX_FLDS)
    CHARACTER(LEN=HE4_MAX_CHAR)              :: msg, loc, geoList

    ! HDF-EOS4 library routines
    INTEGER                                  :: SwInqGFlds

    !---------------------------------
    ! He4SwathGeoFldInfo begins here!
    !---------------------------------

    ! Initialize
    nGeo       = 0
    geoRank(:) = 0
    geoName(:) = ''
    geoType(:) = ''

    ! Get number of geo fields and related info
    nGeo       = SwInqGFlds( sId, geoList, geoRank, typeNum )

    ! Error check
    IF ( nGeo <= 0 ) THEN
       msg = 'Error getting geolocation field information!'
       loc = 'He4SwathGeoFldInfo ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Separate list of field names into an array
    CALL makeCharArrayFromCharList( geoList, ',', geoName )

    ! Get HDF-EOS4 data type names for each data type number
    DO N = 1, nGeo
       geoType(N) = He4DataTypeName( typeNum(N) )
    ENDDO

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nGeo

       DO N = 1, nGeo
          WRITE( 6, 110 ) TRIM( geoName(N) ), geoRank(N), TRIM( geoType(N) )
       ENDDO

100    FORMAT( '===> There are ', i4, ' Geolocation Fields' )
110    FORMAT( '===> ', a25, ' has ', i4 , ' dimensions and is ', a )
    ENDIF

  END SUBROUTINE He4SwathGeoFldInfo
    
!------------------------------------------------------------------------------

  SUBROUTINE He4SwathDataFldInfo( sId, nData, dataRank, dataName, dataType )

    !======================================================================
    ! Subroutine He4SwathDataFieldInfo obtains information about the 
    ! data fields in the HDF-EOS4 swath data structure. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId      (INTEGER  ) : HDF-EOS4 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nData    (INTEGER  ) : Number of data fields
    ! (3) dataRank (INTEGER  ) : Number of dimensions for each data field
    ! (4) dataName (CHARACTER) : Name of each data field
    ! (5) dataType (CHARACTER) : Data type of each data field
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    INTEGER,                     INTENT(OUT) :: ndata
    INTEGER,                     INTENT(OUT) :: dataRank(HE4_MAX_FLDS)
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: dataName(HE4_MAX_FLDS)
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: dataType(HE4_MAX_FLDS)

    ! Local variables
    INTEGER                                  :: N
    INTEGER                                  :: typeNum(HE4_MAX_FLDS)
    CHARACTER(LEN=HE4_MAX_CHAR)              :: msg, loc, dataList

    ! HDF-EOS4 library routines
    INTEGER                                  :: SwInqDFlds

    !---------------------------------
    ! He4SwathGeoFldInfo begins here!
    !---------------------------------

    ! Initialize
    nData       = 0
    dataRank(:) = 0
    dataName(:) = ''
    dataType(:) = ''

    ! Get number of data fields and related info
    nData       = SwInqDFlds( sId, dataList, dataRank, typeNum )

    ! Error check
    IF ( nData <= 0 ) THEN
       msg = 'Error getting data field information!'
       loc = 'He4SwathDataFldInfo ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Separate list of field names into an array
    CALL makeCharArrayFromCharList( dataList, ',', dataName )

    ! Get HDF-EOS4 data type names for each data type number
    DO N = 1, nData
       dataType(N) = He4DataTypeName( typeNum(N) )
    ENDDO

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nData

       DO N = 1, nData
          WRITE( 6, 110 ) TRIM( dataName(N) ), dataRank(N), TRIM( dataType(N) )
       ENDDO

100    FORMAT( '===> There are ', i6, ' Data Fields' )
110    FORMAT( '===> ', a40, ' has ', i4 , ' dimensions and is ', a )
    ENDIF

  END SUBROUTINE He4SwathDataFldInfo
    
!------------------------------------------------------------------------------

  SUBROUTINE He4SwathFldInfo( sId, name, typeName, rank, dims, dimNames )
    
    !======================================================================
    ! Subroutine He4SwathFldInfo obtains information about a particular
    ! data field in the HDF-EOS4 swath data structure. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId      (INTEGER  ) : HDF-EOS4 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nData    (INTEGER  ) : Number of data fields
    ! (3) typeName (CHARACTER) : Name of the data type for this field
    ! (4) rank     (INTEGER  ) : Number of dimensions for each data field
    ! (4) dims     (INTEGER  ) : Integer containing field dimensions
    ! (5) dimNames (CHARACTER) : Array containing names of each dimension
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    CHARACTER(LEN=*),            INTENT(IN)  :: name
    INTEGER,                     INTENT(OUT) :: rank
    INTEGER,                     INTENT(OUT) :: dims(HE4_MAX_DIMS)
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: dimNames(HE4_MAX_DIMS)
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: typeName

    ! Local variables
    INTEGER                                  :: C,   N,   status, typeNum
    CHARACTER(LEN=HE4_MAX_CHAR)              :: msg, loc, dimList

    ! HDF-EOS4 library routines
    INTEGER                                  :: SwFldInfo

    !------------------------------
    ! He4SwathFldInfo begins here!
    !------------------------------

    ! Initialize
    rank     = 0
    dims     = 0
    dimNames = ''

    ! Get number of data fields and related info
    status = SwFldInfo( sId, name, rank, dims, typeNum, dimList )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error getting info about field ' // TRIM( name )
       loc = 'He4SwathFldInfo ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Separate list of dimension names into an array
    CALL makeCharArrayFromCharList( dimList, ',', dimNames )

    ! Get HDF-EOS4 data type name
    typeName = He4DataTypeName( typeNum )

! Comment out for now (bmy, 8/20/07)
!    ! NOTE: Sometimes every other element of DIMS is zero.  
!    ! I don't know why but we can just pack the array to be 
!    ! on the safe side. (bmy, 1/17/06)
!    C = 0
!    DO N = 1, HE4_MAX_DIMS
!       IF ( dims(N) > 0 ) THEN
!          C       = C + 1
!          dims(C) = dims(N)
!          IF ( N > 1 ) dims(N) = 0
!       ENDIF
!    ENDDO

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( name ), rank, TRIM( typeName )
       WRITE( 6, 110 ) dims(1:rank)
100    FORMAT( '===> ', a25 ' has ', i4 , ' dimensions and is ', a )
110    FORMAT( '===> ', 25x,' its dimensions are: ', 10i7 )
    ENDIF

  END SUBROUTINE He4SwathFldInfo

!------------------------------------------------------------------------------

  SUBROUTINE He4SwathFillValue( sId, dataName, dataFill )
    
    !======================================================================
    ! Subroutine He4SwathFillValue reads the missing data "fill" value
    ! for a field contained in the HDF-EOS4 swath data structure. 
    ! (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId      (INTEGER  ) : HDF-EOS4 swath ID 
    ! (2) dataName (CHARACTER) : Name of the field to get fill value for
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------   
    ! (3) dataFill (TYPE (XPLEX)   ) : Fill value for missing data
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(IN)  :: dataName
    TYPE (XPLEX),                      INTENT(OUT) :: dataFill

    ! Local variables
    INTEGER                                  :: status

    ! HDF-EOS4 library routines
    INTEGER                                  :: SwGetFill

    !--------------------------------
    ! He4SwathFillValue begins here!
    !--------------------------------

    ! Get the fill value
    status = SwGetFill( sId, TRIM( dataName ), dataFill )

    ! Set fill value to zero if it does not exist
    IF ( status == FAILURE ) dataFill = 0.0

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( dataName ), dataFill
100    FORMAT( '===> Fill value for ', a, ' is ', es13.6 )
    ENDIF

  END SUBROUTINE He4SwathFillValue
 
!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathAttrs( sId, nAttrs, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4SwathAttrs returns the global attributes associated
    ! with the swath data structure. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nAttrs    (INTEGER  ) : Number of swath attributes
    ! (3) attrName  (CHARACTER) : Array of attribute names
    ! (4) attrValue (CHARACTER) : Array of attribute values
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    INTEGER,                     INTENT(OUT) :: nAttrs
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: attrName(HE4_MAX_ATRS)
    TYPE (XPLEX),                      INTENT(OUT) :: attrValue(HE4_MAX_ATRS)

    ! Local variables
    INTEGER                                  :: N, status, strBufSize
    TYPE (XPLEX)                                   :: value
    CHARACTER(LEN=HE4_MAX_CHAR)              :: attrList
    
    ! HDF_EOS5 library routines
    INTEGER                                  :: SwInqAttrs
    INTEGER                                  :: SwRdAttr

    !----------------------------
    ! He4SwathAttrs begins here!
    !----------------------------

    ! Get list of attribute names
    nAttrs = SwInqAttrs( sId, attrList, strBufSize )

    ! Separate list into array
    CALL makeCharArrayFromCharList( attrList, ',', attrName )

    ! Get the data value for each attribute
    ! For each attribute
    DO N = 1, nAttrs
       status = SwRdAttr( sId, TRIM( attrName(N) ), attrValue )       
    ENDDO
    
    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nAttrs

       DO N = 1, nAttrs
          WRITE( 6, 110 ) TRIM( attrName(N) ), attrValue(N)
       ENDDO

100    FORMAT( '===> There are ', i6, ' Swath Attributes' )
110    FORMAT( '===> ', a20, ' has value ', es13.6 )
    ENDIF

  END SUBROUTINE He4SwathAttrs

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadAttrChar( sId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4SwathAttrChar returns a global attributes of type
    ! CHARACTER associated with the swath data structure. (bmy, 4/8/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (CHARACTER) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    CHARACTER(LEN=*), INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: SwRdAttr
    
    !-----------------------------------
    ! He4SwathReadAttrChar begins here!
    !-----------------------------------

    ! Read attribute
    status = SwRdAttr( sId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4SwathReadAttrChar

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadAttrI2( sId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4SwathAttrI2 returns a global attributes of type
    ! INTEGER*2 associated with the swath data structure. (bmy, 4/8/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (INTEGER*2) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    INTEGER*2,        INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: SwRdAttr
    
    !-----------------------------------
    ! He4SwathReadAttrI2 begins here!
    !-----------------------------------

    ! Read attribute
    status = SwRdAttr( sId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4SwathReadAttrI2

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadAttrI4( sId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4SwathAttrI4 returns a global attributes of type
    ! INTEGER*4 associated with the swath data structure. (bmy, 4/8/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (INTEGER*2) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    INTEGER,          INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: SwRdAttr
    
    !-----------------------------------
    ! He4SwathReadAttrI4 begins here!
    !-----------------------------------

    ! Read attribute
    status = SwRdAttr( sId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4SwathReadAttrI4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadAttrR4( sId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4SwathAttrR4 returns a global attributes of type
    ! TYPE (XPLEX) associated with the swath data structure. (bmy, 4/8/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (INTEGER*2) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    TYPE (XPLEX),           INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: SwRdAttr
    
    !-----------------------------------
    ! He4SwathReadAttrR4 begins here!
    !-----------------------------------

    ! Read attribute
    status = SwRdAttr( sId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4SwathReadAttrR4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadAttrR8( sId, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He4SwathAttrR8 returns a global attributes of type
    ! TYPE (XPLEX) associated with the swath data structure. (bmy, 4/8/08)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS4 swath ID number 
    ! (2) attrName  (CHARACTER) : Name of attribute to read from file
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) attrValue (INTEGER*2) : Value of attribute
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    CHARACTER(LEN=*), INTENT(IN)  :: attrName
    TYPE (XPLEX),           INTENT(OUT) :: attrValue

    ! Local variables
    INTEGER                       :: status
    
    ! HDF4-EOS library routines
    INTEGER                       :: SwRdAttr
    
    !-----------------------------------
    ! He4SwathReadAttrR8 begins here!
    !-----------------------------------

    ! Read attribute
    status = SwRdAttr( sId, TRIM( attrName ), attrValue )       
    
  END SUBROUTINE He4SwathReadAttrR8

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData1dI2( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He4SwathReadData1dI2 reads a 1-D INTEGER*2 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 8/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (INTEGER*2  ) : 1-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX)

    ! Local variables
    INTEGER                       :: status
    INTEGER(HE4_INT)              :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData1dI2 begins here!
    !-----------------------------------

    ! Set up to read data for a given track
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dI2 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData1dI2

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData1dI4( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He4SwathReadData1dI4 reads a 1-D INTEGER*4 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 8/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (TYPE (XPLEX)     ) : 1-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER,          INTENT(OUT) :: fldData(nX)

    ! Local variables
    INTEGER                       :: status
    INTEGER(HE4_INT)              :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData1dI4 begins here!
    !-----------------------------------

    ! Set up to read data for a given track
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dI4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData1dI4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData1dR4( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He4SwathReadData1dR4 reads a 1-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (TYPE (XPLEX)     ) : 1-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX)

    ! Local variables
    INTEGER                       :: status
    INTEGER(HE4_INT)              :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData1dR4 begins here!
    !-----------------------------------

    ! Set up to read data for a given track
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dR4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData1dR4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData1dR8( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He4SwathReadData1dR8 reads a 1-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (TYPE (XPLEX)     ) : 1-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)    :: sId
    INTEGER(HE4_INT), INTENT(IN)    :: nX
    CHARACTER(LEN=*), INTENT(IN)    :: fldName
    TYPE (XPLEX),           INTENT(INOUT) :: fldData(nX)

    ! Local variables
    INTEGER                         :: status    
    INTEGER(HE4_INT)                :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE4_MAX_CHAR)     :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                         :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData1dR8 begins here!
    !-----------------------------------

    ! Set up dimension info
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dR8 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData1dR8

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData2dI2( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He4SwathReadData2dI2 reads a 2-D INTEGER*2 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 8/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (INTEGER*2  ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(2), stride(2), edge(2)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData2dI2 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    edge   = (/ nX, nY /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dI2 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData2dI2

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData2dI4( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He4SwathReadData2dI4 reads a 2-D INTEGER*4 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 8/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (INTEGER    ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER,          INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(2), stride(2), edge(2)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData2dI4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    edge   = (/ nX, nY /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dI4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData2dI4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData2dR4( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He4SwathReadData2dR4 reads a 2-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (TYPE (XPLEX)     ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(2), stride(2), edge(2)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData2dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    edge   = (/ nX, nY /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dR4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData2dR4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData2dR8( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He4SwathReadData2dR8 reads a 2-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (TYPE (XPLEX)     ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status  
    INTEGER(HE4_INT)              :: start(2), stride(2), count(2)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData2dR8 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    count  = (/ nX, nY /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, count, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dR8 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData2dR8

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData3dI2( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He4SwathReadData3dI2 reads a 3-D INTEGER*2 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 8/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 3rd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (INTEGER*2  ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData3I2 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dI2 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData3dI2

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData3dI4( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He4SwathReadData3dI4 reads a 3-D INTEGER*4 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 8/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (INTEGER    ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER,          INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData3dI4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dI4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData3dI4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData3dR4( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He4SwathReadData3dR4 reads a 3-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (TYPE (XPLEX)     ) : 2-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData3dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dR4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData3dR4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData3dR8( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He4SwathReadData3dR8 reads a 3-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 3rd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (TYPE (XPLEX)     ) : 3-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData3dR8 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dR8 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData3dR8

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData4dI2( sId, fldName, nX, nY, nZ, nW, fldData )

    !======================================================================
    ! Routine He4SwathReadData4dI2 reads a 4-D INTEGER*2 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 9/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 3rd dimension of data array
    ! (6 ) nW      (HE4-INTEGER) : 4th dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (7 ) fldData (INTEGER*2  ) : 4-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ, nW
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX,nY,nZ,nW)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(4), stride(4), edge(4)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData4dI2 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0,  0 /) 
    stride = (/  1,  1,  1,  1 /)
    edge   = (/ nX, nY, nZ, nW /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData4dI2 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData4dI2

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData4dI4( sId, fldName, nX, nY, nZ, nW, fldData )

    !======================================================================
    ! Routine He4SwathReadData4dI4 reads a 4-D INTEGER*4 data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 9/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 2nd dimension of data array
    ! (6 ) nW      (HE4-INTEGER) : 4th dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (7 ) fldData (INTEGER    ) : 4-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ, nW
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER,          INTENT(OUT) :: fldData(nX,nY,nZ,nW)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(4), stride(4), edge(4)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData4dI4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0,  0 /) 
    stride = (/  1,  1,  1,  0 /)
    edge   = (/ nX, nY, nZ, nW /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData4dI4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData4dI4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData4dR4( sId, fldName, nX, nY, nZ, nW, fldData )

    !======================================================================
    ! Routine He4SwathReadData4dR4 reads a 4-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 9/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 3rd dimension of data array
    ! (6 ) nW      (HE4-INTEGER) : 4th dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (7 ) fldData (TYPE (XPLEX)     ) : 4-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ, nW
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX,nY,nZ,nW)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(4), stride(4), edge(4)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData4dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0,  0 /) 
    stride = (/  1,  1,  1,  1 /)
    edge   = (/ nX, nY, nZ, nW /)
    
    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData4dR4 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He4SwathReadData4dR4

!-----------------------------------------------------------------------------

  SUBROUTINE He4SwathReadData4dR8( sId, fldName, nX, nY, nZ, nW, fldData )

    !======================================================================
    ! Routine He4SwathReadData4dR8 reads a 4-D TYPE (XPLEX) data block from
    ! an HDF-EOS4 swath data structure.  This routine is included in the
    ! module interface He4SwathReadData. (bmy, 9/20/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS4 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE4-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE4-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE4-INTEGER) : 3rd dimension of data array
    ! (6 ) nW      (HE4-INTEGER) : 3rd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (7 ) fldData (TYPE (XPLEX)     ) : 3-D array w/ data from HDF-EOS4 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE4_INT), INTENT(IN)  :: nX, nY, nZ, nW
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    TYPE (XPLEX),           INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE4_INT)              :: start(4), stride(4), edge(4)
    CHARACTER(LEN=HE4_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS4 library routines
    INTEGER                       :: SwRdFld
    
    !-----------------------------------
    ! He4SwathReadData4dR8 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0,  0 /) 
    stride = (/  1,  1,  1,  1 /)
    edge   = (/ nX, nY, nZ, nW /)

    ! Read data
    status = SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData4dR8 ("He4SwathModule.f90")'
       CALL He4ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He4SwathReadData4dR8

!-----------------------------------------------------------------------------

  SUBROUTINE makeCharArrayFromCharList( list, separator, array )

    !=====================================================================
    ! Subroutine makeCharArrayFromCharList takes a comma-separated word 
    ! list, and places each word into a separate element of a character 
    ! array. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) list      (CHARACTER) : String with comma-separated words
    ! (2) separator (CHARACTER) : String for separator text
    !
    ! Arguments as output:
    ! --------------------------------------------------------------------
    ! (3) array     (CHARACTER) : Array of substrings
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(IN)  :: list
    CHARACTER(LEN=1           ), INTENT(IN)  :: separator
    CHARACTER(LEN=HE4_MAX_CHAR), INTENT(OUT) :: array(:)

    ! local variables
    INTEGER                                  :: P, N, ind(255)
    CHARACTER(LEN=1)                         :: C

    !----------------------------------------
    ! makeCharArrayFromCharList begins here!
    !----------------------------------------

    ! Initialize
    N      = 1
    ind(:) = 0

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

!------------------------------------------------------------------------------

  FUNCTION He4DataTypeName( nType ) RESULT( typeStr )

    !=====================================================================
    ! Subroutine He4DataTypeName returns a descriptive string given a
    ! HDF-EOS4 data type number. (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) nType     (INTEGER) : HDF-EOS number type 
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    INTEGER, INTENT(IN)         :: nType

    ! Local varaibles
    LOGICAL, SAVE               :: FIRST = .TRUE.
    CHARACTER(LEN=HE4_MAX_CHAR) :: typeStr

    !------------------------------
    ! He4DataTypeName begins here!
    !------------------------------

    ! First-time initialization
    IF ( FIRST ) THEN 
       dataTypeName(:)             = ''
       dataTypeName(DFNT_INT16   ) = 'INTEGER*2'
       dataTypeName(DFNT_UINT16  ) = 'Unsigned INTEGER*2'
       dataTypeName(DFNT_INT32   ) = 'INTEGER*4'
       dataTypeName(DFNT_UINT32  ) = 'Unsigned INTEGER*4'
       dataTypeName(DFNT_INT64   ) = 'INTEGER*8'
       dataTypeName(DFNT_UINT64  ) = 'Unsigned INTEGER*8'
       dataTypeName(DFNT_FLOAT32 ) = 'TYPE (XPLEX)'
       dataTypeName(DFNT_FLOAT64 ) = 'TYPE (XPLEX)'
       dataTypeName(DFNT_FLOAT128) = 'REAL*16'
       !dataTypeName(56          ) = 'LOGICAL'
       !dataTypeName(67          ) = 'LOGICAL'
       dataTypeName(DFNT_CHAR    ) = 'CHARACTER'
       dataTypeName(DFNT_CHAR16  ) = 'CHARACTER'
       dataTypeName(DFNT_UCHAR16 ) = 'CHARACTER'
    ENDIF
    
    ! Return value
    typeStr = dataTypeName(nType)

  END FUNCTION He4DataTypeName

!------------------------------------------------------------------------------

END MODULE He4SwathModule









