MODULE AIRSv5_MOD

!*********************************************************************
! Module AIRSv5_MOD contains variables for reading data from AIRS.
! Each array is used to store data read from AIRS file with data for
! 1 day (jaf, 8/20/07)
!
!  Module Variables:
!  ============================================================================
! (1 ) SurfPressure             : surface pressure read from file
! (2 ) COcol                    : columns read from AIRS file
! (3 ) qual                     : quality flag for CO value
! (4 ) dofs                     : degrees of freedom for each observation
! (5 ) Longitude                : vector of longitudes read from AIRS file
! (6 ) Latitude                 : vector of latitudes read from AIRS file
! (7 ) YYYYMMDD                 : date of observation
! (8 ) FILENAME                 : name of input file
! (9 ) iNObs                    : number of observations in one file
! (10) NObss                    : total number of observations in one day
! (11) COlev                    : indices of CO trapeziods
! (12) NTraps                   : number of trapezoids (9)
! (13) Plev                     : AIRS 100 pressure levels
! (14) NLevs                    : number of pressure levels (100) 
! (15) AvgKer                   : AIRS averaging kernel (on 9 trapezoids)
! (16) BotLev                   : Lowest level near surface (see discussion
!                                  of nSurfSup in documentation)
! (17) LocalT                   : Local solar time of granule center
!
!  Module subroutines:
! ===========================================================================
! (1 ) FILE_INFO --> moved to airs_co_obs_mod.f
! (2 ) INFO_AIRS
! (3 ) INIT_READ_AIRS
! (4 ) CLEANUP_AIRS
!
!*********************************************************************
  USE MYTYPE
  USE COMPLEXIFY
  IMPLICIT NONE

  ! Make everything PUBLIC
  PUBLIC

  !=============================================================
  ! MODULE VARIABLES
  !=============================================================

  TYPE (XPLEX), ALLOCATABLE     :: SurfPressure(:)
  TYPE (XPLEX), ALLOCATABLE     :: COcol(:)
  TYPE (XPLEX), ALLOCATABLE     :: COcd(:,:)
  TYPE (XPLEX), ALLOCATABLE     :: LandFrac(:)
  TYPE (XPLEX), ALLOCATABLE     :: H2Ocd(:,:)
  TYPE (XPLEX), ALLOCATABLE     :: Tsurf(:)
  TYPE (XPLEX), ALLOCATABLE     :: AvgKer(:,:,:)
  TYPE (XPLEX), ALLOCATABLE     :: dofs(:)
  TYPE (XPLEX), ALLOCATABLE     :: Longitude(:)
  TYPE (XPLEX), ALLOCATABLE     :: Latitude(:)
  TYPE (XPLEX), ALLOCATABLE     :: Hour(:)
  TYPE (XPLEX), ALLOCATABLE     :: Minute(:)
  CHARACTER(LEN=3), ALLOCATABLE :: DNFlag(:)
  INTEGER, ALLOCATABLE    :: qual(:)
  INTEGER*2, ALLOCATABLE  :: NTraps(:)
  INTEGER*2, ALLOCATABLE  :: LocalT(:)
  INTEGER, ALLOCATABLE    :: iNObs(:)
  INTEGER, ALLOCATABLE    :: BotLev(:)
  INTEGER                 :: Nobss
  INTEGER                 :: NFiles

  ! NTraps0 is hardwired as this does not change for CO.
  ! This is the ideal number of trapezoids. Some retrievals
  ! will not have 9, and this is what is read by NTraps.
  ! This must be changed if you use a different species.
  INTEGER, PARAMETER      :: NTraps0=9
  INTEGER, PARAMETER      :: NLevs=100
  INTEGER*4, ALLOCATABLE  :: COlev(:)
  TYPE (XPLEX)                  :: Plev(NLevs)

  ! Date Information
  INTEGER                 :: YEAR, MONTH, DAY, IDAY
  CHARACTER(LEN=255), ALLOCATABLE  :: FILENAME(:)      

  INTEGER, ALLOCATABLE:: DOMAIN_OBS(:,:)

  !=================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
  !================================================================= 
  CONTAINS

!------------------------------------------------------------------------------

!!$    SUBROUTINE FILE_INFO
!!$!
!!$!******************************************************************************
!!$!  Subroutine FILE_INFO sets the directory and filename given the date. The
!!$!  date is currently hardwired. (jaf, 6/1/07)
!!$!******************************************************************************
!!$
!!$      ! References to F90 modules
!!$      USE TIME_MOD,  ONLY : EXPAND_DATE
!!$      USE FILE_MOD,  ONLY : IOERROR
!!$
!!$      ! Arguments
!!$      CHARACTER(LEN=255)  :: DIR_AIRS
!!$      CHARACTER(LEN=255)  :: FILENAME_IN      
!!$      CHARACTER(LEN=255)  :: file
!!$      CHARACTER(LEN=  8)  :: YYYYMMDDs
!!$      INTEGER             :: IU_FILE, IOS, IOS1, I
!!$      
!!$      !=================================================================
!!$      ! FILE_INFO begins here!
!!$      !=================================================================
!!$      
!!$      ! Set date and corresponding input AIRS filename
!!$      DIR_AIRS = '/san/as04ro/data/obs/airs/'
!!$      !DIR_AIRS = '/san/as04/home/ctm/mak/AIRS/data_airs/'
!!$      !DIR_AIRS = '/as/data-rw/corrections/as/data/airs/'
!!$      FILENAME_IN='input.txt'
!!$
!!$      IU_FILE=15
!!$      OPEN( IU_FILE, FILE=FILENAME_IN, IOSTAT=IOS )
!!$      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'airs:1' )
!!$
!!$      ! zero counters
!!$      Nfiles = 0
!!$
!!$      ! Figure out how many files to read (#lines in the file):
!!$      CALL SYSTEM('wc -l '//trim(FILENAME_IN)//' > tmp.txt')
!!$      
!!$      OPEN( 5, FILE='tmp.txt', IOSTAT=IOS1 )
!!$      IF ( IOS1 /= 0 ) CALL IOERROR( IOS1, 5, 'tmp:1' )
!!$      
!!$      ! Read #lines
!!$      READ( 5, *, IOSTAT=IOS1  ) NFiles
!!$      IF ( IOS1 /= 0 ) CALL IOERROR( IOS1, 5, 'tmp:2' )
!!$   
!!$      ! Close file
!!$      CLOSE( 5 )
!!$
!!$      PRINT*, 'Number of Files: ', NFiles
!!$      ALLOCATE( iNObs (NFiles) )
!!$
!!$      ALLOCATE( FILENAME(NFiles) )
!!$      DO i = 1, NFiles
!!$     
!!$         IF ( IOS < 0 ) EXIT
!!$         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'airs:2' )
!!$         
!!$         READ( IU_FILE,'(a)',IOSTAT=IOS) file
!!$         IF (i .eq. 1) &
!!$              YYYYMMDDs = file(6:9)//file(11:12)//file(14:15)
!!$         FILENAME(i)=trim(DIR_AIRS)//YYYYMMDDs//'/'//trim(file)
!!$         print*, 'filename:', trim(filename(i))
!!$
!!$      ENDDO
!!$
!!$      ! Close file
!!$      CLOSE( IU_FILE )
!!$
!!$      READ (YYYYMMDDs,*) YYYYMMDD
!!$      
!!$      PRINT*,'Date: ',YYYYMMDD
!!$      
!!$      ! Return to calling program
!!$      END SUBROUTINE FILE_INFO
!!$      
!!$!------------------------------------------------------------------------------

      SUBROUTINE INFO_AIRS
!
!******************************************************************************
!  Subroutine INFO_AIRS Info prints info about all VDATA and SDATA fields 
!  contained within the AIRS HDF file. Based on INFO_MOP02 (bmy, 7/3/03, 
!  4/27/05; jaf 8/15/07)
!******************************************************************************

      ! References to F90 modules
      USE He4SwathModule
      USE He4IncludeModule

      INTEGER                     :: fId, sId
      CHARACTER(LEN=HE4_MAX_CHAR) :: swathname = &
           'L2_Support_atmospheric&surface_product'

      ! For dimension information
      INTEGER                     :: nDims
      INTEGER                     :: dims(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: dimNames(HE4_MAX_FLDS)

      ! For swath attributes
      INTEGER                     :: nAttrs
      TYPE (XPLEX)                      :: attrValue(HE4_MAX_ATRS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: attrName(HE4_MAX_ATRS)
      
      ! For geolocation fields
      INTEGER                     :: nGeo
      INTEGER                     :: geoRank(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: geoName(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: geoType(HE4_MAX_FLDS)
      
      ! For data fields
      INTEGER                     :: nData
      INTEGER                     :: dataRank(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: dataName(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: dataType(HE4_MAX_FLDS)

      !=================================================================
      ! INFO_AIRS begins here!
      !=================================================================

      ! Set verbose output
      CALL He4VerboseOutput( .TRUE. )
      
      ! Open HDF-EOS5 swath and get the file ID number
      CALL He4FileOpen( FILENAME(1), fId )
      
      ! Attach to swath and get swath ID number
      CALL He4SwathAttach( fId, swathName, sId )

      ! Get swath attributes
      CALL He4SwathAttrs( sId, nAttrs, attrName, attrValue )

      ! Get swath dimension info
      CALL He4SwathDimInfo( sId, nDims, dims, dimNames )

      ! Get information about geolocation fields
      CALL He4SwathGeoFldInfo( sId, nGeo, geoRank, geoName, geoType )

      ! Get information about data fields
      CALL He4SwathDataFldInfo( sId, nData, dataRank, dataName, dataType )

      !------------------------------------------------------------------------
      ! Cleanup and quit
      !------------------------------------------------------------------------

      ! Detach from swath
      CALL He4SwathDetach( sId )
 
      ! Close HDF-EOS5 file
      CALL He4FileClose( fId )

      ! Return to calling program
    END SUBROUTINE INFO_AIRS

!------------------------------------------------------------------------------

  SUBROUTINE INIT_READ_AIRS

!********************************************************************
! SUBROUTINE INIT_READ_AIRS allocates all module arrays and reads data
! into them from the HDF file. Based on subroutine READ_MOP02.
!*********************************************************************

      ! References to F90 modules
      USE He4ErrorModule
      USE He4SwathModule
      USE He4GridModule
      USE He4IncludeModule
      USE FILE_MOD, ONLY : IOERROR
      USE JULDAY_MOD, ONLY : CALDATE, JULDAY
      USE TIME_MOD, ONLY : YMD_EXTRACT

      ! Arguments
      ! File and Swath Info
      INTEGER  :: fId, sId, as
      CHARACTER(LEN=HE4_MAX_CHAR) :: SWATHNAME = &
           'L2_Support_atmospheric&surface_product'

      ! For dimension information
      INTEGER                     :: nDims
      INTEGER                     :: dims(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: dimNames(HE4_MAX_FLDS)

      ! For data fields
      INTEGER                     :: nData
      INTEGER                     :: dataRank(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: dataName(HE4_MAX_FLDS)
      CHARACTER(LEN=HE4_MAX_CHAR) :: dataType(HE4_MAX_FLDS)

      ! Fill value
      CHARACTER(LEN=HE4_MAX_CHAR) :: name
      TYPE (XPLEX)                      :: dataFill

      ! For data field info
      INTEGER                     :: fldDims(HE4_MAX_DIMS)
      INTEGER                     :: fldRank
      CHARACTER(LEN=HE4_MAX_CHAR) :: fldType
      CHARACTER(LEN=HE4_MAX_CHAR) :: fldDimNames(HE4_MAX_DIMS)

      ! For swath attribute info
      INTEGER                     :: nAttrs
      CHARACTER(LEN=HE4_MAX_CHAR) :: attrName(HE4_MAX_ATRS)      
      INTEGER*2                   :: attrValue(HE4_MAX_ATRS)
      INTEGER                     :: status, strBufSize
      CHARACTER(LEN=HE4_MAX_CHAR) :: attrList
    
      ! HDF_EOS5 library routines
      INTEGER                                  :: SwInqAttrs
      INTEGER                                  :: SwRdAttr

      INTEGER(HE4_INT)            :: nX, nY, nZ, nW
      INTEGER                     :: N, I, J, K
      INTEGER                     :: YYYYMMDD_tmp,HHMMSS,SS
      INTEGER*4                   :: track, xtrack
      TYPE (XPLEX), ALLOCATABLE         :: temp4(:,:)
      TYPE (XPLEX), ALLOCATABLE         :: temp4_3d(:,:,:)
      TYPE (XPLEX), ALLOCATABLE         :: temp4_4d(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE         :: temp8(:,:)
      INTEGER, ALLOCATABLE        :: HH(:,:),MM(:,:)
      INTEGER*2, ALLOCATABLE      :: tempi2(:,:)
      INTEGER*4, ALLOCATABLE      :: tempi4(:,:)
      INTEGER*4, ALLOCATABLE      :: tempi_1d(:)
      CHARACTER(LEN=6)            :: temp_str
      TYPE (XPLEX), ALLOCATABLE         :: Time(:)
      TYPE (XPLEX)                      :: fyr,fday

      ! For reading pressure file
      CHARACTER(LEN=255)          :: Pfile
      INTEGER                     :: IU_FILE, IOS, z, s
      TYPE (XPLEX)                      :: ptemp

      !=================================================================
      ! Init_Read_AIRS begins here!
      !=================================================================


      !----------------------------------------------------------------
      ! Read dimensions and allocate arrays
      !----------------------------------------------------------------

!!$      DO I = 1, NFiles
!!$
!!$      ! Open HDF-EOS5 swath and get the file ID number
!!$      CALL He4FileOpen( FILENAME(I), fId )
!!$      
!!$      ! Attach to swath and get swath ID number
!!$      CALL He4SwathAttach( fId, swathName, sId )
!!$
!!$      ! Get swath dimension info
!!$      CALL He4SwathDimInfo( sId, nDims, dims, dimNames )
!!$
!!$      DO N = 1, nDims
!!$         IF ( TRIM(dimNames(N)) == 'GeoXTrack' ) xtrack = dims(N)
!!$         IF ( TRIM(dimNames(N)) == 'GeoTrack'  )  track = dims(N)
!!$      ENDDO
!!$
!!$      IF (xtrack .ne. 30) PRINT*,FILENAME(I),': xtrack=',xtrack
!!$      IF ( track .ne. 45) PRINT*,FILENAME(I),': track=',track
!!$
!!$      ! Detach from swath
!!$      CALL He4SwathDetach( sId )
!!$      
!!$      ! Close HDF-EOS5 file
!!$      CALL He4FileClose( fId )
!!$
!!$      iNObs(I) = track*xtrack
!!$      ENDDO

      iNObs(:)=1350
      Nobss = SUM(iNObs)
      PRINT*,'### Number of Observations: ',Nobss

      ALLOCATE(   Longitude ( Nobss) )
      ALLOCATE(    Latitude ( Nobss) )
      ALLOCATE(       Hour  ( NObss) )
      ALLOCATE(     Minute  ( NObss) )
      ALLOCATE( SurfPressure( Nobss) )
      ALLOCATE(     DNFlag  ( NObss) )
      ALLOCATE(        COcol( Nobss) )
      ALLOCATE(        Tsurf( Nobss) )
      ALLOCATE(     LandFrac( NObss) )
      ALLOCATE(         qual( Nobss) )
      ALLOCATE(       NTraps( Nobss) )
      ALLOCATE(         dofs( Nobss) )
      ALLOCATE(       BotLev( Nobss) )
      ALLOCATE( AvgKer( NTraps0, NTraps0, Nobss) )
      ALLOCATE(         COcd( NLevs, Nobss ) )
      ALLOCATE(         H2Ocd( NLevs, Nobss ) )
      ALLOCATE(       LocalT( Nobss) )

      Longitude(:) = 0d0
      Latitude(:) = 0d0
      Hour(:) = 0d0
      Minute(:) = 0d0
      SurfPressure(:)=0d0
      DNflag(:) = 'NA'
      COcol(:)=0d0
      Tsurf(:)=0d0
      LandFrac(:)=0d0
      qual(:)=0
      NTraps(:)=0
      dofs(:)=0
      BotLev(:)=0
      AvgKer(:,:,:)=0d0
      LocalT(:)=0
      COcd(:,:)=0d0
      H2Ocd(:,:)=0d0

      ! Note: Trapezoid layer indices and pressure levels do not change,
      !       so these are read only once, from the first file

      ! Open HDF-EOS5 swath and get the file ID number
      CALL He4FileOpen( FILENAME(1), fId )

      ! Attach to swath and get swath ID number
      CALL He4SwathAttach( fId, SWATHNAME, sId )

      !----------------------------------------------------------------
      ! Read trapezoid layer information from HDF-EOS file (1-D INT*4)
      !----------------------------------------------------------------
      ! Field name
      name = 'CO_trapezoid_layers'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*4 variables BEFORE array allocation!
      nX = fldDims(1)

      ! Allocate array
      ALLOCATE( tempi_1d( nX ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'CO_trapezoid_layers' )
      tempi_1d = 0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, tempi_1d )

      ALLOCATE( COlev( nX+1 ) )
      COlev=0
      COlev(1:nX)=tempi_1d
      COlev(nX+1)=NLevs

      DEALLOCATE(tempi_1d)

      !----------------------------------------------------------------
      ! Read pressure level information from HDF-EOS file (1-D INT*4)
      !----------------------------------------------------------------
      ! Field name
      name = 'pressSupp'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*4 variables BEFORE array allocation!
      nX = fldDims(1)

      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, Plev )

      ! Detach from swath
      CALL He4SwathDetach( sId )
 
      ! Close HDF-EOS5 file
      CALL He4FileClose( fId )

      !----------------------------------------------------------------
      ! Loop over files to read all info for one day
      !----------------------------------------------------------------
      DO I = 1, NFiles
      
      ! Open HDF-EOS5 swath and get the file ID number
      CALL He4FileOpen( FILENAME(I), fId )

      ! Attach to swath and get swath ID number
      CALL He4SwathAttach( fId, SWATHNAME, sId )

      !----------------------------------------------------------------
      ! Read local solar time from HDF-EOS file (INTEGER*2)
      ! Read day/night flag from HDF-EOS file (INTEGER*2)
      !----------------------------------------------------------------

      ! Get list of attribute names
      nAttrs = SwInqAttrs( sId, attrList, strBufSize )
      ! Separate list into array
      CALL makeCharArrayFromCharList( attrList, ',', attrName )

      DO N = 1, nAttrs
         status = SwRdAttr( sId, TRIM( attrName(N) ), attrValue(N) )
      ENDDO
      DO N = 1, nAttrs
         IF ( TRIM( attrName(N) ) .EQ. 'LocTimeGranuleCen' ) THEN
            LocalT( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I))) = &
                 attrValue(N)
            EXIT
         ENDIF
      ENDDO
      temp_str=' '
      CALL He4SwathReadAttr( sId, 'DayNightFlag', temp_str)
      DNFlag( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I))) = temp_str
      temp_str=' '
      s = SUM(iNObs(1:I))-iNObs(I)
      DO z = 1, 1350
         IF( DNFlag(s+z) .ne. 'Day') then
            GOTO 222
         ENDIF
      ENDDO
  
      !----------------------------------------------------------------
      ! Read latitude data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------
      print*, 'read other stuff before latitude'
      call flush(6)

      ! Field name
      name = 'Latitude'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp8( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'Latitude' )
      temp8 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp8 )

      Latitude( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp8, (/ iNObs(I) /) )
      DEALLOCATE(temp8)

      !----------------------------------------------------------------
      ! Read longitude data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'Longitude'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp8( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'Longitude' )
      temp8 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp8 )

      Longitude( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp8, (/ iNObs(I) /) )
      DEALLOCATE(temp8)

      !----------------------------------------------------------------
      ! Read time data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'Time'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp8( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'Time' )
      ALLOCATE( HH( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'Time:HH' )
      ALLOCATE( MM( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'Time:MM' )
      temp8 = 0d0
      HH = 0
      MM = 0
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp8 )

     temp8 = JULDAY(1993,1,1d0) + (temp8 / 86400d0)

     DO J = 1, nX
     DO K = 1, nY
        CALL CALDATE( Temp8(J,K), YYYYMMDD_tmp, HHMMSS )
        CALL YMD_EXTRACT(HHMMSS,HH(J,K),MM(J,K),SS)
     ENDDO
     ENDDO

      HOUR( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( HH, (/ iNObs(I) /) )
      MINUTE( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( MM, (/ iNObs(I) /) )
      DEALLOCATE(temp8)
      DEALLOCATE(HH)
      DEALLOCATE(MM)

!      DEALLOCATE(Time)

      !----------------------------------------------------------------
      ! Read surface pressure data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'PSurfStd'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp4( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'PSurfStd')
      temp4 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp4 )

      SurfPressure( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4, (/ iNObs(I) /) )
      DEALLOCATE(temp4)

      !----------------------------------------------------------------
      ! Read CO column data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'CO_total_column'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp4( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'CO_total_column' )
      temp4 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp4 )

      COcol( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4, (/ iNObs(I) /) )
      DEALLOCATE(temp4)

      !----------------------------------------------------------------
      ! Read land fraction data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'landFrac'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp4( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'landFrac' )
      temp4 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp4 )

      LandFrac( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4, (/ iNObs(I) /) )
      DEALLOCATE(temp4)

      !----------------------------------------------------------------
      ! Read surface air temperature data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'TSurfAir'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp4( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'TSurfAir')
      temp4 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp4 )


      Tsurf( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4, (/ iNObs(I) /) )
      DEALLOCATE(temp4)

!!$      !----------------------------------------------------------------
!!$      ! Read CO layer column density data from HDF-EOS file (3-D TYPE (XPLEX))
!!$      !----------------------------------------------------------------
!!$
!!$      ! Field name
!!$      name = 'COCDSup'
!!$      
!!$      ! Get fill value (if necessary to strip out missing data values)
!!$      CALL He4SwathFillValue( sId, name, dataFill )
!!$      
!!$      ! Get field info (array dimensions are in fldDims)
!!$      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
!!$      
!!$      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
!!$      nX = fldDims(1)
!!$      nY = fldDims(2)
!!$      nZ = fldDims(3)
!!$
!!$      ! Allocate array
!!$      ALLOCATE( temp4_3d( nX, nY, nZ ), stat=as )
!!$      IF ( as /= 0 ) CALL He4AllocErr( 'COCDSup' )
!!$      temp4_3d = 0d0
!!$  
!!$      ! Read data from swath
!!$      CALL He4SwathReadData( sId, name, nX, nY, nZ, temp4_3d )
!!$
!!$      COcd( :, SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
!!$           RESHAPE( temp4_3d, (/ NLevs, iNObs(I) /) )
!!$      DEALLOCATE(temp4_3d)
!!$
      !----------------------------------------------------------------
      ! Read H2O layer column density data from HDF-EOS file (3-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'H2OCDSup'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
      nZ = fldDims(3)
  
      ! Allocate array
      ALLOCATE( temp4_3d( nX, nY, nZ ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'H2OCDSup' )
      temp4_3d = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, nZ, temp4_3d )

      H2Ocd( :, SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4_3d, (/ NLevs, iNObs(I) /) )
      DEALLOCATE(temp4_3d)

      !----------------------------------------------------------------
      ! Read averaging kernel data from HDF-EOS file (4-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'CO_ave_kern'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
      nZ = fldDims(3)
      nW = fldDims(4)
  
      ! Allocate array
      ALLOCATE( temp4_4d( nX, nY, nZ, nW), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'CO_ave_kern' )
      temp4_4d(:,:,:,:) = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, nZ, nW, temp4_4d )

      AvgKer( :, :, SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4_4d, (/ NTraps0, NTraps0, iNObs(I)/) )

      DEALLOCATE(temp4_4d)

      !----------------------------------------------------------------
      ! Read degrees of freedom data from HDF-EOS file (2-D TYPE (XPLEX))
      !----------------------------------------------------------------

      ! Field name
      name = 'CO_dof'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate array
      ALLOCATE( temp4( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'CO_dof' )
      temp4 = 0d0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, temp4 )

      dofs( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( temp4, (/ iNObs(I) /) )
      DEALLOCATE(temp4)

      !----------------------------------------------------------------
      ! Read quality data from HDF-EOS file (2-D UNSIGNED INTEGER*4)
      !----------------------------------------------------------------

      ! Field name
      name = 'Qual_CO'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate arrays
      ALLOCATE( tempi2( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'Qual_CO' )
      tempi2 = 0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, tempi2 )

      qual( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( tempi2, (/ iNObs(I) /) )

      DEALLOCATE(tempi2)

      !----------------------------------------------------------------
      ! Read number of AK entries from HDF-EOS file (2-D INTEGER*2)
      !----------------------------------------------------------------

      ! Field name
      name = 'num_CO_Func'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate arrays
      ALLOCATE( tempi2( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'num_CO_Func' )
      tempi2 = 0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, tempi2 )

      NTraps( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( tempi2, (/ iNObs(I) /) )

      DEALLOCATE(tempi2)

      !----------------------------------------------------------------
      ! Read bottom level from HDF-EOS file (2-D INTEGER*4)
      !----------------------------------------------------------------

      ! Field name
      name = 'nSurfSup'
      
      ! Get fill value (if necessary to strip out missing data values)
      CALL He4SwathFillValue( sId, name, dataFill )
      
      ! Get field info (array dimensions are in fldDims)
      CALL He4SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
      
      ! Copy dims to INTEGER*8 variables BEFORE array allocation!
      nX = fldDims(1)
      nY = fldDims(2)
  
      ! Allocate arrays
      ALLOCATE( tempi4( nX, nY ), stat=as )
      IF ( as /= 0 ) CALL He4AllocErr( 'nSurfSup' )
      tempi4 = 0
  
      ! Read data from swath
      CALL He4SwathReadData( sId, name, nX, nY, tempi4 )

      BotLev( SUM(iNObs(1:I))-iNObs(I)+1 : SUM(iNObs(1:I)) ) = &
           RESHAPE( tempi4, (/ iNObs(I) /) )

      DEALLOCATE(tempi4)

      !----------------------------------------------------------------
      ! Cleanup and quit
      !----------------------------------------------------------------
      
222   CONTINUE

      ! Detach from swath
      CALL He4SwathDetach( sId )
      
      ! Close HDF-EOS5 file
      CALL He4FileClose( fId )

      ENDDO
         
      !----------------------------------------------------------------
      ! Echo min & max values
      !----------------------------------------------------------------
      PRINT*, '### Lat min, max: ', MINVAL( Latitude  ), MAXVAL( Latitude  )
      PRINT*, '### Lon min, max: ', MINVAL( Longitude ), MAXVAL( Longitude )
      PRINT*, '### LocalT: ', MINVAL(LocalT)
      PRINT*, '### Hour min, max: ', MINVAL( Hour ), MAXVAL( Hour )
      PRINT*, '### Minute min, max: ', MINVAL( Minute ), MAXVAL( Minute )
      PRINT*, '### SurfPressure min, max:',  MINVAL( SurfPressure ), &
           MAXVAL( SurfPressure )
      PRINT*, '### COcol min, max: ', MINVAL( COcol ), MAXVAL( COcol )
      PRINT*, '### Tsurf min, max: ', MINVAL( Tsurf ), MAXVAL( Tsurf )
      PRINT*, '### LandFrac min, max: ', MINVAL( LandFrac ), MAXVAL( LandFrac )
      PRINT*, '### COcd min, max: ', MINVAL( COcd ), MAXVAL( COcd )
      PRINT*, '### H2Ocd min, max: ', MINVAL( H2Ocd ), MAXVAL( H2Ocd )
!      PRINT*, '### AvgKer min, max: ', MINVAL( AvgKer ), MAXVAL( AvgKer )
!      PRINT*, '### dofs min, max: ', MINVAL( dofs ), MAXVAL( dofs )
!      PRINT*, '### qual min, max: ', MINVAL( qual ), MAXVAL( qual )
!      PRINT*, '### BotLev min, max: ', MINVAL( BotLev ), MAXVAL( BotLev )
!      PRINT*, '### NTraps min, max: ', MINVAL( NTraps ), MAXVAL( NTraps )

!      PRINT*, '### COlev: ', COlev
      PRINT*, '### Plev min, max: ',Plev(1),Plev(NLevs)
!      PRINT*, '### Day Night Flag: ', DNFlag(1350)

      ! Return to calling program
    END SUBROUTINE INIT_READ_AIRS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_AIRS
!
!******************************************************************************
!  Subroutine CLEANUP_AIRS deallocates all module arrays (jaf, 8/15/07)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_AIRS begins here!
      !=================================================================
      IF ( ALLOCATED( iNObs           ) ) DEALLOCATE( iNObs           )
      IF ( ALLOCATED( Filename        ) ) DEALLOCATE( Filename        )
      IF ( ALLOCATED( Latitude        ) ) DEALLOCATE( Latitude        )
      IF ( ALLOCATED( Longitude       ) ) DEALLOCATE( Longitude       )
      IF ( ALLOCATED( Hour            ) ) DEALLOCATE( Hour            )
      IF ( ALLOCATED( Minute          ) ) DEALLOCATE( Minute          )
      IF ( ALLOCATED( SurfPressure    ) ) DEALLOCATE( SurfPressure    )
      IF ( ALLOCATED( DNFlag          ) ) DEALLOCATE( DNFlag          )
      IF ( ALLOCATED( COcol           ) ) DEALLOCATE( COcol           )
      IF ( ALLOCATED( Tsurf           ) ) DEALLOCATE( Tsurf           )
      IF ( ALLOCATED( LandFrac        ) ) DEALLOCATE( LandFrac        )
      IF ( ALLOCATED( COcd            ) ) DEALLOCATE( COcd            )
      IF ( ALLOCATED( H2Ocd           ) ) DEALLOCATE( H2Ocd           )
      IF ( ALLOCATED( AvgKer          ) ) DEALLOCATE( AvgKer          )
      IF ( ALLOCATED( qual            ) ) DEALLOCATE( qual            )
      IF ( ALLOCATED( dofs            ) ) DEALLOCATE( dofs            )
      IF ( ALLOCATED( COlev           ) ) DEALLOCATE( COlev           )
      IF ( ALLOCATED( BotLev          ) ) DEALLOCATE( BotLev          )
      IF ( ALLOCATED( NTraps          ) ) DEALLOCATE( NTraps          )
      IF ( ALLOCATED( LocalT          ) ) DEALLOCATE( LocalT          )
      IF ( ALLOCATED( DOMAIN_OBS      ) ) DEALLOCATE( DOMAIN_OBS      )

      ! Return to calling program
      END SUBROUTINE CLEANUP_AIRS
!------------------------------------------------------------------------------
! End of module
END MODULE AIRSv5_MOD
