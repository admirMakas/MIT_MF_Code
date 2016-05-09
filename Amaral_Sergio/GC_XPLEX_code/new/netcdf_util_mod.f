!$Id: netcdf_util_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: netcdf_util_mod
!
!\\
!\\
! !INTERFACE:
!
      MODULE NETCDF_UTIL_MOD
!
! !USES:
!
      USE NETCDF
      IMPLICIT NONE
      PRIVATE
!     

!
! !PUBLIC MEMBER FUNCTIONS:
!
      ! Functions
      PUBLIC  :: NCDF_GET_VARID  ! Return variable ID
      PUBLIC  :: NCDF_GET_DIMID  ! Return variable ID

      ! Subroutines
      PUBLIC  :: NCDF_OPEN_FOR_READ
      PUBLIC  :: NCDF_CLOSE

      PUBLIC  :: NCDF_VAR_EXIST
      PUBLIC  :: NCDF_DIM_EXIST

      PUBLIC  :: NCDF_GET_VAR   ! Get data from a netCDF file
      
      PUBLIC  :: NCDF_SYNC      ! Push data in buffers to the file ahead of close

      PUBLIC  :: CHECK_NCDF     ! Test for error on netCDF API functions
      PUBLIC  :: NCDF_INIT      ! Set intial variables, etc.     
!
! !PRIVATE MEMBER FUNCTIONS:
!     
! !PUBLIC DATA MEMBERS:
      INTEGER, PUBLIC       :: NCDF_CHAR
      INTEGER, PUBLIC       :: NCDF_INT
      INTEGER, PUBLIC       :: NCDF_REAL

      INTEGER, PUBLIC       :: DEF_LEV ! Deflation level
      
      ! Dimension IDs
      INTEGER, PUBLIC :: LVL_DIMID, LVL_VARID   ! Vertical level dimension
      INTEGER, PUBLIC :: LVLE_DIMID, LVLE_VARID ! Vertical edge dimension
      INTEGER, PUBLIC :: LAT_DIMID, LAT_VARID   ! Latitude dimension
      INTEGER, PUBLIC :: LON_DIMID, LON_VARID   ! Longitude dimension
      INTEGER, PUBLIC :: TRACER_DIMID           ! Tracer dimension
      INTEGER, PUBLIC :: STRLEN_8_DIMID         ! CHAR(LEN=8) dimension
      INTEGER, PUBLIC :: STRLEN_16_DIMID        ! CHAR(LEN=16) dimension
      INTEGER, PUBLIC :: STRLEN_64_DIMID        ! CHAR(LEN=64) dimension
      INTEGER, PUBLIC :: PL_DIMID, PL_VARID     ! Prod or Loss dimension
      INTEGER, PUBLIC :: REC_DIMID, REC_VARID   ! Record dimension (time,unlimited)
      INTEGER, PUBLIC :: TRACER_NAME_VARID

!
! !REMARKS:
!
!  References:
!  ============================================================================
!  (1 )
! !REVISION HISTORY:
!  2 Feb 2011 - L. Murray - Initial version.
!  1 Jun 2011 - L. Murray - Read-only version of the code.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      INTERFACE NCDF_GET_VAR
          ! Compiler chooses which one is appropriate when called
          MODULE PROCEDURE NCDF_GET_VAR_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_1D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_2D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_3D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_4D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_5D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_6D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_7D_CHAR
          MODULE PROCEDURE NCDF_GET_VAR_INT
          MODULE PROCEDURE NCDF_GET_VAR_1D_INT
          MODULE PROCEDURE NCDF_GET_VAR_2D_INT
          MODULE PROCEDURE NCDF_GET_VAR_3D_INT
          MODULE PROCEDURE NCDF_GET_VAR_4D_INT
          MODULE PROCEDURE NCDF_GET_VAR_5D_INT
          MODULE PROCEDURE NCDF_GET_VAR_6D_INT
          MODULE PROCEDURE NCDF_GET_VAR_7D_INT
          MODULE PROCEDURE NCDF_GET_VAR_REAL
          MODULE PROCEDURE NCDF_GET_VAR_1D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_2D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_3D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_4D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_5D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_6D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_7D_REAL
          MODULE PROCEDURE NCDF_GET_VAR_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_1D_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_2D_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_3D_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_4D_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_5D_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_6D_DBLE
          MODULE PROCEDURE NCDF_GET_VAR_7D_DBLE
      END INTERFACE

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_NCDF( STATUS )
      USE ERROR_MOD, ONLY  : GEOS_CHEM_STOP
      IMPLICIT NONE

      INTEGER, INTENT( IN ) :: STATUS
      
      IF ( STATUS /= nf90_noerr ) THEN
         WRITE(6,*) 'Error with netCDF'
         PRINT*,trim( nf90_strerror(status) )
         CALL GEOS_CHEM_STOP
      END IF
      
      END SUBROUTINE CHECK_NCDF
      
!------------------------------------------------------------------------------

      SUBROUTINE NCDF_INIT
      ! Call only once
      IMPLICIT NONE    
      
      ! Set type flags
      NCDF_CHAR = NF90_CHAR
      NCDF_INT  = NF90_INT
      NCDF_REAL = NF90_REAL

      END SUBROUTINE NCDF_INIT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_OPEN_FOR_READ( NCID, FILENAME )
      IMPLICIT NONE
      
      CHARACTER(LEN=*), INTENT( IN  ) :: FILENAME
      INTEGER, INTENT( OUT )          :: NCID
      
      write(6,*) 'netCDF: Opening file for read: ',trim(filename)
      call check_ncdf( nf90_open( trim(filename), nf90_nowrite, ncid ) )

      END SUBROUTINE NCDF_OPEN_FOR_READ

!------------------------------------------------------------------------------

      FUNCTION NCDF_DIM_EXIST( NCID, DIMNAME ) RESULT( DIMEXIST )

      INTEGER,          INTENT( IN  ) :: NCID
      CHARACTER(LEN=*), INTENT( IN  ) :: dimName
      LOGICAL                         :: dimExist

      INTEGER :: dimid, status

      dimExist = .false.

      ! Inquire about the the variable id
      status =  nf90_inq_dimid( ncid, trim(dimname), dimid )      
      if ( status .eq. nf90_NoErr ) dimExist = .true.

      END FUNCTION NCDF_DIM_EXIST

!------------------------------------------------------------------------------

      FUNCTION NCDF_VAR_EXIST( NCID, VARNAME ) RESULT( VAREXIST )

      INTEGER,          INTENT( IN  ) :: NCID
      CHARACTER(LEN=*), INTENT( IN  ) :: VARNAME
      LOGICAL                         :: varExist

      INTEGER :: varid, status

      varExist = .false.

      ! Inquire about the the variable id
      status =  nf90_inq_varid( ncid, trim(varname), varid )      
      if ( status .eq. nf90_NoErr ) varExist = .true.
      ! status .eq. -49 is variable not found

      END FUNCTION NCDF_VAR_EXIST

!------------------------------------------------------------------------------

      FUNCTION NCDF_GET_DIMID( NCID, DIMNAME ) RESULT( DIMID )
      

      INTEGER,          INTENT( IN  ) :: NCID
      CHARACTER(LEN=*), INTENT( IN  ) :: DIMNAME
      INTEGER                         :: DIMID

      ! Get the dimiable id
      call check_ncdf( nf90_inq_dimid( ncid, trim(dimname), dimid ) )

      END FUNCTION NCDF_GET_DIMID

!------------------------------------------------------------------------------

      FUNCTION NCDF_GET_VARID( NCID, VARNAME ) RESULT( VARID )
      

      INTEGER,          INTENT( IN  ) :: NCID
      CHARACTER(LEN=*), INTENT( IN  ) :: VARNAME
      INTEGER                         :: VARID

      ! Get the variable id
      call check_ncdf( nf90_inq_varid( ncid, trim(varname), varid ) )

      END FUNCTION NCDF_GET_VARID

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_1D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_1D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_2D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_2D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_3D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_3D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_4D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_4D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_5D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_5D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_6D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_6D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_7D_CHAR(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      CHARACTER(LEN=*),  INTENT( OUT )   :: ARRAY(:,:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif

      END SUBROUTINE NCDF_GET_VAR_7D_CHAR

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_INT(NCID, VARID, VALUE )
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: VALUE
      
      ! Get the variable
      call check_ncdf( nf90_get_var( ncid, varID, value ) )

      END SUBROUTINE NCDF_GET_VAR_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_1D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_1D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_2D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_2D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_3D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_3D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_4D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_4D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_5D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_5D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_6D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_6D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_7D_INT(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      INTEGER,           INTENT( OUT )   :: ARRAY(:,:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_7D_INT

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_REAL(NCID, VARID, VALUE )
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: VALUE
      
      ! Get the variable
      call check_ncdf( nf90_get_var( ncid, varID, value ) )
      
      END SUBROUTINE NCDF_GET_VAR_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_1D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_1D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_2D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_2D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_3D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_3D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_4D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_4D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_5D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_5D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_6D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_6D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_7D_REAL(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*4,            INTENT( OUT )   :: ARRAY(:,:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_7D_REAL

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_DBLE(NCID, VARID, VALUE )
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: VALUE
      
      ! Get the variable
      call check_ncdf( nf90_get_var( ncid, varID, value ) )
      
      END SUBROUTINE NCDF_GET_VAR_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_1D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_1D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_2D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_2D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_3D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_3D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_4D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_4D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_5D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_5D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_6D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_6D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_GET_VAR_7D_DBLE(NCID, VARID, ARRAY, START, COUNT)
      
      IMPLICIT NONE

      INTEGER,           INTENT( IN  )   :: NCID
      INTEGER,           INTENT( IN  )   :: VARID
      REAL*8,            INTENT( OUT )   :: ARRAY(:,:,:,:,:,:,:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: START(:)
      INTEGER, OPTIONAL, INTENT( IN  )   :: COUNT(:)
      
      ! Get the variable
      if ( present( start ) .and. present( count ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start, count=count ) )
      else if ( present( start ) ) then
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  start=start ) )
      else if ( present( count ) ) then 
         call check_ncdf( nf90_get_var( ncid, varID, array,
     &                                  count=count ) )
      else
         call check_ncdf( nf90_get_var( ncid, varID, array ) )
      endif
    
      END SUBROUTINE NCDF_GET_VAR_7D_DBLE

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_SYNC( NCID )
      

      INTEGER, INTENT( IN ) :: NCID      
      call check_ncdf( nf90_sync( NCID ) )
      
      END SUBROUTINE NCDF_SYNC

!------------------------------------------------------------------------------

      SUBROUTINE NCDF_CLOSE( NCID )
      USE NETCDF

      INTEGER, INTENT( IN ) :: NCID      
      call check_ncdf( nf90_close( NCID ) )
      
      END SUBROUTINE NCDF_CLOSE
      
      END MODULE NETCDF_UTIL_MOD
