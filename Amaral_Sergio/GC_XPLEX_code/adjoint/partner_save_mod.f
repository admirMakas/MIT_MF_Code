      MODULE PARTNER_SAVE_MOD

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "aircraft_nox_mod.f"
      !=================================================================
     

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!=========================================================================

!-------------------------------------------------------------------------

      SUBROUTINE SAVE_EMS_AC_FILE(FILENAME_NC, DATA_NC)
! saving in  netcdf format
      USE NETCDF
#     include "CMN_SIZE"

      INTEGER :: NCID
      INTEGER :: NC_DIM_ID(3)
      INTEGER :: LON_DIM, LAT_DIM, ALT_DIM 
      INTEGER :: TRACER_DIM = 6 ! 6 kinds of aircraft emissions
      INTEGER :: VAR_ID(100)
      INTEGER :: LON_ID=1, LAT_ID=2, ALT_ID=3, TRACER_ID=4
      CHARACTER(LEN=255), INTENT(IN) :: FILENAME_NC
      TYPE (XPLEX), INTENT(IN)    :: DATA_NC(:,:,:,:)
      CHARACTER(LEN=255) :: EMS_AC_NAME(6)
      INTEGER :: N


      LON_DIM=IIPAR
      LAT_DIM=JJPAR
      ALT_DIM=LLPAR
  
      NCID=1
      EMS_AC_NAME(1)='AC_NOx'
      EMS_AC_NAME(2)='AC_HC'
      EMS_AC_NAME(3)='AC_PMNV'
      EMS_AC_NAME(4)='AC_PMFO'
      EMS_AC_NAME(5)='AC_CO'
      EMS_AC_NAME(6)='AC_SOx'

       CALL NC_CHECK(NF90_CREATE(TRIM(FILENAME_NC), NF90_CLOBBER ,NCID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Longitude',LON_DIM,LON_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Latitude',LAT_DIM,LAT_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Altitude',ALT_DIM,ALT_ID))
       NC_DIM_ID = (/LON_ID,LAT_ID,ALT_ID/)
      DO N = 1, TRACER_DIM 
       CALL NC_CHECK(NF90_DEF_VAR(NCID, TRIM(EMS_AC_NAME(N)), 
     &              NF90_DOUBLE, NC_DIM_ID, VAR_ID(N)))
      ENDDO
      CALL NC_CHECK(NF90_ENDDEF(NCID))
      DO N = 1, TRACER_DIM 
      CALL NC_CHECK(NF90_PUT_VAR(NCID,VAR_ID(N),(DATA_NC(:,:,:,N)%r)))
      ENDDO
      CALL NC_CHECK(NF90_CLOSE(NCID))

      END SUBROUTINE SAVE_EMS_AC_FILE
!----------------------------------------------------------------------
      SUBROUTINE SAVE_EMS_AC_FILE_X(FILENAME_NC, DATA_NC)
! saving in  netcdf format
      USE NETCDF
#     include "CMN_SIZE"

      INTEGER :: NCID
      INTEGER :: NC_DIM_ID(3)
      INTEGER :: LON_DIM, LAT_DIM, ALT_DIM
      INTEGER :: TRACER_DIM = 6 ! 6 kinds of aircraft emissions
      INTEGER :: VAR_ID(100)
      INTEGER :: LON_ID=1, LAT_ID=2, ALT_ID=3, TRACER_ID=4
      CHARACTER(LEN=255), INTENT(IN) :: FILENAME_NC
      TYPE (XPLEX), INTENT(IN)    :: DATA_NC(:,:,:,:)
      CHARACTER(LEN=255) :: EMS_AC_NAME(6)
      INTEGER :: N


      LON_DIM=IIPAR
      LAT_DIM=JJPAR
      ALT_DIM=LLPAR

      NCID=1
      EMS_AC_NAME(1)='AC_NOx'
      EMS_AC_NAME(2)='AC_HC'
      EMS_AC_NAME(3)='AC_PMNV'
      EMS_AC_NAME(4)='AC_PMFO'
      EMS_AC_NAME(5)='AC_CO'
      EMS_AC_NAME(6)='AC_SOx'

       CALL NC_CHECK(NF90_CREATE(TRIM(FILENAME_NC), NF90_CLOBBER ,NCID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Longitude',LON_DIM,LON_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Latitude',LAT_DIM,LAT_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Altitude',ALT_DIM,ALT_ID))
       NC_DIM_ID = (/LON_ID,LAT_ID,ALT_ID/)
      DO N = 1, TRACER_DIM
       CALL NC_CHECK(NF90_DEF_VAR(NCID, TRIM(EMS_AC_NAME(N)),
     &              NF90_DOUBLE, NC_DIM_ID, VAR_ID(N)))
      ENDDO
      CALL NC_CHECK(NF90_ENDDEF(NCID))
      DO N = 1, TRACER_DIM
      CALL NC_CHECK(NF90_PUT_VAR(NCID,VAR_ID(N),(DATA_NC(:,:,:,N)%i)))
      ENDDO
      CALL NC_CHECK(NF90_CLOSE(NCID))

      END SUBROUTINE SAVE_EMS_AC_FILE_X
!-----------------------------------------------------------------------------
      SUBROUTINE SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, DATA_NC)
! saving in  netcdf format
      USE NETCDF 
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE COMODE_MOD,     ONLY : JLOP
!      USE TRACER_MOD,    ONLY : N_TRACERS, TRACER_NAME

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h" 

!ITLOOP, IGAS, NAMEGAS

      INTEGER :: NCID
      INTEGER :: NC_DIM_ID(3)
      INTEGER :: LON_DIM, LAT_DIM, ALT_DIM, SPECIE_DIM
      INTEGER :: VAR_ID(200)
      INTEGER :: LON_ID=1, LAT_ID=2, ALT_ID=3, SPECIE_ID=4
      CHARACTER(LEN=255), INTENT(IN) :: FILENAME_NC
      TYPE (XPLEX), INTENT(IN)    :: DATA_NC(:,:)
      TYPE (XPLEX), ALLOCATABLE     :: DATA_NC_4D(:,:,:,:)
      INTEGER :: I,J,L,N,JLOOP, AS
      
      SPECIE_DIM=100
      ALLOCATE( DATA_NC_4D( IIPAR, JJPAR, LLTROP, SPECIE_DIM ), STAT=AS)
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'DATA_NC_4D' )
      DATA_NC_4D = 0d0
      LON_DIM=IIPAR
      LAT_DIM=JJPAR
      ALT_DIM=LLTROP
     
      DO N = 1, SPECIE_DIM
         DO L = 1, LLTROP
            DO J = 1, JJPAR
               DO I = 1, IIPAR
                   JLOOP=JLOP(I,J,L)
                   IF (JLOOP > 0 ) THEN
                      DATA_NC_4D(I,J,L,N)=DATA_NC(JLOOP,N)
                   ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO 

      NCID=1

       CALL NC_CHECK(NF90_CREATE(TRIM(FILENAME_NC), NF90_CLOBBER ,NCID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Longitude',LON_DIM,LON_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Latitude',LAT_DIM,LAT_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Altitude',ALT_DIM,ALT_ID))
       NC_DIM_ID = (/LON_ID,LAT_ID,ALT_ID/)

      DO N = 1, SPECIE_DIM 
       PRINT *,   NAMEGAS(N)
       IF ( TRIM(NAMEGAS(N)) == '' ) CYCLE
       CALL NC_CHECK(NF90_DEF_VAR(NCID, TRIM(NAMEGAS(N)),  
     &                 NF90_DOUBLE, NC_DIM_ID,VAR_ID(N)))
      ENDDO
      CALL NC_CHECK(NF90_ENDDEF(NCID))
      DO N = 1, SPECIE_DIM 
          CALL NC_CHECK(NF90_PUT_VAR(NCID,                            
     &                           VAR_ID(N),(DATA_NC_4D(:,:,:,N)%r)))
      ENDDO
      CALL NC_CHECK(NF90_CLOSE(NCID))
      DEALLOCATE( DATA_NC_4D )
      END SUBROUTINE SAVE_CSPEC_NETCDF_FILE

!----------------------------------------------------------------------

      SUBROUTINE SAVE_NETCDF_FILE(FILENAME_NC, DATA_NC)
! saving in  netcdf format
      USE NETCDF
      USE TRACER_MOD,    ONLY : N_TRACERS, TRACER_NAME
#     include "CMN_SIZE"

      INTEGER :: NCID
      INTEGER :: NC_DIM_ID(3)
      INTEGER :: LON_DIM, LAT_DIM, ALT_DIM, TRACER_DIM
      INTEGER :: VAR_ID(100)
      INTEGER :: LON_ID=1, LAT_ID=2, ALT_ID=3, TRACER_ID=4
      CHARACTER(LEN=255), INTENT(IN) :: FILENAME_NC
      TYPE (XPLEX), INTENT(IN)    :: DATA_NC(:,:,:,:)
      INTEGER :: N

      TRACER_DIM=N_TRACERS

      LON_DIM=IIPAR
      LAT_DIM=JJPAR
      ALT_DIM=LLPAR
  
      NCID=1

       CALL NC_CHECK(NF90_CREATE(TRIM(FILENAME_NC), NF90_CLOBBER ,NCID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Longitude',LON_DIM,LON_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Latitude',LAT_DIM,LAT_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Altitude',ALT_DIM,ALT_ID))
       NC_DIM_ID = (/LON_ID,LAT_ID,ALT_ID/)
      DO N = 1, TRACER_DIM 
       CALL NC_CHECK(NF90_DEF_VAR(NCID, TRIM(TRACER_NAME(N)),  
     &                NF90_DOUBLE, NC_DIM_ID,VAR_ID(N)))
      ENDDO
      CALL NC_CHECK(NF90_ENDDEF(NCID))
      DO N = 1, TRACER_DIM 
      CALL NC_CHECK(NF90_PUT_VAR(NCID,VAR_ID(N),(DATA_NC(:,:,:,N)%r)))
      ENDDO
      CALL NC_CHECK(NF90_CLOSE(NCID))

      END SUBROUTINE SAVE_NETCDF_FILE

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE SAVE_NETCDF_FILE_X(FILENAME_NC, DATA_NC)
! saving in  netcdf format
      USE NETCDF
      USE TRACER_MOD,    ONLY : N_TRACERS, TRACER_NAME
#     include "CMN_SIZE"

      INTEGER :: NCID
      INTEGER :: NC_DIM_ID(3)
      INTEGER :: LON_DIM, LAT_DIM, ALT_DIM, TRACER_DIM
      INTEGER :: VAR_ID(100)
      INTEGER :: LON_ID=1, LAT_ID=2, ALT_ID=3, TRACER_ID=4
      CHARACTER(LEN=255), INTENT(IN) :: FILENAME_NC
      TYPE (XPLEX), INTENT(IN)    :: DATA_NC(:,:,:,:)
      INTEGER :: N

      TRACER_DIM=N_TRACERS

      LON_DIM=IIPAR
      LAT_DIM=JJPAR
      ALT_DIM=LLPAR

      NCID=1

       CALL NC_CHECK(NF90_CREATE(TRIM(FILENAME_NC), NF90_CLOBBER ,NCID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Longitude',LON_DIM,LON_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Latitude',LAT_DIM,LAT_ID))
       CALL NC_CHECK(NF90_DEF_DIM(NCID,'Altitude',ALT_DIM,ALT_ID))
       NC_DIM_ID = (/LON_ID,LAT_ID,ALT_ID/)
      DO N = 1, TRACER_DIM
       CALL NC_CHECK(NF90_DEF_VAR(NCID, TRIM(TRACER_NAME(N)),
     &                NF90_DOUBLE, NC_DIM_ID,VAR_ID(N)))
      ENDDO
      CALL NC_CHECK(NF90_ENDDEF(NCID))
      DO N = 1, TRACER_DIM
      CALL NC_CHECK(NF90_PUT_VAR(NCID,VAR_ID(N),(DATA_NC(:,:,:,N)%i)))
      ENDDO
      CALL NC_CHECK(NF90_CLOSE(NCID))

      END SUBROUTINE SAVE_NETCDF_FILE_X

!----------------------------------------------------------------------
!**********************************************************************
!  Subroutine NC_CHECK checks the netcdf status (jkoo, 03/16/12)
!*********************************************************************
      SUBROUTINE NC_CHECK(status)
      USE NETCDF
         integer, intent ( in) :: status

         if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
      END SUBROUTINE NC_CHECK
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_PARNER_SAVE
!
!******************************************************************************
!  Subroutine CLEANUP_AIRCRAFT deallocates module variables. (srhb, 08/27/09)
!  Added additional allocated variables (jkoo,03/02/09)
!******************************************************************************
!

      ! Return to calling program
      END SUBROUTINE CLEANUP_PARNER_SAVE

!-----------------------------------------------------------------------------
      ! End of module
      END MODULE PARTNER_SAVE_MOD
