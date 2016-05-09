!$Id: scia_ch4_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE SCIA_CH4_MOD
!  
!******************************************************************************
!  Module SCIA_CH4_MOD for SCIAMACHY CH4 observations. 
!  By kjw, added adj32_023 (dkh, 02/12/12) 
!  
!******************************************************************************
!
  
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
 
      ! Parameters
      INTEGER, PARAMETER           :: LLSCIA   = 12
      INTEGER, PARAMETER           :: MAXSCIA  = 50000
      TYPE (XPLEX), PARAMETER            :: ERR_FRAC = 0.015


      ! Record to store data from each TES obs
      TYPE SCIA_CH4_OBS 
         INTEGER                                 :: NYMD
         INTEGER                                 :: NHMS
         INTEGER                                 :: QFLAG
         INTEGER                                 :: TFLAG
         TYPE (XPLEX)                                  :: TIME
         TYPE (XPLEX)                                  :: XCH4
         TYPE (XPLEX),  DIMENSION(LLSCIA)              :: AVGKERNEL
         TYPE (XPLEX),  DIMENSION(LLSCIA)              :: PRESCEN
         TYPE (XPLEX),  DIMENSION(LLSCIA+1)            :: PRESEDGE
         TYPE (XPLEX),  DIMENSION(LLSCIA)              :: PRIOR
         TYPE (XPLEX),  DIMENSION(50)                  :: GCII
         TYPE (XPLEX),  DIMENSION(50)                  :: GCJJ
         TYPE (XPLEX),  DIMENSION(50)                  :: GCFRAC
      ENDTYPE SCIA_CH4_OBS

      TYPE(SCIA_CH4_OBS)                         :: SCIA(MAXSCIA)



      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_SCIA_CH4_OBS( YYYYMMDD, NSCIA )
!
!******************************************************************************
!  Subroutine READ_SCIA_CH4_OBS reads the file and passes back info contained
!  therein. (kjw, 07/20/11) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME   (TYPE (XPLEX)) : SCIA observation filename to read
!
!  Arguments as Output: 
!  ============================================================================
!  (1 ) NSCIA     (INTEGER) : Number of SCIA retrievals for current day 
!
!  Module variable as Output: 
!  ============================================================================
!  (1 ) SCIA_CH4_OBS        : SCIA retrieval for current day 
!
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD,              ONLY : GET_RES_EXT
      USE TIME_MOD,               ONLY : GET_NYMD, GET_NHMS
      USE NETCDF_UTIL_MOD,        ONLY : NCDF_OPEN_FOR_READ
      USE NETCDF_UTIL_MOD,        ONLY : NCDF_GET_VARID
      USE NETCDF_UTIL_MOD,        ONLY : NCDF_GET_VAR
      USE NETCDF_UTIL_MOD,        ONLY : NCDF_CLOSE
      USE ERROR_MOD,              ONLY : ALLOC_ERR
      USE TIME_MOD,               ONLY : EXPAND_DATE

#     include "CMN_SIZE" 

      ! Arguments
      INTEGER,            INTENT(OUT) :: NSCIA
      INTEGER,            INTENT(IN)  :: YYYYMMDD

      ! Information to be stored in module varialbe SCIA_CH4_OBS 
      !TYPE (XPLEX)                          :: XCH4
      !TYPE (XPLEX)                          :: AVG_KERNEL(LLSCIA)
      !TYPE (XPLEX)                          :: PRES(LLSCIA)
      !TYPE (XPLEX)                          :: PRIOR(LLSCIA)
      !INTEGER                         :: QFLAG
      !TYPE (XPLEX)                          :: GCII(50)
      !TYPE (XPLEX)                          :: GCJJ(50)
      !TYPE (XPLEX)                          :: GCfrac(50)

      ! netCDF id's 
      INTEGER                         :: NCID
      INTEGER                         :: nobs_id, yyyymmdd_id, hhmmss_id
      INTEGER                         :: qflag_id, xch4_id, ch4ak_id
      INTEGER                         :: tflag_id
      INTEGER                         :: ch4presedge_id
      INTEGER                         :: ch4prescen_id, ch4prior_id
      INTEGER                         :: gcii_id, gcjj_id, gcfrac_id

      ! Arrays to hold info from NETCDF files
      INTEGER, ALLOCATABLE            :: Xqflag(:)
      INTEGER, ALLOCATABLE            :: Xtflag(:)
      INTEGER, ALLOCATABLE            :: Xnhms(:)
      INTEGER, ALLOCATABLE            :: Xnymd(:)
      TYPE (XPLEX), ALLOCATABLE             :: Xxch4(:)
      TYPE (XPLEX), ALLOCATABLE             :: Xch4ak(:,:)
      TYPE (XPLEX), ALLOCATABLE             :: Xch4prescen(:,:)
      TYPE (XPLEX), ALLOCATABLE             :: Xch4presedge(:,:)
      TYPE (XPLEX), ALLOCATABLE             :: Xch4prior(:,:)
      INTEGER, ALLOCATABLE            :: Xgcii(:,:)
      INTEGER, ALLOCATABLE            :: Xgcjj(:,:)
      TYPE (XPLEX), ALLOCATABLE             :: Xgcfrac(:,:)

      ! Loop indexes, and error handling.
      LOGICAL                         :: file_exist
      INTEGER                         :: NT, NB, AS, NGCFRAC
      INTEGER                         :: HH, MM, SS, NG, LS
      TYPE (XPLEX)                          :: frac
    
      ! Local variables
      CHARACTER(LEN=255)              :: READ_FILENAME



      !=================================================================
      ! READ_SCIA_CH4_OBS begins here!
      !=================================================================

      ! Construct complete filename 
      READ_FILENAME = TRIM( '/home/kjw/scia/data/imapv55/netcdf/ ' ) //
     &                TRIM( 'YYYY/MM/' )  // 
     &                TRIM( 'SCIA_CH4_YYYYMMDD.nc' ) 
      CALL EXPAND_DATE( READ_FILENAME, GET_NYMD(), 0 )

      ! Determine if there are observations today 
      INQUIRE( FILE=READ_FILENAME, exist=file_exist )

      ! If there is no observation file for this day,
      ! Return to calling program
      IF ( .not. file_exist ) THEN
         WRITE(6,*) '    - READ_SCIA_CH4_OBS: file does not exist: ',
     &                              TRIM( READ_FILENAME )
         WRITE(6,*) '                        no observations today.'

         ! Set NSCIA = 0 and Return to calling program
         NSCIA = 0
         RETURN
      ENDIF



      WRITE(6,*) '   - READ_SCIA_CH4_OBS: reading file: ', READ_FILENAME

      ! Open file and assign file id (FID)
      CALL NCDF_OPEN_FOR_READ( NCID, TRIM( READ_FILENAME ) )


      ! Get variable IDs for all variables to be read
      nobs_id      = ncdf_get_varid( NCID, 'Nobs' )
      yyyymmdd_id  = ncdf_get_varid( NCID, 'YYYYMMDD'   )
      hhmmss_id    = ncdf_get_varid( NCID, 'HHMMSS'     )
      qflag_id     = ncdf_get_varid( NCID, 'Qflag'      )
      xch4_id      = ncdf_get_varid( NCID, 'Xch4'       )
      ch4AK_id     = ncdf_get_varid( NCID, 'ch4AK'      )
      ch4presedge_id   = ncdf_get_varid( NCID, 'ch4presedge'   )
      ch4prescen_id   = ncdf_get_varid( NCID, 'ch4prescen'   )
      ch4prior_id  = ncdf_get_varid( NCID, 'ch4prior'  )
      IF ( GET_RES_EXT() .EQ. '4x5' ) THEN
         gcii_id   = ncdf_get_varid( NCID, 'GCII4'  )
         gcjj_id   = ncdf_get_varid( NCID, 'GCJJ4'  )
         gcfrac_id = ncdf_get_varid( NCID, 'GCfrac4'  )
         tflag_id  = ncdf_get_varid( NCID, 'Tflag4' )
         ngcfrac   = 8
      ELSE IF ( GET_RES_EXT() .EQ. '2x25' ) THEN
         gcii_id   = ncdf_get_varid( NCID, 'GCII2'  )
         gcjj_id   = ncdf_get_varid( NCID, 'GCJJ2'  )
         gcfrac_id = ncdf_get_varid( NCID, 'GCfrac2'  )
         tflag_id  = ncdf_get_varid( NCID, 'Tflag2' )
         ngcfrac   = 20
      ELSE IF ( GET_RES_EXT() .EQ. '05x0667' ) THEN
         gcii_id   = ncdf_get_varid( NCID, 'GCII05' )
         gcjj_id   = ncdf_get_varid( NCID, 'GCJJ05' )
         gcfrac_id = ncdf_get_varid( NCID, 'GCfrac05' )
         tflag_id  = ncdf_get_varid( NCID, 'Tflag05' )
         ngcfrac   = 50
      ENDIF


      ! Read Variables from NETCDF data file
         !print*,'     Reading NSCIA'
         ! ---- Number of observations for the day
         CALL NCDF_GET_VAR( NCID, nobs_id,   NSCIA   ) ! integer
         !print*,NSCIA
         !print*,'---------------------------------------------'

         ! ---- SCIAMACHY Quality Flag
         !print*,'     Reading QFLAG'
         ALLOCATE( Xqflag( NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xqflag(:) = -999d0
         CALL NCDF_GET_VAR( NCID, qflag_id,  Xqflag )  ! array of integers
         !print*,XQFLAG(1)
         !print*,'---------------------------------------------'

         ! ---- Tesselation Quality Flag
         !print*,'     Reading TFLAG'
         ALLOCATE( Xtflag( NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xtflag(:) = -999d0
         CALL NCDF_GET_VAR( NCID, tflag_id,  Xtflag )  ! array of integers
         !print*,XTFLAG(1)
         !print*,'---------------------------------------------'

         ! ---- Date of observation (YYYYMMDD)
         !print*,'     Reading YYYYMMDD'
         ALLOCATE( Xnymd( NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xnymd(:) = -999d0
         CALL NCDF_GET_VAR( NCID, yyyymmdd_id, Xnymd ) ! array of integers
         !print*,XNYMD(1)
         !print*,'---------------------------------------------'

         ! ---- Time of observation (HHMMSS)
         !print*,'     Reading HHMMSS'
         ALLOCATE( Xnhms( NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xnhms(:) = -999d0
         CALL NCDF_GET_VAR( NCID, hhmmss_id, Xnhms ) ! array of integers
         !print*,XNHMS(1)
         !print*,'---------------------------------------------'

         ! ---- SCIA CH4 volume mixing ratio [v/v]
         !print*,'     Reading XCH4'
         ALLOCATE( Xxch4( NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xxch4(:) = -999d0
         CALL NCDF_GET_VAR( NCID, xch4_id,   Xxch4   ) ! array of TYPE (XPLEX)
         !print*,XXCH4(1)
         !print*,'---------------------------------------------'

         ! ---- SCIA CH4 Averaging Kernel
         !print*,'     Reading CH4AK'
         ALLOCATE( Xch4ak( LLSCIA, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xch4ak(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, ch4AK_id,  Xch4ak  ) ! array of TYPE (XPLEX) x 12
         !print*,XCH4AK(:,1)
         !print*,'---------------------------------------------'

         ! ---- SCIA Pressure Centers [hPa]
         !print*,'     Reading CH4PRES Centers'
         ALLOCATE( Xch4prescen( LLSCIA, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xch4prescen(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, ch4prescen_id, Xch4prescen )
         ! array of TYPE (XPLEX) x 12
         !print*,XCH4PRESCEN(:,1)
         !print*,'---------------------------------------------'

         ! ---- SCIA Pressure Edges [hPa]
         !print*,'     Reading CH4PRES Edges'
         ALLOCATE( Xch4presedge( LLSCIA+1, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xch4presedge(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, ch4presedge_id, Xch4presedge )
         ! array of TYPE (XPLEX) x 13
         !print*,XCH4PRESEDGE(:,1)
         !print*,'---------------------------------------------'

         ! ---- SCIA CH4 Prior [v/v]
         !print*,'     Reading CH4PRIOR'
         ALLOCATE( Xch4prior( LLSCIA, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         Xch4prior(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, ch4prior_id, Xch4prior) ! array of TYPE (XPLEX) x 12
         !print*,XCH4PRIOR(:,1)
         !print*,'---------------------------------------------'

         ! ---- GEOS-Chem I indices of observation
         !print*,'     Reading GCII'
         ALLOCATE( Xgcii( NGCFRAC, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         XGCII(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, gcii_id, Xgcii ) ! array of integers x NGCFRAC
         !print*,XGCII(:,1)
         !print*,'---------------------------------------------'

         ! ---- GEOS-Chem J indices of observation
         !print*,'     Reading GCJJ'
         ALLOCATE( Xgcjj( NGCFRAC, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         XGCJJ(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, gcjj_id, Xgcjj ) ! array of integers x NGCFRAC
         !print*,XGCJJ(:,1)
         !print*,'---------------------------------------------'

         ! ---- Fraction of observation in each GEOS-Chem grid box
         !print*,'     Reading GCFRAC'
         ALLOCATE( Xgcfrac( NGCFRAC, NSCIA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
         XGCFRAC(:,:) = -999d0
         CALL NCDF_GET_VAR( NCID, gcfrac_id, Xgcfrac ) ! array of TYPE (XPLEX) x NGCFRAC
         !print*,XGCFRAC(:,1)
         !print*,'---------------------------------------------'


      ! Done reading variables from NETCDF. Close file
      CALL NCDF_CLOSE( NCID )

      !NT = 1
      !IF ( NT .EQ. 1 ) THEN
      !   WRITE( 6, * ) 'NSCIA = ',    NSCIA
      !   WRITE( 6, * ) 'QFLAG = ',    XQFLAG(NT)
      !   WRITE( 6, * ) 'NYMD = ',     XNYMD(NT)
      !   WRITE( 6, * ) 'NHMS = ',     XNHMS(NT)
      !   WRITE( 6, * ) 'XCH4 = ',     XXCH4(NT)
      !   WRITE( 6, * ) 'CH4AK = ',    XCH4AK(:,NT)
      !   WRITE( 6, * ) 'CH4PRESCEN = ',  XCH4PRESCEN(:,NT)
      !   WRITE( 6, * ) 'CH4PRESEDGE = ',  XCH4PRESEDGE(:,NT)
      !   WRITE( 6, * ) 'CH4PRIOR = ', XCH4PRIOR(:,NT)
      !   WRITE( 6, * ) 'GCII = ',     XGCII(:,NT)
      !   WRITE( 6, * ) 'GCJJ = ',     XGCJJ(:,NT)
      !   WRITE( 6, * ) 'GCFRAC = ',   XGCFRAC(:,NT)
      !ENDIF
      print*,'Xqflag #good= ',count(Xqflag .gt. 0)

      ! Assign variable output to module variable SCIA_CH4_OBS
      DO NT=1,NSCIA
         ! First initialize variables
         SCIA(NT)%qflag     = -999.
         SCIA(NT)%tflag     = -999.
         SCIA(NT)%nymd      = -999.
         SCIA(NT)%nhms      = -999.
         SCIA(NT)%xch4      = -999.
         DO LS=1,LLSCIA
            SCIA(NT)%avgkernel(LS) = -999.
            SCIA(NT)%prescen(LS)   = -999.
            SCIA(NT)%presedge(LS)  = -999.
            SCIA(NT)%prior(LS)     = -999.
         ENDDO
         SCIA(NT)%presedge(13)     = -999.
         DO NG=1,NGCFRAC
            SCIA(NT)%gcii(NG)      = -999.
            SCIA(NT)%gcjj(NG)      = -999.
            SCIA(NT)%gcfrac(NG)    = -999.
         ENDDO

         ! Place variables into SCIA structure
         SCIA(NT)%qflag     = Xqflag(NT)
         SCIA(NT)%tflag     = Xtflag(NT)
         SCIA(NT)%nymd      = Xnymd(NT)
         SCIA(NT)%nhms      = Xnhms(NT)
         SCIA(NT)%xch4      = Xxch4(NT)
         SCIA(NT)%avgkernel = Xch4ak(:,NT)
         SCIA(NT)%prescen   = Xch4prescen(:,NT)
         SCIA(NT)%presedge  = Xch4presedge(:,NT)
         SCIA(NT)%prior     = Xch4prior(:,NT)
         SCIA(NT)%gcii(1:NGCFRAC)      = Xgcii(1:NGCFRAC,NT)
         SCIA(NT)%gcjj(1:NGCFRAC)      = Xgcjj(1:NGCFRAC,NT)
         SCIA(NT)%gcfrac(1:NGCFRAC)    = Xgcfrac(1:NGCFRAC,NT)
      ENDDO


      ! Calculate fraction of day from NHMS
      DO NT=1,NSCIA
         HH = 0
         MM = 0
         SS = 0
         HH = floor(   SCIA(NT)%nhms                     / 1d4 )
         MM = floor( ( SCIA(NT)%nhms - 1d4*HH          ) / 1d2 )
         SS = floor(   SCIA(NT)%nhms - 1d4*HH - 1d2*MM         )
         frac =        HH / 24d0                 +       
     &                 MM / (24d0 * 60d0)        +
     &                 SS / (24d0 * 60d0 * 60d0)
         SCIA(NT)%TIME = frac
         !IF (NT .eq. 1 ) then
         !   print*,'nhms = ', scia(nt)%nhms
         !   print*,' hh  = ', hh
         !   print*,' mm  = ', mm
         !   print*,' ss  = ', ss
         !   print*,'frac = ', frac
         !endif
      ENDDO


      ! Cleanup allocated arrays
      IF ( ALLOCATED( Xqflag    ) ) DEALLOCATE( Xqflag    )
      IF ( ALLOCATED( Xnymd     ) ) DEALLOCATE( Xnymd     )
      IF ( ALLOCATED( Xnhms     ) ) DEALLOCATE( Xnhms     )
      IF ( ALLOCATED( Xxch4     ) ) DEALLOCATE( Xxch4     )
      IF ( ALLOCATED( Xch4ak    ) ) DEALLOCATE( Xch4ak    )
      IF ( ALLOCATED( Xch4prescen  ) ) DEALLOCATE( Xch4prescen  )
      IF ( ALLOCATED( Xch4presedge  ) ) DEALLOCATE( Xch4presedge  )
      IF ( ALLOCATED( Xch4prior ) ) DEALLOCATE( Xch4prior )
      IF ( ALLOCATED( Xgcii     ) ) DEALLOCATE( Xgcii     )
      IF ( ALLOCATED( Xgcjj     ) ) DEALLOCATE( Xgcjj     )
      IF ( ALLOCATED( Xgcfrac   ) ) DEALLOCATE( Xgcfrac   )


    
      ! Return to calling program
      END SUBROUTINE READ_SCIA_CH4_OBS
!------------------------------------------------------------------------------


      SUBROUTINE CALC_SCIA_CH4_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_SCIA_CH4_FORCE calculates the adjoint forcing from the SCIA
!  CH4 observations and updates the cost function. (kjw, 07/20/11)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (TYPE (XPLEX)) : Cost funciton                        [unitless]
!     
!     
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AD, TROPP
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
      USE TRACER_MOD,         ONLY : STT
      USE TRACER_MOD,         ONLY : TCVV
      USE TRACER_MOD,         ONLY : XNUMOLAIR, XNUMOL
      USE ERROR_MOD,          ONLY : ERROR_STOP

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER                     :: NT, LG, LS, I, J, NB
      INTEGER, SAVE               :: NSCIA
      INTEGER                     :: NTSTART, NTSTOP, nboxes
      TYPE (XPLEX)                      :: GC_PCENTER(LLPAR)
      TYPE (XPLEX)                      :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: thispcen(LLPAR)
      TYPE (XPLEX)                      :: thispedg(LLPAR)
      TYPE (XPLEX)                      :: thisad(LLPAR)
      TYPE (XPLEX)                      :: thisad1
      TYPE (XPLEX)                      :: thisch4(LLPAR)
      TYPE (XPLEX)                      :: GRIDMAP(LLPAR,LLSCIA)
      TYPE (XPLEX)                      :: GC_CH4_onSCIA(LLSCIA)
      TYPE (XPLEX)                      :: molec_air_onSCIA(LLSCIA)
      TYPE (XPLEX)                      :: CH4_PRIOR(LLSCIA)
      TYPE (XPLEX)                      :: frac, frac_total
      TYPE (XPLEX)                      :: thistrop, GC_TROP
      TYPE (XPLEX)                      :: fracreplace
      TYPE (XPLEX)                 :: mass_air, mole_air, molec_air_total
      TYPE (XPLEX)                      :: LHS, RHS, GC_XCH4, GC_CH4
      TYPE (XPLEX)                      :: Sobs, DIFF, FORCE
      TYPE (XPLEX)                      :: thisforce(LLPAR)
      TYPE (XPLEX)                  :: GC_XCH4_ADJ, DIFF_ADJ, GC_CH4_ADJ
      TYPE (XPLEX)                      :: GC_CH4_onSCIA_ADJ(LLSCIA)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: NEW_COST(MAXSCIA)
      TYPE (XPLEX)                      :: TIME_FRAC(MAXSCIA)
      TYPE (XPLEX)                      :: OLD_COST
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      LOGICAL, SAVE               :: DO_FDTEST = .TRUE.
      INTEGER                     :: IOS, ncount
      CHARACTER(LEN=255)          :: FILENAME

      ! Variables for FD testing
      TYPE (XPLEX)                      :: cost_func_pos, cost_func_neg
      TYPE (XPLEX)                      :: cost_func_0
      TYPE (XPLEX)                      :: PERT(LLPAR)
      TYPE (XPLEX)                      :: ADJ_SAVE(LLPAR, 50)
      TYPE (XPLEX)                      :: ADJ(LLPAR, 50)
      TYPE (XPLEX)                      :: FD_CEN(LLPAR, 50)
      TYPE (XPLEX)                      :: FD_POS(LLPAR, 50)
      TYPE (XPLEX)                      :: FD_NEG(LLPAR, 50)


      !=================================================================
      ! CALC_SCIA_CH4_FORCE begins here!
      !=================================================================

      print*, '     - CALC_SCIA_CH4_FORCE '

      NEW_COST(:) = 0d0 


      ! Open files for output
      IF ( FIRST ) THEN
         FILENAME = 'pres.NN.m' 
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_nh3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 102,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'tes_nh3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 103,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'apriori.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 104,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'diff.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 105,      FILE=TRIM( FILENAME    ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'force.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 106,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'nt_ll.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 107,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'adj_nh3_pert.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 108,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'adj_gc_nh3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 109,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'exp_nh3_hat.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 110,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'exp_nh3_hat_dbl.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 111,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         !kjw for testing adjoint of obs operator
         FILENAME = 'test_adjoint_obs.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 116,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FIRST = .FALSE.   ! only open files on first call to 
      ENDIF


      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Check if it is the last hour of a day.
      ! If so, read SCIA CH4 observations for the day
      IF ( GET_NHMS() == 230000 ) THEN 
 
         ! Read the SCIA CH4 file for this day 
         CALL READ_SCIA_CH4_OBS( GET_NYMD(), NSCIA )

         ! If NTES = 0, it means there are no observations today.
         ! Return to calling procedure
         IF ( NSCIA == 0 ) THEN
            WRITE(6,*) '   No SCIA CH4 obs today. Returning 01 ... '
            RETURN
         ENDIF

      ENDIF


      ! If here and NSCIA = 0, there are no more observations today.
      !  There were some, but they've been processed already.
      ! Return to calling procedure
      IF ( NSCIA == 0 ) THEN
         WRITE(6,*) '   No more SCIA CH4 obs today. Returning 02 ... '
         RETURN
      ENDIF
         

      ! Get the range of SCIA retrievals to assimilate in the current hour
      TIME_FRAC(1:NSCIA) = SCIA(1:NSCIA)%TIME
      CALL GET_NT_RANGE( NSCIA, GET_NHMS(), TIME_FRAC, 
     &                   NTSTART, NTSTOP ) 


      ! If no SCIA CH4 observations during this hour, return
      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 
         print*, ' No matching SCIA CH4 obs for this hour'
        RETURN
      ENDIF 


      ! Begin counting number of observations processed in this time step
      ncount = 0


!kjw DO NOT write satellite diagnostic file. It will take up too much space
!  for SCIA assimilations
!      ! Open file for this hour's satellite diagnostics
!      FILENAME = 'diag_sat.YYYYMMDD.hhmm.NN'
!      CALL EXPAND_NAME( FILENAME, N_CALC )
!      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
!      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!      OPEN( IU_FILE,      FILE=TRIM( FILENAME ),  STATUS='UNKNOWN',
!     &       IOSTAT=IOS,  FORM='FORMATTED',       ACCESS='SEQUENTIAL' )
!kjw


! need to update this in order to do i/o with this loop parallel 
!      ! Now do a parallel loop for analyzing data 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, ncount, nboxes, NB, LG, I, J, LS  )
!!$OMP+PRIVATE( LHS, RHS, GC_CH4, GC_XCH4             )
!!$OMP+PRIVATE( frac_total, frac, thistrop, thispcen  )
!!$OMP+PRIVATE( fracreplace, molec_air_onSCIA, molec_air_total  )
!!$OMP+PRIVATE( thispedg, thisad, thisch4, thisforce  )
!!$OMP+PRIVATE( GC_PCENTER, GC_PEDGE, GC_CH4_NATIVE   )
!!$OMP+PRIVATE( GC_TROP, GRIDMAP, GC_CH4_onSCIA       )
!!$OMP+PRIVATE( Sobs, DIFF, FORCE, thisad1 )
!!$OMP+PRIVATE( DIFF_ADJ, GC_XCH4_ADJ, GC_CH4_ADJ     )
!!$OMP+PRIVATE( GC_CH4_onSCIA_ADJ, GC_CH4_NATIVE_ADJ  )

      DO NT  = NTSTART, NTSTOP, -1
!      DO NT = NTSTART,NTSTART-10,-1
!      DO NT = 8776, 8776
         ! Check quality of retrieval 
         IF ( ( SCIA(NT)%QFLAG .ne. 1 ) .OR. 
     &        ( SCIA(NT)%TFLAG .ne. 1 ) ) THEN
            !print*, ' SKIPPING record ', NT
            !print*, ' QFLAG = ', SCIA(NT)%QFLAG
            CYCLE
         ENDIF

         print*, '     - CALC_SCIA_CH4_FORCE: analyzing record ', NT

         ! Count this observation
         ncount = ncount + 1


         ! For safety, initialize these up to LLSCIA 
         CH4_PRIOR(:)     = 0d0
         GC_CH4_NATIVE(:) = 0d0 
         GRIDMAP(:,:)     = 0d0


         ! Get GEOS-Chem pressure and CH4 column corresponding to SCIA
         !    observation. This will not be from a single grid box but 
         !    rather from many as determined by GCII, GCJJ and GCFRAC.
         ! CH4 in [v/v] and pressure in [hPa]

            ! Initialize
            GC_PCENTER(:)    = 0d0
            GC_PEDGE(:)      = 0d0
            frac_total       = 0d0


            ! Determine number of GEOS-Chem boxes covered by the observation
            nboxes = count( SCIA(NT)%GCfrac(:) .gt. 0.0 )

            ! Loop over boxes
            DO NB=1,nboxes

               ! Clear variables to be safe
               I           = 0
               J           = 0
               frac        = 0d0
               thispcen(:) = 0d0
               thispedg(:) = 0d0
               thisad(:)   = 0d0
               thisch4(:)  = 0d0
               thistrop    = 0d0


               ! I and J indices and fractional influence of this box
               I    = SCIA(NT)%GCII(NB)
               J    = SCIA(NT)%GCJJ(NB)
               frac = SCIA(NT)%GCfrac(NB)
               thistrop = TROPP(I,J)

               ! Get column of pressure centers and CH4 values
                  DO LG=1,LLPAR

                     ! Pressure centers [hPa]
                     thispcen(LG) = GET_PCENTER(I,J,LG)

                     ! Pressure edges [hPa]
                     thispedg(LG) = GET_PEDGE(I,J,LG)

                     ! mass per box [kg]
                     thisad(LG)   = AD(I,J,LG)
                     
                  ENDDO
            
                  ! CH4 [kg/box] --> [v/v]
                  !    Numerator   = moles CH4/box
                  !    Denominator = moles air/box
                  thisch4(:) = ( CHK_STT(I,J,:,1 ) * XNUMOL(1) ) / 
     &                         ( thisad(:) * XNUMOLAIR )
            
               ! Add pressure and ch4 columns to total
               GC_PCENTER(:) = GC_PCENTER(:) + thispcen(:) * frac
               GC_PEDGE(:)   = GC_PEDGE(:)   + thispedg(:) * frac
               GC_CH4_NATIVE(:) = GC_CH4_NATIVE(:) + thisch4(:)  * frac
               GC_TROP = thistrop * frac
            
               ! Error checking. sum of fracs should = 1
               frac_total = frac_total + frac
            
            ENDDO
            !print*,'SCIA(NT)%GCFLAG = ',SCIA(NT)%GCFRAC(1:50)
            !print*,'frac_total = ',frac_total

            ! Error checking. sum of fracs should = 1 within reason
            IF ( abs(frac_total-1) .gt. 1d-5 ) THEN
               WRITE( 6, * ) 'ERROR in CALC_SCIA_CH4_FORCE: '
               CALL ERROR_STOP( 'fractions /= 1','GET_GC_PROFILE' )
            ENDIF

         ! Done constructing representative GEOS-Chem profile from 
         !  multiple grid boxes

!         ! dkh debug: compare profiles:
!         print*, ' GC_PCENTER, GC_CH4_NATIVE in [v/v] '
!         WRITE(6,100) (GC_PCENTER(LG), GC_CH4_NATIVE(LG), 
!     &                   LG = LLPAR, 1, -1 )


         ! Get interpolation matrix that maps GEOS-Chem to SCIAMACHY grid
         ! GEOS-Chem grid now in [v/v]
         GRIDMAP(1:LLPAR, 1:LLSCIA) = GET_INTMAP( GC_PEDGE,         
     &                                        SCIA(NT)%PRESEDGE(:)   )

         ! Interpolate GEOS-Chem CH4 column [v/v] to SCIA grid [v/v]
         DO LS = 1, LLSCIA
            GC_CH4_onSCIA(LS) = 0d0 
            DO LG = 1, LLPAR
               GC_CH4_onSCIA(LS) = GC_CH4_onSCIA(LS) 
     &                    + GRIDMAP(LG,LS) * GC_CH4_NATIVE(LG) 
            ENDDO
         ENDDO
!         print*,'GRIDMAP = ',GRIDMAP

         CH4_PRIOR(:) = SCIA(NT)%PRIOR
         !print*, ' SCIA_PRES, GC_CH4_onSCIA [v/v], SCIA_PRIOR  ' 
!         WRITE(6,101) ( SCIA(NT)%PRES(LS), GC_CH4_onSCIA(LS),
!     &                  CH4_PRIOR(LS), LS,    LS = LLSCIA, 1, -1 ) 



         ! Replace GEOS-Chem stratosphere with SCIAMACHY a priori strat
         DO LS=1,LLSCIA

            ! If tropopause pressure less than upper box edge, continue
            IF ( GC_TROP .lt. SCIA(NT)%PRESEDGE(LS+1) )  CONTINUE 

            ! If trop pressure greater than lower box edge,
            !   replace entire box with prior values
            IF ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS) ) THEN
               GC_CH4_onSCIA(LS) = CH4_PRIOR(LS)
            ENDIF

            ! If trop pressure within grid box, replace fraction of value
            IF ( ( GC_TROP .lt. SCIA(NT)%PRESEDGE(LS) )   .AND.  
     &           ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS+1) ) ) THEN
               fracreplace = ( GC_TROP - SCIA(NT)%PRESEDGE(LS+1) ) / 
     &               ( SCIA(NT)%PRESEDGE(LS) - SCIA(NT)%PRESEDGE(LS+1) )
               GC_CH4_onSCIA(LS) = (1-fracreplace) * GC_CH4_onSCIA(LS) + 
     &                             fracreplace * CH4_PRIOR(LS)
            ENDIF
         ENDDO



         ! Convert [v/v] --> [molec/cm2] for application of SCIA AK.
         molec_air_onSCIA(:) = 0d0
         CH4_PRIOR(:) = SCIA(NT)%PRIOR
         DO LS=1,LLSCIA

            ! Get molecules / cm2 of air in each pressure level = molec_air
            !   F=ma where F in one square meter is dPressure
            ! molec/cm2 of air in column
            molec_air_onSCIA(LS) = 
     &           ( SCIA(NT)%PRESEDGE(LS) - SCIA(NT)%PRESEDGE(LS+1) )
     &            * 1d2 / 9.8 * 1d-4 * 1d3 * (1/28.9644) * 6.022d23

            ! CH4 [molec/cm2] = CH4 [v/v] * total_air [molec/cm2]
            CH4_PRIOR(LS)     = CH4_PRIOR(LS)     * molec_air_onSCIA(LS)
            GC_CH4_onSCIA(LS) = GC_CH4_onSCIA(LS) * molec_air_onSCIA(LS)

         ENDDO
         



!         ! dkh debug: compare profiles:
!         print*, ' GC_PCENTER, GC_CH4_NATIVE in [v/v]'
!         WRITE(6,100) (GC_PEDGE(LG), GC_CH4_NATIVE(LG), 
!     &                   LG = LLPAR, 1, -1 )
!         print*, ' SCIA_PRES, GC_CH4_onSCIA [molec/cm2], SCIA_PRIOR  ' 
!         WRITE(6,101) ( SCIA(NT)%PRES(LS), GC_CH4_onSCIA(LS),
!     &                  CH4_PRIOR(LS), LS,    LS = LLSCIA, 1, -1 ) 
! 100  FORMAT(1X,F16.8,1X,E24.12)
!         print*,'total GC on SCIA molec/cm2 = ',SUM(GC_CH4_onSCIA)


!         !--------------------------------------------------------------
!         ! Apply SCIA observation operator
!         !
!         !   x_hat = A ( x_m ) + ( 1 - A ) x_a 
!         !  
!         !  where  
!         !    x_hat = GC modeled column as seen by SCIA [molec/cm2]
!         !    x_a   = SCIA apriori column               [molec/cm2]
!         !    x_m   = GC modeled column                 [molec/cm2]
!         !    A     = SCIA averaging kernel 
!         !--------------------------------------------------------------

         !--------------------------------------------------------------
         ! Apply SCIA observation operator
         !
         !   x_hat = A ( x_m ) + ( 1 - A ) x_a 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by SCIA [v/v]
         !    x_a   = SCIA apriori column               [v/v]
         !    x_m   = GC modeled column                 [v/v]
         !    A     = SCIA averaging kernel 
         !--------------------------------------------------------------


         ! A ( x_m )
         LHS = 0d0
         DO LS = 1, LLSCIA
           LHS = LHS + SCIA(NT)%AVGKERNEL(LS)
     &                       * GC_CH4_onSCIA(LS)
         ENDDO
     
         ! ( 1 - A ) x_a
         RHS = 0d0
         DO LS = 1, LLSCIA
            RHS = RHS + ( ( 1 - SCIA(NT)%AVGKERNEL(LS) )
     &                              * CH4_PRIOR(LS)  )
         ENDDO

         ! x_hat = A ( x_m ) + ( 1 - A ) x_a
         GC_CH4 = RHS + LHS


         ! Convert Units from [molec/cm2] --> [v/v]
         ! Get molecules of air in column using F=ma
         molec_air_total = 0d0
         molec_air_total = SCIA(NT)%PRESEDGE(1) * 1d2 / 9.8 * 1d-4  ! air [kg/cm2]
     &                  * 1d3 * (1/28.9644) * 6.022d23     ! air [molec/cm2]

         ! [molec ch4 / cm2] --> [v/v]   = molec CH4 / molec air in one cm^2
         GC_XCH4 = GC_CH4 / molec_air_total
         !print*,'gc_xch4 = ',gc_xch4


         !--------------------------------------------------------------
         ! Calculate cost function, given S is observation error covariance matrix
         !     Sobs = 1x1 array [ (molec/cm2) ^2 ]
         ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------


         ! Calculate error on this day, 
         !     given in fractional terms as a module variable
         ! Sobs in [v/v]^2
         Sobs = ( ERR_FRAC * 1d-9 * SCIA(NT)%XCH4 ) **2


         ! Calculate difference between modeled and observed profile
         DIFF = GC_XCH4 - 1d-9 * SCIA(NT)%XCH4
         !print*,'NORMAL : DIFF',DIFF


         ! Calculate J(x) = DIFF^T * S_{obs}^{-1} * DIFF 
         NEW_COST(NT) = DIFF**2 / Sobs

         ! Calculate dJ/dx = 2 * DIFF * S_{obs}^{-1}
         FORCE = 0d0
         FORCE = 2 * DIFF / Sobs
         !print*,'NORMAL : FORCE',FORCE

         !print*,'gc_xch4 = ',gc_xch4
         !print*,'SCIA(NT)%XCH4 = ',1d-9 * SCIA(NT)%XCH4
         !print*,'diff =  ',diff
         !print*,'force = ',force
         !print*,'Sobs =  ',Sobs
         ! dkh debug: compare profiles:
         !print*, ' SCIA_PRIOR, XCH4_SCIA, XCH4_GC'
!         WRITE(6,101) (SCIA(NT)%PRIOR(LS), SCIA(NT)%XCH4, GC_XCH4,
!     &             LS,    LS = LLSCIA, 1, -1 )
! 101  FORMAT(1X,E24.16,1X,E24.16,1X,E24.16,1x,i3)
         WRITE(105,101) gc_xch4*1d9, SCIA(NT)%XCH4, 
     &               gc_xch4*1d9-SCIA(NT)%XCH4
 101     FORMAT(F15.8,5X,E15.8,5X,F15.8)

         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

         ! dkh debug
!         print*, 'DIFF , FORCE, Sobs ' 
!         WRITE(6,102) (DIFF, FORCE, Sobs)
! 102  FORMAT(1X,d14.6,1X,d14.6)


         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE


         ! Adjoint of GEOS-Chem - SCIAMACHY difference
         GC_XCH4_ADJ = DIFF_ADJ

         ! Adjoint of unit conversion from [molec/cm2] --> [v/v]
         GC_CH4_ADJ = GC_XCH4_ADJ / molec_air_total
         !print*,'NORMAL : GC_CH4_ADJ',GC_CH4_ADJ


         ! Adjoint of SCIA observation operator
         DO LS=1,LLSCIA
            GC_CH4_onSCIA_ADJ(LS) = SCIA(NT)%AVGKERNEL(LS) *
     &                                   GC_CH4_ADJ
         ENDDO
         !print*,'NORMAL : GC_CH4_ONSCIA_ADJ',GC_CH4_ONSCIA_ADJ

         ! Adjoint of unit conversion [v/v] --> [molec/cm2]
         DO LS=1,LLSCIA
            GC_CH4_onSCIA_ADJ(LS) = GC_CH4_onSCIA_ADJ(LS)
     &                                * molec_air_onSCIA(LS)
         ENDDO
         !print*,'NORMAL : GC_CH4_ONSCIA_ADJ',GC_CH4_ONSCIA_ADJ


         ! Adjoint of replacing GEOS-Chem stratosphere with SCIA prior
         DO LS=1,LLSCIA
            ! If trop pressure within grid box, replace fraction of value
            IF ( ( GC_TROP .lt. SCIA(NT)%PRESEDGE(LS) )   .AND.  
     &           ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS+1) ) ) THEN
               fracreplace = ( GC_TROP - SCIA(NT)%PRESEDGE(LS+1) ) / 
     &               ( SCIA(NT)%PRESEDGE(LS) - SCIA(NT)%PRESEDGE(LS+1) )
               GC_CH4_onSCIA_ADJ(LS) = 
     &                     (1-fracreplace) * GC_CH4_onSCIA_ADJ(LS)
            ENDIF

            ! If trop pressure gt lower grid box boundary, GC_CH4_onSCIA(LS)=0
            IF ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS) ) THEN
               GC_CH4_onSCIA_ADJ(LS) = 0d0
            ENDIF

         ENDDO


         ! Adjoint of interpolation
         DO LG=1,LLPAR
         DO LS=1,LLSCIA
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) + 
     &           GRIDMAP(LG,LS) * GC_CH4_onSCIA_ADJ(LS)
         ENDDO
         ENDDO
         !print*,'NORMAL : GC_CH4_NATIVE_ADJ',GC_CH4_NATIVE_ADJ




         ! Adjoint of GEOS-Chem column averaging
         !   Distribute adjoint forcing across GEOS-Chem grid boxes from 
         !   which the original GEOS-Chem column was calculated.
         ! Adjoint forcing added to adjoint tracer array in this subroutine

            ! Determine number of GEOS-Chem boxes covered by the observation
            nboxes = count( SCIA(NT)%GCfrac(:) .gt. 0.0 )
      
            ! Loop over boxes, placing adjoint variable into the STT_ADJ array
            DO NB=1,nboxes

               ! Clear variables to be safe
               I    = 0
               J    = 0
               frac = 0d0
      
               ! I and J indices and fractional influence of this box
               I    = SCIA(NT)%GCII(NB)
               J    = SCIA(NT)%GCJJ(NB)
               frac = SCIA(NT)%GCfrac(NB)
      
               ! Adjoint of unit conversion from [kg/box] to [v/v]
               DO LG=1,LLPAR

                  ! Get mass in this grid box
                  thisad1 = 0d0
                  thisad1 = AD(I,J,LG)
      
                  ! adjoint of unit conversion
                  thisforce(LG) = GC_CH4_NATIVE_ADJ(LG) * XNUMOL(1) / 
     &                                 ( thisad1 * XNUMOLAIR ) * frac
     
      
                  ! Place adjoint forcing back to adjoint array
                  STT_ADJ(I,J,LG,1) = STT_ADJ(I,J,LG,1) + thisforce(LG) 
      
      
               ENDDO
      
            ENDDO

         ! End distributing adjoint forcing to STT_ADJ array
         !   print*,'thisforce = ',thisforce



! -----------------------------------------------------------------------
!    Use this section to test the adjoint of the TES_CH4 operator by
!          slightly perturbing model [CH4] and recording resultant change
!          in calculated contribution to the cost function.
!
!    This routine will write the following information for each observation
!          to rundir/diagadj/test_adjoint_obs.NN.m
!
!    The adjoint of the observation operator has been tested and validated
!          as of 7/20/10, kjw.
!
!      !IF (( DO_FDTEST ) .AND. ( nboxes .gt. 1 )) THEN
      IF ( DO_FDTEST ) THEN
      WRITE(116,210) '  LG'        , ' box',       '  TROP',
     &               '     GC_PRES',
     &               '      FD_POS', '      FD_NEG', '      FD_CEN',
     &               '         ADJ', '    COST_POS', '    COST_NEG', 
     &               '  FD_POS/ADJ', '  FD_NEG/ADJ', '  FD_CEN/ADJ'
      PERT(:) = 0D0

      COST_FUNC_0 = 0d0
      CALL CALC_SCIA_CH4_FORCE_FD( COST_FUNC_0, PERT, ADJ, NT, NB )
      ADJ_SAVE(:,:) = ADJ(:,:)

      ! Write identifying information to top of satellite diagnostic file
      WRITE(116,212) 'GC_PSURF    ', GC_PEDGE(1)
      WRITE(116,212) 'SCIA PSURF  ', SCIA(NT)%presedge(1)
      WRITE(116,212) 'NEW_COST:   ', NEW_COST(NT)
      WRITE(116,212) 'COST_FUNC_0:', COST_FUNC_0

      ! Determine number of GEOS-Chem boxes covered by the observation
      nboxes = count( SCIA(NT)%GCfrac(:) .gt. 0.0 )

      ! Perform finite difference testing at each vertical level
      ! and for each horizontal grid box in this observation
      DO LG = 1, 47
      DO NB = 1, nboxes

         ! Positive perturbation to GEOS-Chem CH4 columns
         PERT(:) = 0.0
         PERT(LG) = 0.1
         COST_FUNC_pos = 0D0
         CALL CALC_SCIA_CH4_FORCE_FD( COST_FUNC_pos, PERT, ADJ, NT, NB )

         ! Negative perturbation to GEOS-Chem CH4 columns
         PERT(:) = 0.0
         PERT(LG) = -0.1
         COST_FUNC_neg = 0D0
         CALL CALC_SCIA_CH4_FORCE_FD( COST_FUNC_neg, PERT, ADJ, NT, NB )

         ! Calculate dJ/dCH4 from perturbations
         FD_CEN(LG,NB)   = ( COST_FUNC_pos - COST_FUNC_neg ) / 0.2d0
         FD_POS(LG,NB)   = ( COST_FUNC_pos - COST_FUNC_0 )   / 0.1d0
         FD_NEG(LG,NB)   = ( COST_FUNC_0 - COST_FUNC_neg )   / 0.1d0

         ! Write information to satellite diagnostic file
         WRITE(116, 211)  LG, NB,       GC_PCENTER(LG),
     &                    FD_POS(LG,NB),   FD_NEG(LG,NB),
     &                    FD_CEN(LG,NB),   ADJ_SAVE(LG,NB), 
     &                    COST_FUNC_pos, COST_FUNC_neg, 
     &                    FD_POS(LG,NB)/ADJ_SAVE(LG,NB),
     &                    FD_NEG(LG,NB)/ADJ_SAVE(LG,NB),
     &                    FD_CEN(LG,NB)/ADJ_SAVE(LG,NB)
      ENDDO
      ENDDO


      WRITE(116,'(a)') '----------------------------------------------'

 210  FORMAT(A4,2x,A4,A6,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,
     &       A12,2x,A12,2x,A12,2x,A12,2x)
 211  FORMAT(I4,2x,I4,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,
     &        2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6)
 212  FORMAT(A12,F22.6)
 213  FORMAT(A12,I4)
 214  FORMAT(I4,2x,F18.6,2x,F18.6)
! -----------------------------------------------------------------------
      DO_FDTEST = .FALSE.

      ENDIF  ! IF ( DO_FDTEST )




      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function 
      COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))

      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, ' SCIA contribution          = ', COST_FUNC - OLD_COST
      print*, ' # Good Observations analyzed = ', ncount
      print*, ' # Total Observations read    = ', NTSTART-NTSTOP


      ! Return to calling program
      END SUBROUTINE CALC_SCIA_CH4_FORCE

!------------------------------------------------------------------------------






      SUBROUTINE CALC_SCIA_CH4_FORCE_FD( COST_FUNC_A, PERT, ADJ, 
     &                                                NT, boxnum )
!
!******************************************************************************
!  Subroutine CALC_SCIA_CH4_FORCE calculates the adjoint forcing from the SCIA
!  CH4 observations and updates the cost function. (kjw, 07/20/11)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC_A (TYPE (XPLEX)) : Cost funciton (INOUT)                 [unitless]
!  (2 ) PERT        (TYPE (XPLEX)) : Array of perturbations to CH4 column (+/- 0.1, for ex.)
!  (5 ) ADJ         (TYPE (XPLEX)) : Array of adjoint forcings (OUT)
!  (3 ) NT         (INTEGER) : Observation number to process
!  (4 ) NB         (INTEGER) : Box number in which to make perturbation
!     
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AD, TROPP
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
      USE TRACER_MOD,         ONLY : STT
      USE TRACER_MOD,         ONLY : TCVV
      USE TRACER_MOD,         ONLY : XNUMOLAIR, XNUMOL
      USE ERROR_MOD,          ONLY : ERROR_STOP

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC_A
      TYPE (XPLEX), INTENT(OUT)         :: ADJ(LLPAR,50)
      TYPE (XPLEX), INTENT(IN)          :: PERT(LLPAR)
      INTEGER, INTENT(IN)         :: NT
      INTEGER, INTENT(IN)         :: boxnum


      ! Local variables
      INTEGER                     :: LG, LS, I, J, NB
      INTEGER                     :: NSCIA, nboxes
      INTEGER                     :: NTSTART, NTSTOP
      TYPE (XPLEX)                      :: GC_PCENTER(LLPAR)
      TYPE (XPLEX)                      :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: thispcen(LLPAR)
      TYPE (XPLEX)                      :: thispedg(LLPAR)
      TYPE (XPLEX)                      :: thisad(LLPAR)
      TYPE (XPLEX)                      :: thisad1
      TYPE (XPLEX)                      :: thisch4(LLPAR)
      TYPE (XPLEX)                      :: GRIDMAP(LLPAR,LLSCIA)
      TYPE (XPLEX)                      :: GC_CH4_onSCIA(LLSCIA)
      TYPE (XPLEX)                      :: molec_air_onSCIA(LLSCIA)
      TYPE (XPLEX)                      :: CH4_PRIOR(LLSCIA)
      TYPE (XPLEX)                      :: frac, frac_total
      TYPE (XPLEX)                      :: thistrop, GC_TROP
      TYPE (XPLEX)                      :: fracreplace
      TYPE (XPLEX)                 :: mass_air, mole_air, molec_air_total
      TYPE (XPLEX)                      :: LHS, RHS, GC_XCH4, GC_CH4
      TYPE (XPLEX)                      :: Sobs, DIFF, FORCE
      TYPE (XPLEX)                      :: thisforce(LLPAR)
      TYPE (XPLEX)                  :: GC_XCH4_ADJ, DIFF_ADJ, GC_CH4_ADJ
      TYPE (XPLEX)                      :: GC_CH4_onSCIA_ADJ(LLSCIA)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: NEW_COST(MAXSCIA)
      TYPE (XPLEX)                      :: OLD_COST
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME



      !=================================================================
      ! CALC_SCIA_CH4_FORCE_FD begins here!
      !=================================================================


      ! Initialize for safety
      ADJ(:,:)         = 0d0
      GC_CH4_NATIVE(:) = 0d0 
      GRIDMAP(:,:)     = 0d0


      ! Get GEOS-Chem pressure and CH4 column corresponding to SCIA
      !    observation. This will not be from a single grid box but 
      !    rather from many as determined by GCII, GCJJ and GCFRAC.
      ! CH4 in [v/v] and pressure in [hPa]

         ! Initialize
         GC_PCENTER(:)    = 0d0
         GC_PEDGE(:)      = 0d0
         GC_CH4_NATIVE(:) = 0d0
         frac_total       = 0d0


         ! Determine number of GEOS-Chem boxes covered by the observation
         nboxes = count( SCIA(NT)%GCfrac(:) .gt. 0.0 )

         ! Loop over boxes
         DO NB=1,nboxes
         
            ! Clear variables to be safe
            I           = 0
            J           = 0
            frac        = 0d0
            thispcen(:) = 0d0
            thispedg(:) = 0d0
            thisad(:)   = 0d0
            thisch4(:)  = 0d0
            thistrop    = 0d0
         
         
            ! I and J indices and fractional influence of this box
            I    = SCIA(NT)%GCII(NB)
            J    = SCIA(NT)%GCJJ(NB)
            frac = SCIA(NT)%GCfrac(NB)
            thistrop = TROPP(I,J)
         
            ! Get column of pressure centers and CH4 values
            DO LG=1,LLPAR
         
               ! Pressure centers [hPa]
               thispcen(LG) = GET_PCENTER(I,J,LG)
         
               ! Pressure edges [hPa]
               thispedg(LG) = GET_PEDGE(I,J,LG)
         
               ! mass per box [kg]
               thisad(LG)   = AD(I,J,LG)
               
            ENDDO
         
            ! CH4 [kg/box] --> [v/v]
            !    Numerator   = moles CH4/box
            !    Denominator = moles air/box
            ! Only perturb one box given by boxnum input
            IF ( NB .EQ. boxnum ) THEN 
               DO LG=1,LLPAR
               thisch4(LG) = ( CHK_STT(I,J,LG,1 ) * ( 1+PERT(LG) )
     &                       * XNUMOL(1) ) / ( thisad(LG) * XNUMOLAIR )
               ENDDO
            ELSE
               DO LG=1,LLPAR
               thisch4(LG) = ( CHK_STT(I,J,LG,1 ) 
     &                       * XNUMOL(1) ) / ( thisad(LG) * XNUMOLAIR )
               ENDDO
            ENDIF

            ! Add pressure and ch4 columns to total
            GC_PCENTER(:) = GC_PCENTER(:) + thispcen(:) * frac
            GC_PEDGE(:)   = GC_PEDGE(:)   + thispedg(:) * frac
            GC_CH4_NATIVE(:) = GC_CH4_NATIVE(:)     + thisch4(:)  * frac
            GC_TROP = thistrop * frac

            ! Error checking. sum of fracs should = 1
            frac_total = frac_total + frac
         
         ENDDO
         
         
         ! Error checking. sum of fracs should = 1
         IF ( abs(frac_total-1) .gt. 1d-5 ) THEN
            WRITE( 6, * ) 'ERROR in GET_GC_PROFILE: fractions /= 1'
            CALL ERROR_STOP( 'problem','GET_GC_PROFILE' )
         ENDIF

      ! Done getting representative GEOS-Chem profile from many grid boxes



         ! Get interpolation matrix that maps GEOS-Chem to SCIAMACHY grid
         ! GEOS-Chem grid now in [molec/m2]
         GRIDMAP(1:LLPAR, 1:LLSCIA) = GET_INTMAP( GC_PEDGE,         
     &                                        SCIA(NT)%PRESEDGE(:)    )



         ! Interpolate GEOS-Chem CH4 column to SCIA grid  [v/v] --> [v/v]
         DO LS = 1, LLSCIA
            GC_CH4_onSCIA(LS) = 0d0 
            DO LG = 1, LLPAR
               GC_CH4_onSCIA(LS) = GC_CH4_onSCIA(LS) 
     &                 + GRIDMAP(LG,LS) * GC_CH4_NATIVE(LG) 
            ENDDO
         ENDDO


         CH4_PRIOR(:) = SCIA(NT)%PRIOR
      ! Replace GEOS-Chem stratosphere with SCIAMACHY a priori strat
         DO LS=1,LLSCIA

            ! If tropopause pressure less than upper box edge, continue
            IF ( GC_TROP .lt. SCIA(NT)%PRESEDGE(LS+1) )  CONTINUE 

            ! If trop pressure greater than lower box edge,
            !   replace entire box with prior values
            IF ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS) ) THEN
               GC_CH4_onSCIA(LS) = CH4_PRIOR(LS)
            ENDIF

            ! If trop pressure within grid box, replace fraction of value
            IF ( ( GC_TROP .lt. SCIA(NT)%PRESEDGE(LS) )   .AND.  
     &           ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS+1) ) ) THEN
               fracreplace = ( GC_TROP - SCIA(NT)%PRESEDGE(LS+1) ) / 
     &               ( SCIA(NT)%PRESEDGE(LS) - SCIA(NT)%PRESEDGE(LS+1) )
               GC_CH4_onSCIA(LS) = (1-fracreplace) * GC_CH4_onSCIA(LS) + 
     &                             fracreplace * CH4_PRIOR(LS)
            ENDIF
         ENDDO


         ! Convert [v/v] --> [molec/cm2] for application of SCIA AK.
         molec_air_onSCIA(:) = 0d0
         CH4_PRIOR(:) = SCIA(NT)%PRIOR
         DO LS=1,LLSCIA

            ! Get molecules / cm2 of air in each pressure level = molec_air
            !   F=ma where F in one square meter is dPressure
            ! molec/cm2 of air in column
            molec_air_onSCIA(LS) = 
     &              ( SCIA(NT)%PRESEDGE(LS) - SCIA(NT)%PRESEDGE(LS+1) )
     &             * 1d2 / 9.8 * 1d-4 * 1d3 * (1/28.9644) * 6.022d23

            ! CH4 [molec/cm2] = CH4 [v/v] * total_air [molec/cm2]
            CH4_PRIOR(LS)     = CH4_PRIOR(LS)     * molec_air_onSCIA(LS)
            GC_CH4_onSCIA(LS) = GC_CH4_onSCIA(LS) * molec_air_onSCIA(LS)

         ENDDO

        ! dkh debug: compare profiles:
!         print*, ' GC_PRES, GC_CH4_NATIVE '
!         WRITE(6,100) (GC_PRES(LG), GC_CH4_NATIVE(LG), 
!     &                   LG = LLPAR, 1, -1 )
!         print*, ' SCIA_PRES, GC_CH4_onSCIA  ' 
!         WRITE(6,100) (SCIA(NT)%PRES(LS), GC_CH4_onSCIA(LS), 
!     &                       LS = LLSCIA, 1, -1 ) 
! 100  FORMAT(1X,F16.8,1X,F16.8)


         !--------------------------------------------------------------
         ! Apply SCIA observation operator
         !
         !   x_hat = A ( x_m ) + ( 1 - A ) x_a 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by SCIA [molec/cm2]
         !    x_a   = SCIA apriori column               [molec/cm2]
         !    x_m   = GC modeled column                 [molec/cm2]
         !    A     = SCIA averaging kernel 
         !--------------------------------------------------------------

         ! A ( x_m )
         LHS = 0d0
         DO LS = 1, LLSCIA
           LHS = LHS + SCIA(NT)%AVGKERNEL(LS) * GC_CH4_onSCIA(LS)
         ENDDO
     
         ! ( 1 - A ) x_a
         RHS = 0d0
         DO LS = 1, LLSCIA
            RHS = RHS + ( ( 1 - SCIA(NT)%AVGKERNEL(LS) )
     &                              * CH4_PRIOR(LS)  )
         ENDDO

         ! x_hat = A ( x_m ) + ( 1 - A ) x_a
         GC_CH4 = RHS + LHS

         ! Convert Units from [molec/cm2] --> [v/v]
         ! Get molecules of air in column using F=ma
         molec_air_total = 0d0
         molec_air_total = SCIA(NT)%PRESEDGE(1) * 1d2 / 9.8 * 1d-4  ! air [kg/cm2]
     &                  * 1d3 * (1/28.9644) * 6.022d23     ! air [molec/cm2]

         ! [molec ch4 / cm2] --> [v/v]   = molec CH4 / molec air in one cm^2
         GC_XCH4 = GC_CH4 / molec_air_total
!         print*,'FDTEST : GC_XCH4 = ',GC_XCH4


         !--------------------------------------------------------------
         ! Calculate cost function, given S is observation error covariance matrix
         !     Sobs = 1x1 array [ (molec/cm2) ^2 ]
         ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------


         ! Calculate error on this day.
         ! Fractional error = ERR_FRAC, a module variable
         Sobs = ( ERR_FRAC * 1d-9 * SCIA(NT)%XCH4 ) **2


         ! Calculate difference between modeled and observed profile
         DIFF = GC_XCH4 - 1d-9 * SCIA(NT)%XCH4
!         print*,'FDTEST : DIFF = ',DIFF


         ! Calculate J(x) = DIFF^T * S_{obs}^{-1} * DIFF 
         COST_FUNC_A = DIFF**2 / Sobs

         ! Calculate dJ/dx = 2 * DIFF * S_{obs}^{-1}
         FORCE = 0d0
         FORCE = 2 * DIFF / Sobs
         !print*,'FDTEST : FORCE = ',FORCE



         ! dkh debug: compare profiles:
!         print*, ' SCIA_PRIOR, XCH4_SCIA, XCH4_GC'
!         WRITE(6,101) (SCIA(NT)%PRIOR(L), SCIA(NT)%PRIOR(L), GC_XCH4,
!     &             L,    L = LLSCIA, 1, -1 )
! 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1x,i3)

         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

         ! dkh debug
!         print*, 'DIFF , FORCE, Sobs ' 
!         WRITE(6,102) (DIFF, FORCE, Sobs)
! 102  FORMAT(1X,d14.6,1X,d14.6)


         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE


         ! Adjoint of GEOS-Chem - SCIAMACHY difference
         GC_XCH4_ADJ = DIFF_ADJ

         ! Adjoint of unit conversion from [molec/cm2] --> [v/v]
         GC_CH4_ADJ = GC_XCH4_ADJ / molec_air_total
         !print*,'FDTEST : GC_CH4_ADJ = ',GC_CH4_ADJ


         ! Adjoint of SCIA observation operator
         DO LS=1,LLSCIA
            GC_CH4_onSCIA_ADJ(LS) = SCIA(NT)%AVGKERNEL(LS) * 
     &                                   GC_CH4_ADJ
         ENDDO
         !print*,'FDTEST : GC_CH4_ONSCIA_ADJ = ',GC_CH4_ONSCIA_ADJ

         ! Adjoint of unit conversion [v/v] --> [molec/cm2]
         DO LS=1,LLSCIA
            GC_CH4_onSCIA_ADJ(LS) = GC_CH4_onSCIA_ADJ(LS)
     &                                * molec_air_onSCIA(LS)
         ENDDO
         !print*,'FDTEST : GC_CH4_ONSCIA_ADJ = ',GC_CH4_ONSCIA_ADJ

         ! Adjoint of replacing GEOS-Chem stratosphere with SCIA prior
         DO LS=1,LLSCIA
            ! If trop pressure within grid box, replace fraction of value
            IF ( ( GC_TROP .lt. SCIA(NT)%PRESEDGE(LS) )   .AND.  
     &           ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS+1) ) ) THEN
               fracreplace = ( GC_TROP - SCIA(NT)%PRESEDGE(LS+1) ) / 
     &               ( SCIA(NT)%PRESEDGE(LS) - SCIA(NT)%PRESEDGE(LS+1) )
               GC_CH4_onSCIA_ADJ(LS) = 
     &                       (1-fracreplace) * GC_CH4_onSCIA_ADJ(LS)
            ENDIF

            ! If trop pressure gt lower grid box boundary, GC_CH4_onSCIA(LS)=0
            IF ( GC_TROP .gt. SCIA(NT)%PRESEDGE(LS) ) THEN
               GC_CH4_onSCIA_ADJ(LS) = 0d0
            ENDIF

         ENDDO



         ! Adjoint of interpolation
         DO LG=1,LLPAR
         DO LS=1,LLSCIA
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) + 
     &           GRIDMAP(LG,LS) * GC_CH4_onSCIA_ADJ(LS)
         ENDDO
         ENDDO
         !print*,'FDTEST : GC_CH4_NATIVE_ADJ = ',GC_CH4_NATIVE_ADJ



         ! Adjoint of GEOS-Chem column averaging
         !   Distribute adjoint forcing across GEOS-Chem grid boxes from 
         !   which the original GEOS-Chem column was calculated.
         ! Adjoint forcing added to adjoint tracer array in this subroutine

            ! Determine number of GEOS-Chem boxes covered by the observation
            nboxes = count( SCIA(NT)%GCfrac(:) .gt. 0.0 )
      
            ! Loop over boxes, placing adjoint variable into the STT_ADJ array
            DO NB=1,nboxes
      
               ! Clear variables to be safe
               I    = 0
               J    = 0
               frac = 0d0
      
               ! I and J indices and fractional influence of this box
               I    = SCIA(NT)%GCII(NB)
               J    = SCIA(NT)%GCJJ(NB)
               frac = SCIA(NT)%GCfrac(NB)
      
               ! Adjoint of unit conversion from [kg/box] to [v/v]
               DO LG=1,LLPAR
      
                  ! Get mass in this grid box
                  thisad1 = 0d0
                  thisad1 = AD(I,J,LG)
      
                  ! adjoint of unit conversion
                  thisforce(LG) = GC_CH4_NATIVE_ADJ(LG) * XNUMOL(1) / 
     &                                 ( thisad1 * XNUMOLAIR ) * frac  
     

                  ! Calculate adjoint sensitivity for output
                  ADJ(LG,NB) = thisforce(LG) * CHK_STT(I,J,LG,1)

               ENDDO
      
            ENDDO

         ! End distributing adjoint forcing to STT_ADJ array
            !print*,'-------------- FDTEST thisforce ------------------'
            !print*,thisforce



      ! Return to calling program
      END SUBROUTINE CALC_SCIA_CH4_FORCE_FD

!------------------------------------------------------------------------------



      SUBROUTINE GET_NT_RANGE( NSCIA, GCNHMS, TIME_FRAC, 
     &                         NTSTART, NTSTOP )
!
!******************************************************************************
!  Subroutine GET_NT_RANGE 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) GCNYMD (INTEGER) : Current model YYYYMMDD 
!  (2 ) GCNHMS (INTEGER) : Current model HHMMSS  
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) NTSTART (INTEGER) : SCIA record number at which to start
!  (1 ) NTSTOP  (INTEGER) : SCIA record number at which to stop
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TIME_MOD,     ONLY : YMD_EXTRACT

      ! Arguments
      INTEGER, INTENT(IN)   :: NSCIA
      INTEGER, INTENT(IN)   :: GCNHMS
      TYPE (XPLEX),  INTENT(IN)   :: TIME_FRAC(NSCIA)
      INTEGER, INTENT(OUT)  :: NTSTART
      INTEGER, INTENT(OUT)  :: NTSTOP
    
      ! Local variables 
      INTEGER, SAVE         :: NTSAVE
      LOGICAL               :: FOUND_ALL_RECORDS 
      INTEGER               :: NTEST
      INTEGER               :: HH, MM, SS
      TYPE (XPLEX)                :: GC_HH_FRAC
      TYPE (XPLEX)                :: H1_FRAC

      !=================================================================
      ! GET_NT_RANGE begins here!
      !=================================================================


      ! Initialize 
      FOUND_ALL_RECORDS  = .FALSE. 
      NTSTART            = 0
      NTSTOP             = 0

      ! set NTSAVE to NTES every time we start with a new file
      IF ( GCNHMS == 230000 ) NTSAVE = NSCIA


      print*, ' GET_NT_RANGE ', GCNHMS
      print*, ' NTSAVE ', NTSAVE
      print*, ' NSCIA', NSCIA

      CALL YMD_EXTRACT( GCNHMS, HH, MM, SS )


      ! Convert HH from hour to fraction of day 
      GC_HH_FRAC = DCMPLX(HH) / 24d0 
 
      ! one hour as a fraction of day 
      H1_FRAC    = 1d0 / 24d0 

    
      ! All records have been read already 
      IF ( NTSAVE == 0 ) THEN 

         print*, 'All records have been read already '
         RETURN 

      ! No records reached yet
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC < GC_HH_FRAC ) THEN 
           
      
         print*, 'No records reached yet'
         RETURN

      !
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC >=  GC_HH_FRAC ) THEN 
      
         ! Starting record found
         NTSTART = NTSAVE   

         print*, ' Starting : TIME_FRAC(NTSTART) ', 
     &               TIME_FRAC(NTSTART), NTSTART
 
         ! Now search forward to find stopping record
         NTEST = NTSTART

         DO WHILE ( FOUND_ALL_RECORDS == .FALSE. ) 
              
            ! Advance to the next record
            NTEST = NTEST - 1  
           
            ! Stop if we reach the earliest available record 
            IF ( NTEST == 0 ) THEN 
           
               NTSTOP            = NTEST + 1
               FOUND_ALL_RECORDS = .TRUE.

               print*, ' Records found '
               print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP

               ! Reset NTSAVE 
               NTSAVE = NTEST

            ! When the combined test date rounded up to the nearest
            ! half hour is smaller than the current model date, the 
            ! stopping record has been passed.
               !kjw
               ! shouldn't the line below be:
               ! ELSEIF (  TIME_FRAC(NTEST) + H1_FRAC/2d0 <  GC_HH_FRAC ) THEN
               ! (difference is dividing H1_FRAC by 2)
               ! necessary to round to nearest half hour instead of full hour
               !kjw
            ELSEIF (  TIME_FRAC(NTEST) + H1_FRAC <  GC_HH_FRAC ) THEN
          
               print*, ' Testing : TIME_FRAC ', 
     &                  TIME_FRAC(NTEST), NTEST
 
               NTSTOP            = NTEST + 1 
               FOUND_ALL_RECORDS = .TRUE. 

               print*, ' Records found '
               print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP
                  
               ! Reset NTSAVE 
               NTSAVE = NTEST 

            ELSE 
               !print*, ' still looking ', NTEST 
                  
            ENDIF 
                 
         ENDDO 
 
      ELSE

         CALL ERROR_STOP('problem', 'GET_NT_RANGE' ) 

      ENDIF 


      ! Return to calling program
      END SUBROUTINE GET_NT_RANGE

!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( GC_PEDGE, SCIA_PEDGE  )
     &         RESULT      ( M )
!
!******************************************************************************
!  Function GET_INTMAP creates the matrix that places GEOS-Chem column methane 
!     [molec/cm2] onto the 12-level pressure grid used by SCIAMACHY, M.
!         GC[1x47] * M[47x12] = SCIA[1x12]           (kjw, 7/21/11)
!
!  Arguments as Input:
!  ============================================================================
!  (3 ) GC_PEDGE   (TYPE (XPLEX)) : LLPAR bottom pressure edges of GEOS-Chem column
!  (4 ) SCIA_PEDGE (TYPE (XPLEX)) : LLSCIA+1 pressure edges of SCIA column
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) M          (TYPE (XPLEX)) : Interpolation matrix that maps GEOS-Chem to SCIA grid
!     
!  NOTES:
!  (1 ) Based on GET_INTMAP by Daven Henze (I think). See tes_nh3_mod.f for example
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX)             :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)             :: SCIA_PEDGE(LLSCIA+1)
 
      ! Return value 
      TYPE (XPLEX)             :: M(LLPAR,LLSCIA)

      ! Local variables 
      INTEGER  :: LGC, LTM, LS, LG
      TYPE (XPLEX)   :: DIFF, DELTA_SURFP
      TYPE (XPLEX)   :: GUP, GLO, SUP, SLO
      TYPE (XPLEX)   :: column_total(LLSCIA)

      !=================================================================
      ! GET_INTMAP begins here!
      !=================================================================

      ! Initialize output
      M(:,:) = 0D0 

      ! Loop over each pressure level of GEOS-Chem and SCIAMACHY grids
      DO LG=1,LLPAR-1

         ! Get upper and lower pressure edges of GEOS-Chem box
         GUP = GC_PEDGE( LG+1 )
         GLO = GC_PEDGE( LG   )
         
         DO LS=1,LLSCIA

            ! Get top and bottom pressures of SCIA box
            SUP = SCIA_PEDGE( LS+1 )
            SLO = SCIA_PEDGE( LS   )

            ! If both GEOS-Chem edges are within the SCIA box, map value = 1
            IF ( ( GUP .gt. SUP ) .AND. ( GLO .lt. SLO ) ) THEN
               M(LG,LS) = 1
            ENDIF

            ! If both GEOS-Chem stradles a SCIA pressure level, interpolate
            IF ( ( GUP .lt. SUP ) .AND. ( GLO .gt. SUP ) ) THEN
               DIFF       = GLO - GUP
               M(LG,LS+1) = ( SUP - GUP ) / DIFF
               M(LG,LS  ) = ( GLO - SUP ) / DIFF
            ENDIF

         ENDDO
      ENDDO

      ! Add value for uppermost GEOS-Chem grid box
      M(LLPAR,LLSCIA) = 1


      ! Correct for case in which GEOS-Chem pressure is higher than SCIAMACHY
      IF ( GC_PEDGE(1) .GT. SCIA_PEDGE(1) ) THEN


         ! If any part of GEOS-Chem box are under SCIA_PEDGE(1), let 
         !   this GEOS-Chem grid box contribute to the observation because 
         !   SCIA and GEOS-Chem should have same surface pressure. map value = 1
         DO LG=1,LLPAR-1

            ! If GEOS-Chem box entirely below SCIA surface pressure
            IF ( ( GC_PEDGE(LG)   .GT. SCIA_PEDGE(1) ) .AND.    
     &           ( GC_PEDGE(LG+1) .GT. SCIA_PEDGE(1) ) ) THEN 
               M(LG,1) = 1
            ENDIF

            ! If GEOS-Chem box straddles SCIA surface pressure
            IF ( ( GC_PEDGE(LG)   .GT. SCIA_PEDGE(1) ) .AND.    
     &           ( GC_PEDGE(LG+1) .LT. SCIA_PEDGE(1) ) ) THEN 
               DIFF = GC_PEDGE(LG) - GC_PEDGE( LG+1 )
               M(LG,1) = ( SCIA_PEDGE(1) - GC_PEDGE(LG+1) ) / DIFF
            ENDIF
            
         ENDDO
      ENDIF


      ! Correct for case in which GEOS-Chem surface pressure is within 2nd SCIA
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. SCIA_PEDGE(2) ) THEN
         M(1,1) = 1
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 3rd SCIA
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. SCIA_PEDGE(3) ) THEN
         M(1,2) = 1
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 4th SCIA
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. SCIA_PEDGE(4) ) THEN
         M(1,3) = 1
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 5th SCIA
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. SCIA_PEDGE(5) ) THEN
         M(1,4) = 1
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 6th SCIA
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. SCIA_PEDGE(6) ) THEN
         M(1,5) = 1
      ENDIF

      ! Normalize each column of M to 1 so that we are not creating any molecules
      ! when mapping from GEOS-Chem to SCIA grids.

      ! DO NOT do this since we are mapping molc/cm2, not 
      ! Initialize to be safe and calculate column total
      column_total(:) = 0d0
      column_total(:) = SUM( M, DIM=1 )

      ! Normalize columns to column_total
      DO LS=1,LLSCIA
         M(:,LS) = M(:,LS) / column_total(LS)
      ENDDO



      ! Return to calling program
      END FUNCTION GET_INTMAP

!-----------------------------------------------------------------------------



      END MODULE SCIA_CH4_MOD
