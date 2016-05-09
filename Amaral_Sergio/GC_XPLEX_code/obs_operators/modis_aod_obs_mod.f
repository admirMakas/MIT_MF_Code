!$Id: modis_aod_obs_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE MODIS_AOD_OBS_MOD
! 
!******************************************************************************
! Mdoule MODIS_AOD_OBS_MOD contains subroutines necessary to
! 1. READ modis aot observations (Wang et al., 2010), and apriori
! 2. COMPUTE MOIDS-CEOS-Chem difference, cost function, and adj forcing
!
!   (xxu, 7/20/10, 5/17/11)
! Added to standard code  (xxu, dkh, 01/12/12, adj32_011) 
!
!  Module Variables:
!  ============================================================================
! ( 1) NSPECI      (INTEGER)  : # of aerosol species included for adjoint 
! ( 2) NLEV        (INTEGER)  : # of obs layers
! ( 3) MAXMODIS    (INTEGER)  : Max # of obs per day, used for array defining
! ( 4) MODIS       (TYPE   )  : Record data from each MODIS obs
! (  ) IDT_MODIS   (INTEGER)  : Available modis obs tracers' id
!  
!  Module Routines:
!  ============================================================================
! (1) READ_MODIS_AOD_OBS     : Read modis obs from netCDF file
! (2) CALC_MODIS_AOD_FORCE   : Calculates cost function, obs forcing
! (3) CHECK                  : Check status for calling netCDF
! (4) GET_NT_RANGE           : Return the obs range for current hour
! (5) PCENTL()               : Function calculating percentiles
! (6) HPSORT                 : Sort an array by Heapsort method
! 
!  ============================================================================
!  NOTES:
! (1) This module is copied and adapted from Daven Henze's 'tes_o3_mod.f', 
!     which is the operator calculating the adjoint forcing from the TES O3
!     observations. (xxu, 7/20/10)
! (2) Initial code was designed for pixel-based MODIS observations. Here is 
!     modified for those observation aggregated to each grid-box, to aviod 
!     wired single observation values. (xxu,9/9/10)
! (3) Now only consider the troposhere by using ITS_IN_THE_TROP from
!     tropopause_mod.f. (xxu, 6/14/11)
!******************************************************************************
!  
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "aerosol_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CALC_MODIS_AOD_FORCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
 
      ! Parameters
      INTEGER, PARAMETER  :: MAX_MODIS_TRC = 11
      INTEGER, PARAMETER  :: MAX_MODIS_OBS = 500
      INTEGER, PARAMETER  :: NSPECI   = 11
      INTEGER, PARAMETER  :: NLEV     = 47
      INTEGER, PARAMETER  :: MAXMODIS = 500
      TYPE (XPLEX),  PARAMETER  :: OBS_CUTOFF = 1d-3

      ! Variables
      INTEGER             :: IDT_MODIS(MAX_MODIS_TRC)

      ! Record to store data from each MODIS obs
      TYPE MODIS_AOD_OBS
         TYPE (XPLEX)             :: LAT(1)
         TYPE (XPLEX)             :: LON(1)
         TYPE (XPLEX)             :: TIME(1)
         TYPE (XPLEX)           :: SFactor(1)                    ! Ratio of obs to model
         TYPE (XPLEX)           :: AP_MODIS (NLEV,MAX_MODIS_TRC) ! a priori
         TYPE (XPLEX)           :: OBS_MODIS(NLEV,MAX_MODIS_TRC) ! Observations
         TYPE (XPLEX)           :: OER_INV  (NLEV,MAX_MODIS_TRC) ! diagnal elements 
      ENDTYPE MODIS_AOD_OBS

      TYPE(MODIS_AOD_OBS) :: MODIS(MAX_MODIS_OBS)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_MODIS_AOD_OBS
!
!******************************************************************************
!  Subourtine INIT_MODIS_AOD_OBS initialize the MODIS aerosol operator:
!   (1) Assign available modis obs tracers to array IDT_MODIS
!   (2) Check if the used obs are exactly same to those specified by the 
!       adjoint input files
!  (xxu, 6/28/11)
!******************************************************************************
!
      ! Reference to f90 modules
      USE TRACER_MOD,         ONLY : TRACER_NAME
      USE TRACERID_MOD,       ONLY : IDTSO4,  IDTNH4,  IDTNIT
      USE TRACERID_MOD,       ONLY : IDTBCPI, IDTOCPI, IDTBCPO, IDTOCPO
      USE TRACERID_MOD,       ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_THIS_TRACER, NOBS
      USE ERROR_MOD,          ONLY : ERROR_STOP

      ! Local variables
      INTEGER :: T, TT, COUNTER_TRC

      !================================================================
      ! INIT_MODIS_AOD_OBS begins here!
      !================================================================

      ! Assign tracers id to IDT_MODIS
      IDT_MODIS( 1) = IDTSO4
      IDT_MODIS( 2) = IDTNH4
      IDT_MODIS( 3) = IDTNIT
      IDT_MODIS( 4) = IDTBCPI
      IDT_MODIS( 5) = IDTOCPI
      IDT_MODIS( 6) = IDTBCPO
      IDT_MODIS( 7) = IDTOCPO
      IDT_MODIS( 8) = IDTDST1
      IDT_MODIS( 9) = IDTDST2
      IDT_MODIS(10) = IDTDST3
      IDT_MODIS(11) = IDTDST4

      ! Initialize the obs counter
      COUNTER_TRC = 0

      ! Start modis tracer loop
      DO T = 1, MAX_MODIS_TRC
        
         ! Global tracer ID
         TT = IDT_MODIS(T)

         ! Selected tracers to be obs?
         IF ( OBS_THIS_TRACER ( TT ) ) THEN 

            WRITE( 6, 100 ) TT, TRACER_NAME( TT )
            COUNTER_TRC = COUNTER_TRC + 1
     
         ENDIF

      ! Finish modis tracer loop: T 
      ENDDO

      ! Check if the counted obs equals to that specified 
      WRITE( 6, 110 ) COUNTER_TRC
      IF ( COUNTER_TRC /= NOBS )
     &   CALL ERROR_STOP( 'Error: selected modis obs tracer =/ NOBS', 
     &                    'init_modis_aer_obs (modis_aer_obs_mod.f)' )


 100  FORMAT( 3x, 'Used MODIS obs tracers: ', I4, 3x, A6)
 110  FORMAT( 3x, '# of selected MODIS obs: ', I4   )

      ! Return to the calling routine
      END SUBROUTINE INIT_MODIS_AOD_OBS
 
!------------------------------------------------------------------------------
 
      SUBROUTINE READ_MODIS_AOD_OBS( YYYYMMDD, NMODIS )
!
!******************************************************************************
!  Subroutine READ_MODIS_AOD_OBS reads the file and passes back info contained
!  therein. (xxu, 02/19/09) 
! 
!  Based on READ_TES_O3_OBS (dkh, 04/26/10) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD    (INTEGER) : Current year-month-day
!
!  Arguments as Output: 
!  ============================================================================
!  (1 ) NMODIS      (INTEGER) : Number of MODIS retrievals for current day 
!
!  Module variable as Output: 
!  ============================================================================
!  (1 ) MODIS (MODIS_AOD_OBS) : TES retrieval for current day 
!     
!  NOTES:
!  
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,    ONLY : DATA_DIR
      USE NETCDF 
      USE TIME_MOD,         ONLY : EXPAND_DATE

      ! Arguments
      INTEGER, INTENT( IN)   :: YYYYMMDD
      INTEGER, INTENT(OUT)   :: NMODIS    

      ! Local variables 
      INTEGER                :: FID,     NM_ID
      INTEGER                :: TIME_ID, LON_ID,  LAT_ID
      INTEGER                :: AOD1_ID, AOD2_ID, AOD3_ID, EOR_ID
      INTEGER                :: SO4_ID,  NH4_ID,  NIT_ID
      INTEGER                :: BCPI_ID, OCPI_ID, BCPO_ID, OCPO_ID
      INTEGER                :: DST1_ID, DST2_ID, DST3_ID, DST4_ID
!      INTEGER                :: START0(1), COUNT0(1)
!      INTEGER                :: START1(2), COUNT1(2)
!      INTEGER                :: START2(3), COUNT2(3)
      INTEGER                :: N, L, T
      CHARACTER(LEN=5)       :: TMP
      CHARACTER(LEN=255)     :: READ_FILENAME

      TYPE (XPLEX)                   :: EOR_FRAC, EFR_IN
      TYPE (XPLEX), ALLOCATABLE      :: TMP1(:)
      TYPE (XPLEX), ALLOCATABLE      :: TMP2(:,:)
      TYPE (XPLEX)                 :: TMP_OER_INV

      !=================================================================
      ! READ_MODIS_AOD_OBS begins here!
      !=================================================================

      ! filename root 
      READ_FILENAME = TRIM( 'modis_aod_obs_2x25_YYYYMMDD.nc' )

      ! Expand date tokens in filename 
      CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 ) 

      ! Construct complete filename 
      READ_FILENAME = TRIM( DATA_DIR ) // 'MODIS_AOD_OBS_201009/' // 
     &                TRIM( READ_FILENAME )

      ! Print to screen
      WRITE(6,100) TRIM(READ_FILENAME)
 100  FORMAT(' - READ_MODIS_AOD_OBS: reading file: ', A )

      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )
   
      !-------------------------------- 
      ! Get data record IDs
      !-------------------------------- 
      CALL CHECK( NF90_INQ_DIMID( FID, "hhmm",           NM_ID), 101 )
      CALL CHECK( NF90_INQ_VARID( FID, "obserrorfrac",  EOR_ID), 102 )
      CALL CHECK( NF90_INQ_VARID( FID, "dayfrc",       TIME_ID), 103 )
      CALL CHECK( NF90_INQ_VARID( FID, "lon",           LON_ID), 104 )
      CALL CHECK( NF90_INQ_VARID( FID, "lat",           LAT_ID), 105 )
      CALL CHECK( NF90_INQ_VARID( FID, "gc_aod",       AOD1_ID), 106 )
      CALL CHECK( NF90_INQ_VARID( FID, "ret_aod",      AOD2_ID), 107 )
      CALL CHECK( NF90_INQ_VARID( FID, "modis_aod",    AOD3_ID), 108 )
      CALL CHECK( NF90_INQ_VARID( FID, "SO4_GC",        SO4_ID), 109 )
      CALL CHECK( NF90_INQ_VARID( FID, "NH4_GC",        NH4_ID), 110 )
      CALL CHECK( NF90_INQ_VARID( FID, "NIT_GC",        NIT_ID), 111 )
      CALL CHECK( NF90_INQ_VARID( FID, "BCPI_GC",      BCPI_ID), 112 )
      CALL CHECK( NF90_INQ_VARID( FID, "OCPI_GC",      OCPI_ID), 113 )
      CALL CHECK( NF90_INQ_VARID( FID, "BCPO_GC",      BCPO_ID), 114 )
      CALL CHECK( NF90_INQ_VARID( FID, "OCPO_GC",      OCPO_ID), 115 )
      CALL CHECK( NF90_INQ_VARID( FID, "DST1_GC",      DST1_ID), 116 )
      CALL CHECK( NF90_INQ_VARID( FID, "DST2_GC",      DST2_ID), 117 )
      CALL CHECK( NF90_INQ_VARID( FID, "DST3_GC",      DST3_ID), 118 )
      CALL CHECK( NF90_INQ_VARID( FID, "DST4_GC",      DST4_ID), 119 )

      !-------------------------------- 
      ! Read dimensions
      !--------------------------------

      ! READ number of retrievals, NMODIS
      CALL CHECK( NF90_INQUIRE_DIMENSION(FID, NM_ID, TMP, NMODIS), 201)

      ! Print to screen
      WRITE(6,110) NMODIS
 110  FORMAT(' NMODIS = ', I6)

      !-------------------------------- 
      ! Read 0D Data
      !-------------------------------- 

      ! READ observation error fraction
      CALL CHECK( NF90_GET_VAR  ( FID, EOR_ID, EOR_FRAC ), 300)

      ! Print to screen
      WRITE(6,120)  EOR_FRAC
 120  FORMAT(' OBS ERROR FRACTION = ', F10.4)

      !-------------------------------- 
      ! Read 1D Data
      !-------------------------------- 

      ! allocate temporal arrays for 1D data
      ALLOCATE ( TMP1 (NMODIS)      )
      TMP1 = 0

      ! READ latitude 
      CALL CHECK( NF90_GET_VAR ( FID, LAT_ID, TMP1 ), 301)
      MODIS(1:NMODIS)%LAT(1) = TMP1(1:NMODIS)

      ! READ longitude
      CALL CHECK( NF90_GET_VAR ( FID, LON_ID, TMP1 ), 302)
      MODIS(1:NMODIS)%LON(1) = TMP1(1:NMODIS)

      ! READ time
      CALL CHECK( NF90_GET_VAR ( FID, TIME_ID, TMP1 ), 303)
      MODIS(1:NMODIS)%TIME(1) = TMP1(1:NMODIS)

      ! READ AODs and calculate the ratio of obs to model
      CALL CHECK( NF90_GET_VAR ( FID, AOD1_ID, TMP1 ), 304)
      MODIS(1:NMODIS)%SFACTOR(1) = TMP1(1:NMODIS)
      CALL CHECK( NF90_GET_VAR ( FID, AOD2_ID, TMP1 ), 305)
      MODIS(1:NMODIS)%SFACTOR(1) = TMP1(1:NMODIS) 
     &                           / MODIS(1:NMODIS)%SFACTOR(1)

      ! To avoid wired ratios
      DO N = 1, NMODIS
         IF ( MODIS(N)%SFACTOR(1) < 0.5 ) MODIS(N)%SFACTOR(1) = 0.5
         IF ( MODIS(N)%SFACTOR(1) > 10.0 ) MODIS(N)%SFACTOR(1) = 10.0
      ENDDO

      !-------------------------------- 
      ! Read 2D Data
      !--------------------------------

      ! allocate temporal arrays for 2D data 
      ALLOCATE ( TMP2 (NMODIS,NLEV) )
      TMP2 = 0

      ! READ SO4 A Priori
      CALL CHECK( NF90_GET_VAR ( FID, SO4_ID, TMP2 ), 401)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,1) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ NH4 A Priori
      CALL CHECK( NF90_GET_VAR ( FID, NH4_ID, TMP2 ), 402)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,2) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ NIT A Priori
      CALL CHECK( NF90_GET_VAR ( FID, NIT_ID, TMP2 ), 403)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,3) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ BCPI A Priori
      CALL CHECK( NF90_GET_VAR ( FID, BCPI_ID, TMP2 ), 404)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,4) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ OCPI A Priori
      CALL CHECK( NF90_GET_VAR ( FID, OCPI_ID, TMP2 ), 405)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,5) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ BCPO A Priori
      CALL CHECK( NF90_GET_VAR ( FID, BCPO_ID, TMP2 ), 406)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,6) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ OCPO A Priori
      CALL CHECK( NF90_GET_VAR ( FID, OCPO_ID, TMP2 ), 407)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,7) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ DST1 A Priori
      CALL CHECK( NF90_GET_VAR ( FID, DST1_ID, TMP2 ), 408)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,8) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ DST2 A Priori
      CALL CHECK( NF90_GET_VAR ( FID, DST2_ID, TMP2 ), 409)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,9) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ DST3 A Priori
      CALL CHECK( NF90_GET_VAR ( FID, DST3_ID, TMP2 ), 410)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,10) = TMP2(1:NMODIS,L)
      ENDDO

      ! READ DST4 A Priori
      CALL CHECK( NF90_GET_VAR ( FID, DST4_ID, TMP2 ), 411)
      DO L = 1, NLEV
         MODIS(1:NMODIS)%AP_MODIS(L,11) = TMP2(1:NMODIS,L)
      ENDDO

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )

      ! deallocate arrays
      IF ( ALLOCATED(TMP1) )  DEALLOCATE( TMP1 )
      IF ( ALLOCATED(TMP2) )  DEALLOCATE( TMP2 )

      !================================================================
      ! Calculate MODIS obs: obs = apriori * sfactor
      !================================================================
      DO L = 1, NLEV
      DO T = 1, MAX_MODIS_TRC
         MODIS(1:NMODIS)%OBS_MODIS(L,T) = MODIS(1:NMODIS)%AP_MODIS(L,T)
     &                                  * MODIS(1:NMODIS)%SFACTOR(1)
      ENDDO
      ENDDO

      !================================================================
      ! obs error covriance maxtrices and their inverse
      !================================================================

      !----------------------------------------------------------------
      ! NOTE: 
      !  (1) from inverse error fraction to absolute error
      !        OER_INV = err^(-2) = ( obs * err_frac )^(-2)
      !  (2) put a cap on the error
      !        if (obs <= 0.001) OER_INV = 1d0
      !        if (obs >  0.001) OER_INV = MIN(25d0, OER_INV)
      !---------------------------------------------------------------- 

      ! Inverse of error fraction square 
      EFR_IN = 1.d0 / EOR_FRAC / EOR_FRAC

      ! Calculate the obs error inverse matrix
      DO N = 1, NMODIS
      DO L = 1, NLEV
      DO T = 1, MAX_MODIS_TRC

         ! Only consider 
         IF ( MODIS(N)%OBS_MODIS(L,T) >= OBS_CUTOFF ) THEN
            TMP_OER_INV = EFR_IN / MODIS(N)%OBS_MODIS(L,T)
            MODIS(N)%OER_INV(L,T) = MIN( 25d0, TMP_OER_INV )
         ELSE
            MODIS(N)%OER_INV(L,T) = 1d0
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_MODIS_AOD_OBS
!
!------------------------------------------------------------------------------
!
      SUBROUTINE CALC_MODIS_AOD_FORCE( COST_FUNC )

!******************************************************************************
! Subroutine CALC_MODIS_AOD_FORCE calculates the adjoint frocing from the MODIS
! retrieval (Wang et al., 2010) and the cost function. (xxu, 7/20/10)
!
!  Arguments as Input/Output:
!  ============================================================================
! (1 ) COST_FUNC (TYPE (XPLEX))  : Cost funciton                        [unitless]
!  
!  NOTES:
! (1 ) This subroutine is adopted from Daven's "CALC_TES_O3_FORCE", which is 
!      for TES O3. (xxu, 7/20/10)
!
!******************************************************************************

      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC, STT_ADJ, EXPAND_NAME
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_THIS_TRACER
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE DAO_MOD,            ONLY : AD, AIRDEN
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS, GET_TS_CHEM
      USE TRACER_MOD,         ONLY : TCVV, TRACER_NAME
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
 
      ! Local variables
      INTEGER,SAVE                :: N_MODIS_OBS
      INTEGER,SAVE                :: MODIS_TIME(MAX_MODIS_OBS)
      LOGICAL,SAVE                :: FIRST = .TRUE.
      TYPE (XPLEX), SAVE                :: MODIS_DAYFRC(MAX_MODIS_OBS)

      INTEGER                     :: NTSTART, NTSTOP, NT
      INTEGER                     :: IIJJ(2), I, J, L, N
      INTEGER                     :: T, TT, NSP
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME

      TYPE (XPLEX)                      :: NEW_COST(MAX_MODIS_OBS)
      TYPE (XPLEX)                      :: OLD_COST
      TYPE (XPLEX)                      :: H_GC(NLEV, MAX_MODIS_TRC)
      TYPE (XPLEX)                      :: DIFF(NLEV, MAX_MODIS_TRC)
      TYPE (XPLEX)                      :: FORCE(NLEV, MAX_MODIS_TRC)
      TYPE (XPLEX)                      :: DIFF_ADJ(NLEV, MAX_MODIS_TRC)    
     
      !=================================================================
      ! CALC_MODIS_AOD_FORCE begins here!
      !=================================================================

      print*, '     - CALC_MODIS_AOD_FORCE: MODIS AOD forcing '

      ! Initialize
      CALL INIT_MODIS_AOD_OBS

      ! Reset 
      NEW_COST = 0D0

      ! Open files for diagnostic output
      IF ( FIRST ) THEN

         ! Start modis tracer loop
         DO T = 1, MAX_MODIS_TRC
         
            TT = IDT_MODIS(T)
            IF ( OBS_THIS_TRACER( TT ) ) THEN
               FILENAME = 'debug_'//TRIM(TRACER_NAME(TT))//'_ITRNN.m'
               CALL EXPAND_NAME( FILENAME, N_CALC )
               FILENAME = TRIM(DIAGADJ_DIR) // TRIM(FILENAME)
               OPEN( 100+T, FILE=TRIM(FILENAME), STATUS='UNKNOWN',
     &               IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )
            ENDIF
         
         ! Finish modis tracer loop: T      
         ENDDO

      ENDIF  ! FIRST

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Check if it is the last hour of a day 
      IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN

         ! Read the MODIS obs file for this day 
         CALL READ_MODIS_AOD_OBS( GET_NYMD(), N_MODIS_OBS )

         ! Day fraction.
         DO N = 1, N_MODIS_OBS
            MODIS_DAYFRC(1:N) = MODIS(1:N)%TIME(1)
         ENDDO

      ENDIF

      ! Get the range of MODIS retrievals for the current hour
      CALL GET_NT_RANGE( N_MODIS_OBS, GET_NHMS(), 
     &                   MODIS_DAYFRC, NTSTART, NTSTOP )

      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN
         PRINT*, ' No matching MODIS obs for this hour: ', GET_NHMS()
         RETURN
      ENDIF

      PRINT*, ' For hour range: ', GET_NHMS(), MODIS_DAYFRC(NTSTART),
     &       MODIS_DAYFRC(NTSTOP)
      PRINT*, ' found record range: ', NTSTART, NTSTOP

      ! loop for this GC hour
      DO NT  = NTSTART, NTSTOP, -1

         ! Get grid box of current record
         IIJJ  = GET_IJ( TYPE (XPLEX)(MODIS(NT)%LON(1),4), 
     &                   TYPE (XPLEX)(MODIS(NT)%LAT(1),4))
         I     = IIJJ(1)
         J     = IIJJ(2)

         H_GC(:,:) = 0d0
         DIFF(:,:) = 0d0

         ! modeled profile (convert units from kg/box to ppbv)
         DO L = 1, NLEV  
         DO T = 1, MAX_MODIS_TRC

            ! Tracer id
            TT = IDT_MODIS(T)

            ! Check if this is an obs tracer
            IF ( OBS_THIS_TRACER( TT ) ) THEN

               ! GC simulation on the Observation space 
               H_GC(L,T) = CHK_STT(I,J,L,TT) * TCVV(TT) 
     &                   * 1d9 / AD(I,J,L)

               ! Difference of observations from model simulation
               DIFF(L,T) = H_GC(L,T) - MODIS(NT)%OBS_MODIS(L,T)

               ! Adjoint forcing: S_{obs}^{-1} * DIFF
               FORCE(L,T) = MODIS(NT)%OER_INV(L,T) * DIFF(L,T)

               ! Contribution to the cost function: 
               !  1/2 * DIFF^T * S_{obs}^{-1} * DIFF
               NEW_COST(NT) = NEW_COST(NT) + .5d0*DIFF(L,T)*FORCE(L,T)

               ! Now pass the adjoint back to the adjoint tracer array
               DIFF_ADJ(L,T) = FORCE(L,T)
               STT_ADJ(I,J,L,TT) = STT_ADJ(I,J,L,TT) + DIFF_ADJ(L,T)
     &                           * TCVV(TT) * 1d9 / AD(I,J,L)

               ! Debug -xxu
               WRITE(100+T, 110) NT, I, J, L, 
     &                           H_GC(L,T), 
     &                           MODIS(NT)%AP_MODIS(L,T),
     &                           MODIS(NT)%OBS_MODIS(L,T), 
     &                           MODIS(NT)%SFACTOR(1),
     &                           MODIS(NT)%OER_INV(L,T),
     &                           FORCE(L,T), 
     &                           NEW_COST(NT), 
     &                           STT_ADJ(I,J,L,TT)

            ENDIF

         ENDDO
         ENDDO

110      FORMAT(4I6,1P8d10.2)

      ! finish NT loop
      ENDDO ! NT

      ! Update cost function
      COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))
      print*, ' Updated value of COST_FUNC = ', COST_FUNC
      print*, ' MODIS AOD contribution     = ', COST_FUNC - OLD_COST

      ! Reset 
      IF ( FIRST ) FIRST = .FALSE.

      ! debug      
      print*, ' MAX STT_ADJ  = ', MAXVAL(STT_ADJ(:,:,:,:))
      print*, ' MAX in       = ', MAXLOC(STT_ADJ(:,:,:,:))
      print*, ' MAX NEW_COST = ', MAXVAL(NEW_COST)
      print*, ' MAX cost in  = ', MAXLOC(NEW_COST)

      ! Return to calling program
      END SUBROUTINE CALC_MODIS_AOD_FORCE
!
!------------------------------------------------------------------------------
!
      SUBROUTINE CHECK( STATUS, LOCATION )
!
!******************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries routines
!  (dkh, 02/15/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call    
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was made   
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules 
      USE ERROR_MOD,    ONLY  : ERROR_STOP
      USE NETCDF 
 
      ! Arguments
      INTEGER, INTENT(IN)    :: STATUS 
      INTEGER, INTENT(IN)    :: LOCATION
    
      !=================================================================
      ! CHECK begins here!
      !=================================================================

      IF ( STATUS /= NF90_NOERR ) THEN 
        WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
        WRITE(6,*) 'At location = ', LOCATION 
        CALL ERROR_STOP('netCDF error', 'modis_aod_mod')
      ENDIF 

      ! Return to calling program
      END SUBROUTINE CHECK
!
!------------------------------------------------------------------------------
!
      SUBROUTINE GET_NT_RANGE( NTES, HHMMSS, TIME_FRAC, NTSTART, NTSTOP)
!     
!******************************************************************************
!  Subroutine GET_NT_RANGE retuns the range of retrieval records for the 
!  current model hour 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTES   (INTEGER) : Number of TES retrievals in this day 
!  (2 ) HHMMSS (INTEGER) : Current model time 
!  (3 ) TIME_FRAC (TYPE (XPLEX)) : Vector of times (frac-of-day) for the TES retrievals
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) NTSTART (INTEGER) : TES record number at which to start
!  (1 ) NTSTOP  (INTEGER) : TES record number at which to stop
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TIME_MOD,     ONLY : YMD_EXTRACT

      ! Arguments
      INTEGER, INTENT(IN)   :: NTES
      INTEGER, INTENT(IN)   :: HHMMSS
      TYPE (XPLEX),  INTENT(IN)   :: TIME_FRAC(NTES)
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
      IF ( HHMMSS == 230000 ) NTSAVE = NTES

      ! for debug only, need to change back to above line when oneline
      !IF ( HHMMSS == 230000 ) NTSAVE = NTES

      print*, ' GET_NT_RANGE for ', HHMMSS
      print*, ' NTSAVE ', NTSAVE
      print*, ' NTES   ', NTES

      CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )


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
               print*, ' still looking ', NTEST

            ENDIF

         ENDDO

      ELSE

         CALL ERROR_STOP('problem', 'GET_NT_RANGE' )

      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_NT_RANGE
!
!------------------------------------------------------------------------------
!
      END MODULE MODIS_AOD_OBS_MOD
