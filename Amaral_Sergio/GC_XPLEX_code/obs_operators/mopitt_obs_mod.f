      MODULE MOPITT_OBS_MOD
!
!*****************************************************************************
!  Module MOPITT_OBS_MOD contains all the subroutines for the using of MOPITT
!  observation (version 3 and version 4).
!  (zhe 1/19/11)
!  Module Routines:
!  ============================================================================
!  (1 ) READ_MOPITT_FILE       : Read MOPITT hdf file
!  (2 ) CALC_MOPITT_FORCE      : Calculates cost function and STT_ADJ increments
!  (3 ) CALC_AVGKER            : Construct the averging kernel matrix
!  (4 ) BIN_DATA_V4            : Interpolation between different vertical resolutions
!  (5 ) INIT_DOMAIN            : Define the observation window
!  (6 ) CALC_OBS_HOUR          : Calculated hour of morning obs
!  (7 ) ITS_TIME_FOR_MOPITT_OBS: FUNCTION that checks time vs. OBS_HOUR array
!  (8 ) READ_MOP02             : Reads MOPITT data fields from the HDF-EOS file
!  (9 ) APRIORI_MOP02          : Read A priori field for MOPITT version 3
!  (10) INFO_MOP02             : Prints name, dims, type, etc. of MOPITT data fields
!  (11) CLEANUP_MOP02          : Deallocates all module arrays
!  =============================================================================

       USE MYTYPE
       USE COMPLEXIFY
       IMPLICIT NONE

#     include "CMN_SIZE"
#     include "../adjoint/define_adj.h"

      PRIVATE

      PUBLIC OBS_HOUR_MOPITT
      PUBLIC COUNT_TOTAL
      PUBLIC ITS_TIME_FOR_MOPITT_OBS
      PUBLIC READ_MOPITT_FILE
      PUBLIC CALC_MOPITT_FORCE

      !=============================================================================
      ! MODULE VARIABLES
      !=============================================================================

      INTEGER   :: OBS_HOUR_MOPITT(IIPAR,JJPAR)
      INTEGER   :: DOMAIN_OBS(IIPAR,JJPAR)
      TYPE (XPLEX)    :: COUNT_TOTAL

      TYPE (XPLEX)    :: ERR_PERCENT(IIPAR,JJPAR)
      TYPE (XPLEX),  ALLOCATABLE :: A(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: T(:)
      TYPE (XPLEX),  ALLOCATABLE :: XA(:)
      TYPE (XPLEX),  ALLOCATABLE :: AC(:)

      ! MOPITT dimension fields
      INTEGER   :: T_DIM, Z_DIM
      TYPE (XPLEX), ALLOCATABLE :: LATITUDE(:)
      TYPE (XPLEX), ALLOCATABLE :: LONGITUDE(:)
      TYPE (XPLEX), ALLOCATABLE :: PRESSURE(:)
      TYPE (XPLEX), ALLOCATABLE :: SECONDS_IN_DAY(:)
      TYPE (XPLEX), ALLOCATABLE :: MOPITT_GMT(:)
      TYPE (XPLEX), ALLOCATABLE :: TAU(:)

      ! MOPITT data quantities
      TYPE (XPLEX), ALLOCATABLE :: BOTTOM_PRESSURE_V3(:,:)
      TYPE (XPLEX), ALLOCATABLE :: BOTTOM_PRESSURE_V4(:)
      TYPE (XPLEX), ALLOCATABLE :: CO_MIXING_RATIO(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CO_RET_BOT_MIXING_RATIO(:,:)
      TYPE (XPLEX), ALLOCATABLE :: CO_TOTAL_COLUMN(:,:)
      TYPE (XPLEX), ALLOCATABLE :: RET_ERR_COV(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: AVGKER_V4(:,:,:)
      INTEGER, ALLOCATABLE :: CLOUD_DES(:)
      INTEGER, ALLOCATABLE :: SURFACE_INDEX(:)

      ! MOPITT a priori
      INTEGER :: NLEV_AP
      TYPE (XPLEX), ALLOCATABLE :: PLEV_AP(:)
      TYPE (XPLEX), ALLOCATABLE :: CO_MR_AP(:)
      TYPE (XPLEX), ALLOCATABLE :: CH4_MR_AP(:)
      TYPE (XPLEX), ALLOCATABLE :: CO_MR_AP_V4(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CO_MR_AP_BOTTOM(:,:)
      TYPE (XPLEX), ALLOCATABLE :: COV_CO_AP(:,:)

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_MOPITT_FILE( YYYYMMDD, HHMMSS )

!******************************************************************************
!  Subroutine READ_MOPITT_FILE reads the MOPITT hdf file.
!  (mak, 7/12/07, zhe 1/19/11)
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!******************************************************************************

      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE TIME_MOD,  ONLY : EXPAND_DATE, GET_MONTH, GET_YEAR

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, AS
      CHARACTER(LEN=255)  :: DIR_MOPITT
      CHARACTER(LEN=255)  :: DIR_MONTH
      CHARACTER(LEN=255)  :: FILENAMEM
      CHARACTER(LEN=255)  :: FILENAME2
      LOGICAL, SAVE       :: FIRST = .TRUE.

      !=================================================================
      ! READ_MOPITT_FILE begins here!
      !=================================================================

#if defined( MOPITT_V3_CO_OBS )
      DIR_MOPITT = '/users/jk/05/zjiang/mopitt/'
      DIR_MONTH = 'V3/YYYY/MM/'
      FILENAMEM = 'MOP02-YYYYMMDD-L2V5.93.2.val.hdf'
#endif
#if defined( MOPITT_V4_CO_OBS )
      DIR_MOPITT = '/users/jk/05/zjiang/mopitt/'
      DIR_MONTH = 'V4/YYYY/MM/'
      FILENAMEM = 'MOP02-YYYYMMDD-L2V8.0.2.val.hdf'
#endif
#if defined( MOPITT_V5_CO_OBS )
      DIR_MOPITT = '/users/jk/08/zjiang/mopitt/'
      DIR_MONTH = 'V5/YYYY/MM/'
      FILENAMEM = 'MOP02J-YYYYMMDD-L2V10.1.3.beta.hdf'
#endif

      IF ( FIRST ) THEN
         ERR_PERCENT(:,:) = 0.0
         COUNT_TOTAL = 0
         FIRST            = .FALSE.
      ENDIF

      OBS_HOUR_MOPITT(:,:) = -99

      CALL EXPAND_DATE( FILENAMEM, YYYYMMDD, 0 )
      CALL EXPAND_DATE( DIR_MONTH, YYYYMMDD, 0 )

      FILENAME2 = TRIM( DIR_MOPITT ) // TRIM( DIR_MONTH ) // FILENAMEM
      PRINT*, '=== Reading ===:', TRIM( FILENAME2 )

      !CALL INFO_MOP02(FILENAME2)

      CALL READ_MOP02( FILENAME2 )

#if defined( MOPITT_V3_CO_OBS )
      CALL APRIORI_MOP02
#endif

      CALL INIT_DOMAIN

      ! Calculate hour of day when obs should be compared to model
      CALL CALC_OBS_HOUR

      !CALL READ_ERROR_VARIANCE
      !We assume 20% uniform observation error
      ERR_PERCENT(:,:) = 0.2

      END SUBROUTINE READ_MOPITT_FILE
!-------------------------------------------------------------------------------------------------

      SUBROUTINE CALC_MOPITT_FORCE

!******************************************************************************
! CALC_MOPITT_FORCE calculate cost function and STT_ADJ increments
! MOPITT version 3: total column approach
! MOPITT version 4, vertical profile approach
! (zhe 1/19/11)
!******************************************************************************

      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_AP, GET_BP
      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_YEAR
      USE TIME_MOD,     ONLY : GET_HOUR
      USE CHECKPT_MOD,  ONLY : CHK_STT
      USE TRACER_MOD,   ONLY : TCVV
      USE TRACERID_MOD,   ONLY : IDTCO
      USE DAO_MOD,      ONLY : AD
      USE ADJ_ARRAYS_MOD,   ONLY : SET_FORCING, SET_MOP_MOD_DIFF,
     &                  SET_MODEL_BIAS, SET_MODEL, SET_OBS,
     &                  COST_ARRAY, DAY_OF_SIM, IFD, JFD, LFD, NFD,
     &                  COST_FUNC, ADJ_FORCE, STT_ADJ
      USE LOGICAL_ADJ_MOD,  ONLY : LPRINTFD, LDCOSAT
      USE ERROR_MOD,        ONLY : IT_IS_NAN, ERROR_STOP
      USE GRID_MOD,         ONLY : GET_IJ

      ! Local Variables
      INTEGER ::  W, I, J, Z, ZZ, L,LL
      INTEGER ::  LON15, IIJJ(2)
      INTEGER ::  NLEV_RET

      TYPE (XPLEX)  ::  RETLEV(Z_DIM+1)
      TYPE (XPLEX)  ::  P_EDGE(Z_DIM+2), MODEL_COL, MOPITT_COL
      TYPE (XPLEX)  ::  UTC, TAU0
      TYPE (XPLEX)  ::  MODEL_P(LLPAR), MODEL_CO_MR(LLPAR)
      TYPE (XPLEX)  ::  COUNT_GRID(IIPAR,JJPAR)
      TYPE (XPLEX)  ::  COUNT(IIPAR,JJPAR)
      TYPE (XPLEX)  ::  MOP_COL_GRID(IIPAR,JJPAR)
      TYPE (XPLEX)  ::  MODEL_COL_GRID(IIPAR,JJPAR)
      TYPE (XPLEX)  ::  NEW_COST(IIPAR,JJPAR)
      TYPE (XPLEX)  ::  ADJ_F(LLPAR)
      TYPE (XPLEX)  ::  SY
      TYPE (XPLEX)  ::  MODEL_P_EDGE(LLPAR+1)

      TYPE (XPLEX), ALLOCATABLE :: GEOS_RAW(:)
      TYPE (XPLEX), ALLOCATABLE :: MOP_CO(:)
      TYPE (XPLEX), ALLOCATABLE :: DIFF_ADJ(:)
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      TYPE (XPLEX), ALLOCATABLE :: GEOS_CO(:)
      TYPE (XPLEX), ALLOCATABLE :: DIFF_COST(:)
#endif
#if defined( MOPITT_V3_CO_OBS )
      TYPE (XPLEX)  ::  DIFF_COST
#endif

      !=================================================================
      ! CALC_MOPITT_FORCE begins here!
      !=================================================================

      TAU0 = GET_TAU0( GET_MONTH(), GET_DAY(), GET_YEAR() )

      COUNT_GRID(:,:)     = 0d0
      COUNT(:,:)          = 0d0
      MOP_COL_GRID(:,:)   = -999.0
      MODEL_COL_GRID(:,:) = -999.0
      ADJ_FORCE(:,:,:,:)  = 0d0
      NEW_COST(:,:)       = 0d0

      !=================================================================
      ! Loop over MOPITT data
      !=================================================================
      DO W = 1, T_DIM

         ! Compute local time:
         ! Local TIME = GMT + ( LONGITUDE / 15 ) since each hour of time
         ! corresponds to 15 degrees of LONGITUDE on the globe
         LON15 = LONGITUDE(W) / 15.
         UTC   = TAU(W) - TAU0 + LON15
         IF (UTC < 0. )  UTC = UTC + 24
         IF (UTC > 24.)  UTC = UTC - 24


         !Only consider day time MOPITT measurements
         ! am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)
#if defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > 70
     &    .and. LONGITUDE(W) < 150
     &    .and. LATITUDE(W)  > -11
     &    .and. LATITUDE(W)  < 55 ) THEN
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -140
     &    .and. LONGITUDE(W) < -40
     &    .and. LATITUDE(W)  > 10
     &    .and. LATITUDE(W)  < 70 ) THEN
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -126
     &    .and. LONGITUDE(W) < -66
     &    .and. LATITUDE(W)  > 13
     &    .and. LATITUDE(W)  < 57 ) THEN
#else
         IF ( UTC >= 4.5 .and. UTC <= 16.5 ) THEN
#endif

         ! Get grid box
         IIJJ  = GET_IJ( LONGITUDE(W), LATITUDE(W))
         I = IIJJ(1)
         J = IIJJ(2)

         !=================================================================
         ! Data selection
         !=================================================================
         IF( GET_HOUR() == OBS_HOUR_MOPITT(I,J) .and.
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
     &      CLOUD_DES(W) == 2.0 .and.
#endif
     &      CO_TOTAL_COLUMN(1,W) > 5E17 .and.
     &      DOMAIN_OBS(I,J) == 1 ) THEN

            RETLEV(:) = -999.0
            MODEL_COL = 0d0
            MOPITT_COL = 0D0

            ! Create pressure profile
#if defined( MOPITT_V3_CO_OBS )
            RETLEV(1) = BOTTOM_PRESSURE_V3(1,W)
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            RETLEV(1) = BOTTOM_PRESSURE_V4(W)
#endif
            ZZ = 0
            ! Loop over Mopitt levels
            DO Z = 1, Z_DIM
            ! Always start from the bottom pressure,
            ! even if it means skipping a MOPITT pressure level
               IF ( PRESSURE(Z) >= RETLEV(1) ) THEN
                  ZZ = ZZ + 1
                  CYCLE
               ENDIF
               ! Save into profile
               RETLEV(Z+1-ZZ) = PRESSURE(Z)
            ENDDO
            NLEV_RET = Z_DIM+1 - ZZ

#if defined( MOPITT_V3_CO_OBS ) .or. defined( MOPITT_V4_CO_OBS )
            DO L = 2, NLEV_RET
               P_EDGE(L) = ( RETLEV(L-1) + RETLEV(L) ) /2.
            ENDDO
            P_EDGE(1) =  RETLEV(1)
            P_EDGE(NLEV_RET+1) =  RETLEV(NLEV_RET)
#endif
#if defined( MOPITT_V5_CO_OBS )
            DO L = 1, NLEV_RET
               P_EDGE(L) = RETLEV(L)
            ENDDO
            P_EDGE(NLEV_RET+1) = 36
#endif

            ALLOCATE( XA( NLEV_RET ) )
            ALLOCATE( T( NLEV_RET ) )
            ALLOCATE( A( NLEV_RET,NLEV_RET ) )
            ALLOCATE( AC( NLEV_RET ) )
            ALLOCATE( MOP_CO( NLEV_RET ) )
            ALLOCATE( GEOS_RAW( NLEV_RET ) )
            ALLOCATE( DIFF_ADJ( NLEV_RET ) )
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            ALLOCATE( GEOS_CO( NLEV_RET ) )
            ALLOCATE( DIFF_COST( NLEV_RET ) )
#endif

            ! MOPITT CO vertical profile
            MOP_CO(1) = CO_RET_BOT_MIXING_RATIO(1,W)
#if defined( MOPITT_V3_CO_OBS )
            MOP_CO(2:NLEV_RET) = CO_MIXING_RATIO(1,8-NLEV_RET:6,W)
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            MOP_CO(2:NLEV_RET) = CO_MIXING_RATIO(1,11-NLEV_RET:9,W)
#endif
            MOP_CO = MOP_CO * 1E-9

            ! COMPUTE AVERAGING KERNEL
            CALL CALC_AVGKER(NLEV_RET, W, RETLEV, MOP_CO)

            !USE MOPITT SURFACE PRESSURE
            !DO L=1, LLPAR + 1
            !   MODEL_P_EDGE(L) = GET_AP(L) + GET_BP(L) * RETLEV(1)
            !ENDDO

            DO L = 1, LLPAR
               !MOPITT PRESSURE LEVEL
               !MODEL_P(L) = (MODEL_P_EDGE(L) + MODEL_P_EDGE(L+1)) / 2

               ! Get GC pressure levels (mbar)
               MODEL_P(L) = GET_PCENTER(I,J,L)

               ! Obtain archieved forward model results
               ! kg -> v/v
               MODEL_CO_MR(L) = CHK_STT(I,J,L,IDTCO) *
     &                            TCVV(IDTCO) / AD(I,J,L)
            ENDDO

            ! Interplote the model to MOPITT vertical grids
            CALL BIN_DATA_V4(MODEL_P, P_EDGE, MODEL_CO_MR(:),
     &            GEOS_RAW, NLEV_RET, 1)

            !=================================================================
            ! Apply MOPITT observation operator
            !=================================================================

            ! Total Column: C = T * XA  + AC * ( Xm - XA )
            ! Stratosphere Levels are removed
#if defined( MOPITT_V3_CO_OBS )
            DO L = 1, NLEV_RET - 1
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            DO L = 1, NLEV_RET - 2
#endif
#if defined( MOPITT_V3_CO_OBS )
               MODEL_COL = MODEL_COL
     &                   + T(L) * XA(L)
     &                   + AC(L) * (GEOS_RAW(L) -XA(L))
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
               MODEL_COL = MODEL_COL
     &                   + T(L) * XA(L)
     &                   + AC(L) * (LOG10(GEOS_RAW(L))
     &                   - LOG10(XA(L)))
#endif
               MOPITT_COL = MOPITT_COL + T(L) * MOP_CO(L)
            ENDDO


#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            GEOS_CO(:) = 0d0
            ! Smoothed Profile: X_hat = XA  + A * ( Xm - XA )
            DO L = 1, NLEV_RET
               DO LL = 1, NLEV_RET
                  GEOS_CO(L) = GEOS_CO(L)
     &                       + A(L,LL)
     &                       * (LOG10( GEOS_RAW(LL) ) - LOG10( XA(LL) ))
               ENDDO
               GEOS_CO(L) = LOG10( XA(L) ) + GEOS_CO(L)
            ENDDO
#endif

            !=================================================================
            ! COST FUNCTION
            !=================================================================
#if defined( MOPITT_V3_CO_OBS )
            SY = ( ERR_PERCENT(I,J) * MOPITT_COL )**2
            DIFF_COST     = MODEL_COL - MOPITT_COL
            NEW_COST(I,J) = NEW_COST(I,J) + (DIFF_COST ** 2) / SY
            COUNT(I,J) = COUNT(I,J) +1
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            DIFF_COST(:) = 0D0
            SY = ERR_PERCENT(I,J) **2
            DO L = 1, NLEV_RET - 2

               DIFF_COST(L)  = GEOS_CO(L) - LOG10( MOP_CO(L) )
               ! Update to be consistent with merged APCOST routine 
               ! (dkh, 01/18/12, adj32_017) 
               !NEW_COST(I,J) = NEW_COST(I,J)
     &         !              + ( DIFF_COST(L)**2 ) / SY
               NEW_COST(I,J) = NEW_COST(I,J)
     &                       + 0.5d0 * ( DIFF_COST(L)**2 ) / SY
               COUNT(I,J) = COUNT(I,J) + 1

            ENDDO
#endif

            !=================================================================
            ! adjoint operator
            !=================================================================
            DIFF_ADJ(:) = 0D0
#if defined( MOPITT_V3_CO_OBS )
            DO L = 1, NLEV_RET - 1
               DIFF_ADJ(L) = DIFF_COST *  AC(L) / SY
            ENDDO
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            DO L = 1, NLEV_RET
               DO LL = 1, NLEV_RET
                  DIFF_ADJ(L) = DIFF_ADJ(L)
     &                        + A(LL,L) * DIFF_COST(LL) / SY
               ENDDO
               ! fwd code: LOG(GEOS_RAW) - LOG(XA)
               DIFF_ADJ(L) = DIFF_ADJ(L) / GEOS_RAW(L)
            ENDDO
#endif

            CALL BIN_DATA_V4( MODEL_P,  P_EDGE, ADJ_F,
     &                        DIFF_ADJ, NLEV_RET, -1   )

            ! adjoint FORCE
            DO L = 1, LLPAR

               !v/v->kg
               ADJ_FORCE(I,J,L,IDTCO) = ADJ_FORCE(I,J,L,IDTCO)
! Update to be consistent with merged APCOST routine (dkh, 01/18/12, adj32_017) 
!     &                            + 2.0D0 * ADJ_F(L) * TCVV(IDTCO)
     &                            + ADJ_F(L) * TCVV(IDTCO)
     &                            / AD(I,J,L)

            ENDDO

            COUNT_GRID(I,J)     = COUNT_GRID(I,J) + 1.d0
            MOP_COL_GRID(I,J)   = MOP_COL_GRID(I,J) + MOPITT_COL
            MODEL_COL_GRID(I,J) = MODEL_COL_GRID(I,J) + MODEL_COL

            IF ( ALLOCATED( GEOS_RAW ) ) DEALLOCATE( GEOS_RAW )
            IF ( ALLOCATED( MOP_CO   ) ) DEALLOCATE( MOP_CO   )
            IF ( ALLOCATED( DIFF_ADJ ) ) DEALLOCATE( DIFF_ADJ )
            IF ( ALLOCATED( A        ) ) DEALLOCATE( A        )
            IF ( ALLOCATED( AC       ) ) DEALLOCATE( AC       )
            IF ( ALLOCATED( T        ) ) DEALLOCATE( T        )
            IF ( ALLOCATED( XA       ) ) DEALLOCATE( XA       )
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
            IF ( ALLOCATED( GEOS_CO  ) ) DEALLOCATE( GEOS_CO  )
            IF ( ALLOCATED( DIFF_COST) ) DEALLOCATE( DIFF_COST)
#endif

         ENDIF !OBS_HOUR

         ENDIF !local time

      ENDDO  !loop over MOPITT data

      !=================================================================
      ! BIN OUTPUT INFO INTO MODEL GRID BOXES
      !=================================================================
      DO I = 1, IIPAR
         DO J = 1, JJPAR

            IF ( COUNT_GRID(I,J) > 0d0 ) THEN

               !The mean value in the grid
               MOP_COL_GRID(I,J)   = MOP_COL_GRID(I,J)
     &                             / COUNT_GRID(I,J)
               MODEL_COL_GRID(I,J) = MODEL_COL_GRID(I,J)
     &                             / COUNT_GRID(I,J)
               ADJ_FORCE(I,J,:,IDTCO)  = ADJ_FORCE(I,J,:,IDTCO)
     &                             / COUNT_GRID(I,J)
               NEW_COST(I,J)       = NEW_COST(I,J) / COUNT_GRID(I,J)
               COUNT(I,J)          = COUNT(I,J) / COUNT_GRID(I,J)

               !Update adjoint tracer
               STT_ADJ(I,J,:,IDTCO) = STT_ADJ(I,J,:,IDTCO) + 
     &                                   ADJ_FORCE(I,J,:,IDTCO)

               ! Diagnostic stuff: FORCING, MOP_MOD_DIFF, MODEL_BIAS
               IF( LDCOSAT )THEN

                  CALL SET_FORCING( I, J, DAY_OF_SIM,
     &                              ADJ_FORCE(I,J,1,IDTCO) )
                  CALL SET_MOP_MOD_DIFF( I, J, DAY_OF_SIM,
     &               MODEL_COL_GRID(I,J) - MOP_COL_GRID(I,J) )

                  CALL SET_MODEL_BIAS( I, J, DAY_OF_SIM, 1,
     &              ( MODEL_COL_GRID(I,J) - MOP_COL_GRID(I,J) ) /
     &                                 MOP_COL_GRID(I,J)           )
                  CALL SET_MODEL     ( I, J, DAY_OF_SIM, 1,
     &                                 MODEL_COL_GRID(I,J)         )
                  CALL SET_OBS       ( I, J, DAY_OF_SIM, 1,
     &                                 MOP_COL_GRID(I,J)           )

                  COST_ARRAY(I,J,DAY_OF_SIM) =
     &               COST_ARRAY(I,J,DAY_OF_SIM) + NEW_COST(I,J)

               ENDIF

               IF ( IT_IS_NAN( NEW_COST(I,J) ) ) THEN
                  PRINT*, 'I=', I, 'J=', J
                  CALL ERROR_STOP( 'NEW_COST is NaN',
     &                             'CALC_MOPITT_FORCE')
               ENDIF

            ENDIF !COUNT_GRID

         ENDDO
      ENDDO

      IF (LPRINTFD)  THEN
         PRINT*, 'IFD, JFD= ', IFD, JFD
         PRINT*, 'MODEL_STT:', MODEL_COL_GRID(IFD,JFD)
         PRINT*, 'OBS_STT:', MOP_COL_GRID(IFD,JFD)
         PRINT*, 'NEW_COST', NEW_COST(IFD,JFD)
         PRINT*, 'ADJ_FORCE:', ADJ_FORCE(IFD,JFD,:,IDTCO)
         PRINT*, 'STT_ADJ:', STT_ADJ(IFD,JFD,:,IDTCO)
      ENDIF

      ! Update cost function
      PRINT*, 'TOTAL NEW_COST = ', SUM(NEW_COST)
      PRINT*, 'COST_FUNC BEFORE ADDING NEW_COST=', COST_FUNC
      COST_FUNC   = COST_FUNC   + SUM ( NEW_COST )
      COUNT_TOTAL = COUNT_TOTAL + SUM ( COUNT    )
      PRINT*, 'Total observation number:', COUNT_TOTAL

      ! Return to calling program
      END SUBROUTINE CALC_MOPITT_FORCE
!--------------------------------------------------------------------------------------------

      SUBROUTINE CALC_AVGKER( NLEV_RET, W, RETLEV, MOP_CO )

!******************************************************************************
! SUBROUTINE CALC_AVGKER construct the averging kernel matrix
! (zhe 1/19/11)
!******************************************************************************

      INTEGER :: ILEV, JLEV, ILEV2, JLEV2, Z, W
      INTEGER :: NLEV_RET
      TYPE (XPLEX)  :: DELP(NLEV_RET)
      TYPE (XPLEX)  :: RETLEV(NLEV_RET)
      TYPE (XPLEX)  :: MOP_CO(NLEV_RET)

#if defined( MOPITT_V3_CO_OBS )
      INTEGER :: IND0
      INTEGER :: INDRETLEV(7)

      TYPE (XPLEX)  :: COVMAT_RET_INT(7,7)
      TYPE (XPLEX)  :: MAT_UNIT(NLEV_RET, NLEV_RET)
      TYPE (XPLEX)  :: COVMAT_RET(NLEV_RET, NLEV_RET)
      TYPE (XPLEX)  :: SX_I_SA(NLEV_RET, NLEV_RET)

      TYPE (XPLEX)  :: COV_CO_AP_7LEV(NLEV_RET, NLEV_RET)
      TYPE (XPLEX)  :: ICOV_CO_AP_7LEV(NLEV_RET, NLEV_RET)

      LOGICAL      SINGLR       !  O  Set TRUE if A is singular
      LOGICAL      FAIL, MOPF   !  O  Set TRUE if a Fatal Error is detected
      CHARACTER*80 ERRMSG       !  O  Error message written if FAIL is TRUE
#endif

#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      TYPE (XPLEX) :: AVGKER_RET(NLEV_RET, NLEV_RET)
      TYPE (XPLEX), PARAMETER :: log10e = LOG10(2.71828183)
#endif

      !=================================================================
      ! CALC_AVGKER begins here!
      !=================================================================

      A(:,:) = 0d0
      AC(:)  = 0d0

#if defined( MOPITT_V3_CO_OBS )

      DO ILEV = 1, NLEV_RET
      DO JLEV = 1, NLEV_RET
         IF ( ILEV == JLEV ) THEN
            MAT_UNIT(ILEV,JLEV) = 1d0
         ELSE
            MAT_UNIT(ILEV,JLEV) = 0d0
         ENDIF
      ENDDO
      ENDDO

      ! INTERPOLATE FROM PLEV_AP >> RETLEV
      CALL INTERP_AP( PLEV_AP, NLEV_AP,  CO_MR_AP,
     &                RETLEV,  NLEV_RET, XA        )

      ! Interpolate a priori covariance on retrieval levels
      DO ILEV = 1, NLEV_RET
         DO ILEV2 = 1, Nlev_AP-1
            IF ( RETLEV(ILEV) >= PLEV_AP(ILEV2) .and.
     &           RETLEV(ILEV) < PLEV_AP(ILEV2+1)      ) THEN
                 IND0 = ILEV2
            ENDIF
         ENDDO
         INDRETLEV(ILEV) = IND0
      ENDDO

      ! Test on 1st level: must be different from second level!!!
      ! (in this case, Sa singular matrix)
      IF ( INDRETLEV(1) == INDRETLEV(2) ) THEN
         INDRETLEV(1) = INDRETLEV(1) + 1
      ENDIF

      DO ILEV = 1, NLEV_RET
      DO JLEV =1, NLEV_RET
            COV_CO_AP_7LEV(ILEV,JLEV) =
     &         COV_CO_AP( INDRETLEV(ILEV),INDRETLEV(JLEV) )
      ENDDO
      ENDDO

      ! Invert covariance matrix
      CALL GAUSSJ(
     I     NLEV_RET,NLEV_RET,NLEV_RET,
     I     COV_CO_AP_7LEV,
     O     ICOV_CO_AP_7LEV,
     O     SINGLR,FAIL,ERRMSG)
      !IF (SINGLR) STOP 'SINGULAR'
      !If Sa singular: use identity matrix
      IF ( SINGLR ) THEN
         DO ILEV = 1, NLEV_RET
         DO JLEV = 1, NLEV_RET
            IF ( ILEV == JLEV ) THEN
               ICOV_CO_AP_7LEV(ILEV,ILEV) = 1.d0
            ELSE
               ICOV_CO_AP_7LEV(ILEV,ILEV) = 0.d0
            ENDIF
         ENDDO
         ENDDO
      ENDIF

      !Retrieval error covariance matrix:
      !Construct complete covariance matrix (in mixing ratio units)
      COVMAT_RET_INT(1, 2:7) = RET_ERR_COV(1, 1:6, W)
      COVMAT_RET_INT(2, 3:7) = RET_ERR_COV(1, 7:11, W)
      COVMAT_RET_INT(3, 4:7) = RET_ERR_COV(1, 12:15, W)
      COVMAT_RET_INT(4, 5:7) = RET_ERR_COV(1, 16:18, W)
      COVMAT_RET_INT(5, 6:7) = RET_ERR_COV(1, 19:20, W)
      COVMAT_RET_INT(6, 7)   = RET_ERR_COV(1, 21, W)
      DO ILEV = 1, 7
         IF (ILEV == 1) THEN
            COVMAT_RET_INT(ILEV, ILEV) =
     &           ( CO_RET_BOT_MIXING_RATIO(2, W)**2 ) * 1.E-18
         ELSE
            COVMAT_RET_INT(ILEV, ILEV) =
     &           ( CO_MIXING_RATIO(2, ILEV-1, W)**2 ) * 1.E-18
         ENDIF
         DO JLEV = 1, ILEV-1
            COVMAT_RET_INT(ILEV, JLEV) = COVMAT_RET_INT(JLEV, ILEV)
         ENDDO
      ENDDO

      !Remove bad levels from covariance matrix (if psurf lt 850)
      !for a priori and retrieval error covariance matrices
      IF ( NLEV_RET < 7 ) THEN
         DO ILEV = 1, NLEV_RET
            ILEV2 = ILEV + ( 7 - NLEV_RET )
            DO JLEV = 1, NLEV_RET
               JLEV2 = JLEV + ( 7 - NLEV_RET )
               COVMAT_RET(ILEV,JLEV) =
     &              COVMAT_RET_INT(ILEV2,JLEV2)
            ENDDO
         ENDDO
         COVMAT_RET(1,1) = COVMAT_RET_INT(1,1)
      ELSE
         COVMAT_RET = COVMAT_RET_INT
      ENDIF

      ! Compute averaging kernel on retrievals levels:
      ! A = I - Sx Sa^-1
      ! Sx and Sa in vmr
      SX_I_SA = MATMUL( COVMAT_RET, ICOV_CO_AP_7LEV )
      A = MAT_UNIT  - SX_I_SA

      ! Convert to column averaging kernel
      DELP(1)        = ( RETLEV(1) - RETLEV(2) ) / 2.d0
      DELP(NLEV_RET) = 159d0
      DO Z = 2, NLEV_RET-1
         DELP(Z) = ( RETLEV(Z-1) - RETLEV(Z) )  / 2d0
     &           + ( RETLEV(Z)   - RETLEV(Z+1)) / 2d0
      ENDDO
      ! transfer function [v/v -> molec/cm2]
      T = 2.12d+22 * DELP

      AC = MATMUL(T,A)

#endif

#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )

      XA(1) = CO_MR_AP_BOTTOM(1, W)
      XA(2:NLEV_RET) = CO_MR_AP_V4(1,11-NLEV_RET:9,W)
      XA = XA * 1E-9

      !Remove bad levels from averging kernel matrix
      IF ( NLEV_RET < 10 ) THEN
         DO ILEV = 1, NLEV_RET
            ILEV2 = ILEV + ( 10 - NLEV_RET )
            DO JLEV =1, NLEV_RET
               JLEV2 = JLEV + ( 10 - NLEV_RET)
               A(ILEV,JLEV) =
     &              AVGKER_V4(ILEV2,JLEV2,W)
            ENDDO
         ENDDO
      ELSE
         A(:,:) = AVGKER_V4(:,:,W)
      ENDIF
#endif
#if defined( MOPITT_V4_CO_OBS )
      DELP(1) = ( RETLEV(1) - RETLEV(2) ) /2d0
      DELP(2) = DELP(1) + 50D0
      DELP(3:NLEV_RET) = 100D0
#endif
#if defined( MOPITT_V5_CO_OBS )
      DELP(1) = RETLEV(1) - RETLEV(2)
      DELP(2:NLEV_RET-1) = 100D0
      DELP(NLEV_RET) = 74D0
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )

      ! transfer function [v/v -> molec/cm2]
      T = 2.12E+22 * DELP

      ! Convert to column averaging kernel
      DO JLEV = 1, NLEV_RET
         DO ILEV = 1, NLEV_RET
            AC(JLEV) = AC(JLEV) + DELP(ILEV) * MOP_CO(ILEV)
     &               * A(ILEV,JLEV)
         ENDDO
         AC(JLEV) = (2.12E+22 / log10e ) * AC(JLEV)
      ENDDO

#endif

      END SUBROUTINE CALC_AVGKER
!------------------------------------------------------------------------------------

      SUBROUTINE BIN_DATA_V4( P_MODEL,  P_EDGE, DATA_MODEL, DATA_MOP,
     &                        NLEV_RET, FB                            )

!******************************************************************************
!Based on the code from Monika.  (zhe 1/19/11)
!FB = 1 for forward
!FB = -1 for adjoint
!******************************************************************************

      INTEGER :: L, LL, FB
      INTEGER :: NLEV_RET, NB
      TYPE (XPLEX)  :: P_MODEL(LLPAR)
      TYPE (XPLEX)  :: DATA_MODEL(LLPAR), DATA_MOP(NLEV_RET), DATA_TEM
      TYPE (XPLEX)  :: P_EDGE(NLEV_RET+1)

      !=================================================================
      ! BIN_DATA_V4 begins here!
      !=================================================================

      IF (FB > 0) THEN
         
         DO L = 1, NLEV_RET
            DO LL = 1, LLPAR
               IF ( P_MODEL(LL) <= P_EDGE(L) ) THEN
                  DATA_MOP(L) = DATA_MODEL(LL)
                  EXIT
               ENDIF
            ENDDO
         ENDDO
                  
         DO L = 1, NLEV_RET
            NB = 0
            DATA_TEM = 0
            DO LL = 1, LLPAR
               IF ( ( P_MODEL(LL) <= P_EDGE(L)) .and.
     &              ( P_MODEL(LL) > P_EDGE(L+1)) ) THEN
                  DATA_TEM = DATA_TEM + DATA_MODEL(LL)
                  NB = NB + 1
               ENDIF
            ENDDO
            IF (NB > 0) DATA_MOP(L) = DATA_TEM / NB
         ENDDO

      ELSE

         DATA_MODEL(:) = 0.
         DO L = 1, LLPAR
            DO LL = 1, NLEV_RET
               IF ( ( P_MODEL(L) <= P_EDGE(LL)) .and.
     &              ( P_MODEL(L) > P_EDGE(LL+1)) ) THEN
                  DATA_MODEL(L) = DATA_MOP(LL)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


      ! Return to calling program
      END SUBROUTINE BIN_DATA_V4
!-----------------------------------------------------------------------------------

      SUBROUTINE INIT_DOMAIN

!******************************************************************************
!Define the observatio region
!******************************************************************************
#     include "CMN_SIZE"   ! Size parameters

      !local variables
      INTEGER :: I, J

      !=================================================================
      ! INIT_DOMAIN begins here!
      !=================================================================

      DOMAIN_OBS(:,:) = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR

#if   defined( GRID05x0666 )
!     The surrounding region is used as cushion
!     (zhe 11/28/10)
         IF ( J >= 8 .and. J <= JJPAR-7 .and.
     &        I >= 7 .and. I <= IIPAR-6
#elif defined( GRID2x25 )
         IF ( J >= 16 .and. J <= 76   !60S-60N
#elif defined( GRID4x5 )
         IF ( J >= 9 .and. J <= 39    !60S-60N
!         IF ( J >= 21 .and. J <= 37 .and.  !11S-55N
!     &        I >= 51 .and. I <= 67        !70E-150E
#endif
     &       ) DOMAIN_OBS(I,J) = 1d0

      ENDDO
      ENDDO

      PRINT*, sum(DOMAIN_obs), 'MAX observations today'

      END SUBROUTINE INIT_DOMAIN

!-----------------------------------------------------------------------------

      SUBROUTINE CALC_OBS_HOUR

!***************************************************************************
! Subroutine CALC_OBS_HOUR computes an array of hours for each day of obs.
! If there is an obs in a particular gridbox on that day, it assigns the
! hour (0..23). If there isn't, OBS_HOUR stays initialized to -1.
! (mak, 12/14/05)
!***************************************************************************

      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY,
     &                         GET_YEAR, GET_HOUR
      USE GRID_MOD,     ONLY : GET_IJ

#     include "CMN_SIZE"

      TYPE (XPLEX)   :: OBS_HOUR(IIPAR,JJPAR)
      TYPE (XPLEX)   :: TAU0, UTC
      INTEGER  :: W, I, J
      INTEGER  :: LON15, IIJJ(2)
      INTEGER  :: COUNT_GRID(IIPAR,JJPAR)

      !=================================================================
      ! CALC_OBS_HOUR begins here!
      !=================================================================

      ! Get TAU0 from the date (at 0GMT)
      TAU0 = GET_TAU0(GET_MONTH(), GET_DAY(), GET_YEAR())

      OBS_HOUR_MOPITT(:,:) = -1
      OBS_HOUR(:,:)        = 0
      COUNT_GRID(:,:)      = 0

      DO W = 1, T_DIM

         ! Compute local time:
         ! Local TIME = GMT + ( LONGITUDE / 15 ) since each hour of time
         ! corresponds to 15 degrees of LONGITUDE on the globe
         !============================================================
         LON15 = LONGITUDE(W) / 15d0
         UTC   = TAU(W) - TAU0 + LON15
         IF ( UTC < 0d0  )  UTC = UTC + 24
         IF ( UTC > 24d0 )  UTC = UTC - 24

         !Only consider day time MOPITT measurements
         !am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)

#if defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > 70
     &    .and. LONGITUDE(W) < 150
     &    .and. LATITUDE(W)  > -11
     &    .and. LATITUDE(W)  < 55 ) THEN
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -140
     &    .and. LONGITUDE(W) < -40
     &    .and. LATITUDE(W)  > 10
     &    .and. LATITUDE(W)  < 70 ) THEN
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -126
     &    .and. LONGITUDE(W) < -66
     &    .and. LATITUDE(W)  > 13
     &    .and. LATITUDE(W)  < 57 ) THEN
#else
         IF ( UTC >= 4.5 .and. UTC <= 16.5 ) THEN
#endif

         ! Get grid box of current record
            IIJJ  = GET_IJ( LONGITUDE(W), LATITUDE(W))
            I = IIJJ(1)
            J = IIJJ(2)

            ! If there's an obs, calculate the time
            IF ( CO_TOTAL_COLUMN(1,W) > 0d0 ) THEN

               COUNT_GRID(I,J) = COUNT_GRID(I,J) + 1d0
               !Add the time of obs, to be averaged and floored later
               OBS_HOUR(I,J) = OBS_HOUR(I,J) + MOPITT_GMT(W)

            ENDIF
         ENDIF
      ENDDO

      ! average obs_hour on the grid
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( COUNT_GRID(I,J) > 0d0 ) THEN

            OBS_HOUR_MOPITT(I,J) =
     &         FLOOR( OBS_HOUR(I,J) / COUNT_GRID(I,J) )

         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE CALC_OBS_HOUR

!----------------------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_MOPITT_OBS( ) RESULT( FLAG )

!******************************************************************************
!  Function ITS_TIME_FOR_MOPITT_OBS returns TRUE if there are observations
!  available for particular time (hour of a particular day) based on
!  the OBS_HOUR_MOPITT array which holds the hour of obs in each gridbox
!  (computed when file read in mop02_mod.f) (mak, 7/12/07)
!******************************************************************************

      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE

#     include "CMN_SIZE"  ! Size params

      ! Function value
      LOGICAL :: FLAG

      INTEGER :: I,J

      !=================================================================
      ! ITS_TIME_FOR_MOPITT_OBS begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.

      DO J = 1,JJPAR
      DO I = 1,IIPAR
         IF( GET_HOUR()   == OBS_HOUR_MOPITT(I,J) .and.
     &       GET_MINUTE() == 0                          ) THEN

               PRINT*, 'obs_hour was', get_hour(), 'in box', I, J
               FLAG = .TRUE.

               !GOTO 11
               RETURN

         ENDIF
      ENDDO
      ENDDO

      END FUNCTION ITS_TIME_FOR_MOPITT_OBS

!----------------------------------------------------------------------------

      SUBROUTINE READ_MOP02( FILENAME )

!******************************************************************************
!  Subroutine READ_MOP02 allocates all module arrays and reads data into
!  them from the HDF file. (bmy, 7/2/03, zhe 1/19/11)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of MOPITT file to read
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
      USE HdfSdModule
      USE HdfVdModule
      USE ERROR_MOD, ONLY : ALLOC_ERR

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

      ! Local variables
      INTEGER  :: sId, vId, vSize, nDims, dims(4), as
      INTEGER  :: i, year, month, day
      TYPE (XPLEX)   :: TAU0

      !=================================================================
      ! Mop02Read begins here!
      !=================================================================

      ! Deallocate arrays
      CALL CLEANUP_MOP02

      ! Get date from filename (next to the '-' character)
      i    = INDEX( FILENAME, '-' )
      READ( FILENAME(i+1:i+4), '(i4)' ) year
      READ( FILENAME(i+5:i+6), '(i2)' ) month
      READ( FILENAME(i+7:i+8), '(i2)' ) day

      ! Get TAU0 from the date (at 0GMT)
      TAU0 = getTauFromDate( year, month, day )

      ! Open file for HDF-VDATA interface
      CALL vdOpen( FILENAME )

      !=================================================================
      ! VDATA field: Time (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Seconds in Day', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate arrays
      ALLOCATE( SECONDS_IN_DAY( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SECONDS_IN_DAY' )

      ALLOCATE( TAU( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'TAU' )

      ALLOCATE( MOPITT_GMT( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'MOPITT_GMT' )

      ! Read data
      CALL vdGetData( vId, vSize, SECONDS_IN_DAY )

      ! Close field
      CALL vdCloseField( vId )

      ! Compute GMT of MOPITT observations
      MOPITT_GMT = ( DCMPLX( SECONDS_IN_DAY ) / 3600d0 )

      ! Compute TAU values for GAMAP from SECONDS_IN_DAY
      TAU       = MOPITT_GMT + TAU0

      ! Save time dimension in T_DIM
      T_DIM      = vSize

      !=================================================================
      ! VDATA field: LONGITUDE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Longitude', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( LONGITUDE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LONGITUDE' )

      ! Read data
      CALL vdGetData( vId, vSize, LONGITUDE )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: LATITUDE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Latitude', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( LATITUDE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LATITUDE' )

      ! Read data
      CALL vdGetData( vId, vSize, LATITUDE )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: Cloud Description (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Cloud Description', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( CLOUD_DES( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CLOUD_DES' )

      ! Read data
      CALL vdGetData( vId, vSize, CLOUD_DES )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: Surface Index (1-D)
      !=================================================================

      ! Open field for reading
#if defined( MOPITT_V3_CO_OBS )
      CALL vdOpenField( 'Surface Indicator', vId )
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      CALL vdOpenField( 'Surface Index', vId )
#endif
      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( SURFACE_INDEX( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SURFACE_INDEX' )

      ! Read data
      CALL vdGetData( vId, vSize, SURFACE_INDEX )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: PRESSURE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Pressure Grid', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( PRESSURE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'PRESSURE' )

      ! Read data
      CALL vdGetData( vId, vSize, PRESSURE )

      ! Close field
      CALL vdCloseField( vId )

      ! Save PRESSURE dimension in Z_DIM
      Z_DIM = vSize

#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      !=================================================================
      ! VDATA field: Retrieval Bottom Pressure (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Surface Pressure', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( BOTTOM_PRESSURE_V4( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'BOTTOM_PRESSURE_V4' )

      ! Read data
      CALL vdGetData( vId, vSize, BOTTOM_PRESSURE_V4 )

      ! Close field
      CALL vdCloseField( vId )

#endif

      ! Close HDF-VDATA interface
      CALL vdClose( FILENAME )


      ! Open file for HDF-SDATA interface
      CALL sdOpen( FILENAME )

#if defined( MOPITT_V3_CO_OBS )
      !=================================================================
      ! SDATA field: Retrieval Bottom Pressure (2-D)
      !=================================================================
      ! Open field for reading
      CALL sdOpenFieldByName( 'Retrieval Bottom Pressure', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( BOTTOM_PRESSURE_V3( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'BOTTOM_PRESSURE_V3' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), BOTTOM_PRESSURE_V3 )

      ! Close field
      CALL sdCloseField( sId )
#endif

      !=================================================================
      ! SDATA field: CO Mixing Ratio (3-D)
      !=================================================================

      ! Open field
#if defined( MOPITT_V3_CO_OBS )
      CALL sdOpenFieldByName( 'CO Mixing Ratio', sId )
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      CALL sdOpenFieldByName( 'Retrieved CO Mixing Ratio Profile', sId )
#endif

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MIXING_RATIO( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MIXING_RATIO' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), CO_MIXING_RATIO )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: CO Retrieval Bottom Mixing Ratio (3-D)
      !=================================================================

      ! Open field
#if defined( MOPITT_V3_CO_OBS )
      CALL sdOpenFieldByName( 'Retrieval Bottom CO Mixing Ratio', sId )
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      CALL sdOpenFieldByName( 'Retrieved CO Surface Mixing Ratio', sId )
#endif
      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_RET_BOT_MIXING_RATIO( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_RET_BOT_MIXING_RATIO' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_RET_BOT_MIXING_RATIO )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: CO Total Column (2-D)
      !=================================================================

      ! Open field
#if defined( MOPITT_V3_CO_OBS )
      CALL sdOpenFieldByName( 'CO Total Column', sId )
#endif
#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      CALL sdOpenFieldByName( 'Retrieved CO Total Column', sId )
#endif
      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_TOTAL_COLUMN( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_TOTAL_COLUMN' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_TOTAL_COLUMN )

      ! Close field
      CALL sdCloseField( sId )

#if defined( MOPITT_V3_CO_OBS )
      !=================================================================
      ! SDATA field: Retrieval Error Covariance Matrix (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieval Error Covariance Matrix', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( RET_ERR_COV( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'RET_ERR_COV' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), RET_ERR_COV )

      ! Close field
      CALL sdCloseField( sId )
#endif

#if defined( MOPITT_V4_CO_OBS ) .or. defined( MOPITT_V5_CO_OBS )
      !=================================================================
      ! SDATA field: Retrieval Averaging Kernel Matrix (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieval Averaging Kernel Matrix', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( AVGKER_V4( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'AVGKER_V4' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), AVGKER_V4 )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: A Priori CO Mixing Ratio Profile (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'A Priori CO Mixing Ratio Profile', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MR_AP_V4( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP_V4' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), CO_MR_AP_V4 )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: A Priori CO Surface Mixing Ratio (2-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'A Priori CO Surface Mixing Ratio', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MR_AP_BOTTOM( dims(1), dims(2)), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP_BOTTOM' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_MR_AP_BOTTOM )

      ! Close field
      CALL sdCloseField( sId )
#endif

      ! Close file and quit
      CALL sdClose( FILENAME )

      ! Return to calling program
      END SUBROUTINE READ_MOP02

!------------------------------------------------------------------------------------

      SUBROUTINE APRIORI_MOP02

!******************************************************************************
! Calculates averaging kernels
! from a MOPITT retrieved covariance matrix on 7-levels,
! using apriori profile and covariance matrix for v3 MOPITT retrievals
!******************************************************************************

      CHARACTER(LEN=255) AP_FILENAME1
      CHARACTER(LEN=255) AP_FILENAME
      CHARACTER(LEN=80) LINE
      INTEGER :: ILEV,JLEV
      TYPE (XPLEX), ALLOCATABLE :: DATA1(:)

      ! Read CO a priori profile and covariance matrix
      AP_FILENAME1 = 'mopitt_v3_apriori.dat'
      AP_FILENAME = TRIM(AP_FILENAME1)

      OPEN ( FILE = TRIM(AP_FILENAME),
     &     UNIT = 20)
!     &     ,form='unformatted')
      READ(20,101) LINE
      READ(20,101) LINE
      READ(20,*) NLEV_AP
!      print*,'NLEV_AP=',NLEV_AP
      ALLOCATE( PLEV_AP ( NLEV_AP ) )
      ALLOCATE( CO_MR_AP ( NLEV_AP ) )
      ALLOCATE( CH4_MR_AP ( NLEV_AP ) )
      ALLOCATE( COV_CO_AP ( NLEV_AP , NLEV_AP) )
      READ(20,101) LINE
      READ(20,*) PLEV_AP
      !PRINT*, PLEV_AP
      READ(20,101) LINE
      READ(20,*) CO_MR_AP
      READ(20,101) LINE
      READ(20,*) CH4_MR_AP
      READ(20,101) LINE
      DO ILEV= 1, NLEV_AP
         ALLOCATE( DATA1 ( ILEV ) )
         READ(20,101) LINE
         READ(20,*) DATA1
         COV_CO_AP(ILEV,1:ILEV) = DATA1
         COV_CO_AP(1:ILEV,ILEV) = DATA1
         IF ( ALLOCATED( DATA1 ) ) DEALLOCATE( DATA1 )
      ENDDO

      CLOSE(20)

 101  FORMAT(A80)

      END SUBROUTINE APRIORI_MOP02
!-----------------------------------------------------------------------------------------

      SUBROUTINE READ_ERROR_VARIANCE
!
!******************************************************************************
!  Subroutine READ_ERROR_VARIANCE reads observation error from binary punch files
!  (zhe 4/20/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TIME_MOD,   ONLY : GET_TAUb

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_ERROR_VARIANCE begins here!
      !=================================================================

      ! Filename
        FILENAME = TRIM( 'OBS_ERR_' ) // GET_RES_EXT()

      ! Echo some information to the standard output
        WRITE( 6, 110 ) TRIM( FILENAME )
 110    FORMAT( '     - READ_ERROR_VARIANCE: Reading ERR_PERCENT
     &                from: ', a )

      ! Read data from the binary punch file
        CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,
     &           GET_TAUb(),    IGLOB,     JGLOB,
     &           1,  ERR_PERCENT,  QUIET=.TRUE. )

      ! Return to calling program
      END SUBROUTINE READ_ERROR_VARIANCE

!------------------------------------------------------------------------------

      SUBROUTINE INFO_MOP02( FILENAME )
!
!******************************************************************************
!  Subroutine INFO_MOP02 Info prints info about all VDATA and SDATA fields
!  contained within the MOPITT HDF file. (bmy, 7/3/03, 4/27/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of MOPITT file to read
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE HdfSdModule
      USE HdfVdModule

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

      !=================================================================
      ! INFO_MOP02 begins here!
      !=================================================================

      ! Print HDF-VDATA variables
      CALL vdOpen( FILENAME )
      CALL vdPrintInfo
      CALL vdClose( FILENAME )

      ! Print HDF-SDATA variables
      CALL sdOpen( FILENAME )
      CALL sdPrintInfo
      CALL sdClose( FILENAME )

      ! Return to calling program
      END SUBROUTINE INFO_MOP02

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MOP02
!
!******************************************************************************
!  Subroutine CLEANUP_MOP02 deallocates all module arrays (bmy, 4/27/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_MOP02 begins here!
      !=================================================================
      IF ( ALLOCATED( LATITUDE        ) ) DEALLOCATE( LATITUDE        )
      IF ( ALLOCATED( LONGITUDE       ) ) DEALLOCATE( LONGITUDE       )
      IF ( ALLOCATED( PRESSURE        ) ) DEALLOCATE( PRESSURE        )
      IF ( ALLOCATED( CLOUD_DES       ) ) DEALLOCATE( CLOUD_DES       )
      IF ( ALLOCATED( SURFACE_INDEX   ) ) DEALLOCATE( SURFACE_INDEX   )
      IF ( ALLOCATED( TAU             ) ) DEALLOCATE( TAU             )
      IF ( ALLOCATED( SECONDS_IN_DAY  ) ) DEALLOCATE( SECONDS_IN_DAY  )
      IF ( ALLOCATED( MOPITT_GMT      ) ) DEALLOCATE( MOPITT_GMT      )

      IF ( ALLOCATED( BOTTOM_PRESSURE_V3)) THEN
         DEALLOCATE( BOTTOM_PRESSURE_V3)
      ENDIF
      IF ( ALLOCATED( BOTTOM_PRESSURE_V4)) THEN
         DEALLOCATE( BOTTOM_PRESSURE_V4)
      ENDIF
      IF ( ALLOCATED( CO_MIXING_RATIO   )) DEALLOCATE( CO_MIXING_RATIO)
      IF ( ALLOCATED( CO_RET_BOT_MIXING_RATIO)) THEN
         DEALLOCATE( CO_RET_BOT_MIXING_RATIO )
      ENDIF
      IF ( ALLOCATED( CO_TOTAL_COLUMN ) ) DEALLOCATE( CO_TOTAL_COLUMN )
      IF ( ALLOCATED( RET_ERR_COV     ) ) DEALLOCATE( RET_ERR_COV     )
      IF ( ALLOCATED( AVGKER_V4       ) ) DEALLOCATE( AVGKER_V4       )


      IF ( ALLOCATED( PLEV_AP         ) ) DEALLOCATE( PLEV_AP         )
      IF ( ALLOCATED( CO_MR_AP        ) ) DEALLOCATE( CO_MR_AP        )
      IF ( ALLOCATED( CH4_MR_AP       ) ) DEALLOCATE( CH4_MR_AP       )
      IF ( ALLOCATED( CO_MR_AP_V4     ) ) DEALLOCATE( CO_MR_AP_V4     )
      IF ( ALLOCATED( CO_MR_AP_BOTTOM ) ) DEALLOCATE( CO_MR_AP_BOTTOM )
      IF ( ALLOCATED( COV_CO_AP       ) ) DEALLOCATE( COV_CO_AP       )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MOP02

!---------------------------------------------------------------------------------------------------


      END MODULE MOPITT_OBS_MOD
