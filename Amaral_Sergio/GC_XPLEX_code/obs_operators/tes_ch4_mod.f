!$Id: tes_ch4_mod.f,v 1.2 2012/03/01 23:27:52 daven Exp $
      MODULE TES_CH4_MOD
!
!******************************************************************************
!  Module TES_CH4_MOD contains variables and routines which are used to 
!  assimilate real or simulated TES CH4 observations. The module is based on 
!  TES_NH3_MOD (kjw, 7/06/11)
!  Added to adj32_023 (dkh, 02/12/12) 
!******************************************************************************
!


      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
 
      ! Parameters
      INTEGER, PARAMETER           :: MAXLEV = 67
      INTEGER, PARAMETER           :: MAXTES = 2000
      TYPE (XPLEX), PARAMETER            :: ERR_PPB = xplex(40.0,0d0) !Stddev in TES obs
      !LOGICAL                      :: LTES_PSO = .TRUE.

      ! Module Variables
      TYPE (XPLEX)                       :: BIAS_PPB

      ! Record to store data from each TES obs
      TYPE TES_CH4_OBS 
         INTEGER                      :: LTES(1)
         TYPE (XPLEX)                       :: LAT(1)
         TYPE (XPLEX)                       :: LON(1)
         TYPE (XPLEX)                       :: TIME(1)
         TYPE (XPLEX)                       :: ERR(1)
         TYPE (XPLEX)                       :: CH4(MAXLEV)
         TYPE (XPLEX)                       :: GC_CH4(MAXLEV)
         TYPE (XPLEX)                       :: PRES(MAXLEV)
         TYPE (XPLEX)                       :: PRIOR(MAXLEV)
         TYPE (XPLEX)                       :: AVG_KERNEL(MAXLEV,MAXLEV)
         TYPE (XPLEX)                       :: S_OER(MAXLEV,MAXLEV)
         TYPE (XPLEX)                       :: S_OER_INV(MAXLEV,MAXLEV)
      ENDTYPE TES_CH4_OBS  

      TYPE(TES_CH4_OBS)                          :: TES(MAXTES)


      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_TES_CH4_OBS( YYYYMMDD, NTES )
!
!******************************************************************************
!  Subroutine READ_TES_CH4_OBS reads the file and passes back info contained
!  therein. (dkh, 02/19/09) 
! 
!  Based on READ_TES_NH3 OBS (dkh, 04/26/10) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD    INTEGER : Current year-month-day
!
!  Arguments as Output: 
!  ============================================================================
!  (1 ) NTES      (INTEGER) : Number of TES retrievals for current day 
!
!  Module variable as Output: 
!  ============================================================================
!  (1 ) TES    (TES_CH4_OBS) : TES retrieval for current day 
!     
!  NOTES:
!  (1 ) Add calculation of S_OER_INV, though eventually we probably want to 
!        do this offline. (dkh, 05/04/10) 
!  (2 ) Now read data files in BPCH format for better compatibility with
!        the standard GEOS-Chem distribution. (kjw, 06/05/10)
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE DIRECTORY_ADJ_MOD,      ONLY : DIAGADJ_DIR
      USE TIME_MOD,               ONLY : EXPAND_DATE
      USE BPCH2_MOD,              ONLY : READ_BPCH2, GET_TAU0
      USE BPCH2_MOD,              ONLY : OPEN_BPCH2_FOR_READ
      USE TIME_MOD,               ONLY : GET_YEAR, GET_MONTH, GET_DAY
      USE FILE_MOD,               ONLY : IU_FILE
      USE ADJ_ARRAYS_MOD,         ONLY : N_CALC, EXPAND_NAME
      USE LOGICAL_ADJ_MOD,        ONLY : LTES_PSO

      ! From READ_BPCH2
      USE FILE_MOD,               ONLY : IU_FILE, IOERROR


      ! Arguments
      INTEGER,            INTENT(IN)  :: YYYYMMDD
      INTEGER,            INTENT(OUT) :: NTES
    
      ! local variables 
      INTEGER                         :: FID
      INTEGER                         :: LTES
      INTEGER                         :: NT
      INTEGER                         :: YYYY, MM, DD
      INTEGER                         :: START
      CHARACTER(LEN=5)                :: TMP
      CHARACTER(LEN=255)              :: READ_FILENAME
      CHARACTER(LEN=255)              :: FILENAME
      LOGICAL                         :: file_exist

      TYPE (XPLEX), PARAMETER   :: FILL = xplex(-999.0D0,0d0)
      TYPE (XPLEX), PARAMETER  :: TOL  = xplex(1d-04,0d0)
      TYPE (XPLEX)                          :: U(MAXLEV,MAXLEV)
      TYPE (XPLEX)                          :: VT(MAXLEV,MAXLEV)
      TYPE (XPLEX)                          :: S(MAXLEV)
      TYPE (XPLEX)                          :: TMP1
      TYPE (XPLEX)                          :: XTAU
      TYPE (XPLEX)                          :: TEST(MAXLEV,MAXLEV)

      ! From READ_BPCH2
      INTEGER                         :: I, J, L, LL, IOS
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      TYPE (XPLEX)             :: LONRES,    LATRES
      TYPE (XPLEX)             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT     
      CHARACTER(LEN=40)  :: RESERVED
      

      !Arrays in which to read BPCH files
      TYPE (XPLEX)               :: DUMMY_NTES(1)
      TYPE (XPLEX)               :: DUMMY_0D(MAXTES)
      TYPE (XPLEX)               :: DUMMY_1D(MAXTES,MAXLEV,1)
      TYPE (XPLEX)               :: DUMMY_2D(MAXTES,MAXLEV,MAXLEV)


      !=================================================================
      ! READ_TES_CH4_OBS begins here!
      !=================================================================

      ! filename root 
      READ_FILENAME = TRIM( 'tes_ch4_YYYYMMDD.bpch' )
      !READ_FILENAME = TRIM( 'temp_test.bpch' )

      ! Expand date tokens in filename 
      CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 ) 

      ! Construct complete filename 
      READ_FILENAME = TRIM( '/home/kjw/TES/data/V004/bpch/' ) // 
     &                TRIM( READ_FILENAME )

      INQUIRE( FILE=READ_FILENAME, exist=file_exist )


      ! If there is no observation file for this day,
      ! Return to calling program
      IF ( .not. file_exist ) THEN
         WRITE(6,*) '    - READ_TES_CH4_OBS: file does not exist: ',
     &                              TRIM( READ_FILENAME )
         WRITE(6,*) '                        no observations today.'

         ! Set NTES = 0 and Return to calling program
         NTES = 0
         RETURN
      ENDIF


      ! Start Reading Data from BPCH
      WRITE(6,*) '    - READ_TES_CH4_OBS: reading file: ', 
     &                               TRIM( READ_FILENAME )

      ! Read variables from bpch instead of netCDF (kjw, 06/05/10)
      !  1. Open BPCH file for today if it exists. If it doesn't, 
      !       return NTES = 0.
      !  2. Read nTES (tracer=1) from bpch
      !      a. read LTES from bpch, store in TES struct
      !      b. read remaining 0-d data, store in struct
      !      c. read 1-d data, store in TES struct
      !      d. read 2-d data, store in TES struct

      ! READ nTES from BPCH.  Tracer numbers correspond to the following
      ! variables in the TES BPCH files:
      !    Tracer #              Variable
      !       1               targets (# TES obs in file)
      !       2               LTES (# good vertical levels in each obs)
      !       3               Longitude
      !       4               Latitude
      !       5               YYYYMMDD
      !       6               Species
      !       7               Pressure
      !       8               Constraint Vector
      !       9               GEOS-Chem_obs
      !      10               Averaging Kernel
      !      11               Inverse of Observation Error Covar Matrix
      !---------------------------------------------------------------

      ! Tau for the bpch file
      YYYY = INT( floor( YYYYMMDD / 1d4 ) )
      MM   = INT( floor( YYYYMMDD - 1d4*YYYY ) / 1d2 )
      DD   = NINT( YYYYMMDD - 1d4*YYYY - 1d2*MM )
      XTAU = GET_TAU0( MM, DD, YYYY )

      ! Number of TES observations in the file
      WRITE(6,*) '         - Reading: NTES ... '
      print*,'XTAU = ',XTAU
      CALL READ_BPCH2( TRIM(READ_FILENAME), 'IJ-AVG-$',     1,
     &                 XTAU,                   1,     1,
     &                    1,          DUMMY_NTES(1), QUIET=.TRUE. )
      NTES = INT( DUMMY_NTES(1) )
      print*, '         - Found # obs today: NTES = ,', NTES


      !==================================================================
      ! Read data for each TES observation in the current day.
      ! Do NOT use READ_BPCH2 because output dimensions limited size
      !    of global 1x1 grid.
      ! The following lines are modified from READ_BPCH2 (kjw, 07/22/10)
      !
      ! 0-D Data
      !    2. # good vertical levels for each obs.
      !    3. longitude
      !    4. latitude
      !    5. mmdd.frac-of-day
      ! 1-D Data
      !    6. Species (CH4)
      !    7. Pressure
      !    8. Constraint Vector
      !    9. GEOS-Chem_obs
      ! 2-D Data
      !   10. Averaging Kernel
      !   11. Inverse of Observation Error Covariance Matrix
      !==================================================================


      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, READ_FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS ) 
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'tes_ch4_mod:1')

         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'tes_ch4_mod:2' )

         ! Place array into DUMMY_2D
         DUMMY_2D(:,:,:) = 0d0
         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( DUMMY_2D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'tes_ch4_mod:3' )


         ! Test for a match
         IF ( 'IJ-AVG-$' == TRIM( CATEGORY ) .and. XTAU  == ZTAU0 ) THEN

            ! LTES
            IF     ( NTRACER == 2 ) THEN
            WRITE(6,*) '         - Reading: LTES ... '
               TES(1:NTES)%LTES(1) = DUMMY_2D(1:NTES,1,1)

            ! Longitude
            ELSEIF ( NTRACER == 3 ) THEN
            WRITE(6,*) '         - Reading: Longitude ... '
               TES(1:NTES)%LON(1)  = DUMMY_2D(1:NTES,1,1)

            ! Latitude
            ELSEIF ( NTRACER == 4 ) THEN
            WRITE(6,*) '         - Reading: Latitude ... '
               TES(1:NTES)%LAT(1)  = DUMMY_2D(1:NTES,1,1)

            ! MMDD.frac-of-day
            ELSEIF ( NTRACER == 5 ) THEN
            WRITE(6,*) '         - Reading: Frac-of-day ... '
               TES(1:NTES)%TIME(1) = DUMMY_2D(1:NTES,1,1) + 
     &                                             GET_YEAR()*1d4

            ! Species (CH4)
            ELSEIF ( NTRACER == 6 ) THEN
            WRITE(6,*) '         - Reading: CH4 ... '
            DO NT=1,NTES
               LTES  = TES(NT)%LTES(1)
               START = MAXLEV - LTES + 1
               TES(NT)%CH4(1:LTES)  = DUMMY_2D(NT,START:MAXLEV,1)
            ENDDO

            ! Pressure
            ELSEIF ( NTRACER == 7 ) THEN
            WRITE(6,*) '         - Reading: Pressure ... '
            DO NT=1,NTES
               LTES  = TES(NT)%LTES(1)
               START = MAXLEV - LTES + 1
               TES(NT)%PRES(1:LTES) = DUMMY_2D(NT,START:MAXLEV,1)
            ENDDO

            ! Constraint Vector
            ELSEIF ( NTRACER == 8 ) THEN
            WRITE(6,*) '         - Reading: Constraint Vector ... '
            DO NT=1,NTES
               LTES  = TES(NT)%LTES(1)
               START = MAXLEV - LTES + 1
               TES(NT)%PRIOR(1:LTES) = DUMMY_2D(NT,START:MAXLEV,1)
            ENDDO

!            ! Kind of Useless now that LTES_PSO created, kjw 07/25/10
!            ! GEOS-Chem Obs
            ELSEIF ( NTRACER == 9 ) THEN
            WRITE(6,*) '         - Reading: GEOS-Chem Obs ... '
            DO NT=1,NTES
               LTES  = TES(NT)%LTES(1)
               START = MAXLEV - LTES + 1
               TES(NT)%GC_CH4(1:LTES) = DUMMY_2D(NT,START:MAXLEV,1)
            ENDDO

            ! Averaging Kernel
            ELSEIF ( NTRACER == 10) THEN
            WRITE(6,*) '         - Reading: Averaging Kernel ... '
            DO NT=1,NTES
               LTES  = TES(NT)%LTES(1)
               START = MAXLEV - LTES + 1
               TES(NT)%AVG_KERNEL(1:LTES,1:LTES) =  
     &                  DUMMY_2D(NT,START:MAXLEV,START:MAXLEV)
            ENDDO

            ! Inverse of Observation Error Covariance Matrix
            ELSEIF ( NTRACER == 11) THEN
            WRITE(6,*) '         - Reading: S_OER_INV ... '
            DO NT=1,NTES
               LTES  = TES(NT)%LTES(1)
               START = MAXLEV - LTES + 1
               TES(NT)%S_OER_INV(1:LTES,1:LTES) =  
     &                  DUMMY_2D(NT,START:MAXLEV,START:MAXLEV)
            ENDDO


            ENDIF  ! If tracer == # 

         ENDIF  ! If Category and Tau match

      ENDDO

      ! Close today's BPCH file of TES observations
      CLOSE( IU_FILE )


      ! Read Errors and populate TES%GC_CH4 if using pseudo-obs
      IF ( LTES_PSO ) THEN

         ! Make pseudo observations and save in TES%GC_CH4
         CALL MAKE_PSEUDO_OBS( YYYYMMDD, NTES )

      ENDIF

!      ! Save AK and S_OER_INV for one observation.
!      ! The plot to make sure they are correct order
!      WRITE(6,'(a)') ' Writing AK and S_OER_INV files'
!         FILENAME = 'test_ak.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         print*,'FILENAME1 = ',FILENAME
!         FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( IU_FILE,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!      ! Save observation # 600
!         LTES = TES(600)%LTES(1)
!         print*,'LTES of obs # 600 =   ',LTES
!         DO L=1,LTES
!           WRITE(IU_FILE,'(65F16.12)') (TES(600)%AVG_KERNEL(L,LL),
!     &                                             LL=1,LTES)
!         ENDDO
!         CLOSE(IU_FILE)
!
!         FILENAME = 'test_s_obs.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         print*,'FILENAME2 = ',FILENAME
!         FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 189,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!
!      ! Save observation # 600
!         LTES = TES(600)%LTES(1)
!         print*,'LTES of obs # 600 =   ',LTES
!         DO L=1,LTES
!           WRITE(IU_FILE,'(65F16.12)') (TES(600)%S_OER_INV(L,LL),
!     &                                             LL=1,LTES)
!         ENDDO
!      WRITE(6,'(a)') ' Done writing AK and S_OER_INV files'
!
!      CLOSE(IU_FILE)


      ! Check reading against values read from BPCH in IDL
      !print*,'TES(600)%LTES = ',TES(600)%LTES(1)
      !print*,'TES(600)%LON  = ',TES(600)%LON(1)
      !print*,'TES(600)%LAT  = ',TES(600)%LAT(1)
      !print*,'TES(600)%TIME = ',TES(600)%TIME(1)
      !print*,'TES(600)%PRES = ',TES(600)%PRES
      !print*,'TES(600)%CH4  = ',TES(600)%CH4
      !print*,'TES(600)%AK   = ',TES(600)%AVG_KERNEL(1:4,1)
      !print*,'TES(600)%ERR  = ',TES(600)%ERR(1)
      ! Success as of kjw, 07/24/10



      ! Return to calling program
      END SUBROUTINE READ_TES_CH4_OBS
!------------------------------------------------------------------------------


      SUBROUTINE CALC_TES_CH4_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_TES_CH4_FORCE calculates the adjoint forcing from the TES
!  CH4 observations and updates the cost function. (dkh, 02/15/09)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (TYPE (XPLEX)) : Cost funciton                        [unitless]
!     
!     
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!  (1 ) Add more diagnostics.  Now read and write doubled CH4 (dkh, 11/08/09) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AD
      USE DAO_MOD,            ONLY : AIRDEN
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE DAO_MOD,            ONLY : TROPP
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS, GET_YEAR
      USE TIME_MOD,           ONLY : GET_MONTH, GET_DAY, GET_HOUR
      USE TIME_MOD,           ONLY : GET_TS_CHEM, EXPAND_DATE
      USE TRACER_MOD,         ONLY : XNUMOLAIR, XNUMOL
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
!kjw
      USE TRACER_MOD,         ONLY : STT
      USE FILE_MOD,           ONLY : IU_FILE
      USE TIME_MOD,           ONLY : GET_TAUe, GET_TAU
      USE LOGICAL_ADJ_MOD,    ONLY : LTES_PSO
!kjw


#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER                     :: NTSTART, NTSTOP, NT 
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL,     LTES
      INTEGER                     :: JLOOP
      TYPE (XPLEX)                      :: GC_PRES(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4(MAXLEV)
      TYPE (XPLEX)                      :: GC_PSURF
      TYPE (XPLEX)                      :: MAP(LLPAR,MAXLEV)
      TYPE (XPLEX)                      :: CH4_HAT(MAXLEV)
      TYPE (XPLEX)                      :: CH4_HAT_EXP(MAXLEV)
      TYPE (XPLEX)                      :: CH4_PERT(MAXLEV)
      TYPE (XPLEX)                      :: FORCE
      TYPE (XPLEX)                      :: DIFF
      TYPE (XPLEX)                      :: NEW_COST(MAXTES)
      TYPE (XPLEX)                      :: OLD_COST
      TYPE (XPLEX), SAVE                :: TIME_FRAC(MAXTES)
      INTEGER,SAVE                :: NTES
      TYPE (XPLEX)                      :: DOFS
      CHARACTER                   :: F117_STATUS

      !kjw for testing adjoint of tes obs operator
      TYPE (XPLEX)                      :: ADJ(LLPAR)
      TYPE (XPLEX)                      :: ADJ_SAVE(LLPAR)
      TYPE (XPLEX)                      :: PERT(LLPAR)
      TYPE (XPLEX)                      :: FD_CEN(LLPAR)
      TYPE (XPLEX)                      :: FD_POS(LLPAR)
      TYPE (XPLEX)                      :: FD_NEG(LLPAR)
      TYPE (XPLEX)                      :: COST_FUNC_0
      TYPE (XPLEX)                      :: COST_FUNC_1
      TYPE (XPLEX)                      :: COST_FUNC_2
      LOGICAL                     :: ori
      !kjw for testing adjoint of tes obs operator

      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: CH4_HAT_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: CH4_HAT_EXP_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: CH4_PERT_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: GC_CH4_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: DIFF_ADJ
      TYPE (XPLEX)                      :: S_obs_inv
      TYPE (XPLEX)                      :: GC_avg
      TYPE (XPLEX)                      :: OBS_avg
      TYPE (XPLEX)                      :: GC_avg_ADJ
!      TYPE (XPLEX)                      :: Pres_ln(MAXLEV)
!      TYPE (XPLEX)                      :: Pedges_ln(MAXLEV)
!      TYPE (XPLEX)                      :: Pedges(MAXLEV)
!      TYPE (XPLEX)                      :: Pdiff(MAXLEV)
!      TYPE (XPLEX)                      :: Nmolec(MAXLEV)
!      TYPE (XPLEX)                      :: DIFF_onTES(MAXLEV)
!      TYPE (XPLEX)                      :: Totmolec
      TYPE (XPLEX)                      :: OBS_RTVMR, GC_RTVMR
      TYPE (XPLEX)                      :: OBS_RTVMR_ADJ, GC_RTVMR_ADJ
      TYPE (XPLEX)                      :: M_STAR(4,MAXLEV)
      INTEGER                     :: reg
      TYPE (XPLEX)                      :: JJ_this(9)
      TYPE (XPLEX)                      :: Jforce_this(9)
      TYPE (XPLEX)                      :: Jdiff_this(9)

      LOGICAL, SAVE               :: FIRST = .TRUE. 
      LOGICAL, SAVE               :: VERYFIRST = .TRUE.
      LOGICAL, SAVE               :: GOD = .FALSE.
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME



      !=================================================================
      ! CALC_TES_CH4_FORCE begins here!
      !=================================================================

      print*, '     - CALC_TES_CH4_FORCE '


      !------------------GOD RETRIEVAL-------------------
      !  The god retrieval assumes perfect observations in every grid box
      !  at all time steps.  The retrieval should perfectly reproduce
      !  the "true" emission field in identical twin tests.  If it does 
      !  not, then there is a bug in the code. (kjw, 10/07/10)
      !
      !  OUTLINE
      !    A. Get STT in [ppb]
      !    B. Get pseudo-obs in [ppb]
      !    C. Calculate forcing and adjoint variable
      !


      ! If GOD == .TRUE.  THEN do GOD retrieval
      IF ( GOD == .TRUE. ) THEN

         ! Save a value of the cost function
         OLD_COST = COST_FUNC


         ! A. Get STT in [ppb]

         ! B. Get pseudo-obs in [ppb]

         ! C. Calculate forcing and adjoint variable



         ! Return to calling program
         RETURN

      ENDIF  ! End if GOD == .TRUE.

      !--------------------------------------------------


      !kjw for testing
      ori=.TRUE.
      !kjw for testing

      ! Reset 
      NEW_COST = 0D0 
      JJ_this(:)     = 0d0
      Jforce_this(:) = 0d0
      Jdiff_this(:)  = 0d0

      ! Calculate TES vs. GEOS-Chem bais
      IF ( VERYFIRST .AND. (LTES_PSO == .FALSE.) ) CALL CALC_TES_GC_BIAS
      VERYFIRST = .FALSE.

      ! Open files for diagnostic output
      IF ( FIRST ) THEN

         ! Open files for diagnostic output
         FILENAME = 'pres.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
     
         FILENAME = 'gc_ch4.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 102,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
     
         FILENAME = 'tes_ch4.NN.m'
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

         FILENAME = 'ch4_pert_adj.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 108,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_ch4_adj.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 109,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'exp_ch4_hat.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 110,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_press.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 111,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_ch4_native.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 112,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_ch4_on_tes.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 113,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_ch4_on_tes_woStrat.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 115,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_ch4_native_adj.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 114,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         !kjw for testing adjoint of obs operator
         FILENAME = 'test_adjoint_obs.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 116,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      ENDIF

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Check if it is the last hour of a day 
      !kjw. Change this for methane or else we'll never read a file.
      !kjw. For forcing every 1 hour, the following line should work:
      ! IF ( GET_NHMS() == 230000 ) THEN
      IF ( GET_NHMS() == 230000 ) THEN 
 
         ! Read the TES CH4 file for this day 
         CALL READ_TES_CH4_OBS( GET_NYMD(), NTES )

         ! If NTES = 0, it means there are no observations today.
         ! Return to calling procedure
         IF ( NTES == 0 ) THEN
            WRITE(6,*) '   No TES CH4 obs today. Returning 01 ... '
            RETURN
         ENDIF

         ! TIME is YYYYMMDD.frac-of-day. Subtract date and save just time frac
         TIME_FRAC(1:NTES) = TES(1:NTES)%TIME(1) - GET_NYMD()

      ENDIF 

      ! If NTES = 0, it means there are no more observations today.
      ! Return to calling procedure
      IF ( NTES == 0 ) THEN
         WRITE(6,*) '   No TES CH4 obs today. Returning 02 ... '
         RETURN
      ENDIF


      ! Get the range of TES retrievals for the current hour
      CALL GET_NT_RANGE( NTES, GET_NHMS(), TIME_FRAC, NTSTART, NTSTOP ) 

      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 

         print*, ' No matching TES CH4 obs for this hour'
        RETURN
      ENDIF 

      print*, ' for hour range: ', GET_NHMS(), TIME_FRAC(NTSTART), 
     &       TIME_FRAC(NTSTOP)
      print*, ' found record range: ', NTSTART, NTSTOP 


      ! Calculate S_obs_inv for stddev of ERR_PPB [ppb]
      ! Expected difference = ln( 1800 +/- ERR_PPB ) - ln( 1800 )
      !                     = ln( ( 1800 +/- ERR_PPB ) / 1800 )
      !S_obs_inv = 1. / ( LOG( ( 1800. + ERR_PPB ) / 1800. ) )**2
      S_obs_inv = 1. / ( (ERR_PPB*1d-9) ** 2 )

      print*,'kjw debug: calculate S_obs_inv.'
      print*, 'ERR_PPB = ', ERR_PPB
      !print*,'S_obs_inv (should be ~2070) = ',S_obs_inv
      print*,'S_obs_inv (should be ~6.25e14) = ',S_obs_inv



      ! Open file for this hour's satellite diagnostics
      FILENAME = 'diag_sat.YYYYMMDD.hhmm.NN'
      CALL EXPAND_NAME( FILENAME, N_CALC )
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
      OPEN( IU_FILE,      FILE=TRIM( FILENAME ),  STATUS='UNKNOWN',
     &       IOSTAT=IOS,  FORM='FORMATTED',       ACCESS='SEQUENTIAL' )



! need to update this in order to do i/o with this loop parallel 
!!      ! Now do a parallel loop for analyzing data 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, MAP, LTES, IIJJ,  I, J,  L,   LL, JLOOP )
!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_CH4,  DIFF  )
!!$OMP+PRIVATE( GC_CH4_NATIVE, CH4_PERT,   CH4_HAT, FORCE )
!!$OMP+PRIVATE( GC_CH4_NATIVE_ADJ,        GC_CH4_ADJ     )
!!$OMP+PRIVATE( CH4_PERT_ADJ,             CH4_HAT_ADJ    )
!!$OMP+PRIVATE( DIFF_ADJ                                )
      DO NT  = NTSTART, NTSTOP, -1

         !IF ( NT .EQ. 600 ) THEN
         print*, '     - CALC_TES_CH4_FORCE: analyzing record ', NT 

         ! For safety, initialize these up to LLTES 
         GC_CH4(:)       = 0d0 
         MAP(:,:)       = 0d0 
         CH4_HAT_ADJ(:)  = 0d0 
         FORCE          = 0d0 


         ! Copy LTES to make coding a bit cleaner
         LTES = TES(NT)%LTES(1)

         ! Get grid box of current record
         IIJJ  = GET_IJ(XPLX(TES(NT)%LON(1)),XPLX(TES(NT)%LAT(1)))
         I     = IIJJ(1)
         J     = IIJJ(2)

!         ! dkh debug  
!         print*, 'I,J = ', I, J
!         print*,TES(NT)%TIME(1)
!         print*,TES(NT)%LAT(1)
!         print*,TES(NT)%LON(1)

         ! Get GC pressure levels (mbar) 
         DO L = 1, LLPAR
            GC_PRES(L) = GET_PCENTER(I,J,L)
         ENDDO

         ! Get GC surface pressure (mbar) 
         GC_PSURF = GET_PEDGE(I,J,1) 


         ! Calculate the interpolation weight matrix 
         MAP(1:LLPAR,1:LTES) 
     &      = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
     &                    LTES,  TES(NT)%PRES(1:LTES), TES(NT)%PRES(1) )


         ! Get CH4 values at native model resolution
         !DO L = 1, LLPAR


            !kjw. JLOP and CSPEC are variables in comode_mod.f, associated
            !     with smvgear. For getting CH4 values at native model
            !     resolution, write my own thing, it'll probably have to do
            !     with using CH4 from STT when below tropopause, and using
            !     the TES retrieval above the tropopause. Therefore, there
            !     will be no adjoint forcing above the tropopause.
            !kjw.

               !kjw
               ! Get CH4 from restored STT array
               ! Find out units of STT (I think it's [kg/box]).
               !     If so, convert to [v/v] before next step
               GC_CH4_NATIVE(:) = CHK_STT(I,J,:,1)
!               print*,'CHK_STT(14,14,14,1) = ',CHK_STT(14,14,14,1)
!               print*,'    STT(14,14,14,1) = ',STT(14,14,14,1)

         ! Unit Conversion
         DO L=1,LLPAR

            ! Convert from [kg/box] --> [v/v]
            !    Numerator   = moles CH4/box
            !    Denominator = moles air/box
            GC_CH4_NATIVE(L) = (GC_CH4_NATIVE(L)*XNUMOL(1)/6.022d23 ) / 
     &                      ( AD(I,J,L) * XNUMOLAIR / 6.022d23 )
            !!!! Bypass unit conversion
            !!!GC_CH4_NATIVE(:) = GC_CH4_NATIVE(:)
         ENDDO


         ! Interpolate GC CH4 column to TES grid 
         DO LL = 1, LTES
            GC_CH4(LL) = 0d0 
            DO L = 1, LLPAR 
               GC_CH4(LL) = GC_CH4(LL) 
     &                    + MAP(L,LL) * GC_CH4_NATIVE(L) 
            ENDDO
         ENDDO

         IF ( NT == 600 ) THEN
            !print*,'LTES = ',LTES
            !print*,'LLPAR = ',LLPAR
            !print*,'GC_PSURF = ',GC_PSURF

            !WRITE(6,'(a)') 'GEOS-Chem pressure grid'
            !WRITE(6,'(F10.3)') ( GC_PRES(L), L=1,47 )

            !WRITE(6,'(a)') 'TES pressure grid'
            !WRITE(6,'(F10.3)') ( TES(NT)%PRES(L), L=1,65 )

            !WRITE(6,'(a)') '20th row of MAP matrix in CALC_FORCE'
            !WRITE(6,'(5F8.5)') ( MAP(20:24,L) , L=1,65 )
            !print*,'GC_CH4 (observation) = ', TES(NT)%GC_CH4(10)
            !print*,'CH4 (model) = ', GC_CH4(10)


         ENDIF

!         IF ( NT == 600 ) THEN
!         ! dkh debug: compare profiles:
!         print*, ' GC_PRES, GC_native_CH4 [ppb] '
!         WRITE(6,100) (GC_PRES(L), GC_CH4_NATIVE(L)*1d9, 
!     &                   L = LLPAR, 1, -1 )
!         print*, ' TES_PRES, GC_CH4  ' 
!         WRITE(6,100) (TES(NT)%PRES(LL), 
!     &                 GC_CH4(LL)*1d9, LL = LTES, 1, -1 ) 
!         ENDIF
! 100  FORMAT(1X,F16.8,1X,F16.8)


         !--------------------------------------------------------------
         ! Apply TES observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by TES [lnvmr]
         !    x_a   = TES apriori column               [lnvmr]
         !    x_m   = GC modeled column                [lnvmr]
         !    A_k   = TES averaging kernel 
         !--------------------------------------------------------------

         ! x_m - x_a
         DO L = 1, LTES 
           GC_CH4(L)   = MAX(GC_CH4(L), 1d-10)
           CH4_PERT(L) = LOG(GC_CH4(L)) - LOG(TES(NT)%PRIOR(L))
         ENDDO

!         !kjw testing
!         IF ( ori .EQ. .TRUE. ) THEN
!            print*,'CH4_PERT(13) = ',CH4_PERT(13)
!         ENDIF
!         !kjw testing
!!!!!!         ! Bypass !x_m - x_a
!!!!!!         CH4_PERT(:) = GC_CH4(:)
!!!!!!
     
         ! x_a + A_k * ( x_m - x_a )  
         !  AVG_KERNEL indexing may look backwards because BPCH files storing AK 
         !    values use IDL column major indexing.
         DO L = 1, LTES
            CH4_HAT(L)    = 0d0 
            DO LL = 1, LTES
               CH4_HAT(L) = CH4_HAT(L) 
     &                    + TES(NT)%AVG_KERNEL(LL,L) * CH4_PERT(LL) 
            ENDDO
            CH4_HAT(L)    = CH4_HAT(L) + LOG(TES(NT)%PRIOR(L))
         ENDDO
         ! Indexing of Averaging Kernel is seemingly backwards because
         ! TES observation files processed in IDL, which is column major
 
!!!         ! Bypass !x_a + A_k * ( x_m - x_a )  
!!!         CH4_HAT(:) = CH4_PERT(:)
!!!

!         !kjw testing
!         IF ( ori .EQ. .TRUE. ) THEN
!            print*,'CH4_HAT(13) = ',CH4_HAT(13)
!         ENDIF
!         !kjw testing

!         !--------------------------------------------------------------
!         ! Calculate column average ln(vmr) weighted by # density of TES grid
!         !    This operation is self adjoint.
!         ! To get # density [molec / m^2], get pressure differences in grid.
!         !    Get kg/m^2 of air in each box. (F=ma) dP=kg*g
!         !    Get molec/m^2 air  molec/m2 = kg/m2 * XNUMOLAIR
!         ! TES pressure grid linear in ln(pres).  Get pressure at edge of boxes
!         !--------------------------------------------------------------
!         pres_ln(1:LTES)=LOG(TES(NT)%pres(1:LTES))   ! [hPa] --> ln([hPa])
!         Pedges_ln(1)=pres_ln(1)             ! Bottom edge is surface pressure
!         DO L=2,LTES-1
!            Pedges_ln(L) = ( pres_ln(L) + pres_ln(L+1) ) / 2.
!         ENDDO
!         Pedges=EXP(pedges_ln)               ! ln([hPa]) --> [hPa]
!         Pedges(LTES)=0                      ! Top of atmosphere
!         !print*,' kjw debug:  TES pressure edges'
!         !print*,         Pedges
!
!         ! Calculate pressure difference of each LTES-1 boxes
!         DO L=1,LTES-1
!            Pdiff(L) = Pedges(L) - Pedges(L+1)
!         ENDDO
!
!         ! Calculate # molecules air in each LTES-1 obxes
!         Pdiff(:) = 100 * Pdiff(:)              ! [hPa] --> [Pa]
!         Pdiff(:) = Pdiff(:) / 9.8              ! [hPa] --> [kg]
!         Nmolec(:)= Pdiff(:) * XNUMOLAIR        ! [kg]  --> [molec]
!         Totmolec = 0d0
!         DO L=1,LTES-1
!            IF ( TROPP(I,J) < Pedges(L) ) THEN
!               Totmolec = Totmolec + Nmolec(L)
!            ENDIF
!         ENDDO
!
!         IF ( NT .EQ. 600 ) THEN
!            print*,' kjw debug:# molecules in column (should be ~2e29)'
!            print*,         Totmolec
!            print*,' kjw debug: Nmolec'
!            DO L=1,LTES
!            IF ( TROPP(I,J) < Pedges(L) ) THEN
!               print*, Nmolec(L)/Totmolec
!            ENDIF
!            ENDDO
!            print*,'SUM(Nmolec(1:22)/Totmolec)',
!     &                        SUM(Nmolec(1:22)/Totmolec)
!            print*,'SUM(NMOLEC) = ',SUM(Nmolec)
!         ENDIF
!
!         ! Calculate column average ln(vmr) weighted by # density of TES levels
!         ! Only include levels with tropospheric air
!         OBS_avg = 0d0
!         GC_avg  = 0d0
!         DIFF_onTES(:) = 0d0
!         DIFF_onTES(1:LTES) = LOG(TES(NT)%GC_CH4(1:LTES)) - 
!     &                                 CH4_HAT(1:LTES)
!         DO L=1,LTES-1
!            IF ( TROPP(I,J) < Pedges(L) ) THEN
!               OBS_avg = OBS_avg + 
!     &                     LOG(TES(NT)%GC_CH4(L+1)) * Nmolec(L)/Totmolec
!               GC_avg  = GC_avg  + 
!     &                     CH4_HAT(L+1) * Nmolec(L)/Totmolec
!            ENDIF
!         ENDDO

         ! Transform from [ln(vmr)] --> [vmr]
         CH4_HAT_EXP = EXP(CH4_HAT)

         ! Calculate RTVMR for profiles.
         CALL GET_RTVMR( NT, TES(NT)%CH4, OBS_RTVMR, M_STAR )
         CALL GET_RTVMR( NT, CH4_HAT_EXP, GC_RTVMR,  M_STAR )


         ! kjw debug.  Check RTVMR stuff
         !IF ( NT == 600 ) THEN
         !   print*,'Check RTVMR stuff'
         !   print*,'Lat, Lon, PSURF of observation #600 '
         !   print*,TES(NT)%LAT(1),TES(NT)%LON(1),TES(NT)%PRES(1)
         !   print*,'GC_RTVMR = ',GC_RTVMR
         !   print*,'CH4_HAT =  ',CH4_HAT
         !ENDIF



         ! Retrieve RTVMR from TES structure for pseudo-observations
         ! TES(NT)%GC_CH4(67) set during SUBROUTINE MAKE_PSEUDO_OBS
         IF ( LTES_PSO ) THEN
            OBS_RTVMR = 0d0
            OBS_RTVMR = TES(NT)%GC_CH4(67)
         ENDIF

         ! DIFF = model - obs.     units: [ln(vmr)] / m^2
         !IF ( NT == 600 ) THEN
         !   print*,'Error difference [ppb]:  ',TES(NT)%ERR(1)*1d9
         !   print*,'GC_RTVMR = ',GC_RTVMR*1d9
         !   print*,'OBS_RTVMR = ',OBS_RTVMR*1d9
!        !    WRITE(6, '(a)') 'Pseudo-obs         GEOS-Chem '
!        !    WRITE(6, 545)   ( TES(NT)%GC_CH4(L), EXP(CH4_HAT(L)), 
!     &  !                                             L=1,LTES )
         !ENDIF
 545     FORMAT(F16.14,2x,F16.14)


         ! Calculate DOFS for satellite diagnostic file
         DOFS = 0d0
         DO L=1,LTES
            DOFS = DOFS + TES(NT)%AVG_KERNEL(L,L)
         ENDDO

         !--------------------------------------------------------------
         ! Calculate cost function, given S is error on ln(vmr)
         ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------

         ! If using pseudo-obs, do not apply bias
         DIFF = GC_RTVMR - OBS_RTVMR! + BIAS_PPB * 1d-9


         ! Calculate DIFF^T * S_{obs}^{-1} * DIFF
         FORCE = 0d0
         FORCE = 2 * DIFF * S_obs_inv
         NEW_COST(NT) = 0.5d0 * DIFF * FORCE




         ! Write satellite information to file
         ! I,J,LAT,LON,TIME,HOUR,model RTVMR,obs RTVMR,DOFS
         IF ( NT == NTSTART ) THEN
            WRITE(IU_FILE,301) 'I','J','LAT','LON','MONTH','DAY','HOUR',
     &                     'TIME_FRAC','MODEL_RTVMR','OBS_RTVMR','DOFS',
     &                     'ERROR'
         ENDIF
         WRITE(IU_FILE,302) I,J,TES(NT)%LAT(1),TES(NT)%LON(1),
     &                      GET_MONTH(), GET_DAY(), GET_HOUR(),
     &                      TIME_FRAC(NT), 
     &                      1e9*GC_RTVMR, 1e9*OBS_RTVMR, DOFS,
     &                      1e9*TES(NT)%ERR(1)

 301     FORMAT(A4,2x,A4,2x,A8,  2x,A8,  2x,A6,2x,A4,2x,A4,2x,A16,   2x,
     &          A12,  2x,A12,  2x,A7,2x,A8)
 302     FORMAT(I4,2x,I4,2x,F8.3,2x,F8.3,2x,I4,2x,I4,2x,I4,2x,F16.13,2x, 
     &          F12.4,2x,F12.4,2x,F7.3,2x,F8.3)


!         ! Calculate difference between modeled and observed profile
!         ! Eliminate stratospheric forcing in this 
!         DO L = 1, LTES 
!            IF ( TROPP(I,J) < TES(NT)%PRES(L) ) THEN 
!               DIFF(L) = CH4_HAT(L) - LOG( TES(NT)%GC_CH4(L) )
!            ELSE 
!               DIFF(L) = 0d0
!            ENDIF 
!         ENDDO
!
!
!         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF 
!         DO L = 1, LTES
!            FORCE(L)     = 0d0 
!            !FORCE(L)  = FORCE(L) + TES(NT)%S_OER_INV(L,L) * DIFF(L)
!            DO LL = 1, LTES
!               FORCE(L)  = FORCE(L) + TES(NT)%S_OER_INV(L,LL) * DIFF(LL)
!            ENDDO
!            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
!         ENDDO

         !IF ( NT == 600 ) THEN
         !   print*,'DIFF = ', DIFF
         !   print*,'FORCE = ',FORCE
         !   print*,'NEW_COST = ',NEW_COST(NT)
         !ENDIF
!         ! dkh debug: compare profiles:
!         print*, ' CH4_HAT, CH4_TES, CH4_GC [ppb]'
!         WRITE(6,090) ( 1d9 * EXP(CH4_HAT(L)),
!     &                  1d9 * TES(NT)%CH4(L),
!     &                  1d9 * TES(NT)%GC_CH4(L),
!     &                  L,    L = LTES, 1, -1   )
!
!         print*, ' TES_PRIOR, CH4_HAT, CH4_GC  [lnvmr], diag(S^-1)'
!         WRITE(6,101) ( LOG(TES(NT)%PRIOR(L)), CH4_HAT(L), 
!     &                  LOG(TES(NT)%GC_CH4(L)), TES(NT)%S_OER_INV(L,L),
!     &                  L,  L = LTES, 1, -1 )
!         ENDIF
! 090  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,i3)
! 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,d14.6,1x,i3)

         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

      !kjw. We've now calculated:
      !        1) forcing ( 2 * S_{obs}^{-1} * DIFF = FORCE ) units of [lnvmr]^{-1}
      !        2) New contribution to cost function do to diff
      !             DIFF * S_{obs}^{-1} * DIFF
      !     This has all been done on the TES pressure grid
      !
      !     At this point, we need to initialize the adjoint variable: STT_ADJ
      !     Do so by applying the adjoint of all operators used to get
      !               STT --> ln(vmr) for calculating ( F(x)-y )

!         IF ( NT == 600 ) THEN
!         ! dkh debug
!         print*, 'DIFF , FORCE ' 
!         WRITE(6,102) (DIFF(L), FORCE(L), 
!     &       L = LTES, 1, -1 )
!         ENDIF
! 102  FORMAT(1X,d14.6,1X,d14.6)

         ! The adjoint forcing is  2 * S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE

         ! Adjoint of difference
         GC_RTVMR_ADJ = DIFF_ADJ


         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !ADJ_DIFF(:) = 1d0
         !NEW_COST(NT) = ?? SUM(ABS(LOG(CH4_HAT(1:LTES))))
         !print*, ' sumlog =', SUM(ABS(LOG(CH4_HAT(:))))
         !print*, ' sumlog =', ABS(LOG(CH4_HAT(:)))


         ! Adjoint of RTVMR Averaging
         DO L=1,LTES
            CH4_HAT_ADJ(L) = M_STAR(2,L) * GC_RTVMR_ADJ
         ENDDO

         ! kjw debug
         !IF ( NT == 600 ) THEN
         !   print*,'CH4_HAT_ADJ = ',CH4_HAT_ADJ
         !ENDIF

         ! Adjoint of ln(vmr) --> vmr
         DO L=1,LTES
            IF ( CH4_HAT_ADJ(L) /= 0.0 ) THEN
               CH4_HAT_EXP_ADJ(L) = CH4_HAT_ADJ(L) * CH4_HAT_EXP(L)
            ELSE
               CH4_HAT_EXP_ADJ(L) = 0d0
            ENDIF
         ENDDO

         ! kjw debug
         !IF ( NT == 600 ) THEN
         !   print*,'CH4_HAT_EXP_ADJ = ',CH4_HAT_EXP_ADJ
         !ENDIF

!         DO L = 1, LTES-1
!            IF ( TROPP(I,J) < Pedges(L) ) THEN 
!               CH4_HAT_ADJ(L+1) =  GC_RTVMR_ADJ * Nmolec(L)/Totmolec
!            ELSE
!               CH4_HAT_ADJ(L+1) = 0d0
!            ENDIF 
!         ENDDO 


!         ! Adjoint of difference
!         DO L = 1, LTES 
!            IF ( TROPP(I,J) < TES(NT)%PRES(L) ) THEN 
!               CH4_HAT_ADJ(L) =  DIFF_ADJ(L)
!            ELSE
!               CH4_HAT_ADJ(L) = 0d0
!            ENDIF 
!         ENDDO 

         ! adjoint of TES operator
         DO L  = 1, LTES
            CH4_PERT_ADJ(L)    = 0d0
            DO LL = 1, LTES
               CH4_PERT_ADJ(L) = CH4_PERT_ADJ(L) 
     &                        + TES(NT)%AVG_KERNEL(L,LL) 
     &                        * CH4_HAT_EXP_ADJ(LL)
           ENDDO
         ENDDO

         ! Adjoint of x_m - x_a (adjoint of natural log transform)
         DO L = 1, LTES 
           ! fwd code:
           !GC_CH4(L)   = MAX(GC_CH4(L), 1d-10)
           !CH4_PERT(L) = LOG(GC_CH4(L)) - LOG(TES(NT)%PRIOR(L))
           ! adj code:
           IF ( GC_CH4(L) > 1d-10 ) THEN 
              GC_CH4_ADJ(L) = 1d0 / GC_CH4(L) * CH4_PERT_ADJ(L)
           ELSE 
              GC_CH4_ADJ(L) = 1d0 / 1d-10 * CH4_PERT_ADJ(L)
           ENDIF 
         ENDDO

!         ! dkh debug
!         print*, 'CH4_HAT_ADJ, CH4_PERT_ADJ, GC_CH4_ADJ'
!         WRITE(6,103) (CH4_HAT_ADJ(L), CH4_PERT_ADJ(L), GC_CH4_ADJ(L), 
!     &       L = LTES, 1, -1 )
! 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)


         ! adjoint of interpolation 
         DO L  = 1, LLPAR
            GC_CH4_NATIVE_ADJ(L) = 0d0 
            DO LL = 1, LTES
               GC_CH4_NATIVE_ADJ(L) = GC_CH4_NATIVE_ADJ(L)
     &                             + MAP(L,LL) * GC_CH4_ADJ(LL)
            ENDDO
         ENDDO

         ! kjw
         ! Adjoint of interpolation leaves GC_CH4_NATIVE_ADJ with some zeros
         !   in the lower troposphere.  This occurs because the GC pres grid is
         !   finer in lower troposphere than the TES grid. So, when interpolating
         !   from GC --> TES, the contribution from some GC grid boxes to any TES
         !   grid box is zero.  Unfortunately, when we go back from TES to GC grid,
         !   this means that some 


!         WRITE(114,112) ( GC_CH4_NATIVE_ADJ(L),     L=LLPAR,1,-1)
!
!         ! dkh debug
!         print*, 'GC_CH4_NATIVE_ADJ 1 '
!         WRITE(6,104) (GC_CH4_NATIVE_ADJ(L), L = LLPAR, 1, -1 )

         DO L = 1, LLPAR 

            ! Adjoint of unit conversion 
            GC_CH4_NATIVE_ADJ(L) = ( GC_CH4_NATIVE_ADJ(L) * 
     &                                   XNUMOL(1) / 6.022d23 ) / 
     &                     ( AD(I,J,L) * XNUMOLAIR / 6.022d23 )


            ! Just to make sure we're only forcing the troposphere
            IF ( ITS_IN_THE_TROP(I,J,L) ) THEN 
 
               ! Pass adjoint back to adjoint tracer array
               STT_ADJ(I,J,L,1) =
     &            STT_ADJ(I,J,L,1) + GC_CH4_NATIVE_ADJ(L)

            ENDIF 

         ENDDO

         !kjw debug
         !IF ( NT == 600 ) THEN
         !   print*,'GC_CH4_NATIVE_ADJ = ',GC_CH4_NATIVE_ADJ
         !ENDIF
         !kjw debug

!!         ! dkh debug
!         print*, 'GC_CH4_NATIVE_ADJ conv '
!         WRITE(6,104) (GC_CH4_NATIVE_ADJ(L), L = LLPAR, 1, -1 )
! 104  FORMAT(1X,d14.6)
!
!
!         WRITE(101,110) ( TES(NT)%PRES(LL),        LL=LTES,1,-1)
!         WRITE(102,110) ( 1d9 * GC_CH4(LL),        LL=LTES,1,-1)
!         WRITE(103,110) ( 1d9 * TES(NT)%CH4(LL),   LL=LTES,1,-1)
!         WRITE(104,110) ( 1d9 * TES(NT)%PRIOR(LL), LL=LTES,1,-1)
!         WRITE(105,110) ( DIFF(LL),                LL=LTES,1,-1)
!         WRITE(106,112) ( FORCE(LL),               LL=LTES,1,-1)
!         WRITE(107,111) NT, LTES
!         WRITE(108,112) ( CH4_PERT_ADJ(LL),         LL=LTES,1,-1)
!         WRITE(109,112) ( GC_CH4_ADJ(LL),           LL=LTES,1,-1)
!         WRITE(110,110) ( 1d9 * EXP(CH4_HAT(LL)),   LL=LTES,1,-1)
!         WRITE(111,110) ( GC_PRES(L),              L=LLPAR,1,-1)
!         WRITE(112,110) ( 1d9 * GC_CH4_NATIVE(L),   L=LLPAR,1,-1)
!         WRITE(113,110) ( 1d9 * CH4_HAT(LL),        LL=LTES,1,-1)
!         WRITE(115,110) ( 1d9 * TES(NT)%GC_CH4(LL), LL=LTES,1,-1)
 110     FORMAT(F18.6,1X)
 111     FORMAT(i4,1X,i4,1x)
 112     FORMAT(D14.6,1X)


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
!      IF ( NT .EQ. 600 ) THEN
!      WRITE(116,210) '   L'        ,       '  TROP', '     GC_PRES',
!     &               '      FD_POS', '      FD_NEG', '      FD_CEN',
!     &               '         ADJ', '    COST_POS', '    COST_NEG', 
!     &               '  FD_POS/ADJ', '  FD_NEG/ADJ', '  FD_CEN/ADJ'
!      PERT(:) = 1D0
!      CALL CALC_TES_CH4_FORCE_FD( COST_FUNC_0, PERT, ADJ, NT )
!      ori=.FALSE.
!      ADJ_SAVE(:) = ADJ(:)
!      print*, 'dch4:  COST_FUNC_0 = ', COST_FUNC_0
!      WRITE(116,213) 'I           ', I
!      WRITE(116,213) 'J           ', J
!      WRITE(116,213) 'LTES        ',TES(NT)%LTES(1)
!      WRITE(116,212) 'GC_PSURF    ', GC_PSURF
!      WRITE(116,212) 'TES PSURF   ',TES(NT)%PRES(1)
!      WRITE(116,212) 'NEW_COST:   ',NEW_COST(NT)
!      WRITE(116,213) 'NT          ', NT
!      WRITE(116,212) 'COST_FUNC_0:',( COST_FUNC_0 )
!      WRITE(116,212) 'TES(NT).TIME',TES(NT)%TIME(1)
!      DO L = 1, 47
!         PERT(:) = 1D0
!         PERT(L) = 1.1
!         COST_FUNC_1 = 0D0
!         CALL CALC_TES_CH4_FORCE_FD( COST_FUNC_1, PERT, ADJ, NT )
!         PERT(L) = 0.9
!         COST_FUNC_2 = 0D0
!         CALL CALC_TES_CH4_FORCE_FD( COST_FUNC_2, PERT, ADJ, NT )
!         FD_CEN(L)   = ( COST_FUNC_1 - COST_FUNC_2 ) / 0.2d0
!         FD_POS(L)   = ( COST_FUNC_1 - COST_FUNC_0 ) / 0.1d0
!         FD_NEG(L)   = ( COST_FUNC_0 - COST_FUNC_2 ) / 0.1d0
!         WRITE(116, 211)  L, ITS_IN_THE_TROP(I,J,L), GC_PRES(L),
!     &                    FD_POS(L),   FD_NEG(L),
!     &                    FD_CEN(L),   ADJ_SAVE(L), 
!     &                    COST_FUNC_1, COST_FUNC_2, 
!     &                    FD_POS(L)/ADJ_SAVE(L),
!     &                    FD_NEG(L)/ADJ_SAVE(L),
!     &                    FD_CEN(L)/ADJ_SAVE(L)
!      ENDDO
!      WRITE(116,'(a)') '----------------------------------------------'
!
! 210  FORMAT(A4,2x,A6,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,
!     &       A12,2x,A12,2x,A12,2x,A12,2x)
! 211  FORMAT(I4,2x,L6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,
!     &        2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6)
! 212  FORMAT(A12,F22.6)
! 213  FORMAT(A12,I4)
! 214  FORMAT(I4,2x,F18.6,2x,F18.6)
!! -----------------------------------------------------------------------
!      ENDIF  ! IF ( NT .EQ. 600 )


      print*, ' - CALC_TES_CH4_FORCE:  NEW_COST(NT) = ',NEW_COST(NT)


!      ! kjw
!      ! Calculate contribution to cost function for regions in my 9 box grid
!      ! Here, we define regions
!      reg = 0
!      IF ( LTES_PSO == .TRUE. ) THEN
!         IF (( J-1 > JJPAR*2./3. ) .AND. ( I+1 < IIPAR/3.   ))   reg = 1
!         IF (( J-1 > JJPAR*2./3. ) .AND. ( I-1 > IIPAR*2./3.))   reg = 3
!         IF (( J-1 > JJPAR*2./3. ) .AND. ( reg == 0 )       )    reg = 2
!         IF (( J+1 < JJPAR/3.    ) .AND. ( I+1 < IIPAR/3.   ))   reg = 7
!         IF (( J+1 < JJPAR/3.    ) .AND. ( I-1 > IIPAR*2./3.))   reg = 9
!         IF (( J+1 < JJPAR/3.    ) .AND. ( reg == 0 )       )    reg = 8
!         IF (( I+1 < IIPAR/3.    ) .AND. ( reg == 0 )       )    reg = 4
!         IF (( I-1 > IIPAR*2./3. ) .AND. ( reg == 0 )       )    reg = 6
!         IF (( reg == 0 )                                   )    reg = 5
!      ENDIF   ! ENDIF LTES_PSO == .TRUE.
!
!      ! Assign value to proper region
!      JJ_this(reg)     = JJ_this(reg)     + NEW_COST(NT)
!      Jforce_this(reg) = Jforce_this(reg) + FORCE
!      Jdiff_this(reg)  = Jdiff_this(reg)  + DIFF
!
!      JJ(reg)     = JJ(reg)     + NEW_COST(NT)
!      Jforce(reg) = Jforce(reg) + FORCE
!      Jdiff(reg)  = Jdiff(reg)  + DIFF


      ENDDO  ! NT
!!$OMP END PARALLEL DO

      print*, '     - CALC_TES_CH4_FORCE:  finished assimilating ' //
     &               'data this hour.'

      print*,'NEW_COST(NTSTOP:NTSTART) = ', NEW_COST(NTSTOP:NTSTART)
      print*,'SUM(NEW_COST(NTSTOP:NTSTART)) = ', 
     &           SUM(NEW_COST(NTSTOP:NTSTART))
      print*,'NTSTART, NTSTOP = ', NTSTART, NTSTOP

      ! Update cost function 
      COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))

      IF ( FIRST ) FIRST = .FALSE. 


!      ! Print information about JJ_this, Jforce_this, Jdiff_this
!      ! Cost function info by region
!      print*,'EYOEYOEYO, J info below'
!      print*,'Year, month, day, hour, hours before end of simulation'
!      print*,GET_YEAR(), GET_MONTH(), GET_DAY(), GET_HOUR(),
!     &           GET_TAUe() - GET_TAU()
!      WRITE(6,820) 'JJ_this=         ', JJ
!      WRITE(6,820) 'Jforce_this=     ', Jforce
!      WRITE(6,820) 'Jdiff*1d9_this=  ', Jdiff*1d9
! 820  FORMAT(A18, 2x, 9F28.6 )


      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, ' TES contribution this hour = ', COST_FUNC - OLD_COST  


      ! Close Satellite diagnostic file
      CLOSE( IU_FILE )


      ! kjw
      ! Print Information about cost function and forcing according
      ! to region defined by apriori. Print information added during 
      ! this call to CALC_TES_CH4_FORCE
      ! kjw



      ! Return to calling program
      END SUBROUTINE CALC_TES_CH4_FORCE

!------------------------------------------------------------------------------

      SUBROUTINE CALC_TES_CH4_FORCE_FD( COST_FUNC_A, PERT, ADJ, NT )
!
!******************************************************************************
!  Subroutine CALC_TES_CH4_FORCE_FD tests the adjoint of CALC_TES_CH4_FORCE
!  (dkh, 05/05/10) 
!    
!    Can be driven with:
!      PERT(:) = 1D0
!      CALL CALC_TES_CH4_FORCE_FD( COST_FUNC_0, PERT, ADJ )
!      ADJ_SAVE(:) = ADJ(:)
!      print*, 'do3:  COST_FUNC_0 = ', COST_FUNC_0
!      DO L = 1, 30
!         PERT(:) = 1D0
!         PERT(L) = 1.1
!         COST_FUNC = 0D0
!         CALL CALC_TES_CH4_FORCE_FD( COST_FUNC_1, PERT, ADJ )
!         PERT(L) = 0.9
!         COST_FUNC = 0D0
!         CALL CALC_TES_CH4_FORCE_FD( COST_FUNC_2, PERT, ADJ )
!         FD(L)       = ( COST_FUNC_1 - COST_FUNC_2 ) / 0.2d0
!         print*, 'do3:  FD  = ', FD(L), L
!         print*, 'do3:  ADJ = ', ADJ_SAVE(L), L
!         print*, 'do3:  COST = ', COST_FUNC, L
!         print*, 'do3:  FD / ADJ ', FD(L) / ADJ_SAVE(L) , L
!      ENDDO
!
! 
! 
!
!  Arguments as Input/Output:
!  ===========================================================================
!  (1 ) COST_FUNC_A (TYPE (XPLEX)) : Cost funciton                        [unitless]
!     
!     
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!  (1 ) Add more diagnostics.  Now read and write doubled CH4 (dkh, 11/08/09) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
!      USE ADJ_ARRAYS_MOD,     ONLY : CH4_PROF_SAV
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE DAO_MOD,            ONLY : AD
      USE DAO_MOD,            ONLY : AIRDEN
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE DAO_MOD,            ONLY : TROPP
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACER_MOD,         ONLY : XNUMOLAIR, XNUMOL
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
      USE TRACER_MOD,         ONLY : STT
      USE LOGICAL_ADJ_MOD,    ONLY : LTES_PSO


#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC_A
      INTEGER, INTENT(IN)         :: NT
      TYPE (XPLEX), INTENT(IN)          :: PERT(LLPAR) 
      TYPE (XPLEX), INTENT(OUT)         :: ADJ(LLPAR)

      ! Local variables 
      INTEGER                     :: NTSTART, NTSTOP
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL,     LTES
      INTEGER                     :: JLOOP
      TYPE (XPLEX)                      :: GC_PRES(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4(MAXLEV)
      TYPE (XPLEX)                      :: GC_PSURF
      TYPE (XPLEX)                      :: MAP(LLPAR,MAXLEV)
      TYPE (XPLEX)                      :: CH4_HAT(MAXLEV)
      TYPE (XPLEX)                      :: CH4_HAT_EXP(MAXLEV)
      TYPE (XPLEX)                      :: CH4_PERT(MAXLEV)
      TYPE (XPLEX)                      :: FORCE
      TYPE (XPLEX)                      :: DIFF
      TYPE (XPLEX)                      :: NEW_COST!(MAXTES)
      TYPE (XPLEX)                      :: OLD_COST
      !TYPE (XPLEX), SAVE                :: TIME_FRAC!(MAXTES)
      !INTEGER,SAVE                :: NTES 

      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: CH4_HAT_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: CH4_HAT_EXP_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: CH4_PERT_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: GC_CH4_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: DIFF_ADJ
      TYPE (XPLEX)                      :: S_obs_inv
      TYPE (XPLEX)                      :: GC_RTVMR
      TYPE (XPLEX)                      :: OBS_RTVMR
      TYPE (XPLEX)                      :: GC_RTVMR_ADJ
      TYPE (XPLEX)                      :: M_STAR(4,MAXLEV)
      TYPE (XPLEX)                      :: Pres_ln(MAXLEV)
      TYPE (XPLEX)                      :: Pedges_ln(MAXLEV)
      TYPE (XPLEX)                      :: Pedges(MAXLEV)
      TYPE (XPLEX)                      :: Pdiff(MAXLEV)
      TYPE (XPLEX)                      :: Nmolec(MAXLEV)
      TYPE (XPLEX)                      :: Totmolec
      TYPE (XPLEX)                      :: DIFF_onTES(MAXLEV)
   
      !LOGICAL, SAVE               :: FIRST = .TRUE. 
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME


      !=================================================================
      ! CALC_TES_CH4_FORCE_FD begins here!
      !=================================================================

      !print*, '     - CALC_TES_CH4_FORCE_FD '
  
      NEW_COST = 0D0

      ! Calculate S_obs_inv
      !S_obs_inv = 1. / ( LOG( ( 1800. + ERR_PPB ) / 1800. ) )**2
      S_obs_inv = 1. / ( 40.0d-9 ** 2 )


!      ! Open files for output
!      IF ( FIRST ) THEN
!         FILENAME = 'force_adj_stuff.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!     
!
!      ENDIF
!
!      ! Save a value of the cost function first
!      OLD_COST = COST_FUNC_A
!
!      ! Check if it is the last hour of a day 
!      IF ( GET_NHMS() == 230000 ) THEN
! 
!         ! Read the TES CH4 file for this day 
!         CALL READ_TES_CH4_OBS( GET_NYMD(), NTES )
!
!         ! If NTES = 0, it means there are no observations today.
!         ! Return to calling procedure
!         IF ( NTES == 0 ) THEN
!            WRITE(6,*) '   No TES CH4 obs today. Returning 01 ... '
!            RETURN
!         ENDIF
!
!         ! TIME is YYYYMMDD.frac-of-day.  Subtract date and save just time fraction
!         TIME_FRAC(1:NTES) = TES(1:NTES)%TIME(1) - GET_NYMD()
!
!!      ENDIF 
!
!      ! If NTES = 0, it means there are no observations today.
!      ! Return to calling procedure
!      IF ( NTES == 0 ) THEN
!         WRITE(6,*) '   No TES CH4 obs today. Returning 02 ... '
!         RETURN
!      ENDIF
!
!
!      ! Get the range of TES retrievals for the current hour
!      CALL GET_NT_RANGE( NTES, GET_NHMS(), TIME_FRAC, NTSTART, NTSTOP ) 
!
!      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 
!
!         print*, ' No matching TES CH4 obs for this hour'
!        RETURN
!      ENDIF 
!
!      print*, ' for hour range: ', GET_NHMS(), TIME_FRAC(NTSTART), 
!     &       TIME_FRAC(NTSTOP)
!      print*, ' found record range: ', NTSTART, NTSTOP 

!! need to update this in order to do i/o with this loop parallel 
!!!      ! Now do a parallel loop for analyzing data 
!!!$OMP PARALLEL DO
!!!$OMP+DEFAULT( SHARED )
!!!$OMP+PRIVATE( NT, MAP, LTES, IIJJ,  I, J,  L,   LL, JLOOP )
!!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_CH4,  DIFF  )
!!!$OMP+PRIVATE( GC_CH4_NATIVE, CH4_PERT,   CH4_HAT, FORCE )
!!!$OMP+PRIVATE( GC_CH4_NATIVE_ADJ,        GC_CH4_ADJ     )
!!!$OMP+PRIVATE( CH4_PERT_ADJ,             CH4_HAT_ADJ    )
!!!$OMP+PRIVATE( DIFF_ADJ                                )
!      DO NT = NTSTART, NTSTOP, -1 

         !print*, '     - CALC_TES_CH4_FORCE_FD: analyzing record ', NT 

         ! For safety, initialize these up to LLTES 
         GC_CH4(:)       = 0d0 
         MAP(:,:)       = 0d0 
         CH4_HAT_ADJ(:)  = 0d0 
         FORCE       = 0d0 


         ! Copy LTES to make coding a bit cleaner
         LTES = TES(NT)%LTES(1)

         ! Get grid box of current record
         IIJJ  = GET_IJ(XPLX(TES(NT)%LON(1)),XPLX(TES(NT)%LAT(1)))
         I     = IIJJ(1)
         J     = IIJJ(2)

         !print*, 'I,J = ', I, J

         ! Get GC pressure levels (mbar) 
         DO L = 1, LLPAR
            GC_PRES(L) = GET_PCENTER(I,J,L)
         ENDDO

         ! Get GC surface pressure (mbar) 
         GC_PSURF = GET_PEDGE(I,J,1) 


         ! Calculate the interpolation weight matrix 
         MAP(1:LLPAR,1:LTES) 
     &      = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
     &                    LTES,  TES(NT)%PRES(1:LTES), TES(NT)%PRES(1) )


         ! Get CH4 values at native model resolution
         !DO L = 1, LLPAR


            !kjw. JLOP and CSPEC are variables in comode_mod.f, associated
            !     with smvgear. For getting CH4 values at native model
            !     resolution, write my own thing, it'll probably have to do
            !     with using CH4 from STT when below tropopause, and using
            !     the TES retrieval above the tropopause. Therefore, there
            !     will be no adjoint forcing above the tropopause.
            !kjw.

               !kjw
               ! Get CH4 from restored STT array
               ! Find out units of STT (I think it's [kg/box]).
               !     If so, convert to [v/v] before next step
            DO L=1,LLPAR
               GC_CH4_NATIVE(L) = CHK_STT(I,J,L,1) * PERT(L)
            ENDDO

         ! Unit Conversion
         DO L=1,LLPAR

            ! Convert from [kg/box] --> [v/v]
            !    Numerator   = moles CH4/box
            !    Denominator = moles air/box
            GC_CH4_NATIVE(L) = (GC_CH4_NATIVE(L)*XNUMOL(1)/6.022d23 ) / 
     &                      ( AD(I,J,L) * XNUMOLAIR / 6.022d23 )
         ENDDO
            !!!! Bypass unit conversion
!!!            GC_CH4_NATIVE(:) = GC_CH4_NATIVE(:)

         ! Interpolate GC CH4 column to TES grid 
         DO LL = 1, LTES
            GC_CH4(LL) = 0d0 
            DO L = 1, LLPAR 
               GC_CH4(LL) = GC_CH4(LL) 
     &                    + MAP(L,LL) * GC_CH4_NATIVE(L) 
            ENDDO
         ENDDO
!!!         ! Bypass interpolation
!!!         GC_CH4(1:47) = GC_CH4_NATIVE(:)
!!!

!         !kjw testing
!         IF ( ori .EQ. .TRUE. ) THEN
!            print*,'-------------------------'
!            print*,'GC_CH4(13) = ',GC_CH4(13)
!         ENDIF
!         !kjw testing


!         ! dkh debug: compare profiles:
!         print*, ' GC_PRES, GC_native_CH4 [ppb] '
!         WRITE(6,100) (GC_PRES(L), GC_CH4_NATIVE(L)*1d9, 
!     &                   L = LLPAR, 1, -1 )
!         print*, ' TES_PRES, GC_CH4  ' 
!         WRITE(6,100) (TES(NT)%PRES(LL), 
!     &                 GC_CH4(LL)*1d9, LL = LTES, 1, -1 ) 
! 100  FORMAT(1X,F16.8,1X,F16.8)


         !--------------------------------------------------------------
         ! Apply TES observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by TES [lnvmr]
         !    x_a   = TES apriori column               [lnvmr]
         !    x_m   = GC modeled column                [lnvmr]
         !    A_k   = TES averaging kernel 
         !--------------------------------------------------------------
         ! x_m - x_a
         DO L = 1, LTES 
           GC_CH4(L)   = MAX(GC_CH4(L), 1d-10)
           CH4_PERT(L) = LOG(GC_CH4(L)) - LOG(TES(NT)%PRIOR(L))
         ENDDO

!         !kjw testing
!         IF ( ori .EQ. .TRUE. ) THEN
!            print*,'CH4_PERT(13) = ',CH4_PERT(13)
!         ENDIF
!         !kjw testing
!!!         ! Bypass !x_m - x_a
!!!         CH4_PERT(:) = GC_CH4(:)
!!!
     
         ! x_a + A_k * ( x_m - x_a ) 
         !  AVG_KERNEL indexing may look backwards because BPCH files storing AK 
         !    values use IDL column major indexing.
         DO L = 1, LTES
            CH4_HAT(L)    = 0d0 
            DO LL = 1, LTES
               CH4_HAT(L) = CH4_HAT(L) 
     &                    + TES(NT)%AVG_KERNEL(LL,L) * CH4_PERT(LL) 
            ENDDO
            CH4_HAT(L)    = CH4_HAT(L) + LOG(TES(NT)%PRIOR(L))
         ENDDO
!!!         ! Bypass !x_a + A_k * ( x_m - x_a )  
!!!         CH4_HAT(:) = CH4_PERT(:)
!!!


!         !--------------------------------------------------------------
!         ! Calculate column average ln(vmr) weighted by # density of TES grid
!         !    This operation is self adjoint.
!         ! To get # density [molec / m^2], get pressure differences in grid.
!         !    Get kg/m^2 of air in each box. (F=ma) dP=kg*g
!         !    Get molec/m^2 air  molec/m2 = kg/m2 * XNUMOLAIR
!         ! TES pressure grid linear in ln(pres).  Get pressure at edge of boxes
!         !--------------------------------------------------------------
!         pres_ln(1:LTES)=LOG(TES(NT)%pres(1:LTES))   ! [hPa] --> ln([hPa])
!         Pedges_ln(1)=pres_ln(1)             ! Bottom edge is surface pressure
!         DO L=2,LTES-1
!            Pedges_ln(L) = ( pres_ln(L) + pres_ln(L+1) ) / 2.
!         ENDDO
!         Pedges=EXP(pedges_ln)               ! ln([hPa]) --> [hPa]
!         Pedges(LTES)=0                      ! Top of atmosphere
!         !print*,' kjw debug:  TES pressure edges'
!         !print*,         Pedges
!
!         ! Calculate pressure difference of each LTES-1 boxes
!         DO L=1,LTES-1
!            Pdiff(L) = Pedges(L) - Pedges(L+1)
!         ENDDO
!
!         ! Calculate # molecules air in each LTES-1 obxes
!         Pdiff(:) = 100 * Pdiff(:)              ! [hPa] --> [Pa]
!         Pdiff(:) = Pdiff(:) / 9.8              ! [hPa] --> [kg]
!         Nmolec(:)= Pdiff(:) * XNUMOLAIR        ! [kg]  --> [molec]
!         Totmolec = 0d0
!         DO L=1,LTES-1
!            IF ( TROPP(I,J) < Pedges(L) ) THEN
!               Totmolec = Totmolec + Nmolec(L)
!            ENDIF
!         ENDDO
!
!         IF ( NT .EQ. 600 ) THEN
!            !print*,' kjw debug:# molecules in column (should be ~2e29)'
!            !print*,         Totmolec
!         ENDIF
!
!         ! Calculate column average ln(vmr) weighted by # density of TES levels
!         ! Only include levels with tropospheric air
!         OBS_avg = 0d0
!         GC_avg  = 0d0
!         DIFF_onTES(:) = 0d0
!         DIFF_onTES(1:LTES) = LOG(TES(NT)%GC_CH4(1:LTES)) -
!     &                               CH4_HAT(1:LTES)
!         DO L=1,LTES-1
!            IF ( TROPP(I,J) < Pedges(L) ) THEN
!               OBS_avg = OBS_avg + 
!     &                     LOG(TES(NT)%GC_CH4(L+1)) * Nmolec(L)/Totmolec
!               GC_avg  = GC_avg  + 
!     &                     CH4_HAT(L+1) * Nmolec(L)/Totmolec
!            ENDIF
!         ENDDO


         ! Transform from [ln(vmr)] --> [vmr]
         CH4_HAT_EXP = EXP(CH4_HAT)

         ! Calculate RTln(VMR) for profiles.
         CALL GET_RTVMR( NT, TES(NT)%GC_CH4, OBS_RTVMR, M_STAR )
         CALL GET_RTVMR( NT, CH4_HAT_EXP,    GC_RTVMR,  M_STAR )


         ! Retrieve RTVMR from TES structure for pseudo-observations
         IF ( LTES_PSO ) THEN
            OBS_RTVMR = 0d0
            OBS_RTVMR = TES(NT)%GC_CH4(67)
         ENDIF


         !--------------------------------------------------------------
         ! Calculate cost function, given S is error on ln(vmr)
         ! J = [ model - obs ]^T S_{obs}^{-1} [ ln(model - obs ]
         !--------------------------------------------------------------

         ! DIFF = model - obs.     units: [ln(vmr)] / m^2
         DIFF = GC_RTVMR - OBS_RTVMR


         ! Calculate DIFF^T * S_{obs}^{-1} * DIFF
         FORCE = 0d0
         FORCE = 2 * DIFF * S_obs_inv
         NEW_COST = 0.5d0 * DIFF * FORCE


!         ! Calculate difference between modeled and observed profile
!         ! Eliminate stratospheric forcing in this 
!         DO L = 1, LTES 
!            IF ( TROPP(I,J) < TES(NT)%PRES(L) ) THEN 
!               DIFF(L) = CH4_HAT(L) - LOG( TES(NT)%GC_CH4(L) )
!            ELSE 
!               DIFF(L) = 0d0
!            ENDIF 
!         ENDDO 
!         ! Bypass DIFF
!!!!         DIFF(:) = CH4_HAT(:)
!
!
!         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF 
!         DO L = 1, LTES
!            FORCE(L)     = 0d0 
!            DO LL = 1, LTES
!               FORCE(L)  = FORCE(L) + TES(NT)%S_OER_INV(L,LL) *DIFF(LL)
!            ENDDO
!            NEW_COST = NEW_COST + 0.5d0 * DIFF(L) * FORCE(L)
!            !NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
!         ENDDO

!!!         !Bypass this part of adjoint
!!!         DO L=1,LTES
!!!            FORCE(L) = 0.5d0
!!!            NEW_COST=NEW_COST+ 0.5d0 * DIFF(L)
!!!         ENDDO



!         ! dkh debug: compare profiles:
!         print*, ' CH4_HAT, CH4_TES, CH4_GC [ppb]'
!         WRITE(6,090) ( 1d9 * EXP(CH4_HAT(L)),
!     &                  1d9 * TES(NT)%CH4(L),
!     &                  1d9 * TES(NT)%GC_CH4(L),
!     &                  L,    L = LTES, 1, -1   )
!
!         print*, ' TES_PRIOR, CH4_HAT, CH4_GC  [lnvmr], diag(S^-1)'
!         WRITE(6,101) ( LOG(TES(NT)%PRIOR(L)), CH4_HAT(L), 
!     &                  LOG(TES(NT)%GC_CH4(L)), TES(NT)%S_OER_INV(L,L),
!     &                  L,  L = LTES, 1, -1 )
! 090  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,i3)
! 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,d14.6,1x,i3)


!         !--------------------------------------------------------------
!         ! Begin adjoint calculations 
!         !--------------------------------------------------------------

      !kjw. We've now calculated:
      !        1) forcing ( 2 * S_{obs}^{-1} * DIFF = FORCE ) units of [lnvmr]^{-1}
      !        2) New contribution to cost function do to diff
      !             DIFF * S_{obs}^{-1} * DIFF
      !     This has all been done on the TES pressure grid
      !
      !     At this point, we need to initialize the adjoint variable: STT_ADJ
      !     Do so by applying the adjoint of all operators used to get
      !               STT --> ln(vmr) for calculating ( F(x)-y )

         ! dkh debug
!        print*, 'DIFF , FORCE ' 
!         WRITE(6,102) (DIFF(L), FORCE(L), 
!     &       L = LTES, 1, -1 )
! 102  FORMAT(1X,d14.6,1X,d14.6)

         ! The adjoint forcing is  2 * S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE

         ! Adjoint of difference
         GC_RTVMR_ADJ = DIFF_ADJ


         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !ADJ_DIFF(:) = 1d0
         !NEW_COST(NT) = ?? SUM(ABS(LOG(CH4_HAT(1:LTES))))
         !print*, ' sumlog =', SUM(ABS(LOG(CH4_HAT(:))))
         !print*, ' sumlog =', ABS(LOG(CH4_HAT(:)))
 

!         ! Adjoint of Column Averaging
!         DO L = 1, LTES-1
!            IF ( TROPP(I,J) < Pedges(L) ) THEN 
!               CH4_HAT_ADJ(L+1) =  GC_avg_ADJ * Nmolec(L)/Totmolec
!            ELSE
!               CH4_HAT_ADJ(L+1) = 0d0
!            ENDIF 
!         ENDDO 

         ! Adjoint of RTVMR Averaging
         DO L=1,LTES
            CH4_HAT_ADJ(L) = M_STAR(2,L) * GC_RTVMR_ADJ
         ENDDO

         ! Adjoint of ln(vmr) --> vmr
         DO L=1,LTES
            IF ( CH4_HAT_ADJ(L) /= 0.0 ) THEN
               CH4_HAT_EXP_ADJ(L) = CH4_HAT_ADJ(L) * CH4_HAT_EXP(L)
            ELSE
               CH4_HAT_EXP_ADJ(L) = 0d0
            ENDIF
         ENDDO

!         ! Adjoint of difference
!         DO L = 1, LTES 
!            IF ( TROPP(I,J) < TES(NT)%PRES(L) ) THEN 
!               CH4_HAT_ADJ(L) =  DIFF_ADJ(L)
!            ELSE
!               CH4_HAT_ADJ(L) = 0d0
!            ENDIF 
!         ENDDO 
!!!         ! Bypass adjoint of difference
!!!         CH4_HAT_ADJ(:) = DIFF_ADJ(:)



!!!      ! adjoint of TES operator
         DO L  = 1, LTES
            CH4_PERT_ADJ(L)    = 0d0
            DO LL = 1, LTES
               CH4_PERT_ADJ(L) = CH4_PERT_ADJ(L) 
     &                        + TES(NT)%AVG_KERNEL(L,LL) 
     &                        * CH4_HAT_EXP_ADJ(LL)
          ENDDO
         ENDDO
!!!         ! Bypass adjoint of TES operator
!!!         CH4_PERT_ADJ(:) = CH4_HAT_ADJ(:)
!!!

         ! Adjoint of x_m - x_a
         DO L = 1, LTES 
           ! fwd code:
           !GC_CH4(L)   = MAX(GC_CH4(L), 1d-10)
           !CH4_PERT(L) = LOG(GC_CH4(L)) - LOG(TES(NT)%PRIOR(L))
           ! adj code:
           IF ( GC_CH4(L) > 1d-10 ) THEN 
              GC_CH4_ADJ(L) = 1d0 / GC_CH4(L) * CH4_PERT_ADJ(L)
           ELSE 
              GC_CH4_ADJ(L) = 1d0 / 1d-10 * CH4_PERT_ADJ(L)
           ENDIF 
         ENDDO
!!!         ! Bypass adjoint of x_m - x_a
!!!         GC_CH4_ADJ(:) = CH4_PERT_ADJ(:)



         ! dkh debug
!         print*, 'CH4_HAT_ADJ, CH4_PERT_ADJ, GC_CH4_ADJ'
!         WRITE(6,103) (CH4_HAT_ADJ(L), CH4_PERT_ADJ(L), GC_CH4_ADJ(L), 
!     &       L = LTES, 1, -1 )
! 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)


!         ! adjoint of interpolation 
         DO L  = 1, LLPAR
            GC_CH4_NATIVE_ADJ(L) = 0d0 
            DO LL = 1, LTES
               GC_CH4_NATIVE_ADJ(L) = GC_CH4_NATIVE_ADJ(L)
     &                             + MAP(L,LL) * GC_CH4_ADJ(LL)
            ENDDO
         ENDDO
!!!         ! Bypass adjoint of interpolation
!!!         GC_CH4_NATIVE_ADJ(:) = GC_CH4_ADJ(1:47)
!!!

         ! kjw
         ! Adjoint of interpolation leaves GC_CH4_NATIVE_ADJ with some zeros
         !   in the lower troposphere.  This occurs because the GC pres grid is
         !   finer in lower troposphere than the TES grid. So, when interpolating
         !   from GC --> TES, the contribution from some GC grid boxes to any TES
         !   grid box is zero.  Unfortunately, when we go back from TES to GC grid,
         !   this means that some 


!         WRITE(114,112) ( GC_CH4_NATIVE_ADJ(L),     L=LLPAR,1,-1)

         ! dkh debug
         !print*, 'GC_CH4_NATIVE_ADJ 1 '
         ! WRITE(6,104) (GC_CH4_NATIVE_ADJ(L), L = LLPAR, 1, -1 )

         DO L = 1, LLPAR 

            ! Adjoint of unit conversion 
            GC_CH4_NATIVE_ADJ(L) = ( GC_CH4_NATIVE_ADJ(L) * 
     &                                   XNUMOL(1) / 6.022d23 ) / 
     &                     ( AD(I,J,L) * XNUMOLAIR / 6.022d23 )
            !!!! Bypass adjoint of unit conversion
            !!!GC_CH4_NATIVE_ADJ(:) = GC_CH4_NATIVE_ADJ(:)

            ! Just to make sure we're only forcing the troposphere
            IF ( ITS_IN_THE_TROP(I,J,L) ) THEN 
 
               ! Pass adjoint back to adjoint tracer array
               STT_ADJ(I,J,L,1) =
     &            STT_ADJ(I,J,L,1) + GC_CH4_NATIVE_ADJ(L)

               ADJ(L) = GC_CH4_NATIVE_ADJ(L) * CHK_STT(I,J,L,1)
            ELSE
               ADJ(L) = 0
            ENDIF 

         ENDDO


!         ! dkh debug
!         print*, 'GC_CH4_NATIVE_ADJ conv '
!         WRITE(6,104) (GC_CH4_NATIVE_ADJ(L), L = LLPAR, 1, -1 )
! 104  FORMAT(1X,d14.6)
!
!
!         WRITE(101,110) ( TES(NT)%PRES(LL),        LL=LTES,1,-1)
!         WRITE(102,110) ( 1d9 * GC_CH4(LL),         LL=LTES,1,-1)
!         WRITE(103,110) ( 1d9 * TES(NT)%CH4(LL),    LL=LTES,1,-1)
!         WRITE(104,110) ( 1d9 * TES(NT)%PRIOR(LL), LL=LTES,1,-1)
!         WRITE(105,110) ( DIFF(LL),                LL=LTES,1,-1)
!         WRITE(106,112) ( FORCE(LL),               LL=LTES,1,-1)
!         WRITE(107,111) NT, LTES
!         WRITE(108,112) ( CH4_PERT_ADJ(LL),         LL=LTES,1,-1)
!         WRITE(109,112) ( GC_CH4_ADJ(LL),           LL=LTES,1,-1)
!         WRITE(110,110) ( 1d9 * EXP(CH4_HAT(LL)),   LL=LTES,1,-1)
!         WRITE(111,110) ( GC_PRES(L),              L=LLPAR,1,-1)
!         WRITE(112,110) ( 1d9 * GC_CH4_NATIVE(L),   L=LLPAR,1,-1)
!         WRITE(113,110) ( 1d9 * GC_CH4(LL),         LL=LTES,1,-1)
! 110     FORMAT(F18.6,1X)
! 111     FORMAT(i4,1X,i4,1x)
! 112     FORMAT(D14.6,1X)
!
!
!      ENDDO  ! NT
!!!$OMP END PARALLEL DO
!
!      WRITE(116,212) TES(NT)%TIME(1)
! 212  FORMAT(F22.6)

      ! Update cost function 
      ! COST_FUNC = SUM( NEW_COST(NTSTART:NTSTOP))
      COST_FUNC_A = NEW_COST


      ! Return to calling program
      END SUBROUTINE CALC_TES_CH4_FORCE_FD

!!------------------------------------------------------------------------------

      SUBROUTINE MAKE_PSEUDO_OBS( YYYYMMDD, NTES )
!
!******************************************************************************
!  Subroutine MAKE_PSEUDO_OBS populates TES%GC_CH4 with processed GC columns 
!  Processing consists of adding error and applying TES observation operator
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTES     (INTEGER) : Number of TES retrievals in this day 
!  (2 ) YYYYMMDD (INTEGER) : Current model date
!
!  Arguments as Output:
!  ============================================================================
!  (0 )
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE FILE_MOD,           ONLY : IU_FILE, IOERROR
      USE DIRECTORY_MOD,      ONLY : RUN_DIR
      USE DIRECTORY_ADJ_MOD,  ONLY : ADJTMP_DIR, DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH, GET_DAY
      USE BPCH2_MOD,          ONLY : GET_TAU0, READ_BPCH2
      USE TRANSFER_MOD,       ONLY : TRANSFER_2D, TRANSFER_3D
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC, EXPAND_NAME


#     include      "CMN_SIZE"      ! Size params      

      ! Arguments
      INTEGER, INTENT(IN)   :: NTES
      INTEGER, INTENT(IN)   :: YYYYMMDD

      ! Local variables
      CHARACTER(LEN=255)           :: FILENAME
      CHARACTER(LEN=255)           :: FILENAME_ROOT
      CHARACTER(LEN=255)           :: FILENAME_OBS
      INTEGER                      :: IOS, NT, I, J, L 
      INTEGER                      :: IIJJ(2), LL, LTES
      INTEGER                      :: DD, MM, YYYY
      INTEGER                      :: HH, HH_last, HHMMSS
      TYPE (XPLEX)                       :: ARRAY2D(IIPAR,JJPAR,1)
      TYPE (XPLEX)                       :: ARRAY3D(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                       :: MAP(LLPAR,MAXLEV)
      TYPE (XPLEX)                       :: GC_CH4(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                       :: GC_PRES(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                       :: GC_PSURF(IIPAR,JJPAR)
      TYPE (XPLEX)                       :: GC_PRES_this(LLPAR)
      TYPE (XPLEX)                       :: GC_CH4_NATIVE_this(LLPAR)
      TYPE (XPLEX)                       :: GC_CH4_this(MAXLEV)
      TYPE (XPLEX)                       :: CH4_PERT(MAXLEV)
      TYPE (XPLEX)                       :: CH4_HAT(MAXLEV)
      TYPE (XPLEX)                       :: GC_PSURF_this
      TYPE (XPLEX)                       :: XTAU
      TYPE (XPLEX)                       :: day_frac
      TYPE (XPLEX)                       :: GC_PSO_RTVMR
      TYPE (XPLEX)                       :: GC_PSO_RTVMR_werr
      TYPE (XPLEX)                       :: M_STAR(4,MAXLEV)


      !=================================================================
      ! MAKE_PSEUDO_OBS begins here!
      !=================================================================

      ! ----------------------------------------------------------------
      ! Get Today's Error values
      WRITE(6,'(a)') ' MAKE_PSEUDO_OBS - reading random errors'

      ! filename 
      FILENAME = TRIM( 'tes_ch4_random_YYYYMMDD.txt' )

      ! Expand date tokens in filename 
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 9999 ) 

      ! Construct complete filename 
      FILENAME = TRIM('/home/kjw/TES/data/V004/bpch/randoms/')  // 
     &           TRIM( FILENAME )

      ! Open file
      print*,'Opening: ', TRIM(FILENAME)
      OPEN( IU_FILE,      FILE=TRIM( FILENAME ),          
     &      STATUS='OLD', IOSTAT=IOS                    )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_rand_file:1')

      ! Read each line
      DO NT=1,NTES+1
         READ( IU_FILE, '(F16.12)', IOSTAT=IOS ) TES(NT)%ERR(1)
         IF ( NT .EQ. NTES+1 ) THEN
            IF ( IOS < 0 ) THEN
               WRITE(6,'(a)') 'Done reading random errors'
            ELSE
               WRITE(6,'(a)') 'Unexpected end. read_rand_file:2'
            ENDIF
         ENDIF
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      print*,'1TES(300)%ERR(1) = ',TES(300)%ERR(1)

      ! Errors have mean=0, stddev=1. Multiply by ERR_PPB (20-40) / 1d9
      print*,ERR_PPB
      DO NT=1,NTES
         TES(NT)%ERR(1) = TES(NT)%ERR(1) * ERR_PPB / 1.0d9
      ENDDO
            print*,'2TES(300)%ERR(1) = ',TES(300)%ERR(1)


      ! ----------------------------------------------------------------
      ! Get Today's GEOS-Chem columns
      WRITE(6,'(a)') ' MAKE_PSEUDO_OBS - Reading GEOS-Chem observations'

      ! Set HH_last to -1 to guarantee opening file for 1st observation
      HH_last = -1

      ! Loop over each TES observation
      DO NT=1,NTES

         print*, '     - MAKE_PSEUDO_OBS: for record ', NT

         ! Get date and hour of observation (round time)
         YYYY = floor( ( YYYYMMDD * 1d-4 ) )
         MM   = floor( ( YYYYMMDD - 1d4 * YYYY ) * 1d-2 )
         DD   = floor( ( YYYYMMDD - 1d4 * YYYY - 1d2 * MM ) )
         day_frac = TES(NT)%TIME(1) - 1d4 * YYYY -
     &                                1d2 * MM   - DD
         HH = nint( 24. * day_frac )
         IF ( HH == 24 ) HH = 23    ! If last 1/2 hour of day, use HH=23
         HHMMSS = 1d4 * HH

         ! Open new obs file if necessary
         IF ( HH /= HH_last ) THEN

            FILENAME_ROOT = TRIM( RUN_DIR )
            FILENAME_OBS  = 'gctm.obs.YYYYMMDD.hhmm'

            print*,'TES(NT)%TIME(1) = ',TES(NT)%TIME(1)
            print*,'GET_YEAR = ', GET_YEAR()
            print*,'GET_MONTH = ', GET_MONTH()
            print*,'GET_DAY = ', GET_DAY()
            print*,'YYYY = ',YYYY
            print*,'MM = ',MM
            print*,'HH = ',HH
            print*,'HHMMSS = ',HHMMSS
            print*,'FILENAME_OBS = ',FILENAME_OBS
            CALL EXPAND_DATE( FILENAME_OBS, YYYYMMDD, HHMMSS )
            FILENAME_OBS  = TRIM( ADJTMP_DIR   )   //   
     &                      TRIM( FILENAME_OBS )
            print*,'FILENAME_OBS = ', FILENAME_OBS

            ! Get Tau value for BPCH read
            XTAU = GET_TAU0( MM, DD, YYYY, HH )

            ! Get 3D array of GEOS-Chem values
            print*,'Read observations'
            GC_CH4(:,:,:)  = 0d0
            ARRAY3D(:,:,:) = 0d0
            CALL READ_BPCH2( TRIM(FILENAME_OBS), 'IJ-OBS-$', 1,
     &                       XTAU,          IIPAR,     JJPAR,
     &                       LLPAR,         ARRAY3D,   QUIET=.TRUE. )
            !CALL TRANSFER_3D( ARRAY3D(:,:,:), GC_CH4(:,:,:) )
            GC_CH4(:,:,:) = ARRAY3D(:,:,:)

            ! Get 3D array of GEOS-Chem pressure centers
            print*,'Read pressure centers'
            ARRAY3D(:,:,:) = 0d0
            GC_PRES(:,:,:) = 0d0
            CALL READ_BPCH2( TRIM(FILENAME_OBS), 'IJ-OBS-$', 2,
     &                       XTAU,          IIPAR,      JJPAR,
     &                       LLPAR,         ARRAY3D,    QUIET=.TRUE. )
            !CALL TRANSFER_3D( ARRAY3D(:,:,:), GC_PRES(:,:,:) )
            GC_PRES(:,:,:) = ARRAY3D(:,:,:)

            ! Get 2D array of GEOS-Chem surface pressure
            print*,'Read surface pressure'
            GC_PSURF(:,:) = 0d0
            ARRAY2D(:,:,1)  = 0d0
            CALL READ_BPCH2( TRIM(FILENAME_OBS), 'IJ-OBS-$', 3,
     &                       XTAU,          IIPAR,     JJPAR,
     &                       1,             ARRAY2D,     QUIET=.TRUE. )
            !CALL TRANSFER_2D( ARRAY2D(:,:,1), GC_PSURF(:,:) )
            GC_PSURF(:,:) = ARRAY2D(:,:,1)

         ENDIF

         ! RESET a few variables to be safe
         I = 0
         J = 0
         GC_PRES_this(:) = 0d0
         GC_PSURF_this   = 0d0
         CH4_HAT(:)      = 0d0 
         MAP(:,:)        = 0d0
         GC_PSO_RTVMR    = 0d0

         ! Copy LTES to make coding a bit cleaner
         LTES = TES(NT)%LTES(1)

         ! Get grid box of current record
         IIJJ  = GET_IJ(XPLX(TES(NT)%LON(1)),XPLX(TES(NT)%LAT(1)))
         I     = IIJJ(1)
         J     = IIJJ(2)

         ! Get GC pressure levels (mbar) 
         DO L = 1, LLPAR
            GC_PRES_this(L)       = GC_PRES(I,J,L)
            GC_CH4_NATIVE_this(L) = GC_CH4(I,J,L)! + TES(NT)%ERR(1)
         ENDDO


         ! Get GC surface pressure (mbar) 
         GC_PSURF_this = GC_PSURF(I,J) 

         ! Calculate the interpolation weight matrix 
         MAP(1:LLPAR,1:LTES) 
     &      = GET_INTMAP( LLPAR, GC_PRES_this(:),      GC_PSURF_this, 
     &                    LTES,  TES(NT)%PRES(1:LTES), TES(NT)%PRES(1) )


         IF ( NT == 600 ) THEN
            !print*,'LTES = ',LTES
            !print*,'LLPAR = ',LLPAR
            !print*,'GC_PSURF = ',GC_PSURF_this

            !WRITE(6,'(a)') 'GEOS-Chem pressure grid'
            !WRITE(6,'(F10.3)') ( GC_PRES_this(L), L=1,47 )

            !WRITE(6,'(a)') 'TES pressure grid'
            !WRITE(6,'(F10.3)') ( TES(NT)%PRES(L), L=1,65 )

            WRITE(6,'(a)') '20th row of MAP matrix in MAKE_PSEUDO_OBS'
            WRITE(6,'(5F8.5)') ( MAP(20:24,L) , L=65,1,-1)
         ENDIF

         ! Interpolate GC CH4 column to TES grid 
         DO LL = 1, LTES
            GC_CH4_this(LL) = 0d0 
            DO L = 1, LLPAR 
               GC_CH4_this(LL) = GC_CH4_this(LL) 
     &                    + MAP(L,LL) * GC_CH4_NATIVE_this(L) 
            ENDDO
         ENDDO

         !--------------------------------------------------------------
         ! Apply TES observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by TES [lnvmr]
         !    x_a   = TES apriori column               [lnvmr]
         !    x_m   = GC modeled column                [lnvmr]
         !    A_k   = TES averaging kernel 
         !--------------------------------------------------------------

         ! x_m - x_a
         DO L = 1, LTES 
           GC_CH4_this(L) = MAX(GC_CH4_this(L), 1d-10)
           CH4_PERT(L)    = LOG(GC_CH4_this(L)) - LOG(TES(NT)%PRIOR(L))
         ENDDO
     
         ! x_a + A_k * ( x_m - x_a )  
         DO L = 1, LTES
            CH4_HAT(L)    = 0d0 
            DO LL = 1, LTES
               CH4_HAT(L) = CH4_HAT(L) 
     &                    + TES(NT)%AVG_KERNEL(LL,L) * CH4_PERT(LL) 
            ENDDO
            CH4_HAT(L)    = CH4_HAT(L) + LOG(TES(NT)%PRIOR(L))
         ENDDO
         ! Indexing of Averaging Kernel is seemingly backwards because
         ! TES observation files processed in IDL, which is column major


         ! Get RTVMR of GEOS-Chem column w/ TES obs operator applied
         CALL GET_RTVMR( NT, EXP(CH4_HAT), GC_PSO_RTVMR, M_STAR )

         ! Add random error w/ standard deviation = ERR_PPB to RTVMR
         GC_PSO_RTVMR_werr = GC_PSO_RTVMR + TES(NT)%ERR(1)


         ! Place pseudo-observation RTVMR [v/v] in TES structure
         TES(NT)%GC_CH4(67) = GC_PSO_RTVMR_werr

         IF ( NT == 600 ) THEN
            print*,'Error [v/v] = ', TES(NT)%ERR(1)
            print*,'GC_PSO_RTVMR = ',GC_PSO_RTVMR
            print*,'GC_PSO_RTVMR_werr = ',GC_PSO_RTVMR_werr
         ENDIF

         ! Make this hour equal to last hour
         HH_last = HH
         HH      = 0


!         ! Check GEOS-Chem, Error, and CH4_HAT for a given observation
!         ! Success! kjw, 07/25/10
!         IF ( NT == 600 ) THEN
!            ! Write values for one observation to check that it's right
!            FILENAME = 'test_pseudo_obs.NN.m'
!            CALL EXPAND_NAME( FILENAME, N_CALC )
!            FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!            OPEN( IU_FILE,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
!     &          IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!            
!            WRITE(IU_FILE,'(a4,i4)') 'I = ', I
!            WRITE(IU_FILE,'(a4,i4)') 'J = ', J
!            WRITE(IU_FILE,'(a12,F8.3)') 'TES PSURF = ', TES(NT)%PRES(1)
!            WRITE(IU_FILE,'(a12,F8.3)') 'GC  PSURF = ', GC_PSURF_this
!            WRITE(IU_FILE,'(a10,F16.12)') 'ERR_PPB =',1d9*TES(NT)%ERR(1)
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'GEOS-Chem CH4 Native'
!            WRITE(IU_FILE,'(F16.8)') (1d9*GC_CH4(I,J,L), L=1,47)
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'GEOS-Chem CH4 Native w Error'
!            WRITE(IU_FILE,'(F16.8)') (1d9*GC_CH4_native_this(L),L=1,47)
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'GEOS-Chem CH4 on TES'
!            WRITE(IU_FILE,'(F16.8)') (1d9*GC_CH4_this(L), L=1,65)
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'TES a priori'
!            WRITE(IU_FILE,'(F16.8)') ( 1d9*TES(NT)%PRIOR(L), L=1,65 )
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'CH4 HAT'
!            WRITE(IU_FILE,'(F24.12)') ( CH4_HAT(L), L=1,65 )
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'EXP( CH4 HAT )'
!            WRITE(IU_FILE,'(F16.8)') ( 1d9*EXP(CH4_HAT(L)), L=1,65 )
!            WRITE(IU_FILE,'(a)') '-------------------------------------'
!            WRITE(IU_FILE,'(a)') 'GC_CH4 HAT'
!            WRITE(IU_FILE,'(F16.8)') ( 1d9*TES(NT)%GC_CH4(L), L=1,65 )
!
!            CLOSE(IU_FILE)
!         ENDIF


      ENDDO   ! End looping over each observation





      ! Return to calling program
      END SUBROUTINE MAKE_PSEUDO_OBS

!------------------------------------------------------------------------------

      SUBROUTINE CALC_TES_GC_BIAS
!
!******************************************************************************
!  Subroutine CALC_TES_GC_BIAS calculates mean TES bias w.r.t. GEOS-Chem during
!    the entire simulation period. Bias is then stored in module variable BIAS
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE TIME_MOD,           ONLY : GET_TAUb, GET_TAUe
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,           ONLY : GET_TIME_BEHIND_ADJ
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE DIRECTORY_ADJ_MOD,  ONLY : ADJTMP_DIR
      USE BPCH2_MOD,          ONLY : READ_BPCH2, GET_RES_EXT
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE GRID_MOD,           ONLY : GET_IJ
      USE FILE_MOD,           ONLY : IU_FILE,      IOERROR
      USE GRID_MOD,           ONLY : GET_AREA_M2

#     include "CMN_SIZE"

      ! Arguments

      ! Local Variables
      LOGICAL               :: file_exist
      CHARACTER(LEN=255)    :: TES_dir, READ_FILENAME
      CHARACTER(LEN=255)    :: PRS_ROOTNAME, PRS_FILENAME, CHK_FILENAME
      INTEGER               :: NTES, I, J, LTES, NHITS, hh
      INTEGER               :: NT, nday, N, L, ND49_NT, IOS
      INTEGER               :: nymd0, ND49_NTES, LL
      INTEGER               :: IIJJ(2)
      INTEGER               :: MATCHES(MAXTES)
      INTEGER               :: NTSTART, NTSTOP
      TYPE (XPLEX)                :: GC_RTVMRt, OBS_RTVMRt
      TYPE (XPLEX)                :: OBS_RTVMR_today(MAXTES)
      TYPE (XPLEX)                :: GC_RTVMR_today(MAXTES)
      TYPE (XPLEX)                :: OBS_RTVMR_tot(1000)
      TYPE (XPLEX)                :: GC_RTVMR_tot(1000)
      TYPE (XPLEX)                :: NOBS_tot(1000)
      TYPE (XPLEX)                :: ARRAY0(1)
      TYPE (XPLEX)                :: ARRAY1(MAXTES)
      TYPE (XPLEX)                :: ND49_lat(MAXTES)
      TYPE (XPLEX)                :: ND49_lon(MAXTES)
      TYPE (XPLEX)                :: ND49_PSURF(MAXTES)
      TYPE (XPLEX)                :: ARRAY2(MAXTES,LLPAR,1)
      TYPE (XPLEX)                :: ND49_PCEN(MAXTES,LLPAR)
      TYPE (XPLEX)                :: ND49_PEDGE(MAXTES,LLPAR)
      TYPE (XPLEX)                :: ND49_kg_box(MAXTES,LLPAR)
      TYPE (XPLEX)                :: tau0
      TYPE (XPLEX)                :: date0(2)
      TYPE (XPLEX)                :: TIME_FRAC(MAXTES)
      TYPE (XPLEX)                :: TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                :: GC_chk(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                :: GC_ch4_kg(LLPAR)
      TYPE (XPLEX)                :: GC_PRES(LLPAR)
      TYPE (XPLEX)                :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                :: CH4_PERT(MAXLEV)
      TYPE (XPLEX)                :: GC_CH4(MAXLEV)
      TYPE (XPLEX)                :: GC_PSURF
      TYPE (XPLEX)                :: MAP(LLPAR,MAXLEV)
      TYPE (XPLEX)                :: M_STAR(4,MAXLEV)
      TYPE (XPLEX)                :: CH4_HAT(MAXLEV)
      TYPE (XPLEX)                :: CH4_HAT_EXP(MAXLEV)
      TYPE (XPLEX)                :: TAUb, TAUe
      TYPE (XPLEX)                :: BIAS_tot, N_tot

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL,    NV
      INTEGER             :: IJLOOP
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! CALC_TES_GC_BIAS begins here!
      !=================================================================

      ! TES observation root directory
      TES_dir = '/home/kjw/TES/data/V004/bpch/'


      print*,' CALC_TES_GC_BIAS'
      print*,'   current NYMD, NHMS = ',GET_NYMD(), GET_NHMS()


      ! Set TAUb and TAUe
      TAUb = GET_TAUb()
      TAUe = GET_TAUe()

      ! Set some variables for first iteration of loop
      tau0  = TAUe - 24d0
      nday  = 1

      ! Loop over every day in assimilation period
      DO WHILE ( tau0 > TAUb )

         ! Zero arrays to be safe
         GC_RTVMR_today(:) = 0d0
         OBS_RTVMR_today(:) = 0d0

         ! Get NYMD of the day
         date0 = GET_TIME_BEHIND_ADJ( 1380 + (nday-1)*1440 )
         nymd0 = date0(1)
         print*,'nymd0 = ',nymd0
         print*,'tau0 = ',tau0

         ! Get filename of TES observations
         READ_FILENAME = TRIM( 'tes_ch4_YYYYMMDD.bpch' )
         CALL EXPAND_DATE( READ_FILENAME, nymd0, 9999 )
         READ_FILENAME = TRIM( TES_dir ) // TRIM( READ_FILENAME )

         ! Find whether observations exists on this day
         INQUIRE( FILE=READ_FILENAME, exist=file_exist )

         ! If the file exists, proceed.
         IF ( file_exist ) THEN

            ! Read TES_CH4_OBS during the day
            CALL READ_TES_CH4_OBS( nymd0, NTES )

            ! TIME is YYYYMMDD.frac-of-day.
            ! Subtract date and save just time frac
            TIME_FRAC(1:NTES) = TES(1:NTES)%TIME(1) - nymd0

            ! Open ND49 file and find the entries we want
            PRS_FILENAME = 'ND49_trim_' // GET_RES_EXT() //
     &                     '_YYYYMMDD.bpch'
            CALL EXPAND_DATE( PRS_FILENAME, nymd0, 9999 )
            PRS_ROOTNAME = '/home/kjw/GEOS-Chem/runs/ch4/TES/ND49_trim/'
            PRS_FILENAME = TRIM( PRS_ROOTNAME) // TRIM( PRS_FILENAME )


            ! Get # of observations in this BPCH file
            ! Read NTES from ND49 file
            CALL READ_BPCH2( PRS_FILENAME, 'IJ-AVG-$',     1,
     &                               tau0,          1,     1,
     &                                  1,  ARRAY0(1),  QUIET=.TRUE. )
            ND49_NTES = INT( ARRAY0(1) )


            !=================================================================
            ! Open binary punch file and read top-of-file header.
            ! Do some error checking to make sure the file is the right format.
            !=================================================================
            CALL OPEN_BPCH2_FOR_READ( IU_FILE, PRS_FILENAME )
            
            !=================================================================
            ! Read data from the binary punch file 
            !
            ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
            !=================================================================
            DO
               READ( IU_FILE, IOSTAT=IOS ) 
     &              MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
               
               IF ( IOS < 0 ) EXIT
               IF ( IOS > 0 ) CALL IOERROR(IOS,IU_FILE, 'tes_ch4_mod:1')
            
               READ( IU_FILE, IOSTAT=IOS ) 
     &              CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &              NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &              NSKIP
            
               IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_FILE,'tes_ch4_mod:2')
            
               ! Zero Dummy array
               ARRAY2(:,:,:) = 0d0
               READ( IU_FILE, IOSTAT=IOS ) 
     &              ( ( ( ARRAY2(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
            
               IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_FILE,'tes_ch4_mod:3')


               ! Test for a match
               IF ( 'IJ-AVG-$' == TRIM( CATEGORY ) ) THEN

                  ! Longitude
                  IF ( NTRACER == 3 ) THEN
                  WRITE(6,*) '         - Reading: Latitude ... '
                     ND49_lat  = ARRAY2(1:ND49_NTES,1,1)
               
                  ! Latitude
                  ELSEIF ( NTRACER == 4 ) THEN
                  WRITE(6,*) '         - Reading: Longitude ... '
                     ND49_lon  = ARRAY2(1:ND49_NTES,1,1)

               
                  ! Surface Pressure
                  ELSEIF ( NTRACER == 7 ) THEN
                  WRITE(6,*) '         - Reading: PSURF ... '
                     ND49_PSURF  = ARRAY2(1:ND49_NTES,1,1)
               
                  ! Pressure Centers
                  ELSEIF ( NTRACER == 10 ) THEN
                  WRITE(6,*) '         - Reading: Pressure ... '
                     ND49_PCEN(:,:) = ARRAY2(1:ND49_NTES,1:LLPAR,1)

                  ENDIF  ! If tracer == # 
               
               ENDIF  ! If Category match

            ENDDO
            CLOSE( IU_FILE )



            ! Calculate Pressure edges from Pressure centers (approximate)
            DO N=1,ND49_NTES
            ND49_PEDGE(N,1) = ND49_PSURF(N)
            DO L=2,LLPAR
               ND49_PEDGE(N,L) = 0.5d0 * ( ND49_PCEN(N,L) + 
     &                                     ND49_PCEN(N,(L-1)) )
            ENDDO
            ENDDO

            ! Calculate kg air/box from Pressure edges
            ND49_PEDGE(:,:) = ND49_PEDGE(:,:) * 1d2    ! [hPa] --> [Pa]
            DO N=1,ND49_NTES
            DO L=1,LLPAR-1
               ND49_kg_box(N,L) = (ND49_PEDGE(N,L) -ND49_PEDGE(N,(L+1)))
     &                                 / 9.81
            ENDDO
            !ND49_kg_box(N,47) = ND49_PEDGE(N,47) / 9.81
            ND49_kg_box(N,LLPAR) = ND49_PEDGE(N,LLPAR) / 9.81
            ENDDO

            ! Convert 



            ! Associate ND49_NTES information with TES(NT) information
            ! Create array of indices of length = NTES.
            ! It should have values ex. [1, 2, 5, 7, ... ], associating each 
            !   NTES with the matching index of ND49_NTES
            nhits=1    ! nhits counts the # of matches we have.
                        ! nhits should = NTES when these loops finish
            DO NT=1,NTES
            DO ND49_NT=1,ND49_NTES
               IF ( TES(NT)%LAT(1) == ND49_lat(ND49_NT)  .AND.  
     &              TES(NT)%LON(1) == ND49_lon(ND49_NT) ) THEN
                  MATCHES(NT) = ND49_NT
                  nhits=nhits + 1
               ENDIF
            ENDDO
            ENDDO


            ! Loop over every hour in the day
            DO hh=23,0,-1

               ! Get NT range for this hour
               CALL GET_NT_RANGE( NTES, hh*10000, TIME_FRAC, 
     &                                     NTSTART, NTSTOP )

               ! If we have observations during this hour, proceed
               IF ( NTSTART /= 0 .OR. NTSTOP /= 0 ) THEN


                  ! Get GEOS-Chem CH4 values during this hour
                  !------------------------------------------------------
                  CHK_FILENAME = 'gctm.chk.YYYYMMDD.hhmm'
                  CALL EXPAND_DATE( CHK_FILENAME, nymd0, hh*10000 )
                  CHK_FILENAME = TRIM( ADJTMP_DIR   ) // 
     &                           TRIM( CHK_FILENAME )

                  ! Open the binary punch file for input
                  CALL OPEN_BPCH2_FOR_READ( IU_FILE, CHK_FILENAME )
                  READ( IU_FILE, IOSTAT=IOS )
     &                 MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

                  ! IOS < 0 is end-of-file, so exit
                  IF ( IOS < 0 ) EXIT

                  ! IOS > 0 is a real I/O error -- print error message
                  IF ( IOS > 0 ) 
     &                 CALL IOERROR( IOS,IU_FILE,'read_checkpt_file:7' )

                  READ( IU_FILE, IOSTAT=IOS )
     &                 CATEGORY, NTRACER, UNIT, ZTAU0, ZTAU1,  RESERVED,
     &                 NI,       NJ,      NL,   IFIRST,JFIRST, LFIRST,
     &                 NSKIP

                  IF ( IOS /= 0 ) 
     &                 CALL IOERROR(IOS,IU_FILE,'read_checkpt_file:8' )

                  READ( IU_FILE, IOSTAT=IOS )
     &                 ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

                  IF ( IOS /= 0 ) 
     &                 CALL IOERROR( IOS,IU_FILE,'read_checkpt_file:9' )

                  ! Convert from kg/box to [v/v]
                  DO J=1,JJPAR
                  DO I=1,IIPAR
                  DO L=1,LLPAR
                     GC_chk(I,J,L) = TRACER(I,J,L)
                  ENDDO
                  ENDDO
                  ENDDO

                  ! Close file
                  CLOSE( IU_FILE )
                  !------------------------------------------------------


                  !!!!!! GCvmr = GCkg / 16d-3 / ( kg_air(I,J,L) / 29d-3 )

                  ! Loop over observations during this hour
                  DO NT = NTSTART, NTSTOP, -1

                     ! For safety, initialize these up to LLTES 
                     GC_CH4_NATIVE(:) = 0d0
                     GC_CH4(:)      = 0d0
                     MAP(:,:)       = 0d0

                     ! Copy LTES to make coding a bit cleaner
                     LTES = TES(NT)%LTES(1)

                     ! Get grid box of current record
                     IIJJ  = GET_IJ(XPLX(TES(NT)%LON(1)), 
     &                              XPLX(TES(NT)%LAT(1))  )
                     I     = IIJJ(1)
                     J     = IIJJ(2)

                     ! Get GEOS-Chem CH4 in [v/v] from [kg/box]
                     GC_CH4_kg(:) = GC_chk(I,J,:)
                     DO L=1,LLPAR
                     GC_CH4_NATIVE(L) = GC_CH4_kg(L) * 29d-3 / 16d-3 / 
     &                  ( ND49_kg_box(MATCHES(NT),L) * GET_AREA_M2(J) )
                     ENDDO
!                     IF (NT == 600) THEN
!                     print*,'GC_CH4_NATIVE(14) = ', GC_CH4_NATIVE(14)
!                     print*,'GC_CH4_kg(14)     = ', GC_CH4_kg(14)
!                     print*,'ND49_kg_box(14)     = ', 
!     &                             ND49_kg_box(MATCHES(NT),14)
!                     print*,'ND49_PEDGE(14)     = ', 
!     &                             ND49_PEDGE(MATCHES(NT),14)
!                     ENDIF

                     ! Get GEOS-Chem pressure levels
                     GC_PRES(:) = ND49_PCEN(MATCHES(NT),:)
                     GC_PSURF   = ND49_PSURF(MATCHES(NT))

                     ! Calculate the interpolation weight matrix
                     
                     MAP(1:LLPAR,1:LTES) 
     &                   = GET_INTMAP( LLPAR, GC_PRES(:),   GC_PSURF, 
     &                     LTES,  TES(NT)%PRES(1:LTES), TES(NT)%PRES(1))

                     ! Interpolate GC CH4 column to TES grid 
                     DO LL = 1, LTES
                        GC_CH4(LL) = 0d0 
                        DO L = 1, LLPAR 
                           GC_CH4(LL) = GC_CH4(LL) 
     &                                + MAP(L,LL) * GC_CH4_NATIVE(L) 
                        ENDDO
                     ENDDO

                     !--------------------------------------------------------------
                     ! Apply TES observation operator
                     !
                     !   x_hat = x_a + A_k ( x_m - x_a ) 
                     !  
                     !  where  
                     !    x_hat = GC modeled column as seen by TES [lnvmr]
                     !    x_a   = TES apriori column               [lnvmr]
                     !    x_m   = GC modeled column                [lnvmr]
                     !    A_k   = TES averaging kernel 
                     !--------------------------------------------------------------
                     
                     ! x_m - x_a
                     DO L = 1, LTES 
                       GC_CH4(L)   =MAX(GC_CH4(L), 1d-10)
                       CH4_PERT(L) =LOG(GC_CH4(L))-LOG(TES(NT)%PRIOR(L))
                     ENDDO

                     ! x_a + A_k * ( x_m - x_a )  
                     DO L = 1, LTES
                        CH4_HAT(L)    = 0d0 
                        DO LL = 1, LTES
                           CH4_HAT(L) = CH4_HAT(L) 
     &                         + TES(NT)%AVG_KERNEL(LL,L) * CH4_PERT(LL) 
                        ENDDO
                        CH4_HAT(L) = CH4_HAT(L) + LOG(TES(NT)%PRIOR(L))
                     ENDDO

                     ! Transform from [ln(vmr)] --> [ppb]
                     CH4_HAT_EXP = EXP(CH4_HAT)

                     ! Calculate RTVMR for profiles.
                     CALL GET_RTVMR(NT,TES(NT)%CH4, OBS_RTVMRt, M_STAR)
                     CALL GET_RTVMR(NT,CH4_HAT_EXP, GC_RTVMRt,  M_STAR)

                     ! Save RTVMR values
                     GC_RTVMR_today(NT)  = GC_RTVMRt  * 1d9
                     OBS_RTVMR_today(NT) = OBS_RTVMRt * 1d9

                  ENDDO! End looping over each obs during this hour

               ENDIF   ! End if we have obs during this hour
            ENDDO      ! End looping over each hour during the day


         ! If the file does not exist, say so and move to next day
         ELSE
            WRITE(6,*) '    - CALC_TES_GC_BIAS: no files today: ',
     &                              TRIM( READ_FILENAME )
         ENDIF

         ! Average RTVMRs from the day
         OBS_RTVMR_tot(nday) = SUM( OBS_RTVMR_today )
         GC_RTVMR_tot(nday)  = SUM( GC_RTVMR_today  )
         NOBS_tot(nday) = NTES


         ! Increment Time counters
         tau0 = tau0 - 24
         nday = nday + 1
      ENDDO


      ! Calculate mean bias from OBS_RTVMR_tot and GC_RTVMR_tot
      BIAS_tot = SUM( OBS_RTVMR_tot - GC_RTVMR_tot )
      N_tot    = SUM( NOBS_tot      )
      BIAS_PPB = BIAS_tot / N_tot


      print*,'    - CALC_TES_GC_BIAS: '
      print*,'          GC_RTVMR_tot  = ', SUM( GC_RTVMR_tot )/N_tot
      print*,'          OBS_RTVMR_tot = ', SUM( OBS_RTVMR_tot )/N_tot
      print*,'          Total # observations = ', N_tot
      print*,'          Mean Bias [ppb] = ', BIAS_PPB
      print*,'        We hope mean bias ~ 110 ppb'


      ! Return to calling program
      END SUBROUTINE CALC_TES_GC_BIAS

!------------------------------------------------------------------------------


      SUBROUTINE GET_RTVMR( NT, VMR_IN, RTVMR_OUT, M_STAR )
!
!******************************************************************************
!  Subroutine GET_RTVMR returns Representative Tropospheric Volume Mixing Ratio
!  for a given column of ln(vmr).  RTVMR is described in Payne et. al. 2009
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NT  (INTEGER)      : TES observation #
!  (2 ) VMR_IN (REAL)      : CH4 column [ln(vmr)] from which to calculate RTVMR
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) RTVMR  (REAL)      : RTVMR calculated from CH4 column
!  (2 ) M_STAR (REAL)      : Normalized Mapping Matrix
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules

      ! Arguments
      INTEGER, INTENT(IN)    :: NT
      TYPE (XPLEX),  INTENT(IN)    :: VMR_IN(MAXLEV)
      TYPE (XPLEX),  INTENT(OUT)   :: RTVMR_OUT
      TYPE (XPLEX),  INTENT(OUT)   :: M_STAR(4,MAXLEV)
      !INTEGER                :: NT
      !TYPE (XPLEX)                 :: ln_VMR_IN(MAXLEV)
      !TYPE (XPLEX)                 :: RTVMR_OUT

      ! Local Variables
      INTEGER                :: L, LC, LTES
      TYPE (XPLEX)                 :: MAX_AK
      TYPE (XPLEX)                 :: FINE_GRID(MAXLEV)
      TYPE (XPLEX)                 :: COARSE_GRID(4)
      TYPE (XPLEX)                 :: VMR_COARSE(4)
      TYPE (XPLEX)                 :: AK_ROW(MAXLEV)
      TYPE (XPLEX)                 :: temp
      
      LOGICAL                :: FOUND_2nd, FOUND_3rd
      
      

      !=================================================================
      ! GET_RTVMR begins here!
      !=================================================================

      ! Initialize and make necessary variables from arguments

         ! If we've found 2nd and 3rd elements of coarse grid
         FOUND_2nd = .FALSE.
         FOUND_3rd = .FALSE.

         ! To make coding cleaner
         LTES = TES(NT)%LTES(1)

         ! Fine pressure grid
         FINE_GRID(1:LTES) = TES(NT)%PRES(1:LTES)

      ! Construct Coarse Pressure grid
      COARSE_GRID(1) = TES(NT)%PRES(1)        ! Bottom level
      COARSE_GRID(4) = TES(NT)%PRES(LTES)     ! Top level
      AK_ROW(1:LTES) = xplx(SUM(TES(NT)%AVG_KERNEL(1:LTES,1:LTES)%r,2),
     &                      SUM(TES(NT)%AVG_KERNEL(1:LTES,1:LTES)%i,2))
      ! Find max of rows of AK below ~50hPa
      MAX_AK         = MAXVAL( AK_ROW(1:35) )
      IF (NT == 600) THEN
         print*,'--------------------------------------------------'
         print*,'NT = ', NT
         print*,'MAX_AK',MAX_AK
         print*,'AK_ROW',AK_ROW
         print*,'--------------------------------------------------'
      ENDIF

      DO L=LTES,1,-1
         ! First pressure level at which sum of rows of AK > 0.4
         IF ( AK_ROW(L) > 0.4 .AND. TES(NT)%PRES(L) > 30.0  .AND.
     &        FOUND_3rd == .FALSE. ) THEN
            COARSE_GRID(3) = TES(NT)%PRES(L)
            FOUND_3rd      = .TRUE.
         ENDIF
         ! Pressure level at which rows of AK are maximum
         IF ( AK_ROW(L) == MAX_AK  .AND.  FOUND_2nd == .FALSE. ) THEN
            COARSE_GRID(2) = TES(NT)%PRES(L)
            FOUND_2nd      = .TRUE.
         ENDIF
      ENDDO


      ! Now that we have fine and coarse grids, make mapping matrix
      M_STAR = MAKE_RTVMR_MAP( NT, LTES, FINE_GRID, COARSE_GRID )

!      !kjw debug
!      IF ( NT == 600 ) THEN
!         print*,'Checking AK_ROW'
!         print*,AK_ROW
!         print*,'Checking COARSE_GRID'
!         print*,COARSE_GRID
!         print*,'Checking M_STAR ... '
!         print*,'SUM of rows of M_STAR(4,LTES)'
!         print*,SUM(M_STAR,2)
!         print*,'Writing Out M_STAR'
!         WRITE(6,546) (L,M_STAR(1,L),M_STAR(2,L),M_STAR(3,L),
!     &                M_STAR(4,L), L=1,MAXLEV)
! 546     FORMAT(i4, 2x, F10.8, 2x, F10.8, 2x, F10.8, 2x, F10.8)
!      ENDIF



      ! Apply mapping matrix to CH4 column
      DO LC=1,4
         temp = 0d0
         DO L=1,LTES
            temp = temp + M_STAR(LC,L) * VMR_IN(L)
         ENDDO
         VMR_COARSE(LC) = temp
      ENDDO

      ! RTVMR value is 2nd element of the coarse array
      RTVMR_OUT = VMR_COARSE(2)


      ! Return to calling program
      END SUBROUTINE GET_RTVMR

!------------------------------------------------------------------------------

      FUNCTION MAKE_RTVMR_MAP( NT, LTES, FINE_GRID, COARSE_GRID )
     &         RESULT( M_STAR )
!
!******************************************************************************
!  Subroutine MAKE_RTVMR_MAP makes matrix to map 67 element TES grid to 4
!  element RTVMR grid. Adapted from from Mark Shephard and Vivienne
!  Payne's retv_make_map_vhp.pro (acquired by kjw from Vivienne).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NT       (INTEGER)   : # of TES observation
!  (2 ) FINE_GRID   (REAL)   : Fine pressure grid from which to map VMR
!  (3 ) COARSE_GRID (REAL)   : Coarse pressure grid onto which we will map VMR
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) M_STAR      (REAL)   : Normalized mapping matrix   x_coarse = M* x_fine
!          M* is pseudo-inverse of M, which maps coarse to fine grid

!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules

      ! Arguments
      !INTEGER, INTENT(IN)   :: LTES
      !TYPE (XPLEX), INTENT(IN)    :: FINE_GRID(MAXLEV)
      !TYPE (XPLEX), INTENT(IN)    :: COARSE_GRID(4)
      !TYPE (XPLEX), INTENT(OUT)   :: M_STAR(4,MAXLEV)
      INTEGER               :: LTES, NT
      TYPE (XPLEX)                :: FINE_GRID(MAXLEV)
      TYPE (XPLEX)                :: COARSE_GRID(4)
      TYPE (XPLEX)                :: M_STAR(4,MAXLEV)
      
      ! Local Variables
      INTEGER               :: L, LC, IND, K
      INTEGER               :: FINE_INDS(4), FINE_INDS_SIX(6)
      TYPE (XPLEX)                :: MAP_TEMP(4,MAXLEV)
      TYPE (XPLEX)                :: MAP_NORM(4,MAXLEV)
      TYPE (XPLEX)                :: sum_map(4)
      TYPE (XPLEX)                :: xdelta_p, xcoeff

      !=================================================================
      ! MAKE_RTVMR_MAP begins here!
      !=================================================================


      ! Initialize and get required values
      MAP_TEMP(:,:) = 0d0
      MAP_NORM(:,:) = 0d0


      ! Find indices of fine grid which match coarse grid
      FINE_INDS(:) = 0d0
      IND = 1
      DO L=1,LTES
         IF ( FINE_GRID(L) == COARSE_GRID(IND) ) THEN
            FINE_INDS(IND) = L
            IND = IND + 1
         ENDIF
      ENDDO

      ! Make 6-element array of indices
      FINE_INDS_SIX(:)   = 0d0
      FINE_INDS_SIX(1)   = FINE_INDS(1)
      FINE_INDS_SIX(2:5) = FINE_INDS(:)
      FINE_INDS_SIX(6)   = FINE_INDS(4)

         !kjw debug
!      IF ( NT == 600 ) THEN 
!         print*,'Checking FINE_GRID'
!         print*,FINE_GRID
!         print*,'Checking FINE_INDS'
!         print*,FINE_INDS
!         print*,'Checking FINE_INDS_SIX'
!         print*,FINE_INDS_SIX
!      ENDIF


         DO L=1,6
         IF ( FINE_INDS_SIX(L) == 0.0 ) THEN
            print*,'kjw debug: indices of fine grid matches to coarse'
            print*,FINE_INDS_SIX
            print*,'  doh, this is f***ed up.  FINE_INDS(L) = 0. L = ',L
            print*,COARSE_GRID
         ENDIF
         IF ( FINE_INDS_SIX(L) > 67.0 ) THEN
            print*,'kjw debug: indices of fine grid matches to coarse'
            print*,'  doh, this is f***ed up.  FINE_INDS(L) >67. L = ',L
            print*,COARSE_GRID
         ENDIF
         ENDDO

      ! Populate mapping matrix
      K = 1
      DO LC=1,4
         DO L=FINE_INDS_SIX(K),FINE_INDS_SIX(K+2)

            ! Bottom of profile is set a constant perturbation
            IF ( FINE_GRID(L) > COARSE_GRID(LC) .AND. LC == 1 ) THEN
               MAP_TEMP(LC,L) = 1.0d0
            ENDIF

            ! Bottom side of profile
            IF ( LC /= 1 ) THEN
               IF ( FINE_GRID(L) >= COARSE_GRID(LC) .AND. 
     &              FINE_GRID(L) <= COARSE_GRID(LC-1) ) THEN
                  xdelta_p = LOG(COARSE_GRID(LC-1))-LOG(COARSE_GRID(LC))
                  xcoeff   = 1d0 - ( LOG(FINE_GRID(L)) -
     &                               LOG(COARSE_GRID(LC)) ) / xdelta_p
                  MAP_TEMP(LC,L) = xcoeff
               ENDIF
            ENDIF

            ! Top side of profile
            IF ( LC /= 4 ) THEN
               IF ( FINE_GRID(L) <= COARSE_GRID(LC) .AND. 
     &              FINE_GRID(L) >= COARSE_GRID(LC+1) ) THEN
                  xdelta_p = LOG(COARSE_GRID(LC))-LOG(COARSE_GRID(LC+1))
                  xcoeff   = 1d0 - ( -LOG(FINE_GRID(L)) + 
     &                                LOG(COARSE_GRID(LC)) ) / xdelta_p
                  MAP_TEMP(LC,L) = xcoeff
               ENDIF
            ENDIF

            ! Top of profile is set a constant perturbation
            IF ( FINE_GRID(L) < COARSE_GRID(LC) .AND. LC == 4 ) THEN
               MAP_TEMP(LC,L) = 1.0d0
            ENDIF

         ENDDO

         ! Increment Indices between which to fill
         K = K + 1
      ENDDO

!      !kjw debug
!      IF ( NT == 600 ) THEN
!         print*,'Checking M_STAR ... '
!         print*,'SUM of rows of MAP_TEMP(4,LTES)'
!         print*,SUM(MAP_TEMP,2)
!         print*,'Writing Out MAP_TEMP'
!         WRITE(6,547) (L,MAP_TEMP(1,L),MAP_TEMP(2,L),MAP_TEMP(3,L),
!     &                MAP_TEMP(4,L), L=1,MAXLEV)
! 547     FORMAT(i4, 2x, F10.8, 2x, F10.8, 2x, F10.8, 2x, F10.8)
!      ENDIF

      ! Normalize Mapping Matrix
      sum_map(:) = 0d0
      sum_map(:) = xplx(SUM( MAP_TEMP%r, 2 ),SUM( MAP_TEMP%i, 2 ))
      sum_map(:) = 0d0
      DO LC=1,4
         sum_map(LC)    = SUM( MAP_TEMP(LC,:) )
         IF (NT .EQ. 600) THEN 
            !print*,'Sum map',LC
            !print*,sum_map(LC)
         ENDIF
         MAP_NORM(LC,:) = MAP_TEMP(LC,:) / sum_map(LC)
      ENDDO

!      !kjw debug
!      IF ( NT == 600 ) THEN
!         print*,'Checking M_STAR ... '
!         print*,'SUM of rows of MAP_NORM(4,LTES)'
!         print*,SUM(MAP_NORM,2)
!         print*,'Writing Out MAP_NORM'
!         WRITE(6,547) (L,MAP_NORM(1,L),MAP_NORM(2,L),MAP_NORM(3,L),
!     &                MAP_NORM(4,L), L=1,MAXLEV)
!      ENDIF


      ! Assign Map to output variable
      M_STAR = MAP_NORM


      ! Return to calling program
      END FUNCTION MAKE_RTVMR_MAP


!------------------------------------------------------------------------------

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
!  (3 ) TIME_FRAC (REAL) : Vector of times (frac-of-day) for the TES retrievals
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


      print*, ' GET_NT_RANGE for ', HHMMSS
      print*, ' NTSAVE ', NTSAVE
      print*, ' NTES   ', NTES
   
      CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )


      ! Convert HH from hour to fraction of day 
      GC_HH_FRAC = XPLX(HH) / 24d0 
 
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
               print*, ' still looking ', NTEST 
                  
            ENDIF 
                 
         ENDDO 
 
      ELSE

         CALL ERROR_STOP('problem', 'GET_NT_RANGE' ) 

      ENDIF 

      ! Return to calling program
      END SUBROUTINE GET_NT_RANGE

!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP,
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )
     *         RESULT      ( HINTERPZ )
!
!******************************************************************************
!  Function GET_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!     
!  NOTES:
!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod. 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

      ! Arguments
      INTEGER            :: LGC_TOP, LTM_TOP
      TYPE (XPLEX)             :: GC_PRESC(LGC_TOP)
      TYPE (XPLEX)             :: TM_PRESC(LTM_TOP) 
      TYPE (XPLEX)             :: GC_SURFP
      TYPE (XPLEX)             :: TM_SURFP
 
      ! Return value 
      TYPE (XPLEX)             :: HINTERPZ(LGC_TOP, LTM_TOP)

      ! Local variables 
      INTEGER  :: LGC, LTM
      TYPE (XPLEX)   :: DIFF, DELTA_SURFP
      TYPE (XPLEX)   :: LOW, HI

      !=================================================================
      ! GET_HINTERPZ_2 begins here!
      !=================================================================

      HINTERPZ(:,:) = 0D0 
  
!      ! Rescale GC grid according to TM surface pressure
!!         p1_A =     (a1 + b1 (ps_A - PTOP))
!!         p2_A =     (a2 + b2 (ps_A - PTOP))
!!         p1_B =     (a + b (ps_B - PTOP))
!!         p2_B =    *(a + b (ps_B - PTOP))
!!         pc_A = 0.5(a1+a2 +(b1+b2)*(ps_A - PTOP))
!!         pc_B = 0.5(a1+a2 +(b1+b2)*(ps_B - PTOP))
!!         pc_B - pc_A = 0.5(b1_b2)(ps_B-ps_A)
!!         pc_B = 0.5(b1_b2)(ps_B-ps_A) + pc_A
!      DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
!
!      DO LGC = 1, LGC_TOP
!         GC_PRESC(LGC) = ( GET_BP(LGC) + GET_BP(LGC+1))
!     &               * DELTA_SURFP + GC_PRESC(LGC)
!         IF (GC_PRESC(LGC) < 0) THEN 
!            CALL ERROR_STOP( 'highly unlikey', 
!     &                       'read_sciano2_mod.f')
!         ENDIF 
!
!      ENDDO 
      

      ! Loop over each pressure level of TM grid
      DO LTM = 1, LTM_TOP
 
         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, LGC_TOP - 1

            LOW = GC_PRESC(LGC+1)
            HI  = GC_PRESC(LGC)
            IF (LGC == 0) HI = TM_SURFP  !kjw. this line is useless

            ! Linearly interpolate value on the LTM grid 
            IF ( TM_PRESC(LTM) <= HI .and. 
     &           TM_PRESC(LTM)  > LOW) THEN 

               DIFF                = HI - LOW  
               HINTERPZ(LGC+1,LTM) = ( HI - TM_PRESC(LTM)  ) / DIFF
               HINTERPZ(LGC  ,LTM) = ( TM_PRESC(LTM) - LOW ) / DIFF


            ENDIF 
 
            ! dkh debug
            !print*, 'LGC,LTM,HINT', LGC, LTM, HINTERPZ(LGC,LTM)

          ENDDO
       ENDDO

       ! Correct for case where TES pressure is higher than the
       ! highest GC pressure.  In this case, just 1:1 map. 
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             HINTERPZ(:,LTM)   = 0D0
             HINTERPZ(LTM,LTM) = 1D0
          ENDIF
       ENDDO

      ! Return to calling program
      END FUNCTION GET_INTMAP

!!------------------------------------------------------------------------------

!!------------------------------------------------------------------------------
      FUNCTION GET_IJ_2x25( LON, LAT ) RESULT ( IIJJ )

!
!******************************************************************************
!  Subroutine GET_IJ_2x25 returns I and J index from the 2 x 2.5 grid for a 
!  LON, LAT coord. (dkh, 11/08/09) 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LON (TYPE (XPLEX)) : Longitude                          [degrees]
!  (2 ) LAT (TYPE (XPLEX)) : Latitude                           [degrees]
!     
!  Function result
!  ============================================================================
!  (1 ) IIJJ(1) (INTEGER) : Long index                    [none]
!  (2 ) IIJJ(2) (INTEGER) : Lati index                    [none]
!     
!  NOTES:
!
!******************************************************************************
!     
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP

      ! Arguments
      TYPE (XPLEX)    :: LAT, LON
      
      ! Return
      INTEGER :: I, J, IIJJ(2)
      
      ! Local variables 
      TYPE (XPLEX)              :: TLON, TLAT, DLON, DLAT
      TYPE (XPLEX),  PARAMETER  :: DISIZE = xplex(2.5d0,0d0)
      TYPE (XPLEX),  PARAMETER  :: DJSIZE = xplex(2.0d0,0d0)
      INTEGER, PARAMETER  :: IIMAX  = 144
      INTEGER, PARAMETER  :: JJMAX  = 91
      
      
      !=================================================================
      ! GET_IJ_2x25 begins here!
      !=================================================================

      TLON = 180d0 + LON + DISIZE
      TLAT =  90d0 + LAT + DJSIZE
      
      I = TLON / DISIZE
      J = TLAT / DJSIZE

      
      IF ( TLON / DISIZE -XPLX(I)  >= 0.5d0 ) THEN
         I = I + 1
      ENDIF
      
      IF ( TLAT / DJSIZE -XPLX(J)  >= 0.5d0 ) THEN
         J = J + 1
      ENDIF

      
      ! Longitude wraps around
      !IF ( I == 73 ) I = 1 
      IF ( I == ( IIMAX + 1 ) ) I = 1
      
      ! Check for impossible values 
      IF ( I > IIMAX .or. J > JJMAX .or. 
     &     I < 1     .or. J < 1          ) THEN
         CALL ERROR_STOP('Error finding grid box', 'GET_IJ_2x25')
      ENDIF
      
      IIJJ(1) = I
      IIJJ(2) = J
      
      ! Return to calling program
      END FUNCTION GET_IJ_2x25

!!-----------------------------------------------------------------------------
!      SUBROUTINE INIT_TES_CH4
!!
!!*****************************************************************************
!!  Subroutine INIT_TES_CH4 deallocates all module arrays.  (dkh, 02/15/09) 
!!        
!!  NOTES:
!!
!!******************************************************************************
!!     
!      USE ERROR_MOD,  ONLY : ALLOC_ERR
!
!#     include "CMN_SIZE"   ! IIPAR, JJPAR      
!       
!      ! Local variables
!      INTEGER :: AS
!
!      !=================================================================      
!      ! INIT_TES_CH4 begins here
!      !================================================================= 
!
!      ! dkh debug
!      print*, ' INIT_TES_CH4'
!
!      ALLOCATE( CH4_SAVE( LLPAR, MAXTES ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_SAVE' ) 
!      CH4_SAVE = 0d0
!
!
!      TES(	1	)%NYMD 	=	20050704
!      TES(	2	)%NYMD 	=	20050704
!      TES(	3	)%NYMD 	=	20050704
!      TES(	4	)%NYMD 	=	20050704
!      TES(	5	)%NYMD 	=	20050704
!      TES(	6	)%NYMD 	=	20050704
!      TES(	7	)%NYMD 	=	20050704
!      TES(	8	)%NYMD 	=	20050704
!      TES(	9	)%NYMD 	=	20050705
!      TES(	10	)%NYMD 	=	20050705
!      TES(	11	)%NYMD 	=	20050705
!      TES(	12	)%NYMD 	=	20050705
!      TES(	13	)%NYMD 	=	20050705
!      TES(	14	)%NYMD 	=	20050705
!      TES(	15	)%NYMD 	=	20050705
!      TES(	16	)%NYMD 	=	20050705
!      TES(	17	)%NYMD 	=	20050705
!      TES(	18	)%NYMD 	=	20050710
!      TES(	19	)%NYMD 	=	20050710
!      TES(	20	)%NYMD 	=	20050710
!      TES(	21	)%NYMD 	=	20050710
!      TES(	22	)%NYMD 	=	20050710
!      TES(	23	)%NYMD 	=	20050710
!      TES(	24	)%NYMD 	=	20050710
!      TES(	25	)%NYMD 	=	20050710
!      TES(	26	)%NYMD 	=	20050710
!      TES(	27	)%NYMD 	=	20050711
!      TES(	28	)%NYMD 	=	20050711
!      TES(	29	)%NYMD 	=	20050711
!      TES(	30	)%NYMD 	=	20050711
!      TES(	31	)%NYMD 	=	20050712
!      TES(	32	)%NYMD 	=	20050712
!      TES(	33	)%NYMD 	=	20050712
!      TES(	34	)%NYMD 	=	20050712
!      TES(	35	)%NYMD 	=	20050712
!      TES(	36	)%NYMD 	=	20050712
!      TES(	37	)%NYMD 	=	20050712
!      TES(	38	)%NYMD 	=	20050712
!      TES(	39	)%NYMD 	=	20050713
!      TES(	40	)%NYMD 	=	20050713
!      TES(	41	)%NYMD 	=	20050713
!      TES(	42	)%NYMD 	=	20050713
!      TES(	43	)%NYMD 	=	20050713
!      TES(	44	)%NYMD 	=	20050713
!      TES(	45	)%NYMD 	=	20050713
!      TES(	46	)%NYMD 	=	20050713
!      TES(	47	)%NYMD 	=	20050713
!      TES(	48	)%NYMD 	=	20050714
!      TES(	49	)%NYMD 	=	20050714
!      TES(	50	)%NYMD 	=	20050714
!      TES(	51	)%NYMD 	=	20050714
!      TES(	52	)%NYMD 	=	20050714
!      TES(	53	)%NYMD 	=	20050714
!      TES(	54	)%NYMD 	=	20050714
!      TES(	55	)%NYMD 	=	20050714
!      TES(	56	)%NYMD 	=	20050715
!      TES(	57	)%NYMD 	=	20050715
!      TES(	58	)%NYMD 	=	20050715
!      TES(	59	)%NYMD 	=	20050715
!      TES(	60	)%NYMD 	=	20050715
!      TES(	61	)%NYMD 	=	20050715
!      TES(	62	)%NYMD 	=	20050715
!      TES(	63	)%NYMD 	=	20050715
!      TES(	64	)%NYMD 	=	20050715
!      TES(	65	)%NYMD 	=	20050716
!      TES(	66	)%NYMD 	=	20050717
!      TES(	67	)%NYMD 	=	20050717
!      TES(	68	)%NYMD 	=	20050717
!      TES(	69	)%NYMD 	=	20050717
!      TES(	70	)%NYMD 	=	20050717
!      TES(	71	)%NYMD 	=	20050717
!      TES(	72	)%NYMD 	=	20050717
!      TES(	73	)%NYMD 	=	20050717
!      TES(	74	)%NYMD 	=	20050717
!      TES(	75	)%NYMD 	=	20050718
!      TES(	76	)%NYMD 	=	20050718
!      TES(	77	)%NYMD 	=	20050718
!      TES(	78	)%NYMD 	=	20050718
!      TES(	79	)%NYMD 	=	20050719
!      TES(	80	)%NYMD 	=	20050719
!      TES(	81	)%NYMD 	=	20050719
!      TES(	82	)%NYMD 	=	20050719
!      TES(	83	)%NYMD 	=	20050719
!      TES(	84	)%NYMD 	=	20050719
!      TES(	85	)%NYMD 	=	20050719
!      TES(	86	)%NYMD 	=	20050719
!      TES(	87	)%NYMD 	=	20050719
!      
!      TES(	1	)%NHMS 	=	202000
!      TES(	2	)%NHMS 	=	202100
!      TES(	3	)%NHMS 	=	202100
!      TES(	4	)%NHMS 	=	202100
!      TES(	5	)%NHMS 	=	202200
!      TES(	6	)%NHMS 	=	202300
!      TES(	7	)%NHMS 	=	202300
!      TES(	8	)%NHMS 	=	202400
!      TES(	9	)%NHMS 	=	082100
!      TES(	10	)%NHMS 	=	082100
!      TES(	11	)%NHMS 	=	082200
!      TES(	12	)%NHMS 	=	082200
!      TES(	13	)%NHMS 	=	082300
!      TES(	14	)%NHMS 	=	082300
!      TES(	15	)%NHMS 	=	082400
!      TES(	16	)%NHMS 	=	082400
!      TES(	17	)%NHMS 	=	082500
!      TES(	18	)%NHMS 	=	194300
!      TES(	19	)%NHMS 	=	194300
!      TES(	20	)%NHMS 	=	194400
!      TES(	21	)%NHMS 	=	194400
!      TES(	22	)%NHMS 	=	194500
!      TES(	23	)%NHMS 	=	194500
!      TES(	24	)%NHMS 	=	194600
!      TES(	25	)%NHMS 	=	194600
!      TES(	26	)%NHMS 	=	194700
!      TES(	27	)%NHMS 	=	092300
!      TES(	28	)%NHMS 	=	092300
!      TES(	29	)%NHMS 	=	092400
!      TES(	30	)%NHMS 	=	092400
!      TES(	31	)%NHMS 	=	193000
!      TES(	32	)%NHMS 	=	193100
!      TES(	33	)%NHMS 	=	193100
!      TES(	34	)%NHMS 	=	193200
!      TES(	35	)%NHMS 	=	193300
!      TES(	36	)%NHMS 	=	193300
!      TES(	37	)%NHMS 	=	193400
!      TES(	38	)CH4%NHMS 	=	193400
!      TES(	39	)%NHMS 	=	091000
!      TES(	40	)%NHMS 	=	091100
!      TES(	41	)%NHMS 	=	091100
!      TES(	42	)%NHMS 	=	091200
!      TES(	43	)%NHMS 	=	091200
!      TES(	44	)%NHMS 	=	091200
!      TES(	45	)%NHMS 	=	091300
!      TES(	46	)%NHMS 	=	091300
!      TES(	47	)%NHMS 	=	091400
!      TES(	48	)%NHMS 	=	191900
!      TES(	49	)%NHMS 	=	191900
!      TES(	50	)%NHMS 	=	191900
!      TES(	51	)%NHMS 	=	192000
!      TES(	52	)%NHMS 	=	192000
!      TES(	53	)%NHMS 	=	192100
!      TES(	54	)%NHMS 	=	192100
!      TES(	55	)%NHMS 	=	192200
!      TES(	56	)%NHMS 	=	085800
!      TES(	57	)%NHMS 	=	085800
!      TES(	58	)%NHMS 	=	085900
!      TES(	59	)%NHMS 	=	085900
!      TES(	60	)%NHMS 	=	090000
!      TES(	61	)%NHMS 	=	090000
!      TES(	62	)%NHMS 	=	090100
!      TES(	63	)%NHMS 	=	090100
!      TES(	64	)%NHMS 	=	090100
!      TES(	65	)%NHMS 	=	190900
!      TES(	66	)%NHMS 	=	084500
!      TES(	67	)%NHMS 	=	084600
!      TES(	68	)%NHMS 	=	084600
!      TES(	69	)%NHMS 	=	084700
!      TES(	70	)%NHMS 	=	084700
!      TES(	71	)%NHMS 	=	084800
!      TES(	72	)%NHMS 	=	084800
!      TES(	73	)%NHMS 	=	084900
!      TES(	74	)%NHMS 	=	084900
!      TES(	75	)%NHMS 	=	203200
!      TES(	76	)%NHMS 	=	203300
!      TES(	77	)%NHMS 	=	203300
!      TES(	78	)%NHMS 	=	203400
!      TES(	79	)%NHMS 	=	083300
!      TES(	80	)%NHMS 	=	083400
!      TES(	81	)%NHMS 	=	083400
!      TES(	82	)%NHMS 	=	083500
!      TES(	83	)%NHMS 	=	083500
!      TES(	84	)%NHMS 	=	083500
!      TES(	85	)%NHMS 	=	083600
!      TES(	86	)%NHMS 	=	083600
!      TES(	87	)%NHMS 	=	083700
!     
!      TES(	1	)%LAT 	=	31.29
!      TES(	2	)%LAT 	=	33
!      TES(	3	)%LAT 	=	34.64
!      TES(	4	)%LAT 	=	36.2
!      TES(	5	)%LAT 	=	37.91
!      TES(	6	)%LAT 	=	41.1
!      TES(	7	)%LAT 	=	42.8
!      TES(	8	)%LAT 	=	44.43
!      TES(	9	)%LAT 	=	43.54
!      TES(	10	)%LAT 	=	41.84
!      TES(	11	)%LAT 	=	40.2
!      TES(	12	)%LAT 	=	38.65
!      TES(	13	)%LAT 	=	36.94
!      TES(	14	)%LAT 	=	35.3
!      TES(	15	)%LAT 	=	33.74
!      TES(	16	)%LAT 	=	32.03
!      TES(	17	)%LAT 	=	30.39 
!      TES(	18	)%LAT 	=	31.28
!      TES(	19	)%LAT 	=	32.99
!      TES(	20	)%LAT 	=	34.63
!      TES(	21	)%LAT 	=	36.19
!      TES(	22	)%LAT 	=	37.9
!      TES(	23	)%LAT 	=	39.53
!      TES(	24	)%LAT 	=	41.09
!      TES(	25	)%LAT 	=	42.8
!      TES(	26	)%LAT 	=	44.42
!      TES(	27	)%LAT 	=	43.55
!      TES(	28	)%LAT 	=	41.85
!      TES(	29	)%LAT 	=	40.22
!      TES(	30	)%LAT 	=	38.66
!      TES(	31	)%LAT 	=	31.28
!      TES(	32	)%LAT 	=	32.99
!      TES(	33	)%LAT 	=	34.63
!      TES(	34	)%LAT 	=	36.19
!      TES(	35	)%LAT 	=	39.53
!      TES(	36	)%LAT 	=	41.09
!      TES(	37	)%LAT 	=	42.79
!      TES(	38	)%LAT 	=	44.42
!      TES(	39	)%LAT 	=	43.55
!      TES(	40	)%LAT 	=	41.85
!      TES(	41	)%LAT 	=	40.22
!      TES(	42	)%LAT 	=	38.66
!      TES(	43	)%LAT 	=	36.96
!      TES(	44	)%LAT 	=	35.32
!      TES(	45	)%LAT 	=	33.76
!      TES(	46	)%LAT 	=	32.04
!      TES(	47	)%LAT 	=	30.4
!      TES(	48	)%LAT 	=	32.99
!      TES(	49	)%LAT 	=	34.63
!      TES(	50	)%LAT 	=	36.2
!      TES(	51	)%LAT 	=	37.9
!      TES(	52	)%LAT 	=	39.54
!      TES(	53	)%LAT 	=	41.1
!      TES(	54	)%LAT 	=	42.8
!      TES(	55	)%LAT 	=	44.42
!      TES(	56	)%LAT 	=	43.55
!      TES(	57	)%LAT 	=	41.85
!      TES(	58	)%LAT 	=	40.22
!      TES(	59	)%LAT 	=	38.66
!      TES(	60	)%LAT 	=	36.95
!      TES(	61	)%LAT 	=	35.31
!      TES(	62	)%LAT 	=	33.75
!      TES(	63	)%LAT 	=	32.04
!      TES(	64	)%LAT 	=	30.4
!      TES(	65	)%LAT 	=	44.4
!      TES(	66	)%LAT 	=	43.59
!      TES(	67	)%LAT 	=	41.89
!      TES(	68	)%LAT 	=	40.26
!      TES(	69	)%LAT 	=	38.7
!      TES(	70	)%LAT 	=	37
!      TES(	71	)%LAT 	=	35.36
!      TES(	72	)%LAT 	=	33.8
!      TES(	73	)%LAT 	=	32.09
!      TES(	74	)%LAT 	=	30.45
!      TES(	75	)%LAT 	=	31.27
!      TES(	76	)%LAT 	=	32.98
!      TES(	77	)%LAT 	=	34.62
!      TES(	78	)%LAT 	=	36.18
!      TES(	79	)%LAT 	=	43.58
!      TES(	80	)%LAT 	=	41.88
!      TES(	81	)%LAT 	=	40.25
!      TES(	82	)%LAT 	=	38.69
!      TES(	83	)%LAT 	=	36.98
!      TES(	84	)%LAT 	=	35.34
!      TES(	85	)%LAT 	=	33.78
!      TES(	86	)%LAT 	=	32.07
!      TES(	87	)%LAT 	=	30.43
!      
!      TES(	1	)%LON	=	-105.13
!      TES(	2	)%LON	=	-105.6
!      TES(	3	)%LON	=	-106.05
!      TES(	4	)%LON	=	-106.5
!      TES(	5	)%LON	=	-107
!      TES(	6	)%LON	=	-108
!      TES(	7	)%LON	=	-108.57
!      TES(	8	)%LON	=	-109.13
!      TES(	9	)%LON	=	-92.52
!      TES(	10	)%LON	=	-93.09
!      TES(	11	)%LON	=	-93.62
!      TES(	12	)%LON	=	-94.11
!      TES(	13	)%LON	=	-94.62
!      TES(	14	)%LON	=	-95.09
!      TES(	15	)%LON	=	-95.53
!      TES(	16	)%LON	=	-96
!      TES(	17	)%LON	=	-96.44
!      TES(	18	)%LON	=	-95.84
!      TES(	19	)%LON	=	-96.3
!      TES(	20	)%LON	=	-96.76
!      TES(	21	)%LON	=	-97.2
!      TES(	22	)%LON	=	-97.71
!      TES(	23	)%LON	=	-98.21
!      TES(	24	)%LON	=	-98.71
!      TES(	25	)%LON	=	-99.27
!      TES(	26	)%LON	=	-99.83
!      TES(	27	)%LON	=	-107.94
!      TES(	28	)%LON	=	-108.51
!      TES(	29	)%LON	=	-109.04
!      TES(	30	)%LON	=	-109.53
!      TES(	31	)%LON	=	-92.74
!      TES(	32	)%LON	=	-93.2
!      TES(	33	)%LON	=	-93.66
!      TES(	34	)%LON	=	-94.11
!      TES(	35	)%LON	=	-95.11
!      TES(	36	)%LON	=	-95.61
!      TES(	37	)%LON	=	-96.17
!      TES(	38	)%LON	=	-96.73
!      TES(	39	)%LON	=	-104.84
!      TES(	40	)%LON	=	-105.41
!      TES(	41	)%LON	=	-105.94
!      TES(	42	)%LON	=	-106.43
!      TES(	43	)%LON	=	-106.94
!      TES(	44	)%LON	=	-107.42
!      TES(	45	)%LON	=	-107.86
!      TES(	46	)%LON	=	-108.33
!      TES(	47	)%LON	=	-108.76
!      TES(	48	)%LON	=	-90.1
!      TES(	49	)%LON	=	-90.56
!      TES(	50	)%LON	=	-91.01
!      TES(	51	)%LON	=	-91.51
!      TES(	52	)%LON	=	-92.01
!      TES(	53	)%LON	=	-92.51
!      TES(	54	)%LON	=	-93.07
!      TES(	55	)%LON	=	-93.64
!      TES(	56	)%LON	=	-101.74
!      TES(	57	)%LON	=	-102.32
!      TES(	58	)%LON	=	-102.84
!      TES(	59	)%LON	=	-103.33
!      TES(	60	)%LON	=	-103.84
!      TES(	61	)%LON	=	-104.32
!      TES(	62	)%LON	=	-104.76
!      TES(	63	)%LON	=	-105.23
!      TES(	64	)%LON	=	-105.67
!      TES(	65	)%LON	=	-90.54
!      TES(	66	)%LON	=	-98.64
!      TES(	67	)%LON	=	-99.22
!      TES(	68	)%LON	=	-99.75
!      TES(	69	)%LON	=	-100.23
!      TES(	70	)%LON	=	-100.75
!      TES(	71	)%LON	=	-101.22
!      TES(	72	)%LON	=	-101.67
!      TES(	73	)%LON	=	-102.13
!      TES(	74	)%LON	=	-102.57
!      TES(	75	)%LON	=	-108.19
!      TES(	76	)%LON	=	-108.65
!      TES(	77	)%LON	=	-109.11
!      TES(	78	)%LON	=	-109.55
!      TES(	79	)%LON	=	-95.57
!      TES(	80	)%LON	=	-96.14
!      TES(	81	)%LON	=	-96.67
!      TES(	82	)%LON	=	-97.16
!      TES(	83	)%LON	=	-97.67
!      TES(	84	)%LON	=	-98.15
!      TES(	85	)%LON	=	-98.59
!      TES(	86	)%LON	=	-99.06
!      TES(	87	)%LON	=	-99.49
!     
!      TES(	1	)%FILENAME = TRIM('retv_vars.02945_0457_002.cdf')
!      TES(	2	)%FILENAME = TRIM('retv_vars.02945_0457_003.cdf')
!      TES(	3	)%FILENAME = TRIM('retv_vars.02945_0457_004.cdf')
!      TES(	4	)%FILENAME = TRIM('retv_vars.02945_0458_002.cdf')
!      TES(	5	)%FILENAME = TRIM('retv_vars.02945_0458_003.cdf')
!      TES(	6	)%FILENAME = TRIM('retv_vars.02945_0459_002.cdf')
!      TES(	7	)%FILENAME = TRIM('retv_vars.02945_0459_003.cdf')
!      TES(	8	)%FILENAME = TRIM('retv_vars.02945_0459_004.cdf')
!      TES(	9	)%FILENAME = TRIM('retv_vars.02945_0982_002.cdf')
!      TES(	10	)%FILENAME = TRIM('retv_vars.02945_0982_003.cdf')
!      TES(	11	)%FILENAME = TRIM('retv_vars.02945_0982_004.cdf')
!      TES(	12	)%FILENAME = TRIM('retv_vars.02945_0983_002.cdf')
!      TES(	13	)%FILENAME = TRIM('retv_vars.02945_0983_003.cdf')
!      TES(	14	)%FILENAME = TRIM('retv_vars.02945_0983_004.cdf')
!      TES(	15	)%FILENAME = TRIM('retv_vars.02945_0984_002.cdf')
!      TES(	16	)%FILENAME = TRIM('retv_vars.02945_0984_003.cdf')
!      TES(	17	)%FILENAME = TRIM('retv_vars.02945_0984_004.cdf')
!      TES(	18	)%FILENAME = TRIM('retv_vars.02956_0457_002.cdf')
!      TES(	19	)%FILENAME = TRIM('retv_vars.02956_0457_003.cdf')
!      TES(	20	)%FILENAME = TRIM('retv_vars.02956_0457_004.cdf')
!      TES(	21	)%FILENAME = TRIM('retv_vars.02956_0458_002.cdf')
!      TES(	22	)%FILENAME = TRIM('retv_vars.02956_0458_003.cdf')
!      TES(	23	)%FILENAME = TRIM('retv_vars.02956_0458_004.cdf')
!      TES(	24	)%FILENAME = TRIM('retv_vars.02956_0459_002.cdf')
!      TES(	25	)%FILENAME = TRIM('retv_vars.02956_0459_003.cdf')
!      TES(	26	)%FILENAME = TRIM('retv_vars.02956_0459_004.cdf')
!      TES(	27	)%FILENAME = TRIM('retv_vars.02956_1054_002.cdf')
!      TES(	28	)%FILENAME = TRIM('retv_vars.02956_1054_003.cdf')
!      TES(	29	)%FILENAME = TRIM('retv_vars.02956_1054_004.cdf')
!      TES(	30	)%FILENAME = TRIM('retv_vars.02956_1055_002.cdf')
!      TES(	31	)%FILENAME = TRIM('retv_vars.02960_0457_002.cdf')
!      TES(	32	)%FILENAME = TRIM('retv_vars.02960_0457_003.cdf')
!      TES(	33	)%FILENAME = TRIM('retv_vars.02960_0457_004.cdf')
!      TES(	34	)%FILENAME = TRIM('retv_vars.02960_0458_002.cdf')
!      TES(	35	)%FILENAME = TRIM('retv_vars.02960_0458_004.cdf')
!      TES(	36	)%FILENAME = TRIM('retv_vars.02960_0459_002.cdf')
!      TES(	37	)%FILENAME = TRIM('retv_vars.02960_0459_003.cdf')
!      TES(	38	)%FILENAME = TRIM('retv_vars.02960_0459_004.cdf')
!      TES(	39	)%FILENAME = TRIM('retv_vars.02960_1054_002.cdf')
!      TES(	40	)%FILENAME = TRIM('retv_vars.02960_1054_003.cdf')
!      TES(	41	)%FILENAME = TRIM('retv_vars.02960_1054_004.cdf')
!      TES(	42	)%FILENAME = TRIM('retv_vars.02960_1055_002.cdf')
!      TES(	43	)%FILENAME = TRIM('retv_vars.02960_1055_003.cdf')
!      TES(	44	)%FILENAME = TRIM('retv_vars.02960_1055_004.cdf')
!      TES(	45	)%FILENAME = TRIM('retv_vars.02960_1056_002.cdf')
!      TES(	46	)%FILENAME = TRIM('retv_vars.02960_1056_003.cdf')
!      TES(	47	)%FILENAME = TRIM('retv_vars.02960_1056_004.cdf')
!      TES(	48	)%FILENAME = TRIM('retv_vars.02963_0457_003.cdf')
!      TES(	49	)%FILENAME = TRIM('retv_vars.02963_0457_004.cdf')
!      TES(	50	)%FILENAME = TRIM('retv_vars.02963_0458_002.cdf')
!      TES(	51	)%FILENAME = TRIM('retv_vars.02963_0458_003.cdf')
!      TES(	52	)%FILENAME = TRIM('retv_vars.02963_0458_004.cdf')
!      TES(	53	)%FILENAME = TRIM('retv_vars.02963_0459_002.cdf')
!      TES(	54	)%FILENAME = TRIM('retv_vars.02963_0459_003.cdf')
!      TES(	55	)%FILENAME = TRIM('retv_vars.02963_0459_004.cdf')
!      TES(	56	)%FILENAME = TRIM('retv_vars.02963_1054_002.cdf')
!      TES(	57	)%FILENAME = TRIM('retv_vars.02963_1054_003.cdf')
!      TES(	58	)%FILENAME = TRIM('retv_vars.02963_1054_004.cdf')
!      TES(	59	)%FILENAME = TRIM('retv_vars.02963_1055_002.cdf')
!      TES(	60	)%FILENAME = TRIM('retv_vars.02963_1055_003.cdf')
!      TES(	61	)%FILENAME = TRIM('retv_vars.02963_1055_004.cdf')
!      TES(	62	)%FILENAME = TRIM('retv_vars.02963_1056_002.cdf')
!      TES(	63	)%FILENAME = TRIM('retv_vars.02963_1056_003.cdf')
!      TES(	64	)%FILENAME = TRIM('retv_vars.02963_1056_004.cdf')
!      TES(	65	)%FILENAME = TRIM('retv_vars.02967_0459_004.cdf')
!      TES(	66	)%FILENAME = TRIM('retv_vars.02967_1054_002.cdf')
!      TES(	67	)%FILENAME = TRIM('retv_vars.02967_1054_003.cdf')
!      TES(	68	)%FILENAME = TRIM('retv_vars.02967_1054_004.cdf')
!      TES(	69	)%FILENAME = TRIM('retv_vars.02967_1055_002.cdf')
!      TES(	70	)%FILENAME = TRIM('retv_vars.02967_1055_003.cdf')
!      TES(	71	)%FILENAME = TRIM('retv_vars.02967_1055_004.cdf')
!      TES(	72	)%FILENAME = TRIM('retv_vars.02967_1056_002.cdf')
!      TES(	73	)%FILENAME = TRIM('retv_vars.02967_1056_003.cdf')
!      TES(	74	)%FILENAME = TRIM('retv_vars.02967_1056_004.cdf')
!      TES(	75	)%FILENAME = TRIM('retv_vars.02971_0457_002.cdf')
!      TES(	76	)%FILENAME = TRIM('retv_vars.02971_0457_003.cdf')
!      TES(	77	)%FILENAME = TRIM('retv_vars.02971_0457_004.cdf')
!      TES(	78	)%FILENAME = TRIM('retv_vars.02971_0458_002.cdf')
!      TES(	79	)%FILENAME = TRIM('retv_vars.02971_0982_002.cdf')
!      TES(	80	)%FILENAME = TRIM('retv_vars.02971_0982_003.cdf')
!      TES(	81	)%FILENAME = TRIM('retv_vars.02971_0982_004.cdf')
!      TES(	82	)%FILENAME = TRIM('retv_vars.02971_0983_002.cdf')
!      TES(	83	)%FILENAME = TRIM('retv_vars.02971_0983_003.cdf')
!      TES(	84	)%FILENAME = TRIM('retv_vars.02971_0983_004.cdf')
!      TES(	85	)%FILENAME = TRIM('retv_vars.02971_0984_002.cdf')
!      TES(	86	)%FILENAME = TRIM('retv_vars.02971_0984_003.cdf')
!      TES(	87	)%FILENAME = TRIM('retv_vars.02971_0984_004.cdf')
!
!      ! Return to calling program 
!      END SUBROUTINE INIT_TES_CH4
!!------------------------------------------------------------------------------
!
!      SUBROUTINE CLEANUP_TES_CH4
!!
!!*****************************************************************************
!!  Subroutine CLEANUP_TES_CH4 deallocates all module arrays. (dkh, 02/15/09) 
!!        
!!  NOTES:
!!
!!******************************************************************************
!!     
!
!      IF ( ALLOCATED( CH4_SAVE ) )      DEALLOCATE( CH4_SAVE )
!
!
!      ! Return to calling program 
!      END SUBROUTINE CLEANUP_TES_CH4
!!------------------------------------------------------------------------------


!------------------------------------------------------------------------------

!      SUBROUTINE GET_GC_PSEUDO_OBS( NTES )
!!
!!******************************************************************************
!!  Subroutine GET_GC_PSEUDO_OBS replaces TES observatins in TES%CH4 with 
!!  pseudo-observations from a GEOS-Chem run with scaling factors = 1.
!!  The GEOS-Chem profile is mapped to the TES pressure grid and processed with
!!  the TES averaging kernel before being saved in TES%CH4
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!!     
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!!     
!!  NOTES:
!!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod. 
!!
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE ERROR_MOD,     ONLY : ERROR_STOP
!      USE TIME_MOD,      ONLY : GET_NYMD, GET_TIME_BEHIND_ADJ
!      USE PRESSURE_MOD,  ONLY : GET_BP
!
!      ! Arguments
!      INTEGER            :: NTES
! 
!      ! Local variables 
!      INTEGER            :: YYYYMMDD
!      CHARACTER(LEN=255) :: ROOT_FILENAME
!      CHARACTER(LEN=255) :: READ_FILENAME
!
!      !=================================================================
!      ! GET_GC_PSEUDO_OBS begins here!
!      !=================================================================
!
!
!      ! Filename root
!      res_str = 
!      ROOT_FILENAME = TRIM( '/home/kjw/GEOS-Chem/runs' //
!     &                      '/ch4/TES/ND49_' // GET_RES_EXT() )
!      READ_FILENAME = TRIM( 'tsYYYYMMDD.bpch' )
!
!
!      ! Initialize tau_round_old
!      tau_round_old = -1
!
!      ! Loop over all observations
!      DO NT = NTES, 1, -1
!
!         ! Copy LTES to cleanup code
!         LTES = TES(NT)%LTES(1)
!
!         ! Get Tau value for this observation
!         tau_this  = GET_TAU() - 23 + 24d0 * TES(NT)%TIME(1)
!
!         ! Round Tau value to nearest 3 hours to access ND49 files
!         tau_round = 3*NINT( tau_this/3d0 )
!
!         ! If the rounded tau value is different than previous rounded tau,
!         !    we need to read new datablock from ND49 file
!         IF tau_round /= tau_round_old THEN
!            
!            ! If observation occurs after 01:30:00 AM (UTC)
!            IF ( tau_this >= GET_TAU()-23+1.5 ) THEN
!               YYYYMMDD = GET_NYMD()
!            ! If observation occurs in the early morning (UTC)
!            ENDIF ELSE
!               DATE     = GET_TIME_BEHIND_ADJ( 60*24 )
!               YYYYMMDD = DATE(1)
!            ENDIF
!
!            ! Expand date tokens in filename
!            CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 )
!
!            ! Get Filename of GEOS-Chem output to read
!            READ_FILENAME = TRIM( ROOT_FILENAME ) // TRIM( READ_FILENAME )
!
!            WRITE(6,*) '    - READ_GEOS-Chem_CH4_OBS: reading file: ', 
!     &                           READ_FILENAME
!
!            ! Read data from BPCH file
!            CALL READ_BPCH2( READ_FILENAME, 'IJ-AVG-$',   1,      
!     &                        tau_this,      IGLOB,   JGLOB,      
!     &                           LLPAR,      ARRAY,   QUIET=.TRUE.)
!            CALL TRANSFER_3D( ARRAY(:,:,:), GC_CH4_NATIVE_3D(:,:,:) )
!
!         ENDIF
!
!         ! Get GC column [ppb] --> [v/v]
!         GC_CH4_NATIVE(:) = GC_CH4_NATIVE_3D(II,JJ,:) / 1d9
!
!
!         ! Get I,J indices of grid box corresponding to current TES scan
!         IIJJ = GET_IJ( TES(NT)%LON(1), TES(NT)%LAT(1) )
!         II   = IIJJ(1)
!         JJ   = IIJJ(2)
!
!
!         ! Map GEOS-Chem column to TES pressure grid
!
!            ! Reset variables to be safe
!            MAP(:,:)   = 0d0
!            GC_PRES(:) = 0d0
!
!
!            ! Get GC pressure levels (mbar) 
!            DO L = 1, LLPAR
!               GC_PRES(L) = GET_PCENTER(I,J,L)
!            ENDDO
!
!            ! Get GC surface pressure (mbar) 
!            GC_PSURF = GET_PEDGE(I,J,1) 
!
!            ! Calculate the interpolation weight matrix 
!            MAP(1:LLPAR,1:LTES) 
!     &           = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
!     &                         LTES,  TES(NT)%PRES(1:LTES), GC_PSURF  )
!
!            ! Interpolate GC O3 column to TES grid 
!            DO LL = 1, LTES
!               GC_CH4_onTES(LL) = 0d0 
!               DO L = 1, LLPAR 
!                  GC_CH4_onTES(LL) = GC_CH4_onTES(LL) 
!     &                        + MAP(L,LL) * GC_CH4_NATIVE(L) 
!               ENDDO
!            ENDDO
!
!         !--------------------------------------------------------------
!         ! Apply TES observation operator
!         !
!         !   x_hat = x_a + A_k ( x_m - x_a ) 
!         !  
!         !  where  
!         !    x_hat = GC modeled column as seen by TES [lnvmr]
!         !    x_a   = TES apriori column               [lnvmr]
!         !    x_m   = GC modeled column                [lnvmr]
!         !    A_k   = TES averaging kernel 
!         !--------------------------------------------------------------
!
!         ! x_m - x_a
!         DO L = 1, LTES 
!           GC_CH4_onTES(L)   = MAX(GC_CH4_onTES(L), 1d-10)
!           CH4_PERT_onTES(L) = LOG(GC_CH4_onTES(L)) - 
!     &                               LOG(TES(NT)%PRIOR(L))
!         ENDDO
!     
!         ! x_a + A_k * ( x_m - x_a )  
!         DO L = 1, LTES
!            CH4_HAT_onTES(L)    = 0d0 
!            DO LL = 1, LTES
!               CH4_HAT_onTES(L) = CH4_HAT_onTES(L) 
!     &                + TES(NT)%AVG_KERNEL(L,LL) * CH4_PERT_onTES(LL) 
!            ENDDO
!            CH4_HAT_onTES(L)    = CH4_HAT_onTES(L) 
!     &                + LOG(TES(NT)%PRIOR(L))
!         ENDDO
!
!
!         ! Replace stratospheric values with real TES observations
!         !    to prevent adjoint forcing of stratosphere
!         DO L = 1, LTES
!            IF TES(NT)%PRES(L) < TROPP(II,JJ) THEN
!               CH4_HAT_onTES(L) = LOG( TES(NT)%CH4(L) )
!            ENDIF
!         ENDDO
!
!
!         ! Place GEOS-Chem column into the TES_CH4 structure
!         TES(NT)%GC_CH4(:) = CH4_HAT_onTES(:)
!
!
!      ENDDO    ! End looping over each observation
!
!
!
!      END SUBROUTINE GET_GC_PSEUDO_OBS
!
!!--------------------------------------------------------------------------


!      SUBROUTINE SVD(A,N,U,S,VT)
!!
!!******************************************************************************
!!  Subroutine SVD is a driver for the LAPACK SVD routine DGESVD. (dkh, 05/04/10) 
!! 
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) A   (TYPE (XPLEX)) :  N x N matrix to decompose
!!  (2 ) N  (INTEGER) :  N is dimension of A
!!     
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) U   (TYPE (XPLEX)) :  Array of left singular vectors
!!  (2 ) S   (TYPE (XPLEX)) :  Vector of singular values
!!  (3 ) VT  (TYPE (XPLEX)) :  Array of right singular vectors, TRANSPOSED 
!!      
!!     
!!  NOTES:
!!
!!  Copyright (C) 2009-2010 Intel Corporation. All Rights Reserved.
!!  The information and material ("Material") provided below is owned by Intel
!!  Corporation or its suppliers or licensors, and title to such Material remains
!!  with Intel Corporation or its suppliers or licensors. The Material contains
!!  proprietary information of Intel or its suppliers and licensors. The Material
!!  is protected by worldwide copyright laws and treaty provisions. No part of
!!  the Material may be copied, reproduced, published, uploaded, posted,
!!  transmitted, or distributed in any way without Intel's prior express written
!!  permission. No license under any patent, copyright or other intellectual
!!  property rights in the Material is granted to or conferred upon you, either
!!  expressly, by implication, inducement, estoppel or otherwise. Any license
!!  under such intellectual property rights must be express and approved by Intel
!!  in writing.
!!  =============================================================================
!!
!!  DGESVD Example.
!!  ==============
!!
!!  Program computes the singular value decomposition of a general
!!  rectangular matrix A:
!!
!!    8.79   9.93   9.83   5.45   3.16
!!    6.11   6.91   5.04  -0.27   7.98
!!   -9.15  -7.93   4.86   4.85   3.01
!!    9.57   1.64   8.83   0.74   5.80
!!   -3.49   4.02   9.80  10.00   4.27
!!    9.84   0.15  -8.99  -6.02  -5.31
!!
!!  Description.
!!  ============
!!
!!  The routine computes the singular value decomposition (SVD) of a real
!!  m-by-n matrix A, optionally computing the left and/or right singular
!!  vectors. The SVD is written as
!!
!!  A = U*SIGMA*VT
!!
!!  where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
!!  diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
!!  is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
!!  are the singular values of A; they are real and non-negative, and are
!!  returned in descending order. The first min(m, n) columns of U and V are
!!  the left and right singular vectors of A.
!!
!!  Note that the routine returns VT, not V.
!!
!!  Example Program Results.
!!  ========================
!!
!! DGESVD Example Program Results
!!
!! Singular values
!!  27.47  22.64   8.56   5.99   2.01
!!
!! Left singular vectors (stored columnwise)
!!  -0.59   0.26   0.36   0.31   0.23
!!  -0.40   0.24  -0.22  -0.75  -0.36
!!  -0.03  -0.60  -0.45   0.23  -0.31
!!  -0.43   0.24  -0.69   0.33   0.16
!!  -0.47  -0.35   0.39   0.16  -0.52
!!   0.29   0.58  -0.02   0.38  -0.65
!!
!! Right singular vectors (stored rowwise)
!!  -0.25  -0.40  -0.69  -0.37  -0.41
!!   0.81   0.36  -0.25  -0.37  -0.10
!!  -0.26   0.70  -0.22   0.39  -0.49
!!   0.40  -0.45   0.25   0.43  -0.62
!!  -0.22   0.14   0.59  -0.63  -0.44
!!  =============================================================================
!!******************************************************************************
!!
!      ! Arguements 
!      INTEGER,INTENT(IN)     :: N
!      TYPE (XPLEX), INTENT(IN)     :: A(N,N)
!      TYPE (XPLEX), INTENT(OUT)    :: U(N,N)
!      TYPE (XPLEX), INTENT(OUT)    :: S(N)
!      TYPE (XPLEX), INTENT(OUT)    :: VT(N,N)
!
!      ! Local variables 
!      INTEGER, PARAMETER     :: LWMAX = MAXLEV * 35 
!      INTEGER                :: INFO, LWORK
!      TYPE (XPLEX)       :: WORK( LWMAX )
!
!!     .. External Subroutines ..
!      EXTERNAL               :: DGESVD
!
!!     .. Intrinsic Functions ..
!      INTRINSIC              :: INT, MIN
!
!      !=================================================================
!      ! SVD begins here!
!      !=================================================================
!
!!     .. Executable Statements ..
!      !WRITE(*,*)'DGESVD Example Program Results'
!!
!!     Query the optimal workspace.
!!
!      LWORK = -1
!      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
!     $             WORK, LWORK, INFO )
!      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!!
!!     Compute SVD.
!!
!      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
!     $             WORK, LWORK, INFO )
!!
!!     Check for convergence.
!!
!      IF( INFO.GT.0 ) THEN
!         WRITE(*,*)'The algorithm computing SVD failed to converge.'
!         STOP
!      END IF
!
!!  Uncomment the following to print out singlular values, vectors (dkh, 05/04/10) 
!!!
!!!     Print singular values.
!!!
!!      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
!!!
!!!     Print left singular vectors.
!!!
!!      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
!!     $                   N, N, U, N   )
!!!
!!!     Print right singular vectors.
!!!
!!      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
!!     $                   N, N, VT, N    )
!
!      ! Return to calling program
!      END SUBROUTINE SVD
!!------------------------------------------------------------------------------
!      SUBROUTINE DGESVD_EXAMPLE
!
!!     .. Parameters ..
!      INTEGER          M, N
!      PARAMETER        ( M = 6, N = 5 )
!      INTEGER          LDA, LDU, LDVT
!      PARAMETER        ( LDA = M, LDU = M, LDVT = N )
!      INTEGER          LWMAX
!      PARAMETER        ( LWMAX = 1000 )
!!
!!     .. Local Scalars ..
!      INTEGER          INFO, LWORK
!!
!!     .. Local Arrays ..
!      TYPE (XPLEX) A( LDA, N ), U( LDU, M ), VT( LDVT, N ), S( N ),
!     $                 WORK( LWMAX )
!      DATA             A/
!     $  8.79, 6.11,-9.15, 9.57,-3.49, 9.84,
!     $  9.93, 6.91,-7.93, 1.64, 4.02, 0.15,
!     $  9.83, 5.04, 4.86, 8.83, 9.80,-8.99,
!     $  5.45,-0.27, 4.85, 0.74,10.00,-6.02,
!     $  3.16, 7.98, 3.01, 5.80, 4.27,-5.31
!     $                  /
!!
!!     .. External Subroutines ..
!      EXTERNAL         DGESVD
!      !EXTERNAL         PRINT_MATRIX
!!
!!     .. Intrinsic Functions ..
!      INTRINSIC        INT, MIN
!!
!!     .. Executable Statements ..
!      WRITE(*,*)'DGESVD Example Program Results'
!!
!!     Query the optimal workspace.
!!
!      LWORK = -1
!      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
!     $             WORK, LWORK, INFO )
!      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!!
!!     Compute SVD.
!!
!      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
!     $             WORK, LWORK, INFO )
!!
!!     Check for convergence.
!!
!      IF( INFO.GT.0 ) THEN
!         WRITE(*,*)'The algorithm computing SVD failed to converge.'
!         STOP
!      END IF
!!
!!     Print singular values.
!!
!!      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
!!
!!     Print left singular vectors.
!!
!!      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
!!     $                   M, N, U, LDU )
!!
!!     Print right singular vectors.
!!
!!      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
!!     $                   N, N, VT, LDVT )
!!
!!
!!     End of DGESVD Example.
!      END SUBROUTINE DGESVD_EXAMPLE
!------------------------------------------------------------------------------
!
!     Auxiliary routine: printing a matrix.
!
!      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
!      CHARACTER*(*)    DESC
!      INTEGER          M, N, LDA
!      TYPE (XPLEX) A( LDA, * )
!
!      INTEGER          I, J
!
!      WRITE(*,*)
!      WRITE(*,*) DESC
!      DO I = 1, M
!         WRITE(*,9998) ( A( I, J ), J = 1, N )
!      END DO
!
! Change format of output (dkh, 05/04/10) 
! 9998 FORMAT( 11(:,1X,F6.2) )
! 9998 FORMAT( 11(:,1X,E14.8) )
!      RETURN
!
!      END SUBROUTINE PRINT_MATRIX 
!------------------------------------------------------------------------------

      END MODULE TES_CH4_MOD
