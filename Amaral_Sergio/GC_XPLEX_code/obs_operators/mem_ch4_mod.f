!$Id: mem_ch4_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE MEM_CH4_MOD
!  
!******************************************************************************
!  Module MEM_CH4_MOD for CH4 observations. 
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
      INTEGER, PARAMETER         :: LLMEM   = 13
      INTEGER, PARAMETER         :: MAXMEM  = 639059


      ! Record to store information about the new instrument
      TYPE (XPLEX)                     :: AVGKERNEL(    LLMEM, LLMEM )
      TYPE (XPLEX)                     :: OBSERROR(     LLMEM, LLMEM )
      TYPE (XPLEX)                     :: OBSERROR_INV( LLMEM, LLMEM )
      TYPE (XPLEX)                     :: TOTERROR_INV( LLMEM, LLMEM )
      TYPE (XPLEX)                     :: PRESSURE( LLMEM )
      TYPE (XPLEX)                     :: PRESSURE_EDGE( LLMEM )
      TYPE (XPLEX)                     :: RANDNUM( MAXMEM )
      TYPE (XPLEX), ALLOCATABLE        :: CH4_PRIOR(:,:,:)


      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_MEM_INFO
!
!******************************************************************************
!  Subroutine READ_MEM_INFO reads and stores information about the new
!  instrument, specifically AK, pressure levels and error covariance matrices.
!    (kjw, 07/24/11) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME   (CHAR) : MEM filename to read
!
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE FILE_MOD,              ONLY : IOERROR
      USE TIME_MOD,              ONLY : GET_NYMD
      USE BPCH2_MOD,             ONLY : GET_TAU0, READ_BPCH2
      USE BPCH2_MOD,             ONLY : GET_RES_EXT
      USE ERROR_MOD,             ONLY : ALLOC_ERR

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      CHARACTER(LEN=255)              :: FILENAME
    
      ! Local variables
      CHARACTER(LEN=255)              :: READ_FILENAME

      ! netCDF id's 
      INTEGER                         :: NCID, LG, LN
      INTEGER                         :: nobs_id, yyyymmdd_id, hhmmss_id
      INTEGER                         :: qflag_id, xch4_id, ch4ak_id
      INTEGER                         :: ch4pres_id, ch4prior_id
      INTEGER                         :: gcii_id, gcjj_id, gcfrac_id
      LOGICAL, SAVE                   :: LDEBUG = .TRUE.
      TYPE (XPLEX)                          :: XTAU
      TYPE (XPLEX)                     :: DUMMY_PRIOR(IGLOB,JGLOB,LLMEM)

      ! Loop indexes, and error handling.
      INTEGER                         :: IOS, IU_IN, AS



      !=================================================================
      ! READ_MEM_CH4_OBS begins here!
      !=================================================================

      ! Initialize module variabl
      AVGKERNEL(:,:)    = 0d0
      OBSERROR(:,:)     = 0d0
      OBSERROR_INV(:,:) = 0d0
      TOTERROR_INV(:,:) = 0d0
      PRESSURE(:)       = 0d0
      PRESSURE_EDGE(:)  = 0d0
      RANDNUM(:)        = 0d0


      ! Read and store one variable at a time

         ! ------ Averaging Kernel Matrix ------
         ! Filename to read
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/mem/' ) //
     &                   'data/' // TRIM( 'mem_AK.txt' )
         WRITE(6,*) '    - READ_MEM_AK: reading file: ', 
     &                        TRIM(READ_FILENAME)

         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_avg_kernel:1' )

         ! Read File and save info into module variable AVGKERNEL(:,:)
         DO LN=1,LLMEM
            READ( IU_IN, '(13F12.6)', IOSTAT=IOS ) AVGKERNEL(LN,:)
            
            IF ( LDEBUG ) THEN
               WRITE(6,*) 'Avg Kernel, row ',LN
               WRITE(6,'(13F12.6)') AVGKERNEL(LN,:)
            ENDIF

            ! IO status
            IF ( IOS < 0 ) THEN
               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
               WRITE( 6, '(a)' ) 'STOP in READ_MEM_CH4'
            ENDIF
            IF ( IOS > 0 ) THEN
               CALL IOERROR(IOS, IU_IN, 'read_avg_kernel:2')
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IU_IN )
 

         ! ------ Observation Error Covariance Matrix ------
         ! Filename to read
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/mem/' ) //
     &                   'data/' // TRIM( 'mem_obs_error.txt' )
         WRITE(6,*) '    - READ_MEM_OBSERROR: reading file: ', 
     &                        TRIM(READ_FILENAME)


         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_obs_error:1' )

         ! Read File and save info into module variable OBSERROR(:,:)
         DO LN=1,LLMEM
            READ( IU_IN, '(13F18.12)', IOSTAT=IOS ) OBSERROR(LN,:)

            IF ( LDEBUG ) THEN
               WRITE(6,*) 'Obs Error covar, row ',LN
               WRITE(6,'(13F18.12)') OBSERROR(LN,:)
            ENDIF

            ! IO status
            IF ( IOS < 0 ) THEN
               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
               WRITE( 6, '(a)' ) 'STOP in READ_MEM_CH4'
            ENDIF
            IF ( IOS > 0 ) THEN
               CALL IOERROR(IOS, IU_IN, 'read_obs_error:2')
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IU_IN )


         ! ------ Inverse of Observation Error Covariance Matrix ------
         ! Filename to read
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/mem/' ) //
     &                   'data/' // TRIM( 'mem_obs_error_inv.txt' )
         WRITE(6,*) '    - READ_MEM_OBSERROR_INV: reading file: ',
     &                        TRIM(READ_FILENAME)


         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_obs_error:1' )

         ! Read File and save info into module variable OBSERROR_INV(:,:)
         DO LN=1,LLMEM
            READ( IU_IN, '(13F18.6)', IOSTAT=IOS ) OBSERROR_INV(LN,:)

            IF ( LDEBUG ) THEN
               WRITE(6,*) 'Inv Obs Error covar, row ',LN
               WRITE(6,'(13F18.6)') OBSERROR_INV(LN,:)
            ENDIF

            ! IO status
            IF ( IOS < 0 ) THEN
               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
               WRITE( 6, '(a)' ) 'STOP in READ_MEM_CH4'
            ENDIF
            IF ( IOS > 0 ) THEN
               CALL IOERROR(IOS, IU_IN, 'read_obs_error:2')
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IU_IN )


!         ! ------ Total Error Covariance Matrix ------
!         ! Filename to read
!         READ_FILENAME = TRIM( '/home/kjw/new_satellites/mem/' ) //
!     &                 'data/' // TRIM( 'mem_total_error_inv.txt' )
!         WRITE(6,*) '    - READ_MEM_TOTERROR: reading file: ', 
!     &                        TRIM(READ_FILENAME)
!
!
!         ! Open file
!         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
!     &         STATUS='OLD', IOSTAT=IOS )
!         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_tot_error:1' )
!
!         ! Read File and save info into module variable OBSERROR(:,:)
!         DO LN=1,LLMEM
!            READ( IU_IN, '(13F18.12)', IOSTAT=IOS ) TOTERROR_INV(LN,:)
!
!            ! IO status
!            IF ( IOS < 0 ) THEN
!               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
!               WRITE( 6, '(a)' ) 'STOP in READ_MEM_CH4'
!            ENDIF
!            IF ( IOS > 0 ) THEN
!               CALL IOERROR(IOS, IU_IN, 'read_tot_error:2')
!            ENDIF
!         ENDDO
!
!         ! Close file
!         CLOSE( IU_IN )


         ! ------ Pressure Levels ------
         ! Filename to read
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/mem/' ) //
     &                   'data/' // TRIM( 'mem_pressure.txt' )
         WRITE(6,*) '    - READ_MEM_PRESSURE: reading file: ', 
     &                        TRIM(READ_FILENAME)

         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_pressure:1' )

         ! Read File and save info into module variable PRESSURE(:)
         READ( IU_IN, '(13F12.6)', IOSTAT=IOS ) PRESSURE(:)

         ! IO status
         IF ( IOS < 0 ) THEN
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_MEM_CH4'
         ENDIF
         IF ( IOS > 0 ) THEN
            CALL IOERROR(IOS, IU_IN, 'read_pressure:2')
         ENDIF

         ! Close file
         CLOSE( IU_IN )


         ! ------ Pressure Edges ------
         ! By finite difference on log(pressure) grid
         PRESSURE_EDGE(1) = PRESSURE(1)
         PRESSURE_EDGE(LLMEM) = 0.
         DO LN=2,LLMEM-1
            PRESSURE_EDGE(LN) = exp( log(pressure(LN+1)) + 
     &          ( log(PRESSURE(LN)) - log(PRESSURE(LN+1)) ) / 2.  )
         ENDDO


         ! ------ A priori vertical profiles ------
         ALLOCATE( CH4_PRIOR(IGLOB,JGLOB,LLMEM), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_PRIOR' )
         CH4_PRIOR(:,:,:) = 0d0

         FILENAME = '/home/kjw/new_satellites/mem/data/' //  
     &              'mem_prior.' //  GET_RES_EXT() // '.bpch'
         XTAU = GET_TAU0( 1, 1, 1985 )

         WRITE(6,*) '    - READ_CH4_PRIOR: reading file: ', 
     &                        TRIM(FILENAME)

         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,    
     &                    XTAU,      IGLOB,     JGLOB, 
     &                    LLMEM, DUMMY_PRIOR, QUIET=.TRUE.  )
         CH4_PRIOR(:,:,:) = DUMMY_PRIOR(:,:,:)


         ! LDEBUG = FALSE. Only print values first time reading
         LDEBUG = .FALSE.

    
      ! Return to calling program
      END SUBROUTINE READ_MEM_INFO
!------------------------------------------------------------------------------


      SUBROUTINE CALC_MEM_CH4_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_MEM_CH4_FORCE calculates the adjoint forcing from the MEM
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
      USE BPCH2_MOD,          ONLY : GET_RES_EXT, GET_TAU0
      USE BPCH2_MOD,          ONLY : READ_BPCH2, GET_MODELNAME
      USE BPCH2_MOD,          ONLY : BPCH2,     OPEN_BPCH2_FOR_WRITE
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ, CHECK_STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AD, CLDFRC
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR, ADJTMP_DIR
      USE GRID_MOD,           ONLY : GET_YEDGE, GET_AREA_M2
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,           ONLY : GET_DAY, GET_MONTH, GET_YEAR
      USE TIME_MOD,           ONLY : GET_TAU
      USE TIME_MOD,           ONLY : GET_LOCALTIME, EXPAND_DATE
      USE TRACER_MOD,         ONLY : STT
      USE TRACER_MOD,         ONLY : XNUMOL, XNUMOLAIR
      USE TRACER_MOD,         ONLY : TCVV
      USE ERROR_MOD,          ONLY : ERROR_STOP, IT_IS_NAN
      USE LOGICAL_ADJ_MOD,    ONLY : LDCOSAT
      USE FILE_MOD,           ONLY : IU_RST
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET


#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER, SAVE               :: NT   ! # observations processed this day
      INTEGER                     :: LG, LN, LLN, II, JJ, JMIN, OB
      INTEGER                     :: nlev, lind, IU_IN
      INTEGER                     :: nboxes, nobs
      INTEGER                     :: NTSTART, NTSTOP, NTh, NB
      INTEGER, SAVE               :: NTT
      TYPE (XPLEX)              :: GC_CH4_TRUE_ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: DUMMY_TRUE(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: DUMMY_RAND(IGLOB,JGLOB,1)
      TYPE (XPLEX), SAVE                :: RANDOM_GRID(IGLOB,JGLOB)
      TYPE (XPLEX)                      :: GC_PCENTER(LLPAR)
      TYPE (XPLEX)                      :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)                      :: GC_AD(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_OB(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_onMEM(LLMEM)
      TYPE (XPLEX)                      :: GC_CH4_onMEM_OB(LLMEM)
      TYPE (XPLEX)                      :: GRIDMAP(LLPAR,LLMEM)
      TYPE (XPLEX)                  :: OBSERROR_INV_SUPER(LLMEM,LLMEM)
      TYPE (XPLEX)                      :: SIGN(LLMEM,LLMEM)
      TYPE (XPLEX)                      :: OBSERROR_SQRT(LLMEM,LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_EXP(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_OB(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_OB_EXP(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_ADJ(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_EXP_ADJ(LLMEM)
      TYPE (XPLEX)                      :: CH4_PERT(LLMEM)
      TYPE (XPLEX)                      :: CH4_PERT_OB(LLMEM)
      TYPE (XPLEX)                      :: CH4_PERT_ADJ(LLMEM)
      TYPE (XPLEX)                      :: frac, frac_total
      TYPE (XPLEX)                      :: latmin, Jfrac_min, Jfrac
      TYPE (XPLEX)                      :: box_area, cloud_frac
      TYPE (XPLEX)                      :: mass_air, mole_air, mole_ch4
      TYPE (XPLEX)                      :: LHS, RHS, GC_XCH4, XTAU
      TYPE (XPLEX)                      :: PUP, PLO
      TYPE (XPLEX)                      :: XCH4_HAT, XCH4_HAT_OB
      TYPE (XPLEX)                      :: XCH4_HAT_ADJ, XCH4_HAT_OB_ADJ
      TYPE (XPLEX)                      :: SUPER_ERR, S_obs_inv
      TYPE (XPLEX)                      :: SUPER_ERR_EXPECTED
      TYPE (XPLEX)                      :: XWEIGHT(LLMEM)
      TYPE (XPLEX)                      :: DIFF, FORCE
      TYPE (XPLEX)                      :: sumxweight
      TYPE (XPLEX)                      :: DIFF_ADJ
      TYPE (XPLEX)                      :: GC_CH4_onMEM_ADJ(LLMEM)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: NEW_COST(IIPAR*JJPAR*LLPAR)
      TYPE (XPLEX)                      :: OLD_COST
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      LOGICAL, SAVE               :: DO_FDTEST = .TRUE.
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME
      CHARACTER(LEN=255)          :: FILENAME_OBS

      ! Arrays for saving with satellite diagnostic turned on
      TYPE (XPLEX)                      :: hourly_nobs(IIPAR,JJPAR)
      TYPE (XPLEX)                      :: hourly_xch4_sat(IIPAR,JJPAR)
      TYPE (XPLEX)                   :: hourly_xch4_model(IIPAR,JJPAR)
      TYPE (XPLEX)                      :: DATA_FIELD(IIPAR,JJPAR)
      TYPE (XPLEX)                      :: LONRES, LATRES
      INTEGER                     :: TRACER, I0, J0
      INTEGER, PARAMETER          :: HALFPOLAR = 1
      INTEGER, PARAMETER          :: CENTER180 = 1
      CHARACTER(LEN=20)           :: MODELNAME
      CHARACTER(LEN=40)           :: CATEGORY
      CHARACTER(LEN=40)           :: UNIT
      CHARACTER(LEN=40)           :: RESERVED = ''
      CHARACTER(LEN=80)           :: TITLE

      ! Parameters
      TYPE (XPLEX), PARAMETER           :: XCH4_ERR = xplex(8d0,0d0)

      ! Variables for FD testing
      TYPE (XPLEX)                      :: cost_func_pos, cost_func_neg
      TYPE (XPLEX)                      :: cost_func_0
      TYPE (XPLEX)                      :: PERT(LLPAR)
      TYPE (XPLEX)                      :: ADJ_SAVE(LLPAR)
      TYPE (XPLEX)                      :: ADJ(LLPAR)
      TYPE (XPLEX)                      :: FD_CEN(LLPAR)
      TYPE (XPLEX)                      :: FD_POS(LLPAR)
      TYPE (XPLEX)                      :: FD_NEG(LLPAR)
      TYPE (XPLEX)                      :: DOFS


      !=================================================================
      ! CALC_MEM_CH4_FORCE begins here!
      !=================================================================

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


         ! Read CH4 data
         CALL READ_MEM_INFO

         ! Initialize counter for total number of observations processed
         NTT = 0


         FIRST = .FALSE.   ! only open files on first call to 
      ENDIF


!      ! Open file for this hour's satellite diagnostics
!      FILENAME = 'diag_sat.YYYYMMDD.hhmm.NN'
!      CALL EXPAND_NAME( FILENAME, N_CALC )
!      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
!      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!      OPEN( IU_FILE,      FILE=TRIM( FILENAME ),  STATUS='UNKNOWN',
!     &       IOSTAT=IOS,  FORM='FORMATTED',       ACCESS='SEQUENTIAL' )


      ! Check that we haven't added any NaN to the STT_ADJ array
      CALL CHECK_STT_ADJ( 'Start of CALC_MEM_CH4_FORCE' )

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Read "TRUE" state for this time step [kg/box]
      GC_CH4_TRUE_ARRAY(:,:,:) = 0d0
!      FILENAME_OBS = '/home/kjw/GEOS-Chem/gcadj_std/runs/v8-02-01/' //
!     &               'ch4/mem/' // GET_RES_EXT() // '/adjtmp/'  //
!     &               'gctm.obs.YYYYMMDD.hhmm'
      FILENAME_OBS = 'gctm.obs.YYYYMMDD.hhmm'      
      CALL EXPAND_DATE( FILENAME_OBS, GET_NYMD(), GET_NHMS() )
      FILENAME_OBS = TRIM( ADJTMP_DIR   )  //  TRIM( FILENAME_OBS )
      XTAU = GET_TAU()
      CALL READ_BPCH2( TRIM(FILENAME_OBS), 'IJ-OBS-$', 1,
     &                 XTAU,          IIPAR,     JJPAR,
     &                 LLPAR,         DUMMY_TRUE ,  QUIET=.TRUE.)
      GC_CH4_TRUE_ARRAY(:,:,:) = DUMMY_TRUE(:,:,:)

      ! Convert from [kg] --> [v/v]
      DO II=1,IIPAR
      DO JJ=1,JJPAR
      DO LG=1,LLPAR
         GC_CH4_TRUE_ARRAY(II,JJ,LG) = GC_CH4_TRUE_ARRAY(II,JJ,LG)
     &        * ( 1d3 / 16d0 ) / ( AD(II,JJ,LG) * 1d3 / 28.96 )
      ENDDO
      ENDDO
      ENDDO



! need to update this in order to do i/o with this loop parallel 
!!      ! Now do a parallel loop for analyzing data 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, MAP, LLNT, IIJJ,  I, J,  L,   LL    )
!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_CH4,  DIFF  )
!!$OMP+PRIVATE( GC_CH4_NATIVE, CH4_PERT, CH4_HAT, FORCE )
!!$OMP+PRIVATE( ADJ_GC_CH4_NATIVE,       ADJ_GC_CH4     )
!!$OMP+PRIVATE( ADJ_CH4_PERT,            ADJ_CH4_HAT    )
!!$OMP+PRIVATE( ADJ_DIFF                                )

      ! If new day of observations initialize count
      IF ( GET_NHMS() .EQ. 230000 ) THEN
         NT = 0   ! initialize counter of total observations processed today
         NB = 0   ! initialize counter of total boxes processed today

         ! ------ Random Numbers ------
         ! Read error values for this day
         XTAU     = GET_TAU0( GET_MONTH(), GET_DAY(), GET_YEAR() )
         FILENAME = '/home/kjw/new_satellites/mem/data/randnums/'
     &        // 'random.YYYYMMDD.' // GET_RES_EXT() // '.bpch'
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), 0 )
         CALL READ_BPCH2( TRIM(FILENAME), 'IJ-AVG-$', 1,
     &        XTAU,          IGLOB,     JGLOB,
     &        1,             DUMMY_RAND ,  QUIET=.TRUE.)
         RANDOM_GRID(:,:) = DUMMY_RAND(:,:,1)

      ENDIF

      ! Get grid offsets for use with nested grid
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
      

      ! Begin counter for
      NTh = 0    ! number of observations processed this hour
      NB  = 0    ! number of grid boxes processed this hour

      ! Clear satellite diagnostic information to be safe
      IF ( LDCOSAT .EQ. .TRUE. ) THEN 
         hourly_nobs(:,:)       = 0d0
         hourly_xch4_sat(:,:)   = 0d0
         hourly_xch4_model(:,:) = 0d0
      ENDIF

      ! Information for spatial criteria for observations
      latmin  = 40.0

      ! Determine minimum JJ index over which to look for observations
      DO JJ=1, JJPAR-1
         IF ( ( GET_YEDGE(JJ)   .LE. latmin ) .AND. 
     &        ( GET_YEDGE(JJ+1) .GT. latmin ) ) THEN
            JMIN      = JJ
            Jfrac_min = ( GET_YEDGE(JJ+1) - latmin        )  / 
     &                  ( GET_YEDGE(JJ+1) - GET_YEDGE(JJ) )
         ENDIF
      ENDDO

      print*, '    - CALC_MEM_CH4_FORCE  ', GET_NYMD(), GET_NHMS()


      ! Loop over each grid box north of the minimum latitude
      !   1. Determine number of observations in the current grid box
      !   2. Make "super-observation" in current grid box
      !        "super-observation" is one observation with error and
      !        associated error covariance matrix scaled to sqrt(N)
      !        where N is the number of regular observations in box
      DO II = 1,    IIPAR

         ! If not 1400 <= local time < 1500, cycle to next II value
         IF ( ( GET_LOCALTIME( II ) .LT. 14.00 )  .OR.
     &        ( GET_LOCALTIME( II ) .GE. 15.00 ) )  CYCLE

         ! It is 1400-1500 local time, so let's make observations!
         DO JJ = JMIN, JJPAR

            ! For safety, initilize these variables
            nobs       = 0
            cloud_frac = 0.
            box_area   = 0.
            GC_PCENTER(:)    = 0d0
            GC_PEDGE(:)      = 0d0
            GC_AD(:)         = 0d0
            GC_CH4_NATIVE(:) = 0d0
            GC_CH4_onMEM(:)    = 0d0
            GC_CH4_onMEM_OB(:) = 0d0


            ! Fraction of grid box above minimum latitude
            Jfrac = 1.
            IF ( JJ .EQ. JMIN ) Jfrac = Jfrac_min

            ! Determine number of observations in this grid box
            !  # obs = box_area * (1-cloud_fraction) * Jfrac / 100
            !      divide by 100 because each observation takes up 100 km2
            box_area   = GET_AREA_M2( JJ ) * 1d-6 ! [m2] --> [km2]
            cloud_frac = CLDFRC( II, JJ )
            nobs = NINT( ( (1-cloud_frac) * box_area * Jfrac ) / 100. )
            nobs = 10
            IF ( nobs .LT. 1 ) CYCLE 


            ! Get GEOS-Chem pressure and CH4 column corresponding to this grid box.
            ! CH4 in [kg/box] and pressure in [hPa]
            ! Get column of pressure centers and CH4 values
            DO LG=1,LLPAR

               ! Pressure centers [hPa]
               GC_PCENTER(LG) = GET_PCENTER(II,JJ,LG)

               ! Pressure edges [hPa]
               GC_PEDGE(LG)   = GET_PEDGE(II,JJ,LG)

               ! mass per box [kg]
               GC_AD(LG)      = AD(II,JJ,LG)

               ! CH4 values [kg/box] --> [v/v]
               GC_CH4_NATIVE(LG) = ( CHK_STT(II,JJ,LG,1 )
     &                    * XNUMOL(1) ) / ( GC_AD(LG) * XNUMOLAIR )

            ENDDO


            ! Number of vertical levels to use in these observations
            !   Chop off lowermost levels if 
            !     GEOS-Chem surface pressure < MEM pressure levels
            nlev = count( (PRESSURE_EDGE%r) .LT. (GC_PEDGE(1)%r) )
            !IF ( nlev .LT. 13 )  nlev = nlev + 1
            lind = LLMEM + 1 - nlev ! minimum vertical index on MEM grid

            ! Get interpolation matrix that maps GEOS-Chem to MEM grid
            GRIDMAP(1:LLPAR, 1:LLMEM) = 
     &                   GET_INTMAP( GC_PEDGE, PRESSURE_EDGE, nlev )

            ! Get GEOS-Chem column from "truth" run to make pseudo-observations
            GC_CH4_NATIVE_OB(:) = 0d0
            GC_CH4_NATIVE_OB(:) = GC_CH4_TRUE_ARRAY(II,JJ,:)

            ! Interpolate GEOS-Chem CH4 column and observation to MEM grid
            ! Column in [v/v]
            DO LN = lind, LLMEM
               GC_CH4_onMEM(LN) = 0d0 
               GC_CH4_onMEM_OB(LN) = 0d0 
               DO LG = 1, LLPAR
                  GC_CH4_onMEM(LN)    = GC_CH4_onMEM(LN) 
     &                 + GRIDMAP(LG,LN) * GC_CH4_NATIVE(LG) 
                  GC_CH4_onMEM_OB(LN) = GC_CH4_onMEM_OB(LN) 
     &                 + GRIDMAP(LG,LN) * GC_CH4_NATIVE_OB(LG) 
               ENDDO
            ENDDO

            !--------------------------------------------------------------
            ! Apply MEM observation operator
            !
            !   x_hat = x_a + A_k ( x_m - x_a ) 
            !  
            !  where  
            !    x_hat = GC modeled column as seen by MEM [ln(vmr)]
            !    x_a   = MEM apriori column               [ln(vmr)]
            !    x_m   = GC modeled column on MEM grid    [ln(vmr)]
            !    A     = MEM averaging kernel 
            !--------------------------------------------------------------
            
            ! x_m - x_a for model and "observation"
            !    [v/v] --> ln( v/v ) happens here
            DO LN = lind, LLMEM
              GC_CH4_onMEM(LN)   =MAX(GC_CH4_onMEM(LN),   1d-10)
              GC_CH4_onMEM_OB(LN)=MAX(GC_CH4_onMEM_OB(LN),1d-10)
              CH4_PERT(LN)           =LOG( GC_CH4_onMEM(LN) ) - 
     &                                LOG( CH4_PRIOR(II,JJ,LN) )
              CH4_PERT_OB(LN)        =LOG( GC_CH4_onMEM_OB(LN) ) - 
     &                                LOG( CH4_PRIOR(II,JJ,LN) )
            ENDDO

            ! x_a + A_k * ( x_m - x_a ) for model and "observation"
            DO LN = lind, LLMEM
               CH4_HAT(LN)    = 0d0
               CH4_HAT_OB(LN) = 0d0
            
               DO LLN = lind, LLMEM
                  CH4_HAT(LN) = CH4_HAT(LN) 
     &                       + AVGKERNEL(LN,LLN) * CH4_PERT(LLN)
                  CH4_HAT_OB(LN) = CH4_HAT_OB(LN) 
     &                       + AVGKERNEL(LN,LLN) * CH4_PERT_OB(LLN)
               ENDDO
               CH4_HAT(LN)   = CH4_HAT(LN)   +LOG( CH4_PRIOR(II,JJ,LN) )
               CH4_HAT_OB(LN)= CH4_HAT_OB(LN)+LOG( CH4_PRIOR(II,JJ,LN) )
            
            ENDDO


            ! Convert vertical profiles from [ln(vmr)] --> [vmr] before 
            !   calculating XCH4
            CH4_HAT_EXP    = EXP(CH4_HAT)
            CH4_HAT_OB_EXP = EXP(CH4_HAT_OB)


            ! ---- Calculate XCH4 [v/v] from CH4_HAT [v/v] and CH4_HAT_OB [v/v]
            XCH4_HAT    = 0d0
            XCH4_HAT_OB = 0d0
            
            ! Calculate weight of each vertical level on MEM grid for averaging
            !   levels to get XCH4. Weight by # molecules / verical level, which is
            !   proportional to pressure difference between upper and lower bounds
            !   of each box.
            DO LN=lind, LLMEM
            
               ! If ground level, average with same weight as if it were 1st atm level
               IF ( LN .EQ. lind ) THEN
                  PUP = PRESSURE_EDGE(LN+1)
                  PLO = PRESSURE_EDGE(LN  )
               ELSE
                  PUP = PRESSURE_EDGE(LN  )
                  PLO = PRESSURE_EDGE(LN-1)
               ENDIF
            
               Xweight(LN) = PLO - PUP
            ENDDO
            
            !Normalize so that SUM(Xweight) = 1
            sumxweight = SUM( Xweight(:) )
            DO LN=lind,LLMEM
               Xweight(LN) = Xweight(LN) / sumxweight
            ENDDO

            ! Calculate weighted average of CH4_HAT and CH4_HAT_OB
            DO LN=lind, LLMEM
               XCH4_HAT    = XCH4_HAT    + Xweight(LN) * CH4_HAT_EXP(LN)
               XCH4_HAT_OB = XCH4_HAT_OB + 
     &                         Xweight(LN) * CH4_HAT_OB_EXP(LN)
            ENDDO

!            if (( II .eq. 11 ) .AND. (JJ .eq. 39)) then
!               print*,'lind = ',lind
!               DO LN=lind,LLMEM
!                  print*, LN, xweight(LN),
!     &               GC_CH4_onMEM(LN), ch4_hat_exp(LN)
!               ENDDO
!               print*,'---------------------------------------'
!               WRITE(6,'(14F16.8)') 0d0, PRESSURE_EDGE(:)
!               DO LG=1,LLPAR
!                  WRITE(6,'(14F16.8)') GC_PEDGE(LG), GRIDMAP(LG,:)
!               ENDDO
!               print*,'---------------------------------------'
!            endif

            ! Create super observation by adding random error
            ! to XCH4_HAT_OB
            ! SUPER_ERR is 1d-9 * XCH4_ERR[ppb] * N(0,1) / sqrt(nobs)   [v/v]
            !     where 8ppb is expected error on a single XCH4 measurement
            !           N(0,1) is a random number of mean 0, standard deviation 1
            !           nobs is the number of observations merged to form super obs
            ! Expected error of super-observation XCH4
            SUPER_ERR_EXPECTED = 1d-9 * XCH4_ERR / SQRT( XPLX(nobs) )

            ! Multiply expected error of super-observation by 
            !   prescribed random number with mean 0, standard deviation 1
            SUPER_ERR = SUPER_ERR_EXPECTED * RANDOM_GRID( II+I0, JJ+J0 ) 

            ! Add random error to super-observation
            XCH4_HAT_OB = XCH4_HAT_OB + SUPER_ERR    ! add error [v/v]


            !--------------------------------------------------------------
            ! Calculate cost function, given S is observation error 
            ! covariance matrix.
            !     Sobs = 1x1 array [ ln(vmr)^2 ]
            ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
            !--------------------------------------------------------------

            ! Initialize values to be safe
            DIFF        = 0d0
            FORCE       = 0d0
            
            ! Calculate difference between modeled and observed profile
            DIFF        = XCH4_HAT - XCH4_HAT_OB
            
            ! Calculate adjoint forcing: 2 * DIFF^T * S_{obs}^{-1}
            !         and cost function: DIFF^T * S_{obs}^{-1} * DIFF
            ! Inverse observation error covariance matrix of super-obs
            S_obs_inv = 1d0 / (SUPER_ERR_EXPECTED**2)
            FORCE = 2 * DIFF * S_obs_inv
            NEW_COST(NB) = 0.5d0 * DIFF * FORCE

!            print*,'DIFF, XCH4_HAT, XCH4_HAT_OB',
!     &                 DIFF, XCH4_HAT, XCH4_HAT_OB
!            print*,'DIFF, FORCE, S_obs_inv',
!     &                 DIFF, FORCE, S_obs_inv
!            print*,'NB, NEW_COST(NB) = ',NB, NEW_COST(NB) 

         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

         ! Initialize to be safe
         DIFF_ADJ = 0d0

         ! The adjoint forcing is 2 * S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE

         ! Adjoint of GEOS-Chem - Observation difference
         XCH4_HAT_ADJ = DIFF_ADJ


         ! Adjoint of CH4_HAT_EXP --> XCH4_HAT
         DO LN=lind, LLMEM
            CH4_HAT_EXP_ADJ(LN) = XCH4_HAT_ADJ * Xweight(LN)
         ENDDO

         ! Adjoint of CH4_HAT --> CH4_HAT_EXP
         DO LN=lind, LLMEM
            CH4_HAT_ADJ(LN) = CH4_HAT_EXP_ADJ(LN) * CH4_HAT_EXP(LN)
         ENDDO

         ! Adjoint of MEM observation operator
         CH4_PERT_ADJ(:) = 0D0
         DO LN=lind,LLMEM
         DO LLN=lind,LLMEM
            CH4_PERT_ADJ(LN) = CH4_PERT_ADJ(LN) + 
     &              AVGKERNEL(LLN,LN) * CH4_HAT_ADJ(LLN)
         ENDDO
         ENDDO

         ! Adjoint of x_m - x_a
         DO LN = lind, LLMEM 
           ! fwd code:
           !GC_CH4(LN)   = MAX(GC_CH4(LN), 1d-10)
           !CH4_PERT(LN) = LOG(GC_CH4(LN)) - LOG(PRIOR(LN))
           ! adj code:
           IF ( GC_CH4_onMEM(LN) > 1d-10 ) THEN 
              GC_CH4_onMEM_ADJ(LN) = 1d0 / GC_CH4_onMEM(LN) * 
     &                                         CH4_PERT_ADJ(LN)
           ELSE 
              GC_CH4_onMEM_ADJ(LN) = 1d0 / 1d-10 * CH4_PERT_ADJ(LN)
           ENDIF
         ENDDO


         ! Adjoint of interpolation
         DO LN=lind,LLMEM
         DO LG=1,LLPAR
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) + 
     &           GRIDMAP(LG,LN) * GC_CH4_onMEM_ADJ(LN)
         ENDDO
         ENDDO


         ! Adjoint of unit conversion [kg/box] --> [v/v]
         DO LG=1,LLPAR
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG)
     &                 * XNUMOL(1) / ( XNUMOLAIR * GC_AD(LG) )
         ENDDO


         ! Pass adjoing forcing back to adjoint tracer array
         DO LG=1,LLPAR
            STT_ADJ(II,JJ,LG,1) = STT_ADJ(II,JJ,LG,1) + 
     &             GC_CH4_NATIVE_ADJ(LG)
         ENDDO


         ! Update cost function 
         COST_FUNC = COST_FUNC + NEW_COST(NB)
         print*,'--------------------------------'
         print*,'I,J = ',II,JJ
         CALL CHECK_STT_ADJ( 'Inside CALC_MEM_CH4_FORCE' )
         print*,' COST_FUNC, NEW_COST(NB) = ',COST_FUNC, NEW_COST(NB)


         ! Record information for satellite diagnostics
         IF ( LDCOSAT .EQ. .TRUE. ) THEN 
            hourly_nobs(II,JJ)     = hourly_nobs(II,JJ) + nobs
            hourly_xch4_sat(II,JJ) =hourly_xch4_sat(II,JJ) + XCH4_HAT_OB
            hourly_xch4_model(II,JJ)=hourly_xch4_model(II,JJ) + XCH4_HAT
         ENDIF


         ! Increment counters
         NTh = NTh + nobs       ! # obs processed this hour4
         NT  = NT  + nobs       ! # obs processed today
         NTT = NTT + nobs       ! # obs processed total
         NB  = NB  + 1          ! # boxes processed this hour

         ENDDO   ! End looping over each grid box JJ
      ENDDO   ! End looping over each grid box II

      ! Save satellite diagnostic information to file
      IF ( LDCOSAT .EQ. .TRUE. ) THEN
         FILENAME = TRIM( DIAGADJ_DIR ) // 'sat.diagnostic.mem.' //
     &        'YYYYMMDD.hhmm.NN'
         TITLE     = 'Satellite Observation Diagnostic File'
         UNIT      = '[v/v]'
         CATEGORY  = 'IJ-AVG-$'
         MODELNAME = GET_MODELNAME()
         LONRES    = DISIZE
         LATRES    = DJSIZE         

         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
         CALL EXPAND_NAME( FILENAME, N_CALC )
      
         ! Open BPCH file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

         ! Write values to bpch
         TRACER          = 1
         DATA_FIELD(:,:) = hourly_nobs
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,   
     &               HALFPOLAR, CENTER180, CATEGORY,  TRACER,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,  
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         DATA_FIELD )
         TRACER          = 2
         DATA_FIELD(:,:) = hourly_xch4_sat
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,   
     &               HALFPOLAR, CENTER180, CATEGORY,  TRACER,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,  
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         DATA_FIELD )
         TRACER          = 3
         DATA_FIELD(:,:) = hourly_xch4_model
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,   
     &               HALFPOLAR, CENTER180, CATEGORY,  TRACER,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,  
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         DATA_FIELD )

         ! Close file
         CLOSE( IU_RST )

      ENDIF


      ! Check that we haven't added any NaN to the STT_ADJ array
      CALL CHECK_STT_ADJ( 'End of CALC_MEM_CH4_FORCE' )


!!$OMP END PARALLEL DO


! -----------------------------------------------------------------------
!    Use this section to test the adjoint of the MEM_CH4 operator by
!          slightly perturbing model [CH4] and recording resultant change
!          in calculated contribution to the cost function.
!
!    This routine will write the following information for each observation
!          to rundir/diagadj/test_adjoint_obs.NN.m
!
!    The adjoint of the observation operator has been tested and validated
!          as of 7/20/10, kjw.
!
      IF ( DO_FDTEST ) THEN
      WRITE(116,210) '   LG'       ,       '  TROP', '     GC_PRES',
     &               '      FD_POS', '      FD_NEG', '      FD_CEN',
     &               '         ADJ', '    COST_POS', '    COST_NEG', 
     &               '  FD_POS/ADJ', '  FD_NEG/ADJ', '  FD_CEN/ADJ'
      PERT(:) = 0D0

      COST_FUNC_0 = 0d0
      CALL CALC_MEM_CH4_FORCE_FD( COST_FUNC_0, PERT, ADJ )
      ADJ_SAVE(:) = ADJ(:)

      !DO LN=lind,LLMEM
      !   DOFS = DOFS + AVGKERNEL(LN,LN)
      !ENDDO

      ! Write identifying information to top of satellite diagnostic file
      WRITE(116,212) 'COST_FUNC_0: ',( COST_FUNC_0 )


      ! Perform finite difference testing at each vertical level
      DO LG = 1, 47

         ! Positive perturbation to GEOS-Chem CH4 columns
         PERT(:) = 0.0
         PERT(LG) = 0.001
         COST_FUNC_pos = 0D0
         CALL CALC_MEM_CH4_FORCE_FD( COST_FUNC_pos, PERT, ADJ )

         ! Negative perturbation to GEOS-Chem CH4 columns
         PERT(:) = 0.0
         PERT(LG) = -0.001
         COST_FUNC_neg = 0D0
         CALL CALC_MEM_CH4_FORCE_FD( COST_FUNC_neg, PERT, ADJ )

         ! Calculate dJ/dCH4 from perturbations
         FD_CEN(LG) =(COST_FUNC_pos - COST_FUNC_neg) / (2*abs(PERT(LG)))
         FD_POS(LG) = ( COST_FUNC_pos - COST_FUNC_0 )   / abs(PERT(LG))
         FD_NEG(LG) = ( COST_FUNC_0 - COST_FUNC_neg )   / abs(PERT(LG))

         ! Write information to satellite diagnostic file
         WRITE(116, 211)  LG,      GC_PCENTER(LG),
     &                    FD_POS(LG),   FD_NEG(LG),
     &                    FD_CEN(LG),   ADJ_SAVE(LG), 
     &                    COST_FUNC_pos, COST_FUNC_neg, 
     &                    FD_POS(LG)/ADJ_SAVE(LG),
     &                    FD_NEG(LG)/ADJ_SAVE(LG),
     &                    FD_CEN(LG)/ADJ_SAVE(LG)
      ENDDO


      WRITE(116,'(a)') '----------------------------------------------'

 210  FORMAT(A4,2x,A6,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,
     &       A12,2x,A12,2x,A12,2x,A12,2x)
 211  FORMAT(I4,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,
     &        2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6)
 212  FORMAT(A12,F22.6)
 213  FORMAT(A12,I4)
 214  FORMAT(I4,2x,F18.6,2x,F18.6)
! -----------------------------------------------------------------------
         DO_FDTEST = .FALSE.
      ENDIF  ! IF ( DO_FDTEST )



      ! Update cost function 
      !COST_FUNC = COST_FUNC + SUM(NEW_COST(:))

      print*, ' Updated value of COST_FUNC     = ', COST_FUNC 
      print*, ' MEM contribution this hour = ', COST_FUNC - OLD_COST 
      print*, ' # Obs analyzed this hour       = ', NTh
      print*, ' # Obs analyzed today           = ', NT
      print*, ' # Obs analyzed total           = ', NTT



      ! Return to calling program
      END SUBROUTINE CALC_MEM_CH4_FORCE

!------------------------------------------------------------------------------



      SUBROUTINE CALC_MEM_CH4_FORCE_FD( COST_FUNC_A, PERT, ADJ )
!
!******************************************************************************
!  Subroutine CALC_MEM_CH4_FORCE calculates the adjoint forcing from the MEM
!  CH4 observations and updates the cost function. (kjw, 07/20/11)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC_A (TYPE (XPLEX)) : Cost funciton (INOUT)                 [unitless]
!  (2 ) PERT        (TYPE (XPLEX)) : Array of perturbations to CH4 column (+/- 0.1, for ex.)
!  (5 ) ADJ         (TYPE (XPLEX)) : Array of adjoint forcings (OUT)
!     
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD,          ONLY : GET_RES_EXT, GET_TAU0
      USE BPCH2_MOD,          ONLY : READ_BPCH2
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AD, CLDFRC
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR, ADJTMP_DIR
      USE GRID_MOD,           ONLY : GET_YEDGE, GET_AREA_M2
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAU
      USE TIME_MOD,           ONLY : GET_LOCALTIME, EXPAND_DATE
      USE TRACER_MOD,         ONLY : STT
      USE TRACER_MOD,         ONLY : XNUMOL, XNUMOLAIR
      USE TRACER_MOD,         ONLY : TCVV
      USE ERROR_MOD,          ONLY : ERROR_STOP

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC_A
      TYPE (XPLEX), INTENT(OUT)         :: ADJ(LLPAR)
      TYPE (XPLEX), INTENT(IN)          :: PERT(LLPAR)


      ! Local variables
      INTEGER, SAVE               :: NT   ! # observations processed this day
      INTEGER                     :: LG, LN, LLN, II, JJ, JMIN, OB
      INTEGER                     :: nlev, lind, IU_IN
      INTEGER                     :: nboxes, nobs
      INTEGER                     :: NTSTART, NTSTOP, NTh, NB
      INTEGER, SAVE               :: NTT
      TYPE (XPLEX)               :: GC_CH4_TRUE_ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                    :: DUMMY_PRIOR(IIPAR,JJPAR,LLMEM)
      TYPE (XPLEX)                      :: DUMMY_TRUE(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: RANDOM_GRID(IGLOB,JGLOB)
      TYPE (XPLEX)                      :: GC_PCENTER(LLPAR)
      TYPE (XPLEX)                      :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)                      :: GC_AD(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_OB(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_onMEM(LLMEM)
      TYPE (XPLEX)                      :: GC_CH4_onMEM_OB(LLMEM)
      TYPE (XPLEX)                      :: GRIDMAP(LLPAR,LLMEM)
      TYPE (XPLEX)                   :: OBSERROR_INV_SUPER(LLMEM,LLMEM)
      TYPE (XPLEX)                      :: SIGN(LLMEM,LLMEM)
      TYPE (XPLEX)                      :: OBSERROR_SQRT(LLMEM,LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_EXP(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_OB(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_OB_EXP(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_ADJ(LLMEM)
      TYPE (XPLEX)                      :: CH4_HAT_EXP_ADJ(LLMEM)
      TYPE (XPLEX)                      :: CH4_PERT(LLMEM)
      TYPE (XPLEX)                      :: CH4_PERT_OB(LLMEM)
      TYPE (XPLEX)                      :: CH4_PERT_ADJ(LLMEM)
      TYPE (XPLEX)                      :: frac, frac_total
      TYPE (XPLEX)                      :: latmin, Jfrac_min, Jfrac
      TYPE (XPLEX)                      :: box_area, cloud_frac
      TYPE (XPLEX)                      :: mass_air, mole_air, mole_ch4
      TYPE (XPLEX)                      :: LHS, RHS, GC_XCH4, XTAU
      TYPE (XPLEX)                      :: PUP, PLO
      TYPE (XPLEX)                      :: XCH4_HAT, XCH4_HAT_OB
      TYPE (XPLEX)                      :: XCH4_HAT_ADJ, XCH4_HAT_OB_ADJ
      TYPE (XPLEX)                      :: SUPER_ERR, S_obs_inv
      TYPE (XPLEX)                      :: SUPER_ERR_EXPECTED
      TYPE (XPLEX)                      :: XWEIGHT(LLMEM)
      TYPE (XPLEX)                      :: DIFF, FORCE
      TYPE (XPLEX)                      :: sumxweight
      TYPE (XPLEX)                      :: DIFF_ADJ
      TYPE (XPLEX)                      :: GC_CH4_onMEM_ADJ(LLMEM)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: NEW_COST(IIPAR*JJPAR*LLPAR)
      TYPE (XPLEX)                      :: OLD_COST
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      LOGICAL, SAVE               :: DO_FDTEST = .TRUE.
      LOGICAL, SAVE               :: LDEBUG = .FALSE.
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME
      CHARACTER(LEN=255)          :: FILENAME_OBS

      ! Parameters
      TYPE (XPLEX), PARAMETER           :: XCH4_ERR = xplex(8d0,0d0)


      !=================================================================
      ! CALC_MEM_CH4_FORCE_FD begins here!
      !=================================================================

      print*, '     - CALC_MEM_CH4_FORCE_FD '

      NEW_COST(:) = 0d0 


      ! ---- Read "TRUE" state for this time step ----
      GC_CH4_TRUE_ARRAY(:,:,:) = 0d0
!      FILENAME_OBS = '/home/kjw/GEOS-Chem/gcadj_std/runs/v8-02-01/' //
!     &               'ch4/mem/' // GET_RES_EXT() // '/adjtmp/'  //
!     &               'gctm.obs.YYYYMMDD.hhmm'
      FILENAME_OBS = 'gctm.obs.YYYYMMDD.hhmm'
      CALL EXPAND_DATE( FILENAME_OBS, GET_NYMD(), GET_NHMS() )
      FILENAME_OBS = TRIM( ADJTMP_DIR   )  //  TRIM( FILENAME_OBS )
      XTAU = GET_TAU()
      CALL READ_BPCH2( TRIM(FILENAME_OBS), 'IJ-OBS-$', 1,
     &                 XTAU,          IIPAR,     JJPAR,
     &                 LLPAR,         DUMMY_TRUE,   QUIET=.TRUE.)
      GC_CH4_TRUE_ARRAY(:,:,:) = DUMMY_TRUE(:,:,:)

      ! Convert from [kg] --> [v/v]
      DO II=1,IIPAR
      DO JJ=1,JJPAR
      DO LG=1,LLPAR
         GC_CH4_TRUE_ARRAY(II,JJ,LG) = GC_CH4_TRUE_ARRAY(II,JJ,LG)
     &        * ( 1d3 / 16d0 ) / ( AD(II,JJ,LG) * 1d3 / 28.96 )
      ENDDO
      ENDDO
      ENDDO


      ! Select arbitrary II, JJ and NT value
      II=40
      JJ=JJPAR-10
      NB=100
      RANDOM_GRID(:,:)   = 0d0
      RANDOM_GRID(II,JJ) = 1.00

      ! Initialize variables
      GC_PCENTER(:)    = 0d0
      GC_PEDGE(:)      = 0d0
      GC_AD(:)         = 0d0
      GC_CH4_NATIVE(:) = 0d0
      GC_CH4_onMEM(:)    = 0d0
      GC_CH4_onMEM_OB(:) = 0d0
      DIFF          = 0d0
      FORCE         = 0d0


      ! Fraction of grid box above minimum latitude
      Jfrac = 1.
      IF ( JJ .EQ. JMIN ) Jfrac = Jfrac_min

      ! Determine number of observations in this grid box
      !  # obs = box_area * (1-cloud_fraction) * Jfrac / 100
      !      divide by 100 because each observation takes up 100 km2
      box_area   = GET_AREA_M2( JJ ) * 1d-6 ! [m2] --> [km2]
      cloud_frac = CLDFRC( II, JJ )
      nobs = NINT( ( (1-cloud_frac) * box_area * Jfrac ) / 100. )

      ! Get GEOS-Chem pressure and CH4 column corresponding to this grid box.
      ! CH4 in [kg/box] and pressure in [hPa]
      ! Get column of pressure centers and CH4 values
      DO LG=1,LLPAR

         ! Pressure centers [hPa]
         GC_PCENTER(LG) = GET_PCENTER(II,JJ,LG)

         ! Pressure edges [hPa]
         GC_PEDGE(LG)   = GET_PEDGE(II,JJ,LG)

         ! mass per box [kg]
         GC_AD(LG)      = AD(II,JJ,LG)

         ! CH4 values [kg/box] --> [v/v]
         GC_CH4_NATIVE(LG) = ( CHK_STT(II,JJ,LG,1)
     &       * (1+PERT(LG)) * XNUMOL(1) ) / ( GC_AD(LG) * XNUMOLAIR )

      ENDDO


      ! Number of vertical levels to use in these observations
      !   Chop off lowermost levels if 
      !     GEOS-Chem surface pressure < MEM pressure levels
      nlev = count( (PRESSURE_EDGE%r) .LT. (GC_PEDGE(1)%r) )
      IF ( nlev .LT. 13 )  nlev = nlev + 1
      lind = LLMEM + 1 - nlev ! minimum vertical index on MEM grid

      ! Get interpolation matrix that maps GEOS-Chem to MEM grid
      GRIDMAP(1:LLPAR, 1:LLMEM) = 
     &     GET_INTMAP( GC_PEDGE, PRESSURE_EDGE, nlev )

      if ( LDEBUG ) THEN 
         print*,'kjw MAP_GC2MEM, debug'
         print*,'---------------------------------------'
         WRITE(6,'(14F16.8)') 0d0, PRESSURE_EDGE(:)
         DO LG=1,LLPAR
            WRITE(6,'(14F16.8)') GC_PEDGE(LG), GRIDMAP(LG,:)
         ENDDO
         print*,'---------------------------------------'
      endif

      ! Get GEOS-Chem column from "truth" run to make pseudo-observations
      GC_CH4_NATIVE_OB(:) = 0d0
      GC_CH4_NATIVE_OB(:) = GC_CH4_TRUE_ARRAY(II,JJ,:)


      IF ( LDEBUG ) THEN 
         DO LG = 1, LLPAR
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'L, GC_PCENTER, GC_CH4_NATIVE,' //
     &                      'GC_CH4_NATIVE_OB', 
     &                     LG, GC_PCENTER(LG), GC_CH4_NATIVE(LG), 
     &                     GC_CH4_NATIVE_OB(LG)
         ENDDO
      ENDIF
 299        FORMAT(A50,I3,3F30.12)


      ! Interpolate GEOS-Chem CH4 column and observation to MEM grid
      ! Column in [v/v]
      DO LN = lind, LLMEM
         GC_CH4_onMEM(LN) = 0d0 
         GC_CH4_onMEM_OB(LN) = 0d0 
         DO LG = 1, LLPAR
            GC_CH4_onMEM(LN)    = GC_CH4_onMEM(LN) 
     &           + GRIDMAP(LG,LN) * GC_CH4_NATIVE(LG) 
            GC_CH4_onMEM_OB(LN) = GC_CH4_onMEM_OB(LN) 
     &           + GRIDMAP(LG,LN) * GC_CH4_NATIVE_OB(LG) 
         ENDDO
      ENDDO

      IF ( LDEBUG ) THEN 
         DO LN = lind, LLMEM
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'LN, PRESSURE, GC_CH4,GC_CH4_OB', 
     &                     LN, PRESSURE(LN), GC_CH4_onMEM(LN), 
     &                     GC_CH4_onMEM_OB(LN)
         ENDDO
      ENDIF

         !--------------------------------------------------------------
         ! Apply MEM observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by MEM [molec/cm2]
         !    x_a   = MEM apriori column               [molec/cm2]
         !    x_m   = GC modeled column on MEM grid    [molec/cm2]
         !    A     = MEM averaging kernel 
         !--------------------------------------------------------------

         ! x_m - x_a for model and "observation"
         !    [v/v] --> ln( v/v ) happens here
         DO LN = lind, LLMEM
           GC_CH4_onMEM(LN)   =MAX(GC_CH4_onMEM(LN),   1d-10)
           GC_CH4_onMEM_OB(LN)=MAX(GC_CH4_onMEM_OB(LN),1d-10)
           CH4_PERT(LN)           =LOG( GC_CH4_onMEM(LN) ) - 
     &                             LOG( CH4_PRIOR(II,JJ,LN) )
           CH4_PERT_OB(LN)        =LOG( GC_CH4_onMEM_OB(LN) ) - 
     &                             LOG( CH4_PRIOR(II,JJ,LN) )
         ENDDO

         ! x_a + A_k * ( x_m - x_a ) for model and "observation"
         CH4_HAT(:)=CH4_PERT(:)
         DO LN = lind, LLMEM
            CH4_HAT(LN)    = 0d0
            CH4_HAT_OB(LN) = 0d0
         
            DO LLN = lind, LLMEM
               CH4_HAT(LN) = CH4_HAT(LN) 
     &                    + AVGKERNEL(LN,LLN) * CH4_PERT(LLN)
               CH4_HAT_OB(LN) = CH4_HAT_OB(LN) 
     &                    + AVGKERNEL(LN,LLN) * CH4_PERT_OB(LLN)
            ENDDO
            CH4_HAT(LN)    = CH4_HAT(LN)    + LOG( CH4_PRIOR(II,JJ,LN) )
            CH4_HAT_OB(LN) = CH4_HAT_OB(LN) + LOG( CH4_PRIOR(II,JJ,LN) )

         ENDDO

      IF ( LDEBUG ) THEN 
         DO LN = lind, LLMEM
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'LN, PRESSURE, CH4_HAT,CH4_HAT_OB', 
     &                     LN, PRESSURE(LN), exp(CH4_HAT(LN)), 
     &                     exp(CH4_HAT_OB(LN))
            WRITE(6,299) 'LN, CH4_HAT,GC_CH4_onMEM,CH4_PRIOR', 
     &                     LN, exp(CH4_HAT(LN)), 
     &                     GC_CH4_onMEM(LN), CH4_PRIOR(II,JJ,LN)
         ENDDO
      ENDIF


      ! Convert vertical profiles from [ln(vmr)] --> [vmr] before 
      !   calculating XCH4
      CH4_HAT_EXP    = EXP(CH4_HAT)
      CH4_HAT_OB_EXP = EXP(CH4_HAT_OB)

      IF ( LDEBUG ) THEN 
         DO LN = lind, LLMEM
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'CH4_HAT_EXP, CH4_HAT_EXP, CH4_HAT_OB_EXP', 
     &        LN, CH4_HAT_EXP(LN), CH4_HAT_EXP(LN), CH4_HAT_OB_EXP(LN)
         ENDDO
      ENDIF

      ! ---- Calculate XCH4 [v/v] from CH4_HAT [v/v] and CH4_HAT_OB [v/v]
      XCH4_HAT    = 0d0
      XCH4_HAT_OB = 0d0

      ! Calculate weight of each vertical level on MEM grid for averaging
      !   levels to get XCH4. Weight by # molecules / verical level, which is
      !   proportional to pressure difference between upper and lower bounds
      !   of each box.
      DO LN=lind, LLMEM

         ! If ground level, average with same weight as if it were 1st atm level
         IF ( LN .EQ. lind ) THEN
            PUP = PRESSURE_EDGE(LN+1)
            PLO = PRESSURE_EDGE(LN  )
         ELSE
            PUP = PRESSURE_EDGE(LN  )
            PLO = PRESSURE_EDGE(LN-1)
         ENDIF

         Xweight(LN) = PLO - PUP
      ENDDO

      !Normalize so that SUM(Xweight) = 1
      sumxweight = SUM( Xweight(:) )
      DO LN=lind,LLMEM
         Xweight(LN) = Xweight(LN) / sumxweight
      ENDDO


      IF ( LDEBUG ) THEN 
         DO LN=lind,LLMEM
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'Xweight', 
     &             LN,  Xweight(LN),  Xweight(LN),  Xweight(LN) 
         ENDDO
      ENDIF


      ! Calculate weighted average of CH4_HAT and CH4_HAT_OB
      DO LN=lind, LLMEM
         XCH4_HAT    = XCH4_HAT    + Xweight(LN) * CH4_HAT_EXP(LN)
         XCH4_HAT_OB = XCH4_HAT_OB + Xweight(LN) * CH4_HAT_OB_EXP(LN)
      ENDDO

      IF ( LDEBUG ) THEN 
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'XCH4_HAT, XCH4_HAT, XCH4_HAT_OB', 
     &             1,  XCH4_HAT, XCH4_HAT, XCH4_HAT_OB
      ENDIF


      ! Create super observation by adding random error
      ! to XCH4_HAT_OB
      ! SUPER_ERR is 1d-9 * 8ppb * N(0,1) / sqrt(nobs)  [v/v]
      !     where 8ppb is expected error on a single XCH4 measurement
      !           N(0,1) is a random number of mean 0, standard deviation 1
      !           nobs is the number of observations merged to form super obs
      ! Add error of each observation that makes up super-obs. Do this to
      !   preserve error structure across different resolutions.
      SUPER_ERR_EXPECTED = 1d-9 * XCH4_ERR / SQRT( XPLX(nobs) )

      ! Multiply expected error of super-observation by 
      !   prescribed random number with mean 0, standard deviation 1
      SUPER_ERR = SUPER_ERR_EXPECTED * RANDOM_GRID( II, JJ )

      ! Add random error to super-observation
      XCH4_HAT_OB = XCH4_HAT_OB + SUPER_ERR         ! add error [v/v]


      IF ( LDEBUG ) THEN 
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'XCH4_ERR, SUPER_ERR, nobs', 
     &               1, XCH4_ERR, SUPER_ERR, XPLX(nobs)
      ENDIF
      IF ( LDEBUG ) THEN 
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'XCH4_HAT_OB, XCH4_HAT_OB, XCH4_HAT_OB', 
     &             1, XCH4_HAT_OB, XCH4_HAT_OB, XCH4_HAT_OB
      ENDIF


!         ! Add error to create super-observation
!         !   nobs  - number observations in this grid box
!         !   boxno - box number processed during this day
!         ! Magnitude of error in super observation
!         SUPER_ERR = 1d0 / SQRT( REAL(nobs) ) * RANDNUM( NB )
!
!         ! Print information about this grid box to file
!         IF ( SUM(PERT(:)) .EQ. 0. ) THEN
!            WRITE(116,212) 'II = ', REAL(40)
!            WRITE(116,212) 'JJ = ', REAL(JJPAR-10)
!            WRITE(116,212) 'nobs = ', REAL(nobs)
!            WRITE(116,212) 'RANDOM(NB) = ', RANDNUM( NB )
!            WRITE(116,212) 'SUPER_ERR  = ', SUPER_ERR
!         ENDIF
! 212     FORMAT(A12,F22.6)
!
!         ! Calculate sqrt( obserror ) <-- magnitude of error in 1 observation
!         DO LN  = lind, LLMEM
!         DO LLN = lind, LLMEM
!            SIGN(LN,LLN) = OBSERROR(LN,LLN) / ABS( OBSERROR(LN,LLN) )
!            OBSERROR_SQRT(LN,LLN) = SIGN( LN, LLN ) * 
!     &                  SQRT( ABS( OBSERROR(LN,LLN) ) )
!         ENDDO
!         ENDDO
!         print*,'maxval/minval( SIGN ) ', maxval(sign),minval(sign)
!
!         ! Create super observation
!         CH4_HAT_OB_werr(:) = 0d0
!         DO LN  = lind, LLMEM
!            CH4_HAT_OB_werr(LN) = CH4_HAT_OB(LN)
!            DO LLN = lind, LLMEM
!               CH4_HAT_OB_werr(LN) = CH4_HAT_OB_werr(LN) + 
!     &              CH4_HAT_OB(LN) * SUPER_ERR * OBSERROR_SQRT(LN,LLN)
!            ENDDO
!         ENDDO
!
!         IF ( LDEBUG ) THEN 
!            DO LN  = lind, LLMEM
!            DO LLN = lind, LLMEM
!               dummyerr(LN) = CH4_HAT_OB(LN) * OBSERROR_SQRT(LN,LN)
!            ENDDO
!            ENDDO
!            WRITE(6,'(A16,13F18.9)') ,'dummyerr   = ', dummyerr(:)
!            WRITE(6,'(A16,13F18.9)') ,'CH4_HAT_OB = ',exp(CH4_HAT_OB(:))
!            WRITE(6,'(A16,13F18.9)') ,'PERT = ',
!     &                  exp(CH4_HAT_OB(:)+dummyerr(:))
!         ENDIF
!               
!
!
!         ! Scale observation error covariance matrix to nobs
!         DO LN  = lind, LLMEM
!         DO LLN = lind, LLMEM
!            OBSERROR_INV_SUPER(LN,LLN) = 
!     &             OBSERROR_INV(LN,LLN) * REAL(nobs)
!         ENDDO
!         ENDDO
!
!      IF ( LDEBUG ) THEN 
!         DO LN = lind, LLMEM
!            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
!            WRITE(6,299) 'LN, PRESSURE, CH4_HAT_OB, CH4_HAT_OB_werr', 
!     &                     LN, PRESSURE(LN), 
!     &                     exp(CH4_HAT_OB(LN)), exp(CH4_HAT_OB_werr(LN))
!         ENDDO
!      ENDIF
!
         !-------------------------------------------------------------
         ! Calculate cost function, given S is observation error covariance matrix
         !     Sobs = 1x1 array [ ln(vmr)^2 ]
         ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------

         ! Initialize values to be safe
         DIFF        = 0d0
         FORCE       = 0d0

         ! Calculate difference between modeled and observed profile
         DIFF        = XCH4_HAT - XCH4_HAT_OB

         ! Calculate adjoint forcing: 2 * DIFF^T * S_{obs}^{-1}
         !         and cost function: DIFF^T * S_{obs}^{-1} * DIFF
         S_obs_inv   = 1d0 / (SUPER_ERR**2)
         FORCE = 2 * DIFF * S_obs_inv
         NEW_COST(NB) = 0.5d0 * DIFF * FORCE


      IF ( LDEBUG ) THEN 
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'DIFF, FORCE, NEW_COST(NB)', 
     &               1, 1d9*DIFF, 1d9*FORCE, NEW_COST(NB)
      ENDIF

!         ! Initialize values to be safe
!         DIFF(:)         = 0d0
!         FORCE(:)        = 0d0
!
!         ! Calculate difference between modeled and observed profile
!         DO LN = lind, LLMEM
!            DIFF(LN) = CH4_HAT(LN) - CH4_HAT_OB_werr(LN)
!         ENDDO 
!
!            ! Print information about this grid box to file
!         DO LN=lind,LLMEM
!            IF ( LDEBUG ) THEN
!            WRITE(116,213) 'PRESSURE(LN),CH4_HAT(LN),' // 
!     &                     'CH4_HAT_OB(LN),CH4_PRIOR(LN)', 
!     &           PRESSURE( LN ), 1d9 * exp(CH4_HAT(LN)), 
!     &        1d9 * exp(CH4_HAT_OB_werr(LN)), 1d9 * CH4_PRIOR(II,JJ,LN)
!            ENDIF
!         ENDDO
! 213     FORMAT(A60,4F22.6)
!
!
!         ! Calculate adjoint forcing: 2 * DIFF^T * S_{obs}^{-1}
!         !         and cost function: DIFF^T * S_{obs}^{-1} * DIFF
!         DO LN = lind, LLMEM
!            DO LLN = lind, LLMEM
!               FORCE(LN)  = FORCE(LN) + 
!     &              2d0 * OBSERROR_INV_SUPER(LN,LLN) * DIFF(LLN)
!            ENDDO
!            NEW_COST(NB) = NEW_COST(NB) + 0.5*DIFF(LN)*FORCE(LN)
!         ENDDO
!
!
         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------


!         ! The adjoint forcing is 2 * S_{obs}^{-1} * DIFF = FORCE
!         DIFF_ADJ(:) = FORCE(:) 
!
!         ! Adjoint of GEOS-Chem - Observation difference
!         CH4_HAT_ADJ(:) = DIFF_ADJ(:)
!
!         ! Adjoint of adding random error to observation
!         DO LN=lind,LLMEM
!            CH4_HAT_ADJ(LN) = 0d0
!
!            DO LLN=lind,LLMEM
!               CH4_HAT_ADJ(LN) = CH4_HAT_ADJ(LN) + 
!     &           CH4_HAT_ADJ(LLN) * SUPER_ERR * OBSERROR(LLN,LN)
!            ENDDO
!         ENDDO

         ! The adjoint forcing is 2 * S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE

         ! Adjoint of GEOS-Chem - Observation difference
         XCH4_HAT_ADJ = DIFF_ADJ

         ! Adjoint of CH4_HAT_EXP --> XCH4_HAT
         DO LN=lind, LLMEM
            CH4_HAT_EXP_ADJ(LN) = XCH4_HAT_ADJ * Xweight(LN)
         ENDDO

         ! Adjoint of CH4_HAT --> CH4_HAT_EXP
         DO LN=lind, LLMEM
            CH4_HAT_ADJ(LN) = CH4_HAT_EXP_ADJ(LN) * CH4_HAT_EXP(LN)
         ENDDO


         ! Adjoint of MEM observation operator
         CH4_PERT_ADJ(:) = 0D0
         DO LN=lind,LLMEM
         DO LLN=lind,LLMEM
            CH4_PERT_ADJ(LN) = CH4_PERT_ADJ(LN) + 
     &              AVGKERNEL(LLN,LN) * CH4_HAT_ADJ(LLN)
         ENDDO
         ENDDO

         ! Adjoint of x_m - x_a
         DO LN = lind, LLMEM 
           ! fwd code:
           !GC_CH4(LN)   = MAX(GC_CH4(LN), 1d-10)
           !CH4_PERT(LN) = LOG(GC_CH4(LN)) - LOG(PRIOR(LN))
           ! adj code:
           IF ( GC_CH4_onMEM(LN) > 1d-10 ) THEN 
              GC_CH4_onMEM_ADJ(LN) = 1d0 / GC_CH4_onMEM(LN) * 
     &                                         CH4_PERT_ADJ(LN)
           ELSE 
              GC_CH4_onMEM_ADJ(LN) = 1d0 / 1d-10 * CH4_PERT_ADJ(LN)
           ENDIF
         ENDDO

      IF ( LDEBUG ) THEN 
         DO LN=lind,LLMEM
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'GC_CH4_onMEM_ADJ, CH4_PERT_ADJ, CH4_HAT_ADJ', 
     &             LN,  GC_CH4_onMEM_ADJ(LN), CH4_PERT_ADJ(LN), 
     &                  CH4_HAT_ADJ(LN)
         ENDDO
      ENDIF


         ! Adjoint of interpolation
         DO LN=lind,LLMEM
         DO LG=1,LLPAR
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) + 
     &           GRIDMAP(LG,LN) * GC_CH4_onMEM_ADJ(LN)
         ENDDO
         ENDDO


         ! Adjoint of unit conversion
         DO LG=1,LLPAR
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG)
     &                 * XNUMOL(1) / ( XNUMOLAIR * GC_AD(LG) )
         ENDDO


         ! Pass adjoing forcing back to adjoint tracer array
         DO LG=1,LLPAR
            ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) * CHK_STT(II,JJ,LG,1)
         ENDDO

      IF ( LDEBUG ) THEN 
         DO LG=1,LLPAR
            WRITE(6,*)   'kjw DEBUG IN CALC_MEM_FORCE'
            WRITE(6,299) 'GC_CH4_NATIVE_ADJ, ADJ(LG)', 
     &             LN,  GC_CH4_NATIVE_ADJ(LG), ADJ(LG),1
         ENDDO
      ENDIF

         ! Update cost function 
         COST_FUNC_A = COST_FUNC_A + NEW_COST(NB)

         ! Only debug on first pass through routine
         LDEBUG = .FALSE.



      ! Return to calling program
      END SUBROUTINE CALC_MEM_CH4_FORCE_FD


!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( GC_PEDGE, MEM_PEDGE, nlev  )
     &         RESULT      ( M )
!
!******************************************************************************
!  Function GET_INTMAP creates the matrix that places GEOS-Chem column methane 
!     [molec/cm2] onto the 13-level pressure grid used by theoretical instrument, M.
!         GC[1x47] * M[47x13] = MEM[1x13]           (kjw, 7/21/11)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) GC_PEDGE   (TYPE (XPLEX)) : LLPAR bottom pressure edges of GEOS-Chem column
!  (2 ) SCIA_PEDGE (TYPE (XPLEX)) : LLMEM upper pressure edges of MEM column (except 
!                                first entry, which is surface pressure)
!  (3 ) nlev       (TYPE (XPLEX)) : Number of MEM pressure levels to use
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) M          (TYPE (XPLEX)) : Interpolation matrix that maps GEOS-Chem to MEM grid
!     
!  NOTES:
!  (1 ) Based on GET_INTMAP in scia_ch4_mod.f
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX)             :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)             :: MEM_PEDGE(LLMEM)
      INTEGER            :: nlev
 
      ! Return value 
      TYPE (XPLEX)             :: M(LLPAR,LLMEM)

      ! Local variables 
      INTEGER  :: LGC, LTM, LS, LG, LN, LIND
      TYPE (XPLEX)   :: DIFF, DELTA_SURFP
      TYPE (XPLEX)   :: GUP, GLO, NUP, NLO
      TYPE (XPLEX)   :: column_total(LLMEM)
      LOGICAL, SAVE :: LDEBUG = .TRUE.

      !=================================================================
      ! GET_INTMAP begins here!
      !=================================================================

      ! Initialize output
      M(:,:) = 0D0 

      ! Minimum MEM vertical level to use
      lind = LLMEM + 1 - nlev

      ! Loop over each pressure level of GEOS-Chem and MEM grids
      DO LG=1,LLPAR

         ! Get upper and lower pressure edges of GEOS-Chem box
         IF ( LG .EQ. LLPAR ) THEN 
            GUP = 0d0
            GLO = GC_PEDGE( LG   )
         ELSE
            GUP = GC_PEDGE( LG+1 )
            GLO = GC_PEDGE( LG   )
         ENDIF         

         DO LN=lind,LLMEM

            ! Get top and bottom pressures of MEM box
            ! If processing first MEM level, this is surface level, so 
            !    bottom and top of box are same level. Set "bottom" of
            !    MEM box to GEOS-Chem surface pressure so that MEM surface
            !    box avgs GEOS-Chem values between GEOS-Chem surface and 
            !    MEM surface pressures.
            ! GC surface pressure is always > MEM surface pressure because
            !    we chop off lowermost MEM levels if it is not
            IF ( LN .EQ. lind ) THEN 
               NUP = MEM_PEDGE( LN   )
               NLO = GC_PEDGE(  LG   )
            ELSE
               NUP = MEM_PEDGE( LN   )
               NLO = MEM_PEDGE( LN-1 )
            ENDIF

            ! If both GEOS-Chem edges are within the MEM box, map value = 1
            IF ( ( GUP .gt. NUP ) .AND. ( GLO .lt. NLO ) ) THEN
               M(LG,LN) = 1
            ENDIF

            ! If both GEOS-Chem stradles a MEM pressure level, interpolate
            IF ( ( GUP .lt. NUP ) .AND. ( GLO .gt. NUP ) ) THEN
               DIFF       = GLO - GUP
               M(LG,LN+1) = ( NUP - GUP ) / DIFF
               M(LG,LN  ) = ( GLO - NUP ) / DIFF
            ENDIF

         ENDDO
      ENDDO

      ! Add value for uppermost GEOS-Chem grid box
      M(LLPAR,LLMEM) = 1


      ! Correct for case in which GEOS-Chem pressure is higher than MEM
      IF ( GC_PEDGE(1) .GT. MEM_PEDGE(1) ) THEN

         ! If any part of GEOS-Chem box are under MEM_PEDGE(1), let 
         !   this GEOS-Chem grid box contribute to the observation because 
         !   MEM and GEOS-Chem should have same surface pressure. map value = 1
         DO LG=1,LLPAR-1

            ! If GEOS-Chem box entirely below MEM surface pressure
            IF ( ( GC_PEDGE(LG)   .GT. MEM_PEDGE(1) ) .AND.    
     &           ( GC_PEDGE(LG+1) .GT. MEM_PEDGE(1) ) ) THEN 
               M(LG,1) = 1
            ENDIF

            ! If GEOS-Chem box straddles MEM surface pressure
            IF ( ( GC_PEDGE(LG)   .GT. MEM_PEDGE(1) ) .AND.    
     &           ( GC_PEDGE(LG+1) .LT. MEM_PEDGE(1) ) ) THEN 
               DIFF = GC_PEDGE(LG) - GC_PEDGE( LG+1 )
               M(LG,1) = ( MEM_PEDGE(1) - GC_PEDGE(LG+1) ) / DIFF
            ENDIF
            
         ENDDO
      ENDIF


      ! Correct for case in which GEOS-Chem surface pressure is within 2nd MEM
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. MEM_PEDGE(2) ) THEN
         M(1,1) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 3rd MEM
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. MEM_PEDGE(3) ) THEN
         M(1,2) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 4th MEM
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. MEM_PEDGE(4) ) THEN
         M(1,3) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 5th MEM
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. MEM_PEDGE(5) ) THEN
         M(1,4) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 6th MEM
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. MEM_PEDGE(6) ) THEN
         M(1,5) = 0.
      ENDIF

      ! Normalize each column of M to 1 so that we are not creating any molecules
      ! when mapping from GEOS-Chem to MEM grids.

      ! DO NOT do this since we are mapping molc/cm2, not 
      ! Initialize to be safe and calculate column total
      column_total(:) = 0d0
      column_total(:) = xplx(SUM( M%r, 1 ),SUM(M%i,1))

      ! Normalize columns to column_total
      DO LN=1,LLMEM
         IF ( column_total(LN) .EQ. 0. ) CYCLE
         M(:,LN) = M(:,LN) / column_total(LN)
      ENDDO


      !if ( LDEBUG ) THEN 
      !   print*,'kjw GET_INTMAP, debug'
      !   print*,'---------------------------------------'
      !   WRITE(6,'(14F16.8)') 0d0, MEM_PEDGE(:)
      !   DO LG=1,LLPAR
      !      WRITE(6,'(14F16.8)') GC_PEDGE(LG), M(LG,:)
      !   ENDDO
      !   print*,'---------------------------------------'
      !   LDEBUG = .FALSE.
      !endif

      ! Return to calling program
      END FUNCTION GET_INTMAP


!-----------------------------------------------------------------------------



      END MODULE MEM_CH4_MOD
