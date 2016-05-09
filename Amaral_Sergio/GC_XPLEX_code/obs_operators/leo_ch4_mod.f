!$Id: leo_ch4_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE LEO_CH4_MOD
!  
!******************************************************************************
!  Module LEO_CH4_MOD for CH4 observations. 
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
      INTEGER, PARAMETER         :: LLLEO   = 13
      INTEGER, PARAMETER         :: MAXLEO  = 639059


      ! Record to store information about the new instrument
      TYPE (XPLEX)                     :: AVGKERNEL(    LLLEO, LLLEO )
      TYPE (XPLEX)                     :: OBSERROR(     LLLEO, LLLEO )
      TYPE (XPLEX)                     :: OBSERROR_INV( LLLEO, LLLEO )
      TYPE (XPLEX)                     :: TOTERROR_INV( LLLEO, LLLEO )
      TYPE (XPLEX)                     :: PRESSURE( LLLEO )
      TYPE (XPLEX)                     :: PRESSURE_EDGE( LLLEO )
      TYPE (XPLEX)                     :: RANDNUM( MAXLEO )


      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_LEO_INFO
!
!******************************************************************************
!  Subroutine READ_LEO_INFO reads and stores information about the new
!  instrument, specifically AK, pressure levels and error covariance matrices.
!    (kjw, 07/24/11) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME   (CHAR) : LEO filename to read
!
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE FILE_MOD,              ONLY : IOERROR
      USE TIME_MOD,              ONLY : GET_NYMD

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

      ! Loop indexes, and error handling.
      INTEGER                         :: IOS, IU_IN



      !=================================================================
      ! READ_LEO_CH4_OBS begins here!
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
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/leo/' ) //
     &                   'data/' // TRIM( 'leo_AK.txt' )
         WRITE(6,*) '    - READ_LEO_AK: reading file: ', 
     &                        TRIM(READ_FILENAME)


         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_avg_kernel:1' )

         ! Read File and save info into module variable AVGKERNEL(:,:)
         DO LN=1,LLLEO
            READ( IU_IN, '(13F12.6)', IOSTAT=IOS ) AVGKERNEL(LN,:)

            ! IO status
            IF ( IOS < 0 ) THEN
               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
               WRITE( 6, '(a)' ) 'STOP in READ_LEO_CH4'
            ENDIF
            IF ( IOS > 0 ) THEN
               CALL IOERROR(IOS, IU_IN, 'read_avg_kernel:2')
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IU_IN )
 

         ! ------ Observation Error Covariance Matrix ------
         ! Filename to read
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/leo/' ) //
     &                   'data/' // TRIM( 'leo_obs_error.txt' )
         WRITE(6,*) '    - READ_LEO_OBSERROR: reading file: ', 
     &                        TRIM(READ_FILENAME)


         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_obs_error:1' )

         ! Read File and save info into module variable OBSERROR(:,:)
         DO LN=1,LLLEO
            READ( IU_IN, '(13F18.12)', IOSTAT=IOS ) OBSERROR(LN,:)

            ! IO status
            IF ( IOS < 0 ) THEN
               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
               WRITE( 6, '(a)' ) 'STOP in READ_LEO_CH4'
            ENDIF
            IF ( IOS > 0 ) THEN
               CALL IOERROR(IOS, IU_IN, 'read_obs_error:2')
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IU_IN )


         ! ------ Inverse of Observation Error Covariance Matrix ------
         ! Filename to read
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/leo/' ) //
     &                   'data/' // TRIM( 'leo_obs_error_inv.txt' )
         WRITE(6,*) '    - READ_LEO_OBSERROR_INV: reading file: ',
     &                        TRIM(READ_FILENAME)


         ! Open file
         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
     &         STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_obs_error:1' )

         ! Read File and save info into module variable OBSERROR_INV(:,:)
         DO LN=1,LLLEO
            READ( IU_IN, '(13F18.6)', IOSTAT=IOS ) OBSERROR_INV(LN,:)

            ! IO status
            IF ( IOS < 0 ) THEN
               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
               WRITE( 6, '(a)' ) 'STOP in READ_LEO_CH4'
            ENDIF
            IF ( IOS > 0 ) THEN
               CALL IOERROR(IOS, IU_IN, 'read_obs_error:2')
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IU_IN )


!         ! ------ Total Error Covariance Matrix ------
!         ! Filename to read
!         READ_FILENAME = TRIM( '/home/kjw/new_satellites/leo/' ) //
!     &                 'data/' // TRIM( 'leo_total_error_inv.txt' )
!         WRITE(6,*) '    - READ_LEO_TOTERROR: reading file: ', 
!     &                        TRIM(READ_FILENAME)
!
!
!         ! Open file
!         OPEN( IU_IN,        FILE=TRIM( READ_FILENAME ), 
!     &         STATUS='OLD', IOSTAT=IOS )
!         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_IN, 'read_tot_error:1' )
!
!         ! Read File and save info into module variable OBSERROR(:,:)
!         DO LN=1,LLLEO
!            READ( IU_IN, '(13F18.12)', IOSTAT=IOS ) TOTERROR_INV(LN,:)
!
!            ! IO status
!            IF ( IOS < 0 ) THEN
!               WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
!               WRITE( 6, '(a)' ) 'STOP in READ_LEO_CH4'
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
         READ_FILENAME = TRIM( '/home/kjw/new_satellites/leo/' ) //
     &                   'data/' // TRIM( 'leo_pressure.txt' )
         WRITE(6,*) '    - READ_LEO_PRESSURE: reading file: ', 
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
            WRITE( 6, '(a)' ) 'STOP in READ_LEO_CH4'
         ENDIF
         IF ( IOS > 0 ) THEN
            CALL IOERROR(IOS, IU_IN, 'read_pressure:2')
         ENDIF

         ! Close file
         CLOSE( IU_IN )


         ! ------ Pressure Edges ------
         ! By finite difference on log(pressure) grid
         PRESSURE_EDGE(1) = PRESSURE(1)
         PRESSURE_EDGE(LLLEO) = 0.
         DO LN=2,LLLEO-1
            PRESSURE_EDGE(LN) = exp( log(pressure(LN+1)) + 
     &          ( log(PRESSURE(LN)) - log(PRESSURE(LN+1)) ) / 2.  )
         ENDDO

    
      ! Return to calling program
      END SUBROUTINE READ_LEO_INFO
!------------------------------------------------------------------------------


      SUBROUTINE CALC_LEO_CH4_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_LEO_CH4_FORCE calculates the adjoint forcing from the LEO
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
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER, SAVE               :: NT   ! # observations processed this day
      INTEGER                     :: LG, LN, LLN, II, JJ, NB, JMIN, OB
      INTEGER                     :: nlev, lind, IU_IN
      INTEGER                     :: nboxes, nobs
      INTEGER                     :: NTSTART, NTSTOP, NTh
      INTEGER, SAVE               :: NTT
      TYPE (XPLEX)               :: GC_CH4_TRUE_ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: CH4_PRIOR(IIPAR,JJPAR,LLLEO)
      TYPE (XPLEX)                    :: DUMMY_PRIOR(IIPAR,JJPAR,LLLEO)
      TYPE (XPLEX)                      :: DUMMY_TRUE(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: GC_PCENTER(LLPAR)
      TYPE (XPLEX)                      :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)                      :: GC_AD(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_OB(LLPAR)
      TYPE (XPLEX)                      :: thispcen(LLPAR)
      TYPE (XPLEX)                      :: thispedg(LLPAR)
      TYPE (XPLEX)                      :: thisad(LLPAR)
      TYPE (XPLEX)                      :: thisch4(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_onLEO(LLLEO)
      TYPE (XPLEX)                      :: GC_CH4_onLEO_OB(LLLEO)
      TYPE (XPLEX)                      :: GRIDMAP(LLPAR,LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_OB(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_ADJ(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_werr(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_werr_ADJ(LLLEO)
      TYPE (XPLEX)                      :: CH4_PERT(LLLEO)
      TYPE (XPLEX)                      :: CH4_PERT_OB(LLLEO)
      TYPE (XPLEX)                      :: CH4_PERT_ADJ(LLLEO)
      TYPE (XPLEX)                      :: frac, frac_total
      TYPE (XPLEX)                      :: latmin, Jfrac_min, Jfrac
      TYPE (XPLEX)                      :: box_area, cloud_frac
      TYPE (XPLEX)                      :: mass_air, mole_air, mole_ch4
      TYPE (XPLEX)                      :: LHS, RHS, GC_XCH4, XTAU
      TYPE (XPLEX)                      :: DIFF(LLLEO)
      TYPE (XPLEX)                      :: FORCE(LLLEO)
      TYPE (XPLEX)                      :: DIFF_ADJ(LLLEO)
      TYPE (XPLEX)                      :: thisforce(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_onLEO_ADJ(LLLEO)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: NEW_COST(MAXLEO)
      TYPE (XPLEX)                      :: OLD_COST
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      LOGICAL, SAVE               :: DO_FDTEST = .TRUE.
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME
      CHARACTER(LEN=255)          :: FILENAME_OBS

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
      ! CALC_LEO_CH4_FORCE begins here!
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
         CALL READ_LEO_INFO

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



      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Read "TRUE" state for this time step [kg/box]
      GC_CH4_TRUE_ARRAY(:,:,:) = 0d0
      FILENAME_OBS = '/home/kjw/GEOS-Chem/gcadj_std/runs/v8-02-01/' //
     &               'ch4/leo/' // GET_RES_EXT() // '/adjtmp/'  //
     &               'gctm.obs.YYYYMMDD.hhmm'
      CALL EXPAND_DATE( FILENAME_OBS, GET_NYMD(), GET_NHMS() )
      !FILENAME_OBS = TRIM( ADJTMP_DIR   )  //  TRIM( FILENAME_OBS )
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

      ! Read a priori vertical profiles from file
      FILENAME = '/home/kjw/new_satellites/leo/data/' //
     &           'leo_prior.' // GET_RES_EXT() // '.bpch'
      XTAU = GET_TAU0( 1, 1, 1985 )
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,    
     &                 XTAU,      IIPAR,     JJPAR, 
     &                 LLLEO, DUMMY_PRIOR, QUIET=.TRUE.  )
      CH4_PRIOR(:,:,:) = DUMMY_PRIOR(:,:,:)


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
         NT = 0

         ! ------ Random Numbers ------
         ! Open and read random number file. mean = 0, stddev = 1
         FILENAME = '/home/kjw/new_satellites/leo/data/' //
     &              'randnums/random.YYYYMMDD.txt'
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), 0 )
         OPEN( IU_IN, FILE=TRIM( FILENAME ),  STATUS='UNKNOWN',
     &        IOSTAT=IOS,  FORM='FORMATTED',  ACCESS='SEQUENTIAL' )
         DO LG=1,MAXLEO
            READ(IU_IN,'(F13.6)') RANDNUM(LG)
         ENDDO
         CLOSE(IU_IN)

      ENDIF

      ! Begin counter for number of observations processed this hour
      NTh = 0

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

      print*, '    - CALC_LEO_CH4_FORCE  ', GET_NYMD(), GET_NHMS()


      ! Loop over each grid box north of the minimum latitude
      !   1. Determine number of observations in the current grid box
      !   2. Make obseravations
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
            GC_CH4_onLEO(:)    = 0d0
            GC_CH4_onLEO_OB(:) = 0d0


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
               GC_CH4_NATIVE(LG) = ( CHK_STT(II,JJ,LG,1 )
     &                    * XNUMOL(1) ) / ( GC_AD(LG) * XNUMOLAIR )

            ENDDO


            ! Number of vertical levels to use in these observations
            !   Chop off lowermost levels if 
            !     GEOS-Chem surface pressure < LEO pressure levels
            nlev = count( (PRESSURE_EDGE%r) .LT. (GC_PEDGE(1)%r) )
            IF ( nlev .LT. 13 )  nlev = nlev + 1
            lind = LLLEO + 1 - nlev ! minimum vertical index on LEO grid


            ! Get interpolation matrix that maps GEOS-Chem to LEO grid
            GRIDMAP(1:LLPAR, 1:LLLEO) = 
     &                   GET_INTMAP( GC_PEDGE, PRESSURE_EDGE, nlev )

            ! Get GEOS-Chem column from "truth" run to make pseudo-observations
            GC_CH4_NATIVE_OB(:) = 0d0
            GC_CH4_NATIVE_OB(:) = GC_CH4_TRUE_ARRAY(II,JJ,:)

            ! Interpolate GEOS-Chem CH4 column and observation to LEO grid
            ! Column in [v/v]
            DO LN = lind, LLLEO
               GC_CH4_onLEO(LN) = 0d0 
               GC_CH4_onLEO_OB(LN) = 0d0 
               DO LG = 1, LLPAR
                  GC_CH4_onLEO(LN)    = GC_CH4_onLEO(LN) 
     &                 + GRIDMAP(LG,LN) * GC_CH4_NATIVE(LG) 
                  GC_CH4_onLEO_OB(LN) = GC_CH4_onLEO_OB(LN) 
     &                 + GRIDMAP(LG,LN) * GC_CH4_NATIVE_OB(LG) 
               ENDDO
            ENDDO


            !--------------------------------------------------------------
            ! Apply LEO observation operator
            !
            !   x_hat = x_a + A_k ( x_m - x_a ) 
            !  
            !  where  
            !    x_hat = GC modeled column as seen by LEO [molec/cm2]
            !    x_a   = LEO apriori column               [molec/cm2]
            !    x_m   = GC modeled column on LEO grid    [molec/cm2]
            !    A     = LEO averaging kernel 
            !--------------------------------------------------------------
            
            ! x_m - x_a for model and "observation"
            !    [v/v] --> ln( v/v ) happens here
            DO LN = lind, LLLEO
              GC_CH4_onLEO(LN)   =MAX(GC_CH4_onLEO(LN),   1d-10)
              GC_CH4_onLEO_OB(LN)=MAX(GC_CH4_onLEO_OB(LN),1d-10)
              CH4_PERT(LN)           =LOG( GC_CH4_onLEO(LN) ) - 
     &                                LOG( CH4_PRIOR(II,JJ,LN) )
              CH4_PERT_OB(LN)        =LOG( GC_CH4_onLEO_OB(LN) ) - 
     &                                LOG( CH4_PRIOR(II,JJ,LN) )
            ENDDO

            ! x_a + A_k * ( x_m - x_a ) for model and "observation"
            DO LN = lind, LLLEO
               CH4_HAT(LN)    = 0d0
               CH4_HAT_OB(LN) = 0d0
            
               DO LLN = lind, LLLEO
                  CH4_HAT(LN) = CH4_HAT(LN) 
     &                       + AVGKERNEL(LN,LLN) * CH4_PERT(LLN)
                  CH4_HAT_OB(LN) = CH4_HAT_OB(LN) 
     &                       + AVGKERNEL(LN,LLN) * CH4_PERT_OB(LLN)
               ENDDO
               CH4_HAT(LN)   = CH4_HAT(LN)   +LOG( CH4_PRIOR(II,JJ,LN) )
               CH4_HAT_OB(LN)= CH4_HAT_OB(LN)+LOG( CH4_PRIOR(II,JJ,LN) )
            
            ENDDO


            ! Loop over number of observations in this grid box
            DO OB=1,NOBS

               ! Increment number of observations
               NTh = NTh + 1   ! processed this hour
               NT  = NT  + 1   ! processed today
               NTT = NTT + 1   ! processed total

               !print*, '    - CALC_LEO_CH4_FORCE ', OB, ' of ',NOBS


               ! For safety, initialize these up to LLLEO
               CH4_HAT_werr(:) = 0d0
               DIFF(:)         = 0d0
               FORCE(:)        = 0d0
               NEW_COST(:)     = 0d0

               ! Add random error to this observation
               DO LN  = lind, LLLEO

                  CH4_HAT_werr(LN) = CH4_HAT(LN)

                  DO LLN = lind, LLLEO
                     CH4_HAT_werr(LN) = CH4_HAT_werr(LN) + 
     &                    CH4_HAT(LN) * RANDNUM(NT) * OBSERROR(LN,LLN)
                  ENDDO
               ENDDO


               !--------------------------------------------------------------
               ! Calculate cost function, given S is observation error covariance matrix
               !     Sobs = 1x1 array [ (molec/cm2) ^2 ]
               ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
               !--------------------------------------------------------------

               ! Calculate difference between modeled and observed profile
               DO LN = lind, LLLEO
                  DIFF(LN) = CH4_HAT_werr(LN) - CH4_HAT_OB(LN)
               ENDDO 

               ! Calculate adjoint forcing: 2 * DIFF^T * S_{obs}^{-1}
               !         and cost function: DIFF^T * S_{obs}^{-1} * DIFF
               DO LN = lind, LLLEO
                  DO LLN = lind, LLLEO
                     FORCE(LN)  = FORCE(LN) + 
     &                    2d0 * OBSERROR_INV(LN,LLN) * DIFF(LLN)
                  ENDDO
                  NEW_COST(LN) = NEW_COST(LN) + 0.5*DIFF(LN)*FORCE(LN)
               ENDDO



         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

         ! dkh debug
!         print*, 'DIFF , FORCE, Sobs ' 
!         WRITE(6,102) (DIFF, FORCE, Sobs)
! 102  FORMAT(1X,d14.6,1X,d14.6)


         ! The adjoint forcing is 2 * S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ(:) = 2. * FORCE(:)

         ! Adjoint of GEOS-Chem - Observation difference
         CH4_HAT_werr_ADJ(:) = DIFF_ADJ(:)

         ! Adjoint of adding random error to observation
         DO LN=lind,LLLEO
            CH4_HAT_ADJ(LN) = 0d0

            DO LLN=lind,LLLEO
               CH4_HAT_ADJ(LN) = CH4_HAT_ADJ(LN) + 
     &              CH4_HAT_ADJ(LLN) * RANDNUM(NT) * OBSERROR(LLN,LN)
            ENDDO
         ENDDO

         ! Adjoint of LEO observation operator
         DO LN=lind,LLLEO
            CH4_PERT_ADJ(LN) = 0D0

            DO LLN=lind,LLLEO
               CH4_PERT_ADJ(LN) = CH4_PERT_ADJ(LN) + 
     &                 AVGKERNEL(LLN,LN) * CH4_HAT_ADJ(LLN)
            ENDDO
         ENDDO

         ! Adjoint of x_m - x_a
         DO LN = lind, LLLEO 
           ! fwd code:
           !GC_CH4(LN)   = MAX(GC_CH4(LN), 1d-10)
           !CH4_PERT(LN) = LOG(GC_CH4(LN)) - LOG(PRIOR(LN))
           ! adj code:
           IF ( GC_CH4_onLEO(LN) > 1d-10 ) THEN 
              GC_CH4_onLEO_ADJ(LN) = 1d0 / GC_CH4_onLEO(LN) * 
     &                                         CH4_PERT_ADJ(LN)
           ELSE 
              GC_CH4_onLEO_ADJ(LN) = 1d0 / 1d-10 * CH4_PERT_ADJ(LN)
           ENDIF
         ENDDO


         ! Adjoint of interpolation
         DO LN=lind,LLLEO
         DO LG=1,LLPAR
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) + 
     &           GRIDMAP(LG,LN) * GC_CH4_onLEO_ADJ(LN)
         ENDDO
         ENDDO


         ! Adjoint of unit conversion
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
         COST_FUNC = COST_FUNC + SUM(NEW_COST(:))

      ENDDO   ! End looping over each observation in this grid box
      ENDDO   ! End looping over each grid box JJ
      ENDDO   ! End looping over each grid box II

!!$OMP END PARALLEL DO


! -----------------------------------------------------------------------
!    Use this section to test the adjoint of the LEO_CH4 operator by
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
      CALL CALC_LEO_CH4_FORCE_FD( COST_FUNC_0, PERT, ADJ )
      ADJ_SAVE(:) = ADJ(:)

      DO LN=lind,LLLEO
         DOFS = DOFS + AVGKERNEL(LN,LN)
      ENDDO

      ! Write identifying information to top of satellite diagnostic file
      WRITE(116,212) 'COST_FUNC_0:',( COST_FUNC_0 )
      WRITE(116,212) 'RANDOM ERROR',RANDNUM(NT)
      WRITE(116,212) 'DOFS        ',DOFS
      !WRITE(116,*)   (AVGKERNEL(1,LN),LN=1,13)
      !WRITE(116,*)   (OBSERROR(1,LN),LN=1,13)


      ! Perform finite difference testing at each vertical level
      DO LG = 1, 47

         ! Positive perturbation to GEOS-Chem CH4 columns
         PERT(:) = 0.0
         PERT(LG) = 0.001
         COST_FUNC_pos = 0D0
         CALL CALC_LEO_CH4_FORCE_FD( COST_FUNC_pos, PERT, ADJ )

         ! Negative perturbation to GEOS-Chem CH4 columns
         PERT(:) = 0.0
         PERT(LG) = -0.001
         COST_FUNC_neg = 0D0
         CALL CALC_LEO_CH4_FORCE_FD( COST_FUNC_neg, PERT, ADJ )

         ! Calculate dJ/dCH4 from perturbations
         FD_CEN(LG)   = ( COST_FUNC_pos - COST_FUNC_neg ) / 0.2d0
         FD_POS(LG)   = ( COST_FUNC_pos - COST_FUNC_0 )   / 0.1d0
         FD_NEG(LG)   = ( COST_FUNC_0 - COST_FUNC_neg )   / 0.1d0

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
      print*, ' LEO contribution this hour = ', COST_FUNC - OLD_COST 
      print*, ' # Obs analyzed this hour       = ', NTh
      print*, ' # Obs analyzed today           = ', NT
      print*, ' # Obs analyzed total           = ', NTT



      ! Return to calling program
      END SUBROUTINE CALC_LEO_CH4_FORCE

!------------------------------------------------------------------------------



      SUBROUTINE CALC_LEO_CH4_FORCE_FD( COST_FUNC_A, PERT, ADJ )
!
!******************************************************************************
!  Subroutine CALC_LEO_CH4_FORCE calculates the adjoint forcing from the LEO
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
      INTEGER                     :: NT  
      INTEGER                     :: LG, LN, LLN, II, JJ, NB, JMIN, OB
      INTEGER                     :: nlev, lind, IU_IN
      INTEGER                     :: nboxes, nobs
      INTEGER                     :: NTSTART, NTSTOP 
      TYPE (XPLEX)               :: GC_CH4_TRUE_ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: CH4_PRIOR(IIPAR,JJPAR,LLLEO)
      TYPE (XPLEX)                     :: DUMMY_PRIOR(IIPAR,JJPAR,LLLEO)
      TYPE (XPLEX)                      :: DUMMY_TRUE(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                      :: GC_PCENTER(LLPAR)
      TYPE (XPLEX)                      :: GC_PEDGE(LLPAR)
      TYPE (XPLEX)                      :: GC_AD(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_OB(LLPAR)
      TYPE (XPLEX)                      :: thispcen(LLPAR)
      TYPE (XPLEX)                      :: thispedg(LLPAR)
      TYPE (XPLEX)                      :: thisad(LLPAR)
      TYPE (XPLEX)                      :: thisch4(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_onLEO(LLLEO)
      TYPE (XPLEX)                      :: GC_CH4_onLEO_OB(LLLEO)
      TYPE (XPLEX)                      :: GRIDMAP(LLPAR,LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_OB(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_ADJ(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_werr(LLLEO)
      TYPE (XPLEX)                      :: CH4_HAT_werr_ADJ(LLLEO)
      TYPE (XPLEX)                      :: CH4_PERT(LLLEO)
      TYPE (XPLEX)                      :: CH4_PERT_OB(LLLEO)
      TYPE (XPLEX)                      :: CH4_PERT_ADJ(LLLEO)
      TYPE (XPLEX)                      :: frac, frac_total
      TYPE (XPLEX)                      :: latmin, Jfrac_min, Jfrac
      TYPE (XPLEX)                      :: box_area, cloud_frac
      TYPE (XPLEX)                      :: mass_air, mole_air, mole_ch4
      TYPE (XPLEX)                      :: LHS, RHS, GC_XCH4, XTAU
      TYPE (XPLEX)                      :: DIFF(LLLEO)
      TYPE (XPLEX)                      :: FORCE(LLLEO)
      TYPE (XPLEX)                      :: DIFF_ADJ(LLLEO)
      TYPE (XPLEX)                      :: thisforce(LLPAR)
      TYPE (XPLEX)                      :: GC_CH4_onLEO_ADJ(LLLEO)
      TYPE (XPLEX)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: NEW_COST(MAXLEO)
      TYPE (XPLEX)                      :: OLD_COST
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      LOGICAL, SAVE               :: DO_FDTEST = .TRUE.
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME
      CHARACTER(LEN=255)          :: FILENAME_OBS



      !=================================================================
      ! CALC_LEO_CH4_FORCE_FD begins here!
      !=================================================================

      print*, '     - CALC_LEO_CH4_FORCE_FD '

      NEW_COST(:) = 0d0 


      ! Read "TRUE" state for this time step
      GC_CH4_TRUE_ARRAY(:,:,:) = 0d0
!      FILENAME_OBS = '/home/kjw/GEOS-Chem/gcadj_std/runs/v8-02-01/' //
!     &               'ch4/leo/' // GET_RES_EXT() // '/adjtmp/'  //
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

      ! Read a priori vertical profiles from file
      FILENAME = '/home/kjw/new_satellites/leo/data/' //  
     &             'leo_prior.' //  GET_RES_EXT() // '.bpch'
      XTAU = GET_TAU0( 1, 1, 1985 )
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,    
     &                 XTAU,      IGLOB,     JGLOB, 
     &                 LLLEO, DUMMY_PRIOR, QUIET=.TRUE.  )
      CH4_PRIOR(:,:,:) = DUMMY_PRIOR(:,:,:)



      ! Select arbitrary II, JJ and NT value
      II=40
      JJ=JJPAR-10
      NT=100

      ! Initialize variables
      GC_PCENTER(:)    = 0d0
      GC_PEDGE(:)      = 0d0
      GC_AD(:)         = 0d0
      GC_CH4_NATIVE(:) = 0d0
      GC_CH4_onLEO(:)    = 0d0
      GC_CH4_onLEO_OB(:) = 0d0
      CH4_HAT_werr(:)  = 0d0
      DIFF(:)          = 0d0
      FORCE(:)         = 0d0


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
      !     GEOS-Chem surface pressure < LEO pressure levels
      nlev = count( (PRESSURE_EDGE%r) .LT. (GC_PEDGE(1)%r) )
      IF ( nlev .LT. 13 )  nlev = nlev + 1
      lind = LLLEO + 1 - nlev ! minimum vertical index on LEO grid

      ! Get interpolation matrix that maps GEOS-Chem to LEO grid
      GRIDMAP(1:LLPAR, 1:LLLEO) = 
     &     GET_INTMAP( GC_PEDGE, PRESSURE_EDGE, nlev )

      ! Get GEOS-Chem column from "truth" run to make pseudo-observations
      GC_CH4_NATIVE_OB(:) = 0d0
      GC_CH4_NATIVE_OB(:) = GC_CH4_TRUE_ARRAY(II,JJ,:)


      ! Interpolate GEOS-Chem CH4 column and observation to LEO grid
      ! Column in [v/v]
      DO LN = lind, LLLEO
         GC_CH4_onLEO(LN) = 0d0 
         GC_CH4_onLEO_OB(LN) = 0d0 
         DO LG = 1, LLPAR
            GC_CH4_onLEO(LN)    = GC_CH4_onLEO(LN) 
     &           + GRIDMAP(LG,LN) * GC_CH4_NATIVE(LG) 
            GC_CH4_onLEO_OB(LN) = GC_CH4_onLEO_OB(LN) 
     &           + GRIDMAP(LG,LN) * GC_CH4_NATIVE_OB(LG) 
         ENDDO
      ENDDO



         !--------------------------------------------------------------
         ! Apply LEO observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by LEO [molec/cm2]
         !    x_a   = LEO apriori column               [molec/cm2]
         !    x_m   = GC modeled column on LEO grid    [molec/cm2]
         !    A     = LEO averaging kernel 
         !--------------------------------------------------------------

         ! x_m - x_a for model and "observation"
         !    [v/v] --> ln( v/v ) happens here
         DO LN = lind, LLLEO
           GC_CH4_onLEO(LN)   =MAX(GC_CH4_onLEO(LN),   1d-10)
           GC_CH4_onLEO_OB(LN)=MAX(GC_CH4_onLEO_OB(LN),1d-10)
           CH4_PERT(LN)           =LOG( GC_CH4_onLEO(LN) ) - 
     &                             LOG( CH4_PRIOR(II,JJ,LN) )
           CH4_PERT_OB(LN)        =LOG( GC_CH4_onLEO_OB(LN) ) - 
     &                             LOG( CH4_PRIOR(II,JJ,LN) )
         ENDDO

         ! x_a + A_k * ( x_m - x_a ) for model and "observation"
         CH4_HAT(:)=CH4_PERT(:)
         DO LN = lind, LLLEO
            CH4_HAT(LN)    = 0d0
            CH4_HAT_OB(LN) = 0d0
         
            DO LLN = lind, LLLEO
               CH4_HAT(LN) = CH4_HAT(LN) 
     &                    + AVGKERNEL(LN,LLN) * CH4_PERT(LLN)
               CH4_HAT_OB(LN) = CH4_HAT_OB(LN) 
     &                    + AVGKERNEL(LN,LLN) * CH4_PERT_OB(LLN)
            ENDDO
            CH4_HAT(LN)    = CH4_HAT(LN)    + LOG( CH4_PRIOR(II,JJ,LN) )
            CH4_HAT_OB(LN) = CH4_HAT_OB(LN) + LOG( CH4_PRIOR(II,JJ,LN) )
         
         ENDDO


         ! For safety, initialize these up to LLLEO

         ! Add random error to this observation
         CH4_HAT_werr(:) = CH4_HAT(:)
         DO LN  = lind, LLLEO

            CH4_HAT_werr(LN) = CH4_HAT(LN)
            DO LLN = lind, LLLEO
               CH4_HAT_werr(LN) = CH4_HAT_werr(LN) + 
     &              CH4_HAT(LN) * RANDNUM(NT) * OBSERROR(LN,LLN)
            ENDDO
         ENDDO


         !-------------------------------------------------------------
         ! Calculate cost function, given S is observation error covariance matrix
         !     Sobs = 1x1 array [ (molec/cm2) ^2 ]
         ! J = [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------

         ! Calculate difference between modeled and observed profile
         DO LN = lind, LLLEO
            DIFF(LN) = CH4_HAT_werr(LN) - CH4_HAT_OB(LN)
         ENDDO 

         ! Calculate adjoint forcing: 2 * DIFF^T * S_{obs}^{-1}
         !         and cost function: DIFF^T * S_{obs}^{-1} * DIFF
         DO LN = lind, LLLEO
            DO LLN = lind, LLLEO
               FORCE(LN)  = FORCE(LN) + 
     &              2d0 * OBSERROR_INV(LN,LLN) * DIFF(LLN)
            ENDDO
            NEW_COST(LN) = NEW_COST(LN) + 0.5*DIFF(LN)*FORCE(LN)
         ENDDO


         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------


         ! The adjoint forcing is 2 * S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ(:) = FORCE(:) 

         ! Adjoint of GEOS-Chem - Observation difference
         CH4_HAT_werr_ADJ(:) = DIFF_ADJ(:)

         ! Adjoint of adding random error to observation
         DO LN=lind,LLLEO
            CH4_HAT_ADJ(LN) = 0d0

            DO LLN=lind,LLLEO
               CH4_HAT_ADJ(LN) = CH4_HAT_ADJ(LN) + 
     &           CH4_HAT_werr_ADJ(LLN) * RANDNUM(NT) * OBSERROR(LLN,LN)
            ENDDO
         ENDDO
         CH4_HAT_ADJ(:) = CH4_HAT_werr_ADJ(:)

         ! Adjoint of LEO observation operator
         CH4_PERT_ADJ(:) = CH4_HAT_ADJ(:)
         DO LN=lind,LLLEO
            CH4_PERT_ADJ(LN) = 0D0

            DO LLN=lind,LLLEO
               CH4_PERT_ADJ(LN) = CH4_PERT_ADJ(LN) + 
     &                 AVGKERNEL(LLN,LN) * CH4_HAT_ADJ(LLN)
            ENDDO
         ENDDO

         ! Adjoint of x_m - x_a
         DO LN = lind, LLLEO 
           ! fwd code:
           !GC_CH4(LN)   = MAX(GC_CH4(LN), 1d-10)
           !CH4_PERT(LN) = LOG(GC_CH4(LN)) - LOG(PRIOR(LN))
           ! adj code:
           IF ( GC_CH4_onLEO(LN) > 1d-10 ) THEN 
              GC_CH4_onLEO_ADJ(LN) = 1d0 / GC_CH4_onLEO(LN) * 
     &                                         CH4_PERT_ADJ(LN)
           ELSE 
              GC_CH4_onLEO_ADJ(LN) = 1d0 / 1d-10 * CH4_PERT_ADJ(LN)
           ENDIF
         ENDDO


         ! Adjoint of interpolation
         DO LN=lind,LLLEO
         DO LG=1,LLPAR
            GC_CH4_NATIVE_ADJ(LG) = GC_CH4_NATIVE_ADJ(LG) + 
     &           GRIDMAP(LG,LN) * GC_CH4_onLEO_ADJ(LN)
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

         ! Update cost function 
         COST_FUNC_A = COST_FUNC_A + SUM(NEW_COST(:))


      ! Return to calling program
      END SUBROUTINE CALC_LEO_CH4_FORCE_FD


!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( GC_PEDGE, LEO_PEDGE, nlev  )
     &         RESULT      ( M )
!
!******************************************************************************
!  Function GET_INTMAP creates the matrix that places GEOS-Chem column methane 
!     [molec/cm2] onto the 13-level pressure grid used by theoretical instrument, M.
!         GC[1x47] * M[47x13] = LEO[1x13]           (kjw, 7/21/11)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) GC_PEDGE   (TYPE (XPLEX)) : LLPAR bottom pressure edges of GEOS-Chem column
!  (2 ) SCIA_PEDGE (TYPE (XPLEX)) : LLLEO pressure edges of LEO column
!  (3 ) nlev       (TYPE (XPLEX)) : Number of LEO pressure levels to use
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) M          (TYPE (XPLEX)) : Interpolation matrix that maps GEOS-Chem to LEO grid
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
      TYPE (XPLEX)             :: LEO_PEDGE(LLLEO)
      INTEGER            :: nlev
 
      ! Return value 
      TYPE (XPLEX)             :: M(LLPAR,LLLEO)

      ! Local variables 
      INTEGER  :: LGC, LTM, LS, LG, LN, LIND
      TYPE (XPLEX)   :: DIFF, DELTA_SURFP
      TYPE (XPLEX)   :: GUP, GLO, NUP, NLO
      TYPE (XPLEX)   :: column_total(LLLEO)

      !=================================================================
      ! GET_INTMAP begins here!
      !=================================================================

      ! Initialize output
      M(:,:) = 0D0 

      ! Minimum LEO vertical level to use
      lind = LLLEO + 1 - nlev

      ! Loop over each pressure level of GEOS-Chem and LEO grids
      DO LG=1,LLPAR-1

         ! Get upper and lower pressure edges of GEOS-Chem box
         GUP = GC_PEDGE( LG+1 )
         GLO = GC_PEDGE( LG   )
         
         DO LN=lind,LLLEO-1

            ! Get top and bottom pressures of LEO box
            NUP = LEO_PEDGE( LN+1 )
            NLO = LEO_PEDGE( LN   )

            ! If both GEOS-Chem edges are within the LEO box, map value = 1
            IF ( ( GUP .gt. NUP ) .AND. ( GLO .lt. NLO ) ) THEN
               M(LG,LN) = 1
            ENDIF

            ! If both GEOS-Chem stradles a LEO pressure level, interpolate
            IF ( ( GUP .lt. NUP ) .AND. ( GLO .gt. NUP ) ) THEN
               DIFF       = GLO - GUP
               M(LG,LN+1) = ( NUP - GUP ) / DIFF
               M(LG,LN  ) = ( GLO - NUP ) / DIFF
            ENDIF

         ENDDO
      ENDDO

      ! Add value for uppermost GEOS-Chem grid box
      M(LLPAR,LLLEO) = 1


      ! Correct for case in which GEOS-Chem pressure is higher than LEO
      IF ( GC_PEDGE(1) .GT. LEO_PEDGE(1) ) THEN


         ! If any part of GEOS-Chem box are under LEO_PEDGE(1), let 
         !   this GEOS-Chem grid box contribute to the observation because 
         !   LEO and GEOS-Chem should have same surface pressure. map value = 1
         DO LG=1,LLPAR-1

            ! If GEOS-Chem box entirely below LEO surface pressure
            IF ( ( GC_PEDGE(LG)   .GT. LEO_PEDGE(1) ) .AND.    
     &           ( GC_PEDGE(LG+1) .GT. LEO_PEDGE(1) ) ) THEN 
               M(LG,1) = 1
            ENDIF

            ! If GEOS-Chem box straddles LEO surface pressure
            IF ( ( GC_PEDGE(LG)   .GT. LEO_PEDGE(1) ) .AND.    
     &           ( GC_PEDGE(LG+1) .LT. LEO_PEDGE(1) ) ) THEN 
               DIFF = GC_PEDGE(LG) - GC_PEDGE( LG+1 )
               M(LG,1) = ( LEO_PEDGE(1) - GC_PEDGE(LG+1) ) / DIFF
            ENDIF
            
         ENDDO
      ENDIF


      ! Correct for case in which GEOS-Chem surface pressure is within 2nd LEO
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. LEO_PEDGE(2) ) THEN
         M(1,1) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 3rd LEO
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. LEO_PEDGE(3) ) THEN
         M(1,2) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 4th LEO
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. LEO_PEDGE(4) ) THEN
         M(1,3) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 5th LEO
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. LEO_PEDGE(5) ) THEN
         M(1,4) = 0.
      ENDIF
      ! Correct for case in which GEOS-Chem surface pressure is within 6th LEO
      ! pressure level.
      IF ( GC_PEDGE(1) .LT. LEO_PEDGE(6) ) THEN
         M(1,5) = 0.
      ENDIF

      ! Normalize each column of M to 1 so that we are not creating any molecules
      ! when mapping from GEOS-Chem to LEO grids.

      ! DO NOT do this since we are mapping molc/cm2, not 
      ! Initialize to be safe and calculate column total
      column_total(:) = 0d0
      column_total(:) = xplx(SUM( M%r, 1 ),SUM(M%i,1))

      ! Normalize columns to column_total
      DO LN=1,LLLEO
         IF ( column_total(LN) .EQ. 0. ) CYCLE
         M(:,LN) = M(:,LN) / column_total(LN)
      ENDDO

      ! Return to calling program
      END FUNCTION GET_INTMAP


!-----------------------------------------------------------------------------



      END MODULE LEO_CH4_MOD
