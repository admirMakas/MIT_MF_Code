!$Id: gosat_co2_mod.f,v 1.2 2011/02/23 00:08:48 daven Exp $
      MODULE GOSAT_CO2_MOD
  
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
 
      ! Parameters
      INTEGER, PARAMETER           :: MAXLEV = 20
      INTEGER, PARAMETER           :: MAXGOS = 2000


      ! Record to store data from each GOS obs
      TYPE GOS_CO2_OBS 
         INTEGER                      :: LGOS(1)
         TYPE (XPLEX)                       :: LAT(1)
         TYPE (XPLEX)                       :: LON(1)
         TYPE (XPLEX)                       :: TIME(1)
         TYPE (XPLEX)                       :: CO2(MAXLEV)
         TYPE (XPLEX)                       :: PRES(MAXLEV)
         TYPE (XPLEX)                       :: PRIOR(MAXLEV)
         TYPE (XPLEX)                       :: AVG_KERNEL(MAXLEV,MAXLEV)
         TYPE (XPLEX)                       :: S_OER(MAXLEV,MAXLEV)
         TYPE (XPLEX)                       :: S_OER_INV(MAXLEV,MAXLEV)
         INTEGER                      :: QF(1)
      ENDTYPE GOS_CO2_OBS  

      TYPE(GOS_CO2_OBS)                          :: GOS(MAXGOS)

      ! IDTCO2 isn't defined in tracerid_mod because people just assume 
      ! it is one. Define it here for now as a temporary patch. 
      INTEGER, PARAMETER :: IDTCO2   = 1 
      ! Same thing for TCVV(IDTCO2) 
      TYPE (XPLEX),  PARAMETER :: TCVV_CO2 = 28.97d0  / 44d0 

      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_GOS_CO2_OBS( YYYYMMDD, NGOS )
!
!******************************************************************************
!  Subroutine READ_GOS_CO2_OBS reads the file and passes back info contained
!  therein. (dkh, 10/12/10) 
! 
!  Based on READ_TES_NH3 OBS (dkh, 04/26/10) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD    INTEGER : Current year-month-day
!
!  Arguments as Output: 
!  ============================================================================
!  (1 ) NGOS      (INTEGER) : Number of GOS retrievals for current day 
!
!  Module variable as Output: 
!  ============================================================================
!  (1 ) GOS    (GOS_CO2_OBS) : CO2 retrieval for current day 
!     
!  NOTES:
!  (1 ) Add calculation of S_OER_INV, though eventually we probably want to 
!        do this offline. (dkh, 05/04/10) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE NETCDF 
      USE TIME_MOD,               ONLY : EXPAND_DATE



      ! Arguments
      INTEGER,            INTENT(IN)  :: YYYYMMDD
    
      ! local variables 
      INTEGER                         :: FID
      INTEGER                         :: LGOS
      INTEGER                         :: NGOS
      INTEGER                         :: START0(1), COUNT0(1)
      INTEGER                         :: START1(2), COUNT1(2)
      INTEGER                         :: START2(3), COUNT2(3)
      INTEGER                         :: N, J
      INTEGER                         :: NT_ID
      INTEGER                         :: CO2_ID
      INTEGER                         :: PS_ID
      INTEGER                         :: AK_ID
      INTEGER                         :: OE_ID
      INTEGER                         :: AP_ID
      INTEGER                         :: LA_ID
      INTEGER                         :: LO_ID
      INTEGER                         :: DY_ID
      INTEGER                         :: TM_ID
      INTEGER                         :: LV_ID
      INTEGER                         :: OI_ID
      INTEGER                         :: QF_ID
      CHARACTER(LEN=5)                :: TMP
      CHARACTER(LEN=255)              :: READ_FILENAME

      TYPE (XPLEX), PARAMETER               :: FILL = -999.0D0
      TYPE (XPLEX), PARAMETER               :: TOL  = 1d-04
      TYPE (XPLEX)                          :: U(MAXLEV,MAXLEV)
      TYPE (XPLEX)                          :: VT(MAXLEV,MAXLEV)
      TYPE (XPLEX)                          :: S(MAXLEV)
      TYPE (XPLEX)                          :: TMP1
      TYPE (XPLEX)                          :: TEST(MAXLEV,MAXLEV)
      INTEGER                         :: I, II, III

      !=================================================================
      ! READ_GOS_CO2_OBS begins here!
      !=================================================================

      ! filename root 
      READ_FILENAME = TRIM( 'acos-v27-YYYYMMDD.nc' )

      ! Expand date tokens in filename 
      CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 ) 

      ! Construct complete filename 
      READ_FILENAME = TRIM( DATA_DIR ) // TRIM( '../GOSAT_CO2/' ) // 
     &                TRIM( READ_FILENAME )


      WRITE(6,*) '    - READ_GOSAT_CO2_OBS: reading file: ', 
     &   READ_FILENAME

      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )
   
      !-------------------------------- 
      ! Get data record IDs
      !-------------------------------- 
      CALL CHECK( NF90_INQ_DIMID( FID, "nSamples",        NT_ID),  102 )
      CALL CHECK( NF90_INQ_VARID( FID, "species",        CO2_ID ), 102 )
      CALL CHECK( NF90_INQ_VARID( FID, "averagingKernel", AK_ID ), 104 )
      CALL CHECK( NF90_INQ_VARID( FID, "pressure",        PS_ID ), 105 )
!      CALL CHECK( NF90_INQ_VARID( FID, "observationerrorcovariance", 
!     &                                                    OE_ID ), 106 )  
      CALL CHECK( NF90_INQ_VARID( FID, "constraintvector",AP_ID ), 107 )
      CALL CHECK( NF90_INQ_VARID( FID, "latitude",        LA_ID ), 108 )
      CALL CHECK( NF90_INQ_VARID( FID, "longitude",       LO_ID ), 109 )
      CALL CHECK( NF90_INQ_VARID( FID, "date",            DY_ID ), 110 )
      CALL CHECK( NF90_INQ_VARID( FID, "time",            TM_ID ), 111 )
      CALL CHECK( NF90_INQ_VARID( FID, "invCovariance", 
     &                                                    OI_ID ), 112 )  
      CALL CHECK( NF90_INQ_VARID( FID, "lev",             LV_ID ), 113 ) 


      ! READ number of retrievals, NGOS
      CALL CHECK( NF90_INQUIRE_DIMENSION( FID, NT_ID, TMP, NGOS),  202 )

      ! READ quality flag 
      CALL CHECK( NF90_INQ_VARID( FID, "usability",       QF_ID ), 114 ) 

      print*, ' NGOS = ', NGOS

      !-------------------------------- 
      ! Read 0D Data
      !-------------------------------- 

      ! define record size 
      START0 = (/1/)
      COUNT0 = (/1/)

      ! loop over records
      DO N = 1, NGOS
 
         ! Update starting index
         START0(1) = N

         ! READ latitude
         CALL CHECK( NF90_GET_VAR  ( FID, LA_ID,
     &      GOS(N)%LAT,          START0, COUNT0 ), 301 )

         ! READ longitude 
         CALL CHECK( NF90_GET_VAR  ( FID, LO_ID,
     &      GOS(N)%LON,          START0, COUNT0 ), 302 )

         ! READ time
         CALL CHECK( NF90_GET_VAR  ( FID, TM_ID,
     &      GOS(N)%TIME,         START0, COUNT0 ), 303 )

         ! for GOSAT CO2, need to convert from HH.hour_frac to fraction of day
         GOS(N)%TIME(1) = GOS(N)%TIME(1) / 24d0

         ! READ levels
         CALL CHECK( NF90_GET_VAR  ( FID, LV_ID,
     &      GOS(N)%LGOS,         START0, COUNT0 ), 304 )

         ! READ quality flag  aka "usability"
         CALL CHECK( NF90_GET_VAR  ( FID, QF_ID,
     &      GOS(N)%QF,           START0, COUNT0 ), 305 )


      ENDDO

      ! debuggg
      print*, ' passed 0-D'

      print*, ' lev read in = ', GOS(1)%LGOS 

      !-------------------------------- 
      ! Find # of good levels for each
      !-------------------------------- 

      ! define record size 
      START1 = (/1,      1/)
      COUNT1 = (/MAXLEV, 1/)

      ! loop over records
      DO N = 1, NGOS
 
         ! Update starting index
         START1(2) = N
      
         ! READ CO2 column, CO2
         CALL CHECK( NF90_GET_VAR  ( FID, CO2_ID, 
     &      GOS(N)%CO2(1:MAXLEV), START1, COUNT1 ), 401 )

         ! For CO2, we now read the number of good levels in from the data file itself 
!         ! Now determine how many of the levels in CO2 are
!         ! 'good' and how many are just FILL. 
!         J = 1
!         DO WHILE ( J .le. MAXLEV )
! 
!            ! check if the value is good
!            IF ( GOS(N)%CO2(J) > FILL ) THEN
!               
!               ! save the number of good levels as LGOS 
!               GOS(N)%LGOS = MAXLEV - J + 1
!
!               ! and now we can exit the while loop
!               J = MAXLEV + 1
!
!            ! otherwise this level is just filler
!            ELSE 
! 
!               ! so proceed to the next one up
!               J = J + 1 
!
!            ENDIF
!
!         ENDDO

      ENDDO

      ! debuggg
      print*, ' passed 1-D a'
      print*, ' lev calc    = ', GOS(1)%LGOS 

      !-------------------------------- 
      ! Read 1D Data
      !-------------------------------- 

      ! loop over records
      DO N = 1, NGOS

         ! J is number of good levels
         J = GOS(N)%LGOS(1)

         ! define record size 
         ! for CO2, filler values are at the end, not the beginning
         !START1 = (/MAXLEV - J + 1, 1/)
         START1 = (/1, 1/)
         COUNT1 = (/J, 1/)

         ! Update starting index
         START1(2) = N
     
         ! READ CO2 column, CO2
         CALL CHECK( NF90_GET_VAR  ( FID, CO2_ID,
     &      GOS(N)%CO2(1:J),      START1, COUNT1 ), 401 )


         ! READ pressure levels, PRES
         CALL CHECK( NF90_GET_VAR  ( FID, PS_ID,
     &      GOS(N)%PRES(1:J),    START1, COUNT1 ), 402 )

         ! READ apriori CO2 column, PRIOR
         CALL CHECK( NF90_GET_VAR  ( FID, AP_ID, 
     &     GOS(N)%PRIOR(1:J),    START1, COUNT1 ), 403 )


      ENDDO
    
      ! debuggg
      print*, ' passed 1-D b'

    
      !-------------------------------- 
      ! Read 2D Data
      !-------------------------------- 
    
      ! loop over records
      DO N = 1, NGOS
 
         ! J is number of good levels
         J = GOS(N)%LGOS(1)

         ! define record size 
         ! for CO2, filler values are at the end, not the beginning
         !START2 = (/MAXLEV - J + 1, MAXLEV - J + 1, 1/)
         START2 = (/1, 1, 1/)
         COUNT2 = (/J, J, 1/)
 
         ! Update starting index
         START2(3) = N 

         ! READ averaging kernal, AVG_KERNEL
         CALL CHECK( NF90_GET_VAR  ( FID, AK_ID, 
     &      GOS(N)%AVG_KERNEL(1:J,1:J), START2, COUNT2), 501 )

!         ! READ observational error covariance 
!         CALL CHECK( NF90_GET_VAR  ( FID, OE_ID,
!     &      GOS(N)%S_OER(1:J,1:J),      START2, COUNT2), 502 )

         ! READ observational error covariance inverse
         CALL CHECK( NF90_GET_VAR  ( FID, OI_ID,
     &      GOS(N)%S_OER_INV(1:J,1:J),  START2, COUNT2), 503 )

      ENDDO

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )
  
      ! debuggg
      print*, ' passed 2-D '

      print*, ' CO2 prof = ', GOS(1)%CO2(:)
      print*, ' CO2 levels = ', GOS(1)%LGOS(1)
      print*, ' CO2 time   = ', GOS(1)%TIME(1)


!      !-------------------------------- 
!      ! Calculate S_OER_INV
!      !-------------------------------- 
!
!      ! loop over records
!      DO N = 1, NGOS
! 
!         J = GOS(N)%LGOS(1)
!
!         ! add regularization
!         DO II=1,J
!            GOS(N)%S_OER(II,II) = GOS(N)%S_OER(II,II)+ 0.001D0
!         ENDDO
!
!         CALL DGESVD_EXAMPLE
!         
!         CALL SVD( GOS(N)%S_OER(1:J,1:J), J, 
!     &                        U(1:J,1:J), S(1:J), 
!     &                       VT(1:J,1:J)          )
!
!      ! debuggg
!      print*, ' passed SVD 1'
!
!         ! U = S^-1 * U^T  
!         DO I = 1, J
!            DO II = 1, J
!               TEST(I,II) = U(II,I) / S(I)
!            ENDDO
!         ENDDO
!         U    = TEST
!         TEST = 0d0
!   
!      ! debuggg
!      print*, ' passed SVD 2'
!
!         ! S_OER_INV = V * S^-1 * U^T
!         DO I = 1, J
!            DO II = 1, J
!               TMP1 = 0d0
!               DO III = 1, J
!                  TMP1 = TMP1 + VT(III,I) * U(III,II)
!               ENDDO
!               GOS(N)%S_OER_INV(I,II) = TMP1
!            ENDDO
!         ENDDO
!
!      ! debuggg
!      print*, ' passed SVD 3'
!
!         ! TEST: calculate 2-norm of I - S_OER_INV * S_OER
!         DO I = 1, J
!            DO II = 1, J
!            TMP1 = 0d0
!               DO III = 1, J
!                  TMP1 = TMP1 
!     &                 + GOS(N)%S_OER_INV(III,I) * GOS(N)%S_OER(III,II)
!               ENDDO
!               TEST(I,II) = - TMP1
!            ENDDO
!            TEST(I,I) = ( TEST(I,I) + 1 ) ** 2
!         ENDDO
!
!      ! debuggg
!      print*, ' passed SVD 4'
!
!
!         IF ( SUM(TEST(1:J,1:J)) > TOL ) THEN  
!            print*, ' WARNING: inversion error for retv N = ', 
!     &                SUM(TEST(1:J,1:J)), N 
!            print*, '   in GOS obs ', READ_FILENAME 
!         ENDIF 
!
!      ENDDO  ! N
!
!      ! debuggg
!      print*, ' passed SVD f'


      ! Return to calling program
      END SUBROUTINE READ_GOS_CO2_OBS
!------------------------------------------------------------------------------

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
        CALL ERROR_STOP('netCDF error', 'tes_nh3_mod')
      ENDIF 

      ! Return to calling program
      END SUBROUTINE CHECK

!------------------------------------------------------------------------------

      SUBROUTINE CALC_GOS_CO2_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_GOS_CO2_FORCE calculates the adjoint forcing from the GOSAT
!  CO2 observations and updates the cost function. (dkh, 10/12/10) 
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (TYPE (XPLEX)) : Cost funciton                        [unitless]
!     
!     
!  NOTES:
!  
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE DAO_MOD,            ONLY : AD
      USE DAO_MOD,            ONLY : AIRDEN
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACER_MOD,         ONLY : XNUMOLAIR
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP


#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER                     :: NTSTART, NTSTOP, NT 
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL,     LGOS
      INTEGER                     :: JLOOP
      TYPE (XPLEX)                      :: GC_PRES(LLPAR)
      TYPE (XPLEX)                      :: GC_CO2_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_CO2(MAXLEV)
      TYPE (XPLEX)                      :: GC_PSURF
      TYPE (XPLEX)                      :: MAP(LLPAR,MAXLEV)
      TYPE (XPLEX)                      :: CO2_HAT(MAXLEV)
      TYPE (XPLEX)                      :: CO2_PERT(MAXLEV)
      TYPE (XPLEX)                      :: FORCE(MAXLEV)
      TYPE (XPLEX)                      :: DIFF(MAXLEV)
      TYPE (XPLEX)                      :: NEW_COST(MAXGOS)
      TYPE (XPLEX)                      :: OLD_COST
      TYPE (XPLEX), SAVE                :: TIME_FRAC(MAXGOS)
      INTEGER,SAVE                :: NGOS

      TYPE (XPLEX)                      :: GC_CO2_NATIVE_ADJ(LLPAR)
      TYPE (XPLEX)                      :: CO2_HAT_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: CO2_PERT_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: GC_CO2_ADJ(MAXLEV)
      TYPE (XPLEX)                      :: DIFF_ADJ(MAXLEV)
   
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME



      !=================================================================
      ! CALC_GOS_CO2_FORCE begins here!
      !=================================================================

      print*, '     - CALC_GOS_CO2_FORCE '
    
      ! Reset 
      NEW_COST = 0D0 

      ! Open files for diagnostic output
      IF ( FIRST ) THEN
         FILENAME = 'pres.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
     
         FILENAME = 'gc_co2.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 102,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
     
         FILENAME = 'tes_co2.NN.m'
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

         FILENAME = 'co2_pert_adj.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 108,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_co2_adj.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 109,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'exp_co2_hat.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 110,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_press.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 111,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_co2_native.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 112,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_co2_on_tes.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 113,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_co2_native_adj.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 114,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      ENDIF

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Check if it is the last hour of a day 
      IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN 
 
         ! Read the GOS CO2 file for this day 
         CALL READ_GOS_CO2_OBS( GET_NYMD(), NGOS )
 
         ! TIME is YYYYMMDD.frac-of-day.  Subtract date and save just time fraction
         ! Don't need to adjust this for GOSAT CO2, for which TIME is already 
         ! just the time fraction. 
         TIME_FRAC(1:NGOS) = GOS(1:NGOS)%TIME(1)

      ENDIF 

      ! Get the range of GOS retrievals for the current hour
      CALL GET_NT_RANGE( NGOS, GET_NHMS(), TIME_FRAC, NTSTART, NTSTOP ) 

      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 

         print*, ' No matching GOS CO2 obs for this hour'
        RETURN
      ENDIF 

      print*, ' for hour range: ', GET_NHMS(), TIME_FRAC(NTSTART), 
     &       TIME_FRAC(NTSTOP)
      print*, ' found record range: ', NTSTART, NTSTOP 

! need to update this in order to do i/o with this loop parallel 
!!      ! Now do a parallel loop for analyzing data 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, MAP, LGOS, IIJJ,  I, J,  L,   LL, JLOOP )
!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_CO2,  DIFF  )
!!$OMP+PRIVATE( GC_CO2_NATIVE, CO2_PERT,   CO2_HAT, FORCE )
!!$OMP+PRIVATE( GC_CO2_NATIVE_ADJ,        GC_CO2_ADJ     )
!!$OMP+PRIVATE( CO2_PERT_ADJ,             CO2_HAT_ADJ    )
!!$OMP+PRIVATE( DIFF_ADJ                                )
      DO NT  = NTSTART, NTSTOP, -1

         print*, '     - CALC_GOS_CO2_FORCE: analyzing record ', NT 

         ! quality screening
         IF ( GOS(NT)%QF(1) == 0 ) THEN 
            print*, ' BAD QF, skipping record ', NT
            CYCLE
         ENDIF

         ! skip antarctica
         IF ( GOS(NT)%LAT(1) < -60d0 ) THEN 
            print*, ' Skipp all data with latitude < 50 S ', NT
            CYCLE
         ENDIF

         ! For safety, initialize these up to LLGOS
         GC_CO2(:)       = 0d0 
         MAP(:,:)        = 0d0 
         CO2_HAT_ADJ(:)  = 0d0 
         FORCE(:)        = 0d0 


         ! Copy LGOS to make coding a bit cleaner
         LGOS = GOS(NT)%LGOS(1)

         ! Get grid box of current record
         IIJJ  = GET_IJ(DCMPLX(GOS(NT)%LON(1)),DCMPLX(GOS(NT)%LAT(1)))
         I     = IIJJ(1)
         J     = IIJJ(2)

         ! dkh debug  
         print*, 'I,J = ', I, J

         ! Get GC pressure levels (mbar) 
         DO L = 1, LLPAR
            GC_PRES(L) = GET_PCENTER(I,J,L)
         ENDDO

         ! Get GC surface pressure (mbar) 
         GC_PSURF = GET_PEDGE(I,J,1) 

         
         ! Calculate the interpolation weight matrix 
         MAP(1:LLPAR,1:LGOS) 
     &      = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
     &                    LGOS,  GOS(NT)%PRES(1:LGOS), GC_PSURF  )


         ! Get CO2 values at native model resolution
         GC_CO2_NATIVE(:) = CHK_STT(I,J,:,IDTCO2)
 
         ! Convert from kg/box to v/v
         GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2
     &                    / AD(I,J,:) 
 

         ! Interpolate GC CO2 column to TES grid 
         DO LL = 1, LGOS
            GC_CO2(LL) = 0d0 
            DO L = 1, LLPAR 
               GC_CO2(LL) = GC_CO2(LL) 
     &                    + MAP(L,LL) * GC_CO2_NATIVE(L) 
            ENDDO
         ENDDO

         ! dkh debug: compare profiles:
         print*, ' GC_PRES, GC_native_CO2 [ppm] '
         WRITE(6,100) (GC_PRES(L), GC_CO2_NATIVE(L)*1d6, 
     &                   L = LLPAR, 1, -1 )
         print*, ' GOS_PRES, GC_CO2  ' 
         WRITE(6,100) (GOS(NT)%PRES(LL), 
     &                 GC_CO2(LL)*1d6, LL = LGOS, 1, -1 ) 
 100  FORMAT(1X,F16.8,1X,F20.8)


         !--------------------------------------------------------------
         ! Apply GOS observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by TES [vmr]
         !    x_a   = GOS apriori column               [vmr]
         !    x_m   = GC modeled column                [vmr]
         !    A_k   = GOS averaging kernel 
         !--------------------------------------------------------------

         ! x_m - x_a
         DO L = 1, LGOS 
           CO2_PERT(L) = GC_CO2(L) - GOS(NT)%PRIOR(L)
         ENDDO
     
         ! x_a + A_k * ( x_m - x_a )  
         DO L = 1, LGOS
            CO2_HAT(L)    = 0d0 
            DO LL = 1, LGOS
               CO2_HAT(L) = CO2_HAT(L) 
     &                    + GOS(NT)%AVG_KERNEL(L,LL) * CO2_PERT(LL) 
            ENDDO
            CO2_HAT(L)    = CO2_HAT(L) + GOS(NT)%PRIOR(L)
         ENDDO


         !--------------------------------------------------------------
         ! Calculate cost function, given S is error in vmr
         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------

         ! Calculate difference between modeled and observed profile
         DO L = 1, LGOS 
            DIFF(L) = CO2_HAT(L) - GOS(NT)%CO2(L)
         ENDDO 
          
         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF 
         DO L = 1, LGOS
            FORCE(L)     = 0d0 
            DO LL = 1, LGOS
               FORCE(L)  = FORCE(L) + GOS(NT)%S_OER_INV(L,LL) * DIFF(LL)
            ENDDO
            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
         ENDDO

         ! dkh debug: compare profiles:
         print*, ' TES_PRIOR, CO2_HAT, CO2_GOS [ppm], diag(S^-1)'
         WRITE(6,101) ( 1d6 * GOS(NT)%PRIOR(L), 
     &                  1d6 * CO2_HAT(L),
     &                  1d6 * GOS(NT)%CO2(L),
     &                  GOS(NT)%S_OER_INV(L,L),
     &                  L,    L = LGOS, 1, -1   )

 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,d14.6,1x,i3)

         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

         ! dkh debug
         print*, 'DIFF , FORCE ' 
         WRITE(6,102) (DIFF(L), FORCE(L), 
     &       L = LGOS, 1, -1 )
 102  FORMAT(1X,d14.6,1X,d14.6)

         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ(:) = FORCE(:) 

         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !print*, ' FORCE with 1 for sensitivity ' 
         !ADJ_DIFF(:) = 1d0
         !NEW_COST(NT) = ?? SUM(ABS(LOG(O3_HAT(1:LGOS))))
         !print*, ' sumlog =', SUM(ABS(LOG(O3_HAT(:))))
         !print*, ' sumlog =', ABS(LOG(O3_HAT(:)))
 
         ! Adjoint of difference
         DO L = 1, LGOS 
            CO2_HAT_ADJ(L) =  DIFF_ADJ(L)
         ENDDO 

         ! adjoint of TES operator
         DO L  = 1, LGOS
            CO2_PERT_ADJ(L)    = 0d0
            DO LL = 1, LGOS
               CO2_PERT_ADJ(L) = CO2_PERT_ADJ(L) 
     &                         + GOS(NT)%AVG_KERNEL(LL,L) 
     &                         * CO2_HAT_ADJ(LL)
           ENDDO
         ENDDO
 
         ! Adjoint of x_m - x_a
         DO L = 1, LGOS 
           ! fwd code:
           !CO2_PERT(L) = GC_CO2(L) - GOS(NT)%PRIOR(L)
           ! adj code:
           GC_CO2_ADJ(L) = CO2_PERT_ADJ(L)
         ENDDO

         ! dkh debug
         print*, 'CO2_HAT_ADJ, CO2_PERT_ADJ, GC_CO2_ADJ'
         WRITE(6,103) (CO2_HAT_ADJ(L), CO2_PERT_ADJ(L), GC_CO2_ADJ(L), 
     &       L = LGOS, 1, -1 )
 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)

         ! adjoint of interpolation 
         DO L  = 1, LLPAR
            GC_CO2_NATIVE_ADJ(L) = 0d0 
            DO LL = 1, LGOS
               GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)
     &                              + MAP(L,LL) * GC_CO2_ADJ(LL)
            ENDDO
         ENDDO
         
         WRITE(114,112) ( GC_CO2_NATIVE_ADJ(L),     L=LLPAR,1,-1)

         ! Adjoint of unit conversion 
         GC_CO2_NATIVE_ADJ(:) = GC_CO2_NATIVE_ADJ(:) * TCVV_CO2
     &                        / AD(I,J,:) 

           
        ! Pass adjoint back to adjoint tracer array
         STT_ADJ(I,J,:,IDTCO2)  = STT_ADJ(I,J,:,IDTCO2)
     &                          + GC_CO2_NATIVE_ADJ(:)

         ! dkh debug
         print*, 'GC_CO2_NATIVE_ADJ conv '
         WRITE(6,104) (GC_CO2_NATIVE_ADJ(L), L = LLPAR, 1, -1 )
 104  FORMAT(1X,d14.6)


         WRITE(101,110) ( GOS(NT)%PRES(LL),        LL=LGOS,1,-1)
         WRITE(102,110) ( 1d6 * GC_CO2(LL),        LL=LGOS,1,-1)
         WRITE(103,110) ( 1d6 * GOS(NT)%CO2(LL),   LL=LGOS,1,-1)
         WRITE(104,110) ( 1d6 * GOS(NT)%PRIOR(LL), LL=LGOS,1,-1)
         WRITE(105,110) ( DIFF(LL),                LL=LGOS,1,-1)
         WRITE(106,112) ( FORCE(LL),               LL=LGOS,1,-1)
         WRITE(107,111) NT, LGOS
         WRITE(108,112) ( CO2_PERT_ADJ(LL),        LL=LGOS,1,-1)
         WRITE(109,112) ( GC_CO2_ADJ(LL),          LL=LGOS,1,-1)
         WRITE(110,110) ( 1d6 * CO2_HAT(LL),       LL=LGOS,1,-1)
         WRITE(111,110) ( GC_PRES(L),              L=LLPAR,1,-1)
         WRITE(112,110) ( 1d6 * GC_CO2_NATIVE(L),  L=LLPAR,1,-1)
         WRITE(113,110) ( 1d6 * GC_CO2(LL),        LL=LGOS,1,-1)
 110     FORMAT(F18.6,1X)
 111     FORMAT(i4,1X,i4,1x)
 112     FORMAT(D14.6,1X)


      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function 
      COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))

      IF ( FIRST ) FIRST = .FALSE. 

      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, ' GOS contribution           = ', COST_FUNC - OLD_COST  

      ! Return to calling program
      END SUBROUTINE CALC_GOS_CO2_FORCE

!!------------------------------------------------------------------------------
!
!      SUBROUTINE CALC_TES_O3_FORCE_FD( COST_FUNC, PERT, ADJ )
!!
!!******************************************************************************
!!  Subroutine CALC_TES_O3_FORCE_FD tests the adjoint of CALC_TES_O3_FORCE
!!  (dkh, 05/05/10) 
!!    
!!    Can be driven with:
!!      PERT(:) = 1D0
!!      CALL CALC_TES_O3_FORCE_FD( COST_FUNC_0, PERT, ADJ )
!!      ADJ_SAVE(:) = ADJ(:)
!!      print*, 'do3:  COST_FUNC_0 = ', COST_FUNC_0
!!      DO L = 1, 30
!!         PERT(:) = 1D0
!!         PERT(L) = 1.1
!!         COST_FUNC = 0D0
!!         CALL CALC_TES_O3_FORCE_FD( COST_FUNC_1, PERT, ADJ )
!!         PERT(L) = 0.9
!!         COST_FUNC = 0D0
!!         CALL CALC_TES_O3_FORCE_FD( COST_FUNC_2, PERT, ADJ )
!!         FD(L)       = ( COST_FUNC_1 - COST_FUNC_2 ) / 0.2d0
!!         print*, 'do3:  FD  = ', FD(L), L
!!         print*, 'do3:  ADJ = ', ADJ_SAVE(L), L
!!         print*, 'do3:  COST = ', COST_FUNC, L
!!         print*, 'do3:  FD / ADJ ', FD(L) / ADJ_SAVE(L) , L
!!      ENDDO
!!
!! 
!! 
!!
!!  Arguments as Input/Output:
!!  ============================================================================
!!  (1 ) COST_FUNC (TYPE (XPLEX)) : Cost funciton                        [unitless]
!!     
!!     
!!  NOTES:
!!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!!  (1 ) Add more diagnostics.  Now read and write doubled O3 (dkh, 11/08/09) 
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
!      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
!      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
!      USE ADJ_ARRAYS_MOD,     ONLY : O3_PROF_SAV
!      USE CHECKPT_MOD,        ONLY : CHK_STT
!      USE COMODE_MOD,         ONLY : CSPEC, JLOP
!      USE COMODE_MOD,         ONLY : CSPEC_ADJ_FORCE
!      USE DAO_MOD,            ONLY : AD
!      USE DAO_MOD,            ONLY : AIRDEN
!      USE DAO_MOD,            ONLY : BXHEIGHT
!      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
!      USE GRID_MOD,           ONLY : GET_IJ
!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
!      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
!      USE TIME_MOD,           ONLY : GET_TS_CHEM
!      USE TRACER_MOD,         ONLY : XNUMOLAIR
!      USE TRACERID_MOD,       ONLY : IDO3
!      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
!
!
!#     include      "CMN_SIZE"      ! Size params
!
!      ! Arguments
!      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
!  
!      TYPE (XPLEX), INTENT(IN)          :: PERT(LLPAR) 
!      TYPE (XPLEX), INTENT(OUT)         :: ADJ(LLPAR)
!
!      ! Local variables 
!      INTEGER                     :: NTSTART, NTSTOP, NT 
!      INTEGER                     :: IIJJ(2), I,      J
!      INTEGER                     :: L,       LL,     LGOS
!      INTEGER                     :: JLOOP
!      TYPE (XPLEX)                      :: GC_PRES(LLPAR)
!      TYPE (XPLEX)                      :: GC_O3_NATIVE(LLPAR)
!      TYPE (XPLEX)                      :: GC_O3(MAXLEV)
!      TYPE (XPLEX)                      :: GC_PSURF
!      TYPE (XPLEX)                      :: MAP(LLPAR,MAXLEV)
!      TYPE (XPLEX)                      :: O3_HAT(MAXLEV)
!      TYPE (XPLEX)                      :: O3_PERT(MAXLEV)
!      TYPE (XPLEX)                      :: FORCE(MAXLEV)
!      TYPE (XPLEX)                      :: DIFF(MAXLEV)
!      TYPE (XPLEX)                      :: NEW_COST(MAXTES)
!      TYPE (XPLEX)                      :: OLD_COST
!      TYPE (XPLEX), SAVE                :: TIME_FRAC(MAXTES)
!      INTEGER,SAVE                :: NTES 
!
!      TYPE (XPLEX)                      :: GC_O3_NATIVE_ADJ(LLPAR)
!      TYPE (XPLEX)                      :: CO2_HAT_ADJ(MAXLEV)
!      TYPE (XPLEX)                      :: O3_PERT_ADJ(MAXLEV)
!      TYPE (XPLEX)                      :: GC_O3_ADJ(MAXLEV)
!      TYPE (XPLEX)                      :: DIFF_ADJ(MAXLEV)
!   
!      LOGICAL, SAVE               :: FIRST = .TRUE. 
!      INTEGER                     :: IOS
!      CHARACTER(LEN=255)          :: FILENAME
!
!
!
!      !=================================================================
!      ! CALC_TES_O3_FORCE_FD begins here!
!      !=================================================================
!
!      print*, '     - CALC_TES_O3_FORCE_FD '
!  
!      NEW_COST = 0D0
!
!      ! Open files for output
!      IF ( FIRST ) THEN
!         FILENAME = 'pres.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!     
!         FILENAME = 'gc_o3.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 102,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!     
!         FILENAME = 'tes_o3.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 103,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!     
!         FILENAME = 'apriori.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 104,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'diff.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 105,      FILE=TRIM( FILENAME    ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'force.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 106,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'nt_ll.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 107,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'o3_pert_adj.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 108,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'gc_o3_adj.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 109,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'exp_o3_hat.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 110,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'gc_press.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 111,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'gc_o3_native.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 112,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'gc_o3_on_tes.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 113,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'gc_o3_native_adj.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 114,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!
!      ENDIF
!
!      ! Save a value of the cost function first
!      OLD_COST = COST_FUNC
!
!      ! Check if it is the last hour of a day 
!!      IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN 
!      IF ( FIRST ) THEN 
! 
!         ! Read the TES O3 file for this day 
!         CALL READ_TES_O3_OBS( GET_NYMD(), NTES )
! 
!         ! TIME is YYYYMMDD.frac-of-day.  Subtract date and save just time fraction
!         TIME_FRAC(1:NTES) = TES(1:NTES)%TIME(1) - GET_NYMD()
!
!        FIRST = .FALSE. 
!      ENDIF 
!
!!      ! Get the range of TES retrievals for the current hour
!!      CALL GET_NT_RANGE( NTES, GET_NHMS(), TIME_FRAC, NTSTART, NTSTOP ) 
!!
!!      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 
!!
!!         print*, ' No matching TES O3 obs for this hour'
!!        RETURN
!!      ENDIF 
!!
!!      print*, ' for hour range: ', GET_NHMS(), TIME_FRAC(NTSTART), 
!!     &       TIME_FRAC(NTSTOP)
!!      print*, ' found record range: ', NTSTART, NTSTOP 
!
!       NTSTART = 1590
!       NTSTOP  = 1590
!
!! need to update this in order to do i/o with this loop parallel 
!!!      ! Now do a parallel loop for analyzing data 
!!!$OMP PARALLEL DO
!!!$OMP+DEFAULT( SHARED )
!!!$OMP+PRIVATE( NT, MAP, LGOS, IIJJ,  I, J,  L,   LL, JLOOP )
!!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_O3,  DIFF  )
!!!$OMP+PRIVATE( GC_O3_NATIVE, O3_PERT,   O3_HAT, FORCE )
!!!$OMP+PRIVATE( GC_O3_NATIVE_ADJ,        GC_O3_ADJ     )
!!!$OMP+PRIVATE( O3_PERT_ADJ,             CO2_HAT_ADJ    )
!!!$OMP+PRIVATE( DIFF_ADJ                                )
!      DO NT  = NTSTART, NTSTOP, -1 
!
!         print*, '     - CALC_TES_O3_FORCE: analyzing record ', NT 
!
!         ! For safety, initialize these up to LLGOS 
!         GC_O3(:)       = 0d0 
!         MAP(:,:)       = 0d0 
!         CO2_HAT_ADJ(:)  = 0d0 
!         FORCE(:)       = 0d0 
!
!
!         ! Copy LGOS to make coding a bit cleaner
!         LGOS = TES(NT)%LGOS(1)
!
!         ! Get grid box of current record
!         IIJJ  = GET_IJ( REAL(TES(NT)%LON(1),4), REAL(TES(NT)%LAT(1),4))
!         I     = IIJJ(1)
!         J     = IIJJ(2)
!
!         print*, 'I,J = ', I, J
!
!         ! Get GC pressure levels (mbar) 
!         DO L = 1, LLPAR
!            GC_PRES(L) = GET_PCENTER(I,J,L)
!         ENDDO
!
!         ! Get GC surface pressure (mbar) 
!         GC_PSURF = GET_PEDGE(I,J,1) 
!
!         
!         ! Calculate the interpolation weight matrix 
!         MAP(1:LLPAR,1:LGOS) 
!     &      = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
!     &                    LGOS,  TES(NT)%PRES(1:LGOS), GC_PSURF  )
!
!
!         ! Get O3 values at native model resolution
!         DO L = 1, LLPAR
! 
!          
!            ! check if in trop
!            IF ( ITS_IN_THE_TROP(I,J,L) ) THEN 
! 
!               JLOOP = JLOP(I,J,L) 
! 
!               ! get O3 from tropospheric array
!               IF ( JLOOP > 0 ) THEN 
!                   GC_O3_NATIVE(L) = CSPEC(JLOOP,IDO3) * PERT(L)
!    
!               ELSE 
!
!                  ! get O3 from climatology [#/cm2]
!                  GC_O3_NATIVE(L) = O3_PROF_SAV(I,J,L) /
!     &                            ( BXHEIGHT(I,J,L) * 100d0 ) 
!                  !GC_O3_NATIVE(L) = 1d0 
!               ENDIF 
!
!            ELSE 
!
!               ! get O3 from climatology [#/cm2]
!               GC_O3_NATIVE(L) = O3_PROF_SAV(I,J,L) /
!     &                         ( BXHEIGHT(I,J,L) * 100d0 ) 
!               !GC_O3_NATIVE(L) = 1d0 
!                 
!            ENDIF 
!
!            ! Convert from #/cm3 to v/v
!            GC_O3_NATIVE(L) = GC_O3_NATIVE(L) * 1d6 /
!     &                      ( AIRDEN(L,I,J)   * XNUMOLAIR )
! 
!         ENDDO
!
!
!         ! Interpolate GC O3 column to TES grid 
!         DO LL = 1, LGOS
!            GC_O3(LL) = 0d0 
!            DO L = 1, LLPAR 
!               GC_O3(LL) = GC_O3(LL) 
!     &                    + MAP(L,LL) * GC_O3_NATIVE(L) 
!            ENDDO
!         ENDDO
!
!         ! dkh debug: compare profiles:
!         print*, ' GC_PRES, GC_native_O3 [ppb] '
!         WRITE(6,100) (GC_PRES(L), GC_O3_NATIVE(L)*1d9, 
!     &                   L = LLPAR, 1, -1 )
!         print*, ' TES_PRES, GC_O3  ' 
!         WRITE(6,100) (TES(NT)%PRES(LL), 
!     &                 GC_O3(LL)*1d9, LL = LGOS, 1, -1 ) 
! 100  FORMAT(1X,F16.8,1X,F16.8)
!
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
!         DO L = 1, LGOS 
!           GC_O3(L)   = MAX(GC_O3(L), 1d-10)
!           O3_PERT(L) = LOG(GC_O3(L)) - LOG(TES(NT)%PRIOR(L))
!         ENDDO
!     
!         ! x_a + A_k * ( x_m - x_a )  
!         DO L = 1, LGOS
!            O3_HAT(L)    = 0d0 
!            DO LL = 1, LGOS
!               O3_HAT(L) = O3_HAT(L) 
!     &                    + TES(NT)%AVG_KERNEL(L,LL) * O3_PERT(LL) 
!            ENDDO
!            O3_HAT(L)    = O3_HAT(L) + LOG(TES(NT)%PRIOR(L))
!         ENDDO
!
!
!         !--------------------------------------------------------------
!         ! Calculate cost function, given S is error on ln(vmr)
!         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ ln(model - obs ]
!         !--------------------------------------------------------------
!
!         ! Calculate difference between modeled and observed profile
!         DO L = 1, LGOS 
!            IF ( TES(NT)%O3(L) > 0d0 ) THEN 
!               DIFF(L) = O3_HAT(L) - LOG( TES(NT)%O3(L) )
!            ELSE 
!               DIFF(L) = 0d0
!            ENDIF 
!         ENDDO 
!          
!         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF 
!         DO L = 1, LGOS
!            FORCE(L)     = 0d0 
!            DO LL = 1, LGOS
!               FORCE(L)  = FORCE(L) + TES(NT)%S_OER_INV(L,LL) * DIFF(LL)
!            ENDDO
!            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
!         ENDDO
!
!         ! dkh debug: compare profiles:
!         print*, ' TES_PRIOR, O3_HAT, O3_TES [ppb]'
!         WRITE(6,090) ( 1d9 * TES(NT)%PRIOR(L), 
!     &                  1d9 * EXP(O3_HAT(L)),
!     &                  1d9 * TES(NT)%O3(L),
!     &                  L,    L = LGOS, 1, -1   )
!
!         print*, ' TES_PRIOR, O3_HAT, O3_TES  [lnvmr], diag(S^-1)'
!         WRITE(6,101) ( LOG(TES(NT)%PRIOR(L)), O3_HAT(L), 
!     &                  LOG(TES(NT)%O3(L)), TES(NT)%S_OER_INV(L,L),
!     &                  L,  L = LGOS, 1, -1 )
! 090  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,i3)
! 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1X,d14.6,1x,i3)
!
!         !--------------------------------------------------------------
!         ! Begin adjoint calculations 
!         !--------------------------------------------------------------
!
!         ! dkh debug
!         print*, 'DIFF , FORCE ' 
!         WRITE(6,102) (DIFF(L), FORCE(L), 
!     &       L = LGOS, 1, -1 )
! 102  FORMAT(1X,d14.6,1X,d14.6)
!
!         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
!         DIFF_ADJ(:) = FORCE(:) 
!
!         !print*, ' FORCE with 1 for sensitivity ' 
!         !print*, ' FORCE with 1 for sensitivity ' 
!         !print*, ' FORCE with 1 for sensitivity ' 
!         !print*, ' FORCE with 1 for sensitivity ' 
!         !ADJ_DIFF(:) = 1d0
!         !NEW_COST(NT) = ?? SUM(ABS(LOG(O3_HAT(1:LGOS))))
!         !print*, ' sumlog =', SUM(ABS(LOG(O3_HAT(:))))
!         !print*, ' sumlog =', ABS(LOG(O3_HAT(:)))
! 
!         ! Adjoint of difference
!         DO L = 1, LGOS 
!            IF ( TES(NT)%O3(L) > 0d0 ) THEN 
!               CO2_HAT_ADJ(L) =  DIFF_ADJ(L)
!            ENDIF 
!         ENDDO 
!
!         ! adjoint of TES operator
!         DO L  = 1, LGOS
!            O3_PERT_ADJ(L)    = 0d0
!            DO LL = 1, LGOS
!               O3_PERT_ADJ(L) = O3_PERT_ADJ(L) 
!     &                        + TES(NT)%AVG_KERNEL(LL,L) 
!     &                        * CO2_HAT_ADJ(LL)
!           ENDDO
!         ENDDO
! 
!         ! Adjoint of x_m - x_a
!         DO L = 1, LGOS 
!           ! fwd code:
!           !GC_O3(L)   = MAX(GC_O3(L), 1d-10)
!           !O3_PERT(L) = LOG(GC_O3(L)) - LOG(TES(NT)%PRIOR(L))
!           ! adj code:
!           IF ( GC_O3(L) > 1d-10 ) THEN 
!              GC_O3_ADJ(L) = 1d0 / GC_O3(L) * O3_PERT_ADJ(L)
!           ELSE 
!              GC_O3_ADJ(L) = 1d0 / 1d-10 * O3_PERT_ADJ(L)
!           ENDIF 
!         ENDDO
!
!         ! dkh debug
!         print*, 'CO2_HAT_ADJ, O3_PERT_ADJ, GC_O3_ADJ'
!         WRITE(6,103) (CO2_HAT_ADJ(L), O3_PERT_ADJ(L), GC_O3_ADJ(L), 
!     &       L = LGOS, 1, -1 )
! 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)
!
!         ! adjoint of interpolation 
!         DO L  = 1, LLPAR
!            GC_O3_NATIVE_ADJ(L) = 0d0 
!            DO LL = 1, LGOS
!               GC_O3_NATIVE_ADJ(L) = GC_O3_NATIVE_ADJ(L)
!     &                             + MAP(L,LL) * GC_O3_ADJ(LL)
!            ENDDO
!         ENDDO
!         
!         WRITE(114,112) ( GC_O3_NATIVE_ADJ(L),     L=LLPAR,1,-1)
!
!         DO L = 1, LLPAR 
!
!            ! Adjoint of unit conversion 
!            GC_O3_NATIVE_ADJ(L) = GC_O3_NATIVE_ADJ(L) * 1d6 /
!     &              ( AIRDEN(L,I,J)     * XNUMOLAIR )
!
!
!            IF ( ITS_IN_THE_TROP(I,J,L) ) THEN 
! 
!               JLOOP = JLOP(I,J,L) 
!
!               IF ( JLOOP > 0 ) THEN  
! 
!                  ! Pass adjoint back to adjoint tracer array
!                  CSPEC_ADJ_FORCE(JLOOP,IDO3) =
!     &               CSPEC_ADJ_FORCE(JLOOP,IDO3) + GC_O3_NATIVE_ADJ(L)
!
!                  ADJ(L) = GC_O3_NATIVE_ADJ(L) * CSPEC(JLOOP,IDO3)
!
!               ENDIF 
!
!            ENDIF 
!
!         ENDDO
!
!         ! dkh debug
!         print*, 'GC_O3_NATIVE_ADJ conv '
!         WRITE(6,104) (GC_O3_NATIVE_ADJ(L), L = LLPAR, 1, -1 )
! 104  FORMAT(1X,d14.6)
!
!
!         WRITE(101,110) ( TES(NT)%PRES(LL),        LL=LGOS,1,-1)
!         WRITE(102,110) ( 1d9 * GC_O3(LL),         LL=LGOS,1,-1)
!         WRITE(103,110) ( 1d9 * TES(NT)%O3(LL),    LL=LGOS,1,-1)
!         WRITE(104,110) ( 1d9 * TES(NT)%PRIOR(LL), LL=LGOS,1,-1)
!         WRITE(105,110) ( DIFF(LL),                LL=LGOS,1,-1)
!         WRITE(106,112) ( FORCE(LL),               LL=LGOS,1,-1)
!         WRITE(107,111) NT, LGOS
!         WRITE(108,112) ( O3_PERT_ADJ(LL),         LL=LGOS,1,-1)
!         WRITE(109,112) ( GC_O3_ADJ(LL),           LL=LGOS,1,-1)
!         WRITE(110,110) ( 1d9 * EXP(O3_HAT(LL)),   LL=LGOS,1,-1)
!         WRITE(111,110) ( GC_PRES(L),              L=LLPAR,1,-1)
!         WRITE(112,110) ( 1d9 * GC_O3_NATIVE(L),   L=LLPAR,1,-1)
!         WRITE(113,110) ( 1d9 * GC_O3(LL),         LL=LGOS,1,-1)
! 110     FORMAT(F18.6,1X)
! 111     FORMAT(i4,1X,i4,1x)
! 112     FORMAT(D14.6,1X)
!
!
!      ENDDO  ! NT
!!!$OMP END PARALLEL DO
!
!      ! Update cost function 
!      COST_FUNC = SUM(NEW_COST(NTSTOP:NTSTART))
!
!      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
!      print*, ' TES contribution           = ', COST_FUNC - OLD_COST  
!
!      ! Return to calling program
!      END SUBROUTINE CALC_TES_O3_FORCE_FD
!
!!------------------------------------------------------------------------------

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
      IF ( HHMMSS == 230000 ) NTSAVE = NTES - 100
   
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
  


      print*, ' GET_NT_RANGE for ', HHMMSS
      print*, ' NTSAVE ', NTSAVE
      print*, ' NTES   ', NTES
   
      CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )


      ! Convert HH from hour to fraction of day 
      GC_HH_FRAC = DCMPLX(HH) / 24d0 
 
      ! one hour as a fraction of day 
      H1_FRAC    = 1d0 / 24d0 


      ! dkh debug
      print*, ' co2 time frac = ', TIME_FRAC


    
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
            IF (LGC == 0) HI = TM_SURFP

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

       ! Bug fix:  a more general version allows for multiples TES pressure
       ! levels to exist below the lowest GC pressure.  (dm, dkh, 09/30/10) 
       ! OLD code:
       !IF ( TM_PRESC(1) > GC_PRESC(1) ) THEN
       !   HINTERPZ(1,1)         = 1D0 
       !   HINTERPZ(2:LGC_TOP,1) = 0D0 
       !ENDIF
       ! New code:
       ! Loop over each pressure level of TM grid
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             HINTERPZ(1,LTM)         = 1D0
             HINTERPZ(2:LGC_TOP,LTM) = 0D0
          ENDIF
       ENDDO

      ! Return to calling program
      END FUNCTION GET_INTMAP

!!------------------------------------------------------------------------------
!      SUBROUTINE MAKE_O3_FILE(  )
!!
!!******************************************************************************
!!  Subroutine MAKE_O3_FILE saves O3 profiles that correspond to time and
!!  place of TES O3 obs. (dkh, 03/01/09) 
!!
!!  Module variables as Input:
!!  ============================================================================
!!  (1 ) O3_SAVE (TYPE (XPLEX)) : O3 profiles                             [ppmv]
!!     
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE BPCH2_MOD
!      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
!      USE ERROR_MOD,        ONLY : ERROR_STOP
!      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!
!#     include "CMN_SIZE"    ! Size params
!      
!      ! Local variables    
!      INTEGER              :: I, J, I0, J0, L, NT
!      CHARACTER(LEN=120)   :: FILENAME
!      TYPE (XPLEX)               :: DAT(1,LLPAR,MAXTES)
!      INTEGER, PARAMETER   :: IUN = 88 
!      
!      ! For binary punch file, version 2.0
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER, PARAMETER   :: HALFPOLAR = 1
!      INTEGER, PARAMETER   :: CENTER180 = 1
!
!      !=================================================================
!      ! MAKE_O3_FILE begins here!
!      !=================================================================
!      
!      FILENAME = TRIM( 'nh3.bpch' )
!      
!      ! Append data directory prefix
!      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )
!      
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'O3 profile '
!      CATEGORY = 'IJ-AVE-$'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!      UNIT     = 'ppmv'
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the checkpoint file for output -- binary punch format
!      !=================================================================
!
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_O3_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IUN, FILENAME, TITLE )
!
!      ! Temporarily store data in DAT as REAL4
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT ) 
!      DO NT = 1, MAXTES
!
!         DAT(1,:,NT) = REAL(O3_SAVE(:,NT))
!
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IUN,       MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  1,
!     &            UNIT,      1d0,       1d0,       RESERVED,
!     &            1,         LLPAR,     MAXTES,     I0+1,
!     &            J0+1,      1,         DAT )
!
!      ! Close file
!      CLOSE( IUN )        
!
!      print*, ' O3_SAVE sum write = ', SUM(O3_SAVE(:,:))
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_O3_FILE
!
!!------------------------------------------------------------------------------
!      SUBROUTINE READ_O3_FILE(  )
!!
!!******************************************************************************
!!  Subroutine READ_O3_FILE reads the GC modeled O3 profiles that correspond
!!  to the TES O3 times and locations. (dkh, 03/01/09) 
!!
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      USE BPCH2_MOD,         ONLY : READ_BPCH2
!      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
!
!
!#     include "CMN_SIZE"          ! Size parameters
!
!      ! Local variables
!      TYPE (XPLEX)                     :: DAT(1,LLPAR,MAXTES)
!      CHARACTER(LEN=255)         :: FILENAME
!
!      !=================================================================
!      ! READ_USA_MASK begins here!
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( ADJTMP_DIR )           //
!     &           'nh3.bpch'
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_O3_FILE: Reading ', a )
!      
!      
!      ! USA mask is stored in the bpch file as #2
!      CALL READ_BPCH2( FILENAME, 'IJ-AVE-$', 1,
!     &                 1d0,            1,        LLPAR, 
!     &                 MAXTES,    DAT,      QUIET=.TRUE. )
!      
!      ! Cast to TYPE (XPLEX)
!      O3_SAVE(:,:) = DAT(1,:,:)
!      
!      print*, ' O3_SAVE sum read = ', SUM(O3_SAVE(:,:))
!
!      ! Return to calling program
!      END SUBROUTINE READ_O3_FILE
!
!!-----------------------------------------------------------------------------
!      FUNCTION GET_DOUBLED_O3( NYMD, NHMS, LON, LAT ) RESULT( O3_DBL )
!!
!!******************************************************************************
!!  Subroutine GET_DOUBLED_O3 reads and returns the nh3 profiles from 
!!  model run with doubled emissions. (dkh, 11/08/09) 
!!
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      USE BPCH2_MOD,         ONLY : READ_BPCH2
!      USE DIRECTORY_MOD,     ONLY : DATA_DIR
!      USE TIME_MOD,          ONLY : EXPAND_DATE
!      USE TIME_MOD,          ONLY : GET_TAU
!
!
!#     include "CMN_SIZE"          ! Size parameters
!
!      ! Arguments   
!      INTEGER                    :: NYMD, NHMS
!      TYPE (XPLEX)                     :: LON,  LAT
!      
!      ! Function arg 
!      TYPE (XPLEX)                     :: O3_DBL(LLPAR)
!
!      ! Local variables
!      TYPE (XPLEX)                     :: DAT(144,91,20)
!      CHARACTER(LEN=255)         :: FILENAME
!      INTEGER                    :: IIJJ(2)
!
!      !=================================================================
!      ! GET_DOUBLED_O3 begins here!
!      !=================================================================
!
!      ! filename
!      FILENAME = 'nh3.YYYYMMDD.hhmm'
!
!      ! Expand filename
!      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )
!
!      ! Full path to file
!      FILENAME = TRIM( DATA_DIR )           //
!     &           'doubled_nh3/'             // 
!     &           TRIM( FILENAME )           //
!     &           TRIM( '00'     )           
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - GET_DOUBLED_O3: Reading ', a )
!      
!      ! dkh debug
!      print*, ' GET_TAU() = ', GET_TAU()
!      
!      ! Get data
!      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 29, 
!     &                 GET_TAU(), 144,      91, 
!     &                 20,        DAT,      QUIET=.FALSE. )
!      
!      IIJJ = GET_IJ_2x25( LON, LAT )
! 
!      print*, ' found doubled in I/J = ', IIJJ
!
!      ! just the column for the present location, and convert ppb to ppm
!      O3_DBL(1:20)     = REAL(DAT(IIJJ(1),IIJJ(2),:),8) / 1000d0 
!      O3_DBL(21:LLPAR) = 0d0 
!     
!      print*, ' O3_DBL = ', O3_DBL
! 
!      ! Return to calling program
!      END FUNCTION GET_DOUBLED_O3
!
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
      TYPE (XPLEX),  PARAMETER  :: DISIZE = 2.5d0
      TYPE (XPLEX),  PARAMETER  :: DJSIZE = 2.0d0
      INTEGER, PARAMETER  :: IIMAX  = 144
      INTEGER, PARAMETER  :: JJMAX  = 91
      
      
      !=================================================================
      ! GET_IJ_2x25 begins here!
      !=================================================================

      TLON = 180d0 + LON + DISIZE
      TLAT =  90d0 + LAT + DJSIZE
      
      I = TLON / DISIZE
      J = TLAT / DJSIZE

      
      IF ( TLON / DISIZE -DCMPLX(I)  >= 0.5d0 ) THEN
         I = I + 1
      ENDIF
      
      IF ( TLAT / DJSIZE -DCMPLX(J)  >= 0.5d0 ) THEN
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

!------------------------------------------------------------------------------

      SUBROUTINE SVD(A,N,U,S,VT)
!
!******************************************************************************
!  Subroutine SVD is a driver for the LAPACK SVD routine DGESVD. (dkh, 05/04/10) 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) A   (TYPE (XPLEX)) :  N x N matrix to decompose
!  (2 ) N  (INTEGER) :  N is dimension of A
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) U   (TYPE (XPLEX)) :  Array of left singular vectors
!  (2 ) S   (TYPE (XPLEX)) :  Vector of singular values
!  (3 ) VT  (TYPE (XPLEX)) :  Array of right singular vectors, TRANSPOSED 
!      
!     
!  NOTES:
!
*  Copyright (C) 2009-2010 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*  =============================================================================
*
*  DGESVD Example.
*  ==============
*
*  Program computes the singular value decomposition of a general
*  rectangular matrix A:
*
*    8.79   9.93   9.83   5.45   3.16
*    6.11   6.91   5.04  -0.27   7.98
*   -9.15  -7.93   4.86   4.85   3.01
*    9.57   1.64   8.83   0.74   5.80
*   -3.49   4.02   9.80  10.00   4.27
*    9.84   0.15  -8.99  -6.02  -5.31
*
*  Description.
*  ============
*
*  The routine computes the singular value decomposition (SVD) of a real
*  m-by-n matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written as
*
*  A = U*SIGMA*VT
*
*  where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
*  diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
*  is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and are
*  returned in descending order. The first min(m, n) columns of U and V are
*  the left and right singular vectors of A.
*
*  Note that the routine returns VT, not V.
*
*  Example Program Results.
*  ========================
*
* DGESVD Example Program Results
*
* Singular values
*  27.47  22.64   8.56   5.99   2.01
*
* Left singular vectors (stored columnwise)
*  -0.59   0.26   0.36   0.31   0.23
*  -0.40   0.24  -0.22  -0.75  -0.36
*  -0.03  -0.60  -0.45   0.23  -0.31
*  -0.43   0.24  -0.69   0.33   0.16
*  -0.47  -0.35   0.39   0.16  -0.52
*   0.29   0.58  -0.02   0.38  -0.65
*
* Right singular vectors (stored rowwise)
*  -0.25  -0.40  -0.69  -0.37  -0.41
*   0.81   0.36  -0.25  -0.37  -0.10
*  -0.26   0.70  -0.22   0.39  -0.49
*   0.40  -0.45   0.25   0.43  -0.62
*  -0.22   0.14   0.59  -0.63  -0.44
*  =============================================================================
!******************************************************************************
!
      ! Arguements 
      INTEGER,INTENT(IN)     :: N
      TYPE (XPLEX), INTENT(IN)     :: A(N,N)
      TYPE (XPLEX), INTENT(OUT)    :: U(N,N)
      TYPE (XPLEX), INTENT(OUT)    :: S(N)
      TYPE (XPLEX), INTENT(OUT)    :: VT(N,N)

      ! Local variables 
      INTEGER, PARAMETER     :: LWMAX = MAXLEV * 35 
      INTEGER                :: INFO, LWORK
      TYPE (XPLEX)       :: WORK( LWMAX )

*     .. External Subroutines ..
      EXTERNAL               :: DGESVD

*     .. Intrinsic Functions ..
      INTRINSIC              :: INT, MIN

      !=================================================================
      ! SVD begins here!
      !=================================================================

*     .. Executable Statements ..
      !WRITE(*,*)'DGESVD Example Program Results'
*
*     Query the optimal workspace.
*
      print*, ' here 1 '
      LWORK = -1
      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
     $             WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      print*, ' here 2 '
      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
     $             WORK, LWORK, INFO )
*
*     Check for convergence.
*
      print*, ' here 3 '
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF

!  Uncomment the following to print out singlular values, vectors (dkh, 05/04/10) 
!
!     Print singular values.
!
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
!
!     Print left singular vectors.
!
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
     $                   N, N, U, N   )
!
!     Print right singular vectors.
!
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
     $                   N, N, VT, N    )

      ! Return to calling program
      END SUBROUTINE SVD
!------------------------------------------------------------------------------
      SUBROUTINE DGESVD_EXAMPLE

*     .. Parameters ..
      INTEGER          M, N
      PARAMETER        ( M = 6, N = 5 )
      INTEGER          LDA, LDU, LDVT
      PARAMETER        ( LDA = M, LDU = M, LDVT = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      TYPE (XPLEX) A( LDA, N ), U( LDU, M ), VT( LDVT, N ), S( N ),
     $                 WORK( LWMAX )
      DATA             A/
     $  8.79, 6.11,-9.15, 9.57,-3.49, 9.84,
     $  9.93, 6.91,-7.93, 1.64, 4.02, 0.15,
     $  9.83, 5.04, 4.86, 8.83, 9.80,-8.99,
     $  5.45,-0.27, 4.85, 0.74,10.00,-6.02,
     $  3.16, 7.98, 3.01, 5.80, 4.27,-5.31
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DGESVD
      !EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DGESVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF
*
*     Print singular values.
*
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
*
*     Print left singular vectors.
*
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
     $                   M, N, U, LDU )
*
*     Print right singular vectors.
*
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
     $                   N, N, VT, LDVT )

*
*     End of DGESVD Example.
      END SUBROUTINE DGESVD_EXAMPLE
!------------------------------------------------------------------------------
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      TYPE (XPLEX) A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
! Change format of output (dkh, 05/04/10) 
! 9998 FORMAT( 11(:,1X,F6.2) )
 9998 FORMAT( 11(:,1X,E14.8) )
      RETURN

      END SUBROUTINE PRINT_MATRIX 
!------------------------------------------------------------------------------

      END MODULE GOSAT_CO2_MOD
