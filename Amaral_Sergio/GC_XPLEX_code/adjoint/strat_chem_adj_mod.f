!$ID$
!
!  Subroutine STRAT_CHEM_ADJ_MOD performs adjoint of strat chem. 
!
!  Based on forward model routine STRAT_CHEM_MOD.

! !INTERFACE:
!
      MODULE STRAT_CHEM_ADJ_MOD
!
! !USES:
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_STRAT_CHEM_ADJ
      !PUBLIC  :: CLEANUP_STRAT_CHEM
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      !PRIVATE :: INIT_STRAT_CHEM
      !PRIVATE :: GET_RATES
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!     
      ! Scalars            
      !TYPE (XPLEX)               :: DTCHEM
      
      ! Parameters
      !INTEGER, PARAMETER   :: NTR_GMI = 120 ! Number of species
                              ! 118 as output from GMI + NOx + Ox families
      
      !INTEGER, PARAMETER   :: MAX_FM  = 1 ! Max number of species in a fam
      ! Vestigial, as NOx and Ox families pre-processed, but may be useful
      ! for future uses, e.g., ClOx.
      
      ! Arrays 
      !TYPE (XPLEX),  ALLOCATABLE :: PROD(:,:,:,:)
      !TYPE (XPLEX),  ALLOCATABLE :: LOSS(:,:,:,:)
      !INTEGER, ALLOCATABLE :: GMI_TO_GC(:,:)
      !INTEGER, SAVE        :: ncID_strat_rates

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DO_STRAT_CHEM_ADJ
!
! !DESCRIPTION: Function DO\_STRAT\_CHEM is the driver routine for computing
!     the simple linearized stratospheric chemistry scheme for a host of species
!     whose prod/loss rates were determined from the GMI combo model. Ozone is
!     treated using either Linoz or Synoz.
!\\
!\\
! !INTERFACE:
!      
      SUBROUTINE DO_STRAT_CHEM_ADJ
!
! !USES:
!
      USE DAO_MOD,        ONLY : AD, CONVERT_UNITS
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE LOGICAL_MOD,    ONLY : LLINOZ, LPRT
      USE NETCDF_UTIL_MOD
      USE TIME_MOD,       ONLY : GET_MONTH, TIMESTAMP_STRING
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACER_MOD,     ONLY : N_TRACERS, STT, TCVV, TRACER_MW_KG
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL, ITS_IN_THE_STRAT
      ! adj_group (hml, 07/20/11)
      USE STRAT_CHEM_MOD, ONLY : PROD_0, LOSS_0
      USE STRAT_CHEM_MOD, ONLY : PROD, LOSS
      USE STRAT_CHEM_MOD, ONLY : DTCHEM
      USE STRAT_CHEM_MOD, ONLY : GET_RATES
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ 
      USE ADJ_ARRAYS_MOD, ONLY : PROD_SF,     LOSS_SF 
      USE ADJ_ARRAYS_MOD, ONLY : PROD_SF_ADJ, LOSS_SF_ADJ 
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD, ONLY : NSTPL
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE LINOZ_ADJ_MOD,  ONLY : DO_LINOZ_ADJ
      USE CHECKPOINT_MOD, ONLY : READ_BEFSTRAT_CHKFILE
      USE TIME_MOD,       ONLY : GET_NHMS
      USE TIME_MOD,       ONLY : GET_NYMD
      USE TRACER_MOD,     ONLY : STT_STRAT_TMP
      USE LOGICAL_ADJ_MOD,ONLY : LADJ_STRAT

#include "CMN_SIZE"
!
!EOP
!------------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE             :: LASTSEASON = -1
      INTEGER                   :: I, J, L, N, LMIN
      INTEGER                   :: IORD, JORD, KORD
      TYPE (XPLEX)                    :: t, P, k, M0
      TYPE (XPLEX)                    :: P_ADJ, k_ADJ, M0_ADJ
      TYPE (XPLEX)                    :: LOSS_ADJ, PROD_ADJ
      CHARACTER(LEN=16)         :: STAMP
      INTEGER                   :: NHMS
      INTEGER                   :: NYMD


      !===============================
      ! DO_STRAT_CHEM_ADJ begins here!
      !===============================

      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DO_STRAT_CHEM_ADJ: Strat chemistry at ', a )

      !================================================
      ! Determine the rates from disk; merge families
      !================================================

      ! Get the minimum level extent of the tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

      ! Use ITS_A_NEW_MONTH instead, which works for forward and adjoint
      !IF ( GET_MONTH() /= LASTMONTH ) THEN
      IF ( ITS_A_NEW_MONTH() ) THEN 

         WRITE(6,*) 'Getting new strat rates for month: ',GET_MONTH()

         IF ( LPRT ) CALL DEBUG_MSG( '### STRAT_CHEM_ADJ: at GET_RATES')

         ! Read rates for this month
         CALL GET_RATES( GET_MONTH() )

      ENDIF

      IF ( LPRT ) 
     &   CALL DEBUG_MSG( '### STRAT_CHEM_ADJ: at DO_STRAT_CHEM_ADJ' )

      ! READING STT FROM CHECKPOINT FILE (hml, 07/31/11) 
      NHMS     = GET_NHMS()
      NYMD     = GET_NYMD()
      CALL READ_BEFSTRAT_CHKFILE( NYMD, NHMS )

      WRITE(6,*) '-----------------------------------------------------'
      write(6,*) '    Doing strat chem ajdiont (STRAT_CHEM_ADJ_MOD)    '
      WRITE(6,*) '-----------------------------------------------------'

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, k, P, t, M0  )
!$OMP+PRIVATE( k_ADJ,   P_ADJ,   M0_ADJ )
!$OMP+PRIVATE( LOSS_ADJ,     PROD_ADJ   ) 
!$OMP+SCHEDULE( DYNAMIC )
      DO L = LMIN,LLPAR
      DO J = 1,JJPAR
      DO I = 1,IIPAR
      IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

         !when loop over only for NSTPL (hml)
         !DO N=1,N_TRACERS
         DO N = 1 , NSTPL

            ! Include something to expediate skipping past species
            ! that we do not have strat chem for. Prob put tracer on
            ! outermost loop.

            ! Now we will use GMI rate for Ox if LINOZ is off (hml, 10/31/11)  
!            ! Skip Ox -- if we're not using Linoz, we'll use Synoz below
            IF ( ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGOX_SIM() ) .and.
     &           LLINOZ .and. ( N .eq. IDTOx ) ) CYCLE

            ! Use NSTPL instead of N_TRACERS
            !! should we add this?
!            IF ( GMI_TO_GC( N, 1 ) .eq. 0 ) CYCLE

            IF ( LADJ_STRAT ) THEN

               PROD(I,J,L,N) = PROD_0(I,J,L,N) * PROD_SF(I,J,1,N)
               LOSS(I,J,L,N) = LOSS_0(I,J,L,N) * LOSS_SF(I,J,1,N)

            ENDIF

            !===============================
            ! Do chemical production and loss
            !===============================

            ! recalculate forward values to use for adjoint code (hml)
            t = DTCHEM                              ! timestep [s]
            k = LOSS(I,J,L,N)                       ! loss freq [s-1]
            P = PROD(I,J,L,N) * AD(I,J,L) / TCVV(N) ! production term [kg s-1]
            ! Use checkpointed value
            !M0= STT(I,J,L,N)                        ! initial mass [kg]
            M0= STT_STRAT_TMP(I,J,L,N)              ! initial mass [kg]

            ! debug test 
            IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
               print*, N,' STRAT TEST adj: k = ', k 
               print*, N,' STRAT TEST adj: P = ', P
               print*, N,' STRAT TEST adj: M0= ', M0
            ENDIF 

            ! No prod or loss at all
            IF ( k .eq. 0d0 .and. P .eq. 0d0 ) CYCLE

            ! Simple analytic solution to dM/dt = P - kM over [0,t]
            IF ( k .gt. 0d0 ) THEN
              ! fwd code:
              !STT(I,J,L,N) = M0 * exp(-k*t) + (P/k)*(1d0-exp(-k*t))
              ! adj code:
              M0_ADJ = STT_ADJ(I,J,L,N) * exp(-k*t)
              P_ADJ  = STT_ADJ(I,J,L,N) * (1d0 - exp(-k*t))/k
              k_ADJ  = STT_ADJ(I,J,L,N) 
     &               * ( -p/(k**2) + p/(k**2)*exp(-k*t)
     &               + (p*t/k)*exp(-k*t) - t * exp(-k*t) * M0 )
            ELSE
              ! fwd code:
              !STT(I,J,L,N) = M0 + P*t
              ! adj code:
              M0_ADJ = STT_ADJ(I,J,L,N)
              P_ADJ  = STT_ADJ(I,J,L,N) * t 
            ENDIF

            ! fwd code:
            !k = LOSS(I,J,L,N)                       ! loss freq [s-1]
            !P = PROD(I,J,L,N) * AD(I,J,L) / TCVV(N) ! production term [kg s-1]
            !M0 = STT(I,J,L,N)                       ! initial mass [kg]
            ! adj code:
            LOSS_ADJ          = K_ADJ
            PROD_ADJ          = P_ADJ * AD(I,J,L) / TCVV(N) 
            STT_ADJ (I,J,L,N) = M0_ADJ

            IF ( LADJ_STRAT ) THEN

              ! fwd code:
              !PROD(I,J,L,N) = PROD_0(I,J,L,N) * PROD_SF(I,J,1,N)
              !LOSS(I,J,L,N) = LOSS_0(I,J,L,N) * LOSS_SF(I,J,1,N)
              ! adj code:
              PROD_SF_ADJ(I,J,1,N) = PROD_SF_ADJ(I,J,1,N)
     &                             + PROD_0(I,J,L,N) * PROD_ADJ
              LOSS_SF_ADJ(I,J,1,N) = LOSS_SF_ADJ(I,J,1,N)
     &                             + LOSS_0(I,J,L,N) * LOSS_ADJ
            ENDIF

         ENDDO

       ENDIF
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO

      !======================
      ! Stratospheric Ozone !
      !======================
      ! Modified (hml, 10/31/11)
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Put ozone in v/v
         ! fwd code:
         !STT(:,:,:,IDTOX ) = STT(:,:,:,IDTOX) * TCVV( IDTOX ) / AD
         ! adj code: Put ozone back to kg
         STT_ADJ(:,:,:,IDTOX ) = 
     &                         STT_ADJ(:,:,:,IDTOX) * AD / TCVV( IDTOX )
!         STT_ADJ(:,:,:,IDTOX) = 
!     &                       STT_ADJ(:,:,:,IDTOX)* TCVV( IDTOX ) / AD

         IF ( LLINOZ ) THEN
            CALL DO_LINOZ_ADJ    ! Linoz
         ELSE
            ! must use Linoz or strat chem Ox fluxes for the adjoint 
         ENDIF

         ! Now move this into LINOZ (hml, 11/06/11)
         ! Put ozone back to kg
         ! fwd code:
         !STT(:,:,:,IDTOX) = STT(:,:,:,IDTOX) * TCVV( IDTOX ) / AD
         ! adj code: Put ozone in v/v
         STT_ADJ(:,:,:,IDTOX) = 
     &                       STT_ADJ(:,:,:,IDTOX)* TCVV( IDTOX ) / AD
!         STT_ADJ(:,:,:,IDTOX ) = 
!     &                         STT_ADJ(:,:,:,IDTOX) * AD / TCVV( IDTOX )

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         ! fwd code:
         !CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg
         ! adj code:
         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT_ADJ ) ! v/v -> kg

         ! adjoint LINOZ  does not support tagged Ox simulation for now (hml, 10/05/11)
         IF ( LLINOZ ) THEN
            CALL DO_LINOZ_ADJ       ! Linoz
         ELSE
            ! must use Linoz or strat chem Ox fluxes for the adjoint 
         ENDIF

         ! fwd code:
         !CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v
         ! adj code:
         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT_ADJ ) ! kg -> v/v

      ENDIF

      END SUBROUTINE DO_STRAT_CHEM_ADJ
!EOC
!!!!------------------------------------------------------------------------------
!!!!BOP
!!!!
!!!! !IROUTINE: GET_RATES
!!!!
!!!! !DESCRIPTION: Function GET\_RATES reads from disk the chemical production
!!!!  and loss rates for the species of interest
!!!!\\
!!!!\\
!!!! !INTERFACE:
!!!!
!!!      SUBROUTINE GET_RATES( THISMONTH )
!!!!
!!!! !USES:
!!!!
!!!      USE BPCH2_MOD,       ONLY : GET_NAME_EXT, GET_RES_EXT
!!!      USE DIRECTORY_MOD,   ONLY : DATA_DIR
!!!      USE LOGICAL_MOD,     ONLY : LLINOZ
!!!      USE NETCDF_UTIL_MOD
!!!      USE TRACER_MOD,      ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
!!!      USE TRANSFER_MOD,    ONLY : TRANSFER_3D
!!!      USE ADJ_ARRAYS_MOD,  ONLY : PROD   , LOSS   
!!!      USE ADJ_ARRAYS_MOD,  ONLY : PROD_0 , LOSS_0
!!!      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL
!!!
!!!
!!!#     include "CMN_SIZE"
!!!!
!!!! !INPUT PARAMETERS: 
!!!!
!!!      ! Arguments
!!!      INTEGER,INTENT(IN) :: THISMONTH
!!!!
!!!!EOP
!!!!------------------------------------------------------------------------------
!!!!BOC
!!!!
!!!! !LOCAL VARIABLES:
!!!!
!!!      CHARACTER(LEN=255) :: FILENAME
!!!      CHARACTER(LEN=6)   :: SPNAME( NTR_GMI )
!!!      TYPE (XPLEX)             :: ARRAY( IIPAR, JJPAR, LGLOB )
!!!      TYPE (XPLEX)             :: ARRAY2( IIPAR, JJPAR, LLPAR )
!!!      INTEGER            :: N, M, S, F
!!!      INTEGER            :: SPNAME_varID
!!!      INTEGER            :: prod_varID, loss_varID
!!!
!!!      !=================================================================
!!!      ! GET_RATES begins here
!!!      !=================================================================
!!!
!!!      ! In the original schem code, the following species were destroyed 
!!!      ! by photolysis in the stratosphere:
!!!      !  PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O,
!!!      !  N2O5, HNO4, MP
!!!      ! And by reaction with OH for:
!!!      !  ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
!!!      !  PRPE, C3H8, CH2O, C2H6, HNO4, MP
!!!      ! The updated code includes at least all of these, and several more.
!!!
!!!      ! Intialize arrays
!!!      LOSS = 0d0
!!!      PROD = 0d0
!!!
!!!      ! Path to input data
!!!      FILENAME = 'gmi.clim.' //
!!!     &     GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
!!!
!!!      FILENAME = TRIM( DATA_DIR ) // 'strat_chem_201106/' //
!!!     &     TRIM( FILENAME )
!!!
!!!      ! Open the netCDF file containing the rates
!!!      WRITE(6,*) 'Reading in monthly stratospheric prod/loss rates'
!!!      call ncdf_open_for_read( ncID_strat_rates, TRIM(FILENAME) )
!!!
!!!      ! Get the variable IDs for the species, prod and loss rates
!!!      prod_varID = ncdf_get_varid( ncID_strat_rates, 'prod' )
!!!      loss_varID = ncdf_get_varid( ncid_strat_rates, 'loss' )
!!!      
!!!      M = THISMONTH
!!!      
!!!      ! Match to strat chem tracers number (hml, 10/19/11)
!!!      !DO N = 1, N_TRACERS
!!!      DO N = 1, NSTPL
!!!
!!!         F = 1
!!!         
!!!         IF ( GMI_TO_GC( N, F ) .eq. 0 ) CYCLE
!!!         
!!!         S = GMI_TO_GC( N, F )
!!!         
!!!         ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!         ! Read production rate [v/v/s]
!!!         call ncdf_get_var( ncID_strat_rates, prod_varID, array,
!!!     &        start=(/     1,     1,     1,   s,  m  /), 
!!!     &        count=(/ iipar, jjpar, lglob,   1,  1  /)  )
!!!         
!!!         ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to 1:LLPAR
!!!         call transfer_3D( array, array2 )
!!!         
!!!         PROD(:,:,:,N) = PROD(:,:,:,N) + TRACER_COEFF(N,F)*ARRAY2
!!!
!!!         ! Save rates from file to respective arrays (hml, 09/15/11)
!!!         PROD_0(:,:,:,N) = PROD(:,:,:,N)
!!!
!!!         ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!         ! Read loss frequency [s-1]
!!!         call ncdf_get_var( ncID_strat_rates, loss_varID, array,
!!!     &        start=(/     1,     1,     1,   s,  m  /), 
!!!     &        count=(/ iipar, jjpar, lglob,   1,  1  /)  )
!!!         
!!!         ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to 1:LLPAR
!!!         call transfer_3D( array, array2 )
!!!         
!!!         LOSS(:,:,:,N) = LOSS(:,:,:,N) + TRACER_COEFF(N,F)*ARRAY2
!!!
!!!         ! Save rates from file to respective arrays (hml, 09/15/11)
!!!         LOSS_0(:,:,:,N) = LOSS(:,:,:,N)
!!!
!!!      ENDDO
!!!
!!!      call ncdf_close( ncID_strat_rates )
!!!
!!!      END SUBROUTINE GET_RATES
!!!!EOC
!!!!------------------------------------------------------------------------------
!!!!BOP
!!!!
!!!! !IROUTINE: init_strat_chem
!!!!
!!!! !DESCRIPTION: Subroutine INIT\_STRAT\_CHEM allocates all module arrays.  
!!!!  It also opens the necessary rate files.
!!!!\\
!!!!\\
!!!! !INTERFACE:
!!!!      
!!!      SUBROUTINE INIT_STRAT_CHEM
!!!!
!!!! !USES:
!!!!
!!!      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
!!!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!!!      USE ERROR_MOD,     ONLY : ALLOC_ERR
!!!      USE LOGICAL_MOD,   ONLY : LLINOZ
!!!      USE NETCDF_UTIL_MOD
!!!      USE TRACER_MOD,    ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
!!!      USE TIME_MOD,      ONLY : GET_TS_CHEM, EXPAND_DATE
!!!      USE TIME_MOD,      ONLY : GET_NYMDb, GET_NHMSb
!!!      !USE ADJ_ARRAYS_MOD,ONLY : PROD   , LOSS   
!!!
!!!#     include "CMN_SIZE"
!!!! 
!!!!EOP
!!!!------------------------------------------------------------------------------
!!!!BOC
!!!!
!!!! !LOCAL VARIABLES:
!!!!     
!!!      INTEGER :: AS
!!!      CHARACTER(LEN=255) :: FILENAME, FILENAMEOUT
!!!      CHARACTER(LEN=6)   :: SPNAME( NTR_GMI )
!!!      INTEGER            :: spname_varID
!!!      INTEGER            :: N, NN, F
!!!      
!!!      !=================================================================
!!!      ! INIT_STRAT_CHEM begins here!
!!!      !=================================================================
!!!      
!!!      ! Allocate PROD -- array for clim. production rates [v/v/s]
!!!      ALLOCATE( PROD( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
!!!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
!!!      PROD = 0d0
!!!      
!!!      ! Allocate LOSS -- array for clim. loss freq [s-1]
!!!      ALLOCATE( LOSS( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
!!!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS' )
!!!      LOSS = 0d0
!!!      
!!!      ! Allocate GMI_TO_GC -- array for mapping
!!!      ALLOCATE( GMI_TO_GC( N_TRACERS, MAX_FM ) ) 
!!!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GMI_TO_GC' )
!!!      GMI_TO_GC = 0
!!!      
!!!      ! Initialize timestep for chemistry
!!!      DTCHEM = GET_TS_CHEM() * 60d0
!!!      
!!!      ! Determine the mapping for the GMI to the GC variables based on
!!!      ! tracer name, which only needs to be done once per model run.
!!!      
!!!      ! Path to input data
!!!      FILENAME = 'gmi.clim.' //  
!!!     &     GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
!!!     
!!!      FILENAME = TRIM( DATA_DIR ) // 'strat_chem_201106/' //
!!!     &     TRIM( FILENAME )
!!!     
!!!      !write(6,*) 'Opening for read: ',trim(filename)
!!!      
!!!      ! Initialize netCDF (this will be moved to main.f for the standard code)
!!!      call NCDF_INIT
!!!      
!!!      ! Open the input netCDF file
!!!      call ncdf_open_for_read( ncID_strat_rates, trim(filename) )
!!!      
!!!      ! Get the variable IDs for the species names
!!!      spname_varid = ncdf_get_varid( ncID_strat_rates,'species_labels' )
!!!
!!!      ! Get the species names and close file
!!!      call ncdf_get_var( ncID_strat_rates, spname_varid, spname )
!!!      call ncdf_close( ncID_strat_rates )
!!!
!!!      WRITE(6,*) "Linearized stratospheric chemistry performed for:"
!!!
!!!      DO N = 1, N_TRACERS
!!!      DO NN = 1, NTR_GMI
!!!
!!!         ! General case
!!!         IF ( TRIM(TRACER_NAME(N)) .eq. TRIM(SPNAME(NN)) ) THEN
!!!
!!!            IF ( LLINOZ .and. TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
!!!               WRITE(6,*) TRIM(TRACER_NAME(N)) // ' (via Linoz)'
!!!            ELSE
!!!               WRITE(6,*) TRIM(TRACER_NAME(N)) // ' (via GMI rates)'
!!!            ENDIF
!!!
!!!            GMI_TO_GC( N, 1 ) = NN
!!!
!!!         ENDIF
!!!
!!!      ENDDO
!!!      ENDDO
!!!      END SUBROUTINE INIT_STRAT_CHEM
!!!!EOC
!!!!------------------------------------------------------------------------------
!!!!BOP
!!!!
!!!! !IROUTINE: cleanup_strat_chem
!!!!
!!!! !DESCRIPTION: Subroutine CLEANUP\_STRAT\_CHEM deallocates all module 
!!!!  arrays.
!!!!\\
!!!!\\
!!!! !INTERFACE:
!!!!
!!!      SUBROUTINE CLEANUP_STRAT_CHEM
!!!!
!!!! !USES:
!!!      USE NETCDF_UTIL_MOD
!!!!EOP
!!!!------------------------------------------------------------------------------
!!!!BOC
!!!
!!!      IF ( ALLOCATED( PROD ) ) DEALLOCATE( PROD )
!!!      IF ( ALLOCATED( LOSS ) ) DEALLOCATE( LOSS )
!!!
!!!      END SUBROUTINE CLEANUP_STRAT_CHEM
!EOC
      END MODULE STRAT_CHEM_ADJ_MOD
