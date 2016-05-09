!$Id: strat_chem_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: strat_chem_mod
!
! !DESCRIPTION: Module STRAT\_CHEM\_MOD contains variables and routines for 
!  performing a simple linearized chemistry scheme for more realistic
!  upper boundary conditions. Archived 3D monthly climatological production
!  rates and loss frequencies are applied from the GMI combo model.
!
!  In the original schem code (schem.f), only the following species
!  were destroyed by photolysis in the stratosphere:
!    PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O, N2O5, HNO4, MP
!  and by reaction with OH for:
!    ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
!    PRPE, C3H8, CH2O, C2H6, HNO4, MP
!  
!  The updated code includes at least all of these, and many more. The code
!  is flexible enough to automatically apply the rate to any new tracers
!  for future simulations that share the name in tracer_mod with the 
!  GMI name.  (See Documentation).
!
!\\
!\\
! !INTERFACE:
!
      MODULE STRAT_CHEM_MOD
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
      PUBLIC  :: DO_STRAT_CHEM
      PUBLIC  :: CLEANUP_STRAT_CHEM

      ! hml 
      PUBLIC  :: PROD_0
      PUBLIC  :: LOSS_0
      PUBLIC  :: PROD
      PUBLIC  :: LOSS
      PUBLIC  :: DTCHEM
      PUBLIC  :: GET_RATES

! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_STRAT_CHEM
    
      ! hml
      !PRIVATE :: GET_RATES
!
! !ADJOINT GROUP:
!
!
! !PUBLIC DATA MEMBERS:
!
! !REMARKS:
!
!  References:
!  ============================================================================
!  (1 )
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!  22 Oct 2011 - H.-M. Lee - Modified to implement in adjoint. 
!                Now we can calculte strat prod and loss sensitivity. adj32_025
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Scalars
      TYPE (XPLEX)               :: DTCHEM

      ! Parameters
      INTEGER, PARAMETER   :: NTR_GMI = 120 ! Number of species
                              ! 118 as output from GMI + NOx + Ox families

      INTEGER, PARAMETER   :: MAX_FM  = 1 ! Max number of species in a fam
      ! Vestigial, as NOx and Ox families pre-processed, but may be useful
      ! for future uses, e.g., ClOx.

      ! Arrays
      TYPE (XPLEX),  ALLOCATABLE :: PROD(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: LOSS(:,:,:,:)
      INTEGER, ALLOCATABLE :: GMI_TO_GC(:,:)
      INTEGER, SAVE        :: ncID_strat_rates
     
      ! hml 
      TYPE (XPLEX),  ALLOCATABLE :: PROD_0(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: LOSS_0(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DO_STRAT_CHEM
!
! !DESCRIPTION: Function DO\_STRAT\_CHEM is the driver routine for computing
!     the simple linearized stratospheric chemistry scheme for a host of species
!     whose prod/loss rates were determined from the GMI combo model. Ozone is
!     treated using either Linoz or Synoz.
!\\
!\\
! !INTERFACE:
!      
      SUBROUTINE DO_STRAT_CHEM
!
! !USES:
!
      USE DAO_MOD,        ONLY : AD, CONVERT_UNITS
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE LOGICAL_MOD,    ONLY : LLINOZ, LPRT
      USE LINOZ_MOD,      ONLY : DO_LINOZ
      USE NETCDF_UTIL_MOD
      USE UPBDFLX_MOD,    ONLY : UPBDFLX_O3, INIT_UPBDFLX
      USE TIME_MOD,       ONLY : GET_MONTH, TIMESTAMP_STRING
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACER_MOD,     ONLY : N_TRACERS, STT, TCVV, TRACER_MW_KG
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL, ITS_IN_THE_STRAT

      ! adj_group (hml, 07/25/11) 
      USE ADJ_ARRAYS_MOD, ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD, ONLY : NSTPL 
      USE TRACER_MOD,     ONLY : STT_STRAT_TMP 
      USE LOGICAL_ADJ_MOD,ONLY : LADJ
      USE LOGICAL_ADJ_MOD,ONLY : LADJ_STRAT
      USE CHECKPOINT_MOD, ONLY : MAKE_BEFSTRAT_CHKFILE
      USE TIME_MOD,       ONLY : GET_NHMS
      USE TIME_MOD,       ONLY : GET_NYMD
      USE TIME_MOD,       ONLY : GET_TAU


#     include "CMN_SIZE"
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  1 Feb 2011 - L. Murray - Initial version  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE             :: FIRST = .TRUE.
      INTEGER, SAVE             :: LASTMONTH = -999
      INTEGER, SAVE             :: LASTSEASON = -1
      INTEGER                   :: I, J, L, N, LMIN
      INTEGER                   :: IORD, JORD, KORD
      INTEGER                   :: NHMS
      INTEGER                   :: NYMD
      TYPE (XPLEX)                    :: TAU
      TYPE (XPLEX)                    :: t, P, k, M0
      CHARACTER(LEN=16)         :: STAMP


      !===============================
      ! DO_STRAT_CHEM begins here!
      !===============================

      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DO_STRAT_CHEM: Strat chemistry at ', a )

      IF ( FIRST ) THEN

         ! Allocate all module arrays
         CALL INIT_STRAT_CHEM

#if    defined( GEOS_3 )
         ! Initialize some Synoz variables
         IF ( .NOT. ( LLINOZ ) ) THEN
            CALL GET_ORD( IORD, JORD, KORD )
            CALL INIT_UPBDFLX( IORD, JORD, KORD )
         ENDIF
#endif

      ENDIF

      !================================================
      ! Determine the rates from disk; merge families
      !================================================

      ! Get the minimum level extent of the tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()  

      IF ( GET_MONTH() /= LASTMONTH ) THEN

         WRITE(6,*) 'Getting new strat rates for month: ',GET_MONTH()

         IF ( LPRT ) CALL DEBUG_MSG( '### STRAT_CHEM: at GET_RATES' )

         ! Read rates for this month
         CALL GET_RATES( GET_MONTH() )
         
         ! Save month for next iteration
         LASTMONTH = GET_MONTH()
      ENDIF

      ! Set first-time flag to false
      FIRST = .FALSE.    

      IF ( LPRT ) CALL DEBUG_MSG( '### STRAT_CHEM: at DO_STRAT_CHEM' )

      WRITE(6,*) '-----------------------------------------------------'
      write(6,*) '    Doing stratospheric chemistry (STRAT_CHEM_MOD)   '
      WRITE(6,*) '-----------------------------------------------------'

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, k, P, t, M0 )
!$OMP+SCHEDULE( DYNAMIC )
      DO L=LMIN,LLPAR
      DO J=1,JJPAR
      DO I=1,IIPAR
      IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

         ! To loop over only tracers that we have prod and loss rates (hml)
         !DO N=1,N_TRACERS
         DO N=1,NSTPL

            ! Include something to expediate skipping past species
            ! that we do not have strat chem for. Prob put tracer on
            ! outermost loop.

            ! Now we will use GMI rate for Ox if LINOZ is off (hml, 10/31/11)  
!            ! Skip Ox -- if we're not using Linoz, we'll use Synoz below
!            IF ( ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGOX_SIM() ) .and.
!     &           ( N .eq. IDTOx ) ) CYCLE
            IF ( ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGOX_SIM() ) .and.
     &           ( LLINOZ ) .and. ( N .eq. IDTOx ) ) CYCLE

            ! adj_group:  make a version that applies scaling factors 
            ! and use this if the stratosphere adjoint ID #'s are active 
            IF ( LADJ_STRAT ) THEN
  
               ! Check point values of STT
               STT_STRAT_TMP(I,J,L,N) = STT(I,J,L,N)
            
               PROD(I,J,L,N) = PROD_0(I,J,L,N) * PROD_SF(I,J,1,N)
               LOSS(I,J,L,N) = LOSS_0(I,J,L,N) * LOSS_SF(I,J,1,N)

            ENDIF

            !===============================
            ! Do chemical production and loss
            !===============================

            t = DTCHEM                              ! timestep [s]
            k = LOSS(I,J,L,N)                       ! loss freq [s-1]
            P = PROD(I,J,L,N) * AD(I,J,L) / TCVV(N) ! production term [kg s-1]
            M0 = STT(I,J,L,N)                       ! initial mass [kg]

            ! debug test 
            !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
            !   print*, N,' STRAT TEST fwd: k = ', k 
            !   print*, N,' STRAT TEST fwd: P = ', P
            !   print*, N,' STRAT TEST fwd: M0= ', M0
            !ENDIF 

            ! No prod or loss at all
            if ( k .eq. 0d0 .and. P .eq. 0d0 ) cycle

            ! Simple analytic solution to dM/dt = P - kM over [0,t]
            if ( k .gt. 0d0 ) then
               STT(I,J,L,N) = M0 * exp(-k*t) + (P/k)*(1d0-exp(-k*t))
            else
               STT(I,J,L,N) = M0 + P*t
            endif

         ENDDO

      ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Make check point file
      IF ( LADJ ) THEN 
         NHMS     = GET_NHMS()
         NYMD     = GET_NYMD()
         TAU      = GET_TAU()
         CALL MAKE_BEFSTRAT_CHKFILE( NYMD, NHMS, TAU )
      ENDIF 

      !======================
      ! Stratospheric Ozone !
      !======================
      ! Modified (hml, 10/31/11)
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Put ozone in v/v
         STT(:,:,:,IDTOX ) = STT(:,:,:,IDTOX) * TCVV( IDTOX ) / AD
         
         IF ( LLINOZ ) THEN 
            CALL DO_LINOZ       ! Linoz
         ELSE
            ! must use Linoz or strat chem Ox fluxes for the adjoint 
            IF ( .not. LADJ ) THEN
               CALL UPBDFLX_O3     ! Synoz
            ENDIF 
         ENDIF
         
         ! Now move unit conversion into LINOZ (hml, 11/06/11)
         ! Put ozone back to kg
         STT(:,:,:,IDTOX) = STT(:,:,:,IDTOX) * AD / TCVV( IDTOX )

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v
         
         IF ( LLINOZ ) THEN
            CALL DO_LINOZ       ! Linoz
         ELSE
            ! must use Linoz or strat chem Ox fluxes for the adjoint 
            IF ( .not. LADJ ) THEN
               CALL UPBDFLX_O3     ! Synoz
            ENDIF 
         ENDIF
         
         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg
         
      ENDIF

      END SUBROUTINE DO_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_RATES
!
! !DESCRIPTION: Function GET\_RATES reads from disk the chemical production
!  and loss rates for the species of interest
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_RATES( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_NAME_EXT, GET_RES_EXT
      USE DIRECTORY_MOD,   ONLY : DATA_DIR
      USE LOGICAL_MOD,     ONLY : LLINOZ
      USE NETCDF_UTIL_MOD
      USE TRACER_MOD,      ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
      USE TRANSFER_MOD,    ONLY : TRANSFER_3D
      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL


#     include "CMN_SIZE"
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER,INTENT(IN) :: THISMONTH
!
! !REVISION HISTORY: 
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=6)   :: SPNAME( NTR_GMI )
      REAL*8                 :: D_ARRAY(IIPAR,JJPAR,LGLOB)
      TYPE (XPLEX)             :: ARRAY( IIPAR, JJPAR, LGLOB )
      !REAL*8             :: ARRAY( IIPAR, JJPAR, LGLOB )
      TYPE (XPLEX)             :: ARRAY2( IIPAR, JJPAR, LLPAR )
      INTEGER            :: N, M, S, F
      INTEGER            :: SPNAME_varID
      INTEGER            :: prod_varID, loss_varID

      !=================================================================
      ! GET_RATES begins here
      !=================================================================

      ! In the original schem code, the following species were destroyed 
      ! by photolysis in the stratosphere:
      !  PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O,
      !  N2O5, HNO4, MP
      ! And by reaction with OH for:
      !  ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
      !  PRPE, C3H8, CH2O, C2H6, HNO4, MP
      ! The updated code includes at least all of these, and several more.

      ! Initialize arrays
      LOSS = 0d0
      PROD = 0d0
      ! For adjoint (hml)
      !LOSS_0 = 0d0
      !PROD_0 = 0d0

      ! Path to input data
      FILENAME = 'gmi.clim.' // 
     &     GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
      
      FILENAME = TRIM( DATA_DIR ) // 'strat_chem_201106/' //
     &     TRIM( FILENAME )

      ! Open the netCDF file containing the rates
      WRITE(6,*) 'Reading in monthly stratospheric prod/loss rates'
      call ncdf_open_for_read( ncID_strat_rates, TRIM(FILENAME) )

      ! Get the variable IDs for the species, prod and loss rates
      prod_varID = ncdf_get_varid( ncID_strat_rates, 'prod' )
      loss_varID = ncdf_get_varid( ncid_strat_rates, 'loss' )

      M = THISMONTH

      ! Match to strat chem tracers number (hml, 10/19/11)
      !DO N = 1, N_TRACERS
      DO N = 1, NSTPL
      !DO F = 1, 1 !MAX_FM
         F = 1

         IF ( GMI_TO_GC( N, F ) .eq. 0 ) CYCLE

         S = GMI_TO_GC( N, F )

         ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Read production rate [v/v/s]
         call ncdf_get_var( ncID_strat_rates, prod_varID, d_array, 
     &        start=(/     1,     1,     1,   s,  m  /), 
     &        count=(/ iipar, jjpar, lglob,   1,  1  /)  )
         ARRAY(:,:,:) = D_ARRAY(:,:,:)
         ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to 1:LLPAR
         call transfer_3D( array, array2 )

         PROD(:,:,:,N) = PROD(:,:,:,N) + TRACER_COEFF(N,F)*ARRAY2
         
         ! Save rates from file to respective arrays (hml, 09/15/11)
         PROD_0(:,:,:,N) = PROD(:,:,:,N)
 
         ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Read loss frequency [s-1]
         call ncdf_get_var( ncID_strat_rates, loss_varID, d_array,
     &        start=(/     1,     1,     1,   s,  m  /), 
     &        count=(/ iipar, jjpar, lglob,   1,  1  /)  )
         ARRAY(:,:,:) = D_ARRAY(:,:,:)
         ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to 1:LLPAR
         call transfer_3D( array, array2 )

         LOSS(:,:,:,N) = LOSS(:,:,:,N) + TRACER_COEFF(N,F)*ARRAY2

         ! Save rates from file to respective arrays (hml, 09/15/11)
         LOSS_0(:,:,:,N) = LOSS(:,:,:,N)
         
      !ENDDO

      ENDDO

      call ncdf_close( ncID_strat_rates )
      
      END SUBROUTINE GET_RATES
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_strat_chem
!
! !DESCRIPTION: Subroutine INIT\_STRAT\_CHEM allocates all module arrays.  
!  It also opens the necessary rate files.
!\\
!\\
! !INTERFACE:
!      
      SUBROUTINE INIT_STRAT_CHEM
!
! !USES:
!
      USE BPCH2_MOD,   ONLY : GET_NAME_EXT, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LLINOZ
      USE NETCDF_UTIL_MOD
      USE TRACER_MOD,  ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
      USE TIME_MOD,    ONLY : GET_TS_CHEM, EXPAND_DATE
      USE TIME_MOD,    ONLY : GET_NYMDb, GET_NHMSb

#     include "CMN_SIZE"
! 
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      CHARACTER(LEN=255) :: FILENAME, FILENAMEOUT
      CHARACTER(LEN=6)   :: SPNAME( NTR_GMI )
      INTEGER            :: spname_varID
      INTEGER            :: N, NN, F

      !=================================================================
      ! INIT_STRAT_CHEM begins here!
      !=================================================================

      ! Allocate PROD -- array for clim. production rates [v/v/s]
      ALLOCATE( PROD( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
      PROD = 0d0

      ALLOCATE( PROD_0( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD_0' )
      PROD_0 = 0d0

      ! Allocate LOSS -- array for clim. loss freq [s-1]
      ALLOCATE( LOSS( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS' )
      LOSS = 0d0

      ALLOCATE( LOSS_0( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS_0' )
      LOSS_0 = 0d0

      ! Allocate GMI_TO_GC -- array for mapping
      ALLOCATE( GMI_TO_GC( N_TRACERS, MAX_FM ) )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GMI_TO_GC' )
      GMI_TO_GC = 0

      ! Initialize timestep for chemistry
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Determine the mapping for the GMI to the GC variables based on
      ! tracer name, which only needs to be done once per model run.
      
      ! Path to input data
      FILENAME = 'gmi.clim.' // 
     &     GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'

      FILENAME = TRIM( DATA_DIR ) // 'strat_chem_201106/' // 
     &     TRIM( FILENAME )

      !write(6,*) 'Opening for read: ',trim(filename)

      ! Initialize netCDF (this will be moved to main.f for the standard code)
      call NCDF_INIT

      ! Open the input netCDF file
      call ncdf_open_for_read( ncID_strat_rates, trim(filename) )

      ! Get the variable IDs for the species names
      spname_varid = ncdf_get_varid( ncID_strat_rates,'species_labels' )

      ! Get the species names and close file
      call ncdf_get_var( ncID_strat_rates, spname_varid, spname )
      call ncdf_close( ncID_strat_rates )

      WRITE(6,*) "Linearized stratospheric chemistry performed for:"

      DO N = 1, N_TRACERS
      DO NN = 1, NTR_GMI

         ! General case
         IF ( TRIM(TRACER_NAME(N)) .eq. TRIM(SPNAME(NN)) ) THEN

            IF ( LLINOZ .and. TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
               WRITE(6,*) TRIM(TRACER_NAME(N)) // ' (via Linoz)'

            ! Debug, hml
            ELSEIF ( TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
               WRITE(6,*) TRIM(TRACER_NAME(N)) // ' (Ox via GMI)'

            ELSE
               WRITE(6,*) TRIM(TRACER_NAME(N)) // ' (via GMI rates)'
            ENDIF

            GMI_TO_GC( N, 1 ) = NN

         ENDIF              

      ENDDO
      ENDDO

      !DO N=1,N_TRACERS
      !   print*,N,TRACER_NAME(N),TRACER_COEFF(N,:)
      !   print*,'-->',GMI_TO_GC(n,:)
      !ENDDO

      END SUBROUTINE INIT_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_strat_chem
!
! !DESCRIPTION: Subroutine CLEANUP\_STRAT\_CHEM deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_STRAT_CHEM
!
! !USES:
      USE NETCDF_UTIL_MOD
! !REVISION HISTORY: 
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      
      IF ( ALLOCATED( PROD   ) ) DEALLOCATE( PROD   )
      IF ( ALLOCATED( LOSS   ) ) DEALLOCATE( LOSS   )
      IF ( ALLOCATED( PROD_0 ) ) DEALLOCATE( PROD_0 )
      IF ( ALLOCATED( LOSS_0 ) ) DEALLOCATE( LOSS_0 )

      END SUBROUTINE CLEANUP_STRAT_CHEM
!EOC
      END MODULE STRAT_CHEM_MOD
