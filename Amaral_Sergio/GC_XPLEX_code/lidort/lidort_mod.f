!$Id: lidort_mod.f,v 1.5 2012/05/09 22:31:56 nicolas Exp $
      MODULE LIDORT_MOD
!
!******************************************************************************
!  LIDORT_MOD contains routines for calculating aerosol optical properties
!  and there derivaties. (dkh, 01/17/08) 
!
!  Module Variables:
!  ============================================================================
!  (1 )                   (INTEGER)   : 
!
!  NOTES
!  (1 ) Updated from GCv6, renamed LIDORT_MOD from AERO_OPTICAL_MOD. (dkh, 05/18/10) 
!  (2 ) Updated to LIDORT v3p5 (dkh, 07/29/10) 
!******************************************************************************
!
      USE ADJ_ARRAYS_MOD, ONLY : NSPAN
      USE ADJ_ARRAYS_MOD, ONLY : INV_NSPAN

      IMPLICIT NONE

      !  include file for LIDORT Mk 2 threads
#include      "./thread_sourcecode_MkII_F90/LIDORT.PARS_F90"

      PRIVATE 
      PUBLIC :: CALC_RF_FORCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      ! Parameters 
      INTEGER, PARAMETER  :: IU_AOD  = 523
      INTEGER, PARAMETER  :: IU_LID  = 524
      INTEGER, PARAMETER  :: IU_RAD  = 525
      LOGICAL, PARAMETER  :: L_PRINT = .FALSE.

      ! Indicies for MIE TABLE
      INTEGER, PARAMETER  :: NWL_MAX  = 5   ! Maximum number of wavelengths
	  !fha  (6/1/11) added OC
      INTEGER, PARAMETER  :: NSP_MAX  = 3   ! Maximum number of species
      INTEGER, PARAMETER  :: NOP_MAX  = 2   ! Maximum number of optical properties
      INTEGER, PARAMETER  :: NRA_MAX  = 20  ! Maximum number of radii
      !INTEGER, PARAMETER  :: NPM_MAX  = 40  ! Maximum number of phase moments
      INTEGER, PARAMETER  :: NPM_MAX  = 20  ! Maximum number of phase moments
      INTEGER, PARAMETER  :: NPM_ACT  = 5   ! Number of active phase moments


      INTEGER, PARAMETER  :: NSP_SULF = 1   ! mie species index for sulfate
      INTEGER, PARAMETER  :: NSP_BC   = 2   ! mie species index for black carbon
	  !fha  (6/1/11) added OC
      INTEGER, PARAMETER  :: NSP_OC   = 3   ! mie species index for organic carbon
      INTEGER, PARAMETER  :: NOP_BEXT = 1   ! mie property index of extinction coeffs
      INTEGER, PARAMETER  :: NOP_SSA  = 2   ! mie property index of single scattering albedo

      REAL*8,  PARAMETER  :: ABS_FAC  = 1.5d0 ! absorption enhancement.  See Bond and 
                                              ! Bergstrom, 2006; Koch et al., 2009. 
         

      ! Number of TEMIS LER wavelenths
      INTEGER, PARAMETER  :: LLER  = 11

      ! Data
      ! Table of lambdas considered [nm]
      REAL*8              :: LAMBDA_TABLE(NWL_MAX)
      DATA                   LAMBDA_TABLE  / 315d0, 357d0, 500d0,    
     &                                       833d0, 1667d0        / 
      INTEGER, PARAMETER  :: NWL_500  = 3   
 
      ! particle size range. updated to GCv8-03-01
	  ! fha - 6-7-11 added OC range from jv_spec.dat; changed sulfate from 0.25 to 0.23
      REAL*8              :: RMIN(NSP_MAX)
      DATA                   RMIN / 0.07d0,  0.02D0,  0.06D0 /
      REAL*8              :: RMAX(NSP_MAX)
      DATA                   RMAX / 0.23d0,  0.04D0,  0.16D0 /

      ! From integrating the Plank function using T=5796 K (for 6.4d7 W/m2 
      ! emissions from the sun) then assuming a mean distance to earch
      ! of 1.5d11 m. This gives a solar constant of 1367 W / m2. 
      ! note: given the bands here, we're only considering 1254.4 W/m2
      REAL*8                :: SOL_FLUX_TABLE(NWL_MAX)
      DATA                     SOL_FLUX_TABLE  / 9.7635d00,  1.0606d02, 
     &                                           5.0274d02,  4.3814d02, 
     &                                           1.9769d02  /

      ! Parameters 
      REAL*8, PARAMETER :: SMALL = 1d-20

      ! Microphyical Data
      ! size distribution properties
      REAL*8            :: DRY_DIAM(NSP_MAX)
      ! updated to GCv8-03-01
      !DATA                 DRY_DIAM / 0.1d0, 0.02d0 / 
	  !fha - 6-7-11 added OC size (2x RMIN above)
      DATA                 DRY_DIAM / 0.14d0, 0.04d0, 0.12d0 / 
      REAL*8            :: SIGMA(NSP_MAX)
      ! updated to GCv8-03-01
	  !fha - 6-7-11 added OC sigma
      !DATA                 SIGMA    / 2D0,   2D0    /
      DATA                 SIGMA    / 1.6D0,   1.6D0,   1.6D0  /

      ! refractive index.  10^(-7.5) = -3.1622776d-08, 10^(-3.5) = 0.00031622776d0 
      ! From Martin 2004, Table 1
      REAL*8  :: REFRAC_INDEX(NSP_MAX,NWL_MAX,2)
      DATA       REFRAC_INDEX(NSP_SULF,1,:) / 1.46d0, 3.1622776d-08   /
      DATA       REFRAC_INDEX(NSP_SULF,2,:) / 1.46d0, 3.1622776d-08   /
      DATA       REFRAC_INDEX(NSP_SULF,3,:) / 1.46d0, 3.1622776d-08   /
      DATA       REFRAC_INDEX(NSP_SULF,4,:) / 1.40d0, 1.d-07          /
      DATA       REFRAC_INDEX(NSP_SULF,5,:) / 1.40d0, 0.00031622776d0 /
      ! BC props from ??
      DATA       REFRAC_INDEX(NSP_BC,  1,:) / 1.743d0, 0.469d0        /
      DATA       REFRAC_INDEX(NSP_BC,  2,:) / 1.750d0, 0.464d0        /
      DATA       REFRAC_INDEX(NSP_BC,  3,:) / 1.750d0, 0.450d0        /
      DATA       REFRAC_INDEX(NSP_BC,  4,:) / 1.750d0, 0.432d0        /
      DATA       REFRAC_INDEX(NSP_BC,  5,:) / 1.783d0, 0.473d0        /
      !fha 5-26-2011
      ! OC props from d'Almeida, Koepke, Shettle 1991 Atmospheric Aerosols Global climatology and Radiative Characteristics p. 57
      DATA       REFRAC_INDEX(NSP_OC,  1,:) / 1.53d0, -0.007d0        /
      DATA       REFRAC_INDEX(NSP_OC,  2,:) / 1.53d0, -0.00504d0      /
      DATA       REFRAC_INDEX(NSP_OC,  3,:) / 1.53d0, -0.0058d0       /
      DATA       REFRAC_INDEX(NSP_OC,  4,:) / 1.53d0, -0.0105d0       /
      DATA       REFRAC_INDEX(NSP_OC,  5,:) / 1.52d0, -0.01916d0      /
      ! Ion density from Martin et al., 2004 (ACP).  
	  !fha 6-1-11 I added a fake 1d0 for OC here
      REAL*8            :: DENSE(NSP_MAX)
      DATA                 DENSE    / 1.75d0, 1d0, 1d0   /

      ! note: need to update to Daumont-Malicet cross sections

      ! Factor used for averaging results from 1 day using 24 1hr time steps
      ! NSPAN is now defined in input.gcadj
      !INTEGER, PARAMETER  :: NSPAN     = 1
      ! Now reference this from adj_arrays_mod.f
      !REAL*8,  PARAMETER  :: INV_NSPAN = 1.d0 / REAL(NSPAN,8)


      ! Allocatable variables
      REAL*8, ALLOCATABLE :: AOD_SAVE(:,:,:)
      REAL*8, ALLOCATABLE :: AOD_AVE(:,:,:)
      REAL*8, ALLOCATABLE :: LAYER_AOD(:,:,:,:)
      REAL*8, ALLOCATABLE :: LAYER_SSA(:,:,:,:)
      REAL*8, ALLOCATABLE :: LAYER_PHM(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: MASS_ND(:,:)
      REAL*8, ALLOCATABLE :: TEMIS_LER(:,:,:)

      REAL*8, ALLOCATABLE :: GLOB_FLUX_W(:,:,:)
      REAL*8, ALLOCATABLE :: JAC_GLOB_FLUX_W(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: GLOB_FLUX(:,:)
      REAL*8, ALLOCATABLE :: GLOB_AVE_FORCE(:,:)

      REAL*8, ALLOCATABLE :: GLOB_FLUX_W_ADJ(:,:,:)
      REAL*8, ALLOCATABLE :: LAYER_AOD_ADJ(:,:,:,:)
      REAL*8, ALLOCATABLE :: LAYER_SSA_ADJ(:,:,:,:)
      REAL*8, ALLOCATABLE :: LAYER_PHM_ADJ(:,:,:,:,:)

      REAL*8, ALLOCATABLE :: DEBUG_J(:,:)


      !================================================================
      ! Module variables for LIDORT.  These are the arguements of the 
      ! lidort_inputs_basic_master which we only need to load once per
      ! simulation and are independent of the model grid cell. 
      !================================================================

      !  All inputs
      !  ----------

      !  1. CONTROL FLAGS
      !  ================

      !  Basic top-level control
      !   -- Thermal to be reintroduced later 2010
      LOGICAL ::    DO_SOLAR_SOURCES
!      LOGICAL ::    DO_THERMAL_EMISSION

      !  directional control
      LOGICAL ::    DO_UPWELLING
      LOGICAL ::    DO_DNWELLING

      !  Flag for Full Radiance  calculation
      !    If not set, just produce diffuse (Multiple scatter) field
      LOGICAL ::    DO_FULLRAD_MODE

      !  stream angle flag. Normally required for post-processing solutions
      !    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)
      LOGICAL ::    DO_USER_STREAMS

      !  single scatter corrections. (Includes direct beam reflection)
      !    - NADIR    : Plane-parallel line-of-sight
      !    - OUTGOING : Line-of-sight in curved atmosphere
      LOGICAL ::    DO_SSCORR_NADIR         ! May be re-set after Checking
      LOGICAL ::    DO_SSCORR_OUTGOING      ! May be re-set after Checking

      !  Flag for Full-up single scatter calculation. (No MS field)
      !    One of the above two SSCORR flags must be set
      LOGICAL    ::    DO_SSFULL

      !  Flag for performing a complete separate delta-M truncation on the
      !  single scatter corrrection  calculations. **** Use with CAUTION.
      LOGICAL ::    DO_SSCORR_TRUNCATION         ! May be re-set after Checking

      !  Beam particular solution, plane parallel flag
      !    - Not normally required; pseudo-spherical if not set
      LOGICAL ::    DO_PLANE_PARALLEL

      !  Flag for use of BRDF surface
      !    - If not set, default to Lambertian surface
      LOGICAL ::    DO_BRDF_SURFACE

      !  mean value control (1). If set --> Flux output AS WELL AS Intensities
      LOGICAL ::    DO_ADDITIONAL_MVOUT     ! May be re-set after Checking

      !  mean value control (2). If set --> only Flux output (No Intensities)
      !    - DO_USER_STREAMS should be turned off
      LOGICAL ::    DO_MVOUT_ONLY           ! May be re-set after Checking

      !  Beam particular solution: Flag for calculating solar beam paths
      !    ( Chapman factors = slant/vertical path-length ratios)
      !     - This should normally be set. 
      LOGICAL ::    DO_CHAPMAN_FUNCTION     ! May be re-set after Checking

      !  Beam particular solution: Flag for using refraction in solar paths
      !     - This should NOT normally be set. 
      LOGICAL ::    DO_REFRACTIVE_GEOMETRY  ! May be re-set after Checking

      !  Flag for Use of Delta-M scaling
      !    - Should normally be set
      !    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY
      LOGICAL ::    DO_DELTAM_SCALING       ! May be re-set after Checking

      !  double convergence test flag
      LOGICAL ::    DO_DOUBLE_CONVTEST      ! May be re-set after Checking

      !  Performance flags
      !    -- SOLUTION_SAVING gets rid of unneeded RTE computations
      !    -- BVP_TELESCOPING creates reduced Boundary value problems
      !    -- These flags should be used with CAUTION
      !    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers
      LOGICAL ::    DO_SOLUTION_SAVING      ! May be re-set after Checking
      LOGICAL ::    DO_BVP_TELESCOPING      ! May be re-set after Checking

      !  scatterers and phase function control
      !    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
      !    - Isotropic only, if set, phase function is 1
      LOGICAL ::    DO_RAYLEIGH_ONLY        ! May be re-set after Checking
      LOGICAL ::    DO_ISOTROPIC_ONLY       ! May be re-set after Checking

      !  Debug and testing flags
      !   - Normally should not be set
      LOGICAL ::    DO_NO_AZIMUTH           ! May be re-set after Checking
      LOGICAL ::    DO_ALL_FOURIER

      !  Control linearization
      LOGICAL         ::  DO_PROFILE_LINEARIZATION
      LOGICAL         ::  DO_SURFACE_LINEARIZATION

      !  2. CONTROL INTEGERS
      !  ===================
      
      !  Number of discrete ordinate streams
      INTEGER :: NSTREAMS

      !  number of computational layers
      INTEGER :: NLAYERS

      !  Number of fine layers subdividing all computational layers
      !    ( Only required for the outgoing single scattering correction )
      INTEGER :: NFINELAYERS

      !  number of Legendre phase function expansion moments
      INTEGER :: NMOMENTS_INPUT          ! May be re-set after Checking

      !  number of solar beams to be processed
      INTEGER :: NBEAMS

      !  Number of user-defined relative azimuths
      INTEGER :: N_USER_RELAZMS

      !  Number of User-defined viewing zenith angles (0 to 90 degrees)
      INTEGER :: N_USER_STREAMS

      !  Number of User-defined vertical levels for  output
      INTEGER :: N_USER_LEVELS

      !  Linearization control
      !  ---------------------
      
      !  Control for atmospheric linearizations, layer by layer
      LOGICAL         ::  LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER         ::  LAYER_VARY_NUMBER(MAXLAYERS)

      !  Total number of Surface Jacobians
      INTEGER         ::  N_SURFACE_WFS

      !  3. CONTROL NUMBERS
      !  ==================
      
      !  Flux factor ( should be 1 or pi ). Same for all beams.
      REAL(KIND=8) :: FLUX_FACTOR

      !  accuracy for convergence of Fourier series
      REAL(KIND=8) :: LIDORT_ACCURACY

      !  Zenith tolerance (nearness of output zenith cosine to 1.0 )
      !    removed 02 June 2010
      !      REAL(KIND=8) :: ZENITH_TOLERANCE
      
      !  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT
      REAL(KIND=8) :: EARTH_RADIUS

      !  Refractive index parameter
      !  ( Only required for refractive geometry attenuation of the solar beam)
      REAL(KIND=8) :: RFINDEX_PARAMETER

      !  Surface height [km] at which Input geometry is to be specified.
      !    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
      !    -- See special note below
      REAL(KIND=8) :: GEOMETRY_SPECHEIGHT

      !  BOA solar zenith angles (degrees)
      REAL(KIND=8) :: BEAM_SZAS ( MAXBEAMS )
      
      !  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)
      REAL(KIND=8) :: USER_RELAZMS  (MAX_USER_RELAZMS)

      !  User-defined viewing zenith angles input (degrees) 
      REAL(KIND=8) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

      !  User-defined vertical levels for output
      !    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
      !    E.g. For 4.5, this means half way down the 5th layer
      !    E.g. For 0.0, this is output at TOA
      REAL(KIND=8) :: USER_LEVELS   (MAX_USER_LEVELS)


      !  Exception handling 
      !  ------------------
      
      !  Exception handling for File-read Checking. New code, 18 May 2010
      !     Message Length should be at least 120 Characters
      INTEGER         :: STATUS_INPUTREAD
      INTEGER         :: NREADMESSAGES
      CHARACTER*120   :: READMESSAGES(0:MAX_MESSAGES)
      CHARACTER*120   :: READACTIONS (0:MAX_MESSAGES)


      !  ============================================================
      !  ============================================================
      
      !  SPECIAL NOTE on variable GEOMETRY_SPECHEIGHT
      
      !    This is only required when the outgoing sphericity correction is
      !    in operation. Otherwise, the regular pseudo-spherical correction
      !    (wiht or without an exact single-scatter correction) uses the same
      !    set of angles all the up the nadir from the bottom of atmosphere.
      
      !     This height is normally set equal to the height at the lowest
      !     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
      !     In this case, no adjustment to the geometrical inputs is needed
      !     for the outgoing sphericity correction.
      
      !     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
      !     then an adjustment to the geometrical inputs is needed for the
      !     outgoing sphericity correction. This adjustment is internal and
      !     the model output will still be given at the geometrical angles
      !     as specified by the user, even though these angles may not be the
      !     ones at which the calculations were done. This situation will occur
      !     when we are given a BOA geometry but we want to make a calculation
      !     for a reflecting surface (such as a cloud-top) which is above the
      !     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
      !     height HEIGHT_GRID(NLAYERS) = cloud-top.
            
      !     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
      !     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
      !     the calculation will go through without the adjustment. A warning
      !     about this incorrect input choice will be sent to LOGFILE.
      
      !  ============================================================
      !  ============================================================

      !================================================================= 
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!------------------------------------------------------------------------------


      SUBROUTINE CALC_RF_FORCE( COST_FUNC, N_CALC )
!
!******************************************************************************
!  Subroutine CALC_RF_FORCE caculates contribution to cost function from 
!  SW aerosol optical effects. 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL)    : Cost funciton           
!  (2 ) N_CALC    (INTEGER) : Current iteration number 
!     
!     
!  NOTES:
!  (1 ) Updated to work with LIDORT 3p5 and GCv8 adjoint (dkh, 07/29/10) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD,  LFD,  NFD
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : BXHEIGHT,  SUNCOS
      USE ERROR_MOD,      ONLY : IT_IS_NAN, ERROR_STOP
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE CHECKPT_MOD,    ONLY : RP_OUT,    CHK_STT
      USE TIME_MOD,       ONLY : GET_TAUe,  GET_TAUb, GET_TS_CHEM
      USE TIME_MOD,       ONLY : GET_NHMS,  GET_NHMSb
      USE TIME_MOD,       ONLY : GET_NYMD,  GET_NYMDb
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,       ONLY : GET_MONTH
      USE TRACER_MOD,     ONLY : STT
      USE TRACERID_MOD,   ONLY : IDTSO4,    IDTNIT, IDTNH4
      USE TRACERID_MOD,   ONLY : IDTBCPI,   IDTBCPO
	  !fha 6-9-2011 added for debugging output
      USE TRACERID_MOD,   ONLY : IDTOCPI,   IDTOCPO
	  
      USE MIE_MOD,        ONLY : INIT_MIE_MOD

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT

      ! Arguments
      REAL*8, INTENT(INOUT)   :: COST_FUNC
      INTEGER,INTENT(IN)      :: N_CALC
    
      ! Local variables 
      REAL*8                  :: LOCAL_COST(IIPAR,JJPAR,NWL_MAX)
      REAL*8                  :: AOD(NWL_MAX)
      REAL*8                  :: GLOB_FLUX_INST(IIPAR,JJPAR)
      REAL*8                  :: GLOB_FLUX_INST_ADJ(IIPAR,JJPAR)
      INTEGER                 :: I, J, L, IJLOOP
      INTEGER                 :: CLOCK_1, CLOCK_2, CLOCK_RATE
      INTEGER                 :: status_inputread
      LOGICAL, SAVE           :: FIRST = .TRUE. 
      LOGICAL, SAVE           :: NEED_READ_LER = .TRUE.
      !INTEGER, SAVE           :: NCOUNT 
      INTEGER, SAVE           :: NCOUNT = 0 

      integer                 :: n, nf, na, LEN_STRING

      ! Parameters 
      REAL*8, PARAMETER       :: EARTHSA = 510705155749195  ! [m2]

      !=================================================================
      ! CALC_RF_FORCE begins here!
      !=================================================================

      ! dkh debug 
      !print*, 'ddd loc z ', CHK_STT(IFD,JFD,LFD,NFD)
      !print*, 'ddd loc z check', SUM(CHK_STT(IFD,JFD,:,:)) 
      !                          - CHK_STT(IFD,JFD,LFD,NFD)

 
      ! This is now handled in CALC_ADJ_FORCE_FOR_SENSE (dkh, 02/15/11) 
      !NCOUNT = NCOUNT + 1 
      !print*, ' CALC_RF_FORCE: NCOUNT, NSPAN : ', NCOUNT, NSPAN 
      !IF ( NCOUNT > NSPAN ) THEN
      !    print*, ' NCOUNT > NSPAN ' 
      !    RETURN
      !ENDIF 
      ! Actually, still need NCOUNT so that we know when to make the 
      ! aod files (dkh, 03/27/11). 
      NCOUNT = NCOUNT + 1 

 
      !Initialization 
      IF ( FIRST ) THEN 
    
         ! Allocate arrays
         CALL INIT_LIDORT

         CALL INIT_MIE_MOD 

         CALL lidort_inputs_basic_master                               
     &     (  '3p5T_LIDORT_ReadInput.cfg',                            ! INPUT
     &     DO_SOLAR_SOURCES,                                          ! output
     &     DO_UPWELLING,            DO_DNWELLING,                     ! output
     &     DO_FULLRAD_MODE,         DO_USER_STREAMS,                  ! output
     &     DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,               ! output
     &     DO_SSFULL,               DO_SSCORR_TRUNCATION,             ! output
     &     DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,                  ! output
     &     DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                    ! output
     &     DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,           ! output 
     &     DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,               ! output
     &     DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,                ! output
     &     DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,               ! output
     &     DO_ALL_FOURIER,          DO_NO_AZIMUTH,                    ! output
     &     NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,            ! output
     &     NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,     ! output
     &     FLUX_FACTOR, LIDORT_ACCURACY,                              ! output
     &     EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,      ! output
     &     BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_LEVELS,   ! output
     &     STATUS_INPUTREAD, NREADMESSAGES, READMESSAGES, READACTIONS )

         !  Exception handling
         IF ( status_inputread .ne. lidort_success ) THEN
            open(1,file = '3p5T_LIDORT_ReadInput.log', 
     &         status = 'unknown')
            WRITE(1,*)
     &         ' FATAL:   Wrong input from LIDORT input file-read'
            WRITE(1,*)
     &         '  ------ Here are the messages and actions '
            write(1,'(A,I3)')
     &         '    ** Number of messages = ',NREADMESSAGES
            DO N = 1, NREADMESSAGES
               NF = LEN_STRING(READMESSAGES(N))
               NA = LEN_STRING(READACTIONS(N))
               write(1,'(A,I3,A,A)')
     &            'Message # ',N,' : ',READMESSAGES(N)(1:NF)
               write(1,'(A,I3,A,A)')
     &            'Action  # ',N,' : ',READACTIONS(N)(1:NA)
            ENDDO
            close(1)
            STOP
     &         'Read-input fail: Look at file 3p5T_LIDORT_ReadInput.log'
         ENDIF

         FIRST = .FALSE. 

      ENDIF 
     
      IF ( NEED_READ_LER ) THEN

         ! Read TEMIS LER (reflectivities)
         CALL READ_LER_FILE( GET_MONTH() )

         NEED_READ_LER = .FALSE.
      ENDIF
      IF ( ITS_A_NEW_MONTH() ) THEN
      
          NEED_READ_LER = .TRUE.
      
      ENDIF

!      ! Clear
      LOCAL_COST = 0d0  

      CALL SYSTEM_CLOCK(COUNT_RATE=CLOCK_RATE)
      CALL SYSTEM_CLOCK(COUNT=CLOCK_1)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AOD )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      !DO J = JFD, JFD
      !DO I = IFD, IFD
       
         ! Calculate aerosol optical depth in current column
         ! as well as layer-by-layer optical fields (AOD, SSA, PHM)
         CALL CALC_AOD( I, J, AOD )

         ! Safety
         IF ( IT_IS_NAN( SUM(AOD(:)) ) ) THEN
              WRITE(6,*) ' ERROR: AOD is NAN in ', I, J
         ENDIF
 
         ! Calculate contribution to cost function. Scale 
         ! by SCALEF / EARTHSA later. 
         LOCAL_COST(I,J,:) = AOD(:) * GET_AREA_M2( J ) 

         ! Save AOD and contribution to global average AOD
         AOD_SAVE(I,J,:) = AOD_SAVE(I,J,:) + AOD(:) * INV_NSPAN
         AOD_AVE(I,J,:)  = AOD_AVE(I,J,:) 
     &                   + LOCAL_COST(I,J,:) * INV_NSPAN

         ! Diagnostic: average sulfate mass scattering efficiency    
         MASS_ND(I,J) = MASS_ND(I,J) + INV_NSPAN * (
     &                SUM( CHK_STT(I,J,1:LLTROP,IDTSO4) ) )
         
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(COUNT=CLOCK_2)
      WRITE(6,*) ' AOD TIME : ', (CLOCK_2 - CLOCK_1)/
     &                              REAL(CLOCK_RATE)


      CALL SYSTEM_CLOCK(COUNT_RATE=CLOCK_RATE)
      CALL SYSTEM_CLOCK(COUNT=CLOCK_1)

      ! Loop over longitude

!$OMP PARALLEL DO
!$OMP+PRIVATE( I )
!$OMP+PRIVATE( nlayers, nmoments_input )
!$OMP+PRIVATE( DO_PROFILE_LINEARIZATION,DO_SURFACE_LINEARIZATION )
!$OMP+PRIVATE( DO_BRDF_SURFACE, layer_vary_number, layer_vary_flag )
!$OMP+PRIVATE( N_SURFACE_WFS )
      DO I = 1, IIPAR
!      DO I = IFD, IFD

         !Initialize 
         GLOB_FLUX_W(:,I,:)         = 0D0
         JAC_GLOB_FLUX_W(:,:,:,I,:) = 0D0

         ! Call LIDORT routines to calculate radiative fluxes
         ! This routine loops over latitude (in parallel)
         ! and wavelength. 
         CALL LIDORT_DRIVER( I )
 
      ENDDO
!$OMP END PARALLEL DO

      ! Calculate integrated upward radiative flux
      CALL CALC_RADIATIVE_FLUX( GLOB_FLUX_INST )

      ! Calculate integrated upward radiative flux
      CALL CALC_RADIATIVE_FORCE( GLOB_FLUX_INST, COST_FUNC ) 

      ! dkh debug -- 500 nm diagnostics
      print*, 'SUM flux up = ', SUM(GLOB_FLUX_W(NWL_500,:,:))
      WRITE(6,'(A15,I4,2E14.7)') ('SUM flux up = ', L, 
     &  SUM(JAC_GLOB_FLUX_W(2,L,NWL_500,:,:)), 
     &  MAXVAL(JAC_GLOB_FLUX_W(2,L,NWL_500,:,:)), L=1,30)

      ! dkh debug -- 500 nm diagnostics 
      IF ( .not. L_PRINT ) THEN 
         WRITE(*,'(/a/)')'Integrate fluxes at above levels--->'
         WRITE(*,'(3x,1p4e13.5)')(GLOB_FLUX_W(NWL_500,61,32))
         WRITE(*,'(3x,1p4e13.5)')(GLOB_FLUX_W(NWL_500,IFD,JFD))

         WRITE(*,'(/a//6x,a)')
     &   'Upwelling Jacobians w.r.t. layer aerosol extinction TAU -->',
     &   '   ---------- Analytic Jacobian calculations ------------'
         DO L = 1, LLPAR
            WRITE(*,'(i3,2(1x,1p4e14.5))') L, 
!     &      JAC_GLOB_FLUX_W(2,L,NWL_500,61,32)
     &      JAC_GLOB_FLUX_W(2,L,NWL_500,IFD,JFD)
         ENDDO
         WRITE(*,'(/a//6x,a)')
     &   'Upwelling Jacobians w.r.t. layer aerosol single scattering>',
     &   '   ---------- Analytic Jacobian calculations ------------'
         DO L = 1, LLPAR
            WRITE(*,'(i3,2(1x,1p4e14.5))') L, 
!!     &      JAC_GLOB_FLUX_W(2,L,NWL_500,61,32)
     &      JAC_GLOB_FLUX_W(3,L,NWL_500,IFD,JFD)
         ENDDO
         WRITE(*,'(/a//6x,a)')
     &   'Upwelling Jacobians w.r.t. layer aerosol PHM  ------------>',
     &   '   ---------- Analytic Jacobian calculations ------------'
         DO L = 1, LLPAR
            WRITE(*,'(i3,2(1x,1p4e14.5))') L, 
!     &      JAC_GLOB_FLUX_W(2,L,NWL_500,61,32)
     &      JAC_GLOB_FLUX_W(4,L,NWL_500,IFD,JFD)
         ENDDO
         !---------------------------------------------

      ENDIF 

      CALL SYSTEM_CLOCK(COUNT=CLOCK_2)
      WRITE(6,*) ' LIDORT  TIME : ', (CLOCK_2 - CLOCK_1)/
     &                                    REAL(CLOCK_RATE)


      ! dkh debug 
      print*, MAXLOC(AOD_AVE(:,:,:))

      ! Diagnostic printout, see eqn 4 of Martin et al., 2004 ACP
      print*, 'CURRENT AOD STATS: ', 'AOD_ave, Beta_ave, ', 
     &                   ' M_ave (mg SO4 /m2), SCAL'
      print*, 'CURRENT AOD STATS: ', 
     &   SUM(AOD_AVE(:,:,NWL_500))/EARTHSA, 
     &   SUM(AOD_AVE(:,:,NWL_500)) / ( SUM(MASS_ND(:,:)) * 1d03 ),
     &   SUM(MASS_ND(:,:)) * 1d06 / EARTHSA

      print*, ' running BC burden [Tg] = ', 
     &   (SUM(STT(:,:,:,IDTBCPI)) + SUM(STT(:,:,:,IDTBCPO))) * 1d-9

	 ! fha 6-9-2011 added OC total to debug output
      print*, ' running OC burden [Tg] = ', 
     &   (SUM(STT(:,:,:,IDTOCPI)) + SUM(STT(:,:,:,IDTOCPO))) * 1d-9
	 
      ! dkh debug
      print*, ' MAKE_AOD debug ', GET_NYMD(), GET_NYMDb(), GET_NHMS(), 
     &         GET_NHMSb()

      ! Write out AOD file on the last time through 
      IF ( GET_NYMD() == GET_NYMDb() .and. 
     &     GET_NHMS() == GET_NHMSb() .or. 
     &     NCOUNT     == NSPAN             ) THEN 
         CALL MAKE_AOD_FILE( GET_NYMD(), GET_NHMS(),  N_CALC ) 
      ENDIF 
    
      ! fur debugin just forward
      !print*, ' SKIPPING RF ADJOINT CALCULATION *** '
      !print*, ' SKIPPING RF ADJOINT CALCULATION *** '
      !print*, ' SKIPPING RF ADJOINT CALCULATION *** '
      !GOTO 131

      !--------------------------------------------------------------- 
      ! Adjoint code begins here
      !--------------------------------------------------------------- 

      ! fwd code:
      !CALL CALC_RADIATIVE_FLUX( GLOB_FLUX_INST )
      !CALL CALC_RADIATIVE_FORCE( GLOB_FLUX_INST, COST_FUNC ) 
      ! adj code:
      ! note: here GLOB_FLUX_INST_ADJ is output
      CALL ADJ_CALC_RADIATIVE_FORCE( GLOB_FLUX_INST_ADJ )
      ! note: here GLOB_FLUX_INST_ADJ is input
      CALL ADJ_CALC_RADIATIVE_FLUX( GLOB_FLUX_INST_ADJ )

 121  CONTINUE 

      ! dkh debug 
      print*, 'ddd loc a ', STT_ADJ(IFD,JFD,LFD,NFD)

      ! adjoint of LIDORT routine   
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR 
      !DO J = JFD, JFD
      !DO I = IFD, IFD


         ! Get 1D array index
         IJLOOP = ( J - 1 ) * IIPAR  + I

         ! Only bother to do this if the sun is shining 
         ! ISCCP definition of sunset/sunrise.  
         ! Rossow and Duenas, 2004
         ! CYCLE statement leads to OpenMP issues (dkh, 03/27/11) 
         !IF ( SUNCOS(IJLOOP) < 0.0005d0 ) CYCLE 
         IF ( SUNCOS(IJLOOP) >= 0.0005d0 ) THEN
  
            ! fwd code:
            !CALL LIDORT_DRIVER( I, J )
            ! adj code:
            CALL ADJ_LIDORT_DRIVER( I, J )
    
            ! Reset just to be safe
            ! fwd code:
            !!Initialize 
            !GLOB_FLUX_W(:,I,J)         = 0D0
            !JAC_GLOB_FLUX_W(:,:,:,I,J) = 0D0
            ! adj code:
            GLOB_FLUX_W_ADJ(:,I,J)     = 0D0
            JAC_GLOB_FLUX_W(:,:,:,I,J) = 0D0

         ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! dkh debug 
      print*, 'ddd loc b ', STT_ADJ(IFD,JFD,LFD,NFD)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      !DO J = JFD,JFD
      !DO I = IFD,IFD

         ! Get 1D array index
         IJLOOP = ( J - 1 ) * IIPAR  + I
      
         ! Only bother to do this if the sun is shining 
         ! ISCCP definition of sunset/sunrise.  
         ! Rossow and Duenas, 2004
         ! CYCLE leads to OpenMP issues (dkh, 03/27/11) 
         !IF ( SUNCOS(IJLOOP) < 0.0005d0 ) CYCLE
         IF ( SUNCOS(IJLOOP) >= 0.0005d0 ) THEN 
 
            ! Calculate aerosol optical depth in current column
            CALL ADJ_CALC_AOD( I, J )

         ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! dkh debug 
      !print*, 'ddd loc c ', STT_ADJ(IFD,JFD,LFD,NFD)
      !print*, 'ddd loc c HACK zero STT_ADJ'
      !print*, 'ddd loc c HACK zero STT_ADJ'
      !print*, 'ddd loc c HACK zero STT_ADJ'
      !STT_ADJ = 0d0 
      !STT_ADJ(IFD,JFD,LFD,NFD) = 1d0

 131  CONTINUE 

      ! Return to calling program
      END SUBROUTINE CALC_RF_FORCE

!------------------------------------------------------------------------------

      SUBROUTINE CALC_AOD( I, J, AOD )
!
!******************************************************************************
!  Subroutine CALC_AOD calculates the aerosol optical depth in column (AOD)
!  as well as the layer-by-layer value needed for LIDORT input 
!  (dkh, 01/21/08, 07/29/10) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I, J (Integer) : X-Y Grid locations
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) AOD  (REAL*8)  : Aerosol optical depth            [none]
!     
!  Module variables as Output:
!  ============================================================================
!  (1 ) LAYER_AOD  (REAL*8) : Aerosol optical thickness of each layer 
!  (2 ) LAYER_SSA  (REAL*8) : Aerosol single scat albedo of each layer
!  (3 ) LAYER_PHM  (REAL*8) : Aerosol phase momemnt of each layer 
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE CHECKPT_MOD,  ONLY : RP_OUT,  CHK_STT
      USE DAO_MOD,      ONLY : RH,      BXHEIGHT, AIRVOL
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACERID_MOD, ONLY : IDTSO4,  IDTNIT,   IDTNH4
      USE TRACERID_MOD, ONLY : IDTBCPO, IDTBCPI
      USE TRACERID_MOD, ONLY : IDTOCPO, IDTOCPI
	  ! ^fha 6-6-2011: added OC
	  

      ! dkh debug 
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD

#     include "CMN_SIZE"
#     include "CMN"

      ! Arguments
      INTEGER, INTENT(IN)   :: I, J
      REAL*8,  INTENT(OUT)  :: AOD(NWL_MAX)
    
      ! Local variables 
      INTEGER               :: L,    W
      INTEGER               :: NSP
      REAL*8                :: DIAM(NSP_MAX)
      REAL*8                :: WET_DIAM(NSP_MAX)
      REAL*8                :: MASS(NSP_MAX)
      REAL*8                :: NCONC(NSP_MAX)
      REAL*8                :: BEXT, SSA, PHM(NPM_MAX)
      REAL*8                :: RHL
      REAL*8                :: BEXT_TOT
      REAL*8                :: SCAT_TOT
      REAL*8                :: PHM_TOT(NPM_MAX)
 

 
      !=================================================================
      ! CALC_AOD begins here!
      !=================================================================

      ! note: already in a parallel loop over I,J
      !DO L = 1 , LLTROP
      DO L = 1 , LLPAR

         ! Get local RH
         RHL = MIN(99.D0,RH(I,J,L))


         ! Loop over aerosol type
         DO NSP = 1, NSP_MAX 
	  !fha  (6/1/11) added OC
	        IF ( NSP == NSP_BC ) THEN 

               ! Aerosol wet size (um)
               WET_DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, 
     &                                        DRY_DIAM(NSP)/2d0 ) 

               ! Aerosol mass [kg/box]
               MASS(NSP)  = REAL(CHK_STT(I,J,L,IDTBCPI),8)
     &                    + REAL(CHK_STT(I,J,L,IDTBCPO),8)

               ! Total effective diam is volume weighted average
               ! of wet (BCPI) and dry (BCPO) components.  Since
               ! BC assumed to have density of 1, then volume ~ mass
               DIAM(NSP)  = 
     &          (   REAL(CHK_STT(I,J,L,IDTBCPI),8) * WET_DIAM(NSP)
     &            + REAL(CHK_STT(I,J,L,IDTBCPO),8) * DRY_DIAM(NSP) )
     &            / MASS(NSP)

               ! dkh debug
               !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
               !   print*, 'ddd loc md ',  MASS(NSP), DIAM(NSP)
               !ENDIF 

               ! Aerosol number concentration.
               NCONC(NSP) = GET_NCONC( 
     &                      MASS(NSP), 
     &                      0d0, 
!------------------------------------------------------------------
! BUG FIX: since we don't have mass of water, estimate NCONC based 
!  on dry radius (dkh, 01/27/11) 
!     &                      DENSE(NSP), 1.d0,   DIAM(NSP), 
     &                      DENSE(NSP), 1.d0,   DRY_DIAM(NSP), 
!------------------------------------------------------------------
     &                      SIGMA(NSP), AIRVOL(I,J,L)      )

            ELSEIF ( NSP == NSP_OC ) THEN 

               ! Aerosol wet size (um)
               WET_DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, 
     &                                        DRY_DIAM(NSP)/2d0 ) 

               ! Aerosol mass [kg/box]
               MASS(NSP)  = REAL(CHK_STT(I,J,L,IDTOCPI),8)
     &                    + REAL(CHK_STT(I,J,L,IDTOCPO),8)

               ! Total effective diam is volume weighted average
               ! of wet (BCPI) and dry (BCPO) components.  Since
               ! BC assumed to have density of 1, then volume ~ mass
               DIAM(NSP)  = 
     &          (   REAL(CHK_STT(I,J,L,IDTOCPI),8) * WET_DIAM(NSP)
     &            + REAL(CHK_STT(I,J,L,IDTOCPO),8) * DRY_DIAM(NSP) )
     &            / MASS(NSP)

               ! dkh debug
               !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
               !   print*, 'ddd loc md ',  MASS(NSP), DIAM(NSP)
               !ENDIF 

               ! Aerosol number concentration.
               NCONC(NSP) = GET_NCONC( 
     &                      MASS(NSP), 
     &                      0d0, 
!------------------------------------------------------------------
! BUG FIX: since we don't have mass of water, estimate NCONC based 
!  on dry radius (dkh, 01/27/11) 
!     &                      DENSE(NSP), 1.d0,   DIAM(NSP), 
     &                      DENSE(NSP), 1.d0,   DRY_DIAM(NSP), 
!------------------------------------------------------------------
     &                      SIGMA(NSP), AIRVOL(I,J,L)      )

            ELSE

               ! try it this way:
               ! specify dry diam
               !DIAM = 0.1 
               !DIAM(NSP)  = DRY_DIAM(NSP)

                ! Aerosol mass [kg/box]
               MASS(NSP) = REAL(CHK_STT(I,J,L,IDTSO4),8)
     &                   + REAL(CHK_STT(I,J,L,IDTNIT),8)
     &                   + REAL(CHK_STT(I,J,L,IDTNH4),8)

               ! dkh debug
               !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
               !   print*, 'ichk: fwd GET_NCONC sulf ', MASS(NSP), 
               !            DRY_DIAM(NSP)
               !ENDIF 

               ! get NCONC based on dry mass
               NCONC(NSP) = GET_NCONC( 
     &                      MASS(NSP),                     
     &                      0d0,        
     &                      DENSE(NSP), 1.d0, DRY_DIAM(NSP), 
     &                      SIGMA(NSP), AIRVOL(I,J,L)     )

               ! dkh debug
               !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
               !   print*, 'ichk: fwd GET_WETD2 ', MASS(NSP), 
               !            REAL(RP_OUT(I,J,L,4),8), NCONC(NSP)
               !ENDIF 

               
!!               ! get wet diam 
!!               DIAM(NSP) = GET_WETD2(
!!     &                     MASS(NSP),
!!!     &                     REAL(RP_OUT(I,J,L,4),8), 
!!! HACK: for testing SO4
!!     &                     1000d0,
!!     &                     DENSE(NSP), 1.d0, NCONC(NSP),
!!     &                     SIGMA(NSP), AIRVOL(I,J,L)     )

               ! Aerosol wet size (um)
               DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, 
     &                                    DRY_DIAM(NSP)/2d0 ) 

               ! if you feed the above routine zero water, will it
               ! give back the dry diam? yes. 
               !IF ( L == LFD ) THEN
               !  print*, ' DIAM, DRY_DIAM = ', DIAM(NSP), DRY_DIAM(NSP)
               !ENDIF 

            ENDIF    
         ENDDO

!         ! Find total volume (need to convert RP_OUT to kg/box)
!         V_TOT = MASS(NSP_SULF) / DENSE(NSP_SULF)
!     &         + MASS(NSP_BC)   / DENSE(NSP_BC)
!!     &         + REAL(RP_OUT(I,J,L,4),8) 
!     &         + REAL(RP_OUT(I,J,L,4),8) * 1d-9 * AIRVOL(I,J,L)
! 
!         ! Initialize to zero before summing over species
!         LAYER_AOD(I,J,L,:)   = 0d0
!         LAYER_SSA(I,J,L,:)   = 0d0 
!         LAYER_PHM(I,J,L,:,:) = 0d0
!
!!         DO NSP = 1, NSP_MAX   
!         !DO NSP = 2, 2
!         !DO NSP = 1, 1
!
!            ! We are going to combine optical properties from different 
!            ! aerosoly types according to a volume-weighting factor
!            VW = MASS(NSP) / ( DENSE(NSP) * V_TOT )
!              
!            ! For sulfate, add water contribution
!            IF ( NSP == NSP_SULF ) THEN
!               !VW = VW + REAL(RP_OUT(I,J,L,4),8) / ( 1D0 * V_TOT ) 
!               VW = VW + REAL(RP_OUT(I,J,L,4),8) * 1d-9 
!     &            * AIRVOL(I,J,L) / ( 1D0 * V_TOT ) 
!            ENDIF 
!
!            ! fur debugin
!            VW = 1d0
!
!
!               ! dkh debug
!               IF ( I == 1 .and. J == 1 ) THEN 
!                 print*, ' I, J, L, NSP, NCONC, DIAM, MASS, DENSE, VW, 
!     &           RH, dz ',  
!     &           I, J, L, NSP, NCONC(NSP), DIAM(NSP), MASS(NSP),
!     &           DENSE(NSP), VW,
!     &           RHL, 
!     &           BXHEIGHT(I,J,L) 
!               ENDIF 

         !IF ( L == LFD ) THEN
         !   print*, ' HACK ddd loc nd NCONC =', NCONC
         !   print*, ' HACK ddd loc nd DIAM = ', DIAM
         !    DEBUG_J(I,J) = NCONC(2)
         !ENDIF 

         ! Loop over wavelengths
         DO W = 1, NWL_MAX

            BEXT_TOT   = 0d0
            SCAT_TOT   = 0D0 
            PHM_TOT(:) = 0D0

            DO NSP = 1, NSP_MAX

               ! Lookup aerosol optical properties
               CALL MIE_LOOKUP( NCONC(NSP), DIAM(NSP), W,  NSP,    
     &                          BEXT,       SSA,       PHM, I,J,L   )

               !print*, ' mie done ', I, J, L, NSP, W

               !IF ( L == LFD .and. W == 3 .and. NSP == 2 ) THEN 
               !   print*, 'HACK ddd loc bps BEXT = ', BEXT
               !   print*, 'HACK ddd loc bps SSA  = ', SSA
               !   print*, 'HACK ddd loc bps PHM  = ', PHM
               !    DEBUG_J(I,J) = BEXT 
               !ENDIF 

               BEXT_TOT   = BEXT_TOT   + BEXT
               SCAT_TOT   = SCAT_TOT   + BEXT * SSA 
               PHM_TOT(:) = PHM_TOT(:) + PHM(:) * BEXT * SSA 


            ENDDO 

            IF ( BEXT_TOT  < SMALL .or. 
     &           SCAT_TOT  < SMALL      ) THEN
               CALL ERROR_STOP(' Trouble calculating combined opt prop', 
     &                         ' aero_optical_mod.f' )
            ENDIF 

 
            ! Calculate AOD for this layer. BXHEIGHT is [m] while
            ! BEXT is in [1/m]
            LAYER_AOD(I,J,L,W)   = BXHEIGHT(I,J,L) * BEXT_TOT 
    
            ! Layer single scatering albedo
            LAYER_SSA(I,J,L,W)   = SCAT_TOT / BEXT_TOT !/ NCONC_TOT

            ! Layer phase moments
            LAYER_PHM(I,J,L,:,W) = PHM_TOT(:) / SCAT_TOT 
 
            ! dkh debug
            !IF ( I == IFD .and. J == JFD .and. L == LFD 
            !              .and. W == 3                  ) THEN
            !   print*, 'ddd loc w ', LAYER_AOD(I,J,L,W)
            !   print*, 'ddd loc s ', LAYER_SSA(I,J,L,W)
            !
            !ENDIF 

 
         ENDDO ! LAMBDA
       
      ENDDO ! L


      ! Total AOD over all vertical levels
      DO W = 1, NWL_MAX
         AOD(W) = SUM( LAYER_AOD(I,J,:,W) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE CALC_AOD

!------------------------------------------------------------------------------

      FUNCTION GET_WET_DIAM( RH , NSP, DRY_RADIUS ) RESULT( DIAM )
!
!******************************************************************************
!  Subroutine GET_WET_DIAM calculates mean wet radius of aerosol. 
!  Currently based on jv_spec.dat
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) RH         (REAL*8) : Relative humidity                    [%]
!  (2 ) NSP       (INTEGER) : Species index 
!  (3 ) DRY_RADIUS (REAL*8) : Dry mode radius                     [um]
!     
!  Output:
!  ============================================================================
!  (1 ) DIAM       (REAL*8) : Diameter                            [um]
!     
!  NOTES:
!  (1 ) Updated growth factors to GCv8-03-1 (dkh, 09/27/10) 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP

      ! Arguments
      REAL*8, INTENT(IN)  :: RH 
      INTEGER,INTENT(IN)  :: NSP
      REAL*8, INTENT(IN)  :: DRY_RADIUS

      ! Value
      REAL*8              :: DIAM
  
      ! Local variables 
      REAL*8              :: RADIUS
      REAL*8              :: GF

      !=================================================================
      ! GET_WET_DIAM begins here!
      !=================================================================

      SELECT CASE( NSP )
      
         ! Sulfate
         CASE( 1 ) 

            ! GF's udated to GCv8-03-01 (dkh, 09/27/10) 
            IF     ( RH <= 50.d0 ) THEN
               
               ! GF = 1.d0 + ( 1.36 - 1.0 ) / ( 50 - 0 ) * ( RH - 0 )
               GF = 1.d0 + 0.0073d0 * RH

            ELSEIF ( RH <= 70.d0 ) THEN
               
               ! GF = 1.36 + ( 1.52 - 1.36 ) / ( 70 - 50 ) * ( RH - 50 )
               GF = 0.99d0 + 0.0075d0 * RH

            ELSEIF ( RH <= 80.d0 ) THEN
               
               ! GF = 1.52 + ( 1.64 - 1.52 ) / ( 80 - 70 ) * ( RH - 70 )
               GF = 0.68d0 + 0.012d0 * RH

            ELSEIF ( RH <= 90.d0 ) THEN
               
               ! GF = 1.64 + ( 1.87 - 1.64 ) / ( 90 - 80 ) * ( RH - 80 )
               GF = -0.23d0 + 0.023d0 * RH

            ELSEIF ( RH <= 95.d0 ) THEN
               
               ! GF = 1.87 + ( 2.18 - 1.87 ) / ( 95 - 90 ) * ( RH - 90 )
               GF = -3.8d0 + 0.063d0 * RH

            ELSEIF ( RH <= 100.d0 ) THEN
               
               ! GF = 2.18 + ( 3.13 - 2.18 ) / ( 99 - 95 ) * ( RH - 95 )
               GF = -20d0 + 0.24d0 * RH

            ENDIF 

         ! BC
         CASE( 2 )  

            ! GF's udated to GCv8-03-01 (dkh, 09/27/10) 
            IF     ( RH < 70.d0 ) THEN
               
               GF = 1.d0

            ELSEIF ( RH <= 80.d0 ) THEN

               ! GF = 1.d0 + ( 1.2 - 1.0 ) / ( 80 - 70 ) * ( RH - 70 )
               GF = -0.4d0 + 0.02d0 * RH 

            ELSEIF ( RH <= 90.d0 ) THEN

               ! GF = 1.2d0 + ( 1.4 - 1.2 ) / ( 90 - 80 ) * ( RH - 80 )
               GF = -0.4d0 + 0.02d0 * RH 
           
            ELSEIF ( RH <= 95.d0 ) THEN

               ! GF = 1.4d0 + ( 1.5 - 1.4 ) / ( 95 - 90 ) * ( RH - 90 )
               GF = -0.14d0 + 0.017d0 * RH 
           
            ELSEIF ( RH <= 100.d0 ) THEN
 
               ! GF = 1.5d0 + ( 1.9 - 1.5 ) / ( 99 - 95 ) * ( RH - 95 )
               GF = -8.0D0 + 0.1d0 * RH
          
            ENDIF 
			
         ! OC
         CASE( 3 )  

            ! added OC GF's from GCv8-03-01  (fha, 6-1-11)
            ! RH          r-eff         r-eff v r-eff @ RH==0 (GF)
            ! 00          .07           1
            ! 50          .087          1.24
            ! 70          .095          1.36 
            ! 80          .102          1.46
            ! 90          .116          1.66
            ! 95          .133          1.9
            ! 99          .177          2.53
            IF     ( RH <= 50.d0 ) THEN
               
               ! GF = 1.d0 + ( 1.24 - 1.0 ) / ( 50 - 0 ) * ( RH - 0 )
               GF = 1.d0 + 0.0049d0 * RH

            ELSEIF ( RH <= 70.d0 ) THEN
               
               ! GF = 1.24d0 + ( 1.36 - 1.24 ) / ( 70 - 50 ) * ( RH - 0 )
               GF = 0.9571d0 + 0.0057d0 * RH

            ELSEIF ( RH <= 80.d0 ) THEN

               ! GF = 1.36d0 + ( 1.46 - 1.36 ) / ( 80 - 70 ) * ( RH - 70 )
               GF = 0.6571d0 + 0.01d0 * RH 

            ELSEIF ( RH <= 90.d0 ) THEN

               ! GF = 1.46d0 + ( 1.66 - 1.46 ) / ( 90 - 80 ) * ( RH - 80 )
               GF = -0.1429d0 + 0.02d0 * RH 
           
            ELSEIF ( RH <= 95.d0 ) THEN

               ! GF = 1.66d0 + ( 1.9 - 1.66 ) / ( 95 - 90 ) * ( RH - 90 )
               GF = -2.714d0 + 0.04857d0 * RH 
           
            ELSEIF ( RH <= 100.d0 ) THEN
 
               ! GF = 1.9d0 + ( 2.53 - 1.9 ) / ( 99 - 95 ) * ( RH - 95 )
               GF = -13.029D0 + 0.1571d0 * RH
          
            ENDIF 
      
         CASE DEFAULT 
 
            CALL ERROR_STOP( 'bad NSP in GET_WET_DIAM', 
     &                       'aero_optical_mod.f')

      END SELECT 

      RADIUS = DRY_RADIUS * GF

      DIAM = 2D0 * RADIUS 

      ! dkh debug
      !print*, 'RH, DIAM ', RH, DIAM

      ! Return to calling program
      END FUNCTION GET_WET_DIAM

!------------------------------------------------------------------------------

      FUNCTION GET_NCONC( MASS_DRY_kg, MASS_WET, DENSE_DRY, DENSE_WET, 
     &                      D,        SIGMA,  BVOL                  )
     &   RESULT( NCONC )
!
!******************************************************************************
!  Subroutine GET_NCONC calculates the number concentration [#/cm3] that
!  corresponds to the lognormally distributed (D,SIGMA) aerosol mass.
!  (dkh, 01/21/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MASS_DRY_kg  (REAL*8) : Aerosol dry mass concentration             [kg/box]
!  (2 ) MASS_WET  (REAL*8) : Aerosol water concentration                [ug/m3]
!  (3 ) DENSE_DRY (REAL*8) : Dry aerosol density                        [g/cm3]
!  (4 ) DENSE_WET (REAL*8) : Wet aerosol density                        [g/cm3]
!  (5 ) D         (REAL*8) : Assumed aerosol median diameter            [um]
!  (6 ) SIGMA     (REAL*8) : Assumed aerosol standard deviation         [um]
!  (7 ) BVOL      (REAL*8) : Box volume                                 [m3]
!     
!  Output:
!  ============================================================================
!  (1 ) NCONC     (REAL*8) : Aerosol number concentration               [#/cm3]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP

      ! Arguments
      REAL*8, INTENT(IN)     :: MASS_WET
      REAL*8, INTENT(IN)     :: MASS_DRY_kg
      REAL*8, INTENT(IN)     :: DENSE_DRY
      REAL*8, INTENT(IN)     :: DENSE_WET
      REAL*8, INTENT(IN)     :: D
      REAL*8, INTENT(IN)     :: SIGMA
      REAL*8, INTENT(IN)     :: BVOL

      ! Fuction value
      REAL*8                 :: NCONC
        
      ! Local variables 
      REAL*8                 :: DVM
      REAL*8                 :: DENSITY 
      REAL*8                 :: VOLUME
      REAL*8                 :: MASS_TOTAL
      REAL*8                 :: MASS_DRY
 
      ! Parameters
      REAL*8, PARAMETER      :: PI = 3.14159265358979324D0

      !=================================================================
      ! GET_NCONC begins here!
      !=================================================================

      ! Convert mass units ( [kg/box] --> [ug/m3] ) 
      MASS_DRY   = MASS_DRY_kg * 1d9 / BVOL

      ! Get total mass
      MASS_TOTAL = MASS_DRY + MASS_WET
   
      IF ( MASS_TOTAL < 1d-20 ) THEN 
         print*, 'housten, we have a problem' 
         PRINT *, 'Aerosol mass too small: 1'
         NCONC = 0d0  ! <-- this will cause LIDORT to give NANs
!         CALL ERROR_STOP( 'Aersol mass too small: 1', 
!     &                    'aero_optical_mod.f')
         return
      ENDIF 

      ! Average density
      ! Also convert [g/cm3] --> [ug/cm3]
      DENSITY = ( DENSE_DRY * MASS_DRY + DENSE_WET * MASS_WET ) 
     &        / ( MASS_TOTAL ) * 1d6

      ! Total volume [cm3/m3 air]
      VOLUME  = MASS_TOTAL / DENSITY 

      ! Volume mean diameter can be calculated from the log-normal
      ! median diameter and the geometric standard deviation,
      ! see exercise 7.7 of S&P, 1998. 
      ! Also convert units [um] --> [cm]
      DVM = D * EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) * 1D-04

      ! The total number concentration can be calculated from 
      ! the total volume and the volume mean volume  diameter, 
      ! see table 7.2 of S&P, 1998. 
      ! Also convert units [#/m3 air] --> [#/cm3 air]
      NCONC = 6d0 * VOLUME / ( PI * DVM ** 3 ) * 1d-06
 
      ! Return to calling program
      END FUNCTION GET_NCONC

!------------------------------------------------------------------------------
      FUNCTION GET_WETD2( MASS_DRY_kg, MASS_WET, DENSE_DRY, DENSE_WET, 
     &                      NCONC,  SIGMA,    BVOL                 )
     &   RESULT( D )
!
!******************************************************************************
!  Subroutine GET_WETD2 calculates the log-normal median diameter [um] that
!  corresponds to the lognormally distributed (D,SIGMA) wet aerosol mass 
!  with the number concentration calculated from the calculated dry particle 
!  mass and assumed size.
!  (dkh, 02/22/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MASS_DRY_kg (REAL*8) : Aerosol dry mass concentration             [kg/box]
!  (2 ) MASS_WET  (REAL*8) : Aerosol water concentration                [ug/m3]
!  (3 ) DENSE_DRY (REAL*8) : Dry aerosol density                        [g/cm3]
!  (4 ) DENSE_WET (REAL*8) : Wet aerosol density                        [g/cm3]
!  (5 ) NCONC     (REAL*8) : Number concentration                       [#/cm3]
!  (6 ) SIGMA     (REAL*8) : Assumed aerosol standard deviation         [um]
!  (7 ) BVOL      (REAL*8) : Box volume                                 [m3]
!     
!  Output:
!  ============================================================================
!  (1 ) D         (REAL*8) : Aerosol log-normal median diameter         [um]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP

      ! Arguments
      REAL*8, INTENT(IN)     :: MASS_WET
      REAL*8, INTENT(IN)     :: MASS_DRY_kg
      REAL*8, INTENT(IN)     :: DENSE_DRY
      REAL*8, INTENT(IN)     :: DENSE_WET
      REAL*8, INTENT(IN)     :: NCONC
      REAL*8, INTENT(IN)     :: SIGMA
      REAL*8, INTENT(IN)     :: BVOL

      ! Fuction value
      REAL*8                 :: D
        
      ! Local variables 
      REAL*8                 :: DVM
      REAL*8                 :: DENSITY 
      REAL*8                 :: VOLUME
      REAL*8                 :: MASS_TOTAL
      REAL*8                 :: MASS_DRY
 
      ! Parameters
      REAL*8, PARAMETER      :: PI = 3.14159265358979324D0

      !=================================================================
      ! GET_WETD2 begins here!
      !=================================================================

      ! Convert mass units ( [kg/box] --> [ug/m3] ) 
      MASS_DRY   = MASS_DRY_kg * 1d9 / BVOL

      ! Get total mass
      MASS_TOTAL = MASS_DRY + MASS_WET
   
      IF ( MASS_TOTAL < 1d-20 ) THEN 
         print*, 'housten, we have a problem' 
         !NCONC = 0d0  ! <-- this will cause LIDORT to give NANs
         PRINT *, 'Aerosol mass too small: 2'
         RETURN
!         CALL ERROR_STOP( 'Aersol mass too small: 2', 
!     &                    'aero_optical_mod.f')
      ENDIF 
 
      ! Average density
      ! Also convert [g/cm3] --> [ug/cm3]
      DENSITY = ( DENSE_DRY * MASS_DRY + DENSE_WET * MASS_WET ) 
     &        / ( MASS_TOTAL ) * 1d6

      ! Total volume [cm3/m3 air]
      VOLUME  = MASS_TOTAL / DENSITY 

      ! The total number concentration can be calculated from 
      ! the total volume and the volume mean volume  diameter, 
      ! see table 7.2 of S&P, 1998. 
      ! Also convert units [#/m3 air] --> [#/cm3 air]
      DVM = ( 6d0 * VOLUME / ( PI * NCONC ) * 1d-06 ) ** (1.0/3.0)

      ! Volume mean diameter can be calculated from the log-normal
      ! median diameter and the geometric standard deviation,
      ! see exercise 7.7 of S&P, 1998. 
      ! Also convert units [um] --> [cm]
      D = DVM * 1d04 / ( EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) )

 
      ! Return to calling program
      END FUNCTION GET_WETD2

!------------------------------------------------------------------------------

      SUBROUTINE MIE_LOOKUP( NCONC, DIAM, W, NSP, BEXT, SSA, PHM,I,J,L )
!
!******************************************************************************
!  Subroutine MIE_LOOKUP 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [unit]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE MIE_MOD, ONLY : MIE_TABLE 
      USE MIE_MOD, ONLY : MIE_TABLE_PHM

      ! dkh dbug
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD 

!#     include "foo"  

      ! Arguments
      REAL*8,  INTENT(IN)  :: NCONC
      REAL*8,  INTENT(IN)  :: DIAM
      INTEGER, INTENT(IN)  :: W
      INTEGER, INTENT(IN)  :: NSP
      REAL*8,  INTENT(OUT) :: BEXT
      REAL*8,  INTENT(OUT) :: SSA
      REAL*8,  INTENT(OUT) :: PHM(NPM_MAX)
    
      INTEGER, INTENT(IN)  :: L, I, J
      ! Local variables 
      REAL*8               :: N, BEXTT, RAD
      LOGICAL              :: FAIL
      CHARACTER*90         :: MESSAGES(4)
      REAL*8               :: PHM2(0:NPM_MAX)
      REAL*8               :: BEXT0
      REAL*8               :: SSA0

   
      ! Parameters
      ! BUG FIX: need to calculated this separately for BC and SULF, 
      ! as they have different size ranges and spacing in the lookup
      ! table. (dkh, 09/27/10) 
      !REAL*8, PARAMETER    :: RNG = LOG(RMAX/RMIN) / REAL(NRA_MAX)
      REAL*8               :: RNG 

      !=================================================================
      ! MIE_LOOKUP begins here!
      !=================================================================

      ! BUG FIX: need to calculated this separately for BC and SULF, 
      ! as they have different size ranges and spacing in the lookup
      ! table. (dkh, 09/27/10) 
      RNG = LOG(RMAX(NSP)/RMIN(NSP)) / ( REAL(NRA_MAX) - 1d0 )

      RAD = DIAM / 2D0 

      ! Get index in MIE_TABLE that corresponds to current diameter.
      ! Diameters are calculated at NRA_MAX - 1 points between Dmin and 
      ! Dmax with a spacing defined in ln space: 
      !      d_n = d_min * exp( ( n - 1 ) * ln(d_max/d_min) / npts )
      ! The inverse of this formulat gives us the n for any d_n. 
      N   = LOG(RAD/RMIN(NSP)) / RNG + 1

      ! Enforce bounds on this index
      IF ( N < 1.0           ) N = 1.d0
      IF ( N > REAL(NRA_MAX) ) N = REAL(NRA_MAX)

      BEXTT = MIE_TABLE(W,NSP,NOP_BEXT,INT(N))

      ! Table values calculated for N = 1 [#/cm3], so scale up to current
      ! number concentraiton and convert units to 1/m
      BEXT0 = BEXTT * NCONC * 1d-6
     
      IF ( NSP == NSP_SULF ) THEN
         SSA0  = 1d0  ! for sulfate
      ELSE
         SSA0  = MIE_TABLE(W,NSP,NOP_SSA,INT(N))
      ENDIF 


      ! Get phase moment from table
      PHM(1:NPM_MAX) = MIE_TABLE_PHM(W,NSP,INT(N),:)

      ! enhance absorption owing to mixing by factor of 
      ! ABS_FAC = 1.5 (dkh, 02/01/11) 
      IF ( NSP == NSP_BC ) THEN
         BEXT = BEXT0 * ( SSA0 + ABS_FAC * ( 1d0 - SSA0 ) )
         SSA  = SSA0  / ( SSA0 + ABS_FAC * ( 1d0 - SSA0 ) )
      ELSE 
         BEXT = BEXT0 
         SSA  = SSA0 
      ENDIF

     
! Now we calculate these online (dkh, 07/29/10) .
! BUT it is too slow. so go back to lookup table. 
!      ! Call MIE code  
!      CALL GC_FORWARD_MIE 
!     &  ( NPM_MAX, NCONC, DIAM, SIGMA(NSP), REFRAC_INDEX(NSP,W,:),   ! input 
!     &    LAMBDA_TABLE(W) * 1d-03,                                   ! input
!     &    BEXT, SSA, PHM2, FAIL, MESSAGES )                          ! Output
!
!     
!      ! GC_FORWARD_MIE uses the zero element of PHM; just keep 1:NPM_MAX. 
!      PHM(1:NPM_MAX) = PHM2(1:NPM_MAX)
!
!      ! error checking 
!      IF ( FAIL ) THEN 
!        WRITE(*,*) MESSAGES(1)
!         WRITE(*,*) MESSAGES(2)
!         WRITE(*,*) MESSAGES(3)
!         WRITE(*,*) MESSAGES(4)
!      ENDIF 
   
!       !  debug: write results
!       print*, ' !--- MIE DBG ---! '
!       print*, ' NSP, W ', NSP, W 
!       print*, ' NCONC, DIAM ', NCONC, DIAM 
!       print*, ' RTS: BEXT ', BEXT , BEXT2
!       print*, ' RTS: SSA  ', SSA, SSA2
!       DO L = 1, NPM_MAX
!          print*,  'RTS PHM ', PHM(L), PHM2(L)
!       ENDDO
 
      ! dkh debug
      !IF ( W == 3 .and. L == LFD .and. 
      !   I == IFD .and. J == JFD .and. NSP == 2 ) THEN 
      !   print*, ' fwd: BEXTT = ', BEXTT
      !   print*, ' fwd: BEXT0 = ', BEXT0
      !   print*, ' fwd: NCONC = ', NCONC
      !   print*, ' fwd: SSA0  = ', SSA0
      !   print*, ' fwd: N     = ', N, INT(N)
      !ENDIF 

      ! Return to calling program
      END SUBROUTINE MIE_LOOKUP

!------------------------------------------------------------------------------

      SUBROUTINE LIDORT_DRIVER( I )
!
!******************************************************************************
!  Subroutine LIDORT_DRIVER was a modified version of the jactest_3p3_solar2
!  testing routine.  Now it is based on 3p5T_lps_tester.f90. (dkh, 05/18/10) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER)     : Lon index                            [none]
!  (2 ) J (INTEGER)     : Lat index                            [none]
!     
!  Module variables as Output:
!  ============================================================================
!  (1 ) GLOB_FLUX_W     (REAL*8) : Upward reflectance          [%]
!  (2 ) JAC_GLOB_FLUX_W (REAL*8) : Sensitivity of upward reflectance []
!     
!  NOTES:
!  (1 ) This is based upon the jactest_3p3_solar2.f driver provided with 
!       LIDORT by R. Spurr. Original code is lower case. 
!  (2 ) Update for GCv8.  Based on 3p5T_lps_tester.f90 (dkh, 05/18/10) 
!  (3 ) Additional updates and parallelization (dkh, 07/29/10) 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : JFD 
      USE DAO_MOD,         ONLY : SUNCOS

      !  Include files
#     include 'CMN_SIZE'

       ! Arguments 
      INTEGER, INTENT(IN)     :: I 
 
      !---------------------------------------------------------------
      ! LIDORT variables that need to be defined locally because
      ! they will depend on location
      !---------------------------------------------------------------
  

      !  BOA solar zenith angles (degrees)
      REAL(KIND=8) :: BEAM_SZAS ( MAXBEAMS )

      !  Exception handling 
      !  ------------------

      !  Exception handling for LIDORT Input Checking. New code, 18 May 2010
      !     Message Length should be at least 120 Characters

      INTEGER         :: STATUS_INPUTCHECK
      INTEGER         :: NCHECKMESSAGES
      CHARACTER*120   :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER*120   :: CHECKACTIONS (0:MAX_MESSAGES)

      !  Exception handling for LIDORT Model Calculation. New code, 18 May 2010
      INTEGER         :: STATUS_CALCULATION
      CHARACTER*120   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

      !  4. ATMOSPHERIC INPUTS
      !  =====================

      !  PTH inputs
      !  ----------

      !  multilayer Height inputs
      !   (only required for the Chapman function calculations)
      REAL(KIND=8) :: HEIGHT_GRID ( 0:MAXLAYERS )

      !  multilayer atmospheric inputs (P andT)
      !   (only required for the Chapman function calculations, refractive geometry)
      REAL(KIND=8) :: PRESSURE_GRID    ( 0:MAXLAYERS )
      REAL(KIND=8) :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      REAL(KIND=8) :: AIRDENS          ( 1:MAXLAYERS )
      REAL(KIND=8) :: GAS_PROFILE      ( 1:MAXLAYERS )

      !  Number of fine layer gradations 
      !    (only for Chapman function calculations with refractive geometry)
      INTEGER      :: FINEGRID ( MAXLAYERS )

      !  optical property inputs
      !  -----------------------

      !  multilayer optical property (bulk) inputs
      REAL(KIND=8) :: OMEGA_TOTAL_INPUT  ( MAXLAYERS, MAXTHREADS )
      REAL(KIND=8) :: DELTAU_VERT_INPUT  ( MAXLAYERS, MAXTHREADS )

      !  Phase function Legendre-polynomial expansion coefficients
      !   Include all that you require for exact single scatter calculations
      REAL(KIND=8) :: PHASMOMS_TOTAL_INPUT
     &     ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

      !  Optical property linearizations
      !  Layer linearization (bulk property variation) input
      !  Layer linearization (phase function variation) input
      REAL(KIND=8)    
     &   :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(KIND=8)    
     &   :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(KIND=8)    
     &   :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,
     &      MAXLAYERS,MAXTHREADS)

      !  5. Surface inputs (New BRDF arrays, 23 March 2010)
      !  ==================================================

      !  Lambertian Surface control
      REAL(KIND=8) :: LAMBERTIAN_ALBEDO (MAXTHREADS)

      !  Exact (direct bounce) BRDF (same all threads)
      REAL(KIND=8) 
     &   :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, 
     &      MAXBEAMS )

      !  Linearized Exact (direct bounce) BRDF (same all threads)
      REAL(KIND=8)    :: LS_EXACTDB_BRDFUNC 
     &      ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, 
     &        MAXBEAMS )

      !  Fourier components of BRDF, in the following order (same all threads)

      !    incident solar directions,   reflected quadrature streams
      !    incident quadrature streams, reflected quadrature streams
      !    incident solar directions,   reflected user streams
      !    incident quadrature streams, reflected user streams
      REAL(KIND=8) :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, 
     &   MAXBEAMS )
      REAL(KIND=8) :: BRDF_F        ( 0:MAXMOMENTS, MAXSTREAMS, 
     &   MAXSTREAMS )
      REAL(KIND=8) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS,
     &   MAXBEAMS )
      REAL(KIND=8) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, 
     &   MAXSTREAMS )

      !  Linearized Fourier components of BRDF, in the following order (same all threads)

      !    incident solar directions,   reflected quadrature streams
      !    incident quadrature streams, reflected quadrature streams
      !    incident solar directions,   reflected user streams
      !    incident quadrature streams, reflected user streams
      REAL(KIND=8)    :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 
     &   0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(KIND=8)    :: LS_BRDF_F        ( MAX_SURFACEWFS, 
     &   0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(KIND=8)    :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 
     &   0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(KIND=8)    :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 
     &   0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

      !  Subroutine output arguments (Intent(out))
      !  +++++++++++++++++++++++++++++++++++++++++

      !  Main results
      !  ------------

      !  Intensity Results at all angles and optical depths
      REAL(KIND=8) :: INTENSITY 
     &  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, 
     &   MAXTHREADS )

      !  Results for mean-value output
      REAL(KIND=8) :: MEAN_INTENSITY 
     &   (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
      REAL(KIND=8) :: FLUX_INTEGRAL 
     &   (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)

      !  Profile weighting functions
      REAL(kind=8)    :: PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,      
     &   MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

      !  mean intensity weighting functions
      REAL(kind=8)    :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, 
     &   MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
 
      !  flux weighting functions
      REAL(kind=8)    :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS,
     &    MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

      !  Surface weighting functions at user angles
      REAL(kind=8)    :: SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS,
     &                 MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

      !  Flux and mean intensity surface weighting functions
      REAL(kind=8)    :: MINT_SURFACEWF ( MAX_SURFACEWFS, 
     &   MAX_USER_LEVELS,  MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
      REAL(kind=8)    :: FLUX_SURFACEWF ( MAX_SURFACEWFS, 
     &   MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

      !  Bookkeeping output
      !  ------------------

      !  Fourier numbers used
      INTEGER      :: FOURIER_SAVED ( MAXBEAMS, MAXTHREADS ) 

      !  Number of geometries (bookkeeping output)
      INTEGER      :: N_GEOMETRIES


      !  Thread number
      INTEGER  ::          THREAD


      !  Another copy for FD test 
      REAL(KIND=8) :: FLUX_INTEGRAL_BASE
     &   (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)

      ! (dkh, 02/08/08) 
      INTEGER       :: V, N
      INTEGER       :: W, WF


      !  Flag for opening error output file
      LOGICAL         :: OPENFILEFLAG

      !  Help variables

      integer          :: n6,l,ndum,ldum,t,t1,t2,t3,t4,t5,k3,k4,k5
      integer          :: n_totalprofile_wfs,nf,na,LEN_STRING
      integer          :: npert, k2
      double precision :: kd, gaer, waer, taer, parcel, raywt, aerwt
      double precision :: aersca, aerext, molabs, molsca, totsca, totext
      double precision :: omega
      double precision :: molomg(maxlayers),molext(maxlayers)
      double precision :: aermoms(0:maxmoments_input), lamb, eps, epsfac
      double precision :: raymoms(0:2,maxlayers), ratio1, ratio2

      LOGICAL          :: STATUS_PHYSICS

      INTEGER          :: J, K, IJLOOP

      !=================================================================
      ! LIDORT_DRIVER begins here!
      !=================================================================

      ! Set LIDORT number of layers and Legendre moments
      nlayers        = LLPAR
      nmoments_input = NPM_MAX

      eps    = 0d0
      epsfac = 1.0d0 + eps

      !  Initialize linearized inputs
      DO_PROFILE_LINEARIZATION = .TRUE.
      DO_SURFACE_LINEARIZATION = .TRUE.
      N_TOTALPROFILE_WFS       = MAX_ATMOSWFS
      DO_BRDF_SURFACE          = .FALSE.

      !  Initialise
      do n = 1, nlayers
         layer_vary_number(n) = MAX_ATMOSWFS
         layer_vary_flag(n)   = .true.
      enddo

      !  Surface
      N_SuRFACE_WFS   = 1


!!$OMP PARALLEL DO         
!!$OMP+DEFAULT( SHARED )  
!!$OMP+PRIVATE( thread ) 
!!$OMP+PRIVATE( status_physics )
!!$OMP+PRIVATE( STATUS_INPUTCHECK, STATUS_CALCULATION )
!!$OMP+PRIVATE( CHECKACTIONS)
!!$OMP+PRIVATE( NCHECKMESSAGES, CHECKMESSAGES ) 
!!$OMP+PRIVATE( MESSAGE, TRACE_1, TRACE_2, TRACE_3 )       
!!$OMP+PRIVATE( N_GEOMETRIES )
!!$OMP+PRIVATE( AIRDENS, GAS_PROFILE )                                         ! local 
!!$OMP+PRIVATE( FINEGRID, HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID )       ! local
!!$OMP+PRIVATE( BEAM_SZAS )                                                    ! local
!!$OMP+PRIVATE( npert )
!!$OMP+PRIVATE( J )     
!!$OMP+PRIVATE( IJLOOP )
!!$OMP+PRIVATE( W, V, N, WF, K )
      DO THREAD = 1, JJPAR
      !DO THREAD = JFD, JFD

       ! Use THREAD to loop over the latitude index J
       J = THREAD  

       ! First check to see if it is daytime.
       ! Calc 1-D array index
       IJLOOP = ( J - 1 ) * IIPAR + I

       ! Only consider the sunny side
       ! CYCLE leads to OpenMP issues 
       !IF ( SUNCOS(IJLOOP) < 0.0005d0 ) CYCLE 
       IF ( SUNCOS(IJLOOP) >= 0.0005d0 ) THEN

         ! Get LAT / LON / ALT dependent physical properties 
         CALL PREPARE_PHYSICS_GC_IJ( I, J, PRESSURE_GRID, 
     &       TEMPERATURE_GRID, HEIGHT_GRID, AIRDENS, GAS_PROFILE, 
     &       BEAM_SZAS )
  
         IF ( THREAD == 1 ) then 
            npert  = 1
            t1     = THREAD

         ELSEIF ( THREAD == 2 ) then 
            npert  = 21
            t2     = THREAD
            k2     = npert 

         ELSEIF ( THREAD == 3 ) then 
            npert =  23
            t3     = THREAD
            k3     = npert 

         ELSEIF ( THREAD == 4 ) then
            npert =  40
            t4     = THREAD
            k4     = npert 

         ELSEIF ( THREAD == 5 ) then
            !npert =  40
            t5     = THREAD
            k5     = npert 

         ENDIF 

         ! Loop over wavelengths    
         DO W = 1, NWL_MAX

            !  Prepare wavelength dependent physics, abort if failed
            CALL PREPARE_PHYSICS_GC_IJW
     &         ( thread, epsfac, npert, 
     &           status_physics,
     &           AIRDENS, 
     &           GAS_PROFILE, 
     &           LAMBERTIAN_ALBEDO(THREAD),
     &           OMEGA_TOTAL_INPUT(:,THREAD),
     &           DELTAU_VERT_INPUT(:,THREAD), 
     &           PHASMOMS_TOTAL_INPUT(:,:,THREAD),
     &           L_OMEGA_TOTAL_INPUT(:,:,THREAD),
     &           L_DELTAU_VERT_INPUT(:,:,THREAD),
     &           L_PHASMOMS_TOTAL_INPUT(:,:,:,THREAD),
     &           I, J, W )

            IF ( status_physics .EQV. .TRUE. ) STOP'physics_fail---'

            ! Call the actual LIDORT routines 
            CALL LIDORT_LPS_MASTER                                              
     &       ( THREAD,                  DO_SOLAR_SOURCES,                    ! Input
     &         DO_UPWELLING,            DO_DNWELLING,                        ! Input
     &         DO_FULLRAD_MODE,         DO_USER_STREAMS,                     ! Input/Output
     &         DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,                  ! Input/Output
     &         DO_SSFULL,               DO_SSCORR_TRUNCATION,                ! Input/Output
     &         DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,                     ! Input
     &         DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                       ! Input/Output
     &         DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,              ! Input/Output
     &         DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,                  ! Input/Output  
     &         DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,                   ! Input/Output
     &         DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,                  ! Input/Output
     &         DO_ALL_FOURIER,          DO_NO_AZIMUTH,                       ! Input/Output
     &         DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,           ! Input
     &         NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,               ! Input (Nmoment_input re-set)
     &         NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,        ! Input
     &         LAYER_VARY_FLAG,    LAYER_VARY_NUMBER,  N_SURFACE_WFS,        ! Input
     &         FLUX_FACTOR, LIDORT_ACCURACY,                                 ! Input
     &         EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,         ! Input
     &         BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_LEVELS,      ! Input
     &         FINEGRID, HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,       ! Input
     &         DELTAU_VERT_INPUT,    L_DELTAU_VERT_INPUT,                    ! Input
     &         OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,                    ! Input
     &         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                 ! Input
     &         LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC,       ! Input
     &         BRDF_F_0,    BRDF_F,    USER_BRDF_F_0,    USER_BRDF_F,        ! Input
     &         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,     ! Input
     &         INTENSITY, MEAN_INTENSITY, FLUX_INTEGRAL,                     ! Output
     &         N_GEOMETRIES, FOURIER_SAVED,                                  ! Output
     &         PROFILEWF, MINT_PROFILEWF, FLUX_PROFILEWF,                    ! Output
     &         SURFACEWF, MINT_SURFACEWF, FLUX_SURFACEWF,                    ! Output
     &         STATUS_INPUTCHECK,  NCHECKMESSAGES, CHECKMESSAGES,            ! Check 
     &         CHECKACTIONS,                                                 ! Check
     &         STATUS_CALCULATION, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )      ! Check

!            !  Exception handling, write-up (optional)
!            !    Will generate file only if errors or warnings are encountered
!            !      Unit file number  = 35, Filename = '3p5T_LIDORT_Execution.log'
!            CALL LIDORT_STATUS                                                 
!     &       ( THREAD, '3p5T_LIDORT_Execution.log', 35, OPENFILEFLAG,          
!     &         STATUS_INPUTCHECK,  NCHECKMESSAGES, CHECKMESSAGES, 
!     &         CHECKACTIONS,
!     &         STATUS_CALCULATION, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )


            ! Save results in global arrays. 
            ! First loop over beams 
            DO V  = 1, NBEAMS

               ! only save output at N = 1 (i.e., TOA)
               N  = 1

               ! Save flux
               GLOB_FLUX_W(W,I,J) = GLOB_FLUX_W(W,I,J)
     &            + FLUX_INTEGRAL(N,V,UPIDX,THREAD)


               ! Save flux profile weighting functions (derivatives)
               DO WF = 1, MAX_ATMOSWFS
               DO K  = 1, NLAYERS
                  JAC_GLOB_FLUX_W(WF, K,W,I,J)
     &                = JAC_GLOB_FLUX_W(WF,K,W,I,J)
     &                + FLUX_PROFILEWF(WF,K,N,V,UPIDX,THREAD)
               ENDDO
               ENDDO
           ENDDO


         ENDDO ! W 
      
       ENDIF ! SUNCOS 
      ENDDO 
!!$OMP END PARALLEL DO


!      !------------------------------------------------------------
!      ! Repeat for FD calculation
!      !------------------------------------------------------------
!
!      ! set perturbation
!      eps    = 1.0d-03
!      epsfac = 1.0d0 + eps
!
!      ! save a copy of baseline calculations
!      flux_integral_base = flux_integral
!
!!$OMP PARALLEL DO         
!!$OMP+DEFAULT( SHARED )  
!!$OMP+PRIVATE( thread ) 
!!$OMP+PRIVATE( status_physics )
!!$OMP+PRIVATE( STATUS_INPUTCHECK, STATUS_CALCULATION )
!!$OMP+PRIVATE( CHECKACTIONS)
!!$OMP+PRIVATE( NCHECKMESSAGES, CHECKMESSAGES ) 
!!$OMP+PRIVATE( MESSAGE, TRACE_1, TRACE_2, TRACE_3 )       
!!$OMP+PRIVATE( N_GEOMETRIES )
!!$OMP+PRIVATE( AIRDENS, GAS_PROFILE )                                         ! local 
!!$OMP+PRIVATE( FINEGRID, HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID )       ! local
!!$OMP+PRIVATE( BEAM_SZAS )                                                    ! local
!!$OMP+PRIVATE( npert )
!!$OMP+PRIVATE( J )     
!!$OMP+PRIVATE( IJLOOP )
!!$OMP+PRIVATE( W, V, N, WF, K )
!      DO THREAD = 1, 5
!
!         ! Use THREAD to loop over the latitude index J
!         !J = THREAD  
!         J = JFD 
!
!         ! First check to see if it is daytime.
!         ! Calc 1-D array index
!         IJLOOP = ( J - 1 ) * IIPAR + I
!
!         !! Only consider the sunny side
!         IF ( SUNCOS(IJLOOP) < 0.0005d0 ) CYCLE 
!
!         ! Get LAT / LON / ALT dependent physical properties 
!         CALL PREPARE_PHYSICS_GC_IJ( I, J, PRESSURE_GRID, 
!     &       TEMPERATURE_GRID, HEIGHT_GRID, AIRDENS, GAS_PROFILE, 
!     &       BEAM_SZAS )
!
!         IF ( THREAD == 1 ) then 
!            npert  = 1
!            t1     = THREAD
!
!         ELSEIF ( THREAD == 2 ) then 
!            npert  = 21
!            t2     = THREAD
!            k2     = npert 
!
!         ELSEIF ( THREAD == 3 ) then 
!            npert =  23
!            t3     = THREAD
!            k3     = npert 
!
!         ELSEIF ( THREAD == 4 ) then
!            npert =  40
!            t4     = THREAD
!            k4     = npert 
!
!         ELSEIF ( THREAD == 5 ) then
!            npert =  1
!            t5     = THREAD
!
!         ENDIF 
!
!
!         ! Just check a single wavelength
!         W = 1
!
!         !  Prepare physics, abort if failed
!         CALL PREPARE_PHYSICS_GC_IJW
!     &         ( thread, epsfac, npert, 
!     &           status_physics,
!     &           AIRDENS, 
!     &           GAS_PROFILE, 
!     &           LAMBERTIAN_ALBEDO(THREAD),
!     &           OMEGA_TOTAL_INPUT(:,THREAD),
!     &           DELTAU_VERT_INPUT(:,THREAD), 
!     &           PHASMOMS_TOTAL_INPUT(:,:,THREAD),
!     &           L_OMEGA_TOTAL_INPUT(:,:,THREAD),
!     &           L_DELTAU_VERT_INPUT(:,:,THREAD),
!     &           L_PHASMOMS_TOTAL_INPUT(:,:,:,THREAD),
!     &           I, J, W )
!            IF ( status_physics .EQV. .TRUE. ) STOP'physics_fail---'
!
! 
!         ! Call the actual LIDORT routines 
!         CALL LIDORT_LPS_MASTER                                              
!     &       ( THREAD,                  DO_SOLAR_SOURCES,                    ! Input
!     &         DO_UPWELLING,            DO_DNWELLING,                        ! Input
!     &         DO_FULLRAD_MODE,         DO_USER_STREAMS,                     ! Input/Output
!     &         DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,                  ! Input/Output
!     &         DO_SSFULL,               DO_SSCORR_TRUNCATION,                ! Input/Output
!     &         DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,                     ! Input
!     &         DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                       ! Input/Output
!     &         DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,              ! Input/Output
!     &         DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,                  ! Input/Output  
!     &         DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,                   ! Input/Output
!     &         DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,                  ! Input/Output
!     &         DO_ALL_FOURIER,          DO_NO_AZIMUTH,                       ! Input/Output
!     &         DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,           ! Input
!     &         NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,               ! Input (Nmoment_input re-set)
!     &         NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,        ! Input
!     &         LAYER_VARY_FLAG,    LAYER_VARY_NUMBER,  N_SURFACE_WFS,        ! Input
!     &         FLUX_FACTOR, LIDORT_ACCURACY,                                 ! Input
!     &         EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,         ! Input
!     &         BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_LEVELS,      ! Input
!     &         FINEGRID, HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,       ! Input
!     &         DELTAU_VERT_INPUT,    L_DELTAU_VERT_INPUT,                    ! Input
!     &         OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,                    ! Input
!     &         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                 ! Input
!     &         LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC,       ! Input
!     &         BRDF_F_0,    BRDF_F,    USER_BRDF_F_0,    USER_BRDF_F,        ! Input
!     &         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,     ! Input
!     &         INTENSITY, MEAN_INTENSITY, FLUX_INTEGRAL,                     ! Output
!     &         N_GEOMETRIES, FOURIER_SAVED,                                  ! Output
!     &         PROFILEWF, MINT_PROFILEWF, FLUX_PROFILEWF,                    ! Output
!     &         SURFACEWF, MINT_SURFACEWF, FLUX_SURFACEWF,                    ! Output
!     &         STATUS_INPUTCHECK,  NCHECKMESSAGES, CHECKMESSAGES,            ! Check 
!     &         CHECKACTIONS,                                                 ! Check
!     &         STATUS_CALCULATION, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )      ! Check
!
!        !  Exception handling, write-up (optional)
!        !    Will generate file only if errors or warnings are encountered
!        !      Unit file number  = 35, Filename = '3p5T_LIDORT_Execution.log'
!
!         CALL LIDORT_STATUS                                                 
!     & ( THREAD, '3p5T_LIDORT_Execution.log', 35, OPENFILEFLAG,          
!     &   STATUS_INPUTCHECK,  NCHECKMESSAGES, CHECKMESSAGES, 
!     &   CHECKACTIONS,
!     &   STATUS_CALCULATION, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )
!
!      ENDDO 
!!$OMP END PARALLEL DO
!
!
!      ! Display FD vs ADJ ouput 
!      OPEN(36,file = 'Tester_LPS_MT.all', status = 'unknown')
!      write(36,'(/T32,a/T32,a/)')
!     &   'REGULAR FLUX, 3 PROFILE JACOBIANS, 1 SURFACE JACOBIAN',
!     & '====================================================='
!      write(36,'(a,T32,a,a,a/)')' Sun SZA    Level/Output', 
!     & 'Regular F1     Profile AJ2   Profile FD2     Profile AJ3',
!     & '   Profile FD3     Profile AJ4   Profifle FD4',
!     & '    Surface AJ5   Surface FD5'
!      do v = 1, nbeams
!  
!         ! upwelling 
!         do n = 1, n_user_levels
!            write(36,366)v,'Upwelling @',user_levels(n), 
!     &      flux_integral_base(n,v,upidx,t1), 
!     &     flux_profilewf(1,k2,n,v,upidx,t2),
!     &    (flux_integral(n,v,upidx,t2)
!     &       - flux_integral_base(n,v,upidx,t2))/eps,
!     &     flux_profilewf(2,k3,n,v,upidx,t3),
!     &    (flux_integral(n,v,upidx,t3)
!     &       - flux_integral_base(n,v,upidx,t3))/eps,
!     &     flux_profilewf(2,k4,n,v,upidx,t4),
!     &    (flux_integral(n,v,upidx,t4)
!     &       - flux_integral_base(n,v,upidx,t4))/eps,
!     &     LAMBERTIAN_ALBEDO(t5)/epsfac*flux_surfacewf(1,n,v,upidx,t5),
!     &    (flux_integral(n,v,upidx,t5)
!     &       - flux_integral_base(n,v,upidx,t5))/eps      
!         enddo
!     
!!         ! downwelling 
!!         do n = 1, n_user_levels
!!            write(36,366)v,'Dnwelling @',user_levels(n), 
!!     &         flux_integral(n,v,dnidx,t1), 
!!     &    flux_profilewf(1,k2,n,v,dnidx,t2),
!!     &   (flux_integral(n,v,dnidx,t2)
!!     &       - flux_integral_base(n,v,dnidx,t2))/eps,
!!     &    flux_profilewf(2,k3,n,v,dnidx,t3),
!!     &   (flux_integral(n,v,dnidx,t3)
!!     &       - flux_integral_base(n,v,dnidx,t3))/eps,
!!     &    flux_profilewf(2,k4,n,v,dnidx,t4),
!!     &   (flux_integral(n,v,dnidx,t4)
!!     &       - flux_integral_base(n,v,dnidx,t4))/eps,
!!     &     LAMBERTIAN_ALBEDO(t5)/epsfac*flux_surfacewf(1,n,v,dnidx,t5),
!!     &    (flux_integral(n,v,dnidx,t5)
!!     &       - flux_integral_base(n,v,dnidx,t5))/eps      
!!         enddo
!         write(36,*)' '
!      enddo
!
!      close(36)
!366   format(i5,T11,a,f6.2,2x,1pe13.6,4(2x,1p2e14.6))



      !  Close error file if it has been opened
      IF ( OPENFILEFLAG ) then
        CLOSE(35)
        IF ( STATUS_INPUTCHECK .eq. LIDORT_WARNING .and.        
     &      STATUS_CALCULATION .eq. LIDORT_SUCCESS ) then 
          write(*,*)'3p5T.exe: program executed with internal '// 
     &     'defaults, warnings in "3p5T_LIDORT_Execution.log"'
        ELSE
          write(*,*)'3p5T.exe: program Failed to execute, '// 
     &     'Error messages in "3p5T_LIDORT_Execution.log"'
        ENDIF
      ELSE
        !write(*,*)'3p5T.exe: program finished successfully'
      ENDIF


      ! Return to calling program
      END SUBROUTINE LIDORT_DRIVER

!------------------------------------------------------------------------------

      SUBROUTINE CALC_RADIATIVE_FLUX( GLOB_FLUX_INST )

!******************************************************************************
!  Subroutine CALC_RADIATIVE_FLUX calculates the TOA upward radiative flux 
!  [W/m2]. 
!
!  Module variables as Input:
!  ============================================================================
!  (1 ) SUNCOS      (REAL*8) : Cosine of solar zenith angle [rad]
!  (2 ) GLOB_FLUX_W (REAL*8) : Percent of insolation reflected upward [%]
!     
!  Module variables as Output:
!  ============================================================================
!  (1 ) GLOB_FLUX   (REAL*8) : Running total of upward ratiative flux [W/m2/nsteps]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,      ONLY : SUNCOS
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TIME_MOD,     ONLY : GET_NYMD,   GET_NHMS


      ! dkh debug
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD 
#     include "CMN_SIZE"     ! IIPAR, JJPAR

      ! Arguments
      REAL*8, INTENT(INOUT) :: GLOB_FLUX_INST(IIPAR,JJPAR)
    
      ! Local variables 
      INTEGER               :: I, J, W, IJLOOP
      REAL*8                :: TOTAL_FLUX(IIPAR,JJPAR)
      REAL*8                :: GLOB_TEMP(IIPAR,JJPAR)

      ! Parameters 
      REAL*8, PARAMETER     :: EARTHSA = 510705155749195  ! [m2]


      !=================================================================
      ! CALC_RADIATIVE_FLUX begins here!
      !=================================================================

      ! Initialize storage arrays
      TOTAL_FLUX(:,:)      = 0d0
      GLOB_TEMP(:,:)       = 0d0
      GLOB_FLUX_INST(:,:)  = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, W, IJLOOP )    
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      !DO J = JFD, JFD
      !DO I = IFD, IFD

         ! Calc 1-D array index
         IJLOOP = ( J - 1 ) * IIPAR + I

         ! Only consider the sunny side
         IF ( SUNCOS(IJLOOP) > 0.0005d0 ) THEN 

            ! Integrate over wavelengths
            DO W = 1, NWL_MAX

               ! Calculate instantaneous upward flux 
               !  = (% reflected) * (spectral band insolation)
               ! note: zenith angle already accounted for in GLOB_FLUX_W
               GLOB_FLUX_INST(I,J) = GLOB_FLUX_INST(I,J)
     &                             + GLOB_FLUX_W(W,I,J)
     &                             * SOL_FLUX_TABLE(W) 

               ! Track total incoming ratiation
               GLOB_TEMP(I,J)      = GLOB_TEMP(I,J) 
     &                             + SOL_FLUX_TABLE(W) * SUNCOS(IJLOOP)
           
            ENDDO
            
            ! Add current flux to the running total
            GLOB_FLUX(I,J)    = GLOB_FLUX(I,J) 
     &                        + GLOB_FLUX_INST(I,J) * INV_NSPAN
            
            ! Convert from [W/m2] to [W]
            TOTAL_FLUX(I,J)   = GLOB_FLUX(I,J) * GET_AREA_M2(J)
            GLOB_TEMP (I,J)   = GLOB_TEMP(I,J) * GET_AREA_M2(J)

         ENDIF  

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Update count
      !GLOB_FLUX_COUNT = GLOB_FLUX_COUNT + 1    

      ! dkh debug
      print*, ' sum of GLOB_FLUX = ', SUM(TOTAL_FLUX(:,:)) 
     &                                / EARTHSA, !/ GLOB_FLUX_COUNT,
     &                                SUM(GLOB_TEMP(:,:))
     &                                / EARTHSA
!     &                                GLOB_FLUX_COUNT

      ! Save out instantaneous forcing
      CALL MAKE_RAD_FILE( GET_NYMD(), GET_NHMS(), GLOB_FLUX_INST )


      ! Return to calling program
      END SUBROUTINE CALC_RADIATIVE_FLUX

!------------------------------------------------------------------------------

      SUBROUTINE CALC_RADIATIVE_FORCE( GLOB_FLUX_INST, COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_RADIATIVE_FORCE calculates the difference in the 
!  instantaneous upward radiative flux between the current simulation and data
!  from a previous simulation. (dkh, 03/05/08) 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL*8) : Cost function                        [none]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,    ONLY   : CLDFRC
      USE ERROR_MOD,  ONLY   : IT_IS_NAN, ERROR_STOP
      USE GRID_MOD,   ONLY   : GET_AREA_M2
      USE TIME_MOD,   ONLY   : GET_NYMD, GET_NHMS 

      ! dkh debug 
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD 

#     include "CMN_SIZE" 

      ! Arguments
      REAL*8, INTENT(IN)    :: GLOB_FLUX_INST(IIPAR,JJPAR)
      REAL*8, INTENT(INOUT) :: COST_FUNC 

      ! Local variables 
      INTEGER               :: I, J
      REAL*8                :: FORCE(IIPAR,JJPAR)
      REAL*8                :: GLOB_FLUX_NAT(IIPAR,JJPAR)

      ! Parameters 
      REAL*8, PARAMETER     :: EARTHSA = 510705155749195  ! [m2]

      !=================================================================
      ! CALC_RADIATIVE_FORCE begins here!
      !=================================================================

      ! Initialize    
      FORCE(:,:)         = 0d0       
      GLOB_FLUX_NAT(:,:) = 0d0       

      ! Read radiative flux from previous simulation 
      !print*, ' *** SENSITIVITY CALC; USE NAT RAD = 0'
      !print*, ' *** SENSITIVITY CALC; USE NAT RAD = 0'
      !print*, ' *** SENSITIVITY CALC; USE NAT RAD = 0'
      !print*, ' *** SENSITIVITY CALC; USE NAT RAD = 0'
      GLOB_FLUX_NAT(:,:) = 0D0 
      !CALL READ_RAD_FILE( GET_NYMD(), GET_NHMS(), GLOB_FLUX_NAT ) 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J   )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      !DO J = JFD, JFD
      !DO I = IFD, IFD

         !print*, I, ' ', J
         !print*, I, ' ', J, ' ', GLOB_FLUX_INST(I,J),GLOB_AVE_FORCE(I,J)

         ! Forcing is (downward flux) - (downward natural flux). Since these
         ! GLOB_FLUX arrays are actually actualy upward fluxes, subtract 
         ! GLOB_FLUX from GLOB_FLUX_NAT.
         ! Include cloud fraction for clear-sky calculation (dkh, 08/03/11) 
         !FORCE(I,J)          = GLOB_FLUX_NAT(I,J) - GLOB_FLUX_INST(I,J)
         FORCE(I,J)      = ( GLOB_FLUX_NAT(I,J) - GLOB_FLUX_INST(I,J) )
     &                   * CLDFRC(I,J)

         ! Diagnostic: Track global average forcing [W/m2]
         GLOB_AVE_FORCE(I,J) = GLOB_AVE_FORCE(I,J)
     &                       + FORCE(I,J) * INV_NSPAN

         ! Weight it according to local area [W]
         FORCE(I,J)          = FORCE(I,J) * GET_AREA_M2(J)

         !print*, I, ' ', J, ' ', GLOB_FLUX_INST(I,J),GLOB_AVE_FORCE(I,J)
         !print*, I, ' ', J

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
     
      ! Update cost function with global average forcing [W/m2]
      COST_FUNC = COST_FUNC + SUM(FORCE(:,:)) / EARTHSA * INV_NSPAN

      ! dkh debug 
      !! **** HACK
      !COST_FUNC = COST_FUNC + GLOB_FLUX_INST(IFD,JFD) * INV_NSPAN
      !print*, ' *** COST_FUNC ONLY IN JFD, IFD *** '
      !print*, ' *** COST_FUNC ONLY IN JFD, IFD *** '
      !print*, ' *** COST_FUNC ONLY IN JFD, IFD *** '
      !print*, ' *** COST_FUNC ONLY IN JFD, IFD *** '
      !print*, ' *** COST_FUNC ONLY IN JFD, IFD *** '

      ! Error checking and print out
      IF ( IT_IS_NAN( COST_FUNC ) ) THEN
         CALL ERROR_STOP( 'COST_FUNC is NAN', 'CALC_AOD_FORCE')
      ELSE  
         WRITE(6,*) ' CALC_RADIATIVE_FORC: Current COST_FUN = ', 
     &              COST_FUNC
      ENDIF 



      ! Return to calling program
      END SUBROUTINE CALC_RADIATIVE_FORCE

!------------------------------------------------------------------------------
      SUBROUTINE ADJ_CALC_RADIATIVE_FORCE( GLOB_FLUX_INST_ADJ )
!
!******************************************************************************
!  Subroutine ADJ_CALC_RADIATIVE_FORCE is the adjoint of CALC_RADIATIVE_FORCE.
!  (dkh, 04/03/08) 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) GLOB_FLUX_INST_ADJ (REAL*8) :                                    [none]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,    ONLY   : CLDFRC
      USE GRID_MOD,   ONLY   : GET_AREA_M2

#     include "CMN_SIZE" 

      ! Arguments
      REAL*8, INTENT(OUT) :: GLOB_FLUX_INST_ADJ(IIPAR,JJPAR)

      ! Local variables 
      INTEGER               :: J
      INTEGER               :: I

      ! Parameters 
      REAL*8, PARAMETER     :: EARTHSA = 510705155749195  ! [m2]

      !=================================================================
      ! ADJ_CALC_RADIATIVE_FORCE begins here!
      !=================================================================

      ! fwd code:
      !COST_FUNC += SUM( (GLOB_FLUX_NAT(I,J) - GLOB_FLUX_INST(I,J) ) * GET_AREA_M2(J) ) / EARTHSA * INV_NSPAN
      !             * CLDFRC(I,J)

      GLOB_FLUX_INST_ADJ(:,:) = 1D0 / EARTHSA * INV_NSPAN

      ! dkh debug 
      !! **** HACK
      !print*, ' *** HACK '
      !print*, ' *** HACK '
      !print*, ' *** HACK '
      !print*, ' *** HACK '
      !print*, ' force localy '
      !GLOB_FLUX_INST_ADJ(:,:) =  - 1D0 * INV_NSPAN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE(  I, J  )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 
         GLOB_FLUX_INST_ADJ(I,J) = - GLOB_FLUX_INST_ADJ(I,J) 
     &                           * GET_AREA_M2(J)
     &                           * CLDFRC(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ADJ_CALC_RADIATIVE_FORCE

!------------------------------------------------------------------------------

      SUBROUTINE ADJ_CALC_RADIATIVE_FLUX( GLOB_FLUX_INST_ADJ )

!******************************************************************************
!  Subroutine ADJ_CALC_RADIATIVE_FLUX 
!  (dkh, 04/03/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) GLOB_FLUX_INST_ADJ (REAL*8) : 
!     
!  Module variables as Output:
!  ============================================================================
!  (1 ) GLOB_FLUX_W_ADJ  (REAL*8) : 
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,      ONLY : SUNCOS
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD 

#     include "CMN_SIZE"     ! IIPAR, JJPAR

      ! Arguments
      REAL*8, INTENT(IN)    :: GLOB_FLUX_INST_ADJ(IIPAR,JJPAR)
    
      ! Local variables 
      INTEGER               :: I, J, W, IJLOOP

      !=================================================================
      ! ADJ_CALC_RADIATIVE_FLUX begins here!
      !=================================================================

      ! Initialize arrays
      GLOB_FLUX_W_ADJ(:,:,:)  = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, W, IJLOOP )    
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      !DO J = JFD, JFD
      !DO I = IFD, IFD

         ! Calc 1-D array index
         IJLOOP = ( J - 1 ) * IIPAR + I

         ! Only consider the sunny side
         IF ( SUNCOS(IJLOOP) > 0.0005d0 ) THEN 

            ! Integrate over wavelengths
            DO W = 1, NWL_MAX

               ! fwd code:
               !GLOB_FLUX_INST(I,J) = GLOB_FLUX_INST(I,J)
               !                    + GLOB_FLUX_W(W,I,J)
               !                    * SOL_FLUX_TABLE(W) 
               ! adj code:
               GLOB_FLUX_W_ADJ(W,I,J) = GLOB_FLUX_INST_ADJ(I,J) 
     &                                * SOL_FLUX_TABLE(W)

            ENDDO
            
         ENDIF  

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ADJ_CALC_RADIATIVE_FLUX

!------------------------------------------------------------------------------

      SUBROUTINE ADJ_LIDORT_DRIVER( I, J )
!
!******************************************************************************
!  Subroutine ADJ_LIDORT_DRIVER 
!  (dkh, 04/03/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER)     : Lon index                            [none]
!  (2 ) J (INTEGER)     : Lat index                            [none]
!     
!  Module variables as Input:
!  ============================================================================
!  (1 ) GLOB_FLUX_W_ADJ    
!  (2 ) JAC_GLOB_FLUX_W 
!     
!  NOTES:
! 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,      ONLY : IT_IS_NAN, ERROR_STOP
      USE COMODE_MOD,     ONLY : JLOP
      !USE COMODE_MOD,     ONLY : ADJ_CSPEC_FORCE
      USE DAO_MOD,        ONLY : BXHEIGHT
      USE TRACERID_MOD,   ONLY : IDO3

      !  Include files
#     include "CMN_SIZE"       ! LLPAR, LLTROP

       ! Arguments 
      INTEGER, INTENT(IN)    :: I, J

      ! (dkh, 02/08/08) 
      INTEGER                :: W, L, UA, K 
      INTEGER                :: JLOOP
     
      REAL*8                 :: TAU_GAS
      REAL*8                 :: ADJ_GAS_PROFILE

      REAL*8, PARAMETER      :: DU_TO_CM2 = 2.68668d+16

      !=================================================================
      ! ADJ_LIDORT_DRIVER begins here!
      !=================================================================

      ! Loop over wavelengths 
      DO W  = 1, NWL_MAX

         DO K = 1, NLAYERS
     
            L = LLPAR - K + 1
       
! Adjoint of gas absorption.  Leave this out for the now 
!            TAU_GAS = GAS_PROFILE(K) * DU_TO_CM2 * O3_XSEC_TABLE(W)
!
!            ! Trace gas absorption [sw only]
!            ADJ_GAS_PROFILE          = JAC_GLOB_FLUX_W(1,K,W,I,J)
!     &                               * GLOB_FLUX_W_ADJ(W,I,J) 
!     &                               / TAU_GAS
!
!            ! Get 1-D array index for extracting O3 from CSPEC
!            JLOOP = JLOP(I,J,L)
!
!            ! Get O3 from CSPEC within the troposphere
!            IF ( JLOOP > 0 .and. L <= LLTROP ) THEN
!
!               ! fwd code:
!               !GAS_PROFILE(N) = CSPEC(JLOOP,IDO3)
!               !               * BXHEIGHT(I,J,L) * 100d0
!               !               * 3.7219d-17
!
!               ADJ_CSPEC_FORCE(JLOOP,IDO3) = ADJ_CSPEC_FORCE(JLOOP,IDO3)
!     &                                     + ADJ_GAS_PROFILE
!     &                                     * BXHEIGHT(I,J,L) * 100d0
!     &                                     * 3.7219d-17
!            ENDIF


 
            ! Layer aerosol optical depth
            LAYER_AOD_ADJ(I,J,L,W)   = JAC_GLOB_FLUX_W(2,K,W,I,J) 
     &                               * GLOB_FLUX_W_ADJ(W,I,J) 
     &                               / LAYER_AOD(I,J,L,W)

            ! Layer aerosol single scatering albedo 
            LAYER_SSA_ADJ(I,J,L,W)   = JAC_GLOB_FLUX_W(3,K,W,I,J) 
     &                               * GLOB_FLUX_W_ADJ(W,I,J)
     &                               / LAYER_SSA(I,J,L,W)

            ! Layer aerosol phase momemnts          
            LAYER_PHM_ADJ(I,J,L,1,W) = JAC_GLOB_FLUX_W(4,K,W,I,J) 
     &                               * GLOB_FLUX_W_ADJ(W,I,J)
     &                               / LAYER_PHM(I,J,L,1,W)
  
            ! Check if we've actually computed jacobians of all these PHMs
            IF ( MAX_ATMOSWFS == 8 ) THEN 
               LAYER_PHM_ADJ(I,J,L,2:5,W) = JAC_GLOB_FLUX_W(5:8,K,W,I,J) 
     &                                  * GLOB_FLUX_W_ADJ(W,I,J)
     &                                  / LAYER_PHM(I,J,L,2:5,W)
               LAYER_PHM_ADJ(I,J,L,6:NPM_MAX,W) = 0D0
            ELSE 
               LAYER_PHM_ADJ(I,J,L,2:NPM_MAX,W) = 0D0
            ENDIF 

         ENDDO   

      ENDDO 

       ! ********* HACK  works, at least close...
       !LAYER_AOD_ADJ(:,:,:,:)               = 0d0
       !LAYER_SSA_ADJ(:,:,:,:)               = 0d0
       !LAYER_PHM_ADJ(:,:,:,:,:)             = 0d0
       !LAYER_AOD_ADJ(IFD,JFD,6,4)           = 1d0
       !LAYER_SSA_ADJ(IFD,JFD,5,4:4)           = 1d0
       ! LAYER_PHM_ADJ(IFD,JFD,6,1:5,:)       = 1d0

      ! Return to calling program
      END SUBROUTINE ADJ_LIDORT_DRIVER

!------------------------------------------------------------------------------

      SUBROUTINE ADJ_CALC_AOD( I, J )
!
!******************************************************************************
!  Subroutine ADJ_CALC_AOD 
!  (dkh, 04/06/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I, J (Integer) : X-Y Grid locations
!     
!  Output:
!  ============================================================================
!  STT_ADJ and ADJ_FORCE
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE CHECKPT_MOD,    ONLY : RP_OUT,  CHK_STT
      USE DAO_MOD,        ONLY : RH,      BXHEIGHT, AIRVOL
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE TRACERID_MOD,   ONLY : IDTSO4,  IDTNIT,   IDTNH4
      USE TRACERID_MOD,   ONLY : IDTBCPO, IDTBCPI
      USE TRACERID_MOD,   ONLY : IDTOCPO, IDTOCPI
      ! fha 6-6-2011 added OC
	  ! dkh debug 
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(IN)     :: I, J
    
      ! Local variables 
      INTEGER                 :: L,    W
      INTEGER                 :: NSP
      REAL*8                  :: DIAM(NSP_MAX)
      REAL*8                  :: WET_DIAM(NSP_MAX)
      REAL*8                  :: MASS(NSP_MAX)
      REAL*8                  :: NCONC(NSP_MAX)
      REAL*8                  :: BEXT, SSA, PHM(NPM_MAX)
      REAL*8                  :: RHL
      REAL*8                  :: BEXT_TOT
      REAL*8                  :: SCAT_TOT
      REAL*8                  :: PHM_TOT(NPM_MAX)
      !REAL*8            :: NCONC_TOT
 
      ! Local adjoint variables 
      REAL*8                  :: ADJ_DIAM(NSP_MAX)
      REAL*8                  :: ADJ_MASS(NSP_MAX)
      REAL*8                  :: ADJ_NCONC(NSP_MAX)
      REAL*8                  :: ADJ_BEXT, ADJ_SSA, ADJ_PHM(NPM_MAX)
      REAL*8                  :: ADJ_BEXT_TOT
      REAL*8                  :: ADJ_SCAT_TOT
      REAL*8                  :: ADJ_PHM_TOT(NPM_MAX)
      REAL*8                  :: ADJ_MASS_WET
      REAL*8                  :: ADJ_MASS_DRY
      !REAL*8            :: ADJ_NCONC_TOT
      REAL*8                  :: DUM

      ! for testing MIE derivatives 
      !REAL*8 :: DIAM_0,  NCONC_0, PERT
      !REAL*8 :: DIAM_P,  DIAM_N,  NCONC_p,        NCONC_n
      !REAL*8 :: BEXT_1p, BEXT_1n, SSA_1p,         SSA_1n
      !REAL*8 :: BEXT_2p, BEXT_2n, SSA_2p,         SSA_2n
      !REAL*8 :: BEXT_0,  SSA_0,   PHM_0(NPM_MAX)
      !REAL*8 :: PHM_2p(NPM_MAX),  PHM_2n(NPM_MAX)
      !REAL*8 :: ADJ_NCONC_1,      ADJ_NCONC_2
      !REAL*8 :: ADJ_DIAM_1,       ADJ_DIAM_2
      !REAL*8 :: PHM_1p(NPM_MAX),  PHM_1n(NPM_MAX)

      !=================================================================
      ! ADJ_CALC_AOD begins here!
      !=================================================================

             ! dkh debug
!            IF ( I == IFD .and. J == JFD  ) THEN 
!               print*, 'ddd loc w HACK ' 
!               print*, 'ddd loc w HACK ' 
!               print*, 'ddd loc w HACK ' 
!               print*, 'ddd loc w ', LAYER_AOD_ADJ(I,J,LFD,3) 
!               LAYER_AOD_ADJ = 0d0
!               LAYER_SSA_ADJ = 0d0
!               LAYER_PHM_ADJ = 0d0
!              !LAYER_AOD_ADJ(I,J,LFD,3) = 1d0
!              LAYER_SSA_ADJ(I,J,LFD,3) = 1d0
!            ENDIF


      ! note: already in a parallel loop over I,J
      DO L = 1 , LLPAR



         ! Recalculate NCONC, DIAM, MASS

         ! Get local RH, cap at 99%
         RHL = MIN(99.D0,RH(I,J,L))

         ! Loop over aerosol type
         DO NSP = 1, NSP_MAX 
	  !fha  (6/1/11) added OC
            IF ( NSP == NSP_BC ) THEN 

               ! Aerosol wet size (um)
               WET_DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, 
     &                                        DRY_DIAM(NSP)/2d0 ) 

               ! Aerosol mass [kg/box]
               MASS(NSP)  = REAL(CHK_STT(I,J,L,IDTBCPI),8)
     &                    + REAL(CHK_STT(I,J,L,IDTBCPO),8)

               ! Total effective diam is volume weighted average
               ! of wet (BCPI) and dry (BCPO) components.  Since
               ! BC assumed to have density of 1, then volume ~ mass
               DIAM(NSP)  = 
     &            (   REAL(CHK_STT(I,J,L,IDTBCPI),8) * WET_DIAM(NSP)
     &              + REAL(CHK_STT(I,J,L,IDTBCPO),8) * DRY_DIAM(NSP) )
     &            / MASS(NSP)

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: rcl GET_NCONC bc   ', MASS(NSP), 
!     &                     DIAM(NSP)
!               ENDIF 
!
!               IF ( I == IFD .and. J == JFD .and. L == 6   ) THEN
!          print*, 'rcl JJJ MASS = ', MASS(NSP)
!          print*, 'rcl JJJ DIAM = ', DIAM(NSP)
!          print*, 'rcl JJJ WET_DIAM = ', WET_DIAM(NSP)
!          print*, 'rcl JJJ DRY_DIAM = ', DRY_DIAM(NSP)
!          print*, 'rcl JJJ RHL      = ', RHL
!               ENDIF


               ! Aerosol number concentration.
               NCONC(NSP) = GET_NCONC( 
     &                      MASS(NSP), 
     &                      0d0, 
!------------------------------------------------------------------
! BUG FIX: since we don't have mass of water, estimate NCONC based 
!  on dry radius (dkh, 01/27/11) 
!     &                      DENSE(NSP), 1.d0,   DIAM(NSP), 
     &                      DENSE(NSP), 1.d0,   DRY_DIAM(NSP),
!------------------------------------------------------------------
     &                      SIGMA(NSP), AIRVOL(I,J,L)      )
            ELSEIF ( NSP == NSP_OC ) THEN 

               ! Aerosol wet size (um)
               WET_DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, 
     &                                        DRY_DIAM(NSP)/2d0 ) 

               ! Aerosol mass [kg/box]
               MASS(NSP)  = REAL(CHK_STT(I,J,L,IDTOCPI),8)
     &                    + REAL(CHK_STT(I,J,L,IDTOCPO),8)

               ! Total effective diam is volume weighted average
               ! of wet (BCPI) and dry (BCPO) components.  Since
               ! BC assumed to have density of 1, then volume ~ mass
               DIAM(NSP)  = 
     &            (   REAL(CHK_STT(I,J,L,IDTOCPI),8) * WET_DIAM(NSP)
     &              + REAL(CHK_STT(I,J,L,IDTOCPO),8) * DRY_DIAM(NSP) )
     &            / MASS(NSP)

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: rcl GET_NCONC bc   ', MASS(NSP), 
!     &                     DIAM(NSP)
!               ENDIF 
!
!               IF ( I == IFD .and. J == JFD .and. L == 6   ) THEN
!          print*, 'rcl JJJ MASS = ', MASS(NSP)
!          print*, 'rcl JJJ DIAM = ', DIAM(NSP)
!          print*, 'rcl JJJ WET_DIAM = ', WET_DIAM(NSP)
!          print*, 'rcl JJJ DRY_DIAM = ', DRY_DIAM(NSP)
!          print*, 'rcl JJJ RHL      = ', RHL
!               ENDIF


               ! Aerosol number concentration.
               NCONC(NSP) = GET_NCONC( 
     &                      MASS(NSP), 
     &                      0d0, 
!------------------------------------------------------------------
! BUG FIX: since we don't have mass of water, estimate NCONC based 
!  on dry radius (dkh, 01/27/11) 
!     &                      DENSE(NSP), 1.d0,   DIAM(NSP), 
     &                      DENSE(NSP), 1.d0,   DRY_DIAM(NSP),
!------------------------------------------------------------------
     &                      SIGMA(NSP), AIRVOL(I,J,L)      )
			ELSE
               ! try it this way:
               ! specify dry diam
               !DIAM = 0.1 
               !DIAM(NSP)  = DRY_DIAM(NSP)

                ! Aerosol mass [kg/box]
               MASS(NSP) = REAL(CHK_STT(I,J,L,IDTSO4),8)
     &                   + REAL(CHK_STT(I,J,L,IDTNIT),8)
     &                   + REAL(CHK_STT(I,J,L,IDTNH4),8)

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: rcl GET_NCONC sulf ', MASS(NSP), 
!     &                     DRY_DIAM(NSP)
!               ENDIF 

               ! get NCONC based on dry mass
               NCONC(NSP) = GET_NCONC( 
     &                      MASS(NSP),                     
     &                      0d0,        
     &                      DENSE(NSP), 1.d0, DRY_DIAM(NSP), 
     &                      SIGMA(NSP), AIRVOL(I,J,L)     )

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: rcl GET_WETD2 ', MASS(NSP), 
!     &                     REAL(RP_OUT(I,J,L,4),8), NCONC(NSP)
!               ENDIF 

!!               ! get wet diam 
!!               DIAM(NSP) = GET_WETD2(
!!     &                     MASS(NSP),
!!!     &                     REAL(RP_OUT(I,J,L,4),8), 
!!! HACK: for testing SO4
!!     &                     1000d0,
!!     &                     DENSE(NSP), 1.d0, NCONC(NSP),
!!     &                     SIGMA(NSP), AIRVOL(I,J,L)     )
               ! Aerosol wet size (um)
               DIAM(NSP)  = GET_WET_DIAM( RHL, NSP,
     &                                    DRY_DIAM(NSP)/2d0 )


            ENDIF    
         ENDDO

         ! Initialize these, which are summed over wavelength in ADJ_MIE_LOOKUP
         ADJ_NCONC(:) = 0D0 
         ADJ_DIAM(:)  = 0D0 
         ADJ_MASS(:)  = 0D0 

         ! Loop over wavelengths
         DO W = 1, NWL_MAX

            BEXT_TOT   = 0d0
            SCAT_TOT   = 0D0 
            PHM_TOT(:) = 0D0
            !NCONC_TOT  = 0D0
  
            ! Recalculate  BEXT_TOT, SCAT_TOT, PHM_TOT, NCONC_TOT
            DO NSP = 1, NSP_MAX

               ! Lookup aerosol optical properties
               CALL MIE_LOOKUP( NCONC(NSP), DIAM(NSP), W,  NSP,    
     &                          BEXT,       SSA,       PHM, I,J,L   )


               BEXT_TOT   = BEXT_TOT   + BEXT
               SCAT_TOT   = SCAT_TOT   + BEXT * SSA 
               PHM_TOT(:) = PHM_TOT(:) + PHM(:) * BEXT * SSA 
               !NCONC_TOT  = NCONC_TOT  + NCONC(NSP)

               ! dkh debug
               !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN 
               !  !print*, ' I, J, L, W, NSP, NCONC, DIAM, MASS, DENSE, 
               !  !BEXT, SSA in adj',  
               !  !I, J, L, W, NSP, NCONC(NSP), DIAM(NSP), MASS(NSP),
               !  !DENSE(NSP), BEXT, SSA
               !  print*, ' ick: rcl I,J,L,W,NSP ', I, J, L, W, NSP
               !  print*, ' ick: rcl NCONC ', NCONC(NSP)
               !  print*, ' ick: rcl DIAM  ', DIAM(NSP)
               !  print*, ' ick: rcl MASS  ', MASS(NSP)
               !  print*, ' ick: rcl DENSE ', DENSE(NSP)
               !  print*, ' ick: rcl BEXT  ', BEXT
               !  print*, ' ick: rcl SSA   ', SSA
               !ENDIF 

            ENDDO 

            IF ( BEXT_TOT  < SMALL .or. 
     &           SCAT_TOT  < SMALL      ) THEN
               CALL ERROR_STOP(' Trouble calculating combined opt prop', 
     &                         ' aero_optical_mod.f' )
            ENDIF 

            ! dkh debug
            !IF ( I == IFD .and. J == JFD ) THEN
            !   print*, ' jck: rcl I,J,L,W  ', I, J, L, W
            !   print*, ' jck: rcl BEXT_TOT ', BEXT_TOT
            !   print*, ' jck: rcl SCAT_TOT ', SCAT_TOT
            !   print*, ' jck: rcl LAYER_AOD', LAYER_AOD(I,J,L,W)
            !   print*, ' jck: rcl LAYER_SSA', LAYER_SSA(I,J,L,W)
            !   print*, ' jck: rcl RHL      ', RHL
            !ENDIF 


            ! Calculate adjoint of BEXT_TOT, SCAT_TOT, PHM_TOT, NCONC_TOT
            ! fwd code:
            !LAYER_AOD(I,J,L,W)   = BXHEIGHT(I,J,L) * BEXT_TOT 
            !LAYER_SSA(I,J,L,W)   = SCAT_TOT / BEXT_TOT ! / NCONC_TOT
            !LAYER_PHM(I,J,L,:,W) = PHM_TOT(:) / SCAT_TOT 
            ! adj code:
            ADJ_BEXT_TOT   = BXHEIGHT(I,J,L) 
     &                     * LAYER_AOD_ADJ(I,J,L,W) 
!     &                     * 1d-6
     &                     - SCAT_TOT / BEXT_TOT **2 !/ NCONC_TOT
     &                     * LAYER_SSA_ADJ(I,J,L,W)
            ADJ_SCAT_TOT   = 1.0 / BEXT_TOT !/ NCONC_TOT
     &                     * LAYER_SSA_ADJ(I,J,L,W)
     &                     - 1 / SCAT_TOT **2 
     &                     * SUM(PHM_TOT(:) * LAYER_PHM_ADJ(I,J,L,:,W) )
            ADJ_PHM_TOT(:) = LAYER_PHM_ADJ(I,J,L,:,W) / SCAT_TOT
!            ADJ_NCONC_TOT  = - SCAT_TOT  / BEXT_TOT / NCONC_TOT ** 2
!     &                     * LAYER_SSA_ADJ(I,J,L,W)

       !**** HACK
         !IF ( L == LFD .and. W == NWL_500 ) THEN 
         !   ADJ_BEXT_TOT = 0D0 
         !   ADJ_SCAT_TOT = 1D0 
         !   ADJ_PHM_TOT(:)  = 0D0 
         !ELSE
         !   ADJ_BEXT_TOT = 0D0 
         !   ADJ_SCAT_TOT = 0D0 
         !   ADJ_PHM_TOT(:) = 0D0 
         !ENDIF

            !print*, 'ADJ_BEXT_TOT = ', ADJ_BEXT_TOT
       
            DO NSP = 1, NSP_MAX

               ! Recalculate BEXT, SSA, PHM
               ! Lookup aerosol optical properties
               CALL MIE_LOOKUP( NCONC(NSP), DIAM(NSP), W,  NSP,    
     &                          BEXT,       SSA,       PHM, I,J,L   )

               ! dkh debug
               !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
               !  !print*, ' I, J, L, W, NSP, NCONC, DIAM, MASS, DENSE, 
               !  !BEXT, SSA in adj',  
               !  !I, J, L, W, NSP, NCONC(NSP), DIAM(NSP), MASS(NSP),
               !  !DENSE(NSP), BEXT, SSA
               !  print*, ' ick: rc2 I,J,L,W,NSP ', I, J, L, W, NSP
               !  print*, ' ick: rc2 NCONC ', NCONC(NSP)
               !  print*, ' ick: rc2 DIAM  ', DIAM(NSP)
               !  print*, ' ick: rc2 MASS  ', MASS(NSP)
               !  print*, ' ick: rc2 DENSE ', DENSE(NSP)
               !  print*, ' ick: rc2 BEXT  ', BEXT
               !  print*, ' ick: rc2 SSA   ', SSA
               !ENDIF 

               ! Calculate adjoint of BEXT, SSA, PHM, NCONC
               ! fwd code:
               !BEXT_TOT   = BEXT_TOT   + BEXT
               !SCAT_TOT   = SCAT_TOT   + BEXT * SSA 
               !PHM_TOT(:) = PHM_TOT(:) + PHM(:) * BEXT * SSA 
               !!NCONC_TOT  = NCONC_TOT  + NCONC(NSP)
               ! adj code:
               !ADJ_NCONC(NSP) = ADJ_NCONC_TOT
               ADJ_PHM(:)     = ADJ_PHM_TOT(:) * BEXT * SSA 
               ADJ_SSA        = ADJ_SCAT_TOT * BEXT
     &                        + SUM(ADJ_PHM_TOT(:) * PHM(:) ) * BEXT 
               ADJ_BEXT       = ADJ_BEXT_TOT
     &                        + ADJ_SCAT_TOT * SSA
     &                        + SUM(ADJ_PHM_TOT(:) * PHM(:) ) * SSA

               ! dkh debug
               !print*, 'ddd loc bu ', ADJ_SSA, ADJ_BEXT , W, NSP
               !print*, 'ddd loc bu p', ADJ_PHM(:)
               !print*, 'ddd loc bu f1', BEXT, SSA
               !print*, 'ddd loc bu f2', PHM(:)
               !print*, 'ddd loc bu f3', ADJ_PHM_TOT(:)

         !!**** HACK works
!         IF ( L == LFD .and. W == NWL_500 .and. NSP == 2 ) THEN 
!!            print*, ' HACK ddd loc bps ' 
!            ADJ_BEXT   = 1D0 
!            ADJ_PHM(:) = 0D0 
!            ADJ_SSA    = 0D0
!            ADJ_NCONC  = 0D0
!         ELSE 
!            ADJ_BEXT = 0D0 
!            ADJ_PHM(:) = 0D0 
!            ADJ_SSA    = 0D0
!         ENDIF

         
               ! Get ADJ_NCONC, ADJ_DIAM
               CALL ADJ_MIE_LOOKUP( DIAM(NSP), W,  NSP,    
     &                          NCONC(NSP), 
     &                          ADJ_NCONC(NSP), ADJ_DIAM(NSP), 
     &                          ADJ_BEXT,       ADJ_SSA,
     &                          ADJ_PHM, I,J,L   )

               !! dkh debug
               !print*, 'ddd loc bv ', ADJ_NCONC(NSP),ADJ_DIAM(NSP), W 

! Here is some code for testing MIE derivatives.  Should disable OMP 
! loop when calling ADJ_CALC_AOD prior to using this. Also, uncomment
! the block of "testing" variable declarations above (dkh, 03/27/11) 
!----------------------------------------------------------------------
!         IF ( I == IFD .and. L == LFD .and. W == 1 
!     &                 .and. NSP == 2 ) THEN 
!            ! FWD
!            DIAM_0 = DIAM(NSP)
!            NCONC_0 = NCONC(NSP)
!            PERT = 0.1d0 
! 
!            ! Lookup aerosol optical properties
!            CALL MIE_LOOKUP( NCONC(NSP), DIAM(NSP), W,  NSP,    
!     &                       BEXT_0,     SSA_0,     PHM_0, I,J,L   )
!           
!            NCONC_p = NCONC_0  * (1d0 + PERT)
!            NCONC_n = NCONC_0  * (1d0 - PERT)
!            
!            DIAM_p  = DIAM_0 * (1d0 + PERT)
!            DIAM_n  = DIAM_0 * (1d0 - PERT)
!
!            CALL MIE_LOOKUP( NCONC_p,    DIAM_0   , W,  NSP,   
!     &                       BEXT_1p,     SSA_1p,       PHM_1p, I,J,L )
!              
!            CALL MIE_LOOKUP( NCONC_n,    DIAM_0   , W,  NSP,   
!     &                       BEXT_1n,     SSA_1n,       PHM_1n, I,J,L )
!              
!            CALL MIE_LOOKUP( NCONC_0   , DIAM_p,    W,  NSP,   
!     &                       BEXT_2p,     SSA_2p,       PHM_2p, I,J,L  )
!
!            CALL MIE_LOOKUP( NCONC_0   , DIAM_n,    W,  NSP,   
!     &                       BEXT_2n,     SSA_2n,       PHM_2n, I,J,L  )
!
!
!           ADJ_BEXT = 1d0
!           ADJ_SSA  = 0d0
!           ADJ_PHM  = 0d0
!           ADJ_NCONC_1 = 0d0 
!           ADJ_DIAM_1  = 0d0 
!
!           ! Get ADJ_NCONC, ADJ_DIAM
!            CALL ADJ_MIE_LOOKUP( DIAM_0, W,  NSP,    
!     &                          NCONC_0,    
!     &                          ADJ_NCONC_1, ADJ_DIAM_1,
!     &                          ADJ_BEXT,       ADJ_SSA,
!     &                          ADJ_PHM, I,J,L   )
!
!           ADJ_BEXT = 0d0
!           ADJ_SSA  = 1d0
!           ADJ_PHM  = 0d0
!           ADJ_NCONC_2 = 0d0 
!           ADJ_DIAM_2  = 0d0 
!
!           ! Get ADJ_NCONC, ADJ_DIAM
!            CALL ADJ_MIE_LOOKUP( DIAM_0, W,  NSP,    
!     &                          NCONC_0,    
!     &                          ADJ_NCONC_2, ADJ_DIAM_2,
!     &                          ADJ_BEXT,       ADJ_SSA,
!     &                          ADJ_PHM, I,J,L   )
!
!           print*,' dBEXT/dNCONC, dSSA/dNCONC, dBEXT/dDIAM, dSSA/dDIAM',
!     &         I, J, L, W, NSP 
!            WRITE(6,100), ' FD 1st pos: ', ( BEXT_1p - BEXT_0 ) / PERT, 
!     &                                     ( SSA_1p  - SSA_0 )  / PERT,
!     &                                     ( BEXT_2p - BEXT_0)  / PERT,
!     &                                     ( SSA_2p  - SSA_0 )  / PERT 
!            WRITE(6,100), ' FD 1st neg: ', -( BEXT_1n - BEXT_0 ) / PERT, 
!     &                                     -( SSA_1n  - SSA_0 )  / PERT,
!     &                                     -( BEXT_2n - BEXT_0)  / PERT,
!     &                                     -( SSA_2n  - SSA_0 )  / PERT 
!            WRITE(6,100), ' FD 2nd    : ', 
!     &          ( BEXT_1p - BEXT_1n ) / ( 2d0 * PERT ), 
!     &          ( SSA_1p  - SSA_1n )  / ( 2d0 * PERT ),
!     &          ( BEXT_2p - BEXT_2n ) / ( 2d0 * PERT ), 
!     &          ( SSA_2p  - SSA_2n )  / ( 2d0 * PERT )
!            WRITE(6,100), ' ADJ       : ', ADJ_NCONC_1 * NCONC_0, 
!     &                                     ADJ_NCONC_2 * NCONC_0,
!     &                                     ADJ_DIAM_1  * DIAM_0,
!     &                                     ADJ_DIAM_2  * DIAM_0
!            
! 100  FORMAT(A14, 2x, es13.6, 2x, es13.6, 2x, es13.6, 2x, es13.6 )
!
!       ENDIF 
!----------------------------------------------------------------------



!               IF ( I == IFD .and. J == JFD .and.  L == LFD ) THEN 
!                  print*, 'ADJ_NCONC = ', ADJ_NCONC, NSP
!                  print*, 'ADJ_DIAM  = ', ADJ_DIAM , NSP
!               ENDIF 

               ! Reset
               ADJ_PHM(:)     = 0d0
               ADJ_SSA        = 0d0
               ADJ_BEXT       = 0d0

            ENDDO 

            ! reset some variables here
            ! fwd:
            !BEXT_TOT   = 0d0
            !SCAT_TOT   = 0D0 
            !PHM_TOT(:) = 0D0
            !NCONC_TOT  = 0D0
            ! adj:
            ADJ_BEXT_TOT   = 0D0 
            ADJ_SCAT_TOT   = 0D0
            ADJ_PHM_TOT(:) = 0D0
            !ADJ_NCONC_TOT  = 0D0

         ENDDO ! LAMBDA

!         ! dkh debug
!         IF ( I == IFD .and. J == JFD  .and. L == LFD ) THEN 
!             print*, ' ADJ_DIAM, ADJ_NCONC ' , ADJ_DIAM, ADJ_NCONC
!         ENDIF 

!            ! **** HACK  works
!               IF ( L == LFD   ) THEN 
!                  !print*, ' HACK ddd loc nd '
!                  ADJ_NCONC(2) = 1D0
!                  ADJ_DIAM(2) = 0d0
!               ELSE 
!                  ADJ_NCONC(2) = 0D0
!                  ADJ_DIAM(2) = 0d0
!               ENDIF 
!               ADJ_NCONC(1) = 0D0
!              ADJ_DIAM(1)  = 0d0


         ! Calculate adjoint of STT(SO4,NH4,NIT,BC)
         ! Loop over aerosol type
         DO NSP = 1, NSP_MAX 

               ! dkh debug
               !print*, 'ddd loc b xyz loop NSP ', NSP
			!fha  (6/1/11) added OC
            IF ( NSP == NSP_BC ) THEN 

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: ADJ_GET_NCONC bc   ', MASS(NSP), 
!     &                     DIAM(NSP)
!               ENDIF 


               ! Calculate adjoint of MASS, update adjoint of DIAM
               ! fwd:
               !NCONC(NSP) = GET_NCONC( 
               !             MASS(NSP), 
               !             0d0, 
               !             DENSE(NSP), 1.d0,   DIAM(NSP), 
               !             SIGMA(NSP), AIRVOL(I,J,L)      )
               ! adj:
               CALL ADJ_GET_NCONC(
     &                      MASS(NSP), 0d0, 
!----------------------------------------------------------
! BUG FIX: use dry diam for estimating bc NCONC.  Also, 
!  DIAM is not an active variable, so pass 0d0 instead of 
!  ADJ_DIAM.  (dkh, 01/27/11) 
!     &                      DENSE(NSP), 1.d0,   DIAM(NSP), 
!     &                      SIGMA(NSP), AIRVOL(I,J,L),
!     &                      ADJ_DIAM(NSP), ADJ_NCONC(NSP),  
     &                      DENSE(NSP), 1.d0,  DRY_DIAM(NSP), 
     &                      SIGMA(NSP), AIRVOL(I,J,L),
     &                      DUM,           ADJ_NCONC(NSP),  
!----------------------------------------------------------
     &                      ADJ_MASS(NSP) )
 
!! **** HACK works
!       IF ( I == IFD .and. J == JFD .and. L == LFD .and. NSP == 2 ) THEN
!          print*, ' force at ddd loc md ' 
!          ADJ_DIAM(NSP)  = 0D0 
!          ADJ_MASS(NSP)  = 1D0 
!          print*, 'adj JJJ MASS = ', MASS(NSP)
!          print*, 'adj JJJ DIAM = ', DIAM(NSP)
!          print*, 'adj JJJ WET_DIAM = ', WET_DIAM(NSP)
!          print*, 'adj JJJ DRY_DIAM = ', DRY_DIAM(NSP)
!          print*, 'adj JJJ RHL      = ', RHL
!       ELSE 
!          ADJ_DIAM(NSP) = 0D0 
!          ADJ_MASS(NSP) = 0D0 
!       ENDIF 

               ! fwd:
               !DIAM(NSP)  = 
               !   (   REAL(CHK_STT(I,J,L,IDTBCPI),8) * WET_DIAM(NSP)
               !     + REAL(CHK_STT(I,J,L,IDTBCPO),8) * DRY_DIAM(NSP) )
               !   / MASS(NSP)
               ! adj:
               ADJ_MASS(NSP)            =  - DIAM(NSP) / MASS(NSP) 
     &                                  * ADJ_DIAM(NSP) + ADJ_MASS(NSP)

! **** HACK
               STT_ADJ(I,J,L,IDTBCPI) = STT_ADJ(I,J,L,IDTBCPI)
     &                                  + WET_DIAM(NSP) / MASS(NSP) 
     &                                  * ADJ_DIAM(NSP)
               STT_ADJ(I,J,L,IDTBCPO) = STT_ADJ(I,J,L,IDTBCPO)
     &                                  + DRY_DIAM(NSP) / MASS(NSP)
     &                                  * ADJ_DIAM(NSP)

! **** HACK worked
        !IF ( L == 6 ) THEN 
        !   ADJ_MASS(NSP) = 1d0
        !ELSE 
        !   ADJ_MASS(NSP) = 0d0
        !ENDIF 
               ! fwd:
               !MASS(NSP)  = REAL(CHK_STT(I,J,L,IDTBCPI),8)
               !           + REAL(CHK_STT(I,J,L,IDTBCPO),8)
               ! adj: 
! **** HACK
               STT_ADJ(I,J,L,IDTBCPI) = STT_ADJ(I,J,L,IDTBCPI)
     &                                  + ADJ_MASS(NSP)
               STT_ADJ(I,J,L,IDTBCPO) = STT_ADJ(I,J,L,IDTBCPO)
     &                                  + ADJ_MASS(NSP)

               ! fwd:
               !WET_DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, DRY_DIAM(NSP)/2d0 ) 
               ! adj: None, since RH not an active variable, but could 
               ! calculates sensitivity w.r.t. dry diam and RH here. 


            
            ELSEIF ( NSP == NSP_OC  ) THEN 

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: ADJ_GET_NCONC bc   ', MASS(NSP), 
!     &                     DIAM(NSP)
!               ENDIF 


               ! Calculate adjoint of MASS, update adjoint of DIAM
               ! fwd:
               !NCONC(NSP) = GET_NCONC( 
               !             MASS(NSP), 
               !             0d0, 
               !             DENSE(NSP), 1.d0,   DIAM(NSP), 
               !             SIGMA(NSP), AIRVOL(I,J,L)      )
               ! adj:
               CALL ADJ_GET_NCONC(
     &                      MASS(NSP), 0d0, 
!----------------------------------------------------------
! BUG FIX: use dry diam for estimating bc NCONC.  Also, 
!  DIAM is not an active variable, so pass 0d0 instead of 
!  ADJ_DIAM.  (dkh, 01/27/11) 
!     &                      DENSE(NSP), 1.d0,   DIAM(NSP), 
!     &                      SIGMA(NSP), AIRVOL(I,J,L),
!     &                      ADJ_DIAM(NSP), ADJ_NCONC(NSP),  
     &                      DENSE(NSP), 1.d0,  DRY_DIAM(NSP), 
     &                      SIGMA(NSP), AIRVOL(I,J,L),
     &                      DUM,           ADJ_NCONC(NSP),  
!----------------------------------------------------------
     &                      ADJ_MASS(NSP) )
 
!! **** HACK works
!       IF ( I == IFD .and. J == JFD .and. L == LFD .and. NSP == 2 ) THEN
!          print*, ' force at ddd loc md ' 
!          ADJ_DIAM(NSP)  = 0D0 
!          ADJ_MASS(NSP)  = 1D0 
!          print*, 'adj JJJ MASS = ', MASS(NSP)
!          print*, 'adj JJJ DIAM = ', DIAM(NSP)
!          print*, 'adj JJJ WET_DIAM = ', WET_DIAM(NSP)
!          print*, 'adj JJJ DRY_DIAM = ', DRY_DIAM(NSP)
!          print*, 'adj JJJ RHL      = ', RHL
!       ELSE 
!          ADJ_DIAM(NSP) = 0D0 
!          ADJ_MASS(NSP) = 0D0 
!       ENDIF 

               ! fwd:
               !DIAM(NSP)  = 
               !   (   REAL(CHK_STT(I,J,L,IDTBCPI),8) * WET_DIAM(NSP)
               !     + REAL(CHK_STT(I,J,L,IDTBCPO),8) * DRY_DIAM(NSP) )
               !   / MASS(NSP)
               ! adj:
               ADJ_MASS(NSP)            =  - DIAM(NSP) / MASS(NSP) 
     &                                  * ADJ_DIAM(NSP) + ADJ_MASS(NSP)

! **** HACK
               STT_ADJ(I,J,L,IDTOCPI) = STT_ADJ(I,J,L,IDTOCPI)
     &                                  + WET_DIAM(NSP) / MASS(NSP) 
     &                                  * ADJ_DIAM(NSP)
               STT_ADJ(I,J,L,IDTOCPO) = STT_ADJ(I,J,L,IDTOCPO)
     &                                  + DRY_DIAM(NSP) / MASS(NSP)
     &                                  * ADJ_DIAM(NSP)

! **** HACK worked
        !IF ( L == 6 ) THEN 
        !   ADJ_MASS(NSP) = 1d0
        !ELSE 
        !   ADJ_MASS(NSP) = 0d0
        !ENDIF 
               ! fwd:
               !MASS(NSP)  = REAL(CHK_STT(I,J,L,IDTBCPI),8)
               !           + REAL(CHK_STT(I,J,L,IDTBCPO),8)
               ! adj: 
! **** HACK
               STT_ADJ(I,J,L,IDTOCPI) = STT_ADJ(I,J,L,IDTOCPI)
     &                                  + ADJ_MASS(NSP)
               STT_ADJ(I,J,L,IDTOCPO) = STT_ADJ(I,J,L,IDTOCPO)
     &                                  + ADJ_MASS(NSP)

               ! fwd:
               !WET_DIAM(NSP)  = GET_WET_DIAM( RHL, NSP, DRY_DIAM(NSP)/2d0 ) 
               ! adj: None, since RH not an active variable, but could 
               ! calculates sensitivity w.r.t. dry diam and RH here. 


            ELSE

!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: ADJ_GET_WETD2 ', MASS(NSP), 
!     &                     REAL(RP_OUT(I,J,L,4),8), NCONC(NSP)
!               ENDIF 

               ! dkh debug
               !print*, 'ddd loc bx  ', ADJ_NCONC(NSP),
!     &                ADJ_DIAM(NSP)
               !print*, 'ddd loc bx f1', MASS(NSP)
               !print*, 'ddd loc bx f2', REAL(RP_OUT(I,J,L,4),8)
               !print*, 'ddd loc bx f3', DENSE(NSP)
               !print*, 'ddd loc bx f4', NCONC(NSP)
               !print*, 'ddd loc bx f5', SIGMA(NSP)
               !print*, 'ddd loc bx f6', AIRVOL(I,J,L)

               ! fwd:
               !DIAM(NSP) = GET_WETD2(
               !            MASS(NSP),
               !            REAL(RP_OUT(I,J,L,4),8), 
               !            DENSE(NSP), 1.d0, NCONC(NSP),
               !            SIGMA(NSP), AIRVOL(I,J,L)     )
               ! adj: input  is ADJ_NCONC, ADJ_DIAM
               !      output is ADJ_MASS_WET, ADJ_MASS_DRY, ADJ_NCONC
!!               CALL ADJ_GET_WETD2(
!!     &                     MASS(NSP),
!!!     &                     REAL(RP_OUT(I,J,L,4),8), 
!!! HACK: for testing SO4
!!     &                     1000d0,
!!     &                     DENSE(NSP), 1.d0, NCONC(NSP),
!!     &                     SIGMA(NSP), AIRVOL(I,J,L),
!!     &                     ADJ_DIAM(NSP), ADJ_MASS_WET,
!!     &                     ADJ_MASS_DRY,  ADJ_NCONC(NSP) )

! **** HACK
!               ADJ_FORCE(I,J,L,IDADJH2O) = ADJ_MASS_WET
!               ADJ_MASS(NSP) = ADJ_MASS_DRY
               
!               ! dkh debug
!               IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
!                  print*, 'ichk: ADJ_GET_NCONC sulf ', MASS(NSP), 
!     &                     DRY_DIAM(NSP)
!               ENDIF 

               ! **** HACK 
               !IF ( L == LFD ) THEN 
               !   ADJ_NCONC(NSP) = 1D0
               !   ADJ_DIAM(NSP) = 0D0
               !   ADJ_MASS(NSP) = 0D0
               !ELSE 
               !   ADJ_NCONC(NSP) = 0D0
               !   ADJ_DIAM(NSP) = 0D0
               !   ADJ_MASS(NSP) = 0D0
               !ENDIF 

               ! dkh debug
!               print*, 'ddd loc by', ADJ_MASS(NSP), ADJ_NCONC(NSP),
!     &                ADJ_DIAM(NSP)


               ! fwd
               !NCONC(NSP) = GET_NCONC( 
               !             MASS(NSP),                     
               !             0d0,        
               !             DENSE(NSP), 1.d0, DRY_DIAM(NSP), 
               !             SIGMA(NSP), AIRVOL(I,J,L)     )
               ! adj: input is ADJ_DIAM, ADJ_NCONC, ADJ_MASS
               !      output is ADJ_MASS, ADJ_DIAM
               CALL ADJ_GET_NCONC(
     &                      MASS(NSP), 0d0, 
     &                      DENSE(NSP), 1.d0,   DRY_DIAM(NSP), 
     &                      SIGMA(NSP), AIRVOL(I,J,L),
     &                      ADJ_DIAM(NSP), ADJ_NCONC(NSP),  
     &                      ADJ_MASS(NSP) )
 
               ! dkh debug
               !print*, 'ddd loc bz', ADJ_MASS(NSP)

               ! **** HACK 
               !IF ( L == LFD ) THEN 
               !   ADJ_MASS(NSP) = 1D0
               !ELSE 
               !   ADJ_MASS(NSP) = 0D0
               !ENDIF 

                ! Aerosol mass [kg/box]
               !MASS(NSP) = REAL(CHK_STT(I,J,L,IDTSO4),8)
               !          + REAL(CHK_STT(I,J,L,IDTNIT),8)
               !          + REAL(CHK_STT(I,J,L,IDTNH4),8)
               STT_ADJ(I,J,L,IDTSO4) = STT_ADJ(I,J,L,IDTSO4)
     &                                 + ADJ_MASS(NSP)
               STT_ADJ(I,J,L,IDTNIT) = STT_ADJ(I,J,L,IDTNIT)
     &                                 + ADJ_MASS(NSP)
               STT_ADJ(I,J,L,IDTNH4) = STT_ADJ(I,J,L,IDTNH4)
     &                                 + ADJ_MASS(NSP)

            ENDIF    
         ENDDO

      ENDDO ! L


      !print*, ' *** DISABLE ADJOINT H2O FORCING  '
      !print*, ' *** DISABLE ADJOINT H2O FORCING  '
      !print*, ' *** DISABLE ADJOINT H2O FORCING  '
      !print*, ' *** DISABLE ADJOINT H2O FORCING  '

      ! Return to calling program
      END SUBROUTINE ADJ_CALC_AOD

!------------------------------------------------------------------------------
!
      SUBROUTINE ADJ_GET_NCONC( MASS_DRY_kg, MASS_WET, DENSE_DRY, 
     &                       DENSE_WET,   D,        SIGMA,  BVOL,
     &                       ADJ_D,       ADJ_NCONC, ADJ_MASS_DRY )
!
!******************************************************************************
!  Subroutine ADJ_GET_NCONC computes ADJ_MASS_DRY and updates ADJ_D
!  (dkh, 01/21/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MASS_DRY_kg (REAL*8) : Aerosol dry mass concentration             [kg/box]
!  (2 ) MASS_WET  (REAL*8) : Aerosol water concentration                [ug/m3]
!  (3 ) DENSE_DRY (REAL*8) : Dry aerosol density                        [g/cm3]
!  (4 ) DENSE_WET (REAL*8) : Wet aerosol density                        [g/cm3]
!  (5 ) D         (REAL*8) : Assumed aerosol median diameter            [um]
!  (6 ) SIGMA     (REAL*8) : Assumed aerosol standard deviation         [um]
!  (7 ) BVOL      (REAL*8) : Box volume                                 [m3]
!  (8 ) ADJ_D
!  (9 ) ADJ_NCONC
!  (10) ADJ_MASS_DRY
!     
!  Output:
!  ============================================================================
!  (1 ) ADJ_D  (just for diagnostic)
!  (2 ) ADJ_MASS_DRY
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP

      ! Arguments
      REAL*8,  INTENT(IN)    :: MASS_WET
      REAL*8,  INTENT(IN)    :: MASS_DRY_kg
      REAL*8,  INTENT(IN)    :: DENSE_DRY
      REAL*8,  INTENT(IN)    :: DENSE_WET
      REAL*8,  INTENT(IN)    :: D
      REAL*8,  INTENT(IN)    :: SIGMA
      REAL*8,  INTENT(IN)    :: BVOL

      ! Adjoint arguments
      REAL*8,  INTENT(IN)    :: ADJ_NCONC
      REAL*8,  INTENT(INOUT) :: ADJ_MASS_DRY
      REAL*8,  INTENT(INOUT) :: ADJ_D

      ! Local variables 
      REAL*8                 :: DVM
      REAL*8                 :: DENSITY 
      REAL*8                 :: VOLUME
      REAL*8                 :: MASS_TOTAL
      REAL*8                 :: MASS_DRY
      REAL*8                 :: NCONC

      ! Adjoint local variables 
      REAL*8                 :: ADJ_DVM
      REAL*8                 :: ADJ_DENSITY 
      REAL*8                 :: ADJ_VOLUME
      REAL*8                 :: ADJ_MASS_TOTAL
        
      ! Parameters
      REAL*8, PARAMETER      :: PI = 3.14159265358979324D0

      !=================================================================
      ! ADJ_GET_NCONC begins here!
      !=================================================================

      ! Recalculate first

      ! Convert mass units ( [kg/box] --> [ug/m3] ) 
      MASS_DRY   = MASS_DRY_kg * 1d9 / BVOL

      ! Get total mass
      MASS_TOTAL = MASS_DRY + MASS_WET
   
      IF ( MASS_TOTAL < 1d-20 ) THEN 
         print*, 'housten, we have a problem' 
         NCONC = 0d0  ! <-- this will cause LIDORT to give NANs
         PRINT *, 'Aerosol mass too small: 3'
         RETURN
!         CALL ERROR_STOP( 'Aersol mass too small: 3',
!     &                    'aero_optical_mod.f')
      ENDIF 
 
      ! Average density
      ! Also convert [g/cm3] --> [ug/cm3]
      DENSITY = ( DENSE_DRY * MASS_DRY + DENSE_WET * MASS_WET ) 
     &        / ( MASS_TOTAL ) * 1d6

      ! Total volume [cm3/m3 air]
      VOLUME  = MASS_TOTAL / DENSITY 

      ! Volume mean diameter can be calculated from the log-normal
      ! median diameter and the geometric standard deviation,
      ! see exercise 7.7 of S&P, 1998. 
      ! Also convert units [um] --> [cm]
      DVM = D * EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) * 1D-04

      ! fwd:
      !NCONC = 6d0 * VOLUME / ( PI * DVM ** 3 ) * 1d-06
      ! adj:
      ADJ_DVM    = - 18d0 * VOLUME /  ( PI * DVM ** 4 ) * 1d-06 
     &           * ADJ_NCONC
      ADJ_VOLUME = 6D0 / ( PI * DVM ** 3 ) * 1d-06 
     &           * ADJ_NCONC

      ! fwd:
      !DVM = D * EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) * 1D-04
      ! adj: (ADJ_D is an input with previous value)
      ADJ_D  = ADJ_D 
     &       + ADJ_DVM * EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) * 1D-04

      ! fwd:
      !VOLUME  = MASS_TOTAL / DENSITY 
      ! adj:
      ADJ_MASS_TOTAL = ADJ_VOLUME / DENSITY 
      ADJ_DENSITY    = - ADJ_VOLUME * MASS_TOTAL / DENSITY ** 2

      ! fwd:
      !DENSITY = ( DENSE_DRY * MASS_DRY + DENSE_WET * MASS_WET ) 
      !        / ( MASS_TOTAL ) * 1d6
      ! adj: (ADJ_MASS_DRY is an input with previous value)
      ADJ_MASS_DRY   = ADJ_MASS_DRY 
     &               + DENSE_DRY / MASS_TOTAL * 1D6 * ADJ_DENSITY
      ADJ_MASS_TOTAL = ADJ_MASS_TOTAL 
     &               - DENSITY / MASS_TOTAL * ADJ_DENSITY

      ! fwd:
      !MASS_TOTAL = MASS_DRY + MASS_WET
      ! adj:
      ADJ_MASS_DRY = ADJ_MASS_DRY + ADJ_MASS_TOTAL 
   
      ! fwd:
      !MASS_DRY   = MASS_DRY_kg * 1d9 / BVOL
      ! adj:
      ADJ_MASS_DRY   = ADJ_MASS_DRY * 1d9 / BVOL


      ! Return to calling program
      END SUBROUTINE ADJ_GET_NCONC
!------------------------------------------------------------------------------
      
      SUBROUTINE ADJ_GET_WETD2( MASS_DRY_kg, MASS_WET, DENSE_DRY, 
     &                            DENSE_WET, 
     &                      NCONC,  SIGMA,    BVOL, ADJ_D, 
     &                      ADJ_MASS_WET, ADJ_MASS_DRY, ADJ_NCONC) 
!
!******************************************************************************
!  Subroutine ADJ_GET_WETD2 computes ADJ_MASS_DRY, ADJ_MASS_WET and updates 
!  ADJ_NCONC
!  (dkh, 02/22/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MASS_DRY_kg  (REAL*8) : Aerosol dry mass concentration             [kg/box]
!  (2 ) MASS_WET  (REAL*8) : Aerosol water concentration                [ug/m3]
!  (3 ) DENSE_DRY (REAL*8) : Dry aerosol density                        [g/cm3]
!  (4 ) DENSE_WET (REAL*8) : Wet aerosol density                        [g/cm3]
!  (5 ) NCONC     (REAL*8) : Number concentration                       [#/cm3]
!  (6 ) SIGMA     (REAL*8) : Assumed aerosol standard deviation         [um]
!  (7 ) BVOL      (REAL*8) : Box volume                                 [m3]
!  (8 ) ADJ_D
!  (9 ) ADJ_NCONC
!     
!  Output:
!  ============================================================================
!  (1 ) ADJ_MASS_WET
!  (2 ) ADJ_MASS_WET
!  (3 ) ADJ_NCONC
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP

      ! Arguments
      REAL*8, INTENT(IN)     :: MASS_WET
      REAL*8, INTENT(IN)     :: MASS_DRY_kg
      REAL*8, INTENT(IN)     :: DENSE_DRY
      REAL*8, INTENT(IN)     :: DENSE_WET
      REAL*8, INTENT(IN)     :: NCONC
      REAL*8, INTENT(IN)     :: SIGMA
      REAL*8, INTENT(IN)     :: BVOL

      ! Adjoint arguments 
      REAL*8, INTENT(IN)     :: ADJ_D
      REAL*8, INTENT(OUT)    :: ADJ_MASS_WET
      REAL*8, INTENT(OUT)    :: ADJ_MASS_DRY
      REAL*8, INTENT(INOUT)  :: ADJ_NCONC
        
      ! Local variables 
      REAL*8                 :: D
      REAL*8                 :: DVM
      REAL*8                 :: DENSITY 
      REAL*8                 :: VOLUME
      REAL*8                 :: MASS_TOTAL
      REAL*8                 :: MASS_DRY
 
      ! Local adjiont variables 
      REAL*8                 :: ADJ_DVM
      REAL*8                 :: ADJ_DENSITY 
      REAL*8                 :: ADJ_VOLUME
      REAL*8                 :: ADJ_MASS_TOTAL

      ! Parameters
      REAL*8, PARAMETER      :: PI = 3.14159265358979324D0

      !=================================================================
      ! ADJ_GET_WETD2 begins here!
      !=================================================================

      ! Recalculate VOLUME, MASS_TOTAL, DENSITY 
      ! Convert mass units ( [kg/box] --> [ug/m3] ) 
      MASS_DRY   = MASS_DRY_kg * 1d9 / BVOL

      ! Get total mass
      MASS_TOTAL = MASS_DRY + MASS_WET
   
      IF ( MASS_TOTAL < 1d-20 ) THEN 
         print*, 'housten, we have a problem' 
         !NCONC = 0d0  ! <-- this will cause LIDORT to give NANs
         PRINT *, 'Aerosol mass too small: 4'
         RETURN
!         CALL ERROR_STOP( 'Aersol mass too small 4',
!     &                    'aero_optical_mod.f')
      ENDIF 

      ! Trap for small NCONC, as we multiply by 1/(NCONC^2) below
      IF ( NCONC < 1d0 ) THEN 
         !print*, 'warning: NCONC low '
         ADJ_MASS_WET = 0D0
         ADJ_MASS_DRY = 0D0 
         RETURN
      ENDIF   
    
      ! Average density
      ! Also convert [g/cm3] --> [ug/cm3]
      DENSITY = ( DENSE_DRY * MASS_DRY + DENSE_WET * MASS_WET ) 
     &        / ( MASS_TOTAL ) * 1d6

      ! Total volume [cm3/m3 air]
      VOLUME  = MASS_TOTAL / DENSITY 

      ! Don't need to recalc D or DVM for adjoint
      ! fwd:
      !D = DVM * 1d04 / ( EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) )
      ! adj:
      ADJ_DVM = ADJ_D * 1d04 / ( EXP( 1.5d0 * ( LOG( SIGMA ) )**2 ) )

      ! fwd:
      !DVM = ( 6d0 * VOLUME / ( PI * NCONC ) * 1d-06 ) ** (1.0/3.0)
      ! adj:
      ADJ_VOLUME = 6d0 / ( PI * NCONC ) * 1d-06
     &           * 1.d0 / 3d0 * ( 6d0 * VOLUME / ( PI * NCONC ) * 1d-06) 
     &           ** ( - 2.0/3.0) * ADJ_DVM
      ADJ_NCONC  = ( - 6d0 * VOLUME / ( PI * NCONC ** 2 ) * 1d-06 )
     &           * 1.d0 / 3d0 * ( 6d0 * VOLUME / ( PI * NCONC ) * 1d-06)
     &           ** ( - 2.0/3.0) * ADJ_DVM + ADJ_NCONC

      ! fwd:
      !VOLUME  = MASS_TOTAL / DENSITY 
      ! adj:
      ADJ_MASS_TOTAL = ADJ_VOLUME / DENSITY 
      ADJ_DENSITY    = - VOLUME / DENSITY * ADJ_VOLUME

      ! fwd:
      !DENSITY = ( DENSE_DRY * MASS_DRY + DENSE_WET * MASS_WET ) 
      !        / ( MASS_TOTAL ) * 1d6
      ! adj:
      ADJ_MASS_DRY   = DENSE_DRY / MASS_TOTAL * 1D6 * ADJ_DENSITY 
      ADJ_MASS_WET   = DENSE_WET / MASS_TOTAL * 1D6 * ADJ_DENSITY 
      ADJ_MASS_TOTAL = ADJ_MASS_TOTAL 
     &               - DENSITY / MASS_TOTAL * ADJ_DENSITY 

      ! fwd:
      !MASS_TOTAL = MASS_DRY + MASS_WET
      ! adj:
      ADJ_MASS_DRY = ADJ_MASS_DRY + ADJ_MASS_TOTAL 
      ADJ_MASS_WET = ADJ_MASS_WET + ADJ_MASS_TOTAL 
   
      ! fwd:
      !MASS_DRY   = MASS_DRY_kg * 1d9 / BVOL
      ! adj:
      ADJ_MASS_DRY = ADJ_MASS_DRY * 1D9 / BVOL


      ! Return to calling program
      END SUBROUTINE ADJ_GET_WETD2

!------------------------------------------------------------------------------

      SUBROUTINE ADJ_MIE_LOOKUP( DIAM, W, NSP,
     &                           NCONC, 
     &                           ADJ_NCONC, ADJ_DIAM, ADJ_BEXT, ADJ_SSA, 
     &                           ADJ_PHM,I,J,L )
!
!******************************************************************************
!  Subroutine ADJ_MIE_LOOKUP 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DIAM
!  (2 ) W
!  (3 ) NSP
!  (4 ) NCONC
!  (5 ) ADJ_BEXT   
!  (6 ) ADJ_SSA
!  (7 ) ADJ_PHM
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) ADJ_NCONC
!  (2 ) ADJ_DIAM
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE MIE_MOD,       ONLY : MIE_TABLE
      USE MIE_MOD,       ONLY : MIE_TABLE_PHM
      USE MIE_MOD,       ONLY : MIE_TABLE_ADJ

      ! dkh debug
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, LFD 

      ! Arguments
      REAL*8,  INTENT(IN)    :: DIAM
      INTEGER, INTENT(IN)    :: W
      INTEGER, INTENT(IN)    :: NSP
      INTEGER, INTENT(IN)    :: L, I, J
      REAL*8,  INTENT(IN)    :: NCONC 
      REAL*8,  INTENT(INOUT) :: ADJ_NCONC 
      REAL*8,  INTENT(INOUT) :: ADJ_DIAM
      REAL*8,  INTENT(IN)    :: ADJ_BEXT
      REAL*8,  INTENT(IN)    :: ADJ_SSA
      REAL*8,  INTENT(IN)    :: ADJ_PHM(NPM_MAX)

      ! Local variables 
      REAL*8                 :: N, BEXTT, RAD
      REAL*8                 :: ADJ_RAD
      REAL*8                 :: PHM_1(NPM_ACT), PHM_2(NPM_ACT)
      REAL*8                 :: R1, R2
      REAL*8                 :: DPHM_DR(NPM_ACT)
      INTEGER                :: N1, N2
      REAL*8                 :: SSA0
      REAL*8                 :: BEXT0
      REAL*8                 :: ADJ_BEXT0
      REAL*8                 :: ADJ_BEXTT
      REAL*8                 :: ADJ_SSA0
      REAL*8                 :: TEMP
   
      ! Parameters
      !REAL*8, PARAMETER      :: RNG = LOG(RMAX/RMIN) / REAL(NRA_MAX)
      REAL*8                 :: RNG 

      !=================================================================
      ! ADJ_MIE_LOOKUP begins here!
      !=================================================================



      ! BUG FIX: need to calculated this separately for BC and SULF, 
      ! as they have different size ranges and spacing in the lookup
      ! table. (dkh, 09/27/10) 
      RNG = LOG(RMAX(NSP)/RMIN(NSP)) / ( REAL(NRA_MAX) - 1d0 )

      RAD = DIAM / 2D0 

      ! Get index in MIE_TABLE that corresponds to current diameter.
      ! Diameters are calculated at NRA_MAX - 1 points between Dmin and 
      ! Dmax with a spacing defined in ln space: 
      !      d_n = d_min * exp( ( n - 1 ) * ln(d_max/d_min) / npts )
      ! The inverse of this formulat gives us the n for any d_n. 
      N   = LOG(RAD/RMIN(NSP)) / RNG + 1

      ! Enforce bounds on this index
      IF ( N < 1.0           ) N = 1.d0
      IF ( N > REAL(NRA_MAX) ) N = REAL(NRA_MAX)

      ! recalculate BEXTT and BEXT0
      BEXTT = MIE_TABLE(W,NSP,NOP_BEXT,INT(N))
      BEXT0 = BEXTT * NCONC *1d-6
     
      ! recalculate SSA 
      SSA0  = MIE_TABLE(W,NSP,NOP_SSA,INT(N))

      ! Account for absorption enhancement 
      ! fwd:
      !IF ( NSP == NSP_BC ) THEN
      !   BEXT = BEXT0 * ( SSA + ABS_FAC * ( 1d0 - SSA ) )
      !   SSA  = SSA0  / ( SSA + ABS_FAC * ( 1d0 - SSA ) )
      !ENDIF
      ! adj:
      IF ( NSP == NSP_BC ) THEN
 
         ! define this term to make equations simple
         TEMP = SSA0 + ABS_FAC * ( 1d0 - SSA0 )

         !   SSA  = SSA0  / ( SSA0 + ABS_FAC * ( 1d0 - SSA0 ) )
         ADJ_SSA0   = ADJ_SSA * ABS_FAC / ( TEMP ** 2 )
         
         !   BEXT = BEXT0 * ( SSA0 + ABS_FAC * ( 1d0 - SSA0 ) )
         ADJ_SSA0   = ADJ_SSA0
     &              + BEXT0 * ( 1 - ABS_FAC )
     &              * ADJ_BEXT
         ADJ_BEXT0  = ADJ_BEXT * TEMP 

      ELSE 
         ADJ_SSA0  = ADJ_SSA
         ADJ_BEXT0 = ADJ_BEXT
      ENDIF 

      ! Contribution from ADJ_PHM
      ! Use values from table to get the derivative of PHM w.r.t. R 
      ! using finite differencing. 
      IF ( INT(N+1) > NRA_MAX) THEN 
         N1 = INT(N-1)
         N2 = INT(N)
      ELSEIF (INT(N-1) < 1 ) THEN 
         N1 = INT(N)
         N2 = INT(N+1)
      ELSE
         N1 = INT(N-1)
         N2 = INT(N+1)
      ENDIF 

      PHM_1(1:NPM_ACT) = MIE_TABLE_PHM(W,NSP,N1,1:NPM_ACT)
      PHM_2(1:NPM_ACT) = MIE_TABLE_PHM(W,NSP,N2,1:NPM_ACT)
      R1               = RMIN(NSP) * EXP( RNG ) ** N1
      R2               = RMIN(NSP) * EXP( RNG ) ** N2
      DPHM_DR(:)       = ( PHM_2(:) - PHM_1(:) ) / ( R2 - R1 )
      ADJ_RAD          = SUM(DPHM_DR(:) * ADJ_PHM(1:NPM_ACT) )

      ! dkh debug 
      !print*, ' DIAM, N, BEXTT, BEXT = ', DIAM, N, BEXTT, BEXT

      ! Contribution from ADJ_SSA
      ! fwd:
      !IF ( NSP == NSP_SULF ) THEN
      !   SSA  = 1d0  ! for sulfate
      !ELSE
      !   SSA  = MIE_TABLE(W,NSP,NOP_SSA,INT(N))
      !ENDIF 
      ! adj:
      IF ( NSP == NSP_SULF ) THEN
         ADJ_RAD = ADJ_RAD + 0D0  
      ELSE
 
         !-----------------------------------------------------------------
         ! BUG_FIX: use FD sensitivities instead of table.  (dkh, 03/27/11) 
         ! BUG_FIX: use table instead of FD (dkh, 07/01/11) 
         ADJ_RAD = ADJ_RAD + MIE_TABLE_ADJ(W,NSP,NOP_SSA,INT(N)) 
     &           * ADJ_SSA0
         ! OLD: 
         !ADJ_RAD = ADJ_RAD
         !  + (MIE_TABLE(W,NSP,NOP_SSA,N2) - MIE_TABLE(W,NSP,NOP_SSA,N1))
         !  / ( R2 - R1 ) * ADJ_SSA0
         !-----------------------------------------------------------------

      ENDIF 

 
      ! fwd:
      !BEXT0 = BEXTT * NCONC * 1d-6
      ! adj:
      ADJ_BEXTT = NCONC * ADJ_BEXT0 * 1d-6
      ADJ_NCONC = ADJ_NCONC + BEXTT * ADJ_BEXT0 * 1d-6

      ! Contribution from ADJ_BEXT
      ! fwd: 
      ! BEXTT = MIE_TABKE(W,NSP,NOP_BEXT,INT(N))
      ! adj:
      !-----------------------------------------------------------------
      ! BUG_FIX: use FD sensitivities instead of table.  (dkh, 03/27/11) 
      ! BUG_FIX: use table instead of FD. (dkh, 07/01/11) 
      ADJ_RAD = ADJ_RAD + MIE_TABLE_ADJ(W,NSP,NOP_BEXT,INT(N)) 
     &        * ADJ_BEXTT
      !ADJ_RAD = ADJ_RAD 
      !   + (MIE_TABLE(W,NSP,NOP_BEXT,N2) - MIE_TABLE(W,NSP,NOP_BEXT,N1))
      !   / ( R2 - R1 ) * ADJ_BEXTT
      !-----------------------------------------------------------------


      ! fwd:
      !RAD = DIAM / 2D0 
      ! adj:
      ADJ_DIAM = ADJ_DIAM + ADJ_RAD / 2D0

      ! dkh debug
      !ADJ_DIAM = 0d0
      ! this works.  somehow there is a connection between ADJ_DIAM and
      ! ADJ_STT in this situation, but not a connection between STT and DIAM!

      ! dkh debug
      !IF ( W == 3 .and. L == LFD .and. NSP == 2 .and. 
      !    I == IFD .and. J == JFD ) THEN 
      !   print*, ' adj: BEXTT = ', BEXTT
      !   print*, ' adj: BEXT0 = ', BEXT0
      !   print*, ' adj: NCONC = ', NCONC
      !   print*, ' adj: SSA0  = ', SSA0
      !   print*, ' adj: N     = ', N, INT(N)
      !
      !   print*, ' adj: ADJ_BEXT  = ', ADJ_BEXT
      !   print*, ' adj: ADJ_BEXT0 = ', ADJ_BEXT0
      !   print*, ' adj: ADJ_BEXTT = ', ADJ_BEXTT
      !
      !   print*, ' adj: ADJ_NCONC = ', ADJ_NCONC
      !   print*, ' adj: ADJ_DIAM  = ', ADJ_DIAM
      !
      !ENDIF 

      ! Return to calling program
      END SUBROUTINE ADJ_MIE_LOOKUP

!------------------------------------------------------------------------------

      SUBROUTINE PREPARE_PHYSICS_GC_IJ( I, J, PRESSURE_GRID, 
     &  TEMPERATURE_GRID, HEIGHT_GRID, AIRDENS, GAS_PROFILE, BEAM_SZAS )
!
!******************************************************************************
!  Subroutine PREPARE_PHYSICS_GC_IJ sets lidort inputs that depend only 
!  on lat and lon. (dkh, 03/03/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : Lon                                  [unit]
!  (2 ) J (INTEGER) : Lat                                  [unit]
!     
!  Variables as Output:
!  ============================================================================
!  (2 ) PRESSURE_GRID
!  (3 ) AIRDENS
!  (4 ) TEMPERATURE_GRID
!  (5 ) GAS_PROFILE
!  (6 ) NLAYERS   
!  (7 ) NMOMENTS_INPUT	
!  (8 ) HEIGHT_GRID
!
!  NOTES:
!  (1 ) This is based upon the jactest_3p3_solar2.f driver provided with 
!       LIDORT by R. Spurr. Original code is lower case. 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD, ONLY : O3_PROF_SAV
      USE COMODE_MOD,   ONLY : CSPEC,     JLOP
      USE DAO_MOD,      ONLY : AIRDEN,    TS,      T
      USE DAO_MOD,      ONLY : SUNCOS,    BXHEIGHT
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TRACERID_MOD, ONLY : IDO3

#     include "cmn_fj.h"     ! FJ sizes, need for jv_cmn, CMN_SIZE
#     include "jv_cmn.h"     ! DO3
#     include "CMN_GCTM"     ! Rdg0

      ! Arguments
      INTEGER               :: I, J
      REAL(KIND=8)          :: HEIGHT_GRID ( 0:MAXLAYERS )
      REAL(KIND=8)          :: PRESSURE_GRID    ( 0:MAXLAYERS )
      REAL(KIND=8)          :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      REAL(KIND=8)          :: AIRDENS          ( 1:MAXLAYERS )
      REAL(KIND=8)          :: GAS_PROFILE      ( 1:MAXLAYERS )
      REAL(KIND=8)          :: BEAM_SZAS ( MAXBEAMS )

    
      ! Local variables 
      INTEGER               :: L, N, JLOOP, IJLOOP
      REAL*8                :: SS_ZLEVELS( 0: LLPAR )

      ! XNUMOL_AIR: (molecules air / kg air) 
      REAL*8,  PARAMETER    :: XNUMOLAIR  = 6.022d+23 / 0.02897d0
      
      !=================================================================
      ! PREPARE_PHYSICS_GC_IJ begins here!
      !=================================================================

      ! Get model top pressure [hPa].  Can't access directly from CMN_SIZE, 
      ! since we have to use the cmn_fj.h rather than CMN_SIZE. 
      PRESSURE_GRID(0) = GET_PEDGE(I,J,LLPAR+1)

      ! Get layer input.  
      DO L = 1, LLPAR

         ! N starts with 0=TOA.  L = 1 is the ground. 
         N = LLPAR - L + 1 

         ! pressures [hPa] (not used)
         PRESSURE_GRID(N)    = GET_PEDGE(I,J,L)
 
         IF ( L == 1 ) THEN
            SS_ZLEVELS(N) = 0d0 !or altitude?

            ! temperature at surface [k]
            TEMPERATURE_GRID(N) = TS(I,J)
 
         ELSE 
 
            ! temperature at lower edge [k] ( T is mid-layere temp).
            TEMPERATURE_GRID(N) = 0.5d0 * ( T(I,J,L) + T(I,J,L-1) )
 
            ! heights of lower layer boundaries [km]
            SS_ZLEVELS(N)  = SS_ZLEVELS(N+1)
     &                     + Rdg0 * 1d-3 * T(I,J,L-1) 
     &                     * LOG( PRESSURE_GRID(N+1) 
     &                     / PRESSURE_GRID(N) )
         ENDIF 

         ! air density  [molec/cm3], convert from [kg/m3].
         ! NOTE:  AIRDEN index is (L, I, J), not the usual (I, J, L)
         AIRDENS(N)     = AIRDEN(L,I,J) * XNUMOLAIR / 1d6 

         ! Gas VMR [DU]
         IF ( L <= LLTROP ) THEN 

            ! Get 1-D array index for extracting O3 from CSPEC
            JLOOP = JLOP(I,J,L) 

            ! Get O3 from CSPEC within the troposphere
            ! and convert from [#/cm3] --> [DU]
            IF ( JLOOP > 0 )
     &      GAS_PROFILE(N) = CSPEC(JLOOP,IDO3) 
     &                     * BXHEIGHT(I,J,L) * 100d0
     &                     * 3.7219d-17

         ! Otherwise get O3 from FAST-J (scaled climatology)
         ELSE

            ! Convert from [#/cm2] --> [DU]
            GAS_PROFILE(N) = O3_PROF_SAV(I,J,L) * 3.7219d-17

         ENDIF 

         ! dkh debug 
         IF ( I == 61 .and. J == 32 ) THEN 
            print*, 'GAS_PROF = ', GAS_PROFILE(N), L, N
         ENDIF 

      ENDDO

      ! top of atmosphere in [km]
      SS_ZLEVELS(0) =  SS_ZLEVELS(1)
     &              + Rdg0 * 1d-3 * T(I,J,LLPAR)
     &              * LOG( PRESSURE_GRID(1) / GET_PEDGE(I,J,LLPAR+1) )


      ! Set LIDORT number of layers and Legendre moments
      nlayers        = LLPAR
      nmoments_input = NPM_MAX

      ! Set LIDORT height grid for pseudo-spherical treatment
      DO n = 0, nlayers
        height_grid(n) = ss_zlevels(n)
      END DO

      ! note:  probably don't need these error checking and 
      ! calculation of X0 and SX0 -- it's redundant (dkh, 01/19/10) 
      ! Set solar zenith angle and associated quantities
      ! note: this overwrites the values in the lidort input file *.inp
      IJLOOP = ( J - 1 ) * IIPAR + I 
      BEAM_SZAS(1) = ACOS( SUNCOS(IJLOOP) )
!      IF ( BEAM_SZAS(1) .LT. ZERO   .OR.
!     &     BEAM_SZAS(1) .GE. 90.0D0      ) THEN
!          CALL ERROR_STOP('LIDORT: bad input: out-of-range beam',
!     &                    'lidort_mod.f'                   )
!      ENDIF
!      X0(1)  = DCOS( BEAM_SZAS(1) )
!      SX0(1) = DSQRT( ONE - X0(1) * X0(1) )
      
      ! dkh debug
      IF ( L_PRINT ) THEN 
         print*, 'test at cell i61, j32 with zenith angle = ',
     &    BEAM_SZAS(1)
      ENDIF 
 
      ! Return to calling program
      END SUBROUTINE PREPARE_PHYSICS_GC_IJ

!------------------------------------------------------------------------------

      SUBROUTINE PREPARE_PHYSICS_GC_IJW
     & ( TR,                   epsfac,              npert, 
     &   status,               AIRDENS,             GAS_PROFILE, 
     &   LAMBERTIAN_ALBEDO,    OMEGA_TOTAL_INPUT,   DELTAU_VERT_INPUT, 
     &   PHASMOMS_TOTAL_INPUT, L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, 
     &   L_PHASMOMS_TOTAL_INPUT,    I,    J,      W                    )

!
!******************************************************************************
!  Subroutine PREPARE_PHYSICS_GC_IJW sets physics that depends upon lat, lon
!  and wavelength.  (dkh, 03/03/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [unit]
!     
!  NOTES:
!  (1 ) This is based upon the jactest_3p3_solar2.f driver provided with 
!       LIDORT by R. Spurr. Original code is lower case. 
!  (2 ) Updated to 3p5 (dkh, 07/29/10) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE UVALBEDO_MOD,    ONLY : UVALBEDO
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD, NFD 

      !  --------------------
      !  LIDORT Include files
      !  --------------------
#     include "CMN_SIZE"  ! LLPAR

      ! arguments
      LOGICAL          status
      INTEGER               ::   TR
      INTEGER               :: npert
      DOUBLE PRECISION      :: epsfac 
      REAL(KIND=8)          :: AIRDENS          ( 1:MAXLAYERS )
      REAL(KIND=8)          :: GAS_PROFILE      ( 1:MAXLAYERS )

      ! Added for GC implementation
      INTEGER               :: I, J, W
  
      !  -----------------------
      !  LIDORT inputs.  Local, single threaded version.
      !  -----------------------
      REAL(KIND=8) :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )
      REAL(KIND=8) :: DELTAU_VERT_INPUT  ( MAXLAYERS )

      !  Phase function Legendre-polynomial expansion coefficients
      !   Include all that you require for exact single scatter calculations
      REAL(KIND=8) :: PHASMOMS_TOTAL_INPUT
     &     ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

      !  Optical property linearizations
      !  Layer linearization (bulk property variation) input
      !  Layer linearization (phase function variation) input
      REAL(KIND=8)    
     &   :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(KIND=8)    
     &   :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(KIND=8)    
     &   :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,
     &      MAXLAYERS)

      !  Lambertian Surface control
      REAL(KIND=8) :: LAMBERTIAN_ALBEDO 

      ! Local physics storage
      ! ---------------------
      double precision gasvods(NWL_MAX)
      double precision lambdas(NWL_MAX)
      double precision rayleigh_xsec(NWL_MAX)
      double precision rayleigh_depol(NWL_MAX)
      double precision local_albedo(NWL_MAX)
      double precision gas_xsec(NWL_MAX,LLPAR)
      double precision aerosol_deltaus(NWL_MAX,LLPAR)
      double precision aerosol_ssalbs(NWL_MAX,LLPAR)
      double precision 
     &   aerosol_phasmoms(NWL_MAX,LLPAR,0:NPM_MAX)

      ! other local variables
      INTEGER          n, l, nd, ut, uta, k, ll
      double precision aermom,raymom,depol,rayxs,tempdummy,agas,dummy
      double precision scatmom_1, scatmom_2, scatmom_l, deltau, taue_aer
      double precision tg,taua_gas,taue_mol,taus_ray,taus_aer,taus_total
      double precision du_to_cm2, baer, waer, vaer, saer

      double precision gascc, taua_gas0, taua_gas_eps, nx, sum1, sum0

      !----------------------------------------------------------------
      ! Data
      !----------------------------------------------------------------

      ! MOVE TO MODULE HEADER 
      ! O3 absorption cross section [unit]
      ! Data is take from
      !  315nm: Burrows et al., 1999, Fig 4
      !  357nm: Voigt   et al., 2001, Fig 5 (@ 246 K)
      !  500nm: Voigt   et al., 2001, Fig 2 
      !  833nm: Voigt   et al., 2001, Fig 2 (ave)
      ! 1667nm: Bogumil et al., 2003, Fig 1 
      REAL*8              :: O3_XSEC_TABLE(NWL_MAX)
      DATA                   O3_XSEC_TABLE  /  4d-20,      7d-23, 
     &                                         1d-21,      7d-23,
     &                                         0d0  / 

      ! Rayleigh
      ! from Bodhaine et al., 1999, Table 3
      ! Linearly interpolate. Values at 1667 nm are just a guess.
      REAL*8              :: RAY_XSEC_TABLE(NWL_MAX)
      DATA                   RAY_XSEC_TABLE  / 4.5830d-26,  2.700d-26,
     &                                         6.6614d-27,  8.380d-27,
     &                                         0d0  / 
      REAL*8              :: RAY_DEPOL_TABLE(NWL_MAX)
      DATA                   RAY_DEPOL_TABLE / 1.05521d0,   1.05280d0,
     &                                         1.04935d0,   1.04750d0, 
     &                                         1d0  /

      !=================================================================
      ! PREPARE_PHYSICS_GC_IJW begins here!
      !=================================================================

      !  Initialise
      taua_gas_eps = 0.0d0

      !  Conversion, Dobson units to mol/cm/cm
      du_to_cm2 = 2.68668d+16

      ! Wavelengths [nm]
      LAMBDAS(W)        = LAMBDA_TABLE(W)

      !NOTE: the three quantities below are wavelength dependent at 
      ! the very least. 

      ! Local albedo 
      ! from Martin et al., 2004
      ! depnds on season.  should also depend on location?
      ! The GC albedo used for FAST-J depends on season and
      ! location, but is for 340-380 nm
      !LOCAL_ALBEDO(W)   = UVALBEDO(I,J)
      LOCAL_ALBEDO(W)   = TEMIS_LER(I,J,W)

      ! These two also depend slightly on CO2
      ! Rayleigh X-section
      RAYLEIGH_XSEC(W)  = RAY_XSEC_TABLE(W) ! 6.6614d-27
 
      ! Rayleigh depolarization ratio
      RAYLEIGH_DEPOL(W) = RAY_DEPOL_TABLE(W) ! 1.04935

      ! Define the following properties:
      !  Trace gas absorption cross-section, in cm^2/mol.
      !  aerosol extinction optical depth,
      !  aerosol single scattering albedo,
      !  aerosol phase function Legendre expansion coeffs..

      ! Vertical loop.  L is 1 at surface, N is 1 at top. 
      DO N = 1 , LLPAR
           
         L = LLPAR - N + 1 
 
         ! Trace gas absorption cross-section [cm2/mol]
         GAS_XSEC(W,N) =  O3_XSEC_TABLE(W) 
       
         ! Aerosol extinction optical depth 
         AEROSOL_DELTAUS(W,N) = LAYER_AOD(I,J,L,W)

         ! Aerosol single scattering albedoes 
         AEROSOL_SSALBS(W,N)  = LAYER_SSA(I,J,L,W)

         ! Aerosol phase function Legendre expansion coeffs
         AEROSOL_PHASMOMS(W,N,1:NPM_MAX) = 
     &         LAYER_PHM(I,J,L,1:NPM_MAX,W)


      ENDDO

      ! The zeroeth order PHM is always 1.d0
      !note: this is redundant, as the lidort zeroeth order input is hardwired to 1d0
      ! below
      AEROSOL_PHASMOMS(:,:,0) = 1d0

      ! dkh debug
      !IF ( I == 61 .AND. J == 32 .and. W == 3 ) THEN   
      IF ( I == IFD .AND. J == JFD .and. W == 3 ) THEN   
         WRITE(*,'(3i5)'), LLPAR, NPM_MAX, 1
!         WRITE(*,'(1pe15.7)') SS_ZLEVELS(0)
!         WRITE(*,'(1p5e15.7)') (SS_ZLEVELS(N),PRESSURE(N),
!     &       TEMPERATURE(N), AIRDENS(N), 
!     &       GAS_PROFILE(N),N=1,LLPAR)
         WRITE(*,'(2f12.5,1p2e15.6)') LAMBDAS(NWL_500), 
     &      LOCAL_ALBEDO(NWL_500), 
     &      RAYLEIGH_XSEC(NWL_500), RAYLEIGH_DEPOL(NWL_500)
         WRITE(*,'(1p24e15.5)') (GAS_XSEC(NWL_500,N), 
     &      AEROSOL_DELTAUS(NWL_500,N), 
     &      AEROSOL_SSALBS(NWL_500,N), 
     &      AEROSOL_PHASMOMS(NWL_500,N,0:NPM_MAX),
     &      n=1,LLPAR)
         print*, 'ok done' , NPM_MAX
      ENDIF 
 
      ! Surface albedo. (This may be overwritten)
      !  Will only be used if the Lambertian surface flag is set
       lambertian_albedo = local_albedo(w)
       IF ( TR == 5 ) lambertian_albedo = lambertian_albedo * epsfac
      
       ! Rayleigh: Depolarization ratio and second phase function moment
       depol  = rayleigh_depol(w)
       rayxs  = rayleigh_xsec(w)
       raymom  = ( 1.0d0 - depol ) / ( 2.0d0 + depol )

       ! initialise total atmsopheric optical depth for gas
       gasvods(w) = 0.0

       ! start layer loop
       DO n = 1, LLPAR

        ! optical depth due to Rayleigh scattering
        taus_ray = rayxs * airdens(n)

        ! Trace gas layer column amount
        taua_gas = 0.0d0
        gascc = gas_profile(n)

        ! Trace gas absoprtion optical thickness
        tg = gascc * gas_xsec(w,n)
        tg = du_to_cm2 * tg
        taua_gas   = taua_gas + tg
        gasvods(w) = gasvods(w) + tg

        !  Make a finite difference perturbation if flagged
        if ( n == npert .and. TR == 2 ) then 
            taua_gas = taua_gas * epsfac
        endif 



!        ! debug code (ignore)
!        taua_gas0 =  gas_profile(n)*du_to_cm2*gas_xsec(w,n)
!        if (n.eq.npert.and.wq.eq.1)taua_gas_eps = taua_gas - taua_gas0

        ! Trace gas absoprtion + Rayleigh scattering optical depth
        taue_mol = taus_ray + taua_gas

        !  aerosol optical thickness for extinction
        taue_aer = aerosol_deltaus(w,n)

        !   Make a finite difference perturbation if flagged
        if ( n == npert .and. TR == 3 ) taue_aer = taue_aer * epsfac
        if ( n == npert .and. TR == 4 ) taue_aer = taue_aer * epsfac

        ! total optical thickness for layer. Set the LIDORT input
        deltau   = taue_aer + taue_mol

!        if ( n == 21 )print*, 'deltau, taue_aer, taue_mol',
!     &        deltau, taue_aer, taue_mol,TR

        deltau_vert_input(n) = deltau

        ! Derivative factors for computing linearized inputs
        agas = taua_gas / deltau
        vaer = taue_aer / deltau

        ! aerosol scattering optical depth
        taus_aer     = aerosol_ssalbs(w,n) * taue_aer

        ! total scattering optical depth
        taus_total   = taus_ray + taus_aer

        ! total single scattering albedo. Set the LIDORT input
        omega_total_input(n)  = taus_total / deltau

        ! Toggle the SS albedo if near to one
        !  (actually not required with this data set)
        IF (omega_total_input(n) .gt. 0.9999 ) THEN
           omega_total_input(n) = 0.9999
        ENDIF

        ! some additional ratios for computing the linearized inputs
        waer = aerosol_ssalbs(w,n) / omega_total_input(n)
        saer = taus_aer / taus_total

        ! phase function: Zeroth moment  = 1
        phasmoms_total_input(0,n) = 1.0D0

        ! phase function first moment
        scatmom_1 = taus_aer * aerosol_phasmoms(w,n,1)
        phasmoms_total_input(1,n) = scatmom_1 / taus_total

        ! phase function second moment (include Rayleigh):
        aermom    = aerosol_phasmoms(w,n,2)
        scatmom_2 = taus_ray * raymom + taus_aer * aermom
        phasmoms_total_input(2,n) = scatmom_2 / taus_total

        ! phase function moments > 2:
        DO l = 3, NPM_MAX
           aermom              = aerosol_phasmoms(w,n,l)
           scatmom_l           = taus_aer * aermom
           phasmoms_total_input(l,n) = scatmom_l / taus_total
        ENDDO

        ! Linearization inputs if a Linearization Calculation is required
        !if ( dowf ) then
         layer_vary_flag(n)   = .true.
         ! (dkh, 04/04/08) 
         !layer_vary_number(n) = 2
         layer_vary_number(n) = MAX_ATMOSWFS
         l_deltau_vert_input(1,n) = + agas
         l_omega_total_input(1,n) = - agas
         l_deltau_vert_input(2,n) = + vaer
         l_omega_total_input(2,n) = + vaer * ( waer - one )
         l_phasmoms_total_input(1,0,n) = zero
         l_phasmoms_total_input(2,0,n) = zero
         do l = 1, NPM_MAX
          l_phasmoms_total_input(1,l,n) = zero
          baer = aerosol_phasmoms(w,n,l)/phasmoms_total_input(l,n)
          l_phasmoms_total_input(2,l,n) = saer * ( baer - one )
         enddo
        !endif

         ! (dkh, 04/03/08) 
         ! w.r.t. aerosol_ssalb
         l_deltau_vert_input(3,n) = zero
         l_omega_total_input(3,n) = + saer
         l_phasmoms_total_input(3,0,n) = zero
         do l = 1, NPM_MAX
          baer = aerosol_phasmoms(w,n,l)/phasmoms_total_input(l,n)
          l_phasmoms_total_input(3,l,n) = saer * ( baer - one )
         enddo

         ! (dkh, 04/04/08) 
         ! w.r.t. aerosol_phasmoms
         !l_deltau_vert_input(4,n) = zero
         !l_omega_total_input(4,n) = zero
         !l_phasmoms_total_input(4,0,n) = zero
         !do l = 1, NPM_MAX
         ! baer = aerosol_phasmoms(w,n,l)/phasmoms_total_input(l,n)
         ! l_phasmoms_total_input(4,l,n) = saer * baer
         !enddo
         l_deltau_vert_input(4:MAX_ATMOSWFS,n) = zero
         l_omega_total_input(4:MAX_ATMOSWFS,n) = zero
         l_phasmoms_total_input(4:MAX_ATMOSWFS,:,n) = zero
         do ll = 4, MAX_ATMOSWFS
          l = ll - 3
          baer = aerosol_phasmoms(w,n,l)/phasmoms_total_input(l,n)
          l_phasmoms_total_input(ll,l,n) = saer * baer
         enddo


!        ! No linearization required, just zero the inputs
!        if ( .not. dowf ) then
!         layer_vary_flag(n)   = .false.
!         layer_vary_number(n) = 0
!         l_deltau_vert_input(1,n,TR) = zero
!         l_omega_total_input(1,n,TR) = zero
!         l_deltau_vert_input(2,n,TR) = zero
!         l_omega_total_input(2,n,TR) = zero
!         ! (dkh, 04/04/08) 
!         l_deltau_vert_input(3:MAX_ATMOSWFS,n,TR) = zero
!         l_omega_total_input(3:MAX_ATMOSWFS,n,TR) = zero
!         do l = 0, NPM_MAX
!          l_phasmoms_total_input(1,l,n,TR) = zero
!          l_phasmoms_total_input(2,l,n,TR) = zero
!          ! (dkh, 04/04/08) 
!          l_phasmoms_total_input(3:MAX_ATMOSWFS,l,n,TR) = zero
!         enddo
!        endif

       ! End layer loop
       END DO

!!  thermal BB input.
!C   Read pre-prepared BB functions 
!
!      OPEN(45,file='testbb.dat',status='old')
!      read(45,*)dummy,dummy
!      DO n = 0, nlayers
!       READ(45,'(i2,f9.3,1pe18.8)')nd,tempdummy,thermal_bb_input(n)
!      END DO
!      READ(45,'(2x,f9.3,1pe18.8)')tempdummy,surfbb
!      CLOSE(45)

      ! normal return for this routine
      status = .false.
      RETURN

      ! error return for this routine
400   CONTINUE
      status = .TRUE.
      RETURN

      ! Return to calling program
      END SUBROUTINE PREPARE_PHYSICS_GC_IJW

!----------------------------------------------------------------------

      SUBROUTINE MAKE_AOD_FILE( YYYYMMDD, HHMMSS, N_CALC )
!
!******************************************************************************
!  Subroutine MAKE_AOD_FILE saves average AOD (dkh, 01/24/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : Date of average                          [unit]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE DAO_MOD,     ONLY : SUNCOS
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : GET_TAU,     EXPAND_DATE

#     include "CMN_SIZE"    ! Size params

      ! Arguments
      INTEGER              :: YYYYMMDD
      INTEGER              :: HHMMSS
      INTEGER, INTENT(IN)  :: N_CALC

      ! Local variables 
      INTEGER              :: I, J, I0, J0, L, W
      CHARACTER(LEN=120)   :: FILENAME
      REAL*4               :: DAT(IIPAR,JJPAR,LLPAR)
      REAL*4               :: TEMP

      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      !=================================================================
      ! MAKE_AOD_FILE begins here!
      !=================================================================
      DAT(:,:,:) = 0d0 

      FILENAME = TRIM( 'aod.YYYYMMDD.hhmm.NN' )

      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS)

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Append data directory prefix
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'Average AOD     data file '
      CATEGORY = 'IJ-AOD-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UNIT     = 'none'

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_AOD_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_AOD, FILENAME, TITLE )

      !--------------------------------------
      ! write AOD
      !--------------------------------------

      ! Loop over wavelengths
      DO W = 1, NWL_MAX

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )    
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Transfer data to 3D, REAL*4 array
            DAT(I,J,1) = REAL(AOD_SAVE(I,J,W),4)
     
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_AOD,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  W,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         DAT(:,:,1)        )

      ENDDO

      ! Define variables for BINARY PUNCH FILE OUTPUT
      CATEGORY = 'IJ-AOD-$'

      !--------------------------------------
      ! write reflectence
      !--------------------------------------
      UNIT     = '%'

      ! Loop over wavelengths
      DO W = 1, NWL_MAX

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )    
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Transfer data to 3D, REAL*4 array
            DAT(I,J,1) = REAL(GLOB_FLUX_W(W,I,J),4)

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_AOD,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  NWL_MAX+W,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         DAT(:,:,1)        )


      ENDDO 

      !--------------------------------------
      ! write reflectance sensitivities 
      !--------------------------------------
      UNIT     = 'none'

      ! Loop over wavelengths
      DO W = 1, NWL_MAX

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)    
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Transfer data to 3D, REAL*4 array
            DAT(I,J,L) = REAL(JAC_GLOB_FLUX_W(2,L,W,I,J),4)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_AOD,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  2*NWL_MAX+W,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,
     &               J0+1,      1,         DAT               )


      ENDDO  ! W

      !--------------------------------------
      ! write radiative flux
      !--------------------------------------
      UNIT     = 'W/m2'


!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, TEMP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Average
         !IF ( GLOB_FLUX_COUNT > 0 ) THEN 
         !   TEMP    = GLOB_FLUX(I,J) / REAL(GLOB_FLUX_COUNT)
         !ELSE 
         !   TEMP    = 0d0
         !ENDIF 

         ! Transfer data to 3D, REAL*4 array
         !DAT(I,J,1) = REAL(TEMP,4)
         DAT(I,J,1) = REAL(GLOB_FLUX(I,J),4)

         ! For debuggin purposes, I'm curious about the following:
         DAT(I,J,2) = REAL(SUNCOS( ( J - 1 ) * IIPAR + I ),4)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_AOD,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  3*NWL_MAX+1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     2,         I0+1,
     &            J0+1,      1,         DAT(:,:,1:2)           )


      !--------------------------------------
      ! write average radiative force
      !--------------------------------------
      UNIT     = 'W/m2'


!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, TEMP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Average
         !IF ( !GLOB_FLUX_COUNT > 0 ) THEN 
         !   TEMP    = GLOB_AVE_FORCE(I,J) / REAL(GLOB_FLUX_COUNT)
         !ELSE 
         !   TEMP    = 0d0
         !ENDIF 

         ! Transfer data to 3D, REAL*4 array
         DAT(I,J,1) = REAL(GLOB_AVE_FORCE(I,J),4)
         !DAT(I,J,1) = REAL(DEBUG_J(I,J),4)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_AOD,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  3*NWL_MAX+2,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,         I0+1,
     &            J0+1,      1,         DAT(:,:,1)             )


      
      ! Close file
      CLOSE( IU_AOD )

      ! Return to calling program
      END SUBROUTINE MAKE_AOD_FILE
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_RAD_FILE( YYYYMMDD, HHMMSS, GLOB_FLUX_INST )
!
!******************************************************************************
!  Subroutine MAKE_RAD_FILE saves save instantaneous upward SW radiative flux. 
!  (dkh, 03/05/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : Date of average                          
!  (2 ) HHMMSS   (INTEGER) : Hour-minute-second                      
!
!  Arguments as Onput:
!  ============================================================================
!  (1 ) GLOB_FLUX_INST(REAL*8) : Instantaneous upward rad flux        [W/m2]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : GET_TAU,     EXPAND_DATE

#     include "CMN_SIZE"    ! Size params

      ! Arguments
      INTEGER              :: YYYYMMDD, HHMMSS
      REAL*8               :: GLOB_FLUX_INST(IIPAR,JJPAR)

      ! Local variables 
      INTEGER              :: I, J, I0, J0, L
      CHARACTER(LEN=120)   :: FILENAME
      REAL*4               :: DAT(IIPAR,JJPAR,1)

      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      !=================================================================
      ! MAKE_RAD_FILE begins here!
      !=================================================================
      DAT(:,:,:) = 0d0 

      FILENAME = TRIM( 'rad.YYYYMMDD.hhmm' )

      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS)

      ! Append data directory prefix
      FILENAME = TRIM( './nat_rad/' ) // TRIM( FILENAME )

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'Inst. rad  data file '
      CATEGORY = 'RAD_UPFX'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UNIT     = 'W/m2'

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_RAD_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RAD, FILENAME, TITLE )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J   )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Transfer data to 3D, REAL*4 array
         DAT(I,J,1) = REAL(GLOB_FLUX_INST(I,J),4)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_RAD,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,         I0+1,
     &            J0+1,      1,         DAT                   )

      ! Close file
      CLOSE( IU_RAD )

      ! dkh debug
      print*, 'writing GLOB_FLUX_INST = ', SUM(GLOB_FLUX_INST(:,:))

      ! Return to calling program
      END SUBROUTINE MAKE_RAD_FILE
!------------------------------------------------------------------------------

      SUBROUTINE READ_RAD_FILE( YYYYMMDD, HHMMSS, GLOB_FLUX_NAT )
!     
!******************************************************************************
!  Subroutine READ_RAD_FILE reads the contents of aod.* files. (dkh, 03/05/08) 
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD      : Year-Month-Day 
!  (2 ) HHMMSS        :  and Hour-Min-Sec for which to read aod file
!
!  Arguments passed as output:
!  ============================================================================
!  (1 ) GLOB_FLUX_NAT :  Natural upward fluxes
!
!  Notes
!
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY   : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD, ONLY   : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,  ONLY   : IOERROR
      USE TIME_MOD,  ONLY   : EXPAND_DATE

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
      REAL*8, INTENT(OUT) :: GLOB_FLUX_NAT(IIPAR,JJPAR)

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, NTL
      REAL*4              :: TRACER(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_RAD_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      FILENAME = 'rad.YYYYMMDD.hhmm'

      ! Initialize some variables
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open rad file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( FILENAME )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( './nat_rad/' ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_RAD_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RAD, FILENAME )

      !=================================================================
      ! Read checkpointed variables 
      !=================================================================

      READ( IU_RAD, IOSTAT=IOS )
     &  MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RAD,'read_rad_file:1' )

      READ( IU_RAD, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP

      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RAD,'read_rad_file:2')

      READ( IU_RAD, IOSTAT=IOS )
     &     ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RAD,'read_rad_file:3')


      ! Only process rad data (i.e. sw upward flux)
      ! Make sure array dimensions are of global size
      IF ( CATEGORY(1:8) == 'RAD_UPFX' .and. 
     &     NI == IIPAR .and. NJ == JJPAR .and. NL == 1 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            GLOB_FLUX_NAT(I,J) = TRACER(I,J,1)
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE 
        CALL ERROR_STOP('read_rad_file: did not match criteria', 
     &                     'aero_optical_modl.f')
      ENDIF 

      ! dkh debug
      print*, 'reading GLOB_FLUX_NAT = ', SUM(GLOB_FLUX_NAT(:,:))

      ! Close file
      CLOSE( IU_RAD )

      ! Return to calling program
      END SUBROUTINE READ_RAD_FILE

!------------------------------------------------------------------------------

      SUBROUTINE EXPAND_NAME( FILENAME, N_ITRN )
!
!******************************************************************************
!  Subroutine EXPAND_DATE replaces "NN" token within
!  a filename string with the actual values. (bmy, 6/27/02, 12/2/03)
!  (dkh, 9/22/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Filename with tokens to replace
!  (2 ) N_ITRN   (INTEGER  ) : Current iteration number
!
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Modified filename
!
!  NOTES:
!  (1 ) Based on EXPAND_DATE
!
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRREPL
      USE ERROR_MOD,   ONLY : ERROR_STOP

#     include "define.h"

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME
      INTEGER,          INTENT(IN)    :: N_ITRN

      ! Local variables
      CHARACTER(LEN=2)                :: NN_STR

      !=================================================================
      ! EXPAND_NAME begins here!
      !=================================================================

#if   defined( LINUX_PGI )

      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
      ENCODE( 2, '(i2.2)', NN_STR   ) N_ITRN

#else
      
      ! For other platforms, use an F90 internal write (bmy, 9/29/03)
      WRITE( NN_STR,   '(i2.2)' ) N_ITRN

#endif
      
      ! Replace NN token w/ actual value
      CALL STRREPL( FILENAME, 'NN',   NN_STR   )

      
      ! Return to calling program
      END SUBROUTINE EXPAND_NAME

!-----------------------------------------------------------------------


      SUBROUTINE READ_LER_FILE( MM )
!
!******************************************************************************
!  Subroutine READ_LER_FILE reads data file of TEMIS Lambertian Equivalent 
!  reflectivities. (dkh, 10/31/08) 
!
!  Arguments as input:
!  ============================================================================
!  (1 ) MM : Month
!
!  NOTES:
!  (1 ) Based on read_restart. 
! 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE BPCH2_MOD,   ONLY : GET_RES_EXT
      USE DAO_MOD,     ONLY : AD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IU_RST, IOERROR
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE CHARPAK_MOD, ONLY : STRREPL


#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! STT, NSRCX, TCVV, NTRACE, TCMASS, etc..
#     include "CMN_DIAG"   ! TRCOFFSET

      ! Arguments
      INTEGER, INTENT(IN) :: MM

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N
      INTEGER             :: NCOUNT(NNPAR) 
      REAL*4              :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      LOGICAL, SAVE       :: FIRST = .TRUE.

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER(LEN=2)    :: MM_STR
      CHARACTER(LEN=40)   :: INPUT_LER_FILE

      !=================================================================
      ! READ_LER_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_LER_FILE = 'temis_ler.MM.bpch'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_LER_FILE )

      WRITE( MM_STR, '(i2.2)' ) MM
         
      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL STRREPL( FILENAME, 'MM', MM_STR )

      FILENAME = TRIM( DATA_DIR ) // 'temis_ler_200810/' 
     &         // TRIM( FILENAME ) // '.' // GET_RES_EXT()

      ! Echo some input to the screen
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_LER_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      
      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      READ( IU_RST, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 .or. IOS < 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_ler_file:4' )

      READ( IU_RST, IOSTAT=IOS ) 
     &      CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &      NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &      NSKIP
          
      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ler_file:5')

      READ( IU_RST, IOSTAT=IOS ) 
     &    ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ler_file:6')

      !==============================================================
      ! Assign data from the TRACER array to the STT array.
      !==============================================================
  
      ! Save in TEMIS_LER and convert to fraction (data in file is
      ! stored as fraction * 1d3). 

      ! Use values at 335 for 315 
      TEMIS_LER(:,:,1) = TRACER(:,:,1) / 1d3

      ! Interpolate vales at 357 from values at 380 and 335
      TEMIS_LER(:,:,2) = ( TRACER(:,:,2)  - TRACER(:,:,1) )
     &                 / ( 380 - 335 )
     &                 * ( LAMBDA_TABLE(2) - 335 ) / 1d3
     &                 + ( TRACER(:,:,1) ) /  1d3

      ! Interpolate values at 500 from values at 494 and 555
      TEMIS_LER(:,:,3) = ( TRACER(:,:,7)  - TRACER(:,:,6) )
     &                 / ( 555 - 494 )
     &                 * ( LAMBDA_TABLE(3) - 494 ) / 1d3
     &                 + ( TRACER(:,:,6) ) /  1d3

      ! Use values at 758 for 833 and 1667 
      TEMIS_LER(:,:,4) = ( TRACER(:,:,7)  - TRACER(:,:,6) )
     &                 / ( 555 - 494 )
     &                 * ( LAMBDA_TABLE(3) - 494 ) / 1d3
     &                 + ( TRACER(:,:,6) ) /  1d3
      TEMIS_LER(:,:,5) = TEMIS_LER(:,:,4)

      ! Close file
      CLOSE( IU_RST )      


      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_LER_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_LER_FILE

!------------------------------------------------------------------------------

      SUBROUTINE INIT_LIDORT
!
!*****************************************************************************
!  Subroutine INIT_LIDORT deallocates all module arrays.  (dkh, 11/09/07) 
!        
!  NOTES:
!  (1 ) Udate to LIDORT 3p5, change name from INIT_AERO_OPTICAL to INIT_LIDORT
!          (dkh, 07/25/10) 
!  (2 ) Check to make sure JJPAR matches MAXTHREADS
!******************************************************************************
!     
      USE ERROR_MOD,              ONLY : ALLOC_ERR
      USE ERROR_MOD,              ONLY : ERROR_STOP
      
#     include "CMN_SIZE"               ! IIPAR, JJPAR, etc

      !  include file for LIDORT Mk 2 threads
#include      "./thread_sourcecode_MkII_F90/LIDORT.PARS_F90"


      ! Local variables
      INTEGER :: AS

      !=================================================================      
      ! INIT_LIDORT begins here
      !================================================================= 

      ALLOCATE( AOD_SAVE( IIPAR, JJPAR, NWL_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AOD_SAVE' )
      AOD_SAVE = 0d0

      ALLOCATE( AOD_AVE( IIPAR, JJPAR, NWL_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AOD_AVE' )
      AOD_AVE = 0d0

      ALLOCATE( LAYER_AOD( IIPAR, JJPAR, LLPAR, NWL_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAYER_AOD' )
      LAYER_AOD = 0d0

      ALLOCATE( LAYER_SSA( IIPAR, JJPAR, LLPAR, NWL_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAYER_SSA' )
      LAYER_SSA = 0d0

      ALLOCATE( LAYER_PHM( IIPAR, JJPAR, LLPAR, NPM_MAX, NWL_MAX ),
     &    STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAYER_PHM' )
      LAYER_PHM = 0d0

      ALLOCATE( LAYER_AOD_ADJ( IIPAR, JJPAR, LLPAR, NWL_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAYER_AOD_ADJ' )
      LAYER_AOD_ADJ = 0d0

      ALLOCATE( LAYER_SSA_ADJ( IIPAR, JJPAR, LLPAR, NWL_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAYER_SSA_ADJ' )
      LAYER_SSA_ADJ = 0d0

      ALLOCATE( LAYER_PHM_ADJ( IIPAR, JJPAR, LLPAR, NPM_MAX, NWL_MAX ),
     &    STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAYER_PHM_ADJ' )
      LAYER_PHM_ADJ = 0d0

      ALLOCATE( MASS_ND( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASS_ND' )
      MASS_ND = 0d0

      ALLOCATE( GLOB_FLUX_W( NWL_MAX, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GLOB_FLUX_W' )
      GLOB_FLUX_W = 0d0

      ALLOCATE( GLOB_FLUX_W_ADJ( NWL_MAX, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GLOB_FLUX_W_ADJ' )
      GLOB_FLUX_W_ADJ = 0d0

      ALLOCATE( GLOB_FLUX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GLOB_FLUX' )
      GLOB_FLUX = 0d0

      ALLOCATE( 
     &   JAC_GLOB_FLUX_W( MAX_ATMOSWFS, LLPAR, NWL_MAX, IIPAR, JJPAR ),
     &   STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'JAC_GLOB_FLUX_W' )
      JAC_GLOB_FLUX_W = 0d0

      ALLOCATE( GLOB_AVE_FORCE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GLOB_AVE_FORCE' )
      GLOB_AVE_FORCE = 0d0

      ALLOCATE( TEMIS_LER( IIPAR, JJPAR, LLER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TEMIS_LER' )
      TEMIS_LER = 0d0

      ALLOCATE( DEBUG_J( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DEBUG_J' )
      DEBUG_J = 0d0

      ! make sure we've set MAXTHREADS in LIDORT.PAR to match
      ! the current model resolution 
      IF ( MAXTHREADS .ne. JJPAR ) THEN
         CALL ERROR_STOP( 'MAXTHREADS /= JJPAR', 
     &                     'lidort_mod.f' )
      ENDIF 

      ! Return to calling program 
      END SUBROUTINE INIT_LIDORT


!------------------------------------------------------------------------------
      SUBROUTINE CLEANUP_LIDORT
!
!*****************************************************************************
!  Subroutine CLEANUP_LIDORT deallocates all module arrays. (dkh, 11/09/07) 
!        
!  NOTES:
!
!******************************************************************************
!     
      IF ( ALLOCATED( AOD_SAVE        ) )  DEALLOCATE( AOD_SAVE        )
      IF ( ALLOCATED( AOD_AVE         ) )  DEALLOCATE( AOD_AVE         )
      IF ( ALLOCATED( LAYER_AOD       ) )  DEALLOCATE( LAYER_AOD       )
      IF ( ALLOCATED( LAYER_SSA       ) )  DEALLOCATE( LAYER_SSA       )
      IF ( ALLOCATED( LAYER_PHM       ) )  DEALLOCATE( LAYER_PHM       )
      IF ( ALLOCATED( LAYER_AOD_ADJ   ) )  DEALLOCATE( LAYER_AOD_ADJ   )
      IF ( ALLOCATED( LAYER_SSA_ADJ   ) )  DEALLOCATE( LAYER_SSA_ADJ   )
      IF ( ALLOCATED( LAYER_PHM_ADJ   ) )  DEALLOCATE( LAYER_PHM_ADJ   )
      IF ( ALLOCATED( MASS_ND         ) )  DEALLOCATE( MASS_ND         )
      IF ( ALLOCATED( GLOB_FLUX_W     ) )  DEALLOCATE( GLOB_FLUX_W     )
      IF ( ALLOCATED( GLOB_FLUX_W_ADJ ) )  DEALLOCATE( GLOB_FLUX_W_ADJ )
      IF ( ALLOCATED( JAC_GLOB_FLUX_W ) )  DEALLOCATE( JAC_GLOB_FLUX_W )
      IF ( ALLOCATED( GLOB_FLUX       ) )  DEALLOCATE( GLOB_FLUX       )
      IF ( ALLOCATED( GLOB_AVE_FORCE  ) )  DEALLOCATE( GLOB_AVE_FORCE  )
      IF ( ALLOCATED( TEMIS_LER       ) )  DEALLOCATE( TEMIS_LER       )

      ! Return to calling program 
      END SUBROUTINE CLEANUP_LIDORT
!------------------------------------------------------------------------------


      END MODULE LIDORT_MOD 
