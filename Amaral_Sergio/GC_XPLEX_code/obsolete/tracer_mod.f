! $Id: tracer_mod.f,v 1.1.1.1 2009/06/09 21:51:54 daven Exp $
      MODULE TRACER_MOD
!
!******************************************************************************
!  Module TRACER_MOD contains the GEOS-CHEM tracer array STT plus various
!  other related quantities.  TRACER_MOD also contains inquiry functions that
!  can be used to determine the type of GEOS-CHEM simulation.
!  (bmy, 7/20/04, 10/25/05)
!
!  Added two new variables STT_ADJ and PERT for adjoint calculations.
!  (Kumaresh. 01/24/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) SIM_TYPE               : Number denoting simulation type
!  (2 ) N_TRACERS              : Number of GEOS-CHEM tracers
!  (3 ) N_MEMBERS              : Max # of constituents a tracer can have
!  (4 ) ID_TRACER              : Array of tracer numbers
!  (5 ) ID_EMITTED             : Index of which constituent has the emissions
!  (6 ) STT                    : GEOS-CHEM Tracer array [kg] 
!  (6 ) STT_ADJ                : GEOS-CHEM Adjoint tracer array [kg]
!  (6 ) PERT                   : GEOS-CHEM perturbed Tracer array [kg]
!  (7 ) TCVV                   : Molecular weight air / molecular weight tracer
!  (8 ) TRACER_COEFF           : Coefficient of each tracer constituent
!  (9 ) TRACER_MW_G            : Tracer molecular weight [g/mole]
!  (10) TRACER_MW_KG           : Tracer molecular weight [kg/mole]
!  (11) TRACER_N_CONST         : Array of number of constituents per tracer
!  (12) TRACER_NAME            : Array of tracer names
!  (13) TRACER_CONST           : Array of names for tracer constituents
!  (14) SALA_REDGE_um          : Accum mode seasalt radii bin edges [um]
!  (15) SALC_REDGE_um          : Coarse mode seasalt radii bin edges [um]
!  (16) XNUMOL                 : Ratio of (molec/mole) / (kg/mole) = molec/kg
!  (17) XNUMOLAIR              : XNUMOL ratio for air
!
!  Module Routines:
!  ============================================================================
!  (1 ) ITS_A_RnPbBe_SIM       : Returns TRUE if it's a  Rn-Pb-Be  simulation
!  (2 ) ITS_A_CH3I_SIM         : Returns TRUE if it's a  CH3I      simulation
!  (3 ) ITS_A_FULLCHEM_SIM     : Returns TRUE if it's a  fullchem  simulation
!  (4 ) ITS_A_HCN_SIM          : Returns TRUE if it's a  HCN       simulation
!  (5 ) ITS_A_TAGOX_SIM        : Returns TRUE if it's a  Tagged Ox simulation
!  (6 ) ITS_A_TAGCO_SIM        : Returns TRUE if it's a  Tagged CO simulation
!  (7 ) ITS_A_C2H6_SIM         : Returns TRUE if it's a  C2H6      simulation
!  (8 ) ITS_A_CH4_SIM          : Returns TRUE if it's a  CH4       simulation
!  (9 ) ITS_AN_AEROSOL_SIM     : Returns TRUE if it's an aerosol   simulation
!  (10) ITS_A_MERCURY_SIM      : Returns TRUE if it's a  mercury   simulation
!  (11) ITS_A_CO2_SIM          : Returns TRUE if it's a  CO2       simulation
!  (12) ITS_NOT_COPARAM_OR_CH4 : Returns TRUE if it's not CO param or CH4 
!  (13) GET_SIM_NAME           : Returns the name of the current simulation
!  (14) CHECK_STT              : Checks STT array for NaN, Inf, or negatives
!  (15) INIT_TRACER            : Allocates and zeroes all module arrays
!  (16) CLEANUP_TRACER         : Deallocates all module arrays
!
!  Module Routines:
!  ============================================================================
!  (1 ) error_mod.f : Module containing I/O error and NaN check routines 
!
!  NOTES:
!  (1 ) Added function GET_SIM_NAME (bmy, 5/3/05)
!  (2 ) Removed ITS_A_COPARAM_SIM; the CO-OH param is obsolete (bmy, 6/24/05)
!  (3 ) Added ITS_A_CO2_SIM (pns, bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now added XNUMOL, XNUMOLAIR as module variables (bmy, 10/25/05)
!******************************************************************************
!
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER                        :: SIM_TYPE
      INTEGER                        :: N_TRACERS
      INTEGER,           PARAMETER   :: N_MEMBERS = 10
      REAL*8,            PARAMETER   :: XNUMOLAIR =  6.022d+23 / 
     &                                               28.9644d-3

      ! Arrays
      INTEGER,           ALLOCATABLE :: ID_TRACER(:)
      INTEGER,           ALLOCATABLE :: ID_EMITTED(:)
      INTEGER,           ALLOCATABLE :: TRACER_N_CONST(:)
      REAL*8,            ALLOCATABLE :: PERT(:,:,:,:)
      REAL*8,            ALLOCATABLE :: STT_ADJ(:,:,:,:)
      REAL*8,            ALLOCATABLE :: EMIS_ADJ(:,:,:,:)
      REAL*8,            ALLOCATABLE :: DDEP_ADJ(:,:,:,:)
      REAL*8,            ALLOCATABLE :: EMIS_I_ADJ(:,:,:,:)
      REAL*8,            ALLOCATABLE :: F(:,:,:,:)
      REAL*8,            ALLOCATABLE :: STT(:,:,:,:)
      REAL*8,            ALLOCATABLE :: STT_TMP(:,:,:,:)
      REAL*8,            ALLOCATABLE :: TCVV(:)
      REAL*8,            ALLOCATABLE :: TRACER_COEFF(:,:)
      REAL*8,            ALLOCATABLE :: TRACER_MW_G(:)
      REAL*8,            ALLOCATABLE :: TRACER_MW_KG(:)
      REAL*8,            ALLOCATABLE :: XNUMOL(:)
      CHARACTER(LEN=14), ALLOCATABLE :: TRACER_NAME(:)
      CHARACTER(LEN=14), ALLOCATABLE :: TRACER_CONST(:,:)  
      REAL*8,  ALLOCATABLE           :: FP(:,:)  
      INTEGER, ALLOCATABLE           :: IM(:,:)  

      ! Define seasalt radii bin edges [um] here since these
      ! need to be used both in "seasalt_mod.f" and "drydep_mod.f"
      REAL*8                         :: SALA_REDGE_um(2)
      REAL*8                         :: SALC_REDGE_um(2)
      INTEGER                        :: NHMSb, NYMDb
      REAL*8                         :: TAUb

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION ITS_A_RnPbBe_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_RnPbBe_SIM returns TRUE if we are doing a GEOS-CHEM
!  Rn-Pb-Be simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_RnPbBe_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 1 )

      ! Return to calling program
      END FUNCTION ITS_A_RnPbBe_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_CH3I_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_CH3I_SIM returns TRUE if we are doing a GEOS-CHEM
!  CH3I (Methyl Iodide) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CH3I_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 2 )

      ! Return to calling program
      END FUNCTION ITS_A_CH3I_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_FULLCHEM_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_FULLCHEM_SIM returns TRUE if we are doing a GEOS-CHEM
!  full chemistry/aerosol simulation (i.e. via SMVGEAR). (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_FULLCHEM_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 3 )

      ! Return to calling program
      END FUNCTION ITS_A_FULLCHEM_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_HCN_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_HCN_SIM returns TRUE if we are doing a GEOS-CHEM
!  HCN (Hydrogen Cyanide) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_HCN_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 4 )

      ! Return to calling program
      END FUNCTION ITS_A_HCN_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_TAGOX_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_TAGOX_SIM returns TRUE if we are doing a GEOS-CHEM
!  Tagged Ox simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_TAGOX_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 6 )

      ! Return to calling program
      END FUNCTION ITS_A_TAGOX_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_TAGCO_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_TAGCO_SIM returns TRUE if we are doing a GEOS-CHEM
!  CH3I (Methyl Iodide) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_TAGCO_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 7 )

      ! Return to calling program
      END FUNCTION ITS_A_TAGCO_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_C2H6_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_C2H6_SIM returns TRUE if we are doing a GEOS-CHEM
!  C2H6 (Ethane) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_RnPbBe_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 8 )

      ! Return to calling program
      END FUNCTION ITS_A_C2H6_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_CH4_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_CH4_SIM returns TRUE if we are doing a GEOS-CHEM
!  CH4 (Methane) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CH4_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 9 )

      ! Return to calling program
      END FUNCTION ITS_A_CH4_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_AN_AEROSOL_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_SULFATE_SIM returns TRUE if we are doing a GEOS-CHEM
!  offline Sulfate/Carbon/dust/seasalt aerosol simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_AN_AEROSOL_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 10 )

      ! Return to calling program
      END FUNCTION ITS_AN_AEROSOL_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_MERCURY_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_MERCURY_SIM returns TRUE if we are doing a GEOS-CHEM
!  Hg0/Hg2/HgP offline mercury simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_MERCURY_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 11 )

      ! Return to calling program
      END FUNCTION ITS_A_MERCURY_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_CO2_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_CO2_SIM returns TRUE if we are doing a GEOS-CHEM
!  CO2 offline simulation. (pns, bmy, 7/25/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_MERCURY_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 12 )

      ! Return to calling program
      END FUNCTION ITS_A_CO2_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_NOT_COPARAM_OR_CH4() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_NOT_COPARAM_OR_CH4 returns TRUE if we are doing a 
!  GEOS-CHEM simulation other than CO w/ Parameterized OH or CH4. 
!  (bmy, 7/15/04, 6/24/05)
!
!  NOTES:
!  (1 ) The CO-OH param (SIM_TYPE=5) is now obsolete (bmy, 6/24/05)
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_NOT_COPARAM_OR_CH4 begins here!
      !=================================================================
      VALUE = ( SIM_TYPE /= 9 )

      ! Return to calling program
      END FUNCTION ITS_NOT_COPARAM_OR_CH4

!------------------------------------------------------------------------------

      FUNCTION GET_SIM_NAME() RESULT( NAME )
!
!******************************************************************************
!  Function GET_SIM_NAME returns the name (e.g. "NOx-Ox-Hydrocarbon-Aerosol", 
!  "Tagged CO", etc.) of the GEOS-CHEM simulation. (bmy, 5/3/05, 6/24/05)
!
!  NOTES:
!  (1 ) The CO-OH simulation has been removed (bmy, 6/24/05)
!******************************************************************************
!
      ! Function value
      CHARACTER(LEN=40) :: NAME
      
      !=================================================================
      ! GET_SIM_NAME begins here!
      !=================================================================

      ! Pick proper name for each simulation type
      SELECT CASE( SIM_TYPE )
         CASE( 1 ) 
            NAME = 'Rn-Pb-Be'
         CASE( 2 ) 
            NAME = 'CH3I'
         CASE( 3 ) 
            NAME = 'NOx-Ox-Hydrocarbon-Aerosol'
         CASE( 4 )
            NAME = 'HCN'
         CASE( 5 )
            NAME = ''
         CASE( 6 )
            NAME = 'Tagged Ox'
         CASE( 7 )
            NAME = 'Tagged CO'
         CASE( 8 ) 
            NAME = 'Tagged C2H6'
         CASE( 9 )
            NAME = 'CH4'
         CASE( 10 ) 
            NAME = 'Offline Aerosol'
         CASE( 11 ) 
            NAME = 'Mercury'
         CASE DEFAULT
            NAME = 'UNKNOWN'
       END SELECT

      ! Return to calling program
      END FUNCTION GET_SIM_NAME

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_STT( LOCATION )
!
!******************************************************************************
!  Subroutine CHECK_STT checks the STT tracer array for negative values,
!  NaN values, or Infinity values.  If any of these are found, the code
!  will stop with an error message. (bmy, 3/8/01, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1) LOCATION (CHARACTER) : String describing location of error in code
!
!  NOTES:
!  (1 ) CHECK_STT uses the interfaces defined above -- these will do the
!        proper error checking for either SGI or DEC/Compaq platforms.
!        (bmy, 3/8/01)
!  (2 ) Now call GEOS_CHEM_STOP to shutdown safely.  Now use logicals LNAN,
!        LNEG, LINF to flag if we have error conditions, and then stop the
!        run outside of the parallel DO-loop. (bmy, 11/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas (bmy, 3/23/03)
!  (4 ) Moved from "error_mod.f" to "tracer_mod.f" (bmy, 7/15/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD, ONLY : IT_IS_NAN
      USE ERROR_MOD, ONLY : IT_IS_FINITE

#     include "CMN_SIZE"           ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION

      ! Local variables
      LOGICAL                      :: LNEG, LNAN, LINF
      INTEGER                      :: I,    J,    L,   N
      
      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
         IF ( STT(I,J,L,N) < 0d0 ) THEN 
!$OMP CRITICAL
            LNEG = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LNAN = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LINF = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL            

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Stop the run if any of LNEG, LNAN, LINF is true
      !=================================================================
      IF ( LNEG .or. LNAN .or. LINF ) THEN
         WRITE( 6, 120 ) TRIM( LOCATION )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! FORMAT statements
      !=================================================================
 100  FORMAT( 'CHECK_STT: STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )
 120  FORMAT( 'CHECK_STT: STOP at ', a )

      ! Return to calling program
      END SUBROUTINE CHECK_STT

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_STT_ADJ( LOCATION )
!
!******************************************************************************
!  Subroutine CHECK_STT checks the STT tracer array for negative values,
!  NaN values, or Infinity values.  If any of these are found, the code
!  will stop with an error message. (bmy, 3/8/01, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1) LOCATION (CHARACTER) : String describing location of error in code
!
!  NOTES:
!  (1 ) CHECK_STT uses the interfaces defined above -- these will do the
!        proper error checking for either SGI or DEC/Compaq platforms.
!        (bmy, 3/8/01)
!  (2 ) Now call GEOS_CHEM_STOP to shutdown safely.  Now use logicals LNAN,
!        LNEG, LINF to flag if we have error conditions, and then stop the
!        run outside of the parallel DO-loop. (bmy, 11/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas (bmy, 3/23/03)
!  (4 ) Moved from "error_mod.f" to "tracer_mod.f" (bmy, 7/15/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD, ONLY : IT_IS_NAN
      USE ERROR_MOD, ONLY : IT_IS_FINITE

#     include "CMN_SIZE"           ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION

      ! Local variables
      LOGICAL                      :: LNEG, LNAN, LINF
      INTEGER                      :: I,    J,    L,   N
      
      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         !---------------------------
         ! Check for NaN's
         !---------------------------
         IF ( IT_IS_NAN( STT_ADJ(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LNAN = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT_ADJ(I,J,L,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT_ADJ(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LINF = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT_ADJ(I,J,L,N)
!$OMP END CRITICAL            

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Stop the run if any of LNEG, LNAN, LINF is true
      !=================================================================
      IF ( LNAN .or. LINF ) THEN
         WRITE( 6, 120 ) TRIM( LOCATION )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! FORMAT statements
      !=================================================================
 100  FORMAT( 'CHECK_STT: STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )
 120  FORMAT( 'CHECK_STT: STOP at ', a )

      ! Return to calling program
      END SUBROUTINE CHECK_STT_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRACER
!
!******************************************************************************
!  Subroutine INIT_TRACER initializes all module arrays 
!  (bmy, 7/15/04, 10/25/05)
!
!  NOTES:
!  (1 ) Now allocate XNUMOL (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
!      USE gckpp_Global

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! Allocate arrays 
      !=================================================================
      ALLOCATE( ID_TRACER( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_TRACER' )
      ID_TRACER = 0

      ALLOCATE( ID_EMITTED( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_EMITTED' )
      ID_EMITTED = 0

      ALLOCATE( TCVV( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCVV' )
      TCVV = 0d0

      ALLOCATE( TRACER_NAME( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_NAME' )
      TRACER_NAME = ''

      ALLOCATE( TRACER_MW_G( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MW_G' )
      TRACER_MW_G = 0d0

      ALLOCATE( TRACER_MW_KG( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MW_KG' )
      TRACER_MW_KG = 0d0

      ALLOCATE( TRACER_COEFF( N_TRACERS, N_MEMBERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_COEFF' )
      TRACER_COEFF = 0d0

      ALLOCATE( TRACER_CONST( N_TRACERS, N_MEMBERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_CONST' )
      TRACER_CONST = ''

      ALLOCATE( TRACER_N_CONST( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_N_CONST' )
      TRACER_N_CONST = 0

      ALLOCATE( STT( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT' )
      STT = 0d0          

      ALLOCATE( STT_TMP( IIPAR, JJPAR, LLPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_TMP' )
      STT_TMP = 0d0     

      ALLOCATE( PERT( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT' )
      PERT = 0d0

      ALLOCATE( STT_ADJ( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT' )
      STT = 0d0

      ALLOCATE( XNUMOL( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XNUMOL' )
      XNUMOL = 0

c$$$      ALLOCATE( EMIS_ADJ( IIPAR, JJPAR, LLPAR, NCOEFF ), STAT=AS )
c$$$      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_ADJ' )
c$$$      EMIS_ADJ = 0d0
c$$$
c$$$      ALLOCATE( DDEP_ADJ( IIPAR, JJPAR, LLPAR, NCOEFF ), STAT=AS )
c$$$      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DDEP_ADJ' )
c$$$      DDEP_ADJ = 0d0
c$$$
c$$$      ALLOCATE( EMIS_I_ADJ( IIPAR, JJPAR, LLPAR, 40 ), STAT=AS )
c$$$      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_I_ADJ' )
c$$$      EMIS_I_ADJ = 0d0

      ALLOCATE( F( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'F' )
      F = 0d0

      ALLOCATE( FP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FP' )
      FP = 0d0

      ALLOCATE( IM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IM' )
      IM = 0



      END SUBROUTINE INIT_TRACER

!-----------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_TRACER
!
!******************************************************************************
!  Subroutine CLEANUP_TRACER deallocates all module arrays 
!  (bmy, 7/15/04, 10/25/05)
!
!  NOTES:
!  (1 ) Now deallocates XNUMOL (bmy, 10/25/05)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TRACER begins here!
      !=================================================================
      IF ( ALLOCATED( ID_TRACER      ) ) DEALLOCATE( ID_TRACER      )
      IF ( ALLOCATED( ID_EMITTED     ) ) DEALLOCATE( ID_EMITTED     )
      IF ( ALLOCATED( TCVV           ) ) DEALLOCATE( TCVV           )
      IF ( ALLOCATED( TRACER_NAME    ) ) DEALLOCATE( TRACER_NAME    )
      IF ( ALLOCATED( TRACER_COEFF   ) ) DEALLOCATE( TRACER_COEFF   )
      IF ( ALLOCATED( TRACER_CONST   ) ) DEALLOCATE( TRACER_CONST   )
      IF ( ALLOCATED( TRACER_N_CONST ) ) DEALLOCATE( TRACER_N_CONST )
      IF ( ALLOCATED( TRACER_MW_G    ) ) DEALLOCATE( TRACER_MW_G    )
      IF ( ALLOCATED( TRACER_MW_KG   ) ) DEALLOCATE( TRACER_MW_KG   )
      IF ( ALLOCATED( STT            ) ) DEALLOCATE( STT            )
      IF ( ALLOCATED( STT_TMP        ) ) DEALLOCATE( STT_TMP        )
      IF ( ALLOCATED( STT_ADJ        ) ) DEALLOCATE( STT_ADJ        )
      IF ( ALLOCATED( PERT           ) ) DEALLOCATE( PERT           )
      IF ( ALLOCATED( F              ) ) DEALLOCATE( F              )
      IF ( ALLOCATED( XNUMOL         ) ) DEALLOCATE( XNUMOL         )
c$$$      IF ( ALLOCATED( EMIS_ADJ       ) ) DEALLOCATE( EMIS_ADJ       )
c$$$      IF ( ALLOCATED( EMIS_I_ADJ     ) ) DEALLOCATE( EMIS_I_ADJ     )
c$$$      IF ( ALLOCATED( DDEP_ADJ       ) ) DEALLOCATE( DDEP_ADJ       )
      IF ( ALLOCATED( FP             ) ) DEALLOCATE( FP             )
      IF ( ALLOCATED( IM             ) ) DEALLOCATE( IM             )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRACER

!-----------------------------------------------------------------------------

      ! End of module
      END MODULE TRACER_MOD
