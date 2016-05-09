      MODULE AIRCRAFT_MOD
!
!******************************************************************************
!  Module AIRCRAFT_MOD is a replacement for AIRCRAFT_NOX_MOD to include the  
!  monthly averaged FAA aircraft emissions in GEOS-Chem. Species include NOx (with
!  NO/NO2 partitioning, SO2 (with an editable fuel sulfur content, input.geos),
!  SO4 (with an editable S(IV)->S(VI) conversion efficiency), BC based on FOA3
!  below 3000ft AGL and a constant EI above 3000ft AGL (editable), similarly
!  for OC emissions, and hydrocarbon emissions with speciation. 
!  Author:      Steven R.H. Barrett (Cambridge/MIT, 2009)
!  Edited by:   Jamin Koo (MIT, 2011) 
!  Sponsor:     FAA Office of Environment and Energy, ULS Project
!  Notes:
!  (1 ) LTO emissions are defined as <1000 m (not exaclty 3000 ft)
!  (2 ) If LAIRCRAFT in LOGICAL_MOD is off, nothing here gets used, i.e. 
!  (3 ) The module can emit into SMVGEAR below the trop, and into STT above (as a
!        significant fraction of flight is above the trop)
!   (srhb, 8/6/2009)   
!
!  Module Variables:
!  ============================================================================
!  (1 ) CRUISE_FB_MULT   (TYPE (XPLEX) ) : Cruise fuel burn multiplier
!  (2 ) LTO_FB_MULT      (TYPE (XPLEX) ) : LTO fuel burn multiplier
!  (3 ) FSC              (TYPE (XPLEX) ) : Fuel sulfur content (ppm)
!  (4 ) EPSILON          (TYPE (XPLEX) ) : Fuel sulfur -> H2SO4 conversion (%)
!  (5 ) NOX_MULT         (TYPE (XPLEX) ) : NOx multiplier
!  (6 ) HC_MULT          (TYPE (XPLEX) ) : HC multiplier
!  (7 ) BC_MULT          (TYPE (XPLEX) ) : BC multiplier
!  (8 ) OC_MULT          (TYPE (XPLEX) ) : OC multiplier
!  (9 ) CO_MULT          (TYPE (XPLEX) ) : CO multiplier
!  (10) NO2_NOX_CRUISE   (TYPE (XPLEX) ) : % of NOx as NO2 in cruise
!  (11) NO2_NOX_LTO      (TYPE (XPLEX) ) : % of NOx as NO2 in LTO (average)
!  (12) HONO_NOX         (TYPE (XPLEX) ) : % of NOx as HONO on NO2 mass basis
!  (13) EI_BC_CRUISE     (TYPE (XPLEX) ) : EI(BC) in g/kg-fuel for cruise
!  (14) EI_OC_CRUISE     (TYPE (XPLEX) ) : EI(OC) in g/kg-fuel for cruise
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_CURRENT_EMISSIONS: Routine to read NOx emissions from disk
!  (2 ) EMISS_AIRCRAFT        : Routine to emit aircraft pollutants
!  (3 ) INIT_AIRCRAFT         : Routine to allocate/initialize module variables
!  (4 ) CLEANUP_AIRCRAFT      : Routine to deallocate module variables
! 
!  GEOS-CHEM modules referenced by aircraft_mod.f
!  ============================================================================
!  (1 ) error_mod.f    : Module containing NaN and other error check routines
!  (2 ) file_mod.f     : Module containing file unit numbers and error checks
!  (3 ) grid_mod.f     : Module containing horizontal grid information
!  (4 ) pressure_mod.f : Module containing routines to compute P(I,J,L) 
! 
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "aircraft_nox_mod.f"
      !=================================================================
      ! PRIVATE module variables
      PRIVATE :: NAIR

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
               
         ! THC to TOG conversion factor for aircraft emissions
       TYPE (XPLEX), PARAMETER:: THC2TOG= xplex(1.16d0,0d0)
         ! HC speciation is included below
       TYPE (XPLEX), PARAMETER:: ACF_ACET= xplex(0.003693477d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_ALD2= xplex(0.042718224d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_ALK4= xplex(0.213791063d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_C2H6= xplex(0.005214505d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_C3H8= xplex(0.000780871d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_CH2O= xplex(0.123081099d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_PRPE= xplex(0.178041756d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_MACR= xplex(0.005362609d0,0d0)
       TYPE (XPLEX), PARAMETER:: ACF_RCHO= xplex(0.036769436d0,0d0)
      ! Note not all aircraft hydrocarbon species modeled in GEOS-Chem
      ! Switch for LTO emissions
       LOGICAL           :: LLTO_EMIS
      ! Switch for cruise emissions
       LOGICAL           :: LCRUISE_EMIS
   
         ! Cruise fuel burn multiplier
       TYPE (XPLEX)              :: CRUISE_FB_MULT
       TYPE (XPLEX)              :: LTO_FB_MULT
         ! Fuel sulfur content (ppm)
       TYPE (XPLEX)              :: FSC
         ! Fuel sulfur -> H2SO4 conversion efficiency (%)
       TYPE (XPLEX)              :: EPSILON
         ! NOx/HC/BC/OC/CO multiplier to all altitudes
       TYPE (XPLEX)              :: NOX_MULT, HC_MULT, BC_MULT
       TYPE (XPLEX)              :: OC_MULT,  CO_MULT
       TYPE (XPLEX)              :: BGNH3_MULT
         ! Fractions of NO2/HONO of total NOx
       TYPE (XPLEX)         :: NO2_NOX_CRUISE, NO2_NOX_LTO, HONO_NOX
  
         ! Fleet avg. emissions index for black carbon in g/kg for cruise
       TYPE (XPLEX)              :: EI_BC_CRUISE
         ! Fleet avg. emissions index for organic carbon in g/kg for cruise
       TYPE (XPLEX)              :: EI_OC_CRUISE
       INTEGER           :: EXCL_CODE, EXCL_LL_I, EXCL_LL_J
       INTEGER           :: EXCL_UR_I, EXCL_UR_J

      ! Map hour of day to time blocks
       INTEGER, PARAMETER    :: HOUR2HOUR_IDX(24) = (/ 6, 6, !0-1,1-2 
     & 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
     & 6, 6 /) ! Access with hour of day + 1

      ! Define emissions arrays
       TYPE (XPLEX),  ALLOCATABLE :: AC_FB(:,:,:)   ! Fuel burn (kg/s)
       TYPE (XPLEX),  ALLOCATABLE :: AC_CO(:,:,:)   ! CO (kg/s)
       TYPE (XPLEX),  ALLOCATABLE :: AC_HC(:,:,:)   ! HC (kg/s) - need to speciate
       TYPE (XPLEX),  ALLOCATABLE :: AC_NOx(:,:,:)  ! NOx (kg/s)
       TYPE (XPLEX),  ALLOCATABLE :: AC_PMNV(:,:,:) ! Soot (kg/s)
       TYPE (XPLEX),  ALLOCATABLE :: AC_PMFO(:,:,:) ! Organics (kg/s) - only use LTO

       TYPE (XPLEX),  ALLOCATABLE :: EMIS_AC_NOx(:,:,:) ! For SMVGEAR NOx

       CHARACTER(LEN=255)   :: EMISS_DIR      ! Emissions directory

       INTEGER              :: NAIR
 
       ! reading in input.aircraft
       LOGICAL            :: VERBOSE  = .FALSE.
       INTEGER, PARAMETER :: FIRSTCOL = 26
       INTEGER, PARAMETER :: MAXDIM   = 255
       INTEGER            :: IU_GEOS
       LOGICAL            :: LAIRCRAFT
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!=========================================================================

      SUBROUTINE INIT_AIRCRAFT
!
!******************************************************************************
!  Subroutine INIT_AIRCRAFT allocates and initializes module variables.
!  (srhb2, 8/6/09)
!  
!  NOTES:
!  (1 ) It sets all the multipliers, other emission indexes, and fuel sulfur content.
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE CMN_SIZE_MOD
#     include "CMN_SIZE"
#     include "define.h"
      
      ! Local variables
      INTEGER :: AS
      INTEGER :: L

      ! NAIR     - Maximum number for aircraft emissions in SAGE  
#if defined(GEOS_3) && defined(GRIDREDUCED)
	  NAIR = 30
#elif defined(GEOS_5) && defined(GRIDREDUCED)
	  NAIR = 40
#else
#error "AEDT/SAGE aircraft emissions not preprocessed for this condition"
#endif

      CALL READ_AIRCRAFT_MENU
!      ! Default emissions
!      CRUISE_FB_MULT = 1d0
!      LTO_FB_MULT = 1d0
!      NOX_MULT = 1d0
!      HC_MULT = 1d0
!      BC_MULT = 1d0
!      OC_MULT = 1d0
!      CO_MULT = 1d0
!      ! Fractions of NO2/HONO of total NOx
!      NO2_NOX_CRUISE=9 !(%)
!      NO2_NOX_LTO=23 !(%)
!      HONO_NOX=1   !(%)
!      EI_BC_CRUISE=0.03 !(g/kg)
!      EI_OC_CRUISE=0.03 !(g/kg)
!      FSC =  700 
!      EPSILON = 2


      ALLOCATE( AC_FB( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AC_FB' )
      AC_FB = 0d0
  
      ALLOCATE( EMIS_AC_NOx( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_AC_NOx' )
      EMIS_AC_NOx = 0d0
  
      ALLOCATE( AC_CO( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AC_CO' )
      AC_CO = 0d0
  
      ALLOCATE( AC_HC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AC_HC' )
      AC_HC = 0d0
  
      ALLOCATE( AC_NOx( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AC_NOx' )
      AC_NOx = 0d0
  
      ALLOCATE( AC_PMNV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AC_PMNV' )
      AC_PMNV = 0d0
  
      ALLOCATE( AC_PMFO( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AC_PMFO' )
      AC_PMFO = 0d0
  
      ! Return to calling program
      END SUBROUTINE INIT_AIRCRAFT


!============================================================================
! from INPUT_MOD.F
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: split_one_line
!
! !DESCRIPTION: Subroutine SPLIT\_ONE\_LINE reads a line from the input file 
!  (via routine READ\_ONE\_LINE), and separates it into substrings.
!\\
!\\
!  SPLIT\_ONE\_LINE also checks to see if the number of substrings found is 
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori, 
!  you can skip the error check.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, N_EXP, LOCATION ) 
!
! !USES:
!
      USE CHARPAK_MOD, ONLY: STRSPLIT
!
! !INPUT PARAMETERS: 
!
      ! Number of substrings we expect to find
      INTEGER,            INTENT(IN)  :: N_EXP

      ! Name of routine that called SPLIT_ONE_LINE
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION 
!
! !OUTPUT PARAMETERS:
!
      ! Array of substrings (separated by " ")
      CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM)

      ! Number of substrings actually found
      INTEGER,            INTENT(OUT) :: N_SUBSTRS
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                         :: EOF
      CHARACTER(LEN=255)              :: LINE, MSG

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================      

      ! Create error msg
      MSG = 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )

      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = READ_ONE_LINE( EOF )

      ! STOP on End-of-File w/ error msg
      IF ( EOF ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, '(a)' ) 'End of file encountered!' 
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong 
      IF ( N_EXP /= N_SUBSTRS ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
      ENDIF
       
      END SUBROUTINE SPLIT_ONE_LINE
!EOC
!============================================================================
! from INPUT_MOD.F
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_one_line
!
! !DESCRIPTION: Subroutine READ\_ONE\_LINE reads a line from the input file.  
!  If the global variable VERBOSE is set, the line will be printed to stdout.  
!  READ\_ONE\_LINE can trap an unexpected EOF if LOCATION is passed.  
!  Otherwise, it will pass a logical flag back to the calling routine, 
!  where the error trapping will be done.
!\\
!\\
! !INTERFACE:
!
      FUNCTION READ_ONE_LINE( EOF, LOCATION ) RESULT( LINE )
!
! !USES:
!
      USE FILE_MOD,   ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: LOCATION    ! Msg to display
!
! !OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(OUT)          :: EOF         ! Denotes EOF 
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IOS
      CHARACTER(LEN=255) :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file
      READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.

         ! Trap unexpected EOF -- stop w/ error msg if LOCATION is passed
         ! Otherwise, return EOF to the calling program
         IF ( PRESENT( LOCATION ) ) THEN
            MSG = 'READ_ONE_LINE: error at: ' // TRIM( LOCATION )
            WRITE( 6, '(a)' ) MSG
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_ONE_LINE (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            STOP
         ELSE
            RETURN
         ENDIF
      ENDIF

      ! IO Status > 0: true I/O error condition
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_one_line:1' )

      ! Print the line (if necessary)
      IF ( VERBOSE ) WRITE( 6, '(a)' ) TRIM( LINE )

      END FUNCTION READ_ONE_LINE
!EOC
 
 
!------------------------------------------------------------------------------

      SUBROUTINE READ_AIRCRAFT_MENU
!
!******************************************************************************
!  Subroutine READ_AIRCRAFT_MENU reads the AIRCRAFT MENU section of 
!  aircraft configuration file(input.aircraft) (jkoo, 3/3/1011)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IOERROR
          
      ! Local variables
      INTEGER            :: N
      INTEGER            :: IOS
      LOGICAL            :: EOF
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
      CHARACTER(LEN=255) :: FILEAIR
      CHARACTER(LEN=255) :: TITLE

      !=================================================================
      ! READ_AIRCRAFT_MENU begins here!
      !=================================================================

      FILEAIR = 'input.aircraft'

      WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'A I R C R A F T   I N P U T'
      WRITE( 6, 200   ) TRIM( FILEAIR )
 200  FORMAT( 'READ_AIRCRAFT_INPUT_FILE: Reading ', a )

      OPEN( IU_GEOS, FILE=TRIM( FILEAIR ), IOSTAT=IOS )
      IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS, IU_GEOS, 'read_aircraft_file:1' )

      TITLE = READ_ONE_LINE( EOF  )
      TITLE = READ_ONE_LINE( EOF  )
      ! Switch for LTO emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:1' )
      READ( SUBSTRS(1:N), * ) LLTO_EMIS

      ! Switch for cruise emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:2' )
      READ( SUBSTRS(1:N), * ) LCRUISE_EMIS
	  
	  ! Cruise fuel burn multiplier
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:3' )
      READ( SUBSTRS(1:N), * ) CRUISE_FB_MULT%r
      CRUISE_FB_MULT%i=0d0
	  
	  ! LTO fuel burn multiplier
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:4' )
      READ( SUBSTRS(1:N), * ) LTO_FB_MULT%r
      LTO_FB_MULT%i=0d0
	  
	  ! Fuel sulfur content (ppm)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:5' )
      READ( SUBSTRS(1:N), * ) FSC%r
      FSC%i = 0d0
	  
	  ! Fuel sulfur -> H2SO4 conversion efficiency (%)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:6' )
      READ( SUBSTRS(1:N), * ) EPSILON%r
      EPSILON%i = 0d0
	  
	  ! Emissions multipliers
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:7' )
      READ( SUBSTRS(1:N), * ) NOX_MULT%r
      NOX_MULT%i=0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:8' )
      READ( SUBSTRS(1:N), * ) HC_MULT%r
      HC_MULT%i=0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:9' )
      READ( SUBSTRS(1:N), * ) BC_MULT%r
      BC_MULT%i = 0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:10' )
      READ( SUBSTRS(1:N), * ) OC_MULT%r
      OC_MULT%i = 0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:11' )
      READ( SUBSTRS(1:N), * ) CO_MULT%r
      CO_MULT%i = 0d0
	  
	  ! Background (no-aviation) ammonia multiplier
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:12' )
      READ( SUBSTRS(1:N), * ) BGNH3_MULT%r
      BGNH3_MULT%i = 0d0
	  
	  ! NO2/HONO fractions (%) for cruise and LTO (HONO is combined)
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:13' )
      READ( SUBSTRS(1:N), * ) NO2_NOX_CRUISE%r
      NO2_NOX_CRUISE%i = 0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:14' )
      READ( SUBSTRS(1:N), * ) NO2_NOX_LTO%r
      NO2_NOX_LTO%i =0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:15' )
      READ( SUBSTRS(1:N), * ) HONO_NOX%r
      HONO_NOX%i = 0d0
	  
	  ! Fleet avg. EIs for OC and BC for cruise (>1km altitude)
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:16' )
      READ( SUBSTRS(1:N), * ) EI_BC_CRUISE%r
      EI_BC_CRUISE%i =0d0
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:17' )
      READ( SUBSTRS(1:N), * ) EI_OC_CRUISE%r
      EI_OC_CRUISE%i=0d0
	  
	  ! Exclusion code (-1/0/1) and LL and RR corners grid idx
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:18' )
      READ( SUBSTRS(1:N), * ) EXCL_CODE
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_aircraft_menu:19' )
      READ( SUBSTRS(1:N), * ) EXCL_LL_I, EXCL_LL_J
	  CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_aircraft_menu:20' )
      READ( SUBSTRS(1:N), * ) EXCL_UR_I, EXCL_UR_J
	  
	  ! AEDT/SAGE preprocessed emissions directory
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:21' )
      READ( SUBSTRS(1:N), '(a)' ) EMISS_DIR
	  
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aircraft_menu:22' )
      TITLE = READ_ONE_LINE( EOF  )

      CLOSE(IU_GEOS)
      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'AIRCRAFT MENU'
      WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100     ) 'LTO emissions on?           : ', LLTO_EMIS
         WRITE( 6, 100     ) 'Cruise emissions on?        : ', 
     &                        LCRUISE_EMIS
         WRITE( 6, 130     ) 'Cruise fuel burn multiplier : ', 
     &                        CRUISE_FB_MULT
         WRITE( 6, 130     ) 'LTO fuel burn multiplier    : ', 
     &                        LTO_FB_MULT
         WRITE( 6, 130     ) 'Fuel sulfer content (ppm)   : ', FSC
         WRITE( 6, 130     ) 'S(IV)->S(VI) conversion (%) : ', EPSILON
         WRITE( 6, 130     ) 'NOx multiplier              : ', NOX_MULT
         WRITE( 6, 130     ) 'HC multiplier               : ', HC_MULT	
         WRITE( 6, 130     ) 'BC multiplier               : ', BC_MULT	
         WRITE( 6, 130     ) 'OC multiplier               : ', OC_MULT	
         WRITE( 6, 130     ) 'CO multiplier               : ', CO_MULT	
         WRITE( 6, 130     ) 'Background NH3 multiplier   : ', 
     &                        BGNH3_MULT	
         WRITE( 6, 130     ) '% of NOx as NO2 in cruise   : ', 
     &                        NO2_NOX_CRUISE
         WRITE( 6, 130     ) '% of NOx as NO2 in LTO avg. : ', 
     &                        NO2_NOX_LTO	
         WRITE( 6, 130     ) '% of NOx as HONO            : ', HONO_NOX
         WRITE( 6, 130     ) 'EI(BC) in g/kg for cruise   : ',
     &                        EI_BC_CRUISE	
         WRITE( 6, 130     ) 'EI(OC) in g/kg for cruise   : ', 
     &                        EI_OC_CRUISE	
         WRITE( 6, '(  a)' ) 'AEDT/SAGE emissions direct. : ' 
     &                 //        TRIM(EMISS_DIR)
         IF ( EXCL_CODE .EQ. 0 ) THEN
            WRITE( 6, '(  a)' ) 'All aircraft emissions included'
	 ELSE IF ( EXCL_CODE .EQ. -1 ) THEN
	    WRITE( 6, '(  a)' ) 'The following emissions are excluded...'
            WRITE( 6, 120     ) 'LL box I, J: ',
     &                           EXCL_LL_I, EXCL_LL_J
            WRITE( 6, 120     ) 'UR box I, J: ', 
     &                           EXCL_UR_I, EXCL_UR_J
	 ELSE IF ( EXCL_CODE .EQ. 1 ) THEN
			WRITE( 6, '(  a)' ) 'Only the following emissions are included...'
            WRITE( 6, 120     ) 'LL box I, J: ',
     &                           EXCL_LL_I, EXCL_LL_J
            WRITE( 6, 120     ) 'UR box I, J: ', 
     &                           EXCL_UR_I, EXCL_UR_J
         ELSE 
            WRITE (*,* ) 'Wrong emission exclusion selection'
         ENDIF
      WRITE( 6, '(a  )' ) REPEAT( '=', 79 )

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, A   )
 120  FORMAT( A, 2I5 )
 130  FORMAT( A, 2F7.3 )
    
      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================
	  
	  ! None

      ! Return to calling program
      END SUBROUTINE READ_AIRCRAFT_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_CURRENT_EMISSIONS
!
!******************************************************************************
!  Subroutine READ_CURRENT_EMISSIONS reads the aircraft emissions from FAA emissions. 
!  (jkoo, 03/01/11)
!
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD,       ONLY : IU_FILE, IOERROR
      USE TIME_MOD,       ONLY : GET_YEAR, GET_HOUR
      USE TIME_MOD,       ONLY : GET_MONTH, GET_DAY
      USE TIME_MOD,       ONLY : GET_TS_EMIS, GET_HOUR
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1 
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE PRESSURE_MOD, ONLY : GET_PEDGE
!      USE CMN_SIZE_MOD
 
#     include "CMN_SIZE"    ! Size parameters
#     include "define.h"

       CHARACTER(LEN=2)   :: MONTH_STR, DAY_STR
       
       
       ! From FAA emissions database
       INTEGER               :: J, I, K, N, IOS
       INTEGER               :: J0, I0, ii,jj,kk
       INTEGER               :: IREF, JREF, L
       TYPE (XPLEX)                :: FB, CO, HC, NOx, PMNV, PMFO, MTEMP
       TYPE (XPLEX)                :: PLOW, PHIGH, XSUM, PAIR1, PAIR2, TMP

       CHARACTER(LEN=255)    :: FILENAME
       
#if defined(GEOS_3) 
       ! Highest level including LTO emissions
       INTEGER, PARAMETER :: K_LTO_LAST_IDX = 7
       ! Fraction of highest level to consider as LTO
       TYPE (XPLEX), PARAMETER  :: K_LTO_LAST_FRAC = 269.0/328.0
#elif defined(GEOS_4) 
       ! Highest level including LTO emissions
       INTEGER, PARAMETER :: K_LTO_LAST_IDX = 4
       ! Fraction of highest level to consider as LTO
       TYPE (XPLEX), PARAMETER  :: K_LTO_LAST_FRAC = 143.0/692.0
#elif defined(GEOS_5) || defined(MERRA)  
       ! Highest level including LTO emissions
       INTEGER, PARAMETER :: K_LTO_LAST_IDX = 8
       ! Fraction of highest level to consider as LTO
       TYPE (XPLEX),PARAMETER::K_LTO_LAST_FRAC=xplex(55.0/141.0,0d0)
#else
#error "AEDT/SAGE aircraft emissions not preprocessed for the condition"
#endif

#if defined(GEOS_3) && defined(GRID1x1) && defined(NESTED_NA) && defined(GRIDREDUCED)
      FILENAME = TRIM(EMISS_DIR) // 'SAGEULS_GEOS3RG_1x1_NA_'
#elif defined(GEOS_3) && defined(GRID4x5) && defined(GRIDREDUCED)
      FILENAME = TRIM(EMISS_DIR) // 'SAGEULS_GEOS3RG_4x5_'
#elif defined(GEOS_5) && defined(GRID4x5) && defined(GRIDREDUCED)
      FILENAME = TRIM(EMISS_DIR) // 'SAGEULS_GEOS5RG_4x5_'
#elif defined(GEOS_5) && defined(GRID05x0666) && defined(GRIDREDUCED)
      FILENAME = TRIM(EMISS_DIR) // 'SAGEULS_GEOS5RG_05x0666_'
#else
#error "AEDT/SAGE aircraft emissions not preprocessed for the condition"
#endif


      ! Zero emissions arrays
       AC_FB            = 0d0
       AC_CO            = 0d0
       AC_NOx           = 0d0
       AC_PMNV          = 0d0
       AC_PMFO          = 0d0
       AC_HC            = 0d0
       
       ! Construct filename
      IF ( GET_MONTH() < 10 ) THEN
         WRITE( MONTH_STR, '(i1)' )  GET_MONTH()
      ELSE
         WRITE( MONTH_STR, '(i2)' )  GET_MONTH()
      ENDIF
      IF ( GET_DAY() < 10 ) THEN
         WRITE( DAY_STR, '(i1)' )  GET_DAY()
      ELSE
         WRITE( DAY_STR, '(i2)' )  GET_DAY()
      ENDIF

      ! SDE: If working on a leap year, repeat Feb 28th for Feb 29th
      IF (( GET_DAY() .EQ. 29) .AND. ( GET_MONTH() .EQ. 2 )) THEN
        WRITE( DAY_STR, '(i2)' ) 28
      ENDIF


      FILENAME = TRIM(FILENAME) 
     &     // TRIM(MONTH_STR) // '_' // TRIM(DAY_STR)
      FILENAME = TRIM(FILENAME) // '_2006_'
	 
	  SELECT CASE(HOUR2HOUR_IDX(GET_HOUR()+1))
	     CASE(1)
		    FILENAME = TRIM(FILENAME) // '0200-0600.txt'
		 CASE(2)
		    FILENAME = TRIM(FILENAME) // '0600-1000.txt'
		 CASE(3)
		    FILENAME = TRIM(FILENAME) // '1000-1400.txt'
		 CASE(4)
		    FILENAME = TRIM(FILENAME) // '1400-1800.txt'
		 CASE(5)
		    FILENAME = TRIM(FILENAME) // '1800-2200.txt'
		 CASE(6)
		    FILENAME = TRIM(FILENAME) // '2200-0200.txt'
		 CASE DEFAULT
		    PRINT *, 'ERROR in determining aircraft AEDT/SAGE time block'
		    STOP
	  END SELECT
          PRINT *, ''
	  PRINT *,'************************************************'
	  PRINT *,'Reading aircraft emissions: ' // TRIM(FILENAME)


      ! Now open file
       OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS)
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'RCE:1' )
 
       DO 
              !  Read aircraft fuel burn
              !LAT_IDX, LON_IDX, ALT_IDX
              !FB, CO, HC, NOx, PMNV, PMFO
         READ( IU_FILE, *, IOSTAT=IOS ) J, I, K, 
     &   AC_FB(I, J, K)%r, AC_CO(I, J, K)%r, 
     &   AC_HC(I, J, K)%r, AC_NOx(I, J, K)%r, 
     &   AC_PMNV(I, J, K)%r, AC_PMFO(I, J, K)%r
            
         ! IOS < 0 is end of file
         ! IOS > 0 is an I/O error
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'RCE:2' )
       ENDDO ! Reading loop

       ! Get nested-grid offsets
!       I0 = GET_XOFFSET()
!       J0 = GET_YOFFSET(
!       AC_NOx(22,34,1)%i =1e-10!*AC_NOx(22,34,1)%r 
!      where (AC_NOx%r/=0d0) AC_NOx%i = 1d-100
      !AC_CO(:,:,:)%i = 1d-20
      !AC_HC(:,:,:)%i = 1d-20
      !AC_PMNV(:,:,:)%i = 1d-20
      !AC_PMFO(:,:,:)%i = 1d-20
      !AC_FB(:,:,:)%i = 1d-20
      !AC_NOx(:,:,:)%i=1d-100
      ! Correct non-LTO PMNV and PMFO or delete nonLTO/LTO
      ! First BC
       AC_PMNV(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &         AC_FB(:,:,(K_LTO_LAST_IDX+1):LLPAR) *  EI_BC_CRUISE/1000.0d0
      
       MTEMP = (EI_BC_CRUISE/1000.0d0*(1.0d0-K_LTO_LAST_FRAC))
       AC_PMNV(:,:,K_LTO_LAST_IDX) =  AC_PMNV(:,:,K_LTO_LAST_IDX) 
     &                                   * K_LTO_LAST_FRAC
       AC_PMNV(:,:,K_LTO_LAST_IDX) =  AC_PMNV(:,:,K_LTO_LAST_IDX)  
     &                               + AC_FB(:,:,K_LTO_LAST_IDX) * MTEMP
      
      ! Now OC
       AC_PMFO(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &                       AC_FB(:,:,(K_LTO_LAST_IDX+1):LLPAR) 
     &                                *  EI_OC_CRUISE/1000.0d0
      
       MTEMP = (EI_OC_CRUISE/1000.0d0*(1.0d0-K_LTO_LAST_FRAC))
       AC_PMFO(:,:,K_LTO_LAST_IDX) =  AC_PMFO(:,:,K_LTO_LAST_IDX) 
     &                                   * K_LTO_LAST_FRAC
       AC_PMFO(:,:,K_LTO_LAST_IDX) =  AC_PMFO(:,:,K_LTO_LAST_IDX)  
     &                            + AC_FB(:,:,K_LTO_LAST_IDX) * MTEMP


      ! Do we include LTO and/or cruise?
      IF ( .NOT. LLTO_EMIS) THEN
         LTO_FB_MULT = 0.0d0
      ELSE IF ( .NOT. LCRUISE_EMIS ) THEN
         CRUISE_FB_MULT = 0.0d0
      END IF ! Else all emissions included

      ! Put in cruise and LTO FB multipliers
       AC_FB(:,:,1:(K_LTO_LAST_IDX-1)) = AC_FB(:,:,1:(K_LTO_LAST_IDX-1))
     &                                  * LTO_FB_MULT
       AC_FB(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &    AC_FB(:,:,(K_LTO_LAST_IDX+1):LLPAR) * CRUISE_FB_MULT
       MTEMP = K_LTO_LAST_FRAC * LTO_FB_MULT 
     &         + (1.0d0-K_LTO_LAST_FRAC) * CRUISE_FB_MULT
       AC_FB(:,:,K_LTO_LAST_IDX) = AC_FB(:,:,K_LTO_LAST_IDX) * MTEMP
      
       AC_CO(:,:,1:(K_LTO_LAST_IDX-1)) = AC_CO(:,:,1:(K_LTO_LAST_IDX-1))
     &                                  * LTO_FB_MULT
       AC_CO(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &    AC_CO(:,:,(K_LTO_LAST_IDX+1):LLPAR) * CRUISE_FB_MULT
       MTEMP = K_LTO_LAST_FRAC * LTO_FB_MULT 
     &         + (1.0d0-K_LTO_LAST_FRAC) * CRUISE_FB_MULT
       AC_CO(:,:,K_LTO_LAST_IDX) = AC_CO(:,:,K_LTO_LAST_IDX) * MTEMP
      
       AC_HC(:,:,1:(K_LTO_LAST_IDX-1)) = AC_HC(:,:,1:(K_LTO_LAST_IDX-1))
     &                                  * LTO_FB_MULT
       AC_HC(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &    AC_HC(:,:,(K_LTO_LAST_IDX+1):LLPAR) * CRUISE_FB_MULT
       MTEMP = K_LTO_LAST_FRAC * LTO_FB_MULT 
     &         + (1.0d0-K_LTO_LAST_FRAC) * CRUISE_FB_MULT
       AC_HC(:,:,K_LTO_LAST_IDX) = AC_HC(:,:,K_LTO_LAST_IDX) * MTEMP
      
       AC_NOx(:,:,1:(K_LTO_LAST_IDX-1)) = 
     &                   AC_NOx(:,:,1:(K_LTO_LAST_IDX-1)) * LTO_FB_MULT
       AC_NOx(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &    AC_NOx(:,:,(K_LTO_LAST_IDX+1):LLPAR) * CRUISE_FB_MULT
       MTEMP = K_LTO_LAST_FRAC * LTO_FB_MULT 
     &         + (1.0d0-K_LTO_LAST_FRAC) * CRUISE_FB_MULT
       AC_NOx(:,:,K_LTO_LAST_IDX) = AC_NOx(:,:,K_LTO_LAST_IDX) * MTEMP
        
       AC_PMNV(:,:,1:(K_LTO_LAST_IDX-1)) = 
     &                  AC_PMNV(:,:,1:(K_LTO_LAST_IDX-1)) * LTO_FB_MULT
       AC_PMNV(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &    AC_PMNV(:,:,(K_LTO_LAST_IDX+1):LLPAR) * CRUISE_FB_MULT
       MTEMP = K_LTO_LAST_FRAC * LTO_FB_MULT 
     &         + (1.0d0-K_LTO_LAST_FRAC) * CRUISE_FB_MULT
       AC_PMNV(:,:,K_LTO_LAST_IDX) = AC_PMNV(:,:,K_LTO_LAST_IDX) * MTEMP

       AC_PMFO(:,:,1:(K_LTO_LAST_IDX-1)) = 
     &               AC_PMFO(:,:,1:(K_LTO_LAST_IDX-1)) * LTO_FB_MULT
       AC_PMFO(:,:,(K_LTO_LAST_IDX+1):LLPAR) = 
     &    AC_PMFO(:,:,(K_LTO_LAST_IDX+1):LLPAR) * CRUISE_FB_MULT
       MTEMP = K_LTO_LAST_FRAC * LTO_FB_MULT 
     &         + (1.0d0-K_LTO_LAST_FRAC) * CRUISE_FB_MULT
       AC_PMFO(:,:,K_LTO_LAST_IDX) = AC_PMFO(:,:,K_LTO_LAST_IDX) * MTEMP
      
      ! Use overall species multipliers
       AC_NOx(:,:,:)  = AC_NOx(:,:,:)  * NOX_MULT
       AC_HC(:,:,:)   = AC_HC(:,:,:)   * HC_MULT
       AC_PMNV(:,:,:) = AC_PMNV(:,:,:) * BC_MULT
       AC_PMFO(:,:,:) = AC_PMFO(:,:,:) * OC_MULT
       AC_CO(:,:,:)   = AC_CO(:,:,:)   * CO_MULT

	 ! Now do exclusion/inclusion area
      IF ( EXCL_CODE .EQ. 0 ) THEN
		    ! Do nothing - all aircraft emissions included
		 ELSE IF ( EXCL_CODE .EQ. -1 ) THEN
		    ! Exclude specified area
		 AC_NOx(EXCL_LL_I:EXCL_UR_I,EXCL_LL_J:EXCL_UR_J,1:NAIR)
     &          = 0.0d0
		 AC_HC(EXCL_LL_I:EXCL_UR_I,EXCL_LL_J:EXCL_UR_J,1:NAIR)
     &          = 0.0d0
		 AC_PMNV(EXCL_LL_I:EXCL_UR_I,EXCL_LL_J:EXCL_UR_J,1:NAIR)
     &          = 0.0d0
		 AC_PMFO(EXCL_LL_I:EXCL_UR_I,EXCL_LL_J:EXCL_UR_J,1:NAIR)
     &          = 0.0d0
		 AC_CO(EXCL_LL_I:EXCL_UR_I,EXCL_LL_J:EXCL_UR_J,1:NAIR)
     &          = 0.0d0
		 AC_FB(EXCL_LL_I:EXCL_UR_I,EXCL_LL_J:EXCL_UR_J,1:NAIR)
     &          = 0.0d0
	     ELSE IF ( EXCL_CODE .EQ. 1 ) THEN
		    ! Only include specified area
			
			! Delete everything to 'left'
		 AC_NOx(1:(EXCL_LL_I-1),1:JJPAR,1:NAIR)  = 0.0d0
		 AC_HC(1:(EXCL_LL_I-1),1:JJPAR,1:NAIR)  = 0.0d0
		 AC_PMNV(1:(EXCL_LL_I-1),1:JJPAR,1:NAIR)  = 0.0d0
		 AC_PMFO(1:(EXCL_LL_I-1),1:JJPAR,1:NAIR)  = 0.0d0
		 AC_CO(1:(EXCL_LL_I-1),1:JJPAR,1:NAIR)  = 0.0d0
		 AC_FB(1:(EXCL_LL_I-1),1:JJPAR,1:NAIR)  = 0.0d0

			! Delete everything to 'right'
		 AC_NOx((EXCL_UR_I+1):IIPAR,1:JJPAR,1:NAIR) = 0.0d0
		 AC_HC((EXCL_UR_I+1):IIPAR,1:JJPAR,1:NAIR) = 0.0d0
		 AC_PMNV((EXCL_UR_I+1):IIPAR,1:JJPAR,1:NAIR) = 0.0d0
		 AC_PMFO((EXCL_UR_I+1):IIPAR,1:JJPAR,1:NAIR) = 0.0d0
		 AC_CO((EXCL_UR_I+1):IIPAR,1:JJPAR,1:NAIR) = 0.0d0
		 AC_FB((EXCL_UR_I+1):IIPAR,1:JJPAR,1:NAIR) = 0.0d0

			! Delete everything to north
		 AC_NOx(1:IIPAR,(EXCL_UR_J+1):JJPAR,1:NAIR) = 0.0d0
		 AC_HC(1:IIPAR,(EXCL_UR_J+1):JJPAR,1:NAIR) = 0.0d0
		 AC_PMNV(1:IIPAR,(EXCL_UR_J+1):JJPAR,1:NAIR) = 0.0d0
		 AC_PMFO(1:IIPAR,(EXCL_UR_J+1):JJPAR,1:NAIR) = 0.0d0
		 AC_CO(1:IIPAR,(EXCL_UR_J+1):JJPAR,1:NAIR) = 0.0d0
		 AC_FB(1:IIPAR,(EXCL_UR_J+1):JJPAR,1:NAIR) = 0.0d0

			! Delete everything to south
		 AC_NOx(1:IIPAR,1:(EXCL_LL_J-1),1:NAIR) = 0.0d0
		 AC_HC(1:IIPAR,1:(EXCL_LL_J-1),1:NAIR) = 0.0d0
		 AC_PMNV(1:IIPAR,1:(EXCL_LL_J-1),1:NAIR) = 0.0d0
		 AC_PMFO(1:IIPAR,1:(EXCL_LL_J-1),1:NAIR) = 0.0d0
		 AC_CO(1:IIPAR,1:(EXCL_LL_J-1),1:NAIR) = 0.0d0
		 AC_FB(1:IIPAR,1:(EXCL_LL_J-1),1:NAIR) = 0.0d0

			
		 ELSE ! Invalid
	        CALL ERROR_STOP( 
     & 	'Aircraft emissions exlusion flag invalid!',
     &                    'READ_AIRCRAFT_MENU ("input_mod.f")' )
	  ENDIF
       
       ! Check if any emissions <0 
!       AC_FB    =  MAX(AC_FB,0d0)
!       AC_CO    =  MAX(AC_CO,0d0)
!       AC_HC    =  MAX(AC_HC,0d0)
!       AC_NOX   =  MAX(AC_NOX,0d0)
!       AC_PMNV  =  MAX(AC_PMNV,0d0)
!       AC_PMFO  =  MAX(AC_PMFO,0d0)
!       AC_FB    =  MAX(AC_FB,0d0)
       do ii=1,size(AC_FB,1)
         do jj=1,size(AC_FB,2)
           do kk=1,size(AC_FB,3)
              AC_FB(ii,jj,kk) = MAX(AC_FB(ii,jj,kk),0d0)
           enddo
         enddo
       enddo
!       AC_CO    =  MAX(AC_CO,0d0)
       do ii=1,size(AC_CO,1)
         do jj=1,size(AC_CO,2)
           do kk=1,size(AC_CO,3)
              AC_CO(ii,jj,kk) = MAX(AC_CO(ii,jj,kk),0d0)
           enddo
         enddo
       enddo
!       AC_HC    =  MAX(AC_HC,0d0)
       do ii=1,size(AC_HC,1)
         do jj=1,size(AC_HC,2)
           do kk=1,size(AC_HC,3)
              AC_HC(ii,jj,kk) = MAX(AC_HC(ii,jj,kk),0d0)
           enddo
         enddo
       enddo
!       AC_NOX   =  MAX(AC_NOX,0d0)
       do ii=1,size(AC_NOX,1)
         do jj=1,size(AC_NOX,2)
           do kk=1,size(AC_NOX,3)
              AC_NOX(ii,jj,kk) = MAX(AC_NOX(ii,jj,kk),0d0)
           enddo
         enddo
       enddo
!       AC_PMNV  =  MAX(AC_PMNV,0d0)
       do ii=1,size(AC_PMNV,1)
         do jj=1,size(AC_PMNV,2)
           do kk=1,size(AC_PMNV,3)
              AC_PMNV(ii,jj,kk) = MAX(AC_PMNV(ii,jj,kk),0d0)
           enddo
         enddo
       enddo
!       AC_PMFO  =  MAX(AC_PMFO,0d0)
       do ii=1,size(AC_PMFO,1)
         do jj=1,size(AC_PMFO,2)
           do kk=1,size(AC_PMFO,3)
              AC_PMFO(ii,jj,kk) = MAX(AC_PMFO(ii,jj,kk),0d0)
           enddo
         enddo
       enddo

      ! Output total fuelburns etc
	  WRITE (*,*) 'Aircraft fuel burn = ',  
     & SUM(SUM(SUM(AC_FB(:, :, :)%r,3),2),1),
     & SUM(SUM(SUM(AC_FB(:, :, :)%i,3),2),1), ' kg/s'
	  WRITE (*,*) 'Aircraft NOx = ',  
     & SUM(SUM(SUM(AC_NOx(:, :, :)%r,3),2),1), 
     & SUM(SUM(SUM(AC_NOx(:, :, :)%i,3),2),1), ' kg/s'
	  WRITE (*,*) 'Aircraft CO = ',  
     & SUM(SUM(SUM(AC_CO(:, :, :)%r,3),2),1),
     & SUM(SUM(SUM(AC_CO(:, :, :)%i,3),2),1), ' kg/s'
	  WRITE (*,*) 'Aircraft HC= ',  
     & SUM(SUM(SUM(AC_HC(:, :, :)%r,3),2),1), 
     & SUM(SUM(SUM(AC_HC(:, :, :)%i,3),2),1), ' kg/s'
	  WRITE (*,*) 'Aircraft PMNV= ',  
     & SUM(SUM(SUM(AC_PMNV(:, :, :)%r,3),2),1), 
     & SUM(SUM(SUM(AC_PMNV(:, :, :)%i,3),2),1), ' kg/s'
	  WRITE (*,*) 'Aircraft PMFO= ',  
     & SUM(SUM(SUM(AC_PMFO(:, :, :)%r,3),2),1), 
     & SUM(SUM(SUM(AC_PMFO(:, :, :)%i,3),2),1), ' kg/s'

      ! Return to calling program
      END SUBROUTINE READ_CURRENT_EMISSIONS

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_AIRCRAFT
!
!******************************************************************************
!  Subroutine EMISS_AIRCRAFT emits aircraft pollution into SMVGEAR and STT. This
!  is called from DO_EMISSIONS in EMISSIONS_MOD each emissions timestep.
!   (srhb, 8/6/09)
!
!  NOTES: 
!  (1 ) 
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY
      USE TIME_MOD,     ONLY : GET_TS_EMIS, GET_HOUR
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD
      USE LOGICAL_MOD,  ONLY : LVARTROP
      USE TROPOPAUSE_MOD, ONLY: GET_TPAUSE_LEVEL
!      USE CMN_SIZE_MOD
!      USE CMN_MOD
!      USE CMN_DIAG_MOD
       
      TYPE (XPLEX)             :: DTSRCE, TMPMULT
      INTEGER            :: I0, J0, I, J, L, IREF, JREF, TPL
      LOGICAL, SAVE      :: FIRST         = .TRUE.
      INTEGER, SAVE      :: LAST_MONTH    = 0
      INTEGER, SAVE      :: LAST_DAY      = 0
      INTEGER, SAVE      :: LAST_HOUR_IDX = 0

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! PTOP, SIGE, AVP
#     include "CMN_DIAG"  ! Diagnostic switches
      
      !=================================================================
      ! AIREMIS begins here!
      !=================================================================

      
       ! If this is the first time, allocate arrays
       IF ( FIRST ) THEN
          FIRST = .FALSE.

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMAT strings
 100     FORMAT( 'A I R C R A F T   E M I S S I O N S'   )
 110     FORMAT( 'Routines originally by STEVEN BARRETT'    )
 120     FORMAT( 'Edited by Jamin Koo'   )
 130     FORMAT( 'Last Modification Date: 2/24/11'       )

          CALL INIT_AIRCRAFT

       END IF
       
	  ! If not correct hour block and day read in new emissions
	  IF ( HOUR2HOUR_IDX(GET_HOUR()+1) /=  LAST_HOUR_IDX
     & .OR. LAST_MONTH /= GET_MONTH() .OR. LAST_DAY /= GET_DAY()) THEN
	 
	      CALL READ_CURRENT_EMISSIONS
		  LAST_MONTH    = GET_MONTH()
		  LAST_DAY      = GET_DAY()
		  LAST_HOUR_IDX = HOUR2HOUR_IDX(GET_HOUR()+1)
		  
	  END IF ! Otherwise emissions don't need to be read in
	  
      ! Emissions timestep (emissions in kg/timestep)
      DTSRCE  = GET_TS_EMIS() * 60d0

!      PRINT *, 'IDTSO4=', IDTSO4
!	  PRINT *, 'IDTMACR=', IDTMACR
!	  PRINT *, 'IDTRCHO=', IDTRCHO
!	  PRINT *, 'IDTBCPI=', IDTBCPI
!	  PRINT *, 'IDTOCPI=', IDTOCPI
!	  PRINT *, 'IDTNOX=', IDTNOX
!	  PRINT *, 'IDTCO=', IDTCO
!	  PRINT *, 'IDTACET=', IDTACET
!	  PRINT *, 'IDTALD2=', IDTALD2
!	  PRINT *, 'IDTALK4=', IDTALK4
!	  PRINT *, 'IDTC2H6=', IDTC2H6
!	  PRINT *, 'IDTC3H8=', IDTC3H8
!	  PRINT *, 'IDTCH2O=', IDTCH2O
!	  PRINT *, 'IDTPRPE=', IDTPRPE

      ! Do STT-only emissions
       TMPMULT = FSC/1000000.0d0 * EPSILON/100.0d0 
     &                                  * 96.0d0/32.0d0 * DTSRCE
       STT(:,:,1:LLPAR,IDTSO4) = STT(:,:,1:LLPAR,IDTSO4) 
     &                                  + AC_FB(:,:,1:LLPAR) * TMPMULT
!      WRITE (*,*) 'Aircraft SO4= ',  
!     & SUM(SUM(SUM(AC_FB(:,:,:) * TMPMULT/DTSRCE,3),2),1), ' kg/s'

       ! SO2
       TMPMULT = FSC/1000000.0d0 * (1.0d0-EPSILON/100.0d0)
     &                          * 64.0d0/32.0d0 * DTSRCE
       STT(:,:,1:LLPAR,IDTSO2) = STT(:,:,1:LLPAR,IDTSO2) 
     &                                  + AC_FB(:,:,1:LLPAR) * TMPMULT
!       WRITE (*,*) 'Aircraft SO2= ',  
!     & SUM(SUM(SUM(AC_FB(:,:,:) * TMPMULT/DTSRCE,3),2),1), ' kg/s'
      
      ! BC and OC is assumed hydrophobic
       STT(:,:,1:LLPAR,IDTBCPI) = STT(:,:,1:LLPAR,IDTBCPI)
     &                                 + AC_PMNV(:,:,1:LLPAR) * DTSRCE
       STT(:,:,1:LLPAR,IDTOCPI) = STT(:,:,1:LLPAR,IDTOCPI)
     &                                 + AC_PMFO(:,:,1:LLPAR) * DTSRCE

      ! Emitting all emission into STT(Tracer)
      ! NOx fix array
      EMIS_AC_NOx = 0d0

      ! Loop over surface grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

        IF ( LVARTROP ) THEN 
                     ! This is the highest trop box
            TPL = GET_TPAUSE_LEVEL( I, J ) 
        ELSE
                     ! This is the highest trop box
            TPL = GET_TPAUSE_LEVEL( I, J ) - 1
        ENDIF
           
        DO L = 1, LLPAR
          IF ( L > TPL ) THEN
             STT(I,J,L,IDTNOX) = STT(I,J,L,IDTNOX)
     &                        + AC_NOx(I,J,L) * DTSRCE
          ELSE ! In trop
             EMIS_AC_NOx(I,J,L)   = EMIS_AC_NOx(I,J,L) + AC_NOx(I,J,L)
             STT(I,J,L,IDTNOX) = STT(I,J,L,IDTNOX)
     &                        + AC_NOx(I,J,L) * DTSRCE

          ! Start HC emissions -- only in trop as otherwise 
          ! buildup in strat messes up ozone chemistry
             TMPMULT = THC2TOG * ACF_MACR * DTSRCE
             STT(I,J,L,IDTMACR) = STT(I,J,L,IDTMACR) 
     &                                 + AC_HC(I,J,L) * TMPMULT
             TMPMULT = THC2TOG * ACF_RCHO * DTSRCE
             STT(I,J,L,IDTRCHO) = STT(I,J,L,IDTRCHO) 
     &                                 + AC_HC(I,J,L) * TMPMULT

          ! STT CO emissions
             STT(I,J,L,IDTCO) = STT(I,J,L,IDTCO)
     &                           + AC_CO(I,J,L) * DTSRCE
      
          ! STT ACET emissions
             TMPMULT = THC2TOG * ACF_ACET * DTSRCE
             STT(I,J,L,IDTACET) = STT(I,J,L,IDTACET)
     &                           + AC_HC(I,J,L) * TMPMULT
      
          ! STT ALD2 emissions
             TMPMULT = THC2TOG * ACF_ALD2 * DTSRCE
             STT(I,J,L,IDTALD2) = STT(I,J,L,IDTALD2)
     &                           + AC_HC(I,J,L) * TMPMULT
       
          ! STT ALK4 emissions
             TMPMULT = THC2TOG * ACF_ALK4 * DTSRCE
             STT(I,J,L,IDTALK4) = STT(I,J,L,IDTALK4)
     &                           + AC_HC(I,J,L) * TMPMULT
       
          ! STT C2H6 emissions
             TMPMULT = THC2TOG * ACF_C2H6 * DTSRCE
             STT(I,J,L,IDTC2H6) = STT(I,J,L,IDTC2H6)
     &                           + AC_HC(I,J,L) * TMPMULT
      
           ! STT C3H8 emissions
             TMPMULT = THC2TOG * ACF_C3H8 * DTSRCE
             STT(I,J,L,IDTC3H8) = STT(I,J,L,IDTC3H8)
     &                           + AC_HC(I,J,L) * TMPMULT
      
          ! STT CH2O emissions
             TMPMULT = THC2TOG * ACF_CH2O * DTSRCE
             STT(I,J,L,IDTCH2O) = STT(I,J,L,IDTCH2O)
     &                           + AC_HC(I,J,L) * TMPMULT
      
          ! STT PRPE emissions
             TMPMULT = THC2TOG * ACF_PRPE * DTSRCE
             STT(I,J,L,IDTPRPE) = STT(I,J,L,IDTPRPE)
     &                           + AC_HC(I,J,L) * TMPMULT

          ! End HC emissions
                  
               END IF
           END DO
          ! STT NOx emissions
       END DO ! I
       END DO ! J

      ! Return to calling program
      END SUBROUTINE EMISS_AIRCRAFT

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_AIRCRAFT
!
!******************************************************************************
!  Subroutine CLEANUP_AIRCRAFT deallocates module variables. (srhb, 08/27/09)
!  Added additional allocated variables (jkoo,03/02/09)
!******************************************************************************
!
      IF ( ALLOCATED( AC_FB   ) ) DEALLOCATE( AC_FB   )
      IF ( ALLOCATED( AC_CO   ) ) DEALLOCATE( AC_CO   )
      IF ( ALLOCATED( AC_HC   ) ) DEALLOCATE( AC_HC   ) 
      IF ( ALLOCATED( AC_NOx  ) ) DEALLOCATE( AC_NOx  ) 
      IF ( ALLOCATED( AC_PMNV ) ) DEALLOCATE( AC_PMNV ) 
      IF ( ALLOCATED( AC_PMFO ) ) DEALLOCATE( AC_PMFO ) 
      IF ( ALLOCATED( EMIS_AC_NOx ) ) DEALLOCATE( EMIS_AC_NOx ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_AIRCRAFT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE AIRCRAFT_MOD
