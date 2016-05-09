! $Id: vistas_anthro_mod.f,v 1.1.1.1 2009/06/09 21:51:51 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: VISTAS_ANTHRO_MOD
!
! !DESCRIPTION: Module VISTAS\_ANTHRO\_MOD contains variables and routines 
!  to read the  VISTAS anthropogenic emissions. (amv, 11/24/2008)
!\\
!\\
! !INTERFACE: 
!
      MODULE VISTAS_ANTHRO_MOD
! 
! !USES:
!
      USE EPA_NEI_MOD, ONLY : GET_USA_MASK

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_VISTAS_ANTHRO
      PUBLIC :: EMISS_VISTAS_ANTHRO
      PUBLIC :: GET_VISTAS_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      PRIVATE :: INIT_VISTAS_ANTHRO
      PRIVATE :: VISTAS_SCALE_FUTURE
      PRIVATE :: TOTAL_ANTHRO_Tg
!
! !REVISION HISTORY:
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:

      ! Arrays for weekday & weekend emissions
      TYPE (XPLEX),  ALLOCATABLE :: VISTAS_WD_NOx(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: VISTAS_WE_NOx(:,:)

      ! Array for surface area
      TYPE (XPLEX),  ALLOCATABLE :: A_CM2(:)

      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_VISTAS_ANTHRO
!
! !DESCRIPTION: Function GET\_VISTAS\_ANTHRO returns the VISTAS emission for 
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s]. (amv, phs, 1/28/09) 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_VISTAS_ANTHRO( I,        J,           N,  
     &                            WEEKDAY,  MOLEC_CM2_S, KG_S ) 
     &  RESULT( VALUE )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNOx
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N

      ! Return weekday or weekend emissions
      LOGICAL, INTENT(IN)           :: WEEKDAY

      ! OPTIONAL -- return emissions in [molec/cm2/s]
      LOGICAL, INTENT(IN), OPTIONAL :: MOLEC_CM2_S

      ! OPTIONAL -- return emissions in [kg/s]
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S
!
! !RETURN VALUE:
!     
      ! Emissions output
      TYPE (XPLEX)                        :: VALUE    
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS

      !=================================================================
      ! GET_VISTA_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S

      IF ( N == IDTNOx ) THEN

         ! NOx [molec/cm2/s]
         IF ( WEEKDAY ) THEN
            VALUE = VISTAS_WD_NOx(I,J)
         ELSE
            VALUE = VISTAS_WE_NOx(I,J)
         ENDIF

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no VISTAS emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN
            
         ! Convert from [molec/c,2/s] to [kg/s]
         VALUE = VALUE * A_CM2(J) / XNUMOL(N)

      ENDIF

      ! Return to calling program
      END FUNCTION GET_VISTAS_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISS_VISTAS_ANTHRO
!
! !DESCRIPTION: Subroutine EMISS\_VISTAS\_ANTHRO reads the VISTAS emission 
!  fields at 1x1 resolution and regrids them to the current model resolution.
!  (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_VISTAS_ANTHRO
!
! !USES:
! 
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_1x1

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_O3"            ! FSCALYR
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR
      INTEGER                    :: MN, SNo, ScNo
      TYPE (XPLEX)                     :: ARRAY(I1x1,J1x1,1)
      TYPE (XPLEX)                     :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                     :: SC_1x1(I1x1,J1x1)
      TYPE (XPLEX)                     :: TAU2002, TAU
      CHARACTER(LEN=255)         :: FILENAME, VISTAS_DIR
      CHARACTER(LEN=4)           :: SYEAR, SNAME
      CHARACTER(LEN=2)           :: SMN
      CHARACTER(LEN=1)           :: SSMN

      !=================================================================
      ! EMISS_VISTAS_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_VISTAS_ANTHRO
         FIRST = .FALSE.
      ENDIF

      VISTAS_DIR = TRIM( DATA_DIR_1x1 ) // 'VISTAS_200811/'

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

      ! cap maximum scaling year
      IF ( THISYEAR .gt. 2007 ) THEN
         THISYEAR = 2007
      ENDIF

      SNAME = 'NOx'
      SNo = 1
      ScNo = 71
            
      TAU2002 = GET_TAU0( 1, 1, 2002)
      MN = GET_MONTH()

      IF (MN .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'Vistas-' // TRIM(SNAME) // '-'
     &            // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'Vistas-' // TRIM(SNAME) // '-'
     &            // SMN // '.1x1'
      ENDIF

      WRITE( SYEAR, '(i4)' ) THISYEAR


      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_VISTAS_ANTHRO: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                    TAU2002,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX) before regridding
      GEOS_1x1(:,:,1) = ARRAY(:,:,1)

      ! Load ozone season regulation factors
      IF (MN .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &            // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &            // SMN // '.1x1'
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'RATIO-2D', ScNo,
     &                    TAU2002,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Apply Ozone Season Scalars
      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * ARRAY(:,:,1)

      ! Apply Annual Scalar
      IF ( THISYEAR .ne. 2002 ) THEN
         CALL GET_ANNUAL_SCALAR_1x1( ScNo,     2002,
     &                             THISYEAR, SC_1x1 )

         GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)
      ENDIF

      ! Load/Apply weekend/weekday factors
      TAU = GET_TAU0( MN, 1, 1999)
      FILENAME  = TRIM( VISTAS_DIR )
     &            // 'wkend_an_scalar.nei99.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'WD-WE-$', 2,
     &                    TAU,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1 * ARRAY, 
     &                        VISTAS_WE_NOx )

      FILENAME  = TRIM( VISTAS_DIR )
     &            // 'wkday_an_scalar.nei99.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'WD-WE-$', 1,
     &                    TAU,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1 * ARRAY, 
     &                        VISTAS_WD_NOx )

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL VISTAS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR, MN )

      ! Return to calling program
      END SUBROUTINE EMISS_VISTAS_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: VISTAS_SCALE_FUTURE
!
! !DESCRIPTION: Subroutine VISTAS\_SCALE\_FUTURE applies the IPCC future scale 
!  factors to  the VISTAS anthropogenic emissions. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE VISTAS_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff

#     include "CMN_SIZE"             ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! VISTAS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         VISTAS_WE_NOx(I,J)  = VISTAS_WE_NOx(I,J) 
     &                         * GET_FUTURE_SCALE_NOxff( I, J )
         VISTAS_WD_NOx(I,J)  = VISTAS_WD_NOx(I,J) 
     &                         * GET_FUTURE_SCALE_NOxff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE VISTAS_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_ANTHRO_TG
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the 
!  anthropogenic emissions of NOx. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( YEAR, THISMONTH )
!
! !USES:
! 
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : IDTNOX

#     include "CMN_SIZE"   ! Size parameters
!
! !INPUT PARAMETERS:
!
      ! Year and month of data for which to compute totals
      INTEGER, INTENT(IN) :: YEAR, THISMONTH
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,     J
      TYPE (XPLEX)              :: WD_NOX, WE_NOX, F_NOX, A
      CHARACTER(LEN=3)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      WD_NOX  = 0d0
      WE_NOX  = 0d0
      F_NOX   = TRACER_MW_KG(IDTNOX )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Weekday avg emissions
            WD_NOX  = WD_NOX  + VISTAS_WD_NOX (I,J) * A * F_NOX

            ! Weekend avg emissions
            WE_NOX  = WE_NOX  + VISTAS_WE_NOX (I,J) * A * F_NOX

         ENDDO
      ENDDO

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'VISTAS   U S A   E M I S S I O N S', / )

      ! Weekday avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NOx ', THISMONTH, WD_NOX,  '  '
 110  FORMAT( 'Total weekday avg anthro ', a4, ' for 1999/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Weekend avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 120   ) 'NOx ', THISMONTH, WE_NOX,  '  '
 120  FORMAT( 'Total weekend avg anthro ', a4, ' for 1999/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_VISTAS_ANTHRO 
!
! !DESCRIPTION: Subroutine INIT\_VISTAS\_ANTHRO allocates and zeroes all 
!  module arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_VISTAS_ANTHRO
!
! !USES:
! 
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LVISTAS

#     include "CMN_SIZE"    ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_VISTAS_ANTHRO begins here!
      !=================================================================

      ! Return if LVISTAS is false
      IF ( .not. LVISTAS ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( VISTAS_WD_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VISTAS_WD_NOx' )
      VISTAS_WD_NOx = 0d0

      ALLOCATE( VISTAS_WE_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VISTAS_WE_NOx' )
      VISTAS_WE_NOx = 0d0

      !---------------------------------------------------
      ! Pre-store array for grid box surface area in cm2
      !---------------------------------------------------

      ! Allocate array
      ALLOCATE( A_CM2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_CM2' )

      ! Fill array
      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2( J )
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_VISTAS_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_VISTAS_ANTHRO
!
! !DESCRIPTION: Subroutine CLEANUP\_VISTAS\_ANTHRO deallocates all module 
!  arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_VISTAS_ANTHRO
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_STREETS begins here!
      !=================================================================
      IF ( ALLOCATED( A_CM2          ) ) DEALLOCATE( A_CM2          )
      IF ( ALLOCATED( VISTAS_WD_NOx  ) ) DEALLOCATE( VISTAS_WD_NOx  )
      IF ( ALLOCATED( VISTAS_WE_NOx  ) ) DEALLOCATE( VISTAS_WE_NOx  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_VISTAS_ANTHRO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE VISTAS_ANTHRO_MOD
!EOC
