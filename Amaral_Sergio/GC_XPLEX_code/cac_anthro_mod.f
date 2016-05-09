!$Id: cac_anthro_mod.f,v 1.2 2012/03/01 22:00:26 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cac_anthro_mod
!
! !DESCRIPTION: Module CAC\_ANTHRO\_MOD contains variables and routines to 
!  read the  Criteria Air Contaminant Canadian anthropogenic emissions 
!  (amv, phs, 1/28/2009)
!\\
!\\
! !INTERFACE: 
!
      MODULE CAC_ANTHRO_MOD
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
      PUBLIC :: CLEANUP_CAC_ANTHRO
      PUBLIC :: EMISS_CAC_ANTHRO
      PUBLIC :: EMISS_CAC_ANTHRO_05x0666
      PUBLIC :: GET_CANADA_MASK
      PUBLIC :: GET_CAC_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CAC_SCALE_FUTURE
      PRIVATE :: READ_CANADA_MASK
      PRIVATE :: READ_CANADA_MASK_05x0666
      PRIVATE :: INIT_CAC_ANTHRO
      PRIVATE :: TOTAL_ANTHRO_TG
!
! !REVISION HISTORY:
!  28 Jan 2009 - P. Le Sager - Initial Version
!  18 Dec 2009 - Aaron van D - Added EMISS_CAC_ANTHRO_05x0666 routine
!  18 Dec 2009 - Aaron van D - Added READ_CANADA_MASK_05x0666 routine
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!

      ! Arrays for data masks
      INTEGER, ALLOCATABLE :: MASK_CANADA_1x1(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: MASK_CANADA(:,:)
      
      ! Array for surface area
      TYPE (XPLEX),  ALLOCATABLE :: A_CM2(:)

      ! Arrays for emissions
      TYPE (XPLEX),  ALLOCATABLE :: NOx(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: NH3(:,:)
!
! !DEFINED PARAMETERS:
!
      TYPE (XPLEX),PARAMETER::SEC_IN_YEAR=xplex(86400d0*365.25d0,0d0)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_canada_mask
!
! !DESCRIPTION: Function GET\_CANADA\_MASK returns the value of the Canadian 
!  geographic mask at grid box (I,J).  MASK=1 if (I,J) is within Canada, 
!  MASK=0 otherwise. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CANADA_MASK( I, J ) RESULT( THISMASK )
!
! !INPUT PARAMETERS: 
!
      ! Longitude and latitude indices
      INTEGER, INTENT(IN) :: I, J
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      TYPE (XPLEX)              :: THISMASK

      !=================================================================
      ! GET_CANADA_MASK begins here!
      !=================================================================
      THISMASK = MASK_CANADA(I,J)

      END FUNCTION GET_CANADA_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_cac_anthro
!
! !DESCRIPTION: Function GET\_CAC\_ANTHRO returns the Critical Air Contaminants
!  emission for GEOS-Chem grid box (I,J) and tracer N.  Emissions can be 
!  returned in units of [kg/s] or [molec/cm2/s].  (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CAC_ANTHRO( I,    J,     N, 
     &                         MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNOx, IDTCO, IDTSO2, IDTNH3
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N

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
      LOGICAL                       :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_CAC_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      IF ( N == IDTNOx ) THEN

         ! NOx [kg/yr]
         VALUE = NOx(I,J)

      ELSE IF ( N == IDTCO ) THEN

         ! CO [kg/yr]
         VALUE = CO(I,J)

      ELSE IF ( N == IDTSO2 ) THEN

         ! SO2 [kg/yr]
         VALUE = SO2(I,J)

      ELSE IF ( N == IDTNH3 ) THEN

         ! NH3 [kg/yr]
         VALUE = NH3(I,J)

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no CAC emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN
            
         ! Convert from [kg/yr] to [kg/s]
         VALUE = VALUE / SEC_IN_YEAR

      ELSE IF ( DO_MCS ) THEN

         ! Convert NOx from [kg/yr] to [molec/cm2/s]
         VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_YEAR )

      ENDIF

      END FUNCTION GET_CAC_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_cac_anthro
!
! !DESCRIPTION: Subroutine EMISS\_CAC\_ANTHRO reads the Critical Air 
!  Contaminants emission fields at 1x1 resolution and regrids them to the 
!  current model resolution. (amv, phs, 1/28/2009)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_CAC_ANTHRO
!
! !USES:
! 
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : GET_YEAR
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_1x1

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_O3"            ! FSCALYR
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!  (1 ) Emissions are read for a year b/w 2002-2005, and scaled
!        (except NH3) between 1985-2003 if needed (phs, 3/10/08)
!  (2 ) Now accounts for FSCALYR (phs, 3/17/08)
!  18 Dec 2009 - Aaron van D - Use 2005 scale factors for years beyond 2005
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR, SPECIES, SNo, ScNo
      TYPE (XPLEX)                     :: ARRAY(I1x1,J1x1,1)
      TYPE (XPLEX)                     :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                     :: GEOS_1x1_2002(I1x1,J1x1,1)
      TYPE (XPLEX)                     :: GEOS_1x1_2005(I1x1,J1x1,1)
      TYPE (XPLEX)                     :: SC_1x1(I1x1,J1x1)
      TYPE (XPLEX)                     :: TAU2002, TAU2005, TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR, SNAME

      !=================================================================
      ! EMISS_CAC_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_CAC_ANTHRO
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF


      DO SPECIES = 1,4

         IF ( SPECIES .eq. 1 ) THEN
            SNAME = 'NOx'
            SNo = 1
            ScNo = 71
         ELSEIF ( SPECIES .eq. 2 ) THEN
            SNAME = 'CO'
            SNo = 4
            ScNo = 72
         ELSEIF ( SPECIES .eq. 3 ) THEN
            SNAME = 'SOx'
            SNo = 26
            ScNo = 73
         ELSEIF ( SPECIES .eq. 4 ) THEN
            SNAME = 'NH3'
            SNo = 30
            ScNo = 0
         ENDIF
            
         IF ( ( THISYEAR .le. 2002 ) .OR. 
     &        ( THISYEAR .ge. 2005 ) ) THEN
 
            ! TAU values for 2002/2005
            TAU = GET_TAU0( 1, 1, MIN( MAX( THISYEAR, 2002 ), 2005 ) )
            WRITE( SYEAR, '(i4)' ) MIN( MAX( THISYEAR, 2002 ), 2005 )

            ! File name
            FILENAME  = TRIM( DATA_DIR_1x1 ) // 'CAC_200801/CAC' // 
     &                  SYEAR // '-' // TRIM( SNAME ) // '.geos.1x1'

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_CAC_ANTHRO: Reading ', a )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo, 
!phs     &                       TAU,      I1x1,       J1x1-1,     
     &                       TAU,      I1x1,       J1x1,     
     &                       1,        ARRAY,      QUIET=.TRUE. ) 

            ! Cast to TYPE (XPLEX) before regridding
            GEOS_1x1(:,:,1) = ARRAY(:,:,1)

            ! Apply annual scalar factor. Available for 1985-2006,
            ! and NOx, CO and SO2 only.
            IF ( ( THISYEAR .lt. 2002 ) .and. SPECIES .ne. 4 ) THEN

               CALL GET_ANNUAL_SCALAR_1x1( ScNo,     2002, 
     &                                     THISYEAR, SC_1x1 )

               GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

            ELSE IF ((THISYEAR .gt. 2005) .and. SPECIES .ne. 4) THEN

               CALL GET_ANNUAL_SCALAR_1x1( ScNo,     2005,
     &                                     THISYEAR, SC_1x1 )

               GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

            ENDIF

         ELSE

            TAU2002 = GET_TAU0( 1, 1, 2002)
            TAU2005 = GET_TAU0( 1, 1, 2005)

            ! File name for 2002 data
            FILENAME  = TRIM( DATA_DIR_1x1 ) // 'CAC_200801/CAC2002-'
     &                 // TRIM(SNAME) // '.geos.1x1'

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
!phs     &                       TAU2002,   I1x1,      J1x1-1,
     &                       TAU2002,   I1x1,      J1x1,
     &                       1,         ARRAY,     QUIET=.TRUE. )

            ! Cast to TYPE (XPLEX) before regridding
            GEOS_1x1_2002(:,:,1) = ARRAY(:,:,1)

            ! File name for 2005 data
            FILENAME  = TRIM( DATA_DIR_1x1 ) // 'CAC_200801/CAC2005-'
     &                  // TRIM(SNAME) // '.geos.1x1'

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
!phs     &                       TAU2005,   I1x1,      J1x1-1,
     &                       TAU2005,   I1x1,      J1x1,
     &                       1,         ARRAY,     QUIET=.TRUE. )

            ! Cast to TYPE (XPLEX) before regridding
            GEOS_1x1_2005(:,:,1) = ARRAY(:,:,1)

            ! Scale b/w 2002-2005
            GEOS_1x1(:,:,1) = GEOS_1x1_2002(:,:,1) + (THISYEAR-2002.d0)
     &           / 3.d0 * 
     &                   ( GEOS_1x1_2005(:,:,1) - GEOS_1x1_2002(:,:,1) )

         ENDIF

         ! Regrid from GEOS 1x1 --> current model resolution
         IF ( SPECIES .eq. 1 ) THEN

            CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, NOx )

         ELSEIF ( SPECIES .eq. 2 ) THEN

            CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, CO )

         ELSEIF ( SPECIES .eq. 3 ) THEN

            ! Convert SOx to SO2, where SOx is assumed to be 1.4% SO4 and
            ! 98.6% SO2 over NA, based upon Chin et al, 2000, and as
            ! utilized in sulfate_mod.f
            GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * 0.986d0

            CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, SO2 )

         ELSEIF ( SPECIES .eq. 4 ) THEN

            CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, NH3 )

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL CAC_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR )

      END SUBROUTINE EMISS_CAC_ANTHRO
!EOC
!------------------------------------------------------------------------------
!             Dalhousie Atmospheric Compositional Analysis Group              !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_cac_anthro_05x0666
!
! !DESCRIPTION: Subroutine EMISS\_CAC\_ANTHRO\_05x0666 reads the Critical Air
!  Contaminants emission fields at nested NA resolution (1/2 x 2/3)
!  (amv, phs, 11/03/2009)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_CAC_ANTHRO_05x0666
!
! !USES:
!
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE TIME_MOD,          ONLY : GET_YEAR
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_05x0666_NESTED

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_O3"            ! FSCALYR
!
! !REVISION HISTORY:
!  03 Nov 2009 - A. van Donkelaar - Initial Version
!
! !REMARKS:
!  (1 ) Emissions are read for a year b/w 2002-2005, and scaled
!        (except NH3) between 1985-2003 if needed (phs, 3/10/08)
!  (2 ) Now accounts for FSCALYR (phs, 3/17/08)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR, SPECIES, SNo, ScNo
      TYPE (XPLEX)                     :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)                     :: GEOS_05x0666(IIPAR,JJPAR,1)
      TYPE (XPLEX)                   :: GEOS_05x0666_2002(IIPAR,JJPAR,1)
      TYPE (XPLEX)                   :: GEOS_05x0666_2005(IIPAR,JJPAR,1)
      TYPE (XPLEX)                     :: SC_05x0666(IIPAR,JJPAR)
      TYPE (XPLEX)                     :: TAU2002, TAU2005, TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR, SNAME

      !=================================================================
      ! EMISS_CAC_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_CAC_ANTHRO
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF


      DO SPECIES = 1,4

         IF ( SPECIES .eq. 1 ) THEN
            SNAME = 'NOx'
            SNo = 1
            ScNo = 71
         ELSEIF ( SPECIES .eq. 2 ) THEN
            SNAME = 'CO'
            SNo = 4
            ScNo = 72
         ELSEIF ( SPECIES .eq. 3 ) THEN
            SNAME = 'SOx'
            SNo = 26
            ScNo = 73
         ELSEIF ( SPECIES .eq. 4 ) THEN
            SNAME = 'NH3'
            SNo = 30
            ScNo = 0
         ENDIF

         IF ( ( THISYEAR .le. 2002 ) .OR.
     &        ( THISYEAR .ge. 2005 ) ) THEN

            ! TAU values for 2002/2005
            TAU = GET_TAU0( 1, 1, MIN( MAX( THISYEAR, 2002 ), 2005 ) )
            WRITE( SYEAR, '(i4)' ) MIN( MAX( THISYEAR, 2002 ), 2005 )

            ! File name
            FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC' //
     &                  SYEAR // '-' // TRIM( SNAME ) // 
     &                  '.geos.na.1t2x2t3'

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_CAC_ANTHRO_05x0666: Reading ', a )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                       TAU,      IIPAR,       JJPAR,
     &                       1,        ARRAY,      QUIET=.TRUE. )

            ! Cast to TYPE (XPLEX) before regridding
            GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

            ! Apply annual scalar factor. Available for 1985-2006,
            ! and NOx, CO and SO2 only.
            IF ( ( THISYEAR .lt. 2002 ) .and. SPECIES .ne. 4 ) THEN

               CALL GET_ANNUAL_SCALAR_05x0666_NESTED( ScNo, 2002,
     &                                     THISYEAR, SC_05x0666 )

               GEOS_05x0666(:,:,1) = GEOS_05x0666(:,:,1) 
     &                               * SC_05x0666(:,:)

            ELSE IF ((THISYEAR .gt. 2005) .and. SPECIES .ne. 4) THEN

               CALL GET_ANNUAL_SCALAR_05x0666_NESTED( ScNo, 2005,
     &                                     THISYEAR, SC_05x0666 )

              GEOS_05x0666(:,:,1) = GEOS_05x0666(:,:,1) 
     &                               * SC_05x0666(:,:)

            ENDIF

         ELSE

            TAU2002 = GET_TAU0( 1, 1, 2002)
            TAU2005 = GET_TAU0( 1, 1, 2005)

            ! File name for 2002 data
            FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC2002-'
     &                 // TRIM(SNAME) // '.geos.na.1t2x2t3'

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                       TAU2002,   IIPAR,      JJPAR,
     &                       1,         ARRAY,     QUIET=.TRUE. )

            ! Cast to TYPE (XPLEX) before regridding
            GEOS_05x0666_2002(:,:,1) = ARRAY(:,:,1)

            ! File name for 2005 data
            FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC2005-'
     &                  // TRIM(SNAME) // '.geos.na.1t2x2t3'

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                       TAU2005,   IIPAR,    JJPAR,
     &                       1,         ARRAY,     QUIET=.TRUE. )

            ! Cast to TYPE (XPLEX) before regridding
            GEOS_05x0666_2005(:,:,1) = ARRAY(:,:,1)

            ! Scale b/w 2002-2005
            GEOS_05x0666(:,:,1) = GEOS_05x0666_2002(:,:,1) + 
     &              ( THISYEAR - 2002.d0) / 3.d0 *
     &              ( GEOS_05x0666_2005(:,:,1) 
     &                - GEOS_05x0666_2002(:,:,1) )

         ENDIF

         IF ( SPECIES .eq. 1 ) THEN

            NOx(:,:) =  GEOS_05x0666(:,:,1)

         ELSEIF ( SPECIES .eq. 2 ) THEN

            CO(:,:) = GEOS_05x0666(:,:,1)

         ELSEIF ( SPECIES .eq. 3 ) THEN

            ! Convert SOx to SO2, where SOx is assumed to be 1.4% SO4 and
            ! 98.6% SO2 over NA, based upon Chin et al, 2000, and as
            ! utilized in sulfate_mod.f
            SO2(:,:) = GEOS_05x0666(:,:,1) * 0.986d0

         ELSEIF ( SPECIES .eq. 4 ) THEN

            NH3(:,:) = GEOS_05x0666(:,:,1)

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN
         CALL CAC_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR )

      END SUBROUTINE EMISS_CAC_ANTHRO_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cac_scale_future
!
! !DESCRIPTION: Subroutine CAC\_SCALE\_FUTURE applies the IPCC future scale 
!  factors to the Criteria Air Contaminant anthropogenic emissions.
!  (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:

      SUBROUTINE CAC_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

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
      ! STREETS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /yr]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/yr] 
         SO2(I,J)  = SO2(I,J) * GET_FUTURE_SCALE_SO2ff( I, J )

         ! Future NH3 [kg NH3/yr] 
         NH3(I,J)  = NH3(I,J) * GET_FUTURE_SCALE_NH3an( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE CAC_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the 
!  anthropogenic emissions of NOx, CO, SO2 and NH3. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( YEAR )
!
! !USES:
! 
#     include "CMN_SIZE"            ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: YEAR   ! Year of data to compute totals
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,     J
      TYPE (XPLEX)              :: T_NOX, T_CO,  T_SO2,  T_NH3
      CHARACTER(LEN=3)    :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'C. A. C.   C A N A D I A N   E M I S S I O N S', / )


      ! Total NOx [Tg N]
      T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / 46d0 )

      ! Total CO  [Tg CO]
      T_CO  = SUM( CO  ) * 1d-9

      ! Total SO2 [Tg S]
      T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / 64d0 )

      ! Total NH3 [Tg NH3]
      T_NH3 = SUM( NH3 ) * 1d-9

      ! Print totals in [kg]
      WRITE( 6, 110 ) 'NOx ', YEAR, T_NOx,  '[Tg N  ]'
      WRITE( 6, 110 ) 'CO  ', YEAR, T_CO,   '[Tg CO ]'
      WRITE( 6, 110 ) 'SO2 ', YEAR, T_SO2,  '[Tg S  ]'
      WRITE( 6, 110 ) 'NH3 ', YEAR, T_NH3,  '[Tg NH3]'

      ! Format statement
 110  FORMAT( 'C.A.C. Canadian anthro ', a5, 
     &        'for year ', i4, ': ', 2f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      END SUBROUTINE TOTAL_ANTHRO_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_canada_mask
!
! !DESCRIPTION: Subroutine READ\_CANADA\_MASK reads and regrids the Canadian 
!  geographic mask from disk. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CANADA_MASK
!
! !USES:
! 
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1

#     include "CMN_SIZE"       ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)                  :: ARRAY(I1x1,J1x1,1)
      TYPE (XPLEX)                  :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                  :: TAU2000
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_CANADA_MASK begins here!
      !=================================================================

      TAU2000 = GET_TAU0(1,1,2000)

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'CAC_200801/CanadaMask.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CANADA_MASK: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 TAU2000,   I1x1,     J1x1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) before regridding
      GEOS_1x1(:,:,1) = ARRAY(:,:,1)

      ! Save the 1x1 China mask for future use
      MASK_CANADA_1x1(:,:) = GEOS_1x1(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
!      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, MASK_CANADA )

      END SUBROUTINE READ_CANADA_MASK
!EOC
!------------------------------------------------------------------------------
!       Dalhousie University Atmospheric Compositional Analysis Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_canada_mask_05x0666
!
! !DESCRIPTION: Subroutine READ\_CANADA\_MASK\_05x0666 reads the Canadian
!  geographic mask from disk. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CANADA_MASK_05x0666
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1

#     include "CMN_SIZE"       ! Size parameters
!
! !REVISION HISTORY:
!  11 Nov 2009 - A. van Donkelaar - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)                  :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)                  :: GEOS_05x0666(IIPAR,JJPAR,1)
      TYPE (XPLEX)                  :: TAU2000
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_CANADA_MASK begins here!
      !=================================================================

      TAU2000 = GET_TAU0(1,1,2000)

      ! File name
      FILENAME  = TRIM( DATA_DIR ) //
     &            'CAC_200911/CanadaMask.geos.na.1t2x2t3'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CANADA_MASK_05x0666: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 TAU2000,   IIPAR,     JJPAR,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX) before regridding
      MASK_CANADA(:,:) = ARRAY(:,:,1)

      END SUBROUTINE READ_CANADA_MASK_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cac_anthro
!
! !DESCRIPTION: Subroutine INIT\_CAC\_ANTHRO allocates and zeroes all 
!  module arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CAC_ANTHRO
!
! !USES:
! 
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LCAC

#     include "CMN_SIZE"    ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_CAC_ANTHRO begins here!
      !=================================================================

      ! Return if LCAC is false
      IF ( .not. LCAC ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx' )
      NOx = 0d0

      ALLOCATE( CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      ALLOCATE( SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

      ALLOCATE( NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0

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

      !---------------------------------------------------
      ! Read & Regrid masks for CAC emissions
      !---------------------------------------------------

      ALLOCATE( MASK_CANADA_1x1( I1x1, J1x1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CANADA_1x1' )
      MASK_CANADA_1x1 = 0

      ALLOCATE( MASK_CANADA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CANADA' )
      MASK_CANADA = 0d0

      ! Read China & SE Asia masks from disk
      CALL READ_CANADA_MASK

      END SUBROUTINE INIT_CAC_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cac_anthro
!
! !DESCRIPTION: Subroutine CLEANUP\_CAC\_ANTHRO deallocates all module 
!  arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CAC_ANTHRO
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
      IF ( ALLOCATED( MASK_CANADA_1x1) ) DEALLOCATE( MASK_CANADA_1x1) 
      IF ( ALLOCATED( MASK_CANADA    ) ) DEALLOCATE( MASK_CANADA    ) 
      IF ( ALLOCATED( NOx            ) ) DEALLOCATE( NOx            )
      IF ( ALLOCATED( CO             ) ) DEALLOCATE( CO             )
      IF ( ALLOCATED( SO2            ) ) DEALLOCATE( SO2            )
      IF ( ALLOCATED( NH3            ) ) DEALLOCATE( NH3            )

      END SUBROUTINE CLEANUP_CAC_ANTHRO
!EOC
      END MODULE CAC_ANTHRO_MOD

