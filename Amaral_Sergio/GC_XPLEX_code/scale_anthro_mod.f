!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: scale_anthro_mod
!
! !DESCRIPTION: Module SCALE\_ANTHRO\_MOD contains routines to scale 
!  anthropogenic emissions from a base year to a simulation year.
!\\
!\\
! !INTERFACE: 
!
      MODULE SCALE_ANTHRO_MOD
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
      PUBLIC  :: GET_ANNUAL_SCALAR
      PUBLIC  :: GET_ANNUAL_SCALAR_1x1 
      PUBLIC  :: GET_ANNUAL_SCALAR_05x0666_NESTED
      ! add GET_ANNUAL_SCALAR_05x0666_NESTED_CH for backward compatability (dkh, 02/19/11) 
      PUBLIC  :: GET_ANNUAL_SCALAR_05x0666_NESTED_CH

!
! !REVISION HISTORY:
!  28 Jan 2009 - A. v. Donkelaar and P. Le Sager - Initial Version
!
! !REMARKS:
!  (1 ) Add GET_ANNUAL_SCALAR_05x0666_NESTED_CH for nested grid simulations 
!        over China. (tmf, 12/3/09)
!  (2 ) Renamed consistently variables: name depends on relation of variable 
!        to BASE or TARGET year. New data directory to account for updated
!        scale factors for 1985-1989 (phs, 5/7/09)
!  (3 ) Adjusted GET_ANNUAL_SCALAR_05x0666_CH for new scalar format and 
!        renamed to GET_ANNUAL_SCALAR_05x0666 (amv, 10/29/2009)
!  18 Dec 2009 - Aaron van D - Updated scale factors thru 2006  
!  18 Dec 2009 - Aaron van D - Updated routine GET_ANNUAL_SCALAR_05x0666_NESTED
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_annual_scalar
!
! !DESCRIPTION: Subroutine GET\_ANNUAL\_SCALAR returns annual scale 
!  factors to convert B\_YEAR (base year) to T\_YEAR (simulation year), 
!  on the current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ANNUAL_SCALAR( TRACER, B_YEAR, T_YEAR, AS )
!
! !USES:
!
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1

#     include "CMN_SIZE"                         ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: TRACER           ! Tracer number
      INTEGER, INTENT(IN)    :: B_YEAR           ! Base year of emissions
      INTEGER, INTENT(IN)    :: T_YEAR           ! Target year of emissions
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),  INTENT(INOUT) :: AS(IIPAR,JJPAR)  ! Scale factor array
!
! !REVISION HISTORY: 
!  28 Jan 2009 - A. v. Donkelaar and P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)                        :: AS_1x1(I1x1,J1x1)
      TYPE (XPLEX)                        :: AS_1x1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                        :: AS_R8(IIPAR,JJPAR)

      ! Read 1x1 scale factors
      CALL GET_ANNUAL_SCALAR_1x1( TRACER, B_YEAR, T_YEAR, AS_1x1 )

      ! Cast to TYPE (XPLEX)
      AS_1x1x1(:,:,1) = AS_1x1(:,:)

      ! Regrid emissions factors to current model resolution
      CALL DO_REGRID_1x1( 'unitless', AS_1x1x1, AS_R8 )

      ! Cast to TYPE (XPLEX)
      AS(:,:) = AS_R8(:,:)

      END SUBROUTINE GET_ANNUAL_SCALAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_annual_scalar_1x1
!
! !DESCRIPTION: Subroutine GET\_ANNUAL\_SCALAR\_1x1 returns annual scale 
!  factors to convert B\_YEAR (base year) to T\_YEAR (target year), on the 1x1 
!  GEOS-Chem grid.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ANNUAL_SCALAR_1x1( TRACER, B_YEAR, T_YEAR, AS_1x1 )
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE BPCH2_MOD,     ONLY : GET_TAU0, READ_BPCH2

#     include "CMN_SIZE"                           ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: TRACER             ! Tracer number
      INTEGER, INTENT(IN)    :: B_YEAR             ! Base year of emissions
      INTEGER, INTENT(IN)    :: T_YEAR             ! Target year of emissions
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),   INTENT(OUT)  :: AS_1x1(I1x1,J1x1)  ! Scale factor array
!
! !REVISION HISTORY:
!  28 Jan 2009 - A. v. Donkelaar and P. Le Sager - Initial Version
!
! !REMARKS:
!  (1) Scaling factors are for years between 1985 and 2005, on the GEOS-Chem
!       1x1 grid (phs, 3/10/08)
!  18 Dec 2009 - Aaron van D - Updated scale factors through 2006, 
!                              changed to new, directory, reset year limits
!  18 Dec 2009 - Aaron van D - Reformated scale factors to a single file for 
!                              all years, made necessary input changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)                        :: T_1x1(I1x1,J1x1)
      TYPE (XPLEX)                        :: B_1x1(I1x1,J1x1)
      TYPE (XPLEX)                        :: TAU
      CHARACTER(LEN=255)            :: FILENAME,      SCALE_DIR
      CHARACTER(LEN=4)              :: BASE_YYYY_STR, TARG_YYYY_STR
      INTEGER                       :: BASE_YEAR,     TARG_YEAR
      INTEGER                       :: I, J

      !=================================================================
      ! GET_ANNUAL_SCALAR_1x1 begins here!
      !=================================================================

      SCALE_DIR = TRIM( DATA_DIR_1x1 ) // 'anth_scale_factors_200911/'

      ! limit scaling between available years
      BASE_YEAR = MAX( MIN( B_YEAR, 2006 ), 1985 )
      TARG_YEAR = MAX( MIN( T_YEAR, 2006 ), 1985 )

      WRITE( BASE_YYYY_STR, '(i4.4)' ) BASE_YEAR
      WRITE( TARG_YYYY_STR, '(i4.4)' ) TARG_YEAR

      IF ( BASE_YEAR == 2000 ) THEN

         B_1x1(:,:) = 1.d0

      ELSE

         ! Filename
         IF ( TRACER == 71 ) THEN

            ! NOx
            FILENAME = TRIM( SCALE_DIR ) //
     &          'NOx-AnnualScalar.geos.1x1'

         ELSE IF ( TRACER == 72 ) THEN

            ! CO
            FILENAME = TRIM( SCALE_DIR ) //
!     &          'CO-AnnualScalar.geos.1x1'
     &          'CO-AnnualScalar.201105-Fix.geos.1x1'

         ELSE IF ( TRACER == 73 ) THEN

            ! SOx
            FILENAME = TRIM( SCALE_DIR ) //
     &          'SOx-AnnualScalar.geos.1x1'

         ENDIF

         ! Get Tau    
         TAU = GET_TAU0(1,1,BASE_YEAR)

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - GET_ANNUAL_SCALAR_1x1: Reading ', a )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', TRACER,
     &                    TAU,  I1x1,       J1x1,
     &                    1,        B_1x1,      QUIET=.TRUE. )

      ENDIF

      IF ( TARG_YEAR == 2000 ) THEN

         T_1x1(:,:) = 1.d0

      ELSE

         ! Filename
         IF ( TRACER == 71 ) THEN

            ! NOx
            FILENAME = TRIM( SCALE_DIR ) //
     &          'NOx-AnnualScalar.geos.1x1'

         ELSE IF ( TRACER == 72 ) THEN

            ! CO
            FILENAME = TRIM( SCALE_DIR ) //
!     &          'CO-AnnualScalar.geos.1x1'
     &          'CO-AnnualScalar.201105-Fix.geos.1x1'

         ELSE IF ( TRACER == 73 ) THEN

            ! SOx
            FILENAME = TRIM( SCALE_DIR ) //
     &          'SOx-AnnualScalar.geos.1x1'

         ENDIF

         ! Calc Tau
         TAU = GET_TAU0(1,1,TARG_YEAR)

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', TRACER,
     &                    TAU,  I1x1,       J1x1,
     &                    1,        T_1x1,      QUIET=.TRUE. )

      ENDIF

      ! Get scaling and cast as TYPE (XPLEX)
      AS_1x1(:,:) = T_1x1(:,:) / B_1x1(:,:)

      END SUBROUTINE GET_ANNUAL_SCALAR_1x1
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_annual_scalar_05x0666_nested
!
! !DESCRIPTION:  Subroutine GET\_ANNUAL\_SCALAR\_05x0666\_NESTED
!  returns annual scale factors to convert B\_YEAR (base year) to 
!  T\_YEAR (target year), on the 0.5x0.666 GEOS-Chem grid for nested China 
!  domain.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ANNUAL_SCALAR_05x0666_NESTED
     &                     ( TRACER, B_YEAR, T_YEAR, AS )
! !USES:
!
      USE REGRID_1x1_MOD,       ONLY : DO_REGRID_1x1

#     include "CMN_SIZE"             ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: TRACER
      INTEGER, INTENT(IN)  :: B_YEAR
      INTEGER, INTENT(IN)  :: T_YEAR
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),         INTENT(INOUT) :: AS(IIPAR,JJPAR)
!
! !REVISION HISTORY: 
!  28 Jan 2009 - A. v. Donkelaar and P. Le Sager - Initial Version
!  12 Mar 2009 - T-M. Fu     - Initial Version
!  03 Nov 2009 - Aaron van D - rewritten to employ GET_ANNUAL_SCALAR_1x1
!                              and regrid.
!  18 Dec 2009 - Aaron van D - Renamed to GET_ANNUAL_SCALAR_05x0666_NESTED
!  18 Dec 2009 - Aaron van D - Rewrote GET_ANNUAL_SCALAR_05x0666_NESTED to 
!                              retrieve and regrid scale factors by calling 
!                              GET_ANNUAL_SCALAR_1x1 and regridding on fly
!
! !REMARKS:
!  (1) Scaling factors are for years between 1985 and 2005, on the GEOS-Chem
!       0.5x0.666 grid for China domain (tmf, 3/5/09)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! ! LOCAL VARIABLES:
!
      TYPE (XPLEX)                        :: AS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                        :: AS_R8(IIPAR, JJPAR)

      !=================================================================
      ! GET_ANNUAL_SCALAR_05x0666_NESTED begins here!
      !=================================================================

      CALL GET_ANNUAL_SCALAR_1x1( TRACER, B_YEAR, T_YEAR, AS_1x1 )

      CALL DO_REGRID_1x1( 'unitless', AS_1x1, AS_R8 )
      AS(:,:) = AS_R8(:,:)

      END SUBROUTINE GET_ANNUAL_SCALAR_05x0666_NESTED
!EOC
! Keep GET_ANNUAL_SCALAR_05x0666_NESTED_CH here for backwd compatability (dkh, 02/19/11) 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_ANNUAL_SCALAR_05x0666_NESTED_CH
!
! !DESCRIPTION:  Subroutine GET\_ANNUAL\_SCALAR\_05x0666\_NESTED\_CH 
!  returns annual scale factors to convert B\_YEAR (base year) to 
!  T\_YEAR (target year), on the 0.5x0.666 GEOS-Chem grid for nested China 
!  domain. (avd, bmy, phs, 3/10/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ANNUAL_SCALAR_05x0666_NESTED_CH
     &                     ( TRACER, B_YEAR, T_YEAR, AS )
! !USES:
!
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE BPCH2_MOD,            ONLY : GET_TAU0, READ_BPCH2
      USE REGRID_1x1_MOD,       ONLY : DO_REGRID_05x0666

#     include "CMN_SIZE"             ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: TRACER
      INTEGER, INTENT(IN)  :: B_YEAR
      INTEGER, INTENT(IN)  :: T_YEAR
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE (XPLEX),         INTENT(INOUT) :: AS(IIPAR,JJPAR)
!
! !REVISION HISTORY: 
!  12 Mar 2009 - T-M. Fu - Initial Version
!
! !REMARKS:
!  (1) Scaling factors are for years between 1985 and 2005, on the GEOS-Chem
!       0.5x0.666 grid for China domain (tmf, 3/5/09)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! ! LOCAL VARIABLES:
!
      TYPE (XPLEX)                      :: T_05x0666(I05x0666,J05x0666)
      TYPE (XPLEX)                      :: B_05x0666(I05x0666,J05x0666)
      TYPE (XPLEX)                     :: AS_05x0666(I05x0666,J05x0666)
      TYPE (XPLEX)                 :: AS_05x0666x1(I05x0666,J05x0666,1)
      TYPE (XPLEX)                        :: AS_R8(IIPAR, JJPAR)
      TYPE (XPLEX)                        :: TAU2000
      CHARACTER(LEN=255)            :: FILENAME,      SCALE_DIR
      CHARACTER(LEN=4)              :: BASE_YYYY_STR, TARG_YYYY_STR
      INTEGER                       :: BASE_YEAR,     TARG_YEAR
      INTEGER                       :: I, J


      !=================================================================
      ! GET_ANNUAL_SCALAR_05x0666_NESTED_CH begins here!
      !=================================================================

      SCALE_DIR = TRIM( DATA_DIR ) // 'anth_scale_factors_200811/'

      ! limit scaling between available years
      BASE_YEAR = MAX( MIN( B_YEAR, 2005 ), 1985 )
      TARG_YEAR = MAX( MIN( T_YEAR, 2005 ), 1985 )

      WRITE( BASE_YYYY_STR, '(i4.4)' ) BASE_YEAR
      WRITE( TARG_YYYY_STR, '(i4.4)' ) TARG_YEAR

      IF ( BASE_YEAR == 2000 ) THEN

         B_05x0666(:,:) = 1.0

      ELSE

         ! Filename
         IF ( TRACER == 71 ) THEN

            ! NOx
            FILENAME = TRIM( SCALE_DIR ) // 'NOxScalar-' // 
     &                 BASE_YYYY_STR // '-' // '2000.geos.05x0666'

         ELSE IF ( TRACER == 72 ) THEN

            ! CO
            FILENAME = TRIM( SCALE_DIR ) // 'COScalar-' // 
     &                 BASE_YYYY_STR // '-' // '2000.geos.05x0666'

         ELSE IF ( TRACER == 73 ) THEN

            ! SOx
            FILENAME = TRIM( SCALE_DIR ) // 'SOxScalar-' //
     &                 BASE_YYYY_STR // '-' // '2000.geos.05x0666'

         ENDIF

         ! Get Tau    
         TAU2000 = GET_TAU0(1,1,2000)

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - GET_ANNUAL_SCALAR_05x0666_NESTED_CH: Reading ',
     &                   a )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', TRACER,
     &                    TAU2000,  I05x0666,   J05x0666,
     &                    1,        B_05x0666,  QUIET=.TRUE. )

      ENDIF

      IF ( TARG_YEAR == 2000 ) THEN

         T_05x0666(:,:) = 1.0

      ELSE

         ! Filename
         IF ( TRACER == 71 ) THEN

            ! NOx
            FILENAME = TRIM( SCALE_DIR ) // 'NOxScalar-' //
     &                 TARG_YYYY_STR // '-' // '2000.geos.05x0666'

         ELSE IF ( TRACER == 72 ) THEN

            ! CO
            FILENAME = TRIM( SCALE_DIR ) // 'COScalar-' // 
     &                 TARG_YYYY_STR // '-' // '2000.geos.05x0666'

         ELSE IF ( TRACER == 73 ) THEN

            ! SOx
            FILENAME = TRIM( SCALE_DIR ) // 'SOxScalar-' // 
     &                 TARG_YYYY_STR // '-' // '2000.geos.05x0666'

         ENDIF

         ! Calc Tau
         TAU2000 = GET_TAU0(1,1,2000)

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', TRACER,
     &                    TAU2000,  I05x0666,   J05x0666,
     &                    1,        T_05x0666,  QUIET=.TRUE. )

      ENDIF

      ! Get scaling and cast as TYPE (XPLEX)
      AS_05x0666(:,:) = T_05x0666(:,:) / B_05x0666(:,:)

      ! Recast as 3D array
      AS_05x0666x1(:,:,1) = AS_05x0666(:,:)

      ! Regrid emission factors to current model resolution
      CALL DO_REGRID_05x0666( 1, 'unitless', AS_05x0666x1, AS_R8 )

      AS(:,:) = AS_R8(:,:)

      ! Return to calling program
      END SUBROUTINE GET_ANNUAL_SCALAR_05x0666_NESTED_CH
!EOC
!------------------------------------------------------------------------------
      END MODULE SCALE_ANTHRO_MOD


