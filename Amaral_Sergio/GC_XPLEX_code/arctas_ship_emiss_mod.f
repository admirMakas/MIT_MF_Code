! $Id: arctas_ship_emiss_mod.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ARCTAS_SHIP_EMISS_MOD
!
! !DESCRIPTION: Module ARCTAS\_SHIP\_EMISS\_MOD contains variables and 
!  routines to read the Arctas Ship emissions. (phs, 1/28/09)
!\\
!\\
! !INTERFACE: 
!
      MODULE ARCTAS_SHIP_EMISS_MOD
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
      PUBLIC  :: CLEANUP_ARCTAS_SHIP
      PUBLIC  :: EMISS_ARCTAS_SHIP
      PUBLIC  :: GET_ARCTAS_SHIP
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_ARCTAS_SHIP
      PRIVATE :: READ_ARCTAS_SHIP
      PRIVATE :: TOTAL_EMISS_TG
!
! !REVISION HISTORY:
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!  (1) This inventory is based on EDGAR 2000 for NOx, CO, and
!       SO2. But SO2 has been updated by David Street for 2006. BC and OC
!       (from Bond et al, 2004) are also provided. They are a 1996
!       inventory scaled to 2006.  All these emissions were prepared for
!       the ARCTAS 2008 campaign.
!  (2) Only SO2 differs from existing EDGAR/BOND inventories. All other
!       species are disregarded for now, except CO2 that we did not have
!       before.
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!
      ! Arrays
      TYPE (XPLEX),  ALLOCATABLE :: A_CM2(:)

      ! Anthro emiss arrays
      TYPE (XPLEX),  TARGET, ALLOCATABLE :: SO2_SHIP(:,:)
      TYPE (XPLEX),  TARGET, ALLOCATABLE :: CO2_SHIP(:,:)
!
! !DEFINED PARAMETERS:
!
      TYPE (XPLEX),PARAMETER::SEC_IN_YEAR=xplex(86400d0*365.25d0,0d0)

      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_ARCTAS_SHIP
!
! !DESCRIPTION: Function GET\_ARCTAS\_SHIP returns the ARCTAS\_SHIP emission 
!  for GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].  (phs, 1/28/09) 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ARCTAS_SHIP( I, J, N, MOLEC_CM2_S, KG_S ) 
     &     RESULT( VALUE )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : ITS_A_CO2_SIM
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTSO2
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
      ! GET_ARCTAS_SHIP begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.

      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Test for simulation type
      IF ( ITS_A_CO2_SIM() ) THEN
         
         !-------------------
         ! CO2 simulation
         !-------------------
         VALUE  = CO2_SHIP(I,J)

      ELSE

         !-------------------
         ! Other simulations
         !-------------------
         IF ( N == IDTSO2 ) THEN

            ! SO2 [kg/yr]
            VALUE = SO2_SHIP(I,J)

         ELSE

            ! Otherwise return a negative value to indicate
            ! that there are no ARCTAS_SHIP emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN

            ! Convert from [kg/yr] to [kg/s]         
            VALUE = VALUE / SEC_IN_YEAR              

      ELSE IF ( DO_MCS ) THEN

            ! Convert from [kg/yr] to [molec/cm2/s]            
            VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_YEAR )  

      ENDIF

      ! Return to calling program
      END FUNCTION GET_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMISS_ARCTAS_SHIP
!
! !DESCRIPTION: Subroutine EMISS\_ARCTAS\_SHIP reads the ARCTAS\_SHIP emissions
!  from disk. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_ARCTAS_SHIP( YEAR )
!
! !USES:
! 
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1 
      USE TRACER_MOD,     ONLY : ITS_A_CO2_SIM

#     include "CMN_SIZE"       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)     :: YEAR   ! Year of data to read
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE           :: FIRST = .TRUE.
      CHARACTER(LEN=255)      :: FILENAME, DIR

      !=================================================================
      ! EMISS_ARCTAS_SHIP begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_ARCTAS_SHIP
         FIRST = .FALSE.
      ENDIF

 100  FORMAT( '     - EMISS_ARCTAS_SHIP: Reading ', a )

      ! Data directory
      DIR = TRIM( DATA_DIR_1x1 ) // 'ARCTAS_SHIP_2008/'
      

      IF ( ITS_A_CO2_SIM() ) THEN

         !--------------------------
         ! Read CO2 and regrid
         !--------------------------
         FILENAME  = TRIM( DIR ) // 'Arctas_CO2_ship_2008.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_ARCTAS_SHIP( FILENAME, 'CO2-SRCE', 1, CO2_SHIP,
     $        YEAR )
            
      ELSE

         !--------------------------
         ! Read SO2
         !--------------------------
         FILENAME  = TRIM( DIR ) // 'Arctas_SO2_ship_2008.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_ARCTAS_SHIP( FILENAME, 'ANTHSRCE', 26, SO2_SHIP,
     $        YEAR )
         

      ENDIF
     
      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_EMISS_Tg

      ! Return to calling program
      END SUBROUTINE EMISS_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READ_ARCTAS_SHIP
!
! !DESCRIPTION: Subroutine READ\_ARCTAS\_SHIP reads data from one ARCTAS\_SHIP 
!  data file from disk, at GENERIC 1x1 resolution and regrids them to the 
!  current model resolution. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ARCTAS_SHIP( FILENAME, CATEGORY, 
     &                             TRACERN,  ARR,      YEAR )

!
! !USES:
! 
      USE BPCH2_MOD,        ONLY : GET_TAU0,      READ_BPCH2
      USE REGRID_1x1_MOD,   ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR_1x1

#     include "CMN_SIZE"         ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! Year of data to read
      INTEGER,          INTENT(IN)    :: YEAR

      ! Tracer number
      INTEGER,          INTENT(IN)    :: TRACERN

      ! Filename & category under which data is stored in bpch file
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME, CATEGORY
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Array containing output data
      TYPE (XPLEX),           INTENT(INOUT) :: ARR(IIPAR,JJPAR)
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!  (1) Even though the inventory was prepared for Arctas 2008 campaign, CO2 
!       base year is 2000, and SO2 base year is 2006. Input YEAR is used to 
!       scale SO2 into 1985-2004! 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)         :: ARRAY(I1x1,J1x1-1,1)
      TYPE (XPLEX)         :: GEN_1x1(I1x1,J1x1-1)
      TYPE (XPLEX)         :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)         :: SC_1x1(I1x1,J1x1)
      TYPE (XPLEX)         :: TAU2008

      ! TAU0 values for 2008
      TAU2008 = GET_TAU0( 1, 1, 2008 )

      ! Initialize
      SC_1x1 = 1d0

      ! Read data
      CALL READ_BPCH2( FILENAME,  CATEGORY,  TRACERN, 
     &                 TAU2008,   I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)
      
      ! Regrid from GENERIC 1x1 --> GEOS 1x1 
      CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )


      ! Get & apply scaling factor to GEOS 1x1
      IF ( TRACERN == 26 )
     $     CALL GET_ANNUAL_SCALAR_1x1( 73, 2000, YEAR, SC_1x1 )

      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)


      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, ARR )

      END SUBROUTINE READ_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_EMISS_TG
!
! !DESCRIPTION: Subroutine TOTAL\_EMISS\_TG prints the totals for the 
!  anthropogenic or biomass emissions. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_EMISS_TG
!
! !USES:
! 
      USE TRACER_MOD, ONLY : ITS_A_CO2_SIM 

#     include "CMN_SIZE"   ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)     :: T_SO2, T_CO2

      !=================================================================
      ! TOTAL_EMISS_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )  
 100  FORMAT( 'E D G A R   S H I P   E M I S S I O N S', / )


      ! Test for simulation type
      IF ( ITS_A_CO2_SIM() ) THEN

         !-----------------------
         ! CO2 simulation
         !-----------------------

         ! Total CO2 [Tg CO2]
         T_CO2 = SUM( CO2_SHIP ) * 1d-9

         ! Print totals
         WRITE( 6, 110 ) 'CO2 ', 2008, T_CO2,  ' CO2'

      ELSE

         !-----------------------
         ! Other simulations
         !-----------------------

         ! Total SO2 [Tg S]
         T_SO2 = SUM( SO2_SHIP ) * 1d-9 * ( 32d0 / 64d0 )

         ! Print totals in [Tg]
         WRITE( 6, 110 ) 'SO2 ', 2008, T_SO2, '[Tg S  ]'

      ENDIF

      ! Format statement
 110  FORMAT( 'ARCTAS SHIP ', a5, 
     &        'for base year ', i4, ': ', 2f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE TOTAL_EMISS_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_ARCTAS_SHIP
!
! !DESCRIPTION: Subroutine INIT\_ARCTAS\_SHIP allocates and zeroes all module 
!  arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ARCTAS_SHIP
!
! !USES:
! 
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LARCSHIP
      USE TRACER_MOD,  ONLY : ITS_A_CO2_SIM 

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
      ! INIT_ARCTAS_SHIP begins here!
      !=================================================================

      ! Allocate ANTHRO arrays if LARCTAS_SHIP is TRUE
      IF ( LARCSHIP ) THEN
      
         !--------------------------------------------------
         ! Allocate and zero arrays for SHIP emissions
         !--------------------------------------------------
         ! Test for simulation type
         IF ( ITS_A_CO2_SIM() ) THEN

            !-----------------------
            ! CO2 simulation
            !-----------------------
            ALLOCATE( CO2_SHIP( IIPAR, JJPAR ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO2_SHIP' )
            CO2_SHIP = 0d0

         ELSE

            !-----------------------
            ! Other simulations
            !-----------------------
            ALLOCATE( SO2_SHIP( IIPAR, JJPAR ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2_SHIP' )
            SO2_SHIP = 0d0

         ENDIF

         !---------------------------------------------------
         ! Allocate array for grid box surface area in cm2
         !---------------------------------------------------
         ALLOCATE( A_CM2( JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_CM2' )

         ! Fill array
         DO J = 1, JJPAR
            A_CM2(J) = GET_AREA_CM2( J )
         ENDDO

      ENDIF


      ! Return to calling program
      END SUBROUTINE INIT_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_ARCTAS_SHIP
!
! !DESCRIPTION: Subroutine CLEANUP\_ARCTAS\_SHIP deallocates all module 
!  arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ARCTAS_SHIP
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_ARCTAS_SHIP begins here!
      !=================================================================
      IF ( ALLOCATED( A_CM2        ) ) DEALLOCATE( A_CM2        )
      IF ( ALLOCATED( SO2_SHIP     ) ) DEALLOCATE( SO2_SHIP     )
      IF ( ALLOCATED( CO2_SHIP     ) ) DEALLOCATE( CO2_SHIP     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_ARCTAS_SHIP

!------------------------------------------------------------------------------

      ! End of module
      END MODULE ARCTAS_SHIP_EMISS_MOD
!EOC
