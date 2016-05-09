! $Id: bravo_mod.f,v 1.1.1.1 2009/06/09 21:51:53 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: BRAVO_MOD
!
! !DESCRIPTION: \subsection*{Overview}
!  Module BRAVO\_MOD contains variables and routines to read the BRAVO 
!  Mexican anthropogenic emission inventory for NOx, CO, and SO2. 
!  (rjp, kfb, bmy, 6/22/06, 1/30/09)
!
! \subsection*{References}
! \begin{enumerate}
! \item Kuhns, H., M. Green, and Etyemezian, V, \emph{Big Bend Regional 
!       Aerosol and Visibility Observational (BRAVO) Study Emissions 
!       Inventory}, Desert Research Institute, 2003.
! \end{enumerate}
!
! !INTERFACE: 
!
      MODULE BRAVO_MOD
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
      PUBLIC  :: CLEANUP_BRAVO
      PUBLIC  :: EMISS_BRAVO
      PUBLIC  :: GET_BRAVO_MASK
      PUBLIC  :: GET_BRAVO_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      PRIVATE :: BRAVO_SCALE_FUTURE
      PRIVATE :: INIT_BRAVO 
      PRIVATE :: READ_BRAVO_MASK
!
! !REVISION HISTORY:
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  (2 ) Now scale emissions using int-annual scale factors (amv, 08/24/07)
!  (3 ) Now accounts for FSCLYR (phs, 3/17/08)
!  (4 ) Added ProTeX headers (bmy, 1/30/09)
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
! 
      ! Arrays
      TYPE (XPLEX),  ALLOCATABLE :: BRAVO_MASK(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: BRAVO_NOx(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: BRAVO_CO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: BRAVO_SO2(:,:)

      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_BRAVO_MASK
!
! !DESCRIPTION: Function GET\_BRAVO\_MASK returns the value of the Mexico 
!  mask for BRAVO emissions at grid box (I,J).  MASK=1 if (I,J) is in the 
!  BRAVO Mexican region, or MASK=0 otherwise. (rjp, kfb, bmy, 6/22/06)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BRAVO_MASK( I, J ) RESULT( MASK )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I        ! Longitude index
      INTEGER, INTENT(IN) :: J        ! Latitude  index
!
! !RETURN VALUE:
! 
      TYPE (XPLEX)              :: MASK     ! Returns the mask value @ (I,J)
!
! !REVISION HISTORY: 
!  22 Jun 2006 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_BRAVO_MASK begins here!
      !=================================================================
      MASK = BRAVO_MASK(I,J)

      ! Return to calling program
      END FUNCTION GET_BRAVO_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_BRAVO_ANTHRO
!
! !DESCRIPTION: Function GET\_BRAVO\_ANTHRO returns the BRAVO emission 
!  for GEOS-Chem grid box (I,J) and tracer N.  Units are [molec/cm2/s]. 
!  (rjp, kfb, bmy, 6/22/06)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BRAVO_ANTHRO( I, J, N ) RESULT( BRAVO )
!
! !USES:
! 
      USE TRACERID_MOD, ONLY : IDTNOX, IDTCO, IDTSO2
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I       ! Longitude index
      INTEGER, INTENT(IN) :: J       ! Latitude index
      INTEGER, INTENT(IN) :: N       ! Tracer number
!
! RETURN VALUE:
! 
      TYPE (XPLEX)              :: BRAVO   ! Returns emissions at (I,J)
!
! !REVISION HISTORY: 
!  (1 ) added SOx, SOx ship and NH3 emissions, plus optional kg/s output
!       (amv, 06/2008)
!  (2 ) Now returns ship emissions if requested (phs, 6/08)
!  (3 ) Added checks to avoid calling unavailable ship emissions (phs, 6/08)
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_BRAVO_ANTHRO begins here!
      !=================================================================

      ! NOx
      IF ( N  == IDTNOX ) THEN
         BRAVO = BRAVO_NOx(I,J)

      ! CO
      ELSE IF ( N == IDTCO ) THEN
         BRAVO = BRAVO_CO(I,J)

      ! SO2 
      ELSE IF ( N == IDTSO2 ) THEN
         BRAVO = BRAVO_SO2(I,J)

      ! Otherwise return a negative value to indicate
      ! that there are no BRAVO emissions for tracer N
      ELSE
         BRAVO = -1d0

      ENDIF

      ! Return to calling program
      END FUNCTION GET_BRAVO_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISS_BRAVO
!
! !DESCRIPTION: Subroutine EMISS\_BRAVO reads the BRAVO emission fields at 1x1 
!  resolution and regrids them to the current model resolution. 
!  (rjp, kfb, bmy, 6/22/06, 8/9/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_BRAVO
!
! !USES:
! 
      USE BPCH2_MOD,        ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,    ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,      ONLY : LFUTURE
      USE REGRID_1x1_MOD,   ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR_1x1
      USE TIME_MOD,         ONLY : GET_YEAR

#     include "CMN_SIZE"         ! Size parameters
#     include "CMN_O3"           ! 
!
! !REVISION HISTORY: 
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: SCALEYEAR
      TYPE (XPLEX)                  :: ARRAY(I1x1,J1x1-1,1)
      TYPE (XPLEX)                  :: GEN_1x1(I1x1,J1x1-1)
      TYPE (XPLEX)                  :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                  :: SC_1x1(I1x1,J1x1)
      TYPE (XPLEX)                  :: TAU0
      CHARACTER(LEN=255)      :: FILENAME
      
      !=================================================================
      ! EMISS_BRAVO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_BRAVO
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read data from disk
      !=================================================================

      ! Use 1999 for BRAVO emission files (BASE YEAR)
      TAU0  = GET_TAU0( 1, 1, 1999 )
        
      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      !---------------------
      ! Read and regrid NOx
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.NOx.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_BRAVO: Reading ', a )
      
      ! Read NOx [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid NOx [molec/cm2/s] to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Get/Apply annual scalar factor (amv 08/21/2007)
      CALL GET_ANNUAL_SCALAR_1x1( 71, 1999, SCALEYEAR, SC_1x1 )
      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

      ! Regrid NOx [molec/cm2/s] to current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1, BRAVO_NOx )

      !---------------------
      ! Read and regrid CO
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.CO.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
      
      ! Read CO [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 4, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid CO [molec/cm2/s] to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Get/Apply annual scalar factor (amv 08/21/2007)
      CALL GET_ANNUAL_SCALAR_1x1( 72, 1999, SCALEYEAR, SC_1x1 )
      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

      ! Regrid CO [molec/cm2/s] to current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1, BRAVO_CO )

      !---------------------
      ! Read and regrid SO2
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.SO2.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
      
      ! Read SO2 [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 26, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid SO2 [molec/cm2/s] to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Get/Apply annual scalar factor (amv 08/21/2007)
      CALL GET_ANNUAL_SCALAR_1x1( 73, 1999, SCALEYEAR, SC_1x1 )
      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

      ! Regrid SO2 [molec/cm2/s] to current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1, BRAVO_SO2 )

      !=================================================================
      ! Compute IPCC future emissions (if necessary)
      !=================================================================
      IF ( LFUTURE ) THEN 
         CALL BRAVO_SCALE_FUTURE
      ENDIF

      !=================================================================
      ! Print emission totals
      !=================================================================
      CALL TOTAL_ANTHRO_TG( SCALEYEAR )

      ! Return to calling program
      END SUBROUTINE EMISS_BRAVO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BRAVO_SCALE_FUTURE
!
! !DESCRIPTION: Subroutine BRAVO\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the BRAVO anthropogenic emissions. (swu, bmy, 5/30/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BRAVO_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

#     include "CMN_SIZE"             ! Size parameters
!
! !REVISION HISTORY:
!  30 May 2006 - S. Wu & R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! BRAVO_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [molec/cm2/s]
         BRAVO_NOx(I,J) = BRAVO_NOx(I,J)                * 
     &                    GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO [molec/cm2/s]
         BRAVO_CO(I,J)  = BRAVO_CO(I,J)                 *
     &                    GET_FUTURE_SCALE_COff( I, J )

         ! Future ALK4 [atoms C/cm2/s]
         BRAVO_SO2(I,J) = BRAVO_SO2(I,J)                *
     &                    GET_FUTURE_SCALE_SO2ff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE BRAVO_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_ANTHRO_TG
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the amount of BRAVO 
!  anthropogenic emissions that are emitted each year.
!  (rjp, kfb, bmy, 6/26/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( YEAR )
!
! !USES:
! 
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTNOX, IDTCO, IDTSO2

#     include "CMN_SIZE"     ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)   :: YEAR
!
! !REVISION HISTORY: 
!  (1 ) Now YEAR is input to reflect scaling factors applied (phs, 3/17/08) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J
      TYPE (XPLEX)                :: A, B(3), NOx, CO, SO2
      CHARACTER(LEN=3)      :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'B R A V O   M E X I C A N   E M I S S I O N S', /,
     &        'Base year : 1999' )
      
      !----------------
      ! Sum emissions
      !----------------
      
      ! Define conversion factors for kg/molec
      ! (Undefined tracers will be zero)
      B(:) = 0d0
      IF ( IDTNOx > 0 ) B(1) = 1d0 / ( 6.0225d23 / 14d-3 )  ! Tg N
      IF ( IDTCO  > 0 ) B(2) = 1d0 / ( 6.0225d23 / 28d-3 )  ! Tg CO
      IF ( IDTSO2 > 0 ) B(3) = 1d0 / ( 6.0225d23 / 32d-3 )  ! Tg S

      ! Summing variables
      NOX = 0d0   
      CO  = 0d0 
      SO2 = 0d0 

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Convert [molec/cm2/s] to [Tg]
         ! (Multiply by 1d-9 to convert from [kg] to [Tg])
         A = GET_AREA_CM2( J ) * 365.25d0 * 86400d0 * 1d-9 
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Sum emissions (list NOx as Tg N)
            NOX = NOX + ( BRAVO_NOX(I,J) * A * B(1) )
            CO  = CO  + ( BRAVO_CO (I,J) * A * B(2) )
            SO2 = SO2 + ( BRAVO_SO2(I,J) * A * B(3) )
         ENDDO
      ENDDO
 
      !----------------
      ! Print sums
      !----------------

      ! Print totals in [kg]
      WRITE( 6, 110   ) 'NOx ', YEAR, NOx, ' N'
      WRITE( 6, 110   ) 'CO  ', YEAR, CO,  '  '
      WRITE( 6, 110   ) 'SO2 ', YEAR, SO2, ' S'

 110  FORMAT( 'BRAVO anthropogenic ', a4, 
     &        'for year ', i4, ': ', 2f9.4, ' Tg', a2 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_TG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READ_BRAVO_MASK
!
! !DESCRIPTION: Subroutine READ\_BRAVO\_MASK reads the Mexico mask from 
!  disk.  The Mexico mask is the fraction of the grid box (I,J) which lies 
!  w/in the BRAVO Mexican emissions region. (rjp, kfb, bmy, 6/22/06, 8/9/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BRAVO_MASK
!
! !USES:
! 
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1,   DO_REGRID_G2G_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

#     include "CMN_SIZE"       ! Size parameters
!
! !REVISION HISTORY: 
!  (1 ) Now pass UNIT to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)                  :: ARRAY(I1x1,J1x1-1,1)
      TYPE (XPLEX)                  :: GEN_1x1(I1x1,J1x1-1)
      TYPE (XPLEX)                  :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                  :: XTAU
      CHARACTER(LEN=255)      :: FILENAME 

      !=================================================================
      ! READ_BRAVO_MASK begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'BRAVO_200607/BRAVO.MexicoMask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_MEXICO_MASK: Reading ', a )

      ! Get TAU0 for Jan 1985
      XTAU  = GET_TAU0( 1, 1, 1999 )

      ! Mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 XTAU,      I1x1,     J1x1-1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      GEN_1x1(:,:) = ARRAY(:,:,1) 

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )
      
      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, BRAVO_MASK )

      ! Return to calling program
      END SUBROUTINE READ_BRAVO_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_BRAVO
!
! !DESCRIPTION: Subroutine INIT\_BRAVO allocates and zeroes BRAVO module 
!  arrays, and also creates the mask which defines the Mexico region 
!  (rjp, kfb, bmy, 6/26/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_BRAVO
!
! !USES:
! 
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_XMID, GET_YMID
      USE LOGICAL_MOD, ONLY : LBRAVO

#     include "CMN_SIZE"    ! Size parameters
!
! !REVISION HISTORY: 
!  18 Oct 2006 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS

      !=================================================================
      ! INIT_BRAVO begins here!
      !=================================================================

      ! Return if LBRAVO is false
      IF ( .not. LBRAVO ) RETURN
      
      !--------------------------
      ! Allocate and zero arrays
      !--------------------------

      ALLOCATE( BRAVO_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_NOx' )
      BRAVO_NOx = 0d0

      ALLOCATE( BRAVO_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_CO' )
      BRAVO_CO = 0d0

      ALLOCATE( BRAVO_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_SO2' )
      BRAVO_SO2 = 0d0

      !--------------------------
      ! Read Mexico mask
      !--------------------------
     
      ALLOCATE( BRAVO_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_MASK' )
      BRAVO_MASK = 0d0
      
      ! Read the mask
      CALL READ_BRAVO_MASK

      ! Return to calling program
      END SUBROUTINE INIT_BRAVO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_BRAVO
!
! !DESCRIPTION: Subroutine CLEANUP\_BRAVO deallocates all BRAVO module arrays.
!  (rjp, kfb, bmy, 6/26/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_BRAVO
!
! !REVISION HISTORY: 
!  1 Nov 2005 - R. Yantosca - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_BRAVO begins here!
      !=================================================================
      IF ( ALLOCATED( BRAVO_NOx  ) ) DEALLOCATE( BRAVO_NOx  )
      IF ( ALLOCATED( BRAVO_CO   ) ) DEALLOCATE( BRAVO_CO   )
      IF ( ALLOCATED( BRAVO_SO2  ) ) DEALLOCATE( BRAVO_SO2  )
      IF ( ALLOCATED( BRAVO_MASK ) ) DEALLOCATE( BRAVO_MASK )

      ! Return to calling program
      END SUBROUTINE CLEANUP_BRAVO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE BRAVO_MOD
!EOC
