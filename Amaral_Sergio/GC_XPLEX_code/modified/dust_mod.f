! $Id: dust_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE DUST_MOD
!
!******************************************************************************
!  Module DUST_MOD contains routines for computing dust aerosol emissions,
!  chemistry, and optical depths. (rjp, tdf, bmy, 4/14/04, 1/23/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) FRAC_S   (TYPE (XPLEX) ) : Fraction of each size classes (GINOUX only)
!  (2 ) DUSTREFF (TYPE (XPLEX) ) : Dust particle radii [m]
!  (3 ) DUSTDEN  (TYPE (XPLEX) ) : Soil density [kg/m3]
!  (4 ) IDDEP    (INTEGER) : Dust ID flags for drydep
!  (5 ) DRYDST1  (INTEGER) : Index for DST1 in drydep array
!  (6 ) DRYDST2  (INTEGER) : Index for DST2 in drydep array
!  (7 ) DRYDST3  (INTEGER) : Index for DST3 in drydep array
!  (8 ) DRYDST4  (INTEGER) : Index for DST4 in drydep array
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMDUST           : Driver routine for dust chemistry
!  (2 ) DRY_SETTLING       : Routine which performs dust settling
!  (3 ) DRY_DEPOSITION     : Routine which performs dust dry deposition
!  (4 ) EMISSDUST          : Driver routine for dust emission
!  (5 ) SRC_DUST_DEAD      : Dust emissions according to DEAD   source function
!  (6 ) SRC_DUST_GINOUX    : Dust emissions according to GINOUX source function
!  (7 ) RDUST_ONLINE       : Computes dust optical depths (online dust)
!  (8 ) RDUST_OFFLINE      : Computes dust optical depths (monthly mean dust)
!  (9 ) GET_SCALE_GROUP    : Return index of emission group 
!  (10) INIT_DUST          : Allocates & initializes all module variables
!  (11) CLEANUP_DUST       : Deallocates all module variables
!
!  GEOS-CHEM modules referenced by "dust_mod.f"
!  ============================================================================
!  (1 ) dao_mod.f          : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f         : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) directory_mod.f    : Module containing GEOS-CHEM data & met field dirs
!  (4 ) drydep_mod.f       : Module containing GEOS-CHEM drydep routines
!  (5 ) dust_dead_mod.f    : Module containing Zender's DEAD dust routines
!  (6 ) error_mod.f        : Module containing I/O error and NaN check routines
!  (7 ) file_mod.f         : Contains file unit numbers and error checks
!  (8 ) grid_mod.f         : Module containing horizontal grid information
!  (9 ) logical_mod.f      : Module containing GEOS-CHEM logical switches
!  (10) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!  (11) time_mod.f         : Module containing routines for computing time/date
!  (12) tracer_mod.f       : Module containing GEOS-CHEM tracer array STT etc.
!  (13) tracerid_mod.f     : Module containing pointers to tracers & emissions
!  
!  NOTES:
!  (1 ) Bug fix in SRC_DUST_DEAD (bmy, 4/14/04)
!  (2 ) Now references "logical_mod.f", "directory_mod.f", and "tracer_mod.f"
!        Added comments. (bmy, 7/2/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Bug fix in snow height computation (bmy, 11/18/05)
!  (6 ) Now only do drydep if LDRYD=T (bmy, 5/23/06)
!  (7 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (8 ) Updated output print statement in SRC_DUST_DEAD (bmy, 1/23/07)
!  (9 ) Modifications for GEOS-5 (bmy, 1/24/07)
!  (10) Modified to archive only hydrophilic aerosol/aqueous dust surface area
!        (excluding BCPO and OCPO) for aqueous chemistry calculations 
!        Dust surfaces are considered aqueous only when RH > 35% (tmf, 3/6/09)
!  (11) Updated for dust emission adjoint (xxu, dkh, 01/13/12, adj32_011) 
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "dust_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE                :: DRYDST1, DRYDST2,  DRYDST3, DRYDST4 
      ! Make some public for use with adjoint (dkh, 01/10/12, adj32_011) 
      !PRIVATE                :: DUSTDEN, DUSTREFF, FRAC_S,  IDDEP
      PRIVATE                :: FRAC_S

      ! PRIVATE module routines
      PRIVATE                :: DRY_SETTLING
      PRIVATE                :: DRY_DEPOSITION
      PRIVATE                :: SRC_DUST_DEAD
      PRIVATE                :: SRC_DUST_GINOUX

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER                :: DRYDST1, DRYDST2, DRYDST3, DRYDST4
      INTEGER, ALLOCATABLE   :: IDDEP(:)
      TYPE (XPLEX),  ALLOCATABLE   :: FRAC_S(:)
      TYPE (XPLEX),  ALLOCATABLE   :: DUSTREFF(:)
      TYPE (XPLEX),  ALLOCATABLE   :: DUSTDEN(:)
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!-----------------------------------------------------------------------

      SUBROUTINE CHEMDUST
!
!******************************************************************************
!  Subroutine CHEMDUST is the interface between the GEOS-CHEM main program and
!  the dust chemistry routines that mostly calculates dust dry deposition.
!  (tdf, bmy, 3/30/04, 5/23/06)
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" and LDUST from "logical_mod.f"
!        (bmy, 7/20/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now only do dry deposition if LDRYD = T (bmy, 5/23/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LDRYD,   LDUST
      USE DRYDEP_MOD,   ONLY : DEPNAME, NUMDEP
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: N

      !=================================================================
      ! CHEMDUST begins here!
      !=================================================================

      ! Execute on first call only
      IF ( FIRST ) THEN
 
         ! Stop w/ error if dust tracer flags are undefined
         IF ( IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4 == 0 ) THEN
            IF ( LDUST ) THEN 
               CALL ERROR_STOP( 
     &              'LDUST=T but dust tracers are undefined!',
     &              'EMISSDUST ("dust_mod.f")' )
            ENDIF
         ENDIF

         ! Allocate arrays (if necessary)
         CALL INIT_DUST

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'DST1' )
                  DRYDST1 = N
               CASE ( 'DST2' )
                  DRYDST2 = N
               CASE ( 'DST3' )
                  DRYDST3 = N
               CASE ( 'DST4' )
                  DRYDST4 = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! This may lead to out of bounds errors
         IDDEP(1) = DRYDST1
         IDDEP(2) = DRYDST2
         IDDEP(3) = DRYDST3
         IDDEP(4) = DRYDST4
 
         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Do dust settling & deposition
      !=================================================================

      ! Dust settling
      CALL DRY_SETTLING(   STT(:,:,:,IDTDST1:IDTDST4) )

      ! Dust deposition
      IF ( LDRYD ) THEN
         CALL DRY_DEPOSITION( STT(:,:,:,IDTDST1:IDTDST4) )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMDUST

!------------------------------------------------------------------------------

      SUBROUTINE DRY_SETTLING( TC )
!
!******************************************************************************
!  Subroutine DRY_SETTLING computes the dry settling of dust tracers.
!  (tdf, bmy, 3/30/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Dust tracer array 
!
!  NOTES
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!  (2 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
! 
      USE DAO_MOD,      ONLY : T, BXHEIGHT
      USE DIAG_MOD,     ONLY : AD44
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTDST1

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! g0
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! Local variables
      INTEGER               :: I, J, L, N
      TYPE (XPLEX)                :: DT_SETTL, DELZ,  DELZ1
      TYPE (XPLEX)                :: REFF,     DEN,   CONST   
      TYPE (XPLEX)                :: NUM,      LAMDA, FLUX
      TYPE (XPLEX)                :: AREA_CM2, TC0(LLPAR)
      TYPE (XPLEX)                :: TOT1,     TOT2

      ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa      
      TYPE (XPLEX)                :: P 

      ! Diameter of aerosol [um]
      TYPE (XPLEX)                :: Dp

      ! Pressure * DP
      TYPE (XPLEX)                :: PDp 

      ! Temperature (K)    
      TYPE (XPLEX)                :: TEMP        

      ! Slip correction factor
      TYPE (XPLEX)                :: Slip        

      ! Viscosity of air (Pa s)
      TYPE (XPLEX)                :: Visc   

      ! Settling velocity of particle (m/s)
      TYPE (XPLEX)                :: VTS(LLPAR)  
      
      ! Parameters
      TYPE (XPLEX),  PARAMETER    :: C1 = xplex( 0.7674D0,0d0)
      TYPE (XPLEX),  PARAMETER    :: C2 = xplex( 3.079d0,0d0) 
      TYPE (XPLEX),  PARAMETER    :: C3 = xplex( 2.573D-11,0d0)
      TYPE (XPLEX),  PARAMETER    :: C4 = xplex(-1.424d0,0d0)

      !=================================================================
      ! DRY_SETTLING begins here!
      !=================================================================

      ! Dust settling timestep [s]
      DT_SETTL = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,        L,    N,     DEN,  REFF, DP    )
!$OMP+PRIVATE( CONST, AREA_CM2, VTS,  TEMP,  P,    PDP,  SLIP  )
!$OMP+PRIVATE( VISC,  TC0,      DELZ, DELZ1, TOT1, TOT2, FLUX  )

      ! Loop over dust bins
      DO N = 1, NDSTBIN

         ! Initialize
         DEN   = DUSTDEN(N)
         REFF  = DUSTREFF(N)
         DP    = 2D0 * REFF * 1.D6              ! Dp [um] = particle diameter
         CONST = 2D0 * DEN * REFF**2 * G0 / 9D0

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2(J)

            ! Loop over longitudes
            DO I = 1, IIPAR
            
               ! Initialize settling velocity
               DO L = 1, LLPAR
                  VTS(L) = 0d0
               ENDDO

               ! Loop over levels
               DO L = 1, LLPAR

                  ! Get P [kPa], T [K], and P*DP
                  P    = GET_PCENTER(I,J,L) * 0.1d0
                  TEMP = T(I,J,L)
                  PDP  = P * DP

                  !=====================================================
                  ! # air molecule number density
                  ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
                  !
                  ! # gas mean free path
                  ! lamda = 1.d6 / 
                  !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
                  !
                  ! # Slip correction
                  ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * 
                  !      &  exp( -1.1 * Dp / (2. * lamda))) / Dp
                  !=====================================================
                  ! NOTE, Slip correction factor calculations following 
                  !       Seinfeld, pp464 which is thought to be more 
                  !       accurate but more computation required.
                  !=====================================================

                  ! Slip correction factor as function of (P*dp)
                  SLIP = 1d0 + 
     &                   ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP
            
                  !=====================================================
                  ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                  ! which produce slip correction factor with small 
                  ! error compared to the above with less computation.
                  !=====================================================

                  ! Viscosity [Pa s] of air as a function of temp (K)
                  VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

                  ! Settling velocity [m/s]
                  VTS(L) = CONST * SLIP / VISC

               ENDDO

               ! Method is to solve bidiagonal matrix 
               ! which is implicit and first order accurate in Z
               DO L = 1, LLPAR
                  TC0(L) = TC(I,J,L,N)
               ENDDO

               ! We know the boundary condition at the model top
               L           = LLTROP
               DELZ        = BXHEIGHT(I,J,L)
               TC(I,J,L,N) = TC(I,J,L,N) / 
     &                       ( 1.d0 + DT_SETTL * VTS(L) / DELZ )

               DO L = LLTROP-1, 1, -1
                  DELZ        = BXHEIGHT(I,J,L)
                  DELZ1       = BXHEIGHT(I,J,L+1)
                  TC(I,J,L,N) = 1.d0 / 
     &                          ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
     &                 * (TC(I,J,L,N)  + DT_SETTL * VTS(L+1) / DELZ1
     &                 *  TC(I,J,L+1,N) )
               ENDDO

               !========================================================      
               ! ND44: Dry deposition diagnostic [#/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN

                  ! Initialize
                  TOT1 = 0d0
                  TOT2 = 0d0
            
                  ! Compute column totals of TCO(:) and TC(I,J,:,N)
                  DO L = 1, LLPAR
                     TOT1 = TOT1 + TC0(L)
                     TOT2 = TOT2 + TC(I,J,L,N)
                  ENDDO

                  ! Convert dust flux from [kg/s] to [#/cm2/s]
                  FLUX = ( TOT1 - TOT2 ) / DT_SETTL  
                  FLUX = FLUX * XNUMOL(IDTDST1) / AREA_CM2 

                  ! Save in AD44
                  AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRY_SETTLING

!------------------------------------------------------------------------------

      SUBROUTINE DRY_DEPOSITION( TC )
!
!******************************************************************************
!  Subroutine DRY_DEPOSITION computes the loss of dust due to dry deposition
!  at the surface using an implicit method. (tdf, bmy, 3/30/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Dust tracer array   
!
!  NOTES: 
!  (1 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV 
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTDST1

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! local variables
      INTEGER               :: I,   J,   L,      N
      TYPE (XPLEX)                :: OLD, NEW, DTCHEM, FLUX, AREA_CM2

      !=================================================================
      ! DRY_DEPOSITION begins here!
      !=================================================================

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over dust bins
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, OLD, NEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over dust bins
      DO N = 1, NDSTBIN    

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2(J)

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Original dust concentration at surface
               OLD = TC(I,J,1,N)

               ! Dust left after dry deposition
               NEW = OLD * EXP( -DEPSAV(I,J,IDDEP(N)) * DTCHEM  )

               !========================================================
               ! ND44 diagnostic: dust drydep loss [#/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN

                  ! Convert drydep flux from [kg/s] to [#/cm2/s]
                  FLUX     = ( OLD - NEW ) / DTCHEM 
                  FLUX     = FLUX * XNUMOL(IDTDST1) / AREA_CM2 
               
                  ! Store in AD44
                  AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
               ENDIF

               ! Save back into STT
               TC(I,J,1,N) = NEW 
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRY_DEPOSITION

!------------------------------------------------------------------------------

      SUBROUTINE EMISSDUST
!
!******************************************************************************
!  Subroutine EMISSDUST is the driver routine for the dust emission
!  module.  You may call either the GINOUX or the DEAD dust source 
!  function. (tdf, bmy, 3/30/04, 10/3/05)
!
!  NOTES:
!  (1 ) Now reference LDEAD, LDUST, LPRT from "logical_mod.f".  Now reference!
!        STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F(0 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP, DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LDEAD, LDUST, LPRT
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4

#     include "CMN_SIZE"     ! Size parameters
      
      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      
      !=================================================================
      ! EMISSDUST begins here!
      !=================================================================

      ! Execute on first-call only
      IF ( FIRST ) THEN

         ! Return if dust ID flags are not defined
         IF ( IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4 == 0 ) THEN
            IF ( LDUST ) THEN 
               CALL ERROR_STOP( 
     &              'LDUST=T but dust tracers are undefined!',
     &              'EMISSDUST ("dust_mod.f")' )
            ENDIF
         ENDIF
          
         ! Allocate module arrays
         CALL INIT_DUST
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a INIT_DUST' )

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call appropriate emissions routine
      !=================================================================
      IF ( LDEAD ) THEN

         ! Use Zender's DEAD dust source function
         CALL SRC_DUST_DEAD( STT(:,:,:,IDTDST1:IDTDST4) )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_DEAD' )
         
      ELSE

         ! Use Paul Ginoux's dust source function
         CALL SRC_DUST_GINOUX( STT(:,:,:,IDTDST1:IDTDST4) )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_GINOUX')

      ENDIF     
      ! Return to calling program
      END SUBROUTINE EMISSDUST

!------------------------------------------------------------------------------

      SUBROUTINE SRC_DUST_DEAD( TC )
!
!******************************************************************************
!  Subroutine SRC_DUST_DEAD is the DEAD model dust emission scheme, 
!  alternative to Ginoux scheme.  Increments the TC array with emissions 
!  from the DEAD model.  (tdf, bmy, 4/8/04, 1/23/07)
!
!  Input:
!         SRCE_FUNK Source function                               (-)
!         for 1: Sand, 2: Silt, 3: Clay
!         DUSTDEN   Dust density                                  (kg/m3)
!         DUSTREFF  Effective radius                              (um)
!         AD        Air mass for each grid box                    (kg)
!         NTDT      Time step                                     (s)
!         W10M      Velocity at the anemometer level (10meters)   (m/s)
!         GWET      Surface wetness                               (-)
!
!  Parameters used in GEOS-CHEM
!
!  Longitude: IIPAR 
!  Latitude : JJPAR 
!  Levels   : LLPAR  =  20 (GEOS-1), 26 (GEOS-strat), 30 (GEOS-terra)
!  Size bins: NDSTBIN  =   4
!
!  Dust properties used in GOCART
!
!  Size classes: 01-1, 1-1.8, 1.8-3, 3-6 (um)
!  Radius: 0.7, 1.5, 2.5, 4  (um)
!  Density: 2500, 2650, 2650, 2650 (kg/m3)
!
!  NOTES:
!  (1 ) Added OpenMP parallelization, added comments (bmy, 4/8/04)
!  (2 ) Bug fix: DSRC needs to be held PRIVATE (bmy, 4/14/04)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Bug fix: It should be SNOW/1d3 not SNOW*1d3 (tdf, bmy, 11/18/05)
!  (6 ) Updated output statement (bmy, 1/23/07)
!  (7 ) Use SNOMAS (m H2O) for GEOS-5 (bmy, 1/24/07)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : BXHEIGHT,     GWETTOP,   LWI
      USE DAO_MOD,       ONLY : SNOW,         SPHU,      T    
      USE DAO_MOD,       ONLY : TS,           UWND,      VWND
      USE DAO_MOD,       ONLY : SNOMAS
      USE DUST_DEAD_MOD, ONLY : GET_TIME_INVARIANT_DATA, GET_ORO
      USE DUST_DEAD_MOD, ONLY : GET_MONTHLY_DATA,        DST_MBL
      USE DIAG_MOD,      ONLY : AD06
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE GRID_MOD,      ONLY : GET_YMID_R
      USE PRESSURE_MOD,  ONLY : GET_PEDGE,       GET_PCENTER 
      USE TIME_MOD,      ONLY : GET_TS_EMIS,     GET_MONTH   
      USE TIME_MOD,      ONLY : GET_DAY_OF_YEAR, ITS_A_NEW_MONTH
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      ! adj_group: add for emissions scale factors (xxu, 11/02/10, adj32_011)
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EDST1, IDADJ_EDST2
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EDST3, IDADJ_EDST4
      USE ADJ_ARRAYS_MOD,    ONLY : IS_DUST_EMS_ADJ
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND06
#     include "CMN_GCTM"      ! g0

      !----------------
      ! Arguments
      !----------------
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      !-----------------
      ! Local variables
      !-----------------

      ! Scalars
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,      J,      L,       N
      INTEGER                :: M,      IOS,    INC,     LAT_IDX
      INTEGER                :: NDB,    NSTEP
      INTEGER                :: MM
      TYPE (XPLEX)                 :: W10M,   DEN,    DIAM,    U_TS0 
      TYPE (XPLEX)                 :: U_TS,   SRCE_P, Reynol,  YMID_R
      TYPE (XPLEX)                 :: ALPHA,  BETA,   GAMMA,   CW
      TYPE (XPLEX)                 :: DTSRCE, XTAU,   P1,      P2
      TYPE (XPLEX)                 :: DOY             
      CHARACTER(LEN=255)     :: FILENAME

      ! Arrays
      INTEGER                :: OROGRAPHY(IIPAR,JJPAR)
      TYPE (XPLEX)                 :: PSLON(IIPAR)         ! surface pressure
      TYPE (XPLEX)                 :: PTHICK(IIPAR)        ! delta P (L=1)
      TYPE (XPLEX)                 :: PMID(IIPAR)          ! mid layer P (L=1)
      TYPE (XPLEX)                 :: TLON(IIPAR)          ! temperature (L=1)
      TYPE (XPLEX)                 :: THLON(IIPAR)         ! pot. temp. (L=1)
      TYPE (XPLEX)                 :: ULON(IIPAR)          ! U-wind (L=1)
      TYPE (XPLEX)                 :: VLON(IIPAR)          ! V-wind (L=1)
      TYPE (XPLEX)                 :: BHT2(IIPAR)          ! half box height (L=1)
      TYPE (XPLEX)                 :: Q_H2O(IIPAR)         ! specific humidity (L=1)
      TYPE (XPLEX)                 :: ORO(IIPAR)           ! "orography" 
      TYPE (XPLEX)                 :: SNW_HGT_LQD(IIPAR)   ! equivalent snow ht.
      TYPE (XPLEX)                 :: DSRC(IIPAR,NDSTBIN)  ! dust mixing ratio incr.

      !----------------
      ! Parameters
      !----------------
      TYPE (XPLEX), PARAMETER      :: Ch_dust = xplex(9.375d-10,0d0)
      TYPE (XPLEX), PARAMETER      :: G       = xplex(g0%r * 1.D2,0d0)
      TYPE (XPLEX), PARAMETER      :: RHOA    = xplex(1.25D-3,0d0)
      TYPE (XPLEX), PARAMETER      :: CP      = xplex(1004.16d0,0d0)
      TYPE (XPLEX), PARAMETER      :: RGAS= xplex(8314.3d0/28.97d0,0d0)
      TYPE (XPLEX), PARAMETER      :: AKAP    = xplex(RGAS%r/CP%r,0d0)
      TYPE (XPLEX), PARAMETER      :: P1000   = xplex(1000d0,0d0)

      ! External functions
      TYPE (XPLEX),  EXTERNAL       :: SFCWINDSQR

      !=================================================================
      ! SRC_DUST_DEAD begins here!
      !=================================================================      

      ! DTSRCE is the emission timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

      ! DOY is the day of year (0-365 or 0-366)
      DOY    = XPLX( GET_DAY_OF_YEAR() )

      !=================================================================
      ! Read data fields for the DEAD model from disk
      !=================================================================      
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   ) 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       
         ! FORMAT strings
 100     FORMAT( 'D E A D   D U S T   M O B I L I Z A T I O N'         )
 110     FORMAT( 'Routines from DEAD model by Charlie Zender et al'    )
 120     FORMAT( 'Modified for GEOS-CHEM by D. Fairlie and R. Yantosca')
 130     FORMAT( 'Last Modification Date: 1/23/07'                     )

         ! Read fields for DEAD that are time-invariant
         CALL GET_TIME_INVARIANT_DATA
         
         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Read monthly data for DEAD
      IF ( ITS_A_NEW_MONTH() ) THEN
         CALL GET_MONTHLY_DATA
      ENDIF

      ! Determine group (temporal) (xxu, dkh, 01/13/12, adj32_011) 
      MM = GET_SCALE_GROUP()
      ! Print out scaling info
      WRITE(6,*) '- READ / RESCALE DUST EMISSIONS: use SCALE_GROUP ', MM

      !=================================================================
      ! Call dust mobilization scheme
      !=================================================================

      ! Make OROGRAPHY array from GEOS-CHEM LWI
      CALL GET_ORO( OROGRAPHY )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,      P1,   P2,   PTHICK, PMID, TLON        )
!$OMP+PRIVATE( THLON, ULON,   VLON, BHT2, Q_H2O,  ORO,  SNW_HGT_LQD )
!$OMP+PRIVATE( N,     YMID_R, DSRC                                  )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Pressure [hPa] at bottom and top edge of level 1
            P1             = GET_PEDGE(I,J,1)
            P2             = GET_PEDGE(I,J,2)

            ! Pressure thickness of 1st layer [Pa]
            PTHICK(I)      = ( P1 - P2 ) * 100d0

            ! Pressure at midpt of surface layer [Pa]
            PMID(I)        = GET_PCENTER(I,J,1) * 100d0 

            ! Temperature [K] at surface layer
            TLON(I)        = T(I,J,1)

            ! Potential temperature [K] at surface layer
            THLON(I)       = TLON(I) * ( P1000 / PMID(I) )**AKAP

            ! U and V winds at surface [m/s]
            ULON(I)        = UWND(I,J,1)
            VLON(I)        = VWND(I,J,1)
            
            ! Half box height at surface [m]
            BHT2(I)        = BXHEIGHT(I,J,1) / 2.d0

            ! Specific humidity at surface [kg H2O/kg air]
            Q_H2O(I)       = SPHU(I,J,1) / 1000.d0         

            ! Orography at surface
            ! Ocean is 0; land is 1; ice is 2
            ORO(I)         = OROGRAPHY(I,J) 

            ! Snow height [m H2O]
#if   defined( GEOS_5  )
            SNW_HGT_LQD(I) = SNOMAS(I,J)
#else
            SNW_HGT_LQD(I) = SNOW(I,J) / 1000d0
#endif
            ! Dust tracer and increments
            DO N = 1, NDSTBIN
               DSRC(I,N)   = 0.0d0
            ENDDO
         ENDDO 

         !==============================================================
         ! Call dust mobilization driver (DST_MBL) for latitude J
         !==============================================================
         
         ! Latitude in RADIANS
         YMID_R = GET_YMID_R(J)

         ! Call DEAD dust mobilization
         CALL DST_MBL( DOY,    BHT2,  J,     YMID_R, ORO,    
     &                 PTHICK, PMID,  Q_H2O, DSRC,   SNW_HGT_LQD, 
     &                 DTSRCE, TLON,  THLON, VLON,   ULON,
     &                 FIRST,  J )

         ! Update
         DO N = 1, NDSTBIN
         DO I = 1, IIPAR

            ! Include dust adjoint scale factor (xxu, 11/02/10, adj32_011)
            IF ( LADJ_EMS .and. IS_DUST_EMS_ADJ) THEN

               IF(N==1) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST1)
               IF(N==2) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST2)
               IF(N==3) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST3)
               IF(N==4) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST4)

            ENDIF

            ! Add dust emissions into tracer array [kg]
            TC(I,J,1,N) = TC(I,J,1,N) + DSRC(I,N) 
 
            ! ND19 diagnostics [kg]
            IF ( ND06 > 0 ) THEN
               AD06(I,J,N) = AD06(I,J,N) + DSRC(I,N)
            ENDIF
         ENDDO
         ENDDO

      ENDDO  
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRC_DUST_DEAD

!------------------------------------------------------------------------------

      SUBROUTINE SRC_DUST_GINOUX( TC )
!
!******************************************************************************
!  Paul GINOUX dust source function 
!  (Added to GEOS-CHEM, tdf, bmy, 4/8/04, 7/20/04)
!
!  This subroutine updates the surface mixing ratio of dust aerosols for
!  NDSTBIN size bins. The uplifting of dust depends in space on the source 
!  function, and in time and space on the soil moisture and surface
!  wind speed (10 meters). Dust is uplifted if the wind speed is greater
!  than a threshold velocity which is calculated with the formula of
!  Marticorena et al.  (JGR, v.102, p 23277-23287, 1997).
!  To run this subroutine you need the source function which can be
!  obtained by contacting Paul Ginoux at ginoux@rondo.gsfc.nasa.gov
!  If you are not using GEOS DAS met fields, you will most likely need
!  to adapt the adjusting parameter.
! 
!  Contact: Paul Ginoux (ginoux@rondo.gsfc.nasa.gov) 
! 
! 
!  Input:
!    SRCE_FUNK Source function                               (-)
!              for 1: Sand, 2: Silt, 3: Clay
!
!    DUSTDEN   Dust density                                  (kg/m3)
!    DUSTREFF  Effective radius                              (um)
!    AD        Air mass for each grid box                    (kg)
!    NTDT      Time step                                     (s)
!    W10m      Velocity at the anemometer level (10meters)   (m/s)
!    GWET      Surface wetness                               (-)
!       
! 
!  Parameters used in GEOS-CHEM
! 
!  Longitude: IIPAR 
!  Latitude : JJPAR 
!  Levels   : LLPAR  =  20 (GEOS-1), 26 (GEOS-strat), 30 (GEOS-terra)
!  Size bins: NDSTBIN  =   4
! 
!  Dust properties used in GOCART
! 
!  Size classes: 01-1, 1-1.8, 1.8-3, 3-6 (um)
!  Radius: 0.7, 1.5, 2.5, 4  (um)
!  Density: 2500, 2650, 2650, 2650 (kg/m3)
! 
!  References:
!  ============================================================================
!  (1 ) Ginoux, P., M. Chin, I. Tegen, J. Prospero, B. Hoben, O. Dubovik,
!        and S.-J. Lin, "Sources and distributions of dust aerosols simulated
!        with the GOCART model", J. Geophys. Res., 2001
!  (2 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!
!  NOTES:
!  (1 ) Added OpenMP parallelization (bmy, 4/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DAO_MOD,       ONLY : GWETTOP   
      USE DIAG_MOD,      ONLY : AD06
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE GRID_MOD,      ONLY : GET_AREA_M2

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! ND19, LD13 (for now)
#     include "CMN_GCTM"  ! g0

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, N, M, IOS
      INTEGER                :: IPOINT(NDSTBIN) = (/3, 2, 2, 2/)
      TYPE (XPLEX)                 :: ARRAY(IIPAR,JJPAR,3)
      TYPE (XPLEX), SAVE           :: SRCE_FUNC(IIPAR,JJPAR,3)
      TYPE (XPLEX)                 :: W10M,   DEN,  DIAM,   U_TS0, U_TS
      TYPE (XPLEX)                 :: SRCE_P, DSRC, REYNOL, ALPHA, BETA
      TYPE (XPLEX)                 :: GAMMA,  CW,   DTSRCE, AREA_M2
      CHARACTER(LEN=255)     :: FILENAME

      ! Transfer coeff for type natural source  (kg*s2/m5)
      TYPE (XPLEX), PARAMETER      :: CH_DUST  = xplex(9.375d-10,0d0)
      TYPE (XPLEX), PARAMETER      :: G        = xplex(G0%r * 1.d2,0d0)
      TYPE (XPLEX), PARAMETER      :: RHOA     = xplex(1.25d-3,0d0)

      ! External functions
      TYPE (XPLEX),  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRC_DUST_GINOUX begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Read dust source function
      !=================================================================
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   ) 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       
         ! FORMAT strings
 100     FORMAT( 'G I N O U X   D U S T   M O B I L I Z A T I O N'     )
 110     FORMAT( 'Routines originally by Paul Ginoux, GSFC'            )
 120     FORMAT( 'Modified for GEOS-CHEM by D. Fairlie and R. Yantosca')
 130     FORMAT( 'Last Modification Date: 4/6/04'                      )

         ! Filename
         FILENAME = TRIM( DATA_DIR )        // 
     &              'dust_200203/NSP.dust.' // GET_RES_EXT()

         ! Open file
         OPEN( 65, FILE=FILENAME,       STATUS='OLD', 
     &             FORM= 'UNFORMATTED', IOSTAT=IOS )
         IF ( IOS > 0 ) CALL IOERROR( IOS, 65, 'SRC_DUST_GINOUX:1' )

         ! Read data
         READ( 65, IOSTAT=IOS ) ARRAY
         IF ( IOS > 0 ) CALL IOERROR( IOS, 65, 'SRC_DUST_GINOUX:2' )

         ! Close file
         CLOSE( 65 )

         ! Cast to TYPE (XPLEX)
         SRCE_FUNC = ARRAY

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF
        
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,     M,    N,      DEN,   DIAM    )
!$OMP+PRIVATE( REYNOL, ALPHA, BETA, GAMMA,  U_TS0, AREA_M2 )
!$OMP+PRIVATE( CW,     U_TS,  W10M, SRCE_P, DSRC           )

      ! Loop over size bins
      DO N = 1, NDSTBIN

         !==============================================================
         ! Threshold velocity as a function of the dust density and the 
         ! diameter from Bagnold (1941), valid for particles larger 
         ! than 10 um.
         ! 
         ! u_ts0 = 6.5*sqrt(dustden(n)*g0*2.*dustreff(n))
         !
         ! Threshold velocity from Marticorena and Bergametti
         ! Convert units to fit dimensional parameters
         !==============================================================
         DEN    = DUSTDEN(N) * 1.d-3                 ! [g/cm3]
         DIAM   = 2d0 * DUSTREFF(N) * 1.d2           ! [cm in diameter]
         REYNOL = 1331.d0 * DIAM**(1.56d0) + 0.38d0  ! [Reynolds number]
         ALPHA  = DEN * G * DIAM / RHOA
         BETA   = 1d0 + ( 6.d-3 / ( DEN * G * DIAM**(2.5d0) ) )
         GAMMA  = ( 1.928d0 * REYNOL**(0.092d0) ) - 1.d0
         
         !==============================================================
         ! I think the 129.d-5 is to put U_TS in m/sec instead of cm/sec
         ! This is a threshold friction velocity!       from M&B
         ! i.e. Ginoux uses the Gillette and Passi formulation
         ! but has substituted Bagnold's Ut with M&B's U*t.
         ! This appears to be a problem.  (tdf, 4/2/04)
         !==============================================================

         ! [m/s] 
         U_TS0  = 129.d-5 * SQRT( ALPHA ) * SQRT( BETA ) / SQRT( GAMMA )
         M      = IPOINT(N)

         ! Loop over latitudes
         DO J = 1, JJPAR
            
            ! Get grid box surface area [m2]
            AREA_M2 = GET_AREA_M2(J)

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Fraction of emerged surfaces 
               ! (subtract lakes, coastal ocean,...)
               CW = 1.d0

               ! Case of surface dry enough to erode
               IF ( GWETTOP(I,J) < 0.2d0 ) THEN
               
                  U_TS = U_TS0*( 1.2d0 + 
     &                           0.2d0*LOG10( MAX(1.d-3,GWETTOP(I,J))))
                  U_TS = MAX( 0.d0, U_TS )

               ELSE

                  ! Case of wet surface, no erosion
                  U_TS = 100.d0

               ENDIF

               ! 10m wind speed [m/s]
               W10M   = SQRT( SFCWINDSQR(I,J) )

               ! Units are m2
               SRCE_P = FRAC_S(N) * SRCE_FUNC(I,J,M) * AREA_M2

               ! Dust source increment [kg]
               DSRC   = CW     * CH_DUST * SRCE_P * W10M**2
     &                * ( W10M - U_TS )  * DTSRCE  

               ! Not less than zero
               IF ( DSRC < 0.d0 ) DSRC = 0.d0

               ! Dust SOURCE at first model level [kg].
               TC(I,J,1,N) = TC(I,J,1,N) + DSRC 
            
               !========================================================
               ! ND06 diagnostics: dust emissions [kg/timestep]
               !========================================================
               IF ( ND06 > 0 ) THEN
                  AD06(I,J,N) = AD06(I,J,N) + DSRC
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRC_DUST_GINOUX

!------------------------------------------------------------------------------

      SUBROUTINE RDUST_ONLINE( DUST )
!
!******************************************************************************
!  Subroutine RDUST reads global mineral dust concentrations as determined 
!  by P. Ginoux.  Calculates dust optical depth at each level for the
!  FAST-J routine "set_prof.f". (rvm, rjp, tdf, bmy, 4/1/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DUST (TYPE (XPLEX)) : Dust from soils [kg/m3]
!
!  NOTES:
!  (1 ) Bundled into "dust_mod.f" (bmy, 4/1/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f".  Now parallelize over
!        the L-dimension for ND21 diagnostics. (bmy, 7/20/04)
!  (3 ) Archive only hydrophilic aerosol/aqueous dust surface area 
!       (excluding BCPO and OCPO), WTAREA and WERADIUS. (tmf, 3/6/09)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : ERADIUS, IXSAVE, IYSAVE 
      USE COMODE_MOD,    ONLY : IZSAVE,  JLOP,   TAREA
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIAG_MOD,      ONLY : AD21
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE COMODE_MOD,    ONLY : WTAREA, WERADIUS
      USE DAO_MOD,       ONLY : RH

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"   ! LPAR, CMN_SIZE
#     include "jv_cmn.h"   ! ODMDUST, QAA, RAA
#     include "CMN_DIAG"   ! ND21, LD21
#     include "comode.h"   ! NTTLOOP

      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: DUST(IIPAR,JJPAR,LLPAR,NDUST)

      ! Local variables
      INTEGER            :: I, J, JLOOP, L, N
      TYPE (XPLEX)             :: MSDENS(NDUST), XTAU
      ! Added to calculate aqueous dust surface area (WTAREA, WERADIUS)
      ! (tmf, 3/6/09)
      TYPE (XPLEX)             :: XRH
      TYPE (XPLEX)             :: CRITRH      ! Critical RH [%], above which 
                                        !  heteorogeneous chem takes place

      !=================================================================
      ! RDUST_ONLINE begins here!
      !=================================================================

      ! Dust density 
      MSDENS(1) = 2500.0d0
      MSDENS(2) = 2500.0d0
      MSDENS(3) = 2500.0d0
      MSDENS(4) = 2500.0d0
      MSDENS(5) = 2650.0d0
      MSDENS(6) = 2650.0d0
      MSDENS(7) = 2650.0d0

      ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
      CRITRH = 35.0d0   ! [%]

      !=================================================================     
      ! Convert concentration [kg/m3] to optical depth [unitless].
      !
      ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
      !           ( MSDENS * RAA * 1e-6 )
      ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
      !
      !  Units ==> DUST     [ kg/m3    ]
      !            MSDENS   [ kg/m3    ]
      !            RAA      [ um       ]
      !            BXHEIGHT [ m        ]
      !            QAA      [ unitless ]
      !            ODMDUST  [ unitless ]
      !
      ! NOTES: 
      ! (1) Do the calculation at QAA(4,:) (i.e. 999 nm).          
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NDUST
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ODMDUST(I,J,L,N) = 0.75d0        * BXHEIGHT(I,J,L) * 
     &                      DUST(I,J,L,N) * QAA(4,14+N)     / 
     &                     ( MSDENS(N) * RAA(4,14+N) * 1.0D-6 )
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !==============================================================
      ! Calculate Dust Surface Area
      !
      ! Units ==> DUST     [ kg dust/m^3 air    ]
      !           MSDENS   [ kg dust/m^3 dust   ]
      !           RAA      [ um                 ]
      !           TAREA    [ cm^2 dust/cm^3 air ]
      !           ERADIUS  [ cm                 ]
      !
      ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
      !       multiply by 3/radius to convert to surface area in cm2
      !  
      ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
      ! the NDUST FAST-J dust wavelength bins (read into DUST)
      !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N, XRH )
      DO N     = 1, NDUST
      DO JLOOP = 1, NTTLOOP

         ! Compute 3-D grid box indices
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         L = IZSAVE(JLOOP)

         ERADIUS(JLOOP,N) = RAA(4,14+N) * 1.0D-4

         TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                      DUST(I,J,L,N) / MSDENS(N)  

         ! Archive WTAREA and WERADIUS when RH > 35%  (tmf, 6/13/07) 
         ! Get RH  
         XRH     = RH( I, J, L )   ![%]
         WTAREA(JLOOP, N)   = 0.d0
         WERADIUS(JLOOP, N) = 0.d0

         IF ( XRH >= CRITRH ) THEN
            WTAREA(JLOOP, N)   = TAREA(JLOOP, N)
            WERADIUS(JLOOP, N) = ERADIUS(JLOOP, N)
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! ND21 Diagnostic: 
      !
      ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
      ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
      ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
      ! Tracer #4: Dust optical depths at 400 nm (from all size bins)
      ! Tracer #5: Dust surface areas (from all size bins)
      !==============================================================
      IF ( ND21 > 0 ) THEN

         DO N = 1, NDUST

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L )
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               !--------------------------------------
               ! ND21 tracer #4: Dust optical depths
               !--------------------------------------
               AD21(I,J,L,4) = AD21(I,J,L,4) + 
     &              ( ODMDUST(I,J,L,N) * QAA(2,14+N) / QAA(4,14+N) )

               !--------------------------------------
               ! ND21 tracer #5: Dust surface areas
               !--------------------------------------
               IF ( L <= LLTROP ) THEN
               
                  ! Convert 3-D indices to 1-D index
                  ! JLOP is only defined in the tropopause
                  JLOOP = JLOP(I,J,L)
             
                  ! Add to AD21
                  IF ( JLOOP > 0 ) THEN
                     AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                  ENDIF
               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDDO
      ENDIF 

      ! Return to calling program
      END SUBROUTINE RDUST_ONLINE

!------------------------------------------------------------------------------

      SUBROUTINE RDUST_OFFLINE( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine RDUST_OFFLINE reads global mineral dust concentrations as 
!  determined by P. Ginoux.  Calculates dust optical depth at each level for 
!  the FAST-J routine "set_prof.f". (rvm, bmy, 9/30/00, 8/4/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Number of the current month (1-12)
!  (2 ) THISYEAR  (INTEGER) : 4-digit year number (e.g. 1996, 2001)
!
!  NOTES:
!  (1 ) RDUST was patterned after rdaerosol.f (rvm, 9/30/00)
!  (2 ) Don't worry about rewinding the binary file...reading from
!        binary files is pretty fast.  And it's only done once a month.
!  (3 ) Now references punch file utility routines from F90 module
!        "bpch2_mod.f".  Also reference variable DATA_DIR from the
!         header file "CMN_SETUP". (bmy, 9/30/00) 
!  (4 ) Now selects proper GEOS-STRAT dust field for 1996 or 1997.
!        Also need to pass THISYEAR thru the arg list. (rvm, bmy, 11/21/00)
!  (5 ) CONC is now declared as TYPE (XPLEX) (rvm, bmy, 12/15/00)
!  (6 ) Removed obsolete code from 12/15/00 (bmy, 12/21/00)
!  (7 ) CONC(IGLOB,JGLOB,LGLOB,NDUST) is now CONC(IIPAR,JJPAR,LLPAR,NDUST).
!        Now use routine TRANSFER_3D from "transfer_mod.f" to cast from TYPE (XPLEX)
!        to TYPE (XPLEX) and also to convert from {IJL}GLOB to IIPAR,JJPAR,LLPAR 
!        space.  Use 3 arguments in call to GET_TAU0.  Updated comments.
!        (bmy, 9/26/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (9 ) Now reference ERADIUS, IXSAVE, IYSAVE, IZSAVE, TAREA from 
!        "comode_mod.f".  Compute ERADIUS and TAREA for the NDUST dust
!        size bins from FAST-J.  Renamed CONC to DUST to avoid conflicts.
!        Also reference NTTLOOP from "comode.h".  Also added parallel
!        DO-loops.  Also renamed MONTH and YEAR to THISMONTH and THISYEAR
!        to avoid conflicts w/ other variables. (bmy, 11/15/01)
!  (10) Bug fix: Make sure to use 1996 dust data for Dec 1995 for the
!        GEOS-STRAT met field dataset.  Set off CASE statement with an
!        #if defined( GEOS_STRAT ) block. (rvm, bmy, 1/2/02)
!  (11) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (12) Now report dust optical depths in ND21 diagnostic at 400 nm.  Now
!       report dust optical depths as one combined diagnostic field instead 
!        of 7 separate fields.  Now reference JLOP from "comode_mod.f".  
!        Now save aerosol surface areas as tracer #5 of the ND21 diagnostic.  
!        (rvm, bmy, 2/28/02)
!  (13) Remove declaration for TIME, since that is also defined in the
!        header file "comode.h" (bmy, 3/20/02)
!  (14) Now read mineral dust files directly from the DATA_DIR/dust_200203/
!        subdirectory (bmy, 4/2/02)
!  (15) Now reference BXHEIGHT from "dao_mod.f".  Also reference ERROR_STOP
!        from "error_mod.f". (bmy, 10/15/02)
!  (16) Now call READ_BPCH2 with QUIET=TRUE to suppress extra informational
!        output from being printed.  Added cosmetic changes. (bmy, 3/14/03)
!  (17) Since December 1997 dust data does not exist, use November 1997 dust
!        data as a proxy. (bnd, bmy, 6/30/03)
!  (18) Bundled into "dust_mod.f" and renamed to RDUST_OFFLINE. (bmy, 4/1/04)
!  (19) Now references DATA_DIR from "directory_mod.f".  Now parallelize over 
!        the L-dimension for ND21 diagnostic. (bmy, 7/20/04)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (22) Archive only hydrophilic aerosol/aqueous dust surface area 
!       (excluding BCPO and OCPO), WTAREA and WERADIUS. (tmf, 3/6/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE COMODE_MOD,    ONLY : ERADIUS, IXSAVE, IYSAVE
      USE COMODE_MOD,    ONLY : IZSAVE,  JLOP,   TAREA
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIAG_MOD,      ONLY : AD21
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE COMODE_MOD,    ONLY : WTAREA, WERADIUS
      USE DAO_MOD,       ONLY : RH

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"   ! LPAR, CMN_SIZE
#     include "jv_cmn.h"   ! ODMDUST, QAA, RAA
#     include "CMN_DIAG"   ! ND21, LD21
#     include "comode.h"   ! NTTLOOP

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH, THISYEAR

      ! Local variables
      INTEGER             :: I, J, JLOOP, L, N
      INTEGER, SAVE       :: MONTH_LAST = -999
      TYPE (XPLEX)              :: TEMP(IGLOB,JGLOB,LGLOB)
      TYPE (XPLEX)              :: DUST(IIPAR,JJPAR,LLPAR,NDUST)
      TYPE (XPLEX)              :: MSDENS(NDUST), XTAU
      CHARACTER (LEN=255) :: FILENAME

      ! Added to calculate aqueous dust surface area (WTAREA, WERADIUS)
      ! (tmf, 3/6/09)
      TYPE (XPLEX)             :: XRH
      TYPE (XPLEX)             :: CRITRH      ! Critical RH [%], above which 
                                        !  heteorogeneous chem takes place

      !=================================================================
      ! RDUST begins here!
      !
      ! Read aerosol data from the binary punch file during the first 
      ! chemistry timestep and, after that, at the start of each month.
      !=================================================================
      IF ( THISMONTH /= MONTH_LAST ) THEN   
         
         ! Save the current month
         MONTH_LAST = THISMONTH

         ! Get TAU0 value used to index the punch file
         ! Use the "generic" year 1985
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
         
         ! Select proper dust file name for GEOS-1, GEOS-3, or GEOS-4
         FILENAME = TRIM( DATA_DIR ) // 'dust_200203/dust.' //
     &              GET_NAME_EXT()   // '.'                 // 
     &              GET_RES_EXT()

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - RDUST: Reading ', a )

         ! Read aerosol concentrations [kg/m3] for each 
         ! dust type from the binary punch file
         DO N = 1, NDUST 
            CALL READ_BPCH2( FILENAME, 'MDUST-$', N,     XTAU,
     &                       IGLOB,     JGLOB,    LGLOB, TEMP, 
     &                       QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, DUST(:,:,:,N) )
         ENDDO

         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
         !           ( MSDENS * RAA * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         !  Units ==> DUST     [ kg/m3    ]
         !            MSDENS   [ kg/m3    ]
         !            RAA      [ um       ]
         !            BXHEIGHT [ m        ]
         !            QAA      [ unitless ]
         !            ODMDUST  [ unitless ]
         !
         ! NOTES: 
         ! (1) Do the calculation at QAA(4,:) (i.e. 999 nm).          
         !==============================================================
         MSDENS(1) = 2500.0d0
         MSDENS(2) = 2500.0d0
         MSDENS(3) = 2500.0d0
         MSDENS(4) = 2500.0d0
         MSDENS(5) = 2650.0d0
         MSDENS(6) = 2650.0d0
         MSDENS(7) = 2650.0d0

      ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
      CRITRH = 35.0d0   ! [%]

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, NDUST
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ODMDUST(I,J,L,N) = 0.75d0        * BXHEIGHT(I,J,L) * 
     &                         DUST(I,J,L,N) * QAA(4,14+N)     / 
     &                        ( MSDENS(N) * RAA(4,14+N) * 1.0D-6 )
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Echo information
         WRITE( 6, 110 ) 
 110     FORMAT( '     - RDUST: Finished computing optical depths' )

         !==============================================================
         ! Calculate Dust Surface Area
         !
         ! Units ==> DUST     [ kg dust/m^3 air    ]
         !           MSDENS   [ kg dust/m^3 dust   ]
         !           RAA      [ um                 ]
         !           TAREA    [ cm^2 dust/cm^3 air ]
         !           ERADIUS  [ cm                 ]
         !
         ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
         !       multiply by 3/radius to convert to surface area in cm2
         !  
         ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
         ! the NDUST FAST-J dust wavelength bins (read into DUST)
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N, XRH )
         DO N     = 1, NDUST
         DO JLOOP = 1, NTTLOOP

            ! Compute 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ERADIUS(JLOOP,N) = RAA(4,14+N) * 1.0D-4

            TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                         DUST(I,J,L,N) / MSDENS(N)  

         ! Archive WTAREA and WERADIUS when RH > 35% (tmf, 6/13/07) 
         ! Get RH  
         XRH     = RH( I, J, L )   ![%]
         WTAREA(JLOOP, N)   = 0.d0
         WERADIUS(JLOOP, N) = 0.d0

         IF ( XRH >= CRITRH ) THEN
            WTAREA(JLOOP, N)   = TAREA(JLOOP, N)
            WERADIUS(JLOOP, N) = ERADIUS(JLOOP, N)
         ENDIF

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! ND21 Diagnostic: 
         !
         ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
         ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
         ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
         ! Tracer #4: Dust optical depths at 400 nm (from all size bins)
         ! Tracer #5: Dust surface areas (from all size bins)
         !==============================================================
         IF ( ND21 > 0 ) THEN
            DO N = 1, NDUST

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L ) 
               DO L = 1, LD21
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  !--------------------------------------
                  ! ND21 tracer #4: Dust optical depths
                  !--------------------------------------
                  AD21(I,J,L,4) = AD21(I,J,L,4) + 
     &               ( ODMDUST(I,J,L,N) * QAA(2,14+N) / QAA(4,14+N) )

                  !--------------------------------------
                  ! ND21 tracer #5: Dust surface areas
                  !--------------------------------------
                  IF ( L <= LLTROP ) THEN

                     ! Convert 3-D indices to 1-D index
                     ! JLOP is only defined in the tropopause
                     JLOOP = JLOP(I,J,L)
             
                     ! Add to AD21
                     IF ( JLOOP > 0 ) THEN
                        AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                     ENDIF
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDDO
         ENDIF 
      ENDIF

      ! Return to calling program
      END SUBROUTINE RDUST_OFFLINE

!------------------------------------------------------------------------------

      FUNCTION GET_SCALE_GROUP( ) RESULT( CURRENT_GROUP )
!
!********************************************************************************
! Subroutine GET_SCALE_GROUP determines which predifined scaling index corresponds
! to the current time and location  (dkh, 12/02/04)
!
! Added to forward model (xxu, dkh, 01/13/12, adj32_011) 
!
! NOTES
! (1 ) CURRENT_GROUP is currently only a function of TAU
! (2 ) Get rid of I,J as argument. (dkh, 03/28/05)
!
!********************************************************************************

      ! Reference to f90 modules
      USE TIME_MOD,       ONLY : GET_TAU, GET_TAUe, GET_TAUb, GET_MONTH
      USE ADJ_ARRAYS_MOD, ONLY: MMSCL

#     include "CMN_SIZE" ! Size stuff

      ! Arguments
      INTEGER       :: I, J

      ! Local Variables
      TYPE (XPLEX)        :: TOTAL_HR, CURRENT_HR, GROUP_LENGTH
      TYPE (XPLEX)        :: TAU, TAUe, TAUb

      ! Function variable
      INTEGER       :: CURRENT_GROUP
      LOGICAL, SAVE :: MONTHLY = .FALSE.
      INTEGER, SAVE :: MONTH_SAVE
      INTEGER, SAVE :: GROUP_SAVE
      LOGICAL, SAVE :: FIRST = .TRUE.

      !============================================================
      ! GET_SCALE_GROUP begins here!
      !============================================================

      ! Currently there is no spatial grouping

      ! Determine temporal grouping
      IF ( MMSCL == 1 ) THEN
         CURRENT_GROUP = 1
         RETURN
      ENDIF

      IF ( MONTHLY ) THEN
         IF (FIRST) THEN
            MONTH_SAVE = GET_MONTH()
            CURRENT_GROUP = MMSCL
            GROUP_SAVE = MMSCL
            FIRST = .FALSE.
         ENDIF
         IF ( MONTH_SAVE /= GET_MONTH() ) THEN
            MONTH_SAVE = GET_MONTH()
            GROUP_SAVE = GROUP_SAVE - 1
            CURRENT_GROUP = GROUP_SAVE
         ELSE
            CURRENT_GROUP = GROUP_SAVE
         ENDIF

      ELSE
         ! Retrieve time parameters
         TAUe       = GET_TAUe()
         TAUb       = GET_TAUb()
         TAU        = GET_TAU()
         TOTAL_HR   = TAUe - TAUb
         CURRENT_HR = TAU  - TAUb

         ! The last time step always belongs to the last group
         IF ( TAU == TAUe ) THEN
            CURRENT_GROUP = MMSCL
            RETURN
         ELSE

            ! Determine the length of each group
            GROUP_LENGTH = ( TOTAL_HR / MMSCL )

            ! Index is the current time divided by the group length, plus one
            CURRENT_GROUP = ( CURRENT_HR / GROUP_LENGTH ) + 1

         ENDIF

      ENDIF

      ! Return to the calling routine
      END FUNCTION GET_SCALE_GROUP

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DUST
!
!******************************************************************************
!  Subroutine INIT_DUST allocates all module arrays (bmy, 3/30/04)
! 
!  NOTES:
!  (1 ) Now references LDEAD from "logical_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD, ONLY : LDEAD
      USE ERROR_MOD,   ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters
      
      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_DUST begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Drydep flags
      ALLOCATE( IDDEP( NDSTBIN ), STAT=AS )      
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IDDEP' )
      IDDEP = 0

      ! Dust radii
      ALLOCATE( DUSTREFF( NDSTBIN ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DUSTREFF' )
      DUSTREFF(1:NDSTBIN)%r = (/ 0.73d-6, 1.4d-6,  2.4d-6,  4.5d-6  /)
      DUSTREFF(1:NDSTBIN)%i = 0d0
      ! Dust density
      ALLOCATE( DUSTDEN( NDSTBIN ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DUSTREFF' )
      DUSTDEN(1:NDSTBIN)%r = (/ 2500.d0, 2650.d0, 2650.d0, 2650.d0 /)
      DUSTDEN(1:NDSTBIN)%i = 0d0
      ! These only have to be allocated for the Ginoux source function
      IF ( .not. LDEAD ) THEN
         ALLOCATE( FRAC_S( NDSTBIN ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRAC_S' )
         FRAC_S(1:NDSTBIN)%r   = (/ 0.095d0, 0.3d0, 0.3d0, 0.3d0 /)
         FRAC_S(1:NDSTBIN)%i   =0d0
      ENDIF

      ! Reset flag
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_DUST

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DUST
!
!******************************************************************************
!  Subroutine CLEANUP_DUST deallocates all module arrays (bmy, 3/30/04)
! 
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DUST begins here!
      !=================================================================
      IF ( ALLOCATED( IDDEP    ) ) DEALLOCATE( IDDEP    )      
      IF ( ALLOCATED( FRAC_S   ) ) DEALLOCATE( FRAC_S   )
      IF ( ALLOCATED( DUSTREFF ) ) DEALLOCATE( DUSTREFF )
      IF ( ALLOCATED( DUSTDEN  ) ) DEALLOCATE( DUSTDEN  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DUST

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DUST_MOD
