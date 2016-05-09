! $Id: sulfate_adj_mod.f,v 1.7 2012/03/01 22:00:26 daven Exp $
      MODULE SULFATE_ADJ_MOD
!
!******************************************************************************
!  Module SULFATE_MOD contains arrays and routines for performing either a
!  coupled chemistry/aerosol run or an offline sulfate aerosol simulation.
!  Original code taken from Mian Chin's GOCART model and modified accordingly.
!  (rjp, bdf, bmy, 6/22/00, 10/31/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) XNUMOL_OH  (TYPE (XPLEX) ) : Molecules OH  per kg OH          [molec/kg]
!  (2 ) XNUMOL_O3  (TYPE (XPLEX) ) : Molecules O3  per kg O3          [molec/kg]
!  (3 ) XNUMOL_NO3 (TYPE (XPLEX) ) : Molecules NO3 per kg NO3         [molec/kg]
!  (4 ) TCVV_S     (TYPE (XPLEX) ) : Ratio: Molwt air / Molwt S       [unitless]
!  (5 ) DMSo       (TYPE (XPLEX) ) : DMS oceanic emissions            [v/v/timestep]
!  (6 ) DRYH2O2    (INTEGER) : Pointer to H2O2  in DEPVEL array [unitless] 
!  (7 ) DRYSO2     (INTEGER) : Pointer to SO2   in DEPVEL array [unitless]
!  (8 ) DRYSO4     (INTEGER) : Pointer to SO4   in DEPVEL array [unitless]
!  (9 ) DRYSO4s    (INTEGER) : Pointer to SO4s  in DEPVEL array [unitless]
!  (10) DRYMSA     (INTEGER) : Pointer to MSA   in DEPVEL array [unitless]
!  (11) DRYNH3     (INTEGER) : Pointer to NH3   in DEPVEL array [unitless]
!  (12) DRYNH4     (INTEGER) : Pointer to NH4   in DEPVEL array [unitless]
!  (13) DRYNIT     (INTEGER) : Pointer to NIT   in DEPVEL array [unitless]
!  (14) DRYNITs    (INTEGER) : Pointer to NITs  in DEPVEL array [unitless]
!  (15) DRYSO4aq   (INTEGER) : Pointer to SO4aq in DEPVEL array [unitless]
!  (16) DRYAS      (INTEGER) : Pointer to AS    in DEPVEL array [unitless]  
!  (17) DRYAHS     (INTEGER) : Pointer to AHS   in DEPVEL array [unitless]
!  (18) DRYLET     (INTEGER) : Pointer to LET   in DEPVEL array [unitless]
!  (19) DRYNH4aq   (INTEGER) : Pointer to NH4aq in DEPVEL array [unitless]
!  (20) ENH3_an    (TYPE (XPLEX) ) : NH3 anthropogenic emissions      [kg NH3/box/s]
!  (21) ENH3_bb    (TYPE (XPLEX) ) : NH3 biomass emissions            [kg NH3/box/s]
!  (22) ENH3_bf    (TYPE (XPLEX) ) : NH3 biofuel emissions            [kg NH3/box/s]
!  (23) ENH3_na    (TYPE (XPLEX) ) : NH73 natural source emissions    [kg NH3/box/s]
!  (24) ESO2_ac    (TYPE (XPLEX) ) : SO2 aircraft emissions           [kg SO2/box/s]
!  (25) ESO2_an    (TYPE (XPLEX) ) : SO2 anthropogenic emissions      [kg SO2/box/s]
!  (26) ESO2_ev    (TYPE (XPLEX) ) : SO2 eruptive volcanic em.        [kg SO2/box/s]
!  (27) ESO2_nv    (TYPE (XPLEX) ) : SO2 non-eruptive volcanic em.    [kg SO2/box/s]
!  (28) ESO2_bb    (TYPE (XPLEX) ) : SO2 biomass burning emissions    [kg SO2/box/s]
!  (29) ESO2_bf    (TYPE (XPLEX) ) : SO2 biofuel burning emissions    [kg SO2/box/s]
!  (30) ESO2_sh    (TYPE (XPLEX) ) : SO2 ship emissions               [kg SO2/box/s]
!  (31) ESO4_an    (TYPE (XPLEX) ) : SO4 anthropogenic emissions      [kg SO4/box/s]
!  (32) JH2O2      (TYPE (XPLEX) ) : Monthly mean J(H2O2) values      [s-1]
!  (33) O3m        (TYPE (XPLEX) ) : Monthly mean O3 concentration    [v/v]
!  (34) PH2O2m     (TYPE (XPLEX) ) : Monthly mean P(H2O2)             [molec/cm3/s]
!  (35) PMSA_DMS   (TYPE (XPLEX) ) : P(MSA) from DMS                  [v/v/timestep]
!  (36) PSO2_DMS   (TYPE (XPLEX) ) : P(SO2) from DMS                  [v/v/timestep]
!  (37) PSO4_SO2   (TYPE (XPLEX) ) : P(SO4) from SO2                  [v/v/timestep]
!  (38) SSTEMP     (TYPE (XPLEX) ) : Sea surface temperatures         [K]
!  (39) VCLDF      (TYPE (XPLEX) ) : Volume cloud frac. for SO2 aq.   [unitless]
!  (40) NEV        (INTEGER) : Max # of eruptive volcanoes      [unitless]
!  (41) IEV        (INTEGER) : Longitudes of eruptive volcanoes [degrees]  
!  (42) JEV        (INTEGER) : Latitudes of eruptive volcanoes  [degrees ]
!  (43) IHGHT      (INTEGER) : Height of eruptive volcano plume [m]
!  (44) IELVe      (INTEGER) : Elevation of eruptive volcanoes  [m]
!  (45) Eev        (TYPE (XPLEX) ) : SO2 em. from eruptive volcanoes  [kg SO2/box/s]
!  (46) NNV        (INTEGER) : Max # of non-eruptive volcanoes  [unitless]
!  (47) NNVOL      (INTEGER) : Number of non-eruptive volcanoes [unitless]
!  (48) INV        (INTEGER) : Longitude of non-erup volcanoes  [degrees]
!  (49) JNV        (INTEGER) : Latitude of non-erup volcanoes   [degrees]
!  (50) IELVn      (INTEGER) : Elevation of non-erup volcanoes  [m]
!  (51) Env        (INTEGER) : SO2 em. from non-erup volcanoes  [kg SO2/box/s]
!  (52) TCOSZ      (TYPE (XPLEX) ) : Sum of cos(SZA) for offline run  [unitless] 
!  (53) TTDAY      (TYPE (XPLEX) ) : Total daylight length at (I,J)   [minutes]
!  (54) SMALLNUM   (TYPE (XPLEX) ) : Small number - prevent underflow [unitless]
!  (55) COSZM      (TYPE (XPLEX) ) : Array for MAX(cos(SZA)) at (I,J) [unitless]
!  
!  Module Routines:
!  ===========================================================================
!  (1 ) GET_VCLDF         : Computes volume cloud fraction for SO2 chemistry 
!  (2 ) GET_LWC           : Computes liquid water content as a function of T
!  (3 ) CHEMSULFATE       : Driver routine for sulfate/aerosol chemistry
!  (4 ) GRAV_SETTLING     : Routine to compute settling of SO4s and NITs
!  (5 ) CHEM_DMS          : Chemistry routine for DMS tracer
!  (6 ) CHEM_H2O2         : Chemistry routine for H2O2 tracer
!  (7 ) CHEM_SO2          : Chemistry routine for SO2 tracer
!  (8 ) SEASALT_CHEM      : Computes SO2->SO4 and HNO3->nitrate w/in seasalt
!  (9 ) AQCHEM_SO2        : Computes reaction rates for aqueous SO2 chemistry
!  (10) CHEM_SO4          : Chemistry routine for SO4 tracer
!  (11) PHASE_SO4         : Computes phase transition for crystalline tracers 
!  (12) PHASE_RADIATIVE   : Computes radiative forcing for crystalline tracers
!  (13) CHEM_MSA          : Chemistry routine for MSA tracer
!  (14) CHEM_NH3          : Chemistry routine for ammonia tracer
!  (15) CHEM_NH4          : Chemistry routine for ammonium tracer
!  (16) CHEM_NIT          : Chemistry routine for nitrates tracer
!  (17) EMISSSULFATE      : Driver routine for sulfate/aerosol emissions
!  (18) SRCDMS            : Emission routine for DMS tracer
!  (19) SRCSO2            : Emission routine for SO2 tracer
!  (20) SRCSO4            : Emission routine for SO4 tracer
!  (21) SRCNH3            : Emission routine for NH3 tracer
!  (22) GET_OH            : Returns OH for coupled or offline simulations
!  (23) SET_OH            : Resets modified OH in SMVGEAR's CSPEC array
!  (24) GET_NO3           : Returns NO3 for coupled or offline simulations
!  (25) SET_NO3           : Resets modified OH in SMVGEAR's CSPEC array
!  (26) GET_O3            : Returns O3 for coupled or offline simulations
!  (27) READ_NONERUP_VOLC : Reads SO2 emissions from non-eruptive volcanoes
!  (28) READ_ERUP_VOLC    : Reads SO2 emissions from eruptive volcanoes 
!  (29) READ_ANTHRO_SOx   : Reads anthropogenic SO2 and SO4 emissions
!  (30) READ_OCEAN_DMS    : Reads biogenic DMS emissions from oceans
!  (31) READ_SST          : Reads monthly mean sea-surface temperatures
!  (32) READ_BIOFUEL_SO2  : Reads SO2 emissions from biomass burning
!  (33) READ_AIRCRAFT_SO2 : Reads SO2 emissions from aircraft exhaust
!  (34) READ_SHIP_SO2     : Reads SO2 emissions from ship exhaust
!  (35) READ_ANTHRO_NH3   : Reads NH3 emissions from anthropogenic sources
!  (36) READ_NATURAL_NH3  : Reads NH3 emissions from natural sources
!  (37) READ_BIOMASS_NH3  : Reads NH3 biomass burning emissions
!  (38) READ_OXIDANT      : Reads monthly mean O3 and H2O2 for offline run
!  (39) OHNO3TIME         : Computes time arrays for scaling offline OH, NO3
!  (40) INIT_SULFATE      : Allocates & zeroes module arrays
!  (41) CLEANUP_SULFATE   : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by sulfate_mod.f
!  ============================================================================
!  (1 ) biomass_mod.f          : Module w/ routines for biomass burning
!  (2 ) bpch2_mod.f            : Module w/ routines for binary pch file I/O
!  (3 ) bravo_mod.f            : Module w/ routines to read BRAVO emissions
!  (4 ) comode_mod.f           : Module w/ SMVGEAR allocatable arrays
!  (5 ) dao_mod.f              : Module w/ DAO met field arrays
!  (6 ) diag_mod.f             : Module w/ GEOS-Chem diagnostic arrays
!  (7 ) directory_mod.f        : Module w/ GEOS-Chem data & met field dirs
!  (8 ) drydep_mod.f           : Module w/ GEOS-Chem dry deposition routines
!  (9 ) epa_nei_mod.f          : Module w/ routines to read EPA/NEI99 data
!  (10) error_mod.f            : Module w/ NaN, other error check routines
!  (11) file_mod.f             : Module w/ file unit numbers & error checks
!  (12) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (13) grid_mod.f             : Module w/ horizontal grid information
!  (14) global_no3_mod.f       : Module w/ routines to read 3-D NO3 field
!  (15) global_oh_mod.f        : Module w/ routines to read 3-D OH field
!  (16) isoropia_mod.f         : Module w/ ISORROPIA routines for aer thermo eq
!  (17) logical_mod.f          : Module w/ GEOS-Chem logical switches
!  (18) pbl_mix_mod.f          : Module w/ routines for PBL height & mixing
!  (19) pressure_mod.f         : Module w/ routines to compute P(I,J,L)
!  (20) seasalt_mod.f          : Module w/ routines for seasalt chemistry
!  (21) streets_anthro_mod.f   : Module w/ routines for David Streets' emiss
!  (22) time_mod.f             : Module w/ routines to compute time & date
!  (23) tracer_mod.f           : Module w/ GEOS-Chem tracer array STT etc.
!  (24) tracerid_mod.f         : Module w/ pointers to tracers & emissions
!  (25) transfer_mod.f         : Module w/ routines to cast & resize arrays
!  (26) tropopause_mod.f       : Module w/ routines to read ann mean tropopause
!  (27) uvalbedo_mod.f         : Module w/ UV albedo array and reader
!  (28) wetscav_mod.f          : Module w/ routines for wetdep & scavenging
!
!  References:
!  ============================================================================
!  (1 ) Andreae, M.O. & P. Merlet, "Emission of trace gases and aerosols from
!        biomass burning", Global Biogeochem. Cycles, 15, 955-966, 2001.
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!  (4 ) Wanninkhof, R., "Relation between wind speed and gas exchange over
!        the ocean", J. Geophys. Res, 97, 7373-7382, 1992.
!
!  NOTES:
!  (1 ) All module variables are declared PRIVATE (i.e., they can only
!        be seen from within this module (bmy, 6/2/00)
!  (2 ) The routines in "sulfate_mod.f" assume that we are doing chemistry
!        over the global region (e.g. IIPAR=IGLOB, JJPAR=JGLOB). (bmy, 6/8/00)
!  (3 ) Removed obsolete code from DRYDEP_SULFATE (bmy, 12/21/00)
!  (4 ) Removed obsolete commented-out code from module routines (bmy, 4/23/01)
!  (5 ) Now read data files from DATA_DIR/sulfate_sim_200106/ (bmy, 6/19/01)
!  (6 ) Updated comments (bmy, 9/4/01)
!  (7 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Now reference COSSZA from
!        "dao_mod.f". (bmy, 9/27/01)
!  (8 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (9 ) Minor fixes to facilitate compilation on ALPHA (bmy, 11/15/01)
!  (11) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (12) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (13) Now reference "file_mod.f" (bmy, 6/27/02)
!  (14) Now references GET_PEDGE from "pressure_mod.f", which computes P at
!        the bottom edge of grid box (I,J,L).  Also deleted obsolete,
!        commented-out code. (dsa, bdf, bmy, 8/21/02)
!  (15) Added updated code from Rokjin Park and Brendan Field, in order to
!        perform coupled chemistry-aerosol simulations.  Also added parallel
!        DO-loops in several subroutines.  Updated comments, cosmetic
!        changes.  Now reference "error_mod.f" and "wetscav_mod.f".  
!        Now only do chemistry below the tropopause. (rjp, bdf, bmy, 12/6/02)
!  (16) Added ENH3_na array to hold natural source NH3 emissions.  Also now
!        facilitate passing DMS, SO2, SO4, NH3 to SMVGEAR for fullchem
!        simulations.  Added subroutine READ_NATURAL_NH3. (rjp, bmy, 3/23/03)
!  (17) Now references "grid_mod.f" and "time_mod.f".  Also made other minor
!        cosmetic changes. (bmy, 3/27/03)
!  (18) Updated chemistry routines to apply drydep losses throughout the
!        entire PBL. (rjp, bmy, 8/1/03)
!  (19) Now accounts for GEOS-4 PBL being in meters (bmy, 1/15/04)
!  (20) Fix ND44 diag so that we get same results for sp or mp (bmy, 3/24/04)
!  (21) Added COSZM array.  Now use diurnal varying JH2O2 in CHEM_H2O2. 
!        (rjp, bmy, 3/39/04)
!  (22) Added more parallel DO-loops (bmy, 4/14/04)
!  (23) Now add SO2 from ships (bec, bmy, 5/20/04)
!  (24) Now references "directory_mod.f", "logical_mod.f" and "tracer_mod.f".
!        Now removed IJSURF. (bmy, 7/20/04)
!  (25) Can overwrite USA with EPA/NEI99 emissions (rjp, rch, bmy, 11/16/04)
!  (26) Modified for AS, AHS, LET, SO4aq, NH4aq (cas, bmy, 1/11/05)
!  (27) Now also references "pbl_mix_mod.f".  NOTE: Comment out phase 
!        transition  code for now since it is still under development and
!        will take a while to be rewritten. (bmy, 3/15/05)
!  (28) Modified for SO4s, NITs chemistry (bec, 4/13/05)
!  (29) Now reads updated files for SST and offline chemistry.  Now read data
!        for both GCAP and GEOS grids.  Now references "tropopause_mod.f".
!        (bmy, 8/22/05)
!  (30) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (31) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (32) Now read int'annual SST data on GEOS 1x1 grid (bmy, 11/17/05)
!  (33) Bug fix for offline aerosol sim in SEASALT_CHEM (bec, bmy, 3/29/06)
!  (34) Bug fix in INIT_DRYDEP (bmy, 5/23/06)
!  (35) Now references "bravo_mod.f" (rjp, kfb, bmy, 6/26/06)
!  (36) Now references "streets_anthro_mod.f" (yxw, bmy, 8/17/06)
!  (37) Now references "biomass_mod.f" (bmy, 9/27/06)
!  (38) Now prevent seg fault error in READ_BIOFUEL_SO2 (bmy, 11/3/06)
!  (39) Bug fix in SEASALT_CHEM (havala, bec, bmy, 12/8/06)
!  (40) Extra error check for low RH in GRAV_SETTLING (phs, 6/11/08)
!  (41) Now references "cac_anthro_mod.f".  And apply SO2 yearly scale factor
!        to SO2 from GEIA (amv, phs, 3/11/08)  
!  (41) Bug fixes in reading EDGAR data w/ the right tracer number, 
!        when we are doing offline or nonstd simulations (dkh, 10/31/08)
!  (42) Bug fix for AD13_SO2_sh in SRCSO2 (phs, 2/27/09)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "sulfate_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMSULFATE_ADJ
      PUBLIC :: CLEANUP_SULFATE_ADJ
      PUBLIC :: EMISSSULFATE_ADJ
     
      !=================================================================
      ! MODULE VARIABLES (see descriptions listed above)
      !=================================================================

      ! Time variable
      INTEGER              :: ELAPSED_SEC

      ! Logical Flags
      LOGICAL, PARAMETER   :: LENV = .TRUE.
      LOGICAL, PARAMETER   :: LEEV = .TRUE.
      
      ! Parameters
      TYPE (XPLEX),PARAMETER:: XNUMOL_OH  = xplex(6.022d23 / 17d-3,0d0)
      TYPE (XPLEX),PARAMETER:: XNUMOL_O3  = xplex(6.022d23 / 48d-3,0d0)
      TYPE (XPLEX),PARAMETER:: XNUMOL_NO3 = xplex(6.022d23 / 62d-3,0d0)
      TYPE (XPLEX),PARAMETER:: TCVV_S     = xplex(28.97d0  / 32d0,0d0)
      TYPE (XPLEX),PARAMETER:: SMALLNUM   = xplex(1d-20,0d0)

      ! Allocatable arrays
      TYPE (XPLEX),  ALLOCATABLE :: VCLDF(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: ADPSO4_SO2(:,:,:)

      ! Eruptive volcanoes
      INTEGER, PARAMETER   :: NEV=50
      INTEGER              :: NEVOL
      INTEGER, ALLOCATABLE :: IEV(:),   JEV(:)
      INTEGER, ALLOCATABLE :: IDAYs(:), IDAYe(:)
      INTEGER, ALLOCATABLE :: IHGHT(:), IELVe(:)
      TYPE (XPLEX),  ALLOCATABLE :: EEV(:)

      ! Non-eruptive volcanoes 
      INTEGER, PARAMETER   :: NNV=50
      INTEGER              :: NNVOL
      INTEGER, ALLOCATABLE :: INV(:), JNV(:), IELVn(:)
      TYPE (XPLEX),  ALLOCATABLE :: ENV(:)
      
      ! Pointers to drydep species w/in DEPSAV
      INTEGER              :: DRYSO2,  DRYSO4,   DRYMSA,  DRYNH3  
      INTEGER              :: DRYNH4,  DRYNIT,   DRYSO4s, DRYNITs
      INTEGER              :: DRYH2O2, DRYSO4aq, DRYAS,   DRYAHS
      INTEGER              :: DRYLET,  DRYNH4aq

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_VCLDF
!
!******************************************************************************
!  Subroutine GET_VCLDF computes the volume cloud fraction for SO2 chemistry.
!  (rjp, bdf, bmy, 9/23/02)
!
!  References:
!  ============================================================================
!  (1) Sundqvist et al. [1989]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : RH
      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER              :: I,    J,    L
      TYPE (XPLEX)               :: PRES, PSFC, RH2, R0, B0

      ! Parameters
      TYPE(XPLEX),PARAMETER::ZRT=xplex(0.60d0,0d0),ZRS=xplex(0.99d0,0d0)
		
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PSFC, PRES, RH2, R0, B0 )
      DO L = 1, LLTROP
      DO J = 1, JJPAR 
      DO I = 1, IIPAR
	
         ! Surface pressure
         PSFC = GET_PEDGE(I,J,1)

         ! Pressure at the center of the grid box
         PRES = GET_PCENTER(I,J,L)

         ! RH (from "dao_mod.f") is relative humidity [%]
         ! Convert to fraction and store in RH2
         RH2  = RH(I,J,L) * 1.0d-2

         ! Terms from Sundqvist ???
         R0   = ZRT + ( ZRS - ZRT ) * EXP( 1d0 - ( PSFC / PRES )**2.5 )
         B0   = ( RH2 - R0 ) / ( 1d0 - R0 )
	   
         ! Force B0 into the range 0-1
         IF ( RH2 < R0  ) B0 = 0d0
         IF ( B0  > 1d0 ) B0 = 1d0

         ! Volume cloud fraction
         VCLDF(I,J,L) = 1d0 - SQRT( 1d0 - B0 )

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GET_VCLDF

!------------------------------------------------------------------------------

      FUNCTION GET_LWC( T ) RESULT( LWC )
!
!******************************************************************************
!  Function GET_LWC returns the cloud liquid water content at a GEOS-CHEM
!  grid box as a function of temperature. (rjp, bmy, 10/31/02, 1/14/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T (TYPE (XPLEX)) : Temperature value at a GEOS-CHEM grid box [K]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: T

      ! Function value
      TYPE (XPLEX)             :: LWC

      !=================================================================
      ! GET_LWC begins here!
      !=================================================================

      ! Compute Liquid water content in [g/m3]
      IF ( T > 293d0 ) THEN
         LWC = 0.2d0

      ELSE IF ( T >= 280.d0 .AND. T <= 293.d0 ) THEN
         LWC = 0.32d0 - 0.0060d0 * ( T - 273.D0 ) 
 
      ELSE IF ( T >= 248.d0 .AND. T < 280.d0 ) THEN
         LWC = 0.23d0 + 0.0065d0 * ( T - 273.D0 )

      ELSE IF ( T < 248.d0 ) THEN
         LWC = 0.07d0

      ENDIF

      ! Convert from [g/m3] to [m3/m3]
      LWC = LWC * 1.D-6         

      ! Return to calling program
      END FUNCTION GET_LWC

!------------------------------------------------------------------------------
      
      SUBROUTINE CHEMSULFATE_ADJ

!******************************************************************************
!  Subroutine CHEMSULFATE_ADJ is the interface between the GEOS-CHEM main program
!  and the adjoint sulfate chemistry routines.  This version only supports
!  full chemistry.  (dkh, 10/12/05)
!
!  It is attempted to be as similar as possible to the forward
!  routine CHEMSULFATE  (rjp, bdf, bmy, 5/31/00, 3/16/06)
!
!  NOTES:
!  (1 ) Now we call ADJ_CHEM_xxx rather than computing the adjoint chemistry
!        routines directly. (dkh, 10/12/05) 
!  (2 ) Recalculate VCLDF. (dkh, 08/27/06)  
!  (3 ) Updated to GCv8 adjoint (dkh, 09/27/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD,     AIRDEN,  CLDF
      USE DAO_MOD,        ONLY : SUNCOS, CONVERT_UNITS
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH
      USE GLOBAL_NO3_MOD, ONLY : GET_GLOBAL_NO3
      USE LOGICAL_MOD,    ONLY : LCRYST,          LPRT
      USE TIME_MOD,       ONLY : GET_MONTH,       GET_TS_CHEM
      USE TIME_MOD,       ONLY : GET_ELAPSED_SEC, ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : STT,             TCVV 
      USE TRACER_MOD,     ONLY : N_TRACERS,       ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,   ONLY : IDTNITs,         IDTSO4s
 
#     include "CMN_SIZE"     ! Size parameters 

      ! Local variables
      LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE       :: LASTMONTH = -99
      INTEGER             :: I, J, L, N, MONTH
      TYPE (XPLEX)              :: DTCHEM

      ! External functions   
      TYPE (XPLEX),  EXTERNAL   :: BOXVL

      !=================================================================
      ! CHEMSULFATE_ADJ begins here!
      !=================================================================

      ! Get current month
      MONTH = GET_MONTH()

      ! Establish indices w/in DEPSAV array
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_SULFATE_ADJ

         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      ! If it's an offline simulation ...
      IF ( ITS_AN_AEROSOL_SIM() ) THEN
 
         CALL ERROR_STOP('offline aerosol not supported in adjoint',
     &                   'sulfate_adj_mod.f')
!
!         ! Then read monthly data files ...
!         IF ( ITS_A_NEW_MONTH() ) THEN 
!            CALL GET_GLOBAL_OH( MONTH )
!            CALL GET_GLOBAL_NO3( MONTH )
!         ENDIF
!
!         ! And compute time scaling arrays for offline OH, NO3
!         CALL OHNO3TIME
         
      ENDIF

      ! Store NTIME in a shadow variable
      ELAPSED_SEC = GET_ELAPSED_SEC()

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize adjoint module arrays
      !PSO2_DMS = 0d0
      !PMSA_DMS = 0d0
      !PSO4_SO2 = 0d0
      !PSO4_SS  = 0d0
      !PNITs    = 0d0
      ADPSO4_SO2  = 0d0 

      !================================================================= 
      ! Call individual chemistry routines for sulfate/aerosol tracers
      !=================================================================

      ! Redo forward model unit conversion
      
      ! Convert all tracers in STT from [kg] -> [v/v] 
      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CONV UNIT 1 ')

      ! fwd code:
      !! Convert STT from [v/v] -> [kg]
      !CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )
      ! adj code:
      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT_ADJ )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CONV UNIT 2 ')

      ! fwd code:
      !CALL CHEM_NIT
      ! adj code:
      CALL CHEM_NIT_ADJ
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CHEM_NIT_ADJ ' )
 
      ! fwd code:
      !CALL CHEM_NH4
      ! adj code:
      CALL CHEM_NH4_ADJ
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CHEM_NH4_ADJ ' )

      ! fwd code:
      !CALL CHEM_NH3
      ! adj code:
      CALL CHEM_NH3_ADJ
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CHEM_NH3_ADJ ' )

      ! fwd code:
      !CALL CHEM_MSA
      ! adj code:
      CALL CHEM_MSA_ADJ
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CHEM_MSA_ADJ ' )

      ! fwd code:
      !CALL CHEM_SO4
      ! adj code:
      CALL CHEM_SO4_ADJ
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CHEM_SO4_ADJ ' )

      ! SO2 
      ! added in v16 of GCv6 adj (dkh, 08/27/06)  
      CALL GET_VCLDF
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a GET_VCLDF ' )

      ! fwd code:
      !CALL CHEM_SO2
      ! adj code:
      CALL CHEM_SO2_ADJ

      ! Redo forward model unit conversion 
      ! Convert STT from [v/v] -> [kg]
      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CONV UNIT 3 ')


      ! We have already gone thru one chemistry iteration
      FIRSTCHEM = .FALSE. 

      ! For offline runs only ...
      IF ( ITS_AN_AEROSOL_SIM() ) THEN
         
         CALL ERROR_STOP('offline aerosol not supported in adjoint',
     &                   'sulfate_adj_mod.f')

         ! fwd code:
         !! H2O2 (offline only)
         !CALL CHEM_H2O2
         !IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_H2O2' )
         ! adj code:

         ! fwd code:
         !! DMS (offline only)
         !CALL CHEM_DMS
         !IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_DMS' ) 
         ! adj code:
      
      ENDIF


      ! fwd code:
      !! Convert all tracers in STT from [kg] -> [v/v] 
      !CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT )
      ! adj code:
      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT_ADJ )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE_ADJ: a CONV UNIT 4 ')

      ! fwd code:
      !! NITs [kg] gravitational settling 
      !CALL GRAV_SETTLING( STT(:,:,:,IDTNITs), 2 )
      !IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITS' )
      ! adj code:
      print*, ' WARNING: no adjoint of NITs gravitational settling'
         
      ! fwd code:
      !! SO4s [kg] gravitational settling 
      !CALL GRAV_SETTLING( STT(:,:,:,IDTSO4s), 1 )
      !IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4S' )
      ! adj code:
      print*, ' WARNING: no adjoint of SO4s gravitational settling'

      ! Return to calling program
      END SUBROUTINE CHEMSULFATE_ADJ

!!------------------------------------------------------------------------------
!!
!
!      SUBROUTINE GRAV_SETTLING( TC, N )
!!
!!******************************************************************************
!!  Subroutine GRAV_SETTLING performs gravitational settling of sulfate
!!  and nitrate in coarse sea salt (SO4S and NITS).
!!  (bec, rjp, bmy, 4/20/04, 7/20/04, 10/25/05)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) TC (TYPE (XPLEX) ) : Tracer [kg]
!!  (2 ) N  (INTEGER) : N=1 is SO4S; N=2 is NITS
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) TC (TYPE (XPLEX) ) : Contains modified tracer
!!
!!  NOTES:
!!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!!        (bmy, 7/20/04)
!!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!!  (3 ) Now limit relative humidity to [tiny(TYPE (XPLEX)),0.99] range for DLOG
!!         argument (phs, 5/1/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,       ONLY : T, BXHEIGHT, RH
!      USE DIAG_MOD,      ONLY : AD44
!      USE DRYDEP_MOD,    ONLY : DEPSAV
!      USE PRESSURE_MOD,  ONLY : GET_PCENTER
!      USE TRACER_MOD,    ONLY : SALA_REDGE_um,   SALC_REDGE_um,  XNUMOL
!      USE TRACERID_MOD,  ONLY : IDTSO4s,         IDTNITs
!      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
!      USE GRID_MOD,      ONLY : GET_AREA_CM2
!
!#     include "CMN_SIZE"      ! Size parameters
!#     include "CMN_GCTM"      ! g0
!#     include "CMN_DIAG"      ! ND44
!
!      ! Arguments
!      INTEGER, INTENT(IN)    :: N
!      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)
!
!      ! Local variables
!      INTEGER                :: I,      J,     L,        DTCHEM
!      TYPE (XPLEX)                 :: DELZ,   DELZ1, REFF
!      TYPE (XPLEX)                 :: P,      DP,    PDP,      TEMP        
!      TYPE (XPLEX)                 :: CONST,  SLIP,  VISC,     FAC1
!      TYPE (XPLEX)                 :: FAC2,   FLUX,  AREA_CM2, RHB
!      TYPE (XPLEX)                 :: RCM,    RWET,  RATIO_R,  RHO
!      TYPE (XPLEX)                 :: TOT1,   TOT2
!      TYPE (XPLEX)                 :: VTS(LLPAR)  
!      TYPE (XPLEX)                 :: TC0(LLPAR)
!      
!      ! Parameters
!      TYPE (XPLEX),  PARAMETER     :: C1 =  0.7674d0 
!      TYPE (XPLEX),  PARAMETER     :: C2 =  3.079d0 
!      TYPE (XPLEX),  PARAMETER     :: C3 =  2.573d-11
!      TYPE (XPLEX),  PARAMETER     :: C4 = -1.424d0
!      TYPE (XPLEX),  PARAMETER     :: DEN = 2200.0d0 ! [kg/m3] sea-salt density
!
!      ! Arrays
!      INTEGER              :: IDDEP(2)
!      INTEGER              :: IDTRC(2)	
!
!      !=================================================================
!      ! GRAV_SETTLING begins here!
!      !=================================================================
!
!      ! Return if tracers are undefined
!      IF ( IDTSO4s == 0 .and. IDTNITs == 0 ) RETURN
!
!      ! Return if it's the start of the run
!      IF ( GET_ELAPSED_SEC() == 0 ) RETURN
!
!      ! Chemistry timestep [s]
!      DTCHEM = GET_TS_CHEM() * 60d0
!
!      ! Store in IDDEP array
!      IDDEP(1) = DRYSO4s
!      IDDEP(2) = DRYNITs
!
!      ! Tracer array
!      IDTRC(1) = IDTSO4s
!      IDTRC(2) = IDTNITs
!
!      ! Coarse mode
!      REFF = 0.5d-6 * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
!            
!      ! Sea salt radius [cm]
!      RCM  = REFF * 100d0  
!
!      ! Exponential factors
!      FAC1 = C1 * ( RCM**C2 )
!      FAC2 = C3 * ( RCM**C4 )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) 
!!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC, TC0  )
!!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX             )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR       
!
!         ! Initialize 
!         DO L = 1, LLPAR
!            VTS(L) = 0d0
!         ENDDO
!
!         ! Loop over levels
!         DO L = 1, LLPAR
!
!            ! Pressure at center of the level [kPa]
!            P       = GET_PCENTER(I,J,L) * 0.1d0
!
!            ! Temperature [K]
!            TEMP    = T(I,J,L)
!
!            ! Cap RH at 0.99 
!            RHB     = MIN( 0.99d0, RH(I,J,L) * 1d-2 )
!
!            ! Safety check (phs, 5/1/08)
!            RHB     = MAX( TINY(RHB), RHB           )
!
!            ! Aerosol growth with relative humidity in radius [m] 
!            ! (Gerber, 1985)
!            RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHB))+RCM**3.d0)**0.33d0
!
!            ! Ratio dry over wet radii at the cubic power
!            RATIO_R = ( REFF / RWET )**3.d0
!
!            ! Density of the wet aerosol (kg/m3)
!            RHO     = RATIO_R * DEN + ( 1.d0 - RATIO_R ) * 1000.d0
!
!            ! Dp = particle diameter [um]
!            DP      = 2.d0 * RWET * 1.d6        
!
!            ! PdP = P * dP [hPa * um]
!            PDp     = P * Dp
!
!            ! Constant
!            CONST   = 2.d0 * RHO * RWET**2 * g0 / 9.d0
!
!            !===========================================================
!            ! NOTE: Slip correction factor calculations following 
!            ! Seinfeld, pp464 which is thought to be more accurate 
!            ! but more computation required. (rjp, 1/24/02)
!            !
!            ! # air molecule number density
!            ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
!            !
!            ! # gas mean free path
!            ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
!            !
!            ! # Slip correction
!            ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
!            !     &     / (2. * lamda))) / Dp
!            !
!            ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
!            ! which produces slip correction factore with small error
!            ! compared to the above with less computation.
!            !===========================================================  
!          
!            ! Slip correction factor (as function of P*dp)
!            Slip = 1.d0+(15.60d0 + 7.0d0 * EXP(-0.059d0 * PDp)) / PDp
!
!            ! Viscosity [Pa*s] of air as a function of temperature 
!            VISC = 1.458d-6 * (Temp)**(1.5d0) / ( Temp + 110.4d0 )
!
!            ! Settling velocity [m/s]
!            VTS(L) = CONST * Slip / VISC
!         ENDDO
!
!         ! Method is to solve bidiagonal matrix which is
!         ! implicit and first order accurate in z (rjp, 1/24/02)
!
!         ! Save initial tracer concentration in column
!         DO L = 1, LLPAR
!            TC0(L) = TC(I,J,L)
!         ENDDO
!
!         ! We know the boundary condition at the model top
!         L    = LLTROP
!         DELZ = BXHEIGHT(I,J,L)
!
!         TC(I,J,L) = TC(I,J,L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )
!
!         DO L = LLTROP-1, 1, -1
!            DELZ  = BXHEIGHT(I,J,L)
!            DELZ1 = BXHEIGHT(I,J,L+1)
!            TC(I,J,L) = 1.d0 / ( 1.d0 + DTCHEM * VTS(L) / DELZ )
!     &                * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1
!     &                *  TC(I,J,L+1) )
!         ENDDO
!         
!         !==============================================================
!         ! ND44 diagnostic: sea salt loss [molec/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 ) THEN
!
!            ! Initialize
!            TOT1 = 0d0
!            TOT2 = 0d0
!            
!            ! Compute column totals of TCO(:) and TC(I,J,:,N)
!            DO L = 1, LLPAR
!               TOT1 = TOT1 + TC0(L)
!               TOT2 = TOT2 + TC(I,J,L)
!            ENDDO
!
!            ! Surface area [cm2]
!            AREA_CM2 = GET_AREA_CM2( J )
!
!            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
!            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
!            FLUX     = FLUX * XNUMOL(IDTRC(N)) / AREA_CM2 
!   
!            ! Store in AD44 array
!            AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE GRAV_SETTLING
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE CHEM_DMS
!!
!!******************************************************************************
!!  Subroutine CHEM_DMS is the DMS chemistry subroutine from Mian Chin's    
!!  GOCART model, modified for use with the GEOS-CHEM model.
!!  (rjp, bdf, bmy, 5/31/00, 10/25/05)  
!!                                                                           
!!  Module Variables used:                                                     
!!  ============================================================================
!!  (1 ) PSO2_DMS (TYPE (XPLEX) ) : Array for P(SO2) from DMS [v/v]                
!!  (2 ) PMSA_DMS (TYPE (XPLEX) ) : Array for P(MSA) from DMS [v/v]                
!!                                                                            
!!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!!  ============================================================================
!!                                                                           
!!  R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel    
!!         k1 = { 1.7e-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }  
!!         a = 0.75, b = 0.25                                               
!!                                                                           
!!  R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel 
!!         k2 = 1.2e-11*exp(-260/T)                                         
!!                                                                           
!!         DMS_OH = DMS0 * exp(-(r1+r2)* NDT1)                                  
!!         where DMS0 is the DMS concentration at the beginning,            
!!         r1 = k1*[OH], r2 = k2*[OH].                                      
!!                                                                           
!!  R3:    DMS + NO3 ->   SO2 + ...                                         
!!         k3 = 1.9e-13*exp(500/T)                                          
!!                                                                           
!!         DMS = DMS_OH * exp(-r3*NDT1)                                         
!!         where r3 = k3*[NO3].                                             
!!                                                                           
!!  R4:    DMS + X   ->   SO2 + ...                                         
!!         assume to be at the rate of DMS+OH and DMS+NO3 combined.         
!!                                                                           
!!  The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved    
!!  for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They    
!!  are in unit of [v/v/timestep]. 
!!
!!  NOTES: 
!!  (1 ) Now reference AD, AIRDEN, and SUNCOS from "dao_mod.f".  Added 
!!        parallel DO-loops.  Also now extract OH and NO3 from SMVGEAR
!!        for coupled chemistry-aerosol runs. (rjp, bdf, bmy, 9/16/02)
!!  (2 ) Bug fix: remove duplicate definition of RK3 (bmy, 3/23/03)
!!  (3 ) Now use function GET_TS_CHEM from "time_mod.f".  (bmy, 3/27/03)
!!  (4 ) Now reference STT and ITS_A_FULLCHEM_SIM from "tracer_mod.f"
!!        Now replace IJSURF w/ an analytic function. (bmy, 7/20/04)
!!  (5 ) Shift rows 8,9 in AD05 to 9,10 in to make room for P(SO4) from O3 
!!        oxidation in sea-salt aerosols (bec, bmy, 4/13/05)
!!  (6 ) Now remove reference to CMN, it's obsolete.  Now reference 
!!        ITS_IN_THE_STRAT from "tropopause_mod.f". (bmy, 8/22/05)
!!  (7 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      USE DAO_MOD,        ONLY : AD, AIRDEN, SUNCOS, T
!      USE DIAG_MOD,       ONLY : AD05
!      USE DRYDEP_MOD,     ONLY : DEPSAV
!      USE TIME_MOD,       ONLY : GET_TS_CHEM
!      USE TRACER_MOD,     ONLY : STT, ITS_A_FULLCHEM_SIM, XNUMOL
!      USE TRACERID_MOD,   ONLY : IDTDMS
!      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
!
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_GCTM"     ! AIRMW
!#     include "CMN_DIAG"     ! ND05, LD05
!
!      ! Local variables
!      LOGICAL                :: IS_FULLCHEM
!      INTEGER                :: I,   J,    L,      IJLOOP
!      TYPE (XPLEX)                 :: TK,  O2,   RK1,    RK2,    RK3,   F  
!      TYPE (XPLEX)                 :: DMS, DMS0, DMS_OH, DTCHEM, XOH,   XN3 
!      TYPE (XPLEX)                 :: XX,  OH,   OH0,    XNO3,   XNO30, LOH
!      TYPE (XPLEX)                 :: LNO3
!
!      ! Parameters
!      TYPE (XPLEX), PARAMETER      :: FX = 1.0d0
!      TYPE (XPLEX), PARAMETER      :: A  = 0.75d0
!      TYPE (XPLEX), PARAMETER      :: B  = 0.25d0
!
!      ! From D4: only 0.8 efficiency, also some goes to DMSO and lost.  
!      ! So we assume 0.75 efficiency for DMS addtion channel to form     
!      ! products.                                                        
!      TYPE (XPLEX), PARAMETER      :: EFF = 1d0
!      
!      ! External functions   
!      TYPE (XPLEX),  EXTERNAL      :: BOXVL
!      
!      !=================================================================
!      ! CHEM_DMS begins here!
!      !=================================================================
!      IF ( IDTDMS == 0 ) RETURN
!
!      ! Flag for fullchem simulation
!      IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
!
!      ! DTCHEM is the chemistry timestep in seconds
!      DTCHEM = GET_TS_CHEM() * 60d0
!
!      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
!      f  = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
!      
!      !=================================================================
!      ! Do the chemistry over all tropospheric grid boxes!
!      !=================================================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, IJLOOP, J, L, TK, O2, DMS0, OH, XNO3, RK1, RK2 )
!!$OMP+PRIVATE( RK3, DMS_OH, DMS, OH0, XNO30, XOH, XN3, XX, LOH, LNO3  )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO L = 1, LLTROP  
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Skip stratospheric boxes
!         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE
!
!         ! IJLOOP is the 1-D grid box index for SUNCOS
!         IJLOOP = ( (J-1) * IIPAR ) + I
!
!         ! Temperature [K]
!         TK     = T(I,J,L)
!
!         ! Get O2 [molec/cm3], DMS [v/v], OH [molec/cm3], NO3 [molec/cm3]
!         O2     = AIRDEN(L,I,J) * f * 0.21d0
!         DMS0   = STT(I,J,L,IDTDMS)
!         OH     = GET_OH( I, J, L )
!         XNO3   = GET_NO3( I, J, L )
!
!         !==============================================================
!         ! (1) DMS + OH:  RK1 - addition channel  
!         !                RK2 - abstraction channel   
!         !==============================================================
!         RK1 = 0.d0
!         RK2 = 0.d0
!         RK3 = 0.d0
!
!         IF ( OH > 0.d0 ) THEN
!            RK1 = ( 1.7d-42 * EXP( 7810.d0 / TK ) * O2 ) /
!     &            ( 1.d0 + 5.5d-31 * EXP( 7460.d0 / TK ) * O2 ) * OH
!
!            RK2 = 1.2d-11 * EXP( -260.d0 / TK ) * OH 
!         ENDIF
!            
!         !==============================================================
!         ! (2) DMS + NO3 (only happens at night):  
!         !==============================================================
!         IF ( SUNCOS(IJLOOP) <= 0d0 ) THEN
!            RK3 = 1.9d-13 * EXP( 500.d0 / TK ) * XNO3
!         ENDIF
!
!         !==============================================================
!         ! Update DMS concentrations after reaction with OH and NO3, 
!         ! and also account for DMS + X assuming at a rate as 
!         ! (DMS+OH)*Fx in the day and (DMS+NO3)*Fx at night:   
!         ! 
!         ! DMS_OH :  DMS concentration after reaction with OH  
!         ! DMS    :  DMS concentration after reaction with NO3       
!         !           (min(DMS) = 1.0E-32)       
!         !
!         ! NOTE: If we are doing a coupled fullchem/aerosol run, then
!         ! also modify OH and NO3 concentrations after rxn w/ DMS.
!         !==============================================================
!         DMS_OH = DMS0   * EXP( -( RK1 + RK2 ) * Fx * DTCHEM )
!         DMS    = DMS_OH * EXP( -( RK3       ) * Fx * DTCHEM ) 
!         IF ( DMS < SMALLNUM ) DMS = 0d0
!
!         ! Archive initial OH and NO3 for diagnostics
!         OH0    = OH
!         XNO30  = XNO3
!
!         IF ( IS_FULLCHEM ) THEN
!         
!            ! Update OH after rxn w/ DMS (coupled runs only)
!            OH    = OH0 - ( ( DMS0 - DMS_OH ) * AIRDEN(L,I,J) * f )
!            IF ( OH < SMALLNUM ) OH = 0d0
!
!            ! Update NO3 after rxn w/ DMS (coupled runs only)
!            XNO3  = XNO30 - ( ( DMS_OH - DMS ) * AIRDEN(L,I,J) * f )
!            IF ( XNO3 < SMALLNUM ) XNO3 = 0d0
!
!         ENDIF 
!
!         ! Save DMS back to the tracer array
!         STT(I,J,L,IDTDMS) = DMS
!
!         !==============================================================
!         ! Save SO2 and MSA production from DMS oxidation 
!         ! in [mixing ratio/timestep]:    
!         !
!         ! SO2 is formed in DMS+OH addition (0.85) and abstraction 
!         ! (1.0) channels as well as DMS + NO3 reaction.  We also 
!         ! assume that SO2 yield from DMS + X is 1.0.  
!         !
!         ! MSA is formed in DMS + OH addition (0.15) channel. 
!         !==============================================================
!         IF ( ( RK1 + RK2 ) == 0.d0 ) THEN
!            PMSA_DMS(I,J,L) = 0.d0
!         ELSE
!            PMSA_DMS(I,J,L) = ( DMS0 - DMS_OH ) * 
!     &                          B*RK1 / ( ( RK1 + RK2 ) * Fx ) * EFF
!         ENDIF
!
!         PSO2_DMS(I,J,L) =  DMS0 - DMS - PMSA_DMS(I,J,L)
!
!         !==============================================================
!         ! ND05 diagnostic: production and loss  
!         !
!         ! For the offline run, we are reading in monthly mean OH, NO3 
!         ! from disk.  We don't modify these, so LOH = 0 and LNO3 = 0.
!         !==============================================================
!         IF ( ND05 > 0 .and. L <= LD05 ) THEN
!
!            ! P(SO2) from DMS+OH, DMS+NO3, and DMS+X
!            XOH  = ( DMS0   - DMS_OH ) / Fx * AD(I,J,L) / TCVV_S
!            XN3  = ( DMS_OH - DMS    ) / Fx * AD(I,J,L) / TCVV_S
!            XX   = ( ( DMS0 - DMS ) * AD(I,J,L) / TCVV_S ) - XOH - XN3
!        
!            ! Convert L(OH) and L(NO3) from [molec/cm3] to [kg/timestep]
!            LOH  = ( OH0   - OH   ) * BOXVL(I,J,L) / XNUMOL_OH
!            LNO3 = ( XNO30 - XNO3 ) * BOXVL(I,J,L) / XNUMOL_NO3 
!
!            ! Store P(SO2) from DMS + OH [kg S/timestep]
!            AD05(I,J,L,1) = AD05(I,J,L,1) + XOH
!
!            ! Store P(SO2) from DMS + NO3 [kg S/timestep]
!            AD05(I,J,L,2) = AD05(I,J,L,2) + XN3
!
!            ! Store total P(SO2) from DMS [kg S/timestep]
!            AD05(I,J,L,3) = AD05(I,J,L,3) + 
!     &                      ( PSO2_DMS(I,J,L) * AD(I,J,L) / TCVV_S )
!
!            ! Store P(MSA) from DMS [kg S/timestep]
!            AD05(I,J,L,4) = AD05(I,J,L,4) + 
!     &                      ( PMSA_DMS(I,J,L) * AD(I,J,L) / TCVV_S )
!
!            ! Store L(OH) by DMS [kg OH/timestep]
!            AD05(I,J,L,9) = AD05(I,J,L,9) + LOH
!            
!            ! Store L(NO3) by DMS [kg NO3/timestep]
!            AD05(I,J,L,10) = AD05(I,J,L,10) + LNO3
!
!         ENDIF
!
!         !==============================================================
!         ! For a coupled fullchem/aerosol run, save OH [molec/cm3] 
!         ! and NO3 [molec/cm3] back into the CSPEC array of SMVGEAR
!         !==============================================================
!         IF ( IS_FULLCHEM ) THEN
!            CALL SET_OH( I, J, L, OH )
!            CALL SET_NO3( I, J, L, XNO3 )
!         ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE CHEM_DMS
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE CHEM_H2O2
!!
!!******************************************************************************
!!  Subroutine CHEM_H2O2 is the H2O2 chemistry subroutine for offline sulfate
!!  simulations.  For coupled runs, H2O2 chemistry is already computed by
!!  the SMVGEAR module. (rjp, bmy, 11/26/02, 10/25/05)
!!                                                                           
!!  NOTES:
!!  (1 ) Bug fix: need to multiply DXYP by 1d4 for cm2 (bmy, 3/23/03)
!!  (2 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!!        Now use functions GET_MONTH and GET_TS_CHEM from "time_mod.f".
!!        (bmy, 3/27/03)
!!  (3 ) Now references PBLFRAC from "drydep_mod.f".  Now apply dry deposition 
!!        throughout the entire PBL.  Added FREQ variable. (bmy, 8/1/03)
!!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!!        sum into AD44 array.  This preents numerical differences when using
!!        multiple processors. (bmy, 3/24/04)    
!!  (5 ) Now use diurnally-varying JO1D.  Now use new unit conversion for
!!        the ND44 diagnostic. (rjp, bmy, 3/30/04)
!!  (6 ) Now use parallel DO-loop to zero ND44_TMP.  Now uses ITS_A_NEW_MONTH
!!        from time_mod.f. (bmy, 4/14/04)
!!  (7 ) Now reference STT & TCVV from "tracer_mod.f".  Also replace IJSURF
!!        with an analytic function.  Now references DATA_DIR from 
!!        "directory_mod.f". (bmy, 7/20/04)
!!  (8 ) Now suppress output from READ_BPCH with QUIET keyword (bmy, 1/25/05)
!!  (9 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP
!!        from "pbl_mix_mod.f" (bmy, 2/22/05)
!!  (10) Now read offline files from "sulfate_sim_200508/offline".  Now remove 
!!        reference to CMN, it's obsolete.  Now reference ITS_IN_THE_STRAT from 
!!        "tropopause_mod.f". (bmy, 8/22/05)
!!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (12) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,      ONLY : GET_NAME_EXT, GET_RES_EXT
!      USE BPCH2_MOD,      ONLY : GET_TAU0,     READ_BPCH2
!      USE DAO_MOD,        ONLY : AD, AIRDEN, OPTD, SUNCOS, T
!      USE DIAG_MOD,       ONLY : AD44 
!      USE DIRECTORY_MOD,  ONLY : DATA_DIR
!      USE DRYDEP_MOD,     ONLY : DEPSAV
!      USE GRID_MOD,       ONLY : GET_AREA_CM2
!      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP
!      USE TIME_MOD,       ONLY : GET_MONTH, GET_TS_CHEM, ITS_A_NEW_MONTH
!      USE TRACER_MOD,     ONLY : STT,       TCVV,        XNUMOL
!      USE TRACERID_MOD,   ONLY : IDTH2O2
!      USE TRANSFER_MOD,   ONLY : TRANSFER_3D_TROP
!      USE UVALBEDO_MOD,   ONLY : UVALBEDO
!      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
! 
!#     include "cmn_fj.h"       ! IPAR, JPAR, LPAR, CMN_SIZE
!#     include "CMN_GCTM"       ! AIRMW
!#     include "CMN_DIAG"       ! ND44
!      
!      ! Local variables
!      LOGICAL                 :: FIRST     = .TRUE.
!      INTEGER, SAVE           :: LASTMONTH = -99
!      INTEGER                 :: I, J, L, JLOOP
!      TYPE (XPLEX)                  :: ARRAY(IGLOB,JGLOB,LLTROP)
!      TYPE (XPLEX)                  :: ND44_TMP(IIPAR,JJPAR,LLTROP)
!      TYPE (XPLEX)                  :: DT,    Koh,  DH2O2, M,    F ,   XTAU   
!      TYPE (XPLEX)                  :: H2O20, H2O2, ALPHA, FLUX, FREQ, PHOTJ
!      TYPE (XPLEX),  PARAMETER      :: A = 2.9d-12
!      CHARACTER(LEN=255)      :: FILENAME
!
!      !=================================================================
!      ! CHEM_H2O2 begins here!
!      !=================================================================
!      IF ( IDTH2O2 == 0 .or. DRYH2O2 == 0 ) RETURN 
!
!      ! Chemistry timestep [s]
!      DT = GET_TS_CHEM() * 60d0
!
!      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
!      F  = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
!      
!      ! Zero ND44_TMP array
!      IF ( ND44 > 0 ) THEN
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLTROP
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            ND44_TMP(I,J,L) = 0d0
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF
!
!      !=================================================================
!      ! For offline run: read J(H2O2) from disk below
!      !=================================================================
!      IF ( ITS_A_NEW_MONTH() ) THEN
!
!         ! File name to read data 
!         FILENAME = TRIM( DATA_DIR )                       // 
!     &              'sulfate_sim_200508/offline/JH2O2.'    // 
!     &              GET_NAME_EXT() // '.' // GET_RES_EXT()
!           
!         ! Print filename
!         WRITE( 6, 100 ) TRIM( FILENAME )
! 100     FORMAT( '     - CHEM_H2O2: Reading ', a )
!
!         ! Get TAU0 value for this month in "generic" year 1985
!         XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )
!	
!         ! Read J(H2O2) [s-1]  from disk (only up to tropopause)
!         ! limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
!         ! tropopause. This is backward compatibility with 
!         ! offline data set.
!         CALL READ_BPCH2( FILENAME, 'JV-MAP-$', 3,      
!     &     XTAU,        IGLOB,                    JGLOB,      
!     &     LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )
!!     &                    XTAU,      IGLOB,     JGLOB,      
!!     &                    LLTROP,    ARRAY,     QUIET=.TRUE. )
!
!
!
!         ! Cast to TYPE (XPLEX) and resize if necessary
!         CALL TRANSFER_3D_TROP( ARRAY, JH2O2 )
!            
!         ! Reset LASTMONTH
!         !LASTMONTH = GET_MONTH()
!      ENDIF
!
!      !=================================================================
!      ! Loop over tropopsheric grid boxes and do chemistry
!      !=================================================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, M, H2O20, KOH, FREQ, ALPHA, DH2O2, H2O2, FLUX )
!!$OMP+PRIVATE( JLOOP, PHOTJ )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO L  = 1, LLTROP
!      DO J  = 1, JJPAR
!      DO I  = 1, IIPAR
!
!         ! Initialize for safety's sake 
!         FLUX = 0d0
!         FREQ = 0d0
!
!         ! Skip stratospheric boxes
!         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE
!
!         ! Density of air [molec/cm3]
!         M     = AIRDEN(L,I,J) * f  
!
!         ! Initial H2O2 [v/v]
!         H2O20 = STT(I,J,L,IDTH2O2)
!
!         ! Loss frequenty due to OH oxidation [s-1]
!         KOH   = A * EXP( -160.d0 / T(I,J,L) ) * GET_OH(I,J,L)  
!
!         ! H2O2 drydep frequency [1/s].  Account for the fraction
!         ! of grid box (I,J,L) that is located beneath the PBL top.
!         FREQ  = DEPSAV(I,J,DRYH2O2) * GET_FRAC_UNDER_PBLTOP( I, J, L ) 
!
!         ! 1-D grid box index for SUNCOS
!         JLOOP = ( (J-1) * IIPAR ) + I
!
!         ! Impose a diurnal variation of jH2O2 by multiplying COS of 
!         ! solar zenith angle normalized by maximum solar zenith angle 
!         ! because the archived JH2O2 is for local noon time
!         IF ( COSZM(I,J) > 0.d0 ) THEN
!            PHOTJ = JH2O2(I,J,L) * SUNCOS(JLOOP) / COSZM(I,J)
!            PHOTJ = MAX( PHOTJ, 0d0 )
!         ELSE
!            PHOTJ = 0d0
!         ENDIF
!
!         ! Compute loss fraction from OH, photolysis, drydep [unitless].  
!         ALPHA = 1.D0 + ( KOH + PHOTJ + FREQ ) * DT 
!
!         ! Delta H2O2 [v/v]
!         DH2O2 = PH2O2m(I,J,L) * DT / ( ALPHA * M )
!         
!         ! Final H2O2 [v/v]
!         H2O2  = ( H2O20 / ALPHA + DH2O2 )
!         IF ( H2O2 < SMALLNUM ) H2O2 = 0d0
!
!         ! Store final H2O2 in STT
!         STT(I,J,L,IDTH2O2) = H2O2
!
!         !==============================================================
!         ! ND44 diagnostics: H2O2 drydep loss [molec/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN
!
!            ! Convert H2O2 from [v/v] to H2O2 [molec/cm2/s]
!            FLUX = H2O20 * FREQ * DT / ( 1.D0 + FREQ * DT )
!            FLUX = FLUX * AD(I,J,L) / TCVV(IDTH2O2)
!            FLUX = FLUX * XNUMOL(IDTH2O2) / ( GET_AREA_CM2( J ) * DT )
!
!            ! Save dryd flx in ND44_TMP as a placeholder
!            ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
!         ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!      
!      !===============================================================
!      ! ND44: Sum drydep fluxes by level into the AD44 array in
!      ! order to ensure that  we get the same results w/ sp or mp 
!      !===============================================================
!      IF ( ND44 > 0 ) THEN 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!         DO L = 1, LLTROP
!            AD44(I,J,DRYH2O2,1) = AD44(I,J,DRYH2O2,1) + ND44_TMP(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE CHEM_H2O2
!
!!------------------------------------------------------------------------------

      SUBROUTINE CHEM_SO2_ADJ
!
!******************************************************************************
!  Subroutine CHEM_SO2_ADJ is the adjoint of the SO2 chemistry subroutine 
!  (dkh, 10/11/05) 
!                                                                          
!  Based on forward model subroutine CHEM_SO2 (rjp, bmy, 11/26/02, 10/25/05) 
! 
!  NOTES:                                                                   
!  (1 ) Now reference WETSCAV_MOD instead of WETSCAV_ADJ_MOD. (dkh, 10/24/05)  
!  (2 ) BUG FIX Add ADL3S to list of PRIVATE variables. (dkh, 10/08/06)  
!  (3 ) BUG FIX Initialize ADL3S to zero. (dkh, 01/25/07) 
!  (4 ) Updated to GCv8 adj (dkh, 09/27/09) 
!  (5 ) Several bugs fixed (Jamin, Daven) 
!  (6 ) Watch for L3 = 1d0 (dkh, 01/06/12, adj32_005) 
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE CHECKPT_MOD,     ONLY : SO2_CHK
      USE CHECKPT_MOD,     ONLY : H2O2_CHK
      USE DAO_MOD,         ONLY : AD,      AIRDEN, T
      USE DIAG_MOD,        ONLY : AD05,    AD44
      USE DRYDEP_MOD,      ONLY : DEPSAV
      USE DIRECTORY_MOD,   ONLY : DATA_DIR
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE GLOBAL_HNO3_MOD, ONLY : GET_GLOBAL_HNO3
      USE GRID_MOD,        ONLY : GET_AREA_CM2
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD 
      USE PBL_MIX_MOD,     ONLY : GET_FRAC_UNDER_PBLTOP
      USE PRESSURE_MOD,    ONLY : GET_PCENTER
      USE TIME_MOD,        ONLY : GET_TS_CHEM, GET_MONTH
      USE TIME_MOD,        ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,      ONLY : STT,     TCVV,  ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,      ONLY : XNUMOL
      USE TRACERID_MOD,    ONLY : IDTH2O2, IDTSO2
      USE SEASALT_MOD,     ONLY : GET_ALK
      USE WETSCAV_MOD,     ONLY : H2O2s,   SO2s
      USE WETSCAV_ADJ_MOD, ONLY : H2O2s_ADJ,   SO2s_ADJ
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_GCTM"    ! AIRMW

      ! Local variables
      LOGICAL               :: IS_OFFLINE
      INTEGER               :: I,      J,       L,      I1,   I2
      INTEGER               :: II,     NSTEP
      TYPE (XPLEX)                :: K0,     Ki,      KK,     M,    L1
      TYPE (XPLEX)                :: L2,     L3,      Ld,     F,    Fc
      TYPE (XPLEX)                :: RK,     RKT,     DTCHEM, DT_T, TK
      TYPE (XPLEX)                :: F1,     RK1,     RK2,    RK3,  SO20
      TYPE (XPLEX)                :: SO2_cd, H2O20,   O3,     L2S,  L3S
      TYPE (XPLEX)                :: LWC,    KaqH2O2, KaqO3,  PATM, FLUX
      TYPE (XPLEX)                :: ALK,    ALK1,    ALK2,    SO2_ss
      TYPE (XPLEX)                :: Kt1,    Kt2,     AREASS1, AREASS2
      TYPE (XPLEX)                :: PSO4E,  PSO4F,   Kt1N,    Kt2N
      TYPE (XPLEX)                :: AREA_CM2

      ! Parameters
      TYPE (XPLEX),PARAMETER::HPLUS =xplex(3.16227766016837953d-5,0d0)  !pH = 4.5
      TYPE (XPLEX),PARAMETER::MINDAT=xplex(1.d-20,0d0)

C==============================================
C define local TAMC generated variables
C==============================================
      TYPE (XPLEX) adh2o20
      TYPE (XPLEX) adl1
      TYPE (XPLEX) adl2
      TYPE (XPLEX) adl2s
      TYPE (XPLEX) adl3
      TYPE (XPLEX) adl3s
      TYPE (XPLEX) ado3
      TYPE (XPLEX) adso20
      TYPE (XPLEX) adso2_cd


      !=================================================================
      ! CHEM_SO2_ADJ begins here!
      !=================================================================
      IF ( IDTH2O2 == 0 .or. IDTSO2 == 0 .or. DRYSO2 == 0 ) RETURN

      ! Is it an offline simulation?
      IS_OFFLINE = ITS_AN_AEROSOL_SIM()

      ! Read HNO3 for offline simulation
      IF ( IS_OFFLINE ) THEN
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL GET_GLOBAL_HNO3( GET_MONTH() )
         ENDIF
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
      F      = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
      Ki     = 1.5d-12

      
      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO20, H2O20, O3, PATM, TK, K0, M, KK, F1, RK1  )
!$OMP+PRIVATE( RK2, RK, RKT, SO2_cd, L1, Ld, L2, L2S, L3, L3S, FC, LWC )
!$OMP+PRIVATE( KaqH2O2, KaqO3, AREA_CM2, FLUX, ALK, ALK1, ALK2         )
!$OMP+PRIVATE( Kt1, Kt2, AREASS1, AREASS2, SO2_ss, Kt1N, Kt2N          )
!$OMP+PRIVATE( PSO4E, PSO4F                                            )
!$OMP+PRIVATE( ADH2O20, ADL1, ADL2, ADL2S, ADL3, ADO3, ADSO20, ADSO2_CD)
!$OMP+PRIVATE( ADL3S )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP  
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake 
         AREA_CM2 = 0d0
         FLUX     = 0d0
         Ld       = 0d0

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! Initial SO2, H2O2 and O3 [v/v]
         SO20   = SO2_CHK(I,J,L)
         H2O20  = H2O2_CHK(I,J,L)
         O3     = GET_O3(I,J,L)

         ! Initialize local adjoint variables.
         ! Note: O3 influences sulfate chemistry but not vice versa, hence
         !  we initialize ADO3 to zeroo here and update the global array
         !  later. Also, adpso4_so2 got it's initial value during 
         !  ADJ_CHEM_SO4. adso2s(:,:,:) and adh2o2s(:,:,:) get their 
         !  values from adj_wetdep
         ! ( would it be faster to intialize them to zero outside the loop
         !   and then declare them FIRSTPRIVATE ? )
         ADO3     = 0d0
         ADL1     = 0d0
         ADL2     = 0d0
         ADL2S    = 0d0
         ADL3S    = 0d0
         ADL3     = 0d0
         ADSO20   = 0d0
         ADH2O20  = 0D0 
         ADSO2_CD = 0d0

         !---------------------------------------------------------
         ! Initialize parameters not affected by active variables
         !---------------------------------------------------------

         ! PATM  : Atmospheric pressure in atm
         PATM   = GET_PCENTER( I, J, L ) / 1013.25d0

         ! TK : Temperature [K]
         TK     = T(I,J,L)

         IF ( IS_OFFLINE ) THEN

            ! Not yet supported in adjoint:
            !! Gas phase SO4 production is done here in offline run only 
            !! RK1: SO2 + OH(g) [s-1]  (rjp, bmy, 3/23/03)
            !K0  = 3.0d-31 * ( 300.d0 / TK )**3.3d0
            !M   = AIRDEN(L,I,J) * F
            !KK  = K0 * M / Ki
            !F1  = ( 1.d0 + ( LOG10( KK ) )**2 )**( -1 )
            !RK1 = ( K0 * M / ( 1.d0 + KK ) ) * 0.6d0**F1 * GET_OH(I,J,L)
            !  
            CALL ERROR_STOP(' IS_OFFLINE ', 'sulfate_adj_mod.f')

         ELSE 

            ! For online runs, SMVGEAR deals w/ this computation,
            ! so we can simply set RK1 = 0 (rjp, bmy, 3/23/03)
            K0  = 0.d0
            M   = 0.d0
            KK  = 0.d0
            F1  = 0.d0
            RK1 = 0.d0

         ENDIF

         ! SO2 drydep frequency [1/s].  Also accounts for the fraction
         ! of grid box (I,J,L) that is located beneath the PBL top.
         RK2    = DEPSAV(I,J,DRYSO2) * GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! RK: total reaction rate [1/s]
         RK     = ( RK1 + RK2 )
       
         ! RKT: RK * DTCHEM [unitless] (bmy, 6/1/00)
         RKT    =  RK * DTCHEM

         ! Volume cloud fraction (Sundqvist et al 1989) [unitless]
         FC      = VCLDF(I,J,L)

         ! Liquid water content in cloudy area of grid box [m3/m3]
         LWC     = GET_LWC( TK ) * FC

         ! Zero variables
         KaqH2O2 = 0.d0
         KaqO3   = 0.d0
         L2      = 0.d0
         L3      = 0.d0
         L2S     = 0.d0
         L3S     = 0.d0

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
        ! dkh debug
        IF ( LPRINTFD .and. I == IFD .and. 
     &       J == JFD .and. L == LFD ) THEN 
           print*, ' ADJ_CHEM_SO2: rk1, rk2, rk ! ' , rk1, rk2, rk
        ENDIF 

      if (rk .gt. 0.d0) then
        so2_cd = so20*exp(-rkt)
      else
        so2_cd = so20
      endif

      ! Initialize with global arrays above.  "in" and "out" 
      ! were just dummy variables used to create the adjoint code. 
      !adpso4_so2(i,j,l) = adpso4_so2(i,j,l)+adout(3)
      !adout(3) = 0.
      !adso2s(i,j,l) = adso2s(i,j,l)+adout(2)
      !!adout(2) = 0.
      !adh2o2s(i,j,l) = adh2o2s(i,j,l)+adout(1)
      !adout(1) = 0.
      SO2s_ADJ(I,J,L) = SO2s_ADJ(I,J,L) + STT_ADJ(I,J,L,IDTSO2)
      STT_ADJ(I,J,L,IDTSO2) = 0d0
      H2O2s_ADJ(I,J,L) = H2O2s_ADJ(I,J,L) + STT_ADJ(I,J,L,IDTH2O2)
      STT_ADJ(I,J,L,IDTH2O2) = 0d0

      adl1 = adl1+adpso4_so2(i,j,l)
      adl2s = adl2s+adpso4_so2(i,j,l)
      adl3s = adl3s+adpso4_so2(i,j,l)
      adpso4_so2(i,j,l) = 0.d0


      if (fc .gt. 0.d0 .and. so2_cd .gt. mindat .and. tk .gt. 258.) then
        ! Strange, though true, this routine's output, kaqh2o2 and 
        ! kaqo3, don't actually depend upon h2o20 or o3 !
        call aqchem_so2( lwc,tk,patm,so2_cd,h2o20,o3,hplus,kaqh2o2,
     $kaqo3 )

        ! Updates to match modified fwd code (jkoo, dkh, 02/19/11) 
        l2 = exp((so2_cd-h2o20)*kaqh2o2*dtchem)
        l3 = exp((so2_cd-o3)*kaqo3*dtchem)
        !l2s = so2_cd*h2o20*(l2-1.d0)/(so2_cd*l2-h2o20)
        L2S = H2O20 + (H2O20*H2O20 - H2O20*SO2_cd) 
     &                             / ((SO2_cd * L2) - H2O20)
        !l3s = so2_cd*o3*(l3-1.d0)/(so2_cd*l3-o3)
        L3S = O3 + (O3*O3 - O3*SO2_cd) / ((SO2_cd * L3) - O3)


        ! dkh debug
        IF ( LPRINTFD .and. I == IFD .and. 
     &       J == JFD .and. L == LFD ) THEN 
           print*, ' L2 = ', L2
           print*, ' L1 = ', L1
        ENDIF 
        !adh2o20 = adh2o20+adh2o2s(i,j,l)
        !adh2o2s(i,j,l) = 0.d0
        adh2o20 = adh2o20+H2O2s_ADJ(i,j,l)
        H2O2s_ADJ(i,j,l) = 0.d0
        !adso2_cd = adso2_cd+adso2s(i,j,l)
        !adso2s(i,j,l) = 0.d0
        adso2_cd = adso2_cd+SO2s_ADJ(i,j,l)
        SO2s_ADJ(i,j,l) = 0.d0
        adl2s = adl2s-adh2o20*(0.5+sign(0.5d0,h2o20-l2s-mindat))
        adh2o20 = adh2o20*(0.5+sign(0.5d0,h2o20-l2s-mindat))
        adl2s = adl2s-adso2_cd*(0.5+sign(0.5d0,so2_cd-(l2s+l3s)-mindat))
        adl3s = adl3s-adso2_cd*(0.5+sign(0.5d0,so2_cd-(l2s+l3s)-mindat))
        adso2_cd = adso2_cd*(0.5+sign(0.5d0,so2_cd-(l2s+l3s)-mindat))
      ! Now also check to make sure L3 isn't 1d0 (dkh, 07/09/11, adj32_005) 
      IF (1/L3==0 .or. L3 == 1d0 ) THEN
        ! do nothing 
      ELSE 
        adl3 = adl3+adl3s*(so2_cd*o3/(so2_cd*l3-o3)-so2_cd*o3*(l3-1.d0)*
     $so2_cd/((so2_cd*l3-o3)*(so2_cd*l3-o3)))
        ado3 = ado3+adl3s*(so2_cd*(l3-1.d0)/(so2_cd*l3-o3)+so2_cd*o3*
     $(l3-1.d0)/((so2_cd*l3-o3)*(so2_cd*l3-o3)))
        adso2_cd = adso2_cd+adl3s*(o3*(l3-1.d0)/(so2_cd*l3-o3)-so2_cd*
     $o3*(l3-1.d0)*l3/((so2_cd*l3-o3)*(so2_cd*l3-o3)))
        adl3s = 0.d0
      ENDIF 

      ! L2 and L3 can be infinity (jkoo, 10/19/10)
!-------------------------------------------------------------------------
! old version: 
!        adh2o20 = adh2o20+adl2s*(so2_cd*(l2-1.d0)/(so2_cd*l2-h2o20)+
!     $so2_cd*h2o20*(l2-1.d0)/((so2_cd*l2-h2o20)*(so2_cd*l2-h2o20)))
!        adl2 = adl2+adl2s*(so2_cd*h2o20/(so2_cd*l2-h2o20)-so2_cd*h2o20*
!     $(l2-1.d0)*so2_cd/((so2_cd*l2-h2o20)*(so2_cd*l2-h2o20)))
!        adso2_cd = adso2_cd+adl2s*(h2o20*(l2-1.d0)/(so2_cd*l2-h2o20)-
!     $so2_cd*h2o20*(l2-1.d0)*l2/((so2_cd*l2-h2o20)*(so2_cd*l2-h2o20)))
!        adl2s = 0.d0
!        ado3 = ado3-adl3*kaqo3*dtchem*exp((so2_cd-o3)*kaqo3*dtchem)
!        adso2_cd = adso2_cd+adl3*kaqo3*dtchem*exp((so2_cd-o3)*kaqo3*
!     $dtchem)
!        adl3 = 0.d0
!        adh2o20 = adh2o20-adl2*kaqh2o2*dtchem*exp((so2_cd-h2o20)*
!     $kaqh2o2*dtchem)
!        adso2_cd = adso2_cd+adl2*kaqh2o2*dtchem*exp((so2_cd-h2o20)*
!     $kaqh2o2*dtchem)
!        adl2 = 0.d0
!-------------------------------------------------------------------------
      ! Updates
      ! Now also check to make sure L2 isn't 1d0, which can happen, 
      ! rarely, but it causes NaNs (dkh, 07/09/11, adj32_005) 
      !IF (1/L2==0) THEN
      IF (1/L2==0 .or. L2 == 1d0) THEN
        ADH2O20 = ADH2O20+ADL2S
        ADL2 = ADL2
        ADSO2_CD = ADSO2_CD
      ELSE
        adh2o20 = adh2o20+adl2s*(so2_cd*(l2-1.d0)/(so2_cd*l2-h2o20)+
     $so2_cd*h2o20*(l2-1.d0)/((so2_cd*l2-h2o20)*(so2_cd*l2-h2o20)))
        adl2 = adl2+adl2s*(so2_cd*h2o20/(so2_cd*l2-h2o20)-so2_cd*h2o20*
     $(l2-1.d0)*so2_cd/((so2_cd*l2-h2o20)*(so2_cd*l2-h2o20)))
        adso2_cd = adso2_cd+adl2s*(h2o20*(l2-1.d0)/(so2_cd*l2-h2o20)-
     $so2_cd*h2o20*(l2-1.d0)*l2/((so2_cd*l2-h2o20)*(so2_cd*l2-h2o20)))
        adl2s = 0.d0
        adh2o20 = adh2o20-adl2*kaqh2o2*dtchem*L2
        adso2_cd = adso2_cd+adl2*kaqh2o2*dtchem*L2
        adl2 = 0.d0
      ENDIF
      !IF (1/L3==0) THEN
      ! Now also check to make sure L3 isn't 1d0 (dkh, 07/09/11, adj32_005) 
      !IF (1/L3==0) THEN
      IF (1/L3==0 .or. L3 == 1d0 ) THEN
        ! do nothing
      ELSE
        ado3 = ado3-adl3*kaqo3*dtchem*L3
        adso2_cd = adso2_cd+adl3*kaqo3*dtchem*L3
        adl3 = 0.d0
      ENDIF
!-------------------------------------------------------------------------


      else
!        adso2_cd = adso2_cd+adso2s(i,j,l)*(0.5+sign(0.5d0,so2_cd-1.d-32)
!     $)
!        adso2s(i,j,l) = 0.d0
        adso2_cd = adso2_cd+SO2s_ADJ(i,j,l)*
     &(0.5+sign(0.5d0,so2_cd-1.d-32))
        SO2s_ADJ(i,j,l) = 0.d0
        !adh2o20 = adh2o20+adh2o2s(i,j,l)*(0.5+sign(0.5d0,h2o20-1.d-32))
        !adh2o2s(i,j,l) = 0.d0
        adh2o20 = adh2o20+H2O2s_ADJ(i,j,l)
     &            *(0.5+sign(0.5d0,h2o20-1.d-32))
        H2O2s_ADJ(i,j,l) = 0.d0
        ! dkh debug
        IF ( LPRINTFD .and. I == IFD .and. 
     &       J == JFD .and. L == LFD ) THEN 
           print*, ' ADJ_CHEM_SO2: no aqchem ! ' 
        ENDIF 
      endif
      if (rk .gt. 0.d0) then
        adso20 = adso20+adl1*(rk1/rk)
        adso2_cd = adso2_cd-adl1*(rk1/rk)
        adl1 = 0.d0
        adso20 = adso20+adso2_cd*exp(-rkt)
        adso2_cd = 0.d0
      else
        adso20 = adso20+adso2_cd
        adso2_cd = 0.d0
      endif
      
      ! Update global adjoint arrays and reset local variables to zero.
      ! ADJ_STT(IDADJSO2) should be zero at this point, as it will be zeroed
      ! after passing off values to ADSO2s ? -- NO, they both coexist.  H2O2s
      ! and SO2s are just control variables. The actual running concentrations
      ! are kept in STT. 
      STT_ADJ(I,J,L,IDTSO2)   = STT_ADJ(I,J,L,IDTSO2) + ADSO20
      STT_ADJ(I,J,L,IDTH2O2)  = STT_ADJ(I,J,L,IDTH2O2) + ADH2O20
      CALL GET_O3_ADJ( ADO3, I, J, L ) 

      ado3 = 0.d0
      adh2o20 = 0.d0
      adso20 = 0.d0

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_SO2_ADJ

!------------------------------------------------------------------------------
!
!      SUBROUTINE SEASALT_CHEM ( I,      J,     L,   ALK1, ALK2,
!     &                          SO2_cd, Kt1,   Kt2, Kt1N, Kt2N,
!     &                          SO2_ss, PSO4E, PSO4F )
!!
!!******************************************************************************
!!  Function SEASALT_CHEM computes SO4 formed from S(IV) + O3 on seasalt 
!!  aerosols as a function of seasalt alkalinity. (bec, bmy, 4/13/05, 10/7/08)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) I      (INTEGER) :
!!  (2 ) J      (INTEGER) :
!!  (3 ) L      (INTEGER) :
!!  (4 ) O30    (TYPE (XPLEX))  : Initial O3 mixing ratio (v/v]
!!  (5 ) ALK    (TYPE (XPLEX))  : Alkalinity [kg] from seasalt from seasalt_mod
!!  (6 ) SO2_cd (TYPE (XPLEX))  : SO2 mixing ratio [v/v] after gas phase chemistry
!!                           and dry deposition
!!  (7 ) Kt1    (TYPE (XPLEX))  : Rate constant [s-1] for sulfate formation on 
!!                           fine sea-salt aerosols from GET_ALK
!!  (8 ) Kt2    (TYPE (XPLEX))  : Rate constant [s-1] for sulfate formation on 
!!                           coarse sea-salt aerosols from GET_ALK
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (9 ) SO2_ss 	(TYPE (XPLEX)) : SO2 mixing ratio [v/v] updated after SS chemistry
!!  (10) SO4E    (TYPE (XPLEX)) : SO4E (sulfate produced by S(IV)+O3 on fine seasalt)
!!                          mixing ratio [v/v]
!!  (11) SO4F    (TYPE (XPLEX)) : SO4F(sulfate produced by S(IV)+O3 on coarse seasalt)
!!                          mixing ratio [v/v]
!!  (12) O3      (TYPE (XPLEX)) : Updated O3 mixing ratio [v/v] for fullchem runs 
!!                           only. Otherwise O30=O3.
!!
!!  Chemical reactions:
!!  ============================================================================
!!  (R1) SO2 + O3 + ALK => SO4 + O2
!!       Modeled after Chamedies and Stelson, 1992?
!!
!!  NOTES:
!!  (1 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!!  (2 ) Bug fix: now avoid seg fault error if IDTHNO3 is zero, as it would
!!        be for an offline aerosol simulation. (bmy, 3/29/06)
!!  (3 ) Fixed typo in FALK_A_SO2 equation: C_FLUX_C should be C_FLUX_A.
!!        (havala, bec, bmy, 12/8/06)
!!  (4 ) Bug fix for mass balance, replace TITR_HNO3 w/ HNO3_SSC in the
!!        expression for HNO3_ss.  Bug fix: now do equivalent computation 
!!        for GET_GNO3, which is now no longer called because it's in 
!!        "isoropia_mod.f". (bec, bmy, 7/30/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE COMODE_MOD,      ONLY : CSPEC, JLOP, VOLUME
!      USE DAO_MOD,         ONLY : AD, AIRDEN, AIRVOL
!      USE TRACERID_MOD
!      !----------------------------------------------------------------
!      ! DIAGNOSTIC -- leave commented out for now (bec, bmy, 4/13/05)
!      !USE DIAG_MOD,        ONLY : AD09
!      !----------------------------------------------------------------
!      USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
!      USE TIME_MOD,        ONLY : GET_TS_CHEM,        GET_ELAPSED_SEC
!      USE ERROR_MOD,       ONLY : IT_IS_NAN
!      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM, STT
!      USE TRACER_MOD,      ONLY : TCVV,               XNUMOLAIR
!      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
!      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC,    GET_MONTH 
!      USE TIME_MOD,        ONLY : ITS_A_NEW_MONTH
!      USE ISOROPIA_MOD,    ONLY : GET_GNO3
! 
!      ! Add these for GET_GNO3 fix (lyj, bmy, 10/7/08)
!      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
!      USE DAO_MOD,         ONLY : AIRVOL
!
!#     include "CMN_SIZE"  ! Size parameters
!!---------------------------------------------------------------
!! DIAGNOSTIC -- leave commented out for now (bec, bmy, 4/13/05)
!!#     include "CMN_DIAG"  ! ND19
!!---------------------------------------------------------------
!#     include "CMN_GCTM"  ! AIRMW
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: I, J, L 
!      TYPE (XPLEX),  INTENT(IN)  :: SO2_cd, Kt1, Kt2, Kt1N, Kt2N
!      TYPE (XPLEX),  INTENT(IN)  :: ALK1, ALK2
!      TYPE (XPLEX),  INTENT(OUT) :: SO2_ss, PSO4E, PSO4F
!
!      ! Local variables
!      INTEGER              :: JLOOP
!      TYPE (XPLEX)               :: SO2_chem,    DTCHEM
!      TYPE (XPLEX)               :: O3_cspec,    O3_lost
!      TYPE (XPLEX)               :: EQ_1_C,      EQ_2_C
!      TYPE (XPLEX)               :: SO4E,        SO2_new,    SO4F
!      TYPE (XPLEX)               :: SO2_eq,      N_FLUX_A,   N_FLUX_C
!      TYPE (XPLEX)               :: END_ALK,     L5A,        L5C
!      TYPE (XPLEX)               :: EQ1,         EQ2,        TITR_SO2
!      TYPE (XPLEX)               :: TITR_HNO3,   NIT_vv,     NITs_vv
!      TYPE (XPLEX)               :: NIT0,        NITS0
!      TYPE (XPLEX)               :: F_SO2,       FALK_A_SO2, FALK_C_SO2
!      TYPE (XPLEX)               :: EQ_BEG,      F_SO2_A,    F_SO2_C
!      TYPE (XPLEX)               :: ALKA,        ALKC,       TOTAL_ACID_FLUX
!      TYPE (XPLEX)               :: HNO3_EQ,     TOT_FLUX_A, TOT_FLUX_C
!      TYPE (XPLEX)               :: FALK_A_HNO3, HNO3_vv
!      TYPE (XPLEX)               :: FALK_C_HNO3, F_HNO3_A,   F_HNO3_C
!      TYPE (XPLEX)               :: EQ_1_N,      EQ_2_N,     F_HNO3
!      TYPE (XPLEX)               :: HNO3_SSA,    HNO3_SSC,   N_FLUX
!      TYPE (XPLEX)               :: HNO3_EQ_C,   L6A,        L6C   
!      TYPE (XPLEX)               :: C_FLUX_A,    C_FLUX_C,   C_FLUX
!      TYPE (XPLEX)               :: HNO3_ss,     HNO3_kg
!      TYPE (XPLEX),  PARAMETER   :: MINDAT    = 1.0d-20 
!      TYPE (XPLEX),  PARAMETER   :: TCVV_HNO3 = 28.97d0 / 63.0d0 
!
!      !=================================================================
!      ! SEASALT_CHEM begins here!
!      !=================================================================
!
!      ! DTCHEM is the chemistry timestep in seconds
!      DTCHEM = GET_TS_CHEM() * 60d0
!
!      ! Convert SO2 [v/v] to  [eq/gridbox]
!      SO2_eq = ( 2.d0 * SO2_cd * AD(I,J,L) ) / ( TCVV(IDTSO2) * 0.064d0)
!      SO2_eq = MAX( SO2_eq, MINDAT )
!
!      IF ( ITS_A_FULLCHEM_SIM() ) THEN
!
! 	 ! Convert HNO3 [v/v] to [equivalents]
!         HNO3_vv = STT(I,J,L,IDTHNO3)
!         HNO3_eq = HNO3_vv * AD(I,J,L) / ( 28.97d0 / 63d0 ) / 63.d-3
!
!      ELSE
!         
!         !--------------------------------------------------------------------
!         ! Prior to 10/7/08:
!         ! Now that we have switched from ISORROPIA to RPMARES, GET_GNO3
!         ! is no longer defined.  We therefore need to do the equivalent
!         ! computation.  NOTE: This is only an issue for the offline
!         ! aerosol simulation. (lyj, bmy, 10/7/08)
!         ! 
!         ! For more information, please see this Wiki post:
!         ! http://wiki.seas.harvard.edu/geos-chem/index.php/Aerosol_thermodynamical_equilibrium#Bug_in_sulfate_mod.f_caused_by_switch_to_RPMARES
!         ! 
!         !! Get gas-phase HNO3 from ISORROPIA code
!         !CALL GET_GNO3( I, J, L, HNO3_kg )         
!         !--------------------------------------------------------------------
!         
!         ! Get HNO3 in ug/m3, then multiply by volume in m3
!         ! and 1e-6 kg/ug to get HNO3 in kg
!         HNO3_kg = GET_HNO3_UGM3( I, J, L ) * AIRVOL(I,J,L) * 1e-6
!
!	 ! Convert HNO3 [kg] first to [v/v] 
!         HNO3_vv = HNO3_kg * ( 28.97d0 / 63d0 )   / AD(I,J,L)
! 
!         ! Then convert HNO3 [kg] to [equivalents]
!         HNO3_eq = HNO3_kg / 63d-3
!
!      ENDIF
!
!      !-----------
!      ! SO2
!      !-----------
!
!      ! Available flux of SO2 to accum sea salt aerosols [v/v/timestep]
!      L5A      = EXP( -Kt1 * DTCHEM )
!      F_SO2_A  = SO2_cd * ( 1.d0 - L5A )
!      F_SO2_A  = MAX( F_SO2_A, 1.d-32 )
!
!      ! Convert to [eq/timestep] 
!      C_FLUX_A = 2.d0 * F_SO2_A * AD(I,J,L) / TCVV(IDTSO2) / 0.064d0
!
!      ! Available flux of SO2 to coarse sea salt aerosols [v/v/timestep]
!      L5C      = EXP( - Kt2 * DTCHEM )
!      F_SO2_C  = SO2_cd * ( 1.d0 - L5C )
!      F_SO2_C  = MAX( F_SO2_C, 1.0d-32 )
!
!      ! Convert to [eq/timestep] 
!      C_FLUX_C = 2.d0 * F_SO2_C * AD(I,J,L) / TCVV(IDTSO2) / 0.064d0
!
!      ! Total flux of SO2 [v/v/timestep]
!      F_SO2    = F_SO2_A + F_SO2_C 
!
!      ! Total flux of SO2 [eq/timestep]
!      C_FLUX   = C_FLUX_A + C_FLUX_C 
!
!      !-----------
!      ! HNO3
!      !-----------
!
!      ! Available flux of HNO3 to accum sea salt aerosols [v/v/timestep]
!      L6A = EXP( - Kt1N * DTCHEM )
!      F_HNO3_A = HNO3_vv * ( 1.D0 - L6A )
!      F_HNO3_A = MAX( F_HNO3_A, 1.0D-32 )
!
!      ! Convert to [eq/timestep] 
!      N_FLUX_A = F_HNO3_A * AD(I,J,L)/( 28.97d0 / 63d0 )/0.063d0
!
!      ! Available flux of HNO3 to coarse sea salt aerosols [v/v/timestep]
!      L6C = EXP( - Kt2N * DTCHEM )
!      F_HNO3_C = HNO3_vv * ( 1.D0 - L6C )
!      F_HNO3_C = MAX( F_HNO3_C, 1.0D-32 )
!
!      ! convert to [eq/timestep] 
!      N_FLUX_C = F_HNO3_C * AD(I,J,L)/( 28.97d0 / 63d0 )/0.063d0
!
!      ! Total flux of HNO3
!      F_HNO3 = F_HNO3_A + F_HNO3_C ![v/v/timestep]
!      N_FLUX = N_FLUX_A + N_FLUX_C ![eq/timestep]
!
!      !-----------
!      ! Acid
!      !-----------
!
!      ! Total acid flux to accum sea-salt aerosols [eq/box/timestep]
!      TOT_FLUX_A = C_FLUX_A + N_FLUX_A 
!      TOT_FLUX_A = MAX( TOT_FLUX_A, MINDAT )
!
!      ! Total acid flux to coarse sea-salt aerosols [eq/box/timestep]
!      TOT_FLUX_C = C_FLUX_C + N_FLUX_C 
!      TOT_FLUX_C = MAX( TOT_FLUX_C, MINDAT )
!
!      ! Total  acid flux to sea salt aerosols
!      TOTAL_ACID_FLUX = TOT_FLUX_A + TOT_FLUX_C
!
!      ! Total available alkalinity [eq]
!      EQ1 = ALK1 * 0.07d0
!      EQ2 = ALK2 * 0.07d0
!
!      !----------------------------------------------------------------------
!      ! NOTE: This was a sensitivity simulation, keep for future reference
!      !       cf Alexander et al 2005 (bec, bmy, 4/13/05)
!      !! Total available alkalinity [eq] doubled for Sievering run
!      !EQ1 = ALK1 * 0.14d0
!      !EQ2 = ALK2 * 0.14d0
!      !----------------------------------------------------------------------
!
!      !----------------------------------------------------------------------
!      ! DIAGNOSTIC -- leave uncommented for now (bec, bmy, 4/13/05)
!      !! Write out beginning alkalinity [eq SO4]
!      !EQ_BEG = EQ1 + EQ2
!      !IF ( ND09 > 0 ) AD09(I,J,L,1) = AD09(I,J,L,1) + EQ_BEG
!      !----------------------------------------------------------------------
!
!      IF ( TOT_FLUX_A > EQ1 ) THEN
!
!	 ! Fraction of alkalinity available for each acid
!         FALK_A_SO2  = C_FLUX_A / TOT_FLUX_A
!	 FALK_A_HNO3 = N_FLUX_A / TOT_FLUX_A
!         FALK_A_SO2  = MAX( FALK_A_SO2, MINDAT )
!	 FALK_A_HNO3 = MAX( FALK_A_HNO3, MINDAT )
!
!      ELSE
!
!	 FALK_A_SO2  = 1.0d0
!	 FALK_A_HNO3 = 1.0d0
!
!      ENDIF
!      
!      IF ( TOT_FLUX_C > EQ2 ) THEN
!
!         ! Fraction of flkalinity available for each acid
!	 FALK_C_SO2  = C_FLUX_C/TOT_FLUX_C
!	 FALK_C_HNO3 = N_FLUX_C/TOT_FLUX_C
!         FALK_C_SO2  = MAX( FALK_C_SO2, MINDAT )
!	 FALK_C_HNO3 = MAX( FALK_C_HNO3, MINDAT )
!
!      ELSE
!
!	 FALK_C_SO2  = 1.0d0
!	 FALK_C_HNO3 = 1.0d0
!
!      ENDIF
!
!      ! Alkalinity available for S(IV) --> S(VI)
!      EQ_1_C       = EQ1 * FALK_A_SO2
!      EQ_1_C       = MAX( EQ_1_C, MINDAT )
!      EQ_1_N       = EQ1 * FALK_A_HNO3
!      EQ_1_N       = MAX( EQ_1_N, MINDAT )
!                  
!      EQ_2_C       = EQ2 * FALK_C_SO2
!      EQ_2_C       = MAX( EQ_2_C, MINDAT )
!      EQ_2_N       = EQ2 * FALK_C_HNO3
!      EQ_2_N       = MAX( EQ_2_N, MINDAT )
!
!      !-----------------
!      ! Fine Seasalt
!      !-----------------
!
!      ! don't produce more SO4 than available ALK or SO2
!      SO4E         = MIN( C_FLUX_A, EQ_1_C, SO2_eq ) 
!      SO4E         = MAX( SO4E, MINDAT )
!
!      ! Update SO2 concentration [eq/box] 
!      SO2_new      = SO2_eq - SO4E
!      SO2_new      = MAX( SO2_new, MINDAT )
!
!      !-----------------
!      ! Coarse Seasalt
!      !-----------------     
!      IF ( SO2_new > MINDAT ) THEN
!
! 	 ! don't produce more SO4 than available ALK or SO2
!	 SO4F      = MIN( C_FLUX_C, SO2_new, EQ_2_C ) 
!	 SO4F      = MAX( SO4F, MINDAT )
!
!	 !Update SO2 concentration [eq] 
!	 SO2_chem  = SO2_new - SO4F
!	 SO2_chem  = MAX( SO2_chem, MINDAT )
!      ELSE
!	 SO4F      = MINDAT
!	 SO2_chem  = MINDAT
!      ENDIF
!
!      ! Alkalinity titrated by S(IV) --> S(VI) [eq]
!      TITR_SO2     = SO4E + SO4F
!
!      !-------------------------------------------------------------------
!      ! DIAGNOSTIC -- leave uncommented for now
!      !! write out in diagnostic
!      !IF ( ND09 > 0 ) AD09(I,J,L,2) = AD09(I,J,L,2) + TITR_SO2
!      !-------------------------------------------------------------------
!
!      !Modified SO2 [eq] converted back to [v/v]
!      SO2_ss       = SO2_chem * 0.064 * TCVV(IDTSO2)/AD(I,J,L)/2.0d0
!      SO2_ss       = MAX( SO2_ss, MINDAT )
!
!      !SO4E produced converted from [eq/timestep] to [v/v/timestep] 
!      PSO4E        = SO4E * 0.096 * TCVV(IDTSO4)/AD(I,J,L)/2.0d0
!
!      !SO4F produced converted from [eq/timestep] to [v/v/timestep] 
!      PSO4F        = SO4F * 0.096 * TCVV(IDTSO4S)/AD(I,J,L)/2.0d0
!
!      ! Alkalinity titrated by HNO3
!      HNO3_SSA     = MIN(N_FLUX_A, HNO3_EQ, EQ_1_N)
!      HNO3_SSA     = MAX(HNO3_SSA, MINDAT)
!      HNO3_EQ_C    = HNO3_EQ - HNO3_SSA
!      HNO3_EQ_C    = MAX(HNO3_EQ_C, MINDAT)
!      HNO3_SSC     = MIN(N_FLUX_C, HNO3_EQ_C, EQ_2_N)
!      HNO3_SSC     = MAX(HNO3_SSC, MINDAT)
!      TITR_HNO3    = HNO3_SSA + HNO3_SSC
!
!      !----------------------------------------------------------------------
!      ! DIAGNOSTIC -- leave commented out for now
!      ! !write out alkalinity titrated by HNO3(g)
!      !IF ( ND09 > 0 ) AD09(I,J,L,3) = AD09(I,J,L,3) + TITR_HNO3
!      !----------------------------------------------------------------------
!
!      ! HNO3 lost [eq/timestep] converted back to [v/v/timestep]
!      IF ( IDTHNO3 > 0 ) THEN
!
!         ! Coupled sim: IDTHNO3 is defined, so use it
!         HNO3_ss = HNO3_SSC * 0.063 * TCVV(IDTHNO3)/AD(I,J,L)
!         STT(I,J,L,IDTHNO3) = MAX( HNO3_vv - HNO3_ss, MINDAT )
!
!      ELSE
!
!         ! Offline aerosol sim: IDTHNO3 isn't defined, use TCVV_HNO3
!         HNO3_ss = TITR_HNO3 * 0.063 * TCVV_HNO3 / AD(I,J,L)
!
!      ENDIF
!
!      ! NITS produced converted from [eq/timestep] to [v/v/timestep] 
!      PNITs(I,J,L) = HNO3_SSC * 0.063 * TCVV(IDTNITS)/AD(I,J,L)
!	 
!      ! Modified accum alkalinity 
!      ALKA         = EQ1 - (SO4E + HNO3_SSA)
!      ALKA         = MAX( ALKA, MINDAT )
!
!      !------------------------------------------------------------------------
!      ! Uncomment this if you want to transport alkalinity (bec, bmy, 4/13/05)
!      ![eq] --> [kg] --> [v/v] use this only if transporting alkalinity
!      !ALKAvv = (ALKA * TCVV(IDTSAL1))/(7.0d-2 * AD(I,J,L))
!      !ALKAvv = MAX( ALKAvv, MINDAT )
!      !------------------------------------------------------------------------
!
!      ! Modified accum alkalinity 
!      ALKC         = EQ2 - (SO4F + HNO3_SSC)
!      ALKC         = MAX( ALKC, MINDAT )
!      
!      !------------------------------------------------------------------------
!      ! Uncomment this if you want to transport alkalinity (bec, bmy, 4/13/05)
!      !! [eq] --> [kg] --> [v/v] use this only if transporting alkalinity
!      !ALKCvv = (ALKC * TCVV(IDTSAL2))/(7.0d-2 * AD(I,J,L))
!      !ALKCvv = MAX(ALKCvv, MINDAT)
!      !------------------------------------------------------------------------
!
!      !------------------------------------------------------------------------
!      ! DIAGNOSTIC -- leave commented out for now (bec, bmy, 4/13/05)
!      !! write out ending alkalinity
!      !END_ALK = ALKA + ALKC
!      !IF ( ND09 > 0 ) AD09(I,J,L,4) = AD09(I,J,L,4) + END_ALK
!      !------------------------------------------------------------------------
!
!      ! Return to calling program
!      END SUBROUTINE SEASALT_CHEM
!
!!------------------------------------------------------------------------------

      SUBROUTINE AQCHEM_SO2( LWC, T,     P,       SO2, H2O2, 
     &                       O3,  Hplus, KaqH2O2, KaqO3 ) 
!
!******************************************************************************
!  Function AQCHEM_SO2 computes the reaction rates for aqueous SO2 chemistry.
!  (rjp, bmy, 10/31/02, 12/12/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LWC     (TYPE (XPLEX)) : Liquid water content [m3/m3] = 1.E-6*L [g/m3]
!  (2 ) T       (TYPE (XPLEX)) : Temperature [K]
!  (3 ) P       (TYPE (XPLEX)) : Pressure [atm]
!  (4 ) SO2     (TYPE (XPLEX)) : SO2  mixing ratio [v/v]
!  (5 ) H2O2    (TYPE (XPLEX)) : H2O2 mixing ratio [v/v]
!  (6 ) O3      (TYPE (XPLEX)) : O3   mixing ratio [v/v]
!  (7 ) HPLUS   (TYPE (XPLEX)) : Concentration of H+ ion (i.e. the pH) [v/v]
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) KaqH2O2 (TYPE (XPLEX)) : Reaction rate for H2O2
!  (9 ) KaqO3   (TYPE (XPLEX)) : Reaction rate for O3
!
!  Chemical Reactions:
!  ============================================================================
!  (R1) HSO3- + H2O2(aq) + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]   
!
!      d[S(VI)]/dt = k[H+][H2O2(aq)][HSO3-]/(1 + K[H+]) 
!      [Seinfeld and Pandis, 1998, page 366]
!
!  (R2) SO2(aq) + O3(aq) =>                                        
!       HSO3-   + O3(aq) =>  
!       SO3--   + O3(aq) =>
!       [Jacob, 1986; Jacobson, 1999]
!
!       d[S(VI)]/dt = (k0[SO2(aq)] + k1[HSO3-] + K2[SO3--])[O3(aq)]
!       [Seinfeld and Pandis, 1998, page 363]
!
!  Reaction rates can be given as
!       Ra     = k [H2O2(ag)] [S(IV)]  [mole/liter*s]  OR
!       Krate  = Ra LWC R T / P        [1/s]
!
!  Where:
!       LWC = Liquid water content(g/m3)*10-6 [m3(water)/m3(gas)]
!       R   = 0.08205  (atm L / mol-K), Universal gas const.
!       T   = Temperature (K)
!       P   = Pressure (atm)
!
!  Procedure:
!  ============================================================================
!  (a ) Given [SO2] which is assumed to be total SO2 (gas+liquid) in 
!        equilibrium between gas and liquid phase. 
!
!  (b ) We can compute SO2(g) using Henry's law 
!          P(so2(g)) = Xg * [SO2]
!          Xg = 1/(1 + Faq), Fraction of SO2 in gas
!       where: 
!          Faq   = Kheff * R * T * LWC, 
!          KHeff = Effective Henry's constant
!
!  (c ) Then Calculate Aquous phase, S[IV] concentrations
!        S[IV] = Kheff * P(so2(g) in atm) [M]
!
!  (d ) The exact same procedure is applied to calculate H2O2(aq)
!
!  NOTES:
!  (1 ) Updated by Rokjin Park (rjp, bmy, 12/12/02)
!******************************************************************************
!
      ! Arguments
      TYPE (XPLEX), INTENT(IN)  :: LWC, T, P, SO2, H2O2, O3, HPLUS
      TYPE (XPLEX), INTENT(OUT) :: KaqH2O2, KaqO3

      ! Local variables
      TYPE (XPLEX), PARAMETER   :: R = xplex(0.08205d0,0d0) 
      TYPE (XPLEX)              :: KH2O2,  RA,     KS1, KS2,    HCSO2
      TYPE (XPLEX)              :: FHCSO2, XSO2G,  SIV, HSO3,   XSO2AQ
      TYPE (XPLEX)              :: XHSO3,  XSO3,   KH1, HCH2O2, FHCH2O2
      TYPE (XPLEX)              :: XH2O2G, H2O2aq, KO0, KO1,    KO2
      TYPE (XPLEX)              :: HCO3,   XO3g,   O3aq

      !=================================================================
      ! AQCHEM_SO2 begins here!
      !
      ! Aqueous reaction rate
      ! HSO3- + H2O2 + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
      !=================================================================

      ! [Jacob, 1986]
      KH2O2 = 6.31d14 * EXP( -4.76d3 / T )  

!      ! [Jacobson, 1999]
!      KH2O2 = 7.45d07 * EXP( -15.96d0 * ( (298.15/T) - 1.) ) / 
!     &                  ( 1.d0 + 13.d0 * Hplus)

      !=================================================================
      ! Equilibrium reaction of SO2-H2O
      !    SO2 + H2O = SO2(aq)        (s0)
      !    SO2(ag)   = HSO3- + H+     (s1)
      !    HSO3-     = SO3-- + H+     (s2)
      !
      ! Reaction constant for Aqueous chemistry -- No big difference 
      ! between Jacob and Jacobson, choose one of them.
      !
      ! Reaction rate dependent on Temperature is given
      !   H = A exp ( B (T./T - 1) ) 
      !
      ! For equilibrium reactions of SO2:
      !            As1      Bs1   As2      Bs2  
      !  Seinfeld  1.30d-2  7.02  6.60d-8  3.76   [1998]
      !  Jacob     1.30d-2  6.75  6.31d-8  5.05   [1986]
      !  Jacobson  1.71d-2  7.04  5.99d-8  3.74   [1996]
      !=================================================================
      Ks1    = 1.30d-2 * EXP( 6.75d0 * ( 298.15d0 / T - 1.d0 ) )
      Ks2    = 6.31d-8 * EXP( 5.05d0 * ( 298.15d0 / T - 1.d0 ) )

      ! SIV Fraction
      XSO2aq = 1.d0/(1.d0 + Ks1/Hplus + Ks1*Ks2/(Hplus*Hplus))
      XHSO3  = 1.d0/(1.d0 + Hplus/Ks1 + Ks2/Hplus)
      XSO3   = 1.d0/(1.d0 + Hplus/Ks2 + Hplus*Hplus/(Ks1*Ks2))

      ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
      HCSO2  = 1.22d0 * EXP( 10.55d0 * ( 298.15d0 / T - 1.d0) )         
      FHCSO2 = HCSO2 * (1.d0 + (Ks1/Hplus) + (Ks1*Ks2 / (Hplus*Hplus)))
      
      XSO2g  = 1.d0 / ( 1.d0 + ( FHCSO2 * R * T * LWC ) )
      SIV    = FHCSO2 * XSO2g * SO2 * P
!      HSO3   = Ks1 * HCSO2 * XSO2g * SO2 * P

      !=================================================================
      ! H2O2 equilibrium reaction
      ! H2O2 + H2O = H2O2.H2O
      ! H2O2.H2O   = HO2- + H+   1)
      !
      ! Reaction rate dependent on Temperature is given
      !   H = A exp ( B (T./T - 1) ) 
      !
      ! For equilibrium reactions of SO2
      !            Ah1       Bh1
      !  Jacob     1.58E-12  -12.49  [1986]
      !  Jacobson  2.20E-12  -12.52  [1996]
      !=================================================================
      Kh1 = 2.20d-12 * EXP( -12.52d0 * ( 298.15d0 / T - 1.d0 ) )

      ! Henry's constant [mol/l-atm] and Effective Henry's constant for H2O2
      ! [Seinfeld and Pandis, 1998]
      ! HCH2O2  = 7.45D4 * EXP( 24.48d0 * ( 298.15d0 / T - 1.d0) ) 

      ! [Jacobson,1999]
      HCH2O2  = 7.45D4 * EXP( 22.21d0 * (298.15d0 / T - 1.d0) )
      FHCH2O2 = HCH2O2 * (1.d0 + (Kh1 / Hplus))

      XH2O2g  = 1.d0 / ( 1.d0 + ( FHCH2O2 * R * T * LWC ) )
!      H2O2aq  = FHCH2O2 * XH2O2g * H2O2 * P

      ! Conversion rate from SO2 to SO4 via reaction with H2O2
      KaqH2O2  = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2g
     &         * P * LWC * R * T            ! [v/v/s]

      !=================================================================
      !  Aqueous reactions of SO2 with O3
      !  SO2(aq) + O3 =>                       (0)
      !  HSO3-   + O3 => SO4-- + H+ + O2       (1)
      !  SO3--   + O3 => SO4-- + O2            (2)
      !
      ! NOTE
      ! [Jacob, 1986]
      !    KO1  = 3.49E12 * EXP( -4.83E3 / T )  
      !    KO2  = 7.32E14 * EXP( -4.03E3 / T )
      !
      ! [Jacobson, 1999]
      !    KO0  = 2.40E+4
      !    KO1  = 3.70E+5 * EXP( -18.56 * ((298.15/T) - 1.))
      !    KO2  = 1.50E+9 * EXP( -17.72 * ((298.15/T) - 1.))
      !
      ! Rate constants from Jacobson is larger than those of Jacob
      ! and results in faster conversion from S(IV) to S(VI)
      ! We choose Jacob 1) 2) and Jacobson 0) here
      !=================================================================
      KO0 = 2.40d+4
      KO1 = 3.49d12 * EXP( -4.83d3 / T )  
      KO2 = 7.32d14 * EXP( -4.03d3 / T )

      !=================================================================
      ! H2O2 equilibrium reaction
      ! O3 + H2O = O3.H2O
      !  HCO3  = 1.13E-2 * EXP( 8.51 * (298.15/T -1.) ), S & P
      !  HCO3  = 1.13E-2 * EXP( 7.72 * (298.15/T -1.) ), Jacobson
      !=================================================================

      ! Calculate Henry's Law constant for atmospheric temperature
      HCO3  = 1.13d-2 * EXP( 8.51d0 * ( 298.15d0 / T - 1.d0 ) )

      XO3g  = 1.d0 / ( 1.d0 + ( HCO3 * R * T * LWC ) )
!      O3aq  = HCO3 * XO3g * O3 * P
      
      ! Conversion rate from SO2 to SO4 via reaction with O3
      KaqO3 = (KO0*XSO2AQ + KO1*XHSO3 + KO2*XSO3) * FHCSO2 * XSO2g
     &      * P * HCO3 * XO3g * LWC * R * T   ! [v/v/s]

      ! Return to calling program
      END SUBROUTINE AQCHEM_SO2

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_SO4_ADJ
!
!******************************************************************************
!  Subroutine CHEM_SO4_ADJ is the adjoint SO4 chemistry subroutine for the 
!  fwd routine CHEM_SO4.  (dkh, 10/13/05)
!
!  * Does not yet include seasalt sulfate * 
!
!  See original routine for notes.  (rjp, bdf, cas, bmy, 5/31/00, 5/23/06) 
!
!  Module Variables Used:
!  ============================================================================
!  (1 ) ADPSO4_SO2 (TYPE (XPLEX) ) : Array for adjoint of P(SO4) from SO2 
!                                                                           
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!  The Only production is from SO2 oxidation (save in CHEM_SO2), and the only  
!  loss is dry depsition here.  Wet deposition will be treated in "wetdep.f".
!                                                                          
!  SO4 = SO4_0 * exp(-kt) + PSO4_SO2/kt * (1.-exp(-kt))                    
!    where k = dry deposition.                                             
!                      
!  NOTES:              
!  
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,    ONLY : LCRYST, LSSALT
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT,    TCVV,     XNUMOL
      USE TRACERID_MOD,   ONLY : IDTSO4, IDTSO4s,  IDTAS,   IDTAHS 
      USE TRACERID_MOD,   ONLY : IDTLET, IDTSO4aq, IDTNH4aq
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER               :: I,      J,      L,        N,     N_ND44
      TYPE (XPLEX)             :: AS,     AS0,    AHS,      AHS0,  LET   
      TYPE (XPLEX)           :: LET0,   SO4,    SO40,     SO4aq, SO4aq0
      TYPE (XPLEX)            :: SO4s,   SO40s,  RKT,      RKTs,  E_RKT
      TYPE (XPLEX)                :: E_RKTs, DTCHEM, AREA_CM2, FLUX
      TYPE (XPLEX)                :: ADJ_SO4, ADJ_SO40

      !=================================================================
      ! CHEM_SO4_ADJ begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTSO4 == 0 .or. IDTSO4s == 0 ) RETURN
      IF ( DRYSO4 == 0 .or. DRYSO4s == 0 ) RETURN

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,    L,   AREA_CM2, RKT,   RKTs, E_RKT )
!$OMP+PRIVATE( E_RKTs, FLUX, SO4, SO4s,     SO40,  SO40s       )
!$OMP+PRIVATE( ADJ_SO40, ADJ_SO4 )
!$OMP+SCHEDULE( DYNAMIC ) 
      DO L = 1, LLTROP 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake
         AREA_CM2 = 0d0
         RKT      = 0d0
         RKTs     = 0d0
         E_RKT    = 0d0
         E_RKTs   = 0d0
         FLUX     = 0d0
         SO4      = 0d0
         SO4s     = 0d0


         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         !==============================================================
         ! Initial concentrations before chemistry
         !==============================================================

         ! Initial ADJ_SO4
         ADJ_SO4 = STT_ADJ(I,J,L,IDTSO4)

         ! SO4 drydep frequency [1/s].  Also accounts for the fraction
         ! of each vertical level that is located below the PBL top
         RKT  = DEPSAV(I,J,DRYSO4)  * GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! RKT > 0 denotes that we have SO4 drydep occurring
         IF ( RKT > 0d0 ) THEN
            
            !-----------------------------------------------------------
            ! CASE 1: SO4 production from SO2 and SO4 loss by drydep
            !-----------------------------------------------------------

            ! Fraction of SO4 lost to drydep [unitless]
            RKT   = RKT * DTCHEM
            
            ! Pre-compute exponential term for use below
            E_RKT = EXP( -RKT ) 

            ! Adjoint of SO4 change due to deposition
            ADJ_SO40 =  ADJ_SO4 * E_RKT

            ! Adjoint of SO4 change due to SO2 production
            ADPSO4_SO2(I,J,L) = ADJ_SO4 / RKT
     &                            * ( 1.d0 - E_RKT )

         ELSE

            !-----------------------------------------------------------
            ! CASE 2: Production of SO4 from SO2; no SO4 drydep loss
            !-----------------------------------------------------------

            ! Adjoint of SO4 change due to SO2 production
            ADJ_SO40            = ADJ_SO4
            ADPSO4_SO2(I,J,L)   = ADJ_SO4

         ENDIF

         ! Update global array
         STT_ADJ(I,J,L,IDTSO4) = ADJ_SO40

         
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_SO4_ADJ

!!------------------------------------------------------------------------------
!
!      !SUBROUTINE PHASE_SO4
!      !
!      ! *** Currently under development ***
!      !
!      !END SUBROUTINE PHASE_SO4
!
!!------------------------------------------------------------------------------
!
!      !SUBROUTINE PHASE_RADIATIVE
!      !
!      ! *** Currently under development ***
!      !
!      !END SUBROUTINE PHASE_RADIATIVE
!
!!------------------------------------------------------------------------------

      SUBROUTINE CHEM_MSA_ADJ
!
!******************************************************************************
!  Subroutine ADJ_CHEM_MSA is the adjoint MSA chemistry subroutine for the 
!  fwd routine CHEM_MSA. See original routine for notes.  (dkh, 10/13/05)
!
!  See original routine (rjp, bdf, bmy, 5/31/00, 10/25/05) for notes. 
!
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!                                                                          
!  MSA = MSA_0 * exp(-kt)
!    where k = dry deposition.                                             
!                      
!  NOTES:  
!  (1 ) Updated to GCv8 adjoint (dkh, 09/28/09) 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD,   ONLY : IDTMSA

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN_GCTM"       ! AIRMW
#     include "CMN_DIAG"       ! ND44

      ! Local variables
      INTEGER                 :: I,      J,    L,        PBL_MAX
      TYPE (XPLEX)                  :: DTCHEM, MSA0, MSA,      RK       
      TYPE (XPLEX)               :: RKT,    FLUX, AREA_CM2, F_UNDER_TOP
      TYPE (XPLEX)                  :: ADJ_MSA, ADJ_MSA0

      !=================================================================
      ! CHEM_MSA_ADJ begins here!
      !=================================================================
      IF ( IDTMSA == 0 .or. DRYMSA == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, MSA0, RKT, MSA, AREA_CM2, FLUX )
!$OMP+PRIVATE( ADJ_MSA0, ADJ_MSA )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake
         RKT      = 0d0

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Initial ADJ_MSA
         ADJ_MSA = STT_ADJ(I,J,L,IDTMSA)

         ! Only apply drydep loss to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0 ) THEN
         
            ! MSA drydep frequency [1/s].  Also accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            RKT = DEPSAV(I,J,DRYMSA) * F_UNDER_TOP

            ! RKT > 0 denotes that we have drydep occurring
            IF ( RKT > 0.d0 ) THEN

               ! Fraction of MSA lost to drydep [unitless]
               RKT = RKT * DTCHEM

               ! Adjoint of MSA change due to deposition
               ADJ_MSA0 =  ADJ_MSA * EXP( -RKT )

            ELSE

               ADJ_MSA0            = ADJ_MSA

            ENDIF

            ! Update global array
            STT_ADJ(I,J,L,IDTMSA) = ADJ_MSA0

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_MSA_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH3_ADJ
!
!******************************************************************************
!  Subroutine ADJ_CHEM_NH3 is the adjoint NH3 chemistry subroutine for the 
!  fwd routine CHEM_NH3. (dkh, 10/13/05)
!
!  See originat routine (rjp, bdf, bmy, 1/2/02, 10/25/05)  for notes. 
!                                                                          
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!                                                                          
!  NH3 = NH3_0 * exp(-kt)
!    where k = dry deposition.                                             
!                      
!  NOTES:  
!  (1 ) Updated to GCv8 adjoint (dkh, 09/28/09) 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD,   ONLY : IDTNH3

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER :: I,      J,        L,    PBL_MAX
      TYPE (XPLEX)  :: DTCHEM, NH30,     NH3
      TYPE (XPLEX)  :: FREQ,   AREA_CM2, FLUX, F_UNDER_TOP
      TYPE (XPLEX)  :: ADJ_NH3, ADJ_NH30

      !=================================================================
      ! CHEM_NH3_ADJ begins here!
      !=================================================================
      IF ( IDTNH3 == 0 .or. DRYNH3 == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, FREQ, NH30, NH3, AREA_CM2, FLUX )
!$OMP+PRIVATE( ADJ_NH30, ADJ_NH3 )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN

            ! NH3 drydep frequency [1/s].  Also accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH3) * F_UNDER_TOP

            ! Only compute drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Initial ADJ_NH3
               ADJ_NH3 = STT_ADJ(I,J,L,IDTNH3)

               ADJ_NH30 = ADJ_NH3 * EXP( -FREQ * DTCHEM ) 

               STT_ADJ(I,J,L,IDTNH3) = ADJ_NH30

            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      ! Return to calling program
      END SUBROUTINE CHEM_NH3_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH4_ADJ
!
!******************************************************************************
!  Subroutine CHEM_NH4_ADJ is the adjoint NH4 chemistry subroutine for the 
!  fwd routine CHEM_NH4. (dkh, 10/13/05)
!
!  See original (rjp, bdf, bmy, 1/2/02, 10/25/05)  for notes. 
!                                                                          
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!                                                                          
!  NH4 = NH4_0 * exp(-kt)
!    where k = dry deposition.                                             
!                      
!  NOTES:              
!  (1 ) Updated to GCv8 adjoint 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD,   ONLY : IDTNH4

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER :: I,      J,    L,        PBL_MAX
      TYPE (XPLEX)  :: DTCHEM, NH4,  NH40
      TYPE (XPLEX)  :: FREQ,   FLUX, AREA_CM2, F_UNDER_TOP
      TYPE (XPLEX)  :: ADJ_NH4, ADJ_NH40

      !=================================================================
      ! CHEM_NH4_ADJ begins here!
      !=================================================================
      IF ( IDTNH4 == 0 .or. DRYNH4 == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, FREQ, NH40, NH4, AREA_CM2, FLUX )
!$OMP+PRIVATE( ADJ_NH40, ADJ_NH4 )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )       

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN

            ! NH4 drydep frequency [1/s].  Also accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH4) * F_UNDER_TOP

            ! Only apply drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Initial ADJ_NH4
               ADJ_NH4 = STT_ADJ(I,J,L,IDTNH4)
 
               ! Adjoint of NH4 change due to deposition
               ADJ_NH40 = ADJ_NH4 * EXP( -FREQ * DTCHEM ) 

               ! Update global array
               STT_ADJ(I,J,L,IDTNH4) = ADJ_NH40

            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      ! Return to calling program
      END SUBROUTINE CHEM_NH4_ADJ

!!------------------------------------------------------------------------------
!
!      SUBROUTINE CHEM_NH4aq
!!
!!******************************************************************************
!!  Subroutine CHEM_NH4aq removes NH4aq from the surface via dry deposition.
!!  (cas, bmy, 1/6/05, 10/25/05)  
!!                                                                          
!!  Reaction List:
!!  ============================================================================
!!  (1 ) NH4aq = NH4_0aq * EXP( -dt )  where d = dry deposition rate [s-1]
!!        
!!  NOTES:
!!  (1 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!!        "pbl_mix_mod.f".  Also reference GET_PBL_MAX_L from "pbl_mix_mod.f"
!!        Vertical DO-loops can run up to PBL_MAX and not LLTROP. (bmy, 2/22/05)
!!  (31) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,      ONLY : AD
!      USE DIAG_MOD,     ONLY : AD44
!      USE DRYDEP_MOD,   ONLY : DEPSAV
!      USE GRID_MOD,     ONLY : GET_AREA_CM2
!      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
!      USE TIME_MOD,     ONLY : GET_TS_CHEM
!      USE TRACER_MOD,   ONLY : STT, TCVV, XNUMOL
!      USE TRACERID_MOD, ONLY : IDTNH4aq
!
!#     include "CMN_SIZE"  ! Size parameters
!#     include "CMN_DIAG"  ! ND44
!
!      ! Local variables
!      INTEGER :: I,      J,     L,        PBL_MAX
!      TYPE (XPLEX)  :: DTCHEM, NH4aq, NH4aq0
!      TYPE (XPLEX)  :: FREQ,   FLUX,  AREA_CM2, F_UNDER_TOP
!      TYPE (XPLEX)  :: T44(IIPAR,JJPAR,LLTROP)
!
!      !=================================================================
!      ! CHEM_NH4 begins here!
!      !=================================================================
!      IF ( IDTNH4aq == 0 .or. DRYNH4aq == 0 ) RETURN
!
!      ! DTCHEM is the chemistry interval in seconds
!      DTCHEM  = GET_TS_CHEM() * 60d0 
!
!      ! Model level where maximum PBL height occurs 
!      PBL_MAX = GET_PBL_MAX_L()      
!
!      ! Zero T44 array
!      IF ( ND44 > 0 ) THEN
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLTROP 
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            T44(I,J,L) = 0d0
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,      J,     L,        F_UNDER_TOP, FREQ  )
!!$OMP+PRIVATE( NH4aq0, NH4aq, AREA_CM2, FLUX               )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO L = 1, PBL_MAX
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
!         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )     
!
!         ! Only apply drydep to boxes w/in the PBL
!         IF ( F_UNDER_TOP > 0d0 ) THEN
!
!            ! NH4 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
!            ! of each grid box (I,J,L) that is located beneath the PBL top
!            FREQ = DEPSAV(I,J,DRYNH4aq) * F_UNDER_TOP
!
!            ! Only apply drydep loss if FREQ is nonzero
!            IF ( FREQ > 0d0 ) THEN
!
!               ! Initial NH4 [v/v]
!               NH4aq0 = STT(I,J,L,IDTNH4aq)
!         
!               ! Amount of NH4 lost to drydep [v/v]
!               NH4aq = NH4aq0 * ( 1d0 - EXP( -FREQ * DTCHEM ) )
!
!               ! Prevent underflow condition
!               IF ( NH4aq < SMALLNUM ) NH4aq = 0d0
!
!               ! Subtract NH4 lost to drydep from initial NH4 [v/v]
!               STT(I,J,L,IDTNH4aq) = NH4aq0 - NH4aq
!
!               !========================================================
!               ! ND44 diagnostic: Drydep flux of NH4 [molec/cm2/s]
!               !========================================================
!               IF ( ND44 > 0 .and. NH4aq > 0d0 ) THEN
!         
!                  ! Surface area [cm2]
!                  AREA_CM2 = GET_AREA_CM2( J )
!                  
!                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
!                  FLUX = NH4aq * AD(I,J,L)        / TCVV(IDTNH4aq)
!                  FLUX = FLUX  * XNUMOL(IDTNH4aq) / AREA_CM2 / DTCHEM
!
!                  ! Store dryd flx in ND44_TMP as a placeholder
!                  T44(I,J,L) = T44(I,J,L) + FLUX
!               ENDIF
!            ENDIF
!         ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      !===============================================================
!      ! ND44: Sum drydep fluxes by level into the AD44 array in
!      ! order to ensure that  we get the same results w/ sp or mp 
!      !===============================================================
!      IF ( ND44 > 0 ) THEN 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!         DO L = 1, PBL_MAX
!            AD44(I,J,DRYNH4aq,1) = AD44(I,J,DRYNH4aq,1) + T44(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE CHEM_NH4aq
!
!!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NIT_ADJ
!
!******************************************************************************
!  Subroutine ADJ_CHEM_NIT is the adjoint NIT chemistry subroutine for the 
!  fwd routine CHEM_NIT. See original routine for notes.  (dkh, 10/13/05)
!
!  * Does not deal with NITs *
!
!  See original (rjp, bdf, bmy, 1/2/02, 5/23/06) for notes. 
!                                                                          
!  Reaction List:
!  ============================================================================
!  (1 ) NIT = NIT_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
!  NOTES:
!  (1 ) Updates to GCv8 adjoint (dkh, 09/28/09) 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,    ONLY : LSSALT
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD,   ONLY : IDTNIT, IDTNITs

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER :: I,        J,     L,    N,     N_ND44, PBL_MAX
      TYPE (XPLEX)  :: DTCHEM,   NIT,   NITs, NIT0,  NIT0s,  E_RKT
      TYPE (XPLEX)  :: E_RKTs,   FLUX,  FREQ, FREQs, RKT,    RKTs   
      TYPE (XPLEX)  :: AREA_CM2, F_UNDER_TOP
      TYPE (XPLEX)  :: ADJ_NIT, ADJ_NIT0

      !=================================================================
      ! CHEM_NIT_ADJ begins here!
      !=================================================================

      ! Return if tracers are not defined
      !IF ( IDTNIT == 0 .or. IDTNITs == 0 ) RETURN
      !IF ( DRYNIT == 0 .or. DRYNITs == 0 ) RETURN
      IF ( IDTNIT == 0 ) RETURN
      IF ( DRYNIT == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()      

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()      
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,        NIT0,        NIT0s, NIT   )
!$OMP+PRIVATE( NITs, FREQ,   FREQs,    F_UNDER_TOP, RKT,   E_RKT ) 
!$OMP+PRIVATE( RKTs, E_RKTs, AREA_CM2, FLUX                      )
!$OMP+PRIVATE( ADJ_NIT0, ADJ_NIT )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize variables
         NIT   = 0d0
         NITs  = 0d0
         FREQ  = 0d0
         FREQs = 0d0

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )     

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN 

            !===========================================================
            ! NIT chemistry
            !===========================================================

            ! NIT drydep frequency [1/s].  Also accounts for the fraction
            ! of each vertical level that is located below the PBL top
            FREQ  = DEPSAV(I,J,DRYNIT)  * F_UNDER_TOP

            ! If there is drydep ...
            IF ( FREQ > 0d0 ) THEN

               ! Fraction of NIT lost to drydep [unitless] (bec, 12/15/04)
               RKT  = FREQ  * DTCHEM

               ! Pre-compute the exponential term
               E_RKT = EXP( -RKT )

               ! Initial ADJ_NIT
               ADJ_NIT = STT_ADJ(I,J,L,IDTNIT)

               ! Adjoint of NIT change due to deposition
               ADJ_NIT0 = ADJ_NIT * E_RKT

               ! Update global array
               STT_ADJ(I,J,L,IDTNIT) = ADJ_NIT0


            ENDIF

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      ! Return to calling program
      END SUBROUTINE CHEM_NIT_ADJ

!-----------------------------------------------------------------------------

      SUBROUTINE EMISSSULFATE_ADJ
!
!******************************************************************************
!  Subroutine EMISSSULFATE_ADJ is the interface between the adjoint 
!  of GEOS-CHEM and the adjoint of the sulfate emission routines.
!  (dkh, 04/19/06)  
!
!  Based on forward routine EMISSSULFATE (bmy, 6/7/00, 10/3/05)
! 
!  NOTES:
!  (1 ) Update to GCv8 (dkh, 11/03/09) 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE LOGICAL_MOD,       ONLY : LSHIPSO2,   LPRT,   LBIOMASS !(win,5/1/09)
      USE TIME_MOD,          ONLY : GET_SEASON, GET_MONTH
      USE TIME_MOD,          ONLY : GET_YEAR,   ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : STT,        ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,      ONLY : IDTNITs,    IDTSO4s
      USE TRACERID_MOD,      ONLY : IDTDMS,     IDTSO2 
      USE TRACERID_MOD,      ONLY : IDTSO4,     IDTNH3
      USE GFED2_BIOMASS_MOD, ONLY : GFED2_IS_NEW
      USE LOGICAL_MOD,       ONLY : LANTHRO
  
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: FIRSTEMISS = .TRUE. 
      INTEGER            :: NSEASON, MONTH, YEAR   

      !=================================================================
      ! EMISSSULFATE_ADJ begins here!
      !=================================================================

      ! Get the season and month
      NSEASON = GET_SEASON()
      MONTH   = GET_MONTH()
      YEAR    = GET_YEAR()

      ! fwd code:
      !!=================================================================
      !! Add emissions into the STT tracer array
      !!=================================================================
      !IF ( IDTDMS /= 0 ) CALL SRCDMS( STT(:,:,:,IDTDMS)          )
      !IF ( IDTSO2 /= 0 ) CALL SRCSO2( STT(:,:,:,IDTSO2), NSEASON )
      !IF ( IDTSO4 /= 0 ) CALL SRCSO4( STT(:,:,:,IDTSO4)          )
      !IF ( IDTNH3 /= 0 ) CALL SRCNH3( STT(:,:,:,IDTNH3)          )
      ! adj code:
      IF ( IDTNH3 /= 0 ) CALL SRCNH3_ADJ( STT_ADJ(:,:,:,IDTNH3)   )
      ! add these later as needed 
      !IF ( IDTSO4 /= 0 ) CALL SRCSO4_ADJ( STT_ADJ(:,:,:,IDTSO4)   )
      IF ( IDTSO2 /= 0 ) CALL SRCSO2_ADJ( STT_ADJ(:,:,:,IDTSO2),NSEASON)
      !IF ( IDTDMS /= 0 ) CALL SRCDMS_ADJ( STT_ADJ(:,:,:,IDTDMS)   )

      ! If we wanted to pass sensitivities all the way back to specific
      ! emissions inventories, we could carry on like so...
!      IF(  ( GFED2_IS_NEW() .or. ITS_A_NEW_MONTH_ADJ() ) .AND. 
!     &     ( LBIOMASS )  ) THEN 
!
!         ! fwd code:
!         !CALL GET_BIOMASS_NH3
!         ! adj code:
!         print*, ' need to add GET_BIOMASS_NH3_ADJ '
!         !CALL GET_BIOMASS_NH3_ADJ
!         !IF ( LPRT ) CALL DEBUG_MSG('### EMISSSULFATE: GET_BM_NH3_ADJ')
!
!         ! fwd code:
!         !CALL GET_BIOMASS_SO2
!         ! adj code:
!         print*, ' need to add GET_BIOMASS_SO2_ADJ'
!         !CALL GET_BIOMASS_SO2_ADJ
!         !IF ( LPRT ) CALL DEBUG_MSG('### EMISSSULFATE: GET_BM_SO2_ADJ')
!         
!      ENDIF
!
!      !=================================================================
!      ! If this is a new month, read in the monthly mean quantities
!      !=================================================================
!      IF ( ITS_A_NEW_MONTH() ) THEN 
!
!         ! fwd code:
!         !! Read oxidants for the offline simulation only
!         !IF ( ITS_AN_AEROSOL_SIM() ) CALL READ_OXIDANT( MONTH )
!         ! adj code: not supported 
!
!         ! fwd code: 
!         !! Also read ship exhaust SO2 if necessary 
!         !CALL READ_SHIP_SO2( MONTH )
!         ! adj code: 
!         CALL READ_SHIP_SO2_ADJ( MONTH )
!
!         
!         ! Add LANTHRO switch to turn off anthropogenic emissions.
!         ! (ccc, 4/15/09)
!         IF ( LANTHRO ) THEN
! 
!            ! fwd code: 
!            !CALL READ_AIRCRAFT_SO2( MONTH )
!            !CALL READ_ANTHRO_SOx( MONTH, NSEASON )
!            !CALL READ_ANTHRO_NH3( MONTH )
!            ! adj code: 
!            CALL READ_ANTHRO_NH3_ADJ( MONTH )
!            CALL READ_ANTHRO_SOx_ADJ( MONTH, NSEASON )
!            CALL READ_AIRCRAFT_SO2_ADJ( MONTH )
!
!         ENDIF
!
!         ! fwd code: 
!         !! Read monthly mean data
!         !CALL READ_SST( MONTH, YEAR )
!         !CALL READ_OCEAN_DMS( MONTH )
!         !CALL READ_BIOFUEL_SO2( MONTH )
!         !CALL READ_BIOFUEL_NH3( MONTH )
!         !CALL READ_NATURAL_NH3( MONTH )
!         ! adj code: 
!         CALL READ_NATURAL_NH3_ADJ( MONTH )
!         CALL READ_BIOFUEL_NH3_ADJ( MONTH )
!         CALL READ_BIOFUEL_SO2_ADJ( MONTH )
!         CALL READ_OCEAN_DMS_ADJ( MONTH )
!         CALL READ_SST_ADJ( MONTH, YEAR )
!
!
!      ENDIF


      ! Return to calling program
      END SUBROUTINE EMISSSULFATE_ADJ

!-----------------------------------------------------------------------------

      SUBROUTINE SRCNH3_ADJ( TC_ADJ )
!
!******************************************************************************
!  Subroutine ADJ_SRCNH3 computes adjoint of emission from adjoint of tracer
!  (dkh, 05/05/06) 
! 
!  Based on forward routine SRCNH3 (rjp, bmy, 12/17/01, 2/22/05)
!
!  Arguments as Input
!  ============================================================================
!  (1 ) TC_ADJ (TYPE (XPLEX) ) : Array for adjoint of NH3 tracer mass in kg
!
!  NOTES:
!  (1 ) Update to GCv8 
!  (2 ) Now adjoints are w.r.t. NH3an, ENH3_bb, NH3bf and ENH3_na
!  (3 ) Add LEMS_ABS option (dkh, 02/17/11) 
!  (4 ) Add support for LNEI05 (dkh, 02/19/11) 
!******************************************************************************
!
      ! References to F90 modules
      USE CAC_ANTHRO_MOD, ONLY : GET_CANADA_MASK
      USE CAC_ANTHRO_MOD, ONLY : GET_CAC_ANTHRO
      USE DIAG_MOD,       ONLY : AD13_NH3_an, AD13_NH3_bb
      USE DIAG_MOD,       ONLY : AD13_NH3_bf, AD13_NH3_na
      USE DAO_MOD,        ONLY : PBL
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE EPA_NEI_MOD,    ONLY : GET_EPA_ANTHRO, GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,    ONLY : GET_USA_MASK
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE LOGICAL_MOD,    ONLY : LNEI99, LCAC, LNEI05
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,       ONLY : GET_DAY_OF_WEEK, GET_TS_EMIS
      USE TRACER_MOD,     ONLY : XNUMOL
      USE TRACERID_MOD,   ONLY : IDTNH3

      USE ADJ_ARRAYS_MOD, ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ENH3_an
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ENH3_na
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ENH3_bb
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ENH3_bf
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD
      USE ADJ_ARRAYS_MOD, ONLY : EMS_ADJ
      USE LOGICAL_ADJ_MOD,ONLY : LPRINTFD
      USE LOGICAL_ADJ_MOD,ONLY : LEMS_ABS
      USE SULFATE_MOD,    ONLY : ENH3_an
      USE SULFATE_MOD,    ONLY : ENH3_bb
      USE SULFATE_MOD,    ONLY : ENH3_bf
      USE SULFATE_MOD,    ONLY : ENH3_na

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND13
#     include "CMN_GCTM"     ! SCALE_HEIGHT
      
      ! Argumetns
      TYPE (XPLEX),  INTENT(IN)    :: TC_ADJ(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL                :: WEEKDAY
      INTEGER                :: I, J, L,  K, NTOP, DAY_NUM
      TYPE (XPLEX)                 :: FEMIS,    DTSRCE
      TYPE (XPLEX)                 :: AREA_CM2, EPA_AN,  EPA_BF
      TYPE (XPLEX)                 :: CAC_AN
      TYPE (XPLEX)                 :: NH3an(IIPAR,JJPAR)
      TYPE (XPLEX)                 :: NH3bf(IIPAR,JJPAR)
      TYPE (XPLEX)                 :: TNH3_ADJ

      !=================================================================
      ! SRCNH3_ADJ begins here!
      !=================================================================


      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0
 
      ! Get current day of the week
      DAY_NUM = GET_DAY_OF_WEEK()

      ! Is it a weekday?
      WEEKDAY = ( DAY_NUM > 0 .and. DAY_NUM < 6 )

      !=================================================================
      ! Overwrite USA with EPA/NEI NH3 emissions (if necessary)
      ! Store emissions into local arrays NH3an, NH3bf
      !=================================================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, AREA_CM2, EPA_AN, EPA_BF, CAC_AN )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, EPA_AN, EPA_BF, CAC_AN, L )
      DO J = 1, JJPAR
         
         !-------------------------------------------------------------------
         ! NOTE: There seems to be some problems with the EPA/NEI NH3 
         ! emissions.  Therefore we will use the existing emissions for NH3 
         ! until further notice.  Comment out the lines below until 
         ! further notice.  (bmy, 11/17/04) 
         !! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )
         !-------------------------------------------------------------------

         DO I = 1, IIPAR

            !-----------------------------------------------------------------
            ! NOTE: There seems to be some problems with the EPA/NEI NH3 
            ! emissions.  Therefore we will use the existing emissions for 
            ! NH3 until further notice.  Comment out the lines below until 
            ! further notice.  (bmy, 11/17/04) 
            !
            !! If we are using EPA/NEI99 emissions ...
            !IF ( LNEI99 ) THEN
            !
            !   ! If we are over the USA ...
            !   IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
            !
            !      ! Read NH3 anthro emissions in [molec NH3/cm2/s]
            !      EPA_AN     = GET_EPA_ANTHRO(  I, J, IDTNH3, WEEKDAY )
            !      EPA_BF     = GET_EPA_BIOFUEL( I, J, IDTNH3, WEEKDAY )
            !
            !      ! Convert from [molec NH3/cm2/s] to [kg NH3/box/sec]
            !      NH3an(I,J) = EPA_AN * AREA_CM2 / XNUMOL(IDTNH3)
            !      NH3bf(I,J) = EPA_BF * AREA_CM2 / XNUMOL(IDTNH3)
            !
            !   ELSE
            !
            !      ! If we are not over the USA, just use the regular 
            !      ! emissions in NH3_an and NH3bf (bmy, 11/16/04)
            !      NH3an(I,J) = ENH3_an(I,J)
            !      NH3bf(I,J) = ENH3_bf(I,J)
            !
            !   ENDIF
            !
            !ELSE
            !-----------------------------------------------------------------

               ! If we are not using the EPA/NEI emissions, just copy the 
               ! regular ENH3_an and ENH3_bf to local arrays. (bmy, 11/16/04)
               NH3an(I,J) = ENH3_an(I,J)
               NH3bf(I,J) = ENH3_bf(I,J)
            
            !-----------------------------------------------------------------
            ! NOTE: There seems to be some problems with the EPA/NEI NH3 
            ! emissions.  Therefore we will use the existing emissions for 
            ! NH3 until further notice.  Comment out the lines below until 
            ! further notice.  (bmy, 11/17/04) 
            !ENDIF
            !-----------------------------------------------------------------

            ! If we are using CAC emissions and over CANADA
            IF ( LCAC ) THEN
            IF ( GET_CANADA_MASK( I, J ) > 0d0 ) THEN

               ! Read NH3 anthro emissions in [molec NH3/cm2/s]
               CAC_AN = GET_CAC_ANTHRO( I, J, IDTNH3, 
     &                                  MOLEC_CM2_S=.TRUE.)

               ! Convert from [molec NH3/cm2/c] to [kg NH3/box/sec]
               NH3an(I,J) = CAC_AN * AREA_CM2 / XNUMOL(IDTNH3)

            ENDIF
            ENDIF

            ! If we are using NEI 2005 over North America
            IF ( LNEI05 ) THEN
            IF ( NEI05_MASK( I, J ) > 0d0 ) THEN

               ! For 2D NH3an (dkh, 02/19/11) 
               NH3an(I,J) = 0d0 

               DO L = 1, NOXLEVELS

                  ! Read NH3 anthro emissions in [molec NH3/cm2/s]                  
                  EPA_AN = GET_NEI2005_ANTHRO( I, J, L, IDTNH3,
     &                                   WEEKDAY, MOLEC_CM2_S=.TRUE.)

                  ! Convert from [molec NH3/cm2/c] to [kg NH3/box/sec]
                  ! Keep NH3an 2D for now. (dkh, 02/19/11) 
                  !NH3an(I,J,L) = EPA_AN * AREA_CM2 / XNUMOL(IDTNH3)
                  NH3an(I,J) = NH3an(I,J) 
     &                       + EPA_AN * AREA_CM2 / XNUMOL(IDTNH3)

               ENDDO

            ENDIF
            ENDIF

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Partition NH3 emissions into the STT tracer array
      !=================================================================

      ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, FEMIS )
!$OMP+PRIVATE( TNH3_ADJ )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
          
         ! Layer where the PBL top happens
         NTOP   = CEILING( (GET_PBL_TOP_L( I, J )) )

         !==============================================================
         ! Add NH3 emissions [kg NH3/box] into the tracer array
         ! Partition total NH3 throughout the entire boundary layer
         !==============================================================

         TNH3_ADJ = 0d0 

         ! Loop over all levels in the boundary layer
         DO L = 1, NTOP

            ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
            FEMIS     = GET_FRAC_OF_PBL( I, J, L )

            ! fwd code: 
            !! Add NH3 emissions into tracer array [kg NH3/timestep]
            !TC(I,J,L) = TC(I,J,L) + ( TNH3 * FEMIS * DTSRCE )
            ! adj code:
            TNH3_ADJ = TNH3_ADJ + ( TC_ADJ(I,J,L) * FEMIS * DTSRCE )

         ENDDO

         ! fwd code:
         !! Sum all types of NH3 emission [kg/box/s]
         !TNH3   = NH3an(I,J) + ENH3_bb(I,J) + 
     &   !         NH3bf(I,J) + ENH3_na(I,J)
         ! adj code:
         EMS_SF_ADJ(I,J,1,IDADJ_ENH3_an) = 
     &      EMS_SF_ADJ(I,J,1,IDADJ_ENH3_an) + TNH3_ADJ * NH3an(I,J)
         EMS_SF_ADJ(I,J,1,IDADJ_ENH3_bb) = 
     &     EMS_SF_ADJ(I,J,1,IDADJ_ENH3_bb)  + TNH3_ADJ * ENH3_bb(I,J)
         EMS_SF_ADJ(I,J,1,IDADJ_ENH3_bf) = 
     &     EMS_SF_ADJ(I,J,1,IDADJ_ENH3_bf)  + TNH3_ADJ * NH3bf(I,J)
         EMS_SF_ADJ(I,J,1,IDADJ_ENH3_na) = 
     &     EMS_SF_ADJ(I,J,1,IDADJ_ENH3_na)  + TNH3_ADJ * ENH3_na(I,J)
        
         ! dkh debug
         IF ( I == IFD .and. J == JFD .AND. LPRINTFD ) THEN
            print*, ' SRCNH3 adj : NH3an  = ', NH3an(I,J)
            print*, ' SRCNH3 adj : ENH3_bb= ', ENH3_bb(I,J)
            print*, ' SRCNH3 adj : NH3bf  = ', NH3bf(I,J)
            print*, ' SRCNH3 adj : ENH3_na= ', ENH3_na(I,J)
            print*, ' SRCNH3 adj : TNH3_ADJ = ', TNH3_ADJ
            print*, ' EMS_SF_ADJ(ENH3_an) = ', 
     &         EMS_SF_ADJ(I,J,1,IDADJ_ENH3_an)  
         ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Optional diagnostic -- add sensitivity w.r.t to absolute emissions (dkh, 02/17/11) 
      IF ( LEMS_ABS ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, FEMIS )
!$OMP+PRIVATE( TNH3_ADJ )
!$OMP+SCHEDULE( DYNAMIC )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
          
            ! Layer where the PBL top happens
            NTOP   = CEILING( (GET_PBL_TOP_L( I, J )) )

            !==============================================================
            ! Add NH3 emissions [kg NH3/box] into the tracer array
            ! Partition total NH3 throughout the entire boundary layer
            !==============================================================
   
            TNH3_ADJ = 0d0 
   
            ! Loop over all levels in the boundary layer
            DO L = 1, NTOP
   
               ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
               FEMIS     = GET_FRAC_OF_PBL( I, J, L )
   
               ! fwd code: 
               !! Add NH3 emissions into tracer array [kg NH3/timestep]
               !TC(I,J,L) = TC(I,J,L) + ( TNH3 * FEMIS * DTSRCE )
               ! adj code:
               TNH3_ADJ = TNH3_ADJ + ( TC_ADJ(I,J,L) * FEMIS * DTSRCE )

            ENDDO

            ! fwd code:
            !! Sum all types of NH3 emission [kg/box/s]
            !TNH3   = NH3an(I,J) + ENH3_bb(I,J) + 
     &      !         NH3bf(I,J) + ENH3_na(I,J)
            ! adj code: also convert to J / (kg/box/timestep)
            EMS_ADJ(I,J,1,IDADJ_ENH3_an) 
     &         = EMS_ADJ(I,J,1,IDADJ_ENH3_an) + TNH3_ADJ / DTSRCE
            EMS_ADJ(I,J,1,IDADJ_ENH3_bb)
     &         = EMS_ADJ(I,J,1,IDADJ_ENH3_bb) + TNH3_ADJ / DTSRCE
            EMS_ADJ(I,J,1,IDADJ_ENH3_bf)
     &         = EMS_ADJ(I,J,1,IDADJ_ENH3_bf) + TNH3_ADJ / DTSRCE
            EMS_ADJ(I,J,1,IDADJ_ENH3_na)
     &         = EMS_ADJ(I,J,1,IDADJ_ENH3_na) + TNH3_ADJ / DTSRCE
        

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF 

      ! Return to calling program
      END SUBROUTINE SRCNH3_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE SRCSO2_ADJ( TC_ADJ, NSEASON )
!
!******************************************************************************
!  Subroutine SRCSO2_ADJ computes the adjoint of SO2 emisson scaling factors. 
!  (dkh, 05/04/06, 02/04/10) 
!
!  Based on forward model code SRCSO2  (rjp, bdf, bmy, 6/2/00, 2/27/09)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  (1 ) NSEASON (INTEGER) : Season number: 1=DJF; 2=MAM; 3=JJA; 4=SON
!  (2 ) TC_ADJ  (TYPE (XPLEX) ) : adoint of SO2 tracer mass
! 
!  NOTES:
!  (1 ) Doesn't deal with aircraft or volcanoe source yet. 
!  (2 ) Updated to GCv8 (dkh, 02/04/10) 
!  (3 ) Add support for LEMS_ABS 
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE BRAVO_MOD,      ONLY : GET_BRAVO_ANTHRO, GET_BRAVO_MASK
      USE CAC_ANTHRO_MOD, ONLY : GET_CANADA_MASK,  GET_CAC_ANTHRO
      USE DIAG_MOD,       ONLY : AD13_SO2_an,      AD13_SO2_ac
      USE DIAG_MOD,       ONLY : AD13_SO2_bb,      AD13_SO2_nv
      USE DIAG_MOD,       ONLY : AD13_SO2_ev,      AD13_SO2_bf
      USE DIAG_MOD,       ONLY : AD13_SO2_sh
      USE DAO_MOD,        ONLY : BXHEIGHT, PBL
      USE EPA_NEI_MOD,    ONLY : GET_EPA_ANTHRO,   GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,    ONLY : GET_USA_MASK
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,    ONLY : LBRAVO, LNEI99,   LSHIPSO2
      USE LOGICAL_MOD,    ONLY : LCAC,   LNEI05
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_OF_PBL,  GET_PBL_TOP_L
      USE TIME_MOD,       ONLY : GET_TS_EMIS,      GET_DAY_OF_YEAR 
      USE TIME_MOD,       ONLY : GET_DAY_OF_WEEK
      USE TRACER_MOD,     ONLY : XNUMOL
      USE TRACERID_MOD,   ONLY : IDTSO2

      ! adj_group
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ESO2_an1,   IDADJ_ESO2_an2
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ESO2_bb,    IDADJ_ESO2_bf
      USE ADJ_ARRAYS_MOD, ONLY : IDADJ_ESO2_sh
      USE ADJ_ARRAYS_MOD, ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD, ONLY : EMS_ADJ
      USE LOGICAL_ADJ_MOD,ONLY : LEMS_ABS
      USE SULFATE_MOD,    ONLY : ESO2_an
      USE SULFATE_MOD,    ONLY : ESO2_bb
      USE SULFATE_MOD,    ONLY : ESO2_bf
      USE SULFATE_MOD,    ONLY : ESO2_sh

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND13, LD13 (for now)
#     include "CMN_GCTM"     ! SCALE_HEIGHT

      ! Arguments
      INTEGER, INTENT(IN)    :: NSEASON
      TYPE (XPLEX),  INTENT(INOUT) :: TC_ADJ(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL                :: WEEKDAY
      INTEGER                :: I, J, K, L, LV1, LV2, NTOP, JDAY
      INTEGER                :: DAY_NUM
      TYPE (XPLEX)                 :: ZH(0:LLPAR), DZ(LLPAR), SO2(LLPAR)
      TYPE (XPLEX)                 :: DTSRCE,      HGHT,      SO2SRC
      TYPE (XPLEX)                 :: SLAB,        SLAB1
      TYPE (XPLEX)                 :: TSO2,        FEMIS
      TYPE (XPLEX)                 :: AREA_CM2,    AN,        BF
      TYPE (XPLEX)                 :: SO2an(IIPAR,JJPAR,2)
      TYPE (XPLEX)                 :: SO2bf(IIPAR,JJPAR)

      ! for adjoint 
      TYPE (XPLEX)                 :: SO2SRC_ADJ
      TYPE (XPLEX)                 :: SO2_ADJ(LLPAR)
      TYPE (XPLEX)                 :: TSO2_ADJ

      ! Ratio of molecular weights: S/SO2
      TYPE (XPLEX),  PARAMETER     :: S_SO2 = xplex(32d0 / 64d0,0d0)

      !=================================================================
      ! SRCSO2_ADJ begins here!
      !================================================================

      ! DTSRCE is the emission timestep in seconds
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! JDAY is the day of year (0-365 or 0-366)
      JDAY    = GET_DAY_OF_YEAR()

      ! Get current day of the week
      DAY_NUM = GET_DAY_OF_WEEK()

      ! Is it a weekday?
      WEEKDAY = ( DAY_NUM > 0 .and. DAY_NUM < 6 )

      !=================================================================
      ! First we recalculate fwd model emissions 
      !=================================================================
!  Adjoint of LENV not supported yet 
!      !=================================================================
!      ! SO2 emissions from non-eruptive volcanoes [kg SO2/box/s].
!      ! Assume that emission only occurs at the crater altitude.
!      !=================================================================
!      IF ( LENV ) THEN
!
!         ! Initialize
!         ESO2_nv = 0.d0
!
!         ! Loop thru each non-erupting volcano
!         DO K = 1, NNVOL
!
!            ! Elevation of volcano crater
!            HGHT  = DBLE( IELVn(k) )
!
!            ! Altitude of crater from the ground
!            ZH(0) = 0.d0
!
!            ! Loop over levels
!            DO L = 1, LLPAR
!
!               ! Thickness of layer [m] w/ crater
!               DZ(L) = BXHEIGHT(INV(K),JNV(K),L)
!
!               ! Increment altitude
!               ZH(L) = ZH(L-1) + DZ(L)
!
!               ! If we are at the crater altitude, add emissions and exit
!               IF ( ZH(L-1) <= HGHT .AND. ZH(L) > HGHT ) THEN 
!                  ESO2_nv(INV(K),JNV(K),L) = 
!     &                 ESO2_nv(INV(K),JNV(K),L) + Env(K)
!                  EXIT
!               ENDIF
!            ENDDO
!         ENDDO
!      ENDIF  
!
!      !=================================================================
!      ! Calculate eruptive volcanic emission of SO2.
!      !=================================================================
!      IF ( LEEV ) THEN
!
!         ! Initialize
!         ESO2_ev = 0.D0
!
!         ! Loop thru each erupting volcano
!         DO K = 1, NEVOL
!
!            ! Test to see if the volcano is erupting
!            IF ( JDAY < IDAYS(K) .OR. JDAY > IDAYe(K) ) GOTO 20
!
!            !===========================================================
!            ! Define a slab at the top 1/3 of the volcano plume.
!            !===========================================================
!            HGHT  = DBLE( IHGHT(K) )
!
!            ! slab bottom height
!            SLAB1 = HGHT - ( HGHT - DBLE ( IELVe(K) ) ) / 3.d0 
!
!            ! Slab thickness
!            SLAB  = HGHT - SLAB1 
!            ZH(0) = 0.d0 
!        
!            ! Loop over each level
!            DO L = 1, LLPAR
!
!               ! DZ is the thickness of level L [m]
!               DZ(L) = BXHEIGHT(IEV(K),JEV(K),L)
!
!               ! ZH is the height of the top edge of 
!               ! level L, measured from the ground up [m]
!               ZH(L) = ZH(L-1) + DZ(L) 
!
!               ! max model erup.height
!               IF ( L == LLPAR .AND. HGHT > ZH(L) ) THEN 
!                  LV2 = LLPAR
!                  !HGHT = ZH(L)
!                  !SLAB1 = SLAB1 - ( HGHT - ZH(L) )
!               ENDIF
!
!               !========================================================
!               ! If Zh(l) <= bottom of the slab, go to next level.
!               !========================================================
!               IF ( ZH(L) <= SLAB1 ) GOTO 22
!
!               !========================================================
!               ! If the slab is only in current level: CASE 1
!               !       ---------------------------------- ZH(L)
!               ! HGHT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               ! SLAB1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!               !       ---------------------------------- ZH(L-1)
!               !========================================================
!               IF ( ZH(L-1) <= SLAB1 .AND. ZH(L) >= HGHT ) THEN
!                  LV1   = L
!                  LV2   = L
!                  DZ(L) = SLAB
!
!               !========================================================
!               ! The slab extends more then one level.  Find the 
!               ! lowest (lv1) and the highest (lv2) levels:
!               !       --------------------------------- ZH(L)
!               ! HGHT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               !       --------------------------------- ZH(L-1)
!               ! 
!               !       --------------------------------- ZH(L)
!               ! SLAB1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               !       --------------------------------- ZH(L-1)
!               !========================================================
!               ELSE IF (ZH(L-1) <= SLAB1 .AND. ZH(L) > SLAB1)  THEN
!                  LV1   = L
!                  DZ(L) = ZH(L) - SLAB1
!                  
!               ELSE IF (ZH(L-1) < HGHT  .AND. ZH(L) > HGHT )  THEN
!                  LV2   = L
!                  DZ(L) = HGHT - ZH(L-1)
!                  EXIT          ! do 20 
!        
!               ENDIF
!               
!               ! Go to next level
! 22            CONTINUE         
!            ENDDO
!
!            !===========================================================
!            ! Calculate SO2 emission in the levels between LV1 and LV2.  
!            ! Convert Eev from [kg SO2/box/event] to [kg SO2/box/s].  
!            ! ESO2_ev is distributed evenly with altitude among the slab.
!            !===========================================================
!            DO L = LV1, LV2
!               ESO2_ev(IEV(K),JEV(K),L) = ESO2_ev(IEV(K),JEV(K),L) + 
!     &              EEV(K) / ( (IDAYe(K)-IDAYs(K)+1) * 24.d0 * 3600.d0 )
!     &              * DZ(L) / SLAB
!            ENDDO
!
!            ! Go to next volcano
! 20         CONTINUE 
!         ENDDO
!
!      ENDIF  

      !=================================================================
      ! Overwrite USA    w/ EPA/NEI99 (anthro+biofuel) SO2 emissions 
      ! Overwrite MEXICO w/ BRAVO     (anthro only   ) SO2 emissions
      ! Overwrite CANADA w/ CAC       (anthro only   ) SO2 emissions
      !-----------------------------------------------------------------
      ! Note that we:
      ! Overwrite ASIA   w/ STREETS and EUROPE w/ EMEP
      !  in READ_ANTHRO_SOx.
      !
      ! In both cases, SO4 is a fraction of provided SO2 (except for
      ! EPA). It is done in READ_ANTHRO_SOX in the 1st case, and in
      ! SRCSO4 for inventories dealt with here. EPA is the only one to
      ! provide direct SO4, which is why we deal with it here, even
      ! though it does not have to be like that. Historical.
      !     
      ! So, since we have EPA here, we have to deal with BRAVO and
      ! CAC here to deal with their overlaping mask.
      !     
      ! (amv, phs, 3/12/08, 8/24/09)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, AREA_CM2, AN, BF )
!$OMP+PRIVATE( I, J, L, AREA_CM2, AN, BF )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Initialize
         AN       = 0d0
         BF       = 0d0

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2(J)

         DO I = 1, IIPAR
!------- Prior to 3/8/07 (phs)         
            !!-----------------------------------------------------------
            !! If we are using EPA/NEI99 (anthro + biofuel) ...
            !!-----------------------------------------------------------
            !IF ( LNEI99 ) THEN
            !
            !   ! If we are over the USA ...
            !   IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
            !
            !      ! Read EPA/NEI99 SO2 emissions in [molec/cm2/s]
            !      AN           = GET_EPA_ANTHRO( I, J, IDTSO2, WEEKDAY )
            !      BF           = GET_EPA_BIOFUEL(I, J, IDTSO2, WEEKDAY )
            !   
            !      ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
            !      ! Place all anthro SO2 into surface layer
            !      SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
            !      SO2an(I,J,2) = 0d0
            !   
            !      ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
            !      SO2bf(I,J)   = BF * AREA_CM2 / XNUMOL(IDTSO2)
            !
            !   ELSE
            !
            !      ! If we are not over the USA, then just use the regular 
            !      ! emissions from ESO2_an and ESO2_bf (bmy, 11/16/04)
            !      SO2an(I,J,1) = ESO2_an(I,J,1)
            !      SO2an(I,J,2) = ESO2_an(I,J,2)
            !      SO2bf(I,J)   = ESO2_bf(I,J)
            !
            !   ENDIF
            !
            !ELSE
            ! 
            !   ! If we are not using EPA/NEI99 emissions, then just copy 
            !   ! ESO2_an and ESO2_bf into local arrays (bmy, 11/16/04)
            !   SO2an(I,J,1) = ESO2_an(I,J,1)
            !   SO2an(I,J,2) = ESO2_an(I,J,2)
            !   SO2bf(I,J)   = ESO2_bf(I,J)
            !
            !ENDIF
            !
            !!-----------------------------------------------------------
            !! If we are using BRAVO emissions over Mexico ...
            !!-----------------------------------------------------------
            !IF ( LBRAVO ) THEN
            !
            !   ! If we are over Mexico ...
            !   IF ( GET_BRAVO_MASK( I, J ) > 0d0 ) THEN
            !
            !      ! Read BRAVO SO2 emissions in [molec/cm2/s]
            !      AN           = GET_BRAVO_ANTHRO( I, J, IDTSO2 )
            !   
            !      ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
            !      ! Place all anthro SO2 into surface layer
            !      SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
            !      SO2an(I,J,2) = 0d0
            !
            !   ELSE
            !
            !      ! If we are not over MEXICO, then just use 
            !      ! the regular emissions from ESO2_an
            !      SO2an(I,J,1) = ESO2_an(I,J,1)
            !      SO2an(I,J,2) = ESO2_an(I,J,2)
            !
            !   ENDIF
            !
            !ELSE
            ! 
            !   ! If we are not using BRAVO emissions, then just copy 
            !   ! ESO2_an and ESO2_bf into local arrays
            !   SO2an(I,J,1) = ESO2_an(I,J,1)
            !   SO2an(I,J,2) = ESO2_an(I,J,2)
            !   SO2bf(I,J)   = ESO2_bf(I,J)
            !
            !ENDIF

            !-----------------------------------------------------------
            ! Default SO2 from GEIA or EDGAR (w/ optional STREETS for 
            ! ASIA, and EMEP for Europe)
            !-----------------------------------------------------------
            SO2an(I,J,1) = ESO2_an(I,J,1)
            SO2an(I,J,2) = ESO2_an(I,J,2)
            SO2bf(I,J)   = ESO2_bf(I,J)

            !-----------------------------------------------------------
            ! If we are using EPA/NEI99 over the USA (anthro + biofuel)
            !-----------------------------------------------------------
            IF ( LNEI99 ) THEN
            IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
            
               ! Read EPA/NEI99 SO2 emissions in [molec/cm2/s]
               AN           = GET_EPA_ANTHRO( I, J, IDTSO2, WEEKDAY )
               BF           = GET_EPA_BIOFUEL(I, J, IDTSO2, WEEKDAY )
               
               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer
               SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
               SO2an(I,J,2) = 0d0
               
               ! Convert biofuel SO2 from [molec/cm2/s] to [kg/box/s] 
               SO2bf(I,J)   = BF * AREA_CM2 / XNUMOL(IDTSO2)

            ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are using BRAVO emissions over Mexico (anthro)
            !-----------------------------------------------------------
            IF ( LBRAVO )  THEN 
            IF ( GET_BRAVO_MASK( I, J ) > 0d0 )  THEN
            
               ! Read BRAVO SO2 emissions in [molec/cm2/s]
               AN           = GET_BRAVO_ANTHRO( I, J, IDTSO2 )
               
               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer.
               ! Add to USA emissions if on the border. 
               IF ( LNEI99 .and. GET_USA_MASK( I, J) > 0d0 ) THEN

                  SO2an(I,J,1) = SO2an(I,J,1) + 
     &                           AN * AREA_CM2 / XNUMOL(IDTSO2)
               ELSE
                  SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
               ENDIF

               SO2an(I,J,2) = 0d0

            ENDIF
            ENDIF


            !-----------------------------------------------------------
            ! If we are using CAC emissions over Canada ...
            !-----------------------------------------------------------
            IF ( LCAC ) THEN
            IF ( GET_CANADA_MASK( I, J ) > 0d0 ) THEN

               ! Read CAC SO2 emissions in [molec/cm2/s]
               AN           = GET_CAC_ANTHRO( I, J, IDTSO2,
     &                            MOLEC_CM2_s=.TRUE. )

               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer.
               ! Add to USA emissions if on the border. 
               IF ( LNEI99 .and. GET_USA_MASK( I, J) > 0d0 ) THEN

                  SO2an(I,J,1) = SO2an(I,J,1) + 
     &                           AN * AREA_CM2 / XNUMOL(IDTSO2)
               ELSE
                  SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
               ENDIF

               SO2an(I,J,2) = 0d0

            ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are using EPA/NEI 2005 over USA.
            ! Must be called after CAC and BRAVO to simply overwrite
            ! where they overlap
            ! Modify for use with 2-level SO2an (dkh, 02/19/11) 
            !-----------------------------------------------------------  
            IF ( LNEI05 ) THEN            
            IF ( NEI05_MASK( I, J ) > 0d0 ) THEN
            
               ! Read USA SO2 emissions in [molec/cm2/s]
               ! Level L=1 (dkh, 02/19/11) 
               AN = GET_NEI2005_ANTHRO( I, J, 1, IDTSO2, WEEKDAY,
     &                         MOLEC_CM2_s=.TRUE. )

               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s]
               SO2an(I,J,1) =  AN * AREA_CM2 / XNUMOL(IDTSO2)
               SO2an(I,J,2) =  0d0

               ! Read USA SO2 emissions in [molec/cm2/s]
               DO L = 2, NOXLEVELS
                  AN = GET_NEI2005_ANTHRO( I, J, L, IDTSO2, WEEKDAY,
     &                            MOLEC_CM2_s=.TRUE. )

                  ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s]
                  !SO2an(I,J,L) =  AN * AREA_CM2 / XNUMOL(IDTSO2)
                  SO2an(I,J,2) = SO2an(I,J,2) 
     &                         + AN * AREA_CM2 / XNUMOL(IDTSO2)
               ENDDO
                   

            ENDIF
            ENDIF


         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Adjoint of Add SO2 emissions into model levels
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, SO2, TSO2, FEMIS, SO2SRC )
!$OMP+PRIVATE( SO2SRC_ADJ, SO2_ADJ,      TSO2_ADJ      )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Top of the boundary layer
         NTOP = CEILING( (GET_PBL_TOP_L( I, J )) ) 

         DO L = 1, LLPAR
            
            ! Add SO2 to TC array [kg/box/timestep]
            ! fwd code:
            !TC(I,J,L) = TC(I,J,L) + ( SO2SRC * DTSRCE )
            ! adj code:
            SO2SRC_ADJ = TC_ADJ(I,J,L) * DTSRCE 

            ! SO2 emissions [kg/box/s]
            ! fwd code:
            !SO2SRC = SO2(L)         + ESO2_ac(I,J,L) + 
            !         ESO2_nv(I,J,L) + ESO2_ev(I,J,L) 
            ! adj code:
            SO2_ADJ(L) = SO2SRC_ADJ 
            ! adj code not implemented yet
            ! ESO2_ac_ADJ(I,J,L) = ... etc

         ENDDO

         !===============================================================
         ! Partition the total anthro and biomass SO2 emissions thru
         ! the entire boundary layer (if PBL top is higher than level 2)
         !===============================================================
         IF ( NTOP > 2 ) THEN

            TSO2_ADJ = 0d0

            ! Loop thru levels in the PBL
            DO L  = 1, NTOP
                 
               ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
               FEMIS  = GET_FRAC_OF_PBL( I, J, L )
                 
               ! Partition total SO2 into level K
               ! fwd_code:
               !SO2(L) = FEMIS * TSO2
               ! adj_code:
               TSO2_ADJ = TSO2_ADJ + SO2_ADJ(L) * FEMIS


            ENDDO

            ! Also add SO2 from ship exhaust if necessary (bec, bmy, 5/20/04)
            ! fwd code:
            !TSO2 = TSO2 + ESO2_sh(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_sh)
            ! adj code:
            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_sh) = 
     &         EMS_SF_ADJ(I,J,1,IDADJ_ESO2_sh) + ESO2_sh(I,J) * TSO2_ADJ
     
            ! Sum of anthro (surface + 100m), biomass, biofuel SO2 at (I,J)
            ! fwd code:
            !TSO2  =  SO2an(I,J,1) * EMS_SF(I,J,1,IDADJ_ESO2_an1)
            !       + SO2an(I,J,2) * EMS_SF(I,J,1,IDADJ_ESO2_an2)
            !       + ESO2_bb(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_bb)
            !       + SO2bf(I,J)   * EMS_SF(I,J,1,IDADJ_ESO2_bf)
            ! adj code:
            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an1) = 
     &           EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an1) 
     &         + SO2an(I,J,1) * TSO2_ADJ

            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an2) = 
     &           EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an2) 
     &         + SO2an(I,J,2) * TSO2_ADJ

            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bb) =  
     &           EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bb) 
     &         + ESO2_bb(I,J) * TSO2_ADJ

            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bf) =  
     &           EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bf) 
     &         + SO2bf(I,J) * TSO2_ADJ


         !===============================================================
         ! If PBL top occurs lower than or close to the top of level 2,
         ! then then surface SO2 goes into level 1 and the smokestack
         ! stack SO2 goes into level 2.
         !===============================================================
         ELSE
 

            ! Also add ship exhaust SO2 into surface if necessary 
            ! (bec, bmy, 5/20/04)
            ! fwd code:
            !SO2(1) = SO2(1) 
            !       + ESO2_sh(I,J) * EMS_SF(I,J,IDADJ_ESO2_sh)
            ! adj code:
            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_sh) =
     &         EMS_SF_ADJ(I,J,1,IDADJ_ESO2_sh) +
     &         ESO2_sh(I,J) * SO2_ADJ(1)

            ! fwd code:       
            !SO2(1) = SO2an(I,J,1) * EMS_SF(I,J,1,IDADJ_ESO2_an1)
            !       + ESO2_bb(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_bb)  
            !       + SO2bf(I,J)   * EMS_SF(I,J,1,IDADJ_ESO2_bf)
            !SO2(2) = SO2an(I,J,2) * EMS_SF(I,J,1,IDADJ_ESO2_an2)
            ! adj code:       
            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an1) =
     &         EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an1) +
     &         SO2an(I,J,1) * SO2_ADJ(1)

            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bb) =
     &         EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bb) +
     &         ESO2_bb(I,J) * SO2_ADJ(1)

            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bf) =
     &         EMS_SF_ADJ(I,J,1,IDADJ_ESO2_bf) +
     &         SO2bf(I,J) * SO2_ADJ(1)

            EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an2) =
     &         EMS_SF_ADJ(I,J,1,IDADJ_ESO2_an2) +
     &         SO2an(I,J,2) * SO2_ADJ(2)


         ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Optional diagnostic -- sensitivity w.r.t. absolute emissions (dkh, 02/17/11) 
      IF ( LEMS_ABS ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, SO2, TSO2, FEMIS, SO2SRC )
!$OMP+PRIVATE( SO2SRC_ADJ, SO2_ADJ,      TSO2_ADJ      )
!$OMP+SCHEDULE( DYNAMIC )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Top of the boundary layer
            NTOP = CEILING( (GET_PBL_TOP_L( I, J )) ) 
   
            DO L = 1, LLPAR
               
               ! Add SO2 to TC array [kg/box/timestep]
               ! fwd code:
               !TC(I,J,L) = TC(I,J,L) + ( SO2SRC * DTSRCE )
               ! adj code:
               SO2SRC_ADJ = TC_ADJ(I,J,L) * DTSRCE 
   
               ! SO2 emissions [kg/box/s]
               ! fwd code:
               !SO2SRC = SO2(L)         + ESO2_ac(I,J,L) + 
               !         ESO2_nv(I,J,L) + ESO2_ev(I,J,L) 
               ! adj code:
               SO2_ADJ(L) = SO2SRC_ADJ 
               ! adj code not implemented yet
               ! ESO2_ac_ADJ(I,J,L) = ... etc
   
            ENDDO
   
            !===============================================================
            ! Partition the total anthro and biomass SO2 emissions thru
            ! the entire boundary layer (if PBL top is higher than level 2)
            !===============================================================
            IF ( NTOP > 2 ) THEN
   
               TSO2_ADJ = 0d0
   
               ! Loop thru levels in the PBL
               DO L  = 1, NTOP
                    
                  ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
                  FEMIS  = GET_FRAC_OF_PBL( I, J, L )
                    
                  ! Partition total SO2 into level K
                  ! fwd_code:
                  !SO2(L) = FEMIS * TSO2
                  ! adj_code:
                  TSO2_ADJ = TSO2_ADJ + SO2_ADJ(L) * FEMIS
   
   
               ENDDO

               ! Also add SO2 from ship exhaust if necessary (bec, bmy, 5/20/04)
               ! fwd code:
               !TSO2 = TSO2 + ESO2_sh(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_sh)
               ! adj code: also convert to ( J / kg / box / timestep )
               EMS_ADJ(I,J,1,IDADJ_ESO2_sh) 
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_sh) 
     &            + TSO2_ADJ / DTSRCE
     
               ! Sum of anthro (surface + 100m), biomass, biofuel SO2 at (I,J)
               ! fwd code:
               !TSO2  =  SO2an(I,J,1) * EMS_SF(I,J,1,IDADJ_ESO2_an1)
               !       + SO2an(I,J,2) * EMS_SF(I,J,1,IDADJ_ESO2_an2)
               !       + ESO2_bb(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_bb)
               !       + SO2bf(I,J)   * EMS_SF(I,J,1,IDADJ_ESO2_bf)
               ! adj code: also convert to J / ( kg/box/timestep )
               EMS_ADJ(I,J,1,IDADJ_ESO2_an1) 
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_an1) 
     &            + TSO2_ADJ / DTSRCE

               EMS_ADJ(I,J,1,IDADJ_ESO2_an2) 
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_an2) 
     &            + TSO2_ADJ / DTSRCE

               EMS_ADJ(I,J,1,IDADJ_ESO2_bb) 
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_bb) 
     &            + TSO2_ADJ / DTSRCE

               EMS_ADJ(I,J,1,IDADJ_ESO2_bf)
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_bf) 
     &            + TSO2_ADJ / DTSRCE


            !===============================================================
            ! If PBL top occurs lower than or close to the top of level 2,
            ! then then surface SO2 goes into level 1 and the smokestack
            ! stack SO2 goes into level 2.
            !===============================================================
            ELSE
 

               ! Also add ship exhaust SO2 into surface if necessary 
               ! (bec, bmy, 5/20/04)
               ! fwd code:
               !SO2(1) = SO2(1) 
               !       + ESO2_sh(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_sh)
               ! adj code: also convert to J / (kg/box/timestep)
               EMS_ADJ(I,J,1,IDADJ_ESO2_sh)
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_sh)
     &            + SO2_ADJ(1) / DTSRCE

               ! fwd code:       
               !SO2(1) = SO2an(I,J,1) * EMS_SF(I,J,1,IDADJ_ESO2_an1)
               !       + ESO2_bb(I,J) * EMS_SF(I,J,1,IDADJ_ESO2_bb)  
               !       + SO2bf(I,J)   * EMS_SF(I,J,1,IDADJ_ESO2_bf)
               !SO2(2) = SO2an(I,J,2) * EMS_SF(I,J,1,IDADJ_ESO2_an2)
               ! adj code: also convert to kg/box/timestep
               EMS_ADJ(I,J,1,IDADJ_ESO2_an1)
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_an1)
     &            + SO2_ADJ(1) / DTSRCE

               EMS_ADJ(I,J,1,IDADJ_ESO2_bb)
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_bb)
     &            + SO2_ADJ(1) / DTSRCE

               EMS_ADJ(I,J,1,IDADJ_ESO2_bf)
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_bf)
     &            + SO2_ADJ(1) / DTSRCE

               EMS_ADJ(I,J,1,IDADJ_ESO2_an2)
     &            = EMS_ADJ(I,J,1,IDADJ_ESO2_an2)
     &            + SO2_ADJ(2) / DTSRCE


            ENDIF 

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF  

      ! Return to calling program
      END SUBROUTINE SRCSO2_ADJ

!------------------------------------------------------------------------------
!
!      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )
!!
!!******************************************************************************
!!  Function GET_OH returns OH from SMVGEAR's CSPEC array (for coupled runs)
!!  or monthly mean OH (for offline runs).  Imposes a diurnal variation on
!!  OH for offline simulations. (bmy, 12/16/02, 7/20/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!!
!!  NOTES:
!!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!!  (3 ) Now reference ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM from 
!!        "tracer_mod.f".  Also replace IJSURF w/ an analytic function. 
!!        (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE COMODE_MOD,    ONLY : CSPEC, JLOP
!      USE DAO_MOD,       ONLY : SUNCOS
!      USE ERROR_MOD,     ONLY : ERROR_STOP
!      USE GLOBAL_OH_MOD, ONLY : OH
!      USE TIME_MOD,      ONLY : GET_TS_CHEM
!      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
!      USE TRACERID_MOD,  ONLY : IDOH
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: I, J, L
!
!      ! Local variables
!      INTEGER             :: JLOOP
!      TYPE (XPLEX)              :: OH_MOLEC_CM3
! 
!      !=================================================================
!      ! GET_OH begins here!
!      !=================================================================
!      IF ( ITS_A_FULLCHEM_SIM() ) THEN
!
!         !---------------------
!         ! Coupled simulation
!         !---------------------
!
!         ! JLOOP = SMVGEAR 1-D grid box index
!         JLOOP = JLOP(I,J,L)
!
!         ! Take OH from the SMVGEAR array CSPEC
!         ! OH is defined only in the troposphere
!         IF ( JLOOP > 0 ) THEN
!            OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
!         ELSE
!            OH_MOLEC_CM3 = 0d0
!         ENDIF
!
!      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
!
!         !---------------------
!         ! Offline simulation
!         !---------------------
!
!         ! 1-D grid box index for SUNCOS
!         JLOOP = ( (J-1) * IIPAR ) + I
!
!         ! Test for sunlight...
!         IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN
!
!            ! Impose a diurnal variation on OH during the day
!            OH_MOLEC_CM3 = OH(I,J,L)                      *           
!     &                     ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
!     &                     ( 1440d0        / GET_TS_CHEM() )
!
!            ! Make sure OH is not negative
!            OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
!               
!         ELSE
!
!            ! At night, OH goes to zero
!            OH_MOLEC_CM3 = 0d0
!
!         ENDIF
!
!      ELSE
!
!         !---------------------
!         ! Invalid simulation
!         !---------------------         
!         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')
!
!      ENDIF
!
!      ! Return to calling program
!      END FUNCTION GET_OH
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE SET_OH( I, J, L, OH ) 
!!
!!******************************************************************************
!!  Function SET_OH saves the modified OH value back to SMVGEAR's CSPEC array
!!  for coupled sulfate/aerosol simulations. (bmy, 12/16/02)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!!  (4  ) OH        (TYPE (XPLEX) ) : OH at grid box (I,J,L) to be saved into CSPEC
!!
!!  NOTES:
!!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 12/16/02)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE COMODE_MOD,   ONLY : CSPEC, JLOP
!      USE TRACERID_MOD, ONLY : IDOH
! 
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: I, J, L
!      TYPE (XPLEX),  INTENT(IN) :: OH
!
!      ! Local variables
!      INTEGER             :: JLOOP
!
!      !=================================================================
!      ! SET_OH begins here!
!      !=================================================================
!
!      ! JLOOP = SMVGEAR 1-D grid box index
!      JLOOP = JLOP(I,J,L) 
!
!      ! Replace OH into CSPEC(troposphere only)
!      IF ( JLOOP > 0 ) THEN
!         CSPEC(JLOOP,IDOH) = OH
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE SET_OH
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_NO3( I, J, L ) RESULT( NO3_MOLEC_CM3 ) 
!!
!!******************************************************************************
!!  Function GET_NO3 returns NO3 from SMVGEAR's CSPEC array (for coupled runs)
!!  or monthly mean OH (for offline runs).  For offline runs, the concentration
!!  of NO3 is set to zero during the day. (rjp, bmy, 12/16/02)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!!
!!  NOTES:
!!  (1 ) Now references ERROR_STOP from "error_mod.f".  We also assume that
!!        SETTRACE has been called to define IDNO3.  Now also set NO3 to 
!!        zero during the day. (rjp, bmy, 12/16/02)
!!  (2 ) Now reference ITS_A_FULLCHEM_SIM and ITS_AN_AEROSOL_SIM from 
!!        "tracer_mod.f".  Also remove reference to CMN.   Also replace
!!        IJSURF with an analytic function. (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE COMODE_MOD,     ONLY : CSPEC, JLOP
!      USE DAO_MOD,        ONLY : AD,    SUNCOS
!      USE ERROR_MOD,      ONLY : ERROR_STOP
!      USE GLOBAL_NO3_MOD, ONLY : NO3
!      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
!      USE TRACERID_MOD,   ONLY : IDNO3
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: I, J, L
!
!      ! Local variables
!      INTEGER             :: JLOOP
!      TYPE (XPLEX)              :: NO3_MOLEC_CM3
! 
!      ! External functions
!      TYPE (XPLEX),  EXTERNAL   :: BOXVL
!
!      !=================================================================
!      ! GET_NO3 begins here!
!      !=================================================================
!      IF ( ITS_A_FULLCHEM_SIM() ) THEN
!
!         !--------------------
!         ! Coupled simulation
!         !--------------------
!            
!         ! 1-D SMVGEAR grid box index
!         JLOOP = JLOP(I,J,L)
!
!         ! Take NO3 from the SMVGEAR array CSPEC
!         ! NO3 is defined only in the troposphere
!         IF ( JLOOP > 0 ) THEN
!            NO3_MOLEC_CM3 = CSPEC(JLOOP,IDNO3)
!         ELSE
!            NO3_MOLEC_CM3 = 0d0
!         ENDIF
!
!      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
!
!         !==============================================================  
!         ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
!         ! in [v/v].  Convert these to [molec/cm3] as follows:
!         !
!         !  vol NO3   moles NO3    kg air     kg NO3/mole NO3
!         !  ------- = --------- * -------- * ---------------- =  kg NO3 
!         !  vol air   moles air      1        kg air/mole air
!         !
!         ! And then we convert [kg NO3] to [molec NO3/cm3] by:
!         !  
!         !  kg NO3   molec NO3   mole NO3     1     molec NO3
!         !  ------ * --------- * -------- * ----- = --------- 
!         !     1     mole NO3     kg NO3     cm3       cm3
!         !          ^                    ^
!         !          |____________________|  
!         !            this is XNUMOL_NO3
!         !
!         ! If at nighttime, use the monthly mean NO3 concentration from
!         ! the NO3 array of "global_no3_mod.f".  If during the daytime,
!         ! set the NO3 concentration to zero.  We don't have to relax to 
!         ! the monthly mean concentration every 3 hours (as for HNO3) 
!         ! since NO3 has a very short lifetime. (rjp, bmy, 12/16/02) 
!         !==============================================================
!
!         ! 1-D grid box index for SUNCOS
!         JLOOP = ( (J-1) * IIPAR ) + I
!
!         ! Test if daylight
!         IF ( SUNCOS(JLOOP) > 0d0 ) THEN
!
!            ! NO3 goes to zero during the day
!            NO3_MOLEC_CM3 = 0d0
!              
!         ELSE
!
!            ! At night: Get NO3 [v/v] and convert it to [kg]
!            NO3_MOLEC_CM3 = NO3(I,J,L) * AD(I,J,L) * ( 62d0/28.97d0 ) 
!               
!            ! Convert NO3 from [kg] to [molec/cm3]
!            NO3_MOLEC_CM3 = NO3_MOLEC_CM3 * XNUMOL_NO3 / BOXVL(I,J,L)
!                  
!         ENDIF
!            
!         ! Make sure NO3 is not negative
!         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0d0 )
!
!      ELSE
!
!         !--------------------
!         ! Invalid simulation
!         !--------------------
!         CALL ERROR_STOP( 'Invalid NSRCX!','GET_NO3 (sulfate_mod.f)')
!
!      ENDIF
!
!      ! Return to calling program
!      END FUNCTION GET_NO3
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE SET_NO3( I, J, L, NO3 ) 
!!
!!******************************************************************************
!!  Function SET_NO3 saves the modified NO3 value back to SMVGEAR's CSPEC array
!!  for coupled lfate/aerosol simulations. (rjp, bmy, 12/16/02, 7/20/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!!  (4  ) NO3     (TYPE (XPLEX) ) : OH at grid box (I,J,L) to be saved into CSPEC
!!
!!  NOTES:
!!  (1 ) We assume SETTRACE has been called to define IDNO3. (bmy, 12/16/02)
!!  (2 ) Remove references to "error_mod.f" and CMN (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE COMODE_MOD,   ONLY : CSPEC, JLOP
!      USE TRACERID_MOD, ONLY : IDNO3
!      
!#     include "CMN_SIZE"  ! Size parameters 
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: I, J, L
!      TYPE (XPLEX),  INTENT(IN) :: NO3
!
!      ! Local variables
!      INTEGER             :: JLOOP
!
!      !=================================================================
!      ! SET_NO3 begins here!
!      !=================================================================
!
!      ! 1-D grid box index for CSPEC
!      JLOOP = JLOP(I,J,L) 
!
!      ! Replace OH into CSPEC (troposphere only)
!      IF ( JLOOP > 0 ) THEN
!         CSPEC(JLOOP,IDNO3) = NO3
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE SET_NO3
!      
!!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_VV )
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Now remove reference to CMN, it's obsolete.  (bmy, 8/22/05)
!  (4 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP, VOLUME
      USE DAO_MOD,      ONLY : AIRDEN
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM  
      USE TRACER_MOD,   ONLY : XNUMOLAIR
      USE TRACERID_MOD, ONLY : IDO3

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)   :: I, J, L

      ! Local variables
      INTEGER               :: JLOOP
      TYPE (XPLEX)                :: O3_VV
 
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Get O3 from CSPEC [molec/cm3] and convert it to [v/v]
         ! O3 data will only be defined below the tropopause
         IF ( JLOOP  > 0 ) THEN
            O3_VV = ( CSPEC(JLOOP,IDO3) * 1d6       ) / 
     &              ( AIRDEN(L,I,J)     * XNUMOLAIR )
         ELSE
            O3_VV = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
         
      !   !--------------------
      !   ! Offline simulation
      !   !--------------------
      !
      !   ! Get O3 [v/v] for this gridbox & month
      !   ! O3 data will only be defined below the tropopause
      !   IF ( L <= LLTROP ) THEN
      !      O3_VV = O3m(I,J,L)
      !   ELSE
      !      O3_VV = 0d0
      !   ENDIF
      !
      !ELSE

         !--------------------
         ! Invalid simulation
         !--------------------         
         CALL ERROR_STOP( 'Invalid NSRCX!', 
     $                    'GET_O3 (sulfate_adj_mod.f)')

      ENDIF

      ! Return to calling program
      END FUNCTION GET_O3

!!------------------------------------------------------------------------------
!      
!      SUBROUTINE READ_NONERUP_VOLC
!!
!!******************************************************************************
!!  Subroutine READ_NONERUP_VOLC reads SO2 emissions from non-eruptive
!!  volcanoes. (rjp, bdf, bmy, 9/19/02, 10/3/05)
!!
!!  NOTES:
!!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (3 ) Now read files from "sulfate_sim_200508/" (bmy, 7/28/05)
!!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!******************************************************************************
!! 
!      ! References to F90 modules
!      USE BPCH2_MOD,     ONLY : GET_RES_EXT, GET_TAU0, READ_BPCH2
!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Local variables
!      INTEGER              :: I, IOS, J, K, L
!      TYPE (XPLEX)               :: EE
!      CHARACTER(LEN=255)   :: FILENAME
!
!      !=================================================================
!      ! READ_NONERUP_VOLC begins here!
!      !=================================================================
!
!      ! Initialize
!      K        = 1
!      Env      = 0.d0
!      FILENAME = TRIM( DATA_DIR )                  // 
!     &           'sulfate_sim_200508/volcano.con.' // GET_RES_EXT()
!
!      !=================================================================
!      ! Read NON-eruptive volcanic SO2 emission (GEIA) into Env.  
!      ! Convert Env from [Mg SO2/box/day] to [kg SO2/box/s].
!      !=================================================================
!
!      ! Fancy output
!      WRITE( 6, 100 ) TRIM( FILENAME ) 
! 100  FORMAT( '     - READ_NONERUP_VOLC: Reading ', a )
!     
!      ! Open file 
!      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
!      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_nonerup_volc:1' )
!
!      ! Read header lines
!      DO L = 1, 2
!         READ( IU_FILE, '(a)', IOSTAT=IOS )             
!         IF ( IOS > 0 ) THEN
!            CALL IOERROR( IOS, IU_FILE, 'read_nonerup_volc:2' )
!         ENDIF
!      ENDDO
!
!      ! Read data values
!      DO 
!         READ( IU_FILE, '(49x,i4,e11.3,1x,2i4)', IOSTAT=IOS ) 
!     &        IELVn(k), EE, INV(K), JNV(k) 
!         
!         ! Check for EOF
!         IF ( IOS < 0 ) EXIT
!
!         ! Trap I/O error
!         IF ( IOS > 0 ) THEN
!            CALL IOERROR( IOS, IU_FILE, 'read_nonerup_volc:3' )
!         ENDIF
!
!         ! Unit conversion: [Mg SO2/box/day] -> [kg SO2/box/s]
!         Env(k) = EE * 1000.d0 / ( 24.d0 * 3600.d0 )
!
!         ! Increment counter
!         K = K + 1 
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_FILE )
!
!      ! NNVOL = Number of non-eruptive volcanoes
!      NNVOL = K - 1
!
!      ! Return to calling program
!      END SUBROUTINE READ_NONERUP_VOLC
!
!!------------------------------------------------------------------------------
!      
!      SUBROUTINE READ_ERUP_VOLC
!!
!!***************************************************************************** 
!!  Subroutine READ_ERUP_VOLC reads SO2 emissions from eruptive
!!  volcanoes. (rjp, bdf, bmy, 9/19/02, 10/3/05)
!!
!!  NOTES:
!!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (3 ) Now read files from "sulfate_sim_200508/" (bmy, 7/28/05)
!!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!*****************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,     ONLY : GET_RES_EXT, GET_TAU0, READ_BPCH2
!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
!      
!#     include "CMN_SIZE"   ! Size parameters
! 
!      ! Local variables
!      INTEGER              :: I, IOS, IUNIT, J, K, L, M
!      TYPE (XPLEX)               :: A, B, Fe, X, EE
!      CHARACTER(LEN=255)   :: FILENAME
!
!      !==================================================================
!      ! READ_ERUP_VOLC begins here
!      !==================================================================
!
!      ! Initialize
!      K        = 1
!      Eev(:)   = 0.d0
!      FILENAME = TRIM( DATA_DIR )                        // 
!     &           'sulfate_sim_200508/volcano.erup.1990.' // 
!     &           GET_RES_EXT()
!      
!      !==================================================================
!      ! Read eruptive volcanic SO2 emission (based on Smithsonian data 
!      ! base, SO2 emission and cloud height are a function of VEI.  
!      ! Data are over-written if TOMS observations are available.  
!      ! Also define a slab with a thickness of 1/3 of the cloud column, 
!      ! and SO2 are emitted uniformely within the slab.  
!      !
!      ! Convert Ee from [kton SO2] to [kg SO2/box] and store in Eev.
!      ! ESO2_ev(i,j,l) in [kg SO2/box/s] will be calculated in SRCSO2.  
!      !==================================================================
!      
!      ! Fancy output
!      WRITE( 6, 100 ) TRIM( FILENAME ) 
! 100  FORMAT( '     - READ_ERUP_VOLC: Reading ', a )
!   
!      ! Open file 
!      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
!      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_erup_volc:1' )
!
!      ! Read header lines
!      DO L = 1, 2
!         READ( IU_FILE, '(a)', IOSTAT=IOS )
!         IF ( IOS > 0 ) THEN
!            CALL IOERROR( IOS, IU_FILE, 'read_erup_volc:2' )
!         ENDIF
!      ENDDO
!
!         ! Read data values
!      DO 
!         READ( IU_FILE, '(47x,3i6,6x,i6,es11.3,1x,2i4)', IOSTAT=IOS )
!     &        IELVe(K), IDAYs(K), IDAYe(K), IHGHT(K), 
!     &        Ee,       IEV(K),   JEV(K)
!         
!         ! Check for EOF
!         IF ( IOS < 0 ) EXIT
!
!          ! Trap I/O error
!         IF ( IOS > 0 ) THEN
!            CALL IOERROR( IOS, IU_FILE, 'sulfate_readyr:6' )
!         ENDIF
!
!         ! Unit conversion: [kton SO2/box/event] -> [kg SO2/box/event]
!         Eev(k) = Ee * 1.d6
!
!         ! Increment count
!         K = K + 1
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_FILE )
!
!      ! NEVOL = Number of eruptive volcanoes
!      NEVOL = K - 1
!
!      ! Return to calling program
!      END SUBROUTINE READ_ERUP_VOLC
!      
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_ANTHRO_SOx( THISMONTH, NSEASON )
!!
!!******************************************************************************
!!  Suborutine READ_ANTHRO_SOx reads the anthropogenic SOx from disk, 
!!  and partitions it into anthropogenic SO2 and SO4. 
!!  (rjp, bdf, bmy, 9/20/02, 10/31/08)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!  
!!  NOTES:
!!  (1 ) Now use functions GET_XMID and GET_YMID to compute lon and lat
!!        centers of grid box (I,J).  Now replace DXYP(JREF)*1d4 with routine
!!        GET_AREA_CM2 of "grid_mod.f".  Now use functions GET_MONTH and
!!        GET_YEAR of time_mod.f".  Now call READ_BPCH2 with QUIET=.TRUE. 
!!        (bmy, 3/27/03)
!!  (2 ) Now references DATA_DIR from "directory_mod.f".  Also removed
!!        reference to CMN, it's not needed. (bmy, 7/20/04)
!!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both
!!        GCAP and GEOS grids (bmy, 8/16/05)
!!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!!  (6 ) Now computes future SOx emissions (swu, bmy, 5/30/06)
!!  (7 ) Now can read either EDGAR or GEIA emissions (avd, bmy, 7/14/06)
!!  (8 ) Now overwrite David Streets' SO2, if necessary (yxw, bmy, 8/14/06)
!!  (9 ) Now accounts for FSCLYR (phs, 3/17/08)
!!  (9 ) Bug fix: Using tracer #30 in the call to GET_STREETS_ANTHRO can cause
!!        problems when adding or removing species.  Replace w/ IDTNH3.
!!        (dkh, 10/31/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR
!      USE EDGAR_MOD,            ONLY : GET_EDGAR_ANTH_SO2
!      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO
!      USE EMEP_MOD,             ONLY : GET_EUROPE_MASK
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
!      USE GRID_MOD,             ONLY : GET_XMID,        GET_YMID
!      USE GRID_MOD,             ONLY : GET_AREA_CM2
!      USE LOGICAL_MOD,          ONLY : LFUTURE,         LEDGARSOx
!      USE LOGICAL_MOD,          ONLY : LSTREETS,        LEMEP
!      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
!      USE STREETS_ANTHRO_MOD,   ONLY : GET_STREETS_ANTHRO
!      USE TIME_MOD,             ONLY : GET_YEAR
!      USE TRACER_MOD,           ONLY : XNUMOL
!      USE TRACERID_MOD,         ONLY : IDTSO2, IDTSO4
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!      USE SCALE_ANTHRO_MOD,     ONLY : GET_ANNUAL_SCALAR
!
!#     include "CMN_SIZE"             ! Size parameters
!#     include "CMN_O3"               ! FSCALYR
!                                    
!      ! Arguments                   
!      INTEGER, INTENT(IN)           :: THISMONTH, NSEASON
!                                    
!      ! Local variables             
!      INTEGER                       :: I, J, L, IX, JX, IOS
!      INTEGER, SAVE                 :: LASTYEAR = -99
!      INTEGER                       :: SCALEYEAR
!      TYPE (XPLEX)                        :: E_SOx(IGLOB,JGLOB,2)
!      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                        :: XTAU, Fe, X, Y, AREA_CM2
!      TYPE (XPLEX)                        :: EDG_SO2
!      CHARACTER(LEN=4)              :: SYEAR
!      CHARACTER(LEN=255)            :: FILENAME
!
!      !=================================================================
!      ! READ_ANTHRO_SOx begins here!
!      !=================================================================
!         
!      IF ( FSCALYR < 0 ) THEN
!         SCALEYEAR = GET_YEAR()
!      ELSE
!         SCALEYEAR = FSCALYR
!      ENDIF
!
!
!      IF ( LEDGARSOx ) THEN
!
!         !==============================================================
!         ! Use EDGAR SOx emissions
!         !
!         ! Partition SOx into SO2 and SO4, according to the following
!         ! fractions (cf Chin et al, 2000):
!         ! 
!         ! Europe     [ 36N-78N,  12.5W-60.0E ]:  5.0% of SOx is SO4
!         !                                       95.0% of SOx is SO2   
!         !
!         ! N. America [ 26N-74N, 167.5W-52.5W ]:  1.4% of SOx is SO4
!         !                                       98.6% of SOx is SO2
!         !
!         ! Everywhere else:                       3.0% of SOx is SO4
!         !                                       97.0% of SOx is SO2
!         !==============================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, X, Y, EDG_SO2, Fe ) 
!
!         ! Loop over latitudes
!         DO J = 1, JJPAR
!
!            ! Latitude [degrees]
!            Y = GET_YMID( J )
!
!            ! Loop over longitudes
!            DO I = 1, IIPAR
!
!               ! Longitude [degrees]
!               X  = GET_XMID( I )
!
!               ! Get EDGAR SO2 emissions [kg/s]
!               ! NOTE: Future emissions are already applied!
!               EDG_SO2 = GET_EDGAR_ANTH_SO2( I, J, KG_S=.TRUE. ) 
!
!               ! If we are using David Streets' emissions ...
!               IF ( LSTREETS ) THEN
!
!                  ! If we are over the SE Asia region ...
!                  IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN
!
!                     ! Overwrite EDGAR SO2 w/ David Streets' [kg SO2/s]
!                     EDG_SO2 = GET_STREETS_ANTHRO( I,      J, 
!     &                                             IDTSO2, KG_S=.TRUE. )
!                     
!                     ! Streets 2006 includes biofuels.
!                     IF ( SCALEYEAR >= 2006 ) ESO2_bf(I,J) = 0d0
!
!                  ENDIF
!               ENDIF
!
!               ! If we are using EMEP over Europe...
!               IF ( LEMEP ) THEN
!
!                  IF (GET_EUROPE_MASK(I,J) > 0d0) THEN
!
!!-----------------------------------------------------------------------------
!! Prior to 11/14/08:
!! BUG FIX: Using tracer #26 in the call to GET_EMEP_ANTHRO can cause 
!! problems when adding or removing species.  Replace w/ IDTSO2. 
!! (phs, 11/14/08) 
!!                     EDG_SO2 = GET_EMEP_ANTHRO( I, J, 26, KG_S=.TRUE. )
!!-----------------------------------------------------------------------------
!                     EDG_SO2 = GET_EMEP_ANTHRO( I,      J,
!     $                                          IDTSO2, KG_S=.TRUE. )
!
!                  ENDIF
!
!               ENDIF
!
!               ! Compute SO4/SOx fraction for EUROPE
!               IF ( ( X >= -12.5 .and. X <= 60.0 )  .and.
!     &              ( Y >=  36.0 .and. Y <= 78.0 ) ) THEN
!                  Fe = 0.05d0
!
!               ! Compute SO4/SOx fraction for NORTH AMERICA
!               ELSE IF ( ( X >= -167.5 .and. X <= -52.5 )  .and.
!     &                   ( Y >=   26.0 .and. Y <=  74.0 ) ) THEN
!                  Fe = 0.014d0
!
!               ! Compute SO4/SOx fraction for EVERYWHERE ELSE
!               ELSE
!                  Fe = 0.03d0
!
!               ENDIF
!
!               ! Store SO2 emission [kg SO2/s]
!               ESO2_an(I,J,1) = EDG_SO2
!               ESO4_an(I,J,2) = 0d0
!
!               ! Compute SO4 from SO2 [kg SO4/s]
!               ESO4_an(I,J,1) = EDG_SO2 * Fe / ( 1.d0 - Fe )
!               ESO4_an(I,J,2) = 0d0
!            ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!      ELSE
!
!         !==============================================================
!         ! Use GEIA SOx emissions
!         !==============================================================
!
!         ! Define filename
!         FILENAME = TRIM( DATA_DIR )                          //
!     &              'fossil_200104/merge_nobiofuels.'         //
!     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT() 
!     
!         ! Echo output
!         WRITE( 6, 100 ) TRIM( FILENAME )
! 100     FORMAT( '     - READ_ANTHRO_SOx: Reading ', a )
!
!         ! Pick the right TAU value for the given season
!         ! Seasons: 1=DJF, 2=MAM, 3=JJA, 4=SON
!         SELECT CASE ( NSEASON )
!            CASE ( 1 )
!               XTAU = -744d0
!            CASE ( 2 )
!               XTAU = 1416d0
!            CASE ( 3 )
!               XTAU = 3624d0
!            CASE ( 4 )
!               XTAU = 5832d0
!         END SELECT
!
!         ! Read anthropogenic SOx [molec/cm2/s] 
!         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 27, 
!     &                    XTAU,      IGLOB,     JGLOB,     
!     &                    2,         E_SOx,     QUIET=.TRUE. )
!
!         !=================================================================
!         ! Read in yearly SO2 scale factors here
!         ! (For now we only have 1998, deal w/ other years later)
!         !=================================================================
!!-- prior to 3/11/08
!         !IF ( LASTYEAR < 0 ) THEN
!         !
!         !   ! put in SOX scale year here (hardwired to 1998 for now)
!         !   SYEAR    = '1998'
!         !   FILENAME = TRIM( DATA_DIR )                    // 
!         !&                 'sulfate_sim_200508/scalefoss.SOx.' //  
!         !&                 GET_RES_EXT()  // '.' // SYEAR
!         !
!         !   ! Echo output
!         !   WRITE( 6, 100 ) TRIM( FILENAME )
!         !
!         !   ! Get TAU value (use Jan 1, 1998 for scale factors)
!         !   XTAU = GET_TAU0( 1, 1, 1998 )
!         !
!         !   ! Read anthropogenic SOx [molec/cm2/s] 
!         !   CALL READ_BPCH2( FILENAME, 'SCALFOSS', 3, 
!         !&                       XTAU,      IGLOB,     JGLOB,     
!         !&                       1,         ARRAY,     QUIET=.TRUE. )
!         !
!         !   ! Cast from TYPE (XPLEX) to COMPLEX*16
!         !   CALL TRANSFER_2D( ARRAY(:,:,1), SOx_SCALE )
!         !
!         !   ! Reset LASTYEAR
!         !   LASTYEAR = GET_YEAR()
!         !ENDIF
!
!
!         ! Get annual scalar factor (amv, 08/24/07)
!         CALL GET_ANNUAL_SCALAR( 73, 1985, SCALEYEAR, SOx_SCALE )
!
!         !==============================================================
!         ! Partition SOx into SO2 and SO4, according to the following
!         ! fractions (cf Chin et al, 2000):
!         ! 
!         ! Europe     [ 36N-78N,  12.5W-60.0E ]:  5.0% of SOx is SO4
!         !                                       95.0% of SOx is SO2   
!         !
!         ! N. America [ 26N-74N, 167.5W-52.5W ]:  1.4% of SOx is SO4
!         !                                       98.6% of SOx is SO2
!         !
!         ! Everywhere else:                       3.0% of SOx is SO4
!         !                                       97.0% of SOx is SO2
!         !==============================================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, AREA_CM2, Y, X, Fe )
!         DO L = 1, 2
!         DO J = 1, JJPAR
!
!            ! Grid box surface area [cm2]
!            AREA_CM2 = GET_AREA_CM2( J )
!
!            ! Latitude [degrees]
!            Y = GET_YMID( J )
!
!            DO I = 1, IIPAR
!
!               ! Longitude [degrees]  
!               X = GET_XMID( I )
!
!               ! First scale SOx to the given fossil fuel year
!               E_SOx(I,J,L) = E_SOx(I,J,L) * SOx_SCALE(I,J)
!            
!               ! Compute future SOx emissions (if necessary)
!               IF ( LFUTURE ) THEN
!                  E_SOx(I,J,L) = E_SOx(I,J,L)                  * 
!     &                           GET_FUTURE_SCALE_SO2ff( I, J )
!               ENDIF
!
!               ! Compute SO4/SOx fraction for EUROPE
!               IF ( ( X >= -12.5 .and. X <= 60.0 )  .and. 
!     &              ( Y >=  36.0 .and. Y <= 78.0 ) ) THEN
!                  Fe = 0.05d0
!
!               ! Compute SO4/SOx fraction for NORTH AMERICA
!               ELSE IF ( ( X >= -167.5 .and. X <= -52.5 )  .and.   
!     &                   ( Y >=   26.0 .and. Y <=  74.0 ) ) THEN
!                  Fe = 0.014d0
! 
!               ! Compute SO4/SOx fraction for EVERYWHERE ELSE
!               ELSE
!                  Fe = 0.03d0
!             
!               ENDIF
!         
!               ! Compute SO2 (tracer #2) from SOx
!               ! Convert from [molec SOx/cm2/s] to [kg SO2/box/s]
!               ESO2_an(I,J,L) = E_SOx(I,J,L) * ( 1.d0 - Fe ) * 
!     &                          AREA_CM2 / XNUMOL(IDTSO2)            
!
!               ! If we are using David Streets' emissions
!               ! Remember: those include BF if Year is GE 2006
!               IF ( LSTREETS ) THEN
!
!                  ! If we are over the SE Asia region
!                  IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN
!
!                     ! Overwrite GEIA SO2 w/ David Streets' SO2 [kg SO2/s]
!                     ESO2_an(I,J,1) = GET_STREETS_ANTHRO( I, J, IDTSO2, 
!     &                                                    KG_S=.TRUE. )
!
!                     ! Zero 2nd level of emissions
!                     ESO2_an(I,J,2) = 0d0
!
!                  ENDIF
!               ENDIF
!
!               IF ( LEMEP ) THEN
!
!                  IF (GET_EUROPE_MASK(I,J) > 0d0 ) THEN
!
!                     ESO2_an(I,J,1) = GET_EMEP_ANTHRO( I, J, IDTSO2,
!     &                                                 KG_S=.TRUE. ) 
!
!                     ESO2_an(I,J,2) = 0d0
! 
!                  ENDIF
!
!               ENDIF
!
!!--- prior 6/23/08
!!      Now calculate SO4 from SO2, since SOx not available with STREETS and EMEP
!!               ! Compute SO4 (tracer #3) from SOx
!!               ! Convert from [molec SOx/cm2/s] to [kg SO4/box/s]
!!               ESO4_an(I,J,L) = E_SOx(I,J,L) * Fe *
!!     &                          AREA_CM2 / XNUMOL(IDTSO4)
!
!                ESO4_an(I,J,L) = ESO2_an(I,J,L) * Fe / (1.d0-Fe)
!
!            ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE READ_ANTHRO_SOx
! 
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_OCEAN_DMS( THISMONTH )
!!
!!***************************************************************************** 
!!  Subroutine READ_OCEAN_DMS reads seawater concentrations of DMS (nmol/L).
!!  (rjp, bdf, bmy, 9/20/02, 10/3/05)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (4 ) Now read files from "sulfate_sim_200508/".  Now read data for both
!!        GCAP and GEOS grids (bmy, 8/16/05) 
!!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!*****************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"  ! Size parameters 
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: THISMONTH
!      
!      ! Local variables
!      TYPE (XPLEX)              :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)              :: XTAU
!      CHARACTER(LEN=255)  :: FILENAME
!
!      !==================================================================
!      ! READ_OCEAN_DMS begins here!
!      !==================================================================
!
!      ! File name
!      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'sulfate_sim_200508/DMS_seawater.'       //
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!
!      ! Echo output
!      WRITE( 6, 100 ) TRIM( FILENAME )  
! 100  FORMAT( '     - READ_OCEAN_DMS: Reading ', a ) 
!
!      ! TAU value (use generic year 1985)
!      XTAU = GET_TAU0( THISMONTH, 1, 1985 )
!
!      ! Read ocean DMS [nmol/L]
!      CALL READ_BPCH2( FILENAME, 'BIOGSRCE',    25, 
!     &                 XTAU,      IGLOB,        JGLOB,      
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. ) 
!
!      ! Cast from TYPE (XPLEX) to COMPLEX*16 
!      CALL TRANSFER_2D( ARRAY(:,:,1), DMSo )
!      
!      ! Return to calling program
!      END SUBROUTINE READ_OCEAN_DMS
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_SST( THISMONTH, THISYEAR )
!!
!!***************************************************************************** 
!!  Subroutine READ_SST reads monthly mean sea surface temperatures.
!!  (rjp, bdf, bmy, 9/18/02, 11/17/05)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!  (2 ) THISYEAR  (INTEGER) : Current 4-digit year 
!!
!!  NOTES:
!!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (4 ) Now use interannual SST data from NOAA if present; otherwise use
!!        climatological SST data.  Now read data for both GCAP and GEOS 
!!        grids (bmy, 8/16/05) 
!!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (6 ) Now read int'annual SST data on the GEOS 1x1 grid (bmy, 11/17/05)
!!*****************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,  ONLY : DATA_DIR,        DATA_DIR_1x1
!      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
!      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"       ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)     :: THISMONTH, THISYEAR
!      
!      ! Local variables
!      TYPE (XPLEX)                  :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                  :: ARRAY1(I1x1,J1x1,1)
!      TYPE (XPLEX)                  :: XTAU
!      CHARACTER(LEN=4)        :: SYEAR
!      CHARACTER(LEN=255)      :: FILENAME
!
!      !==================================================================
!      ! READ_SST begins here!
!      !==================================================================
!
!      IF ( THISYEAR >= 1985 .and. THISYEAR <= 2004 ) THEN 
!         
!         !------------------------------------
!         ! Use interannual SST data from NOAA
!         ! Data exists for 1985 - 2004,
!         ! Add other years as necessary
!         !------------------------------------
!
!         ! Make a string for THISYEAR
!         WRITE( SYEAR, '(i4)' ) THISYEAR
!
!         ! File name
!         FILENAME = TRIM( DATA_DIR_1x1 )       // 
!     &              'SST_200508/SST.geos.1x1.' // SYEAR
!
!         ! Echo output
!         WRITE( 6, 100 ) TRIM( FILENAME )  
! 100     FORMAT( '     - READ_SST: Reading ', a ) 
!
!         ! TAU value (use this year)
!         XTAU = GET_TAU0( THISMONTH, 1, THISYEAR )
!
!         ! Read sea surface temperature [K]
!         CALL READ_BPCH2( FILENAME, 'GMAO-2D',      69, 
!     &                    XTAU,      I1x1,          J1x1,     
!     &                    1,         ARRAY1(:,:,1), QUIET=.TRUE. ) 
!
!         ! Regrid from 1x1 and cast to TYPE (XPLEX)
!         CALL DO_REGRID_1x1( 'K', ARRAY1, SSTEMP )
!
!      ELSE
!
!         !-------------------------------
!         ! Use climatological SST data 
!         !-------------------------------
!
!         ! File name
!         FILENAME = TRIM( DATA_DIR )          // 
!     &              'sulfate_sim_200508/SST.' // GET_NAME_EXT_2D() //
!     &              '.'                       // GET_RES_EXT()
!
!         ! Echo output
!         WRITE( 6, 100 ) TRIM( FILENAME )  
!
!         ! TAU value (use generic year 1985)
!         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
!
!         ! Read sea surface temperature [K]
!         CALL READ_BPCH2( FILENAME, 'DAO-FLDS',    5, 
!     &                    XTAU,      IGLOB,        JGLOB,     
!     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. ) 
!
!         ! Cast from TYPE (XPLEX) to COMPLEX*16 
!         CALL TRANSFER_2D( ARRAY(:,:,1), SSTEMP )
!
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE READ_SST
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_BIOFUEL_SO2( THISMONTH )
!!
!!******************************************************************************
!!  Subroutine READ_BIOFUEL_SO2 reads monthly mean biomass burning 
!!  emissions for SO2.  SOx is read in, and converted to SO2. 
!!  (rjp, bdf, bmy, phs, 1/16/03, 12/23/08)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!!  (2 ) Modified molar ratio of biomass burning SO2 per CO.  Added SO2
!!        emission from biofuel burning. (rjp, bmy, 1/16/03)
!!  (3 ) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!!        Now replace MONTH with the argument THISMONTH.  Now call READ_BPCH2
!!        with QUIET=.TRUE. (bmy, 3/27/03)
!!  (4 ) Now references DATA_DIR from "directory_mod.f".  Also removed
!!        references to CMN and CMN_SETUP. (bmy, 7/20/04)
!!  (5 ) Now can read either seasonal or interannual biomass burning emissions.
!!        Now references routines from both "logical_mod.f" and "time_mod.f".
!!        Now reads SO2 biomass emissions directly rather than computing
!!        it by mole fraction from CO. (rjp, bmy, 1/11/05)
!!  (6 ) Now read data for both GCAP and GEOS grids (bmy, 8/16/05) 
!!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (8 ) Now computes future biomass emissions, if necessary (swu, bmy, 5/30/06)
!!  (9 ) Now only read the biofuel, we have moved the biomass-reading code
!!        to "gc_biomass_mod.f" for compatibility with GFED2 biomass emissions
!!        (bmy, 9/27/06)
!!  (10) Now prevent seg fault if BIOMASS emissions are turned off.
!!        (bmy, 10/3/06)
!!  (11) Renamed READ_BIOFUEL_SO2, and move all biomass code to GET_BIOMASS_SO2
!!        to account for several GFED2 products (yc, phs, 12/23/08)   
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBSO2
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2bf
!      USE GRID_MOD,             ONLY : GET_AREA_CM2
!      USE LOGICAL_MOD,          ONLY : LBIOMASS,        LFUTURE
!      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
!      USE TRACER_MOD,           ONLY : XNUMOL
!      USE TRACERID_MOD,         ONLY : IDTSO2
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!      
!#     include "CMN_SIZE"             ! Size parameters
!                                    
!      ! Arguments                   
!      INTEGER, INTENT(IN)           :: THISMONTH
!                                    
!      ! Local variables             
!      INTEGER                       :: I, J, THISYEAR
!      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                        :: BIOCO(IIPAR,JJPAR)
!!-- prior 12/23/08            
!!      TYPE (XPLEX)                        :: CONV, XTAU
!      TYPE (XPLEX)                        :: XTAU
!      CHARACTER(LEN=4  )            :: CYEAR
!      CHARACTER(LEN=255)            :: FILENAME
!                                    
!      ! Days per month              
!      TYPE (XPLEX) :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0, 31d0, 30d0, 
!     &                         31d0, 31d0, 30d0, 31d0, 30d0, 31d0/)
!
!      !=================================================================
!      ! READ_BIOFUEL_SO2 begins here!
!      !=================================================================
!
!      !=================================================================
!      ! Compute biofuel SO2 from biofuel CO.  Use a molar 
!      ! ratio of 0.0015 moles SO2/mole CO from biofuel burning. 
!      ! (Table 2, [Andreae and Merlet, 2001])
!      !=================================================================
!      
!      ! File name for biofuel burning file
!      FILENAME = TRIM( DATA_DIR )          // 
!     &           'biofuel_200202/biofuel.' // GET_NAME_EXT_2D() //
!     &           '.'                       // GET_RES_EXT()
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_BIOFUEL_SO2: Reading ', a )
!
!      ! Get TAU0 value (use generic year 1985)
!      XTAU = GET_TAU0( 1, 1, 1985 )
!
!      ! Read Biofuel burning of CO [kg/yr]
!      CALL READ_BPCH2( FILENAME, 'BIOFSRCE',    4, 
!     &                 XTAU,      IGLOB,        JGLOB,     
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE.  ) 
!
!      ! Cast from TYPE (XPLEX) to COMPLEX*16
!      CALL TRANSFER_2D( ARRAY(:,:,1), BIOCO )
!
!      !=================================================================
!      ! Unit conversion to [kg SO2/s]
!      !=================================================================
!  
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!----------------------------------------------------------------------
!! Prior to 1/28/09:
!! Need to remove CONV from the PRIVATE statement (bmy, 1/28/09)
!!!$OMP+PRIVATE( I, J, CONV )
!!----------------------------------------------------------------------
!!$OMP+PRIVATE( I, J )
!      ! Loop over longitudes
!      DO J = 1, JJPAR
!
!!-- prior 12/23/08            
!!         ! Conversion factor for [cm2 * kg/molec]
!!         CONV = GET_AREA_CM2( J ) / XNUMOL(IDTSO2)
!
!         ! Loop over latitudes
!         DO I = 1, IIPAR
!
!!-- prior 12/23/08            
!!            ! Convert biomass SO2 from [molec SO2/cm2/s] -> [kg SO2/s] 
!!            ! NOTE: Future scale has already been applied by this point
!!            IF ( LBIOMASS ) THEN
!!               ESO2_bb(I,J) = BIOMASS(I,J,IDBSO2) * CONV
!!            ELSE
!!               ESO2_bb(I,J) = 0d0
!!            ENDIF
!
!            ! Convert biofuel SO2 from [kg CO/yr] to [kg SO2/s]
!            ESO2_bf(I,J) = ( BIOCO(I,J) * 64d-3 * 0.0015d0 /
!     &                     ( 28d-3 * 86400.d0 * 365.25d0 ) )
!
!            ! Apply future emissions to biofuel SO2, if necessary
!            IF ( LFUTURE ) THEN
!               ESO2_bf(I,J) = ESO2_bf(I,J) * 
!     &                        GET_FUTURE_SCALE_SO2bf( I, J )
!            ENDIF
!         ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE READ_BIOFUEL_SO2
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE GET_BIOMASS_SO2
!!
!!******************************************************************************
!!  Subroutine GET_BIOMASS_SO2 retrieve monthly/8-day/3hr biomass burning 
!!  emissions for SO2.  (yc, phs, 12/23/08)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) 
!!
!!  NOTES:
!!  (1 ) Extracted from old module subroutine READ_BIOMASS_SO2
!!        (yc, phs, 12/23/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBSO2
!      USE GRID_MOD,             ONLY : GET_AREA_CM2
!      USE TRACER_MOD,           ONLY : XNUMOL
!      USE TRACERID_MOD,         ONLY : IDTSO2
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!      
!#     include "CMN_SIZE"             ! Size parameters
!                                    
!      ! Local variables             
!      INTEGER                       :: I, J
!      TYPE (XPLEX)                        :: CONV
!                                    
!      !=================================================================
!      ! GET_BIOMASS_SO2 begins here!
!      !=================================================================
!      ! Unit conversion to [kg SO2/s]
!  
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, CONV )
!
!      ! Loop over longitudes
!      DO J = 1, JJPAR
!
!         ! Conversion factor for [cm2 * kg/molec]
!         CONV = GET_AREA_CM2( J ) / XNUMOL(IDTSO2)
!
!         ! Loop over latitudes
!         DO I = 1, IIPAR
!
!            ! Convert biomass SO2 from [molec SO2/cm2/s] -> [kg SO2/s] 
!            ! NOTE: Future scale has already been applied by this point
!            ESO2_bb(I,J) = BIOMASS(I,J,IDBSO2) * CONV
!
!         ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE GET_BIOMASS_SO2
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_AIRCRAFT_SO2( THISMONTH )
!!
!!******************************************************************************
!!  Subroutine READ_AIRCRAFT_SO2 reads monthly mean aircraft fuel emissions 
!!  and converts them to SO2 emissions. (rjp, bdf, bmy, 9/18/02, 10/3/05)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!!        GCAP and GEOS grids (bmy, 8/16/05)
!!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      USE BPCH2_MOD,     ONLY : GET_RES_EXT, GET_TAU0, READ_BPCH2
!      USE DAO_MOD,       ONLY : BXHEIGHT
!      USE DIRECTORY_MOD, ONLY : DATA_DIR 
!      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: THISMONTH
!
!      ! Local variables
!      INTEGER             :: I, IOS, J, K, L
!      TYPE (XPLEX)              :: ACSO2(IGLOB,JGLOB,20)
!      TYPE (XPLEX)              :: FAC, FUEL, DZ(LLPAR), ZH(0:LLPAR)
!      CHARACTER(LEN=255)  :: FILENAME
!
!      ! Month names
!      CHARACTER(LEN=3)    :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr', 
!     &                                      'may', 'jun', 'jul', 'aug',
!     &                                      'sep', 'oct', 'nov', 'dec'/)
!
!      !=================================================================
!      ! READ_AIRCRAFT_SO2 begins here!
!      !=================================================================
!      
!      ! Zero arrays
!      ESO2_ac = 0d0
!      ACSO2   = 0d0
!      
!      ! File name
!      FILENAME = TRIM( DATA_DIR )               // 
!     &           'sulfate_sim_200508/aircraft.' // GET_RES_EXT() //
!     &           '.1992.'                       // CMONTH(THISMONTH)
!
!      ! Echo output
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_AIRCRAFT_SO2: Reading ', a )     
!
!      !=================================================================
!      ! Read aircraft emissions.  These are fuel burned in [kg/box/day],
!      ! from AEAP for 1992.  SO2 emission is calculated by assuming    
!      ! an emission index EI of 1.0, i.e., 1g of SO2 emitted per kg    
!      ! of fuel burned.  It is also assumed that there is no diurnal   
!      ! variation of emission rate. Convert to [kg SO2/box/s]. 
!      !=================================================================
!
!      ! Open file 
!      OPEN( IU_FILE, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
!      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_aircraft_so2:1' )
!
!      ! Read header line
!      READ( IU_FILE, '(/)', IOSTAT=IOS )
!      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_aircraft_so2:2' )
!      
!      ! Read data values until an EOF is found
!      DO 
!         READ( IU_FILE, '(3i4,e11.3)', IOSTAT=IOS ) I, J, L, FUEL
!
!         ! EOF encountered
!         IF ( IOS < 0 ) EXIT
!
!         ! I/O error condition
!         IF ( IOS > 0 ) THEN
!            CALL IOERROR( IOS, IU_FILE, 'read_aircraft_so2:3' )
!         ENDIF
!
!         ! Unit conversion: [kg Fuel/box/day] -> [kg SO2/box/s]
!         ! Assuming an emission index of 1.0, 
!         ! 1 g SO2 / kg fuel burned [Weisenstein et al., 1996]
!         ACSO2(I,J,L+1) = 1.d-3 * FUEL / ( 24.d0 * 3600d0 )
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_FILE )
!
!      !=================================================================
!      ! Interpolate from the 1-km grid to the given GEOS-CHEM grid
!      ! NOTE: we need to account for window grids (bmy, 9/20/02)
!      !=================================================================
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! ACSO2 is the aircraft SO2 on the 1-km vertical grid
!         FUEL = SUM( ACSO2(I,J,:) )
!         IF ( FUEL < 1d-20 ) CYCLE
!
!         ! There are 20 1-km levels
!         DO K = 1, 20
!
!            ! Initialize
!            ZH(0) = 0.d0
!
!            ! Loop over levels
!            DO L = 1, LLPAR
!
!               ! Altitude of top edge of level L, from ground [km]
!               ZH(L) = ZH(L-1) + ( BXHEIGHT(I,J,L) * 1d-3 )
!               
!               IF ( ZH(L-1) > DBLE(K)   ) EXIT
!               IF ( ZH(L  ) < DBLE(K-1) ) CYCLE
!               
!               IF ( ZH(L) < DBLE(K) ) THEN
!                  FAC            = ZH(L) - MAX( ZH(L-1), DBLE(K-1) )
!                  ESO2_ac(I,J,L) = ESO2_ac(I,J,L) + ACSO2(I,J,K) * FAC
!               ELSE
!                  FAC            = DBLE(K) - MAX( ZH(L-1), DBLE(K-1) )
!                  ESO2_ac(I,J,L) = ESO2_ac(I,J,L) + ACSO2(I,J,K) * FAC
!                  EXIT
!               ENDIF		     
!            ENDDO
!         ENDDO     
!      ENDDO
!      ENDDO
!
!      ! Return to calling program
!      END SUBROUTINE READ_AIRCRAFT_SO2
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_SHIP_SO2( THISMONTH )
!!
!!******************************************************************************
!!  Subroutine READ_SHIP_SO2 reads in ship SO2 emissions, from either Corbett
!!  et al or EDGAR inventories. (bec, qli, 10/01/03, 7/14/06)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!!
!!  NOTES:
!!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (2 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!!        GCAP and GEOS grids. (bmy, 8/16/05)
!!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!!  (5 ) Now get EDGAR ship SO2 emissions if necessary.  Also apply future
!!        emissions scale factors to the default Corbett et al ship emissions.
!!        (avd, bmy, 7/14/06)
!!  (6 ) Now references GET_ARCTAS_HIP from 'arctas_ship_emiss_mod.f" and
!!        GET_EMEP_ANTHRO to get ARCTAS and EMEP SO2 ship emissions (phs, 12/5/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE ARCTAS_SHIP_EMISS_MOD,ONLY : GET_ARCTAS_SHIP
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR 
!      USE EDGAR_MOD,            ONLY : GET_EDGAR_SHIP_SO2
!      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO, GET_EUROPE_MASK
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
!      USE GRID_MOD,             ONLY : GET_AREA_CM2
!      USE LOGICAL_MOD,          ONLY : LEDGARSHIP,      LFUTURE, 
!     &                                 LARCSHIP,        LSHIPSO2,
!     $                                 LEMEPSHIP
!      USE TRACER_MOD,           ONLY : XNUMOL
!      USE TRACERID_MOD,         ONLY : IDTSO2
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"             ! Size parameters 
!
!      ! Arguments
!      INTEGER, INTENT(IN)           :: THISMONTH
!
!      ! Local variables
!      INTEGER                       :: I, J
!      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                        :: SHIPSO2(IIPAR,JJPAR)
!      TYPE (XPLEX)                        :: XTAU, AREA_CM2
!      CHARACTER (LEN=255)           :: FILENAME
!
!      !=================================================================
!      ! READ_SHIP_SO2 begins here!
!      !=================================================================
!
!      ! Reset
!      ESO2_sh = 0D0
!
!
!      ! Test for EDAGR last, since this is default inventory by design.
!      ! So we can still use EDGAR SHIP to get ship-NOX and CO, and
!      ! overwrite ship-SO2 with ARCTAS or Colbert (phs, 12/5/08)
!
!      !-----------------------------------------------------------
!      ! Use ARCTAS SHIP emissions (EDGAR 2006 update) 
!      !-----------------------------------------------------------
!      IF ( LARCSHIP ) THEN
!            
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!         ! Read ARCTAS SO2 emissions in [kg SO2/BOX/s]
!         ESO2_sh(I,J) = GET_ARCTAS_SHIP( I, J, IDTSO2,  KG_S=.TRUE. )
!               
!            IF ( LEMEPSHIP ) THEN
!               IF ( GET_EUROPE_MASK(I,J) > 0d0 )
!     $           ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2, 
!     &                                          KG_S=.TRUE.,
!     $                                          SHIP=.TRUE.)
!            ENDIF
!
!         ENDDO 
!         ENDDO
!!$OMP END PARALLEL DO
!
!      !----------------------------------------
!      ! Or Corbett et al ship SO2 emissions
!      !----------------------------------------
!      ELSE IF ( LSHIPSO2 ) THEN
!
!         ! Filename
!         FILENAME = TRIM( DATA_DIR )           // 
!     &           'sulfate_sim_200508/shipSOx.' // GET_NAME_EXT_2D() //
!     &           '.'                           // GET_RES_EXT()
!
!         ! Echo some information to the standard output
!         WRITE( 6, 110 ) TRIM( FILENAME )
! 110     FORMAT( '     - READ_SHIP_SO2 ', a )
!      
!         ! TAU value at the beginning of this month
!         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
!      
!         ! Read in this month's ship SO2 emissions [molec SO2/cm2/s]
!         CALL READ_BPCH2( FILENAME, 'SO2-SHIP',     26,  
!     &                    XTAU,      IIPAR,         JJPAR, 
!     &                    1,         ARRAY(:,:,1),  QUIET=.TRUE. )
!
!         ! Cast from TYPE (XPLEX) to COMPLEX*16
!         CALL TRANSFER_2D( ARRAY(:,:,1), SHIPSO2 )
!
!         ! Loop over latitudes
!         DO J = 1, JJPAR
!
!            ! Grid box surface area [cm2]
!            AREA_CM2 = GET_AREA_CM2( J )
!
!            ! Loop over longitudes
!            DO I = 1, IIPAR
!            
!               ! Convert [molec SO2/cm2/s] to [kg SO2/box/s]
!               ESO2_sh(I,J) = SHIPSO2(I,J) * AREA_CM2 / XNUMOL(IDTSO2)
!
!               ! Apply future emissions (if necessary)
!               IF ( LFUTURE ) THEN
!                   ESO2_sh(I,J) = ESO2_sh(I,J) *
!     &                            GET_FUTURE_SCALE_SO2ff( I, J ) 
!               ENDIF
!
!               IF ( LEMEPSHIP ) THEN
!                  IF ( GET_EUROPE_MASK(I,J) > 0d0 )
!     $                 ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2,
!     $                                         KG_S=.TRUE., SHIP=.TRUE.)
!               ENDIF               
!
!            ENDDO 
!         ENDDO
!
!      !-----------------------------------------------------------
!      ! Test for EDGAR ship emissions
!      !-----------------------------------------------------------
!      ELSE IF ( LEDGARSHIP ) THEN 
!
!         !----------------------------------------
!         ! Use EDGAR ship SO2 emissions
!         !----------------------------------------
!
!         ! Get EDGAR ship SO2 [kg SO2/box/s]
!         ! NOTE: Future emissions have already been applied! 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ESO2_sh(I,J) = GET_EDGAR_SHIP_SO2( I, J, KG_S=.TRUE. )
!
!            IF ( LEMEPSHIP ) THEN
!               IF ( GET_EUROPE_MASK(I,J) > 0d0 )
!     $           ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2, 
!     &                                          KG_S=.TRUE.,
!     $                                          SHIP=.TRUE.)
!            ENDIF
!         ENDDO 
!         ENDDO
!!$OMP END PARALLEL DO
!         
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE READ_SHIP_SO2
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_ANTHRO_NH3( THISMONTH )
!!
!!******************************************************************************
!!  Subroutine READ_ANTHRO_NH3 reads the monthly mean anthropogenic 
!!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!!  (rjp, bdf, bmy, 9/20/02, 10/31/08)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Renamed from NH3_READ to READ_ANTHRO_NH3.  Also updated comments,
!!        made cosmetic changes. (bmy, 9/20/02)
!!  (2 ) Changed filename to NH3_anthsrce.geos.*.  Also now reads data under
!!        category name "NH3-ANTH". (rjp, bmy, 3/23/03)
!!  (3 ) Now reads from NH3emis.monthly.geos.* files.  Now call READ_BPCH2
!!        with QUIET=.TRUE. (bmy, 3/27/03)
!!  (4 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (5 ) Now read files from "sulfate_sim_200508/". Now read data for both 
!!        GCAP and GEOS grids. (bmy, 8/16/05)
!!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (6 ) Now compute future emissions, if necessary (swu, bmy, 5/30/06)
!!  (7 ) Now overwrite w/ David Streets' NH3, if necessary (yxw, bmy, 8/17/06)
!!  (8 ) Bug fix: Using tracer #30 in the call to GET_STREETS_ANTHRO can cause
!!        problems when adding or removing species.  Replace w/ IDTNH3.
!!        (dkh, 10/31/08)
!!  (9 ) Now check if NH3 Streets is available (phs, 12/10/08)      
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR
!      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO
!      USe EMEP_MOD,             ONLY : GET_EUROPE_MASK
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an
!      USE LOGICAL_MOD,          ONLY : LFUTURE,         LSTREETS
!      USE LOGICAL_MOD,          ONLY : LEMEP
!      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
!      USE STREETS_ANTHRO_MOD,   ONLY : GET_STREETS_ANTHRO
!      USE TRACERID_MOD,         ONLY : IDTNH3
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"             ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)           :: THISMONTH
!
!      ! Local variables
!      LOGICAL                       :: WEEKDAY
!      INTEGER                       :: I, J, DAY_NUM
!      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                        :: AREA_CM2, EPA_NEI, XTAU
!      TYPE (XPLEX)              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
!     &                                      31d0, 30d0, 31d0, 31d0, 
!     &                                      30d0, 31d0, 30d0, 31d0 /)
!      CHARACTER(LEN=255)            :: FILENAME
!      TYPE (XPLEX)                        :: STREETS
!
!      !=================================================================
!      ! READ_ANTHRO_NH3 begins here!
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( DATA_DIR )                         //
!     &           'sulfate_sim_200508/NH3_anthsrce.'       //
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!
!      ! Echo output
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_ANTHRO_NH3: Reading ', a )
!      
!      ! Get TAU value (use year 1990, the year of the data!)
!      XTAU = GET_TAU0( THISMONTH, 1, 1990 )
!	
!      ! Read 1990 NH3 emissions [kg N/box/mon]
!      CALL READ_BPCH2( FILENAME, 'NH3-ANTH',    29,  
!     &                 XTAU,      IGLOB,        JGLOB,       
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
!
!      ! Cast from TYPE (XPLEX) to COMPLEX*16
!      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_an )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Convert from [kg N/box/mon] to [kg NH3/box/s]
!         ENH3_an(I,J) = ENH3_an(I,J) * ( 17.d0 / 14.d0 ) 
!     &                / ( NMDAY(THISMONTH) * 86400.d0 )
!
!         ! Compute future NH3an emissions, if necessary
!         ! Moved here since Streets and EMEP should have already
!         ! applied FUTURE scale factors if needed 
!         IF ( LFUTURE ) THEN
!            ENH3_an(I,J) = ENH3_an(I,J) * GET_FUTURE_SCALE_NH3an( I, J )
!         ENDIF
!
!         ! If we are using David Streets' emissions ...
!         IF ( LSTREETS ) THEN
!
!            ! If we are over the SE Asia region ...
!            IF ( GET_SE_ASIA_MASK( I, J ) > 0d0  ) THEN
!
!               ! Overwrite with David Streets emissions [kg NH3/s]
!!------------------------------------------------------------------------------
!! Prior to 12/10/08:
!! Now check first that NH3 is available (phs, 12/10/08)               
!!               ENH3_an(I,J) = GET_STREETS_ANTHRO( I,      J, 
!!     &                                            IDTNH3, KG_S=.TRUE.)
!               STREETS = GET_STREETS_ANTHRO( I,      J, 
!     &                                       IDTNH3, KG_S=.TRUE.)
!
!               IF ( .not. ( STREETS < 0d0 ) )
!     $              ENH3_an(I,J) = STREETS
!               
!            ENDIF
!         ENDIF
!
!         IF ( LEMEP ) THEN
!            IF ( GET_EUROPE_MASK(I,J) > 0d0) THEN
!               ENH3_an(I,J) = GET_EMEP_ANTHRO(I,J,IDTNH3,KG_S=.TRUE.)
!            ENDIF
!         ENDIF
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE READ_ANTHRO_NH3
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_NATURAL_NH3( THISMONTH )
!!
!!******************************************************************************
!!  Subroutine READ_NATURAL_NH3 reads the monthly mean natural 
!!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!!  (rjp, bdf, bmy, 9/20/02, 10/3/05)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Updated FORMAT string.  Now also call READ_BPCH2 with QUIET=.TRUE.
!!        (bmy, 4/8/03)
!!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!!        GCAP and GEOS grids. (bmy, 8/16/05)
!!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: THISMONTH
!
!      ! Local variables
!      TYPE (XPLEX)              :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)              :: XTAU
!      TYPE (XPLEX)              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
!     &                                      31d0, 30d0, 31d0, 31d0, 
!     &                                      30d0, 31d0, 30d0, 31d0 /)
!      CHARACTER(LEN=255)  :: FILENAME
!
!      !=================================================================
!      ! READ_NATURAL_NH3 begins here!
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( DATA_DIR )                         //
!     &           'sulfate_sim_200508/NH3_natusrce.'       //
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!
!      ! Echo output
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_NATURAL_NH3: Reading ', a )
!      
!      ! Get TAU value (use year 1990, the year of the data!)
!      XTAU = GET_TAU0( THISMONTH, 1, 1990 )
!	
!      ! Read 1990 NH3 emissions [kg N/box/mon]
!      CALL READ_BPCH2( FILENAME, 'NH3-NATU',    29,  
!     &                 XTAU,      IGLOB,        JGLOB,       
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
!
!      ! Cast from TYPE (XPLEX) to COMPLEX*16
!      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_na )
!
!      ! Convert from [kg N/box/mon] to [kg NH3/box/s]
!      ENH3_na = ENH3_na * ( 17.d0 / 14.d0 ) /
!     &          ( NMDAY(THISMONTH) * 86400.d0 ) 
! 
!      ! Return to calling program
!      END SUBROUTINE READ_NATURAL_NH3
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_BIOFUEL_NH3( THISMONTH ) 
!!
!!******************************************************************************
!!  Subroutine READ_BIOFUEL_NH3 reads the monthly mean biomass NH3 
!!  and biofuel emissions from disk and converts to [kg NH3/box/s]. 
!!  (rjp, bdf, bmy, phs, 9/20/02, 12/23/08)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Renamed from NH3_READ to READ_BIOMASS_NH3.  Also updated comments,
!!        made cosmetic changes.  Now reads in both biomass and biofuel
!!        emissions. (rjp, bmy, 12/13/02)
!!  (2 ) Now replace DXYP(J+J0) with routine GET_AREA_M2 of "grid_mod.f"
!!        Now use function GET_YEAR from "time_mod.f".  Replace MONTH with 
!!        THISMONTH when referencing the NMDAY variable.  Now call READ_BPCH2
!!        with QUIET=.TRUE. (bmy, 3/27/03)
!!  (3 ) If using interannual biomass emissions, substitute seasonal emissions
!!        for years where internannual emissions do not exist.  Now also
!!        reference GET_TAU from "time_mod.f" (bmy, 5/15/03)
!!  (4 ) Now use ENCODE statement for PGI/F90 on Linux (bmy, 9/29/03)
!!  (5 ) Changed cpp switch name from LINUX to LINUX_PGI (bmy, 12/2/03)
!!  (6 ) Now references DATA_DIR from "directory_mod.f".  Now references LBBSEA
!!        from "logical_mod.f".  Removed references to CMN and CMN_SETUP.
!!        (bmy, 7/20/04)
!!  (7 ) Now can read either seasonal or interannual biomass burning emissions.
!!        Now references routines from both and "time_mod.f".  Now reads SO2 
!!        biomass emissions directly rather than computing it by mole fraction 
!!        from CO. (rjp, bmy, 1/11/05)
!!  (8 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!!        GCAP and GEOS grids. (bmy, 8/16/05)
!!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (10) Now only read the biofuel, we have moved the biomass-reading code to 
!!        "gc_biomass_mod.f" for compatibility with GFED2 biomass emissions
!!        (bmy, 9/27/06)
!!  (11) Prevent seg fault error when LBIOMASS=F (bmy, 11/3/06)
!!  (12) Renamed READ_BIOFUEL_NH3, and move all biomass code to GET_BIOMASS_NH3
!!        to account for several GFED2 products (yc, phs, 12/23/08)   
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBNH3
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3bf
!      USE LOGICAL_MOD,          ONLY : LBIOMASS,        LFUTURE
!      USE GRID_MOD,             ONLY : GET_AREA_CM2
!      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
!      USE TRACER_MOD,           ONLY : XNUMOL
!      USE TRACERID_MOD,         ONLY : IDTNH3
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"             ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)           :: THISMONTH
!
!      ! Local variables
!      INTEGER                       :: I, J, THISYEAR
!      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                        :: XTAU, DMASS!, CONV
!      TYPE (XPLEX)                        :: NMDAY(12) = (/31d0, 28d0, 31d0, 
!     &                                               30d0, 31d0, 30d0, 
!     &                                               31d0, 31d0, 30d0, 
!     &                                               31d0, 30d0, 31d0/)
!      CHARACTER(LEN=4  )            :: CYEAR
!      CHARACTER(LEN=255)            :: FILENAME
!
!      !=================================================================
!      ! READ_BIOFUEL_NH3 begins here!
!      !=================================================================
!
!      !=================================================================
!      ! Read NH3 biofuel emissions
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'sulfate_sim_200508/NH3_biofuel.'        // 
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!   
!      ! Echo output
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_BIOFUEL_NH3: Reading ', a )
!
!      ! Get TAU0 value for 1998
!      XTAU  = GET_TAU0( THISMONTH, 1, 1998 )
!
!      ! Read NH3 biofuel data [kg NH3/box/month]
!      CALL READ_BPCH2( FILENAME, 'BIOFSRCE',    29, 
!     &                 XTAU,      IGLOB,        JGLOB,       
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
! 
!      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize if necesary
!      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_bf )
!
!      ! Store NH3 in ENH3_bf array [kg NH3/box/s]
!      ENH3_bf = ENH3_bf / ( NMDAY(THISMONTH) * 86400.d0 )
!
!      !=================================================================
!      ! Convert units and apply IPCC future emissions (if necessary)
!      !=================================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!------------------------------------------------------------------
!! Prior to 1/28/09:
!! Need to remove CONV from the PRIVATE statement (bmy, 1/28/09)
!!!$OMP+PRIVATE( I, J, CONV )
!!------------------------------------------------------------------
!!$OMP+PRIVATE( I, J )
!
!      ! Loop over latitudes
!      DO J = 1, JJPAR
! 
!!-- prior 12/23/08            
!!         ! Conversion factor for [cm2 * kg/molec]
!!         CONV = GET_AREA_CM2( J ) / XNUMOL(IDTNH3)
!
!         ! Loop over longitudes
!         DO I = 1, IIPAR
!
!!-- prior 12/23/08            
!!            ! Convert biomass NH3 from [molec NH3/cm2/s] -> [kg NH3/s]
!!            ! NOTE: Future scale is applied by this point (if necessary)
!!            IF ( LBIOMASS ) THEN
!!               ENH3_bb(I,J) = BIOMASS(I,J,IDBNH3) * CONV
!!            ELSE
!!               ENH3_bb(I,J) = 0d0
!!            ENDIF
!
!            ! Scale biofuel NH3 to IPCC future scenario (if necessary)
!            IF ( LFUTURE ) THEN
!               ENH3_bf(I,J) = ENH3_bf(I,J) * 
!     &                        GET_FUTURE_SCALE_NH3bf( I, J )
!            ENDIF
!         ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!      
!      ! Return to calling program
!      END SUBROUTINE READ_BIOFUEL_NH3
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE GET_BIOMASS_NH3
!!
!!******************************************************************************
!!  Subroutine GET_BIOMASS_NH3 retrieve the monthly/8days/3hr mean biomass NH3 
!!  (yc, phs, 12/23/08)
!!
!!  Arguments as input:
!!  ===========================================================================
!!  (1 ) 
!!
!!  NOTES:
!!  (1 ) Extracted from old module subroutine READ_BIOMASS_NH3
!!        (yc, phs, 12/23/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBNH3
!      USE GRID_MOD,             ONLY : GET_AREA_CM2
!      USE TRACER_MOD,           ONLY : XNUMOL
!      USE TRACERID_MOD,         ONLY : IDTNH3
!
!#     include "CMN_SIZE"             ! Size parameters
!
!      ! Local variables
!      INTEGER                       :: I, J
!      TYPE (XPLEX)                        :: CONV
!      
!      !=================================================================
!      ! READ_BIOMASSBURN_NH3 begins here!
!      !=================================================================
!      ! Convert units
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, CONV )
!
!      ! Loop over latitudes
!      DO J = 1, JJPAR
! 
!         ! Conversion factor for [cm2 * kg/molec]
!         CONV = GET_AREA_CM2( J ) / XNUMOL(IDTNH3)
!
!         ! Loop over longitudes
!         DO I = 1, IIPAR
!
!            ! Convert biomass NH3 from [molec NH3/cm2/s] -> [kg NH3/s]
!            ! NOTE: Future scale is applied by this point (if necessary)
!            ENH3_bb(I,J) = BIOMASS(I,J,IDBNH3) * CONV
!
!         ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!      
!      ! Return to calling program
!      END SUBROUTINE GET_BIOMASS_NH3
!      
!!------------------------------------------------------------------------------
!
!      SUBROUTINE READ_OXIDANT( MONTH )
!!
!!******************************************************************************
!!  Subroutine READ_OXIDANT reads in monthly mean H2O2 and O3 fields for the
!!  offline sulfate + aerosol simulation. (rjp, bdf, bmy, 11/1/02, 10/3/05)
!!
!!  Arguments as input:
!!  ============================================================================
!!  (1 ) MONTH    (INTEGER  ) : Emission timestep in minutes
!!
!!  NOTES:
!!  (1 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (3 ) Now read files from "sulfate_sim_200508/offline/".  Now read data
!!        for both GEOS and GCAP grids (bmy, 8/16/05)
!!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!******************************************************************************
!!  
!      ! References to F90 modules
!      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE TRANSFER_MOD,  ONLY : TRANSFER_3D_TROP
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: MONTH
!
!      ! Local variables 
!      INTEGER             :: I, J, L, K      
!      TYPE (XPLEX)              :: ARRAY(IGLOB,JGLOB,LLTROP)
!      TYPE (XPLEX)              :: XTAU
!      CHARACTER(LEN=255)  :: FILENAME
!
!      !=================================================================
!      ! READ_OXIDANT begins here !
!      !
!      ! Oxidant fields were computed for 1998 using coupled aerosol
!      ! and gas chemistry GEOS-CHEM by Brendan Field (bdf, 5/23/02).  
!      ! Bob Yantosca has regridded these fields to all GEOS-CHEM grids.  
!      ! Data is saved from the surface to the tropopause. 
!      !=================================================================
!
!      ! Use generic year 1985
!      XTAU = GET_TAU0( MONTH, 1, 1985 )
!
!      !=================================================================
!      ! Read monthly mean PH2O2 (from HO2 + HO2 = H2O2) [molec/cm3/s]
!      !=================================================================
!      FILENAME = TRIM( DATA_DIR )                      // 
!     &           'sulfate_sim_200508/offline/PH2O2.'   // 
!     &           GET_NAME_EXT() //  '.' // GET_RES_EXT()
!
!      ! Echo filename
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_OXIDANT: Reading ', a ) 
!
!      ! Read data
!      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 5,     
!      ! limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
!      ! tropopause. This is backward compatibility with 
!      ! offline data set.
!     &     XTAU,        IGLOB,                    JGLOB,      
!     &     LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )
!!     &                 XTAU,      IGLOB,     JGLOB,     
!!     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )
!
!      ! Cast to TYPE (XPLEX) and resize if necessary
!      CALL TRANSFER_3D_TROP( ARRAY, PH2O2m )
!            
!      !=================================================================
!      ! Read monthly mean O3 [v/v]
!      !=================================================================
!      FILENAME = TRIM( DATA_DIR )                      // 
!     &           'sulfate_sim_200508/offline/O3.'      //
!     &           GET_NAME_EXT() // '.' // GET_RES_EXT()
!
!      ! Echo filename
!      WRITE( 6, 100 ) TRIM( FILENAME )
!
!      ! Read data
!      ! limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
!      ! tropopause. This is backward compatibility with 
!      ! offline data set.
!      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 51,     
!     &     XTAU,        IGLOB,                    JGLOB,      
!     &     LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )
!!     &                 XTAU,      IGLOB,     JGLOB,     
!!     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )
!
!      ! Cast to TYPE (XPLEX) and resize if necessary
!      CALL TRANSFER_3D_TROP( ARRAY, O3m ) 
!      
!      ! Return to calling program
!      END SUBROUTINE READ_OXIDANT
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE OHNO3TIME
!!
!!******************************************************************************
!!  Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!!  angle over a 24 hour day, as well as the total length of daylight. 
!!  This is needed to scale the offline OH and NO3 concentrations.
!!  (rjp, bmy, 12/16/02, 3/30/04)
!!
!!  NOTES:
!!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!!        values.  Figure this out later (rjp, bmy, 1/10/03)
!!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!!        "time_mod.f". (bmy, 3/27/03)
!!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!!        the COSZM array. (rjp, bmy, 3/30/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
!      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
!      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
!#     include "CMN_SIZE"  ! Size parameters
!#     include "CMN_GCTM"
!
!      ! Local variables
!      LOGICAL, SAVE       :: FIRST = .TRUE.
!      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
!      TYPE (XPLEX)              :: A0, A1, A2, A3, B1, B2, B3
!      TYPE (XPLEX)              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
!      TYPE (XPLEX)              :: SUNTMP(MAXIJ)
!      
!      !=================================================================
!      ! OHNO3TIME begins here!
!      !=================================================================
!
!      !  Solar declination angle (low precision formula, good enough for us):
!      A0 = 0.006918
!      A1 = 0.399912
!      A2 = 0.006758
!      A3 = 0.002697
!      B1 = 0.070257
!      B2 = 0.000907
!      B3 = 0.000148
!      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.
!
!      DEC = A0 - A1*cos(  R) + B1*sin(  R)
!     &         - A2*cos(2*R) + B2*sin(2*R)
!     &         - A3*cos(3*R) + B3*sin(3*R)
!
!      LHR0 = int(float( GET_NHMSb() )/10000.)
!
!      ! Only do the following at the start of a new day
!      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
!      
!         ! Zero arrays
!         TTDAY(:,:) = 0d0
!         TCOSZ(:,:) = 0d0
!         COSZM(:,:) = 0d0
!
!         ! NDYSTEP is # of chemistry time steps in this day
!         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         
!
!         ! NT is the elapsed time [s] since the beginning of the run
!         NT = GET_ELAPSED_SEC()
!
!         ! Loop forward through NDYSTEP "fake" timesteps for this day 
!         DO N = 1, NDYSTEP
!            
!            ! Zero SUNTMP array
!            SUNTMP(:) = 0d0
!
!            ! IJLOOP is the 1-D loop index for SUNCOS
!            IJLOOP = 0
!
!            ! Loop over surface grid boxes
!            DO J = 1, JJPAR
!
!               ! Grid box latitude center [radians]
!               YMID_R = GET_YMID_R( J )
!
!            DO I = 1, IIPAR
!
!               ! Increment IJLOOP
!               IJLOOP = IJLOOP + 1
!               TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
!         
!               DO WHILE (TIMLOC .lt. 0)
!                  TIMLOC = TIMLOC + 24.0
!               ENDDO
!
!               DO WHILE (TIMLOC .gt. 24.0)
!                  TIMLOC = TIMLOC - 24.0
!               ENDDO
!
!               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180
!
!            !===========================================================
!            ! The cosine of the solar zenith angle (SZA) is given by:
!            !     
!            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
!            !                   
!            ! where LAT = the latitude angle, 
!            !       DEC = the solar declination angle,  
!            !       AHR = the hour angle, all in radians. 
!            !
!            ! If SUNCOS < 0, then the sun is below the horizon, and 
!            ! therefore does not contribute to any solar heating.  
!            !===========================================================
!
!               ! Compute Cos(SZA)
!               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
!     &                          cos(YMID_R) * cos(DEC) * cos(AHR)
!
!               ! TCOSZ is the sum of SUNTMP at location (I,J)
!               ! Do not include negative values of SUNTMP
!               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )
!
!               ! COSZM is the peak value of SUMTMP during a day at (I,J)
!               ! (rjp, bmy, 3/30/04)
!               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )
!
!               ! TTDAY is the total daylight time at location (I,J)
!               IF ( SUNTMP(IJLOOP) > 0d0 ) THEN
!                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
!               ENDIF
!            ENDDO
!            ENDDO
!
!            !### Debug
!            !PRINT*, '### IN OHNO3TIME'
!            !PRINT*, '### N       : ', N
!            !PRINT*, '### NDYSTEP : ', NDYSTEP
!            !PRINT*, '### NT      : ', NT
!            !PRINT*, '### JDAY    : ', JDAY
!            !PRINT*, '### RLAT    : ', RLAT
!            !PRINT*, '### XMID    : ', XMID
!            !PRINT*, '### SUNTMP  : ', SUNTMP
!            !PRINT*, '### TCOSZ   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )
!            !PRINT*, '### TTDAY   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )
!
!            ! Increment elapsed time [sec]
!            NT = NT + ( GET_TS_CHEM() * 60 )             
!         ENDDO
!
!         ! Reset first-time flag
!         FIRST = .FALSE.
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE OHNO3TIME
!
!!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE GET_O3_ADJ( ADO3, I, J, L )
!
!******************************************************************************
!  Subroutine GET_O3_ADJ is the adjoint of GET_O3. 
!  (dkh, 10/12/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1  )           (TYPE (XPLEX))  : Local adjoint of ozone conc
!  (2-4) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) The only reason we can just add the new value  of ADO3 to the pre-existing
!        adjoint of ozone is that ADO3 is not involved in any equations in 
!        ADJ_CHEM_SO2 other than ado3 = ado3 + ...  
!        Otherwise we would have to initialize ado3 with ADCSPEC, then replace
!        ADCSPEC directly here. (dkh, 10/24/05) 
!  (2 ) Updated to GCv8 adjoint (dkh, 09/28/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP, VOLUME, CSPEC_ADJ
      USE DAO_MOD,      ONLY : AIRDEN
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,   ONLY : XNUMOLAIR
      USE TRACERID_MOD, ONLY : IDO3

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      TYPE (XPLEX), INTENT(IN)  :: ADO3

      ! Local variables
      INTEGER             :: JLOOP

      !=================================================================
      ! GET_O3_ADJ begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Add additional sensitivity gleaned from sulfate 
         ! chemistry to previous adjoint ozone value. 
         IF ( JLOOP  > 0 ) THEN
            CSPEC_ADJ(JLOOP,IDO3) = CSPEC_ADJ(JLOOP,IDO3) +
     &       ( ADO3 * 1d6 ) / ( AIRDEN(L,I,J) * XNUMOLAIR )
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         CALL ERROR_STOP( 'Invalid NSRCX!',
     &     'GET_O3_ADJ (sulfate_adj_mod.f)')

      ENDIF 


      ! Return to calling program
      END SUBROUTINE GET_O3_ADJ

!------------------------------------------------------------------------------

          
      SUBROUTINE INIT_SULFATE_ADJ
!
!******************************************************************************
!  Subroutine INIT_ADJ_SULFATE initializes and zeros all allocatable arrays
!  declared in "adj_sulfate_mod.f" (dkh, 10/12/05)
!
!  NOTES:
!  (1 ) Updated for GCv8 adjoint (dkh, 09/28/09) 
!
!******************************************************************************
!


      ! Reference to f90 modules 
      USE DRYDEP_MOD,    ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE LOGICAL_MOD,   ONLY : LDRYD

         
#     include "CMN_SIZE"      ! Size parameters
            
      ! Local variables 
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS, I, J, N, IJLOOP

      !=================================================================
      ! INIT_SULFATE_ADJ begins here!
      !=================================================================

      ! Return if we have already initialized
      IF ( IS_INIT ) RETURN

      ALLOCATE( VCLDF( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VCLDF' )
      VCLDF = 0d0
      
      ALLOCATE( ADPSO4_SO2( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ADPSO4_SO2' )
      ADPSO4_SO2 = 0d0


      ! Initialize flags
      DRYH2O2  = 0
      DRYSO2   = 0
      DRYSO4   = 0
      DRYSO4s  = 0
      DRYMSA   = 0
      DRYNH3   = 0
      DRYNH4   = 0
      DRYNIT   = 0
      DRYSO4s  = 0
      DRYAS    = 0
      DRYAHS   = 0
      DRYLET   = 0
      DRYSO4aq = 0
      DRYNH4aq = 0  

      IF ( LDRYD ) THEN
         
         ! Locate position of each tracer in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'H2O2'   )
                  DRYH2O2  = N
               CASE ( 'SO2'   )
                  DRYSO2   = N
               CASE ( 'SO4'   )
                  DRYSO4   = N
               CASE ( 'SO4S'   )
                  DRYSO4s  = N
               CASE ( 'MSA'   )
                  DRYMSA   = N
               CASE ( 'NH3'   )
                  DRYNH3   = N
               CASE ( 'NH4'   )
                  DRYNH4   = N
               CASE ( 'NIT'   )
                  DRYNIT   = N
               CASE ( 'NITS'   )
                  DRYNITs  = N
               CASE ( 'AS'    )
                  DRYAS    = N
               CASE ( 'AHS'   )
                  DRYAHS   = N
               CASE ( 'LET'   )
                  DRYLET   = N
               CASE ( 'SO4aq' )
                  DRYSO4aq = N
               CASE ( 'NH4aq' )
                  DRYNH4aq = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

      ENDIF

      ! Reset IS_INIT so we do not allocate arrays again
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_SULFATE_ADJ

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SULFATE_ADJ
!
!******************************************************************************
!  Subroutine CLEANUP_SULFATE_ADJ deallocates all previously allocated arrays 
!  for sulfate emissions -- call at the end of the run (bmy, 6/1/00, 5/3/06)
! 
!  NOTES:
!  (1 ) Updated for GCv8 adjoint (dkh, 09/28/09) 
!******************************************************************************
! 
      !=================================================================
      ! CLEANUP_SULFATE_ADJ begins here!
      !=================================================================
      IF ( ALLOCATED( ADPSO4_SO2 ) ) DEALLOCATE( ADPSO4_SO2 )
      IF ( ALLOCATED( VCLDF      ) ) DEALLOCATE( VCLDF      )

      ! Return to calling program
      END SUBROUTINE CLEANUP_SULFATE_ADJ

!------------------------------------------------------------------------------

      END MODULE SULFATE_ADJ_MOD
