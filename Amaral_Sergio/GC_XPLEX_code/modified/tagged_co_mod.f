! $Id: tagged_co_mod.f,v 1.5 2012/03/01 22:00:27 daven Exp $
      MODULE TAGGED_CO_MOD
!
!******************************************************************************
!  Module TAGGED_CO_MOD contains variables and routines used for the 
!  geographically tagged CO simulation. (bmy, 7/28/00, 2/14/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BB_REGION   : Array for geographic regions for biomass burning CO
!  (2 ) FF_REGION   : Array for geographic regions for fossil fuel CO
!  (3 ) SUMISOPCO   : Array for production of CO from Isoprene
!  (4 ) SUMMONOCO   : Array for production of CO from Monoterpenes
!  (5 ) SUMCH3OHCO  : Array for production of CO from CH3OH (methanol)
!  (6 ) SUMACETCO   : Array for production of CO from Acetone
!  (7 ) EMACET      : Array for hold monthly mean acetone emissions
!  (8 ) CO_PRODS    : Array for P(CO) from CH4      in the stratosphere
!  (9 ) CO_LOSSS    : Array for L(CO) from CO  + OH in the stratosphere
!  (10) FMOL_CO     : molecular weight of CO
!  (11) XNUMOL_CO   : molec CO / kg CO
!  (12) FMOL_CO     : molecular weight of ISOP
!  (13) XNUMOL_CO   : molec ISOP / kg ISOP
!  (14) FMOL_MONO   : molecular weight of MONOTERPENES
!  (15) XNUMOL_MONO : molec MONOTERPENES / kg MONOTERPENES
!
!  Module Routines:
!  ============================================================================
!  (1 ) DEFINE_FF_CO_REGIONS : Defines geographic regions for fossil fuel CO
!  (2 ) DEFINE_BB_CO_REGIONS : Defines geographic regions for biomass burn. CO
!  (3 ) EMISS_TAGGED_CO      : Emits CO into geographically "tagged" tracers
!  (4 ) CHEM_TAGGED_CO       : Does chemistry for "tagged" CO tracers
!  (5 ) GET_ALPHA_ISOP       : Returns CO yield from isoprene as f(NOx)
!  (6 ) READ_PCO_LCO_STRAT   : Reads data into CO_PRODS and CO_LOSSS
!  (7 ) GET_PCO_LCO_STRAT    : Extracts data from CO_PRODS and CO_LOSSS
!  (8 ) READ_ACETONE         : Reads biog acetone and acetone from grasslands
!  (9 ) INIT_TAGGED_CO       : Allocates and initializes module arrays
!  (10) CLEANUP_TAGGED_CO    : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by tagged_co_mod.f
!  ============================================================================
!  (1 ) biofuel_mod.f    : Module w/ routines to read biofuel emissions
!  (2 ) biomass_mod.f    : Module w/ routines to read biomass emissions
!  (3 ) bpch2_mod.f      : Module w/ routines for binary punch file I/O
!  (4 ) dao_mod.f        : Module w/ arrays for DAO met fields
!  (5 ) diag_mod.f       : Module w/ GEOS-CHEM diagnostic arrays
!  (6 ) diag_pl_mod.f    : Module w/ routines for prod & loss diag's
!  (7 ) directory_mod.f  : Module w/ GEOS-CHEM data & met field dirs
!  (8 ) error_mod.f      : Module w/ I/O error and NaN check routines
!  (9 ) geia_mod         : Module w/ routines to read anthro emissions
!  (10) global_oh_mod.f  : Module w/ routines to read 3-D OH field
!  (11) global_nox_mod.f : Module w/ routines to read 3-D NOx field
!  (12) grid_mod.f       : Module w/ horizontal grid information
!  (13) logical_mod.f    : Module w/ GEOS-CHEM logical switches
!  (14) pressure_mod.f   : Module w/ routines to compute P(I,J,L)
!  (15) time_mod.f       : Module w/ routines for computing time & date
!  (16) tracer_mod.f     : Module w/ GEOS-CHEM tracer array STT etc.
!  (17) tropopause_mod.f : Module w/ routines to read ann mean tropopause
!  (18) logical_adj_mod.f: Module w/ adj logical flags
! 
!  Tagged CO Tracers (you can modify these as needs be!)
!  ============================================================================
!  (1 ) Total CO
!  (2 ) CO from North American fossil fuel 
!  (3 ) CO from European fossil fuel 
!  (4 ) CO from Asian fossil fuel 
!  (5 ) CO from fossil fuel from everywhere else
!  (6 ) CO from South American biomass burning 
!  (7 ) CO from African biomass burning 
!  (8 ) CO from Southeast Asian biomass burning
!  (9 ) CO from Oceania biomass burning
!  (10) CO from European biomass burning
!  (11) CO from North American biomass burning
!  (12) CO chemically produced from Methane
!  (13) CO from Biofuel burning (whole world)
!  (14) CO chemically produced from Isoprene 
!  (15) CO chemically produced from Monoterpenes
!  (16) CO chemically produced from Methanol (CH3OH)
!  (17) CO chemically produced from Acetone
!
!  NOTES:
!  (1 ) Removed obsolete code from CHEM_TAGGED_CO (bmy, 12/21/00)
!  (2 ) Added CO sources from oxidation of biofuel VOC's, biomass burning
!        VOC's, fossil fuel VOC's, and natural VOC's (bnd, bmy, 1/2/01)
!  (3 ) Added chemical P(CO) from CH3OH and MONOTERPENES (bnd, bmy, 1/2/01)
!  (4 ) Now cap SCALEYEAR at 1997 in "emiss_tagged_co" (bnd, bmy, 4/6/01) 
!  (5 ) Removed obsolete commented-out code (bmy, 4/23/01)
!  (6 ) Added new module variables SUMACETCO, EMACET, CO_PRODS, CO_LOSSS,
!        ISOP96, MONO96, and MEOH96.  Also added new module routines 
!        GET_ALPHA_ISOP, READ_PCO_LCO_STRAT, GET_PCO_LCO_STRAT, 
!        READ_ACETONE, and READ_BIOG_FOR_GEOS3. (bnd, bmy, 6/14/01)
!  (7 ) Now read files from DATA_DIR/tagged_CO_200106/ (bmy, 6/19/01) 
!  (8 ) Removed ISOP96, MONO96, and CH3OH96 since we now use the new GEOS-3
!        fields and no longer have to correct for the surface temperature.
!        (bmy, 8/21/01)
!  (9 ) Bug fix: don't call GLOBAL_NOX_MOD in routine CHEM_TAGGED_CO 
!        unless a logical switch is set (bmy, 8/28/01)
!  (10) Updated comments (bmy, 9/6/01)
!  (11) Deleted obsolete code for 1998 GEOS-3 fix.  Also now archive ND46
!        diagnostic as [atoms C/cm2/s] (bmy, 9/13/01)
!  (12) Bug fix in CHEM_TAGGED_CO: now save CO sources/sinks into the
!        Total CO tracer (N=1). (qli, bmy, 9/21/01)
!  (13) Resize arrays of (IGLOB,JGLOB) to (IIPAR,JJPAR) (bmy, 9/28/01)
!  (14) Removed obsolete code from 9/28/01 (bmy, 10/23/01)
!  (15) Updated comments (bmy, 2/15/02)
!  (16) Removed double definition of SUMCH3OHCO (bmy, 3/20/02)
!  (17) Now use P(I,J) + PTOP instead of PS(I,J) (bmy, 4/11/02)
!  (18) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)  
!  (19) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/21/02)
!  (20) Now reference AD, BXHEIGHT, T and SUNCOS from "dao_mod.f".  Also 
!        removed obsolete code from various routines.  Now references
!        ERROR_STOP from "error_mod.f". (bmy, 10/15/02)
!  (21) Now references "grid_mod.f" and the new "time_mod.f". (bmy, 2/3/03)
!  (22) Bug fix for NTAU in EMISS_TAGGED_CO.  Bug fix for FILENAME in routine
!        READ_PCO_LCO_STRAT. (ave, bnd, bmy, 6/3/03)
!  (23) Updated arg list in call to EMISOP in EMISS_TAGGED_CO (bmy, 12/9/03)
!  (24) Now references "directory_mod.f", "logical_mod.f", "tracer_mod.f".
!        Now remove IJLOOP_CO. (bmy, 7/20/04)
!  (25) Fixed bug in CHEM_TAGGED_CO (bmy, 3/7/05)
!  (26) Now reads data from both GEOS and GCAP grids.  Now also references
!        "tropopause_mod.f". (bmy, 8/16/05)
!  (27) Now modified for new "biomass_mod.f" (bmy, 4/5/06)
!  (28) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  (29) Routines GET_ALPHA_ISOP, GET_PCO_LCO_STRAT, READ_PCO_LCO_STRAT,
!        READ_ACETONE and INIT_TAGGED_CO are public now. Variable
!        EMACET is public now.  (phs, 9/18/07)
!  (30) Revert emiss to v8-01-01 If using adjoint, also set LSPLIT to FALSE
!       (adj_group, 6/8/09)
!******************************************************************************
!
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tagged_co_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE BB_REGION, FF_REGION,   SUMCH3OHCO
      PRIVATE SUMISOPCO, SUMMONOCO,   SUMACETCO
      PRIVATE CO_PRODS,  CO_LOSSS,    FMOL_CO,   XNUMOL_CO
      PRIVATE FMOL_ISOP, XNUMOL_ISOP, FMOL_MONO, XNUMOL_MONO      
      
      ! PRIVATE module routines
      PRIVATE DEFINE_FF_CO_REGIONS, DEFINE_BB_CO_REGIONS 
      PRIVATE GET_SCALE_GROUP

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER, ALLOCATABLE :: BB_REGION(:,:)
      INTEGER, ALLOCATABLE :: FF_REGION(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMCH3OHCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMISOPCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMMONOCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMACETCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMACET(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CO_PRODS(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CO_LOSSS(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: COSF(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CO_CH4(:,:,:,:)

      ! FMOL_CO     - kg CO / mole CO
      ! XNUMOL_CO   - molecules CO / kg CO
      TYPE (XPLEX),PARAMETER::FMOL_CO  = xplex(28d-3,0d0)
      TYPE (XPLEX),PARAMETER::XNUMOL_CO= xplex(6.022d+23/FMOL_CO%r,0d0)

      ! FMOL_ISOP   - kg ISOP / mole ISOP
      ! XNUMOL_ISOP - molecules CO / kg ISOP
      TYPE (XPLEX),PARAMETER::FMOL_ISOP   = xplex(12d-3,0d0)
      TYPE(XPLEX),PARAMETER::XNUMOL_ISOP=xplex(6.022d+23/
     & FMOL_ISOP%r,0d0)

      ! FMOL_MONO   - kg MONO / mole MONO
      ! XNUMOL_MONO - molecules MONO / kg MONO
      TYPE (XPLEX),PARAMETER::FMOL_MONO  =xplex(12d-3,0d0)
      TYPE (XPLEX),PARAMETER::XNUMOL_MONO=xplex(6.022d+23/
     & FMOL_MONO%r,0d0)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_FF_CO_REGIONS( REGION )
!
!******************************************************************************
!  Subroutine DEFINE_FF_CO_REGIONS defines the geographic regions
!  for fossil fuel emissions for the Tagged CO simulation. 
!  (bnd, bmy, 7/21/00, 2/3/03)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) REGION (INTEGER) : Array of Fossil Fuel CO regions 
!
!  NOTES:
!  (1 ) This subroutine only has to be at the beginning of the run, 
!        since the regions don't change with time.
!  (2 ) Regions are hardwired for now -- change if necessary
!  (3 ) Now use regions from bnd simulation (bmy, 6/14/01)
!  (4 ) REGION is now of size (IIPAR,JJPAR) (bmy, 9/28/01)
!  (5 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (6 ) Now uses routines GET_XMID and GET_YMID from "grid_mod.f" to compute
!        grid box lon & lat.  Remove XMID, YMID from arg list. (bmy, 2/3/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(OUT) :: REGION(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: I, J
      TYPE (XPLEX)               :: X, Y

      !=================================================================
      ! DEFINE_FF_CO_REGIONS begins here!
      !
      ! Loop over entire globe -- multiprocessor
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, X, Y )
      DO J = 1, JJPAR
         Y = GET_YMID( J )         ! Latitude [degrees]

         DO I = 1, IIPAR
            X = GET_XMID( I )      ! Longitude [degrees]
         
            ! Region #2 -- North American FF CO 
            IF ( ( X >= -172.5 .and. X < -17.5 )   .and. 
     &           ( Y >=   24.0 .and. Y <  88.0 ) ) THEN
               REGION(I,J) = 2
            
            ! Region #3 -- European FF CO (1st sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 72.5 )  .and. 
     &                ( Y >=  36.0 .and. Y < 45.0 ) ) THEN
               REGION(I,J) = 3

            ! Region #3 -- European FF CO (2nd sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 172.5 )  .and. 
     &                ( Y >=  45.0 .and. Y <  88.0 ) ) THEN
               REGION(I,J) = 3

            ! Region #4 -- Asian FF CO (1st sub-box) 
            ELSE IF ( ( X >= 70.0 .and. X < 152.5 )  .and.
     &                ( Y >=  8.0 .and. Y <  45.0 ) ) THEN
               REGION(I,J) = 4
   
            ! Region #5 -- FF CO from everywhere else
            ELSE
               REGION(I,J) = 5
               
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEFINE_FF_CO_REGIONS

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_BB_CO_REGIONS( REGION )
!
!******************************************************************************
!  Subroutine DEFINE_BB_CO_REGIONS defines the geographic regions 
!  for biomass burning emissions for the Tagged CO simulation. 
!  (bnd, bmy, 7/21/00, 2/4/03)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) REGION (INTEGER) : Array of Fossil Fuel CO regions 
!
!  NOTES:
!  (1 ) This subroutine only has to be at the beginning of the run, 
!        since the regions don't change with time.
!  (2 ) Regions are hardwired for now -- change if necessary
!  (3 ) Now use regions from bnd simulation (bmy, 6/14/01)
!  (4 ) REGION is now of size (IIPAR,JJPAR) (bmy, 9/28/01)
!  (5 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (6 ) Now uses routines GET_XMID and GET_YMID from "grid_mod.f" to compute
!        grid box lon & lat.  Remove XMID, YMID from arg list. (bmy, 2/3/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(OUT) :: REGION(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: I, J
      TYPE (XPLEX)               :: X, Y

      !=================================================================
      ! DEFINE_BB_CO_REGIONS begins here!
      ! Loop over entire globe
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, X, Y )
      DO J = 1, JJPAR
         Y = GET_YMID( J )      ! Latitude [degrees]

         DO I = 1, IIPAR
            X = GET_XMID( I )   ! Longitude [degrees]

            ! Region #6 -- South American BB CO 
            IF ( ( X >= -112.5 .and. X < -32.5 )  .and.
     &           ( Y >=  -56.0 .and. Y <  24.0 ) ) THEN
               REGION(I,J) = 6

            ! Region #7 -- African BB CO 
            ELSE IF ( ( X >= -17.5 .and. X < 70.0 )  .and.
     &                ( Y >= -48.0 .and. Y < 36.0 ) ) THEN
               REGION(I,J) = 7

            ! Region #8 -- SE Asian BB CO (1st sub-box) 
            ELSE IF ( ( X >= 70.0 .and. X < 152.5 )  .and.
     &                ( Y >=  8.0 .and. Y <  45.0 ) ) THEN
               REGION(I,J) = 8

            ! Region #9 -- Oceania BB CO 
            ELSE IF ( ( X >= 70.0 .and. X < 170.0 )  .and.
     &                ( Y <  8.0 ) ) THEN
               REGION(I,J) = 9

            ! Region #10 -- European BB CO (1st sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 72.5 )  .and.
     &                ( Y >=  36.0 .and. Y < 45.0 ) ) THEN
               REGION(I,J) = 10

            ! Region #10 -- European BB CO (2nd sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 172.5 )  .and.
     &                ( Y >=  45.0 .and. Y <  88.0 ) ) THEN
               REGION(I,J) = 10 

            ! Region #11 -- North America BB CO  
            ! (rest of world = North America)
            ELSE 
               REGION(I,J) = 11

            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEFINE_BB_CO_REGIONS

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_TAGGED_CO
!
!******************************************************************************
!  Subroutine EMISS_TAGGED_CO reads in CO emissions for the Tagged CO run
!  (bey, bdf, bnd, bmy, 7/21/00, 6/30/08)
!
!  NOTES:
!  (1 ) Adapted from "emiss_co.f" (bmy, 7/2100)
!  (2 ) Add multiprocessor commands for DO-loops (bmy, 7/21/00)
!  (3 ) Added references to F90 modules "biomass_mod.f" and "biofuel_mod.f".
!        Also, TWOODIJ is now called BIOFUEL.  Finally, BURNEMIS is now 
!        referenced with IREF = I + I0 and JREF = J + J0. (bmy, 9/12/00)
!  (4 ) Now reference FF_REGION with IREF & JREF.  Removed obsolete code
!        from 7/00 and 9/00.  Also define NFAM in "ndxx_setup.f" (bmy, 10/6/00)
!  (5 ) Added variables EMMO and references to functions EMMONOT, EMCH3OH -- 
!        these are for CO production by VOC oxidation. (bnd, bmy, 1/2/01)
!  (6 ) Merged ISOP and MONOTERPENE sources into the same DO-loop, for
!        computational efficiency (bmy, 1/3/01)
!  (7 ) Make sure that SCALEYEAR does not go higher than 1997.  1997 is
!        the last year we have FF scale factor data. (bnd, bmy, 4/6/01)
!  (8 ) Now call READ_ACETONE to read acetone emissions into EMACET.  
!        Also call READ_BIOG_FOR_GEOS3 to read acetone emissions into 
!        ISOP96, MONO96, CH3OH96 for GEOS-3 simulations. (bnd, bmy, 6/14/01)
!  (9 ) Remove GEOS-3 fix for biogenic fields (ISOP96, MONO96, CH3OH96),
!        since we now use updated met fields w/o the surface temperature
!        problem (bmy, 8/21/01)
!  (10) Now archive ND46 as atoms C/cm2/s here.  Also deleted obsolete code
!        for the 1998 GEOS-3 fix. (bmy, 9/13/01)
!  (11) BXHEIGHT is now dimensioned IIPAR,JJPAR,LLPAR.  BIOFUEL(:,IREF,JREF)
!        is now BIOFUEL(:,I,J).  BB_REGION(IREF,JREF) is now BB_REGION(I,J).
!        FF_REGION(IREF,JREF) is now FF_REGION(I,J). (bmy, 9/28/01)
!  (12) Removed obsolete code from 9/13/01 and 9/28/01 (bmy, 10/22/01)
!  (13) Now reference SUNCOS "dao_mod.f".  Remove SUNCOS and BXHEIGHT from 
!        the arg list -- BXHEIGHT isn't used!.  Now make FIRSTEMISS a local
!        SAVEd variable.  Now do not let SCALEYEAR exceed 1998.  Now
!        reference IDBCO, IDBFCO from "tracerid_mod.f". (bmy, 1/13/03)
!  (14) Remove XMID, YLMID from call to routines DEFINE_BB_CO_REGIONS and 
!        DEFINE_FF_CO_REGIONS.   Now replace DXYP(JREF)*1d4 with routine
!        GET_AREA_CM2 of "grid_mod.f".  Removed IMONTH from call to BIOBURN.
!        Now uses functions GET_MONTH, GET_YEAR, GET_ELAPSED_MIN, and 
!        GET_TS_EMIS from the new "time_mod.f".  Now only passes I in the
!        call to GET_IHOUR. Now use functions GET_XOFFSET, GET_YOFFSET from
!        "grid_mod.f".  I0 and J0 are now local variables. (bmy, 2/11/03)
!  (15) Bug fix: NTAU should be the integer value of TAU (bmy, 6/10/03)
!  (16) Now pass I, J to EMISOP (bmy, 12/9/03)
!  (17) Now reference LSPLIT, LANTHRO, LBIOMASS, & LBIOFUEL from 
!        "logical_mod.f".  Now reference STT from "tracer_mod.f".  Now replace
!        IJLOOP_CO w/ an analytic function. (bmy, 7/20/04)
!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (19) Now modified for the new "biomass_mod.f" (bmy, 4/5/06)
!  (20) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  (21) Now take CO emissions from the EMISRR array, which is archived by
!        calling EMISSDR.  Remove duplicate scaling and other operations which
!        are now done in EMISSDR. Remove references to BIOFUEL_BURN and
!        all routines from GEIA_MOD. (jaf, mak, bmy, 2/14/08)
!  (22) Bug fix: Now use IDECO to be consistent (phs, bmy, 6/30/08)
!  (23) Now distribute CO biomass burning source through the PBL
!        (yc, phs, 12/23/08)
!  (24) For adjoint run, revert to v8-01-01 and set LSPLIT To .false.
!       (adj_group, 6/08/09)
!  (25) For adjoint run, now distribute biomass emiss throughout the BL
!       (dkh, 02/04/11) 
!******************************************************************************
!
      ! References to F90 modules
      USE BIOFUEL_MOD,  ONLY : BIOFUEL
      USE BIOMASS_MOD,  ONLY : BIOMASS,       IDBCO
      USE DAO_MOD,      ONLY : SUNCOS
      USE DIAG_MOD,     ONLY : AD29,          AD46
      USE GRID_MOD,     ONLY : GET_XOFFSET,   GET_YOFFSET, GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LSPLIT,        LANTHRO
      USE LOGICAL_MOD,  ONLY : LBIOMASS,      LBIOFUEL
      USE PBL_MIX_MOD,  ONLY : GET_PBL_MAX_L, GET_FRAC_OF_PBL
      USE TIME_MOD,     ONLY : GET_MONTH, GET_TAU
      USE TIME_MOD,     ONLY : GET_YEAR, GET_TS_EMIS  
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDBFCO,        IDECO
      USE LOGICAL_ADJ_MOD, ONLY : LADJ, LADJ_EMS
      USE ADJ_ARRAYS_MOD,  ONLY : ADCOEMS, EMS_SF
      USE ADJ_ARRAYS_MOD,       ONLY : IFD, JFD
      USE LOGICAL_ADJ_MOD,      ONLY : LPRINTFD
      USE BPCH2_MOD
      
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! FSCALYR, SCNR89, TODH, EMISTCO, EMISRR
#     include "CMN_DIAG"     ! Diagnostic arrays & switches

      ! Local variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE.
      INTEGER                :: I, J, L, N, I0, J0, M
      INTEGER                :: AS, IREF, JREF, IJLOOP, PBL_MAX
      INTEGER, SAVE          :: LASTMONTH = -999

      ! For now these are defined in CMN_O3
      !TYPE (XPLEX)                 :: EMISTCO(IGLOB,JGLOB)
      !TYPE (XPLEX)                 :: FLIQCO2(IGLOB,JGLOB)

      TYPE (XPLEX)                 :: TMMP, EMXX, EMX, EMMO, F_OF_PBL 
      TYPE (XPLEX)                 :: EMAC, E_CO, DTSRCE, AREA_CM2
      TYPE (XPLEX)                 :: XTAU, ED_CO
      
      ! External functions
      TYPE (XPLEX), EXTERNAL       :: XLTMMP,  EMISOP, BOXVL
      TYPE (XPLEX), EXTERNAL       :: EMMONOT, EMCH3OH

      !=================================================================
      ! EMISS_TAGGED_CO begins here!
      !
      ! Do the following only on the first call to EMISS_TAGGED_CO...
      !=================================================================
      IF ( FIRSTEMISS ) THEN 

         ! Allocate all module arrays
         CALL INIT_TAGGED_CO

         ! Define geographic regions for both fossil fuel & biomass burning
         CALL DEFINE_FF_CO_REGIONS( FF_REGION )         
         CALL DEFINE_BB_CO_REGIONS( BB_REGION )

         ! Set first-time flag to false
         FIRSTEMISS = .FALSE.
         
         XTAU = GET_TAU()
         ! Obsolete (zhej, dkh, 01/17/12, adj32_017) 
         !CALL GET_OPT_CO( XTAU )  !(zhe 11/28/10)
         
      ENDIF

      !=================================================================
      ! Once a month, read acetone from disk.  For GEOS-3, also read 
      ! P(CO) from ISOPRENE, MONOTERPENES, and METHANOL from 1996
      !=================================================================
      IF ( GET_MONTH() /= LASTMONTH ) THEN
      
         ! Read acetone for this month
         CALL READ_ACETONE( GET_MONTH() )
      
         ! Save month for next iteration
         LASTMONTH = GET_MONTH()
    
      ENDIF
      
      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()
      ! ===============================================================
      ! If we're doing adjoint run
      IF ( LADJ .AND. LADJ_EMS ) THEN 
         LSPLIT = .FALSE.
         print*, 'For an adjoint run, revert to emissions from v8-01-01'

          ! Determine group (temporal)
          M = GET_SCALE_GROUP()
          ! Print out scaling info
          WRITE(6,*) '- READ / RESCALE EMISSIONS: use SCALE_GROUP ',  M 


         !=================================================================
         ! Process Anthropogenic (Fossil Fuel) CO emissions
         !
         ! Anthropogenic emissions are enhanced by 18.5% below.  This
         ! accounts for production of CO from oxidation of certain VOC's, 
         ! which are not explicitly carried by GEOS-CHEM as anthropogenic
         ! species.  This needs to be done here since a different scale
         ! factor is used for the full chemistry run.  Also update the
         ! ND29 diagnostic below, in order to archive the correct 
         ! emissions. (bmy, 6/14/01)
         !
         ! NOTES: 
         ! (1) Anthro CO emissions come from the GEIA/Piccot inventory.
         !     (bmy, 1/2/01)
         ! (2) Need to save ND29 diagnostics here (bmy, 1/2/01)
         !=================================================================
         ED_CO = 0

         ! COSF now obsolete (zhej, dkh, 01/17/12, adj32_017) 
         PBL_MAX = GET_PBL_MAX_L()

         IF ( LANTHRO ) THEN

! F_OF_PBL and ED_CO thread private (zhej, dkh, 01/17/12, adj32_017) 
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N, IREF, JREF, E_CO, F_OF_PBL, ED_CO )

         DO J = 1, JJPAR
            JREF = J + J0

            DO I = 1, IIPAR
               IREF = I + I0

               ! E_CO = Fossil Fuel CO emissions in [molec CO/s]
               ! EMISRR is archived in "emissdr.f" (jaf, mak, bmy, 2/14/08)
               !-------------------------------------------------------------
               ! Prior to 6/30/08:
               ! Now use IDECO to be consistent (bmy, 6/30/08)
               !E_CO = EMISRR(IREF,JREF,1)
               !-------------------------------------------------------------
               E_CO = EMISRR(IREF,JREF,IDECO)

               ! Convert from [molec CO/s] to [kg CO] 
               ! (jaf, mak, bmy, 2/14/08)
               E_CO = E_CO * ( DTSRCE / XNUMOL_CO )
               
               IF ( ND29 > 0 ) THEN  !(zhe 11/28/10)
                   AD29(I,J,6) = AD29(I,J,6) + E_CO
               ENDIF
               
               !----------------------------------------------------------------------
               ! COSF now obsolete (zhej, dkh, 01/17/12, adj32_017)
               ! OLD:
               !E_CO = E_CO * COSF(I,J,1)  !(zhe 11/28/10)
               !
               !! Add FF CO to Tracer #1 -- total CO [kg CO]
               !STT(I,J,1,1) = STT(I,J,1,1) + 
               !            EMS_SF(I,J,M,ADCOEMS) * E_CO
               !! Split FF CO into geographic regions -- Tracers #2 - #5
               !IF ( LSPLIT ) THEN
               !   N            = FF_REGION(I,J)
               !   STT(I,J,1,N) = STT(I,J,1,N) + E_CO
               !ENDIF
               ! NEW:
               DO L = 1, PBL_MAX

                  F_OF_PBL = GET_FRAC_OF_PBL(I,J,L)

                  STT(I,J,L,1) = STT(I,J,L,1) + E_CO * F_OF_PBL
     &                         * EMS_SF(I,J,M,ADCOEMS)

                  !Split FF CO into geographic regions -- Tracers #6 - #11
                  IF ( LSPLIT ) THEN
                     N            = BB_REGION(I,J)
                     STT(I,J,L,N) = STT(I,J,L,N) + E_CO * F_OF_PBL
                  ENDIF

               ENDDO
               !----------------------------------------------------------------------
               

               IF(I == IFD .AND. J == JFD .and. LPRINTFD) THEN
                 PRINT*, 'CO_FF=', E_CO
                 ED_CO = ED_CO + E_CO
               ENDIF

            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
      !=================================================================
      ! Process biomass burning CO emissions, stored in array
      ! BIOMASS(:,:,IDBCO) which has units of [molec/cm2/s]
      !
      ! The default Duncan et al 2001 biomass burning emissions are 
      ! enhanced by 16.4% within the routines of "gc_biomass_mod.f".  
      ! This accounts for production of CO from oxidation of certain 
      ! VOC's, which are not explicitly carried by GEOS-Chem as biomass 
      ! burning species.  The scaling needs to be done in "biomass_mod.f"
      ! so that the diagnostics will archive the correct emissions. 
      ! (bmy, 6/8/01, 9/27/06)
      !
      ! GFED2 CO biomass burning emissions are not scaled any further.
      !
      ! NOTES:
      ! (1) Some forest fires generate strong convection columns.  
      !      However, we release biomass burning emissions only into 
      !      the surface layer. (bnd, bmy, 1/3/01)
      ! (2) ND29 diagnostics are saved within routine BIOBURN.
      !      (bmy, 1/3/01)
      ! (3) Earlier versions of the adjoint put all the biomass burning
      !      emissions in level 1.  Now we distribute them throughout
      !      the boundary layer just like standard model. (dkh, 02/04/11) 
      !=================================================================
      ! (zhej, dkh, 01/17/12, adj32_017) 
      !PBL_MAX = GET_PBL_MAX_L()

      IF ( LBIOMASS ) THEN

! Add ED_CO as thread PRIVATE (zhej, dkh, 01/17/12, adj32_017) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( E_CO, I, J, L, F_OF_PBL, N, AREA_CM2, ED_CO )



         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert [molec CO/cm2/s] to [kg CO] and store in E_CO
            E_CO = ( BIOMASS(I,J,IDBCO) / XNUMOL_CO ) *
     &             ( AREA_CM2           * DTSRCE    )
            
            IF ( ND29 > 0 ) THEN  !(zhe 11/28/10)
                AD29(I,J,6) = AD29(I,J,6) + E_CO
            ENDIF
            
            ! COSF obsolete (zhej, dkh, 01/17/12, adj32_017) 
            !E_CO = E_CO * COSF(I,J,1)  !(zhe 11/28/10)

!---------------------------------------------------------------------------
! OLD: just emis into one level
!            ! Add biomass burning to tracer #1 -- total CO [kg CO]
!            STT(I,J,1,1) = STT(I,J,1,1) + 
!     &                EMS_SF(I,J,M,ADCOEMS) * E_CO
!               ! Split BB CO into geographic regions -- Tracers #6 - #11
! NEW: throughout the BL (dkh, 02/04/11) 
            DO L = 1, PBL_MAX

               F_OF_PBL = GET_FRAC_OF_PBL (I,J,L)

               ! Add biomass burning to tracer #1 -- total CO [kg CO]
!-- prior to 12/23/08
!               STT(I,J,1,1) = STT(I,J,1,1) + E_CO
               STT(I,J,L,1) = STT(I,J,L,1) + E_CO * F_OF_PBL 
     &                      * EMS_SF(I,J,M,ADCOEMS)

! LSPLIT = FALSE for adjoint 
!               ! Split BB CO into geographic regions -- Tracers #6 - #11
!               IF ( LSPLIT ) THEN
!                  N            = BB_REGION(I,J)
!!-- prior to 12/23/08
!!                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO
!                  STT(I,J,L,N) = STT(I,J,L,N) + E_CO * F_OF_PBL                  
!               ENDIF

               ! Reinstate (zhej, dkh, 01/17/12, adj32_017) 
               !Split BB CO into geographic regions -- Tracers #6 - #11
               IF ( LSPLIT ) THEN
                  N            = BB_REGION(I,J)
                  STT(I,J,L,N) = STT(I,J,L,N) + E_CO * F_OF_PBL
               ENDIF

            ENDDO
 
!---------------------------------------------------------------------------
            
            IF(I == IFD .AND. J == JFD .and. LPRINTFD) THEN
               PRINT*, 'CO_BB=', E_CO
               ED_CO = ED_CO + E_CO
            ENDIF

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Process biofuel (formerly wood burning) CO emissions
      ! stored in BIOFUEL(IDBCO,IREF,JREF) in [molec/cm3/s]
      !
      ! Biofuel burning emissions are enhanced by 18.9% within the
      ! routines of "biofuel_mod.f".  This accounts for production 
      ! of CO from oxidation of certain VOC's, which are not explicitly
      ! carried by GEOS-CHEM as biofuel burning species.  The scaling
      ! needs to be done in "biofuel_mod.f" so that the diagnostics 
      ! will archive the correct emissions. (bmy, 6/8/01)
      !
      ! NOTES:
      ! (1 ) ND29 diagnostics are saved within routine BIOFUEL_BURN.
      !       (bmy, 1/2/01)
      ! (2 ) Now use IDBFCO to index the proper element of the 
      !       biofuel burning array (bmy, 6/8/01).
      ! (3 ) Now add biofuel burning to tagged tracer #13 (bmy, 6/14/01)
      !=================================================================
      IF ( LBIOFUEL ) THEN
! make F_OF_PBL and ED_CO thread private (zhej, dkh, 01/17/12, adj32_017) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( E_CO, I, J, N, F_OF_PBL, ED_CO )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Convert biofuel CO from [molec CO/cm3/s] to [kg CO]
            E_CO = BIOFUEL(IDBFCO,I,J) / XNUMOL_CO * 
     &             BOXVL(I,J,1)        * DTSRCE

            IF ( ND29 > 0 ) THEN  !(zhe 11/28/10)
                AD29(I,J,6) = AD29(I,J,6) + E_CO
            ENDIF
            
            !----------------------------------------------------------------------
            ! COSF obsolete (zhej, dkh, 01/17/12, adj32_017) 
            ! OLD:
            !E_CO = E_CO * COSF(I,J,1) !(zhe 11/28/10)
            ! 
            !! Add biofuel CO burning to tracer #1 -- total CO [kg CO]
            !STT(I,J,1,1) = STT(I,J,1,1) + 
     &      !           EMS_SF(I,J,M,ADCOEMS) * E_CO
            !! Add biofuel CO burning [kg CO] as Tracer #13 (bmy, 6/14/01)
            !!
            !! %%% NOTE: If you are modifying the tagged CO simulation, 
            !! %%% and your simulation has less than 13 tracers, then 
            !! %%% then comment out this IF block.  If you don't you can
            !! %%% get an array-out-of-bounds error (bmy, 6/11/08)
            !IF ( LSPLIT ) THEN
            !   STT(I,J,1,13) = STT(I,J,1,13) + E_CO
            !ENDIF
            ! NEW:
            DO L = 1, PBL_MAX

               F_OF_PBL = GET_FRAC_OF_PBL(I,J,L)

               STT(I,J,L,1) = STT(I,J,L,1) + E_CO * F_OF_PBL
     &                      * EMS_SF(I,J,M,ADCOEMS)

               !Split BF CO into geographic regions -- Tracers #6 - #11
               IF ( LSPLIT ) THEN
                  N            = BB_REGION(I,J)
                  STT(I,J,L,N) = STT(I,J,L,N) + E_CO * F_OF_PBL
               ENDIF

            ENDDO
            !----------------------------------------------------------------------

            IF(I == IFD .AND. J == JFD .and. LPRINTFD) THEN
               PRINT*, 'CO_BF=', E_CO
               ED_CO = ED_CO + E_CO
            ENDIF

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
      
      IF(LPRINTFD) THEN
        PRINT*, 'STT=', STT(IFD,JFD,1,1)
        PRINT*, 'TOTAL DIRECT CO=', ED_CO
      ENDIF

      ELSE ! LADJ and/or LADJ_EMS is false
      !=================================================================
      ! Process Anthropogenic (Fossil Fuel) CO emissions
      !
      ! Anthropogenic emissions are enhanced by 18.5% below.  This
      ! accounts for production of CO from oxidation of certain VOC's, 
      ! which are not explicitly carried by GEOS-CHEM as anthropogenic
      ! species.  This needs to be done here since a different scale
      ! factor is used for the full chemistry run.  Also update the
      ! ND29 diagnostic below, in order to archive the correct 
      ! emissions. (bmy, 6/14/01)
      !
      ! NOTES: 
      ! (1) Anthro CO emissions come from the GEIA/Piccot inventory.
      !     (bmy, 1/2/01)
      ! (2) Need to save ND29 diagnostics here (bmy, 1/2/01)
      !=================================================================
      IF ( LANTHRO ) THEN

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( E_CO, I, IREF, J, JREF, N )
         DO J = 1, JJPAR
            JREF = J + J0

            DO I = 1, IIPAR
               IREF = I + I0

               ! E_CO = Fossil Fuel CO emissions in [molec CO/s]
               ! EMISRR is archived in "emissdr.f" (jaf, mak, bmy, 2/14/08)
               E_CO = EMISRR(IREF,JREF,IDECO)

               ! Convert from [molec CO/s] to [kg CO] 
               ! (jaf, mak, bmy, 2/14/08)
               E_CO = E_CO * ( DTSRCE / XNUMOL_CO )  

               ! Add FF CO to Tracer #1 -- total CO [kg CO]
               STT(I,J,1,1) = STT(I,J,1,1) + E_CO

               ! Split FF CO into geographic regions -- Tracers #2 - #5
               IF ( LSPLIT ) THEN
                  N            = FF_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO
               ENDIF
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Process biomass burning CO emissions, stored in array
      ! BIOMASS(:,:,IDBCO) which has units of [molec/cm2/s]
      !
      ! The default Duncan et al 2001 biomass burning emissions are 
      ! enhanced by 16.4% within the routines of "gc_biomass_mod.f".  
      ! This accounts for production of CO from oxidation of certain 
      ! VOC's, which are not explicitly carried by GEOS-Chem as biomass 
      ! burning species.  The scaling needs to be done in "biomass_mod.f"
      ! so that the diagnostics will archive the correct emissions. 
      ! (bmy, 6/8/01, 9/27/06)
      !
      ! GFED2 CO biomass burning emissions are not scaled any further.
      !
      ! NOTES:
      ! (1) Some forest fires generate strong convection columns.  
      !      However, we release biomass burning emissions only into 
      !      the surface layer. (bnd, bmy, 1/3/01)
      ! (2) ND29 diagnostics are saved within routine BIOBURN.
      !      (bmy, 1/3/01)
      !=================================================================
      PBL_MAX = GET_PBL_MAX_L()

      IF ( LBIOMASS ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( E_CO, I, J, L, F_OF_PBL, N, AREA_CM2 )

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert [molec CO/cm2/s] to [kg CO] and store in E_CO
            E_CO = ( BIOMASS(I,J,IDBCO) / XNUMOL_CO ) *
     &             ( AREA_CM2           * DTSRCE    ) 

            DO L = 1, PBL_MAX
           
               F_OF_PBL = GET_FRAC_OF_PBL (I,J,L) 

               ! Add biomass burning to tracer #1 -- total CO [kg CO]
!-- prior to 12/23/08
!               STT(I,J,1,1) = STT(I,J,1,1) + E_CO
               STT(I,J,L,1) = STT(I,J,L,1) + E_CO * F_OF_PBL
            
               ! Split BB CO into geographic regions -- Tracers #6 - #11
               IF ( LSPLIT ) THEN
                  N            = BB_REGION(I,J)
!-- prior to 12/23/08
!                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO
                  STT(I,J,L,N) = STT(I,J,L,N) + E_CO * F_OF_PBL                  
               ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Process biofuel (formerly wood burning) CO emissions
      ! stored in BIOFUEL(IDBCO,IREF,JREF) in [molec/cm3/s]
      !
      ! Biofuel burning emissions are enhanced by 18.9% within the
      ! routines of "biofuel_mod.f".  This accounts for production 
      ! of CO from oxidation of certain VOC's, which are not explicitly
      ! carried by GEOS-CHEM as biofuel burning species.  The scaling
      ! needs to be done in "biofuel_mod.f" so that the diagnostics 
      ! will archive the correct emissions. (bmy, 6/8/01)
      !
      ! NOTES:
      ! (1 ) ND29 diagnostics are saved within routine BIOFUEL_BURN.
      !       (bmy, 1/2/01)
      ! (2 ) Now use IDBFCO to index the proper element of the 
      !       biofuel burning array (bmy, 6/8/01).
      ! (3 ) Now add biofuel burning to tagged tracer #13 (bmy, 6/14/01)
      !=================================================================
      IF ( LBIOFUEL ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( E_CO, I, J, N )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Convert biofuel CO from [molec CO/cm3/s] to [kg CO]
            E_CO = BIOFUEL(IDBFCO,I,J) / XNUMOL_CO * 
     &             BOXVL(I,J,1)        * DTSRCE

            ! Add biofuel CO burning to tracer #1 -- total CO [kg CO]
            STT(I,J,1,1) = STT(I,J,1,1) + E_CO

            ! Add biofuel CO burning [kg CO] as Tracer #13 (bmy, 6/14/01)
            !
            ! %%% NOTE: If you are modifying the tagged CO simulation, 
            ! %%% and your simulation has less than 13 tracers, then 
            ! %%% then comment out this IF block.  If you don't you can
            ! %%% get an array-out-of-bounds error (bmy, 6/11/08)
            IF ( LSPLIT ) THEN
               STT(I,J,1,13) = STT(I,J,1,13) + E_CO
            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ENDIF ! LADJ true or false
      !=================================================================
      ! Process emissions of ISOPRENE, MONOTERPENES, METHANOL
      ! and ACETONE -- save into summing arrays for later use
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, IJLOOP, TMMP, EMXX, EMMO, EMAC )
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

         DO I = 1, IIPAR

            ! 1-D array index 
            IJLOOP = ( (J-1) * IIPAR ) + I

            !===========================================================
            ! The CO yields from ISOP, MONOTERPENES, and CH3OH will be
            ! computed in subroutine CHEM_TAGGED_CO.  P(CO) from CH3OH
            ! will be scaled to isoprene emissions within subroutine
            ! CHEM_TAGGED_CO (bnd, bmy, 6/14/01)
            !===========================================================

            ! Surface air temperature [K]
            TMMP = XLTMMP(I,J,IJLOOP)
      
            ! ISOP and MONOTERPENE emissions in [atoms C/box/time step] 
            ! SUNCOS is COSINE( solar zenith angle ) 
            EMXX = EMISOP( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL_ISOP ) 
            EMMO = EMMONOT( IJLOOP, TMMP, XNUMOL_MONO )

            ! Store ISOP and MONOTERPENE emissions [atoms C/box/time step]
            ! for later use in the subroutine CHEM_TAGGED_CO
            SUMISOPCO(I,J) = SUMISOPCO(I,J) + EMXX
            SUMMONOCO(I,J) = SUMMONOCO(I,J) + EMMO

            ! ND46 -- save biogenic emissions [atoms C/cm2/s] here 
            IF ( ND46 > 0 ) THEN

               ! Isoprene 
               AD46(I,J,1) = AD46(I,J,1) + ( EMXX / AREA_CM2 / DTSRCE )

               ! Monoterpenes 
               AD46(I,J,4) = AD46(I,J,4) + ( EMMO / AREA_CM2 / DTSRCE )

            ENDIF

            !===========================================================
            ! For GEOS-1, GEOS-STRAT, GEOS-3, extract acetone emission
            ! fluxes the EMACET array for the current month
            !===========================================================

            ! EMAC = [atoms C/box/s] from acetone
            EMAC = EMACET( I, J )

            ! Sum acetone loss for use in chemco_decay
            ! Units = [atoms C/box/timestep]
            SUMACETCO(I,J) = SUMACETCO(I,J) + (EMAC * DTSRCE * AREA_CM2) 
            
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMISS_TAGGED_CO

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_TAGGED_CO
!
!******************************************************************************
!  Subroutine CHEM_TAGGED_CO performs CO chemistry on geographically
!  "tagged" CO tracers.  Loss is via reaction with OH. 
!  (qli, bnd, bdf, bmy, 10/19/99, 5/31/06)
!
!  NOTES:
!  (1 ) Now do chemistry all the way to the model top. 
!  (2 ) Use monthly mean OH fields for oxidation.
!  (3 ) Now reference the monthly mean OH array and the routine which reads
!        it from disk in "global_oh_mod.f" (bmy, 7/28/00)
!  (4 ) Removed obsolete code from 10/6/00 (bmy, 12/21/00)
!  (5 ) Added P(CO) from CH3OH and MONOTERPENES.  Also account for the 
!        variation of CH4 conc. w/ latitude and year (bnd, bmy, 1/2/01)
!  (6 ) Removed obsolete commented-out code (bmy, 4/23/01)
!  (7 ) Updated with new stuff from bnd.  Also references module routines
!        & arrays from "global_nox_mod.f" and "error_mod.f".  Removed
!        THISMONTH as an argument since we can use the MONTH value
!        from the "CMN" include file. (bmy, 6/14/01)
!  (8 ) Remove GEOS-3 fix for biogenic fields (ISOP96, MONO96, CH3OH96),
!        since we now use updated met fields w/o the surface temperature
!        problem (bmy, 8/21/01)
!  (9 ) Now only call GLOBAL_NOX_MOD to define the BNOX array, if switch
!        ALPHA_ISOP_FROM_NOX is set.  The NOx concentrations used to compute
!        the CO yield from isoprene are currently only at 4x5 (bmy, 8/28/01)
!  (10) Bug fix: now make sure to add CO production and loss into the 
!        STT(:,:,:,1), which is the Total CO tracer. (qli, bmy, 9/21/01)
!  (11) Updated comments.  Bug fix: multiply CO_OH (after using it to update 
!        tagged tracers) by GCO (the initial value of STT in molec/cm3) to 
!        convert it to an amount of CO lost by OH [molec/cm3]. (bmy, 2/19/02)
!  (12) Removed PS as an argument; use P(I,J) + PTOP instead of PS, in order
!        to ensure that we use P and AD computed from the same preesure
!        by AIRQNT. (bmy, 4/11/02)
!  (13) Now use GET_PCENTER from "pressure_mod.f" to compute the pressure
!        at the midpoint of box (I,J,L).  Also deleted obsolete, commented-out
!        code. (dsa, bdf, bmy, 8/21/02)
!  (14) Now reference AD and T from "dao_mod.f".  Now make FIRSTCHEM a local
!        SAVEd variable. (bmy, 11/15/02)
!  (15) Now replace YLMID(J) with routine GET_YMID of "grid_mod.f".  Now
!        uses functions GET_TS_CHEM, GET_MONTH, GET_YEAR from the new
!        "time_mod.f".  (bmy, 2/10/03) 
!  (16) Now reference STT & N_TRACERS from "tracer_mod.f".  Now references
!        LSPLIT from "logical_mod.f".  Now references AD65 from 
!        "diag_pl_mod.f".  Updated comments. (bmy, 7/20/04)
!  (17) Bug fix: re-insert ELSE between (1a-1) and (1a-2); it appears to have
!        been mistakenly deleted. (bmy, 3/7/05)
!  (18) Now references ITS_IN_THE_STRAT from "tropopause_mod.f".  Now remove
!        reference to "CMN", it's obsolete. (bmy, 8/22/05)
!  (19) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (20) Remove reference to "global_ch4_mod.f" (bmy, 5/31/06)
!  (21) Add BIOGSRC and adj scaling (adj_group, 6/08/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD, AIRVOL, T, DELP
      USE DIAG_PL_MOD,    ONLY : AD65
      USE ERROR_MOD,      ONLY : CHECK_VALUE
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH,   OH
      USE GLOBAL_NOX_MOD, ONLY : GET_GLOBAL_NOX,  BNOX
      USE GRID_MOD,       ONLY : GET_YMID
      USE LOGICAL_MOD,    ONLY : LSPLIT
      USE PRESSURE_MOD,   ONLY : GET_PCENTER
      USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_MONTH, GET_YEAR
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE TRACER_MOD,     ONLY : N_TRACERS,       STT
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
      USE ADJ_ARRAYS_MOD, ONLY : ADCOVOX, EMS_SF, NNEMS, ADCOEMS
      USE LOGICAL_ADJ_MOD,ONLY : LADJ, LADJ_EMS, LPRINTFD, LDEVOC  !(zhe 11/28/10)
      USE DIAG_MOD,       ONLY : AD29
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD
      USE BPCH2_MOD
      
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND65

      ! Local variables
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, N, MONTH, M
      TYPE (XPLEX)                 :: ALPHA_CH4, ALPHA_ISOP, ALPHA_MONO
      TYPE (XPLEX)                 :: DTCHEM,    GCO,        PCO   
      TYPE (XPLEX)                 :: STTCO,     KRATE,      CH4
      TYPE (XPLEX)                 :: CO_CH4,    CO_ISOP,    CO_MONO
      TYPE (XPLEX)                 :: CO_CH3OH,  CO_OH,      CO_ACET
      TYPE (XPLEX)                 :: CH4RATE,   DENS,       ALPHA_ACET
      TYPE (XPLEX)                 :: CORATE,  YMID,  CO_VOC

      ! For saving CH4 latitudinal gradient
      TYPE (XPLEX),  SAVE          :: A3090S, A0030S, A0030N, A3090N

      ! External functions
      TYPE (XPLEX),  EXTERNAL      :: BOXVL

      ! WTAIR = molecular weight of air (g/mole)
      TYPE (XPLEX),  PARAMETER     :: WTAIR = xplex(28.966d0,0d0)
      
      ! Switch to scale yield of isoprene from NOx concentration or not
      LOGICAL, PARAMETER     :: ALPHA_ISOP_FROM_NOX = .FALSE.
      CHARACTER(LEN=255)     :: FILENAME
      ! Obsolete (zhej, dkh, 01/17/12, adj32_017) 
      !INTEGER, SAVE          :: TIME_STEPS = 1

      !=================================================================
      ! CHEM_TAGGED_CO begins here!
      !
      ! Do the following on the first calla to CHEM_TAGGED_CO...
      !=================================================================
      IF ( FIRSTCHEM ) THEN
        
          FIRSTCHEM = .FALSE.
         
      ENDIF

      IF ( LADJ .AND. LADJ_EMS .AND. NNEMS .GT. 1) THEN
         M = GET_SCALE_GROUP()
         ! compute temporal scaling group
         WRITE(6,*) 
     &      '- READ/RESCALE background source: use SCALE_GROUP ',M 
      ENDIF 

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Read in OH, NOx, P(CO), and L(CO) fields for the current month
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN
         
         ! Get current month
         MONTH = GET_MONTH()
         
         ! Global OH 
         CALL GET_GLOBAL_OH( MONTH )

         ! Global NOx -- need this to determine 
         ! ALPHA_ISOP which is a function of NOx
         IF ( ALPHA_ISOP_FROM_NOX ) CALL GET_GLOBAL_NOX( MONTH )
            
         ! Read in the loss/production of CO in the stratosphere.
         CALL READ_PCO_LCO_STRAT( MONTH )
      ENDIF

      !=================================================================
      ! Get the yearly and latitudinal gradients for CH4
      ! This only needs to be called once per year
      !
      ! NOTE: If you are going to run w/ future emissions you must 
      !       pass the future emissions year to GET_GLOBAL_CH4.  
      !       See the modification that was made in "chemdr.f" 
      !       (bmy, 1/24/08)
      !=================================================================
      IF ( ITS_A_NEW_YEAR() ) THEN 
         CALL GET_GLOBAL_CH4( GET_YEAR(), .TRUE., 
     &                        A3090S, A0030S, A0030N, A3090N )
      ENDIF

      !=================================================================
      ! Do tagged CO chemistry -- Put everything within a large 
      ! DO-loop over all grid boxes to facilitate parallelization
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED                                             ) 
!$OMP+PRIVATE( I, J, L, N, STTCO, GCO, DENS, CH4RATE, CO_CH4, CH4 )
!$OMP+PRIVATE( ALPHA_CH4, KRATE, CO_ISOP, CO_CH3OH, ALPHA_ISOP    )
!$OMP+PRIVATE( CO_MONO, ALPHA_MONO, CO_ACET, ALPHA_ACET, CORATE   )
!$OMP+PRIVATE( CO_OH, PCO, YMID, CO_VOC                           )
      DO L = 1, LLPAR
      DO J = 1, JJPAR

         ! Latitude of grid box
         YMID = GET_YMID( J )

      DO I = 1, IIPAR

         !==============================================================
         ! (0) Define useful quantities
         !==============================================================

         ! STTCO [molec CO/cm3/kg CO] converts [kg CO] --> [molec CO/cm3]
         ! kg CO/box * box/cm3 * mole/0.028 kg CO * Avog.#/mole
         STTCO = 1d0 / AIRVOL(I,J,L) / 1d6 / FMOL_CO * 6.023d23

         ! GCO is CO concentration in [molec CO/cm3]
         GCO   = STT(I,J,L,1) * STTCO

         ! DENS is the number density of air [molec air/cm3]
         DENS  = AD(I,J,L) * 1000.d0 / BOXVL(I,J,L) * 6.023d23  / WTAIR

         !==============================================================
         ! (1a) Production of CO by reaction with CH4
         !==============================================================

         ! Initialize
         CO_CH4 = 0d0

         ! Test level for stratosphere or troposphere
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

            !===========================================================
            ! (1a-1) Production of CO from CH4 in the stratosphere 
            !===========================================================

            ! Call GET_PCO_LCO_STRAT to get the P(CO) rate from CH4
            CH4RATE = GET_PCO_LCO_STRAT( .TRUE., I, J, L )
            
            ! Convert units of CH4RATE from [v/v/s] to [molec CO/cm3]
            CO_CH4 = CH4RATE * DTCHEM * DENS

         ELSE
            
            !===========================================================
            ! (1a-2) Production of CO from CH4 in the troposphere 
            !===========================================================

            ! CH4 concentration [ppbv] for the given latitude band 
            ! (bmy, 1/2/01)
            CH4 = A3090S
            IF ( YMID >= -30.0 .and. YMID < 0.0  ) CH4 = A0030S
            IF ( YMID >=   0.0 .and. YMID < 30.0 ) CH4 = A0030N
            IF ( YMID >=  30.0                   ) CH4 = A3090N
            
            ! Convert CH4 from [ppbv] to [molec CH4/cm3]
            CH4 = CH4 * 1d-9 * DENS

            ! Yield of CO from CH4: estimated to be 95-100% (acf)
            ALPHA_CH4 = 1d0

            ! Calculate updated rate constant [s-1] (bnd, bmy, 1/2/01)
            KRATE = 2.45D-12 * EXP( -1775.d0 / T(I,J,L) ) 
 
            ! Production of CO from CH4 = alpha * k * [CH4] * [OH] * dt 
            ! Units are [molec CO/cm3]
            CO_CH4 = ALPHA_CH4 * KRATE * CH4 * OH(I,J,L) * DTCHEM 

         ENDIF

         ! Check CO_CH4 for NaN or Infinity
         CALL CHECK_VALUE( CO_CH4,  (/ I, J, L, 0 /),         
     &                    'CO_CH4', 'STOP at tagged_co_mod:1' )

         !==============================================================
         ! (1b) Production of CO from ISOPRENE and METHANOL (CH3OH)
         !==============================================================

         ! Initialize
         CO_ISOP  = 0d0
         CO_CH3OH = 0d0

         ! Isoprene is emitted only into the surface layer 
         IF ( L == 1 ) THEN 

            !===========================================================
            ! Yield of CO from ISOP: 30%, from Miyoshi et al., 1994.
            ! They estimate globally 105 Tg C/yr of CO is produced 
            ! from isoprene oxidation. 
            !-----------------------------------------------------------
            ! We need to scale the Isoprene flux to get the CH3OH 
            ! (methanol) flux.  Currently, the annual isoprene flux in 
            ! GEOS-CHEM is ~ 397 Tg C.
            !
            ! Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH 
            ! oxidation based on Singh et al. 2000 [JGR 105, 3795-3805] 
            ! who estimate a global methanol source of 122 Tg yr-1, of 
            ! which most (75 Tg yr-1) is "primary biogenic".  He also 
            ! recommends for now that the CO flux from CH3OH oxidation 
            ! be scaled to monthly mean isoprene flux.
            !
            ! To get CO from METHANOL oxidation, we must therefore 
            ! multiply the ISOPRENE flux by the following scale factor:
            !  ( 100 Tg CO / 397 Tg C ) * ( 12 g C/mole / 28 g CO/mole )
            !-----------------------------------------------------------
            ! We now call GET_ALPHA_ISOP to get the yield factor of
            ! CO produced from isoprene, as a function of NOx, or
            ! as a constant. (bnd, bmy, 6/14/01)
            !===========================================================

            ! Get CO yield from isoprene
            IF ( ALPHA_ISOP_FROM_NOX ) THEN
               ALPHA_ISOP = GET_ALPHA_ISOP( .TRUE., BNOX(I,J,L) )         
            ELSE
               ALPHA_ISOP = GET_ALPHA_ISOP( .FALSE. )
            ENDIF

            ! P(CO) from Isoprene Flux = ALPHA_ISOP * Flux(ISOP)
            ! Convert from [molec ISOP/box] to [molec CO/cm3]
            !
            ! Units of SUMISOPCO are [atoms C/box/time step]. 
            ! Division by 5 is necessary to convert to 
            ! [molec ISOP/box/timestep].
            !
            ! Units of ALPHA_ISOP are [molec CO/molec ISOP]
            ! Units of CO_ISOP are [molec CO/cm3]
            CO_ISOP = SUMISOPCO(I,J) / BOXVL(I,J,L) / 5.d0 * ALPHA_ISOP

            ! P(CO) from CH3OH is scaled to Isoprene Flux (see above)
            ! Units are [molec CO/cm3]
            CO_CH3OH = ( SUMISOPCO(I,J) / BOXVL(I,J,L) ) * 
     &                 ( 100d0          / 397d0        ) * 
     &                 ( 12d0           / 28d0         )
 
            ! Zero SUMISOPCO and SUMCH3OHCO for the next emission step
            SUMISOPCO(I,J)  = 0d0
            SUMCH3OHCO(I,J) = 0d0

            ! Check CO_ISOP for NaN or Infinity
            CALL CHECK_VALUE( CO_ISOP,  (/ I, J, L, 0 /),
     &                       'CO_ISOP', 'STOP at tagged_co_mod:2' )

            ! Check CO_CH4 for NaN or Infinity
            CALL CHECK_VALUE( CO_CH3OH,  (/ I, J, L, 0 /),
     &                       'CO_CH3OH', 'STOP at tagged_co_mod:3' )

         ENDIF

         !==============================================================
         ! (1c) Production of CO from MONOTERPENE oxidation
         !==============================================================

         ! Initialize
         CO_MONO = 0.d0

         ! Monoterpenes are emitted only into the surface layer
         IF ( L == 1 ) THEN

            !===========================================================
            ! Assume the production of CO from monoterpenes is 
            ! instantaneous even though the lifetime of intermediate 
            ! species may be on the order of hours or days.  This 
            ! assumption will likely cause CO from monoterpene 
            ! oxidation to be too high in the box in which the 
            ! monoterpene is emitted.
            !-----------------------------------------------------------
            ! The CO yield here is taken from:
            !   Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
            !   Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 
            !     (1998)
            !
            ! Hatakeyama:  "The ultimate yield of CO from the 
            !   tropospheric oxidation of terpenes (including both O3 
            !   and OH reactions) was estimated to be 20% on the carbon 
            !   number basis."  They studied ALPHA- & BETA-pinene.
            !
            ! Vinckier  :  "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
            !-----------------------------------------------------------
            ! Calculate source of CO per time step from monoterpene 
            ! flux (assume lifetime very short) using the C number basis:
            !
            !   CO [molec CO/cm3] = Flux [atoms C from MONO/box] /
            !                       Grid Box Volume [cm^-3]       *
            !                       ALPHA_MONO 
            !
            ! where ALPHA_MONO = 0.2 as explained above.
            !===========================================================

            ! Yield of CO from MONOTERPENES: 20% (see above)
            ALPHA_MONO = 0.20d0

            ! P(CO) from Monoterpene Flux =  alpha * Flux(Mono)
            ! Units are [molec CO/cm3]
            CO_MONO = ( SUMMONOCO(I,J) / BOXVL(I,J,L) ) * ALPHA_MONO

            ! Zero SUMMONOCO for the next emission step
            SUMMONOCO(I,J) = 0d0

            ! Check CO_MONO for NaN or Infinity
            CALL CHECK_VALUE( CO_MONO,  (/ I, J, L, 0 /),
     &                       'CO_MONO', 'STOP at tagged_co_mod:4' )

         ENDIF

         !==============================================================
         ! (1d) Production of CO from oxidation of ACETONE 
         !
         ! ALPHA_ACET = 2/3 to get a yield for CO.  This accounts 
         ! for acetone loss from reaction with OH And photolysis.  
         ! The acetone sources taken into account are:
         ! 
         ! (a) Primary emissions of acetone from biogenic sources
         ! (b) Secondary production of acetone from monoterpene 
         !      oxidation
         ! (c) Secondary production of acetone from ALK4 and 
         !      propane oxidation
         ! (d) Direct emissions of acetone from biomass burning and 
         !      fossil fuels
         ! (e) direct emissions from ocean
         !
         ! Calculate source of CO per time step from biogenic acetone 
         ! # molec CO/cc = ALPHA * ACET Emission Rate * dt
         !==============================================================

         ! Initialize
         CO_ACET = 0.d0

         ! Biogenic acetone sources are emitted only into the surface layer
         IF ( L == 1 ) THEN

            ! Yield of CO from ACETONE: 2/3 (see above)
            ALPHA_ACET = 2.D0 / 3.D0
            
            ! Units are [molec CO/cc]
            CO_ACET = SUMACETCO(I,J) / BOXVL(I,J,L) * ALPHA_ACET

            ! Zero SUMACETCO for the next emission step
            SUMACETCO(I,J) = 0d0

            ! Check CO_ACET for NaN or Infinity
            CALL CHECK_VALUE( CO_ACET,  (/I, J, L, 0 /),
     &                       'CO_ACET', 'STOP at tagged_co_mod:5' )

         ENDIF

         !==============================================================
         ! (1e) Add production of CO into the following tagged tracers:
         !
         ! (a) Tracer #12: CO produced from CH4
         ! (b) Tracer #14: CO produced from ISOPRENE
         ! (c) Tracer #15: CO produced from MONOTERPENES
         ! (d) Tracer #16: CO produced from METHANOL (CH3OH)
         ! (e) Tracer #17: CO produced from ACETONE
         !
         ! %%% NOTE: If you are modifying the tagged CO simulation, 
         ! %%% and your simulation has less than 12 tracers, then 
         ! %%% then comment out this section.  If you don't you can
         ! %%% get an array-out-of-bounds error (bmy, 6/11/08)
         !==============================================================
         IF ( LSPLIT ) THEN
            STT(I,J,L,12) = STT(I,J,L,12) + CO_CH4   / STTCO
            STT(I,J,L,14) = STT(I,J,L,14) + CO_ISOP  / STTCO 
            STT(I,J,L,15) = STT(I,J,L,15) + CO_MONO  / STTCO
            STT(I,J,L,16) = STT(I,J,L,16) + CO_CH3OH / STTCO
            STT(I,J,L,17) = STT(I,J,L,17) + CO_ACET  / STTCO
         ENDIF

         !==============================================================
         ! (2a) Loss of CO due to chemical reaction w/ OH
         !==============================================================

         ! Select out tropospheric or stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

            !===========================================================
            ! (2a-1) Stratospheric loss of CO due to chemical rxn w/ OH
            !===========================================================

            ! Get the L(CO) rate in the stratosphere in [s-1]
            CORATE = GET_PCO_LCO_STRAT( .FALSE., I, J, L )

            ! CO_OH is the fraction of CO lost to OH [unitless]
            CO_OH = CORATE * DTCHEM

            ! Check CO_ACET for NaN or Infinity
            CALL CHECK_VALUE( CO_OH,  (/ I, J, L, 0 /),
     &                       'CO_OH', 'STOP at tagged_co_mod:6' )

            ! Handle strat loss by OH for regional CO tracers
            IF ( LSPLIT ) THEN

               ! Loop over regional CO tracers
               DO N = 2, N_TRACERS

                  ! Loss
                  STT(I,J,L,N) = STT(I,J,L,N) * ( 1d0 - CO_OH )

                  ! STT shouldn't be less than zero
                  IF ( STT(I,J,L,N) < 0d0 ) STT(I,J,L,N) = 0d0

                  ! Error check
                  CALL CHECK_VALUE( STT(I,J,L,N), (/ I, J, L, N /),
     &                              'STT','STOP at tagged_co_mod.f:7' )

                  ! ND65 diagnostic -- loss of CO by OH for "tagged" tracers
                  IF ( ND65 > 0 .and. L <= LD65 ) THEN
                     AD65(I,J,L,N) = AD65(I,J,L,N) +
     &                               ( CORATE * STT(I,J,L,N) * STTCO )
                  ENDIF
               ENDDO
            ENDIF

            ! CO_OH above is just the fraction of CO lost by OH.  Here
            ! we multiply it by GCO (the initial value of STT in molec/cm3)
            ! to convert it to an amount of CO lost by OH [molec/cm3]
            ! (bmy, 2/19/02)
            CO_OH = GCO * CO_OH

         ELSE
   
            !===========================================================
            ! (2a-2) Tropospheric loss of CO due to chemical rxn w/ OH
            !
            !  DECAY RATE
            !  The decay rate (KRATE) is calculated by:
            !
            !     OH + CO -> products (JPL '97)
            !     k = (1 + 0.6Patm) * 1.5E-13
            !
            !  KRATE has units of [ molec^2 CO / cm6 / s ]^-1, 
            !  since this is a 2-body reaction.
            !===========================================================

            ! Pressure at the center of sigma level L,
            ! expressed as a fraction of surface pressure in [Pa]
            PCO = GET_PCENTER(I,J,L) / 1.01325d3

            ! Decay rate
            KRATE = ( 1.d0 + ( 0.6d0 * PCO ) ) * 1.5d-13

            ! CO_OH = Tropospheric loss of CO by OH [molec/cm3]
            CO_OH = KRATE * GCO * OH(I,J,L) * DTCHEM

            ! Handle trop loss by OH for regional CO tracers
            IF ( LSPLIT ) THEN

               ! Loop over regional CO tracers
               DO N = 2, N_TRACERS

                  ! Use tropospheric rate constant 
                  STT(I,J,L,N) = STT(I,J,L,N) *
     &                 ( 1d0 - KRATE * OH(I,J,L) * DTCHEM )

                  ! Error check
                  CALL CHECK_VALUE( STT(I,J,L,N), (/ I, J, L, N /),
     &                              'STT','STOP at tagged_co_mod.f:8' )

                  ! ND65 diagnostic -- loss of CO by OH for "tagged" tracers
                  IF ( ND65 > 0 .and. L <= LD65 ) THEN
                     AD65(I,J,L,N) = AD65(I,J,L,N) + 
     &                    ( KRATE * OH(I,J,L)*STT(I,J,L,N) * STTCO )
                  ENDIF
               ENDDO
            ENDIF

         ENDIF
         
         !==============================================================
         ! Save the total chemical production from various sources
         ! into the total CO tracer STT(I,J,L,1)
         !==============================================================

         ! GCO is the total CO before chemistry was applied [molec CO/cm3]
         ! Add to GCO the sources and sinks listed above
         CO_VOC = CO_MONO + CO_ACET + CO_CH3OH + CO_ISOP
         
         IF ( LADJ .AND. LADJ_EMS .AND. NNEMS .GT. 1) THEN
            IF ( LDEVOC ) THEN  !(zhe 11/28/10)
            !----------------------------------------------------------
            ! COSF now obsolete (zhe, dkh, 01/17/12, adj32_017) 
            ! OLD:
            ! GCO = GCO + EMS_SF(I,J,M,ADCOEMS) * CO_VOC * COSF(I,J,1) 
            !+ EMS_SF(I,J,M,ADCOVOX) * CO_CH4 * COSF(I,J,2) - CO_OH
            !ELSE    
            !  GCO = GCO + EMS_SF(I,J,M,ADCOVOX) * (CO_CH4 + CO_VOC) * 
            !   COSF(I,J,2) - CO_OH
            ! NEW:
                GCO = GCO + EMS_SF(I,J,M,ADCOEMS) * CO_VOC
     &                    + EMS_SF(I,J,M,ADCOVOX) * CO_CH4 - CO_OH
             ELSE
                GCO = GCO + EMS_SF(I,J,M,ADCOVOX) * (CO_CH4 + CO_VOC)
     &                                                       - CO_OH
            !----------------------------------------------------------
            ENDIF
         ELSE
            GCO = GCO + CO_CH4 + CO_VOC - CO_OH
         ENDIF
         
         IF ( ND29 > 0 ) THEN  !(zhe 11/28/10)
            AD29(I,J,7) = AD29(I,J,7) + CO_VOC / STTCO
            AD29(I,J,8) = AD29(I,J,8) + CO_CH4 / STTCO
         ENDIF
         
         ! Removed (zhej, dkh, 01/17/12, adj32_017) 
         !Save the emission for the backward run
         !(zhe 11/28/10)
         !CALL SET_CO_CH4(I, J, L, TIME_STEPS, CO_CH4/STTCO)

         ! Convert net CO from [molec CO/cm3] to [kg] and store in STT
         STT(I,J,L,1) = GCO / STTCO
         
         IF(I == IFD .AND. J == JFD .AND. L==1 .and. LPRINTFD) THEN
          PRINT*, 'STTCO:', STTCO
          PRINT*, 'CO_OH/STTCO=', CO_OH/STTCO
          PRINT*, 'CO_VOC:', CO_VOC / STTCO
          PRINT*, 'CO_CH4:', CO_CH4 / STTCO
         ENDIF

         !==============================================================
         ! Archive ND65 diagnostics -- Production & Loss of CO
         ! Also save P(CO) from CH3OH and MONOTERPENES (bmy, 1/2/01)
         !==============================================================
         IF ( ND65 > 0 .and. L <= LD65 ) THEN

            ! Loss of CO by OH (global) [s-1]
            N             = 1
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_OH / DTCHEM )

            ! Production of CO from Isoprene [molec CO/cm3/s]
            N             = N_TRACERS + 1
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_ISOP / DTCHEM )

            ! Production of CO from CH4 [molec CO/cm3/s]
            N             = N_TRACERS + 2
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_CH4 / DTCHEM )

            ! Production of CO from CH3OH [molec CO/cm3/s]            
            N             = N_TRACERS + 3
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_CH3OH / DTCHEM )

            ! Production of CO from MONO [molec CO/cm3/s]
            N             = N_TRACERS + 4
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_MONO / DTCHEM )    

            ! Production of CO from ACET [molec CO/cm3/s]
            N             = N_TRACERS + 5 
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_ACET / DTCHEM )
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(LPRINTFD) THEN
          PRINT*, 'STT:', STT(IFD,JFD,1,1)
      ENDIF

      ! Obsolete (zhej, dkh, 01/17/12, adj32_017) 
      ! TIME_STEPS = TIME_STEPS + 1  !ZHE

      ! Return to calling program
      END SUBROUTINE CHEM_TAGGED_CO

!------------------------------------------------------------------------------

      FUNCTION GET_ALPHA_ISOP( FROM_NOX, NOX ) RESULT( ALPHA_ISOP )
!
!******************************************************************************
!  Function GET_ALPHA_ISOP returns the CO yield from Isoprene (ALPHA_ISOP) 
!  either as a function of NOx or as a constant. (bnd, bmy, 6/13/01. 7/20/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FROM_NOX (LOGICAL) : If =T, will take ALPHA_ISOP as f(NOx)
!                            If =F, will set ALPHA_ISOP to a constant
!  (2 ) NOX      (TYPE (XPLEX) ) : NOx concentration in ppbv
! 
!  NOTES:
!  (1 ) Now make NOx an optional argument (bmy, 8/28/01)
!  (2 ) Now reference ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (3 ) Updated comments (bmy, 7/20/04)
!******************************************************************************
!      
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      LOGICAL, INTENT(IN)           :: FROM_NOX
      TYPE (XPLEX),  INTENT(IN), OPTIONAL :: NOX

      ! Function Value
      TYPE (XPLEX)                        :: ALPHA_ISOP

      !=================================================================
      ! GET_ALPHA_ISOP begins here!
      !
      ! If FROM_NOX = T, then we are taking the CO yield from Isoprene 
      ! as a function of NOx.  Otherwise, we set it to a constant.
      !=================================================================
      IF ( FROM_NOX ) THEN

         ! Make sure NOX is passed if FROM_NOX = TRUE
         IF ( PRESENT( NOX ) ) THEN

            ! The CO yield here is taken from acf's calculation 
            ! (from group meeting (1-20-99).  Assumming linearity, 
            ! ALPHA=0.8*[NOx]+0.6, with an upper and lower limit of 2.1 
            ! and 0.8, respectively. [NOx] is concentration in ppbv !!
            ALPHA_ISOP = ( 0.8d0 * NOX ) + 0.6D0

            ! Force lower & upper limits
            IF ( NOX < 0.5d0 ) ALPHA_ISOP = 0.8d0
            IF ( NOX > 1.8d0 ) ALPHA_ISOP = 2.1d0
         
         ELSE

            ! Stop w/ error message
            CALL ERROR_STOP( 'NOx argument not passed!', 
     &                       'GET_ALPHA_ISOP (tagged_co_mod.f)' )
           
         ENDIF
   
      ELSE

         ! Otherwise, use a 30% yield from Miyoshi et al., 1994.
         ! They estimate globally 105 Tg C/yr of CO is produced
         ! from isoprene oxidation.
         ! ALPHA_ISOP = (0.3 molec CO/atoms C) x (5 atoms C/molec ISOP)
         ALPHA_ISOP = 1.5d0

      ENDIF
      
      ! Return to calling program
      END FUNCTION GET_ALPHA_ISOP

!------------------------------------------------------------------------------
 
      SUBROUTINE READ_PCO_LCO_STRAT( THISMONTH ) 
!
!******************************************************************************
!  Subroutine READ_PCO_LCO_STRAT reads production and destruction
!  rates for CO in the stratosphere. (bnd, bmy, 9/13/00, 10/3/05)
!
!  NOTES:
!  (1 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (2 ) Added to module "tagged_co_mod.f" (bmy, 6/13/01)
!  (3 ) ARRAY needs to be of size (IGLOB,JGLOB).  Use TRANSFER_ZONAL from
!        "transfer_mod.f" to cast data from TYPE (XPLEX) to COMPLEX*16, and also to
!        copy data to an array of size (JJPAR,LLPAR).  Use 3 arguments (M/D/Y)
!        in call to GET_TAU0.  Use JGLOB,LGLOB in call to READ_BPCH2. 
!        (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (5 ) Updated comments (bmy, 2/15/02)
!  (6 ) Update FILENAME so that it looks in the "pco_lco_200203" subdirectory
!        of DATA_DIR. (bnd, bmy, 6/30/03)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE TRANSFER_MOD,  ONLY : TRANSFER_ZONAL

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      CHARACTER(LEN=255)  :: FILENAME
      TYPE (XPLEX)              :: ARRAY(1,JGLOB,LGLOB) 
      TYPE (XPLEX)              :: XTAU

      !=================================================================
      ! READ_PCO_LCO_STRAT begins here!
      !=================================================================

      ! Initialize some variables
      ARRAY    = 0e0
      CO_PRODS = 0d0
      CO_LOSSS = 0d0

      ! TAU value at the beginning of this month
      ! Use "generic" year 1985
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      !=================================================================
      ! Read in CO production rates [v/v/s]
      !=================================================================

      ! Construct filename
      FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COprod.' //
     &           GET_NAME_EXT()   // '.' // GET_RES_EXT()
!     &           'geos4'          // '.' // GET_RES_EXT()
!      print*, 'dkh hack: use GEOS_4 P/L fields'
!      print*, 'dkh hack: use GEOS_4 P/L fields'
!      print*, 'dkh hack: use GEOS_4 P/L fields'
!      print*, 'dkh hack: use GEOS_4 P/L fields'
      
      ! Read P(CO)
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,     
     &                 XTAU,      1,         JGLOB,     
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (JJPAR,LLPAR)
      CALL TRANSFER_ZONAL( ARRAY(1,:,:), CO_PRODS )

      !=================================================================
      ! Read in CO loss rates [s^-1]
      !=================================================================

      ! Construct filename
      FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COloss.' //
     &           GET_NAME_EXT()   // '.' // GET_RES_EXT()
!     &           'geos4'     // '.' // GET_RES_EXT()
!      print*, 'dkh hack: use GEOS_4 P/L fields'
!      print*, 'dkh hack: use GEOS_4 P/L fields'
!      print*, 'dkh hack: use GEOS_4 P/L fields'
!      print*, 'dkh hack: use GEOS_4 P/L fields'

      ! Read L(CO)
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 10,    
     &                 XTAU,      1,         JGLOB,     
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (JJPAR,LLPAR)
      CALL TRANSFER_ZONAL( ARRAY(1,:,:), CO_LOSSS )

      ! Return to calling program
      END SUBROUTINE READ_PCO_LCO_STRAT

!------------------------------------------------------------------------------

      FUNCTION GET_PCO_LCO_STRAT( IS_PROD, I, J, L ) RESULT( RATE )
!
!******************************************************************************
!  Function GET_CO_STRAT_RATE uses production and loss rates for CO to 
!  calculate net production of CO in the stratosphere.  The purpose of this 
!  SR is to prevent high CO concentrations from building up in the 
!  stratosphere; in these layers only transport is simulated (i.e., no 
!  chemistry).  For a long simulation, a buildup of high concentrations 
!  could occur causing the stratosphere to become a significant source of CO. 
!  (bnd, bmy, 6/13/01, 7/20/04)
!
!  Arguments as Input: 
!  ============================================================================
!  (1  ) IS_PROD (LOGICAL) : If =T, then return CO production rate [v/v/s]
!                            If =F, then return CO loss rate [s^-1]
!  (2-5) I, J, L (INTEGER) : Grid box indices (lon,lat,alt)
!  
!  Arguments as Output: 
!  ============================================================================
!  (6  ) RATE    (TYPE (XPLEX) ) : CO production [v/v/s] or loss rate [s-1]
!
!  NOTES:
!  (1 ) Production (mixing ratio/sec) and loss (1/sec) rates provided by 
!        Dylan Jones.  Only production by CH4+OH and destruction by CO+OH 
!        are considered.
!  (2 ) The annual mean tropopause is stored in the LPAUSE array 
!       (from header file "CMN").  LPAUSE is defined such that: 
!       Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!               LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!  (3 ) LPAUSE_MIN = minimun tropopause height.  Start L-loop from the 
!        lowest stratospheric level!
!  (4 ) Added to module "tagged_co_mod.f" (bmy, 6/18/01)
!  (5 ) Updated comments (bmy, 2/19/02)
!  (6 ) Removed reference to CMN, it's not needed (bmy, 7/20/04)
!******************************************************************************
!
#     include "CMN_SIZE"       ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN)  :: IS_PROD
      INTEGER, INTENT(IN)  :: I, J, L

      ! Function return value
      TYPE (XPLEX)               :: RATE
 
      !=================================================================
      ! GET_PCO_LCO_STRAT begins here!
      !
      ! Pick P(CO) or L(CO) depending on IS_PROD
      !=================================================================
      IF ( IS_PROD ) THEN
         RATE = CO_PRODS(J,L)   ! P(CO) from CH4 + OH in [v/v/s]
      ELSE
         RATE = CO_LOSSS(J,L)   ! L(CO) from CO + OH  in [s^-1]
      ENDIF
    
      ! Return to calling program
      END FUNCTION GET_PCO_LCO_STRAT

!-----------------------------------------------------------------------------

      SUBROUTINE READ_ACETONE( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_ACETONE reads in biogenic acetone emissions from
!  a binary punch file. (bdf, bnd, bmy, 6/19/01, 10/3/05)
!
!  Arguments as Input
!  =========================================================================== 
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Eliminate variables that aren't used anymore.  Updated comments
!        and made some cosmetic changes. (bnd, bmy, 6/14/01)
!  (2 ) Added to "tagged_co_mod.f" (bmy, 6/14/01)
!  (3 ) Now read acetone file from DATA_DIR/tagged_CO_200106 (bmy, 6/19/01)
!  (4 ) ARRAY needs to be of size (IGLOB,JGLOB).  Use TRANSFER_2D from 
!        "transfer_mod.f" to cast data from TYPE (XPLEX) to COMPLEX*16, and also to 
!        copy data to an array of size (IIPAR,JJPAR).  Use 3 arguments (M/D/Y)
!        in call to GET_TAU0.  Use IGLOB,JGLOB in call to READ_BPCH2.
!        Added array TEMP(IIPAR,JJPAR).  Updated comments. (bmy, 9/28/01)
!  (5 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (6 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (7 ) Now reads data from both GEOS and GCAP grids (bmy, 8/16/05)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments 
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      TYPE (XPLEX)              :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)              :: TEMP(IIPAR,JJPAR)
      TYPE (XPLEX)              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME
     
      !=================================================================
      ! READ_ACETONE begins here!
      !=================================================================

      ! Name of file with acetone data
      FILENAME = TRIM( DATA_DIR )            // 
     &           'tagged_CO_200106/acetone.' // GET_NAME_EXT_2D() //
     &           '.'                         // GET_RES_EXT()

      ! Echo to stdout
      WRITE( 6, '(a)' ) 'READING ', TRIM( FILENAME )
 
      !  Initialize some variables      
      EMACET   = 0d0
       
      ! TAU value at the beginning of this month for 1994
      XTAU = GET_TAU0( THISMONTH, 1, 1994 )

      !=================================================================
      ! Read direct biogenic acetone emissions [atoms C/cm2/s]
      !=================================================================

      ! Initialize ARRAY
      ARRAY = 0e0

      ! Read data
      CALL READ_BPCH2( FILENAME, 'EMISACET', 6, XTAU,
     &                 IIPAR,     JJPAR,     1, ARRAY )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), EMACET )

      !=================================================================
      ! Read acetone from grasslands [atoms C/cm2/s]
      !=================================================================

      ! Initialize ARRAY
      ARRAY = 0e0

      ! Read data
      CALL READ_BPCH2( FILENAME, 'EMISACET', 7, XTAU,
     &                 IIPAR,     JJPAR,     1, ARRAY )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and add
      CALL TRANSFER_2D( ARRAY(:,:,1), TEMP )

      ! Add acetone from grasslands to direct biogenic emissions
      EMACET(:,:) = EMACET(:,:) + TEMP(:,:)

      ! Return to calling program
      END SUBROUTINE READ_ACETONE

!------------------------------------------------------------------------------
      
      FUNCTION GET_SCALE_GROUP( ) RESULT( CURRENT_GROUP )
!
!********************************************************************************
! Subroutine GET_SCALE_GROUP determines which predifined scaling index corresponds
! to the current time and location  (dkh, 12/02/04)
!
! NOTES
! (1 ) CURRENT_GROUP is currently only a function of TAU
! (2 ) Get rid of I,J as argument. (dkh, 03/28/05)
!
!********************************************************************************

      ! Reference to f90 modules
      USE TIME_MOD, ONLY : GET_TAU, GET_TAUe, GET_TAUb, GET_MONTH
      USE ADJ_ARRAYS_MOD, ONLY: MMSCL

#     include "CMN_SIZE" ! Size stuff

      ! Arguments
      INTEGER      :: I, J

      ! Local Variables
      TYPE (XPLEX)       :: TOTAL_HR, CURRENT_HR, GROUP_LENGTH
      TYPE (XPLEX)       :: TAU, TAUe, TAUb

      ! Function variable
      INTEGER      :: CURRENT_GROUP
      LOGICAL, SAVE :: MONTHLY = .TRUE.
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
            GROUP_LENGTH =( TOTAL_HR / MMSCL )

            ! Index is the current time divided by the group length, plus one
            CURRENT_GROUP = ( CURRENT_HR / GROUP_LENGTH ) + 1

         ENDIF

      ENDIF

      END FUNCTION GET_SCALE_GROUP

! Obsolete (zhej, dkh, 01/17/12, adj32_017) 
!!-----------------------------------------------------------------------------------------
!
!      SUBROUTINE GET_OPT_CO( XTAU )
!!
!!******************************************************************************
!!  Subroutine GET_OPT_CO reads CO source scaling factors from binary punch files
!!  (zhe zhe 1/19/11)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD
!      USE ADJ_ARRAYS_MOD,   ONLY : NOR
!
!      USE MYTYPE
!      USE COMPLEXIFY
!      IMPLICIT NONE
!
!#     include "CMN_SIZE"    ! Size parameters
!
!      ! Arguments
!      TYPE (XPLEX), INTENT(IN)  :: XTAU
!      INTEGER             :: I, J
!
!      ! Local variables
!      CHARACTER(LEN=255)   :: FILENAME
!
!      !=================================================================
!      ! GET_OPT_CO begins here!
!      !=================================================================
!#IF DEFINED(NESTED_CH )
!      ! Filename
!        FILENAME = TRIM( 'COSF_' ) // GET_RES_EXT()
!
!      ! Echo some information to the standard output
!        WRITE( 6, 110 ) TRIM( FILENAME )
! 110    FORMAT( '     - GET_OPT_CO: Reading COSF from: ', a )
!      
!      ! Read data from the binary punch file
!        CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,     
!     &           XTAU,      IGLOB,     JGLOB,      
!     &           2,     COSF,     QUIET=.TRUE. )
!        PRINT*, 'XTAU=', XTAU
!          
!      ! Nested observation region
!        COSF(NOR(1):NOR(2),NOR(3):NOR(4),:) = 1.0
!
!#ELSE
!          PRINT*, 'Using the original emission inventory'
!#ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE GET_OPT_CO
!!------------------------------------------------------------------------------
!
!      SUBROUTINE SET_CO_CH4( I, J, L, N, VALUE ) 
!!
!!******************************************************************************
!!  Subroutine SET_CO_CH4 save the CO_CH4 in the forward run (zhe 11/28/10)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) I     (INTEGER) : Index for 1st dimension of array
!!  (2 ) J     (INTEGER) : Index for 1st dimension of array
!!  (3 ) L     (INTEGER) : Index for 1st dimension of array
!!  (4 ) N     (INTEGER) : Index for time steps
!!  (5 ) VALUE (TYPE (XPLEX) ) : Value to store in (I,J,N)th element of array
!!  NOTES:
!!******************************************************************************
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: I, J, L, N
!      TYPE (XPLEX),  INTENT(IN) :: VALUE
!      INTEGER             :: H
!
!      !=================================================================
!      ! SET_CO_CH4 begins here!
!      !=================================================================
!      H = MOD(N,24)
!      IF (H .EQ. 0) H = 24
!      CO_CH4(I, J, L, H) = VALUE
!
!      ! Return to calling program
!      END SUBROUTINE SET_CO_CH4
!
!!-----------------------------------------------------------------------------
!
!      FUNCTION GET_CO_CH4( I, J, L, N) RESULT( VALUE )
!!
!!******************************************************************************
!!  Subroutine GET_CO_VOC gets a value for CO_VOC (zhe 11/28/10)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) I     (INTEGER) : Index for 1st dimension of array
!!  (2 ) J     (INTEGER) : Index for 1st dimension of array
!!  (3 ) L     (INTEGER) : Index for 1st dimension of array
!!  (4 ) N     (INTEGER) : Index for time steps
!!  (6 ) VALUE (TYPE (XPLEX) ) : Value to store in (I,J,N)th element of array
!!
!!  NOTES:
!!*****************************************************************************
!      ! Arguments
!      INTEGER, INTENT(IN) :: I, J, L, N
!
!      ! Function value
!      TYPE (XPLEX)              :: VALUE
!      INTEGER             :: H
!
!      !=================================================================
!      ! GET_CO_VOCH4 begins here!
!      !=================================================================
!      
!      H = MOD(N,24)
!      IF (H .EQ. 0) H = 24
!      VALUE = CO_CH4(I, J, L, H)
!
!      ! Return to calling program
!      END FUNCTION GET_CO_CH4
!      
!!-----------------------------------------------------------------------------------------

      SUBROUTINE INIT_TAGGED_CO
!
!******************************************************************************
!  Subroutine INIT_TAGGED_CO allocates memory to module arrays.
!  (bmy, 7/19/00, 9/18/07)
!
!  NOTES:
!  (1 ) Added ISOP96, MONO96, CH3OH96 for GEOS-3 (bnd, bmy, 6/14/01)
!  (2 ) Removed ISOP96, MONO96, CH3OH96 for GEOS-3, since the new GEOS-3
!        fields make these no longer necessary (bmy, 8/21/09)
!  (3 ) Now allocate BB_REGION, FF_REGION as (IIPAR,JJPAR) (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (5 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (6 ) Now remove IJLOOP_CO (bmy, 7/20/04)
!  (7 ) Now public. Now references ITS_A_H2HD_SIM from "tracer_mod.f". 
!        Allocate needed variables if H2/HD simulation (phs, 9/18/07)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE TRACER_MOD, ONLY : ITS_A_H2HD_SIM

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS, I, J, IJLOOP
      
      !=================================================================
      ! INIT_TAGGED_CO begins here!
      !=================================================================

      ! Allocate SUMACETCO -- array for CO from isoprene
      ALLOCATE( EMACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMACET' )         
      EMACET = 0d0 

      ! Allocate and initialize IJLOOP_CO -- 1-D array index
      ALLOCATE( CO_PRODS( JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_PRODS' )    
      CO_PRODS = 0d0

      ! Allocate and initialize IJLOOP_CO -- 1-D array index
      ALLOCATE( CO_LOSSS( JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_LOSSS' )    
      CO_LOSSS = 0d0

      ! Allocate only what is needed
      IF ( ITS_A_H2HD_SIM() ) RETURN

      ! Allocate BB_REGION -- array for biomass burning regions
      ALLOCATE( BB_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BB_REGION' )         
      BB_REGION = 0

      ! Allocate FF_REGION -- array for fossil fuel regions
      ALLOCATE( FF_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FF_REGION' )
      FF_REGION = 0
      
      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMISOPCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMISOPCO' )         
      SUMISOPCO = 0d0 

      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMMONOCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMMONOCO' )         
      SUMMONOCO = 0d0 

      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMCH3OHCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMCH3OHCO' )         
      SUMCH3OHCO = 0d0 

      ! Allocate SUMACETCO -- array for CO from isoprene
      ALLOCATE( SUMACETCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMACETCO' )         
      SUMACETCO = 0d0
      
      ! Obsolete (zhej, dkh, 01/17/12, adj32_017) 
      !ALLOCATE( COSF( IIPAR, JJPAR, 2 ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSF' )
      !COSF(:,:,:) = 1D0
      !
      !ALLOCATE( CO_CH4( IIPAR, JJPAR, LLPAR, 24), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_CH4' )
      !CO_CH4 = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_TAGGED_CO

!------------------------------------------------------------------------------
  
      SUBROUTINE CLEANUP_TAGGED_CO
!
!******************************************************************************
!  Subroutine CLEANUP_TAGGED_CO deallocates memory from previously
!  allocated module arrays (bmy, 7/19/00, 7/20/04)
!
!  NOTES:
!  (1 ) Added ISOP96, MONO96, CH3OH96 for GEOS-3 (bnd, bmy, 6/14/01)
!  (2 ) Removed ISOP96, MONO96, CH3OH96 for GEOS-3, since the new GEOS-3
!        fields make these no longer necessary (bmy, 8/21/09)
!  (3 ) Now remove IJLOOP_CO (bmy, 7/20/04)
!******************************************************************************
!
      IF ( ALLOCATED( BB_REGION  ) ) DEALLOCATE( BB_REGION  )
      IF ( ALLOCATED( FF_REGION  ) ) DEALLOCATE( FF_REGION  )
      IF ( ALLOCATED( SUMISOPCO  ) ) DEALLOCATE( SUMISOPCO  )
      IF ( ALLOCATED( SUMMONOCO  ) ) DEALLOCATE( SUMMONOCO  )
      IF ( ALLOCATED( SUMCH3OHCO ) ) DEALLOCATE( SUMCH3OHCO )
      IF ( ALLOCATED( SUMACETCO  ) ) DEALLOCATE( SUMACETCO  )
      IF ( ALLOCATED( EMACET     ) ) DEALLOCATE( EMACET     )
      IF ( ALLOCATED( CO_PRODS   ) ) DEALLOCATE( CO_PRODS   )
      IF ( ALLOCATED( CO_LOSSS   ) ) DEALLOCATE( CO_LOSSS   )

      ! Obsolete (zhej, dkh, 01/17/12, adj32_017) 
      !IF ( ALLOCATED( COSF       ) ) DEALLOCATE( COSF       )
      !IF ( ALLOCATED( CO_CH4     ) ) DEALLOCATE( CO_CH4     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TAGGED_CO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE TAGGED_CO_MOD
