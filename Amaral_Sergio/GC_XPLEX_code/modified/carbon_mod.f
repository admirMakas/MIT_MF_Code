! $Id: carbon_mod.f,v 1.4 2011/02/23 00:08:48 daven Exp $
      MODULE CARBON_MOD
!
!******************************************************************************
!  Module CARBON_MOD contains arrays and routines for performing a 
!  carbonaceous aerosol simulation.  Original code taken from Mian Chin's 
!  GOCART model and modified accordingly. (rjp, bmy, 4/2/04, 2/19/09)
!
!  4 Aerosol species : Organic and Black carbon 
!                    : hydrophilic (soluble) and hydrophobic of each
!  For secondary organic aerosol (SOA) simulation orginal code developed
!  by Chung and Seinfeld [2002] and Hong Liao from John Seinfeld's group 
!  at Caltech was taken and further modified accordingly (rjp, bmy, 7/15/04)
!  This simulation introduces additional following species:
!     ALPH, LIMO, ALCO, SOA1, SOA2, SOA3, SOA4, SOG1, SOG2, SOG3, SOG4
!
!  References for BC/OC emissions:
! --------------------------------
!   Bond, T.C. et al.: Historical emissions of black and organic carbon
! aerosol fromenergy-related combustion, 1850-2000, Global Biogeochem. Cycles,
! 21 GB2018, doi: 10.1029/2006GB002840, 2007.
!
!      Base Year is 2000. More at http://www.hiwater.org
!
!  Module Variables:
!  ============================================================================
!  (1 ) ANTH_BLKC        (TYPE (XPLEX) ) : BC anthropogenic emissions       [kg C ]
!  (2 ) ANTH_ORGC        (TYPE (XPLEX) ) : OC anthropogenic emissions       [kg C ]
!  (3 ) APROD            (TYPE (XPLEX) ) : Aerosol mass ratio               
!  (4 ) BCCONV           (TYPE (XPLEX) ) : Hydrophilic BC from Hydrophobic  [kg C ]
!  (5 ) BIOB_BLKC        (TYPE (XPLEX) ) : BC biomass emissions             [kg C ]
!  (6 ) BIOB_ORGC        (TYPE (XPLEX) ) : OC biomass emissions             [kg C ]
!  (7 ) BIOF_BLKC        (TYPE (XPLEX) ) : BC biofuel emissions             [kg C ]
!  (8 ) BIOF_ORGC        (TYPE (XPLEX) ) : OC biofuel emissions             [kg C ]
!  (9 ) BIOG_ALPH        (TYPE (XPLEX) ) : A-PINENE biogenic emissions      [kg]
!  (10) BIOG_LIMO        (TYPE (XPLEX) ) : LIMONENE biogenic emissions      [kg]
!  (11) BIOG_ALCO        (TYPE (XPLEX) ) : ALCOHOL biogenic emissions       [kg]
!  (12) BIOG_TERP        (TYPE (XPLEX) ) : TERPENE biogenic emissions       [kg]
!  (13) BIOG_SESQ        (TYPE (XPLEX) ) : SESQTERPENE biogenic emissions   [kg]
!  (14) DIUR_ORVC        (TYPE (XPLEX) ) : Diurnal variation of NVOC        [kg C]
!  (15) DRYBCPI          (INTEGER) : Index for BCPI in drydep array 
!  (16) DRYOCPI          (INTEGER) : Index for OCPI in drydep array 
!  (17) DRYBCPO          (INTEGER) : Index for BCPO in drydep array 
!  (18) DRYOCPO          (INTEGER) : Index for OCPO in drydep array 
!  (19) DRYALPH          (INTEGER) : Index for ALPH in drydep array 
!  (20) DRYLIMO          (INTEGER) : Index for LIMO in drydep array 
!  (21) DRYALCO          (INTEGER) : Index for ALCO in drydep array 
!  (22) DRYSOG1          (INTEGER) : Index for SOG1 in drydep array 
!  (23) DRYSOG2          (INTEGER) : Index for SOG2 in drydep array 
!  (24) DRYSOG3          (INTEGER) : Index for SOG3 in drydep array 
!  (25) DRYSOA1          (INTEGER) : Index for SOA1 in drydep array 
!  (26) DRYSOA2          (INTEGER) : Index for SOA2 in drydep array 
!  (27) DRYSOA3          (INTEGER) : Index for SOA3 in drydep array 
!  (28) EF_BLKC          (TYPE (XPLEX) ) : Emission factors for BC          [kg/kg]
!  (29) EF_ORGC          (TYPE (XPLEX) ) : Emission factors for OC          [kg/kg]
!  (30) GEIA_ORVC        (TYPE (XPLEX) ) : NVOC emissions from GEIA         [kg C ]
!  (31) I1_NA            (TYPE (XPLEX) ) : Starting lon index for N. America region
!  (32) I2_NA            (TYPE (XPLEX) ) : Ending   lon index for N. America region
!  (33) J1_NA            (TYPE (XPLEX) ) : Starting lat index for N. America region
!  (34) J2_NA            (TYPE (XPLEX) ) : Ending   lat index for N. America region
!  (35) GPROD            (TYPE (XPLEX) ) : Gas mass ratio                   
!  (36) MHC              (INTEGER) : Number of VOC clasees
!  (37) NDAYS            (INTEGER) : Array w/ number of days per month
!  (38) NPROD            (INTEGER) : Number of products by oxdation
!  (39) OCCONV           (TYPE (XPLEX) ) : Hydrophilic OC from Hydrophobic  [kg C ]
!  (40) ORVC_SESQ        (TYPE (XPLEX) ) : SESQTERPENE concentration        [kg]
!  (41) ORVC_TERP        (TYPE (XPLEX) ) : MONOTERPENES concentration       [kg]
!  (42) TCOSZ            (TYPE (XPLEX) ) : Summing array for SUNCOS
!  (43) TERP_ORGC        (TYPE (XPLEX) ) : Lumped terpene emissions         [kg C ]
!  (44) SMALLNUM         (TYPE (XPLEX) ) : A small positive number
!  (45) USE_BOND_BIOBURN (LOGICAL) : Flag to use annual biomass emiss.
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMCARBON         : Driver program for carbon aerosol chemistry
!  (2 ) CHEM_BCPO          : Chemistry routine for hydrophobic BC (aka EC)
!  (3 ) CHEM_BCPI          : Chemistry routine for hydrophilic BC (aka EC)
!  (4 ) CHEM_OCPO          : Chemistry routine for hydrophobic OC
!  (5 ) CHEM_OCPI          : Chemistry routine for hydrophilic OC
!  (6 ) SOA_CHEMISTRY      : Driver routine for SOA chemistry
!  (7 ) SOA_EQUIL          : Determines mass of SOA
!  (8 ) ZEROIN             : Finds root of an equation by bisection method
!  (9 ) SOA_PARA           : Gets SOA yield parameters
!  (10) CHEM_NVOC          : Computes oxidation of HC by O3, OH, NO3
!  (11) SOA_PARTITION      : Partitions mass of SOA gas & aerosol tracers
!  (12) SOA_LUMP           : Returns organic gas & aerosol back to STT
!  (13) SOA_DEPO           : Performs dry deposition of SOA tracers & species
!  (14) EMISSCARBON        : Driver routine for carbon aerosol emissions
!  (15) BIOGENIC_OC        : Computes biogenic OC [each time step]
!  (16) ANTHRO_CARB_TBOND  : Computes anthropogenic OC/EC [annual data]
!  (17) ANTHRO_CARB_COOKE  : Computes anthropogenic OC/EC [monthly data]
!  (18) BIOMASS_CARB_TBOND : Computes biomass burning OC/EC [annual data]
!  (19) BIOMASS_CARB_GEOS  : Computes biomass burning OC/EC [monthly data]
!  (20) EMITHIGH           : Computes complete mixing of emission within PBL
!  (21) OHNO3TIME          : Computes the sum of the cosine of SZA
!  (22) GET_OH             : Returns monthly-mean OH conc.  at grid box (I,J,L)
!  (23) GET_NO3            : Returns monthly-mean O3 conc.  at grid box (I,J,L)
!  (24) GET_O3             : Returns monthly-mean NO3 conc. at grid box (I,J,L)
!  (25) GET_DOH            : Returns ISOP lost to rxn w/ OH at grid box (I,J,L)
!  (26) INIT_CARBON        : Allocates and zeroes all module arrays
!  (27) CLEANUP_CARBON     : Deallocates all module arrays
!
!  NOTE: Choose either (16) or (17) for ANTHROPOGENIC emission
!        Choose either (18) or (19) for BIOMASS BURNING emission.
!
!  GEOS-CHEM modules referenced by carbon_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f    : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) drydep_mod.f       : Module w/ routines for dry deposition
!  (6 ) error_mod.f        : Module w/ I/O error and NaN check routines
!  (7 ) global_no3_mod.f   : Module w/ routines to read 3-D NO3 field
!  (8 ) global_oh_mod.f    : Module w/ routines to read 3-D OH  field
!  (9 ) global_o3_mod.f    : Module w/ routines to read 3-D O3  field
!  (10) grid_mod.f         : Module w/ horizontal grid information
!  (11) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (12) megan_mod.f        : Module w/ routines to read MEGAN biogenic emiss
!  (13) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (14) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (15) time_mod.f         : Module w/ routines for computing time & date
!  (16) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc. 
!  (17) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!  (18) transfer_mod.f     : Module w/ routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Added code from the Caltech group for SOA chemistry (rjp, bmy, 7/15/04)
!  (2 ) Now references "directory_mod.f", "logical_mod.f", "tracer_mod.f".
!        (bmy, 7/20/04)
!  (3 ) Now read data from carbon_200411/ subdir of DATA_DIR.  Also added
!        some extra debug output.  Now read T. Bond yearly emissions as 
!        default, but overwrite N. America with the monthly Cooke/RJP 
!        emissions.  Added module variables I1_NA, I2_NA, J1_NA, J2_NA.
!        (rjp, bmy, 12/1/04)
!  (4 ) Now can read seasonal or interannual BCPO, OCPO biomass emissions.
!        Also parallelize loop in OHNO3TIME. (rjp, bmy, 1/18/05)
!  (5 ) Now references "pbl_mix_mod.f".  Bug fix: now make sure only to save
!        up to LD07 levels for the ND07 diagnostic in SOA_LUMP. (bmy, 3/4/05)
!  (6 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Now references "megan_mod.f".  Also now references XNUMOL and 
!        XNUMOLAIR from "tracer_mod.f" (tmf, bmy, 10/25/05)
!  (9 ) Bug fix for GCAP in BIOGENIC_OC (bmy, 4/11/06)
!  (10) Updated for SOA production from ISOP (dkh, bmy, 5/22/06)
!  (11) Updated for IPCC future emission scale factors.  Also added function
!        GET_DOH to return ISOP that has reacted w/ OH. (swu, dkh, bmy, 6/1/06)
!  (12) Now add SOG condensation onto SO4, NH4, NIT (rjp, bmy, 8/3/06)
!  (13) Minor fix for 20 carbon tracers. (phs, 9/14/06)
!  (14) Now remove reading of biomass emissions from "carbon_mod.f", since
!        they are better done in gc_biomass_mod.f.  This will allow us to
!        standardize treatment of GFED2 or default BB emissions.  Also applied
!        a typo fix in SOA_LUMP. (tmf, bmy, 10/16/06)
!  (15) Prevent seg fault error in BIOMASS_CARB_GEOS (bmy, 11/3/06)
!  (16) Corrected typos in SOA_LUMP.  Now also save GPROD and APROD to disk
!        for each new diagnostic interval. (dkh, tmv, havala, bmy, 2/6/07)
!  (17) Modifications for 0.5 x 0.666 nested grids (yxw, dan, bmy, 11/6/08)
!  (18) Now account for various GFED2 products (yc, phs, 12/23/08) 
!  (19) Now add future scaling to BIOMASS_CARB_GEOS (hotp, swu, 2/19/09)
!  (20) Added SOA production from dicarbonyls (tmf, 3/2/09)
!  (21) Bugfix: cleanup ORVC_TERP and ORVC_SESQ (tmf, 3/2/09)
!  (22) Replace USE_MONTHLY_BIOB with USE_BOND_BIOBURN, since this hardwired
!        flag is a switc b/w annual Bond biomass burning emissions, and default
!        GC source, which can be monthly/8 days/3hr.
!        Implement changes for reading new Bond files (eml, phs, 5/18/09)
!  (23) Add option for non-local PBL scheme (lin, 06/09/08)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "carbon_mod.f"
      !=================================================================

      ! Declare everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMCARBON
      PUBLIC :: CLEANUP_CARBON
      PUBLIC :: EMISSCARBON
      PUBLIC :: WRITE_GPROD_APROD

      ! adj_group: make these public (dkh, 09/09/09) 
      PUBLIC :: DRYBCPI,   DRYOCPI,   DRYBCPO,   DRYOCPO
      ! ... and these (dkh, 11/11/09) 
      PUBLIC :: ANTH_ORGC, ANTH_BLKC, BIOB_ORGC, BIOB_BLKC
      PUBLIC :: BIOF_ORGC, BIOF_BLKC

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      !LOGICAL             :: USE_MONTHLY_BIOB = .TRUE.
      LOGICAL             :: USE_BOND_BIOBURN = .FALSE.
      INTEGER             :: DRYBCPI, DRYOCPI, DRYBCPO, DRYOCPO
      INTEGER             :: DRYALPH, DRYLIMO, DRYALCO
      INTEGER             :: DRYSOG1, DRYSOG2, DRYSOG3, DRYSOG4
      INTEGER             :: DRYSOA1, DRYSOA2, DRYSOA3, DRYSOA4
      INTEGER             :: I1_NA,   J1_NA
      INTEGER             :: I2_NA,   J2_NA
      INTEGER             :: DRYSOAG, DRYSOAM

      ! Parameters
      INTEGER, PARAMETER  :: MHC      = 6
      INTEGER, PARAMETER  :: NPROD    = 3  
      TYPE (XPLEX),  PARAMETER  :: SMALLNUM = xplex(1d-20,0d0)

      ! Arrays
      TYPE (XPLEX), ALLOCATABLE :: ANTH_BLKC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: ANTH_ORGC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOB_BLKC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOB_ORGC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOF_BLKC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOF_ORGC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: EF_BLKC(:,:)
      TYPE (XPLEX), ALLOCATABLE :: EF_ORGC(:,:)
      TYPE (XPLEX), ALLOCATABLE :: TERP_ORGC(:,:)
      TYPE (XPLEX), ALLOCATABLE :: BCCONV(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: OCCONV(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOG_ALPH(:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOG_LIMO(:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOG_ALCO(:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOG_TERP(:,:)
      TYPE (XPLEX), ALLOCATABLE :: BIOG_SESQ(:,:)    
      TYPE (XPLEX), ALLOCATABLE :: DIUR_ORVC(:,:)
      TYPE (XPLEX), ALLOCATABLE :: GEIA_ORVC(:,:)
      TYPE (XPLEX), ALLOCATABLE :: TCOSZ(:,:)
      TYPE (XPLEX), ALLOCATABLE :: ORVC_SESQ(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: ORVC_TERP(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: GPROD(:,:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: APROD(:,:,:,:,:)
      ! Cloud fraction - for cloud droplet uptake of dicarbonyls 
      ! (tmf, 12/07/07) 
      TYPE (XPLEX), ALLOCATABLE :: VCLDF(:,:,:)

      ! Days per month (based on 1998)
      INTEGER             :: NDAYS(12) = (/ 31, 28, 31, 30, 31, 30, 
     &                                      31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMCARBON
!
!******************************************************************************
!  Subroutine CHEMCARBON is the interface between the GEOS-CHEM main program 
!  and the carbon aerosol chemistry routines that calculates dry deposition, 
!  chemical conversion between hydrophilic and hydrophobic, and SOA production.
!  (rjp, bmy, 4/1/04, 9/14/06)
!
!  NOTES:
!  (1 ) Added code from the Caltech group for SOA chemistry.  Also now 
!        reference "global_oh_mod.f", "global_o3_mod.f", "global_no3_mod.f".
!        (rjp, bmy, 7/8/04)
!  (2 ) Now reference LSOA and LEMIS from CMN_SETUP.  Now only call OHNO3TIME
!        if it hasn't been done before w/in EMISSCARBON. (rjp, bmy, 7/15/04)
!  (3 ) Now reference LSOA, LEMIS, LPRT from "logical_mod.f".  Now reference
!        STT and ITS_AN_AEROSOL_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now updated for SOA production from ISOP. (dkh, bmy, 6/1/06)
!  (6 ) Bug fix for aerosol sim w/ 20 tracers (phs, 9/14/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,     ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH
      USE GLOBAL_NO3_MOD, ONLY : GET_GLOBAL_NO3
      USE GLOBAL_O3_MOD,  ONLY : GET_GLOBAL_O3
      USE LOGICAL_MOD,    ONLY : LSOA, LEMIS, LPRT
      USE TIME_MOD,       ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : STT, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,   ONLY : IDTBCPI, IDTBCPO, IDTOCPI
      USE TRACERID_MOD,   ONLY : IDTOCPO, IDTSOG4, IDTSOA4
      USE TRACERID_MOD,   ONLY : IDTSOAG, IDTSOAM

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      LOGICAL, SAVE           :: FIRSTCHEM = .TRUE.
      INTEGER                 :: N, THISMONTH

      !=================================================================
      ! CHEMCARBON begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_CARBON

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'BCPI' )
                  DRYBCPI = N
               CASE ( 'OCPI' )
                  DRYOCPI = N
               CASE ( 'BCPO' )
                  DRYBCPO = N
               CASE ( 'OCPO' )
                  DRYOCPO = N
               CASE ( 'ALPH' )
                  DRYALPH = N
               CASE ( 'LIMO' )
                  DRYLIMO = N
               CASE ( 'ALCO' )
                  DRYALCO = N
               CASE ( 'SOG1' )
                  DRYSOG1 = N
               CASE ( 'SOG2' )
                  DRYSOG2 = N
               CASE ( 'SOG3' )
                  DRYSOG3 = N
               CASE ( 'SOG4' )
                  DRYSOG4 = N
               CASE ( 'SOA1' )
                  DRYSOA1 = N
               CASE ( 'SOA2' )
                  DRYSOA2 = N
               CASE ( 'SOA3' )
                  DRYSOA3 = N
               CASE ( 'SOA4' )
                  DRYSOA4 = N
               CASE ( 'SOAG' )
                  DRYSOAG = N
               CASE ( 'SOAM' )
                  DRYSOAM = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! Zero SOG4 and SOA4 (SOA from ISOP in gas & aerosol form)
         ! for offline aerosol simulations.  Eventually we should have
         ! archived isoprene oxidation fields available for offline
         ! simulations but for now we just set them to zero. 
         ! (dkh, bmy, 6/1/06)
         IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! temp fix for aerosol w/ 20 tracers simulation (phs)
            IF ( IDTSOG4 .NE. 0 ) THEN   
               STT(:,:,:,IDTSOG4) = 0d0
               STT(:,:,:,IDTSOA4) = 0d0
            ENDIF
         ENDIF
         
         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      !=================================================================
      ! Do chemistry for carbon aerosol tracers 
      !=================================================================

      ! Chemistry for hydrophobic BC
      IF ( IDTBCPO > 0 ) THEN
         CALL CHEM_BCPO( STT(:,:,:,IDTBCPO) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPO' )
      ENDIF

      ! Chemistry for hydrophilic BC
      IF ( IDTBCPI > 0 ) THEN
         CALL CHEM_BCPI( STT(:,:,:,IDTBCPI) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPI' )
      ENDIF

      ! Chemistry for hydrophobic OC
      IF ( IDTOCPO > 0 ) THEN
         CALL CHEM_OCPO( STT(:,:,:,IDTOCPO) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPO' )
      ENDIF

      ! Chemistry for hydrophilic OC
      IF ( IDTOCPI > 0 ) THEN 
         CALL CHEM_OCPI( STT(:,:,:,IDTOCPI) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPI' )
      ENDIF

      !=================================================================
      ! Do chemistry for secondary organic aerosols 
      !=================================================================
      IF ( LSOA ) THEN

         ! Read offline OH, NO3, O3 fields from disk
         IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! Current month
            THISMONTH = GET_MONTH()

            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_OH(  THISMONTH )
               CALL GET_GLOBAL_NO3( THISMONTH )
               CALL GET_GLOBAL_O3(  THISMONTH )
            ENDIF

            ! Compute time scaling arrays for offline OH, NO3
            ! but only if it hasn't been done in EMISSCARBON
            IF ( LSOA .and. ( .not. LEMIS ) ) THEN
               CALL OHNO3TIME
               IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARB: a OHNO3TIME' )
            ENDIF
         ENDIF

         ! Compute SOA chemistry
         ! NOTE: This is SOA production from the reversible mechanism only 
         ! (tmf, 12/07/07) 
         CALL SOA_CHEMISTRY
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a SOA_CHEM' )

         ! If SOAG and SOAM are declared, switch on 
         !    SOA production from dicarbonyls (tmf, 12/07/07) 
         IF ( IDTSOAG > 0 ) THEN

            ! Get grid box cloud fraction 
            ! (tmf, 2/26/07)
            CALL GET_VCLDF

            ! Cloud uptake
            CALL SOAG_CLOUD
            IF ( LPRT ) 
     &       CALL DEBUG_MSG('### CHEMCARBON: a SOAG_CLOUD')        

            ! Aqueous aerosol uptake
            CALL SOAG_LIGGIO_DIFF
            IF ( LPRT ) 
     &       CALL DEBUG_MSG('### CHEMCARBON: a SOAG_LIGGIO_DIFF')        

         ENDIF

         IF ( IDTSOAM > 0 ) THEN
         
            ! Get grid box cloud fraction 
            ! (tmf, 2/26/07)
            CALL GET_VCLDF

            ! Cloud uptake
            CALL SOAM_CLOUD
            IF ( LPRT ) 
     &       CALL DEBUG_MSG('### CHEMCARBON: a SOAM_CLOUD')        

            ! Aqueous aerosol uptake
            CALL SOAM_LIGGIO_DIFF
            IF ( LPRT ) 
     &       CALL DEBUG_MSG( '### CHEMCARBON: a SOAM_LIGGIO_DIFF' )        


           
         ENDIF   


      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMCARBON

!-----------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPO( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPO converts hydrophobic BC to hydrophilic BC and
!  calculates the dry deposition of hydrophobic BC. (rjp, bmy, 4/1/04,10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophobic BC tracer 
!
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!  (2 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP 
!        from "pbl_mix_mod.f" (bmy, 2/17/05)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_BC 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTBCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,       J,   L
      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)     
      TYPE (XPLEX)                 :: DTCHEM, FLUX, KBC, FREQ
      TYPE (XPLEX)              :: TC0,    CNEW, RKT, AREA_CM2, BL_FRAC
      TYPE (XPLEX),  PARAMETER     :: BC_LIFE = xplex(1.15D0,0d0)

      !=================================================================
      ! CHEM_BCPO begins here!
      !=================================================================

      ! Return if BCPO isn't defined
      IF ( IDTBCPO == 0 .or. DRYBCPO == 0 ) RETURN

      ! Initialize
      KBC    = 1.D0 / ( 86400d0 * BC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BCCONV(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !
      ! Hydrophobic(2) --> Hydrophilic(1) ,  k  = 1.0e-5          
      ! Both aerosols are dry-deposited,     kd = Dvel/DELZ (sec-1)      
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial BC mass [kg]
         TC0  = TC(I,J,L)

         ! Zero drydep freq
         FREQ = 0d0

         ! Fraction of grid box underneath PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! BC drydep frequency [1/s] -- BL_FRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYBCPO) * BL_FRAC

         ENDIF

         ! Amount of BCPO left after chemistry and drydep [kg]
         RKT  = ( KBC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Amount of BCPO converted to BCPI [kg/timestep]
         BCCONV(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )

         !==============================================================
         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

             ! Surface area [cm2]
             AREA_CM2 = GET_AREA_CM2( J )

             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]  
             ! XNUMOL is the ratio [molec tracer/kg tracer]   
             FLUX     = TC0 - CNEW - BCCONV(I,J,L) 
             FLUX     = FLUX * XNUMOL(IDTBCPO) / ( DTCHEM * AREA_CM2 )

             ! Store in ND44_TMP as a placeholder
             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF

         !==============================================================
         ! ND07 diagnostic: H-philic BC from H_phobic BC [kg/timestep]
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
             AD07_BC(I,J,L) = AD07_BC(I,J,L) + BCCONV(I,J,L)
         ENDIF

         ! Store new concentration back into tracer array
         TC(I,J,L) = CNEW
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYBCPO,1) = AD44(I,J,DRYBCPO,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPO

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPI( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPI calculates dry deposition of hydrophilic BC.
!  (rjp, bmy, 4/1/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!  (2 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f" (bmy, 2/17/05)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTBCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, L
      TYPE (XPLEX)                 :: DTCHEM, FLUX, BL_FRAC, AREA_CM2
      TYPE (XPLEX)                 :: TC0,    CNEW, CCV,     FREQ
      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_BCPI begins here!
      !=================================================================

      ! Return if BCPI isn't defined
      IF ( IDTBCPI == 0 .or. DRYBCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for ND44 diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, FREQ, BL_FRAC, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic BC [kg]
         TC0 = TC(I,J,L)

         ! H-philic BC that used to be H-phobic BC [kg]
         CCV = BCCONV(I,J,L)
         
         ! Fraction of grid box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L ) 

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency
            FREQ = DEPSAV(I,J,DRYBCPI) * BL_FRAC
            
            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT)) 
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! Comment out for now
            !CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !     + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep flux [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 .and. FREQ > 0d0 ) THEN
  
               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [molec/cm2/s]
               FLUX = ( TC0 + CCV - CNEW ) 
               FLUX = FLUX * XNUMOL(IDTBCPI) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, omit the exponential to save on clock cycles
            CNEW = TC0 + CCV

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Save new concentration of H-philic IC in tracer array
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Zero out the BCCONV array for the next iteration
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BCCONV(I,J,L) = 0.d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYBCPI,1) = AD44(I,J,DRYBCPI,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPI

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_OCPO( TC )
!
!******************************************************************************
!  Subroutine CHEM_OCPO converts hydrophobic OC to hydrophilic OC and
!  calculates the dry deposition of hydrophobic OC. (rjp, bmy, 4/1/04,10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophobic OC tracer [kg]
!
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!  (2 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f" (bmy, 2/17/05)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_OC 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTOCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                 :: DTCHEM, FLUX, KOC,  BL_FRAC
      TYPE (XPLEX)                 :: TC0,    FREQ, CNEW, RKT, AREA_CM2
      TYPE (XPLEX),  PARAMETER     :: OC_LIFE = xplex(1.15D0,0d0)

      !=================================================================
      ! CHEM_OCPO begins here!
      !=================================================================

      ! Return if OCPO isn't defined
      IF ( IDTOCPO == 0 .or. DRYOCPO == 0 ) RETURN

      ! Initialize
      KOC    = 1.D0 / ( 86400d0 * OC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCCONV(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !    Hydrophobic --> Hydrophilic,  k  = 1.0e-5          
      !    Aerosols are dry-deposited,   kd = DEPSAV (sec-1)      
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial OC [kg]
         TC0  = TC(I,J,L)

         ! Zero drydep freq 
         FREQ = 0d0

         ! Fraction of grid box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! OC drydep frequency [1/s] -- BL_FRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYOCPO) * BL_FRAC

         ENDIF

         ! Amount of OCPO left after chemistry and drydep [kg]
         RKT  = ( KOC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Amount of OCPO converted to OCPI [kg/timestep]
         OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

         !==============================================================
         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

             ! Surface area [cm2]
             AREA_CM2 = GET_AREA_CM2( J )

             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
             ! XNUMOL is the ratio [molec tracer/kg tracer]     
             FLUX     = TC0 - CNEW - OCCONV(I,J,L)
             FLUX     = FLUX * XNUMOL(IDTOCPO) / ( DTCHEM * AREA_CM2 )

             ! Store in ND44_TMP as a placeholder
             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF

         !==============================================================
         ! ND07 diagnostic: H-Philic OC from H-phobic [kg/timestep]
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
            AD07_OC(I,J,L) = AD07_OC(I,J,L) + OCCONV(I,J,L) 
         ENDIF

         ! Store modified OC concentration back in tracer array
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYOCPO,1) = AD44(I,J,DRYOCPO,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF   

      ! Return to calling program
      END SUBROUTINE CHEM_OCPO

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_OCPI( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPI calculates dry deposition of hydrophilic OC.
!  (rjp, bmy, 4/1/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!  (2 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f" (bmy, 2/17/05)
!  (3 ) Bug fix: add BL_FRAC to the PRIVATE list (mak, bmy, 10/3/05)
!  (4 ) Now refrerences XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTOCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      TYPE (XPLEX)                 :: DTCHEM, FLUX, BL_FRAC
      TYPE (XPLEX)                 :: TC0, CNEW, CCV, FREQ, AREA_CM2
      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_OCPI begins here!
      !=================================================================
      IF ( IDTOCPI == 0 .or. DRYOCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for drydep diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, BL_FRAC, FREQ, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic OC [kg]
         TC0 = TC(I,J,L)

         ! H-philic OC that used to be H-phobic OC [kg]
         CCV = OCCONV(I,J,L)

         ! Fraction of box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DRYOCPI) * BL_FRAC

            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT))
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !       + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
               FLUX = ( TC0 + CCV - CNEW ) 
               FLUX = FLUX * XNUMOL(IDTOCPI) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            CNEW = TC0 + CCV

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Store modified concentration back in tracer array [kg]
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Zero OCCONV array for next timestep
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCCONV(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYOCPI,1) = AD44(I,J,DRYOCPI,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF    

      ! Return to calling program
      END SUBROUTINE CHEM_OCPI

!-------------------------------------------------------------------------------

      SUBROUTINE SOAG_LIGGIO_DIFF
!
!******************************************************************************
!  Subroutine SOAG_LIGGIO_DIFF produces SOA on aqueous aerosol surfaces
!   from GLYX following the uptake model used for N2O5, and the gamma 
!   from Liggio et al. [2005]. (tmf, 5/30/06)
!
!  Procedure:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  (1 ) SOAG (SOA product of GLYX is produced at existing hydrophilic aerosol
!        surface.
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : WTAREA, WERADIUS 
      USE COMODE_MOD,   ONLY : AIRDENS, JLOP
      USE DAO_MOD,      ONLY : AIRVOL, T, RH

      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE DIAG_MOD,     ONLY : AD07_SOAGM
      USE TIME_MOD,     ONLY : GET_TS_CHEM, GET_MONTH
      USE TRACER_MOD,   ONLY : STT 
      USE TRACERID_MOD, ONLY : IDTGLYX, IDTSOAG
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44, ND07, LD07
#     include "comode.h"     ! AD, WTAIR, other SMVGEAR variables


      ! Local variables
      INTEGER              :: I, J, L, JLOOP, N
      TYPE (XPLEX)               :: XTEMP       ! Temperature [K]
      TYPE (XPLEX)               :: XSQTEMP     ! SQRT of Temperature
      TYPE (XPLEX)               :: XAD         ! Air density [molec/cm3]
      TYPE (XPLEX)               :: XRADIUS     ! particle radius [cm]
      TYPE (XPLEX)               :: XDFKG       ! Gas phase diffusion coeff [cm2/s]
      TYPE (XPLEX)               :: XARSL1K     ! 1st order k
      TYPE (XPLEX)               :: XRH         ! Relative humidity [%]
      TYPE (XPLEX)               :: XAIRM3      ! Air volume in grid box [m3]
      TYPE (XPLEX)               :: XGASM       ! Gas mass at grid box before uptake [kg]
      TYPE (XPLEX)               :: XGASC       ! Gas concentration at grid box before uptake [molec/cm3]
      TYPE (XPLEX)               :: XWAREA      ! Wet aerosol surface area at grid box [cm^2 wet sfc area of aerosol cm^-3 air]
      TYPE (XPLEX)               :: XUPTK0      ! Potential uptake of gas by aerosol in grid box by aerosol type [molec/cm3]
      TYPE (XPLEX)               :: XUPTK1      ! Potential uptake of gas by aerosol in grid box by aerosol type [kg]
      TYPE (XPLEX)               :: XUPTKSUM    ! Potential uptake of gas by aerosol in grid box [kg]
      TYPE (XPLEX)               :: XUPTK       ! Actual uptake of gas by aerosol in grid box [kg]
                                          !  XUPTK <= STT( I, J, L, IDTGLYX )
      TYPE (XPLEX)               :: XGAMMA      ! Uptake coefficient 

      ! Local variables not changing 
      TYPE (XPLEX)               :: DTCHEM      ! Chemistry time step [s]
      TYPE (XPLEX)               :: XMW         ! Molecular weight of gas [g/mole]
      TYPE (XPLEX)               :: XSQMW       ! Square root of molecular weight [g/mole]
      TYPE (XPLEX)               :: CRITRH      ! Critical RH [%], above which 
                                          !  heteorogeneous chem takes place
      TYPE (XPLEX)               :: XNAVO       ! Avogadro number

      !=================================================================
      ! SOAG_LIGGIO_DIFF begins here!
      !=================================================================

      ! Get chemistry time step
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Molecular weight of GLYX [g/mole]
      XMW   = 58.d0
      XSQMW = SQRT( XMW )

      ! Critical RH, above which heteorogeneous chem takes place
      CRITRH = 35.0d0   ! [%]

      ! Avogadro number
      XNAVO = 6.022d23

      ! Uptake coefficient from Liggio et al. [2005b]
      XGAMMA = 2.9d-3

      !=================================================================
      ! Loop over grid boxes
      !=================================================================
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

            ! Get 1-D index
            JLOOP   = JLOP( I, J, L )

            ! Get RH  
            XRH     = RH( I, J, L )   ![%]

            ! initialize for safety
            XUPTK0   = 0d0
            XUPTK1   = 0d0
            XUPTKSUM = 0d0
            XUPTK    = 0d0

            ! Get T
            XTEMP   = T( I, J, L )
            XSQTEMP = SQRT( XTEMP )

            ! Get air density  [molec/cm3]
            XAD     = AIRDENS( JLOOP )

            ! Get air volumne [m3]
            XAIRM3  = AIRVOL( I, J, L )

            ! Get gas mass at grid box [kg]
            XGASM   = STT( I, J, L, IDTGLYX )

            ! Get gas concentration at grid box [molec/cm3]
            XGASC   = XGASM / (XMW*1.d-3) * XNAVO / (XAIRM3*1.d6) 

            !---------------------------------------
            ! Gas phase diffusion coeff [cm2/s]           
            !---------------------------------------
            XDFKG = 9.45D17 / XAD * XSQTEMP * 
     &              SQRT( 3.472D-2 + (1.D0/XMW) )


            !========================================================
            ! Calculate heteorogeneous uptake only if the grid box
            !  relative humidity XRH is >= critical relative humidity CRITRH
            !========================================================
            IF ( XRH >= CRITRH ) THEN

               ! Loop over sulfate and other aerosols
               DO N = 1, NDUST + NAER

                  !---------------------------------------
                  ! Total available wet aerosol area 
                  !  archived in 'aerosol_mod.f.glyx'
                  !  XWAREA [ cm^2 wet sfc area of aerosol cm^-3 air ]
                  !---------------------------------------
                  XWAREA  = WTAREA( JLOOP, N) 

                  IF ( XWAREA > 0D0 ) THEN 

                     ! Get particle radius [cm]
                     XRADIUS = WERADIUS( JLOOP, N )   

                     !---------------------------------------
                     ! First order rate constant
                     !---------------------------------------
                     XARSL1K = XWAREA / 
     &               (XRADIUS/XDFKG + 2.749064D-4*XSQMW/XGAMMA/XSQTEMP)

                     !---------------------------------------
                     ! Calculate potential uptake: Liggio et al. (2005b) Eq (3)
                     !   
                     !   d( organic carbon conc ) / dt = 
                     !      XARSL1K * XGASC 
                     !---------------------------------------
                     XUPTK0 = XARSL1K * XGASC * DTCHEM
                     XUPTK1 = XUPTK0 / XNAVO*(XMW*1.d-3)*(XAIRM3*1.d6)
                     XUPTKSUM = XUPTKSUM + XUPTK1

	            ENDIF
           
               ENDDO

                  ! However, the mass of gas being absorbed by aerosol 
                  !  cannot exceed the original amount of gas XGASM
                  XUPTK  = MIN( XUPTKSUM, XGASM )
            
                  ! Update GLYX in the STT array
                  STT( I, J, L, IDTGLYX ) = STT( I, J, L, IDTGLYX ) -
     &                                      XUPTK

                  ! Update SOAG in the STT array
                  STT( I, J, L, IDTSOAG ) = STT( I, J, L, IDTSOAG ) + 
     &                                      XUPTK

            ENDIF

         !==============================================================
         ! ND07 diagnostic: SOAG from GLYX [kg/timestep] on aerosol
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
            AD07_SOAGM(I,J,L,1) = AD07_SOAGM(I,J,L,1) + XUPTK
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      !=================================================================       
      ! Calculate dry-deposition
      !=================================================================
      CALL SOA_DEPO( STT(:,:,:,IDTSOAG), DRYSOAG, IDTSOAG )

      ! Return to calling program 
      END SUBROUTINE SOAG_LIGGIO_DIFF

!------------------------------------------------------------------------------

      SUBROUTINE SOAM_LIGGIO_DIFF
!
!******************************************************************************
!  Subroutine SOAG_LIGGIO_DIFF produces SOA on aqueous aerosol surfaces
!   from GLYX following the uptake model used for N2O5, and the gamma 
!   from Liggio et al. [2005]. (tmf, 5/30/06)
!
!  Procedure:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  (1 ) SOAM (SOA product of MGLY) is produced at existing hydrophilic aerosol
!        surface.
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : WTAREA, WERADIUS
      USE COMODE_MOD,   ONLY : AIRDENS, JLOP
      USE DAO_MOD,      ONLY : AIRVOL, T, RH

      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE DIAG_MOD,     ONLY : AD07_SOAGM
      USE TIME_MOD,     ONLY : GET_TS_CHEM, GET_MONTH
      USE TRACER_MOD,   ONLY : STT 
      USE TRACERID_MOD, ONLY : IDTMGLY, IDTSOAM
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44, ND07, LD07
#     include "comode.h"     ! AD, WTAIR, other SMVGEAR variables


      ! Local variables
      INTEGER              :: I, J, L, JLOOP, N
      TYPE (XPLEX)               :: XTEMP       ! Temperature [K]
      TYPE (XPLEX)               :: XSQTEMP     ! SQRT of Temperature
      TYPE (XPLEX)               :: XAD         ! Air density [molec/cm3]
      TYPE (XPLEX)               :: XRADIUS     ! particle radius [cm]
      TYPE (XPLEX)               :: XDFKG       ! Gas phase diffusion coeff [cm2/s]
      TYPE (XPLEX)               :: XARSL1K     ! 1st order k
      TYPE (XPLEX)               :: XRH         ! Relative humidity [%]
      TYPE (XPLEX)               :: XAIRM3      ! Air volume in grid box [m3]
      TYPE (XPLEX)               :: XGASM       ! Gas mass at grid box before uptake [kg]
      TYPE (XPLEX)               :: XGASC       ! Gas concentration at grid box before uptake [molec/cm3]
      TYPE (XPLEX)               :: XWAREA      ! Wet aerosol surface area at grid box [cm^2 wet sfc area of aerosol cm^-3 air]
      TYPE (XPLEX)               :: XUPTK0      ! Potential uptake of gas by aerosol in grid box by aerosol type [molec/cm3]
      TYPE (XPLEX)               :: XUPTK1      ! Potential uptake of gas by aerosol in grid box by aerosol type [kg]
      TYPE (XPLEX)               :: XUPTKSUM    ! Potential uptake of gas by aerosol in grid box [kg]
      TYPE (XPLEX)               :: XUPTK       ! Actual uptake of gas by aerosol in grid box [kg]
                                          !  XUPTK <= STT( I, J, L, IDTGLYX )
      TYPE (XPLEX)               :: XGAMMA      ! Uptake coefficient 

      ! Local variables not changing 
      TYPE (XPLEX)               :: DTCHEM      ! Chemistry time step [s]
      TYPE (XPLEX)               :: XMW         ! Molecular weight of gas [g/mole]
      TYPE (XPLEX)               :: XSQMW       ! Square root of molecular weight [g/mole]
      TYPE (XPLEX)               :: CRITRH      ! Critical RH [%], above which 
                                          !  heteorogeneous chem takes place
      TYPE (XPLEX)               :: XNAVO       ! Avogadro number

      !=================================================================
      ! SOAG_LIGGIO_DIFF begins here!
      !=================================================================

      ! Get chemistry time step
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Molecular weight of MGLY [g/mole]
      XMW   = 72.d0
      XSQMW = SQRT( XMW )

      ! Critical RH, above which heteorogeneous chem takes place
      CRITRH = 35.0d0   ! [%]

      ! Avogadro number
      XNAVO = 6.022d23

      ! Uptake coefficient from Liggio et al. [2005b]
      XGAMMA = 2.9d-3

      ! Create RH field -- relative humidity (in dao_mod.f)
!      CALL MAKE_RH

      !=================================================================
      ! Loop over grid boxes
      !=================================================================
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

            ! Get 1-D index
            JLOOP   = JLOP( I, J, L )

            ! Get RH  
            XRH     = RH( I, J, L )   ![%]

            ! initialize for safety
            XUPTK0   = 0d0
            XUPTK1   = 0d0
            XUPTKSUM = 0d0
            XUPTK    = 0d0

            ! Get T
            XTEMP   = T( I, J, L )
            XSQTEMP = SQRT( XTEMP )

            ! Get air density  [molec/cm3]
            XAD     = AIRDENS( JLOOP )

            ! Get air volumne [m3]
            XAIRM3  = AIRVOL( I, J, L )

            ! Get gas mass at grid box [kg]
            XGASM   = STT( I, J, L, IDTMGLY )

            ! Get gas concentration at grid box [molec/cm3]
            XGASC   = XGASM / (XMW*1.d-3) * XNAVO / (XAIRM3*1.d6) 

            !---------------------------------------
            ! Gas phase diffusion coeff [cm2/s]           
            !---------------------------------------
            XDFKG = 9.45D17 / XAD * XSQTEMP * 
     &              SQRT( 3.472D-2 + (1.D0/XMW) )


            !========================================================
            ! Calculate heteorogeneous uptake only if the grid box
            !  relative humidity XRH is >= critical relative humidity CRITRH
            !========================================================
            IF ( XRH >= CRITRH ) THEN

               ! Loop over sulfate and other aerosols
               DO N = 1, NDUST + NAER

                  !---------------------------------------
                  ! Total available wet aerosol area 
                  !  archived in 'aerosol_mod.f.glyx'
                  !  XWAREA [ cm^2 wet sfc area of aerosol cm^-3 air ]
                  !---------------------------------------
                  XWAREA  = WTAREA( JLOOP, N) 

                  IF ( XWAREA > 0D0 ) THEN 

                     ! Get particle radius [cm]
                     XRADIUS = WERADIUS( JLOOP, N )   

                     !---------------------------------------
                     ! First order rate constant
                     !---------------------------------------
                     XARSL1K = XWAREA / 
     &               (XRADIUS/XDFKG + 2.749064D-4*XSQMW/XGAMMA/XSQTEMP)

                     !---------------------------------------
                     ! Calculate potential uptake: Liggio et al. (2005b) Eq (3)
                     !   
                     !   d( organic carbon conc ) / dt = 
                     !      XARSL1K * XGASC 
                     !---------------------------------------
                     XUPTK0 = XARSL1K * XGASC * DTCHEM
                     XUPTK1 = XUPTK0 / XNAVO*(XMW*1.d-3)*(XAIRM3*1.d6)
                     XUPTKSUM = XUPTKSUM + XUPTK1

	            ENDIF
           
               ENDDO

                  ! However, the mass of gas being absorbed by aerosol 
                  !  cannot exceed the original amount of gas XGASM
                  XUPTK  = MIN( XUPTKSUM, XGASM )
            
                  ! Update MGLY in the STT array
                  STT( I, J, L, IDTMGLY ) = STT( I, J, L, IDTMGLY ) -
     &                                      XUPTK

                  ! Update SOAM in the STT array
                  STT( I, J, L, IDTSOAM ) = STT( I, J, L, IDTSOAM ) + 
     &                                      XUPTK

            ENDIF

         !==============================================================
         ! ND07 diagnostic: SOAM from MGLY [kg/timestep] on aerosol
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
            AD07_SOAGM(I,J,L,2) = AD07_SOAGM(I,J,L,2) + XUPTK
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      !=================================================================       
      ! Calculate dry-deposition
      !=================================================================
      CALL SOA_DEPO( STT(:,:,:,IDTSOAM), DRYSOAM, IDTSOAM )

      ! Return to calling program 
      END SUBROUTINE SOAM_LIGGIO_DIFF


!------------------------------------------------------------------------------

      SUBROUTINE SOA_CHEMISTRY
!
!******************************************************************************
!  Subroutine SOA_CHEMISTRY performs SOA formation.  This code is from the
!  Caltech group (Hong Liao, Serena Chung, et al) and was modified for 
!  GEOS-CHEM. (rjp, bmy, 7/8/04, 8/3/06)
!
!  Procedure:
!  ============================================================================
!  (1 ) Read in NO3, O3, OH in CHEM_SOA
!  (2 ) Scales these fields using OHNO3TIME in sulfate_mod.f (see GET_OH)
!  (3 ) Calculate reaction rates (Serena's OCHEMPARAETER)
!  (4 ) CALCULATE DELHC
!  (5 ) get T0M gas products
!  (6 ) equilibrium calculation
!
!  There are total of 42 tracers considered in this routine:
!
!  4 classes of primary carbonaceous aerosols:  
!     BCPI = Hydrophilic black carbon
!     OCPI = Hydrophilic organic carbon
!     BCPO = Hydrophobic black carbon
!     OCPO = Hydrophobic organic carbon
!
!  6 reactive biogenic hydrocarbon groups (NVOC):
!     ALPH = a-pinene, b-pinene, sabinene, carene, terpenoid ketones
!     LIMO = limonene
!     TERP = a-terpinene, r-terpinene, terpinolene
!     ALCO = myrcene, terpenoid alcohols, ocimene
!     SESQ = sesquiterpenes
!     ISOP = Isoprene
!
!  NOTE: TERP and SESQ are not tracers because of their high reactivity
!
!  32 organic oxidation products by O3+OH and NO3: 
!     6 ( 3 gases + 3 aerosols ) from each of first four NVOC  = 24 
!     4 ( 2 gases + 2 aerosols ) from sesquiterpenes oxidation = 4
!     4 ( 2 gases + 2 aerosols ) from isoprene oxidation       = 4
!
!  NOTE: We aggregate these into 6 tracers according to HC classes
!     SOG1 = lump of gas products of first three (ALPH+LIMO+TERP) HC oxidation.
!     SOG2 = gas product of ALCO oxidation
!     SOG3 = gas product of SESQ oxidation 
!     SOG4 = gas product of ISOP oxidation
!     SOA1 = lump of aerosol products of first 3 (ALPH+LIMO+TERP) HC oxidation.
!     SOA2 = aerosol product of ALCO oxidation
!     SOA3 = aerosol product of SESQ oxidation 
!     SOA4 = aerosol product of ISOP oxidation
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now modified for SOG4, SOA4 -- products of oxidation by isoprene.
!        (dkh, bmy, 6/1/06)
!  (3 ) Now consider SOG condensation onto SO4, NH4, NIT aerosols (if SO4,
!        NH4, NIT are defined as tracers). (rjp, bmy, 8/3/06) 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE DAO_MOD,      ONLY : T, AD, AIRVOL, SUNCOS
      USE DIAG_MOD,     ONLY : AD07_HC 
      USE TRACER_MOD,   ONLY : STT 
      USE TRACERID_MOD, ONLY : IDTALCO, IDTALPH, IDTLIMO, IDTOCPI
      USE TRACERID_MOD, ONLY : IDTOCPO, IDTSOA1, IDTSOA2, IDTSOA3
      USE TRACERID_MOD, ONLY : IDTSOA4, IDTSOG1, IDTSOG2, IDTSOG3
      USE TRACERID_MOD, ONLY : IDTSOG4, IDTSO4,  IDTNH4,  IDTNIT
      USE TIME_MOD,     ONLY : GET_TS_CHEM,      GET_MONTH
      USE TIME_MOD,     ONLY : ITS_TIME_FOR_BPCH
      USE LOGICAL_MOD,  ONLY : LDICARB

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I,        J,        L,        N    
      INTEGER               :: JHC,      IPR
      TYPE (XPLEX)                :: RTEMP,    VOL,      FAC,      MPOC 
      TYPE (XPLEX)                :: MNEW,     MSOA_OLD, MPRODUCT, CAIR
      TYPE (XPLEX)                :: LOWER,    UPPER,    TOL,      VALUE
      TYPE (XPLEX)                :: KO3(MHC), KOH(MHC), KNO3(MHC)
      TYPE (XPLEX)                :: ALPHA(NPROD,MHC),   KOM(NPROD,MHC)
      TYPE (XPLEX)                :: GM0(NPROD,MHC),     AM0(NPROD,MHC)
      TYPE (XPLEX)             :: ORG_AER(NPROD,MHC), ORG_GAS(NPROD,MHC)

      !=================================================================
      ! SOA_CHEMISTRY begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,        J,        L,     JHC,   IPR,   GM0,  AM0  )
!$OMP+PRIVATE( VOL,      FAC,      RTEMP, KO3,   KOH,   KNO3, CAIR )
!$OMP+PRIVATE( MPRODUCT, MSOA_OLD, VALUE, UPPER, LOWER, MNEW, TOL  )
!$OMP+PRIVATE( ORG_AER,  ORG_GAS,  ALPHA, KOM,   MPOC              )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Volume of grid box [m3]
         VOL    = AIRVOL(I,J,L)    

         ! conversion factor from kg to ug/m3
         FAC    = 1.D9 / VOL       

         ! air conc. in kg/m3
         CAIR   = AD(I,J,L) / VOL  

         ! Temperature [K]
         RTEMP  = T(I,J,L)         

         ! Get SOA yield parameters
         CALL SOA_PARA( I, J, L, RTEMP, KO3, KOH, KNO3, ALPHA, KOM )

         ! Partition mass of gas & aerosol tracers 
         ! according to 5 VOC classes & 3 oxidants
         CALL SOA_PARTITION( I, J, L, GM0, AM0 )

         ! Compute oxidation of hydrocarbons by O3, OH, NO3
         CALL CHEM_NVOC( I, J, L, KO3, KOH, KNO3, ALPHA, GM0 )

         !==============================================================
         ! Equilibrium calculation between GAS (SOG) and Aerosol (SOA)
         !==============================================================

         ! Total VOC oxidation products (gas and aerosol) [kg]
         MPRODUCT = 0.D0   

         ! Initialize aerosol-only total [kg]
         MSOA_OLD = 0.D0  
     
         ! Compute total VOC and aerosol-only total
         DO JHC = 1, MHC
         DO IPR = 1, NPROD
            MPRODUCT = MPRODUCT + GM0(IPR,JHC) + AM0(IPR,JHC)
            MSOA_OLD = MSOA_OLD + AM0(IPR,JHC)
         ENDDO
         ENDDO

         ! No need to proceed because there is no
         ! products that need to be re-equilibrated.
         IF ( ( MPRODUCT / AD(I,J,L) ) <= 29.D-18 ) CYCLE

         ! Individual SOA's; units: [ug/m3]
         DO JHC = 1, MHC
         DO IPR = 1, NPROD
            ORG_GAS(IPR,JHC) = GM0(IPR,JHC) * FAC
            ORG_AER(IPR,JHC) = AM0(IPR,JHC) * FAC
         ENDDO
         ENDDO
                
         IF ( IDTSO4 > 0 .and. IDTNH4 > 0 .and. IDTNIT > 0 ) THEN

            !-----------------------------------------------------------
            ! Condensation of SOG onto seed aerosols (rjp, bmy, 8/3/06)
            !
            ! If SO4, NH4, and NIT are defined tracers in this run
            ! (i.e. as in full-chemistry simulations), then:
            !  
            ! (1) Compute the condensation of SOG on seed aerosols 
            !     including preexisting OC and inorganic aerosols
            !
            ! (2) Use higher ratio (2.1) of molecular weight of 
            !     organic mass per carbon mass accounting for non-carbon 
            !     components attached to OC [Turpin and Lim, 2001]
            !
            ! NOTE...The remaining question is:
            !
            !     Do we have to include other types of aerosols such 
            !     as dust, sea salt, and aerosol water which should be 
            !     available in ISORROPIA?
            !
            ! I think that aerosol water is valuable information for 
            ! this and also AOD computation so at least we need to 
            ! archive it if including it as a tracer is too expensive.
            ! (Rokjin Park, 8/3/06)
            !-----------------------------------------------------------

            !-----------------------------------------------------------
            ! The standard code reversibly-partitions SOA mass onto all 
            ! aqueous aerosols.
            ! However, if the dicarbonyl SOA formation pathway is included, 
            ! this would lead to an overestimate of SOA mass compare to 
            ! observations during ICARTT. So LDICARB was introduced to 
            ! change to partitioning only onto pre-existing organic aerosols
            ! when adding the dicarbonyl SOA formation pathway. (tmf, 3/06/09)
            !-----------------------------------------------------------
 
	    IF ( .not. LDICARB ) THEN
               ! First compute SOG condensation onto OC aerosol
               MPOC = ( STT(I,J,L,IDTOCPI) + STT(I,J,L,IDTOCPO) ) * FAC
               MPOC = MPOC * 2.1d0
            ENDIF

            ! Then compute SOG condensation onto SO4, NH4, NIT aerosols
            MPOC = MPOC + ( STT(I,J,L,IDTSO4) + 
     &                      STT(I,J,L,IDTNH4) + 
     &                      STT(I,J,L,IDTNIT) ) * FAC

         ELSE

            !-----------------------------------------------------------
            ! If SO4, NIT, and NH3 are all undefined for this run, 
            ! then we will omit SOG condensation onto these aerosols 
            ! (at least for the time being).
            !
            ! We carry carbon mass only in STT array and here multiply 
            ! by the ratio (2.1) of molecular weight of organic mass
            ! per carbon mass.  This ratio accounts for non-carbon 
            ! components attached to OC [Turpin and Lim, 2001].
            !
            ! %%% NOTE: Ideally we should read in offline monthly    %%%
            ! %%% mean concentrations of SO4, NH4, NIT, and then use %%%
            ! %%% these to compute SOG condensation.  Implement this %%%
            ! %%% later.                                             %%%
            ! %%%                                                    %%%
            ! %%% For the time being, GEOS-Chem users are STRONGLY   %%%
            ! %%% encouraged to include sulfate aerosol and 2dary    %%%
            ! %%% organic aerosol tracers in offline carbon aerosol  %%%
            ! %%% simulations.  Without these extra aerosol tracers, %%%
            ! %%% the carbon aerosol simulation alone may under-     %%%
            ! %%% estimate SOA formation.                            %%%
            !
            ! (rjp, bmy, 8/3/06)
            !-----------------------------------------------------------
 
            ! Compute SOG condensation onto OC aerosol
            MPOC = ( STT(I,J,L,IDTOCPI) + STT(I,J,L,IDTOCPO) ) * FAC
            MPOC = MPOC * 2.1d0

         ENDIF

         !==============================================================
         ! Solve for MNEW by solving for SOA=0
         !==============================================================
         IF ( ( MPOC / ( CAIR*1.D9 ) ) <= 2.1D-18 ) THEN
            VALUE = 0.D0
            UPPER = 0.D0

            DO JHC = 1, MHC
            DO IPR = 1, NPROD
               VALUE = VALUE + KOM(IPR,JHC) *
     &                 (ORG_GAS(IPR,JHC) + ORG_AER(IPR,JHC))

               UPPER = UPPER + ORG_GAS(IPR,JHC) + ORG_AER(IPR,JHC)
            ENDDO
            ENDDO

            IF ( VALUE <= 1.D0 ) THEN
               MNEW  = 0.D0
            ELSE
               LOWER = 1.D-18 * ( CAIR * 1.D9 )
               TOL   = 1.D-18 
               MNEW  = ZEROIN(LOWER,UPPER,TOL,MPOC,ORG_AER,ORG_GAS,KOM)
            ENDIF

         ELSE

            UPPER = MPOC
            DO JHC = 1, MHC
            DO IPR = 1, NPROD
               UPPER = UPPER + ORG_GAS(IPR,JHC) + ORG_AER(IPR,JHC)
            ENDDO
            ENDDO

            LOWER = MPOC
            TOL   = 1.D-9*MPOC
            MNEW  = ZEROIN(LOWER,UPPER,TOL,MPOC,ORG_AER,ORG_GAS,KOM)

         ENDIF

         !==============================================================
         ! Equilibrium partitioning into new gas and aerosol 
         ! concentrations for individual contributions of SOA
         !==============================================================
         IF ( MNEW > 0.D0 ) THEN
             DO JHC = 1, MHC
             DO IPR = 1, NPROD

                ORG_AER(IPR,JHC) = KOM(IPR,JHC)*MNEW /
     &                            (1.D0 + KOM(IPR,JHC) * MNEW ) *
     &                            (ORG_AER(IPR,JHC) + ORG_GAS(IPR,JHC))

                IF ( KOM(IPR,JHC).NE.0D0 ) THEN
                    ORG_GAS(IPR,JHC) = ORG_AER(IPR,JHC) * 1.D8 /
     &                                 ( KOM(IPR,JHC) * MNEW * 1.D8 )
                ELSE
                   ORG_GAS(IPR,JHC) = 0.D0
                ENDIF

             ENDDO
             ENDDO

             ! STORE PRODUCT INTO T0M 
             DO JHC = 1, MHC
             DO IPR = 1, NPROD
                GM0(IPR,JHC) = ORG_GAS(IPR,JHC) / FAC
                AM0(IPR,JHC) = ORG_AER(IPR,JHC) / FAC
             ENDDO
             ENDDO

         !==============================================================
         ! Mnew=0.D0, all SOA evaporates to the gas-phase
         !==============================================================
         ELSE

             DO JHC = 1, MHC
             DO IPR = 1, NPROD
                GM0(IPR,JHC) = GM0(IPR,JHC) + AM0(IPR,JHC)
                AM0(IPR,JHC) = 1.D-18 * AD(I,J,L)
             ENDDO
             ENDDO

         ENDIF

         ! Lump SOA
         CALL SOA_LUMP( I, J, L, GM0, AM0 )
       
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================       
      ! Calculate dry-deposition
      !=================================================================
      CALL SOA_DEPO( STT(:,:,:,IDTALPH), DRYALPH, IDTALPH )
      CALL SOA_DEPO( STT(:,:,:,IDTLIMO), DRYLIMO, IDTLIMO )
      CALL SOA_DEPO( STT(:,:,:,IDTALCO), DRYALCO, IDTALCO )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG1), DRYSOG1, IDTSOG1 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG2), DRYSOG2, IDTSOG2 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG3), DRYSOG3, IDTSOG3 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG4), DRYSOG4, IDTSOG4 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA1), DRYSOA1, IDTSOA1 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA2), DRYSOA2, IDTSOA2 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA3), DRYSOA3, IDTSOA3 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA4), DRYSOA4, IDTSOA4 )

      ! Return to calling program
      END SUBROUTINE SOA_CHEMISTRY

!------------------------------------------------------------------------------

      FUNCTION SOA_EQUIL( MASS, MPOC, AEROSOL, GAS, KOM ) 
     &         RESULT( SOA_MASS )
!
!******************************************************************************
!  Subroutine SOA_EQUIL solves SOAeqn=0 to determine Mnew (= mass)
!  See Eqn (27) on page 70 of notes.  Originally written by Serena Chung at
!  Caltech, and modified for inclusion into GEOS-CHEM. (rjp, bmy, 7/8/04)
!
!  This version does NOT assume that the gas and aerosol phases are in 
!  equilibrium before chemistry; therefore, gas phase concentrations are 
!  needed explicitly.  The gas and aerosol phases are assumed to be in 
!  equilibrium after chemistry.
! 
!  Note: Unlike FUNCTION SOA, this function assumes no reactions.  It only 
!  considers the partitioning of existing products of VOC oxidation.
!
!  HC_JHC + OXID_IOXID - > 
!    alpha(1,IOXID,JHC) [SOAprod_gas(1,IOXID,JHC)+SOAprod(1,IOXID,JHC)]+
!    alpha(2,IOXID,JHC) [SOAprod_gas(2,IOXID,JHC)+SOAprod(2,IOXID,JHC)]
!
!  SOAprod_gas(IPR,IOXID,JHC) <--> SOAprod(IPR,IOXID,JHC)   
!                                           (aerosol phase)
!
!  w/ equilibrium partitioning:
!
!                                   SOAprod(IPR,IOXID,JHC)
!    SOAprod_gas(IPR,IOXID,JHC) = ------------------------
!                                     Kom(IPR,IOXID,JHC)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MASS    (TYPE (XPLEX)) : Pre-existing aerosol mass                [ug/m3]
!  (2 ) MPOC    (TYPE (XPLEX)) : POA Mass                                 [ug/m3]
!  (3 ) AEROSOL (TYPE (XPLEX)) : Aerosol concentration                    [ug/m3]
!  (4 ) GAS     (TYPE (XPLEX)) : Gas-phase concentration                  [ug/m3]
!  (5 ) KOM     (TYPE (XPLEX)) : Equilibrium gas-aerosol partition coeff. [m3/ug]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: MASS, MPOC         
      TYPE (XPLEX), INTENT(IN) :: AEROSOL(NPROD,MHC)
      TYPE (XPLEX), INTENT(IN) :: GAS(NPROD,MHC)
      TYPE (XPLEX), INTENT(IN) :: KOM(NPROD,MHC)

      ! Local variables
      INTEGER            :: JHC,   IPR
      TYPE (XPLEX)             :: VALUE, SOA_MASS

      !=================================================================
      ! SOA_EQUIL begins here!
      !=================================================================

      ! Equation (39) on page 139 of notes:
      VALUE = 0.D0

      ! 
      DO JHC = 1, MHC
      DO IPR = 1, NPROD
         VALUE = VALUE + KOM(IPR,JHC)                        / 
     &                   ( 1.D0 + KOM(IPR,JHC) * MASS      ) * 
     &                   ( GAS(IPR,JHC) + AEROSOL(IPR,JHC) )
      ENDDO                
      ENDDO   

      ! Compute SOA mass
      SOA_MASS = VALUE + ( 1.D5 * MPOC ) / ( 1.D5 * MASS ) - 1.0D0

      ! Return to calling program
      END FUNCTION SOA_EQUIL

!------------------------------------------------------------------------------

      FUNCTION ZEROIN(AX,BX,TOL,MPOC,AEROSOL,GAS,KOM) RESULT( MNEW )
!
!******************************************************************************
! NOTE: This function may be problematic -- it uses GOTO's, which are not
! good for parallelization. (bmy, 7/8/04)
!
! shc I got this code from http://www.netlib.org
!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result ( .ge. 0.0d0)
!
!
!  output..
!
!  zeroin abcissa approximating a zero of  f  in the interval ax,bx
!
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!  NOTES:
!  (1 ) Change dabs to ABS and dsign to SIGN, in order to avoid conflicts
!        with intrinsic function names on the PGI compiler. (bmy, 12/2/04)
!******************************************************************************
!
      TYPE (XPLEX), intent(in) :: ax,bx,tol
      TYPE (XPLEX), INTENT(IN) :: Mpoc      
      TYPE (XPLEX), INTENT(IN) :: Aerosol(NPROD,MHC), Gas(NPROD,MHC)
      TYPE (XPLEX), INTENT(IN) :: Kom(NPROD,MHC)

      !local variables
      TYPE (XPLEX)             :: MNEW
      TYPE (XPLEX)             :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
c
c  compute eps, the relative machine precision
c
      eps = 1.0d0
   10 eps = eps/2.0d0
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d0) go to 10
c
c initialization
c
      a  = ax
      b  = bx
      fa = SOA_equil( A, MPOC, Aerosol, GAS, Kom )
      fb = SOA_equil( B, MPOC, Aerosol, GAS, Kom ) 
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d

   30 if (ABS(fc) .ge. ABS(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*ABS(b) + 0.5d0*tol
      xm = 0.5D0*(c - b)
      if (ABS(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (ABS(e) .lt. tol1) go to 70
      if (ABS(fa) .le. ABS(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = ABS(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - ABS(tol1*q))) go to 70
      if (p .ge. ABS(0.5d0*e*q)) go to 70

      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (ABS(d) .gt. tol1) b = b + d
      if (ABS(d) .le. tol1) b = b + SIGN(tol1, xm)

      fb = SOA_equil( B, MPOC, Aerosol, GAS, Kom ) 
      if ((fb*(fc/ABS(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 MNEW = b

      ! Return to calling program
      END FUNCTION ZEROIN

!------------------------------------------------------------------------------

      FUNCTION RTBIS( X1,   X2,      XACC, 
     &                MPOC, AEROSOL, GAS, KOM ) RESULT( ROOT )
!
!******************************************************************************
!  Function RTBIS finds the root of the function SOA_EQUIL via the bisection
!  method.  Original algorithm from "Numerical Recipes" by Press et al, 
!  Cambridge UP, 1986.  Modified for inclusion into GEOS-CHEM. (bmy, 7/8/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X1      (TYPE (XPLEX)) : Endpoint #1 
!  (2 ) X2      (TYPE (XPLEX)) : Endpoint #2
!  (3 ) XACC    (TYPE (XPLEX)) : Desired accuracy of solution
!  (4 ) MPOC    (TYPE (XPLEX)) : POA Mass                                 [ug/m3]
!  (5 ) AEROSOL (TYPE (XPLEX)) : Aerosol concentration                    [ug/m3]
!  (6 ) GAS     (TYPE (XPLEX)) : Gas-phase concentration                  [ug/m3]
!  (7 ) KOM     (TYPE (XPLEX)) : Equilibrium gas-aerosol partition coeff. [m3/ug]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: X1, X2, XACC, MPOC
      TYPE (XPLEX), INTENT(IN) :: AEROSOL(NPROD,MHC)
      TYPE (XPLEX), INTENT(IN) :: GAS(NPROD,MHC)
      TYPE (XPLEX), INTENT(IN) :: KOM(NPROD,MHC)

      ! Local variables
      INTEGER, PARAMETER :: JMAX = 100
      INTEGER            :: J
      TYPE (XPLEX)             :: ROOT, DX, F, FMID, XMID

      !=================================================================
      ! RTBIS begins here!
      !=================================================================

      ! Compute value of function SOA_EQUIL at endpoints
      FMID = SOA_EQUIL( X2, MPOC, AEROSOL, GAS, KOM )
      F    = SOA_EQUIL( X1, MPOC, AEROSOL, GAS, KOM )

      ! Test if we are bracketing a root
      IF ( F * FMID >= 0d0 ) THEN
         CALL ERROR_STOP( 'Root must be bracketed!', 
     &                    'RTBIS ("carbon_mod.f")' )
      ENDIF

      ! Set initial root and interval
      IF ( F < 0d0 ) THEN
         ROOT = X1
         DX   = X2 - X1
      ELSE
         ROOT = X2
         DX   = X1 - X2
      ENDIF

      ! Loop until max iteration count
      DO J = 1, JMAX

         ! Halve the existing interval
         DX   = DX * 0.5D0

         ! Compute midpoint of new interval
         XMID = ROOT + DX

         ! Compute value of function SOA_EQUIL at new midpoint
         FMID = SOA_EQUIL( XMID, MPOC, AEROSOL, GAS, KOM )

         ! We have found the root!
         IF ( FMID <= 0D0 ) ROOT = XMID

         ! We have reached the tolerance, so return
         IF ( ABS( DX ) < XACC .OR. FMID == 0.D0 ) RETURN
      ENDDO

      ! Stop with error condition
      CALL ERROR_STOP( 'Too many bisections!', 
     &                 'RTBIS ("carbon_mod.f")' )

      ! Return to calling program
      END FUNCTION RTBIS

!------------------------------------------------------------------------------

      SUBROUTINE SOA_PARA( II,  JJ,  LL,   TEMP, 
     &                     KO3, KOH, KNO3, RALPHA, KOM )
!
!******************************************************************************
!  Subroutine SOA_PARA gves mass-based stoichiometric coefficients for semi-
!  volatile products from the oxidation of hydrocarbons.  It calculates 
!  secondary organic aerosol yield parameters.  Temperature effects are
!  included.  Original code from the CALTECH group and modified for inclusion
!  to GEOS-CHEM. (rjp, bmy, 7/8/04, 6/30/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) II     (INTEGER) : GEOS-Chem longitude index 
!  (2 ) JJ     (INTEGER) : GEOS-Chem latitude index
!  (3 ) LL     (INTEGER) : GEOS-Chem altitude index
!  (4 ) TEMP   (TYPE (XPLEX) ) : Temperature [k]
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) KO3    (TYPE (XPLEX) ) : Rxn rate for HC oxidation by O3        [cm3/molec/s]
!  (6 ) KOH    (TYPE (XPLEX) ) : Rxn rate for HC oxidation by OH        [cm3/molec/s]
!  (7 ) KNO3   (TYPE (XPLEX) ) : Rxn rate for HC oxidation by NO3       [cm3/molec/s]
!  (8 ) RALPHA (TYPE (XPLEX) ) : Mass-based stoichiometric coefficients  [unitless]
!  (9 ) KOM    (TYPE (XPLEX) ) : Equilibrium gas-aerosol partition coeff [m3/ug]
!
!  References:
!  ============================================================================
!  PHOTO-OXIDATION RATE CONSTANTS OF ORGANICS come from:
!  (1 ) Atkinson, el al., Int. J. Chem.Kinet., 27: 941-955 (1995)
!  (2 ) Shu and Atkinson, JGR 100: 7275-7281 (1995)
!  (3 ) Atkinson, J. Phys. Chem. Ref. Data 26: 215-290 (1997)   
!  (4 ) Some are reproduced in Table 1 of Griffin, et al., JGR 104: 3555-3567
!  (5 ) Chung and Seinfeld (2002)
!
!  ACTIVATION ENERGIES come from:
!  (6 ) Atkinson, R. (1994) Gas-Phase Tropospheric Chemistry of Organic 
!        Compounds.  J. Phys. Chem. Ref. Data, Monograph No.2, 1-216. 
!  (7 ) They are also reproduced in Tables B.9 and B.10 of Seinfeld and 
!        Pandis (1988).
!  
!  PARAMETERS FOR ISOPRENE:
!  (8 ) Kroll et al., GRL, 109, L18808 (2005)
!  (9 ) Kroll et al., Environ Sci Tech, in press (2006)
!  (10) Henze and Seinfeld, GRL, submitted (2006)
!
!  NOTES:
!  (1 ) Now use temporary variables TMP1, TMP2, TMP3 to pre-store the values
!        of exponential terms outside of DO-loops (bmy, 7/8/04)
!  (2 ) Add parameters for isoprene.  Now include grid cell location in
!        subroutine arguments.  Define a reference temperature at 295.
!        Now use ITS_IN_THE_TROP to determine if we are in a tropospheric
!        grid box.  Now pass II, JJ, LL via the argument list.
!        (dkh, bmy, 5/22/06)
!  (3 ) Corrected confusing documentation. (clh, bmy, 6/30/08)
!******************************************************************************
!
      ! References to F90 modules
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP

      ! Arguments
      INTEGER, INTENT(IN)  :: II, JJ, LL
      TYPE (XPLEX),  INTENT(IN)  :: TEMP
      TYPE (XPLEX),  INTENT(OUT) :: KO3(MHC), KOH(MHC), KNO3(MHC)
      TYPE (XPLEX),  INTENT(OUT) :: RALPHA(NPROD,MHC), KOM(NPROD,MHC)

      ! Local variables
      INTEGER              :: IPR,  JHC,  J
      TYPE (XPLEX)               :: TMP1, TMP2, TMP3, OVER

      ! Activation Energy/R [K] for O3, OH, NO3 (see Refs #6-7)
      TYPE (XPLEX), PARAMETER    :: ACT_O3     =  xplex(732.0d0 ,0d0)  
      TYPE (XPLEX), PARAMETER    :: ACT_OH     = xplex(-400.0d0,0d0)   
      TYPE (XPLEX), PARAMETER    :: ACT_NO3    = xplex(-490.0d0,0d0)   

      ! Heat of vaporization (from CRC Handbook of Chemistry & Physics)
      TYPE (XPLEX), PARAMETER    :: HEAT_VAPOR = xplex(5.d3,0d0)     

      ! Reciprocal reference temperatures at 298K and 310K
      TYPE (XPLEX), PARAMETER    :: REF295     = xplex(1d0 / 295d0,0d0)
      TYPE (XPLEX), PARAMETER    :: REF298     = xplex(1d0 / 298d0,0d0)
      TYPE (XPLEX), PARAMETER    :: REF310     = xplex(1d0 / 310d0,0d0)

      !=================================================================
      ! SOA_PARA begins here!
      !=================================================================

      ! Photo-oxidation rates of O3 [cm3/molec/s] (See Refs #1-4)
      KO3(1) = 56.15d-18
      KO3(2) = 200.d-18
      KO3(3) = 7707.d-18
      KO3(4) = 422.5d-18
      KO3(5) = ( 11600.D0 + 11700.D0 ) / 2.D0 * 1.D-18

      ! Photo-oxidation rates of OH [cm3/molec/s] (See Refs #1-4)
      KOH(1) = 84.4d-12
      KOH(2) = 171.d-12
      KOH(3) = 255.d-12
      KOH(4) = 199.d-12
      KOH(5) = ( 197.d0 + 293.d0 ) / 2.d0 * 1.d-12

      ! Photo-oxidation rate of NO3 [cm3/molec/s] (See Refs #1-4)
      KNO3(1) = 6.95d-12
      KNO3(2) = 12.2d-12
      KNO3(3) = 88.7d-12
      KNO3(4) = 14.7d-12
      KNO3(5) = ( 19.d0 + 35.d0 ) / 2.d0 * 1.d-12

      !=================================================================
      ! Temperature Adjustments of KO3, KOH, KNO3
      !=================================================================

      ! Reciprocal temperature [1/K]
      OVER = 1.0d0 / TEMP

      ! Compute the exponentials once outside the DO loop
      TMP1 = EXP( ACT_O3  * ( REF298 - OVER ) )
      TMP2 = EXP( ACT_OH  * ( REF298 - OVER ) )
      TMP3 = EXP( ACT_NO3 * ( REF298 - OVER ) )

      ! Multiply photo-oxidation rates by exponential of temperature
      DO JHC = 1, MHC 
         KO3(JHC)  = KO3(JHC)  * TMP1
         KOH(JHC)  = KOH(JHC)  * TMP2
         KNO3(JHC) = KNO3(JHC) * TMP3
      ENDDO

      ! If we are in the troposphere, then calculate ISOP oxidation rates
      ! If we are in the stratosphere, set ISOP oxidation rates to zero
      ! (dkh, bmy, 5/22/06)
      IF ( ITS_IN_THE_TROP( II, JJ, LL ) ) THEN
         KO3(6)  = 1.05d-14 * EXP(-2000.d0 / TEMP )
         KOH(6)  = 2.70d-11 * EXP( 390.d0 / TEMP )
         KNO3(6) = 3.03d-12 * EXP(-446.d0 / TEMP )
      ELSE
         KO3(6)  = 0d0
         KOH(6)  = 0d0
         KNO3(6) = 0d0
      ENDIF

      !=================================================================
      ! SOA YIELD PARAMETERS
      ! 
      ! Aerosol yield parameters for photooxidation of biogenic organics
      ! The data (except for C7-C10 n-carbonyls, aromatics, and higher 
      ! ketones are from: 
      !
      ! (7) Tables 1 and 2 of Griffin, et al., Geophys. Res. Lett. 
      !      26: (17)2721-2724 (1999)
      !
      ! These parameters neglect contributions of the photooxidation 
      ! by NO3. 
      !
      ! For the aromatics, the data are from
      ! (8) Odum, et al., Science 276: 96-99 (1997).
      !
      ! Isoprene (dkh, bmy, 5/22/06)
      ! Unlike the other species, we consider oxidation by purely OH. 
      ! CHEM_NVOC has been adjusted accordingly. There's probably 
      ! significant SOA formed from NO3 oxidation, but we don't know 
      ! enough to include that yet.  Data for the high NOX and low NOX 
      ! parameters are given in Kroll 05 and Kroll 06, respectively.  
      ! The paramters for low NOX are given in Table 1 of Henze 06.
      !=================================================================

      ! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
      RALPHA(1,1) = 0.067d0            
      RALPHA(2,1) = 0.35425d0

      ! LIMONENE
      RALPHA(1,2) = 0.239d0
      RALPHA(2,2) = 0.363d0

      ! Average of TERPINENES and TERPINOLENE
      RALPHA(1,3) = 0.0685d0
      RALPHA(2,3) = 0.2005d0

      ! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
      RALPHA(1,4) = 0.06675d0
      RALPHA(2,4) = 0.135d0

      ! Average of BETA-CARYOPHYLLENE and and ALPHA-HUMULENE
      RALPHA(1,5) = 1.0d0
      RALPHA(2,5) = 0.0d0

      ! Using BETA-PINENE for all species for NO3 oxidation
      ! Data from Table 4 of Griffin, et al., JGR 104 (D3): 3555-3567 (1999)
      RALPHA(3,:) = 1.d0           

      ! Here we define some alphas for isoprene (dkh, bmy, 5/22/06)

      ! high NOX  [Kroll et al, 2005]
      !RALPHA(1,6) = 0.264d0
      !RALPHA(2,6) = 0.0173d0
      !RALPHA(3,6) = 0d0

      ! low NOX   [Kroll et al, 2006; Henze and Seinfeld, 2006]
      RALPHA(1,6) = 0.232d0
      RALPHA(2,6) = 0.0288d0
      RALPHA(3,6) = 0d0

      !=================================================================
      ! Equilibrium gas-particle partition coefficients of 
      ! semi-volatile compounds [ug-1 m**3]
      !=================================================================

      ! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
      KOM(1,1) = 0.1835d0
      KOM(2,1) = 0.004275d0

      ! LIMONENE
      KOM(1,2) = 0.055d0
      KOM(2,2) = 0.0053d0

      ! Average of TERPINENES and TERPINOLENE
      KOM(1,3) = 0.133d0
      KOM(2,3) = 0.0035d0

      ! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
      KOM(1,4) = 0.22375d0
      KOM(2,4) = 0.0082d0

      ! Average of BETA-CARYOPHYLLENE and and ALPHA-HUMULENE
      KOM(1,5) = ( 0.04160d0 + 0.0501d0 ) / 2.d0
      KOM(2,5) = 0.0d0

      ! NOT APPLICABLE -- using BETA-PINENE for all species
      ! Data from Table 4 of Griffin, et al., JGR 104 (D3): 3555-3567 (1999)
      KOM(3,:) = 0.0163d0

      ! Again, for isoprene we only consider two products, 
      ! both from OH oxidation. (dkh, bmy, 5/22/06)

      ! High NOX
      !KOM(1,6) = 0.00115d0
      !KOM(2,6) = 1.52d0
      !KOM(3,6) = 0d0

      ! Low NOX
      KOM(1,6) = 0.00862d0
      KOM(2,6) = 1.62d0
      KOM(3,6) = 0d0
                           
      !=================================================================
      ! Temperature Adjustments of KOM
      !=================================================================

      ! Reciprocal temperature [1/K]
      OVER = 1.0D0 / TEMP

      ! Divide TEMP by 310K outside the DO loop
      TMP1 = ( TEMP / 310.D0 )

      ! Compute the heat-of-vaporization exponential term outside the DO loop
      TMP2 = EXP( HEAT_VAPOR * ( OVER - REF310 ) )

      ! Multiply KOM by the temperature and heat-of-vaporization terms
      DO JHC = 1, 5 
      DO IPR = 1, 3
         KOM(IPR,JHC) = KOM(IPR,JHC) * TMP1 * TMP2
      ENDDO
      ENDDO

      !--------------------------------------------------------
      ! For isoprene products, reference temperature is 295 K. 
      ! (dkh, bmy, 5/22/06)
      !--------------------------------------------------------

      ! Divide TEMP by 295K outside the DO loop
      TMP1 = ( TEMP / 295.D0 )

      ! Compute the heat-of-vaporization exponential term outside the DO loop
      TMP2 = EXP( HEAT_VAPOR * ( OVER - REF295 ) )

      ! Multiply KOM by the temperature and heat-of-vaporization terms
      JHC = 6
      DO IPR = 1, 3
         KOM(IPR,JHC) = KOM(IPR,JHC) * TMP1 * TMP2
      ENDDO

      ! Return to calling program
      END SUBROUTINE SOA_PARA

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NVOC( I, J, L, KO3, KOH, KNO3, ALPHA, GM0 )
!
!******************************************************************************
!  Subroutine CHEM_NVOC computes the oxidation of Hydrocarbon by O3, OH, and 
!  NO3.  This comes from the Caltech group (Hong Liao, Serena Chung, et al)
!  and was incorporated into GEOS-CHEM. (rjp, bmy, 7/6/04, 6/1/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER) : GEOS-Chem longitude index
!  (2 ) J      (INTEGER) : GEOS-Chem latitude index
!  (3 ) L      (INTEGER) : GEOS-Chem altitude index
!  (4 ) KO3    (TYPE (XPLEX) ) : Rxn rate for HC oxidation by O3        [cm3/molec/s]
!  (5 ) KOH    (TYPE (XPLEX) ) : Rxn rate for HC oxidation by OH        [cm3/molec/s]
!  (6 ) KNO3   (TYPE (XPLEX) ) : Rxn rate for HC oxidation by NO3       [cm3/molec/s]
!  (7 ) ALPHA  (TYPE (XPLEX) ) : Mass-based stoichiometric coefficients [unitless]  
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) GM0   (TYPE (XPLEX) ) : Gas mass for each HC and its oxidation product [kg]
!  
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Updated for SOA from isoprene.  Now calls GET_DOH. (dkh, bmy, 6/1/06)
!******************************************************************************
!
      ! References to F90 kmodules
      USE TRACER_MOD,    ONLY : STT
      USE TRACERID_MOD,  ONLY : IDTALCO,     IDTALPH,  IDTLIMO
      USE TIME_MOD,      ONLY : GET_TS_CHEM, GET_MONTH

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J, L
      TYPE (XPLEX),  INTENT(IN)    :: KO3(MHC), KOH(MHC), KNO3(MHC)
      TYPE (XPLEX),  INTENT(IN)    :: ALPHA(NPROD,MHC)
      TYPE (XPLEX),  INTENT(INOUT) :: GM0(NPROD,MHC)

      ! Local variables
      INTEGER                :: JHC, IPR
      TYPE (XPLEX)             :: DELHC(NPROD), CHANGE(MHC), NMVOC(MHC)
      TYPE (XPLEX)                 :: OHMC, TTNO3, TTO3, DTCHEM, RK
      TYPE (XPLEX)                 :: OVER, DO3, DOH, DNO3

      !=================================================================
      ! CHEM_NVOC begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM  = GET_TS_CHEM() * 60d0 

      ! Get offline OH, NO3, O3 concentrations [molec/cm3]
      OHMC    = GET_OH(  I, J, L ) 
      TTNO3   = GET_NO3( I, J, L )
      TTO3    = GET_O3(  I, J, L ) 

      ! 6 classes NVOC concentrations followed by 4 primary species
      NMVOC(1) = STT(I,J,L,IDTALPH)
      NMVOC(2) = STT(I,J,L,IDTLIMO)
      NMVOC(3) = ORVC_TERP(I,J,L)
      NMVOC(4) = STT(I,J,L,IDTALCO)
      NMVOC(5) = ORVC_SESQ(I,J,L)

      ! Initialize DELHC so that the values from the previous
      ! time step are not carried over.
      DELHC(:) = 0.D0

      !=================================================================
      ! Change in NVOC concentration due to photooxidation [kg]
      !=================================================================
      DO JHC = 1, MHC

         IF ( JHC <= 5 ) THEN

            !-------------------------------
            ! Individual gas-phase products
            !-------------------------------

            RK          = KO3(JHC)*TTO3 + KOH(JHC)*OHMC
     &                  + KNO3(JHC)*TTNO3
            CHANGE(JHC) = NMVOC(JHC) * ( 1.D0 - EXP( -RK * DTCHEM ) )

            ! In case that the biogenic hydrocarbon is the limiting reactant
            IF ( CHANGE(JHC) >= NMVOC(JHC) ) CHANGE(JHC) = NMVOC(JHC)
      
            ! NMVOC concentration after oxidation reactions
            NMVOC(JHC) = NMVOC(JHC) - CHANGE(JHC)

            IF( CHANGE(JHC) > 1.D-16 ) THEN
               OVER  = 1.D0 / RK
               DO3   = CHANGE(JHC) * KO3(JHC)  * TTO3  * OVER ![kg]
               DOH   = CHANGE(JHC) * KOH(JHC)  * OHMC  * OVER ![kg]
               DNO3  = CHANGE(JHC) * KNO3(JHC) * TTNO3 * OVER ![kg]
            ELSE
               DO3   = 0.D0
               DOH   = 0.D0
               DNO3  = 0.D0
            ENDIF
                  
            ! VOC change by photooxidation of O3 and OH [kg]
            DELHC(1) =  DO3 + DOH 
            DELHC(2) =  DO3 + DOH 
    
            ! VOC change by photooxidation of NO3 [kg]
            DELHC(3) = DNO3
          
            ! Lump OH and O3 oxidation for HC 1-5
            DO IPR = 1, NPROD
               GM0(IPR,JHC) = GM0(IPR,JHC) + ALPHA(IPR,JHC) * DELHC(IPR)
            ENDDO
      
         ELSE

            !-------------------------------
            ! SOA from ISOPRENE
            !-------------------------------

            ! Get ISOP lost to rxn with OH [kg]
            DOH = GET_DOH( I, J, L )

            ! Consider only OH oxidation for isoprene.  Also convert 
            ! from mass of carbon to mass of isoprene. (dkh, bmy, 5/22/06)
            DO IPR = 1, 3
               GM0(IPR,JHC) = GM0(IPR,JHC) + ALPHA(IPR,JHC) * DOH
     &                                     * 68d0           / 60d0 
            ENDDO
         ENDIF

      ENDDO                     

      !=================================================================
      ! Store Hydrocarbon remaining after oxidation rxn back into STT
      !=================================================================
      STT(I,J,L,IDTALPH) = MAX( NMVOC(1), 1.D-32 )
      STT(I,J,L,IDTLIMO) = MAX( NMVOC(2), 1.D-32 )
      ORVC_TERP(I,J,L)   = MAX( NMVOC(3), 1.D-32 )
      STT(I,J,L,IDTALCO) = MAX( NMVOC(4), 1.D-32 )
      ORVC_SESQ(I,J,L)   = MAX( NMVOC(5), 1.D-32 )

      ! Return to calling program
      END SUBROUTINE CHEM_NVOC

!------------------------------------------------------------------------------

      SUBROUTINE SOA_PARTITION( I, J, L, GM0, AM0 )
!
!******************************************************************************
!  Subroutine SOA_PARTITION partitions the mass of gas and aerosol 
!  tracers according to five Hydrocarbon species and three oxidants.
!  (rjp, bmy, 7/7/04, 5/22/06)
!
!  NOTE: GPROD and APROD are mass ratios of individual oxidation 
!        products of gas/aerosol to the sum of all. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I   (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J   (INTEGER) : GEOS-CHEM latitude index
!  (3 ) L   (INTEGER) : GEOS-CHEM altitude index
! 
!  Arguments as Output:
!  ============================================================================
!  (4 ) GM0 (TYPE (XPLEX) ) : Gas mass for each HC and its oxidation product     [kg]
!  (5 ) AM0 (TYPE (XPLEX) ) : Aerosol mass for each HC and its oxidation product [kg]
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Updated for SOG4, SOA4 (bmy, 5/22/06)
!******************************************************************************
!
      ! Refrences to F90 modules
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSOA1, IDTSOA2, IDTSOA3, IDTSOA4
      USE TRACERID_MOD, ONLY : IDTSOG1, IDTSOG2, IDTSOG3, IDTSOG4

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L
      TYPE (XPLEX),  INTENT(OUT) :: GM0(NPROD,MHC), AM0(NPROD,MHC)

      ! Local variables
      INTEGER              :: JHC, IPR

      !=================================================================
      ! SOA_PARTITION begins here!
      !=================================================================

      ! Initialize
      DO JHC = 1, 3
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = 0D0
         AM0(IPR,JHC) = 0D0
      ENDDO
      ENDDO

      ! Partition the lump of first three HC (ALPH + LIMO + TERP)
      ! oxidation products.  These are grouped together because they
      ! have the same molecular weight.
      DO JHC = 1, 3
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG1) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA1) * APROD(I,J,L,IPR,JHC)
      ENDDO                       
      ENDDO

      ! Alcohol
      JHC = 4
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG2) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA2) * APROD(I,J,L,IPR,JHC)
      ENDDO

      ! Sesqterpene
      JHC = 5
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG3) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA3) * APROD(I,J,L,IPR,JHC)
      ENDDO

      ! Isoprene
      JHC = 6
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG4) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA4) * APROD(I,J,L,IPR,JHC)
      ENDDO
      
      ! Return to calling program
      END SUBROUTINE SOA_PARTITION

!------------------------------------------------------------------------------

      SUBROUTINE SOA_LUMP( I, J, L, GM0, AM0 )
!
!******************************************************************************
!  Subroutine SOA_LUMP returns the organic gas and aerosol back to the
!  STT array.  (rjp, bmy, 7/7/04, 2/6/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I   (INTEGER) : Longitude index
!  (2 ) J   (INTEGER) : Latitude index
!  (3 ) L   (INTEGER) : Altitude index
!  (4 ) GM0 (TYPE (XPLEX) ) : Gas mass for each HC and its oxidation product     [kg]
!  (5 ) AM0 (TYPE (XPLEX) ) : Aerosol mass for each HC and its oxidation product [kg]
! 
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Bug fix: make sure L <= LD07 before saving into AD07 array, or else
!        we will get an out-of-bounds error. (bmy, 3/4/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Updated for SOG4, SOA4 (dkh, bmy, 5/22/06)
!  (5 ) Typo fix: GPROD should be APROD in a couple places (tmf, bmy, 10/16/06)
!  (6 ) Bug fix: For SOA4, GPROD and APROD should have default values of 0.5,
!        instead of 1.0 (dkh, bmy, 2/6/07)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD07_HC 
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSOA1, IDTSOA2, IDTSOA3, IDTSOA4
      USE TRACERID_MOD, ONLY : IDTSOG1, IDTSOG2, IDTSOG3, IDTSOG4

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! ND44, ND07, LD07
      
      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L
      TYPE (XPLEX),  INTENT(IN)  :: GM0(NPROD,MHC), AM0(NPROD,MHC)

      ! Local variables
      INTEGER              :: JHC, IPR
      TYPE (XPLEX)               :: GASMASS, AERMASS

      !=================================================================
      ! SOA_LUMP begins here!
      !=================================================================

      ! Initialize
      GASMASS = 0D0
      AERMASS = 0D0

      ! Compute total gas & aerosol mass
      DO JHC = 1, 3
      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)               
      ENDDO                       
      ENDDO

      !----------------------------
      ! SOA1 net Production (kg)
      !----------------------------
      IF ( ND07 > 0 .and. L <= LD07 ) THEN
         AD07_HC(I,J,L,1) = AD07_HC(I,J,L,1)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA1) )
      ENDIF

      STT(I,J,L,IDTSOG1) = MAX( GASMASS, 1.D-32 )
      STT(I,J,L,IDTSOA1) = MAX( AERMASS, 1.D-32 )

      IF ( STT(I,J,L,IDTSOG1) > 1.0D-6 ) THEN
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG1)
         ENDDO                       
         ENDDO
      ELSE
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC)= 0.111111111d0
         ENDDO
         ENDDO   
      ENDIF

      IF ( STT(I,J,L,IDTSOA1) > 1.0D-6 ) THEN
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA1)
         ENDDO                       
         ENDDO
      ELSE
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = 0.111111111d0
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! Lump of products of fourth Hydrocarbon class (ALCOHOL)
      !=================================================================

      JHC     = 4
      GASMASS = 0D0
      AERMASS = 0D0
      
      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)               
      ENDDO

      !---------------------------
      ! SOA2 net Production (kg)
      !---------------------------
      IF ( ND07 > 0 .and. L <= LD07 ) THEN
         AD07_HC(I,J,L,2) = AD07_HC(I,J,L,2)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA2) )
      ENDIF

      STT(I,J,L,IDTSOG2) = MAX(GASMASS, 1.D-32)
      STT(I,J,L,IDTSOA2) = MAX(AERMASS, 1.D-32)

      IF ( STT(I,J,L,IDTSOG2) > 1.0D-6 ) THEN
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG2)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = 0.333333333d0
         ENDDO
      ENDIF
      
      IF ( STT(I,J,L,IDTSOA2) > 1.0D-6 ) THEN
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA2)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = 0.333333333d0
         ENDDO
      ENDIF

      !=================================================================
      ! Lump of products of fifth Hydrocarbon class (SESQTERPINE)
      !=================================================================
      JHC     = 5
      GASMASS = 0D0
      AERMASS = 0D0

      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)               
      ENDDO

      !---------------------------
      ! SOA3 net Production (kg)
      !---------------------------
      IF ( ND07 > 0 .and. L <= LD07 ) THEN
         AD07_HC(I,J,L,3) = AD07_HC(I,J,L,3)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA3) )
      ENDIF

      STT(I,J,L,IDTSOG3) = MAX(GASMASS, 1.D-32)
      STT(I,J,L,IDTSOA3) = MAX(AERMASS, 1.D-32)

      IF ( STT(I,J,L,IDTSOG3) > 1.0D-6 ) THEN
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG3)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = 0.5D0
         ENDDO
      ENDIF

      IF ( STT(I,J,L,IDTSOA3) > 1.0D-6 ) THEN
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA3)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = 0.5D0
         ENDDO
      ENDIF

      ! make sure there is no second oxidation product 
      ! for SESQTERPENE by OH + O3
      GPROD(I,J,L,2,JHC) = 0.D0
      APROD(I,J,L,2,JHC) = 0.D0


      !=================================================================
      ! Lump of products of sixth Hydrocarbon class (ISOP) 
      ! (dkh, bmy, 5/22/06)
      !=================================================================
      JHC     = 6
      GASMASS = 0D0
      AERMASS = 0D0

      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)
      ENDDO

      !---------------------------
      ! SOA4 net Production (kg)
      !---------------------------

      IF ( ND07 > 0 .and. L <= LD07 ) THEN
         AD07_HC(I,J,L,4) = AD07_HC(I,J,L,4)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA4) )
      ENDIF

      STT(I,J,L,IDTSOG4) = MAX( GASMASS, 1.D-32 )
      STT(I,J,L,IDTSOA4) = MAX( AERMASS, 1.D-32 )

      IF ( STT(I,J,L,IDTSOG4) > 1.0D-6 ) THEN
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG4)
         ENDDO
      ELSE
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = 0.5d0   ! only one product
         ENDDO
      ENDIF

      IF ( STT(I,J,L,IDTSOA4) > 1.0D-6 ) THEN
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA4)
         ENDDO
      ELSE
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = 0.5d0
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE SOA_LUMP

!------------------------------------------------------------------------------

      SUBROUTINE SOA_DEPO( TC, DEPID, TRID )
!
!******************************************************************************
!  Subroutine SOA_DEPO computes dry-deposition of a particular SOA species.
!  (rjp, bmy, 7/8/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC    (TYPE (XPLEX) ) : Array of SOA tracer 
!  (2 ) DEPID (INTEGER) : Dry deposition ID # (from DEPVEL) 
!  (3 ) TRID  (INTEGER) : GEOS-CHEM tracer number 
! 
!  NOTES:
!  (1 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!  (2 ) Replace PBLFRAC from "drydep_mod.f" with  GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f" (bmy, 2/17/05)
!  (3 ) Bug fix: Add BL_FRAC to the PRIVATE list (mak, bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : XNUMOL

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)
      INTEGER, INTENT(IN)    :: DEPID, TRID

      ! Local variable
      INTEGER                :: I, J, L
      TYPE (XPLEX)                 :: DTCHEM, FLUX, BL_FRAC
      TYPE (XPLEX)                 :: TC0, CNEW, FREQ, AREA_CM2
      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! SOA_DEPO begins here!
      !=================================================================

      ! Return if tracer ID or tracer ID is undefined
      IF ( TRID == 0 .OR. DEPID == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for drydep diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, BL_FRAC, FREQ, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial SOA [kg]
         TC0 = TC(I,J,L)

         ! Fraction of box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DEPID) * BL_FRAC

            ! Amount of SOA[G] left after drydep [kg]
            CNEW = TC0 * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
               FLUX = ( TC0 - CNEW ) 
               FLUX = FLUX * XNUMOL(TRID) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            CNEW = TC0

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Store modified concentration back in tracer array [kg]
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DEPID,1) = AD44(I,J,DEPID,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF    

      ! Return to calling program
      END SUBROUTINE SOA_DEPO

!------------------------------------------------------------------------------

      SUBROUTINE EMISSCARBON
!
!******************************************************************************
!  Subroutine EMISSCARBON is the interface between the GEOS-CHEM model
!  and the CARBONACEOUS AEROSOL emissions (rjp, bmy, 1/24/02, 9/25/06)
!
!  NOTES:
!  (1 ) Now references LSOA from "CMN_SETUP".  Also now call OHNO3TIME since
!        biogenic emissions also have a diurnal variation. (rjp, bmy, 7/15/04)
!  (2 ) Now references LSOA and LPRT from "logical_mod.f".  Now references
!        STT from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Bug fix: removed "," from FORMAT 111.  Also added extra DEBUG_MSG
!        output after calling emissions routines. (bmy, 11/19/04)
!  (4 ) Now always call ANTHRO_CARB_TBOND and ANTHRO_CARB_COOKE.  This will
!        read the T. Bond et al [2004] emissions but overwrite the North
!        America region with monthly-mean emissions from Cooke et al [1999] 
!        with imposed seasonality from R. Park [2003].  (bmy, 12/1/04)
!  (5 ) Now remove THISMONTH from the arg list to BIOMASS_CARB_GEOS 
!        (bmy, 9/25/06)
!  (6 ) Now check that GFED2 has been updated if we do not use the annual
!        Bond Biomass emission (phs, yc, 12/18/08)
!  (7 ) adj_group: Now add scaling factors for adjoint (dkh, 11/10/09) 
!  (8 ) Now reads monthly (eml, phs, 5/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,          ONLY : AD07
      USE DAO_MOD,           ONLY : PBL
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE LOGICAL_MOD,       ONLY : LSOA,      LPRT, LCOOKE
      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : STT
      USE GFED2_BIOMASS_MOD, ONLY : GFED2_IS_NEW
      !USE TRACERID_MOD
     
      ! adj_group: add for emissions scale factors
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EBCPI_an, IDADJ_EBCPO_an
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EOCPI_an, IDADJ_EOCPO_an
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EBCPI_bb, IDADJ_EBCPO_bb
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EOCPI_bb, IDADJ_EOCPO_bb
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EBCPI_bf, IDADJ_EBCPO_bf
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EOCPI_bf, IDADJ_EOCPO_bf
      USE ADJ_ARRAYS_MOD,    ONLY : IS_CARB_EMS_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! ND07

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J, MONTH, N
      TYPE (XPLEX)               :: BCSRC(IIPAR,JJPAR,2)
      TYPE (XPLEX)               :: OCSRC(IIPAR,JJPAR,2)

      !=================================================================
      ! EMISSCARBON begins here!
      !
      ! Read carbonaceous aerosols from disk and compute hydrophilic 
      ! and hydrophobic fractions. NOTE, CARBON AEROSOLS HAVE TO BE 
      ! ORDERED AS Hydrophilic(BC[1], OC[2]) Hydrophobic(BC[3], OC[4]).
      !=================================================================      

      !--------------------------
      ! Read time-invariant data
      !--------------------------
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100 )

         ! Echo info about ANTHRO emissionls
         WRITE( 6, 110 )
         IF ( LCOOKE ) WRITE( 6, 111 )
         WRITE( 6, 112 )
         WRITE( 6, 113 )

         ! Monthly or annual BIOMASS emissions?
         IF ( USE_BOND_BIOBURN ) THEN
            WRITE( 6, 120 )
         ELSE
            WRITE( 6, 130 )
         ENDIF
         
         ! Write spacer
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMAT strings
 100     FORMAT( 'C A R B O N   A E R O S O L   E M I S S I O N S'    )
 110     FORMAT( 'w/ ANTHROPOGENIC emissions from Bond et al [2007]'  )
 111     FORMAT( 'w/ North American emissions from Cooke et al [1999]')
 112     FORMAT( 'w/ North American emissions having imposed'         )
 113     FORMAT( '   seasonality following Park et al [2003]'         )
 120     FORMAT( 'w/ BIOMASS emissions from Bond et al [2004]'        )
 130     FORMAT( 'w/ BIOMASS emissions from GEOS-CHEM inventory'      )


         ! Initialize arrays
         CALL INIT_CARBON       
         
!----------------------
! prior to 5/18/09 - now reads monthly data below (phs)         
!         ! Read annual mean anthro emissions from T. Bond [2004]
!         CALL ANTHRO_CARB_TBOND
!         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a A_CRB_TBOND' )
!---------------------

         ! Read annual mean biomass emissions if necessary
         IF ( USE_BOND_BIOBURN ) THEN
            CALL BIOMASS_CARB_TBOND
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a B_CRB_TBOND' )
         ENDIF

         ! Reset flag
         FIRST = .FALSE.
      ENDIF

      ! Compute time scaling arrays which are used both for
      ! biogenic emission and offline OH (rjp, bmy, 7/15/04)
      IF ( LSOA ) THEN
         CALL OHNO3TIME
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: after OHNO3TIME' )
      ENDIF
      
      !------------------------------
      ! Read monthly-mean ANTHRO data
      !------------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN
      
         ! Current month
         MONTH = GET_MONTH()

         ! Read in monthly emissions from T Bond [2007], with imposed
         ! seasonality over North America by R. Park [2003]
         CALL ANTHRO_CARB_TBOND( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a A_CRB_TBOND' )

         IF ( LCOOKE ) THEN
         ! Overwrite the T. Bond [2004] emissions over North America
         ! with monthly mean anthro emissions from Cooke et al [1999] 
         ! having imposed seasonality by R. Park [2003]
            CALL ANTHRO_CARB_COOKE( MONTH )
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a A_CRB_COOKE' )
         ENDIF

!-----------------------------------
! Prior to 12/18/08
!         ! Read monthly mean biomass emissions
!         IF ( .not. USE_BOND_BIOBURN ) THEN
!            CALL BIOMASS_CARB_GEOS
!            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a B_CRB_COOKE' )
!         ENDIF
!-----------------------------------

      ENDIF

      !-----------------------------------
      ! Read monthly/8-day/3-hr mean biomass emissions
      !-----------------------------------
      IF ( .not. USE_BOND_BIOBURN ) THEN

         IF ( GFED2_IS_NEW() .or. ITS_A_NEW_MONTH() ) THEN 
            CALL BIOMASS_CARB_GEOS
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a BB_CRB_GEOS' )
         ENDIF
            
      ENDIF

           
      !--------------------------
      ! Compute biogenic OC
      !--------------------------
      CALL BIOGENIC_OC
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: after BIOGENIC_OC' )

      !=================================================================
      ! Sum up BC and OC sources. 
      ! N=1 is HYDROPHILIC; N=2 is HYDROPHOBIC.
      !
      ! COMMENT: Maybe someday we'll want to play with the different 
      ! emission height for different source type.  For example the
      ! carbon from biomass burning could be emitted to the higher 
      ! altitude due to the thermal bouyancy and shallow convection.
      ! The current setting to use EMITHIGH seems rather inefficient 
      ! but robust for sensitivity studies for emission height 
      ! variation on carbon concentrations, so please keep using the 
      ! current setup until we decide otherwise. (rjp, 4/2/02)
      !=================================================================

      ! adj_group:  make a version that applies scaling factors 
      ! and use this if the emissions adjoint ID #'s are active 
      IF ( LADJ_EMS .and. IS_CARB_EMS_ADJ ) THEN 
        
         IF ( MMSCL /= 1 ) CALL ERROR_STOP('MMSCL /= 1','carbon_mod.f')

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Total HYDROPHILIC BC source [kg]
           BCSRC(I,J,1)= ANTH_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_an)
     &                 + BIOF_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bf)
     &                 + BIOB_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bb)

            ! Total HYDROPHOBIC BC source [kg]
           BCSRC(I,J,2)= ANTH_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_an)
     &                 + BIOF_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bf)
     &                 + BIOB_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bb)

            IF ( LSOA ) THEN

               CALL ERROR_STOP('LSOA not supported for adj', 
     &                         'carbon_mod.f')

               !! Total HYDROPHILIC OC source [kg]
               !! (Don't use archived TERP_ORGC if LSOA=T)
               !OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
     &         !               BIOF_ORGC(I,J,1) + 
     &         !               BIOB_ORGC(I,J,1)

            ELSE

               ! Total HYDROPHILIC OC source [kg]
               ! (Use archived TERP_ORGC for if LSOA=F)
               OCSRC(I,J,1)
     &            = ANTH_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_an)
     &            + BIOF_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bf)
     &            + BIOB_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bb)
     &            + TERP_ORGC(I,J)

            ENDIF

            ! Total HYDROPHOBIC OC source [kg]
            OCSRC(I,J,2) 
     &         = ANTH_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_an)
     &         + BIOF_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bf)
     &         + BIOB_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bb)
         ENDDO
         ENDDO

!$OMP END PARALLEL DO

      ELSE 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Total HYDROPHILIC BC source [kg]
            BCSRC(I,J,1) = ANTH_BLKC(I,J,1) + 
     &                     BIOF_BLKC(I,J,1) + 
     &                     BIOB_BLKC(I,J,1)   

            ! Total HYDROPHOBIC BC source [kg]
            BCSRC(I,J,2) = ANTH_BLKC(I,J,2) +
     &                     BIOF_BLKC(I,J,2) +
     &                     BIOB_BLKC(I,J,2)  

            IF ( LSOA ) THEN

               ! Total HYDROPHILIC OC source [kg]
               ! (Don't use archived TERP_ORGC if LSOA=T)
               OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
     &                        BIOF_ORGC(I,J,1) + 
     &                        BIOB_ORGC(I,J,1)

            ELSE

               ! Total HYDROPHILIC OC source [kg]
               ! (Use archived TERP_ORGC for if LSOA=F)
               OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
     &                        BIOF_ORGC(I,J,1) + 
     &                        BIOB_ORGC(I,J,1) + 
     &                        TERP_ORGC(I,J)

            ENDIF

            ! Total HYDROPHOBIC OC source [kg]
            OCSRC(I,J,2) = ANTH_ORGC(I,J,2) + 
     &                     BIOF_ORGC(I,J,2) + 
     &                     BIOB_ORGC(I,J,2) 
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF  ! EMS_ADJ
 
      ! Sum up all carbon tracers throughout the boundary layer
      CALL EMITHIGH( BCSRC, OCSRC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARB: after EMITHIGH' )

      !=================================================================
      ! ND07 diagnostic: Carbon aerosol emissions [kg/timestep]
      !=================================================================
      ! adj_group:  include scaling factors if they are activated 
      IF ( ND07 > 0 .and. IS_CARB_EMS_ADJ ) THEN

         IF ( MMSCL /= 1 ) CALL ERROR_STOP('MMSCL /= 1','carbon_mod.f')

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Anthropogenic BC source
            AD07(I,J,1) = AD07(I,J,1)  
     &         + ANTH_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_an)
     &         + ANTH_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_an)
            
            ! Biogenic BC source
            AD07(I,J,2) = AD07(I,J,2)        
     &         + BIOB_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bb)
     &         + BIOB_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bb)

            ! Biofuel BC source
            AD07(I,J,3) = AD07(I,J,3)        
     &         + BIOF_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bf)
     &         + BIOF_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bf)

            ! Anthropogenic OC source
            AD07(I,J,4) = AD07(I,J,4)        
     &         + ANTH_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_an)
     &         + ANTH_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_an)

            ! Biomass OC source
            AD07(I,J,5) = AD07(I,J,5)        
     &         + BIOB_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bb)
     &         + BIOB_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bb)

            ! Biofuel OC source
            AD07(I,J,6) = AD07(I,J,6)         
     &         + BIOF_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bf)
     &         + BIOF_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bf)

            ! Terpene source
            AD07(I,J,7) = AD07(I,J,7)        + TERP_ORGC(I,J)

            IF ( LSOA ) THEN

               CALL ERROR_STOP('LSOA not supported for adj', 
     &                         'carbon_mod.f')
!               ! ALPHA-PINENE
!               AD07(I,J,8)  = AD07(I,J,8)    + BIOG_ALPH(I,J)
!
!               ! LIMONENE
!               AD07(I,J,9)  = AD07(I,J,9)    + BIOG_LIMO(I,J)
!
!               ! TERPENE
!               AD07(I,J,10) = AD07(I,J,10)   + BIOG_TERP(I,J)
!
!               ! ALCOHOL
!               AD07(I,J,11) = AD07(I,J,11)   + BIOG_ALCO(I,J)
!
!               ! SESQTERPENE
!               AD07(I,J,12) = AD07(I,J,12)   + BIOG_SESQ(I,J)

            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARB: after ND07' )

      ELSEIF ( ND07 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Anthropogenic BC source
            AD07(I,J,1) = AD07(I,J,1)        +
     &                    ( ANTH_BLKC(I,J,1) + ANTH_BLKC(I,J,2) )
            
            ! Biogenic BC source
            AD07(I,J,2) = AD07(I,J,2)        +
     &                    ( BIOB_BLKC(I,J,1) + BIOB_BLKC(I,J,2) )

            ! Biofuel BC source
            AD07(I,J,3) = AD07(I,J,3)        +
     &                    ( BIOF_BLKC(I,J,1) + BIOF_BLKC(I,J,2) )

            ! Anthropogenic OC source
            AD07(I,J,4) = AD07(I,J,4)        +
     &                    ( ANTH_ORGC(I,J,1) + ANTH_ORGC(I,J,2) )

            ! Biomass OC source
            AD07(I,J,5) = AD07(I,J,5)        +
     &                    ( BIOB_ORGC(I,J,1) + BIOB_ORGC(I,J,2) )

            ! Biofuel OC source
            AD07(I,J,6) = AD07(I,J,6)        + 
     &                    ( BIOF_ORGC(I,J,1) + BIOF_ORGC(I,J,2) )

            ! Terpene source
            AD07(I,J,7) = AD07(I,J,7)        + TERP_ORGC(I,J)

            IF ( LSOA ) THEN

               ! ALPHA-PINENE
               AD07(I,J,8)  = AD07(I,J,8)    + BIOG_ALPH(I,J)

               ! LIMONENE
               AD07(I,J,9)  = AD07(I,J,9)    + BIOG_LIMO(I,J)

               ! TERPENE
               AD07(I,J,10) = AD07(I,J,10)   + BIOG_TERP(I,J)

               ! ALCOHOL
               AD07(I,J,11) = AD07(I,J,11)   + BIOG_ALCO(I,J)

               ! SESQTERPENE
               AD07(I,J,12) = AD07(I,J,12)   + BIOG_SESQ(I,J)

            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARB: after ND07' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE EMISSCARBON

!------------------------------------------------------------------------------

      SUBROUTINE BIOGENIC_OC
!
!******************************************************************************
!  Subroutine BIOGENIC_OC emits secondary organic carbon aerosols.
!  Also modified for SOA tracers. (rjp, bmy, 4/1/04, 1/24/08)
!
!  Terpene emissions as a source of OC:  TERP.GEIA90.a1.2x2.5.*
!  Assuming 10% yield of OC(hydrophilic) from terpene emission.
!
!  NOTES:
!  (1 ) Now separate computation for FULLCHEM and OFFLINE runs (bmy, 7/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f".  Now references LSOA
!        from "logical_mod.f". (bmy, 7/20/04)
!  (3 ) Now reads data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!  (4 ) Now can use MEGAN biogenic emissions (tmf, bmy, 10/20/05)
!  (5 ) For GCAP, need to use GET_NAME_EXT_2D in NVOC file name (bmy, 4/11/06)
!  (6 ) Bug fix: add MEGAN emissions to TERP_ORGC when SOA emissions are
!        turned on (dkh, bmy, 1/24/08)
!  (7 ) Change LMEGAN switch to LMEGANMONO switch (ccc, 3/2/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D,  GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,         READ_BPCH2
      USE DAO_MOD,       ONLY : SUNCOS
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LMEGANMONO,       LSOA
      USE MEGAN_MOD,     ONLY : GET_EMMONOT_MEGAN
      USE TIME_MOD,      ONLY : GET_MONTH,        GET_TS_CHEM
      USE TIME_MOD,      ONLY : GET_TS_EMIS,      ITS_A_NEW_MONTH
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, IJLOOP, THISMONTH
      TYPE (XPLEX)                 :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                 :: CONVERT(NVEGTYPE)
      TYPE (XPLEX)                 :: GMONOT(NVEGTYPE)
      TYPE (XPLEX)                 :: FD2D(IIPAR,JJPAR)
      TYPE (XPLEX)                 :: TMMP, EMMO, VALUE
      TYPE (XPLEX)                 :: XTAU, STEPS_PER_MON
      TYPE (XPLEX), PARAMETER:: FC1 = xplex(136.2364D0 / 120.11D0,0d0)
      TYPE (XPLEX), PARAMETER:: FC2 = xplex(154.2516D0 / 120.11D0,0d0)
      TYPE (XPLEX), PARAMETER:: FC3 = xplex(204.3546D0 / 180.165D0,0d0)
      TYPE (XPLEX), PARAMETER:: FC4 = xplex(152.D0     / 120.11D0,0d0)
      CHARACTER(LEN=255)     :: FILENAME

      ! Fraction of yield of OC (hydrophilic) from terpene emission
      TYPE (XPLEX), PARAMETER      :: FBIOG = xplex(1.0d-1,0d0)

      ! External functions
      TYPE (XPLEX),  EXTERNAL      :: XLTMMP
      TYPE (XPLEX),  EXTERNAL      :: EMMONOT

      !=================================================================
      ! BIOGENIC_OC begins here!
      !=================================================================

      ! Get ISOPRENE baseline emissions (first-time only)
      IF ( FIRST ) THEN
         CALL RDISOPT( CONVERT )
         CALL RDMONOT( GMONOT  )
         CALL SETBASE( CONVERT, GMONOT )

         ! Resety first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! If secondary organic aerosols are turned off ...
      ! Compute biogenic organic carbon as 0.1 * MONOTERPENES
      !=================================================================
      IF ( .not. LSOA ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP, TMMP, EMMO )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! 1-D loop index
            IJLOOP         = ( (J-1) * IIPAR ) + I

             ! Surface temperature [K]
            TMMP           = XLTMMP(I,J,IJLOOP)

            ! Get monoterpenes from MEGAN or GEIA [kg C/box]
            IF ( LMEGANMONO ) THEN
               EMMO = GET_EMMONOT_MEGAN( I, J, TMMP,xplex(1d0,0d0) )
            ELSE
               EMMO = EMMONOT( IJLOOP, TMMP, xplex(1d0,0d0) )
            ENDIF

            ! Fraction of EMMO that converts into OC [kg/box/timestep]
            TERP_ORGC(I,J) = EMMO * FBIOG
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! If secondary organic aerosols are turned on ...
      ! Then use CALTECH algorithm 
      !=================================================================
      ELSE 

         ! Current month
         THISMONTH     = GET_MONTH()

         ! Number of emission timesteps per month
         STEPS_PER_MON = ( ( 1440 * NDAYS(THISMONTH) ) / GET_TS_EMIS() )

         !-----------------------------------------
         ! Read data from disk if it's a new month
         !-----------------------------------------
         IF ( ITS_A_NEW_MONTH() ) THEN
         
            ! Get TAU0 value to index the punch file
            XTAU  = GET_TAU0( THISMONTH, 1, 1990 )

            ! Filename for carbon aerosol from fossil fuel use
            FILENAME = TRIM( DATA_DIR )      //
     &                 'carbon_200411/NVOC.' // GET_NAME_EXT_2D() //
     &                 '.'                   // GET_RES_EXT()

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - BIOGENIC_OC: Reading ', a )

            ! Read NVOC emission in kg/month
            CALL READ_BPCH2( FILENAME, 'NVOCSRCE', 35, 
     &                       XTAU,      IGLOB,     JGLOB,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast to TYPE (XPLEX) and resize
            CALL TRANSFER_2D( ARRAY(:,:,1), GEIA_ORVC )

            ! from kgC/month to kgC/timestep
            GEIA_ORVC(:,:) = GEIA_ORVC(:,:) / STEPS_PER_MON
         ENDIF

         !------------------------------
         ! Get TERP_ORGC and DIUR_ORVC
         !------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP, TMMP, EMMO )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! 1-D loop index
            IJLOOP         = ( (J-1) * IIPAR ) + I

            ! Surface temperature [K]
            TMMP           = XLTMMP(I,J,IJLOOP)

            ! Monoterpene emission [kg C/box/timestep]
            IF ( LMEGANMONO ) THEN
               TERP_ORGC(I,J)=GET_EMMONOT_MEGAN(I,J,TMMP,xplex(1d0,0d0))
            ELSE
               TERP_ORGC(I,J) = EMMONOT( IJLOOP, TMMP,xplex(1d0,0d0) )
            ENDIF

            !---------------------------------------------------
            ! Impose a diurnal variation on NVOC during the day
            !---------------------------------------------------
            IF ( SUNCOS(IJLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN
               DIUR_ORVC(I,J) = GEIA_ORVC(I,J)*
     &                          ( SUNCOS(IJLOOP) / TCOSZ(I,J) ) *
     &                          ( 1440d0        / GET_TS_CHEM() )

               ! Make sure ORVC is not negative
               DIUR_ORVC(I,J) = MAX( DIUR_ORVC(I,J), 0d0 )

            ELSE

               ! At night, ORVC goes to zero
               DIUR_ORVC(I,J) = 0d0

            ENDIF

            !===========================================================
            ! For SOA (Hong Liao, 02/13/04)
            ! 
            ! The input emission data files have units of 
            ! [kg C/box/timestep]
            !
            ! The variable scale is used to convert to the relevant 
            ! units as follows:
            ! (1) Convert from kg C/step to kg compound/step
            ! (2) Multiply by the fraction of monoterpenes that 
            !      contributes to the particular species of interest.
            !
            ! The fraction of monoterpenes is from Table 4 of Griffin
            !  et al., Geophys. Res. Lett. 26 (17): 2721-2724 (1999)
            !===========================================================

            ! ALPHA-PINENE (0.35)
            ! BETA-PINENE lumped with ALPHA-PINENE (0.23)
            ! SABINENE lumped with ALPHA-PINENE (0.05)
            ! D3-CARENE lumped with ALPHA-PINENE (0.04)
            BIOG_ALPH(I,J) = TERP_ORGC(I,J) * FC1 * 0.67D0

            ! TERPENOID KETONE is lumped with SABINENE
            ! Then SABINENE is lumped with ALPHA-PINENE
            BIOG_ALPH(I,J) = BIOG_ALPH(I,J)
     &                     + ( DIUR_ORVC(I,J) * FC4 * 0.04D0 ) !using campher

            ! LIMONENE
            BIOG_LIMO(I,J) = TERP_ORGC(I,J) * FC1 * 0.23D0

            ! TERPINENE is lumped with TERPINOLENE
            BIOG_TERP(I,J) = TERP_ORGC(I,J) * FC1 * 0.03D0

            ! MYRCENE is lumped with TERPENOID ALCOHOL (0.05)
            ! OCIMENE lumped with TERPENOID ALCOHOL (0.02)
            BIOG_ALCO(I,J) = TERP_ORGC(I,J) * FC1 * 0.07D0

            ! Other reactive volatile organic carbon emissions
            BIOG_ALCO(I,J) = BIOG_ALCO(I,J)
     &                     + ( DIUR_ORVC(I,J) * FC2 * 0.09D0 ) !using LINALOOL

            ! We do not transport SESQ (C15H24) 
            ! because its chemical lifetime is short (reactive)
            BIOG_SESQ(I,J) = DIUR_ORVC(I,J) * FC3 * 0.05D0

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
      END SUBROUTINE BIOGENIC_OC

!------------------------------------------------------------------------------
! v8-02-01, now use v8-02-02 version 
!
!      SUBROUTINE ANTHRO_CARB_TBOND
!!
!!******************************************************************************
!!  Subroutine ANTHRO_CARB_TBOND computes annual mean anthropogenic and 
!!  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!!  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!!  fractions. (rjp, bmy, 4/2/04, 5/30/06)
!!
!!  Emissions data comes from the Bond et al [2004] inventory and has units
!!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!!
!!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!!
!!  NOTES:
!!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!!  (3 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!!  (4 ) Now compute future emissions of BC,OC if necessary. (swu, bmy, 5/30/06)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCbf
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCff
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCbf
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCff
!      USE LOGICAL_MOD,          ONLY : LFUTURE
!      USE TIME_MOD,             ONLY : GET_TS_EMIS
!      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"             ! Size parameters
!
!      ! Local variables
!      INTEGER                       :: I, J
!      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)                        :: XTAU, STEPS_PER_YR, FUT_SCL
!      TYPE (XPLEX)                        :: FD2D(IIPAR,JJPAR)
!      CHARACTER(LEN=255)            :: FILENAME
!
!      ! Hydrophilic fraction of BLACK CARBON (aka ELEMENTAL CARBON)
!      TYPE (XPLEX), PARAMETER             :: FHB = 0.2d0
!
!      ! Hydrophilic fraction of ORGANIC CARBON 
!      TYPE (XPLEX), PARAMETER             :: FHO = 0.5d0 
!
!      !=================================================================
!      ! ANTHRO_CARB_TBOND begins here!
!      !=================================================================
!
!      ! Number of emission timesteps per year
!      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )
!
!      ! Get TAU0 value to index the punch file
!      XTAU         = GET_TAU0( 1, 1, 2001 )
!
!      !=================================================================
!      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
!      ! anthropogenic sources as tracer #34 in [kg C/year].  
!      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
!      !=================================================================
!
!      ! Filename for carbon aerosol from fossil fuel use
!      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'carbon_200411/BCOC_TBond_fossil.'       // 
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - ANTHRO_CARB_TBOND: Reading ', a )
!
!      ! Read BLCK emission
!      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
!     &                 XTAU,      IGLOB,     JGLOB,     
!     &                 1,         ARRAY,     QUIET=.TRUE. ) 
!
!      ! Cast to TYPE (XPLEX) and resize
!      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, FUT_SCL )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!         
!         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
!         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR
!         
!         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
!         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR
!        
!         ! Compute future emissions of BLACK CARBON (if necessary)
!         IF ( LFUTURE ) THEN
!            FUT_SCL          = GET_FUTURE_SCALE_BCff( I, J )
!            ANTH_BLKC(I,J,1) = ANTH_BLKC(I,J,1) * FUT_SCL
!            ANTH_BLKC(I,J,2) = ANTH_BLKC(I,J,2) * FUT_SCL
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      !=================================================================
!      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
!      ! in [kg C/year].  Then Convert to [kg C/timestep] and store in 
!      ! ANTH_ORGC.
!      !=================================================================
!      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
!     &                 XTAU,      IGLOB,     JGLOB,     
!     &                 1,         ARRAY,     QUIET=.TRUE. ) 
!
!      ! Cast to TYPE (XPLEX) and resize
!      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, FUT_SCL )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
!         ANTH_ORGC(I,J,1) =          FHO *   FD2D(I,J) / STEPS_PER_YR
!
!         ! Hydrophobic ORGANIC CARBON from anthropogenics [kgC/timestep]
!         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR
!
!         ! Compute future emissions of ORGANIC CARBON (if necessary)
!         IF ( LFUTURE ) THEN
!            FUT_SCL          = GET_FUTURE_SCALE_OCff( I, J )
!            ANTH_ORGC(I,J,1) = ANTH_ORGC(I,J,1) * FUT_SCL
!            ANTH_ORGC(I,J,2) = ANTH_ORGC(I,J,2) * FUT_SCL
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      !=================================================================
!      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
!      ! combustion as tracer #34 in [kg C/year].  Then convert to 
!      ! [kg C/timestep] and store in BIOF_BLKC.
!      !=================================================================
!
!      ! Filename
!      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'carbon_200411/BCOC_TBond_biofuel.'      // 
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
!
!      ! Read data
!      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
!     &                 XTAU,      IGLOB,     JGLOB,
!     &                 1,         ARRAY,     QUIET=.TRUE. ) 
!
!      ! Cast to TYPE (XPLEX) and resize
!      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, FUT_SCL )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Hydrophilic BLACK CARBON from biofuels [kg C /timestep]
!         BIOF_BLKC(I,J,1) =          FHB *   FD2D(I,J) / STEPS_PER_YR
!         
!         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
!         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR
!
!         ! Compute future emissions of BLACK CARBON (if necessary)
!         IF ( LFUTURE ) THEN
!            FUT_SCL          = GET_FUTURE_SCALE_BCbf( I, J )
!            BIOF_BLKC(I,J,1) = BIOF_BLKC(I,J,1) * FUT_SCL
!            BIOF_BLKC(I,J,2) = BIOF_BLKC(I,J,2) * FUT_SCL
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      !=================================================================
!      ! Read ORGANIC CARBON from biofuel combustion as tracer #35 in
!      ! [kg C/year].  Convert to [kg C/timestep] and store in BIOF_BLKC.
!      !=================================================================
!      
!      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
!     &                 XTAU,      IGLOB,     JGLOB,     
!     &                 1,         ARRAY,     QUIET=.TRUE. )
!
!      ! Cast to TYPE (XPLEX) and resize
!      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!         
!         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
!         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR
!
!         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
!         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR
!
!         ! Compute future emissions of BLACK CARBON (if necessary)
!         IF ( LFUTURE ) THEN
!            FUT_SCL          = GET_FUTURE_SCALE_OCbf( I, J )
!            BIOF_ORGC(I,J,1) = BIOF_ORGC(I,J,1) * FUT_SCL
!            BIOF_ORGC(I,J,2) = BIOF_ORGC(I,J,2) * FUT_SCL
!         ENDIF          
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE ANTHRO_CARB_TBOND
!
!!------------------------------------------------------------------------------

      SUBROUTINE ANTHRO_CARB_TBOND( THISMONTH )
!
!******************************************************************************
!  Subroutine ANTHRO_CARB_TBOND reads monthly mean anthropogenic and biofuel
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (eml 4/17/09, rjp, bmy, 4/2/04, 5/30/06)
!
!  Emissions data comes from Bond et al [GBC, 2007] inventory and has units
!  of [kg C/yr], which is converted to [kg C/timestep] below. Seasonality is
!  applied over the US as in Park [2003].
!
!-----------------------------------
!     Prior to 12/18/08
!     !  OLD:
!     !  Subroutine ANTHRO_CARB_TBOND computes annual mean anthropogenic and 
!     !  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!     !  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!     !  fractions. (rjp, bmy, 4/2/04, 5/30/06)
!     !
!     !  Emissions data comes from the Bond et al [2004] inventory and has units
!     !  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!     !
!----------------------------------
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!  (3 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now compute future emissions of BC,OC if necessary. (swu, bmy, 5/30/06)
!  (5 ) Now reads in monthly data from Bond et al [2007] (eml, 4/17/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCbf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCbf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCff
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE TIME_MOD,             ONLY : GET_TS_EMIS
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

#     include "CMN_SIZE"             ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: THISMONTH

      ! Local variables
      INTEGER                       :: I, J
      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                        :: XTAU, STEPS_PER_YR, FUT_SCL
      TYPE (XPLEX)                        :: STEPS_PER_MON
      TYPE (XPLEX)                        :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)            :: FILENAME

      ! Hydrophilic fraction of BLACK CARBON (aka ELEMENTAL CARBON)
      TYPE (XPLEX), PARAMETER             :: FHB = xplex(0.2d0,0d0)

      ! Hydrophilic fraction of ORGANIC CARBON 
      TYPE (XPLEX), PARAMETER             :: FHO = xplex(0.5d0,0d0) 

      !=================================================================
      ! ANTHRO_CARB_TBOND begins here!
      !=================================================================

!-----------------------------------
!     Prior to 12/18/08
!      ! Number of emission timesteps per year
!      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )
!
!      ! Get TAU0 value to index the punch file
!      XTAU         = GET_TAU0( 1, 1, 2001 )
!----------------------------------
      
      ! Number of emission timesteps per month
      STEPS_PER_MON = ( ( 1440 * NDAYS ( THISMONTH ) ) / GET_TS_EMIS() )

      ! Get TAU0 value to index the punch file
      XTAU          = GET_TAU0( THISMONTH, 1, 2000 )   

      
      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
      ! anthropogenic sources as tracer #34 in [kg C/month].  
      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
      !=================================================================
      ! Filename for carbon aerosol from fossil fuel use
      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'carbon_200905/BCOC_TBond_fossil.2000.'  // 
     &           'carbon_200909/BCOC_TBond_fossil.2000.'  // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
      

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - ANTHRO_CARB_TBOND: Reading ', a )

      ! Read BLCK emission
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FUT_SCL )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
!-----------------------------------
! Prior to 12/18/08
!        ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR
         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON
         
         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
!-----------------------------------
! Prior to 12/18/08
!        ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR
         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_BCff( I, J )
            ANTH_BLKC(I,J,1) = ANTH_BLKC(I,J,1) * FUT_SCL
            ANTH_BLKC(I,J,2) = ANTH_BLKC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
      ! in [kg C/month].  Then Convert to [kg C/timestep] and store in 
      ! ANTH_ORGC.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FUT_SCL )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
!-----------------------------------
! Prior to 12/18/08
!        ANTH_ORGC(I,J,1) =          FHO *   FD2D(I,J) / STEPS_PER_YR
         ANTH_ORGC(I,J,1) =          FHO *   FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from anthropogenics [kgC/timestep]
!-----------------------------------
! Prior to 12/18/08
!        ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR
         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of ORGANIC CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_OCff( I, J )
            ANTH_ORGC(I,J,1) = ANTH_ORGC(I,J,1) * FUT_SCL
            ANTH_ORGC(I,J,2) = ANTH_ORGC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
      ! combustion as tracer #34 in [kg C/year].  Then convert to 
      ! [kg C/timestep] and store in BIOF_BLKC.
      !=================================================================
      ! Filename
      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'carbon_200905/BCOC_TBond_biofuel.2000.' // 
     &           'carbon_200909/BCOC_TBond_biofuel.2000.' // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()


      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FUT_SCL )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic BLACK CARBON from biofuels [kg C /timestep]
!-----------------------------------
! Prior to 12/18/08
!        BIOF_BLKC(I,J,1) =          FHB *   FD2D(I,J) / STEPS_PER_YR
         BIOF_BLKC(I,J,1) =          FHB *   FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
!-----------------------------------
! Prior to 12/18/08
!        BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR
         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON


         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_BCbf( I, J )
            BIOF_BLKC(I,J,1) = BIOF_BLKC(I,J,1) * FUT_SCL
            BIOF_BLKC(I,J,2) = BIOF_BLKC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from biofuel combustion as tracer #35 in
      ! [kg C/year].  Convert to [kg C/timestep] and store in BIOF_BLKC.
      !=================================================================
      
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
!-----------------------------------
! Prior to 12/18/08
!        BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR
         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
!-----------------------------------
! Prior to 12/18/08
!        BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR
         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_OCbf( I, J )
            BIOF_ORGC(I,J,1) = BIOF_ORGC(I,J,1) * FUT_SCL
            BIOF_ORGC(I,J,2) = BIOF_ORGC(I,J,2) * FUT_SCL
         ENDIF          
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ANTHRO_CARB_TBOND

!------------------------------------------------------------------------------


      SUBROUTINE ANTHRO_CARB_COOKE( THISMONTH )
!
!******************************************************************************
!  Subroutine ANTHRO_CARB_COOKE computes monthly mean anthropogenic and 
!  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!  fractions. (rjp, bmy, 4/2/04, 5/30/06)
!
!  Emissions data comes from the Cooke et al. [1999] inventory and 
!  seasonality imposed by Park et al. [2003].  The data has units of 
!  [kg C/month].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR.  Now only apply
!        Cooke/RJP emissions over the North American region (i.e. the region
!        bounded by indices I1_NA, J1_NA, I2_NA, J2_NA).  (rjp, bmy, 12/1/04)
!  (3 ) Now can read data from both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now compute future emissions of BC,OC if necessary. (swu, bmy, 5/30/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCbf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCbf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCff
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE TIME_MOD,             ONLY : GET_TS_EMIS
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

#     include "CMN_SIZE"             ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: THISMONTH

      ! Local variables
      INTEGER                       :: I, J
      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                       :: XTAU, STEPS_PER_MON, FUT_SCL
      TYPE (XPLEX)                        :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)            :: FILENAME

      ! Hydrophilic fraction of BLACK CARBON aerosol
      TYPE (XPLEX), PARAMETER             :: FHB = xplex(0.2d0,0d0)

      ! Hydrophilic fraction of ORGANIC CARBON aerosol
      TYPE (XPLEX), PARAMETER             :: FHO = xplex(0.5d0,0d0)

      !=================================================================
      ! ANTHRO_CARB_COOKE begins here!
      !=================================================================

      ! Return if we are running on a nested grid (e.g. China) which
      ! does not cover the North America region (rjp, bmy, 12/1/04)
      IF ( I1_NA + J1_NA + I2_NA + J2_NA == 0 ) RETURN

      ! Number of emission timesteps per month
      STEPS_PER_MON = ( ( 1440 * NDAYS( THISMONTH ) ) / GET_TS_EMIS() )
      
      ! Get TAU0 value to index the punch file
      XTAU = GET_TAU0( THISMONTH, 1, 1998 )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
      ! anthropogenic sources as tracer #34 in [kg C/month].  
      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
      !
      ! The ANTH_BLKC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of ANTH_BLKC over North America below.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )               //
     &           'carbon_200411/BCOC_anthsrce.' // GET_NAME_EXT_2D() // 
     &           '.'                            // GET_RES_EXT()       

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - ANTHRO_CARB_COOKE: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA

         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_BCff( I, J )
            ANTH_BLKC(I,J,1) = ANTH_BLKC(I,J,1) * FUT_SCL
            ANTH_BLKC(I,J,2) = ANTH_BLKC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO

      !=================================================================
      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
      ! in [kg C/month].  Then Convert to [kg C/timestep] and store in 
      ! ANTH_ORGC.
      ! 
      ! The ANTH_ORGC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of ANTH_ORGC over North America below.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA
         
         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,1) = FHO * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of ORGANIC CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_OCff( I, J )
            ANTH_ORGC(I,J,1) = ANTH_ORGC(I,J,1) * FUT_SCL
            ANTH_ORGC(I,J,2) = ANTH_ORGC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
      ! combustion over Canada and the US as tracer #34 in [kg C/year].  
      ! Then convert to [kg C/timestep] and store in BIOF_BLKC.
      !
      ! Seasonality has been imposed using the heating degree approach 
      ! for year 1998 [Park et al., 2003].
      !
      ! The BIOF_BLKC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of BIOF_BLKC over North America below.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )              //
     &           'carbon_200411/BCOC_biofuel.' // GET_NAME_EXT_2D() //
     &           '.'                           // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA
         
         ! Hydrophilic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_BCbf( I, J )
            BIOF_BLKC(I,J,1) = BIOF_BLKC(I,J,1) * FUT_SCL
            BIOF_BLKC(I,J,2) = BIOF_BLKC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO

      !=================================================================
      ! Read ORGANIC CARBON emission from biofuel combustion over 
      ! Canada and the US as tracer #35 in [kg C/year].  Then convert 
      ! to [kg C/timestep] and store in BIOF_ORGC.
      !
      ! Seasonality has been imposed using the heating degree approach 
      ! for year 1998 [Park et al., 2003].
      !
      ! The BIOF_ORGC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of BIOF_ORGC over North America below.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA

         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON

         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_OCbf( I, J )
            BIOF_ORGC(I,J,1) = BIOF_ORGC(I,J,1) * FUT_SCL
            BIOF_ORGC(I,J,2) = BIOF_ORGC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ANTHRO_CARB_COOKE

!------------------------------------------------------------------------------

      SUBROUTINE BIOMASS_CARB_TBOND
!
!******************************************************************************
!  Subroutine BIOMASS_CARB_TBOND computes annual mean biomass burning 
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (rjp, bmy, 4/2/04, 5/30/06)
!
!  Emissions data comes from the Bond et al [2004] inventory and has units
!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!  (3 ) Now can read from both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now compute future emissions of BC,OC if necessary (swu, bmy, 5/30/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCbb
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCbb
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE TIME_MOD,             ONLY : GET_TS_EMIS
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

#     include "CMN_SIZE"             ! Size parameters 

      ! Local variables
      INTEGER                       :: I, J
      TYPE (XPLEX)                        :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                        :: XTAU, STEPS_PER_YR, FUT_SCL
      TYPE (XPLEX)                        :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)            :: FILENAME

      ! Hydrophilic fraction of carbonaceous aerosols
      TYPE (XPLEX), PARAMETER             :: FHB = xplex(0.2d0,0d0)
      TYPE (XPLEX), PARAMETER             :: FHO = xplex(0.5d0,0d0)

      !=================================================================
      ! BIOMASS_CARB_TBOND begins here!
      !=================================================================

      ! Number of emission timesteps per year
      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )

      ! Filename containing biomass emissions
      FILENAME = TRIM( DATA_DIR )                         //
     &           'carbon_200411/BCOC_TBond_biomass.'      // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()


      ! Get TAU0 value to index the punch file
      XTAU = GET_TAU0( 1, 1, 2001 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - BIOMASS_CARB_TBOND: Reading ', a )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from  
      ! biomass burning as tracer #34 in [kg C/year].  Then 
      ! convert to [kg C/timestep] and store in BIOB_BLKC.
      !=================================================================  
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FUT_SCL )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR

         ! Compute future emissions of BLACK CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_BCbb( I, J )
            BIOB_BLKC(I,J,1) = BIOB_BLKC(I,J,1) * FUT_SCL
            BIOB_BLKC(I,J,2) = BIOB_BLKC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from biomass burning as tracer #35 in 
      ! [kg C/year].  Then convert to [kg C/timestep] and store in 
      ! BIOF_BLKC.
      !=================================================================  
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX) and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FUT_SCL )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

         ! Compute future emissions of ORGANIC CARBON (if necessary)
         IF ( LFUTURE ) THEN
            FUT_SCL          = GET_FUTURE_SCALE_OCbb( I, J )
            BIOB_ORGC(I,J,1) = BIOB_ORGC(I,J,1) * FUT_SCL
            BIOB_ORGC(I,J,2) = BIOB_ORGC(I,J,2) * FUT_SCL
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE BIOMASS_CARB_TBOND

!------------------------------------------------------------------------------

      SUBROUTINE BIOMASS_CARB_GEOS
!
!******************************************************************************
!  Subroutine BIOMASS_CARB_GEOS computes monthly mean biomass burning 
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (rjp, bmy, 4/2/04, 2/19/09)
!
!  Emissions are contained in the BIOMASS array of "biomass_mod.f", and will 
!  contain biomass emissions from either the Duncan et al [2001] inventory or 
!  the GFED2 inventory, depending on the option selected at runtime startup.  
!  BIOMASS has units of [atoms C/cm3/s].  Units will be converted to
!  [kg C/timestep] below. 
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f".  Also removed CMN,
!        it's obsolete. (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!  (3 ) Now read BCPO, OCPO biomass burning data directly from files instead
!        of computing from emission factors. (rjp, bmy, 1/11/05)
!  (4 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (5 ) Now compute future emissions of BC,OC if necessary (swu, bmy, 5/30/06)
!  (6 ) Now get biomass emissions from the BIOMASS array of "biomass_mod.f",
!        which will contain either GFED2 or default emissions.  Also move
!        file-reading code to gc_biomass_mod.f. (bmy, 9/25/06)
!  (7 ) Prevent seg fault error when LBIOMASS=F (bmy, 11/3/06)
!  (8 ) Now apply future emissions if necessary (hotp, swu, bmy, 2/19/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,          ONLY : BIOMASS,  IDBBC, IDBOC
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,          ONLY : LBIOMASS, LFUTURE
      USE TIME_MOD,             ONLY : GET_TS_EMIS
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTBCPO,  IDTOCPO
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCbb
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCbb

#     include "CMN_SIZE"             ! Size parameters

      !-------------------
      ! Local variables
      !-------------------

      ! Hydrophilic fraction of BLACK CARBON
      TYPE (XPLEX), PARAMETER     :: FHB = xplex(0.2d0,0d0)

      ! Hydrophilic fraction of ORGANIC CARBON
      TYPE (XPLEX), PARAMETER     :: FHO = xplex(0.5d0,0d0) 
                                    
      INTEGER               :: I,       J
      TYPE (XPLEX)                :: A_CM2,   BIOBC,   BIOOC
      TYPE (XPLEX)                :: CONV_BC, CONV_OC, DTSRCE, FUT_SCL
      CHARACTER(LEN=255)    :: BC_FILE, OC_FILE

      !=================================================================
      ! BIOMASS_CARB_GEOS begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE  = 60d0 * GET_TS_EMIS()

      ! Conversion factor for [s * kg/molec]
      CONV_BC = DTSRCE / XNUMOL(IDTBCPO)
      CONV_OC = DTSRCE / XNUMOL(IDTOCPO)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, BIOBC, BIOOC, A_CM2, FUT_SCL )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box area [cm2]
         A_CM2 = GET_AREA_CM2( J )
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert [molec/cm2/s] --> [kg C/timestep]
            IF ( LBIOMASS ) THEN
               BIOBC         = BIOMASS(I,J,IDBBC) * A_CM2 * CONV_BC
               BIOOC         = BIOMASS(I,J,IDBOC) * A_CM2 * CONV_OC
            ELSE
               BIOBC         = 0d0
               BIOOC         = 0d0
            ENDIF

            ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
            BIOB_BLKC(I,J,1) =          FHB   * BIOBC
         
            ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
            BIOB_BLKC(I,J,2) = ( 1.D0 - FHB ) * BIOBC

            ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
            BIOB_ORGC(I,J,1) =          FHO   * BIOOC

            ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
            BIOB_ORGC(I,J,2) = ( 1.D0 - FHO ) * BIOOC

            ! Apply future emissions for GCAP (if necessary)
            IF ( LFUTURE ) THEN

               ! Compute future emissions of ORGANIC CARBON
               FUT_SCL          = GET_FUTURE_SCALE_OCbb( I, J )
               BIOB_ORGC(I,J,1) = BIOB_ORGC(I,J,1) * FUT_SCL
               BIOB_ORGC(I,J,2) = BIOB_ORGC(I,J,2) * FUT_SCL

               ! Compute future emissions of BLACK CARBON
               FUT_SCL          = GET_FUTURE_SCALE_BCbb( I, J )
               BIOB_BLKC(I,J,1) = BIOB_BLKC(I,J,1) * FUT_SCL
               BIOB_BLKC(I,J,2) = BIOB_BLKC(I,J,2) * FUT_SCL
            ENDIF 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO  

        ! Return to calling program
        END SUBROUTINE BIOMASS_CARB_GEOS

!------------------------------------------------------------------------------

      SUBROUTINE EMITHIGH( BCSRC, OCSRC )
!
!******************************************************************************
!  Subroutine EMITHIGH mixes tracer completely from the surface to the PBL
!  top. (rjp, bmy, 4/2/04, 2/17/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) BCSRC (TYPE (XPLEX)) : Array which holds Total BC (H-phobic & H-philic)
!  (2 ) OCSRC (TYPE (XPLEX)) : Array which holds Total OC (H-phobic & H-philic)
!
!  NOTES:
!  (1 ) Now also mix ALPH, LIMO, ALCO tracers (rjp, bmy, 7/8/04)
!  (2 ) Now reference STT from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Remove references to "dao_mod.f", "pressure_mod.f", and "error_mod.f".
!        Rewrote for computational expediency using routines from
!        "pbl_mix_mod.f".  (bmy, 2/17/05)
!******************************************************************************
!
      ! References to F90 modules
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL,  GET_PBL_MAX_L
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTBCPI, IDTBCPO, IDTOCPI, IDTOCPO
      USE TRACERID_MOD, ONLY : IDTALPH, IDTLIMO, IDTALCO
      USE ERROR_MOD
#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: BCSRC(IIPAR,JJPAR,2)
      TYPE (XPLEX), INTENT(IN) :: OCSRC(IIPAR,JJPAR,2)

      ! Local variables
      LOGICAL            :: IS_BCPO, IS_OCPO, IS_BCPI, IS_OCPI
      LOGICAL            :: IS_ALPH, IS_LIMO, IS_ALCO
      INTEGER            :: I,       J,       L,       PBL_MAX
      TYPE (XPLEX)             :: F_OF_PBL

      !=================================================================
      ! EMITHIGH begins here!
      !=================================================================

      ! Define logical flags for expediency
      IS_BCPI = ( IDTBCPI > 0 )
      IS_OCPI = ( IDTOCPI > 0 ) 
      IS_BCPO = ( IDTBCPO > 0 )
      IS_OCPO = ( IDTOCPO > 0 )
      IS_ALPH = ( IDTALPH > 0 )
      IS_LIMO = ( IDTLIMO > 0 )
      IS_ALCO = ( IDTALCO > 0 )

      ! Maximum extent of PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L()

      !=================================================================
      ! Partition emissions throughout the boundary layer
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_OF_PBL )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
         F_OF_PBL = GET_FRAC_OF_PBL( I, J, L )

         ! Hydrophilic BLACK CARBON
         IF ( IS_BCPI ) THEN
            STT(I,J,L,IDTBCPI) = STT(I,J,L,IDTBCPI) + 
     &                           ( F_OF_PBL * BCSRC(I,J,1) )
         ENDIF

         ! Hydrophilic ORGANIC CARBON
         IF ( IS_OCPI ) THEN
            STT(I,J,L,IDTOCPI) = STT(I,J,L,IDTOCPI) + 
     &                           ( F_OF_PBL * OCSRC(I,J,1) )
         ENDIF
            
         ! Hydrophobic BLACK CARBON
         IF ( IS_BCPO ) THEN
            STT(I,J,L,IDTBCPO) = STT(I,J,L,IDTBCPO) + 
     &                           ( F_OF_PBL * BCSRC(I,J,2) )
         ENDIF

         ! Hydrophobic ORGANIC CARBON
         IF ( IS_OCPO ) THEN
            STT(I,J,L,IDTOCPO) = STT(I,J,L,IDTOCPO) + 
     &                           ( F_OF_PBL * OCSRC(I,J,2) )
         ENDIF

         ! ALPHA-PINENE
         IF ( IS_ALPH ) THEN
            STT(I,J,L,IDTALPH) = STT(I,J,L,IDTALPH) + 
     &                           ( F_OF_PBL * BIOG_ALPH(I,J) )
         ENDIF

         ! LIMONENE
         IF ( IS_LIMO ) THEN
            STT(I,J,L,IDTLIMO) = STT(I,J,L,IDTLIMO) + 
     &                           ( F_OF_PBL * BIOG_LIMO(I,J) )

            ORVC_TERP(I,J,L)   = ORVC_TERP(I,J,L) + 
     &                           ( F_OF_PBL * BIOG_TERP(I,J) )
         ENDIF

         ! ALCOHOL and SESQTERPENE (not a tracer)
         IF ( IS_ALCO ) THEN
            STT(I,J,L,IDTALCO) = STT(I,J,L,IDTALCO) + 
     &                           ( F_OF_PBL * BIOG_ALCO(I,J) )
               
            ORVC_SESQ(I,J,L)   = ORVC_SESQ(I,J,L) + 
     &                           ( F_OF_PBL * BIOG_SESQ(I,J) )
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMITHIGH

!------------------------------------------------------------------------------

      SUBROUTINE OHNO3TIME
!
!******************************************************************************
!  Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 1/18/05)
!
!  NOTES:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (bmy, 1/18/05)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      TYPE (XPLEX)              :: A0, A1, A2, A3, B1, B2, B3
      TYPE (XPLEX)              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      TYPE (XPLEX)              :: SUNTMP(MAXIJ)
      
      !=================================================================
      ! OHNO3TIME begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * dble( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(dble( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TCOSZ(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( J )

            DO I = 1, IIPAR

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC =(LHR0)+dble(NT)/3600.0 + GET_XMID(I)/15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!------------------------------------------------------------------------------

      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_OH returns OH from SMVGEAR's CSPEC array (for coupled runs)
!  or monthly mean OH (for offline runs).  Imposes a diurnal variation on
!  OH for offline simulations. (bmy, 7/9/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP
      USE DAO_MOD,       ONLY : SUNCOS
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,  ONLY : IDOH

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      TYPE (XPLEX)              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !---------------------
         ! Coupled simulation
         !---------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Take OH from the SMVGEAR array CSPEC
         ! OH is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
         ELSE
            OH_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !---------------------
         ! Offline simulation
         !---------------------

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

            ! Impose a diurnal variation on OH during the day
            OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                     ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                     ( 1440d0        / GET_TS_CHEM() )

            ! Make sure OH is not negative
            OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
         ELSE

            ! At night, OH goes to zero
            OH_MOLEC_CM3 = 0d0

         ENDIF

      ELSE

         !---------------------
         ! Invalid sim type!
         !---------------------        
         CALL ERROR_STOP( 'Invalid Simulation Type!', 
     &                    'GET_OH ("carbon_mod.f")' )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_OH

!------------------------------------------------------------------------------

      FUNCTION GET_NO3( I, J, L ) RESULT( NO3_MOLEC_CM3 ) 
!
!******************************************************************************
!  Function GET_NO3 returns NO3 from SMVGEAR's CSPEC array (for coupled runs)
!  or monthly mean OH (for offline runs).  For offline runs, the concentration
!  of NO3 is set to zero during the day. (rjp, bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) Now references ERROR_STOP from "error_mod.f".  We also assume that
!        SETTRACE has been called to define IDNO3.  Now also set NO3 to 
!        zero during the day. (rjp, bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,     ONLY : CSPEC, JLOP
      USE DAO_MOD,        ONLY : AD,    SUNCOS
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GLOBAL_NO3_MOD, ONLY : NO3
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,   ONLY : IDNO3

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      TYPE (XPLEX)              :: NO3_MOLEC_CM3
      TYPE (XPLEX),PARAMETER::XNUMOL_NO3=xplex(6.022d23 / 62d-3,0d0)
 
      ! External functions
      TYPE (XPLEX),  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_NO3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !----------------------
         ! Fullchem simulation
         !----------------------

         ! 1-D SMVGEAR grid box index
         JLOOP = JLOP(I,J,L)

         ! Take NO3 from the SMVGEAR array CSPEC
         ! NO3 is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            NO3_MOLEC_CM3 = CSPEC(JLOOP,IDNO3)
         ELSE
            NO3_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !==============================================================  
         ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
         ! in [v/v].  Convert these to [molec/cm3] as follows:
         !
         !  vol NO3   moles NO3    kg air     kg NO3/mole NO3
         !  ------- = --------- * -------- * ---------------- =  kg NO3 
         !  vol air   moles air      1        kg air/mole air
         !
         ! And then we convert [kg NO3] to [molec NO3/cm3] by:
         !  
         !  kg NO3   molec NO3   mole NO3     1     molec NO3
         !  ------ * --------- * -------- * ----- = --------- 
         !     1     mole NO3     kg NO3     cm3       cm3
         !          ^                    ^
         !          |____________________|  
         !            this is XNUMOL_NO3
         !
         ! If at nighttime, use the monthly mean NO3 concentration from
         ! the NO3 array of "global_no3_mod.f".  If during the daytime,
         ! set the NO3 concentration to zero.  We don't have to relax to 
         ! the monthly mean concentration every 3 hours (as for HNO3) 
         ! since NO3 has a very short lifetime. (rjp, bmy, 12/16/02) 
         !==============================================================

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test if daylight
         IF ( SUNCOS(JLOOP) > 0d0 ) THEN

            ! NO3 goes to zero during the day
            NO3_MOLEC_CM3 = 0d0
              
         ELSE

            ! At night: Get NO3 [v/v] and convert it to [kg]
            NO3_MOLEC_CM3 = NO3(I,J,L) * AD(I,J,L) * ( 62d0/28.97d0 ) 
               
            ! Convert NO3 from [kg] to [molec/cm3]
            NO3_MOLEC_CM3 = NO3_MOLEC_CM3 * XNUMOL_NO3 / BOXVL(I,J,L)
                  
         ENDIF
            
         ! Make sure NO3 is not negative
         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0d0 )

      ELSE

         !----------------------
         ! Invalid sim type!
         !----------------------       
         CALL ERROR_STOP( 'Invalid Simulation Type!',
     &                    'GET_NO3 ("carbon_mod.f")' )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_NO3

!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Now reference XNUMOLAIR from "tracer_mod.f" (bmy, 10/20/05)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP, VOLUME
      USE DAO_MOD,       ONLY : SUNCOS, AD
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GLOBAL_O3_MOD, ONLY : O3
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,    ONLY : XNUMOLAIR
      USE TRACERID_MOD,  ONLY : IDO3

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      TYPE (XPLEX)              :: O3_MOLEC_CM3
      TYPE (XPLEX),  PARAMETER::XNUMOL_O3 =xplex(6.022d23 / 48d-3,0d0)
 
      ! External functions
      TYPE (XPLEX),  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_O3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Get O3 from CSPEC [molec/cm3]
         ! O3 data will only be defined below the tropopause
         IF ( JLOOP  > 0 ) THEN
            O3_MOLEC_CM3 = CSPEC(JLOOP,IDO3)
         ELSE
            O3_MOLEC_CM3 = 0d0
         ENDIF
         
      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !--------------------
         ! Offline simulation
         !--------------------

         ! Get O3 [v/v] for this gridbox & month
         ! O3 data will only be defined below the tropopause
         IF ( L <= LLTROP ) THEN

            ! Get O3 [v/v] and convert it to [kg]
            O3_MOLEC_CM3 = O3(I,J,L) * AD(I,J,L) * ( 48d0/28.97d0 )
               
            ! Convert O3 from [kg] to [molec/cm3]
            O3_MOLEC_CM3 = O3_MOLEC_CM3 * XNUMOL_O3 / BOXVL(I,J,L)
         ELSE
            O3_MOLEC_CM3 = 0d0
         ENDIF

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

            ! Impose a diurnal variation on OH during the day
            O3_MOLEC_CM3 = O3_MOLEC_CM3                     *        
     &                     ( SUNCOS(JLOOP) / TCOSZ(I,J) )   *
     &                     ( 1440d0        / GET_TS_CHEM() )

            ! Make sure OH is not negative
            O3_MOLEC_CM3 = MAX( O3_MOLEC_CM3, 0d0 )

         ELSE
            O3_MOLEC_CM3 = 0d0
         ENDIF

      ELSE

         !--------------------
         ! Invalid sim type!
         !--------------------
         CALL ERROR_STOP( 'Invalid Simulation Type!', 
     &                     'GET_O3 ("carbon_mod.f")' )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------

      FUNCTION GET_DOH( I, J, L ) RESULT( DOH )
!
!******************************************************************************
!  Function GET_DOH returns the amount of isoprene [kg] that has reacted with 
!  OH during the last chemistry time step. (dkh, bmy, 6/01/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP, VOLUME
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,    ONLY : XNUMOL,             TRACER_COEFF
      USE TRACERID_MOD,  ONLY : IDTISOP

#     include "CMN_SIZE"      ! Size parameters
#     include "comode.h"      ! ILISOPOH   (dkh, 05/29/06)  

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J, L

      ! Function returns 
      TYPE (XPLEX)                 :: DOH

      ! Local variables
      INTEGER                :: JLOOP

      !=================================================================
      ! GET_DOH begins here!
      !=================================================================

      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------

         ! Get 1-D index from current 3-D position
         JLOOP = JLOP(I,J,L)

         ! Test if we are in the troposphere
         IF ( JLOOP > 0 ) THEN 
 
            !-----------------------------------------------------------
            ! Get DOH from CSPEC [molec/cm3] and 
            ! convert to [kg C isop / box]
            ! 
            ! CSPEC(JLOOP,ILISOPOH)    : molec isop (lost to OH) / cm3
            !  XNUMOL(IDTISOP)         : atom C / kg C isop
            !  TRACER_COEFF(IDTISOP,1) : atom C / molec isop
            !  VOLUME                  : cm3 / box
            !-----------------------------------------------------------
            DOH = CSPEC(JLOOP,ILISOPOH)   * 
     &            VOLUME(JLOOP)           * 
     &            TRACER_COEFF(IDTISOP,1) / 
     &            XNUMOL(IDTISOP) 
 
         ELSE

            ! Otherwise set DOH=0
            DOH = 0d0
 
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !--------------------
         ! Offline simulation
         !--------------------   

         ! ISOP from OH not is yet supported for
         ! offline aerosol simulations, set DOH=0
         DOH = 0d0

      ELSE

         !--------------------
         ! Invalid sim type!
         !--------------------
         CALL ERROR_STOP( 'Invalid simulation type!', 
     &                    'GET_DOH ("carbon_mod.f")' )

      ENDIF 

      ! Return to calling program
      END FUNCTION GET_DOH

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
!  (1 ) Copied from 'sulfate_mod.f' for cloud uptake of GLYX and MGLY (tmf, 2/26/07)
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
!  (1 ) Copied from 'sulfate_mod.f' for cloud uptake of GLYX and MGLY (tmf, 2/26/07)
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

      ! Return to calling program
      END FUNCTION GET_LWC

!------------------------------------------------------------------------------

      SUBROUTINE SOAG_CLOUD
!
!******************************************************************************
!  Subroutine SOAG_CLOUD produces SOAG from GLYX during a cloud event.
!  Mimics the SO2 -> SO4 process from 'sulfate_mod.f'.  (tmf, 2/26/07)
!
!  Procedure:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  (1 ) SOAG (SOA product of GLYX is produced at existing hydrophilic aerosol
!        surface. (tmf, 2/26/07)
!  (2 ) Assume marine and continental cloud droplet size (tmf, 2/26/07)
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DAO_MOD,         ONLY : AD, T, AIRVOL
      USE DAO_MOD,         ONLY : IS_LAND        ! return true if sfc grid box is land
      USE DIAG_MOD,        ONLY : AD07_SOAGM
      USE TIME_MOD,        ONLY : GET_TS_CHEM
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT
      USE TRACER_MOD,      ONLY : STT
      USE TRACERID_MOD,    ONLY : IDTGLYX, IDTSOAG

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Local variables
      INTEGER   :: I, J, L
      TYPE (XPLEX)    :: DTCHEM      ! Chemistry time step [s]
      TYPE (XPLEX)    :: XAIRM       ! Air mass in grid box [kg/box]
      TYPE (XPLEX)    :: XAIRM3      ! Air volume in grid box [m3]
      TYPE (XPLEX)    :: XGASM       ! Gas mass at grid box before uptake [kg]
      TYPE (XPLEX)    :: XGASC       ! Gas concentration at grid box before uptake [molec/cm3]
      TYPE (XPLEX)    :: XGASMIX     ! Gas mixing ratio [v/v]

      TYPE (XPLEX)    :: XCLDR    ! cloud droplet radius [cm]
      TYPE (XPLEX)    :: XDF      ! gas-phase diffusivity [cm2/s]
      TYPE (XPLEX)    :: XMS      ! Mean molecular speed [cm/s] 
      TYPE (XPLEX)    :: XKT      ! phase-transfer coefficient [1/s]
      TYPE (XPLEX)    :: XTEMP    ! Temperature [K]
      TYPE (XPLEX)    :: XDELTAC  ! Potential maximum change of gas concentration due to cloud chemistry [molecules/cm3]
      TYPE (XPLEX)    :: XUPTKMAX    ! Potential maximum uptake of gas by cloud in grid box [kg]
      TYPE (XPLEX)    :: XUPTK       ! Actual uptake of gas by cloud in grid box [kg]
                                          !  XUPTK <= STT( I, J, L, IDTGLYX )
      TYPE (XPLEX)    :: FC          ! Cloud fraction by volume [unitless]
      TYPE (XPLEX)    :: LWC         ! Liquid water content [g/m3]

      ! Parameters
      TYPE (XPLEX), PARAMETER :: XCLDR_CONT =  xplex(6.d-4,0d0)  ! Cloud droplet radius in continental warm clouds [cm]
      TYPE (XPLEX), PARAMETER :: XCLDR_MARI = xplex(10.d-4,0d0)  ! Cloud droplet radius in marine warm clouds [cm]
      TYPE (XPLEX), PARAMETER :: XMW = xplex(58.d0,0d0)    ! Molecular weight of glyoxal [g/mole]
      TYPE (XPLEX), PARAMETER :: XNAVO = xplex(6.023d23,0d0)    ! Avogadro's number
      TYPE (XPLEX), PARAMETER :: MINDAT = xplex(1.d-20,0d0)   ! Minimum GLYX mixing ratio to calculate cloud uptake
      TYPE (XPLEX), PARAMETER :: XGAMMA = xplex(2.9d-3,0d0)   ! Uptake coefficient (Assume XGAMMA = 2.9d-3 following Liggio et al., 2005)

      !=================================================================
      ! SOAG_CLOUD
      !=================================================================

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over tropospheric grid boxes
      DO L = 1, LLTROP  
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! initialize for safety
         XUPTKMAX = 0d0
         XUPTK    = 0d0

         ! Get temperature
         XTEMP   = T( I, J, L )

         ! Get air mass  [kg/box]
         XAIRM   = AD( I, J, L  )

         ! Get air volumne [m3]
         XAIRM3  = AIRVOL( I, J, L )

         ! Get gas mass at grid box [kg]
         XGASM   = STT( I, J, L, IDTGLYX )

         ! Get gas concentration at grid box [molec/cm3]
         XGASC   = XGASM / (XMW*1.d-3) * XNAVO / (XAIRM3*1.d6) 

         ! GET gas mixing ratio [v/v]
         XGASMIX = XGASM / XMW / ( XAIRM / 28.97d0 )

         ! Volume cloud fraction (Sundqvist et al 1989) [unitless]
         FC      = VCLDF(I,J,L)

         ! Liquid water content in cloudy area of grid box [g/m3]
         LWC     = GET_LWC( XTEMP ) * FC

         !==============================================================
         ! If (1) there is cloud, (2) there is GLYX present, and 
         ! (3) the T > -15 C, then compute cloud uptake
         !==============================================================
         IF ( ( FC     > 0.d0   )  .AND. 
     &        ( XGASMIX > MINDAT )  .AND. 
     &        ( XTEMP   > 258.0  ) ) THEN

            IF ( IS_LAND(I,J) ) THEN
               XCLDR = XCLDR_CONT     ! Continental cloud droplet radius  [m]
            ELSE
               XCLDR = XCLDR_MARI     ! Marine cloud droplet radius  [m]
            ENDIF

            !---------------------------------------
            ! Gas phase diffusivity [cm2/s]       [Lim et al., 2005 Eq. (4)]
            !---------------------------------------
            XDF = 1.9d0 * (XMW**(-0.667))

            !---------------------------------------
            ! Mean molecular speed [cm/s]         [Lim et al., 2005 Eq. (5)]
            !  XMS = SQRT( ( 8 * Boltzmann const * Temperature * N_Avogadro  ) / 
            !              ( pi * molecular weight [g/mole] ) )
            !      = SQRT( 2.117d8 * Temperature / molecular weight )
            !---------------------------------------
            XMS = SQRT( 2.117d8 * XTEMP / XMW )

            !---------------------------------------
            ! Phase transfer coeff [1/s]          [Lim et al., 2005 Eq. (3)] 
            ! XGAMMA = ALPHA, XGAMMA = 2.9d-3 following Liggio et al., 2005
            !---------------------------------------
            XKT = 1.d0 / ( ( XCLDR * XCLDR / 3.d0 / XDF ) + 
     &                     ( 4.d0 * XCLDR / 3.d0 / XMS / XGAMMA ) )

            !---------------------------------------
            ! Maximum potential change in concentration [molecules/cm3]   [Lim et al., 2005 Eq. (1)]
            !---------------------------------------
            XDELTAC = LWC * XKT * XGASC * DTCHEM

            !---------------------------------------
            ! Maximum potential uptake of gas mass [kg/box]
            !---------------------------------------
            XUPTKMAX = XDELTAC * 1.d6 / XNAVO * XMW * 1.d-3 * XAIRM3

            !---------------------------------------
            ! However, the mass of gas being absorbed by aerosol 
            !  cannot exceed the original amount of gas XGASM
            !---------------------------------------
            XUPTK  = MIN( XUPTKMAX, XGASM )
            
            ! Update GLYX in the STT array
            STT( I, J, L, IDTGLYX ) = STT( I, J, L, IDTGLYX ) -
     &                                XUPTK

            ! Update SOAG in the STT array
            STT( I, J, L, IDTSOAG ) = STT( I, J, L, IDTSOAG ) + 
     &                                XUPTK

            !==============================================================
            ! ND07 diagnostic: SOAG from GLYX in cloud [kg/timestep]
            !==============================================================
            IF ( ND07 > 0 .and. L <= LD07 ) THEN
               AD07_SOAGM(I,J,L,3) = AD07_SOAGM(I,J,L,3) + XUPTK
            ENDIF


         ENDIF    ! End of IN CLOUD criteria

      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program 
      END SUBROUTINE SOAG_CLOUD

!------------------------------------------------------------------------------

      SUBROUTINE SOAM_CLOUD
!
!******************************************************************************
!  Subroutine SOAM_CLOUD produces SOAM from MGLY during a cloud event.
!  Mimics the SO2 -> SO4 process from 'sulfate_mod.f'.  (tmf, 2/26/07)
!
!  Procedure:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  (1 ) SOAM (SOA product of MGLY is produced at existing hydrophilic aerosol
!        surface. (tmf, 2/26/07)
!  (2 ) Assume typical marine and continental cloud droplet size (tmf, 2/26/07)
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DAO_MOD,         ONLY : AD, T, AIRVOL
      USE DAO_MOD,         ONLY : IS_LAND        ! return true if sfc grid box is land
      USE DIAG_MOD,        ONLY : AD07_SOAGM
      USE TIME_MOD,        ONLY : GET_TS_CHEM
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT
      USE TRACER_MOD,      ONLY : STT
      USE TRACERID_MOD,    ONLY : IDTMGLY, IDTSOAM

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Local variables
      INTEGER   :: I, J, L
      TYPE (XPLEX)    :: DTCHEM      ! Chemistry time step [s]
      TYPE (XPLEX)    :: XAIRM       ! Air mass in grid box [kg/box]
      TYPE (XPLEX)    :: XAIRM3      ! Air volume in grid box [m3]
      TYPE (XPLEX)    :: XGASM       ! Gas mass at grid box before uptake [kg]
      TYPE (XPLEX)    :: XGASC       ! Gas concentration at grid box before uptake [molec/cm3]
      TYPE (XPLEX)    :: XGASMIX     ! Gas mixing ratio [v/v]

      TYPE (XPLEX)    :: XCLDR       ! cloud droplet radius [cm]
      TYPE (XPLEX)    :: XDF         ! gas-phase diffusivity [cm2/s]
      TYPE (XPLEX)    :: XMS         ! Mean molecular speed [cm/s] 
      TYPE (XPLEX)    :: XKT         ! phase-transfer coefficient [1/s]
      TYPE (XPLEX)    :: XTEMP       ! Temperature [K]
      TYPE (XPLEX)    :: XDELTAC     ! Potential maximum change of gas concentration due to cloud chemistry [molecules/cm3]
      TYPE (XPLEX)    :: XUPTKMAX    ! Potential maximum uptake of gas by cloud in grid box [kg]
      TYPE (XPLEX)    :: XUPTK       ! Actual uptake of gas by cloud in grid box [kg]
                                          !  XUPTK <= STT( I, J, L, IDTGLYX )
      TYPE (XPLEX)    :: FC          ! Cloud fraction by volume [unitless]
      TYPE (XPLEX)    :: LWC         ! Liquid water content [g/m3]

      ! Parameters
      TYPE (XPLEX), PARAMETER :: XCLDR_CONT =  xplex(6.d-4,0d0)  ! Cloud droplet radius in continental warm clouds [cm]
      TYPE (XPLEX), PARAMETER :: XCLDR_MARI = xplex(10.d-4,0d0)  ! Cloud droplet radius in marine warm clouds [cm]
      TYPE (XPLEX), PARAMETER :: XMW = xplex(72.d0,0d0)    ! Molecular weight of methylglyoxal [g/mole]
      TYPE (XPLEX), PARAMETER :: XNAVO = xplex(6.023d23,0d0)    ! Avogadro's number
      TYPE (XPLEX), PARAMETER :: MINDAT = xplex(1.d-20,0d0)   ! Minimum GLYX mixing ratio to calculate cloud uptake
      TYPE (XPLEX), PARAMETER :: XGAMMA = xplex(2.9d-3,0d0)   ! Uptake coefficient (Assume XGAMMA = 2.9d-3 following Liggio et al., 2005)

      !=================================================================
      ! SOAG_CLOUD
      !=================================================================

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over tropospheric grid boxes
      DO L = 1, LLTROP  
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! initialize for safety
         XUPTKMAX = 0d0
         XUPTK    = 0d0

         ! Get temperature
         XTEMP   = T( I, J, L )

         ! Get air mass  [kg/box]
         XAIRM   = AD( I, J, L  )

         ! Get air volumne [m3]
         XAIRM3  = AIRVOL( I, J, L )

         ! Get gas mass at grid box [kg]
         XGASM   = STT( I, J, L, IDTMGLY )

         ! Get gas concentration at grid box [molec/cm3]
         XGASC   = XGASM / (XMW*1.d-3) * XNAVO / (XAIRM3*1.d6) 

         ! GET gas mixing ratio [v/v]
         XGASMIX = XGASM / XMW / ( XAIRM / 28.97d0 )

         ! Volume cloud fraction (Sundqvist et al 1989) [unitless]
         FC      = VCLDF(I,J,L)

         ! Liquid water content in cloudy area of grid box [g/m3]
         LWC     = GET_LWC( XTEMP ) * FC

         !==============================================================
         ! If (1) there is cloud, (2) there is MGLY present, and 
         ! (3) the T > -15 C, then compute cloud uptake
         !==============================================================
         IF ( ( FC     > 0.d0   )  .AND. 
     &        ( XGASMIX > MINDAT )  .AND. 
     &        ( XTEMP   > 258.0  ) ) THEN

            IF ( IS_LAND(I,J) ) THEN
               XCLDR = XCLDR_CONT     ! Continental cloud droplet radius  [m]
            ELSE
               XCLDR = XCLDR_MARI     ! Marine cloud droplet radius  [m]
            ENDIF

            !---------------------------------------
            ! Gas phase diffusivity [cm2/s]       [Lim et al., 2005 Eq. (4)]
            !---------------------------------------
            XDF = 1.9d0 * (XMW**(-0.667))

            !---------------------------------------
            ! Mean molecular speed [cm/s]         [Lim et al., 2005 Eq. (5)]
            !  XMS = SQRT( ( 8 * Boltzmann const * Temperature * N_Avogadro  ) / 
            !              ( pi * molecular weight [g/mole] ) )
            !      = SQRT( 2.117d8 * Temperature / molecular weight )
            !---------------------------------------
            XMS = SQRT( 2.117d8 * XTEMP / XMW )

            !---------------------------------------
            ! Phase transfer coeff [1/s]          [Lim et al., 2005 Eq. (3)] 
            ! XGAMMA = ALPHA, XGAMMA = 2.9d-3 following Liggio et al., 2005
            !---------------------------------------
            XKT = 1.d0 / ( ( XCLDR * XCLDR / 3.d0 / XDF ) + 
     &                     ( 4.d0 * XCLDR / 3.d0 / XMS / XGAMMA ) )

            !---------------------------------------
            ! Maximum potential change in concentration [molecules/cm3]   [Lim et al., 2005 Eq. (1)]
            !---------------------------------------
            XDELTAC = LWC * XKT * XGASC * DTCHEM

            !---------------------------------------
            ! Maximum potential uptake of gas mass [kg/box]
            !---------------------------------------
            XUPTKMAX = XDELTAC * 1.d6 / XNAVO * XMW * 1.d-3 * XAIRM3

            !---------------------------------------
            ! However, the mass of gas being absorbed by aerosol 
            !  cannot exceed the original amount of gas XGASM
            !---------------------------------------
            XUPTK  = MIN( XUPTKMAX, XGASM )
            
            ! Update MGLY in the STT array
            STT( I, J, L, IDTMGLY ) = STT( I, J, L, IDTMGLY ) -
     &                                XUPTK

            ! Update SOAM in the STT array
            STT( I, J, L, IDTSOAM ) = STT( I, J, L, IDTSOAM ) + 
     &                                XUPTK

            !==============================================================
            ! ND07 diagnostic: SOAM from MGLY in cloud [kg/timestep]
            !==============================================================
            IF ( ND07 > 0 .and. L <= LD07 ) THEN
               AD07_SOAGM(I,J,L,4) = AD07_SOAGM(I,J,L,4) + XUPTK
            ENDIF


         ENDIF    ! End of IN CLOUD criteria

      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program 
      END SUBROUTINE SOAM_CLOUD

! <<<<<
!------------------------------------------------------------------------------


      SUBROUTINE WRITE_GPROD_APROD( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine WRITE_GPROD_APROD writes the SOA quantities GPROD and APROD to 
!  disk at the start of a new diagnostic interval. (tmf, havala, bmy, 2/6/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : YYYY/MM/DD value at which to write file
!  (2 ) HHMMSS   (INTEGER) : hh:mm:ss value at which to write file
!  (3 ) TAU      (TYPE (XPLEX) ) : TAU value corresponding to YYYYMMDD, HHMMSS
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,  ONLY : BPCH2,         GET_MODELNAME
      USE BPCH2_MOD,  ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE FILE_MOD,   ONLY : IU_FILE
      USE GRID_MOD,   ONLY : GET_XOFFSET,   GET_YOFFSET
      USE TIME_MOD,   ONLY : EXPAND_DATE

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
      TYPE (XPLEX),  INTENT(IN) :: TAU

      ! Local variables
      INTEGER             :: HALFPOLAR
      INTEGER, PARAMETER  :: CENTER180 = 1
      INTEGER             :: I0, IPR, J0, JHC, N
      TYPE (XPLEX)              :: LONRES,   LATRES
      TYPE (XPLEX)              :: ARRAY(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED = ''
      CHARACTER(LEN=80)   :: TITLE 
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! WRITE_GPROD_APROD begins here
      !=================================================================
      
      ! Define variables for binary punch file
      FILENAME  = 'restart_gprod_aprod.YYYYMMDDhh'
      TITLE     = 'GEOS-Chem SOA restart file: GPROD & APROD'
      UNIT      = 'kg/kg'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. )
      
      ! Replace date & time tokens in FILENAME
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )

      ! Loop over VOC classes and products
      DO JHC = 1, MHC 
      DO IPR = 1, NPROD

         ! Tracer number
         N = ( JHC - 1 ) * NPROD + IPR

         !----------------------
         ! Write GPROD to file
         !----------------------
         
         ! Initialize
         CATEGORY = 'IJ-GPROD'
         ARRAY    = GPROD(:,:,:,IPR,JHC)

         ! Write to disk
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      TAU,       TAU,      RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,    I0+1,            
     &               J0+1,      1,         ARRAY )

         !----------------------
         ! Write APROD to file 
         !----------------------

         ! Initialize
         CATEGORY = 'IJ-APROD'
         ARRAY    = APROD(:,:,:,IPR,JHC)

         ! Write to disk
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      TAU,       TAU,      RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,    I0+1,            
     &               J0+1,      1,         ARRAY )

      ENDDO
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE WRITE_GPROD_APROD

!------------------------------------------------------------------------------

      SUBROUTINE READ_GPROD_APROD( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine READ_GPROD_APROD writes the SOA quantities GPROD and APROD from 
!  disk at the start of a new diagnostic interval. (tmf, havala, bmy, 2/6/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : YYYY/MM/DD value at which to write file
!  (2 ) HHMMSS   (INTEGER) : hh:mm:ss value at which to write file
!  (3 ) TAU      (TYPE (XPLEX) ) : TAU value corresponding to YYYYMMDD, HHMMSS
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,  ONLY : READ_BPCH2
      USE FILE_MOD,   ONLY : IU_FILE
      USE TIME_MOD,   ONLY : EXPAND_DATE

#     include "CMN_SIZE"   ! Size parameters      

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
      TYPE (XPLEX),  INTENT(IN) :: TAU

      ! Local variables
      INTEGER             :: IPR, JHC, N
      TYPE (XPLEX)              :: ARRAY(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_GPROD_APROD begins here!
      !=================================================================

      ! File name
      FILENAME = 'restart_gprod_aprod.YYYYMMDDhh'

      ! Replace date & time tokens in FILENAME
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Loop over VOC classes and products
      DO JHC = 1, MHC 
      DO IPR = 1, NPROD

         ! Tracer number
         N = ( JHC - 1 ) * NPROD + IPR

         !-----------------------
         ! Read GPROD from file
         !-----------------------

         ! Read data
         CALL READ_BPCH2( FILENAME, 'IJ-GPROD', N, 
     &                    TAU,       IIPAR,     JJPAR,     
     &                    LLPAR,     ARRAY,     QUIET=.TRUE. ) 

         ! Cast to TYPE (XPLEX)
         GPROD(:,:,:,IPR,JHC) = ARRAY

         !-----------------------
         ! Read APROD from file
         !-----------------------

         ! Read data
         CALL READ_BPCH2( FILENAME, 'IJ-APROD', N, 
     &                    TAU,       IIPAR,     JJPAR,     
     &                    LLPAR,     ARRAY,     QUIET=.TRUE. ) 

         ! Cast to TYPE (XPLEX)
         APROD(:,:,:,IPR,JHC) = ARRAY

      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_GPROD_APROD

!------------------------------------------------------------------------------

      SUBROUTINE INIT_CARBON
!
!******************************************************************************
!  Subroutine INIT_CARBON initializes all module arrays. 
!  (rjp, bmy, 4/1/04, 11/6/08)
!
!  NOTES:
!  (1 ) Also added arrays for secondary organic aerosols (rjp, bmy, 7/8/04)
!  (2 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!  (3 ) Now reference LSOA from "logical_mod.f" not CMN_SETUP.  Now call
!        GET_BOUNDING_BOX from "grid_mod.f" to compute the indices I1_NA,
!        I2_NA, J1_NA, J2_NA which define the N. America region. (bmy, 12/1/04)
!  (4 ) Now call READ_GPROD_APROD to read GPROD & APROD from disk. 
!        (tmf, havala, bmy, 2/6/07)
!  (5 ) Now set I1_NA, I2_NA, J1_NA, J2_NA appropriately for both 1 x 1 and
!        0.5 x 0.666 nested grids (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR, ERROR_STOP
      USE GRID_MOD,    ONLY : GET_BOUNDING_BOX
      USE LOGICAL_MOD, ONLY : LCHEM,     LSOA 
      USE TIME_MOD,    ONLY : GET_NYMDb, GET_NHMSb, GET_TAUb

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS, INDICES(4), YYYYMMDD, HHMMSS
      TYPE (XPLEX)        :: COORDS(4), TAU

      !=================================================================
      ! INIT_CARBON begins here!
      !=================================================================
      
      ! Return if we already allocated arrays
      IF ( IS_INIT ) RETURN

      ALLOCATE( ANTH_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ANTH_BLKC' )
      ANTH_BLKC = 0d0

      ALLOCATE( ANTH_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ANTH_ORGC' )
      ANTH_ORGC = 0d0

      ALLOCATE( BIOB_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOB_BLKC' )
      BIOB_BLKC = 0d0

      ALLOCATE( BIOB_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOB_ORGC' )
      BIOB_ORGC = 0d0

      ALLOCATE( BIOF_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOF_BLKC' )
      BIOF_BLKC = 0d0

      ALLOCATE( BIOF_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOF_ORGC' )
      BIOF_ORGC = 0d0

      ALLOCATE( TERP_ORGC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TERP_ORGC' )
      TERP_ORGC = 0d0

      ALLOCATE( BCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONV' )
      BCCONV = 0d0

      ALLOCATE( OCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV' )
      OCCONV = 0d0

      !=================================================================
      ! These only have to be allocated if we are
      ! reading in monthly/8-day/3-hr mean
      !=================================================================
      IF ( .not. USE_BOND_BIOBURN ) THEN

         ALLOCATE( EF_BLKC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EF_BLKC' )
         EF_BLKC = 0d0

         ALLOCATE( EF_ORGC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EF_ORGC' )
         EF_ORGC = 0d0

      ENDIF

      !=================================================================
      ! SOA arrays only have to be allocated if LSOA = T
      !=================================================================
      IF ( LSOA ) THEN

         ALLOCATE( BIOG_ALPH( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_ALPH' )
         BIOG_ALPH = 0d0

         ALLOCATE( BIOG_LIMO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_LIMO' )
         BIOG_LIMO = 0d0

         ALLOCATE( BIOG_ALCO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_ALCO' )
         BIOG_ALCO = 0d0

         ALLOCATE( BIOG_TERP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_TERP' )
         BIOG_TERP = 0d0

         ALLOCATE( BIOG_SESQ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_SESQ' )
         BIOG_SESQ = 0d0

         ALLOCATE( DIUR_ORVC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIUR_ORVC' )
         DIUR_ORVC = 0d0

         ALLOCATE( GEIA_ORVC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEIA_ORVC' )
         GEIA_ORVC = 0d0

         ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
         TCOSZ = 0d0

         ALLOCATE( ORVC_TERP( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ORVC_TERP' )
         ORVC_TERP = 0d0

         ALLOCATE( ORVC_SESQ( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ORVC_SESQ' )
         ORVC_SESQ = 0d0

         ALLOCATE( VCLDF( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'VCLDF' )
         VCLDF = 0d0

         ALLOCATE( GPROD( IIPAR, JJPAR, LLPAR, NPROD, MHC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GPROD' )
         GPROD = 0D0

         ALLOCATE( APROD( IIPAR, JJPAR, LLPAR, NPROD, MHC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'APROD' )
         APROD = 0D0

         !--------------------------------------------------------------
         ! Read GPROD and APROD from disk from the last simulation
         ! NOTE: do this here after GPROD, APROD are allocated!
         !--------------------------------------------------------------

         IF ( LCHEM ) THEN 

            ! Get time values at start of the run
            YYYYMMDD = GET_NYMDb()
            HHMMSS   = GET_NHMSb()
            TAU      = GET_TAUb()

            ! Read GPROD, APROD from restart file
            CALL READ_GPROD_APROD( YYYYMMDD, HHMMSS, TAU )
         ENDIF
      ENDIF

      !=================================================================
      ! Compute indices which define the N. America region so that we 
      ! can overwrite T. Bond emissions w/ Cooke/RJP emissions 
      !=================================================================

#if   defined( NESTED_NA )
      ! For 1x1 N. America nested grid: set indices to grid extent
      I1_NA = 1
      J1_NA = 1
      I2_NA = IIPAR
      J2_NA = JJPAR

#elif defined( NESTED_CH )
      ! For 1x1 China nested grid: we don't cover N. America region
      ! Setting these to zero will turn off Cooke/RJP emissions
      I1_NA = 0
      J1_NA = 0
      I2_NA = 0
      J2_NA = 0

#else

      ! Definition of the N. American bounding box
      ! with LL corner (10N,165W) and UR corner (90N,40W)
      !            Lon_LL  Lat_LL  Lon_UR  Lat_UR
      COORDS%r = (/ -165d0,  10d0,  -40d0,   90d0  /)
      COORDS%i = 0d0
      ! Get the indices corresponding to the lon/lat values in COORDS
      CALL GET_BOUNDING_BOX( COORDS, INDICES )

      ! Copy values from INDEX array to scalars
      I1_NA = INDICES(1)
      J1_NA = INDICES(2)
      I2_NA = INDICES(3)
      J2_NA = INDICES(4)

#endif
              
      ! Reset IS_INIT before exiting
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_CARBON

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_CARBON
!
!******************************************************************************
!  Subroutine CLEANUP_CARBON deallocates all module arrays 
!  (rjp, bmy, 4/1/04, 7/8/04)
!
!  NOTES:
!  (1 ) Now deallocate arrays for secondary organic aerosols (rjp, bmy, 7/8/04)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_CARBON begins here!
      !=================================================================
      IF ( ALLOCATED( ANTH_BLKC ) ) DEALLOCATE( ANTH_BLKC )
      IF ( ALLOCATED( ANTH_ORGC ) ) DEALLOCATE( ANTH_ORGC )
      IF ( ALLOCATED( BIOB_BLKC ) ) DEALLOCATE( BIOB_BLKC )
      IF ( ALLOCATED( BIOB_ORGC ) ) DEALLOCATE( BIOB_ORGC )
      IF ( ALLOCATED( BIOF_BLKC ) ) DEALLOCATE( BIOF_BLKC )
      IF ( ALLOCATED( BIOF_ORGC ) ) DEALLOCATE( BIOF_ORGC )
      IF ( ALLOCATED( TERP_ORGC ) ) DEALLOCATE( TERP_ORGC )
      IF ( ALLOCATED( BCCONV    ) ) DEALLOCATE( BCCONV    )
      IF ( ALLOCATED( OCCONV    ) ) DEALLOCATE( OCCONV    )
      IF ( ALLOCATED( EF_BLKC   ) ) DEALLOCATE( EF_BLKC   )
      IF ( ALLOCATED( EF_ORGC   ) ) DEALLOCATE( EF_ORGC   )
      IF ( ALLOCATED( BIOG_ALPH ) ) DEALLOCATE( BIOG_ALPH )
      IF ( ALLOCATED( BIOG_LIMO ) ) DEALLOCATE( BIOG_LIMO )
      IF ( ALLOCATED( BIOG_ALCO ) ) DEALLOCATE( BIOG_ALCO )
      IF ( ALLOCATED( BIOG_TERP ) ) DEALLOCATE( BIOG_TERP )
      IF ( ALLOCATED( BIOG_SESQ ) ) DEALLOCATE( BIOG_SESQ )
      IF ( ALLOCATED( DIUR_ORVC ) ) DEALLOCATE( DIUR_ORVC )
      IF ( ALLOCATED( GEIA_ORVC ) ) DEALLOCATE( GEIA_ORVC )
      IF ( ALLOCATED( TCOSZ     ) ) DEALLOCATE( TCOSZ     )
      IF ( ALLOCATED( GPROD     ) ) DEALLOCATE( GPROD     )
      IF ( ALLOCATED( APROD     ) ) DEALLOCATE( APROD     )
      IF ( ALLOCATED( ORVC_TERP ) ) DEALLOCATE( ORVC_TERP )
      IF ( ALLOCATED( ORVC_SESQ ) ) DEALLOCATE( ORVC_SESQ )

      IF ( ALLOCATED( VCLDF     ) ) DEALLOCATE( VCLDF     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_CARBON

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CARBON_MOD
