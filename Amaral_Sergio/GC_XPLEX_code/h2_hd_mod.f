! $Id: h2_hd_mod.f,v 1.1.1.1 2009/06/09 21:51:51 daven Exp $
      MODULE H2_HD_MOD
!
!******************************************************************************
!  Module H2_HD_MOD contains variables and routines used for the 
!  geographically tagged H2-HD simulation. (lyj, hup, phs, bmy, 9/18/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) SUMISOPCO   : Array for production of CO from Isoprene
!  (2 ) SUMMONOCO   : Array for production of CO from Monoterpenes
!  (3 ) SUMCH3OHCO  : Array for production of CO from CH3OH (methanol)
!  (4 ) SUMACETCO   : Array for production of CO from Acetone
!* (5 ) EMACET      : Array for hold monthly mean acetone emissions
!  (8 ) FMOL_H2     : molecular weight of H2
!  (9 ) XNUMOL_H2   : molec H2 / kg H2
!  (10) FMOL_ISOP   : molecular weight of ISOP
!  (11) XNUMOL_ISOP : molec ISOP / kg ISOP
!  (12) FMOL_MONO   : molecular weight of MONOTERPENES
!  (13) XNUMOL_MONO : molec MONOTERPENES / kg MONOTERPENES
!  (14) EMOCEAN     : Array for hold monthly mean ocean H2 emissions
!  (15) H2CO_YIELD  : Array for photochemical yield of H2 vs CO
!
! * = shared w/ tagged_co_mod.f where it belongs.
!
!  Module Routines:
!  ============================================================================
!  (1 ) EMISS_H2_HD      : Emissions of H2 and HD
!  (2 ) CHEM_H2_HD       : Does chemistry for H2 and HD tracers
!  (3 ) INIT_H2_HD       : Allocates and initializes module arrays
!  (4 ) CLEANUP_H2_HD    : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by h2_hd_mod.f
!  ============================================================================
!  (1 ) biofuel_mod.f    : Module w/ routines to read biofuel emissions
!  (2 ) biomass_mod.f    : Module w/ routines to read biomass emissions
!  (3 ) bpch2_mod.f      : Module w/ routines for binary punch file I/O
!  (4 ) dao_mod.f        : Module w/ arrays for DAO met fields
!  (5 ) diag_mod.f       : Module w/ GEOS-CHEM diagnostic arrays
!  (6 ) directory_mod.f  : Module w/ GEOS-CHEM data & met field dirs
!  (7 ) error_mod.f      : Module w/ I/O error and NaN check routines
!  (8 ) geia_mod         : Module w/ routines to read anthro emissions
!  (9 ) global_oh_mod.f  : Module w/ routines to read 3-D OH field
!  (10) global_nox_mod.f : Module w/ routines to read 3-D NOx field
!  (11) global_ch4_mod.f : Module containing routines to read 3-D CH4 field
!  (12) global_o1d_mod.f : Module containing routines to read 3-D O1D field
!  (13) grid_mod.f       : Module w/ horizontal grid information
!  (14) logical_mod.f    : Module w/ GEOS-CHEM logical switches
!  (15) pressure_mod.f   : Module w/ routines to compute P(I,J,L)
!  (16) tagged_co_mod.f  : Module w/ CO arrays and routines
!  (16) time_mod.f       : Module w/ routines for computing time & date
!  (17) tracer_mod.f     : Module w/ GEOS-CHEM tracer array STT etc.
!  (18) tropopause_mod.f : Module w/ routines to read ann mean tropopause
! 
!  Tracers 
!  ============================================================================
!  (1 ) Total H2
!  (2)  Total HD
!
!  NOTES:
!  (1 ) Based on "tagged_co_mod.f" (lyj, bmy, phs, 5/10/07)
!******************************************************************************
!
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "h2_hd_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE SUMCH3OHCO
      PRIVATE SUMISOPCO, SUMMONOCO,   SUMACETCO
      PRIVATE FMOL_H2,   XNUMOL_H2
      PRIVATE FMOL_HD,   XNUMOL_HD,   EMOCEAN,   H2CO_YIELD
      PRIVATE FMOL_ISOP, XNUMOL_ISOP, FMOL_MONO, XNUMOL_MONO
      
      ! PRIVATE module routines
      PRIVATE INIT_H2_HD,       READ_OCEAN_H2
      PRIVATE READ_H2YIELD

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      TYPE (XPLEX),  ALLOCATABLE :: SUMCH3OHCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMISOPCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMMONOCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: SUMACETCO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMOCEAN(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: H2CO_YIELD(:,:,:)

      ! FMOL_H2     = kg H2 / mole H2
      ! XNUMOL_H2   = molecules H2 / kg H2
      TYPE (XPLEX),  PARAMETER   :: FMOL_H2     = xplex(2d-3,0d0)
      TYPE (XPLEX),  PARAMETER   :: XNUMOL_H2   = xplex(6.022d+23/
     & FMOL_H2%r,0d0)
      
      ! FMOL_HD     = kg HD / mole HD
      ! XNUMOL_HD   = molecules HD / kg HD
      TYPE (XPLEX),  PARAMETER   :: FMOL_HD     = xplex(3d-3,0d0)
      TYPE (XPLEX),  PARAMETER   :: XNUMOL_HD   = xplex(6.022d+23/
     & FMOL_HD%r,0d0)

      ! FMOL_ISOP   - kg ISOP / mole ISOP
      ! XNUMOL_ISOP - molecules CO / kg ISOP
      TYPE (XPLEX),  PARAMETER   :: FMOL_ISOP   = xplex(12d-3,0d0)
      TYPE (XPLEX),  PARAMETER   :: XNUMOL_ISOP = xplex(6.022d+23/
     & FMOL_ISOP%r,0d0)

      ! FMOL_MONO   - kg MONO / mole MONO
      ! XNUMOL_MONO - molecules MONO / kg MONO
      TYPE (XPLEX),  PARAMETER   :: FMOL_MONO   = xplex(12d-3,0d0)
      TYPE (XPLEX),  PARAMETER   :: XNUMOL_MONO = xplex(6.022d+23/
     & FMOL_MONO%r,0d0)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE EMISS_H2_HD
!
!******************************************************************************
!  Subroutine EMISS_H2_HD reads in emissions for the H2/HD simulation.
!  (lyj, phs, bmy, 9/18/07)
!
!  NOTES:
!    (1 ) Now references GET_ANNUAL_SCALAR (phs, 3/11/08)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOFUEL_MOD,   ONLY : BIOFUEL,     BIOFUEL_BURN
      USE BIOMASS_MOD,   ONLY : BIOMASS,     IDBCO
      USE DAO_MOD,       ONLY : SUNCOS,      BXHEIGHT
      USE DIAG_MOD,      ONLY : AD29,        AD46,          AD10em
      USE GEIA_MOD,      ONLY : GET_IHOUR,   GET_DAY_INDEX, READ_GEIA
      USE GEIA_MOD,      ONLY : READ_LIQCO2, READ_TOTCO2,   READ_TODX
      USE GRID_MOD,      ONLY : GET_XOFFSET, GET_YOFFSET,   GET_AREA_CM2
      USE LOGICAL_MOD,      ONLY : LANTHRO,     LGFED2BB
      USE LOGICAL_MOD,      ONLY : LBIOMASS,    LBIOFUEL,      LNEI99
      USE LOGICAL_MOD,      ONLY : LSTREETS,    LEDGAR,        LBRAVO
      USE TIME_MOD,         ONLY : GET_MONTH,   GET_TAU 
      USE TIME_MOD,         ONLY : GET_YEAR,    GET_TS_EMIS  
      USE TRACER_MOD,       ONLY : STT
      USE TRACERID_MOD,     ONLY : IDBFCO,      IDTH2,         IDTHD
      USE TAGGED_CO_MOD,    ONLY : INIT_TAGGED_CO,  READ_ACETONE, EMACET
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! FSCALYR, SCNR89, TODH, EMISTCO
#     include "CMN_DIAG"     ! Diagnostic arrays & switches

      ! Local variables
      INTEGER                :: I, J, L, N, I0, J0
      INTEGER                :: AS, IREF, JREF, IJLOOP
      INTEGER                :: SCALEYEAR, IHOUR, NTAU, MONTH

      ! SAVED variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE.
      INTEGER, SAVE          :: LASTYEAR = -999, LASTMONTH = -999

      ! For now these are defined in CMN_O3
      !TYPE (XPLEX)                 :: EMISTCO(IGLOB,JGLOB)
      !TYPE (XPLEX)                 :: FLIQCO2(IGLOB,JGLOB)

      TYPE (XPLEX)                 :: TMMP, EMXX,   EMX,  EMMO,   EMME   
      TYPE (XPLEX)            :: EMAC, SFAC89, E_CO, E_H2, E_HD, DTSRCE 
      TYPE (XPLEX)                 :: AREA_CM2, EMOC
      TYPE (XPLEX)                 :: CONVERT(NVEGTYPE) 
      TYPE (XPLEX)                 :: GMONOT(NVEGTYPE)

      ! External functions
      TYPE (XPLEX), EXTERNAL       :: XLTMMP,  EMISOP, BOXVL
      TYPE (XPLEX), EXTERNAL       :: EMMONOT, EMCH3OH

      !=================================================================
      ! EMISS_H2_HD begins here!
      !
      ! Do the following only on the first call to EMISS_H2_HD...
      !=================================================================
      IF ( FIRSTEMISS ) THEN

         ! Read polynomial coeffs' for isoprene emissions
         CALL RDLIGHT

         ! Read conversion tables for isoprene & monoterpene emissions
         ! Also read acetone emissions (bnd, bmy, 6/8/01)
         CALL RDISOPT( CONVERT )
         CALL RDMONOT( GMONOT  ) 

         ! Set the base level of isoprene & monoterpene emissions
         CALL SETBASE( CONVERT, GMONOT )

         ! Read time-of-day and day-of-week scale factors for GEIA emissions
         CALL READ_TODX( TODN, TODH, TODB, SCNR89 )

         ! Read anthropogenic CO emissions from GEIA
         CALL READ_GEIA( E_CO=EMISTCO )

         ! Allocate all module arrays
         CALL INIT_H2_HD
         CALL INIT_TAGGED_CO

         !! Define geographic regions for both fossil fuel & biomass burning
         !CALL DEFINE_FF_CO_REGIONS( FF_REGION )         
         !CALL DEFINE_BB_CO_REGIONS( BB_REGION )

         ! Set first-time flag to false
         FIRSTEMISS = .FALSE.
      ENDIF

      MONTH = GET_MONTH()

      !=================================================================
      ! Once a month, read acetone from disk.  For GEOS-3, also read 
      ! P(CO) from ISOPRENE, MONOTERPENES, and METHANOL from 1996.
      ! Also read ocean H2 source and the relative H2/CO 
      ! photochemical yield.
      !=================================================================
      IF ( MONTH /= LASTMONTH ) THEN
      
         ! Read acetone for this month
         CALL READ_ACETONE( MONTH )
      
         ! Read ocean emissions
         CALL READ_OCEAN_H2( MONTH )

         ! Read H2/CO photochemical yield
         CALL READ_H2YIELD( MONTH )

         ! Save month for next iteration
         LASTMONTH = MONTH
      ENDIF

      !=================================================================
      ! If FSCALYR < 0 then use this year (JYEAR) for scaling the
      ! fossil fuel emissions.  Otherwise, use the value of FSCALYR 
      ! as specified in 'input.ctm'.
      !
      ! Modified to use new scaling factor (phs, 3/11/08)
      !=================================================================
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
!------------------
! prior to 3/11/08
!         ! Cap SCALEYEAR at 1998 for now (bmy, 1/13/03) 
!         IF ( SCALEYEAR > 1998 ) SCALEYEAR = 1998
!------------------
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      IF ( SCALEYEAR /= LASTYEAR ) THEN
!------------------
! prior to 3/11/08
!         CALL READ_LIQCO2( SCALEYEAR, FLIQCO2 )
!-----------------
         ! now use updated scalars (phs, 3/11/08)
         CALL GET_ANNUAL_SCALAR( 72, 1985, SCALEYEAR, FLIQCO2 )
         LASTYEAR = SCALEYEAR
      ENDIF

      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      !=================================================================
      ! Process Anthropogenic (Fossil Fuel) H2 emissions based on
      ! CO emissions scaled by H2/CO emission ratio of 0.042 H2/CO 
      ! from the GEIA/Piccot inventory (Hauglustaine and Ehhalt, 2002).
      !
      ! Anthropogenic emissions are enhanced by 18.5% below.  This
      ! accounts for production of H2 from oxidation of certain VOC's,
      ! which are not explicitly carried by GEOS-CHEM as anthropogenic
      ! species.  This needs to be done here since a different scale
      ! factor for the full chemistry run is used. The 18.5% is further
      ! multiplied by the relative H2/CO photochemical yield.  Also update
      ! the ND29 diagnostic below, in order to archive the correct 
      ! emissions. (bmy, 6/14/01)
      !
      ! For HD, we include an isotopic signature of -196 permil from
      ! Quay and Gerst (2001).
      !
      ! NOTES: 
      ! (1) Anthro CO emissions come from the GEIA/Piccot inventory.
      !     (bmy, 1/2/01)
      ! (2) Need to save ND29 diagnostics here (bmy, 1/2/01)
      !=================================================================
      IF ( LANTHRO ) THEN

         ! NTAU is just the integral value of TAU (ave, bmy, 6/10/03)
         NTAU = GET_TAU()

         ! SFAC89 is the Weekday/Saturday/Sunday scale factor
         SFAC89 = SCNR89( 2, GET_DAY_INDEX( NTAU ) ) 


!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( E_CO, E_H2, E_HD, AREA_CM2, I, IHOUR, IREF, J, JREF, N )
         DO J = 1, JJPAR
            JREF = J + J0

            ! Grid box surface areas [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            DO I = 1, IIPAR
               IREF = I + I0

               ! E_CO is FF CO emissions in [molec CO/cm^2/s]
               ! Scale E_CO by the day-of-week scale factor SFAC89
               E_CO = EMISTCO(IREF,JREF) * SFAC89

               ! Scale E_CO by the time-of-day scale factor TODH
               ! IHOUR is the index for the time-of-day scale factor TODH
               IHOUR = GET_IHOUR( I )
               E_CO  = E_CO * TODH(IHOUR)

               ! Scale E_CO by FLIQCO2, which reflects the yearly 
               ! increase in FF CO emissions for each country
               E_CO = E_CO * FLIQCO2(IREF,JREF)


               !%%%  Need to also overwrite emissions with EDGAR etc %%% 
               IF ( LEDGAR ) THEN
                  ! etc
               ENDIF
      
               IF ( LSTREETS ) THEN
                  ! etc
               ENDIF
      
               IF ( LBRAVO ) THEN
                  ! etc
               ENDIF
      
               IF ( LNEI99 ) THEN
                  ! etc
               ENDIF


               ! ND29 diagnostic -- store Fossil Fuel CO [molec/cm2/s] 
               IF ( ND29 > 0 ) THEN
                  AD29(I,J,1) = AD29(I,J,1) + E_CO * 1.185d0
               ENDIF

               ! To obtain H2 FF emissions, scale E_CO by Anthropogenic ratio 
               ! H2/CO of 0.042 from Hauglustaine and Ehhalt, 2002
               ! Convert from mass to molecules
	       ! 0.042 gH2/gCO x (28 mol/g)/(2 mol/g) = 0.588 molec H2/CO 
               ! Enhance H2 by 18.5% to account for oxidation
               ! from Anthropogenic VOC's (bnd, bmy, 6/8/01) and
               ! multiply with H2/CO photochemical yield
               ! (hup, jaegle 4/20/2007)
               E_H2 = E_CO * ( 0.185d0 * H2CO_YIELD(I,J,1) + 0.588d0 )
	       
	       ! Calculate FF emissions for HD, using the -196 permil 
               ! signature measured by (Gerst & Quay, 2001).
               ! Convert permil to D/H ratio
               ! (D/H)ff = (delD/1000d0+1d0)*vsmow =  1.2523d-4 
               ! with vsmow = 155.76d-6 [ Vienna Mean Standard Ocean
               ! Water] and delD = -196 permil
               ! We further need to multiply by 2 to get the DH/H2
               ! ratio (we count the hydrogens vs the deuteriums).
               ! (DH/H2)ff = 1.2523d-4 * 2.d0
	       ! (hup, jaegle, 11/16/2003) 
	       E_HD = E_H2 * 1.2523d-4 * 2.d0

               ! ND10 diagnostic -- store Fossil Fuel H2,HD [molec/cm2/s]
               IF ( ND10 > 0 ) THEN
                  AD10em(I,J,1) = AD10em(I,J,1) + E_H2
               ENDIF
	       
               ! Convert E_H2 from [molec H2/cm2/s] to [kg H2]
               E_H2 = E_H2 * ( AREA_CM2 * DTSRCE / XNUMOL_H2 )
	       
	       ! Convert E_HD from [molec HD/cm2/s] to [kg HD]
	       E_HD = E_HD * ( AREA_CM2 * DTSRCE / XNUMOL_HD )

               ! Add FF H2 to Tracer #1 -- total H2 [kg H2]
	       ! Add FF HD to Tracer #2 -- total HD [kg HD]
               STT(I,J,1,IDTH2) = STT(I,J,1,IDTH2) + E_H2
	       STT(I,J,1,IDTHD) = STT(I,J,1,IDTHD) + E_HD

            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Process biomass burning CO emissions, stored in array
      ! BIOMASS(:,:,IDBCO) which has units of [molec/cm2/s]
      ! To obtain the H2 biomass burning emissions, the CO emissions
      ! are scaled by a molar ratio of 0.29 molH2/molCO based on
      ! Andreae and Merlet (2001, Table 2). This value is obtained
      ! by taking the mean of Savannah, Grassland, Tropical Forest
      ! and Extratropical Forest H2/CO ratios weighted by Dry
      ! Matter Burned.
      !
      ! The default Duncan et al 2001 biomass burning emissions are
      ! enhanced by 11.% within here to account for the VOC.
      !
      ! For HD, we include an isotopic signature of -290 permil from
      ! Quay and Gerst (2001).
      !
      ! GFED2 CO biomass burning emissions are not scaled in a full
      ! chemistry, meaning that VOC oxidation is already included.
      ! Then formula below for CO->H2 is modified to account for 
      ! that assumption (phs).
      !
      ! NOTES:
      ! (1) Some forest fires generate strong convection columns.  
      !      However, we release biomass burning emissions only into 
      !      the surface layer. (bnd, bmy, 1/3/01)
      !=================================================================
      IF ( LBIOMASS ) THEN

!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( E_H2, E_HD, I, J, N, AREA_CM2 )

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Loop over longitudes
            DO I = 1, IIPAR

            ! Convert [molec CO/cm2/s] to [molec H2/cm2/s]
            ! Scale by 0.29 mol H2/mol CO and increase by 
            ! 0.11*H2CO_YIELD to account for voc oxidation.
            IF ( LGFED2BB ) THEN
               BIOMASS(I,J,IDBCO) = BIOMASS(I,J,IDBCO) *
     &              ( H2CO_YIELD(I,J,1) + ( 0.29d0 / 0.11d0 ) )
            ELSE
               BIOMASS(I,J,IDBCO) = BIOMASS(I,J,IDBCO) *
     &              ( ( 0.11d0 * H2CO_YIELD(I,J,1) ) + 0.29d0 )
            ENDIF


            ! ND10 diagnostic -- store biomass burning of H2 [molec H2/cm2/s]
            IF ( ND10 > 0 ) AD10em(I,J,2) = AD10em(I,J,2) + 
     &                                      BIOMASS(I,J,IDBCO)

            ! Convert [molec H2/cm2/s] to [kg H2] and store in E_H2.
            E_H2 = ( BIOMASS(I,J,IDBCO) / XNUMOL_H2 ) *
     &             ( AREA_CM2           * DTSRCE    )

           ! Convert [molec H2/cm3/s] to [kg HD] and store in E_HD.
	   ! Scale E_HD by biomass burning ratio 1.1d-4 molecules HD/H2 
	   ! accd isotopic signature of -293permil (Gerst & Quay, 2001)
           ! change to -290 permil (as in Gerst and Quay) - jaegle
           ! add missing factor of 2
	   ! Calculate BF emissions for HD, using the -290 permil 
           ! isotopic signature measured by (Gerst & Quay, 2001).
           ! Convert permil to D/H ratio
           ! (D/H)ff = (delD/1000d0+1d0)*vsmow =  1.1059d-4
           ! with vsmow = 155.76d-6 [ Vienna Mean Standard Ocean
           ! Water] and delD = -290 permil
           ! We further need to multiply by 2 to get the DH/H2
           ! ratio (we count the hydrogens vs the deuteriums).
           ! (DH/H2)ff = 1.1059d-4 * 2.d0
           
            E_HD = ( ( BIOMASS(I,J,IDBCO) * 1.1059d-4 * 2.d0) 
     &                / XNUMOL_HD ) * ( AREA_CM2 * DTSRCE )

            ! Add H2 HD biomass burning to corresponding tracers -
            ! - total H2/HD [kg H2/HD]
            STT(I,J,1,IDTH2) = STT(I,J,1,IDTH2) + E_H2
            STT(I,J,1,IDTHD) = STT(I,J,1,IDTHD) + E_HD

            ENDDO
         ENDDO
!!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Process biofuel (formerly wood burning) CO emissions
      ! stored in BIOFUEL(IDBCO,IREF,JREF) in [molec/cm3/s]
      !
      ! Biofuel burning emissions must be enhanced by 18.9% 
      ! to account for CO production from oxidation of certain VOC's, 
      ! which are not explicitly carried by GEOS-CHEM as biofuel burning
      ! species.
      ! In case of H2/HD simulation,the scaling is done here.
      !
      ! Scale CO emissions by 0.322 molH2/molCO from 
      ! Andreae and Merlet [2001].
      ! Scaled by the relative H2 to CO photochemical yield.
      !
      ! For HD, we include an isotopic signature of -290 permil from
      ! Quay and Gerst (2001) - we assume that biofuel isotopic
      ! signature is the same as biomass burning.
      !
      ! NOTES:
      ! (1 ) Use IDBFCO to index the proper element of the 
      !       biofuel burning array (bmy, 6/8/01).
      !=================================================================
      IF ( LBIOFUEL ) THEN
         CALL BIOFUEL_BURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( E_H2, E_HD, I, J, N )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Convert from [molec CO/cm3/s] to [kg H2]
            ! BIOFUEL(IDBFCO,I,J) contains biofuel CO emissions.
            ! Scale by 0.322d0 H2/CO from Andreae and Merlet
            ! and enhance by 18.9% * yield h2/co from photochemical
            ! oxidation (lyj, 2/10/07)
            E_H2 = BIOFUEL(IDBFCO,I,J)  / XNUMOL_H2 *
     &             BOXVL(I,J,1)        * DTSRCE  *
     &             ( ( 0.189d0 * H2CO_YIELD(I,J,1) ) + 0.322d0 )


            ! Calculate the HD emissions using a -290 permil
            ! isotopic signature (see BB emissions above).
            ! jaegle (4/20/07)

            E_HD = ( BIOFUEL(IDBFCO,I,J) * 1.1059d-4 ) * 
     &             ( ( 0.189d0 * H2CO_YIELD(I,J,1) ) + 0.322d0 ) /
     &             XNUMOL_HD * ( BOXVL(I,J,1) * DTSRCE ) * 2.d0

            ! ND10 -- store Biofuel Fuel H2 [molec/cm2/s]
            IF ( ND10 > 0 ) THEN
               AD10em(I,J,3) = AD10em(I,J,3) + ( BIOFUEL(IDBFCO,I,J) *
     &                 ( ( 0.189d0 * H2CO_YIELD(I,J,1) ) + 0.322d0 ) *
     &                       BXHEIGHT(I,J,1)*100d0 )
            ENDIF 

            ! Add H2 HD biomass burning to corresponding tracers -
            ! - total H2/HD [kg H2/HD]
            STT(I,J,1,IDTH2) = STT(I,J,1,IDTH2) + E_H2
            STT(I,J,1,IDTHD) = STT(I,J,1,IDTHD) + E_HD

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Process emissions of ISOPRENE, MONOTERPENES, METHANOL
      ! and ACETONE -- save into summing arrays for later use
      ! Also process ocean emissions for H2.
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, IJLOOP, TMMP, EMXX, EMMO, EMAC, EMOC )
!$OMP+PRIVATE( E_H2, E_HD )
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

         DO I = 1, IIPAR

            ! 1-D array index 
            IJLOOP = ( (J-1) * IIPAR ) + I

            !===========================================================
            ! The CO and H2 yields from ISOP, MONOTERPENES, and CH3OH will be
            ! computed in subroutine CHEM_H2_HD.  P(CO) from CH3OH
            ! will be scaled to isoprene emissions within subroutine
            ! CHEM_H2_HD 
            !===========================================================

            ! Surface air temperature [K]
            TMMP = XLTMMP(I,J,IJLOOP)
      
            ! ISOP and MONOTERPENE emissions in [atoms C/box/time step] 
            ! SUNCOS is COSINE( solar zenith angle ) 
            EMXX = EMISOP( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL_ISOP ) 
            EMMO = EMMONOT( IJLOOP, TMMP, XNUMOL_MONO )

            ! Store ISOP and MONOTERPENE emissions [atoms C/box/time step]
            ! for later use in the subroutine CHEM_H2_HD
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

	    !===========================================================
            ! For GEOS-4 (?) and GEOS-3, extract ocean emissions
            ! fluxes the EMOCEAN array for the current month
            !===========================================================

            ! EMOC = [molec/cm2/s] Ocean emissions of H2 from N-fixation
            EMOC = EMOCEAN( I, J )
	    
	    ! Calculate HD ocean source (molec HD/cm2/s)
            ! assume delD=-628 permil
            ! add missing factor of 2 (jaegle 12/23/05)
	    
	    E_HD = EMOC * 5.79427d-5 * 2.d0
	    
            ! diag 10 
            IF ( ND10 > 0 ) THEN
                 AD10em(I,J,4) = AD10em(I,J,4) + EMOC
                 AD10em(I,J,5) = AD10em(I,J,5) + E_HD	  
            ENDIF

            ! Convert ocean H2 source from [molec H2/cm2/s] to [kg H2]
            E_H2 =  EMOC  * ( AREA_CM2 * DTSRCE / XNUMOL_H2 ) 
            
	    ! Scale E_HD by ocean isotopic signature ratio 
	    ! of -628 permil (Rice & Quay, 2007)
            ! (DH/H2)ocean = 5.79427d-5 * 2
            
            E_HD =  ( EMOC * 5.79427d-5 * 2d0) * 
     &              ( AREA_CM2 * DTSRCE / XNUMOL_HD )  
            
            ! Add H2 HD biomass burning to corresponding tracers -
            ! - total H2/HD [kg H2/HD]
            STT(I,J,1,IDTH2) = STT(I,J,1,IDTH2) + E_H2
            STT(I,J,1,IDTHD) = STT(I,J,1,IDTHD) + E_HD

            
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMISS_H2_HD

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_H2_HD
!
!******************************************************************************
!  Subroutine CHEM_H2_HD performs H2 and HD chemistry.  Chemical production is
!  by oxidation of BVOC and CH4.  Loss is via reaction with OH and uptake by 
!  soils. In the stratosphere, H2 is also lost by reaction with O(1D).  For 
!  HD, we include the fractionation from photochemical oxidation (162 permil),
!  and loss by OH and soil uptake. (lyj, hup, phs, 9/18/07)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD, AIRVOL, T
      USE DIAG_MOD,       ONLY : AD10
      USE ERROR_MOD,      ONLY : CHECK_VALUE
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH,   OH
      USE GLOBAL_O1D_MOD, ONLY : GET_GLOBAL_O1D, O1D
      USE GLOBAL_NOX_MOD, ONLY : GET_GLOBAL_NOX,  BNOX
      USE GRID_MOD,       ONLY : GET_YMID, GET_AREA_M2, GET_AREA_CM2
      USE PRESSURE_MOD,   ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,       ONLY : GET_TS_CHEM,     GET_MONTH, GET_YEAR
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE TRACER_MOD,     ONLY : N_TRACERS,       STT
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
      USE TRACERID_MOD,   ONLY : IDTH2, IDTHD
      USE TAGGED_CO_MOD,  ONLY : GET_ALPHA_ISOP,    READ_PCO_LCO_STRAT
      USE TAGGED_CO_MOD,  ONLY : GET_PCO_LCO_STRAT

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND
#     include "CMN_DIAG"     ! ND65

      ! Local variables
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, N, MONTH
      TYPE (XPLEX)                 :: ALPHA_CH4, ALPHA_ISOP, ALPHA_MONO
      TYPE (XPLEX)                 :: DTCHEM,    GH2,        PCO
      TYPE (XPLEX)                 :: GHD,       STTHD
      TYPE (XPLEX)                 :: STTH2,     KRATE,      CH4
      TYPE (XPLEX)                 :: H2_CH4,    H2_ISOP,    H2_MONO
      TYPE (XPLEX)                 :: H2_CH3OH,  H2_OH,      H2_ACET
      TYPE (XPLEX)                 :: HD_OH,     O1DRATE,    HD_RATE
      TYPE (XPLEX)                 :: HD_O1D,    H2_O1D,     HD_CH4
      TYPE (XPLEX)                 :: CH4RATE,   DENS,       ALPHA_ACET
      TYPE (XPLEX)                 :: CORATE,    YMID,       DPHOTO
      TYPE (XPLEX)                 :: KTEST(IIPAR, JJPAR)
      TYPE (XPLEX)                 :: AREA_CM2, SVEL, FLUX
      TYPE (XPLEX)                 :: THIK, P2, DRYF, TESTVAR
      TYPE (XPLEX)                 :: FRACLOST, H2_RATE

      ! For saving CH4 latitudinal gradient
      TYPE (XPLEX),  SAVE          :: A3090S, A0030S, A0030N, A3090N

      ! External functions
      TYPE (XPLEX),  EXTERNAL      :: BOXVL

      ! WTAIR = molecular weight of air (g/mole)
      TYPE (XPLEX),  PARAMETER     :: WTAIR = xplex(28.966d0,0d0)
      
      ! Switch to scale yield of isoprene from NOx concentration or not
      LOGICAL, PARAMETER     :: ALPHA_ISOP_FROM_NOX = .FALSE.

      !=================================================================
      ! CHEM_H2_HD begins here!
      !
      ! Do the following on the first calla to CHEM_H2_HD...
      !=================================================================
      IF ( FIRSTCHEM ) THEN
         FIRSTCHEM = .FALSE.
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0
      
      !=================================================================
      ! Read in OH, O1D, NOx, P(CO), and L(CO) fields for the current month
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN
         
         ! Get current month
         MONTH = GET_MONTH()

         ! Global OH 
         CALL GET_GLOBAL_OH( MONTH )

	 ! Global O1D
	 CALL GET_GLOBAL_O1D( MONTH )

         ! Global NOx -- need this to determine 
         ! ALPHA_ISOP which is a function of NOx
         IF ( ALPHA_ISOP_FROM_NOX ) CALL GET_GLOBAL_NOX( MONTH )
            
         ! Read in the loss/production of CO in the stratosphere.
         CALL READ_PCO_LCO_STRAT( MONTH )
      ENDIF

      !=================================================================
      ! Get the yearly and latitudinal gradients for CH4
      ! This only needs to be called once per year
      !=================================================================
      IF ( ITS_A_NEW_YEAR() ) THEN 
         CALL GET_GLOBAL_CH4( GET_YEAR(), .TRUE., 
     &                        A3090S, A0030S, A0030N, A3090N )
      ENDIF

      !=================================================================
      ! Do H2 and HD chemistry 
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED                                             ) 
!$OMP+PRIVATE( I, J, L, N, STTH2, GH2, DENS, CH4RATE, H2_CH4, CH4 )
!$OMP+PRIVATE( ALPHA_CH4, KRATE, H2_ISOP, H2_CH3OH, ALPHA_ISOP    )
!$OMP+PRIVATE( H2_MONO, ALPHA_MONO, H2_ACET, ALPHA_ACET, CORATE   )
!$OMP+PRIVATE( H2_OH, YMID, DPHOTO, SVEL, FLUX, AREA_CM2    )
!$OMP+PRIVATE( GHD, STTHD, HD_OH, HD_CH4, HD_O1D, H2_RATE, HD_RATE)
!$OMP+PRIVATE( H2_O1D, DRYF, P2, THIK, FRACLOST  )
      
      DO L = 1, LLPAR
      DO J = 1, JJPAR

         ! Latitude of grid box
         YMID = GET_YMID( J )

      DO I = 1, IIPAR

         !==============================================================
         ! (0) Define useful quantities
         !==============================================================

         ! STTH2 [molec H2/cm3/kg H2] converts [kg H2] --> [molec H2/cm3]
         ! kg H2/box * box/cm3 * mole/0.002 kg H2 * Avog.#/mole
         STTH2 = 1d0 / AIRVOL(I,J,L) / 1d6 / FMOL_H2 * 6.023d23
	 
	 STTHD = 1d0 / AIRVOL(I,J,L) / 1d6 / FMOL_HD * 6.023d23 
                  
         ! GH2 is H2 concentration in [molec H2/cm3]
         GH2   = STT(I,J,L,IDTH2) * STTH2

	 ! GHD is HD concentration in [molec HD/cm3]
	 GHD   = STT(I,J,L,IDTHD) * STTHD

         ! DENS is the number density of air [molec air/cm3]
         DENS  = AD(I,J,L) * 1000.d0 / BOXVL(I,J,L) * 6.023d23  / WTAIR

         !==============================================================
         ! (1a) Production of H2 and HD by methane oxidation
         !==============================================================

         ! Initialize
         H2_CH4 = 0d0

         ! Test level for stratosphere or troposphere
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

            !===========================================================
            ! (1a-1) Production of H2 from CH4 in the stratosphere
            !        This is based on the CO production rates in the
            !        stratosphere, scaled by a 0.43 H2/CO yield.
            !        For HD, we deal with the stratosphere in
            !        upbdflx_mod.f, using an adapted version of
            !        SYNOZ.
            !===========================================================

            ! Call GET_PCO_LCO_STRAT to get the P(CO) rate from CH4
            CH4RATE = GET_PCO_LCO_STRAT( .TRUE., I, J, L )
            
            ! Convert units of CH4RATE from [v/v/s] to [molec H2/cm3]
            ! In the stratosphere, we use a H2/CO yield of 0.43
	    H2_CH4 = CH4RATE * DTCHEM * DENS * 0.43d0 

            ! no HD from CH4
            HD_CH4 = 0d0

            ! Set stratospheric signature of photochemical oxidation to
            ! zero (lyj, 01/13/06)
            DPHOTO = 0d0

         ELSE
            
            !===========================================================
            ! (1a-2) Production of H2 from CH4 in the troposphere
	    ! From the photocissoc coef of formaldehyde one finds that
	    ! the channel leading to H2 + CO as dissoc products contrib
	    ! roughly 50% of the yield of CO accd to Novelli, 1999
	    ! (hup, 5/27/2004)
            !  We assume an isotopic fractionation signature from
            !  photochemical oxidation fo 162 permil. See
            !  Price et al. [2007]
            !===========================================================

            ! CH4 concentration [ppbv] for the given latitude band 
            ! (bmy, 1/2/01)
            CH4 = A3090S
            IF ( YMID >= -30.0 .and. YMID < 0.0  ) CH4 = A0030S
            IF ( YMID >=   0.0 .and. YMID < 30.0 ) CH4 = A0030N
            IF ( YMID >=  30.0                   ) CH4 = A3090N
            
            ! Convert CH4 from [ppbv] to [molec CH4/cm3]
            CH4       = CH4 * 1d-9 * DENS

            ! Yield of CO from CH4: estimated to be 95-100% (acf)
            ALPHA_CH4 = 1d0

            ! Calculate updated rate constant [s-1] (bnd, bmy, 1/2/01)
            KRATE     = 2.45D-12 * EXP( -1775.d0 / T(I,J,L) ) 
 
            ! Production of CO from CH4 = alpha * k * [CH4] * [OH] * dt
            ! Scale this by the H2 yield relative to CO
            ! Units are [molec H2/cm3]
            H2_CH4    = ALPHA_CH4 *  H2CO_YIELD(I,J,L) * KRATE * 
     &                  CH4 *        OH(I,J,L)         * DTCHEM
	    
            ! HD Photochemical Signature in troposphere of
            ! delD(hv)=162 permil (hup, 8/14/2006)
            ! This means H/D(hv)=1.8099311d-4
            ! and H2/HD(hv) = 1.8099311d-4 * 2.d0
	    DPHOTO = 1.8099311d-4 * 2.d0
	    
	    ! Calculated CH4 oxidation source of HD in troposphere
	    HD_CH4 = H2_CH4 * DPHOTO

         ENDIF
	 
         ! Check H2_CH4 for NaN or Infinity
         CALL CHECK_VALUE( H2_CH4,  (/ I, J, L, 0 /),
     *                    'H2_CH4', 'STOP at h2_hd_mod:1' )


         !==============================================================
         ! (1b) Production of H2 from ISOPRENE and METHANOL (CH3OH)
         !==============================================================

         ! Initialize
         H2_ISOP  = 0d0
         H2_CH3OH = 0d0

         ! Isoprene is emitted only into the surface layer 
         IF ( L == 1 ) THEN 

            !===========================================================
            ! Yield of CO from ISOP: 30%, from Miyoshi et al., 1994.
            ! They estimate globally 105 Tg C/yr of CO is produced
            ! from isoprene oxidation.
            !
            ! Increased this factor from 30% to 50% (bnd, bmy, 1/3/01)
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
            !-----------------------------------------------------------
            ! To obtain H2 sources, we scale all of these BVOC emissions
            ! by the H2/CO photochemical yield. (hup,jaegle, 04/20/07)
            !===========================================================

            ! Get CO yield from isoprene
            IF ( ALPHA_ISOP_FROM_NOX ) THEN
               ALPHA_ISOP = GET_ALPHA_ISOP( .TRUE., BNOX(I,J,L) )         
            ELSE
               ALPHA_ISOP = GET_ALPHA_ISOP( .FALSE. )
            ENDIF

            ! scale ALPHA_ISOP by the relative H2/CO yield
            ALPHA_ISOP = ALPHA_ISOP * H2CO_YIELD(I,J,L)

            ! P(CO) from Isoprene Flux = ALPHA_ISOP * Flux(ISOP)
            ! Convert from [molec ISOP/box] to [molec CO/cm3]
            ! Also account for the fact that ISOP has 5 carbons
            H2_ISOP = SUMISOPCO(I,J) / BOXVL(I,J,L) / 5.d0 * ALPHA_ISOP
	    
            ! P(CO) from CH3OH is scaled to Isoprene Flux (see above)
	    ! Units are [molec CO/cm3]
            ! For H2, scale by the relative CO/H2 molar weights and
            ! the H2/CO yields.
            H2_CH3OH = ( SUMISOPCO(I,J) / BOXVL(I,J,L) ) *
     &                 ( 100d0          / 397d0        ) *
     &                 ( H2CO_YIELD(I,J,L)             ) *
     &                 ( 12d0           / 28d0         ) 
 
            ! Zero SUMISOPCO and SUMCH3OHCO for the next emission step
            SUMISOPCO(I,J)  = 0d0
            SUMCH3OHCO(I,J) = 0d0

            ! Check H2_ISOP for NaN or Infinity
            CALL CHECK_VALUE( H2_ISOP,  (/ I, J, L, 0 /),
     &                       'H2_ISOP', 'STOP at h2_hd_mod:2' )

            ! Check H2_CH4 for NaN or Infinity
            CALL CHECK_VALUE( H2_CH3OH,  (/ I, J, L, 0 /),
     &                       'H2_CH3OH', 'STOP at h2_hd_mod:3' )

         ENDIF

         !==============================================================
         ! (1c) Production of H2 from MONOTERPENE oxidation
         !==============================================================

         ! Initialize
         H2_MONO = 0.d0

         ! Monoterpenes are emitted only into the surface layer
         IF ( L == 1 ) THEN

            !===========================================================
            ! Assume the production of H2 from monoterpenes is 
            ! instantaneous even though the lifetime of intermediate
            ! species may be on the order of hours or days.  This
            ! assumption will likely cause H2 from monoterpene
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
            !-----------------------------------------------------------
            ! For H2, scale by the H2 to CO photochemical yield.
            !===========================================================

            ! Yield of CO from MONOTERPENES: 20% (see above)
            ALPHA_MONO = 0.20d0

            ! P(CO) from Monoterpene Flux =  alpha * Flux(Mono)
            ! Units are [molec H2/cm3]. Scale by the
            ! H2/CO photochemical yield.
            H2_MONO = ( SUMMONOCO(I,J) / BOXVL(I,J,L)      ) * 
     &                ( ALPHA_MONO     * H2CO_YIELD(I,J,L) )

            ! Zero SUMMONOCO for the next emission step
            SUMMONOCO(I,J) = 0d0

            ! Check H2_MONO for NaN or Infinity
            CALL CHECK_VALUE( H2_MONO,  (/ I, J, L, 0 /),
     &                       'H2_MONO', 'STOP at h2_hd_mod:4' )

         ENDIF

         !==============================================================
         ! (1d) Production of H2 from oxidation of ACETONE 
         !
         ! ALPHA_ACET = 2/3 to get a yield for CO.  This accounts
         ! for acetonWRITE (6, *) 'DTC 2', DTC, AREA_CM2e loss 
	 ! from reaction with OH And photolysis.
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
         !-----------------------------------------------------------
         ! For H2, scale by the H2 to CO photochemical yield.
         !==============================================================

         ! Initialize
         H2_ACET = 0.d0

         ! Biogenic acetone sources are emitted only into the surface layer
         IF ( L == 1 ) THEN

            ! Yield of CO from ACETONE: 2/3 (see above)
            ALPHA_ACET = 2.D0 / 3.D0
            
            ! Units are [molec H2/cc]. Scale by the H2/CO yield.
            H2_ACET = SUMACETCO(I,J) / BOXVL(I,J,L) * 
     &                ALPHA_ACET * H2CO_YIELD (I,J,L)

            ! Zero SUMACETCO for the next emission step
            SUMACETCO(I,J) = 0d0

            ! Check H2_ACET for NaN or Infinity
            CALL CHECK_VALUE( H2_ACET,  (/I, J, L, 0 /),
     &                       'H2_ACET', 'STOP at h2_hd_mod:5' )


         ENDIF


         !==============================================================
         ! (2a) Loss of H2 and HD due to chemical reaction w/ OH and O1D
         !==============================================================

         ! Initialize
         H2_OH  = 0.d0
         HD_OH  = 0.d0
         H2_O1D = 0.d0
         HD_O1D = 0.d0

         ! Select out tropospheric or stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

            !===========================================================
            ! (2a-1) Stratospheric loss H2 HD due to chem rxn w/OH & O1D
            !===========================================================

            ! Get the L(CO) rate in the stratosphere in [s-1]
            CORATE = GET_PCO_LCO_STRAT( .FALSE., I, J, L )
	    
	    ! Rate constants for H2+OH and HD+OH from JPL, 2004 (hup, 06/28/04)
            H2_RATE = 5.5D-12 * EXP( -2000.d0 / T(I,J,L) )
            HD_RATE = 5.0D-12 * EXP( -2130.d0 / T(I,J,L) )
	    
            ! H2_OH = Stratospheric loss of H2 by OH [molec/cm3]
	    ! alpha=0.52 for HD-OH in strat 230K (Rockmann et al., 2003)
	    ! alpha changes with temperature so just use the explicit 
	    ! rate constant for HD + OH rxn(hup, 5/27/2005)
            H2_OH = H2_RATE * GH2 * OH(I,J,L) * DTCHEM

            ! For now we overwrite this with zero, as the
            ! stratospheric HD is dealt with simply in
            ! updflux_mod.f  - in the future this could be improved
            ! by including more explicitely the isotopic enrichment
            ! of CH4 in the stratosphere (jaegle, 4/20/07)
	    HD_OH = 0d0
	    
	    ! Stratospheric loss of H2 and HD by O1D. Decay rate is same 
	    !(1.1D-10 cm3 molecule-1 s-1) in stratosphere.(hup, 4/26/2005)    
	    O1DRATE = 1.1D-10
	    
	    H2_O1D = O1DRATE * GH2 * O1D(I,J,L) * DTCHEM

            ! Do not deal with this for HD.
	    HD_O1D = 0d0
	    
	    
            ! Check values for NaN or Infinity
            CALL CHECK_VALUE( H2_OH,  (/ I, J, L, 0 /),
     &                      'H2_OH', 'STOP at h2_hd_mod:6' )
     
     	    CALL CHECK_VALUE( H2_O1D,  (/ I, J, L, 0 /),
     &                      'H2_O1D', 'STOP at h2_hd_mod:6' )

         ELSE

            !===========================================================
            ! (2b-2) Tropospheric loss of H2 due to chemical rxn w/ OH
            !
            !  DECAY RATE
            !  The decay rate (H2_RATE) is calculated by:
	    !
	    !     OH + H2 -> products
	    !     use: H2_RATE = 5.5D-12 * EXP( -2000.d0 / T )
            !     HD + H2 -> products
            !     use: HD_RATE = 5.0D-12 * EXP(-2130.d0 / T )
	    !
            !  H2_RATE has units of cm^3/molec sec
	    !
            !     (hup, 05/28/04)
            !===========================================================

	    ! Loss for H2 and HD by reaction with OH accd JPL, 2004 
	    ! (hup, 06/28/04)
            H2_RATE = 5.5D-12 * EXP( -2000.d0 / T(I,J,L) )
	    HD_RATE = 5.0D-12 * EXP( -2130.d0 / T(I,J,L) )

            ! H2_OH = Tropospheric loss of H2 by OH [molec/cm3]
            H2_OH   = H2_RATE * GH2 * OH(I,J,L) * DTCHEM
	    
	    ! HD_OH = Tropospheric loss of HD by OH [molec/cm3]
            ! Calculated explicitely using the HD rate constant.
            ! This results in a ~0.6 fractionation effect at 300K.
	    HD_OH   = HD_RATE * GHD * OH(I,J,L) * DTCHEM

	    ! HD_O1D and H2_O1D are both zero in the troposphere
	    H2_O1D  = 0d0
            HD_O1D  = 0d0   

         ENDIF
         
         !==============================================================
         ! Save the total chemical production from various sources
         ! into the total H2 tracer STT(I,J,L,1)
         !==============================================================

         ! GH2 is the total H2 before chemistry 
	 ! is applied [molec H2/cm3]. Add production from
         ! oxidation of CH4, monoterpenes, acetone, methanol,
         ! isoprene and remove loss from reactions with OH and O1D
         GH2 = GH2 + H2_CH4   + H2_MONO + H2_ACET +
     &               H2_CH3OH + H2_ISOP - H2_OH   - H2_O1D
    
         ! For HD, do the same, scaling the BVOC terms by the 
         ! photochemical enrichement term (DPHOTO=162 permil)
         ! For methane, this is already done.
         GHD = GHD + HD_CH4 + 
     &         ( (H2_MONO + H2_ACET + H2_CH3OH + H2_ISOP ) * DPHOTO )
     &         - HD_OH - HD_O1D
                         
         ! Convert net H2 from [molec H2/cm3] to [kg] and store in STT
         STT(I,J,L,IDTH2) = GH2 / STTH2
	 
	 ! Convert net HD from [molec H2/cm3] to [kg] and store in STT
         STT(I,J,L,IDTHD) = GHD / STTHD

         !==============================================================
         ! Archive ND10 diagnostics -- Production & Loss of H2
         !==============================================================
         IF ( ND10 > 0 .and. L <= LD10 ) THEN

            ! Loss of H2 by OH [molec CO/cm3/s]
            N             = 1
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_OH / DTCHEM )

            ! Production of H2 from Isoprene [molec CO/cm3/s]
            N             = 2
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_ISOP / DTCHEM )

            ! Production of H2 from CH4 [molec CO/cm3/s]
            N             = 3
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_CH4 / DTCHEM )

            ! Production of from CH3OH [molec CO/cm3/s]
            N             = 4
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_CH3OH / DTCHEM )

            ! Production of H2 from MONO [molec CO/cm3/s]
            N             = 5
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_MONO / DTCHEM )

            ! Production of H2 from ACET [molec CO/cm3/s]
            N             = 6
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_ACET / DTCHEM )
	    
	    ! Loss of H2 by O1D in the stratosphere [molec CO/cm3/s]
	    N             = 7
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( H2_O1D / DTCHEM )
	    
            !Loss of HD by OH [molec CO/cm3/s]
            N             = 8
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( HD_OH / DTCHEM )

            ! Production of HD from Isoprene [molec CO/cm3/s]
            N             = 9
            AD10(I,J,L,N) = AD10(I,J,L,N)+
     &                         ((H2_ISOP *DPHOTO) / DTCHEM)

            ! Production of HD from CH4 [molec CO/cm3/s]
            N             = 10
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( HD_CH4 / DTCHEM ) 

            ! Production of HD from CH3OH [molec CO/cm3/s]
            N             = 11
            AD10(I,J,L,N) = AD10(I,J,L,N)+
     &                        ((H2_CH3OH*DPHOTO) / DTCHEM) 

            ! Production of HD from MONO [molec CO/cm3/s]
            N             = 12
            AD10(I,J,L,N) = AD10(I,J,L,N) +
     &                        ((H2_MONO * DPHOTO) / DTCHEM)

            ! Production of HD from ACET [molec CO/cm3/s]
            N             = 13
            AD10(I,J,L,N) = AD10(I,J,L,N) +
     &                          ((H2_ACET * DPHOTO) / DTCHEM)  
	    
	    ! Loss of HD by O1D in the stratosphere [molec CO/cm3/s]
            N             = 14
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( HD_O1D / DTCHEM )
	    
	    ! Alpha (ratio of OH k rates kHD/kH2)
            N             = 15
            AD10(I,J,L,N) = AD10(I,J,L,N) + ( HD_RATE /H2_RATE )
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_H2_HD

!-----------------------------------------------------------------------------

      SUBROUTINE READ_OCEAN_H2( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_OCEAN_H2 reads in oceanic H2 emissions from nitrogen 
!  fixation. (hup, lyj, phs, bmy, 9/18/07)
!
!  Ocean H2 emissions are based on the N2 oceanic fixation rates
!  determined by Curtis Deutsch (University of Washington) by
!  assimilating observed nutrient distributions in the oceans:
!  "Spatial coupling of nitrogen inputs and losses in the ocean",
!  Deutsch et al., Nature 445, 163-167 (2007).
!  
!  The oceanic N2 fixation rates are read in and then scaled to
!  obtain a total ocean H2 source of 6 TgH2/yr. This source is
!  assumed to be constant and does not vary annually.
!
!  Arguments as Input
!  =========================================================================== 
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments 
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      TYPE (XPLEX)              :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)              :: TEMP(IIPAR,JJPAR)
      TYPE (XPLEX)              :: XTAU, YMID
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_OCEAN_H2 begins here!
      !=================================================================
 
      ! Name of file with annual ocean H2 emissions
      FILENAME = TRIM( DATA_DIR )                   // 
     &           'hydrogen_200704/Ocean_H2_annual.' // GET_NAME_EXT() // 
     &           '.'                                // GET_RES_EXT()

      ! Echo to stdout
      WRITE( 6, '(8x,a)' ) 'Reading ', TRIM( FILENAME )

      !  Initialize some variables      
      EMOCEAN   = 0d0
       
      ! TAU value at the beginning of this month for 1994
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      !=================================================================
      ! Read ocean emissions
      !=================================================================

      ! Initialize ARRAY
      ARRAY = 0e0
      
      ! Read data!
       CALL READ_BPCH2( FILENAME, 'H2FIX',    1,  
     &                 XTAU,      IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.FALSE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), EMOCEAN )
      
      ! Return to calling program
      END SUBROUTINE READ_OCEAN_H2

!------------------------------------------------------------------------------

      SUBROUTINE READ_H2YIELD( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_H2YIELD reads in the relative H2/CO yield from photochemical
!  production. This has been archived monthly (PH2/PCO using the PRODLOSS
!  diagnostic and turning H2 on as an active species) from a full chemistry
!  simulation at 4x5, v7-03-03, year 2001, GEOS-3 met fields.
!  (lyj, hup, phs, bmy, 9/18/07)
!
!  Arguments as Input
!  =========================================================================== 
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE GRID_MOD,      ONLY : GET_YMID

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments 
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      TYPE (XPLEX)              :: ARRAY(IGLOB,JGLOB,LGLOB)
      TYPE (XPLEX)              :: TEMP(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: XTAU!, YMID
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_H2YIELD begins here!
      !=================================================================
 
      ! File with H2/CO yields
      FILENAME = TRIM( DATA_DIR )             //
     &           'hydrogen_200704/H2COyield.' // GET_NAME_EXT() // 
     &           '.'                          // GET_RES_EXT()

      ! Echo to stdout
      WRITE( 6, '(a)' ) 'READING ', TRIM( FILENAME )

      ! TAU value at the beginning of this month for 1985
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      !=================================================================
      ! Read Monthly H2/CO relative yields
      !=================================================================

      ! Initialize ARRAY and Yield
      ARRAY = 0e0
      
      ! Read data!
       CALL READ_BPCH2( FILENAME, 'PORL-L=$',    1,  
     &                  XTAU,      IGLOB,        JGLOB,      
     &                  LGLOB,     ARRAY )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_3D( ARRAY, H2CO_YIELD )
      
      ! Return to calling program
      END SUBROUTINE READ_H2YIELD

!------------------------------------------------------------------------------

      SUBROUTINE INIT_H2_HD
!
!******************************************************************************
!  Subroutine INIT_H2_HD allocates memory to module arrays.
!  (lyj, hyp, phs, bmy, 9/18/07)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS, I, J, IJLOOP
      
      !=================================================================
      ! INIT_H2_HD begins here!
      !=================================================================

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

      ! Allocate array for H2 from ocean n2 fixation
      ALLOCATE( EMOCEAN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMOCEAN' )         
      EMOCEAN = 0d0 
      
      ! Allocate array for H2 yield from photoch. production
      ALLOCATE( H2CO_YIELD( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2CO_YIELD' )         
      H2CO_YIELD = 0d0 
      
      ! Return to calling program
      END SUBROUTINE INIT_H2_HD

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_H2_HD
!
!******************************************************************************
!  Subroutine CLEANUP_H2_HD deallocates memory from previously
!  allocated module arrays. (lyj, hup, phs, bmy, 9/18/07)
!
!  NOTES:
!******************************************************************************
!
      IF ( ALLOCATED( SUMISOPCO  ) ) DEALLOCATE( SUMISOPCO  )
      IF ( ALLOCATED( SUMMONOCO  ) ) DEALLOCATE( SUMMONOCO  )
      IF ( ALLOCATED( SUMCH3OHCO ) ) DEALLOCATE( SUMCH3OHCO )
      IF ( ALLOCATED( SUMACETCO  ) ) DEALLOCATE( SUMACETCO  )
      IF ( ALLOCATED( EMOCEAN    ) ) DEALLOCATE( EMOCEAN    )

      ! Return to calling program
      END SUBROUTINE CLEANUP_H2_HD

!------------------------------------------------------------------------------

      ! End of module
      END MODULE H2_HD_MOD
