! $Id: emissions_adj_mod.f,v 1.7 2012/03/04 18:37:57 daven Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION: Module EMISSIONS\_MOD is used to call the proper emissions 
!  subroutines for the various GEOS-CHEM simulations. (bmy, 2/11/03, 2/14/08)
!\\
!\\
! !INTERFACE: 
!
      MODULE EMISSIONS_ADJ_MOD
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
      PUBLIC :: DO_EMISSIONS_ADJ
!
! !REVISION HISTORY:
!  (1 ) Now references DEBUG_MSG from "error_mod.f"
!  (2 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (3 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/2/04)
!  (4 ) Now references "seasalt_mod.f" (rjp, bmy, bec, 4/20/04)
!  (5 ) Now references "logical_mod" & "tracer_mod.f" (bmy, 7/20/04)
!  (6 ) Now references "epa_nei_mod.f" and "time_mod.f" (bmy, 11/5/04)
!  (7 ) Now references "emissions_mod.f" (bmy, 12/7/04)
!  (8 ) Now calls EMISSSULFATE if LCRYST=T.  Also read EPA/NEI emissions for 
!        the offline aerosol simulation. (bmy, 1/11/05)
!  (9 ) Remove code for the obsolete CO-OH param simulation (bmy, 6/24/05)
!  (10) Now references "co2_mod.f" (pns, bmy, 7/25/05)
!  (11) Now references "emep_mod.f" (bdf, bmy, 10/1/05)
!  (12) Now references "gfed2_biomass_mod.f" (bmy, 3/30/06)
!  (13) Now references "bravo_mod.f" (rjp, kfb, bmy, 6/26/06)
!  (14) Now references "edgar_mod.f" (avd, bmy, 7/6/06)
!  (15) Now references "streets_anthro_mod.f" (yxw, bmy, 8/18/06)
!  (16) Now references "h2_hd_mod.f" (lyj, phs, 9/18/07)
!  (17) Now calls EMISSDR for tagged CO simulation (jaf, mak, bmy, 2/14/08)
!  (18) Now references "cac_anthro_mod.f" (amv, phs, 03/11/08)
!  (19) Now references "vistas_anthro_mod.f" (amv, 12/02/08)
!  (20) Bug fixe : add specific calls for Streets for the grid 0.5x0.666.
!        (dan, ccc, 3/11/09)
!  (21) Add support for CH4 (kjw, dkh, 02/12/12, adj32_023) 
!EOP


      ! add for adjoint  (dkh, 06/02/08) 
      TYPE (XPLEX), ALLOCATABLE  :: BURNEMIS_orig(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: BIOFUEL_orig(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: EMISRR_orig(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: EMISRRB_orig(:,:,:)


      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------



      SUBROUTINE DO_EMISSIONS_ADJ
!
!******************************************************************************
!  Subroutine ADJ_DO_EMISSIONS is the driver routine which calls the appropriate
!  adjoint emissions subroutine for the various GEOS-CHEM simulations. 
!  Currently only supported for NSRCS = 3
!
!  NOTES
!  ( 1) Continued updates to v8 and CO simulation (mak, 7/1/09)
!  ( 2) The approach here is that we read emissions in the adj the same way
!       we do in the forward. I think we plan to replace that with storage
!       in EMS_orig, but until then, we'll recompute/reread. (mak,7/1/09)
!  (3 ) Now check that adjoint emissions ID #'s defined before calling
!       fullchem adjoint emissions routines (dkh, 11/11/09) 
!  (4 ) Now call EMISSCO2_ADJ. (dkh, 05/06/10) 
!  (5 ) Now include dust emissions adjoint (xxu, dkh, 01/09/12, adj32_011) 
!******************************************************************************
!
      ! References to F90 modules
      ! these two not yet ready (mak, 7/1/09)
      USE ADJ_ARRAYS_MOD,         ONLY : IS_CARB_EMS_ADJ
      USE ADJ_ARRAYS_MOD,         ONLY : IS_SULF_EMS_ADJ
      USE ADJ_ARRAYS_MOD,         ONLY : IS_DUST_EMS_ADJ
      USE CARBON_ADJ_MOD,         ONLY : EMISSCARBON_ADJ      
      USE DUST_ADJ_MOD,           ONLY : EMISSDUST_ADJ
      USE ERROR_MOD,              ONLY : ERROR_STOP, DEBUG_MSG
      USE SULFATE_ADJ_MOD,        ONLY : EMISSSULFATE_ADJ

      ! from EMISSIONS_MOD (mak, 7/1/09)
      USE BIOMASS_MOD,            ONLY : NBIOMAX
      USE BIOMASS_MOD,            ONLY : COMPUTE_BIOMASS_EMISSIONS
      USE ARCTAS_SHIP_EMISS_MOD,  ONLY : EMISS_ARCTAS_SHIP
      USE BRAVO_MOD,              ONLY : EMISS_BRAVO
      USE C2H6_MOD,               ONLY : EMISSC2H6
      USE CAC_ANTHRO_MOD,         ONLY : EMISS_CAC_ANTHRO
      USE CAC_ANTHRO_MOD,         ONLY : EMISS_CAC_ANTHRO_05x0666
      USE CARBON_MOD,             ONLY : EMISSCARBON
      USE CH3I_MOD,               ONLY : EMISSCH3I
      USE CO2_ADJ_MOD,            ONLY : EMISSCO2_ADJ
      USE DUST_MOD,               ONLY : EMISSDUST
      USE EDGAR_MOD,              ONLY : EMISS_EDGAR
      USE EMEP_MOD,               ONLY : EMISS_EMEP
      USE EMEP_MOD,               ONLY : EMISS_EMEP_05x0666
      USE EPA_NEI_MOD,            ONLY : EMISS_EPA_NEI
      USE GLOBAL_CH4_MOD,         ONLY : EMISSCH4
      USE GLOBAL_CH4_ADJ_MOD,     ONLY : EMISSCH4_ADJ
      USE H2_HD_MOD,              ONLY : EMISS_H2_HD
      USE HCN_CH3CN_MOD,          ONLY : EMISS_HCN_CH3CN
      USE Kr85_MOD,               ONLY : EMISSKr85
      USE LOGICAL_MOD                
      USE RnPbBe_MOD,             ONLY : EMISSRnPbBe
      USE SEASALT_MOD,            ONLY : EMISSSEASALT
      USE STREETS_ANTHRO_MOD,     ONLY : EMISS_STREETS_ANTHRO
      USE STREETS_ANTHRO_MOD,     ONLY : EMISS_STREETS_ANTHRO_05x0666
      USE NEI2005_ANTHRO_MOD,     ONLY : EMISS_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD,     ONLY : EMISS_NEI2005_ANTHRO_05x0666
      USE SULFATE_MOD,            ONLY : EMISSSULFATE 
      USE TIME_MOD,               ONLY : GET_MONTH,       GET_YEAR
      USE TIME_MOD,               ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE TRACER_MOD                 
      USE TAGGED_CO_ADJ_MOD,      ONLY : EMISS_TAGGED_CO_ADJ
      USE VISTAS_ANTHRO_MOD,      ONLY : EMISS_VISTAS_ANTHRO

      USE AIRCRAFT_ADJ_MOD,       ONLY : EMISS_AIRCRAFT_ADJ ! jkoo
      USE AIRCRAFT_MOD,           ONLY : LAIRCRAFT
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"                 ! FSCALYR

      INTEGER                     :: MONTH, YEAR

      !=================================================================
      ! DO_EMISSIONS_ADJ begins here!
      !=================================================================

      ! Get year and month
      MONTH = GET_MONTH()

      ! check if emissions year differs from met field year
      IF ( FSCALYR < 0 ) THEN
         YEAR = GET_YEAR()
      ELSE
         YEAR = FSCALYR
      ENDIF

      ! Get biomass burning emissions for use below
      IF ( LBIOMASS ) THEN
         CALL COMPUTE_BIOMASS_EMISSIONS( GET_YEAR(), MONTH )
      ENDIF

      IF ( ITS_A_FULLCHEM_SIM() ) THEN 

            ! haven't made these routines yet, but this is where they would go...
            !IF ( LSSALT ) CALL ADJ_EMISSSEASALT
  
            ! Add support for dust adjoint (xxu, dkh, 01/09/12, adj32_011) 
            IF ( LDUST  .and. IS_DUST_EMS_ADJ ) CALL EMISSDUST_ADJ

            ! Adjoint of carbon emissions 
            IF ( LCARB .and. IS_CARB_EMS_ADJ ) CALL EMISSCARBON_ADJ
 
            ! Adjoint of sulfate emissions (dkh, 11/04/09) 
            IF ( LSULF .and. IS_SULF_EMS_ADJ ) CALL EMISSSULFATE_ADJ

            ! Adjoint of gas-phase emissions is in setemis_adj.f 
 
            ! Adjoint of Aircraft emissions (jkoo,3/9/10)
            IF (  LAIRCRAFT ) THEN
                  CALL EMISS_AIRCRAFT_ADJ
            ENDIF

 
      ! (yhmao, dkh, 01/13/12, adj32_013) 
      ELSE IF (ITS_AN_AEROSOL_SIM()) THEN 

            IF ( LCARB .and. IS_CARB_EMS_ADJ ) CALL EMISSCARBON_ADJ 

            IF ( LDUST  .and. IS_DUST_EMS_ADJ ) CALL EMISSDUST_ADJ

      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN

         !--------------------
         ! Tagged CO
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (pdk, phs, 3/17/09)
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Read CAC emissions
         ! Now support nested (zhej, dkh, 01/16/12, adj32_015) 
         IF ( LCAC .and. ITS_A_NEW_YEAR() ) THEN
#if   defined( GRID05x0666 )      
            CALL EMISS_CAC_ANTHRO_05x0666
#else 
            CALL EMISS_CAC_ANTHRO
#endif
         ENDIF

         ! Read EDGAR emissions once per month
!----------------
! prior to 3/11/08
!         IF ( LEDGAR .and. ITS_A_NEW_MONTH() ) THEN
!----------------
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR( YEAR, MONTH )
         ENDIF

         ! Read EPA (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) CALL EMISS_EPA_NEI

         ! Now support nested (zhej, dkh, 01/16/12, adj32_015) 
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if    defined( GRID05x0666 )
            CALL EMISS_EPA_NEI                  ! Use NEI99 biofuel, nested
            CALL EMISS_NEI2005_ANTHRO_05x0666   ! Use NEI05 anthro,  global
#else
            CALL EMISS_EPA_NEI                  ! Use NEI99 biofuel, nested
            CALL EMISS_NEI2005_ANTHRO           ! Use NEI05 anthro,  global
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR()  ) CALL EMISS_BRAVO

         ! Read EMEP (Europe) emissions once per year (adj32_015)
         IF ( LEMEP  .and. ITS_A_NEW_MONTH()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
#else
            CALL EMISS_EMEP
#endif
         ENDIF

         ! Now call EMISSDR for Tagged CO fossil fuel emissions, 
         ! so that we get the same emissions for Tagged CO as 
         ! we do for the full-chemistry (jaf, mak, bmy, 2/14/08)
         CALL EMISSDR 
         
         ! Emit tagged CO
         CALL EMISS_TAGGED_CO_ADJ
   

      ! Add support for CH4 simulation (dkh, 02/12/12, adj32_023) 
      ELSE IF ( ITS_A_CH4_SIM() ) THEN

         CALL EMISSCH4_ADJ

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

          ! don't have anything for tag OX emissions adjoint yet
          print*, ' warning: emissions adj for tagged OX not supported'

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         ! Emit CO2
         CALL EMISSCO2_ADJ

      ELSE 
      !============= we could add other simulation mode later !cs

      !              ....................

         CALL ERROR_STOP(' Other values of NSRCX not supported yet',
     &        ' ADJ_DO_EMISSIONS')
         
      ENDIF 


      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG ( '### ADJ_DO_EMISSIONS: a EMISSIONS' )

      ! Return to calling program
      END SUBROUTINE DO_EMISSIONS_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE INIT_EMISSIONS_ADJ
!     
!******************************************************************************
!  Subroutine INIT_EMISSIONS initializes all module arrays (dkh, 06/01/06)  
!
!  NOTES:
!  ( 1) Replace NBIOTRCE with NBIOMAX in v8 update (I think that's equivalent)
!       (mak, 7/1/09)
! 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR, ERROR_STOP
      USE BIOMASS_MOD, ONLY : NBIOMAX
      USE BIOFUEL_MOD, ONLY : NBFTRACE

#     include "CMN_SIZE" ! Size parameters
#     include "comode.h" ! ITLOOP


      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_EMISSIONS     begins here!
      !=================================================================

      ! Return if we already allocated arrays
      IF ( IS_INIT ) RETURN

      ALLOCATE( BIOFUEL_orig( NBFTRACE, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOFUEL_orig' )
      BIOFUEL_orig = 0d0

      ALLOCATE( BURNEMIS_orig( NBIOMAX, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BURNEMIS_orig' )
      BURNEMIS_orig = 0d0

      ! fix (dkh, 05/04/09) 
      !ALLOCATE( EMISRR_orig( IIPAR, JJPAR, 2:NEMPARA+NEMPARB ), STAT=AS)
      ALLOCATE( EMISRR_orig( IIPAR, JJPAR, 1:NEMPARA+NEMPARB ), STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISRR_orig' )
      EMISRR_orig = 0d0

      ! fix (dkh, 05/04/09) 
      !ALLOCATE(EMISRRB_orig( IIPAR, JJPAR, 2:NEMPARA+NEMPARB ), STAT=AS)
      ALLOCATE(EMISRRB_orig( IIPAR, JJPAR, 1:NEMPARA+NEMPARB ), STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISRRB_orig' )
      EMISRRB_orig = 0d0

      ! Reset IS_INIT
      IS_INIT = .TRUE. 
      
      ! Return to calling progam 
      END SUBROUTINE INIT_EMISSIONS_ADJ
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_EMISSIONS
!     
!******************************************************************************
!  Subroutine CLEANUP_EMISSIONS deallocates all module arrays  
!   (dkh, 06/01/06)  
!
!  NOTES:
! 
!******************************************************************************
!     
      !=================================================================
      ! CLEANUP_EMISSIONS begins here!
      !=================================================================
      IF ( ALLOCATED( BIOFUEL_orig   ) ) DEALLOCATE( BIOFUEL_orig   )
      IF ( ALLOCATED( BURNEMIS_orig  ) ) DEALLOCATE( BURNEMIS_orig  )
      IF ( ALLOCATED( EMISRR_orig    ) ) DEALLOCATE( EMISRR_orig    )
      IF ( ALLOCATED( EMISRRB_orig   ) ) DEALLOCATE( EMISRRB_orig   )
      
      ! Return to calling program
      END SUBROUTINE CLEANUP_EMISSIONS


      END MODULE EMISSIONS_ADJ_MOD
!EOC
