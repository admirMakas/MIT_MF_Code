!$Id: co2_mod.f,v 1.5 2012/03/01 22:00:26 daven Exp $
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: co2_mod
!     
! !DESCRIPTION: Module CO2\_MOD contains variables and routines used for the 
!  CO2 simulation.  A tagged CO2 simulation capability has now been added.
!\\   
!\\   
!  References:
!
!  \begin{itemize}
!  \item Andres, R.J, G. Marland, I. Fung, and E. Matthews, \emph{A 1x1 
!        distribution of carbon dioxide emissions from fossil fuel 
!        consumption and cement manufacture}, \underline{Glob. Biogeochem. 
!        Cycles}, \textbf{10}, 419-429, 1996.
!  \item Corbett and Koehler (2003) \emph{Updated emissions from ocean 
!        shipping}, \underline{J. Geophys. Res.}, \textbf{108}, D20, 4650.
!  \item Corbett and Koehler (2004) \emph{Considering alternative input 
!        parameters in an activity-based ship fuel consumption and emissions 
!        model: Reply ...} \underline{J. Geophys. Res.}, D23303.
!  \item Endresen et al. (2007) \emph{A historical reconstruction of ships 
!        fuel consumption and emissions}, \underline{J. Geophys. Res.} 
!        \textbf{112}, D12301.
!  \item Kim et al. (2005) \emph{System for assessing Aviation's Global 
!        Emissions (SAGE) Version 1.5 global Aviation Emissions Inventories 
!        for 2000-2004}
!  \item Kim et al. (2007) \emph{System for assessing Aviation's Global 
!        Emissions (SAGE) Part 1: Model description and inventory results} 
!  \item LeQuere et al. (2009) \emph{Trends in the sources and sinks of carbon 
!        dioxide}, \underline{Nature Geoscience}, doi:10.1038/ngeo689.
!  \item Olsen and Randerson (2004), \emph{Differences between surface and 
!        column atmospheric CO2 and implications for carbon cycle research}, 
!        \underline{J. Geophys. Res.}, \textbf{109}, D02301,
!  \item Potter et al. (1993), \emph{Terrestrial Ecosystem Production: 
!        A process model based on global satellite and surface data}, 
!        \underline{Glob. Biogeochem. Cycles}, \textbf{7}(4), 811-841.
!  \item Randerson, J.T, M.V. Thompson, T.J.Conway, I.Y. Fung, and C.B. Field,
!        \emph{The contribution of terrestrial sources and sinks to trends 
!        in the seasonal cycle of atmospheric carbon dioxide}, 
!        \underline{Glob. Biogeochem. Cycles},\textbf{11}, 535-560, 1997.
!  \item Suntharalingam et al. (2005) \emph{Infulence of reduced carbon 
!        emissions and oxidation on the distribution of atmospheric CO2: 
!        Implications for inversion analysis}, BGC, 19, GB4003.
!  \item Takahashi, T, R. Feely, R. Weiss, R. Wanninkof, D. Chipman, 
!        S. Sutherland, and T. Takahashi (1997), \emph{Global air-sea flux 
!        of CO2: An estimate based on measurements of sea-air pCO2 difference},
!        \underline{Proceedings of the National Academy of Sciences}, 
!        \textbf{94}, 8292-8299.
!  \item Takahashi, T, et al. (2009), \emph{Climatological mean and decadal 
!        change in surface ocean pCO2, and net sea-air CO2 flux over the 
!        global oceans}, \textbf{Deep-Sea Research II}, 
!        doi:10.1016/jdsr2/2008.12.009.
!  \item Yevich, R. and J. A. Logan, \emph{An assesment of biofuel use and 
!        burning of agricultural waste in the developing world}, 
!        \underline{Glob. Biogeochem. Cycles}, \textbf{17}, 1095, 
!        doi:10.1029/2002GB001952, 2003.
!  \item Sausen, R. and Schumann, U. "Estimates of the Climate Response to
!        Aircraft CO2 and NOx Emissions Scenarios", Climate Change, 
!        44: 27-58, 2000
!  \item Wilkersen, J.T. et al. \emph{Analysis of emission data from global 
!        commercial Aviation: 2004 and 2006}, \underline{Atmos. chem. Phys. 
!        Disc.}, \textbf{10}, 2945-2983, 2010.
!  \end{itemize}
!
! !INTERFACE: 
!     
      MODULE CO2_MOD
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
      PUBLIC  :: CLEANUP_CO2
      PUBLIC  :: EMISSCO2

! adj_group:  make these public as well (dkh, 03/07/11) 
      PUBLIC :: READ_BBIO_DIURNALCYCLE
      PUBLIC :: READ_BBIO_DAILYAVERAGE
      PUBLIC :: READ_FOSSILCO2
      PUBLIC :: READ_OCEANCO2
      PUBLIC :: READ_SHIPCO2_EDGAR
      PUBLIC :: READ_SHIPCO2_ICOADS
      PUBLIC :: READ_AVIATION_CO2
      PUBLIC :: READ_CHEMCO2
      PUBLIC :: CHEM_SURF
      PUBLIC :: CHEMCO2
      PUBLIC :: EMFOSSCO2
      PUBLIC :: EMOCCO2
      PUBLIC :: EMBIOCO2
      PUBLIC :: EMBIOFUELCO2
      PUBLIC :: EMBIONETCO2
      PUBLIC :: EMSHIPCO2
      PUBLIC :: EMPLANECO2
      PUBLIC :: EMIS_SUB
      PUBLIC :: XNUMOL_CO2


! adj group: some of these are now public (dkh, 03/07/11) 
! !PRIVATE MEMBER FUNCTIONS:
!     
      !PRIVATE :: READ_CHEMCO2   
      !PRIVATE :: READ_FOSSILCO2
      !PRIVATE :: CHEM_SURF  
      PRIVATE :: AVIATION_DOM_CORR  
      !PRIVATE :: READ_OCEANCO2
      PRIVATE :: READ_ANNUAL_BIOFUELCO2
      !PRIVATE :: READ_SHIPCO2_EDGAR
      !PRIVATE :: READ_SHIPCO2_ICOADS
      !PRIVATE :: READ_AVIATION_CO2
      PRIVATE :: READ_ANNUAL_BIONET_CO2 
      !PRIVATE :: READ_BBIO_DAILYAVERAGE
      !PRIVATE :: READ_BBIO_DIURNALCYCLE
      PRIVATE :: TOTAL_BIOMASS_TG
      PRIVATE :: DEF_BIOSPH_CO2_REGIONS_F
      PRIVATE :: DEF_OCEAN_CO2_REGIONS_F
      PRIVATE :: DEF_FOSSIL_CO2_REGIONS_F
      PRIVATE :: INIT_CO2 
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%   BUYER BEWARE! Tagged CO2 tracers only work for 2 x 2.5 grid!   %%%
!  %%%   Someone will have to make this more general later on...        %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                             .
! !REVISION HISTORY:
!  16 Aug 2005 - P. Suntharalingam - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now references biomass_mod.f (bmy, 9/27/06)
!  (3 ) Tagged CO2 capability developed (dbj)
!  (4 ) Implemented monthly and annual fossil fuel inventories 
!        (R.Nassar 2009-03-10)
!  (5 ) Implemented CO2 emissions from shipping and aviation (R.Nassar 2010)
!  (6 ) Implemented monthly CO2 chemical production and surface correction 
!        (R.Nassar 2010)   
!  25 Feb 2011 - R. Nassar  - Now read updated CDIAC CO2 emissions data
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, ALLOCATABLE :: FOSSIL_REGION(:,:)
      INTEGER, ALLOCATABLE :: BIOSPH_REGION(:,:)
      INTEGER, ALLOCATABLE :: OCEAN_REGION(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMFOSSCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMOCCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMBIOCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMBIOBRNCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMBIOFUELCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMBIONETCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMSHIPCO2(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMPLANECO2(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CHEMCO2(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EMIS_SUB(:,:)
!
! !DEFINED PARAMETERS:
!
      ! FMOL_CO2     - kg CO2 / mole CO2 
      TYPE (XPLEX),  PARAMETER   :: FMOL_CO2   = xplex(44d-3,0d0)

      ! XNUMOL_CO2   - molecules CO2 / kg CO2 
      TYPE (XPLEX),PARAMETER::XNUMOL_CO2=xplex(6.022d+23/FMOL_CO2%r,0d0)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissco2
!
! !DESCRIPTION: Subroutine EMISSCO2 is the driver routine for CO2 emissions. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSCO2
!
! !USES:
!
      USE BIOMASS_MOD,  ONLY : BIOMASS
      USE DIAG04_MOD,   ONLY : AD04, ND04
      USE DIAG04_MOD,   ONLY : AD04_plane,    AD04_chem
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_DAY,       GET_DAY_OF_YEAR
      USE TIME_MOD,     ONLY : GET_HOUR,      GET_MONTH
      USE TIME_MOD,     ONLY : GET_YEAR,      GET_TS_CHEM, GET_TS_EMIS 
      USE TIME_MOD,     ONLY : ITS_A_NEW_DAY, ITS_A_NEW_MONTH
      USE TRACER_MOD,   ONLY : N_TRACERS,     STT

      ! adj_group: IDBCO2 still in biomass_mod.f (dkh, 03/07/11) 
      !USE TRACERID_MOD, ONLY : IDBCO2
      USE BIOMASS_MOD,  ONLY : IDBCO2

      USE LOGICAL_MOD,  ONLY : LGENFF,      LANNFF,   LMONFF, LSTREETS
      USE LOGICAL_MOD,  ONLY : LSEASBB,     LGFED2BB, L8DAYBB, LBIOFUEL
      USE LOGICAL_MOD,  ONLY : LBIODAILY,   LBIODIURNAL
      USE LOGICAL_MOD,  ONLY : LBIONETORIG, LBIONETCLIM
      USE LOGICAL_MOD,  ONLY : LOCN1997,    LOCN2009ANN, LOCN2009MON
      USE LOGICAL_MOD,  ONLY : LSHIPEDG,    LSHIPICO,    LPLANE
      USE LOGICAL_MOD,  ONLY : LBIOSPHTAG,  LFOSSILTAG,  LFFBKGRD
      USE LOGICAL_MOD,  ONLY : LSHIPTAG,    LPLANETAG
      USE LOGICAL_MOD,  ONLY : LSHIPSCALE,  LPLANESCALE
      USE LOGICAL_MOD,  ONLY : LCHEMCO2

      ! adj_group (dkh, 03/07/11) 
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2ff,  IDADJ_ECO2ocn
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2bal, IDADJ_ECO2bb
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2bf,  IDADJ_ECO2nte
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2shp, IDADJ_ECO2pln
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2che, IDADJ_ECO2sur
      USE ADJ_ARRAYS_MOD,  ONLY : GET_SCALE_GROUP
      USE ADJ_ARRAYS_MOD,  ONLY : EMS_SF
      
      ! dkh debug 
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD


      
#     include "CMN_SIZE"    ! Size parameters
!
! !REMARKS:
!  The initial condition for CO2 has to be at least 50 ppm or higher or else
!  the balanced biosphere fluxes will make STT negative. (pns, bmy, 8/16/05)
! 
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) We now get CO2 biomass emissions from biomass_mod.f.  This allows us 
!        to use either GFED2 or default Duncan et al biomass emissions. 
!        (bmy, 9/27/06)
!  (3 ) Tagged tracer capability added. This requires the editable region 
!        files Regions_land.dat and Regions_ocean.dat in the run directory
!        (rnassar,dbj, 2009)
!  (4 ) New tracers for emissions from international and domestic shipping,
!        international and domestic aviation, and the chemical CO2 source
!        from the oxidation of CO, CH4, and other organics (rnassar,dbj, 2009)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRST  = .TRUE.

      ! Local variables
      INTEGER                :: I,     IJLOOP, J,    L,     N,   NN
      INTEGER                :: DAY,   DOY,    HOUR, MONTH, YEAR   
      TYPE (XPLEX)                 :: A_CM2, DTSRCE, E_CO2
      TYPE (XPLEX)                 :: biomass_sum, bionet_sum
      TYPE (XPLEX), SAVE           :: CHEMSRC(IIPAR,JJPAR,LLPAR)  ! dbj

      ! External functions     
      TYPE (XPLEX), EXTERNAL       :: BOXVL    ! dbj

      ! adj_group 
      INTEGER                :: M

      !=================================================================
      ! EMISSCO2 begins here!
      !=================================================================
      ! First-time initialization

      IF ( FIRST ) THEN 

         ! Allocate arrays and read annual-mean data
         CALL INIT_CO2

         ! Set up tagged regions for balanced biosphere & ocean
         IF ( LBIOSPHTAG ) THEN
            CALL DEF_BIOSPH_CO2_REGIONS_F( BIOSPH_REGION )
            CALL DEF_OCEAN_CO2_REGIONS_F( OCEAN_REGION )
         ENDIF

         ! Set up tagged regions for fossil fuel
         IF ( LFOSSILTAG ) THEN
            CALL DEF_FOSSIL_CO2_REGIONS_F( FOSSIL_REGION )
         ENDIF

         ! Set first-time flag to false
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read in monthly, daily or variable emissions fields
      !=================================================================      

      ! Emission timestep
!      DTSRCE = 60d0 * GET_TS_CHEM() !Line from the orginal code
      DTSRCE = 60d0 * GET_TS_EMIS()
      
      ! Time variables
      DAY    = GET_DAY()
      DOY    = GET_DAY_OF_YEAR()
      HOUR   = GET_HOUR()
      MONTH  = GET_MONTH()
      YEAR   = GET_YEAR()

      ! adj_group
      M      = GET_SCALE_GROUP()

!------------------------------------------------------------------------------
!      ! Read monthly-mean biomass burning emissions
!      IF ( LBIOBRNCO2 .and. ITS_A_NEW_MONTH() ) THEN 
!         CALL READ_MONTH_BIOBRN_CO2( MONTH, YEAR )
!      ENDIF
!  This requires a subroutine called READ_MONTH_BIOBRN_CO2 
!  GFEDv2 biomass burning emissions are a better choice !Ray Nassar
!
!------------------------------------------------------------------------------
!  At present, biomass burning emissions are dealt with in the following way:
!
!  1)  main.f calls do_emissions in emissions.f
!  2)  do_emissions calls compute_biomass_emissions in biomass_mod.f
!  3a) compute_biomass_emissions calls gfed2_compute_biomass 
!      in gfed2_biomass_mod.f
!               ** OR **
!  3b) compute_biomass_emissions calls gc_read_biomass_co2 in gc_biomass_mod.f
!------------------------------------------------------------------------------

!      ! Check if Balanced Biosphere emissions are required  
!      IF ( LBIOCO2 ) THEN  
!         ! If LUSECASANEP is TRUE ...
!         IF ( LUSECASANEP ) THEN

      !----------------------------------------------------------------
      ! Read in 3-hourly or daily balanced biosphere data
      !----------------------------------------------------------------
      IF ( LBIODIURNAL ) THEN

         write(*,*) '*** USING DIURNAL CASA NEP ***'

         ! ... then use 3-hourly NEP emissions for Bal Bio ...
         IF ( MOD( HOUR, 3 ) == 0 ) THEN
            CALL READ_BBIO_DIURNALCYCLE( MONTH, DAY, HOUR, DOY )
         ENDIF
         
      ELSEIF ( LBIODAILY ) THEN

         ! ... otherwise use constant daily emissions of NEP for Bal Bio
         IF ( ITS_A_NEW_DAY() ) THEN
            CALL READ_BBIO_DAILYAVERAGE( MONTH, DAY, DOY ) 
         ENDIF
         
      ENDIF

      !----------------------------------------------------------------
      ! Fluxes with "possible" monthly variability are called below
      ! In some cases the annual file is just called at the start 
      ! of the month
      !----------------------------------------------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN 
	
         ! Fossil fuel emissions
         IF ( LMONFF .OR. LANNFF .OR. LGENFF ) THEN
            CALL READ_FOSSILCO2 
         ENDIF
                               
         ! Oceanic exchange
         IF ( LOCN1997 .OR. LOCN2009ANN .OR. LOCN2009MON ) THEN
            CALL READ_OCEANCO2
         ENDIF

         ! Ship emissions from EDGAR
         IF ( LSHIPEDG ) CALL READ_SHIPCO2_EDGAR

         ! Ship emissions from ICOADS
         IF ( LSHIPICO ) CALL READ_SHIPCO2_ICOADS
 
         ! Aircraft CO2 emissions
         IF ( LPLANE   ) CALL READ_AVIATION_CO2
		
         ! Get chemical source   ! dbj
         IF ( LCHEMCO2 ) THEN
            CALL READ_CHEMCO2
            CALL CHEM_SURF
            CHEMSRC = CHEMCO2
         ENDIF
      ENDIF
	
      !=================================================================
      ! Process emissions and save diagnostics
      !=================================================================

      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, A_CM2, E_CO2 )  ! dbj
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

      ! Loop over longitudes
      DO I = 1, IIPAR

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from fossil fuel emissions 
         !-------------------------------------------
         IF ( LGENFF .or. LANNFF .or. LMONFF ) THEN

            ! Fossil fuel emissions of CO2 [molec/cm2/s]
            E_CO2          = EMFOSSCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2ff > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ff)
            ENDIF
            
            ! ND04 diag: Fossil Fuel CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,1) = AD04(I,J,1) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

	     ! Add to Tracer #2: Fossil CO2 [kg]
            IF ( N_TRACERS > 1 ) THEN
               STT(I,J,1,2)   = STT(I,J,1,2) + E_CO2

               ! Split Fossil Fuel CO2 into geographic regions
               IF ( LFOSSILTAG ) THEN
                  N            = FOSSIL_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO2
               ENDIF
            ENDIF

	     ! Add to Tracer #12: Background with Fossil CO2 [kg]
            IF ( LFFBKGRD ) THEN
               STT(I,J,1,12)   = STT(I,J,1,12) + E_CO2
            ENDIF

	   ! Note: To define the background as including fossil fuels here 
           ! (instead of during data processing) for tagged-CO2 runs which 
           ! will be used for estimating natural bio/ocean fluxes etc., and 
           ! accepting the FF inventory
           ! use STT(I,J,1,12) = STT(I,J,1,12) + E_CO2.
	   ! Shipping and aviation emissions can be included in a similar way.
	   
         ENDIF

         !-------------------------------------------
         ! #3: CO2 from ocean exchange
         !-------------------------------------------
         IF ( LOCN1997 .or. LOCN2009ANN .or. LOCN2009MON ) THEN

            ! Ocean CO2 emissions in [molec/cm2/s]
            E_CO2          = EMOCCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            ! Now move this here so ND04 is correct (dkh, 02/08/12, adj32_018) 
            IF ( LADJ .and. IDADJ_ECO2ocn > 0 ) THEN

               ! dkh debug
               IF ( I == IFD .and. J == JFD .and. LPRINTFD ) THEN
                  print*, ' ECO2onc fwd = ', E_CO2
                  print*, ' ECO2onc SF fwd = ', 
     &                     EMS_SF(I,J,M,IDADJ_ECO2ocn)
               ENDIF
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ocn)
            ENDIF

            ! ND04 diag: Ocean CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,2) = AD04(I,J,2) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2


            ! Add to Tracer #1: Total CO2 [kg]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #3: Ocean CO2 [kg]
            IF ( N_TRACERS > 2 ) THEN
               STT(I,J,1,3)   = STT(I,J,1,3) + E_CO2   

               ! Split ocean CO2 into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = OCEAN_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO2
               ENDIF
            ENDIF
         ENDIF

         !-------------------------------------------
         ! #4: CO2 from balanced biosphere emissions
         !-------------------------------------------
         IF ( LBIODAILY .OR. LBIODIURNAL ) THEN

            ! Balanced biosphere CO2 [molec/cm2/s]
            E_CO2         = EMBIOCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            ! Now move this here so ND04 is correct (dkh, 02/08/12, adj32_018) 
            IF ( LADJ .and. IDADJ_ECO2bal > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2bal)
            ENDIF

            ! ND04 diag: Bal Bio CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,3) = AD04(I,J,3) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg CO2]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 

            ! Add to Tracer #1 -- total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #4 -- Bal Bio CO2 [kg CO2]
            IF ( N_TRACERS > 3 ) THEN
               STT(I,J,1,4)   = STT(I,J,1,4) + E_CO2
                
               	! Split biospheric CO2 exchange into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = BIOSPH_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO2
               ENDIF
            ENDIF	
         ENDIF

         !-------------------------------------------
         ! #5: CO2 from biomass burning emissions
         !-------------------------------------------
         IF ( LSEASBB .OR. LGFED2BB .OR. L8DAYBB ) THEN 

            ! Biomass burning emissions [molec/cm2/s]
            E_CO2          = BIOMASS(I,J,IDBCO2)
            !E_CO2          = EMBIOBRNCO2(I,J) 
            !This was from older versions, see note above

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2bb > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2bb)
            ENDIF

            ! ND04 diag: Biomass burning CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,4) = AD04(I,J,4) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 

            ! Add to Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            
            ! Add to Tracer #5: Biomass burning CO2 [kg CO2]
            IF ( N_TRACERS > 4 ) THEN
               STT(I,J,1,5)   = STT(I,J,1,5) + E_CO2

               ! Split Bioburn CO2 into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = BIOSPH_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO2
               ENDIF
            ENDIF	
         ENDIF

         !-------------------------------------------
         ! #6: CO2 from biofuel emissions
         !-------------------------------------------
         IF ( LBIOFUEL ) THEN

           ! Biofuel CO2 emissions [molec/cm2/s]
            E_CO2          = EMBIOFUELCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2bf > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2bf)
            ENDIF

            ! ND04 diag: Biofuel CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,5) = AD04(I,J,5) + E_CO2
            ENDIF

            ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #6: Biofuel CO2 [kg CO2]
            IF (N_TRACERS > 5) THEN
               STT(I,J,1,6)   = STT(I,J,1,6) + E_CO2

               ! Split BF CO2 into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = BIOSPH_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO2
               ENDIF
            ENDIF
         ENDIF

         !-------------------------------------------
         ! #7: CO2 from net terrestrial exchange
         !-------------------------------------------
         IF ( LBIONETORIG .OR. LBIONETCLIM ) THEN

            ! CO2 from net terrestrial exchange [molec/cm2/s]
            E_CO2          = EMBIONETCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2nte > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2nte)
            ENDIF

            ! ND04 diag: net terrestrial exchange [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,6) = AD04(I,J,6) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #7: Net Terr exchange CO2 [kg]
            IF ( N_TRACERS > 6 ) THEN
               STT(I,J,1,7)   = STT(I,J,1,7) + E_CO2

               ! Split Net Terr Exch CO2 into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = BIOSPH_REGION(I,J)
                  STT(I,J,1,N) = STT(I,J,1,N) + E_CO2
               ENDIF
            ENDIF	
         ENDIF

         !-------------------------------------------
         ! #8: CO2 from ship emissions
         !-------------------------------------------
         IF ( LSHIPEDG .OR. LSHIPICO ) THEN

            ! Ship CO2 emissions [molec/cm2/s]
            E_CO2          = EMSHIPCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2shp > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2shp)
            ENDIF

            ! ND04 diag: Ship CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,7) = AD04(I,J,7) + E_CO2
            ENDIF

            ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #8: Ship CO2 [kg CO2]
            IF (N_TRACERS > 7) THEN
               STT(I,J,1,8)   = STT(I,J,1,8) + E_CO2
 
              	! Tagged tracer for global ship emissions
               IF ( LSHIPTAG ) THEN
                  STT(I,J,1,53) = STT(I,J,1,53) + E_CO2
               ENDIF
            ENDIF

	   ! Add to Tracer #12: Background with Fossil CO2 [kg]
           !-------------------------------------------
	   !IF ( LFFBKGRD ) THEN
           !   STT(I,J,1,12)   = STT(I,J,1,12) + E_CO2
	   !ENDIF
           !-------------------------------------------
	   ! Uncomment to include ship CO2 emissions in the background

         ENDIF
	   
         !-------------------------------------------
         ! #9: CO2 from aircraft emissions
         !-------------------------------------------
         IF ( LPLANE ) THEN
            DO L = 1, LLPAR

               ! Aircraft CO2 emissions (3-D) [molec/cm3/s]
               E_CO2          = EMPLANECO2(I,J,L)

               ! adj_group: apply scaling factors (dkh, 04/25/10) 
               IF ( LADJ .and. IDADJ_ECO2pln > 0 ) THEN
                  E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2pln)
               ENDIF

               ! ND04 diag: Aircraft CO2 [molec/cm3/s] 
               IF ( ND04 > 0 ) THEN
                  AD04_plane(I,J,L) = AD04_plane(I,J,L) + E_CO2
               ENDIF

               ! Convert E_CO2 from [molec CO2/cm3/s] to [kg]
               E_CO2    = E_CO2 * BOXVL(I,J,L) * DTSRCE / XNUMOL_CO2

               ! Add to Tracer #1: Total CO2 [kg CO2]
               STT(I,J,L,1)   = STT(I,J,L,1) + E_CO2

               ! Add to Tracer #9: Aircraft CO2 [kg CO2]
               IF ( N_TRACERS > 8 ) THEN
                  STT(I,J,L,9) = STT(I,J,L,9) + E_CO2
 
              	  ! Tagged tracer for global ship emissions
                  IF ( LPLANETAG ) THEN
                     STT(I,J,L,54) = STT(I,J,L,54) + E_CO2
                  ENDIF
               ENDIF
	
	      ! Add to Tracer #12: Background with Fossil CO2 [kg]
              !-------------------------------------------
	      !IF (LFFBKGRD) THEN
	      !   STT(I,J,L,12)   = STT(I,J,L,12) + E_CO2
	      !ENDIF
              !-------------------------------------------
	      ! Uncomment to include aviation CO2 emissions in the background

            ENDDO
         ENDIF
 
         !-------------------------------------------
         ! #10 CO2 production from CO oxidation
         !-------------------------------------------
         IF ( LCHEMCO2 ) THEN
            DO L = 1, LLPAR
         
               E_CO2 = CHEMSRC(I,J,L)

               ! adj_group: apply scaling factors (dkh, 04/25/10) 
               IF ( LADJ .and. IDADJ_ECO2che > 0 ) THEN
                  E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2che)
               ENDIF


               ! ND04 diag: CO2 chem source from CO loss (3-D) [molec/cm3/s] 
               IF ( ND04 > 0 ) THEN
                  AD04_chem(I,J,L) = AD04_chem(I,J,L) + E_CO2
               ENDIF

               !  Convert from [molec/cm3/s] to [kg]
               E_CO2  = E_CO2 * BOXVL(I,J,L) * DTSRCE / XNUMOL_CO2

               ! Add to Tracer #1: Total CO2 [kg CO2]
               STT(I,J,L,1)   = STT(I,J,L,1) + E_CO2

               ! Add to Tracer #10: Chemical Source of CO2 [kg CO2]
               IF (N_TRACERS > 9) THEN
                  STT(I,J,L,10)   = STT(I,J,L,10) + E_CO2
               ENDIF
		
            ENDDO
         ENDIF

         !-------------------------------------------
         ! #11 CO2 surface correction for CO oxidation
         !-------------------------------------------
         IF ( LCHEMCO2 ) THEN
            
            E_CO2 = EMIS_SUB(I,J) ! EMIS_SUB is positive, but is subtracted

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            ! Now move this here so ND04 is correct (dkh, 02/08/12, adj32_018) 
            IF ( LADJ .and. IDADJ_ECO2sur > 0 ) THEN
               E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2sur)
            ENDIF

            ! ND04 diag: CO2 chem source surface correction [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,10) = AD04(I,J,10) - E_CO2 ! SUBTRACT
            ENDIF

            ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Subtract from Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) - E_CO2

            ! Subtract from Tracer #11: Chem Source Surf Correction [kg CO2]
            IF ( N_TRACERS > 10 ) THEN
               STT(I,J,1,11)   = STT(I,J,1,11) - E_CO2
            ENDIF
		
         ENDIF

         !-------------------------------------------
         ! #12: Background CO2
         !-------------------------------------------
         ! Background CO2 without fossil fuels is obtained by setting 
         ! tracer 12 in the restart file, equal to tracer 1 at the start of 
         ! a run.  No sources or sinks (chemical or surface) act on this 
         ! tracer if LFFBKGRD == FALSE, it is simply the advected intial CO2.
	 ! To include fossil fuels in the background, see the first loop.

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE EMISSCO2 
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_chemco2
!
! !DESCRIPTION: Reads the chemical source of CO2 [molec/cm3/s] from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CHEMCO2
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE BPCH2_MOD,     ONLY : GET_MODELNAME, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,   READ_BPCH2
      USE TIME_MOD,      ONLY : GET_MONTH,  GET_YEAR

#     include "CMN_SIZE"      ! Size parameters
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: YEAR, MONTH   
      TYPE (XPLEX)                 :: ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                 :: TAU
      CHARACTER(LEN=4)       :: YEAR_STR
      CHARACTER(LEN=255)     :: FILENAME

      YEAR = GET_YEAR()
      if (YEAR < 2000) then
        YEAR = 2000
      endif
      MONTH  = GET_MONTH()
      TAU = GET_TAU0( MONTH, 1, YEAR )

	WRITE( YEAR_STR, '(i4)' ) YEAR 

      FILENAME = TRIM( DATA_DIR )                     // 
     &           'CO2_201003/ChemSrc/CO2_prod_rates_' // 
     &           TRIM( YEAR_STR )        // '.'       // 
     &           TRIM( GET_MODELNAME() ) // '.'       // 
     &           GET_RES_EXT()

      ARRAY = 0.0e0

      Print*,'Reading CO2 Chem. production rates from file = ',
     &   trim(filename)

      !=================================================================
      ! Read chemical source of CO2 [molec/cm3/s]
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 4, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 LLPAR,     ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      ! assume that incoming data is on the same vertical grid
      CHEMCO2 = (ARRAY)

      END SUBROUTINE READ_CHEMCO2   
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_fossilco2
!
! !DESCRIPTION: Subroutine READ\_FOSSILCO2 reads in fossil fuel CO2 
!  emissions from a bpch file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_FOSSILCO2
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,      ONLY : GET_YEAR, GET_MONTH
      USE LOGICAL_MOD,   ONLY : LGENFF, LANNFF, LMONFF, LCHEMCO2, LPLANE

#     include "CMN_SIZE"      ! Size parameters
!
! !REMARKS:
!  Original data provided by Robert Andres (CDIAC), personal communication
!                                                                             .
!  If GENFF=T, then annual data for 1995 are read (but tau is for 1985)
!  If ANNFF=T, then annual data for a given year (1985-2006) are read
!  If MONFF=T, then annual data for a given month (198501-200612) are read
!                                                                             .
!  ANNFF and MONFF for 2007-2009 were developed based on scaling using 
!  preliminary data on the CDIAC website for 2007-2008 and LeQuere et al. 
!  (2009) for 2009 
!                                                                             .
!      -- Ray Nassar 2010-03-10
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!  25 Feb 2011 - R. Nassar           - Now point to annual_v2010 and
!                                      monthly_v2010 directories, which
!                                      contain updated CO2 data from CDIAC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, YEAR, MONTH   
      TYPE (XPLEX)                 :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                 :: TAU, GLOB_SCL_FAC
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=4)       :: YEAR_STR
      CHARACTER(LEN=2)       :: MONTH_STR

      !=================================================================
      ! READ_FOSSILCO2 begins here!
      !=================================================================

      ! Time variables
      YEAR   = GET_YEAR()
      MONTH  = GET_MONTH()

      IF ( YEAR > 2006 ) THEN
         WRITE( 6, 90 ) 
 90      FORMAT('Fossil Fuel CO2 data for 2007-2009 are preliminary!')
      ENDIF
      IF ( YEAR > 2009 ) THEN
         YEAR = 2009
         WRITE( 6, 95 ) 
 95      FORMAT( 'YEAR > 2009; Using Fossil CO2 emissions for 2009!')
      ENDIF

      WRITE( YEAR_STR,  '(i4)'   ) YEAR 
      WRITE( MONTH_STR, '(i2.2)' ) MONTH
      
      IF ( LMONFF ) THEN
         LGENFF = .FALSE.
         LANNFF = .FALSE.
      ENDIF

      IF ( LANNFF ) THEN
         LGENFF = .FALSE.
      ENDIF

      !=================================================================
      ! Filename and tau for fossil fuel CO2 data 
      !=================================================================
      IF ( LGENFF ) THEN

         TAU      = GET_TAU0( 1, 1, 1985 )
         FILENAME = TRIM( DATA_DIR )        // 
     &              'CO2_200508/fossil95_CO2.' // GET_NAME_EXT_2D() //
     &              '.'                        // GET_RES_EXT()

         WRITE( 6, 100 ) TRIM( FILENAME )

      ELSE IF ( LANNFF ) THEN

         TAU      = GET_TAU0( 1, 1, YEAR )
         FILENAME = TRIM( DATA_DIR )                                // 
!------------------------------------------------------------------------------
! Prior to 2/25/11: 
! Now use updated CO2 annual emissions from CDIAC (cf Bob Andres)
! (rnassar, bmy, 2/25/11)
!     &              'CO2_201003/fossilfuel_andres/annual/ff.'       // 
!------------------------------------------------------------------------------
     &              'CO2_201003/fossilfuel_andres/annual_v2010/ff.' // 
     &              YEAR_STR          // '.'                        //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         WRITE( 6, 110 ) TRIM( FILENAME )

      ELSE IF ( LMONFF ) THEN

         TAU      = GET_TAU0( MONTH, 1, YEAR )
         FILENAME = TRIM( DATA_DIR )                                //
!------------------------------------------------------------------------------
! Prior to 2/25/11: 
! Now use updated CO2 annual emissions from CDIAC (cf Bob Andres)
! (rnassar, bmy, 2/25/11)                          // 
!     &              'CO2_201003/fossilfuel_andres/monthly/ff.'      //  
!------------------------------------------------------------------------------
     &              'CO2_201003/fossilfuel_andres/monthly_v2010/ff.'//
     &              YEAR_STR          // MONTH_STR // '.'           //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         WRITE( 6, 120 ) TRIM( FILENAME )
     
      ENDIF
     
      ! FORMATS
 100  FORMAT( '     - READ_GENERIC_FOSSCO2: Reading ', a ) 
 110  FORMAT( '     - READ_ANNUAL_FOSSCO2: Reading ',  a ) 
 120  FORMAT( '     - READ_MONTHLY_FOSSCO2: Reading ', a ) 

      !=================================================================
      ! Read fossil fuel CO2 [molec/cm2/s]
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      CALL TRANSFER_2D( ARRAY(:,:,1), EMFOSSCO2 )

      IF ( LPLANE ) THEN
         CALL AVIATION_DOM_CORR( EMFOSSCO2 )
      ENDIF

      END SUBROUTINE READ_FOSSILCO2
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_surf
!
! !DESCRIPTION: This subroutine reads the fossil fuel distribution from 
!  file to be used for part of the spatial distribution of the CO2 surface 
!  correction, based on a value of 4.89%, similar to that used in 
!  Suntharalingam et al. (2005).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_SURF
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_TAU0,   READ_BPCH2
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE LOGICAL_MOD,    ONLY : LGENFF, LANNFF, LMONFF
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      USE TIME_MOD,       ONLY : GET_YEAR,GET_MONTH
      USE GRID_MOD,       ONLY : GET_AREA_CM2

#     include "CMN_SIZE"       ! Size parameters
!
! !REMARKS:
!  Methane source distribution are read for the same purpose from 2004 data
!  provided by Kevin Wecht.
!                                                                             .
!  Monoterpenes and Isoprene are read and treated as representative NMVOCs.
!                                                                             .
!     -- Ray Nassar 2010-03-27
!
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version 
!  25 Feb 2011 - R. Nassar           - Now point to annual_v2010 and
!                                      monthly_v2010 directories, which
!                                      contain updated CO2 data from CDIAC
!EOP!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=4)        :: YEAR_STR
      CHARACTER(LEN=2)        :: MONTH_STR

      INTEGER                 :: I, J, YEAR, MONTH
      INTEGER, PARAMETER      :: NDAYS2004(12) = (/31,29,31,30,31,30,
     &                                             31,31,30,31,30,31/)

      TYPE (XPLEX)                  :: ARRAY(IIPAR,JJPAR,1)

      TYPE (XPLEX)                 :: TAU, MONFAC, NMHCFAC, A_CM2(JJPAR)
      
      TYPE (XPLEX)                  :: FOSS_MASS(IIPAR,JJPAR),  FOSS_SUM
      TYPE (XPLEX)                  :: LIVE_MASS(IIPAR,JJPAR),  LIVE_SUM
      TYPE (XPLEX)                 :: WASTE_MASS(IIPAR,JJPAR), WASTE_SUM
      TYPE (XPLEX)                  :: RICE_MASS(IIPAR,JJPAR),  RICE_SUM
      TYPE (XPLEX)                  :: WET_MASS(IIPAR,JJPAR),   WET_SUM
      TYPE (XPLEX)                 :: OTHER_MASS(IIPAR,JJPAR), OTHER_SUM
      TYPE (XPLEX)                  :: ISO_MASS(IIPAR,JJPAR),   ISO_SUM
      TYPE (XPLEX)                  :: MONO_MASS(IIPAR,JJPAR),  MONO_SUM
      TYPE (XPLEX)                  :: TOT_MASS(IIPAR,JJPAR),   TOT_SUM
      
      TYPE (XPLEX)                  :: FOSSIL_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: LIVE_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: WASTE_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: RICE_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: WET_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: OTHER_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: ISO_CORR(IIPAR,JJPAR)
      TYPE (XPLEX)                  :: MONO_CORR(IIPAR,JJPAR)

      TYPE (XPLEX), PARAMETER:: PERCENT_CORRECTION = xplex(4.89d0,0d0)

      !For # molecules <--> mass in kg
      TYPE (XPLEX), PARAMETER:: CH4FAC =xplex( 6.022d23/16d-3,0d0)
      TYPE (XPLEX), PARAMETER       :: CFAC= xplex(6.022d23/12d-3,0d0)

      !-----------------------------------------------------------------
      ! Get month and year
      !-----------------------------------------------------------------
      MONTH  = GET_MONTH()
      YEAR   = GET_YEAR()

      WRITE( YEAR_STR,  '(i4)'   ) YEAR 
      WRITE( MONTH_STR, '(i2.2)' ) MONTH 

      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2(J)
      ENDDO

      !-----------------------------------------------------------------
      ! Read Generic or annual or monthly fossil fuel emissions file
      !-----------------------------------------------------------------
      IF ( LMONFF ) THEN
         LGENFF = .FALSE.
         LANNFF = .FALSE.
      ENDIF

      IF ( LANNFF ) THEN
         LGENFF = .FALSE.
      ENDIF

      IF ( LGENFF ) THEN

         TAU      = GET_TAU0( 1, 1, 1985 )
         FILENAME = TRIM( DATA_DIR )                                 // 
     &              'CO2_200508/fossil95_CO2.'                       // 
     &               GET_NAME_EXT_2D() // '.'      // GET_RES_EXT()

      ELSE IF ( LANNFF ) THEN

         TAU      = GET_TAU0( 1, 1, YEAR )
         FILENAME = TRIM( DATA_DIR )                                 // 
!------------------------------------------------------------------------------
! Prior to 2/25/11: 
! Now use updated CO2 annual emissions from CDIAC (cf Bob Andres)
! (rnassar, bmy, 2/25/11)                   
!     &              'CO2_201003/fossilfuel_andres/annual/ff.'        // 
!------------------------------------------------------------------------------
     &              'CO2_201003/fossilfuel_andres/annual_v2010/ff.'  // 
     &              YEAR_STR          // '.'                         //
     &              GET_NAME_EXT_2D() // '.'       //  GET_RES_EXT()

      ELSE IF ( LMONFF ) THEN

         TAU      = GET_TAU0( MONTH, 1, YEAR )
         FILENAME = TRIM( DATA_DIR )                                 // 
!------------------------------------------------------------------------------
! Prior to 2/25/11: 
! Now use updated CO2 annual emissions from CDIAC (cf Bob Andres)
! (rnassar, bmy, 2/25/11)    
!     &              'CO2_201003/fossilfuel_andres/monthly/ff.'      //  
!------------------------------------------------------------------------------
     &              'CO2_201003/fossilfuel_andres/monthly_v2010/ff.' //  
     &              YEAR_STR          // MONTH_STR // '.'            //
     &              GET_NAME_EXT_2D() // '.'       // GET_RES_EXT()

      ENDIF
	
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         FOSSIL_corr(I,J) = (PERCENT_CORRECTION/100d0)*ARRAY(I,J,1)
      ENDDO
      ENDDO

      !-----------------------------------------------------------------
      ! Read Monthly CH4 emissions
      !-----------------------------------------------------------------
      TAU      = GET_TAU0( MONTH, 1, 2004 )

      FILENAME = TRIM( DATA_DIR )                 // 
     &           'CO2_201003/ChemSrc/CH4_source.' // 
     &           GET_NAME_EXT_2D() // '.'         //
     &           GET_RES_EXT()

      ! %%% Livestock %%%
      WRITE( 6, 40 ) TRIM( FILENAME )
 40   FORMAT( '     - READ_LIVESTOCK: Reading ', a ) 

      CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 4, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  LIVE_MASS)

      ! %%% Waste %%%
      WRITE( 6, 50 ) TRIM( FILENAME )
 50   FORMAT( '     - READ_WASTE: Reading ', a ) 

      CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 5, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  WASTE_MASS)

      ! %%% Rice %%%
      WRITE( 6, 60 ) TRIM( FILENAME )
 60   FORMAT( '     - READ_RICE: Reading ', a ) 

      CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 7, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  RICE_MASS)

      ! %%% Wetlands %%%
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_WETLANDS: Reading ', a ) 

      CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 10, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  WET_MASS)

      ! %%% Natural %%%
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - READ_OTHER_NATURAL: Reading ', a ) 

      CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 12, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  OTHER_MASS)

      !-----------------------------------------------------------------
      ! Print raw monthly totals in Tg CH4
      !-----------------------------------------------------------------
      LIVE_SUM  = sum(LIVE_MASS(:,:)*1d-9)
      WASTE_SUM = sum(WASTE_MASS(:,:)*1d-9)
      RICE_SUM  = sum(RICE_MASS(:,:)*1d-9)
      WET_SUM   = sum(WET_MASS(:,:)*1d-9)
      OTHER_SUM = sum(OTHER_MASS(:,:)*1d-9)

      write(6,200) " GLOBAL LIVESTOCK    ", LIVE_SUM, "  TgCH4/month"
      write(6,200) " GLOBAL WASTE        ", WASTE_SUM,"  TgCH4/month"
      write(6,200) " GLOBAL RICE         ", RICE_SUM, "  TgCH4/month"
      write(6,200) " GLOBAL WETLANDS     ", WET_SUM,  "  TgCH4/month"
      write(6,200) " GLOBAL OTHER NATURAL", OTHER_SUM,"  TgCH4/month"

      !-----------------------------------------------------------------
      ! Convert kg/gridbox/month to molecules/cm2/s
      !-----------------------------------------------------------------
      MONFAC = ndays2004(month)*86400d0
 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         LIVE_CORR(I,J)      = LIVE_MASS(I,J)*CH4FAC/MONFAC/A_CM2(J)
         WASTE_CORR(I,J)     = WASTE_MASS(I,J)*CH4FAC/MONFAC/A_CM2(J)
         RICE_CORR(I,J)      = RICE_MASS(I,J)*CH4FAC/MONFAC/A_CM2(J)
         WET_CORR(I,J)       = WET_MASS(I,J)*CH4FAC/MONFAC/A_CM2(J)
         OTHER_CORR(I,J)     = OTHER_MASS(I,J)*CH4FAC/MONFAC/A_CM2(J)
      ENDDO
      ENDDO

      !-----------------------------------------------------------------
      ! Read Monthly Isoprene emissions
      !-----------------------------------------------------------------
      TAU = GET_TAU0( MONTH, 1, 2004 )

      FILENAME = TRIM( DATA_DIR )                            //
     &           'CO2_201003/ChemSrc/Isoprene-2004.'         //
     &           GET_NAME_EXT_2D()  // '.' // GET_RES_EXT()

      WRITE( 6, 150 ) TRIM( FILENAME )
 150  FORMAT( '     - READ_ISOPRENE: Reading ', a ) 

      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  ISO_corr)

      !-----------------------------------------------------------------
      ! Read Monthly Monoterpene emissions
      !-----------------------------------------------------------------

      FILENAME = TRIM( DATA_DIR )                            // 
     &           'CO2_201003/ChemSrc/Monoterpene-2004.'      //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      WRITE( 6, 160 ) TRIM( FILENAME )
 160  FORMAT( '     - READ_MONOTERPENE: Reading ', a ) 

      ! NOTE: use same TAU0 as for isoprene
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 4, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
      CALL TRANSFER_2D( ARRAY(:,:,1),  MONO_corr)

      !-----------------------------------------------------------------
      ! Take the sum of all surface corrections. 
      !-----------------------------------------------------------------
      ! NMHCFAC is a scale factor which combines the CO yield from 
      ! monoterpenes and isoprenes (~0.2) but increase it to use their
      ! spatial distribution as a proxy for other NMHCs.
      !-----------------------------------------------------------------

      NMHCFAC = 0.333d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         EMIS_SUB(I,J) = FOSSIL_corr(I,J)     
     &                 + LIVE_corr(I,J) 
     &                 + WASTE_corr(I,J)     
     &                 + RICE_corr(I,J)      
     &                 + WET_corr(I,J)      
     &                 + OTHER_corr(I,J)  
     &                 + ISO_corr(I,J)*NMHCFAC
     &                 + MONO_corr(I,J)*NMHCFAC
      ENDDO
      ENDDO

      !-----------------------------------------------------------------

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         FOSS_MASS(I,J) = FOSSIL_CORR(I,J)*A_CM2(J)*MONFAC/CFAC
         ISO_MASS(I,J)  = ISO_CORR(I,J)*A_CM2(J)*NMHCFAC*MONFAC/CFAC
         MONO_MASS(I,J) = MONO_CORR(I,J)*A_CM2(J)*NMHCFAC*MONFAC/CFAC
         TOT_MASS(I,J)  = EMIS_SUB(I,J)*A_CM2(J)*MONFAC/CFAC
      ENDDO
      ENDDO

      FOSS_SUM = sum(FOSS_MASS(:,:)*1d-9)
      ISO_SUM  = sum(ISO_MASS(:,:)*1d-9)
      MONO_SUM = sum(MONO_MASS(:,:)*1d-9)
      TOT_SUM  = sum(TOT_MASS(:,:)*1d-9)
	 
	write(6,200) " GLOBAL FOSS CORR       ", FOSS_SUM, " TgC/month"
	write(6,200) " GLOBAL ISOPRENE        ", ISO_SUM,  " TgC/month"
	write(6,200) " GLOBAL MONTERPENE      ", MONO_SUM, " TgC/month"
	write(6,200) " GLOBAL TOTAL SURF CORR ", TOT_SUM,  " TgC/month"

 200  FORMAT( A, F9.5, A ) 

      END SUBROUTINE CHEM_SURF  
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aviation_dom_corr
!
! !DESCRIPTION: This subroutine downscales national fossil fuels emissions 
!  for the CO2 which is atttibuted to domestic aviation based on Kim et al. 
!  (2005,2007).  It should only be used when the aviation emissions are 
!  turned on since these emissions will instead be emitted throughout the 
!  troposphere.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AVIATION_DOM_CORR( EMFOSS )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_TAU0,   READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LGENFF
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,      ONLY : GET_YEAR, ITS_A_LEAPYEAR
      USE GRID_MOD,      ONLY : GET_AREA_CM2

#     include "CMN_SIZE"      ! Size parameters
!
! !INPUT PARAMETERS:
!
      TYPE (XPLEX), INTENT(INOUT)  :: EMFOSS(IIPAR,JJPAR)  ! Fuel to be scaled

! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version 
!  25 Feb 2011 - R. Nassar           - Now point to annual_v2010 and
!                                      monthly_v2010 directories, which
!                                      contain updated CO2 data from CDIAC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=4)   :: YEAR_STR
      INTEGER            :: I, J, N, YEAR, NDAYS
      TYPE (XPLEX)             :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)             :: REGION(IIPAR,JJPAR)
      TYPE (XPLEX)             :: TAU, U_CONV
      TYPE (XPLEX)             :: DOM_COR(8), REG_SUM(8), AV_SCLFAC(8)
      TYPE (XPLEX)             :: ANN_FOSS(IIPAR,JJPAR), A_CM2(JJPAR)
      TYPE (XPLEX)             :: ANN_FOSS_NOAREA(IIPAR,JJPAR)

      !-----------------------------------------------------------------
      ! Read Generic or annual fossil fuel emissions file
      !-----------------------------------------------------------------
      YEAR   = GET_YEAR()
      WRITE( YEAR_STR, '(i4)' ) YEAR 

      print*, " Staring Domestic Aviation Correction"

      IF ( LGENFF ) THEN
         TAU      = GET_TAU0( 1, 1, 1985 )
         FILENAME = TRIM( DATA_DIR )                               // 
     &              'CO2_200508/fossil95_CO2.'                     // 
     &              GET_NAME_EXT_2D() // '.'   // GET_RES_EXT()
      ELSE
         TAU      = GET_TAU0( 1, 1, YEAR )
         FILENAME = TRIM( DATA_DIR )                               // 
!------------------------------------------------------------------------------
! Prior to 2/25/11: 
! Now use updated CO2 annual emissions from CDIAC (cf Bob Andres)
! (rnassar, bmy, 2/25/11)        
!     &              'CO2_201003/fossilfuel_andres/annual/ff.'      //  
!------------------------------------------------------------------------------
     &              'CO2_201003/fossilfuel_andres/annual_v2010/ff.'//  
     &              YEAR_STR          // '.'                       //
     &              GET_NAME_EXT_2D() // '.'   //  GET_RES_EXT()
      ENDIF
	
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16
      CALL TRANSFER_2D( ARRAY(:,:,1), ANN_FOSS )

      !-----------------------------------------------------------------
      ! Calculate box areas, covert to area-independent emission mass
      ! i.e. molecules/cm2/s --> Tg/gridbox/year
      !-----------------------------------------------------------------
      NDAYS = 365
      IF ( ITS_A_LEAPYEAR() ) NDAYS = 366

      U_CONV = NDAYS*86400*(12./6.022D23)*1D-12

      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2(J)
         ANN_FOSS_NOAREA(:,J) =  ANN_FOSS(:,J)*A_CM2(J)*U_CONV
      ENDDO

      !-----------------------------------------------------------------
      ! Read region/continent file 
      !-----------------------------------------------------------------
      TAU      = GET_TAU0( 1, 1, 2004 )
      FILENAME = TRIM( DATA_DIR )                          //  
     &           'CO2_201003/Aviation_Regions.'            // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         REGION,    QUIET=.TRUE. )

      !-----------------------------------------------------------------
      ! Sum fossil fuel CO2 (as C) in each region
      !-----------------------------------------------------------------
      REG_SUM(:) = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO N = 1, 8
         IF ( NINT(REGION(I,J) ) == N ) THEN
            REG_SUM(N) =  REG_SUM(N) + ANN_FOSS_NOAREA(I,J)
         ENDIF	
      ENDDO
      ENDDO
      ENDDO

      !-----------------------------------------------------------------
      ! Mean domestic CO2 from aviation for 2000-2004 in TgC/yr 
      ! from Kim et al. (2005). See report or bpch for exact regions.
      !-----------------------------------------------------------------
      DOM_COR(1) =  49.6   ! North America                             
      DOM_COR(2) =   2.8        ! South America                             
      DOM_COR(3) =   1.8   ! Eastern Europe                          
      DOM_COR(4) =  12.3   ! Western Europe     
      DOM_COR(5) =  16.1  ! Asia     
      DOM_COR(6) =   1.0  ! Africa                   
      DOM_COR(7) =   2.2  ! MiddleEast 
      DOM_COR(8) =   2.0  ! Oceania                     

      !-----------------------------------------------------------------
      ! Calculate aviation scale factors then apply them
      !-----------------------------------------------------------------
      DO N = 1,8
         AV_SCLFAC(N) = (REG_SUM(N) - DOM_COR(N)) / REG_SUM(N)
      ENDDO

      print*, " Scaling down EMFOSS "

      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO N = 1, 8
	 IF (NINT(REGION(I,J))==N) EMFOSS(I,J) = EMFOSS(I,J)*AV_SCLFAC(N)
      ENDDO
      ENDDO
      ENDDO
	
      print*, 
     & " Domestic Aviation CO2 subtracted from land fossil fuel CO2"

      END SUBROUTINE AVIATION_DOM_CORR  
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_oceanco2
!
! !DESCRIPTION: Subroutine READ\_OCEANCO2 reads in either
!
! \begin{itemize}
! \item Annual mean oceanic CO2 exchange from Takahashi 1997
! \item Annual mean oceanic CO2 exchange from Takahashi 2009
! \item Aonthly mean oceanic CO2 exchange from Takahashi 2009
! \end{itemize}
!
!  from a binary punch file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_OCEANCO2
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_MONTH
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE LOGICAL_MOD,   ONLY : LOCN1997, LOCN2009ANN, LOCN2009MON

#     include "CMN_SIZE"      ! Size parameters
!
! !REMARKS:
!  See References Above
!
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!  25 Feb 2011 - R. Nassar           - Now point to annual_v2010 and
!                                      monthly_v2010 directories, which
!                                      contain updated CO2 data from CDIAC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, MONTH
      TYPE (XPLEX)                 :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=2)       :: MONTH_STR

      !=================================================================
      ! READ_OCEANCO2 begins here!
      !=================================================================

      MONTH = GET_MONTH()
      WRITE( MONTH_STR, '(i2.2)' ) MONTH
      
      IF ( LOCN1997 ) THEN

         ! Start of "generic" year 1985 
         TAU = GET_TAU0( 1, 1, 1985 ) 

         FILENAME = TRIM( DATA_DIR )                             // 
     &              'CO2_200508/ocean_CO2.'                      //
     &               GET_NAME_EXT_2D() //  '.' // GET_RES_EXT()
     
      ELSE IF ( LOCN2009ANN ) THEN
	
         ! Start of "generic" year 2000 
         TAU      = GET_TAU0( 1, 1, 2000 ) 

         FILENAME = TRIM( DATA_DIR )                             // 
     &              'CO2_201003/Ocean/Taka2009_OceanCO2_annual.' // 
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

	
      ELSE IF ( LOCN2009MON ) THEN

         ! Start of "generic" month in 2000 
         TAU      = GET_TAU0( MONTH, 1, 2000 ) 

         FILENAME = TRIM( DATA_DIR  )                            // 
     &              'CO2_201003/Ocean/Taka2009_OceanCO2_'        // 
     &              MONTH_STR         // '.'                     // 
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ENDIF

      WRITE( 6, 100 ) TRIM( FILENAME )
	
      ! Read ocean CO2 data [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 2, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX) and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMOCCO2 )

      ! Return to calling program
 100  FORMAT( '     - READ_OCEANCO2: Reading ', a ) 

      END SUBROUTINE READ_OCEANCO2
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_annual_biofuelco2
!
! !DESCRIPTION: Subroutine READ\_ANNUAL\_BIOFUELCO2 reads in annual mean 
!  biofuel CO2 emissions from a binary punch file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ANNUAL_BIOFUELCO2
!
! !USES:
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
!
! !REMARKS:
!  References:
!  (1 ) Yevich and Logan 2001 gridded (1x1) dataset in combination with 
!        emission factors for CO2 per kg drymatter burned
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J
      ! adj_group: (zhej, dkh, 01/17/12, adj32_015) 
      !TYPE (XPLEX)                 :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                 :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_BIOFUELCO2 begins here!
      ! Use 1985 emissions or 1995 scaled values
      ! "Burn in fields" not included (this is already in GFED)
      !=================================================================

      FILENAME = TRIM( DATA_DIR )                 // 
     &          'CO2_201003/biofuel/biofuel_CO2.' // 
     &           GET_NAME_EXT_2D() // '.'         // 
     &           GET_RES_EXT()     // '-1995'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_BIOFUELCO2: Reading ', a ) 

!     TAU = GET_TAU0( 1, 1, 1985 )
      TAU = GET_TAU0( 1, 1, 1995 )

      ! Read biofuel CO2 emissions [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 5, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOFUELCO2 )

      END SUBROUTINE READ_ANNUAL_BIOFUELCO2
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_shipco2_edgar
!
! !DESCRIPTION: Subroutine READ\_SHIPCO2\_EDGAR reads in annual mean ship CO2 
!  emissions from a binary punch file.  Scaling is based on Endresen et al. 
!  (2007).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SHIPCO2_EDGAR
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR, DATA_DIR_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1 
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE TIME_MOD,       ONLY : GET_YEAR
	
#     include "CMN_SIZE"       ! Size parameters
!
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, N, YEAR
      TYPE (XPLEX)		 :: A_CM2(JJPAR)
      TYPE (XPLEX)             :: ARRAY(360,180,1)
      TYPE (XPLEX)             :: GEN_1x1(360,180,1)
      TYPE (XPLEX)             :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)             :: GEOS_GRID(IIPAR,JJPAR,1)
      TYPE (XPLEX)             :: TAU
      TYPE (XPLEX),PARAMETER::SEC_IN_YEAR=xplex(86400d0 * 365.25d0,0d0)
      TYPE (XPLEX),  PARAMETER :: GlobSTot = xplex(0.117 ,0d0)
                             ! EDGAR global total value in PgC/yr
      TYPE (XPLEX)             :: GlobSTotNew(25) 
                             ! For scaling to a specified global total PgC/yr

      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_SHIPCO2_EDGAR begins here!
      !=================================================================

      YEAR = GET_YEAR()

      ! Fill array
      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2( J )
      ENDDO

      ! File contaning EDGAR ship CO2 data
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &     'ARCTAS_SHIP_2008/Arctas_CO2_ship_2008.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_SHIP_CO2_EDGAR: Reading ', a ) 

      ! TAU value for start of 2008
      TAU = GET_TAU0( 1, 1, 2008 )

      ! Read CO2 shipping emissions [kg/yr]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       360,     180,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX) before regridding
      GEN_1x1(:,:,1) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 --> GEOS 1x1 
      CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1 )

      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 1, 'kg/yr', GEOS_1x1, GEOS_GRID)

      ! Convert units kg CO2 / yr --> molecules/cm2/s 
      DO J = 1, JJPAR
      	EMSHIPCO2(:,J) = GEOS_GRID(:,J,1) * 6.022D23 
     &            / ( 44.0d-3 * SEC_IN_YEAR * A_CM2(J))
      ENDDO

      !-----------------------------------------------------------------
      ! Global Ship Totals for the years 1985 to 2009
      !-----------------------------------------------------------------
      ! These were based on a linear fit to 1985 to 2002 
      ! values from Endresen et al. (2007)
      ! The 2007 value was used for 2009 similar to the case for 
      ! overall fossil fuel use related to the recession
      ! We have not attempted to extrapolate beyond 2009
      !-----------------------------------------------------------------
      ! Values are essentially all international (domestic negligible)
      !-----------------------------------------------------------------

      GlobSTotNew(1:25)%r = (/ 0.122, 0.125, 0.128, 0.132, 0.135,
     &                       0.138, 0.141, 0.144, 0.148, 0.151,
     &                       0.154, 0.157, 0.160, 0.163, 0.167,
     &                       0.170, 0.173, 0.176, 0.179, 0.182,
     &                       0.186, 0.189, 0.192, 0.195, 0.192 /)
      GlobSTotNEW(1:25)%i = 0d0  
      IF (YEAR <= 1985) THEN
         n = 1
      ELSEIF (YEAR >= 2009) THEN
         n = 25
      ELSE
         n = YEAR - 1985 + 1
      ENDIF

      ! Apply scaling to obtain a specified yearly global total
      DO J = 1, JJPAR
         EMSHIPCO2(:,J) = EMSHIPCO2(:,J)*(GlobSTotNew(n)/GlobSTot)
      ENDDO

      END SUBROUTINE READ_SHIPCO2_EDGAR
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_shipco2_icoads
!
! !DESCRIPTION: Subroutine READ\_SHIPCO2\_ICOADS reads in ICOADS monthly 
!  ship CO2 emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SHIPCO2_ICOADS
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR, DATA_DIR_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      USE FILE_MOD,       ONLY : IU_FILE, IOERROR
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1 
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE TIME_MOD,       ONLY : GET_YEAR, GET_MONTH
	
#     include "CMN_SIZE"       ! Size parameters
!
! !REMARKS:
!  This subroutine reads from bpch files at GEOS 1x1 (half-polar) resolution
!  although the original data are provided as 0.1 deg x 0.1 deg. Regridding to 
!  the current resolution is carried out in the code.
!                                                                             .
!  References:
!  (1) Corbett and Koehler (2003) "Updated emissions from ocean shipping", 
!       JGR 108, D20, 4650.
!  (2) Corbett and Koehler (2004) "Considering alternative input parameters in 
!       an activity-based ship fuel consumption and emissions model: Reply ..."
!       JGR, 109, D23303.
!  (3) Endresen et al. (2007) "A historical reconstruction of ships fuel 
!       consumption and emissions", JGR, 112, D12301.
!                                                                             .
!  NOTE: The Corbett website values do not sum to the values in any Corbett
!  et al. or Wang (2008) papers. It is not clear if this relates to the 
!  ongoing dispute between Corbett et al.(2003,2004) and Endresen et al.
!  (2003,2004,2007)
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, MONTH, YEAR, IOS, N
      TYPE (XPLEX)		     :: A_CM2(JJPAR)
      ! adj_group: (zhej, dkh, 01/17/12, adj32_015) 
      !TYPE (XPLEX)                 :: ARRAY(360,181,1)
      TYPE (XPLEX)                 :: ARRAY(I1x1,J1x1,1)
      TYPE (XPLEX)                 :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                 :: GEOS_GRID(IIPAR,JJPAR,1)
      TYPE (XPLEX)                 :: TAU
      TYPE (XPLEX),PARAMETER::SEC_IN_YEAR=xplex(86400d0*365.25d0,0d0)
      TYPE (XPLEX),  PARAMETER     :: GlobSTot = xplex(0.1760,0d0) 
                                !ICOADS global total value in PgC/yr
      TYPE (XPLEX)                 :: GlobSTotNew(25) 
                                !For scaling to a specified global total PgC/yr

      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=2)       :: MONTH_STR
	
      !=================================================================
      ! READ_SHIPCO2_ICOADS begins here!
      !=================================================================

      ! Fill array
      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2( J )
      ENDDO

      YEAR = GET_YEAR()
      MONTH = GET_MONTH()
 
      WRITE( MONTH_STR, '(i2.2)' ) MONTH
      
      FILENAME = TRIM( DATA_DIR_1x1 )           // 
     &           'CO2_201003/ship_ICOADS/co2_'  // 
     &           MONTH_STR // '.geos.1x1'
	
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_MONTHLY_SHIP_CO2_ICOADS: Reading ', a ) 

      ! TAU value for start of month in 2004
      TAU = GET_TAU0( MONTH, 1, 2004 )

      ! Read CO2 shipping emissions
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 7, 
     &                 TAU,       360,       181,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX) before regridding
      GEOS_1x1(:,:,1) = ARRAY(:,:,1)

      ! Regrid to current model resolution
      CALL DO_REGRID_1x1( 1, 'molec/cm2/s', GEOS_1x1, GEOS_GRID )

      DO J = 1, JJPAR
         EMSHIPCO2(:,J) = GEOS_GRID(:,J,1) 
      ENDDO

      !-----------------------------------------------------------------
      ! Global Ship Totals for the years 1985 to 2009
      !-----------------------------------------------------------------
      ! These were based on a linear fit to 1985 to 2002 
      ! values from Endresen et al. (2007)
      ! The 2007 value was used for 2009 similar to the case for 
      ! overall fossil fuel use related to the recession
      ! We have not attempted to extrapolate beyond 2009
      !-----------------------------------------------------------------
      ! Values are essentially all international (domestic negligible)
      !-----------------------------------------------------------------

      GlobSTotNew(1:25)%r = (/ 0.122, 0.125,0.128, 0.132, 0.135,
     &                       0.138, 0.141,0.144, 0.148, 0.151,
     &                       0.154, 0.157,0.160, 0.163, 0.167,
     &                       0.170, 0.173,0.176, 0.179, 0.182,
     &                       0.186, 0.189,0.192, 0.195, 0.192 /)
      GlobSTotNew(1:25)%i = 0d0
      IF ( YEAR <= 1985 ) THEN
         n = 1
      ELSEIF ( YEAR >= 2009 ) THEN
         n = 25
      ELSE
         n = YEAR - 1985 + 1
      ENDIF
      
      ! Apply scaling to obtain a specified yearly global total
      DO J = 1, JJPAR
         EMSHIPCO2(:,J) = EMSHIPCO2(:,J)*(GlobSTotNew(n)/GlobSTot)
      ENDDO

      END SUBROUTINE READ_SHIPCO2_ICOADS
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aviation_co2
!
! !DESCRIPTION: Subroutine READ\_AVIATION\_CO2 reads monthly mean aircraft 
!  fuel emissions and converts them to CO2 emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_AVIATION_CO2
!
! !USES:
!
      ! Reference to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT, GET_TAU0, READ_BPCH2
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE TIME_MOD,      ONLY : GET_MONTH, GET_YEAR

#     include "CMN_SIZE"      ! Size parameters
!
! !REMARKS:
!  This is a modified version of READ_AIRCRAFT_SO2 from:
!                                            rjp, bdf, bmy, 9/18/02, 10/3/05
!                                                                             .
!  The sulfate data are based on an inventory by the Atmospheric Effects of 
!  Aviation Project (AEAP) for the year 1992.  
!                                                                             .
!  CO2 emission factor of 3155 g/kg fuel was taken from
!  (1) Kim et al. (2005) System for assessing Aviation's Global Emissions 
!       (SAGE) Federal Aviation Administration Office of Environment and 
!       Energy Version 1.5 (FAA-EE-2005-02), Global Aviation Emissions 
!       Inventories for 2000 through 2004.
!  (2) Kim et al. (2007) System for assessing Aviation's Global Emissions
!       (SAGE) Part 1: Model description and inventory results 
! 
! !REVISION HISTORY: 
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Reading of GlobPTot values from input.geos has not yet been implemented
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, IOS, J, K, L, THISMONTH, N, YEAR
      TYPE (XPLEX)              :: ACCO2(IGLOB,JGLOB,20)
      TYPE (XPLEX)              :: FAC, FUEL, DZ(LLPAR), ZH(0:LLPAR)
      TYPE (XPLEX), PARAMETER   :: TINY = xplex(1d-20,0d0)
      TYPE (XPLEX), EXTERNAL    :: BOXVL

      ! For scaling to a specified global total PgC/yr
      TYPE (XPLEX), PARAMETER   :: GlobPTot    = xplex(0.12,0d0) !1992 AEAP estimate
      TYPE (XPLEX)              :: GlobPTotNew(22)
	
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=3)    :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr', 
     &                                      'may', 'jun', 'jul', 'aug',
     &                                      'sep', 'oct', 'nov', 'dec'/)

      !=================================================================
      ! READ_AVIATION_CO2 begins here!
      !=================================================================
      
      ! Zero arrays
      EMPLANECO2 = 0d0
      ACCO2   = 0d0
      
	THISMONTH = GET_MONTH()

      ! File name
      FILENAME = TRIM( DATA_DIR )               // 
     &           'sulfate_sim_200508/aircraft.' // GET_RES_EXT() //
     &           '.1992.'                       // CMONTH(THISMONTH)

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_AVIATION_CO2: Reading ', a )     

      !=================================================================
      ! Read aircraft emissions.  These are fuel burned in [kg/box/day],
      ! from AEAP for 1992.  CO2 emission is calculated by assuming    
      ! an emission index EI of 3155, i.e., 3155 g of CO2 emitted per kg    
      ! of fuel burned.  It is also assumed that there is no diurnal   
      ! variation of emission rate. Convert to [kg CO2/box/s]. 
      !=================================================================

      ! Open file 
      OPEN( IU_FILE, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_aviation_co2:1' )

      ! Read header line
      READ( IU_FILE, '(/)', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_aviation_co2:2' )
      
      ! Read data values until an EOF is found
      DO 
         READ( IU_FILE, '(3i4,e11.3)', IOSTAT=IOS ) I, J, L, FUEL

         ! EOF encountered
         IF ( IOS < 0 ) EXIT

         ! I/O error condition
         IF ( IOS > 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'read_aviation_co2:3' )
         ENDIF

         ! Unit conversion: [kg Fuel/box/day] -> [kg CO2/box/s]
         ! Assuming an emission index of 1.0, 
         ! 3155 g CO2 / kg fuel burned [Kim et al., 2005, 2007]
         ACCO2(I,J,L+1) = 3.155 * FUEL / ( 24.d0 * 3600d0 )

      ENDDO

      ! Close file
      CLOSE( IU_FILE )
	
      !=================================================================
      ! Interpolate from the 1-km grid to the given GEOS-CHEM grid
      ! NOTE: we need to account for window grids (bmy, 9/20/02)
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! ACCO2 is the aircraft CO2 on the 1-km vertical grid
         FUEL = SUM( ACCO2(I,J,:) )
         IF ( FUEL < TINY ) CYCLE

         ! There are 20 1-km levels
         DO K = 1, 20

            ! Initialize
            ZH(0) = 0.d0

            ! Loop over levels
            DO L = 1, LLPAR

               ! Altitude of top edge of level L, from ground [km]
               ZH(L) = ZH(L-1) + ( BXHEIGHT(I,J,L) * 1d-3 )
               
               IF ( ZH(L-1) > XPLX(K)   ) EXIT
               IF ( ZH(L  ) < XPLX(K-1) ) CYCLE
               
               IF ( ZH(L) < XPLX(K) ) THEN
                  FAC            = ZH(L) - MAX( ZH(L-1), XPLX(K-1) )
                  EMPLANECO2(I,J,L) = EMPLANECO2(I,J,L)+ACCO2(I,J,K)*FAC
               ELSE
                  FAC            = XPLX(K)-MAX( ZH(L-1), XPLX(K-1) )
                  EMPLANECO2(I,J,L) = EMPLANECO2(I,J,L)+ACCO2(I,J,K)*FAC
                  EXIT
               ENDIF		     
            ENDDO
         ENDDO     

      ENDDO
      ENDDO

      YEAR = GET_YEAR()

      !-----------------------------------------------------------------
      ! Global Aviation Totals for 1985-2006
      !-----------------------------------------------------------------
      ! 1985-1995 values come directly from Sausen & Schumann (2000)
      ! 1996-1999 value of 0.155 assumed for bridging
      ! 2000-2005 values come directly from Kim et al. (2007)
      ! 2006 value comes from Wilkerson et al. (2010)
      ! We have not attempted to extrapolate outside of this range
      !-----------------------------------------------------------------
      ! A correction for domestic fuel use is carried out in the 
      ! fossil fuel subroutine.
      !-----------------------------------------------------------------

	GlobPTotNew(1:22)%r = (/
     &       0.1234, 0.1299, 0.1356, 0.1414, 0.1465,
     &       0.1469, 0.1434, 0.1420, 0.1441, 0.1500, 0.1543,  ! 1985-1995
     &       0.1550, 0.1550, 0.1550, 0.1550,                  ! 1996-1999
     &       0.1560, 0.1462, 0.1470, 0.1519, 0.1620, 0.1748,  ! 2000-2005
     &       0.16225 /)                                        ! 2006
        GlobPTotNew(1:22)%i = 0d0
	IF ( YEAR <= 1985 ) THEN
           n = 1
        ELSEIF ( YEAR >= 2006 ) THEN
           n = 22
	ELSE
           n = YEAR - 2000 + 1
	ENDIF

      !=================================================================
      ! Convert units from kg/box/s to molecules/cm3/s
      !
      ! Notes:
      ! 1) box volume from BOXVL is in cm3
      ! 2) mass is for CO2 (not C) with emission factor of 3155 g/kg  
      ! 3) optional global scaling to specified value
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, LLPAR
         EMPLANECO2(I,J,L) = EMPLANECO2(I,J,L) / BOXVL(I,J,L)
     &                     * 6.022d23 / 44.01d-3

	 ! Apply scaling to obtain a specified global Total
         EMPLANECO2(I,J,L) = EMPLANECO2(I,J,L)*(GlobPTotNew(n)/GlobPTot)
      ENDDO     
      ENDDO
      ENDDO

      END SUBROUTINE READ_AVIATION_CO2
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_annual_bionet_co2
!
! !DESCRIPTION: Subroutine READ\_ANNUAL\_BIONET\_CO2 reads in annual mean 
!  values of for Net Terrestrial exchange from a binary punch file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ANNUAL_BIONET_CO2
!
! !USES:
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR, DATA_DIR_1x1
      USE FILE_MOD,       ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      USE LOGICAL_MOD,    ONLY : LBIONETORIG,   LBIONETCLIM
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1 

#     include "CMN_SIZE"       ! Size parameters
!
! !REMARKS:
!  The two choices are:
!  (1 ) Old Net Terrestrial Exchange for Year 2000 from David Baker 
!        (pers. comm.) from undocumented Transcom 3 inversion results
!  (2 ) New Baker et al [2006] Transcom 3 climatology 1991-2000 minus
!       GFEDv2 climatology 1997-2007.
!                                                                             .
!  References:
!  (1 ) Baker et al. (2006), Transcom3 inversion intercomparison: Impact of 
!        Transport model errors on the interannual vaiability of regional CO2
!        fluxes, 1988-2003, Glob. Biogeochem. Cycles, 20, GB1002.
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, IOS
      TYPE (XPLEX)                 :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)                 :: ARRAY_1x1(360,180,1)
      TYPE (XPLEX)                 :: GEN_1x1(360,180,1)
      TYPE (XPLEX)                 :: GEOS_1x1(I1x1,J1x1,1)
      TYPE (XPLEX)                 :: GEOS_GRID(IIPAR,JJPAR,1)
      TYPE (XPLEX)		     :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_BIONET_CO2 begins here!
      !=================================================================
      
      ! Initialize ARRAY
      ARRAY = 0e0

      !------------------------------------
      ! Read original Bionet data
      !------------------------------------

      ! Filename
      IF ( LBIONETORIG ) THEN

         FILENAME = TRIM( DATA_DIR )                //
     &              'CO2_200508/net_terr_exch_CO2.' // 
     &               GET_NAME_EXT_2D() // '.'       // 
     &               GET_RES_EXT()     // '.txt'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Open file
         OPEN( IU_FILE,          FILE=TRIM( FILENAME ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
         IF ( IOS > 0 ) CALL IOERROR(IOS,IU_FILE, 'read_ann_bionet:1')      

         ! Read data
         READ( IU_FILE, '(7e13.6)', IOSTAT=IOS )
     &        ( ( ARRAY(I,J,1), I=1,IIPAR), J=1,JJPAR )
         IF ( IOS > 0 ) CALL IOERROR(IOS,IU_FILE, 'read_ann_bionet:2')

         ! Close file
         CLOSE( IU_FILE )

         ! Cast to TYPE (XPLEX) and resize if necessary
         CALL TRANSFER_2D( ARRAY(:,:,1), EMBIONETCO2 )

      ENDIF

      !------------------------------------
      ! Read climatological Bionet data
      !------------------------------------       
      IF ( LBIONETCLIM ) THEN

         ! TAU value for start of "generic" year 2000
         TAU = GET_TAU0( 1, 1, 2000 )

         ! Filename
         FILENAME = TRIM( DATA_DIR )                           //
     &              'CO2_201003/Net_terrestrial_exch_5.29Pg.'  //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT() 
	
         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read Net Terrestrial CO2 Exchange [molec/cm2/s]
         CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 6, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

         ! Cast to TYPE (XPLEX) and resize if necessary
         CALL TRANSFER_2D( ARRAY(:,:,1), EMBIONETCO2 )

!-----------------------------------------------------------------------
! Commented block of code below is nearly equivalent to the block of code 
! above but regrids during the model run from generic 1x1. This approach  
! is better for working with resolutions finer than 2x2.5, but at present,
! requires a 2-step regridding which is not exaclty equivalent.
!-----------------------------------------------------------------------
!
!         FILENAME = TRIM( DATA_DIR_1x1 )             //
!     &      'CO2_201003/Net_terrestrial_exch_1x1_4.47Pg.bpch'
!	
!         ! Echo info
!         WRITE( 6, 100 ) TRIM( FILENAME )
!
!         ! Read Net Terrestrial CO2 Exchange [molec/cm2/s]
!         CALL READ_BPCH2( FILENAME, 'CO2-SRCE',  6, 
!     &                 TAU,       360,           180,     
!     &                 1,         ARRAY_1x1,     QUIET=.TRUE. )
!
!         ! Cast to TYPE (XPLEX) before regridding
!         GEN_1x1(:,:,1) = ARRAY_1x1(:,:,1)
!
!         ! Regrid from GENERIC 1x1 --> GEOS 1x1 
!         CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1 )
!
!         ! Regrid from GEOS 1x1 --> current model resolution
!         CALL DO_REGRID_1x1( 1, 'molec/cm2/s', GEOS_1x1, GEOS_GRID)
!
!	   EMBIONETCO2(:,:) = GEOS_GRID(:,:,1)
!
!-----------------------------------------------------------------------

      ENDIF

 100  FORMAT( '     - READ_ANNUAL_BIONETCO2: Read ', a ) 

      END SUBROUTINE READ_ANNUAL_BIONET_CO2 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_bbio_dailyaverage
!
! !DESCRIPTION: Subroutine READ\_DAILY\_BBIO\_CO2 reads in daily values for 
!  balanced biospheric exchange from a binary punch file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BBIO_DAILYAVERAGE( MONTH, DAY, DOY ) 
!
! !USES:
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,      ONLY : GET_YEAR,        ITS_A_LEAPYEAR

#     include "CMN_SIZE"      ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: MONTH   ! Current month (1-12)
      INTEGER, INTENT(IN) :: DAY     ! Current day (1-31)
      INTEGER, INTENT(IN) :: DOY     ! Current day of year (0-366)
!
! !REMARKS:
!  Data Source: CASA gridded (1x1) dataset for from M. Thompson
!  Monthly values interpolated to daily values : 365 daily files 
!  NB : These files DO NOT have the diurnal cycle in daily emissions
!  See routine ' ' to read in files with diurnal cycle imposed
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Added fixes for leapyears
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, YEAR      ! dbj
      TYPE (XPLEX)                 :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                 :: TAU
      CHARACTER(LEN=3  )     :: SDOY
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_BBIO_DAILYAVERAGE begins here!
      !=================================================================

      YEAR   = GET_YEAR()    ! dbj

      ! Create string for day of year
      ! ----------------------------------------------------------------
      ! Since 2000 is a leap-year, we need to add one to get the correct 
      ! day-of-year for non-leap years, if the date is after Feb 28.
      !                                                     ! Ray Nassar
      ! ----------------------------------------------------------------
      
      IF( ITS_A_LEAPYEAR() .OR. DOY <= 59 ) THEN
         WRITE( SDOY, '(i3.3)' ) DOY 
      ELSE
         WRITE( SDOY, '(i3.3)' ) DOY + 1 
      ENDIF

      ! Get TAU value corresponding to DOY in year 2000
      TAU = GET_TAU0( MONTH, DAY, 2000 )

!      write(*,*) 'BB day ave DOY, SDOY, Tau = ',DOY, SDOY, Tau
      !-----------------------------------------------------------------

!      ! Make a string from DOY
!      WRITE( SDOY, '(i3.3)' ) DOY 

      ! Name of file with Balanced Bio CO2 data
      FILENAME = TRIM( DATA_DIR )                      //
     &           'CO2_200508/BBIO_DAILYAVG/CO2.daily.' //
     &           GET_NAME_EXT_2D() // '.'              //
     &           GET_RES_EXT()     // '.'              // SDOY

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DAILYAVERAGE: Reading ', a ) 

!      ! Get TAU value corresponding to DOY in year 2000
!      TAU = GET_TAU0( MONTH, DAY, 2000 )

      ! Read balanced biosphere CO2 [molec/cm2/s] from disk
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 3, 
     &                 TAU,       IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )

      END SUBROUTINE READ_BBIO_DAILYAVERAGE
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_bbio_diurnalcycle
!
! !DESCRIPTION: Subroutine READ\_BBIO\_DIURNALCYCLE reads CASA daily Net 
!  Ecosystem Production (NEP) fluxes but with a diurnal cycle imposed.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BBIO_DIURNALCYCLE( MONTH, DAY, HOUR, DOY )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D 
      USE TIME_MOD,      ONLY : GET_YEAR, ITS_A_LEAPYEAR

#     include "CMN_SIZE"      ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: MONTH   ! Current month (1-12)
      INTEGER, INTENT(IN) :: DAY     ! Current day (1-31)
      INTEGER, INTENT(IN) :: HOUR    ! Current hour (0-23)
      INTEGER, INTENT(IN) :: DOY     ! Current day of year (0-365)
!
! !REMARKS:
!  References
!  (1 ) Olsen and Randerson (2004), Differences between surface and column 
!        atmospheric CO2 and implications for carbon cycle research, J. 
!        Geophys. Res., 109, D02301,
!  (2 ) Potter et al. (1993), terrestrial Ecosystem Production: A process
!        model based on global satellite and surface data, Glob. Biogeochem.
!        Cycles, 7(4), 811-841.
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Added fixes for leapyears 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: YEAR
      TYPE (XPLEX)                 :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)                 :: TAU
      CHARACTER(LEN=3 )      :: SDOY
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_BBIO_DIURNALCYCLE begins here!
      !=================================================================
      
      YEAR   = GET_YEAR()    ! dbj

      ! Create string for day of year
      !-----------------------------------------------------------------
      ! Since 1985 is NOT a leap-year, we need to account for years that
      ! are leap years if the date is Feb 29 or after by subtracting one 
      !                                                     ! Ray Nassar
      !-----------------------------------------------------------------
      
      IF( ITS_A_LEAPYEAR() .AND. DOY >= 60 ) THEN
         WRITE( SDOY, '(i3.3)' ) DOY-1 
      ELSE
         WRITE( SDOY, '(i3.3)' ) DOY 
      ENDIF

      ! Get TAU of this month & day in "generic" year 1985
      IF( MONTH == 2 .AND. DAY == 29) THEN
         TAU = GET_TAU0( 2,     28,  1985, HOUR )
      ELSE
         TAU = GET_TAU0( MONTH, DAY, 1985, HOUR )
      ENDIF

!      ! Create string for day of year
!      WRITE( SDOY, '(i3.3)' ) DOY

      ! File name
      IF (SDOY == '189') THEN
         FILENAME = TRIM( DATA_DIR )               //
     &              'CO2_200508/BBIO_DIURNAL/nep.' // 
     &              GET_NAME_EXT_2D() // '.'       //  
     &              GET_RES_EXT()     // '.'       // 
     &              SDOY              // '.orig'
      ELSE
         FILENAME = TRIM( DATA_DIR )  //
     &              'CO2_200508/BBIO_DIURNAL/nep.' // 
     &              GET_NAME_EXT_2D() // '.'       //  
     &              GET_RES_EXT()     // '.'       // 
     &              SDOY
      ENDIF

      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DIURNALCYCLE: Reading ', a )

      !-----------------------------------------------------------------
      ! Read Net Ecosytem Productivity [molec CO/cm2/s] from disk
      ! The CASA fluxes use atmospheric convention: 
      ! positive = into atm; negative = into biosphere
      !-----------------------------------------------------------------
      CALL READ_BPCH2( FILENAME, 'GLOB-NPP', 2, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
       
      ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )
         
      END SUBROUTINE READ_BBIO_DIURNALCYCLE
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_biomass_tg
!
! !DESCRIPTION: Subroutine TOTAL\_BIOMASS\_Tg prints the amount of biomass 
!  burning emissions that are emitted each month in Tg or Tg 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_BIOMASS_Tg( BBARRAY, MOLWT, NAME )
!
! !USES:
!
      USE GRID_MOD, ONLY : GET_AREA_CM2

#     include "CMN_SIZE"  ! Size parameters
!
! !INPUT PARAMETERS:
!
      TYPE (XPLEX),           INTENT(IN) :: MOLWT                 ! Mol wt [kg/mole]
      CHARACTER(LEN=*), INTENT(IN) :: NAME                  ! Species name 
      TYPE (XPLEX),           INTENT(IN) :: BBARRAY(IIPAR,JJPAR)  ! BB Emissions
                                                            ! [molec/cm2/month]
!
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: I, J
      TYPE (XPLEX)                       :: TOTAL, A_CM2
      CHARACTER(LEN=6)             :: UNIT

      !=================================================================
      ! TOTAL_BIOMASS_TG begins here!
      !=================================================================

      ! Initialize summing variable
      TOTAL = 0d0

      ! Convert from [molec  /cm2/month] to [kg  /month]
      ! or      from [molec C/cm2/month] to [kg C/month]
      DO J = 1, JJPAR
         A_CM2 = GET_AREA_CM2( J ) ! Grid box surface area [cm2]

         DO I = 1, IIPAR
            TOTAL = TOTAL + BBARRAY(I,J) * A_CM2 * ( MOLWT / 6.023d23 )
         ENDDO
      ENDDO
     
      ! Convert from kg --> Tg
      TOTAL = TOTAL * 1d-9

      ! Define unit string
      IF ( NAME == 'NOx' .or. NAME == 'CO' .or. NAME == 'CH2O' ) THEN
         UNIT = '[Tg  ]'
      ELSE
         UNIT = '[Tg C]'
      ENDIF

      ! Write totals
      WRITE( 6, 100 ) NAME, TOTAL, UNIT
 100  FORMAT( 'Sum Biomass ', a4, 1x, ': ', f9.3, 1x, a9  )

      END SUBROUTINE TOTAL_BIOMASS_TG
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: def_biosph_co2_regions_f
!
! !DESCRIPTION: Subroutine DEF\_BIOSPH\_CO2\_REGIONS defines the land 
!  biospheric and ocean CO2 exchange regions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DEF_BIOSPH_CO2_REGIONS_F( REGION )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)   :: REGION(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%   BUYER BEWARE! Tagged CO2 tracers only work for 2 x 2.5 grid!   %%%
!  %%%   Someone will have to make this more general later on...        %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, IOS
      INTEGER                :: TMP(IIPAR,JJPAR)
      INTEGER                :: LAND_REG(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: LANDFILE
      CHARACTER(LEN=144)     :: ROW
      CHARACTER(LEN=1)       :: CHAR1(IIPAR,JJPAR)

      !=================================================================
      ! Reading LAND BIOSPHERE REGIONS
      !=================================================================
      
      LANDFILE  = 'Regions_land.dat'

      WRITE(*,*) ' '
 100  FORMAT( '     - READ_REGIONS: Reading ', a ) 
      WRITE( 6, 100 ) TRIM( LANDFILE )

      ! Initialize ARRAY
      LAND_REG = 0

      ! Open file
      OPEN( IU_FILE, FILE = TRIM( LANDFILE ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )      

      ! Read data
      DO J = 1, JJPAR
         IF (IIPAR ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
         IF (IIPAR == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
         WRITE (*,'(A)') ROW

         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

         DO I = 1, IIPAR
            CHAR1(I,J) = ROW(I:I)
            IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
            READ (CHAR1(I,J),'(I1)') TMP(I,J)
         ENDDO
      ENDDO

       ! Close file
      CLOSE( IU_FILE )

	! Flip array in the North-South Direction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         LAND_REG(I,J) = TMP(I,JJPAR-J+1)
      ENDDO
      ENDDO
      WRITE(*,*) ' '

      !=================================================================
      ! Loop over entire globe -- multiprocessor
      !=================================================================

      DO J = 1, JJPAR
      DO I = 1, IIPAR
!-----------------------------------------------------------------------
         ! Tracer #13 -- Canadian Tundra       
         IF (LAND_REG(I,J) == 1 .and. I >= 5 .and. I <= 60) THEN
            REGION(I,J) = 13
!----------------------------------------------------------------------
         ! Tracer #14 -- NA Boreal Forest                 
         ELSE IF (LAND_REG(I,J) == 2 .and. I <= 60) THEN
            REGION(I,J) = 14
!-----------------------------------------------------------------------
         ! Tracer #15 -- Western US/Mexico                 
         ELSE IF (LAND_REG(I,J) == 3 .and. I <= 60) THEN
            REGION(I,J) = 15
!-----------------------------------------------------------------------
         ! Tracer #16 -- Central NA Agricultural  
         ELSE IF (LAND_REG(I,J) == 4 .and. I <= 60) THEN
            REGION(I,J) = 16
!-----------------------------------------------------------------------      
         ! Tracer #17 -- NA Mixed Forest 
         ELSE IF (LAND_REG(I,J) == 5 .and. I <= 60) THEN
            REGION(I,J) = 17
!-----------------------------------------------------------------------      
         ! Tracer #18 -- Central America and Caribbean
         ELSE IF (LAND_REG(I,J) == 6 .and. I <= 60) THEN
            REGION(I,J) = 18
!-----------------------------------------------------------------------      
         ! Tracer #19 -- SA Tropical Rain Forest 
         ELSE IF (LAND_REG(I,J) == 7 .and. I <= 60) THEN
            REGION(I,J) = 19
!-----------------------------------------------------------------------      
         ! Tracer #20 -- SA Coast and Mountains
         ELSE IF (LAND_REG(I,J) == 8 .and. I <= 60) THEN
            REGION(I,J) = 20
!-----------------------------------------------------------------------      
         ! Tracer #21 -- SA Wooded Grasslands
         ELSE IF (LAND_REG(I,J) == 9 .and. I <= 60) THEN
            REGION(I,J) = 21
!-----------------------------------------------------------------------      
         ! Tracer #22 -- Eurasian Tundra
         ELSE IF (LAND_REG(I,J) == 1 .and. (I>60 .or. I<=5)) THEN
            REGION(I,J) = 22
!-----------------------------------------------------------------------      
         ! Tracer #23 -- Eurasian Boreal Coniferous Forest
         ELSE IF (LAND_REG(I,J) == 2 .and. I > 60 .and. J > 65) THEN
            REGION(I,J) = 23
!-----------------------------------------------------------------------      
         ! Tracer #24 -- Eurasian Boreal Deciduous Forest
         ELSE IF (LAND_REG(I,J) == 5 .and. I > 60 .and. J > 65) THEN
            REGION(I,J) = 24
!-----------------------------------------------------------------------      
         ! Tracer #25 -- South and Central Europe
         ELSE IF (LAND_REG(I,J) == 6 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 25
!-----------------------------------------------------------------------      
         ! Tracer #26 -- Central Asian Grasslands
         ELSE IF (LAND_REG(I,J) == 4 .and. I > 60 .and. J > 46) THEN
            REGION(I,J) = 26
!-----------------------------------------------------------------------      
         ! Tracer #27 -- Central Asian Desert
         ELSE IF (LAND_REG(I,J) == 8 .and. I >100 .and. I <117) THEN
            REGION(I,J) = 27
!-----------------------------------------------------------------------      
         ! Tracer #28 -- East Asia Mainland
         ELSE IF (LAND_REG(I,J) == 3 .and. I > 100) THEN
            REGION(I,J) = 28
!-----------------------------------------------------------------------      
         ! Tracer #29 -- Japan
         ELSE IF (LAND_REG(I,J) == 9 .and. I > 100) THEN
            REGION(I,J) = 29
!-----------------------------------------------------------------------      
         ! Tracer #30 -- Northern African Desert
         ELSE IF (LAND_REG(I,J) == 8 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 30
!-----------------------------------------------------------------------      
         ! Tracer #31 -- Northern Africa Grasslands
         ELSE IF (LAND_REG(I,J) == 3 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 31
!-----------------------------------------------------------------------      
         ! Tracer #32 -- Africa Tropical Forest
         ELSE IF (LAND_REG(I,J) == 7 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 32
!-----------------------------------------------------------------------      
         ! Tracer #33 -- Southern Africa Grasslands 
         ELSE IF (LAND_REG(I,J) == 4 .and. I > 60 .and. J < 50) THEN
            REGION(I,J) = 33
!-----------------------------------------------------------------------      
         ! Tracer #34 -- Southern African Desert
         ELSE IF (LAND_REG(I,J) == 9 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 34
!-----------------------------------------------------------------------      
         ! Tracer #35 -- Middle East
         ELSE IF (LAND_REG(I,J) == 2 .and. J > 40 .and. J < 65) THEN
            REGION(I,J) = 35
!-----------------------------------------------------------------------      
         ! Tracer #36 -- India and bordering countries
         ELSE IF (LAND_REG(I,J) == 5 .and. I > 60 .and. J < 65) THEN
            REGION(I,J) = 36
!-----------------------------------------------------------------------      
         ! Tracer #37 -- Maritime Asia (Indonesia, Malaysia, New Guinea, etc.)
         ELSE IF (LAND_REG(I,J) == 7 .and. I > 100) THEN
            REGION(I,J) = 37
!-----------------------------------------------------------------------      
         ! Tracer #38 -- Australian Forest/Grassland
         ELSE IF (LAND_REG(I,J) == 6 .and. I > 100) THEN
            REGION(I,J) = 38
!-----------------------------------------------------------------------      
         ! Tracer #39 -- Australian Desert
         ELSE IF (LAND_REG(I,J) == 8 .and. I >116 .and. J < 46) THEN
            REGION(I,J) = 39
!-----------------------------------------------------------------------      
         ! Tracer #40 -- New Zealand
         ELSE IF (LAND_REG(I,J) == 2 .and. I > 120) THEN
            REGION(I,J) = 40
!-----------------------------------------------------------------------      
         ! Tracer #52 -- CO2 from everywhere else (Remote Islands & Ice Caps)
         ELSE
            REGION(I,J) = 52
!-----------------------------------------------------------------------      
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE DEF_BIOSPH_CO2_REGIONS_F
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: def_ocean_co2_regions_f
!
! !DESCRIPTION: Subroutine DEF\_OCEAN\_CO2\_REGIONS defines CO2 regions 
!  for ocean exchange.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DEF_OCEAN_CO2_REGIONS_F( REGION )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)  :: REGION(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%   BUYER BEWARE! Tagged CO2 tracers only work for 2 x 2.5 grid!   %%%
!  %%%   Someone will have to make this more general later on...        %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER               :: I, J, IOS
      INTEGER               :: TMP(IIPAR,JJPAR), OCEAN_REG(IIPAR,JJPAR)
      CHARACTER(LEN=255)    :: OCEANFILE
      CHARACTER(LEN=144)    :: ROW
      CHARACTER(LEN=1)      :: CHAR1(IIPAR,JJPAR)

      !=================================================================
      ! DEF_CO2_OCEAN_REGIONS begins here!
      !=================================================================
      
      OCEANFILE = 'Regions_ocean.dat'

      WRITE( 6, 100 ) TRIM( OCEANFILE )
 100  FORMAT( '     - READ_REGIONS: Reading ', a ) 
	WRITE(*,*) ' '

      ! Initialize ARRAYS
      OCEAN_REG = 0

      ! Open file
      OPEN( IU_FILE, FILE = TRIM( OCEANFILE ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )      

      ! Read data
      DO J = 1, JJPAR
         IF (IIPAR ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
         IF (IIPAR == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
         WRITE (*,'(A)') ROW

         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

         DO I = 1, IIPAR
            CHAR1(I,J) = ROW(I:I)
            IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
            READ (CHAR1(I,J),'(I1)') TMP(I,J)
         ENDDO
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

	! Flip array in the North-South Direction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCEAN_REG(I,J) = TMP(I,JJPAR-J+1)
      ENDDO
      ENDDO
      WRITE(*,*) ' '

      !=================================================================
      ! Loop over entire globe -- multiprocessor
      !=================================================================

      DO J = 1, JJPAR
      DO I = 1, IIPAR
!-----------------------------------------------------------------------
         ! Tracer #41 -- Arctic Ocean       
         IF       (OCEAN_REG(I,J) == 5 .and. J > 60) THEN
            REGION(I,J) = 41
!-----------------------------------------------------------------------
         ! Tracer #42 -- North Pacific                 
         ELSE IF  (OCEAN_REG(I,J) == 1) THEN
            REGION(I,J) = 42
!-----------------------------------------------------------------------
         ! Region #43 -- Tropical West Pacific                 
         ELSE IF  (OCEAN_REG(I,J) == 2) THEN
            REGION(I,J) = 43
!-----------------------------------------------------------------------
         ! Tracer #44 -- Tropical East Pacific 
         ELSE IF  (OCEAN_REG(I,J) == 3) THEN
            REGION(I,J) = 44
!-----------------------------------------------------------------------      
         ! Tracer #45-- South Pacific 
         ELSE IF  (OCEAN_REG(I,J) == 4) THEN
            REGION(I,J) = 45
!-----------------------------------------------------------------------      
         ! Tracer #46 -- North Atlantic
         ELSE IF  (OCEAN_REG(I,J) == 6 .and. J > 45) THEN
            REGION(I,J) = 46
!-----------------------------------------------------------------------      
         ! Tracer #47 -- Tropical Atlantic
         ELSE IF  (OCEAN_REG(I,J) == 7) THEN
            REGION(I,J) = 47
!-----------------------------------------------------------------------      
         ! Tracer #48 -- South Atlantic
         ELSE IF  (OCEAN_REG(I,J) == 8) THEN
            REGION(I,J) = 48
!-----------------------------------------------------------------------      
         ! Tracer #49 -- Tropical Indian Ocean
         ELSE IF  (OCEAN_REG(I,J) == 5 .and. J < 60) THEN
            REGION(I,J) = 49
!-----------------------------------------------------------------------      
         ! Tracer #50 -- Southern Indian Ocean
         ELSE IF  (OCEAN_REG(I,J) == 6 .and. J < 45) THEN
            REGION(I,J) = 50
!-----------------------------------------------------------------------      
         ! Tracer #51 -- Southern (Antacrtic) Ocean
         ELSE IF  (OCEAN_REG(I,J) == 9) THEN
            REGION(I,J) = 51
!-----------------------------------------------------------------------      
         ! Tracer #52 -- CO2 from everywhere else (Remote Islands & Ice Caps)
         ELSE
            REGION(I,J) = 52
!-----------------------------------------------------------------------      
		ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE DEF_OCEAN_CO2_REGIONS_F
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: def_fossil_co2_regions_f
!
! !DESCRIPTION:  Subroutine DEF\_FOSSIL\_CO2\_REGIONS defines CO2 regions 
!  for anthropogenic emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DEF_FOSSIL_CO2_REGIONS_F( REGION )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)   :: REGION(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%   BUYER BEWARE! Tagged CO2 tracers only work for 2 x 2.5 grid!   %%%
!  %%%   Someone will have to make this more general later on...        %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, IOS
      INTEGER                :: TMP(IIPAR,JJPAR), REG_CODE(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=144)     :: ROW
      CHARACTER(LEN=1)       :: CHAR1(IIPAR,JJPAR)

      !=================================================================
      ! DEF_CO2_FOSSIL_REGIONS begins here!
      !=================================================================
      
      FILENAME  = 'Regions_land.dat'

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_REGIONS: Reading ', a ) 

      ! Initialize ARRAYS
      REG_CODE = 0

      ! Open file
      OPEN( IU_FILE, FILE = TRIM( FILENAME ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )      

      ! Read data
      DO J = 1, JJPAR
         IF (IIPAR ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
         IF (IIPAR == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
         WRITE (*,'(A)') ROW
                
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

         DO I = 1, IIPAR
            CHAR1(I,J) = ROW(I:I)
            IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
            READ (CHAR1(I,J),'(I1)') TMP(I,J)
         ENDDO
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Flip array in the North-South Direction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         REG_CODE(I,J) = TMP(I,JJPAR-J+1)
      ENDDO
      ENDDO

      !=================================================================
      ! Loop over entire globe -- multiprocessor
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
!-----------------------------------------------------------------------
         ! Tracer #13 -- Canadian Tundra       
         IF (REG_CODE(I,J) == 1 .and. I >= 5 .and. I <= 60) THEN
            REGION(I,J) = 13
!-----------------------------------------------------------------------
         ! Tracer #14 -- NA Boreal Forest                 
         ELSE IF (REG_CODE(I,J) == 2 .and. I <= 60) THEN
            REGION(I,J) = 14
!-----------------------------------------------------------------------
         ! Tracer #15 -- Western US/Mexico                 
         ELSE IF (REG_CODE(I,J) == 3 .and. I <= 60) THEN
            REGION(I,J) = 15
!-----------------------------------------------------------------------
         ! Tracer #16 -- Central NA Agricultural  
         ELSE IF (REG_CODE(I,J) == 4 .and. I <= 60) THEN
            REGION(I,J) = 16
!-----------------------------------------------------------------------      
         ! Tracer #17 -- NA Mixed Forest 
         ELSE IF (REG_CODE(I,J) == 5 .and. I <= 60) THEN
            REGION(I,J) = 17
!-----------------------------------------------------------------------      
         ! Tracer #18 -- Central America and Caribbean
         ELSE IF (REG_CODE(I,J) == 6 .and. I <= 60) THEN
            REGION(I,J) = 18
!-----------------------------------------------------------------------      
         ! Tracer #19 -- SA Tropical Rain Forest 
         ELSE IF (REG_CODE(I,J) == 7 .and. I <= 60) THEN
            REGION(I,J) = 19
!-----------------------------------------------------------------------      
         ! Tracer #20 -- SA Coast and Mountains
         ELSE IF (REG_CODE(I,J) == 8 .and. I <= 60) THEN
            REGION(I,J) = 20
!-----------------------------------------------------------------------      
         ! Tracer #21 -- SA Wooded Grasslands
         ELSE IF (REG_CODE(I,J) == 9 .and. I <= 60) THEN
            REGION(I,J) = 21
!-----------------------------------------------------------------------      
         ! Tracer #22 -- Eurasian Tundra
         ELSE IF (REG_CODE(I,J) == 1 .and. (I>60 .or. I<=5)) THEN
            REGION(I,J) = 22
!-----------------------------------------------------------------------      
         ! Tracer #23 -- Eurasian Boreal Coniferous Forest
         ELSE IF (REG_CODE(I,J) == 2 .and. I > 60 .and. J > 55) THEN
            REGION(I,J) = 23
!-----------------------------------------------------------------------      
         ! Tracer #24 -- Eurasian Boreal Deciduous Forest
         ELSE IF (REG_CODE(I,J) == 5 .and. I > 60 .and. J > 64) THEN
            REGION(I,J) = 24
!-----------------------------------------------------------------------      
         ! Tracer #25 -- South and Central Europe
         ELSE IF (REG_CODE(I,J) == 6 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 25
!-----------------------------------------------------------------------      
         ! Tracer #26 -- Central Asian Grasslands
         ELSE IF (REG_CODE(I,J) == 4 .and. I > 60 .and. J > 46) THEN
            REGION(I,J) = 26
!-----------------------------------------------------------------------      
         ! Tracer #27 -- Central Asian Desert
         ELSE IF (REG_CODE(I,J) == 8 .and. I >100 .and. I <117) THEN
            REGION(I,J) = 27
!-----------------------------------------------------------------------      
         ! Tracer #28 -- East Asia Mainland
         ELSE IF (REG_CODE(I,J) == 3 .and. I > 100) THEN
            REGION(I,J) = 28
!-----------------------------------------------------------------------      
         ! Tracer #29 -- Japan
         ELSE IF (REG_CODE(I,J) == 9 .and. I > 100) THEN
            REGION(I,J) = 29
!-----------------------------------------------------------------------      
         ! Tracer #30 -- Northern African Desert
         ELSE IF (REG_CODE(I,J) == 8 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 30
!-----------------------------------------------------------------------      
         ! Tracer #31 -- Northern Africa Grasslands
         ELSE IF (REG_CODE(I,J) == 3 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 31
!-----------------------------------------------------------------------      
         ! Tracer #32 -- Africa Tropical Forest
         ELSE IF (REG_CODE(I,J) == 7 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 32
!-----------------------------------------------------------------------      
         ! Tracer #33 -- Southern Africa Grasslands 
         ELSE IF (REG_CODE(I,J) == 4 .and. I > 60 .and. J < 50) THEN
            REGION(I,J) = 33
!-----------------------------------------------------------------------      
         ! Tracer #34 -- Southern African Desert
         ELSE IF (REG_CODE(I,J) == 9 .and. I > 60 .and. I <100) THEN
            REGION(I,J) = 34
!-----------------------------------------------------------------------      
         ! Tracer #35 -- Middle East
         ELSE IF (REG_CODE(I,J) == 2 .and. J > 40 .and. J < 60) THEN
            REGION(I,J) = 35
!-----------------------------------------------------------------------      
         ! Tracer #36 -- India and bordering countries
         ELSE IF (REG_CODE(I,J) == 5 .and. I > 60 .and. J < 64) THEN
            REGION(I,J) = 36
!-----------------------------------------------------------------------
         ! Tracer #37 -- Maritime Asia (Indonesia, Malaysia, New Guinea, etc.)
         ELSE IF (REG_CODE(I,J) == 7 .and. I > 100) THEN
            REGION(I,J) = 37
!-----------------------------------------------------------------------      
         ! Tracer #38 -- Australian Forest/Grassland
         ELSE IF (REG_CODE(I,J) == 6 .and. I > 100) THEN
            REGION(I,J) = 38
!-----------------------------------------------------------------------      
         ! Tracer #39 -- Australian Desert
         ELSE IF (REG_CODE(I,J) == 8 .and. I > 116 .and. J <45) THEN
            REGION(I,J) = 39
!-----------------------------------------------------------------------      
         ! Tracer #40 -- New Zealand
         ELSE IF (REG_CODE(I,J) == 2 .and. I > 120) THEN
            REGION(I,J) = 40
!-----------------------------------------------------------------------      
         ! Tracer #52 -- CO2 from everywhere else (Remote Islands & Ice Caps)
         ELSE
            REGION(I,J) = 52
!-----------------------------------------------------------------------      
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE DEF_FOSSIL_CO2_REGIONS_F
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_co2
!
! !DESCRIPTION: Subroutine INIT\_CO2 allocates memory to module arrays and 
!  reads in annual mean emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CO2 
!
! !USES:
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LGENFF,  LANNFF,   LMONFF, LSTREETS
      USE LOGICAL_MOD, ONLY : LSEASBB, LGFED2BB, L8DAYBB, LBIOFUEL
      USE LOGICAL_MOD, ONLY : LBIODAILY,   LBIODIURNAL
      USE LOGICAL_MOD, ONLY : LBIONETORIG, LBIONETCLIM
      USE LOGICAL_MOD, ONLY : LOCN1997,    LOCN2009ANN, LOCN2009MON
      USE LOGICAL_MOD, ONLY : LFFBKGRD
      USE LOGICAL_MOD, ONLY : LSHIPEDG,    LSHIPICO,    LPLANE
      USE LOGICAL_MOD, ONLY : LBIOSPHTAG,  LFOSSILTAG
      USE LOGICAL_MOD, ONLY : LSHIPTAG,    LPLANETAG
      USE TRACER_MOD,  ONLY : N_TRACERS

#     include "CMN_SIZE"
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_CO2 begins here!
      !=================================================================
           
      ! Exit if we have already intialized 
      IF ( IS_INIT ) RETURN
	
      ! Array for Fossil fuel CO2
      ALLOCATE( EMFOSSCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMFOSSCO2' )         
      EMFOSSCO2 = 0d0 

      ! Array for CO2 from ocean exchange
      ALLOCATE( EMOCCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMOCCO2' )         
      EMOCCO2 = 0d0 

      ! Array for Balanced Bio CO2
      ALLOCATE( EMBIOCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOCO2' )         
      EMBIOCO2 = 0d0 

      ! Array for Biomass burning CO2
      ALLOCATE( EMBIOBRNCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOBRNCO2' )
      EMBIOBRNCO2 = 0d0

      ! Array for Biofuel CO2
      ALLOCATE( EMBIOFUELCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOFUELCO2' )         
      EMBIOFUELCO2 = 0d0 

      ! Array for NET BIO CO2
      ALLOCATE( EMBIONETCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIONETCO2' )
      EMBIONETCO2  = 0d0

      ! Array for Ship CO2
      ALLOCATE( EMSHIPCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMSHIPCO2' )         
      EMSHIPCO2 = 0d0 

      ! Array for Aircraft CO2
      ALLOCATE( EMPLANECO2( IIPAR, JJPAR, LLPAR), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMPLANECO2' )         
      EMPLANECO2 = 0d0 

      ! Array for chemical source of CO2   ! dbj
      ALLOCATE( CHEMCO2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHEMCO2' )
      CHEMCO2 = 0d0

      ! Array for compensation to chemical source !dbj
      ALLOCATE( EMIS_SUB( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_SUB' )
      EMIS_SUB = 0d0

      ! Array for Fossil Fuel regions
      ALLOCATE( FOSSIL_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FOSSIL_REGION' )
      FOSSIL_REGION = 0

      ! Array for Biospheric regions
      ALLOCATE( BIOSPH_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOSPH_REGION' )
      BIOSPH_REGION = 0

      ! Array for Ocean Regions
      ALLOCATE( OCEAN_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCEAN_REGION' )
      OCEAN_REGION = 0

      !=================================================================
      ! Read in annual mean emissions
      !=================================================================

      ! Biofuel emissions
      IF (LBIOFUEL) CALL READ_ANNUAL_BIOFUELCO2

      ! Net terrestrial exchange
      IF (LBIONETORIG .OR. LBIONETCLIM) CALL READ_ANNUAL_BIONET_CO2

      ! Reset IS_INIT flag
      IS_INIT = .TRUE.
      
      END SUBROUTINE INIT_CO2 
!EOC
!------------------------------------------------------------------------------
!                        University of Toronto and                            !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_co2
!
! !DESCRIPTION: Subroutine CLEANUP\_CO2 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CO2 
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_CO2 begins here!
      !=================================================================
      IF ( ALLOCATED( EMFOSSCO2    ) ) DEALLOCATE( EMFOSSCO2    )
      IF ( ALLOCATED( EMOCCO2      ) ) DEALLOCATE( EMOCCO2      )
      IF ( ALLOCATED( EMBIOCO2     ) ) DEALLOCATE( EMBIOCO2     )
      IF ( ALLOCATED( EMBIOBRNCO2  ) ) DEALLOCATE( EMBIOBRNCO2  )
      IF ( ALLOCATED( EMBIOFUELCO2 ) ) DEALLOCATE( EMBIOFUELCO2 )
      IF ( ALLOCATED( EMBIONETCO2  ) ) DEALLOCATE( EMBIONETCO2  )
      IF ( ALLOCATED( EMSHIPCO2    ) ) DEALLOCATE( EMSHIPCO2    )
      IF ( ALLOCATED( EMPLANECO2   ) ) DEALLOCATE( EMPLANECO2   )
      IF ( ALLOCATED( CHEMCO2      ) ) DEALLOCATE( CHEMCO2      )
      IF ( ALLOCATED( EMIS_SUB     ) ) DEALLOCATE( EMIS_SUB     )

      END SUBROUTINE CLEANUP_CO2
!EOC
      END MODULE CO2_MOD
