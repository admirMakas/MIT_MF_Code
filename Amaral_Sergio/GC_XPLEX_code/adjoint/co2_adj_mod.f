! $Id: co2_adj_mod.f,v 1.3 2012/03/01 22:00:26 daven Exp $
      MODULE CO2_ADJ_MOD
!
!******************************************************************************
!  Module CO2_ADJ_MOD contains variables and routines used for the CO2 
!  adjoint simulation. (dkh, 04/25/10) 
! 
!  Based on the forward module  (pns, bmy, 8/16/05, 9/27/06) 
!
!  Module Variables:
!  ============================================================================
!
!  Module Procedures:
!  ============================================================================
!  (1 ) EMISSCO2_ADJ           : Adjoint of emits CO2 into individual tracers
!
!  GEOS-CHEM modules referenced by "co2_mod.f"
!  ============================================================================
!  (1 ) biomass_mod.f          : Module w/ routines for biomass burning
!  (2 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (3 ) diag04_mod.f           : Module w/ routines for CO2 diagnostics
!  (4 ) directory_mod.f        : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (6 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (7 ) grid_mod.f             : Module w/ horizontal grid information
!  (8 ) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (9 ) time_mod.f             : Module w/ routines for computing time & date
!  (10) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.
!  (11) transfer_mod.f         : Module w/ routines to cast & resize arrays 
!  (12) dao_mod.f              : Module w/ routines for working with DAO met fields 
!  (13) regrid_1x1_mod.f       : Modele w/ routines for regridding to and from 1x1
!
!  NOTES: 
!  (1 ) See forward model module for complete documentation  
!******************************************************************************
!
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "co2_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 

      ! ... except these routines
      PUBLIC :: EMISSCO2_ADJ

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!-----------------------------------------------------------------------------
   
      SUBROUTINE EMISSCO2_ADJ
!
!******************************************************************************
!  Subroutine EMISSCO2_ADJ is the adjoint routine for CO2 emissions. (dkh, 04/25/10) 
!
!  Based on forward model code.  (pns, bmy, 8/16/05, 9/27/06)
!
!  NOTES:
!
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,   ONLY : BIOMASS,       IDBCO2
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TIME_MOD,      ONLY : GET_DAY,       GET_DAY_OF_YEAR
      USE TIME_MOD,      ONLY : GET_HOUR,      GET_MONTH
      USE TIME_MOD,      ONLY : GET_YEAR,      GET_TS_CHEM, GET_TS_EMIS 
      USE TIME_MOD,      ONLY : ITS_A_NEW_DAY, ITS_A_NEW_MONTH
      USE TRACER_MOD,    ONLY : N_TRACERS,     STT

      USE LOGICAL_MOD, ONLY : LGENFF, LANNFF, LMONFF, LSTREETS
      USE LOGICAL_MOD, ONLY : LSEASBB, LGFED2BB, L8DAYBB, LBIOFUEL
      USE LOGICAL_MOD, ONLY : LBIODAILY, LBIODIURNAL
      USE LOGICAL_MOD, ONLY : LBIONETORIG, LBIONETCLIM
      USE LOGICAL_MOD, ONLY : LOCN1997, LOCN2009ANN, LOCN2009MON
      USE LOGICAL_MOD, ONLY : LSHIPEDG, LSHIPICO, LPLANE
      USE LOGICAL_MOD, ONLY : LBIOSPHTAG, LFOSSILTAG, LFFBKGRD
      USE LOGICAL_MOD, ONLY : LSHIPTAG, LPLANETAG
      USE LOGICAL_MOD, ONLY : LSHIPSCALE, LPLANESCALE
      USE LOGICAL_MOD, ONLY : LCHEMCO2
 
      ! adj_group
      USE CO2_MOD,         ONLY : READ_BBIO_DIURNALCYCLE
      USE CO2_MOD,         ONLY : READ_BBIO_DAILYAVERAGE
      USE CO2_MOD,         ONLY : READ_FOSSILCO2
      USE CO2_MOD,         ONLY : READ_OCEANCO2
      USE CO2_MOD,         ONLY : READ_SHIPCO2_EDGAR
      USE CO2_MOD,         ONLY : READ_SHIPCO2_ICOADS
      USE CO2_MOD,         ONLY : READ_AVIATION_CO2
      USE CO2_MOD,         ONLY : READ_CHEMCO2
      USE CO2_MOD,         ONLY : CHEM_SURF
      USE CO2_MOD,         ONLY : CHEMCO2
      USE CO2_MOD,         ONLY : EMFOSSCO2
      USE CO2_MOD,         ONLY : EMOCCO2
      USE CO2_MOD,         ONLY : EMBIOCO2
      USE CO2_MOD,         ONLY : EMBIOFUELCO2
      USE CO2_MOD,         ONLY : EMBIONETCO2
      USE CO2_MOD,         ONLY : EMSHIPCO2
      USE CO2_MOD,         ONLY : EMPLANECO2
      USE CO2_MOD,         ONLY : EMIS_SUB
      USE CO2_MOD,         ONLY : XNUMOL_CO2
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2ff,  IDADJ_ECO2ocn
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2bal, IDADJ_ECO2bb
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2bf,  IDADJ_ECO2nte
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2shp, IDADJ_ECO2pln
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ECO2che, IDADJ_ECO2sur
      USE ADJ_ARRAYS_MOD,  ONLY : GET_SCALE_GROUP
      USE ADJ_ARRAYS_MOD,  ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE TIME_MOD,        ONLY : ITS_A_NEW_DAY_ADJ
      USE TIME_MOD,        ONLY : GET_NHMSe
	
      ! dkh debug
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      LOGICAL, SAVE          :: FIRST  = .TRUE.

      ! Local variables
      INTEGER                :: I,     IJLOOP, J,    L,     N,   NN
      INTEGER                :: M
      INTEGER                :: DAY,   DOY,    HOUR, MONTH, YEAR   
      TYPE (XPLEX)                 :: A_CM2, DTSRCE, E_CO2
      TYPE (XPLEX)                 :: E_CO2_ADJ
      TYPE (XPLEX)                 :: biomass_sum, bionet_sum
      TYPE (XPLEX), SAVE           :: CHEMSRC(IIPAR,JJPAR,LLPAR)  ! dbj

      ! External functions     
      TYPE (XPLEX), EXTERNAL       :: BOXVL    ! dbj

      !=================================================================
      ! EMISSCO2_ADJ begins here!
      !=================================================================
      IF ( FIRST ) THEN
         CHEMSRC = CHEMCO2
         FIRST   = .FALSE. 
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

      M      = GET_SCALE_GROUP()
 
      ! adjust HOUR for adjoint
      IF (     MOD(GET_NHMSe(),3) == 0 ) THEN 
         HOUR = HOUR - 2
         IF ( HOUR < 0 ) HOUR = HOUR + 24
      ELSEIF ( MOD(GET_NHMSe(),2) == 0 ) THEN 
         HOUR = HOUR - 1
         IF ( HOUR < 0 ) HOUR = HOUR + 24
      ELSEIF ( MOD(GET_NHMSe(),1) == 0 ) THEN 
         HOUR = HOUR
         IF ( HOUR < 0 ) HOUR = HOUR + 24
      ENDIF 
         
!--------------------------------------------------------------------------------------
!      ! Read monthly-mean biomass burning emissions
!      IF ( LBIOBRNCO2 .and. ITS_A_NEW_MONTH() ) THEN 
!         CALL READ_MONTH_BIOBRN_CO2( MONTH, YEAR )
!      ENDIF
!  This requires a subroutine called READ_MONTH_BIOBRN_CO2 
!  GFEDv2 biomass burning emissions are a better choice !Ray Nassar
!
!--------------------------------------------------------------------------------------
!  At present, biomass burning emissions are dealt with in the following way:
!
!  1) main.f calls do_emissions in emissions.f
!  2) do_emissions calls compute_biomass_emissions in biomass_mod.f
!  3a) compute_biomass_emissions calls gfed2_compute_biomass in gfed2_biomass_mod.f
!               ** OR **
!  3b) compute_biomass_emissions calls gc_read_biomass_co2 in gc_biomass_mod.f
!--------------------------------------------------------------------------------------

!      ! Check if Balanced Biosphere emissions are required  
!      IF ( LBIOCO2 ) THEN  
!         ! If LUSECASANEP is TRUE ...
!         IF ( LUSECASANEP ) THEN

      IF ( LBIODIURNAL ) THEN

         write(*,*) '*** USING DIURNAL CASA NEP ***'

         ! ... then use 3-hourly NEP emissions for Bal Bio ...
         IF ( MOD( HOUR, 3 ) == 0 ) THEN
               CALL READ_BBIO_DIURNALCYCLE( MONTH, DAY, HOUR, DOY )
         ENDIF

         ! dkh debug
         !print*, 'CO2dbg adj HOUR ', HOUR, MOD( HOUR, 3 )

      ELSEIF ( LBIODAILY ) THEN

         ! ... otherwise use constant daily emissions of NEP for Bal Bio
         IF ( ITS_A_NEW_DAY_ADJ() ) THEN
            CALL READ_BBIO_DAILYAVERAGE( MONTH, DAY, DOY ) 
         ENDIF

      ENDIF

!-----------------------------------------------------------------------
! Fluxes with "possible" monthly variability are called below
! In some cases the annual file is just called at the start of the month
!-----------------------------------------------------------------------

      IF (ITS_A_NEW_MONTH() ) THEN 
	
      ! Fossil fuel emissions
        IF (LMONFF .OR. LANNFF .OR. LGENFF) CALL READ_FOSSILCO2 

      ! Oceanic exchange
      	IF (LOCN1997 .OR. LOCN2009ANN .OR. LOCN2009MON) THEN
		CALL READ_OCEANCO2
	ENDIF

      ! Ship emissions from EDGAR
      	IF (LSHIPEDG) CALL READ_SHIPCO2_EDGAR

      ! Ship emissions from ICOADS
      	IF (LSHIPICO) CALL READ_SHIPCO2_ICOADS

      ! Aircraft CO2 emissions
      	IF (LPLANE) CALL READ_AVIATION_CO2
		
      ! Get chemical source   ! dbj
      	IF (LCHEMCO2) THEN
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
!$OMP+PRIVATE( E_CO2_ADJ             )
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

      ! Loop over longitudes
      DO I = 1, IIPAR

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from fossil fuel emissions 
         !-------------------------------------------
         IF ( LGENFF .OR. LANNFF .OR. LMONFF) THEN

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Fossil fuel emissions of CO2 [molec/cm2/s]
            E_CO2          = EMFOSSCO2(I,J)
           
            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2ff > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ff)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) = 
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF

         !-------------------------------------------
         ! #3: CO2 from ocean exchange
         !-------------------------------------------
         IF (LOCN1997 .OR. LOCN2009ANN .OR. LOCN2009MON) THEN

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Ocean CO2 emissions in [molec/cm2/s]
            E_CO2          = EMOCCO2(I,J)

               ! dkh debug
               IF ( I == IFD .and. J == JFD ) THEN
               print*, ' ECO2onc adj = ', E_CO2
               print*, ' E_CO2_ADJ   = ', E_CO2_ADJ
               print*, ' STT_ADJ(CO2)= ', STT_ADJ(I,J,1,1)
               ENDIF

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2ocn > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ocn)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2ocn) = 
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2ocn) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF

         !-------------------------------------------
         ! #4: CO2 from balanced biosphere emissions
         !-------------------------------------------
         IF ( LBIODAILY .OR. LBIODIURNAL ) THEN

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Balanced biosphere CO2 [molec/cm2/s]
            E_CO2         = EMBIOCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2bal > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2bal)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2bal) = 
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2bal) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF

         !-------------------------------------------
         ! #5: CO2 from biomass burning emissions
         !-------------------------------------------
         IF ( LSEASBB .OR. LGFED2BB .OR. L8DAYBB ) THEN 

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Biomass burning emissions [molec/cm2/s]
            E_CO2          = BIOMASS(I,J,IDBCO2)
            !E_CO2          = EMBIOBRNCO2(I,J) !This was from older versions, see note above

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2bb > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2bb)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2bb) = 
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2bb) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF

         !-------------------------------------------
         ! #6: CO2 from biofuel emissions
         !-------------------------------------------
         IF ( LBIOFUEL ) THEN

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

           ! Biofuel CO2 emissions [molec/cm2/s]
            E_CO2          = EMBIOFUELCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2bf > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2bf)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2bf) = 
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2bf) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF

         !-------------------------------------------
         ! #7: CO2 from net terrestrial exchange
         !-------------------------------------------
         IF ( LBIONETORIG .OR. LBIONETCLIM ) THEN

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            ! CO2 from net terrestrial exchange [molec/cm2/s]
            E_CO2          = EMBIONETCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2nte > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2nte)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2nte) =
     &           EMS_SF_ADJ(I,J,M,IDADJ_ECO2nte) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF

         !-------------------------------------------
         ! #8: CO2 from ship emissions
         !-------------------------------------------
         IF ( LSHIPEDG .OR. LSHIPICO ) THEN

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,1,1)

            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Ship CO2 emissions [molec/cm2/s]
            E_CO2          = EMSHIPCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2shp > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2shp)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2shp) =
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2shp) + E_CO2_ADJ * E_CO2
            ENDIF 
              
         ENDIF 

         !-------------------------------------------
         ! #9: CO2 from aircraft emissions
         !-------------------------------------------
         IF ( LPLANE ) THEN
	   DO L = 1, LLPAR

            ! fwd code:
            !STT(I,J,L,1)   = STT(I,J,L,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,L,1)

            !--------------------------------------------------------
            ! BUG FIX: unit conversion (dkh, 02/08/12, adj32_018) 
            ! OLD CODE: 
            !! fwd code:
            !!E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            !! adj code:
            !!E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2
            ! NEW CODE: 
            ! fwd code:
            !E_CO2    = E_CO2 * BOXVL(I,J,L) * DTSRCE / XNUMOL_CO2
            ! adj code:
             E_CO2_ADJ      = E_CO2_ADJ * BOXVL(I,J,L) 
     &                      * DTSRCE    / XNUMOL_CO2
            !--------------------------------------------------------



            ! Aircraft CO2 emissions (3-D) [molec/cm3/s]
            E_CO2          = EMPLANECO2(I,J,L)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2pln > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2pln)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2pln) =
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2pln) + E_CO2_ADJ * E_CO2
            ENDIF 
              
           ENDDO
         ENDIF
 
         !-------------------------------------------
         ! #10 CO2 production from CO oxidation
         !-------------------------------------------
         IF ( LCHEMCO2 ) THEN
	   DO L = 1, LLPAR
         
            ! fwd code:
            !STT(I,J,L,1)   = STT(I,J,L,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = STT_ADJ(I,J,L,1)

            !--------------------------------------------------------
            ! BUG FIX: unit conversion (dkh, 02/08/12, adj32_018) 
            ! OLD CODE: 
            !! fwd code:
            !!E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            !! adj code:
            !E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2
            ! NEW CODE: 
            ! fwd code:
            !E_CO2  = E_CO2 * BOXVL(I,J,L) * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * BOXVL(I,J,L) 
     &                     * DTSRCE    / XNUMOL_CO2
            !--------------------------------------------------------


            E_CO2 = CHEMSRC(I,J,L)

            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( LADJ .and. IDADJ_ECO2che > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2che)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2che) =
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2che) + E_CO2_ADJ * E_CO2
            ENDIF 
           ENDDO 
         ENDIF 

         !-------------------------------------------
         ! #11 CO2 surface correction for CO oxidation
         !-------------------------------------------
         IF ( LCHEMCO2 ) THEN
      
            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) - E_CO2
            ! adj code:
            !E_CO2_ADJ = STT_ADJ(I,J,1,1)
            E_CO2_ADJ = - STT_ADJ(I,J,1,1)
            
            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2
   
            E_CO2 = EMIS_SUB(I,J)  ! EMIS_SUB is positive, but is subtracted
            
            ! adj_group: apply scaling factors (dkh, 04/25/10) 
            IF ( IDADJ_ECO2sur > 0 ) THEN 
               ! fwd code: 
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2sur)
               ! adj code:
               EMS_SF_ADJ(I,J,M,IDADJ_ECO2sur) =
     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2sur) + E_CO2_ADJ * E_CO2
            ENDIF


          ENDIF 

            ! dkh debug
          IF ( I == IFD .and. J == JFD .and. LPRINTFD ) THEN 
               print*, 'CO2dbg adj E_CO2 ff ', EMFOSSCO2(I,J)
               print*, 'CO2dbg adj E_CO2 ocn', EMOCCO2(I,J)
               print*, 'CO2dbg adj E_CO2 bal', EMBIOCO2(I,J)
               print*, 'CO2dbg adj E_CO2 bb ', BIOMASS(I,J,IDBCO2)
               print*, 'CO2dbg adj E_CO2 bf ', EMBIOFUELCO2(I,J)
               print*, 'CO2dbg adj E_CO2 nte', EMBIONETCO2(I,J)
               print*, 'CO2dbg adj E_CO2 shp', EMSHIPCO2(I,J)
               print*, 'CO2dbg adj E_CO2 pln', EMPLANECO2(I,J,LFD)
               print*, 'CO2dbg adj E_CO2 che', CHEMSRC(I,J,LFD)
               print*, 'CO2dbg adj E_CO2 sur', EMIS_SUB(I,J)
           ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

            ! dkh debug
               print*, 'CO2dbg adj E_CO2 ff ', SUM(EMFOSSCO2(:,:))
               print*, 'CO2dbg adj E_CO2 ocn', SUM(EMOCCO2(:,:))
               print*, 'CO2dbg adj E_CO2 bal', SUM(EMBIOCO2(:,:))
               print*, 'CO2dbg adj E_CO2 bb ', SUM(BIOMASS(:,:,IDBCO2))
               print*, 'CO2dbg adj E_CO2 bf ', SUM(EMBIOFUELCO2(:,:))
               print*, 'CO2dbg adj E_CO2 nte', SUM(EMBIONETCO2(:,:))
               print*, 'CO2dbg adj E_CO2 shp', SUM(EMSHIPCO2(:,:))
               print*, 'CO2dbg adj E_CO2 pln', SUM(EMPLANECO2(:,:,:))
               print*, 'CO2dbg adj E_CO2 che', SUM(CHEMSRC(:,:,:))
               print*, 'CO2dbg adj E_CO2 sur', SUM(EMIS_SUB(:,:))


      ! Return to calling program
      END SUBROUTINE EMISSCO2_ADJ

!-----------------------------------------------------------------------------


      ! End of module
      END MODULE CO2_ADJ_MOD
