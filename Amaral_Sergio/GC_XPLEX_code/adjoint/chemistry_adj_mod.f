! $Id: chemistry_adj_mod.f,v 1.9 2012/03/01 22:00:26 daven Exp $
      MODULE CHEMISTRY_ADJ_MOD
!
!******************************************************************************
!  Module CHEMISTRY_MOD is used to call the proper chemistry subroutine
!  for the various GEOS-CHEM simulations. (bmy, 4/14/03, 4/2/08)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CHEMISTRY       : Driver which calls various chemistry routines
!
!  GEOS-CHEM modules referenced by chemistry_mod.f
!  ============================================================================
!  (1 ) acetone_mod.f      : Module w/ routines for ACET chemistry
!  (2 ) c2h6_mod.f         : Module w/ routines for C2H6 chemistry
!  (3 ) carbon_mod.f       : Module w/ routines for carbon arsl chem.
!  (4 ) ch3i_mod.f         : Module w/ routines for CH3I chemistry
!  (5 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (6 ) diag_pl_mod.f      : Module w/ routines to save P(Ox), L(Ox)
!  (7 ) drydep_mod.f       : Module w/ GEOS-CHEM drydep routines
!  (8 ) dust_mod.f         : Module w/ routines for dust arsl chem.
!  (9 ) error_mod.f        : Module w/ NaN and error checks
!  (10) global_ch4_mod.f   : Module w/ routines for CH4 chemistry
!  (11) hcn_ch3cn_mod.f    : Module w/ routines for HCN and CH3CN chemistry
!  (12) Kr85_mod.f         : Module w/ routines for Kr85 chemistry
!  (13) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (14) RnPbBe_mod.f       : Module w/ routines for Rn-Pb-Be chemistry
!  (15) rpmares_mod.f      : Module w/ routines for arsl phase equilib.
!  (16) seasalt_mod.f      : Module w/ routines for seasalt chemistry
!  (17) sulfate_mod.f      : Module w/ routines for sulfate chemistry
!  (18) tagged_co_mod.f    : Module w/ routines for Tagged CO chemistry
!  (19) tagged_ox_mod.f    : Module w/ routines for Tagged Ox chemistry
!  (20) time_mod.f         : Module w/ routines to compute time & date
!  (21) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc. 
!  (22) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!
!  NOTES:
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Now references "tagged_ox_mod.f"(bmy, 8/18/03)
!  (4 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (5 ) Bug fix: Now also call OPTDEPTH for GEOS-4 (bmy, 1/27/04)
!  (6 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/5/04)
!  (7 ) Now references "seasalt_mod.f" (rjp, bec, bmy, 4/20/04)
!  (8 ) Now references "logical_mod.f", "tracer_mod.f", "diag20_mod.f", and
!        "diag65_mod.f", and "aerosol_mod." (bmy, 7/20/04)
!  (9 ) Now references "mercury_mod.f" (bmy, 12/7/04)
!  (10) Updated for SO4s, NITs chemistry (bec, bmy, 4/13/05)
!  (11) Now call CHEM_HCN_CH3CN from "hcn_ch3cn_mod.f".  Also remove all
!        references to the obsolete CO-OH param simulation. (xyp, bmy, 6/24/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Updated for SOA from isoprene (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) For now, replace use RPMARES instead of ISORROPIA. (bmy, 4/2/08)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_CHEMISTRY_ADJ
!
!******************************************************************************
!  Subroutine DO_CHEMISTRY is the driver routine which calls the appropriate
!  chemistry subroutine for the various GEOS-CHEM simulations. 
!  (bmy, 2/11/03, 9/18/07)
!
!  NOTES:
!  (1 ) Now reference DELP, T from "dao_mod.f" since we need to pass this
!        to OPTDEPTH for GEOS-1 or GEOS-STRAT met fields (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Removed call to CHEMO3, it's obsolete.  Now calls CHEM_TAGGED_OX !
!        from "tagged_ox_mod.f" when NSRCX==6.  Now calls Kr85 chemistry if 
!        NSRCX == 12 (jsw, bmy, 8/20/03)
!  (4 ) Bug fix: added GEOS-4 to the #if block in the call to OPTDEPTH.
!        (bmy, 1/27/04)
!  (5 ) Now calls CHEMCARBON and CHEMDUST to do carbon aerosol & dust 
!        aerosol chemistry (rjp, tdf, bmy, 4/2/04)
!  (6 ) Now calls CHEMSEASALT to do seasalt aerosol chemistry 
!        (rjp, bec, bmy, 4/20/04)
!  (7 ) Now references "logical_mod.f" & "tracer_mod.f".  Now references
!        AEROSOL_CONC, AEROSOL_RURALBOX, and RDAER from "aerosol_mod.f".  
!        Now includes "CMN_DIAG" and "comode.h".  Also call READER, READCHEM, 
!        and INPHOT to initialize the FAST-J arrays so that we can save out !
!        AOD's to the ND21 diagnostic for offline runs. (bmy, 7/20/04)
!  (8 ) Now call routine CHEMMERCURY from "mercury_mod.f" for an offline
!        Hg0/Hg2/HgP simulation. (eck, bmy, 12/7/04)
!  (9 ) Now do not call DO_RPMARES if we are doing an offline aerosol run
!        with crystalline sulfur & aqueous tracers (cas, bmy, 1/7/05)
!  (10) Now use ISOROPIA for aer thermodyn equilibrium if we have seasalt 
!        tracers defined, or RPMARES if not.  Now call CHEMSEASALT before
!        CHEMSULFATE.  Now do aerosol thermodynamic equilibrium before
!        aerosol chemistry for offline aerosol runs.  Now also reference 
!        CLDF from "dao_mod.f" (bec, bmy, 4/20/05)
!  (11) Now modified for GCAP met fields.  Now call CHEM_HCN_CH3CN from 
!        "hcn_ch3cn_mod.f".  Also remove allreferences to the obsolete 
!         CO-OH param simulation. (xyp, bmy, 6/23/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Removed ISOP_PRIOR as a local variable (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) Now use DRYFLXH2HD and CHEM_H2_HD for H2/HD sim (lyj, phs, 9/18/07)
!  (17) Bug fix: now hardwired to use RPMARES since ISORROPIA can return very
!        unphysical values at low RH.  Wait for ISORROPIA II. (bmy, 4/2/08)
!  (18) Updated to support offline BC aerosol (yhmao, dkh, 01/13/12, adj32_013) 
!  (19) Add support for CH4 (kjw, dkh, 02/12/12, adj32_023) 
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,     ONLY : OCEAN_SINK_ACET
      USE AEROSOL_MOD,     ONLY : AEROSOL_CONC, AEROSOL_RURALBOX
      USE AEROSOL_MOD,     ONLY : RDAER,        SOILDUST
      USE CARBON_MOD,      ONLY : CHEMCARBON
      USE C2H6_MOD,        ONLY : CHEMC2H6
      USE CH3I_MOD,        ONLY : CHEMCH3I
      USE DAO_MOD,         ONLY : CLDF,    DELP
      USE DAO_MOD,         ONLY : OPTDEP,  OPTD,   T
      USE DRYDEP_MOD,      ONLY : DRYFLX, DRYFLXRnPbBe, DRYFLXH2HD
      USE DUST_MOD,        ONLY : CHEMDUST, RDUST_ONLINE
      USE DUST_ADJ_MOD,    ONLY : CHEMDUST_ADJ
      USE ERROR_MOD,       ONLY : DEBUG_MSG
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
      USE H2_HD_MOD,       ONLY : CHEM_H2_HD
      USE HCN_CH3CN_MOD,   ONLY : CHEM_HCN_CH3CN
      USE ISOROPIA_MOD,    ONLY : DO_ISOROPIA
      USE Kr85_MOD,        ONLY : CHEMKr85
      USE LOGICAL_MOD,     ONLY : LCARB, LCHEM,  LCRYST, LDUST
      USE LOGICAL_MOD,     ONLY : LPRT,  LSSALT, LSULF,  LSOA
      USE MERCURY_MOD,     ONLY : CHEMMERCURY
      USE OPTDEPTH_MOD,    ONLY : OPTDEPTH
      USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
      USE RPMARES_MOD,     ONLY : DO_RPMARES
      USE SEASALT_MOD,     ONLY : CHEMSEASALT
      USE SULFATE_MOD,     ONLY : CHEMSULFATE
      USE TAGGED_OX_MOD,   ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
      USE TIME_MOD,        ONLY : ITS_TIME_FOR_CHEM
      USE TRACER_MOD,      ONLY : N_TRACERS,       STT  
      USE TRACER_MOD,      ONLY : ITS_A_C2H6_SIM
      USE TRACER_MOD,      ONLY : ITS_A_CH3I_SIM
      USE TRACER_MOD,      ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,      ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,      ONLY : ITS_A_HCN_SIM
      USE TRACER_MOD,      ONLY : ITS_A_MERCURY_SIM
      USE TRACER_MOD,      ONLY : ITS_A_RnPbBe_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGOX_SIM
      USE TRACER_MOD,      ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,      ONLY : ITS_NOT_COPARAM_OR_CH4
      USE TRACERID_MOD,    ONLY : IDTACET, IDTISOP

      ! adjoint modules:
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, LFD, NFD 
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ
      USE CARBON_ADJ_MOD,    ONLY : CHEMCARBON_ADJ
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_CHEM
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM
      USE LOGICAL_ADJ_MOD,   ONLY : LPRINTFD
      USE RPMARES_ADJ_MOD,   ONLY : DO_RPMARES_ADJ
      USE RPMARES_MOD,       ONLY : RECOMP_RPMARES
      USE SULFATE_ADJ_MOD,   ONLY : CHEMSULFATE_ADJ
      USE TAGGED_CO_ADJ_MOD, ONLY : CHEM_TAGGED_CO_ADJ
      ! lzh 12/08/2009  add adjoint for tagged ox simulation
      USE TAGGED_OX_ADJ_MOD, ONLY : CHEM_TAGGED_OX_ADJ
      USE GLOBAL_CH4_ADJ_MOD,ONLY : CHEMCH4_ADJ


#     include "CMN_SIZE"        ! Size parameters
#     include "CMN_DIAG"        ! NDxx flags
#     include "comode.h"        ! NPHOT

      ! Local variables
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: N_TROP

      !=================================================================
      ! DO_CHEMISTRY_ADJ begins here!
      !=================================================================

      ! Compute optical depths (except for CH4 simulation)
      IF ( .not. ITS_A_CH4_SIM() ) THEN
         CALL OPTDEPTH( LLPAR, CLDF, OPTDEP, OPTD )
      ENDIF

      !=================================================================
      ! If LADJ_CHEM=T then call the adjoint chemistry subroutines
      !=================================================================
      IF ( LADJ_CHEM ) THEN 

         !---------------------------------
         ! NOx-Ox-HC (w/ or w/o aerosols) 
         !---------------------------------
         IF ( ITS_A_FULLCHEM_SIM() ) THEN 

            ! Adjoint of remove acetone ocean sink (it is self-adjoint)
            IF ( IDTACET /= 0 ) THEN
               CALL OCEAN_SINK_ACET( STT_ADJ(:,:,1,IDTACET) )
            ENDIF

            ! Do carbonaceous aerosol chemistry
            IF ( LCARB ) CALL CHEMCARBON_ADJ

            ! Also do sulfate chemistry
            IF ( LSULF ) THEN

               ! Do aerosol thermodynamic equilibrium
               !------------------------------------------------------------
               ! Prior to 4/2/08:
               ! Bug fix: ISORROPIA can return very unphysical values when
               ! RH is very low.  We will replace the current version of
               ! ISORROPIA with ISORROPIA II.  In the meantime, we shall
               ! use RPMARES to do the ATE computations. (bmy, 4/2/08)
               !IF ( LSSALT ) THEN
               !
               !   ! ISOROPIA takes Na+, Cl- into account
               !   CALL DO_ISOROPIA
               !
               !ELSE

                  ! RPMARES does not take Na+, Cl- into account
                  IF ( LAERO_THERM ) THEN 
 
                     ! Recalculate intermediate values 
                     CALL RECOMP_RPMARES

                     ! Diagnostic
                     IF ( LPRINTFD ) THEN
                        WRITE(6,*) 'Before RPMARES_ADJ: STT_ADJ(FD) = ',
     &                              STT_ADJ(IFD,JFD,LFD,NFD)
                     ENDIF

                     ! Call adjoint aerosol thermodynamics routine
                     CALL DO_RPMARES_ADJ

              	  ENDIF 

               !ENDIF
               !------------------------------------------------------------

               ! Do sulfate chemistry
               CALL CHEMSULFATE_ADJ
               
            ENDIF

            ! Call SMVGEAR routines
            CALL CHEMDR_ADJ

            ! Do seasalt aerosol chemistry
            IF ( LSSALT ) print*, ' ADJ of CHEMSEASALT not supported'
!            IF ( LSSALT ) CALL CHEMSEASALT

            ! Do dust aerosol chemistry
            IF ( LDUST ) CALL CHEMDUST_ADJ

            ! ND44 drydep fluxes
!            CALL DRYFLX     

            ! ND43 chemical production
!            CALL DIAGOH

      !---------------------------------
      ! Offline aerosol simulation
      !---------------------------------
      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            ! Define loop index and other SMVGEAR arrays
            ! N_TROP, the # of trop boxes, is returned
            !CALL AEROSOL_RURALBOX( N_TROP )


            ! Initialize FAST-J quantities for computing AOD's
            IF ( FIRST ) THEN
               CALL READER( FIRST )
 
               ! turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
               !CALL READCHEM
               !CALL INPHOT( LLTROP, NPHOT )

               ! Reset NCS with NCSURBAN
               NCS     = NCSURBAN

               ! Reset NTLOOP and NTTLOOP after call to READER
               ! with the actual # of boxes w/in the ann mean trop
               NTLOOP  = N_TROP
               NTTLOOP = N_TROP

               ! Reset first-time flag
               FIRST = .FALSE.
            ENDIF

            !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            ! Compute aerosol & dust concentrations [kg/m3]
            ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
            !CALL AEROSOL_CONC

            !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            ! Compute AOD's and surface areas
            !CALL RDAER

            !*** AEROSOL THERMODYNAMIC EQUILIBRIUM ***
            !-------------------------------------------------------------
            ! Prior to 4/2/08:
            ! Bug fix: ISORROPIA can return very unphysical values when
            ! RH is very low.  We will replace the current version of
            ! ISORROPIA with ISORROPIA II.  In the meantime, we shall
            ! use RPMARES to do the ATE computations. (bmy, 4/2/08)
            !IF ( LSSALT ) THEN
            !
            !   ! ISOROPIA takes Na+, Cl- into account
            !   CALL DO_ISOROPIA
            !
            !ELSE

               !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
               ! RPMARES does not take Na+, Cl- into account
               ! (skip for crystalline & aqueous offline run)
               !IF ( .not. LCRYST ) CALL DO_RPMARES

            !ENDIF
            !-------------------------------------------------------------

            !*** SEASALT AEROSOLS ***
            !IF ( LSSALT ) CALL CHEMSEASALT
            IF ( LSSALT ) 
     &         CALL ERROR_STOP( ' need CHEMSEASALT_ADJ ',
     &                          ' chemistry_adj_mod.f'   )

            !*** SULFATE AEROSOLS ***
            IF ( LSULF .or. LCRYST ) THEN

               ! Do sulfate chemistry
               !CALL CHEMSULFATE
               CALL ERROR_STOP( ' need CHEMSULFATE_ADJ ',
     &                          ' chemistry_adj_mod.f'   )

            ENDIF
               
            !*** CARBON AND 2NDARY ORGANIC AEROSOLS ***
            ! (yhmao, dkh, 01/13/12, adj32_013) 
            IF ( LCARB ) CALL CHEMCARBON_ADJ

            !*** MINERAL DUST AEROSOLS ***
            IF ( LDUST ) THEN 

               ! Do dust aerosol chemsitry 
               ! Adjoint now supported (dkh, 01/13/12, adj32_011) 
               CALL CHEMDUST_ADJ

               ! Compute dust OD's & surface areas
               !CALL RDUST_ONLINE( SOILDUST )
            ENDIF

      !---------------------------------
      ! Rn-Pb-Be
      !---------------------------------                 
      ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN
         CALL ERROR_STOP('Simulation not supported: 2 ', 
     &                   'chemistry_adj_mod.f')

         CALL CHEMRnPbBe 
         CALL DRYFLXRnPbBe
                  
      !---------------------------------
      ! CH3I
      !---------------------------------
      ELSE IF ( ITS_A_CH3I_SIM() ) THEN
         CALL ERROR_STOP('Simulation not supported: 3 ', 
     &                   'chemistry_adj_mod.f')

         CALL CHEMCH3I

      !---------------------------------            
      ! HCN
      !---------------------------------
      ELSE IF ( ITS_A_HCN_SIM() ) THEN
         CALL ERROR_STOP('Simulation not supported: 4 ', 
     &                   'chemistry_adj_mod.f')
         CALL CHEM_HCN_CH3CN( N_TRACERS, STT )

      !---------------------------------
      ! Tagged O3
      !---------------------------------
      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN 
 
         ! lzh 12/08/2009 add tagged ox adjoint
         CALL CHEM_TAGGED_OX_ADJ
 
      !---------------------------------
      ! Tagged CO
      !---------------------------------
      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
         !mak debug
         print*, 'its tag CO chemistry adj'

         CALL CHEM_TAGGED_CO_ADJ

      !---------------------------------
      ! C2H6
      !---------------------------------
      ELSE IF ( ITS_A_C2H6_SIM() ) THEN
         CALL ERROR_STOP('Simulation not supported: 6 ', 
     &                  'chemistry_adj_mod.f')
         CALL CHEMC2H6

      !---------------------------------
      ! CH4 now supported (adj32_023) 
      !---------------------------------
      ELSE IF ( ITS_A_CH4_SIM() ) THEN

         CALL CHEMCH4_ADJ

      !---------------------------------
      ! Mercury
      !---------------------------------
      ELSE IF ( ITS_A_MERCURY_SIM() ) THEN
         CALL ERROR_STOP('Simulation not supported: 8 ', 
     &                  'chemistry_adj_mod.f')

        ! Do Hg chemistry
        CALL CHEMMERCURY
          
      !---------------------------------
      ! Offline H2/HD
      !---------------------------------
      ELSE IF ( ITS_A_H2HD_SIM() ) THEN
        CALL ERROR_STOP('Simulation not supported: 9 ', 
     &                  'chemistry_adj_mod.f')
        CALL CHEM_H2_HD
        CALL DRYFLXH2HD
 
!-----------------------------------------------------------------------------
! Prior to 7/19/04:
! Fully install Kr85 run later (bmy, 7/19/04)
!            !---------------------------------
!            ! Kr85   
!            !---------------------------------
!            CASE ( 12 )
!               CALL CHEMKr85
!-----------------------------------------------------------------------------
         ENDIF

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CHEMISTRY_ADJ' )
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE DO_CHEMISTRY_ADJ

!------------------------------------------------------------------------------
! Use GCKPP_ADJ_DRIVER for solving chemistry in both directions (dkh, 07/31/09) 
! mak, comment out for now, while testing tagged CO (mak, 6/20/09)
!      SUBROUTINE GCKPP_DRIVER_ADJ( ) 
!
!******************************************************************************
!     Driver routine to perform adjoint integration of the full KPP chemistry 
!     mechanism. Based on Daven Henze's GCKPP_DRIVER.    (Kumaresh, 01/24/2008)
!******************************************************************************
!
      ! Reference to f90 modules
c$$$      USE COMODE_MOD,       ONLY : JLOP, CSPEC, IXSAVE, CSPEC_FOR_KPP,
c$$$     &                             IYSAVE, IZSAVE, R_KPP, HSAVE_KPP,
c$$$     &                             CSPEC_ADJ, CSPEC_ADJ_FOR_KPP, 
c$$$     &                             EMIS_RATE
c$$$      USE TRACER_MOD,       ONLY : DDEP_ADJ, EMIS_ADJ, EMIS_I_ADJ
c$$$      USE TIME_MOD,         ONLY : GET_TS_CHEM, GET_LOCALTIME
c$$$      USE GCKPP_UTIL,       ONLY : Shuffle_kpp2user,INIT_KPP
c$$$      USE GCKPP_Initialize, ONLY : Initialize
c$$$      USE GCKPP_Rates,      ONLY : UPDATE_RCONST
c$$$      USE GCKPP_Monitor,    ONLY : SPC_NAMES
c$$$      USE ERROR_MOD,        ONLY : ERROR_STOP
c$$$      USE LOGICAL_MOD,      ONLY : LEMIS, LDRYD
c$$$      USE GCKPP_Global,     ONLY : SMAL2, VAR, VAR_ADJ, V_CSPEC,
c$$$     &                             V_CSPEC_ADJ, VAR_R_ADJ, RCONST
c$$$      USE gckpp_Function
c$$$      USE gckpp_Model
c$$$
c$$$      USE GCKPP_adj_Initialize,    ONLY : Initialize_adj
c$$$      USE GCKPP_adj_Integrator_em, ONLY : INTEGRATE_em_adj, NIERR,
c$$$     &                                    Nhnew, Nhexit
c$$$      USE GCKPP_adj_Integrator,    ONLY : INTEGRATE_adj
c$$$      
c$$$      ! Local variables
c$$$      TYPE (XPLEX)        :: T, TIN, TOUT
c$$$      INTEGER       :: ICNTRL(20)
c$$$      complex\(kind\=dp\) :: RCNTRL(20)
c$$$      INTEGER       :: ISTATUS(20)
c$$$      INTEGER       :: I, J, L, N, JJLOOP
c$$$      INTEGER       :: IH, JH, LH 
c$$$      INTEGER       :: TID, OMP_GET_THREAD_NUM 
c$$$      complex\(kind\=dp\) :: RSTATE(20)
c$$$      LOGICAL, SAVE :: FIRST = .TRUE. 
c$$$
c$$$      INTEGER, PARAMETER :: NADJ = NVAR
c$$$      complex\(kind\=dp\), DIMENSION(NVAR,NADJ) :: ATOL_adj, RTOL_adj
c$$$
c$$$!~~~> Tests
c$$$      complex\(kind\=dp\) :: VAR0(NVAR), VAR1(NVAR), VAR2(NVAR),fd,ad
c$$$
c$$$!~~~  > Output variables     
c$$$      complex\(kind\=dp\) :: Vdot(NVAR)
c$$$
c$$$      !=================================================================
c$$$
c$$$      STEPMIN = 0.0d0
c$$$      STEPMAX = 0.0d0
c$$$      
c$$$      DO i=1,NVAR
c$$$         RTOL(i) = 1.0d-3
c$$$         ATOL(i) = 1.0d-2
c$$$      END DO   
c$$$
c$$$      DO j=1,NADJ
c$$$        DO i=1,NVAR
c$$$          RTOL_adj(i,j) = 0!1.0d-4
c$$$          ATOL_adj(i,j) = 0!1.0d-10
c$$$        END DO
c$$$      END DO
c$$$    
c$$$!     -------------
c$$$      CALL INIT_KPP
c$$$!     -------------
c$$$     
c$$$      ! Set parameters to default. See comments in RosenbrockADJ for
c$$$      ! a list of the defaults.
c$$$      ICNTRL(:) = 0
c$$$      RCNTRL(:) = 0.d0
c$$$
c$$$      ! Change some parameters from the default to new values
c$$$      ICNTRL(1) = 1    ! Autonomous
c$$$      ICNTRL(2) = 0    ! Nonautonomous
c$$$
c$$$      ! Select Integrator
c$$$      !    ICNTRL(3)  -> selection of a particular Rosenbrock method
c$$$      !        = 0 :  default method is Rodas3
c$$$      !        = 1 :  method is  Ros2
c$$$      !        = 2 :  method is  Ros3 
c$$$      !        = 3 :  method is  Ros4 
c$$$      !        = 4 :  method is  Rodas3
c$$$      !        = 5:   method is  Rodas4
c$$$      ICNTRL(3) = 4    
c$$$
c$$$      ICNTRL(7) = 2           ! 1 = No adjoint, 2 = discrete adjoint
c$$$      
c$$$      IF(FIRST)THEN
c$$$
c$$$         
c$$$         RSTATE(2) = 0d0
c$$$                                ! reset FIRST flag 
c$$$         FIRST = .FALSE. 
c$$$
c$$$      ENDIF
c$$$
c$$$      ! GET TS_CHEM and convert it to seconds. 
c$$$      DT = GET_TS_CHEM() * 60d0
c$$$
c$$$      ! Set time parameters. 
c$$$      T = 0d0
c$$$      TIN = T
c$$$      TOUT = T + DT   
c$$$
c$$$      !=================================================================
c$$$      ! Solve Chemistry
c$$$      !=================================================================
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( JJLOOP, I, J, L, N, RSTATE, ISTATUS )
c$$$!$OMP+FIRSTPRIVATE( RCNTRL, ICNTRL )
c$$$!$OMP+COPYIN( TIME )
c$$$!$OMP+SCHEDULE( DYNAMIC )
c$$$      DO JJLOOP = 1,NTT
c$$$         
c$$$         JLOOP = JJLOOP
c$$$         ! Get 3D coords from SMVGEAR's 1D coords
c$$$         I = IXSAVE(JJLOOP)
c$$$         J = IYSAVE(JJLOOP)
c$$$         L = IZSAVE(JJLOOP)
c$$$      
c$$$         DO N =1, NVAR
c$$$            V_CSPEC(N) = CSPEC_FOR_KPP(JLOOP,N)
c$$$            !V_CSPEC_ADJ(N) = CSPEC_ADJ_FOR_KPP(JLOOP,N)
c$$$            V_CSPEC_ADJ(N) = CSPEC_ADJ_(JLOOP,N)
c$$$         END DO
c$$$
c$$$         ! Pass tracer concentrations from CSPEC_FOR_KPP to KPP working vectors VAR, FIX.
c$$$         ! This also initializes the constant rate constants.
c$$$         CALL Initialize()
c$$$        
c$$$         CALL Initialize_adj()
c$$$
c$$$         RCNTRL(3) = HSAVE_KPP(I,J,L)
c$$$
c$$$         ! Recalculate rate constants
c$$$         CALL Update_RCONST()    !*******************!
c$$$         
c$$$         !------switch---------
c$$$         IF(LEMIS.or.LDRYD)THEN
c$$$            CALL INTEGRATE_EM_ADJ(1, VAR, VAR_ADJ, VAR_R_ADJ, TIN, TOUT,
c$$$     &           ATOL_adj, RTOL_adj, ICNTRL, RCNTRL, ISTATUS, RSTATE)
c$$$         ELSE
c$$$            CALL INTEGRATE_ADJ(1, VAR, VAR_ADJ, TIN, TOUT,ATOL_adj, 
c$$$     &           RTOL_adj, ICNTRL, RCNTRL, ISTATUS, RSTATE)
c$$$         ENDIF
c$$$         !--------------------
c$$$
c$$$         IF ( ISTATUS(20) < 0 ) THEN !**************!
c$$$            rcntrl(3)  = 0d0
c$$$            CALL Initialize( )  ! v2.1 
c$$$            CALL Initialize_adj( )            
c$$$            CALL Update_RCONST()
c$$$            !------switch---------
c$$$            IF(LEMIS.or.LDRYD)THEN
c$$$               CALL INTEGRATE_EM_ADJ(1, VAR, VAR_ADJ, VAR_R_ADJ, TIN, 
c$$$     &              TOUT, ATOL_adj, RTOL_adj, ICNTRL, RCNTRL, ISTATUS, 
c$$$     &              RSTATE)
c$$$            ELSE
c$$$               CALL INTEGRATE_ADJ(1, VAR, VAR_ADJ, TIN, TOUT,ATOL_adj, 
c$$$     &              RTOL_adj, ICNTRL, RCNTRL, ISTATUS, RSTATE)
c$$$            ENDIF
c$$$            !---------------------
c$$$            IF ( ISTATUS(20) < 0 ) THEN 
c$$$              print*, 'failed twice !!! '
c$$$              CALL ERROR_STOP('IERR < 0 ', 'INTEGRATE_ADJ')
c$$$            ENDIF
c$$$         ENDIF 
c$$$
c$$$         ! Set negative values to SMAL2
c$$$         DO N = 1, NVAR
c$$$            VAR(N) = MAX(VAR(N),SMAL2)
c$$$         ENDDO
c$$$
c$$$         CALL Shuffle_kpp2user(VAR_ADJ,V_CSPEC_ADJ)
c$$$         CALL Shuffle_kpp2user(VAR,V_CSPEC)
c$$$
c$$$         DO N =1, NVAR
c$$$            CSPEC(JLOOP,N)     = V_CSPEC(N)
c$$$            CSPEC_ADJ(JLOOP,N) = V_CSPEC_ADJ(N)
c$$$         END DO
c$$$
c$$$        !------switch---------
c$$$         IF(LEMIS.or.LDRYD)THEN
c$$$        !==================================
c$$$        !   Scaled Emission Adjoints for NO, NO2, CO, ALK4
c$$$        !   ISOP,  ACET, PRPE, C3H8, C2H6, MEK, ALD2, CH2O
c$$$        !----------------------------------
c$$$         DO N =1, 12           !232-243 emission variables
c$$$            EMIS_ADJ(I,J,L,N) = EMIS_ADJ(I,J,L,N) 
c$$$     &                        + VAR_R_ADJ(N)*RCONST(N+231)
c$$$         END DO
c$$$        !----------------------------------
c$$$
c$$$        !==================================
c$$$        !   Drydeposition Rate Adjoints
c$$$        !----------------------------------
c$$$         DO N =13, NCOEFF      !244-253 drydep variables
c$$$            DDEP_ADJ(I,J,L,N) = DDEP_ADJ(I,J,L,N) 
c$$$     &                        + VAR_R_ADJ(N)*RCONST(N+231)
c$$$         END DO
c$$$        !----------------------------------
c$$$         
c$$$        !==================================
c$$$        !  Scaled Individual Source Emissions
c$$$        !----------------------------------
c$$$         DO N =1, 3           !1-3 NOx (1-Anthro, 2-Soil, 3-Aircraft/Lightning)
c$$$            EMIS_I_ADJ(I,J,L,N) = EMIS_I_ADJ(I,J,L,N) 
c$$$     &                        + VAR_R_ADJ(1)*EMIS_RATE(JLOOP,N)
c$$$         END DO
c$$$         DO N=4, 13           !4-13 Anthropogenic (except NOx)
c$$$            EMIS_I_ADJ(I,J,L,N) = EMIS_I_ADJ(I,J,L,N) 
c$$$     &                        + VAR_R_ADJ(N-2)*EMIS_RATE(JLOOP,N)
c$$$         END DO
c$$$         DO N=14, 24          !14-24 Biomass Burning
c$$$            EMIS_I_ADJ(I,J,L,N) = EMIS_I_ADJ(I,J,L,N) 
c$$$     &                        + VAR_R_ADJ(N-13)*EMIS_RATE(JLOOP,N)
c$$$         END DO
c$$$         DO N=25, 35          !25-35 Biofuel Burning
c$$$            EMIS_I_ADJ(I,J,L,N) = EMIS_I_ADJ(I,J,L,N) 
c$$$     &                        + VAR_R_ADJ(N-24)*EMIS_RATE(JLOOP,N)
c$$$         END DO
c$$$        !----------------------------------
c$$$         ENDIF
c$$$
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE GCKPP_DRIVER_ADJ
!------------------------------------------------------------------------------
      ! End of module
      END MODULE CHEMISTRY_ADJ_MOD
