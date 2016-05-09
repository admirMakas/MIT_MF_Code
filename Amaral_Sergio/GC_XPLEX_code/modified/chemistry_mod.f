! $Id: chemistry_mod.f,v 1.15 2012/03/04 18:40:03 daven Exp $
      MODULE CHEMISTRY_MOD
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
      
      SUBROUTINE DO_CHEMISTRY
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
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,     ONLY : OCEAN_SINK_ACET
      USE AEROSOL_MOD,     ONLY : AEROSOL_CONC, AEROSOL_RURALBOX
      USE AEROSOL_MOD,     ONLY : RDAER,        SOILDUST
      USE C2H6_MOD,        ONLY : CHEMC2H6
      USE CARBON_MOD,      ONLY : CHEMCARBON
      USE CH3I_MOD,        ONLY : CHEMCH3I
      USE DAO_MOD,         ONLY : CLDF,    DELP
      USE DAO_MOD,         ONLY : OPTDEP,  OPTD,   T
      USE DRYDEP_MOD,      ONLY : DRYFLX, DRYFLXRnPbBe, DRYFLXH2HD
      USE DUST_MOD,        ONLY : CHEMDUST, RDUST_ONLINE
      USE ERROR_MOD,       ONLY : DEBUG_MSG,ERROR_STOP,GEOS_CHEM_STOP
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
      USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
      USE TAGGED_OX_MOD,   ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
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
      ! ADJ_GROUP
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      USE LOGICAL_ADJ_MOD, ONLY : LAERO_THERM
      USE TRACERID_MOD
      
#     include "CMN_SIZE"        ! Size parameters
#     include "CMN_DIAG"        ! NDxx flags
#     include "comode.h"        ! NPHOT

      ! Local variables
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: N_TROP,i,j,k,n,l

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

      ! Compute optical depths (except for CH4 simulation)
      IF ( .not. ITS_A_CH4_SIM() ) THEN
         CALL OPTDEPTH( LLPAR, CLDF, OPTDEP, OPTD )
      ENDIF
      !=================================================================
      ! If LCHEM=T then call the chemistry subroutines
      !=================================================================
      IF ( LCHEM ) THEN 

         !---------------------------------
         ! NOx-Ox-HC (w/ or w/o aerosols) 
         !---------------------------------
         IF ( ITS_A_FULLCHEM_SIM() ) THEN 

            ! Call SMVGEAR routines
            CALL CHEMDR
!             do n=1,size(STT,4)
!      do l=1,size(STT,3)
!      do j=1,size(STT,2)
!      do i=1,size(STT,1)
!      if (abs(STT(i,j,l,n)%i)>1e-210) then
!          print*,'STT CHEMDR',STT(i,j,l,n),n,i,j,l
!          CALL GEOS_CHEM_STOP
!      endif
!      enddo
!      enddo
!      enddo
!      enddo
            ! Do seasalt aerosol chemistry
            IF ( LSSALT ) THEN
               IF ( LADJ ) WRITE(6,*) 'WARNING: NEED SSALT ADJOINT'
               !CALL CHEMSEASALT
            END IF

            ! Also do sulfate chemistry
            IF ( LSULF ) THEN

               ! Do sulfate chemistry
               CALL CHEMSULFATE
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
                  ! adj_group: implement flag for aerosol thermo
                  IF ( LAERO_THERM ) THEN 
                     CALL DO_RPMARES
                  ENDIF 

               !ENDIF
               !------------------------------------------------------------
               
            ENDIF

            ! Do carbonaceous aerosol chemistry
            IF ( LCARB ) CALL CHEMCARBON  

            ! Do dust aerosol chemistry
            IF ( LDUST ) THEN
               CALL CHEMDUST
            END IF

            ! ND44 drydep fluxes
            CALL DRYFLX     
            ! ND43 chemical production
            CALL DIAGOH
            ! Remove acetone ocean sink
            IF ( IDTACET /= 0 ) THEN
               CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
            ENDIF
         !---------------------------------
         ! Offline aerosol simulation
         ! Now enabled for OC, BC and dust (yhmao dkh, 01/13/12, adj32_013) 
         !---------------------------------
         ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            ! Define loop index and other SMVGEAR arrays
            ! N_TROP, the # of trop boxes, is returned
            !CALL AEROSOL_RURALBOX( N_TROP )

            !  turn off for BC and DUST only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            !! Initialize FAST-J quantities for computing AOD's
            !IF ( FIRST ) THEN
            ! 
            !   !CALL READCHEM
            !   !CALL INPHOT( LLTROP, NPHOT )
            !
            !   ! Reset NCS with NCSURBAN
            !   NCS     = NCSURBAN
            !
            !   ! Reset NTLOOP and NTTLOOP after call to READER
            !   ! with the actual # of boxes w/in the ann mean trop
            !   NTLOOP  = N_TROP
            !   NTTLOOP = N_TROP
            !
            !   ! Reset first-time flag
            !   FIRST = .FALSE.
            !ENDIF

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

               ! RPMARES does not take Na+, Cl- into account
               ! (skip for crystalline & aqueous offline run)
               !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
               !IF ( .not. LCRYST ) CALL DO_RPMARES

            !ENDIF
            !-------------------------------------------------------------

            !*** SEASALT AEROSOLS ***
            !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            !IF ( LSSALT ) CALL CHEMSEASALT

            !*** SULFATE AEROSOLS ***
            !  turn off for BC only aerosol sim (yhmao, dkh, 01/13/12, adj32_013) 
            !IF ( LSULF .or. LCRYST ) THEN
            !
            !   ! Do sulfate chemistry
            !   CALL CHEMSULFATE
            !
            !ENDIF
               
            !*** CARBON AND 2NDARY ORGANIC AEROSOLS ***
            IF ( LCARB ) CALL CHEMCARBON
            !*** MINERAL DUST AEROSOLS ***
            IF ( LDUST ) THEN 

               ! Do dust aerosol chemsitry
               CALL CHEMDUST
               ! Compute dust OD's & surface areas
               ! Turn off for dust only offline aerosol sim (dkh, 03/04/12, adj32_013) 
               !CALL RDUST_ONLINE( SOILDUST )

            ENDIF

         !---------------------------------
         ! Rn-Pb-Be
         !---------------------------------                 
         ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN
            CALL CHEMRnPbBe 
            CALL DRYFLXRnPbBe
                  
         !---------------------------------
         ! CH3I
         !---------------------------------
         ELSE IF ( ITS_A_CH3I_SIM() ) THEN
            CALL CHEMCH3I
         !---------------------------------            
         ! HCN
         !---------------------------------
         ELSE IF ( ITS_A_HCN_SIM() ) THEN
            CALL CHEM_HCN_CH3CN( N_TRACERS, STT )
         !---------------------------------
         ! Tagged O3
         !---------------------------------
         ELSE IF ( ITS_A_TAGOX_SIM() ) THEN 
            CALL CHEM_TAGGED_OX
         !---------------------------------
         ! Tagged CO
         !---------------------------------
         ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
            CALL CHEM_TAGGED_CO
         !---------------------------------
         ! C2H6
         !---------------------------------
         ELSE IF ( ITS_A_C2H6_SIM() ) THEN
            CALL CHEMC2H6
         !---------------------------------
         ! CH4
         !---------------------------------
         ELSE IF ( ITS_A_CH4_SIM() ) THEN

            CALL CHEMCH4
         !---------------------------------
         ! Mercury
         !---------------------------------
         ELSE IF ( ITS_A_MERCURY_SIM() ) THEN

            ! Do Hg chemistry
            CALL CHEMMERCURY
         !---------------------------------
         ! Offline H2/HD
         !---------------------------------
         ELSE IF ( ITS_A_H2HD_SIM() ) THEN
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
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE DO_CHEMISTRY

!------------------------------------------------------------------------------

      !SUBROUTINE GCKPP_DRIVER( ) 
      SUBROUTINE GCKPP_ADJ_DRIVER( DIRECTION  ) 
!
!******************************************************************************
!  Subroutine GCKPP_DRIVER is the driver routine which performs the 
!  adjoint integration of the full chemistry mechanism.  (dkh, 07/18/05)
! 
!
!  NOTES:
!  
!   (1 ) Add argument DIRECTION so that we can use the same driver for fwd 
!         calculation and for recomputation during the bwd run. (dkh, 08/26/05)  
!   (2 ) Parallelize the loop over grid cells. (dkh, 10/07/05)  
!   (3 ) Call ADJ_CALCRATE to calculate adjoint of emissions and dep of 
!         species for which these processes are included in full chemistry
!         mechanism (dkh, 06/04/06). 
!   (4 ) Add FIRST and now call SET_SMV2KPP. (dkh, 06/06/06)  
!   (5 ) Update for KPP v2.2.  Many of the controls are indexed differently
!         since v2.1.  (dkh, 07/23/06) 
!   (6 ) Updated to GCv7 (ks, 01/08)
!   (7 ) Updated for GCv8 (dkh, ks, mak, cs, 06/09)
!   (8 ) Updated to run in parallel, fix format (dkh, 07/28/09) 
!   (9 ) Reuse internal time step (H) from the previous integration (dkh, 07/28/09) 
!   (10) Now only redo chemistry calculation in forward mode (dkh, 07/08/11, adj32_004) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,       ONLY : IFD, JFD, LFD 
      ! Now use HSAVE in comode (dkh, 04/25/10) 
      !USE CHECKPT_MOD,          ONLY : HSAVE
      USE COMODE_MOD,           ONLY : JLOP, CSPEC, IXSAVE
      USE COMODE_MOD,           ONLY : CSPEC_FOR_KPP,
     &                                 IYSAVE, IZSAVE, R_KPP
      USE COMODE_MOD,           ONLY : CSPEC_ADJ
      USE COMODE_MOD,           ONLY : HSAVE
      USE ERROR_MOD,            ONLY : ERROR_STOP,GEOS_CHEM_STOP
      USE GCKPP_ADJ_UTIL,       ONLY : Shuffle_kpp2user
      USE GCKPP_ADJ_UTIL,       ONLY : INIT_KPP
      USE GCKPP_ADJ_Initialize, ONLY : Initialize
      USE GCKPP_ADJ_Initialize, ONLY : Initialize_adj
      USE GCKPP_ADJ_Global,     ONLY : SMALL2,  VAR, VAR_ADJ, V_CSPEC,
     &                                 V_CSPEC_ADJ, VAR_R_ADJ, RCONST    
      USE GCKPP_ADJ_Global,     ONLY : ATOL
      USE GCKPP_ADJ_Global,     ONLY : RTOL
      USE GCKPP_ADJ_Global,     ONLY : NTT
      !USE GCKPP_ADJ_Global,     ONLY : TIME
      USE GCKPP_ADJ_Global,     ONLY : JLOOP
      USE GCKPP_ADJ_Global,     ONLY : DT
      USE GCKPP_ADJ_Global,     ONLY : STEPMIN
      USE GCKPP_ADJ_Global,     ONLY : STEPMAX
      USE GCKPP_ADJ_Global,     ONLY : RTOLS
      USE GCKPP_ADJ_Global,     ONLY : IND
      USE GCKPP_ADJ_Global,     ONLY : NCOEFF
      USE GCKPP_ADJ_Global,     ONLY : JCOEFF
      USE GCKPP_ADJ_Integrator, ONLY : INTEGRATE_ADJ, NIERR
      USE GCKPP_ADJ_Integrator, ONLY : NHnew, Nhexit
      USE GCKPP_ADJ_Parameters, ONLY : NVAR
      USE GCKPP_ADJ_Parameters, ONLY : ind_CO2, ind_DRYDEP
      USE GCKPP_ADJ_Parameters, ONLY : ind_LISOPOH
      USE GCKPP_ADJ_Rates,      ONLY : UPDATE_RCONST
      USE GCKPP_ADJ_Monitor,    ONLY : SPC_NAMES
      USE LOGICAL_ADJ_MOD,      ONLY : LPRINTFD
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ_EMS
      USE TIME_MOD,             ONLY : GET_TS_CHEM, GET_LOCALTIME
      USE TIME_MOD,             ONLY : GET_NYMD,    GET_NYMDb
      USE TIME_MOD,             ONLY : GET_NHMS,    GET_NHMSb
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ_ONLY !jkoo
      USE TRACER_MOD,           ONLY : STT
#     include "CMN_SIZE"             ! Size params
#     include "comode.h"             ! NMTRATE

      ! Argument 
      INTEGER                       :: DIRECTION 

      ! Local variables
      TYPE (XPLEX)                        :: T, TIN, TOUT
      INTEGER                       :: ICNTRL(20)
      TYPE (XPLEX)                        :: RCNTRL(20)
      INTEGER                       :: ISTATUS(20)
      INTEGER                       :: I, J, L, N, JJLOOP
      INTEGER                       :: IH, JH, LH 
      INTEGER                       :: TID, OMP_GET_THREAD_NUM 
      INTEGER                       :: ICOEFF, NNREAC
      TYPE (XPLEX)                        :: RSTATE(20)
      TYPE (XPLEX)                        :: RRATE_ADJ(NMTRATE)
      LOGICAL, SAVE                 :: FIRST = .TRUE. 

      ! SAFETY KLUDGE: hardwire this for now.  Change when parser is fixed (dkh, 01/27/10)
      INTEGER, PARAMETER            :: NVAR_GC = 87

      !=================================================================
      ! GCKPP_ADJ_DRIVER begins here!
      !=================================================================
      WRITE(6,*) '      - GCKPP_ADJ_DRIVER, DIRECTION =  ', DIRECTION

      ! Init VAR_ADJ
      VAR_ADJ(:) = 0.d0
      VAR(:)=0d0
      STEPMIN = 0.0d0
      STEPMAX = 0.0d0
      
      RTOLS   = 1d-1 
      DO I=1,NVAR
         RTOL(I) = RTOLS
         ATOL(I) = 1.0d-2
      END DO

      ! Set parameters to default. See comments in RosenbrockADJ for
      ! a list of the defaults.
      ICNTRL(:) = 0
      RCNTRL(:) = 0.d0

      ! Change some parameters from the default to new values
      ICNTRL(1) = 1    ! Autonomous
      ICNTRL(2) = 0    ! Nonautonomous

      ! Select Integrator
      !    ICNTRL(3)  -> selection of a particular Rosenbrock method
      !        = 0 :  default method is Rodas3
      !        = 1 :  method is  Ros2
      !        = 2 :  method is  Ros3 
      !        = 3 :  method is  Ros4 
      !        = 4 :  method is  Rodas3
      !        = 5:   method is  Rodas4
      ICNTRL(3) = 4   

      ! No adjoint needed for fwd run
      IF ( DIRECTION > 0 ) THEN

         ICNTRL(7) = 1    ! No adjoint

      ! Pick an adjoint method for bwd run  
      ELSEIF (  DIRECTION < 0 ) THEN

         !ICNTRL(7) = 1    ! No adjoint
         ICNTRL(7) = 2    ! Discrete adjoint
         !ICNTRL(7) = 3    ! Continuous adjoint

      ENDIF
 
      IF ( FIRST ) THEN
         
  
         CALL INIT_KPP 

         RSTATE(Nhexit) = 0d0
                               
         FIRST          = .FALSE. 

         ! SAFETY KLUDGE: to overcome lacking in the KPP <--> SMVGEAR 
         ! mapping, we've implemented some hardwired code below. 
         ! Here we are just checking to see if species indices are the
         ! same as the hardwired values.  Hopefully remove this
         ! once the parser is fixed. (dkh, 01/27/10) 
         IF ( ind_CO2     .ne. 13 .or. ind_DRYDEP .ne. 14 .OR. 
     &        ind_LISOPOH .ne. 15                              ) THEN
            CALL ERROR_STOP('Need to adjust hardwired mapping',
     &                      'GCKPP_ADJ_DRIVER')
         ENDIF


      ENDIF

      ! GET TS_CHEM and convert it to seconds. 
      DT = GET_TS_CHEM() * 60d0

      ! Set time parameters. 
      T = 0d0
      TIN = T
      TOUT = T + DT
      
      !=================================================================
      ! Solve Chemistry
      !=================================================================
      
      ! dkh debug
      IF ( LPRINTFD ) THEN 
         print*, ' NTT = ', NTT
         print*, ' IND = ', IND(:)
      ENDIF 

      ! OLD CODE: we don't use common blocks anymore.  
      !! Apparently this needs to be set in order to use THREADPRIVATE
      !! common blocks. 
      !CALL OMP_SET_DYNAMIC(.FALSE.)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JJLOOP, I, J, L, N, RSTATE, ISTATUS )
!$OMP+PRIVATE( ICOEFF, NNREAC, RRATE_ADJ           )
!$OMP+FIRSTPRIVATE( RCNTRL, ICNTRL )
!!$OMP+COPYIN( TIME )
!$OMP+SCHEDULE( DYNAMIC )
      DO JJLOOP = 1, NTT
!         print*,'JJLOOP/NTT',JJLOOP,NTT
         JLOOP = JJLOOP

         ! Get 3D coords from SMVGEAR's 1D coords
         I = IXSAVE(JJLOOP)
         J = IYSAVE(JJLOOP)
         L = IZSAVE(JJLOOP)

         DO N =1, NVAR
            V_CSPEC(N) = CSPEC_FOR_KPP(JLOOP,N)
         END DO 


         ! Pass tracer concentrations from CSPEC_FOR_KPP to KPP working vectors VAR, FIX.
         ! This also initializes the constant rate constants.
         CALL Initialize()
         ! Init sensitivity vector 
         IF ( DIRECTION < 0 ) THEN

             DO N =1, NVAR
               V_CSPEC_ADJ(N) = CSPEC_ADJ(JLOOP,N)
            END DO 

            ! Update from ks
            !CALL INIT_Y_ADJ( Y_ADJ(1:NVAR,NJ), ADCSPEC(JJLOOP,1:NVAR) )
            CALL Initialize_adj()

         ENDIF

         ! SAFETY KLUDGE:
         ! These are inactive species. Since they don't get inlcuded 
         ! in the mapping, have to set them here manually (dkh, 08/06/09) 
         VAR(13:14) = 0d0
         VAR(15)    = 1d-99

         ! Recalculate rate constants
         CALL Update_RCONST()    
         ! Error check
         IF ( MAXVAL(RCONST) .eq. 0.d0 ) THEN
             print*, jloop
             print*, I, J, L
             CALL ERROR_STOP('Rate contants are zero','chemistry_mod.f')
         ENDIF


         ! Start with the last internal time step from the previous
         ! integartion
         RCNTRL(Nhnew) = HSAVE(I,J,L)

         IF ( LPRINTFD .and. I==IFD .and. J==JFD .and. L==LFD ) THEN 
            print*, ' JLOOP in CHEM = ', JLOOP
            print*, ' V_CSPEC= ', V_CSPEC
            print*, ' VAR    = ', VAR
            print*, ' RCONST = ', RCONST
            print*, ' R(EMIS) = ', RCONST(229:242)
            print*, ' R(DRYD) = ', RCONST(243:252)
            print*, ' HSTART= ', RCNTRL(Nhnew)
            IF ( DIRECTION < 0 ) print*, ' VAR_ADJ = ', VAR_ADJ 
         
            ! Inscpect reaction rates and species concentrations 
            ! for consistancy in fwd and adj chemistry (dkh, 08/06/09) 
            IF ( .not. LADJ_ONLY ) THEN
               CALL CINSPECT( RCONST, VAR, DIRECTION ) 
            ENDIF
         ENDIF 

         
!         do i=1,size(VAR,1)
!         if (isnan(VAR(i)).or.isnan(TIN).or.isnan(TOUT).or.
!     &          exponent(VAR(i)%i)>0d0) then
!            print*,'VAR is NaN bef int_adj',VAR(i),i
!            CALL GEOS_CHEM_STOP
!         endif
!         enddo
                  
         CALL INTEGRATE_ADJ( 1, VAR, VAR_ADJ, VAR_R_ADJ, TIN, TOUT, 
     &        ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATE)
!         do i=1,size(VAR,1)
!         if (isnan(VAR(i)).or.exponent(VAR(i)%i)>0d0) then
!            print*,'VAR is aft int_adj',VAR(i),i
!            CALL GEOS_CHEM_STOP
!         endif
!         enddo
     
         ! Sometimes the integration will fail the first time.  
         ! That's OK.  Just reset and try again w/ smaller step size. 
         IF ( ISTATUS(NIERR) < 0 ) THEN 
            print*, 'IERR = ', RSTATE(20)
            print*, 'RSTAT = ', RSTATE
            print*, 'ISTAT = ', ISTATUS
            print*, 'RCNTRL = ', RCNTRL
            print*, 'ICNTRL = ', ICNTRL
            print*, 'JLOOP, I, J, L ', JLOOP, I, J, L
            print*, ' REDO THIS CELL WITH H_START = 0 '
            rcntrl(3)  = 0d0
            CALL Initialize( )  ! v2.1 
            VAR(15) = 1d-99
            IF ( DIRECTION < 0 ) CALL Initialize_adj()
            CALL Update_RCONST()

            ! Now only repeat this for forward run (dkh, xxu, 07/08/11, adj32_004) 
            IF ( DIRECTION > 0 ) THEN 

               CALL INTEGRATE_ADJ(1, VAR, VAR_ADJ, VAR_R_ADJ, TIN, TOUT, 
     &                      ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATE)

 
               ! There seems to be an issue with random failing in cells
               ! I = J = x, L = 1, where x = 1:8, only on the first 
               ! time through.  An OpenMP issues? Move along if fails
               ! twice in these cells during first call to chemistry. 
               IF ( ISTATUS(NIERR) < 0 ) THEN 
                  print*, 'failed twice !!! '
                  IF ( I < 9 .and. I == J .and. L == 1 .and. 
     &                GET_NYMD() == GET_NYMDb() .and. 
     &                GET_NHMS() == GET_NHMSb() ) THEN
                     print*, 'just move on !!! '
                  ELSE 
                     CALL ERROR_STOP('IERR < 0 ', 'INTEGRATE_ADJ')
                  ENDIF 

               ENDIF

            ! Now reset adjoint values and move on. (dkh, 07/08/11, adj32_004) 
            ELSE 

               WRITE(6,*) ' JUST RESET ADJOINT VALUES AND MOVE ON '

               DO N =1, NVAR
                  V_CSPEC_ADJ(N) = CSPEC_ADJ(JLOOP,N)
               END DO
               CALL Initialize_adj()
               IF ( LADJ_EMS ) THEN
                  DO ICOEFF = 1, NCOEFF
                    NNREAC = JCOEFF(ICOEFF)
                    VAR_R_ADJ(ICOEFF) = RRATE_ADJ(IND(NNREAC))
                  ENDDO
               ENDIF 
            ENDIF 

         ENDIF 
         ! Set negative values to SMALL2
         DO N = 1, NVAR
            VAR(N) = MAX(VAR(N),SMALL2)
         ENDDO
         !print*,'VAR aft SMALL2',VAR
         IF ( LPRINTFD .and. I==IFD .and. J==JFD .and. L==LFD ) THEN 
            print*, ' VAR after  = ', VAR
            print*, ' VAR_ADJ after = ', VAR_ADJ 
            print*, ' HSAVE         = ', RSTATE(3)
         ENDIF 
 
         ! Save last internal time step, Hexit, during fwd integration. 
         !IF ( DIRECTION > 0 ) HSAVE(I,J,L) = RSTAT(Nhexit)
         ! It seems to be faster if we use hlast rather than hexit. (dkh, 07/23/06) 
         IF ( DIRECTION > 0 ) HSAVE(I,J,L) = RSTATE(3)

         ! SAFETY KLUDGE:  the KPP <--> SMVGEAR mapping isn't 1:1.  There
         ! are species included in KPP but not in SMVGEAR active list,
         ! so only map back up to NVAR_GC.  (dkh, 01/27/10) 

         ! Map the KPP results back to CSPEC_FORKPP
         ! This is only necessary for DIRECTION > 0, but is included for 
         ! debugging in both directions for the moment. 
         CALL SHUFFLE_KPP2USER( VAR, V_CSPEC )  
         !DO N =1, NVAR
         DO N =1, NVAR_GC
            CSPEC_FOR_KPP(JLOOP,N)     = V_CSPEC(N)
         END DO

         ! Map results back to CSPEC_ADJ
         IF ( DIRECTION < 0 ) THEN 
            CALL SHUFFLE_KPP2USER( VAR_ADJ, V_CSPEC_ADJ )  
            !DO N =1, NVAR
            DO N =1, NVAR_GC
               CSPEC_ADJ(JLOOP,N) = V_CSPEC_ADJ(N)
            END DO

            ! Map the KPP rate adjoints back to RRATE_ADJ
            IF ( LADJ_EMS ) THEN 
               DO ICOEFF = 1, NCOEFF
                 NNREAC = JCOEFF(ICOEFF)
                 RRATE_ADJ(IND(NNREAC)) = VAR_R_ADJ(ICOEFF)
               ENDDO 

               ! Calculate the adjoint of emissions and drydep rates
               CALL CALCRATE_ADJ( RRATE_ADJ, I, J, L )
            ENDIF 

         ENDIF 
      ENDDO
!$OMP END PARALLEL DO
      ! Need to update this routine before implementing 
      !! Calculate ARR: save the reaction rate adjoints 
      !IF ( DIRECTION < 0 ) CALL SAVE_INST_ARR

      ! Return to calling program
      END SUBROUTINE GCKPP_ADJ_DRIVER

!------------------------------------------------------------------------------

      SUBROUTINE CINSPECT( RCONST, VAR, DIRECTION )
!
!******************************************************************************
!  Subroutine CINSPECT save reaction rates and species concentrations
!  during the forward run and checks to make sure these match values
!  during the adjoint run. (dkh, 08/06/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : LONGITUDE center of grid box [degrees]
!
!  NOTES:
!  (1 ) Now make VAR_FD and RCONST_FD allocatable and move them to 
!         adj_arrays_mod.f (dkh, 02/23/12, adj32_026) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,       ONLY : IFD, JFD, LFD 
      USE ADJ_ARRAYS_MOD,       ONLY : VAR_FD
      USE ADJ_ARRAYS_MOD,       ONLY : RCONST_FD
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE GCKPP_ADJ_PARAMETERS, ONLY : NVAR, NREACT

      ! Arguements 
      TYPE (XPLEX)                        :: RCONST(NREACT)
      TYPE (XPLEX)                        :: VAR(NVAR)
      INTEGER                       :: DIRECTION

      ! Local variables
      INTEGER                       :: R, N
      INTEGER, SAVE                 :: NCHEM = 1
      TYPE (XPLEX)                        :: X
      !TYPE (XPLEX),  SAVE                 :: VAR_FD(NVAR,1000)
      !TYPE (XPLEX),  SAVE                 :: RCONST_FD(NREACT,1000)
      LOGICAL, SAVE                 :: PASS = .TRUE. 
     
      ! Parameters
      !TYPE (XPLEX), PARAMETER             :: EPS = 0.000005d0
      ! for debugging:
      !TYPE (XPLEX), PARAMETER             :: EPS = 0.005000d0
      TYPE (XPLEX), PARAMETER             :: EPS = xplex(0.000500d0,0d0)

      !=================================================================
      ! CINSPECT begins here!
      !=================================================================

      ! No longer needed, as NCHEM_MAX dynamic (dkh, 02/23/12, adj32_026) 
      !IF ( NCHEM > 1000 ) THEN 
      !   CALL ERROR_STOP( 'Need to boost NCHEM max', 'chemistry_mod.f')
      !ENDIF 

      IF ( DIRECTION > 0 ) THEN

         ! Save from forward run 
         RCONST_FD(:,NCHEM) = RCONST(:)
         VAR_FD   (:,NCHEM) = VAR(:)
         NCHEM              = NCHEM + 1 

      ELSE
      
         NCHEM = NCHEM - 1

         ! Check species concentration against saved values from forward
         DO N = 1, NVAR
            
            IF ( ABS( VAR(N) - VAR_FD(N,NCHEM) ) > 
     &           ABS( VAR(N) + VAR_FD(N,NCHEM) ) * EPS ) THEN   
               ! Don't care if VAR_FD = 1.000398d-99 and VAR = 1.000000d-99
               IF ( ABS( VAR(N) + VAR_FD(N,NCHEM) ) > 1d-98 ) THEN 
                  WRITE(6,*) ' Error in VAR ', N
                  WRITE(6,*) '  - Forward  value  =', VAR_FD(N,NCHEM)
                  WRITE(6,*) '  - Backward value  =', VAR(N)
                  PASS   = .FALSE.
               ENDIF 
            ENDIF 

         ENDDO

         ! Check reaction rates against saved values from forward
         DO R = 1, NREACT

            IF ( ABS( RCONST(R) - RCONST_FD(R,NCHEM) ) > 
     &           ABS( RCONST(R) + RCONST_FD(R,NCHEM) ) * EPS ) THEN
               WRITE(6,*) ' Error in  RCONST ', R
               WRITE(6,*) ' - Forward  value  =', RCONST_FD(R,NCHEM)
               WRITE(6,*) ' - Backward value  =', RCONST(R)
               PASS   = .FALSE.
            ENDIF

         ENDDO

         IF ( .not. PASS ) THEN 
            CALL ERROR_STOP('Stop at CINSPECT', 'chemistry_mod.f')
         ENDIF 

      ENDIF

      ! Return to calling program
      END SUBROUTINE CINSPECT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CHEMISTRY_MOD
