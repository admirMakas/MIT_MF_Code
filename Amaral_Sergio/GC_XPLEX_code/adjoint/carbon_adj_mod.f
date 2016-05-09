! $Id: carbon_adj_mod.f,v 1.4 2011/02/23 00:08:47 daven Exp $
      MODULE CARBON_ADJ_MOD
!
!******************************************************************************
!  Module CARBON_ADJ_MOD contains arrays and routines for performing an offline 
!  carbonaceous aerosol adjoint simulation.  Original code taken from forward 
!  routines in CARBON_MOD and modified accordingly.  (dkh, 03/01/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BCCONV_ADJ     (TYPE (XPLEX)) : Adjoint of BCCONV
!  (2 ) OCCONV_ADJ     (TYPE (XPLEX)) : Adjoint of OCCONV
!
!  Module Routines:
!  ============================================================================
!  (1 ) ADJ_CHEMCARBON     : Driver program for adjoint carbon aerosol chemistry
!  (2 ) ADJ_CHEM_BCPO      : Chemistry routine for hydrophobic BC (aka EC)
!  (3 ) ADJ_CHEM_BCPI      : Chemistry routine for hydrophilic BC (aka EC)
!  (4 ) ADJ_CHEM_OCPO      : Chemistry routine for hydrophobic OC
!  (5 ) ADJ_CHEM_OCPI      : Chemistry routine for hydrophilic OC
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
!  (1 ) See original forward module for all notes. 
!  (2 ) Change BCCONV and OCCONV to ADJ_BCCONV and ADJ_OCCONV. (dkh, 03/22/07)
!  (3 ) Updated to GCv8 (dkh, 09/09/09) 
!
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
      PUBLIC :: CHEMCARBON_ADJ
      PUBLIC :: EMISSCARBON_ADJ
      PUBLIC :: CLEANUP_CARBON_ADJ

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

! Comment out module variables from forward routine that we don't use 
! for the adjoint.  (dkh, 09/09/09) 
!      ! Scalars
!      LOGICAL             :: USE_MONTHLY_BIOB = .TRUE.
!      INTEGER             :: DRYBCPI, DRYOCPI, DRYBCPO, DRYOCPO
!      INTEGER             :: DRYALPH, DRYLIMO, DRYALCO
!      INTEGER             :: DRYSOG1, DRYSOG2, DRYSOG3, DRYSOG4
!      INTEGER             :: DRYSOA1, DRYSOA2, DRYSOA3, DRYSOA4
!      INTEGER             :: I1_NA,   J1_NA
!      INTEGER             :: I2_NA,   J2_NA
!      INTEGER             :: DRYSOAG, DRYSOAM
!      
!      ! Parameters
!      INTEGER, PARAMETER  :: MHC      = 6
!      INTEGER, PARAMETER  :: NPROD    = 3  
!      TYPE (XPLEX),  PARAMETER  :: SMALLNUM = 1d-20
!
!      ! Arrays
!      TYPE (XPLEX), ALLOCATABLE :: ANTH_BLKC(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: ANTH_ORGC(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOB_BLKC(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOB_ORGC(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOF_BLKC(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOF_ORGC(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: EF_BLKC(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: EF_ORGC(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: TERP_ORGC(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BCCONV(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: OCCONV(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOG_ALPH(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOG_LIMO(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOG_ALCO(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOG_TERP(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: BIOG_SESQ(:,:)    
!      TYPE (XPLEX), ALLOCATABLE :: DIUR_ORVC(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: GEIA_ORVC(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: TCOSZ(:,:)
!      TYPE (XPLEX), ALLOCATABLE :: ORVC_SESQ(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: ORVC_TERP(:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: GPROD(:,:,:,:,:)
!      TYPE (XPLEX), ALLOCATABLE :: APROD(:,:,:,:,:)
!      ! Cloud fraction - for cloud droplet uptake of dicarbonyls 
!      ! (tmf, 12/07/07) 
!      TYPE (XPLEX), ALLOCATABLE :: VCLDF(:,:,:)

      TYPE (XPLEX), ALLOCATABLE :: BCCONV_ADJ(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: OCCONV_ADJ(:,:,:)

      ! Days per month (based on 1998)
      INTEGER             :: NDAYS(12) = (/ 31, 28, 31, 30, 31, 30, 
     &                                      31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMCARBON_ADJ
!
!******************************************************************************
!  Subroutine CHEMCARBON is the interface between the GEOS-CHEM main 
!  program and the adjoint carbon aerosol chemistry routines that calculates
!  dry deposition and chemical conversion between hydrophilic and 
!  hydrophobic.
!
!  NOTES:
!  (1 ) Based on CHEMCARBON from forward model. (rjp, bmy, 4/1/04, 9/14/06)
!      The only differences are:
!        i.  Use STT_ADJ instead of STT
!        ii. Call CHEM_xxxx_ADJ rather than CHEM_xxxx
!
!  NOTES:
!  (1 ) See forword module for all notes. 
!  (2 ) Updated to GCv8 (dkh, 09/09/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DRYDEP_MOD,     ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE ERROR_MOD,      ONLY : ERROR_STOP
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
      ! CHEMCARBON_ADJ begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_CARBON_ADJ

         ! Don't need to repeat the rest of this (dkh, 09/09/09) 
!         ! Find drydep species in DEPSAV
!         DO N = 1, NUMDEP
!            SELECT CASE ( TRIM( DEPNAME(N) ) )
!               CASE ( 'BCPI' )
!                  DRYBCPI = N
!               CASE ( 'OCPI' )
!                  DRYOCPI = N
!               CASE ( 'BCPO' )
!                  DRYBCPO = N
!               CASE ( 'OCPO' )
!                  DRYOCPO = N
!               CASE ( 'ALPH' )
!                  DRYALPH = N
!               CASE ( 'LIMO' )
!                  DRYLIMO = N
!               CASE ( 'ALCO' )
!                  DRYALCO = N
!               CASE ( 'SOG1' )
!                  DRYSOG1 = N
!               CASE ( 'SOG2' )
!                  DRYSOG2 = N
!               CASE ( 'SOG3' )
!                  DRYSOG3 = N
!               CASE ( 'SOG4' )
!                  DRYSOG4 = N
!               CASE ( 'SOA1' )
!                  DRYSOA1 = N
!               CASE ( 'SOA2' )
!                  DRYSOA2 = N
!               CASE ( 'SOA3' )
!                  DRYSOA3 = N
!               CASE ( 'SOA4' )
!                  DRYSOA4 = N
!               CASE ( 'SOAG' )
!                  DRYSOAG = N
!               CASE ( 'SOAM' )
!                  DRYSOAM = N
!               CASE DEFAULT
!                  ! Nothing
!            END SELECT        
!         ENDDO
!
!         ! Zero SOG4 and SOA4 (SOA from ISOP in gas & aerosol form)
!         ! for offline aerosol simulations.  Eventually we should have
!         ! archived isoprene oxidation fields available for offline
!         ! simulations but for now we just set them to zero. 
!         ! (dkh, bmy, 6/1/06)
!         IF ( ITS_AN_AEROSOL_SIM() ) THEN
!
!            ! temp fix for aerosol w/ 20 tracers simulation (phs)
!            IF ( IDTSOG4 .NE. 0 ) THEN   
!               STT(:,:,:,IDTSOG4) = 0d0
!               STT(:,:,:,IDTSOA4) = 0d0
!            ENDIF
!         ENDIF
         
         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      !=================================================================
      ! Do chemistry for carbon aerosol tracers 
      !=================================================================

      ! Chemistry for hydrophilic OC
      IF ( IDTOCPI > 0 ) THEN 
         CALL CHEM_OCPI_ADJ( STT_ADJ(:,:,:,IDTOCPI) )
         IF ( LPRT ) 
     &      CALL DEBUG_MSG( '### CHEMCARBON_ADJ: a CHEM_OCPI_ADJ' )
      ENDIF

      ! Chemistry for hydrophobic OC
      IF ( IDTOCPO > 0 ) THEN
         CALL CHEM_OCPO_ADJ( STT_ADJ(:,:,:,IDTOCPO) )
         IF ( LPRT ) 
     &      CALL DEBUG_MSG( '### CHEMCARBON_ADJ: a CHEM_OCPO_ADJ' )
      ENDIF

      ! Chemistry for hydrophilic BC
      IF ( IDTBCPI > 0 ) THEN
         CALL CHEM_BCPI_ADJ( STT_ADJ(:,:,:,IDTBCPI) )
         IF ( LPRT ) 
     &      CALL DEBUG_MSG( '### CHEMCARBON_ADJ: a CHEM_BCPI_ADJ' )
      ENDIF


      ! Chemistry for hydrophobic BC
      IF ( IDTBCPO > 0 ) THEN
         CALL CHEM_BCPO_ADJ( STT_ADJ(:,:,:,IDTBCPO) )
         IF ( LPRT ) 
     &      CALL DEBUG_MSG( '### CHEMCARBON_ADJ: a CHEM_BCPO_ADJ' )
      ENDIF


      !=================================================================
      ! Do chemistry for secondary organic aerosols 
      !=================================================================
      IF ( LSOA ) THEN

         CALL ERROR_STOP('SOA not supported yet for adjoint', 
     &                   'carbon_adj_mod.f') 

!         ! Read offline OH, NO3, O3 fields from disk
!         IF ( ITS_AN_AEROSOL_SIM() ) THEN
!
!            ! Current month
!            THISMONTH = GET_MONTH()
!
!            IF ( ITS_A_NEW_MONTH() ) THEN
!               CALL GET_GLOBAL_OH(  THISMONTH )
!               CALL GET_GLOBAL_NO3( THISMONTH )
!               CALL GET_GLOBAL_O3(  THISMONTH )
!            ENDIF
!
!            ! Compute time scaling arrays for offline OH, NO3
!            ! but only if it hasn't been done in EMISSCARBON
!            IF ( LSOA .and. ( .not. LEMIS ) ) THEN
!               CALL OHNO3TIME
!               IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARB: a OHNO3TIME' )
!            ENDIF
!         ENDIF
!
!         ! Compute SOA chemistry
!         ! NOTE: This is SOA production from the reversible mechanism only 
!         ! (tmf, 12/07/07) 
!         CALL SOA_CHEMISTRY
!         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a SOA_CHEM' )
!
!         ! If SOAG and SOAM are declared, switch on 
!         !    SOA production from dicarbonyls (tmf, 12/07/07) 
!         IF ( IDTSOAG > 0 ) THEN
!
!            ! Get grid box cloud fraction 
!            ! (tmf, 2/26/07)
!            CALL GET_VCLDF
!
!            ! Cloud uptake
!            CALL SOAG_CLOUD
!            IF ( LPRT ) 
!     &       CALL DEBUG_MSG('### CHEMCARBON: a SOAG_CLOUD')        
!
!            ! Aqueous aerosol uptake
!            CALL SOAG_LIGGIO_DIFF
!            IF ( LPRT ) 
!     &       CALL DEBUG_MSG('### CHEMCARBON: a SOAG_LIGGIO_DIFF')        
!
!         ENDIF
!
!         IF ( IDTSOAM > 0 ) THEN
!         
!            ! Get grid box cloud fraction 
!            ! (tmf, 2/26/07)
!            CALL GET_VCLDF
!
!            ! Cloud uptake
!            CALL SOAM_CLOUD
!            IF ( LPRT ) 
!     &       CALL DEBUG_MSG('### CHEMCARBON: a SOAM_CLOUD')        
!
!            ! Aqueous aerosol uptake
!            CALL SOAM_LIGGIO_DIFF
!            IF ( LPRT ) 
!     &       CALL DEBUG_MSG( '### CHEMCARBON: a SOAM_LIGGIO_DIFF' )        
!
!
!           
!         ENDIF   


      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMCARBON_ADJ

!-----------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPO_ADJ( TC_ADJ )
!
!******************************************************************************
!  Subroutine ADJ_CHEM_BCPO converts hydrophobic BC to hydrophilic BC and
!  calculates the dry deposition of hydrophobic BC adjoints. (dkh, 03/02/05)
!
!  Based on forward model by (rjp, bmy, 4/1/04,10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophobic BC tracer 
!
!  NOTES:
!  (1 ) See forward model. 
!  (2 ) Based on CHEM_OPCI from forward model. The only differences are:
!        i.   Check if ABS( CNEW ) < SMALLNUM   
!        ii.  Include STT_ADJ
!        iii. Take out ND44 stuff
!  (3 ) Updated to include adjoint of OCPO --> OCPI.  Comment out diagnostics
!        from forward model.  (dkh, 03/22/07) 
!  (4 ) Updated to GCv8 (dkh, 09/09/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE CARBON_MOD,   ONLY : DRYBCPO
      USE DAO_MOD,      ONLY : AD
      !USE DIAG_MOD,     ONLY : AD44, AD07_BC 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTBCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC_ADJ(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,       J,   L
      !TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)     
      TYPE (XPLEX)                 :: DTCHEM, FLUX, KBC, FREQ, BL_FRAC
      TYPE (XPLEX)                 :: TC0_ADJ, CNEW_ADJ, RKT, AREA_CM2
      TYPE (XPLEX),  PARAMETER     :: BC_LIFE = xplex(1.15D0,0d0)

      !=================================================================
      ! CHEM_BCPO_ADJ begins here!
      !=================================================================

      ! Return if BCPO isn't defined
      IF ( IDTBCPO == 0 .or. DRYBCPO == 0 ) RETURN

      ! Initialize
      KBC    = 1.D0 / ( 86400d0 * BC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

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
!$OMP+PRIVATE( I, J, L, TC0_ADJ, FREQ, BL_FRAC, RKT, CNEW_ADJ )
!$OMP+PRIVATE( AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Store new concentration back into tracer array
         CNEW_ADJ = TC_ADJ(I,J,L)

         ! Zero drydep freq
         FREQ = 0d0

         ! Fraction of box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )


         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! BC drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYBCPO) * BL_FRAC

         ENDIF

         ! Prevent underflow condition
         !IF ( ABS( CNEW ) < SMALLNUM ) CNEW = 0d0

         ! Amount of BCPO converted to BCPI [kg/timestep]
         ! fwd code:
         !BCCONV(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )
         TC0_ADJ  =   BCCONV_ADJ(I,J,L) * KBC / ( KBC + FREQ )
         ! CNEW_ADJ is calculated as:
         ! CNEW_ADJ = CNEW_ADJ - BCCONV_ADJ(I,J,L) * KBC / ( KBC + FREQ )
         ! same thing, but faster:
         CNEW_ADJ = CNEW_ADJ - TC0_ADJ

         ! Amount of BCPO left after chemistry and drydep [kg]
         RKT  = ( KBC + FREQ ) * DTCHEM
         ! fwd code:
         !CNEW = TC0 * EXP( -RKT )
         TC0_ADJ = TC0_ADJ + CNEW_ADJ * EXP( -RKT )

!         !==============================================================
!         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN
!
!             ! Surface area [cm2]
!             AREA_CM2 = GET_AREA_CM2( J )
!
!             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]  
!             ! XNUMOL is the ratio [molec tracer/kg tracer]   
!             FLUX     = TC0 - CNEW - BCCONV(I,J,L) 
!             FLUX     = FLUX * XNUMOL(IDTBCPO) / ( DTCHEM * AREA_CM2 )
!
!             ! Store in ND44_TMP as a placeholder
!             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
!         ENDIF
!
!         !==============================================================
!         ! ND07 diagnostic: H-philic BC from H_phobic BC [kg/timestep]
!         !==============================================================
!         IF ( ND07 > 0 .and. L <= LD07 ) THEN
!             AD07_BC(I,J,L) = AD07_BC(I,J,L) + BCCONV(I,J,L)
!         ENDIF

         ! Initial BC mass [kg]
         ! fwd code:
         !TC0  = TC(I,J,L)
         TC_ADJ(I,J,L) = TC0_ADJ

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

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
!         DO L = 1, LLPAR
!            AD44(I,J,DRYBCPO,1) = AD44(I,J,DRYBCPO,1) + ND44_TMP(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF      

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! fwd code:
         !BCCONV(I,J,L) = 0d0
         BCCONV_ADJ(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
!         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_BCPO_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPI_ADJ( TC_ADJ )
!
!******************************************************************************
!  Subroutine CHEM_BCPI_ADJ calculates dry deposition of hydrophilic BC adjoint
!  (dkh, 03/02/05)
!
!  Based on forward code by (rjp, bmy, 4/1/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophilic BC adjoint
! 
!  NOTES:
!  (1 ) Based on CHEM_BCPI from forward model.  The only differences are:
!        i.   Check if ABS( CNEW ) < SMALLNUM
!        ii.  Return if IDADJBCPI is 0
!        iii. Include ADJ_STT
!  (2 ) Updated to include adjoint of BCPO --> BCPI.  Comment out diagnostics
!        from forward model.  (dkh, 03/22/07) 
!  (3 ) Updated to GCv8 (dkh, 09/09/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE CARBON_MOD,   ONLY : DRYBCPI
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP      
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTBCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC_ADJ(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, L
      TYPE (XPLEX)             :: DTCHEM, FLUX, BL_FRAC, AREA_CM2, FREQ
      TYPE (XPLEX)             :: TC0_ADJ, CNEW_ADJ, CCV_ADJ
!      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_BCPI_ADJ begins here!
      !=================================================================

      ! Return if BCPI isn't defined
      IF ( IDTBCPI == 0 .or. DRYBCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Zero out the BCCONV_ADJ array for the next iteration
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! fwd code:
         !BCCONV(I,J,L) = 0.d0
         BCCONV_ADJ(I,J,L) = 0.d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!      ! Initialize for ND44 diagnostic
!      IF ( ND44 > 0 ) THEN
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            ND44_TMP(I,J,L) = 0d0
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0_ADJ, CCV_ADJ, FREQ, BL_FRAC)
!$OMP+PRIVATE( CNEW_ADJ, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Save new concentration of H-philic IC in tracer array
         CNEW_ADJ = TC_ADJ(I,J,L) 

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

            ! note -- this was already commented out of fwd code 
            ! Comment out for now
            !CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !     + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            ! fwd code:
            !CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )
            TC0_ADJ = CNEW_ADJ * EXP( -FREQ * DTCHEM )
            ! adjoint for CCV_ADJ is:
            ! CCV_ADJ = CNEW_ADJ * EXP( -FREQ * DTCHEM )
            ! or, same but faster:
            CCV_ADJ = TC0_ADJ

!            !===========================================================
!            ! ND44 diagnostic: drydep flux [atoms C/cm2/s]
!            !===========================================================
!            IF ( ND44 > 0 .and. FREQ > 0d0 ) THEN
!  
!               ! Surface area [cm2]
!               AREA_CM2 = GET_AREA_CM2( J )
!
!               ! Convert drydep loss from [kg/timestep] to [molec/cm2/s]
!               FLUX = ( TC0 + CCV - CNEW ) 
!               FLUX = FLUX * XNUMOL(IDTBCPI) / ( AREA_CM2 * DTCHEM )
!             
!               ! Store in ND44_TMP as a placeholder
!               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
!            ENDIF

         ELSE

            ! Otherwise, omit the exponential to save on clock cycles
            ! fwd code:
            !CNEW = TC0 + CCV
            TC0_ADJ = CNEW_ADJ
            CCV_ADJ = CNEW_ADJ

         ENDIF
      
         ! Prevent underflow condition
         !IF ( ABS( CNEW ) < SMALLNUM ) CNEW = 0d0

         ! H-philic BC that used to be H-phobic BC [kg]
         ! fwd code:
         !CCV = BCCONV(I,J,L)
         BCCONV_ADJ(I,J,L) = CCV_ADJ
         
         ! Initial H-philic BC [kg]
         ! fwd code:
         !TC0 = TC(I,J,L)
         TC_ADJ(I,J,L) = TC0_ADJ

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  


!      !=================================================================
!      ! ND44: Sum drydep fluxes by level into the AD44 array in
!      ! order to ensure that  we get the same results w/ sp or mp 
!      !=================================================================
!      IF ( ND44 > 0 ) THEN 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!         DO L = 1, LLPAR
!            AD44(I,J,DRYBCPI,1) = AD44(I,J,DRYBCPI,1) + ND44_TMP(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPI_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_OCPI_ADJ( TC_ADJ )
!
!******************************************************************************
!  Subroutine CHEM_OCPI_ADJ calculates dry deposition of hydrophilic OC adjoint
!  (dkh, 03/02/05)
!
!  Based on forward code by (rjp, bmy, 4/1/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!  (1 ) Based on CHEM_OPCI from forward model. The only differences are:
!        i.   Check if ABS( CNEW ) < SMALLNUM   
!        ii.  Return if IDADJOPCI is 0
!        iii. Include ADJ_STT
!  (2 ) Updated to include adjoint of OCPO --> OCPI.  Comment out diagnostics
!        from forward model.  (dkh, 03/22/07) 
!  (4 ) Updated to GCv8 (dkh, 09/09/09) 
!  (5 ) BUG FIX: now declare BL_FRAC thread private (dkh, 07/30/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE CARBON_MOD,   ONLY : DRYOCPI
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTOCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! ND44

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC_ADJ(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      TYPE (XPLEX)                 :: DTCHEM, FLUX, BL_FRAC, AREA_CM2
      TYPE (XPLEX)                 :: TC0_ADJ, CNEW_ADJ, CCV_ADJ, FREQ
!      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_OCPI_ADJ begins here!
      !=================================================================
      IF ( IDTOCPI == 0 .or. DRYOCPI == 0  ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!      ! Initialize for drydep diagnostic
!      IF ( ND44 > 0 ) THEN
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            ND44_TMP(I,J,L) = 0d0
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF

      !=================================================================
      ! Zero OCCONV_ADJ array for next timestep
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! fwd code:
         !OCCONV(I,J,L) = 0d0
         OCCONV_ADJ(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! BUG FIX: BL_FRAC needs to be thread private (dkh, 07/30/10) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0_ADJ, CCV_ADJ, FREQ, CNEW_ADJ )
!$OMP+PRIVATE( AREA_CM2, FLUX )
!$OMP+PRIVATE( BL_FRAC )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Store modified concentration back in tracer array [kg]
         ! fwd code: 
         !TC(I,J,L) = CNEW
         CNEW_ADJ = TC_ADJ(I,J,L)

         ! dkh -- don't take adjoint of this.  It would require
         ! recalculation of fwd CNEW -- probably not worth while. 
         ! Prevent underflow condition
         !IF ( ABS( CNEW ) < SMALLNUM ) CNEW = 0d0

         ! Fraction of grid box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Recalculate drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DRYOCPI) * BL_FRAC

            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT))
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! dkh -- this was already commented out of fwd code
            ! CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !       + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]D
            ! fwd code:
            !CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )
            TC0_ADJ = CNEW_ADJ * EXP( -FREQ * DTCHEM )
            ! adjoint code for CCV is:
            ! CCV_ADJ = CNEW_ADJ * EXP( -FREQ * DTCHEM )
            ! same thing, except faster:
            CCV_ADJ = TC0_ADJ

!            !===========================================================
!            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
!            !===========================================================
!            IF ( ND44 > 0 ) THEN
!
!               ! Surface area [cm2]
!               AREA_CM2 = GET_AREA_CM2( J )
!
!               ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
!               FLUX = ( TC0 + CCV - CNEW ) 
!               FLUX = FLUX * XNUMOL(IDTOCPI) / ( AREA_CM2 * DTCHEM )
!             
!               ! Store in ND44_TMP as a placeholder
!               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
!            ENDIF

         ELSE

            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            ! fwd code:
            !CNEW = TC0 + CCV
            TC0_ADJ = CNEW_ADJ
            CCV_ADJ = CNEW_ADJ

         ENDIF
      
         ! Initial H-philic OC [kg]
         ! fwd code:
         !TC0 = TC(I,J,L)
         TC_ADJ(I,J,L) = TC0_ADJ

         ! H-philic OC that used to be H-phobic OC [kg]
         ! fwd code:
         !CCV = OCCONV(I,J,L)
         OCCONV_ADJ(I,J,L) = CCV_ADJ

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

!      !=================================================================
!      ! ND44: Sum drydep fluxes by level into the AD44 array in
!      ! order to ensure that  we get the same results w/ sp or mp 
!      !=================================================================
!      IF ( ND44 > 0 ) THEN 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!         DO L = 1, LLPAR
!            AD44(I,J,DRYOCPI,1) = AD44(I,J,DRYOCPI,1) + ND44_TMP(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF    

      ! Return to calling program
      END SUBROUTINE CHEM_OCPI_ADJ


!------------------------------------------------------------------------------

      SUBROUTINE CHEM_OCPO_ADJ( TC_ADJ )
!
!******************************************************************************
!  Subroutine CHEM_OCPO_ADJ converts adjoint of hydrophobic OC to hydrophilic OC and
!  calculates the dry deposition of hydrophobic OC. (dkh, 03/02/05)
!
!  Based on forward model by (rjp, bmy, 4/1/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Array of hydrophobic OC tracer [kg]
!
!  NOTES:
!  (1 ) Based on CHEM_OCPO from forward model. The only differences are:
!        i.   Check if ABS( CNEW ) < SMALLNUM
!        ii.  Return if IDADJOCPO is 0
!        iii. Include ADJ_STT
!  (2 ) Updated to include adjoint of OCPO --> OCPI.  Comment out diagnostics
!        from forward model.  (dkh, 03/22/07) 
!  (3 ) Updated to GCv8 (dkh, 09/09/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE CARBON_MOD,   ONLY : DRYOCPO
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_OC 
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTOCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT) :: TC_ADJ(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
!      TYPE (XPLEX)                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                 :: DTCHEM, FLUX, KOC,  BL_FRAC
      TYPE (XPLEX)           :: TC0_ADJ,   FREQ, CNEW_ADJ, RKT, AREA_CM2
      TYPE (XPLEX),  PARAMETER     :: OC_LIFE = xplex(1.15D0,0d0)

      !=================================================================
      ! CHEM_OCPO_ADJ begins here!
      !=================================================================

      ! Return if OCPO isn't defined
      IF ( IDTOCPO == 0 .or. DRYOCPO == 0 ) RETURN

      ! Initialize
      KOC    = 1.D0 / ( 86400d0 * OC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0


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
!$OMP+PRIVATE( I, J, L, TC0_ADJ, FREQ, BL_FRAC, RKT )
!$OMP+PRIVATE( CNEW_ADJ, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Store modified OC concentration back in tracer array
         ! fwd code: 
         !TC(I,J,L) = CNEW
         CNEW_ADJ = TC_ADJ(I,J,L)

         ! Zero drydep freq 
         FREQ = 0d0

         ! Fraction of box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! OC drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYOCPO) * BL_FRAC

         ENDIF

         ! Amount of OCPO converted to OCPI [kg/timestep]
         ! fwd code:
         !OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )
         TC0_ADJ  =   OCCONV_ADJ(I,J,L) * KOC / ( KOC + FREQ )
         ! adjoint code is:
         ! CNEW_ADJ = CNEW_ADJ - OCCONV_ADJ(I,J,L) * KOC / ( KOC + FREQ )
         ! same thing, except faster:
         CNEW_ADJ = CNEW_ADJ - TC0_ADJ

         ! Amount of OCPO left after chemistry and drydep [kg]
         RKT  = ( KOC + FREQ ) * DTCHEM
         ! fwd code:
         !CNEW = TC0 * EXP( -RKT )
         TC0_ADJ = TC0_ADJ + CNEW_ADJ * EXP( -RKT )

         ! dkh -- don't take adjoint of this
         ! Prevent underflow condition
         !IF ( ABS( CNEW ) < SMALLNUM ) CNEW = 0d0

!         !==============================================================
!         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN
!
!             ! Surface area [cm2]
!             AREA_CM2 = GET_AREA_CM2( J )
!
!             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
!             ! XNUMOL is the ratio [molec tracer/kg tracer]     
!             FLUX     = TC0 - CNEW - OCCONV(I,J,L)
!             FLUX     = FLUX * XNUMOL(IDTOCPO) / ( DTCHEM * AREA_CM2 )
!
!             ! Store in ND44_TMP as a placeholder
!             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
!         ENDIF
!
!         !==============================================================
!         ! ND07 diagnostic: H-Philic OC from H-phobic [kg/timestep]
!         !==============================================================
!         IF ( ND07 > 0 .and. L <= LD07 ) THEN
!            AD07_OC(I,J,L) = AD07_OC(I,J,L) + OCCONV(I,J,L) 
!         ENDIF

         ! Initial OC [kg]
         ! fwd code:
         !TC0  = TC(I,J,L)
         TC_ADJ(I,J,L) = TC0_ADJ

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

!      !=================================================================
!      ! ND44: Sum drydep fluxes by level into the AD44 array in
!      ! order to ensure that  we get the same results w/ sp or mp 
!      !=================================================================
!      IF ( ND44 > 0 ) THEN 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!         DO L = 1, LLPAR
!            AD44(I,J,DRYOCPO,1) = AD44(I,J,DRYOCPO,1) + ND44_TMP(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!      ENDIF   

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! fwd code:
         !OCCONV(I,J,L) = 0d0
         OCCONV_ADJ(I,J,L) = 0d0

!         ! Initialize for drydep diagnostic
!         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_OCPO_ADJ
!------------------------------------------------------------------------------

      SUBROUTINE EMISSCARBON_ADJ
!
!******************************************************************************
!  Subroutine EMISSCARBON_ADJ is the adjoint of EMISSCARBON.   (dkh, 04/26/06)  
  
!  It is based on the forward model subroutine EMISSCARBON which is the interface 
!  between the GEOS-CHEM modeland the CARBONACEOUS AEROSOL emissions 
!  (rjp, bmy, 1/24/02, 9/25/06)
!
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 11/10/09) 
!  (2 ) Add LEMS_ABS option (dkh, 02/17/11) 
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,          ONLY : AD07
      USE DAO_MOD,           ONLY : PBL
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_MOD,       ONLY : LSOA,      LPRT
      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : STT
      USE GFED2_BIOMASS_MOD, ONLY : GFED2_IS_NEW
      !USE TRACERID_MOD

      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EBCPI_an, IDADJ_EBCPO_an
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EOCPI_an, IDADJ_EOCPO_an
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EBCPI_bb, IDADJ_EBCPO_bb
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EOCPI_bb, IDADJ_EOCPO_bb
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EBCPI_bf, IDADJ_EBCPO_bf
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EOCPI_bf, IDADJ_EOCPO_bf
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_ADJ
      USE CARBON_MOD,        ONLY : ANTH_ORGC,      ANTH_BLKC
      USE CARBON_MOD,        ONLY : BIOB_ORGC,      BIOB_BLKC
      USE CARBON_MOD,        ONLY : BIOF_ORGC,      BIOF_BLKC
      USE LOGICAL_ADJ_MOD,   ONLY : LEMS_ABS 

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! ND07

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J, MONTH, N
      TYPE (XPLEX)               :: BCSRC_ADJ(IIPAR,JJPAR,2)
      TYPE (XPLEX)               :: OCSRC_ADJ(IIPAR,JJPAR,2)

      !=================================================================
      ! EMISSCARBON_ADJ begins here!
      !=================================================================      

      ! fwd code:
      !CALL EMITHIGH( BCSRC, OCSRC )
      ! adj code: 
      CALL EMITHIGH_ADJ( BCSRC_ADJ, OCSRC_ADJ )
      IF ( LPRT ) 
     &    CALL DEBUG_MSG( '### EMISCARB_ADJ: after EMITHIGH_ADJ' )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! fwd code: 
         !   ! Total HYDROPHILIC BC source [kg]
         !  BCSRC(I,J,1)= ANTH_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_an)
         !              + BIOF_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bf)
         !              + BIOB_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bb)
         ! adj code:
         EMS_SF_ADJ(I,J,1,IDADJ_EBCPI_an) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EBCPI_an) 
     &    + ANTH_BLKC(I,J,1) * BCSRC_ADJ(I,J,1)

         EMS_SF_ADJ(I,J,1,IDADJ_EBCPI_bf) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EBCPI_bf) 
     &    + BIOF_BLKC(I,J,1) * BCSRC_ADJ(I,J,1)

         EMS_SF_ADJ(I,J,1,IDADJ_EBCPI_bb) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EBCPI_bb) 
     &    + BIOB_BLKC(I,J,1) * BCSRC_ADJ(I,J,1)

         ! fwd code:
         !   ! Total HYDROPHOBIC BC source [kg]
         !  BCSRC(I,J,2)= ANTH_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPI_an)
         !              + BIOF_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bf)
         !              + BIOB_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bb)
         ! adj code:
         EMS_SF_ADJ(I,J,1,IDADJ_EBCPO_an) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EBCPO_an) 
     &    + ANTH_BLKC(I,J,2) * BCSRC_ADJ(I,J,2)

         EMS_SF_ADJ(I,J,1,IDADJ_EBCPO_bf) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EBCPO_bf) 
     &    + BIOF_BLKC(I,J,2) * BCSRC_ADJ(I,J,2)

         EMS_SF_ADJ(I,J,1,IDADJ_EBCPO_bb) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EBCPO_bb) 
     &    + BIOB_BLKC(I,J,2) * BCSRC_ADJ(I,J,2)

         IF ( LSOA ) THEN

             CALL ERROR_STOP('LSOA not supported yet', 
     &                       'carbon_adj_mod.f')
!            ! Total HYDROPHILIC OC source [kg]
!            ! (Don't use archived TERP_ORGC if LSOA=T)
!            OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
!     &                     BIOF_ORGC(I,J,1) + 
!     &                     BIOB_ORGC(I,J,1)

         ELSE

            ! fwd code:
            !  ! Total HYDROPHILIC OC source [kg]
              ! (Use archived TERP_ORGC for if LSOA=F)
            !  OCSRC(I,J,1)
            !     = ANTH_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_an)
            !     + BIOF_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bf)
            !     + BIOB_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bb)
            !     + TERP_ORGC(I,J)
            ! adj code:
            EMS_SF_ADJ(I,J,1,IDADJ_EOCPI_an) 
     &         = EMS_SF_ADJ(I,J,1,IDADJ_EOCPI_an) 
     &         + ANTH_ORGC(I,J,1) * OCSRC_ADJ(I,J,1)

            EMS_SF_ADJ(I,J,1,IDADJ_EOCPI_bf) 
     &         = EMS_SF_ADJ(I,J,1,IDADJ_EOCPI_bf) 
     &         + BIOF_ORGC(I,J,1) * OCSRC_ADJ(I,J,1)

            EMS_SF_ADJ(I,J,1,IDADJ_EOCPI_bb) 
     &         = EMS_SF_ADJ(I,J,1,IDADJ_EOCPI_bb) 
     &         + BIOB_ORGC(I,J,1) * OCSRC_ADJ(I,J,1)

         ENDIF

         ! fwd code:
         !    ! Total HYDROPHOBIC OC source [kg]
         !   OCSRC(I,J,2)       
         !      = ANTH_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_an)
         !      + BIOF_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bf)
         !      + BIOB_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bb)
         ! adj code:
         EMS_SF_ADJ(I,J,1,IDADJ_EOCPO_an) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EOCPO_an) 
     &    + ANTH_ORGC(I,J,2) * OCSRC_ADJ(I,J,2)

         EMS_SF_ADJ(I,J,1,IDADJ_EOCPO_bf) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EOCPO_bf) 
     &    + BIOF_ORGC(I,J,2) * OCSRC_ADJ(I,J,2)

         EMS_SF_ADJ(I,J,1,IDADJ_EOCPO_bb) 
     &    = EMS_SF_ADJ(I,J,1,IDADJ_EOCPO_bb) 
     &    + BIOB_ORGC(I,J,2) * OCSRC_ADJ(I,J,2)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Optional diagnostic -- also save out the emissions adjoints (dkh, 02/17/11) 
      ! (absolute sensitivities per emissions rather than per scaling factor)
      IF ( LEMS_ABS ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! fwd code: 
            !   ! Total HYDROPHILIC BC source [kg]
            !  BCSRC(I,J,1)= ANTH_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_an)
            !              + BIOF_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bf)
            !              + BIOB_BLKC(I,J,1) * EMS_SF(I,J,1,IDADJ_EBCPI_bb)
            ! adj code:
            EMS_ADJ(I,J,1,IDADJ_EBCPI_an) 
     &       = EMS_ADJ(I,J,1,IDADJ_EBCPI_an) + BCSRC_ADJ(I,J,1)

            EMS_ADJ(I,J,1,IDADJ_EBCPI_bf) 
     &       = EMS_ADJ(I,J,1,IDADJ_EBCPI_bf) + BCSRC_ADJ(I,J,1)

            EMS_ADJ(I,J,1,IDADJ_EBCPI_bb) 
     &       = EMS_ADJ(I,J,1,IDADJ_EBCPI_bb) + BCSRC_ADJ(I,J,1)

            ! fwd code:
            !   ! Total HYDROPHOBIC BC source [kg]
            !  BCSRC(I,J,2)= ANTH_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPI_an)
            !              + BIOF_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bf)
            !              + BIOB_BLKC(I,J,2) * EMS_SF(I,J,1,IDADJ_EBCPO_bb)
            ! adj code:
            EMS_ADJ(I,J,1,IDADJ_EBCPO_an) 
     &       = EMS_ADJ(I,J,1,IDADJ_EBCPO_an) + BCSRC_ADJ(I,J,2)

            EMS_ADJ(I,J,1,IDADJ_EBCPO_bf) 
     &       = EMS_ADJ(I,J,1,IDADJ_EBCPO_bf) + BCSRC_ADJ(I,J,2)

            EMS_ADJ(I,J,1,IDADJ_EBCPO_bb) 
     &       = EMS_ADJ(I,J,1,IDADJ_EBCPO_bb) + BCSRC_ADJ(I,J,2)

            IF ( LSOA ) THEN

                CALL ERROR_STOP('LSOA not supported yet', 
     &                       'carbon_adj_mod.f')
!               ! Total HYDROPHILIC OC source [kg]
!               ! (Don't use archived TERP_ORGC if LSOA=T)
!               OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
!     &                        BIOF_ORGC(I,J,1) + 
!     &                        BIOB_ORGC(I,J,1)

            ELSE

               ! fwd code:
               !  ! Total HYDROPHILIC OC source [kg]
               ! (Use archived TERP_ORGC for if LSOA=F)
               !  OCSRC(I,J,1)
               !     = ANTH_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_an)
               !     + BIOF_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bf)
               !     + BIOB_ORGC(I,J,1) * EMS_SF(I,J,1,IDADJ_EOCPI_bb)
               !     + TERP_ORGC(I,J)
               ! adj code:
               EMS_ADJ(I,J,1,IDADJ_EOCPI_an) 
     &            = EMS_ADJ(I,J,1,IDADJ_EOCPI_an) + OCSRC_ADJ(I,J,1)

               EMS_ADJ(I,J,1,IDADJ_EOCPI_bf) 
     &            = EMS_ADJ(I,J,1,IDADJ_EOCPI_bf) + OCSRC_ADJ(I,J,1)

               EMS_ADJ(I,J,1,IDADJ_EOCPI_bb) 
     &            = EMS_ADJ(I,J,1,IDADJ_EOCPI_bb) + OCSRC_ADJ(I,J,1)

            ENDIF

            ! fwd code:
            !    ! Total HYDROPHOBIC OC source [kg]
            !   OCSRC(I,J,2)       
            !      = ANTH_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_an)
            !      + BIOF_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bf)
            !      + BIOB_ORGC(I,J,2) * EMS_SF(I,J,1,IDADJ_EOCPO_bb)
            ! adj code:
            EMS_ADJ(I,J,1,IDADJ_EOCPO_an) 
     &       = EMS_ADJ(I,J,1,IDADJ_EOCPO_an) + OCSRC_ADJ(I,J,2)

            EMS_ADJ(I,J,1,IDADJ_EOCPO_bf) 
     &       = EMS_ADJ(I,J,1,IDADJ_EOCPO_bf) + OCSRC_ADJ(I,J,2)

            EMS_ADJ(I,J,1,IDADJ_EOCPO_bb) 
     &       = EMS_ADJ(I,J,1,IDADJ_EOCPO_bb) + OCSRC_ADJ(I,J,2)

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF 

      ! Return to calling program
      END SUBROUTINE EMISSCARBON_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE EMITHIGH_ADJ( BCSRC_ADJ, OCSRC_ADJ )
!
!******************************************************************************
!  Subroutine EMITHIGH_ADJ is the adjoint of EMITHIGH (dkh, 04/26/06)
!  
!  Based on forward routine EMITHIGHT that mixes tracer completely from the 
!  surface to the PBL top. (rjp, bmy, 4/2/04, 2/17/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) BCSRC (TYPE (XPLEX)) : Array which holds Total BC (H-phobic & H-philic)
!  (2 ) OCSRC (TYPE (XPLEX)) : Array which holds Total OC (H-phobic & H-philic)
!
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 11/11/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_OF_PBL,  GET_PBL_MAX_L
      USE TRACER_MOD,     ONLY : STT
      USE TRACERID_MOD,   ONLY : IDTBCPI, IDTBCPO, IDTOCPI, IDTOCPO
      USE TRACERID_MOD,   ONLY : IDTALPH, IDTLIMO, IDTALCO

#     include "CMN_SIZE"       ! Size parameters

      ! Arguments
      TYPE (XPLEX), INTENT(OUT)     :: BCSRC_ADJ(IIPAR,JJPAR,2)
      TYPE (XPLEX), INTENT(OUT)     :: OCSRC_ADJ(IIPAR,JJPAR,2)

      ! Local variables
      LOGICAL                 :: IS_BCPO, IS_OCPO, IS_BCPI, IS_OCPI
      LOGICAL                 :: IS_ALPH, IS_LIMO, IS_ALCO
      INTEGER                 :: I,       J,       L,       PBL_MAX
      TYPE (XPLEX)                  :: F_OF_PBL

      !=================================================================
      ! EMITHIGH_ADJ begins here!
      !=================================================================

      ! initialize
      BCSRC_ADJ = 0d0 
      OCSRC_ADJ = 0d0 
      
      ! Define logical flags for expediency
      IS_BCPI = ( IDTBCPI > 0 )
      IS_OCPI = ( IDTOCPI > 0 ) 
      IS_BCPO = ( IDTBCPO > 0 )
      IS_OCPO = ( IDTOCPO > 0 )
      IF ( IDTALPH > 0 ) 
     &   CALL ERROR_STOP( 'ALPH not supported', 'carbon_adj_mod')
      IF ( IDTLIMO > 0 ) 
     &   CALL ERROR_STOP( 'LIMO not supported', 'carbon_adj_mod')
      IF ( IDTALCO > 0 ) 
     &   CALL ERROR_STOP( 'ALCO not supported', 'carbon_adj_mod')

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
            ! fwd code:
            !STT(I,J,L,IDTBCPI) = STT(I,J,L,IDTBCPI) + 
            !                     ( F_OF_PBL * BCSRC(I,J,1) )
            ! adj code: 
            BCSRC_ADJ(I,J,1) = BCSRC_ADJ(I,J,1) 
     &                       + F_OF_PBL * STT_ADJ(I,J,L,IDTBCPI) 
            
         ENDIF

         ! Hydrophilic ORGANIC CARBON
         IF ( IS_OCPI ) THEN
            ! fwd code:
            !STT(I,J,L,IDTOCPI) = STT(I,J,L,IDTOCPI) + 
            !                     ( F_OF_PBL * OCSRC(I,J,1) )
            ! adj code: 
            OCSRC_ADJ(I,J,1) = OCSRC_ADJ(I,J,1)
     &                       + F_OF_PBL * STT_ADJ(I,J,L,IDTOCPI) 
         ENDIF
            
         ! Hydrophobic BLACK CARBON
         IF ( IS_BCPO ) THEN
            ! fwd code:
            !STT(I,J,L,IDTBCPO) = STT(I,J,L,IDTBCPO) + 
            !                     ( F_OF_PBL * BCSRC(I,J,2) )
            ! adj code:
            BCSRC_ADJ(I,J,2) = BCSRC_ADJ(I,J,2)
     &                       + F_OF_PBL * STT_ADJ(I,J,L,IDTBCPO) 
         ENDIF

         ! Hydrophobic ORGANIC CARBON
         IF ( IS_OCPO ) THEN
            ! fwd code:
            !STT(I,J,L,IDTOCPO) = STT(I,J,L,IDTOCPO) + 
            !                     ( F_OF_PBL * OCSRC(I,J,2) )
            ! adj code:
            OCSRC_ADJ(I,J,2) = OCSRC_ADJ(I,J,2)
     &                       + F_OF_PBL * STT_ADJ(I,J,L,IDTOCPO) 
         ENDIF

         ! remaining species not yet included in adjoint 
!         ! ALPHA-PINENE
!         IF ( IS_ALPH ) THEN
!            STT(I,J,L,IDTALPH) = STT(I,J,L,IDTALPH) + 
!     &                           ( F_OF_PBL * BIOG_ALPH(I,J) )
!         ENDIF
!
!         ! LIMONENE
!         IF ( IS_LIMO ) THEN
!            STT(I,J,L,IDTLIMO) = STT(I,J,L,IDTLIMO) + 
!     &                           ( F_OF_PBL * BIOG_LIMO(I,J) )
!
!            ORVC_TERP(I,J,L)   = ORVC_TERP(I,J,L) + 
!     &                           ( F_OF_PBL * BIOG_TERP(I,J) )
!         ENDIF
!
!         ! ALCOHOL and SESQTERPENE (not a tracer)
!         IF ( IS_ALCO ) THEN
!            STT(I,J,L,IDTALCO) = STT(I,J,L,IDTALCO) + 
!     &                           ( F_OF_PBL * BIOG_ALCO(I,J) )
!               
!            ORVC_SESQ(I,J,L)   = ORVC_SESQ(I,J,L) + 
!     &                           ( F_OF_PBL * BIOG_SESQ(I,J) )
!         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMITHIGH_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE INIT_CARBON_ADJ
!     
!******************************************************************************
!  Subroutine INIT_CARBON_ADJ initializes all module arrays (rjp, bmy, 4/1/04)
!     
!  NOTES:
!******************************************************************************
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" ! Size parameters
      
      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS
      
      !=================================================================
      ! INIT_CARBON_ADJ begins here!
      !=================================================================
      
      ! Return if we already allocated arrays
      IF ( IS_INIT ) RETURN

      ALLOCATE( BCCONV_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONV_ADJ' )
      BCCONV_ADJ = 0d0
      
      ALLOCATE( OCCONV_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV_ADJ' )
      OCCONV_ADJ = 0d0

      ! Reset IS_INIT
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_CARBON_ADJ

!------------------------------------------------------------------------------
      SUBROUTINE CLEANUP_CARBON_ADJ
!
!******************************************************************************
!  Subroutine CLEANUP_CARBON_ADJ deallocates all module arrays (rjp, bmy, 4/1/04)
!        
!  NOTES:
!******************************************************************************
!        
      !=================================================================
      ! CLEANUP_CARBON_ADJ begins here!
      !=================================================================
      IF ( ALLOCATED( BCCONV_ADJ ) ) DEALLOCATE( BCCONV_ADJ )
      IF ( ALLOCATED( OCCONV_ADJ ) ) DEALLOCATE( OCCONV_ADJ )

      ! Return to calling program
      END SUBROUTINE CLEANUP_CARBON_ADJ

!------------------------------------------------------------------------------


      END MODULE CARBON_ADJ_MOD 
