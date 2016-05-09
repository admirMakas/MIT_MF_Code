! $Id: convection_adj_mod.f,v 1.5 2010/04/25 17:18:58 daven Exp $
      MODULE CONVECTION_ADJ_MOD
!
!******************************************************************************
!  Module CONVECTION_MOD contains routines which select the proper convection
!  code for GEOS-3, GEOS-4, GEOS-5, or GCAP met field data sets. 
!  (bmy, 6/28/03, 1/31/08)
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CONVECTION       : Wrapper routine, chooses correct convection code
!  (2 ) DO_GEOS4_CONVECT    : Calls GEOS-4 convection routines 
!  (3 ) DO_GCAP_CONVECT     : Calls GCAP convection routines
!  (4 ) NFCLDMX             : Convection routine for GEOS-3 and GEOS-5 met
!
!  GEOS-CHEM modules referenced by convection_mod.f
!  ============================================================================
!  (1 ) dao_mod.f           : Module w/ containing arrays for DAO met fields   
!  (2 ) diag_mod.f          : Module w/ GEOS-Chem diagnostic arrays
!  (3 ) fvdas_convect_mod.f : Module w/ convection code for fvDAS met fields
!  (4 ) grid_mod.f          : Module w/ horizontal grid information
!  (5 ) logical_mod.f       : Module w/ GEOS-Chem logical switches
!  (6 ) ocean_mercury_mod.f : Module w/ routines for Hg(0) ocean flux
!  (7 ) pressure_mod.f      : Module w/ routines to compute P(I,J,L)
!  (8 ) time_mod.f          : Module w/ routines for computing time
!  (9 ) tracer_mod.f        : Module w/ GEOS-Chem tracer array STT etc
!  (10) tracerid_mod.f      : Module w/ GEOS-Chem tracer ID flags etc
!  (11) wetscav_mod.f       : Module w/ routines for wetdep/scavenging
!
!  NOTES:
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also now references
!        "error_mod.f".  Now make F in routine NFCLDMX a 4-D array to avoid
!        memory problems on the Altix. (bmy, 1/27/04)
!  (2 ) Bug fix: Now pass NTRACE elements of TCVV to FVDAS_CONVECT in routine 
!        DO_CONVECTION (bmy, 2/23/04)  
!  (3 ) Now references "logical_mod.f" and "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Now also references "ocean_mercury_mod.f" and "tracerid_mod.f" 
!        (sas, bmy, 1/19/05)
!  (5 ) Now added routines DO_GEOS4_CONVECT and DO_GCAP_CONVECT by breaking 
!        off code from DO_CONVECTION, in order to implement GCAP convection
!        in a much cleaner way. (swu, bmy, 5/25/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Shut off scavenging in shallow convection for GCAP (swu, bmy, 11/1/05)
!  (8 ) Modified for tagged Hg simulation (cdh, bmy, 1/6/06)
!  (9 ) Bug fix: now only call ADD_Hg2_WD if LDYNOCEAN=T (phs, 2/8/07)
!  (10) Fix for GEOS-5 met fields in routine NFCLDMX (swu, 8/15/07)
!  (11) Resize DTCSUM array in NFCLDMX to save memory (bmy, 1/31/08)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "convection_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: DO_CONVECTION_ADJ

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_CONVECTION_ADJ
!
!******************************************************************************
!  Subroutine DO_CONVECTION_ADJ calls the adjoint of the appropriate 
!  convection driver program for different met field data sets. 
!  Based on forward code (swu, bmy, 5/25/05, 2/8/07). (ks,mak,dkh, 08/25/09) 
!
!  NOTES:
!  (1 ) Updated for GCv8 adjoint (dkh, 08/25/09) 
!        
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE DAO_MOD,         ONLY : CLDMAS,    CMFMC, DTRAIN
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE TRACER_MOD,      ONLY : N_TRACERS, TCVV,  STT
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM
      USE WETSCAV_MOD,     ONLY : H2O2s, SO2s
      USE WETSCAV_MOD,     ONLY : RESTORE_CONV


#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER :: I, J, L, N


#if   defined( GCAP ) 

      !-------------------------
      ! GCAP met fields
      !-------------------------

      ! Call GEOS-4 driver routine
      !CALL DO_GCAP_CONVECT
      CALL ERROR_STOP( 'GCAP not supported for adjoint', 
     &                 'convection_adj_mod.f'            )

#elif defined( GEOS_4 )

      !-------------------------
      ! GEOS-4 met fields
      !-------------------------

      ! Call GEOS-4 driver routine
      CALL DO_GEOS4_CONVECT_ADJ

#elif defined( GEOS_5 )

      !-------------------------
      ! GEOS-5 met fields
      !-------------------------

      ! Restore checkpted values of H2O2s and SO2s (dkh, 11/22/05) 
      IF ( ITS_A_FULLCHEM_SIM() ) THEN 
 
         CALL RESTORE_CONV

         IF ( LPRINTFD ) THEN 
            WRITE(6,*) ' H2O2s before conv adj = ', H2O2s(IFD,JFD,LFD)
            WRITE(6,*) '  SO2s before conv adj = ', SO2s(IFD,JFD,LFD)
         ENDIF 
     
      ENDIF 

      ! Call the S-J Lin convection routine for GEOS-1, GEOS-S, GEOS-3
      CALL NFCLDMX_ADJ( N_TRACERS, TCVV, CMFMC(:,:,2:LLPAR+1), DTRAIN,
     &                  STT _ADJ)

#elif defined( GEOS_3 )

      ! Restore checkpted values of H2O2s and SO2s (dkh, 11/22/05) 
      IF ( ITS_A_FULLCHEM_SIM() ) THEN 

         CALL RESTORE_CONV

         IF ( LPRINTFD ) THEN 
            WRITE(6,*) ' H2O2s before convection = ', H2O2s(IFD,JFD,LFD)
            WRITE(6,*) '  SO2s before convection = ', SO2s(IFD,JFD,LFD)
         ENDIF 
      
      ENDIF 

      !-------------------------
      ! GEOS-3 met fields
      !-------------------------

      ! Call the S-J Lin convection routine for GEOS-1, GEOS-S, GEOS-3
      CALL NFCLDMX_ADJ( N_TRACERS, TCVV, CLDMAS, DTRAIN, STT_ADJ )

#endif

      ! Return to calling program
      END SUBROUTINE DO_CONVECTION_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE DO_GEOS4_CONVECT_ADJ
!
!******************************************************************************
!  Subroutine DO_GEOS4_CONVECT_ADJ is the adjooint of the GEOS4 convection. 
!  Based on DO_GEOS4_CONVECT (swu, bmy, 5/25/05, 10/3/05) with adjoint 
!  updated to GCv8 (ks, mak, dkh, 08/25/09) 
!
!  NOTES:
!  (1 ) Updated to GCv8
! 
!*****************************************************************************
!     
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ
      USE CHECKPT_MOD,       ONLY : CHK_STT_CON
      USE DAO_MOD,           ONLY : HKETA, HKBETA, ZMEU, ZMMU, ZMMD
      USE DIAG_MOD,          ONLY : AD37
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FVDAS_CONVECT_ADJ_MOD, ONLY : FVDAS_CONVECT_ADJ
      USE LOGICAL_MOD,       ONLY : LPRT
      USE PRESSURE_MOD,      ONLY : GET_PEDGE
      USE TIME_MOD,          ONLY : GET_TS_CONV
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS, TCVV
      USE WETSCAV_MOD,       ONLY : COMPUTE_F
      USE WETSCAV_MOD,       ONLY : RESTORE_CONV

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! ND37, LD37 

      ! Local variables 
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, ISOL, J, L, L2, N, NSTEP  
      INTEGER                    :: INDEXSOL(N_TRACERS) 
      INTEGER                    :: CONVDT    
      TYPE (XPLEX)                     :: F(IIPAR,JJPAR,LLPAR,N_TRACERS)
      TYPE (XPLEX)                     :: RPDEL(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: DP(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: P1, P2, TDT   

      !=================================================================
      ! DO_GEOS4_CONVECT_ADJ begins here!
      !=================================================================
      
      ! Convection timestep [s]
      CONVDT = GET_TS_CONV() * 60d0 
       
      ! NSTEP is the # of internal convection timesteps.  According to
      ! notes in the old convection code, 300s works well. (swu, 12/12/03)
      NSTEP  = CONVDT / 300    
      NSTEP  = MAX( NSTEP, 1 ) 

      ! TIMESTEP*2; will be divided by 2 before passing to CONVTRAN 
      TDT    = XPLX( CONVDT ) * 2.0D0 / XPLX( NSTEP )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV_ADJ: a INIT_FV' )

      !=================================================================
      ! Before calling convection, compute the fraction of insoluble
      ! tracer (Finsoluble) lost in updrafts.  Finsoluble = 1-Fsoluble.
      !=================================================================
      ! Need this too for full chemistry. (dkh, 10/01/08) 
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
         CALL RESTORE_CONV 
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ISOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, N_TRACERS

         ! Get fraction of tracer scavenged and the soluble tracer 
         ! index (ISOL). For non-soluble tracers, F=0 and ISOL=0.
         CALL COMPUTE_F( N, F(:,:,:,N), ISOL ) 
         
         ! Store ISOL in an array for later use
         INDEXSOL(N) = ISOL

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ! GEOS-4 convection routines need the insoluble fraction
            F(I,J,L,N) = 1d0 - F(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
       
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV_ADJ: a COMPUTE_F' )

      !=================================================================
      ! Compute pressure thickness arrays DP and RPDEL
      ! These arrays are indexed from atm top --> surface
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, L2, P1, P2 )
      DO L = 1, LLPAR

         ! L2 runs from the atm top down to the surface
         L2 = LLPAR - L + 1

         ! Loop over surface grid boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
            ! Pressure at bottom and top edges of grid box [hPa]
            P1 = GET_PEDGE(I,J,L)
            P2 = GET_PEDGE(I,J,L+1)

            ! DP = Pressure difference between top & bottom edges [Pa]
            DP(I,J,L2) = ( P1 - P2 ) * 100.0d0

            ! RPDEL = reciprocal of DP [1/hPa]
            RPDEL(I,J,L2) = 100.0d0 / DP(I,J,L2) 
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV_ADJ: a DP, RPDEL' )
   
      !=================================================================
      ! Flip arrays in the vertical and call FVDAS_CONVECT
      !=================================================================

      ! Call the fvDAS convection routines (originally from NCAR!)
      CALL FVDAS_CONVECT_ADJ( TDT,      
     &                    N_TRACERS, 
     &                    CHK_STT_CON(:,:,LLPAR:1:-1,:),    
     &                    RPDEL,         
     &                    HKETA (:,:,LLPAR:1:-1  ),
     &                    HKBETA(:,:,LLPAR:1:-1  ), 
     &                    ZMMU  (:,:,LLPAR:1:-1  ),    
     &                    ZMMD  (:,:,LLPAR:1:-1  ),  
     &                    ZMEU  (:,:,LLPAR:1:-1  ),  
     &                    DP,     
     &                    NSTEP,    
     &                    F     (:,:,LLPAR:1:-1,:),         
     &                    TCVV,   
     &                    INDEXSOL,STT_ADJ(:,:,LLPAR:1:-1,:) )


      !### Debug! 
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV_ADJ: a FVDAS_CONVECT')

      ! Return to calling program
      END SUBROUTINE DO_GEOS4_CONVECT_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE NFCLDMX_ADJ( NC, TCVV, CLDMAS, DTRN, Q )
!
!******************************************************************************
! Subroutine ADJ_NFDCLDMX is based on the original NFCLDMX code, where the
!  loop over the tracers has been extracted, sent to TAMC, and reinserted.
! (dkh, 02/22/05)
!
!  Arguments as input:
!  ==========================================================================
!  (1 ) NC     : TOTAL number of tracers (soluble + insoluble)  [unitless]
!  (2 ) TCVV   : MW air (g/mol) / MW of tracer (g/mol)          [unitless]
!  (3 ) CLDMAS : Cloud mass flux (at upper edges of each level) [kg/m2/s]
!  (4 ) DTRN   : Detrainment mass flux                          [kg/m2/s]
!
!  Arguments as Input/Output:
!  ============================================================================
!  (5 )  Q     : Tracer concentration                           [v/v]
!
!  NOTES:
!  (1 ) See orignial NFCLDMX for references,descriptions and notes.
!  (2 ) TAMC code and varialbes are lowercase.
!  (3 ) Use COMPUTE_ADJ_F from WETSCAV_ADJ_MOD rather than COMPUTE_F from 
!        WETSCAV_MOD	
!  (4 ) Get rid of excess array element copying (dkh, 03/01/05)
!  (5 ) Leave out ( Q  + DELQ > 0 ) condition, as we don't need to force 
!        the adjoints to be positive definite.
!  (6 ) Add support for carbon, dust, ss.  (dkh, 03/05/05)
!  (7 ) Now include CMN_ADJ to allow for printout. (dkh, 03/14/05)
!  (8 ) Rebuild adjoing so that can loop easily over I,J (dkh, 03/22/05)
!  (9 ) Now reference WETSCAV_MOD instead of WETSCAV_ADJ_MOD. (dkh, 10/24/05)  
!  (10) Updated to GCv8 (dkh, 08/25/09) 
!  (11) BUG FIX: Now correctly reset adjoints for GEOS_5 (dkh, 04/21/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,           ONLY : AD  !,   CLDMAS, DTRN=>DTRAIN
      USE DIAG_MOD,          ONLY : AD37, AD38,   CONVFLUP
      USE GRID_MOD,          ONLY : GET_AREA_M2
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN
      USE OCEAN_MERCURY_MOD, ONLY : ADD_Hg2_WD
      USE PRESSURE_MOD,      ONLY : GET_BP, GET_PEDGE
      USE TIME_MOD,          ONLY : GET_TS_CONV
      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM
      USE TRACERID_MOD,      ONLY : IS_Hg2
      USE WETSCAV_MOD,       ONLY : COMPUTE_F

! DISABLE this for now.  It needs to be further validated. (dkh, 10/12/08) 
!      !>>>
!      ! Now include adjoint of F (dkh, 10/03/08) 
!      USE WETSCAV_MOD,       ONLY : QC_SO2
!      USE WETSCAV_MOD,       ONLY : ADJ_COMPUTE_F
!      USE WETSCAV_MOD,       ONLY : ADJ_F
!      USE WETSCAV_MOD,       ONLY : RESTORE_CONV
!      USE TRACERID_MOD,      ONLY : IDTSO2
!      !<<<
      USE LOGICAL_ADJ_MOD,   ONLY : LPRINTFD
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, LFD, NFD

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches & arrays

      ! Arguments
      INTEGER, INTENT(IN)    :: NC
      TYPE (XPLEX),  INTENT(IN)    :: CLDMAS(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX),  INTENT(IN)    :: DTRN(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX),  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NC)
      TYPE (XPLEX),  INTENT(IN)    :: TCVV(NC)


      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      LOGICAL, SAVE          :: IS_Hg = .TRUE.
      INTEGER                :: I, J, K, KTOP, L, N, NDT
      INTEGER                :: IC, ISTEP, JUMP, JS, JN, NS
      INTEGER                :: IMR, JNP, NLAY
      TYPE (XPLEX),  SAVE          :: DSIG(LLPAR)
      TYPE (XPLEX)                 :: SDT, CMOUT, ENTRN, DQ, AREA_M2
      TYPE (XPLEX)                 :: T0, T1, T2, T3, T4, TSUM, DELQ
      TYPE (XPLEX)                 :: DTCSUM(IIPAR,JJPAR,LLPAR,NC)

      ! F is the fraction of tracer lost to wet scavenging in updrafts
      TYPE (XPLEX)                 :: F(IIPAR,JJPAR,LLPAR,NC)

      ! Local Work arrays (Comment out those that are superfluous for adj)
      TYPE (XPLEX)                 :: BMASS(IIPAR,JJPAR,LLPAR)
      !TYPE (XPLEX)                 :: QB(IIPAR,JJPAR)
      !TYPE (XPLEX)                 :: MB(IIPAR,JJPAR)
      !TYPE (XPLEX)                 :: QC(IIPAR,JJPAR)

      ! TINY = a very small number
      TYPE (XPLEX), PARAMETER      :: TINY = xplex(1d-14,0d0)

      ! ISOL is an index for the diagnostic arrays
      INTEGER                :: ISOL

      ! QC_PRES and QC_SCAV are the amounts of tracer
      ! preserved against and lost to wet scavenging
      ! Not needed for adjoint
      !TYPE (XPLEX)                 :: QC_PRES, QC_SCAV

      ! DNS is the TYPE (XPLEX) value for NS
      TYPE (XPLEX)                 :: DNS

      ! Amt of Hg2 scavenged out of the column (sas, bmy, 1/19/05)
      TYPE (XPLEX)                 :: WET_Hg2

      !>>>
      ! Now include adjoint of F (dkh, 10/03/08) 
      TYPE (XPLEX)                 :: F_SO2(IIPAR,JJPAR,LLPAR)
      !<<<

C==============================================
C define arguments (comment out those already defined)
C==============================================
      TYPE (XPLEX) adq_in(llpar)
      TYPE (XPLEX) adq_out(llpar)
      TYPE (XPLEX) vbmass(llpar)
      TYPE (XPLEX) vcldmas(llpar)
      !TYPE (XPLEX) dsig(llpar)
      TYPE (XPLEX) vdtrn(llpar)
      TYPE (XPLEX) vf(llpar)
      !integer ktop
      !integer ns
      !TYPE (XPLEX) sdt

C==============================================
C define local variables (comment out those already defined)
C==============================================
      TYPE (XPLEX) addelq
      TYPE (XPLEX) adq(llpar)
      TYPE (XPLEX) adqb
      TYPE (XPLEX) adqc
      TYPE (XPLEX) adqc_pres
      TYPE (XPLEX) adt1
      TYPE (XPLEX) adt2
      TYPE (XPLEX) adt3
      TYPE (XPLEX) adt4
      TYPE (XPLEX) adtsum
      !TYPE (XPLEX) cmout
      !TYPE (XPLEX) entrn
      integer ip1
      !integer istep
      !integer k
      TYPE (XPLEX) mb

      !=================================================================
      ! ADJ_NFCLDMX begins here!
      !=================================================================



      ! First-time initialization
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'N F C L D M X  -- by S-J Lin'
         WRITE( 6, '(a)' ) 'Modified for GEOS-CHEM by Bob Yantosca'
         WRITE( 6, '(a)' ) 'Last Modification Date: 1/27/04'
         WRITE( 6, '(a)' ) 'Adjoint constucted with TAMC: dkh, 03/01/05'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

#if   !defined( GEOS_5 ) 
         ! NOTE: We don't need to do this for GEOS-5 (bmy, 6/27/07)
         ! DSIG is the sigma-level thickness (NOTE: this assumes that
         ! we are using a pure-sigma grid.  Use new routine for fvDAS.)
         DO L = 1, LLPAR
            DSIG(L) = GET_BP(L) - GET_BP(L+1)
         ENDDO
#endif


         ! Reset first time flag
         FIRST = .FALSE.
      ENDIF

      ! Define dimensions
      IMR  = IIPAR
      JNP  = JJPAR
      NLAY = LLPAR

      ! Convection timestep [s]
      NDT  = GET_TS_CONV() * 60d0

      !=================================================================
      ! Define active convective region, from J = JS(outh) to
      ! J = JN(orth), and to level K = KTOP.
      !
      ! Polar regions are too cold to have moist convection.
      ! (Dry convection should be done elsewhere.)
      !
      ! We initialize the ND14 diagnostic each time we start a new
      ! time step loop.  Only initialize DTCSUM array if the ND14
      ! diagnostic is turned on.  This saves a quite a bit of time.
      ! (bmy, 12/15/99)
      !=================================================================
      IF ( ND14 > 0 ) DTCSUM = 0d0

      KTOP = NLAY - 1
      JUMP = (JNP-1) / 20
      JS   = 1 + JUMP
      JN   = JNP - JS + 1

      !=================================================================
      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.
      !=================================================================
      NS  = NDT / 300
      NS  = MAX(NS,1)
      SDT = XPLX(NDT) / XPLX(NS)
      DNS = XPLX( NS )

!=============================================================================
!  BMASS has units of kg/m^2 and is equivalent to AD(I,J,L) / AREA_M2
!
!   Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!  -------------+---------+--------+-------+----+--------  =  -----
!               | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!  This is done to keep BMASS in the same units as CLDMAS * SDT
!
!  We can parallelize over levels here.  The only quantities that need to
!  be held local are the loop counters (I, IC, J, JREF, K). (bmy, 5/2/00)
!
!  Now use routine GET_AREA_M2 from "grid_mod.f" to get surface area of
!  grid boxes in m2. (bmy, 2/4/03)
!=============================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_M2, K )
!$OMP+SCHEDULE( DYNAMIC )
      DO K = 1, NLAY
         DO J = 1, JJPAR
            AREA_M2 = GET_AREA_M2( J )
            DO I = 1, IMR
               BMASS(I,J,K) = AD(I,J,K) / AREA_M2
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! (1)  T r a c e r   L o o p
      !
      ! We now parallelize over tracers, since tracers are independent
      ! of each other.  The parallel loop only takes effect if you
      ! compile with the f90 "-mp" switch.  Otherwise the compiler will
      ! interpret the parallel-processing directives as comments, and
      ! the loop will execute on a single thread.
      !
      ! The following types of quantities must be held local for
      ! parallelization:
      ! (1) Loop counters ( I, IC, ISTEP, J, K )
      ! (2) Scalars that are assigned values inside the tracer loop:
      !     ( CMOUT, DELQ, ENTRN, ISOL, QC_PRES, etc. )
      ! (3) Arrays independent of tracer ( F, MB, QB, QC )
      !=================================================================
      !>>>
      ! Now include adjoint of F (dkh, 10/03/08) 
      ! OLD:
      !DO IC = 1, NC
      !   CALL COMPUTE_ADJ_F( IC, F(:,:,:,IC), ISOL )
      !ENDDO
      ! NEW:
      DO IC = 1, NC
         CALL COMPUTE_F( IC, F(:,:,:,IC), ISOL )
      ENDDO
!
! DISABLE this for now.  It needs to be further validated. (dkh, 10/12/08) 
!      F_SO2(:,:,:) = F(:,:,:,IDTSO2)
!
!      !<<<
 
      IF ( LPRINTFD ) THEN
         WRITE(165,*) ' Convection variables ',
     &                ' AD(FD) =      ', AD(IFD,JFD,LFD),
     &                ' CLDMAS =      ', CLDMAS(IFD,JFD,LFD),
     &                ' DTRN   =      ', DTRN(IFD,JFD,LFD),
     &                ' GET_BP =      ', GET_BP(LFD),
     &                ' GET_AREA_M2 = ', GET_AREA_M2(JFD),
     &                ' F =           ', F(IFD,JFD,LFD,NFD)
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( CMOUT, DELQ, ENTRN, I, IC, ISOL, ISTEP, J, K           )
!$OMP+PRIVATE( MB, T0, T1, T2, T3, T4, TSUM                           )
!$OMP+PRIVATE( WET_Hg2                                                )
!$OMP+PRIVATE( addelq, adqb, adqc, adqc_pres, adt1, adt2, adt3, adt4  )
!$OMP+PRIVATE( adtsum, ip1, adq_out, vdtrn, vbmass, vcldmas, vf       )
!$OMP+PRIVATE( adq_in, adq                                            )
!$OMP+SCHEDULE( DYNAMIC )
      DO IC = 1, NC
      DO J  = JS, JN
      DO I  = 1, IMR

      adq_out(:)  = Q     (I,J,:,IC)
      vdtrn  (:)  = DTRN  (I,J,:)
      vbmass (:)  = BMASS (I,J,:)
      vcldmas(:)  = CLDMAS(I,J,:)
      vf     (:)  = F     (I,J,:,IC)  

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      addelq = 0.
      do ip1 = 1, llpar
        adq(ip1) = 0.
      end do
      adqb = 0.
      adqc = 0.
      adqc_pres = 0.
      adt1 = 0.
      adt2 = 0.
      adt3 = 0.
      adt4 = 0.
      adtsum = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      adq(:) = adq(:)+adq_out(:)
      adq_in(:) = 0d0

      adq_out(:) = 0.
 
      ! IF ( LPRINTFD .and. i == IFD .and. j == JFD .and.
      !&    ic == STT2ADJ(NFD) ) THEN
      !     print*, 'adq = ', adq
      ! ENDIF

      do istep = ns, 1, -1
        do k = ktop, 3, -1
          if (vcldmas(k-1) .gt. tiny) then
            cmout = vcldmas(k)+vdtrn(k)
            entrn = cmout-vcldmas(k-1)
            addelq = addelq+adq(k)

            ! note: need to implement CONVECTION_FLOW_CHK
            ! fwd code:
            !IF ( Q(I,J,K,IC) + DELQ < 0.0d0 ) THEN
            !    DELQ = -Q(I,J,K,IC)
            !ENDIF
            ! adj code:
            !IF ( CONVECTION_FLOW_CHK(I,J,ISTEP,1) ) THEN 
            !    ADQ(K) = -ADDELQ
            !ENDIF 

            adtsum = adtsum+addelq*(sdt/vbmass(k))
            addelq = 0.
            adt1 = adt1+adtsum
            adt2 = adt2+adtsum
            adt3 = adt3+adtsum
            adt4 = adt4+adtsum
            adtsum = 0.
            adq(k) = adq(k)-adt4*vcldmas(k-1)
            adt4 = 0.
            adq(k+1) = adq(k+1)+adt3*vcldmas(k)
            adt3 = 0.
            adqc = adqc-adt2*vcldmas(k)
            adt2 = 0.
            adqc_pres = adqc_pres+adt1*vcldmas(k-1)
            adt1 = 0.
            if (entrn .ge. 0) then
              adq(k) = adq(k)+adqc*(entrn/cmout)
              adqc_pres = adqc_pres+adqc*(vcldmas(k-1)/cmout)
              adqc = 0.
            endif
            adqc = adqc+adqc_pres*(1.d0-vf(k))
! DISABLE this for now.  It needs to be further validated. (dkh, 10/12/08) 
!            !>>>
!            ! Now include adjoint of F(SO2) (dkh, 10/03/08) 
!            ! fwd code: 
!            !QC_PRES  = QC(I,J) * ( 1d0 - F(I,J,K,IC) )
!            ! adj code:
!            IF ( IC == IDTSO2 ) THEN
!               ADJ_F(I,J,K) = ADJ_F(I,J,K)
!     &                      - QC_SO2(I,J,K,ISTEP) * ADQC_PRES
!            ENDIF
!            !<<<
            adqc_pres = 0.
          else

#if   defined( GEOS_5 ) 
            IF ( CLDMAS(I,J,K) > TINY ) THEN

               ! fwd code:
               !Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ
               ! adj code:
               ADDELQ = ADQ(K) 

               ! note: need to implement CONVECTION_FLOW_CHK
               ! fwd code:
               !IF ( Q(I,J,K,IC) + DELQ < 0.0d0 ) THEN
               !    DELQ = -Q(I,J,K,IC)
               !ENDIF
               ! adj code:
               !IF ( CONVECTION_FLOW_CHK(I,J,ISTEP,2) ) THEN 
               !    ADQ(K) = -ADDELQ
               !ENDIF 

               ! fwd code:
               !DELQ = ( SDT / BMASS(I,J,K) ) * (T2 + T3)
               ! adj code:
               ADT2 = ( SDT / VBMASS(K) ) * ADDELQ
               ADT3 = ( SDT / VBMASS(K) ) * ADDELQ
               ! BUG FIX: make sure to reset ADDELQ (dkh, 04/21/10)
               ADDELQ = 0d0

               ! fwd code:
               !T3   =  CLDMAS(I,J,K  ) * Q (I,J,K+1,IC)
               ! adj code:
               ADQ(K+1) = ADQ(K+1) + VCLDMAS(K) * ADT3
              ! BUG FIX: make sure to reset ADT3 (dkh, 04/21/10)
               ADT3 = 0d0

               ! fwd code:
               !T2   = -CLDMAS(I,J,K  ) * QC(I,J)
               ! adj code:
               ADQC = ADQC - VCLDMAS(K) * ADT2
              ! BUG FIX: make sure to reset ADT2 (dkh, 04/21/10)
               ADT2 = 0d0


            ENDIF
#endif

            adq(k) = adq(k)+adqc
            adqc = 0.
          endif
        end do

        !   IF ( LPRINTFD .and. i == IFD .and. j == JFD .and.
        !&    ic == STT2ADJ(NFD) ) THEN
        !     print*, 'adq = ', adq
        !   ENDIF

        if (vcldmas(2) .gt. tiny) then
          mb = vbmass(1)+vbmass(2)
          adqc = adqc+adq(1)
          adq(1) = 0.
          adqc = adqc+adq(2)
          adq(2) = 0.
          adq(3) = adq(3)+adqc*(vcldmas(2)*sdt/(mb+vcldmas(2)*sdt))
          adqb = adqb+adqc*(mb/(mb+vcldmas(2)*sdt))
          adqc = 0. 
#if   defined ( GEOS_5 )
          ! for GEOS-5 (dkh, 08/25/09) 
          adq(2) = adq(2)+adqb*(( GET_PEDGE(I,J,2) - GET_PEDGE(I,J,3) )
     &                      /( GET_PEDGE(I,J,1) - GET_PEDGE(I,J,3) ) )
          adq(1) = adq(1)+adqb*(( GET_PEDGE(I,J,1) - GET_PEDGE(I,J,2) )
     &                      /( GET_PEDGE(I,J,1) - GET_PEDGE(I,J,3) ) )
#else     
          ! for GEOS-3
          adq(2) = adq(2)+adqb*(dsig(2)/(dsig(1)+dsig(2)))
          adq(1) = adq(1)+adqb*(dsig(1)/(dsig(1)+dsig(2)))
#endif 
          adqb = 0.
        else
          adq(3) = adq(3)+adqc
          adqc = 0.
        endif
      end do
      adq_in(:) = adq_in(:)+adq(:)
      adq(:) = 0.

      Q(I,J,:,IC) = adq_in(:)


      ENDDO           !I
      ENDDO           !J
      ENDDO           !IC
!$OMP END PARALLEL DO

! DISABLE this for now.  It needs to be further validated. (dkh, 10/12/08) 
!      !>>>  
!      ! Now include adjoint of F(SO2) (dkh, 10/03/08) 
!      ! Restore H2O2s and SO2s to their pre-convection values
!      CALL RESTORE_CONV
!            
!      ! fwd code: 
!      !CALL COMPUTE_F( IC, F(:,:,:,IC), ISOL )
!      ! adj code: 
!      CALL ADJ_COMPUTE_F( F_SO2(:,:,:) )
!      !<<<

      ! Return to calling program
      END SUBROUTINE NFCLDMX_ADJ
!------------------------------------------------------------------------------

      END MODULE CONVECTION_ADJ_MOD
