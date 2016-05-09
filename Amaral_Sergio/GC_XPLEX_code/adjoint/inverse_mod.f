!$Id: inverse_mod.f,v 1.20 2012/03/04 19:34:15 daven Exp $
      MODULE INVERSE_MOD
!
!*****************************************************************************
!  Module INVERSE_MOD contains all the subroutines that used to be in
!  inverse.f. While having these routines in the top most program file worked
!  on SGI, it didn't work on Linux, so had to move all to a module.
!  (dkh, 02/05)!
!  Module Variables:
!  ============================================================================
!  (1 ) COST_FUNC   (TYPE (XPLEX))      : Value of cost function
!  (2 ) N_CALC      (INTEGER)     : Optimization iteration number
!  (3 ) N_CALC_STOP (INTEGER)     : Maximum optimization iteration number
!  (4 ) F           (DOUBLE)      : For optimization routine
!  (5 ) X           (DOUBLE, ALLOC): Vector of active varialbes
!  (6 ) GRADNT      (DOUBLE, ALLOC): Vector of adjoint gradients
!  (7 ) XP          (DOUBLE, ALLOC): Vector of active strat prod varialbes
!  (8 ) GRADNT_P    (DOUBLE, ALLOC): Vector of strat prod adjoint gradients
!  (9 ) XL          (DOUBLE, ALLOC): Vector of active strat loss varialbes
!  (10) GRADNT_L    (DOUBLE, ALLOC): Vector of adjoint strat loss gradients
!
!  Module Routines
!  ============================================================================
!  (1 ) SET_SF              : Initializes ICS_SF and EMS_SF
!  (2 ) SET_LOG_SF          : Initializes ICS_SF and EMS_SF for log scaling
!  (3 ) GET_X_FROM_SF       : Turns SF array into a vector X for optimization
!  (4 ) GET_SF_FROM X       : Turns vector X into array SF after optimization
!  (5 ) GET_GRADNT_FROM_ADJ : Turns ADJ_STT array into vector GRADNT for opt.
!  (6 ) MAKE_GDT_FILE       : Save GRADNT values at iteration N_CALC to adjtmp/*gdt*
!  (7 ) READ_GDT_FILE       : Reads saved GRADNT values from previous iterations
!  (8 ) MAKE_SF_FILE        : Saves SF at iteration N_CALC to adjtmp/*sf*
!  (9 ) READ_SF_FILE        : Reads saved SF from previous iterations
!  (10) EXPAND_NAME         : Adds iteration number to file names
!  (11) DISPLAY_STUFF       : Echo various things at each iteration
!  (12) SET_SF_FORFD        : Set the scaling factors for finite difference test.
!  (13) MAKE_CFN_FILE       : Save cost function to cnf.* file
!  (14) READ_CFN_FILE       : Read cost function from cnf.* file
!  (15) SET_OPT_RANGE       : Set range of parameters to optimize
!  (16) CALC_NOPT           : Set range of parameters to optimize
!  (17) ITER_CONDITION      : Write out iteration diagnostics to gctm.iteration
!  (18) MAYBE_DO_GEOS_CHEM_ADJ: For FDTEST determine if need to call adjoint 
!  (19) INIT_INVERSE        : Initialize allocatable arrays
!  (20) CLEANUP_INVERSE     : Deallocatte arrays
!
!  Modules referenced by "inverse_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f       : Module containing routines for binary pch file I/O
!  (2 ) charpak_mod.f     : Module containing string handling routines
!  (3 ) error_mod.f       : Module containing NaN and other error check routines
!  (4 ) file_mod.f        : Module containing file unit numbers & error checks
!  (5 ) grid_mod.f        : Module containing horizontal grid information
!  (6 ) restart_mod.f     : Module containing CHECK_DIMENSIONS
!  (7 ) time_mod.f        : Module containing routines to compute time & date
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added subroutine INIT_REGIONAL_ICS (dkh, 02/12/05)
!  (3 ) Now use IDADJxxx (03/03/05)
!  (4 ) Don't zero the adjoints of NO3, NIT, and NH4 
!  (5 ) Now save EMS_ICS from reference run to EMS_ICS_orig, a mod variable
!        Also update MAKE_GDT and MAKE_ICS to handle all emissions.
!        (dkh, 03/29/05)
!  (6 ) Remove all duplicate declarations of N_CALC and N_CALC_STOP. Now this
!        is always treated as a module variable. (dkh, 02/15/06)  
!  (7 ) Update MAKE_ICS_FILE to support writing initial NOx emisions. (dkh, 08/27/06)  
!  (8 ) Bug fix: change N to 1 in TRACER(I,J,1) while writing scaled 
!       emissions. (dkh, 10/26/06)  
!  (9 ) BUG FIX: make ADJ_STT_FD allocatable. (dkh, 03/21/07) 
!  (10) Update to support LOG_OPT pre-processor option. 
!  (11) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!*****************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
 
#     include "define_adj.h"    ! obs operators
      
      !====================================================================
      ! MODULE VARIABLES  ( those that used to be program variables )
      !====================================================================
      !TYPE (XPLEX), ALLOCATABLE  :: EMS_ICS_orig(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: X(:)
      TYPE (XPLEX), ALLOCATABLE  :: GRADNT(:)

      !For strat prod & loss SF (hml, 08/11/14)
      TYPE (XPLEX), ALLOCATABLE  :: XP(:)
      TYPE (XPLEX), ALLOCATABLE  :: GRADNT_P(:)
      TYPE (XPLEX), ALLOCATABLE  :: XL(:)
      TYPE (XPLEX), ALLOCATABLE  :: GRADNT_L(:)

      !====================================================================
      ! MODULE ROUTINES 
      !====================================================================
      CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE SET_SF 
!
!*****************************************************************************
!  Subroutine SET_SF sets the intial conditions used for a GEOS_CHEM run
!  (dkh, 9/16/04).  
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Switch to using IDADJxxx  (dkh, 03/03/05)
!  (4 ) Rename to SET_SF, replace CMN_ADJ with adjoint_array_mod 
!        (dkh, ks, mak, cs  06/07/09) 
!  (5 ) Now get first guesses from input.gcadj file (mak, 9/23/09)
!  (6 ) Now use ICS_SF_DEFAULT and ICS_SF_DEFAULT instad of ICS_SF_tmp
!        and EMS_SF_tmp. (dkh, 02/09/11) 
!  (7 ) Now support strat fluxes LADJ_STRAT and add flags to avoid accessing 
!        unallocated arrays (hml, dkh, 02/20/12, adj32_025) 
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : ICS_SF, ICS_SF0
      USE ADJ_ARRAYS_MOD,     ONLY : EMS_SF, EMS_SF0
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC_STOP
      USE ADJ_ARRAYS_MOD,     ONLY : IDADJ_ENH3_an
      USE ADJ_ARRAYS_MOD,     ONLY : ADCOEMS
      USE ADJ_ARRAYS_MOD,     ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD,     ONLY : NFD, MFD, EMSFD
      USE ADJ_ARRAYS_MOD,     ONLY : ICSFD
      !USE ADJ_ARRAYS_MOD,     ONLY : ICS_SF_tmp, EMS_SF_tmp
      USE ADJ_ARRAYS_MOD,     ONLY : ICS_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,     ONLY : EMS_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,     ONLY : MMSCL
      USE ADJ_ARRAYS_MOD,     ONLY : PROD_SF,   PROD_SF0
      USE ADJ_ARRAYS_MOD,     ONLY : LOSS_SF,   LOSS_SF0
      USE ADJ_ARRAYS_MOD,     ONLY : PROD_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,     ONLY : LOSS_SF_DEFAULT
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,    ONLY : LICS,   LADJ_EMS
      USE LOGICAL_ADJ_MOD,    ONLY : LADJ_STRAT
      USE TRACERID_MOD,       ONLY : IDTOX,  IDTNOX

#     include "CMN_SIZE"           ! Size params
#     include "define_adj.h"       ! obs operators

      ! local variables 
      INTEGER                    :: I
      INTEGER                    :: J
      INTEGER                    :: L
      INTEGER                    :: M
 
      !=================================================================
      ! SET_SF begins here!
      !=================================================================


      ! Set to defaults or user defined values 
      IF ( N_CALC_STOP .EQ. 0) THEN

         ! Set default scaling factors to 1d0 everywhere for reference run
         ! (perfect model generating pseudo observations)
         ICS_SF(:,:,:,:)   = 1.d0

         IF ( LADJ_EMS )  EMS_SF(:,:,:,:) = 1.d0
            
         IF ( LADJ_STRAT ) THEN
            PROD_SF(:,:,:,:) = 1.d0
            LOSS_SF(:,:,:,:) = 1.d0
         ENDIF

      ELSE
 
         ! Now define defaults for all in input.gcadj (dkh, 02/09/11) 
         !EMS_SF(:,:,:,:) = 1.d0
         !ICS_SF(:,:,:,:) = 1.d0
         !! otherwise, use values from input.gcadj file for ICSFD and EMSFD
         !EMS_SF(:,:,:,EMSFD) = EMS_SF_tmp
         !ICS_SF(:,:,:,ICSFD) = ICS_SF_tmp
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L ) 
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ICS_SF (I,J,L,:)  = ICS_SF_DEFAULT (:)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         IF ( LADJ_EMS ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M ) 
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               EMS_SF(I,J,M,:)   = EMS_SF_DEFAULT(:)
 
               IF ( LADJ_STRAT ) THEN
                  PROD_SF(I,J,M,:)  = PROD_SF_DEFAULT(:)
                  LOSS_SF(I,J,M,:)  = LOSS_SF_DEFAULT(:)
               ENDIF 

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF

      ! the following options for PSEUDO_OBS should not become obsolete
      ! We don't even have to remember to change initial SF between pseudo obs
      ! and perturbed run, since the adjustment will be made automaticlally above
      ! the following #if statement can be removed.
      ! one thing we can add is NFD selection to EMS_SF, so that we don't have to 
      ! perturb all emissions, but only one. not sure if this would be good or 
      ! would just complicate things...
      ! (mak, 9/23/09)

#if defined ( PSEUDO_OBS )

      ! Make the initial guess for iteration N_CALC == 1 
      ! BUG FIX:  make sure this happens every time the optimization 
      !  cycles through N_CALC = 1 as well. (mak, dkh, 09/08/09) 
      !IF ( N_CALC == 1 ) THEN 
!      IF ( N_CALC == 1 
!     &     .or. ( N_CALC == 0 .and. N_CALC_STOP > 1 ) ) THEN 
      IF ( N_CALC == 1 .or.  
     &           ( N_CALC == 0 .and. N_CALC_STOP > 0 ) ) THEN 
 
         ! For control parameters = initial conditions 
         IF ( LICS ) THEN

            ! Now enforce defaults for all set in input.gcadj (dkh, 02/09/11) 
            !! BUG FIX: enforce defualt scaling factors before using SF_tmp 
            !! (dkh, 07/30/10) 
            !ICS_SF(:,:,:,:) = 1.d0
            !
            !print*, 'set ICS_SF to', ICS_SF_tmp
            !! Start with an initial guess for ICS_SF that is wrong
            !! Let's set the default to perturb everything to avoid
            !! hardwiring (mak, 6/18/09)
            !! now this is done via input.gcadj file (mak, 9/23/09)
            !ICS_SF(:,:,:,ICSFD) = ICS_SF_tmp !0.5d0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L ) 
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ICS_SF(I,J,L,:) = ICS_SF_DEFAULT(:)

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

  
         ELSEIF ( LADJ_EMS ) THEN 
         
            ! Now enforce defaults for all set in input.gcadj (dkh, 02/09/11) 
            !!! BUG FIX: enforce defualt scaling factors before using SF_tmp 
            !! (dkh, 07/30/10) 
            !EMS_SF(:,:,:,:) = 1.d0
            !
            !! Start with an initial guess for EMS_SF that is wrong
            !EMS_SF(:,:,1,EMSFD) = EMS_SF_tmp 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M ) 
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               EMS_SF(I,J,M,:) = EMS_SF_DEFAULT(:)
               
               IF ( LADJ_STRAT ) THEN
                  PROD_SF(I,J,M,:)  = PROD_SF_DEFAULT(:)
                  LOSS_SF(I,J,M,:)  = LOSS_SF_DEFAULT(:)
               ENDIF

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF 
     
      ENDIF 

#endif 

      ! Save a copy of the initial guess of the scaling factors 
      ! for use later in calculating the a priori penalty term 
      ICS_SF0 (:,:,:,:) = ICS_SF (:,:,:,:)

      ! Add flags (hml, 02/23/12)
      IF ( LADJ_EMS )  EMS_SF0 (:,:,:,:) = EMS_SF (:,:,:,:)

      IF ( LADJ_STRAT ) THEN
         PROD_SF0(:,:,:,:) = PROD_SF(:,:,:,:)
         LOSS_SF0(:,:,:,:) = LOSS_SF(:,:,:,:)
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_SF
!-----------------------------------------------------------------------------

      SUBROUTINE SET_LOG_SF 
!
!*****************************************************************************
!  Subroutine SET_LOG_SF sets the intial conditions used for a GEOS_CHEM run
!  (dkh, 9/16/04).  
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Switch to using IDADJxxx  (dkh, 03/03/05)
!  (4 ) Rename to SET_LOG_SF, replace CMN_ADJ with adjoint_array_mod 
!        (dkh, ks, mak, cs  06/07/09) 
!  (5 ) Now use ICS_SF_DEFAULT and ICS_SF_DEFAULT instad of ICS_SF_tmp
!        and EMS_SF_tmp. (dkh, 02/09/11) 
!  (6 ) Add flags to avoid accessing unallocated arrays 
!        (hml, dkh, 02/27/12, adj32_025) 
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : ICS_SF, ICS_SF0
      USE ADJ_ARRAYS_MOD,     ONLY : EMS_SF, EMS_SF0
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : IDADJ_ENH3_an
      USE ADJ_ARRAYS_MOD,     ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD,     ONLY : NFD, MFD, EMSFD
      USE ADJ_ARRAYS_MOD,     ONLY : ICSFD
      USE ADJ_ARRAYS_MOD,     ONLY : IDADJ_ENH3_an
      !USE ADJ_ARRAYS_MOD,     ONLY : ICS_SF_tmp, EMS_SF_tmp
      USE ADJ_ARRAYS_MOD,     ONLY : ICS_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,     ONLY : EMS_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC_STOP
      USE ADJ_ARRAYS_MOD,     ONLY : MMSCL
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,    ONLY : LICS,   LADJ_EMS
      USE LOGICAL_ADJ_MOD,    ONLY : LADJ_STRAT 
      USE TRACERID_MOD,       ONLY : IDTOX

#     include "CMN_SIZE"           ! Size params

      ! Internal varaibles
      INTEGER                    :: I
      INTEGER                    :: J
      INTEGER                    :: L
      INTEGER                    :: M


      !=================================================================
      ! SET_LOG_SF begins here!
      !=================================================================
  
      IF ( LADJ_STRAT ) THEN 
         CALL ERROR_STOP(' LADJ_STRAT not yet implemented for LOG_OPT',
     &                   ' subroutine SET_LOG_SF, inverse_mod.f ' )
      ENDIF 
 

      ! Set to defaults or user defined values 
      ! Add flags (hml, 02/23/12)
      IF ( N_CALC_STOP .EQ. 0) THEN
         ! Set default scaling factors to 0d0 everywhere for reference run
         ! (perfect model generating pseudo observations)
         ICS_SF(:,:,:,:) = 0.d0
         IF ( LADJ_EMS) EMS_SF(:,:,:,:) = 0.d0
      ELSE
 
   
         ! Now define defaults for all in input.gcadj (dkh, 02/09/11) 
         !EMS_SF(:,:,:,:) = 0.d0
         !ICS_SF(:,:,:,:) = 0.d0
         !! otherwise, use values from input.gcadj file for ICSFD and EMSFD
         !EMS_SF(:,:,:,EMSFD) = LOG(EMS_SF_tmp)
         !ICS_SF(:,:,:,ICSFD) = LOG(ICS_SF_tmp)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L ) 
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ICS_SF(I,J,L,:) = LOG(ICS_SF_DEFAULT(:))

         ENDDO
         ENDDO
         ENDDO

         ! Add flags (hml, 02/23/12)
         IF ( LADJ_EMS ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M ) 
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               EMS_SF(I,J,M,:) = LOG(EMS_SF_DEFAULT(:))

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF


#if defined ( PSEUDO_OBS )

      ! BUG FIX:  make sure this happens every time the optimization 
      !  cycles through N_CALC = 1 as well. (mak, dkh, 09/08/09) 
      !IF ( N_CALC == 1 ) THEN 
!      IF ( N_CALC == 1 
!     &     .or. ( N_CALC == 0 .and. N_CALC_STOP > 1 ) ) THEN 
      IF ( N_CALC == 1 .or.  
     &           ( N_CALC == 0 .and. N_CALC_STOP > 0 ) ) THEN

 
         ! For control parameters = initial conditions 
         IF ( LICS ) THEN

            ! Now define defaults for all in input.gcadj (dkh, 02/09/11) 
            !! BUG FIX: enforce defualt scaling factors before using SF_tmp 
            !! (dkh, 07/30/10) 
            !ICS_SF(:,:,:,:) = 0.d0
            !
            !! Start with an initial guess for ICS_SF that is wrong
            !ICS_SF(:,:,:,ICSFD) = LOG(ICS_SF_tmp)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L ) 
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ICS_SF(I,J,L,:) = LOG(ICS_SF_DEFAULT(:))

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

  
         ENDIF 
         IF ( LADJ_EMS ) THEN 

            ! Now define defaults for all in input.gcadj (dkh, 02/09/11) 
            !! BUG FIX: enforce defualt scaling factors before using SF_tmp 
            !!! (dkh, 07/30/10) 
            !EMS_SF(:,:,:,:) = 0.d0
            !
            !! Start with an initial guess for EMS_SF that is wrong
            !EMS_SF(:,:,1,EMSFD) = LOG(EMS_SF_tmp)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M ) 
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               EMS_SF(I,J,M,:) = LOG(EMS_SF_DEFAULT(:))

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF 
     
      ENDIF 

#endif 

      ! Save a copy of the initial guess of ICS_SF for regularization
      ICS_SF0(:,:,:,:) = ICS_SF(:,:,:,:)

      ! Save a copy of the initial guess of EMS_SF for regularization
      ! Add flags (hml, 02/23/12)
       IF ( LADJ_EMS ) EMS_SF0(:,:,:,:) = EMS_SF(:,:,:,:)


      ! Return to calling program
      END SUBROUTINE SET_LOG_SF

!-----------------------------------------------------------------------------

      SUBROUTINE GET_X_FROM_SF
!
!*****************************************************************************
!  Subroutine GET_X_FROM_ICS compiles the vector X of initial conditions from
!   the array STT_IC.  (dkh, 9/16/04)
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Rename to SET_LOG_SF, replace CMN_ADJ with adjoint_array_mod 
!        (dkh, ks, mak, cs  06/07/09) 
!  (4 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : ICS_SF, EMS_SF, MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, NFD
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LICS, LADJ_EMS
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT  
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS

#     include   "CMN_SIZE"        ! Size params

      ! Local variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM

      !=================================================================
      ! GET_X_FROM_SF begins here!
      !=================================================================
      IF ( LICS ) THEN 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, I_DUM)
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                   + (  IIPAR * JJPAR * ( L - 1 )  )
     &                   + (  IIPAR * JJPAR * LLPAR * ( N - 1 )  )

            ! Load X from active tracer concentrations
            X(I_DUM) = ICS_SF(I,J,L,N)

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSEIF ( LADJ_EMS ) THEN 

         IF ( ITS_A_TAGCO_SIM() ) THEN 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
            DO N = 1, 1
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                      + (  IIPAR * JJPAR * ( M - 1 )  )
     &                      + (  IIPAR * JJPAR * MMSCL * ( N - 1 )  )

               ! Load X from active tracer concentrations
               X(I_DUM) = EMS_SF(I,J,M,N)

            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            IF ( NNEMS == 2 ) THEN
               N = 2
               print*, IIPAR*JJPAR*MMSCL,'adding backgnd component to X'
               X(IIPAR*JJPAR*MMSCL+1) = EMS_SF(1,1,1,N)
            ENDIF 

         ELSE 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
            DO N = 1, NNEMS
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                      + (  IIPAR * JJPAR * ( M - 1 )  )
     &                      + (  IIPAR * JJPAR * MMSCL * ( N - 1 )  )

               ! Load X from active tracer concentrations
               X(I_DUM) = EMS_SF(I,J,M,N)

            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            IF ( LADJ_STRAT ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
               DO N = 1, NSTPL
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  
                  I_DUM    = I + (  IIPAR * ( J - 1)  ) 
     &                         + (  IIPAR * JJPAR * ( M - 1 )  )
     &                         + (  IIPAR * JJPAR * MMSCL * ( N - 1 )  )

                  ! Load X from active tracer concentrations
                  XP(I_DUM) = PROD_SF(I,J,M,N)
                  XL(I_DUM) = LOSS_SF(I,J,M,N)
                  IF ( I == IFD.and.J == JFD.and.N == NFD )THEN
                     print*, 'inverse_0: I_DUM = ' ,
     &                        I_DUM
                     print*, 'inverse_0: XL(I_DUM) = ' ,
     &                        XL(I_DUM)
                     print*, 'inverse_0: LOSS_SF = ' ,
     &                        LOSS_SF(I,J,M,N)
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               X( IIPAR*JJPAR*MMSCL*NNEMS + 1         :
     &            IIPAR*JJPAR*MMSCL*(NSTPL+NNEMS)) = XP(:)
               X( IIPAR*JJPAR*MMSCL*(NSTPL+NNEMS) + 1 :
     &            IIPAR*JJPAR*MMSCL*(2*NSTPL+NNEMS)) = XL(:)
       
            ENDIF 

         ENDIF  

      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_X_FROM_SF

!-----------------------------------------------------------------------------

      SUBROUTINE GET_SF_FROM_X
!
!*****************************************************************************
!  Subroutine GET_SF_FROM_X compiles the array of scaling factors  from
!   the vector X.  (dkh, 9/16/04)
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Rename to SET_LOG_SF, replace CMN_ADJ with adjoint_array_mod 
!        (dkh, ks, mak, cs  06/07/09) 
!  (4 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : ICS_SF, EMS_SF, NNEMS, MMSCL
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, NFD
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LICS, LADJ_EMS 
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS

#     include   "CMN_SIZE"        ! Size params

      ! Local Variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM

      !=================================================================
      ! GET_SF_FROM_X begins here!
      !=================================================================
      IF ( LICS ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, I_DUM)
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            I_DUM        = I + (  IIPAR * ( J - 1)  )
     &                       + (  IIPAR * JJPAR * ( L - 1 )  )
     &                       + (  IIPAR * JJPAR * LLPAR * ( N - 1 )  )

            ! Update the tracer concentrations from X
            ICS_SF(I,J,L,N) = X(I_DUM)

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSEIF ( LADJ_EMS ) THEN 

         IF ( ITS_A_TAGCO_SIM() ) THEN 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
            DO N = 1, 1
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               I_DUM        = I + (  IIPAR * ( J - 1)  )
     &                          + (  IIPAR * JJPAR * ( M - 1 )  )
     &                          + (  IIPAR * JJPAR * MMSCL * ( N - 1 ) )

               ! Update the tracer concentrations from X
               EMS_SF(I,J,M,N) = X(I_DUM)
   
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            IF ( NNEMS == 2 ) THEN
               N = 2
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Update the tracer concentrations from X
                  EMS_SF(I,J,M,N) = X(IIPAR*JJPAR*MMSCL+1)

               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO      
 
            ENDIF

         ELSE 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
            DO N = 1, NNEMS
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               I_DUM        = I + (  IIPAR * ( J - 1)  )
     &                          + (  IIPAR * JJPAR * ( M - 1 )  )
     &                          + (  IIPAR * JJPAR * MMSCL * ( N - 1 ) )

               ! Update the tracer concentrations from X
               EMS_SF(I,J,M,N) = X(I_DUM)
   
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! For strat prod and loss (hml)
            IF ( LADJ_STRAT ) THEN

               XP(:) = X(IIPAR*JJPAR*MMSCL*NNEMS+1:
     &               IIPAR*JJPAR*MMSCL*(NNEMS+NSTPL))
               XL(:) = X(IIPAR*JJPAR*MMSCL*(NNEMS+NSTPL)+1:
     &               IIPAR*JJPAR*MMSCL*(NNEMS+2*NSTPL))

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
               DO N = 1, NSTPL
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  I_DUM     = I + (  IIPAR * ( J - 1)  )
     &                          + (  IIPAR * JJPAR * ( M - 1 )  )
     &                          + (  IIPAR * JJPAR * MMSCL * ( N - 1 ) )


                  ! Update the tracer concentrations from X
                  PROD_SF(I,J,M,N) = XP(I_DUM)
                  LOSS_SF(I,J,M,N) = XL(I_DUM)
                  ! Debug (hml, 10/15/11)
                  IF ( I == IFD.and.J == JFD.and.N == NFD )THEN
                     print*, 'inverse_1: I_DUM = ' ,
     &                        I_DUM
                     print*, 'inverse_1: XL(I_DUM) = ' ,
     &                        XL(I_DUM)
                     print*, 'inverse_1: LOSS_SF = ' ,
     &                        LOSS_SF(I,J,M,N)
                  ENDIF

               ENDDO
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF 
         ENDIF 
      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_SF_FROM_X

!-----------------------------------------------------------------------------

      SUBROUTINE GET_GRADNT_FROM_ADJ
!
!*****************************************************************************
!  Subroutine GET_GRADNT_FROM_ADJ compiles the gradient vector from the array
!   of adjoint values.  (dkh, 9/16/04)
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Don't zero the NIT, NH4 and NO3 gradnts (dkh, 03/03/05)
!  (4 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL, NNEMS 
      USE ADJ_ARRAYS_MOD,    ONLY : ICS_SF_ADJ, EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_SF_ADJ, LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, LFD, NFD
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LICS, LADJ_EMS 
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS

#     include   "CMN_SIZE"  ! Size params

      ! Local Variables
      INTEGER :: I, J, L, M, N
      INTEGER :: I_DUM
      INTEGER :: I_DUM_TMP

      !=================================================================
      ! GET_GRADNT_FROM_ADJ begins here!
      !=================================================================

      IF ( LICS ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, I_DUM)
          DO N = 1, N_TRACERS
          DO L = 1, LLPAR
          DO J = 1, JJPAR
          DO I = 1, IIPAR

             I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                   + (  IIPAR * JJPAR * ( L - 1 )  )
     &                   + (  IIPAR * JJPAR * LLPAR * ( N - 1 )  )

             GRADNT(I_DUM) =  ICS_SF_ADJ(I,J,L,N)

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSEIF( LADJ_EMS ) THEN 

         IF ( ITS_A_TAGCO_SIM() ) THEN 
 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
            DO N = 1, 1
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                      + (  IIPAR * JJPAR * ( M - 1 )  )
     &                      + (  IIPAR * JJPAR * MMSCL * ( N - 1 )  )

               GRADNT(I_DUM) =  EMS_SF_ADJ(I,J,M,N)
   
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            IF( NNEMS == 2 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
               DO N = 2, 2 !NNEMS=2, but get zonal average for bkg
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  I_DUM    = (IIPAR*JJPAR*MMSCL) + 1
 
                  ! KLUDGE:  Ask MAK about this. 
                  ! sum zonally
                  GRADNT(I_DUM) =  GRADNT(I_DUM) + EMS_SF_ADJ(I,J,M,N)

               ENDDO
               ENDDO
               ENDDO

                  ! Update to include CH4 oxidation (zhej, 01/16/12, adj32_017) 
                  ! OLD:
                  !! KLUDGE:  Ask MAK about this. 
                  !! average zonally and per layer
                  !GRADNT(I_DUM) = GRADNT(I_DUM)
     &            !              / ( IIPAR * JJPAR * LLPAR * MMSCL)
                  ! NEW:
                  GRADNT(I_DUM) = GRADNT(I_DUM) /
     &                           ( IIPAR * JJPAR * MMSCL)

               ENDDO
!$OMP END PARALLEL DO

            ENDIF 

         ELSE 

            IF ( .NOT. LADJ_STRAT ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
               DO N = 1, NNEMS
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                         + (  IIPAR * JJPAR * ( M - 1 )  )
     &                         + (  IIPAR * JJPAR * MMSCL * ( N - 1 )  )

                  GRADNT(I_DUM) =  EMS_SF_ADJ(I,J,M,N)
   
               ENDDO
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO


            ELSEIF ( LADJ_STRAT ) THEN

               ! For strat prod & loss (hml, 08/29/11)
               I_DUM_TMP = IIPAR * JJPAR * MMSCL * NNEMS

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, N, I_DUM)
               DO N = 1, NSTPL
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  I_DUM    = I + (  IIPAR * ( J - 1)  )
     &                         + (  IIPAR * JJPAR * ( M - 1 )  )
     &                         + (  IIPAR * JJPAR * MMSCL * ( N - 1 )  )

                  GRADNT_P(I_DUM) = PROD_SF_ADJ(I,J,M,N)
                  GRADNT_L(I_DUM) = LOSS_SF_ADJ(I,J,M,N)
                 ! Debug (hml, 10/15/11)
                 IF ( I == IFD.and.J == JFD.and.N == NFD )THEN
                     print*, 'inverse_2: GRADNT_L = ' ,
     &                        GRADNT_L(I_DUM)
                 ENDIF

               ENDDO
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               GRADNT( I_DUM_TMP + 1 :
     &               I_DUM_TMP + IIPAR * JJPAR * MMSCL * NSTPL )
     &               = GRADNT_P(:)
               GRADNT( I_DUM_TMP + IIPAR * JJPAR * MMSCL * NSTPL + 1 :
     &               IIPAR * JJPAR * MMSCL * 2 * NSTPL )
     &               = GRADNT_L(:)
               ! Debug (hml, 10/15/11)
               print*, 'inverse_3: GRADNT = ' ,
     &                  GRADNT( I_DUM_TMP + IFD * JFD * 2*NFD )

            ENDIF
         ENDIF

      ENDIF 

      ! Return to calling program
      END SUBROUTINE GET_GRADNT_FROM_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_GDT_FILE( )
!
!******************************************************************************
!  Subroutine MAKE_GDT_FILE creates a binary file of ADJ_xxx
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) ICS_SF_ADJ : Array of adjoint gradients to be written
!  (3 ) EMS_SF_ADJ : Array of adjoint gradients to be written
!
!  NOTES:
!  (1 ) Just like MAKE_OBS_FILE except
!       - write to .adj. file
!  (2 ) Changed name to MAKE_GDT_FILE.  Now the .adj. files are trajectories,
!       and the .gdt. files are final gradients  (dkh, 10/03/04)
!  (3 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (4 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (5 ) Now use CATEGORY = 'IJ-GDE-$' for 'EMISSIONS' case. (dkh, 03/29/05)
!  (6 ) No longer pass COST_FUNC in the header; use cnf.* files. (dkh, 02/13/06)  
!  (7 ) Rename everything, replace CMN_ADJ, move nonessential stuff
!       to diagnostic files  (dkh, ks, mak, cs  06/07/09) 
!  (8 ) Add normalized gradients IJ-GDEN$ (dkh, 05/06/10) 
!  (9 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : ICS_SF_ADJ, EMS_SF_ADJ 
      USE ADJ_ARRAYS_MOD,      ONLY : NNEMS, MMSCL 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,      ONLY : COST_FUNC
      USE ADJ_ARRAYS_MOD,      ONLY : PROD_SF_ADJ, LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NSTPL
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE GRID_MOD,            ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_ADJ_MOD,     ONLY : LICS, LADJ_EMS 
      USE LOGICAL_ADJ_MOD,     ONLY : LADJ_STRAT
      USE LOGICAL_MOD,         ONLY : LPRT
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,          ONLY : N_TRACERS

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NEMIS(NCS)

      ! Local Variables
      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)               :: EMS_3D(IIPAR,JJPAR,MMSCL)
      CHARACTER(LEN=255)   :: FILENAME
      TYPE (XPLEX)               :: PROD_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)               :: LOSS_3D(IIPAR,JJPAR,MMSCL)
     
      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_GDT_FILE begins here!
      !=================================================================

      ! Clear intermediate arrays
      EMS_3D (:,:,:) = 0d0
      PROD_3D(:,:,:) = 0d0
      LOSS_3D(:,:,:) = 0d0

      ! Hardwire output file for now
      OUTPUT_GDT_FILE = 'gctm.gdt.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM GDT File: ' //
     &           'Final gradient values '
      UNIT     = 'none' 
      CATEGORY = 'IJ-GDT-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_GDT_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      IF ( LICS ) THEN 
         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================
         DO N = 1, N_TRACERS

            !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
              TRACER(I,J,L) = ICS_SF_ADJ (I,J,L,N)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &                  J0+1,      1,         TRACER )
         ENDDO

      ENDIF 
      IF ( LADJ_EMS ) THEN 

         ! Reset CATEGORY as labeling in gamap is different 
         CATEGORY = 'IJ-GDE-$'

         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================
         DO N = 1, NNEMS

            !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               EMS_3D(I,J,M) = XPLX(EMS_SF_ADJ(I,J,M,N))
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
     &                  J0+1,      1,         EMS_3D )
         ENDDO

! Reset CATEGORY as labeling in gamap is different 
         CATEGORY = 'IJ-GDEN$'
         UNIT     = 'none'

         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================
         DO N = 1, NNEMS

            !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               EMS_3D(I,J,M) = XPLX(EMS_SF_ADJ(I,J,M,N))
     &                       / COST_FUNC
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
            
            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
     &                  J0+1,      1,         EMS_3D )
         ENDDO

         ! Strat prod and loss (hml)
         IF ( LADJ_STRAT ) THEN


            !==============================================================
            ! Write each observed quantity to the observation file
            !==============================================================
            DO N = 1, NSTPL
               
               !Temporarily store quantities in the PROD_3D, LOSS_3D array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                 PROD_3D(I,J,M) = XPLX(PROD_SF_ADJ(I,J,M,N))
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               CATEGORY = 'IJ-GDP-$'
               UNIT     = 'J'
               CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY,  N,
     &                     UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                     IIPAR,     JJPAR,     MMSCL,         I0+1,
     &                     J0+1,      1,         PROD_3D )

            ENDDO

            ! Strat loss
            DO N = 1, NSTPL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                 LOSS_3D(I,J,M) = XPLX(LOSS_SF_ADJ(I,J,M,N))
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               CATEGORY = 'IJ-GDL-$'
               UNIT     = 'J'
               CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY,  N,
     &                     UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                     IIPAR,     JJPAR,     MMSCL,         I0+1,
     &                     J0+1,      1,         LOSS_3D )


            ENDDO
         ENDIF
      ENDIF

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_GDT_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_GDT_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_GDT_FILE ( )
!
!******************************************************************************
!  Subroutine READ_GDT_FILE reads the gctm.gdt file into ADJ_xxx
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Notes
!  (1 ) now called GDT instead of ADJ
!  (2 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (3 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (4 ) Now use CATEGORY = 'IJ-GDE-$' for EMISSIONS case. (dkh, 03/29/05)
!  (5 ) No longer pass COST_FUNC in the header; use cnf.* files. (dkh, 02/13/06)  
!  (6 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : ICS_SF_ADJ, EMS_SF_ADJ 
      USE ADJ_ARRAYS_MOD,      ONLY : NNEMS, MMSCL 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,      ONLY : PROD_SF_ADJ,LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NSTPL
      USE BPCH2_MOD,           ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST, IOERROR
      USE LOGICAL_ADJ_MOD,     ONLY : LICS, LADJ_EMS 
      USE LOGICAL_ADJ_MOD,     ONLY : LADJ_STRAT
      USE LOGICAL_MOD,         ONLY : LPRT
      USE RESTART_MOD,         ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,            ONLY : EXPAND_DATE
      USE TRACER_MOD,          ONLY : N_TRACERS

#     include "CMN_SIZE"   ! Size parameters

      ! Local Variables
      INTEGER             :: I, IOS, J, L, M, N
      INTEGER             :: NCOUNT(NNPAR)
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: EMS_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)              :: PROD_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)              :: LOSS_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16              :: D_EMS_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16              :: D_PROD_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16              :: D_LOSS_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER(LEN=20)   :: INPUT_GDT_FILE

      !=================================================================
      ! READ_GDT_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_GDT_FILE = 'gctm.gdt.NN'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) =0e0
      !=================================================================
      ! Open gradient file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_GDT_FILE )

      ! Replace NN tokens in FILENAME w/ actual values
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add OPTDATA_DIR prefix to FILENAME
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'G D T   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_GDT_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )


      IF ( LICS ) THEN 
         !=================================================================
         ! Read adjoints -- store in the TRACER array
         !=================================================================
         DO N = 1, N_TRACERS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES , D_LATRES, HALFPOLAR, CENTER180
            LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES) 
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
            ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
            TRACER%r = dble(D_TRACER)
            TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:6')

            !==============================================================
            ! Assign data from the TRACER array to the ADJ_STT array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-GDT-$' ) THEN

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  ICS_SF_ADJ(I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

      ENDIF 
      IF ( LADJ_EMS ) THEN 

         !=================================================================
         ! Read adjoints -- store in the TRACER array
         !=================================================================
         DO N = 1, NNEMS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
            LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
            ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_EMS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
            EMS_3D%r = dble(D_EMS_3D)
            EMS_3D%i = dimag(D_EMS_3D)         
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:7')

            !==============================================================
            ! Assign data from the TRACER array to the ADJ_STT array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-GDE-$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  EMS_SF_ADJ(I,J,M,N) = EMS_3D(I,J,M)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

         IF ( LADJ_STRAT ) THEN
            !==============================================================
            ! Read adjoints -- store in the TRACER array
            !==============================================================
            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, D_LONRES , D_LATRES, HALFPOLAR, CENTER180
               LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a real I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR
     &                             ( IOS,IU_RST,'read_gdt_file:8' )

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, D_ZTAU0,D_ZTAU1,RESERVED,
     &              NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &              NSKIP
               ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
               IF ( IOS /= 0 ) CALL IOERROR
     &                              ( IOS,IU_RST,'read_gdt_file:9')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( D_PROD_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
               IF ( IOS /= 0 ) CALL IOERROR
     &                              ( IOS,IU_RST,'read_gdt_file:10')
             PROD_3D%r = dble(D_PROD_3D)
            PROD_3D%i = dimag(D_PROD_3D)

               !===========================================================
               ! Assign data from the TRACER array to the ADJ_STT array.
               !===========================================================

               ! Only process observation data (i.e. aerosol and precursors)
               IF ( CATEGORY(1:8) == 'IJ-GDP-$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     PROD_SF_ADJ(I,J,M,N) = PROD_3D(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO

            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, D_LONRES , D_LATRES, HALFPOLAR, CENTER180
              LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a real I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR
     &                             ( IOS,IU_RST,'read_gdt_file:8b')

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, D_ZTAU0,D_ZTAU1,RESERVED,
     &              NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &              NSKIP
             ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
               IF ( IOS /= 0 ) CALL IOERROR
     &                              ( IOS,IU_RST,'read_gdt_file:9b')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( D_LOSS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
               LOSS_3D%r = dble(D_LOSS_3D)
               LOSS_3D%i = dimag(D_LOSS_3D)
               IF ( IOS /= 0 ) CALL IOERROR
     &                              (IOS,IU_RST,'read_gdt_file:10b')

               ! Only process observation data (i.e. aerosol and precursors)
               IF ( CATEGORY(1:8) == 'IJ-GDL-$' ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     LOSS_SF_ADJ(I,J,M,N) = LOSS_3D(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_GDT_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_GDT_FILE

! needs to be updated 
!-----------------------------------------------------------------------
!
!      SUBROUTINE MAKE_GDT_DIAG_FILE( )
!!
!!******************************************************************************
!!  Subroutine MAKE_GDT_DIAG_FILE creates a binary file of daignostics
!!  relatied to the adjoint gradients.  (dkh, 06/07/09) 
!!  (dkh, 9/17/04)
!!
!!  Module Variable as Input:
!!  ============================================================================
!!  (1 ) N_CALC       : Current iteration number
!!  (2 ) ICS_SF_ADJ   : Array of adjoint gradients to be written
!!  (3 ) EMS_SF_ADJ   : Array of adjoint gradients to be written
!!  (4 ) ADJ_BURNEMIS : Array of biomass burning sensitivities
!!  (5 ) ADJ_BIOFUEL  : Array of biofuel sensitivities 
!!  (6 ) ADJ_EMISRR   : 
!!  (7 ) ADJ_EMISRRB  : 
!!
!!  NOTES:
!! 
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE ADJ_ARRAYS_MOD,     ONLY : ADJ_BURNEMIS
!      USE ADJ_ARRAYS_MOD,     ONLY : ADJ_BIOFUEL
!      USE ADJ_ARRAYS_MOD,     ONLY : ADJ_EMISRR
!      USE ADJ_ARRAYS_MOD,     ONLY : ADJ_EMISRRB
!      USE BPCH2_MOD
!      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJTMP_DIR 
!      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
!      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
!      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
!      USE BIOMASS_MOD,       ONLY : NBIOTRCE
!      USE BIOFUEL_MOD,       ONLY : NBFTRACE
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "CMN"        ! LPRT
!#     include "comode.h"   ! NEMIS(NCS)
!
!      ! Local Variables
!      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
!      INTEGER              :: YYYY, MM, DD,  HH, SS
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: EMS_3D(IIPAR,JJPAR,MMSCL)
!      CHARACTER(LEN=255)   :: FILENAME
!     
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER, PARAMETER   :: HALFPOLAR = 1
!      INTEGER, PARAMETER   :: CENTER180 = 1
!
!      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE
!
!      !=================================================================
!      ! MAKE_GDT_FILE begins here!
!      !=================================================================
!
!      ! Clear intermediate arrays
!      EMS_3D(:,:,:) = 0d0
!
!      ! Hardwire output file for now
!      OUTPUT_GDT_FILE = 'gctm.gdt.diag.NN'
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM GDT File: ' //
!     &           'Gradient diagnostics  '
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the adjoint file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output observation file name into a local variable
!      FILENAME = TRIM( OUTPUT_GDT_FILE )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, N_CALC )
!
!      ! Add the OPTDATA_DIR prefix to the file name
!      FILENAME = TRIM( DIAGADJTMP_DIR ) //  TRIM( FILENAME )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_GDT_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      !=================================================================
!      ! Normalized sensitivies 
!      !=================================================================
!
!      ! Reset CATEGORY as labeling in gamap is different 
!      CATEGORY = 'IJ-GDEN$'
!      UNIT     = '%'
!
!      !=================================================================
!      ! Write each observed quantity to the observation file
!      !=================================================================
!      DO N = 1, NNEMS
!
!         ! Temporarily store quantities in the TRACER array
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M )
!         DO M = 1, MMSCL
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            EMS_3D(I,J,M) = REAL(ADJ_EMS(I,J,M,N)) / COST_FUNC
!     &                    * 100d0 
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N + NNEMS,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &               J0+1,      1,         EMS_3D )
!
!      ENDDO
!
!
!      !=================================================================
!      ! Normalized VOC sensitivies - EMISRR (anthro hydrocarbons)
!      !=================================================================
!      CATEGORY = 'DEMISRR' 
!      DO N = 1, NEMIS(NCS)
!
!         ! Temporarily store quantities in the TRACER array
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            EMS_3D(I,J,1) = REAL(ADJ_EMISRR(I,J,N)) / COST_FUNC
!     &                     * 100d0 
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         ! dkh debug 
!         print*, ' ADJ EMISRR = ', maxval(adj_emisrr(:,:,n)), n
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &               J0+1,      1,         EMS_3D )
!      ENDDO
!
! 
!      !=================================================================
!      ! Normalized VOC sensitivies - EMISRRB (biogenic hydrocarbons)
!      !=================================================================
!      CATEGORY = 'DEMISRRB' 
!      DO N = 1, NEMIS(NCS)
!
!         ! Temporarily store quantities in the TRACER array
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            EMS_3D(I,J,1) = REAL(ADJ_EMISRRB(I,J,N)) / COST_FUNC
!     &                    * 100d0 
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         ! dkh debug 
!         print*, ' ADJ EMISRRB = ', maxval(adj_emisrrb(:,:,n)), n
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &               J0+1,      1,         EMS_3D )
!      ENDDO
!
!      !=================================================================
!      ! Normalized VOC sensitivies - BOIFUEL 
!      !=================================================================
!      CATEGORY = 'DBIOFUEL'
!      DO N = 1, NBFTRACE
!
!         ! Temporarily store quantities in the TRACER array
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            EMS_3D(I,J,1) = REAL(ADJ_BIOFUEL(I,J,N)) / COST_FUNC
!     &                    * 100d0 
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &               J0+1,      1,         EMS_3D )
!
!         ! dkh debug 
!         print*, ' ADJ BIOFUEL= ', maxval(ADJ_BIOFUEL(:,:,n)), n
!
!      ENDDO
!
!
!      !=================================================================
!      ! Normalized VOC sensitivies - BURNEMIS 
!      !=================================================================
!      CATEGORY = 'DBURNEMIS' 
!      DO N = 1, NBIOTRCE
!
!         ! Temporarily store quantities in the TRACER array
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            EMS_3D(I,J,1) = REAL(ADJ_BURNEMIS(I,J,N)) / COST_FUNC
!     &                    * 100d0 
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &               J0+1,      1,         EMS_3D )
!
!        ! dkh debug 
!         print*, ' ADJ BURNEMIS= ', maxval(ADJ_BURNEMIS(:,:,n)), n
!
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_RST )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_GDT_DIAG_FILE: wrote file' )
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_GDT_DIAG_FILE
!
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_SF_FILE ( )
!
!******************************************************************************
!  Subroutine MAKE_SF_FILE creates a binary file of STT_IC or EMS_ICS
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!  (2 ) ICS_SF    : Initial conditions scaling factors 
!  (3 ) EMS_SF    : Emissions scaling factors 
!
!  NOTES:
!  (1 ) Just like MAKE_ADJ_FILE except
!       - write to .ics. file
!  (2 ) Add support for ACTIVE_VARS == 'EMISSIONS' case (dkh, 11/27/04)
!  (3 ) Add support for ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (4 ) Change UNIT to unitless and change title to Scale factors (dkh, 03/06/05)
!  (5 ) Change output for ACTIVE_VARS == 'EMISSIONS' case.
!        Now use label IJ-EMS-$, and update gamap code accordingly. 
!        First write the scaling factors, in consecutive species. Temporal 
!         varations in the emissions, if any, will be in the L direction.
!        Next, write out the optimized emissions themselves.
!        Finally, write out the difference between orig and optimized emissions.
!        (dkh, 03/28/05)
!  (6 )  Use EMS_orig instead of ESO4_an_orig so that we can loop over N.
!  (7 )  Update to add support for writing NOx emissions. (dkh, 08/27/06)  
!  (8 )  Only write the value of the scaling facotr in locations where the 
!         actual emission is greater than zero.  Also include the current 
!         scale emissions themselves in every *ics* file.  (dkh, 09/22/06)  
!  (9 )  Add suppport for LOG_OPT
!  (10)  Standardize units for saving emissions. (dkh, 06/16/07) 
!  (11)  Add option to print prior and posterior emissions totals. (dkh, 06/16/07) 
!  (12)  Change names, replace CMN_ADJ. (dkh, ks, mak, cs  06/08/09) 
!  (13)  Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC, ICS_SF, EMS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC, ICS_SF, EMS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE BPCH2_MOD
      USE DIRECTORY_MOD,     ONLY : TEMP_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE GRID_MOD,          ONLY : GET_AREA_CM2
      USE LOGICAL_ADJ_MOD,   ONLY : LICS, LADJ_EMS
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT
      USE LOGICAL_MOD,       ONLY : LPRT
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TIME_MOD,          ONLY : GET_TS_CHEM
      USE TIME_MOD,          ONLY : GET_TS_EMIS
      USE TIME_MOD,          ONLY : GET_TAUb, GET_TAUe
      USE TRACER_MOD,        ONLY : N_TRACERS


#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_O3"     ! EMISRN
#     include "comode.h"   ! NEMIS(NCS)

      ! Local Variables
      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      INTEGER              :: NOFFSET
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)   :: FILENAME
      TYPE (XPLEX)               :: TEMP 
      TYPE (XPLEX)               :: NEMIS_DT
      TYPE (XPLEX)               :: USA_MASK(IIPAR,JJPAR)
      TYPE (XPLEX)               :: EMS_TOTAL(NNEMS)
      TYPE (XPLEX)               :: EMS_PERCENT(IIPAR,JJPAR,NNEMS)
      LOGICAL, PARAMETER   :: LPRINT_TOTAL = .TRUE. 

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: OUTPUT_SF_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      ! Parameters
      TYPE (XPLEX), PARAMETER    :: SEC_PER_YEAR = xplex(3.1536d7,0d0)
      TYPE (XPLEX), PARAMETER    :: MIN_PER_YEAR = xplex(5.2560d5,0d0)
      TYPE (XPLEX), PARAMETER    :: TG_PER_KG    = xplex(1d-09,0d0)

      !=================================================================
      ! MAKE_SF_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_SF_FILE = 'gctm.sf.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM SF File: ' //
     &           'Scale Factors'
      UNIT     = 'unitless'
      CATEGORY = 'IJ-ICS-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_SF_FILE )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add OPTDATA_DIR prefix to FILENAME
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_SF_FILE:  Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      IF ( LICS ) THEN

         CATEGORY = 'IJ-ICS-$'

         !=================================================================
         ! Write each observed quantity to the ics file
         !=================================================================
         DO N = 1, N_TRACERS

            !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               TRACER(I,J,L) = ICS_SF(I,J,L,N)
            ENDDO
            ENDDO
            ENDDO

!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &                  J0+1,      1,         TRACER )
         ENDDO

      ENDIF 

      IF ( LADJ_EMS ) THEN 

         CATEGORY = 'IJ-EMS-$'
         UNIT     = 'unitless'

         !=================================================================
         ! Write each observed quantity to the ics file
         !=================================================================
         DO N = 1, NNEMS
         
            !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M , TEMP )
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               TRACER(I,J,M) = EMS_SF(I,J,M,N)

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
     &                  J0+1,      1,         TRACER )
         ENDDO

         IF ( LADJ_STRAT ) THEN

            !==============================================================
            ! Write each observed quantity to the ics file
            !==============================================================
            DO N = 1, NSTPL

               !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M ) 
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  TRACER(I,J,M) = PROD_SF(I,J,M,N)
               ENDDO
               ENDDO
               ENDDO

!$OMP END PARALLEL DO

               CATEGORY = 'IJ-STRP$'

               CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY,  N,
     &                     UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                     IIPAR,     JJPAR,     MMSCL,     I0+1,
     &                     J0+1,      1,         TRACER )
            ENDDO
            
            DO N = 1, NSTPL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M ) 
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  TRACER(I,J,M) = LOSS_SF(I,J,M,N)
               ENDDO
               ENDDO
               ENDDO

!$OMP END PARALLEL DO
               
               CATEGORY = 'IJ-STRL$'
               
               CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY,  N,
     &                     UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                     IIPAR,     JJPAR,     MMSCL,     I0+1,
     &                     J0+1,      1,         TRACER )
            
            ENDDO
         ENDIF
      ENDIF 


      !### Debug

      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SF_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_SF_FILE
! needs to be updated:
!!------------------------------------------------------------------------------
!
!      SUBROUTINE MAKE_SF_DIAG_FILE ( )
!!
!!******************************************************************************
!!  Subroutine MAKE_SF_DIAG_FILE creates a binary file of diagnostics
!!  related to scaling factor values. (dkh, 06/08/09) 
!!
!!  Module Variable as Input:
!!  ============================================================================
!!  (1 ) N_CALC    : Current iteration number
!!  (2 ) ICS_SF    : Initial conditions scaling factors 
!!  (3 ) EMS_SF    : Emissions scaling factors 
!!
!!  NOTES:
!!  (1)  Split this off from MAKE_ICS_FILE  (dkh, ks, mak, cs  06/08/09) 
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF, ICS_SF
!      USE BPCH2_MOD
!      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
!      USE DIRECTORY_MOD,     ONLY : TEMP_DIR
!      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJTMP_DIR 
!      USE FILE_MOD,  ONLY : IU_RST,      IOERROR
!      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET, GET_AREA_CM2
!      USE TIME_MOD,  ONLY : EXPAND_DATE, GET_TAU
!      USE TIME_MOD,  ONLY : GET_TS_CHEM
!      USE TIME_MOD,  ONLY : GET_TS_EMIS
!      USE TIME_MOD,  ONLY : GET_TAUb, GET_TAUe
!      USE TRACERID_MOD, ONLY : IDTNH3, IDTNOX, IDTBCPI, IDTSO2
!      USE SULFATE_MOD, ONLY : EMS_orig
!      USE LIGHTNING_NOX_MOD, ONLY : EMS_orig_li
!      USE EMISSIONS_MOD, ONLY : BIOFUEL_orig
!      USE EMISSIONS_MOD, ONLY : BURNEMIS_orig
!      USE EMISSIONS_MOD, ONLY : EMISRR_orig
!      USE EMISSIONS_MOD, ONLY : EMISRRB_orig
!      USE BIOMASS_MOD, ONLY   : NBIOTRCE
!      USE BIOFUEL_MOD, ONLY   : NBFTRACE
!      USE DAO_MOD,     ONLY   : BXHEIGHT
!
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "CMN"        ! LPRT, LLIGHTNOX
!#     include "CMN_O3"     ! EMISRN
!#     include "comode.h"   ! NEMIS(NCS)
!
!      ! Local Variables
!      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
!      INTEGER              :: YYYY, MM, DD,  HH, SS
!      INTEGER              :: NOFFSET
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: TRACER_VOC(IIPAR,JJPAR,20)
!      TYPE (XPLEX)               :: TRACER_US(IIPAR,JJPAR,LLPAR)
!      CHARACTER(LEN=255)   :: FILENAME
!      TYPE (XPLEX)               :: TEMP 
!      TYPE (XPLEX)               :: NEMIS_DT
!      TYPE (XPLEX)               :: USA_MASK(IIPAR,JJPAR)
!      TYPE (XPLEX)               :: EMS_TOTAL(NNEMS)
!      TYPE (XPLEX)               :: EMS_PERCENT(IIPAR,JJPAR,NNEMS)
!      LOGICAL, PARAMETER   :: LPRINT_TOTAL = .TRUE. 
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER, PARAMETER   :: HALFPOLAR = 1
!      INTEGER, PARAMETER   :: CENTER180 = 1
!
!      CHARACTER(LEN=20)    :: OUTPUT_SF_DIAG_FILE
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE
!
!      ! Parameters
!      TYPE (XPLEX), PARAMETER    :: SEC_PER_YEAR = 3.1536d7
!      TYPE (XPLEX), PARAMETER    :: MIN_PER_YEAR = 5.2560d5
!      TYPE (XPLEX), PARAMETER    :: TG_PER_KG    = 1d-09
!
!      !=================================================================
!      ! MAKE_SF_DIAG_FILE begins here!
!      !=================================================================
!
!      ! Hardwire output file for now
!      OUTPUT_SF_DIAG_FILE = 'gctm.sf.diag.NN'
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM SF  File: ' //
!     &           'Scale Factors Diagnostics'
!      UNIT     = 'unitless'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the adjoint file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output observation file name into a local variable
!      FILENAME = TRIM( OUTPUT_SF_DIAG_FILE )
!
!      ! Replace NN token w/ actual value
!      CALL EXPAND_NAME( FILENAME, N_CALC )
!
!      ! Add OPTDATA_DIR prefix to FILENAME
!      FILENAME = TRIM( DIAGADJTMP_DIR ) // TRIM( FILENAME )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_SF_DIAG_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      IF ( NEMS ) THEN 
!
!
!         ! Also write the actual emissions.
!         ! Go ahead and include this every time.
!         CATEGORY = 'IJ-EM0-$'
!         UNIT     = 'molecule/cm2/s'
!    
!         ! emdt / sim = hr / sim * min / hr * emdt / min
!         NEMIS_DT = ( GET_TAUe() - GET_TAUb() ) * 60d0 / GET_TS_EMIS()
!            
!         DO N = 1, NNEMS
!
!            ! Compile TRACER [ molec / cm2 / s ] 
!            ! Get the actual emission in the current cell
!            ! Original emissions are in EMS_orig, but in a variety of unit
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!                  
!               ! lightning NOx
!               IF ( N == IDADJEMS_ENOxli ) THEN 
! 
!                  ! Add to prevent allocation segfault (dkh, 10/10/08) 
!                  IF ( LLIGHTNOX ) THEN 
!   
!                     ! molec NOx / cm2 / s total sim  -> molec NOx / cm2 / s
!                     TRACER(I,J,N) = EMS_orig_li(I,J)
!     &                             / NEMIS_DT                 ! number of emissions
!
!                  ELSE 
!                      TRACER(I,J,N)  = 0D0
!                  ENDIF 
!
! 
!               ! soil NOx
!               ELSEIF ( N == IDADJEMS_ENOxso ) THEN 
!   
!                  ! molec NOx / cm2 / s total -> molec / cm2 / s 
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          / NEMIS_DT                 ! number of emissions
!
! 
!               ! BC / OC
!               ELSEIF ( N == IDADJEMS_BCan .or. N == IDADJEMS_OCan .or. 
!     &                  N == IDADJEMS_BCbb .or. N == IDADJEMS_OCbb .or.  
!     &                  N == IDADJEMS_BCbf .or. N == IDADJEMS_OCbf     )
!     &            THEN 
!
!                  ! Convert from kg / yr to molec C / cm2 / s
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          * XNUMOL(IDTBCPI)
!     &                          / GET_AREA_CM2(J)
!     &                          / SEC_PER_YEAR
!
!
!               ! Anth NOx
!               ELSEIF ( N == IDADJEMS_ENOx1 .or. N == IDADJEMS_ENOx2 )
!     &            THEN
!
!                  ! Convert from kg / box / emdt to molec / cm2 / s
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          * XNUMOL(IDTNOX)
!     &                          / GET_AREA_CM2(J)
!     &                          / ( GET_TS_EMIS() * 60.d0 )  ! seconds per emdt
!
!               ! NH3
!               ELSEIF ( N == IDADJEMS_ENH3_an .or.
!     &                  N == IDADJEMS_ENH3_na .or. 
!     &                  N == IDADJEMS_ENH3_bb .or. 
!     &                  N == IDADJEMS_ENH3_bf )
!     &            THEN
!
!                  ! Convert from kg NH3 / box / s to molec / cm2 / s
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          * XNUMOL(IDTNH3)
!     &                          / GET_AREA_CM2(J)
!
!               ! SO2
!               ELSEIF ( N == IDADJEMS_ESO2_bb .or.
!     &                  N == IDADJEMS_ESO2_bf .or. 
!     &                  N == IDADJEMS_ESO2_sh )
!     &            THEN
!
!                  ! Convert from kg SO2 / box / s to molec / cm2 / s
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          * XNUMOL(IDTSO2) 
!     &                          / GET_AREA_CM2(J)
!
!               ! Volcano SO2 emissions (dkh, cklee 09/14/08) 
!               ELSEIF ( N == IDADJEMS_ESO2_ev .or.    !(added,cklee)
!     &                  N == IDADJEMS_ESO2_nv )       !(added,cklee)
!     &            THEN 
!                  ! Convert from kg SO2 / box / s total to molec / cm2 / s
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          * XNUMOL(IDTSO2) 
!     &                          / GET_AREA_CM2(J)
!     &                          / NEMIS_DT
!
!               ! Anth SOx
!               ELSEIF ( N == IDADJEMS_ESOx1 .or. 
!     &                  N == IDADJEMS_ESOx2      ) THEN 
!
!                  ! it's already in molec SOx / cm2 / s
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!
!               ELSE 
!
!                  CALL ERROR_STOP('undefined emissions', 
!     &                           'inverse_mod.f')
!               ENDIF
!
!#if defined ( LOG_OPT )
!               ! Apply current scaling
!               TRACER(I,J,N) = TRACER(I,J,N) * EXP(EMS_ICS(I,J,1,N))
!
!#else
!               ! Apply current scaling
!               TRACER(I,J,N) = TRACER(I,J,N) * EMS_ICS(I,J,1,N)
!#endif 
!
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!               
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  50+N,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     1,         I0+1,
!     &                  J0+1,      1,         TRACER(:,:,N) )
!            
!         ENDDO
!
!
!         ! Also write the normalized emissions
!         CALL READ_USA_MASK( USA_MASK )
!         CATEGORY = 'IJ-EMP-$'
!         UNIT     = '%'
!
!         EMS_TOTAL(:)       = 0d0
!         EMS_PERCENT(:,:,:) = 0d0
!
!         DO N = 1, NNEMS
!
!            ! Note: not in parallel, would need another tmp array for that
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!
!               IF ( USA_MASK(I,J) > 0d0 ) THEN
!               
!                  EMS_TOTAL(N) = EMS_TOTAL(N) 
!     &                         + TRACER(I,J,N) * GET_AREA_CM2(J)
!               ENDIF 
!
!            ENDDO
!            ENDDO
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!
!               IF ( EMS_TOTAL(N) == 0d0 .or.
!     &              USA_MASK(I,J) == 0d0 ) THEN
!
!                  ! Not sure what to store as "actual emission" for lightning NOx
!                  EMS_PERCENT(I,J,N) = 0d0
!
!               ELSE
!
!                  ! emissions percentages 
!                  EMS_PERCENT(I,J,N) = TRACER(I,J,N) * GET_AREA_CM2(J)
!     &                               / EMS_TOTAL(N) * 100d0
!
!               ENDIF
!
!
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  50+N+NNEMS,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     1,         I0+1,
!     &                  J0+1,      1,         EMS_PERCENT(:,:,N) )
!
!            ! dkh debug
!            print*, 'EMS_PERCENT total = ', SUM(EMS_PERCENT(:,:,N)), N
!
!         ENDDO
! 
!         !NOFFSET = 0
!
!         ! VOC emissions -- anth hydrocarbons (EMISRR)
!         CATEGORY  = 'EMISRR' 
!         print*, 'make_ics db: nemis = ', NEMIS(NCS)
!         DO N = 1, NEMIS(NCS)
!
!            ! Compile TRACER [ molec / cm2 / s ] 
!            ! Get the actual emission in the current cell
!            ! Original emissions are in EMS_orig, but in a variety of unit
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!                  
!               ! molec C / box / s total sim  -> molec C / cm2 / s
!               TRACER_VOC(I,J,N) = EMISRR_orig(I,J,N)
!     &                          / GET_AREA_CM2(J)
!     &                          / NEMIS_DT
!
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!               
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N, 
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     1,         I0+1,
!     &                  J0+1,      1,         TRACER_VOC(:,:,N) )
!            
!            ! dkh debug
!            print*, 'max EMISRR = ', MAXVAL(EMISRR_orig(:,:,N)), N
!
!         ENDDO
!        
!         ! VOC emissions -- biogenic hydrocarbons (EMISRRB)
!         CATEGORY  = 'EMISRRB' 
!         DO N = 1, NEMIS(NCS)
!
!            ! Compile TRACER [ molec / cm2 / s ] 
!            ! Get the actual emission in the current cell
!            ! Original emissions are in EMS_orig, but in a variety of unit
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!                  
!               ! molec C / box / s total sim  -> molec C / cm2 / s
!               TRACER_VOC(I,J,N) = EMISRRB_orig(I,J,N)
!     &                          / GET_AREA_CM2(J)
!     &                          / NEMIS_DT
!
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!               
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N, 
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     1,         I0+1,
!     &                  J0+1,      1,         TRACER_VOC(:,:,N) )
!            
!            ! dkh debug
!            print*, 'max EMISRRB = ', MAXVAL(EMISRRB_orig(:,:,N)), N
!
!         ENDDO
!         !NOFFSET = NOFFSET + NEMIS(NCS)
!
!         ! VOC emissions - BIOFUEL
!         CATEGORY  = 'BIOFUEL'
!         DO N = 1, NBFTRACE
!
!            ! Compile TRACER [ molec / cm2 / s ] 
!            ! Get the actual emission in the current cell
!            ! Original emissions are in EMS_orig, but in a variety of unit
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!                  
!               ! molec C / cm3 / s total sim  -> molec C / cm2 / s
!               TRACER_VOC(I,J,N) = BIOFUEL_orig(N,I,J)
!     &                           * BXHEIGHT(I,J,1) * 100d0
!     &                           / NEMIS_DT
!
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!               
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N, 
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     1,         I0+1,
!     &                  J0+1,      1,         TRACER_VOC(:,:,N) )
!            
!            ! dkh debug
!            print*, 'max BIOFUEL = ', MAXVAL(BIOFUEL_orig(N,:,:)), N
!
!         ENDDO
!
!         !NOFFSET = NOFFSET + NBFTRACE
!
!         ! VOC emissions - BURNEMIS
!         CATEGORY  = 'BURNEMIS'
!         DO N = 1, NBIOTRCE
!
!            ! Compile TRACER [ molec / cm2 / s ] 
!            ! Get the actual emission in the current cell
!            ! Original emissions are in EMS_orig, but in a variety of unit
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!                  
!               ! molec C / cm3 / s total sim  -> molec C / cm2 / s
!               TRACER_VOC(I,J,N) = BURNEMIS_orig(N,I,J)
!     &                           * BXHEIGHT(I,J,1) * 100d0
!     &                           / NEMIS_DT
!
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!               
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     1,         I0+1,
!     &                  J0+1,      1,         TRACER_VOC(:,:,N) )
!            
!            ! dkh debug
!            print*, 'max BURNEMIS= ', MAXVAL(BURNEMIS_orig(N,:,:)), N
!
!         ENDDO
!     
!      ENDIF
!
!      ! Close file
!      CLOSE( IU_RST )
!
!      IF ( LPRINT_TOTAL ) THEN
!         ! print out scaled emissions totals
!         CALL READ_USA_MASK( USA_MASK )
!
!         ! Tracer is now going to be in units of Tg X / yr / box
!         TRACER   =0d0
!         TRACER_US=0d0
!
!         IF ( NNEMS > LLPAR ) CALL ERROR_STOP('baddd','inverse_mod')
!
!         DO N = 1, NNEMS
!
!            ! Units of emission for NOx from EMISRN are different
!            ! Units of carbon emission also different.  Skip em
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!
!               ! Get the actual emission in the current cell
!
!               ! lightning NOx
!               IF ( N == IDADJEMS_ENOxli ) THEN 
!                IF ( LLIGHTNOX ) THEN
!
!                  ! molec NOx / cm2 / s  -> Tg N / yr
!                  TRACER(I,J,N) = EMS_orig_li(I,J)
!     &                          / NEMIS_DT                     ! number of emissions
!     &                          * SEC_PER_YEAR                 ! s/yr
!     &                          * GET_AREA_CM2(J)              ! cm^2
!     &                          / XNUMOL(IDTNOX)               ! molec / kg of NO2
!     &                          * TG_PER_KG                    ! Tg / kg
!     &                          * 14.d0 / 46.d0                ! g N / g NO2
!
!                  TRACER_US(I,J,N) = TRACER(I,J,N) * USA_MASK(I,J)
!                 ELSE
!                  TRACER(I,J,N) = 0d0
!                  TRACER_US(I,J,N) = 0d0
!                 ENDIF
!
!               ! soil NOx
!               ELSEIF ( N == IDADJEMS_ENOxso ) THEN
!
!                  ! Not sure what to store as "actual emission" for soil NOx
!                  !TRACER(I,J,N) = 0d0
!                  ! molec NOx / cm2 / s total  -> Tg N / yr
!                  TRACER(I,J,N) = EMS_orig(I,J,N)
!     &                          / NEMIS_DT
!     &                          * SEC_PER_YEAR                 ! s/yr
!     &                          * GET_AREA_CM2(J)              ! cm^2
!     &                          / XNUMOL(IDTNOX)               ! molec / g of NO2
!     &                          * TG_PER_KG                    ! Tg / kg
!     &                          * 14.d0 / 46.d0                ! g N / g NO2
!
!                  TRACER_US(I,J,N) = TRACER(I,J,N) * USA_MASK(I,J)
!
!
!               ! BC / OC
!               ELSEIF ( N == IDADJEMS_BCan .or. N == IDADJEMS_OCan .or.
!     &                  N == IDADJEMS_BCbb .or. N == IDADJEMS_OCbb .or.  
!     &                  N == IDADJEMS_BCbf .or. N == IDADJEMS_OCbf     )
!     &            THEN
!
!                  ! Convert from kg C / yr to Tg C / year
!                  TRACER(I,J,N)    = EMS_orig(I,J,N) * TG_PER_KG
!                  TRACER_US(I,J,N) = EMS_orig(I,J,N) * TG_PER_KG 
!     &                             * USA_MASK(I,J)
!
!               ! Anth NOx
!               ELSEIF ( N == IDADJEMS_ENOx1 .or. N == IDADJEMS_ENOx2 )
!     &            THEN
!
!                  ! Convert from kg NOx / emdt  to Tg N / year
!                  TRACER(I,J,N)    = EMS_orig(I,J,N) * TG_PER_KG
!     &                             * MIN_PER_YEAR / GET_TS_EMIS()
!     &                             * 14d0 / 46d0 
!                  TRACER_US(I,J,N) = EMS_orig(I,J,N) * TG_PER_KG
!     &                             * MIN_PER_YEAR / GET_TS_EMIS()
!     &                             * 14d0 / 46d0 
!     &                             * USA_MASK(I,J)
!
!               ! SO2
!               ELSEIF ( N == IDADJEMS_ESO2_bb .or.
!     &                  N == IDADJEMS_ESO2_bf .or.
!     &                  N == IDADJEMS_ESO2_sh )
!
!     &            THEN
!
!                  ! Convert from kg SO2 / box / s to Tg S / year
!                  TRACER(I,J,N)    = EMS_orig(I,J,N) 
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 0.5d0
!                  TRACER_US(I,J,N) = EMS_orig(I,J,N) 
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 0.5d0 * USA_MASK(I,J)
!
!               ! Volcano SO2 emissions (dkh, cklee 09/14/08) 
!               ELSEIF ( N == IDADJEMS_ESO2_ev .or.      
!     &                  N == IDADJEMS_ESO2_nv )        
!     &            THEN 
!                  ! Convert from kg SO2 / box / s total to Tg S / year
!                  TRACER(I,J,N)    = EMS_orig(I,J,N)
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 0.5d0
!     &                             / NEMIS_DT
!                  TRACER_US(I,J,N) = EMS_orig(I,J,N)
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 0.5d0 * USA_MASK(I,J)
!     &                             / NEMIS_DT
!
!               ! NH3
!               ELSEIF ( N == IDADJEMS_ENH3_an .or.
!     &                  N == IDADJEMS_ENH3_na .or.
!     &                  N == IDADJEMS_ENH3_bb .or.
!     &                  N == IDADJEMS_ENH3_bf )
!     &            THEN
!
!                  ! Convert from kg NH3 / box / s to Tg N / year
!                  TRACER(I,J,N)    = EMS_orig(I,J,N) 
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 14d0 / 17d0
!                  TRACER_US(I,J,N) = EMS_orig(I,J,N) 
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 14d0 / 17d0
!     &                             * USA_MASK(I,J)
!
!               ! Anth SOx
!               ELSEIF ( N == IDADJEMS_ESOx1   .or.
!     &                  N == IDADJEMS_ESOx2 )
!     &            THEN
!
!                  ! Convert from molec SOx / cm2 / s to Tg S / year
!                  TRACER(I,J,N)    = EMS_orig(I,J,N) * GET_AREA_CM2(J)
!     &                             / XNUMOL(IDTSO2)
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 0.5d0
!                  TRACER_US(I,J,N) = EMS_orig(I,J,N) * GET_AREA_CM2(J)
!     &                             / XNUMOL(IDTSO2)
!     &                             * SEC_PER_YEAR * TG_PER_KG
!     &                             * 0.5d0
!     &                             * USA_MASK(I,J)
!
!               ELSE
!
!                  CALL ERROR_STOP('undefined emissions',
!     &                           'inverse_mod.f')
!
!               ENDIF
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!          ENDDO
!
!          print*, 'PRIOR EMISSIONS'
!          print*, 'TOTAL SOx1      [Tg S/y] = ', SUM(TRACER(:,:,1))
!          print*, 'TOTAL SOx2      [Tg S/y] = ', SUM(TRACER(:,:,2))
!          print*, 'TOTAL SO2_sh    [Tg S/y] = ', SUM(TRACER(:,:,3))
!          print*, 'TOTAL SO2_bb    [Tg S/y] = ', SUM(TRACER(:,:,4))
!          print*, 'TOTAL SO2_bf    [Tg S/y] = ', SUM(TRACER(:,:,5))
!          print*, 'TOTAL NH3_bb    [Tg N/y] = ', SUM(TRACER(:,:,6))
!          print*, 'TOTAL NH3_bf    [Tg N/y] = ', SUM(TRACER(:,:,7))
!          print*, 'TOTAL NH3_an    [Tg N/y] = ', SUM(TRACER(:,:,8))
!          print*, 'TOTAL NH3_na    [Tg N/y] = ', SUM(TRACER(:,:,9))
!          print*, 'TOTAL BCan      [Tg C/y] = ', SUM(TRACER(:,:,10))
!          print*, 'TOTAL OCan      [Tg C/y] = ', SUM(TRACER(:,:,11))
!          print*, 'TOTAL BCbf      [Tg C/y] = ', SUM(TRACER(:,:,12))
!          print*, 'TOTAL OCbf      [Tg C/y] = ', SUM(TRACER(:,:,13))
!          print*, 'TOTAL BCbb      [Tg C/y] = ', SUM(TRACER(:,:,14))
!          print*, 'TOTAL OCbb      [Tg C/y] = ', SUM(TRACER(:,:,15))
!          print*, 'TOTAL NOx1      [Tg N/y] = ', SUM(TRACER(:,:,16))
!          print*, 'TOTAL NOx2      [Tg N/y] = ', SUM(TRACER(:,:,17))
!          print*, 'TOTAL NOx_li    [Tg N/y] = ', SUM(TRACER(:,:,18))
!          print*, 'TOTAL NOx_so    [Tg N/y] = ', SUM(TRACER(:,:,19))
!          print*, 'TOTAL SO2_ev    [Tg S/y] = ', SUM(TRACER(:,:,20))
!          print*, 'TOTAL SO2_nv    [Tg S/y] = ', SUM(TRACER(:,:,21))
!          print*, 'TOTAL US SOx1   [Tg S/y] = ', SUM(TRACER_US(:,:,1))
!          print*, 'TOTAL US SOx2   [Tg S/y] = ', SUM(TRACER_US(:,:,2))
!          print*, 'TOTAL US SO2_sh [Tg S/y] = ', SUM(TRACER_US(:,:,3))
!          print*, 'TOTAL US SO2_bb [Tg S/y] = ', SUM(TRACER_US(:,:,4))
!          print*, 'TOTAL US SO2_bf [Tg S/y] = ', SUM(TRACER_US(:,:,5))
!          print*, 'TOTAL US NH3_bb [Tg N/y] = ', SUM(TRACER_US(:,:,6))
!          print*, 'TOTAL US NH3_bf [Tg N/y] = ', SUM(TRACER_US(:,:,7))
!          print*, 'TOTAL US NH3_an [Tg N/y] = ', SUM(TRACER_US(:,:,8))
!          print*, 'TOTAL US NH3_na [Tg N/y] = ', SUM(TRACER_US(:,:,9))
!          print*, 'TOTAL US BCan   [Tg C/y] = ', SUM(TRACER_US(:,:,10))
!          print*, 'TOTAL US OCan   [Tg C/y] = ', SUM(TRACER_US(:,:,11))
!          print*, 'TOTAL US BCbf   [Tg C/y] = ', SUM(TRACER_US(:,:,12))
!          print*, 'TOTAL US OCbf   [Tg C/y] = ', SUM(TRACER_US(:,:,13))
!          print*, 'TOTAL US BCbb   [Tg C/y] = ', SUM(TRACER_US(:,:,14))
!          print*, 'TOTAL US OCbb   [Tg C/y] = ', SUM(TRACER_US(:,:,15))
!          print*, 'TOTAL US NOx1   [Tg N/y] = ', SUM(TRACER_US(:,:,16))
!          print*, 'TOTAL US NOx2   [Tg N/y] = ', SUM(TRACER_US(:,:,17))
!          print*, 'TOTAL US NOx_li [Tg N/y] = ', SUM(TRACER_US(:,:,18))
!          print*, 'TOTAL US NOx_so [Tg N/y] = ', SUM(TRACER_US(:,:,19))
!          print*, 'TOTAL US SO2_ev [Tg N/y] = ', SUM(TRACER_US(:,:,20))
!          print*, 'TOTAL US SO2_nv [Tg N/y] = ', SUM(TRACER_US(:,:,21))
!
!
!          DO N = 1, NNEMS
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J ) 
!               DO J = 1, JJPAR
!               DO I = 1, IIPAR
!#if defined ( LOG_OPT )
!               ! Apply current scaling
!               TRACER(I,J,N) = TRACER(I,J,N) * EXP(EMS_ICS(I,J,1,N))
!               TRACER_US(I,J,N) = TRACER_US(I,J,N)
!     &                          * EXP(EMS_ICS(I,J,1,N))
!
!#else
!               ! Apply current scaling
!               TRACER(I,J,N) = TRACER(I,J,N) * EMS_ICS(I,J,1,N)
!               TRACER_US(I,J,N) = TRACER_US(I,J,N) * EMS_ICS(I,J,1,N)
!#endif 
!
!               ENDDO
!               ENDDO
!!$OMP END PARALLEL DO
!
!          ENDDO
!
!          print*, 'POSTERIOR EMISSIONS'
!          print*, 'TOTAL SOx1      [Tg S/y] = ', SUM(TRACER(:,:,1))
!          print*, 'TOTAL SOx2      [Tg S/y] = ', SUM(TRACER(:,:,2))
!          print*, 'TOTAL SO2_sh    [Tg S/y] = ', SUM(TRACER(:,:,3))
!          print*, 'TOTAL SO2_bb    [Tg S/y] = ', SUM(TRACER(:,:,4))
!          print*, 'TOTAL SO2_bf    [Tg S/y] = ', SUM(TRACER(:,:,5))
!          print*, 'TOTAL NH3_bb    [Tg N/y] = ', SUM(TRACER(:,:,6))
!          print*, 'TOTAL NH3_bf    [Tg N/y] = ', SUM(TRACER(:,:,7))
!          print*, 'TOTAL NH3_an    [Tg N/y] = ', SUM(TRACER(:,:,8))
!          print*, 'TOTAL NH3_na    [Tg N/y] = ', SUM(TRACER(:,:,9))
!          print*, 'TOTAL BCan      [Tg C/y] = ', SUM(TRACER(:,:,10))
!          print*, 'TOTAL OCan      [Tg C/y] = ', SUM(TRACER(:,:,11))
!          print*, 'TOTAL BCbf      [Tg C/y] = ', SUM(TRACER(:,:,12))
!          print*, 'TOTAL OCbf      [Tg C/y] = ', SUM(TRACER(:,:,13))
!          print*, 'TOTAL BCbb      [Tg C/y] = ', SUM(TRACER(:,:,14))
!          print*, 'TOTAL OCbb      [Tg C/y] = ', SUM(TRACER(:,:,15))
!          print*, 'TOTAL NOx1      [Tg N/y] = ', SUM(TRACER(:,:,16))
!          print*, 'TOTAL NOx2      [Tg N/y] = ', SUM(TRACER(:,:,17))
!          print*, 'TOTAL NOx_li    [Tg N/y] = ', SUM(TRACER(:,:,18))
!          print*, 'TOTAL NOx_so    [Tg N/y] = ', SUM(TRACER(:,:,19))
!          print*, 'TOTAL SO2_ev    [Tg S/y] = ', SUM(TRACER(:,:,20))
!          print*, 'TOTAL SO2_nv    [Tg S/y] = ', SUM(TRACER(:,:,21))
!          print*, 'TOTAL US SOx1   [Tg S/y] = ', SUM(TRACER_US(:,:,1))
!          print*, 'TOTAL US SOx2   [Tg S/y] = ', SUM(TRACER_US(:,:,2))
!          print*, 'TOTAL US SO2_sh [Tg S/y] = ', SUM(TRACER_US(:,:,3))
!          print*, 'TOTAL US SO2_bb [Tg S/y] = ', SUM(TRACER_US(:,:,4))
!          print*, 'TOTAL US SO2_bf [Tg S/y] = ', SUM(TRACER_US(:,:,5))
!          print*, 'TOTAL US NH3_bb [Tg N/y] = ', SUM(TRACER_US(:,:,6))
!          print*, 'TOTAL US NH3_bf [Tg N/y] = ', SUM(TRACER_US(:,:,7))
!          print*, 'TOTAL US NH3_an [Tg N/y] = ', SUM(TRACER_US(:,:,8))
!          print*, 'TOTAL US NH3_na [Tg N/y] = ', SUM(TRACER_US(:,:,9))
!          print*, 'TOTAL US BCan   [Tg C/y] = ', SUM(TRACER_US(:,:,10))
!          print*, 'TOTAL US OCan   [Tg C/y] = ', SUM(TRACER_US(:,:,11))
!          print*, 'TOTAL US BCbf   [Tg C/y] = ', SUM(TRACER_US(:,:,12))
!          print*, 'TOTAL US OCbf   [Tg C/y] = ', SUM(TRACER_US(:,:,13))
!          print*, 'TOTAL US BCbb   [Tg C/y] = ', SUM(TRACER_US(:,:,14))
!          print*, 'TOTAL US OCbb   [Tg C/y] = ', SUM(TRACER_US(:,:,15))
!          print*, 'TOTAL US NOx1   [Tg N/y] = ', SUM(TRACER_US(:,:,16))
!          print*, 'TOTAL US NOx2   [Tg N/y] = ', SUM(TRACER_US(:,:,17))
!          print*, 'TOTAL US NOx_li [Tg N/y] = ', SUM(TRACER_US(:,:,18))
!          print*, 'TOTAL US NOx_so [Tg N/y] = ', SUM(TRACER_US(:,:,19))
!          print*, 'TOTAL US SO2_ev [Tg N/y] = ', SUM(TRACER_US(:,:,20))
!          print*, 'TOTAL US SO2_nv [Tg N/y] = ', SUM(TRACER_US(:,:,21))
! 
! 
!      ENDIF 
!
!
!      !### Debug
!
!      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SF_DIAG_FILE: wrote file' )
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_SF_DIAG_FILE
!
!!------------------------------------------------------------------------------
!
      SUBROUTINE READ_SF_FILE ( )
!
!******************************************************************************
!  Subroutine READ_SF_FILE reads the gctm.sf.* file into ICS_SF or EMS_SF
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Notes
!  (1 ) Add support for ACTIVE_VARS == 'EMISSIONS' case (dkh, 11/27/04)
!  (2 ) Add support for ACTIVE_VARS == 'FDTEST' case (dkh, 02/17/05)
!  (3 ) Now use CATEGORY = 'IJ-EMS-$' for ACTIVE_VARS == 'EMISSIONS' case.
!       (dkh, 03/28/05)
!  (4 ) Change name from ICS to SF, replace CMN_ADJ (dkh, ks, mak, cs  06/08/09)
!  (5 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : EMS_SF, ICS_SF
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : NNEMS, MMSCL
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,     ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD,     ONLY : NSTPL
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD,  ONLY : OPTDATA_DIR
      USE ERROR_MOD,          ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,           ONLY : IU_RST, IOERROR
      USE LOGICAL_ADJ_MOD,    ONLY : LICS, LADJ_EMS
      USE LOGICAL_ADJ_MOD,    ONLY : LADJ_STRAT
      USE LOGICAL_MOD,        ONLY : LPRT
      USE RESTART_MOD,        ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TRACER_MOD,         ONLY : N_TRACERS

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! LPRT

      ! Local Variables
      INTEGER             :: I, IOS, J, L, M, N
      INTEGER             :: NCOUNT(NNPAR)
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER(LEN=20)   :: INPUT_SF_FILE

      !=================================================================
      ! READ_SF_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_SF_FILE = 'gctm.sf.NN'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:)=0e0
      !=================================================================
      ! Open SF file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_SF_FILE )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add OPTDATA_DIR prefix to FILENAME
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      ! can hardwire this to read a specific file from another run:
      !FILENAME = TRIM( 'opt_ics/ADJv27fi04r10/gctm.ics.16' )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'S F   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_SF_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      IF ( LICS ) THEN 

         !=================================================================
         ! Read initial conditions -- store in the TRACER array
         !=================================================================
         DO N = 1, N_TRACERS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
           LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
            ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
            TRACER%r = dble(D_TRACER)
            TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:6')

            !==============================================================
            ! Assign data from the TRACER array to the xxx_IC array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-ICS-$' ) THEN

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  ICS_SF (I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

             ENDIF
         ENDDO

      ENDIF 

      IF ( LADJ_EMS ) THEN 

         !=================================================================
         ! Read emission scale factors -- store in the TRACER array
         !=================================================================
         DO N = 1, NNEMS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
           LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
            ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
            TRACER%r = dble(D_TRACER)
            TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:6')

            !==============================================================
            ! Assign data from the TRACER array to the xxx_IC array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-EMS-$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  EMS_SF(I,J,M,N) = TRACER(I,J,M)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

             ENDIF
         ENDDO

         ! Strat prod and loss (hml)
         IF ( LADJ_STRAT ) THEN

            !=================================================================
            ! Read strat prod & loss scale factors -- store in the TRACER array
            !=================================================================
            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
             LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a real I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,
     &                                      'read_strat_file:4' )

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, D_ZTAU0,D_ZTAU1,RESERVED,
     &              NI,       NJ,       NL,    IFIRST, JFIRST, LFIRST,
     &              NSKIP
             ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:5')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
               TRACER%r = dble(D_TRACER)
            TRACER%i = dimag(D_TRACER)
               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:6')

               !==============================================================
               ! Assign data from the TRACER array to the xxx_STR array.
               !==============================================================

               ! Only process observation data (i.e. aerosol and precursors)

               IF ( CATEGORY(1:8) == 'IJ-STRP$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     PROD_SF(I,J,M,N) = TRACER(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO

            !=================================================================
            ! Read strat prod & loss scale factors -- store in the TRACER array
            !=================================================================
            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
              LONRES%r = dble(D_LONRES)
            LATRES%r = dble(D_LATRES)
            LONRES%i = dimag(D_LONRES)
            LATRES%i = dimag(D_LATRES)
               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a real I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,
     &                                      'read_strat_file:4' )

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, D_ZTAU0,D_ZTAU1,RESERVED,
     &              NI,       NJ,       NL,    IFIRST, JFIRST, LFIRST,
     &              NSKIP
              ZTAU0%r = dble(D_ZTAU0)
            ZTAU1%r = dble(D_ZTAU1)
            ZTAU0%i = dimag(D_ZTAU0)
            ZTAU1%i = dimag(D_ZTAU1)
               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:5')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
             TRACER%r = dble(D_TRACER)
            TRACER%i = dimag(D_TRACER)
               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:6')

               IF ( CATEGORY(1:8) == 'IJ-STRL$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     LOSS_SF(I,J,M,N) = TRACER(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
                ENDIF
             ENDDO
         ENDIF
      ENDIF

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_SF_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_SF_FILE

!-----------------------------------------------------------------------

      SUBROUTINE MAKE_SAT_DIAG_FILE ( type)
!
!******************************************************************************
!  Subroutine MAKE_DIAG_FILE creates a binary file of a diagnostic array 
!  calculated in CALC_ADJ_FORCING in adjoint_mod.f
!  (mak, 02/09/06, 2/17/06, zhe 08/29/10)
!
!  ============================================================================
!  (1 ) MODEL_BIAS
!  NOTES:
!  (1 ) Just like MAKE_ADJ_FILE except
!       - write to .force. file
!  (2 ) Add support for ACTIVE_VARS == 'EMISSIONS' case (dkh, 11/27/04)
!  (3 ) Add support for ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (4 ) Change UNIT to unitless and change title to Scale factors (dkh, 03/06/05)
!  (5 ) Change output for ACTIVE_VARS == 'EMISSIONS' case.
!        Now use label IJ-EMS-$, and update gamap code accordingly. 
!        First write the scaling factors, in consecutive species. Temporal 
!         varations in the emissions, if any, will be in the L direction.
!        Next, write out the optimized emissions themselves.
!        Finally, write out the difference between orig and optimized emissions.
!        (dkh, 03/28/05)
!  (6 )  Use EMS_orig instead of ESO4_an_orig so that we can loop over N.
!  (7 )  Move EMS_org declaration to CMN_ADJ, (mak)
!  (8 )  Updated to v8, adj_group, 6/09/09, (mak, 6/22/09)
!  (9 )  Bug fixed, the flog SDFLAG is added, zhe 8/29/10
!  (10)  Update MOPITT obs operators (zhe, 1/19/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG,   ERROR_STOP
      USE FILE_MOD,           ONLY : IU_RST,      IOERROR
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,           ONLY : EXPAND_DATE, GET_TAU, GET_CT_EMIS
      USE ADJ_ARRAYS_MOD,     ONLY : GET_MODEL_BIAS, GET_FORCING, 
     &                               GET_MODEL,   GET_OBS, COST_ARRAY, 
     &                               COST_ARRAY,  GET_DOFS,
     &                               OBS_COUNT,   GET_EMS_ORIG, 
     &                               N_CALC, SAT, DAYS, MMSCL
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
#if defined ( MOPITT_V3_CO_OBS) || defined (MOPITT_V4_CO_OBS )  
      USE MOPITT_OBS_MOD,     ONLY : OBS_HOUR_MOPITT   !(zhe 1/19/11)
#endif
#if defined(AIRS_CO_OBS)
      USE AIRS_CO_OBS_MOD,    ONLY : OBS_HOUR_AIRS_CO
#endif
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE LOGICAL_MOD,        ONLY : LPRT
      USE LOGICAL_ADJ_MOD,    ONLY : LHMOD, LHOBS, LMODBIAS, LOBS_COUNT,
     &                               LDOFS

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! TRCOFFSET, TINDEX

      ! Arguments
      integer, intent(in)  :: type ! type of diag file
      INTEGER              :: NN

      ! Local Variables
      INTEGER              :: I, I0, IOS, J,  J0, L, M, N,H,s
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,DAYS,sat)
      TYPE (XPLEX), ALLOCATABLE  :: TRACER_EMS(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: TRACER_COST(:,:,:)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: OUTPUT_ICS_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE
      LOGICAL              :: SDFLAG


! INPUTS:
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      !=================================================================
      ! MAKE_SAT_DIAG_FILE begins here!
      !=================================================================

      SDFLAG = .FALSE.

      ! Hardwire output file for now
      IF( TYPE == 1 .AND. LHMOD ) THEN

         OUTPUT_ICS_FILE = 'gctm.model.NN'
         UNIT = 'molec/cm2'
         CATEGORY = 'IJ-AVG-$'
         SDFLAG = .TRUE.
 
      ELSEIF( TYPE == 2 .AND. LHOBS ) THEN

         OUTPUT_ICS_FILE = 'gctm.obs.NN'
         UNIT     = 'molec/cm2'
         CATEGORY = 'IJ-AVG-$'
         TITLE    = 'GEOS-CHEM observation file: '
         SDFLAG = .TRUE.

      ELSEIF( TYPE ==3 .AND. LMODBIAS ) THEN

         OUTPUT_ICS_FILE = 'gctm.modelbias.NN'
         UNIT     = '%'
         CATEGORY = 'IJ-AVG-$'
         TITLE    = 'GEOS-CHEM model bias File: ' //
     &              'model - obs bias'
         SDFLAG = .TRUE.


c$$$     IF( type == 4 ) THEN
c$$$
c$$$         OUTPUT_ICS_FILE = 'gctm.emsorig'
c$$$         TITLE    = 'GEOS-CHEM emissions file: '

      ELSEIF( TYPE == 5 .AND. LOBS_COUNT ) THEN

         OUTPUT_ICS_FILE = 'gctm.costf.NN'
         TITLE    = 'GEOS-CHEM cost file: '
         SDFLAG = .TRUE.

      ELSEIF( TYPE == 6 .AND. LDOFS ) THEN

         OUTPUT_ICS_FILE = 'gctm.dofs.NN'
         TITLE    = 'Degrees of Freedom of Signal for sats: '
         UNIT     = 'unitless'
         CATEGORY = 'IJ-DOF-$'
         SDFLAG = .TRUE.

      ! (zhe, dkh, 02/04/11) 
      ELSEIF( TYPE == 7 ) THEN

         OUTPUT_ICS_FILE = 'gctm.forcing.NN'
         TITLE    = 'Adjoint forcing: '
         UNIT     = 'unitless'
         CATEGORY = 'IJ-AVG-$'
         SDFLAG = .TRUE.

      ENDIF

      IF (SDFLAG) THEN

      ! zero TRACER array, for clarity
      TRACER(:,:,:,:) = 0d0

      ! Define variables for BINARY PUNCH FILE OUTPUT

      ! now passed in
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_ICS_FILE )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add OPT_DATA_DIR prefix to FILENAME
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_SAT_DIAG_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each observed quantity to the ics file
      !=================================================================
      !Temporarily store quantities in the TRACER array

      ! Loop over number of satellites
      DO s = 1, sat

            IF( TYPE == 1 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, DAYS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! average over all days, add all days in IDL
               TRACER(I,J,L,s) = GET_MODEL(I,J,L,s)
            ENDDO
            ENDDO
            ENDDO
           
!$OMP END PARALLEL DO

         ELSEIF( TYPE == 2 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
            DO L = 1, DAYS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! average over all days, add all days in IDL
               TRACER(I,J,L,s) = GET_OBS(I,J,L,s)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
            print*, 'obs#:', s
            print*, 'min obs:',minval(tracer(:,:,:,s))
            print*, 'max obs:',maxval(tracer(:,:,:,s))

         ELSEIF( TYPE == 3 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, DAYS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! average over all days, add all days in IDL
               TRACER(I,J,L,s) = GET_MODEL_BIAS(I,J,L,s)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ! (zhe, dkh, 02/04/11) 
         ELSEIF( TYPE == 7 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, DAYS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! average over all days, add all days in IDL
               TRACER(I,J,L,s) = GET_FORCING(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO



         ENDIF ! TYPE

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  s,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     DAYS,     I0+1,
     &                  J0+1,      1,         TRACER(:,:,:,s) )

      ENDDO   ! s = 1,SAT


      IF (TYPE .EQ. 6 ) THEN

         DO s = 1, sat

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, DAYS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! average over all days, add all days in IDL
               TRACER(I,J,L,s) = GET_DOFS(I,J,L,s)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &           HALFPOLAR, CENTER180, CATEGORY,  s,
     &           UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &           IIPAR,     JJPAR,     DAYS,     I0+1,
     &           J0+1,      1,         TRACER(:,:,:,s) )

         ENDDO                  ! s = 1,SAT
      ENDIF                     !TYPE == 6

      ! Comment for now and later decide if we want it (mak,6/22/09)
c$$$      IF( TYPE == 4 ) THEN
c$$$         ALLOCATE (TRACER_EMS(IIPAR,JJPAR,MMSCL))
c$$$         TRACER_EMS = 0e0
c$$$
c$$$        ! The following taken from ND29
c$$$         UNIT = 'kg/box/h'
c$$$         CATEGORY ='CO--SRCE'
c$$$         NN = TINDEX(29,1)
c$$$
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I, J, M)
c$$$         DO M = 1, MMSCL
c$$$         DO J = 1, JJPAR
c$$$         DO I = 1, IIPAR
c$$$            ! average over all days, add all days in IDL
c$$$            TRACER_EMS(I,J,M) = GET_EMS_ORIG(I,J,M)*EMS_ICS(I,J,M,1)
c$$$     &           /DBLE(GET_CT_EMIS() )
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO 
c$$$!$OMP END PARALLEL DO
c$$$
c$$$         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
c$$$     &        HALFPOLAR, CENTER180, CATEGORY,   NN,
c$$$     &        UNIT,      GET_TAU(), GET_TAU(), RESERVED,
c$$$     &        IIPAR,     JJPAR,     MMSCL,     I0+1,
c$$$     &        J0+1,      1,         TRACER_EMS )
c$$$         
c$$$
c$$$         IF( ALLOCATED( TRACER_EMS ) ) DEALLOCATE(TRACER_EMS)
c$$$
c$$$      ENDIF ! TYPE=4
      
      IF( TYPE == 5 ) THEN
         
         !ALLOCATE (TRACER_COST(IFDSIZE, JFDSIZE, LFDSIZE ))
         ALLOCATE (TRACER_COST(IIPAR, JJPAR, DAYS ))
         TRACER_COST = 0e0

         ! The following taken from ND29
         UNIT = 'unitless'
         CATEGORY ='COSTF'
         NN = 8301

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
         DO L = 1, DAYS !LFDSIZE
         DO J = 1, JJPAR !JFDSIZE
         DO I = 1, IIPAR !IFDSIZE
            ! COST_ARRAY
            TRACER_COST(I,J,L) = COST_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO      
!$OMP END PARALLEL DO

         !print*, 'min/max of COST_ARRAY going to file:'
         !print*, minval(TRACER_COST), maxval(TRACER_COST)

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,   NN,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,   JJPAR,   DAYS,     I0+1,
     &                  J0+1,      1,         TRACER_COST )

         PRINT*, 'FINISHED STORING COSTF'

         IF( ALLOCATED( TRACER_COST ) ) DEALLOCATE(TRACER_COST)

         ALLOCATE (TRACER_COST(IIPAR,JJPAR,1))
         TRACER_COST = 0e0

         ! The following taken from ND29
         UNIT = 'unitless'
         CATEGORY ='OBSCT'
         NN = 8401

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
         DO L = 1, 1
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ! OBS_COUNT ARRAY
            TRACER_COST(I,J,1) = OBS_COUNT(I,J)
         ENDDO
         ENDDO
         ENDDO      
!$OMP END PARALLEL DO

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,   NN,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     1,     I0+1,
     &                  J0+1,      1,         TRACER_COST )

         print*, 'finished saving OBSCT, tot obs#:',sum(obs_count)

#if defined (MOPITT_V3_CO_OBS) || defined (MOPITT_V4_CO_OBS)
         ! (MOPITT version 4, zhe 1/19/11)

         ! store OBS_HOUR, but note that it's only for the last day of sim
         CATEGORY ='OBSHR'
         NN = 8501
         TRACER_COST(:,:,:) = 0e0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
         DO L = 1, 1
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ! OBS_COUNT ARRAY
            TRACER_COST(I,J,1) = OBS_HOUR_MOPITT(I,J)
         ENDDO
         ENDDO
         ENDDO      
!$OMP END PARALLEL DO

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,   NN,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     1,     I0+1,
     &                  J0+1,      1,         TRACER_COST )

         PRINT*, 'FINISHED OBS HOUR MOPITT'
#endif

#if defined(AIRS_CO_OBS)
         ! store OBS_HOUR, but note that it's only for the last day of sim
         CATEGORY ='OBSHR'
         NN = 8502
         TRACER_COST(:,:,:) = 0e0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
         DO L = 1, 1
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ! OBS_COUNT ARRAY
            TRACER_COST(I,J,1) = OBS_HOUR_AIRS_CO(I,J)
         ENDDO
         ENDDO
         ENDDO      
!$OMP END PARALLEL DO

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,   NN,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     1,     I0+1,
     &                  J0+1,      1,         TRACER_COST )

         PRINT*, 'FINISHED storing OBS_HOUR_AIRS_CO'
#endif
 
         IF( ALLOCATED( TRACER_COST ) ) DEALLOCATE(TRACER_COST)

      ENDIF ! TYPE=5

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF( LPRT ) CALL DEBUG_MSG( '### MAKE_SAT_DIAG_FILE: wrote file')

      ENDIF !SDFLAG


      ! Return to calling program
      END SUBROUTINE MAKE_SAT_DIAG_FILE

!------------------------------------------------------------------------------
! Now move this in adj_arrays_mod.f (dkh, 10/15/09) 
!
!      SUBROUTINE EXPAND_NAME( FILENAME, N_ITRN )
!!
!!******************************************************************************
!!  Subroutine EXPAND_DATE replaces "NN" token within
!!  a filename string with the actual values. (bmy, 6/27/02, 12/2/03)
!!  (dkh, 9/22/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) FILENAME (CHARACTER) : Filename with tokens to replace
!!  (2 ) N_ITRN   (INTEGER  ) : Current iteration number
!!
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) FILENAME (CHARACTER) : Modified filename
!!
!!  NOTES:
!!  (1 ) Based on EXPAND_DATE
!!
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE CHARPAK_MOD, ONLY : STRREPL
!      USE ERROR_MOD,   ONLY : ERROR_STOP
!
!#     include "define.h"
!
!      ! Arguments
!      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME
!      INTEGER,          INTENT(IN)    :: N_ITRN
!
!      ! Local variables
!      CHARACTER(LEN=2)                :: NN_STR
!
!      !=================================================================
!      ! EXPAND_NAME begins here!
!      !=================================================================
!
!#if   defined( LINUX_PGI )
!
!      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
!      ENCODE( 2, '(i2.2)', NN_STR   ) N_ITRN
!
!#else
!
!      ! For other platforms, use an F90 internal write (bmy, 9/29/03)
!      WRITE( NN_STR,   '(i2.2)' ) N_ITRN
!
!#endif
!
!      ! Replace NN token w/ actual value
!      CALL STRREPL( FILENAME, 'NN',   NN_STR   )
!
!
!      ! Return to calling program
!      END SUBROUTINE EXPAND_NAME
!
!!-----------------------------------------------------------------------------

      SUBROUTINE DISPLAY_STUFF( LOCATION )
!
!********************************************************************************
! Subroutine DISPLAY_STUFF writes output to the screen during optimization
! (dkh, 11/28/04)
!
! NOTES
! (1 ) Rearragne the structure so that LOCATION is outermost selection, then 
!       ACTIVE_VARS == xx is subselection.  Add support for LOCATION 4 ( final
!       iteration ).  dkh, 02/17/05
! (2 ) Update to v8 and new interface/var names (mak, 6/19/09)
! (3 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!
!********************************************************************************
!
      ! References to f90 modules
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE ADJ_ARRAYS_MOD,  ONLY : COST_FUNC
      USE ADJ_ARRAYS_MOD,  ONLY : COST_FUNC_SAV
      USE ADJ_ARRAYS_MOD,  ONLY : FD_DIFF
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, MFD, LFD, NFD, EMSFD
      USE ADJ_ARRAYS_MOD,  ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : ICS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ_FD
      USE ADJ_ARRAYS_MOD,  ONLY : N_CALC, N_CALC_STOP
      USE ADJ_ARRAYS_MOD,  ONLY : ICS_SF, EMS_SF
      USE ADJ_ARRAYS_MOD,  ONLY : NNEMS, ICSFD
      USE ADJ_ARRAYS_MOD,  ONLY : STRFD
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL
      USE LOGICAL_ADJ_MOD, ONLY : LFDTEST, LADJ_EMS, LICS
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_STRAT
      USE LOGICAL_ADJ_MOD, ONLY : LFD_SPOT
      USE TRACER_MOD,      ONLY : N_TRACERS
      
#     include "CMN_SIZE"        ! Size params

      ! Argument
      INTEGER                  :: LOCATION

      ! Local variables
      INTEGER                  :: I_DUM
      INTEGER                  :: N
      TYPE (XPLEX)                   :: FINAL_ADJ_GRAD
      TYPE (XPLEX)                   :: FINAL_FD_GRAD
 
      !============================================================
      ! DISPLAY_STUFF starts here!
      !============================================================


      SELECT CASE ( LOCATION )

         ! Read/Write an iteration
         CASE( 1 )

            IF ( LICS ) THEN

               WRITE(6,*) ' ICS_SF(1,1,1,:) is ', ICS_SF(1,1,1,:)
               WRITE(6,*) ' ICS_SF(:,:,1,1) range is ',
     &                    MINVAL(ICS_SF(:,:,1,1) ),
     &                    ' to ', MAXVAL(ICS_SF(:,:,1,1) )

               WRITE(6,*) ' ICS_SF(1,1,:,1) range is ',
     &                    MINVAL(ICS_SF(1,1,:,1) ),
     &                    ' to ', MAXVAL(ICS_SF(1,1,:,1) )

            ELSEIF( LADJ_EMS ) THEN
  
               ! Nothing

            ELSEIF( LFDTEST ) THEN

               IF (LICS) THEN
                  WRITE(6,*) ' ICS_SF(FD) is ',ICS_SF(IFD,JFD,LFD,ICSFD)
               ENDIF
               IF (LADJ_EMS) THEN
                  WRITE(6,*) ' EMS_SF(FD) is ',EMS_SF(IFD,JFD,LFD,EMSFD)

                  ! Strat prod and loss (hml)
                  IF (LADJ_STRAT) THEN
                     WRITE(6,*) ' PROD_SF(FD) is '
     &                            ,PROD_SF(IFD,JFD,LFD,STRFD)
                     WRITE(6,*) ' LOSS_SF(FD) is '
     &                            ,LOSS_SF(IFD,JFD,LFD,STRFD)
                  ENDIF

               ENDIF
                  
            ELSE
              CALL ERROR_STOP( 'ACTIVE_VARS not defined!',
     &                   'DISPLAY_STUFF' )
            ENDIF 

         ! After loading gradient
         CASE( 2 )

            IF ( LICS .AND. LADJ_EMS ) THEN 
               
               WRITE(6,*) ' ICS_SF(:,:,1,1) range is ',
     &              MINVAL( ICS_SF(:,:,1,1) ), ' to ',
     &              MAXVAL( ICS_SF(:,:,1,1) )
               
               
               WRITE(6,*) ' EMS_SF(:,:,1,1) range is ',
     &              MINVAL( EMS_SF(:,:,1,1) ), ' to ',
     &              MAXVAL( EMS_SF(:,:,1,1) )
               
               print*, ' GRADNT range is',
     &              MINVAL( GRADNT ), ' to ',
     &              MAXVAL( GRADNT )

            ELSEIF ( LFDTEST .AND. LICS ) THEN

               ! for now, the I_DUM calculation is only supported for LICS, 
               ! not LADJ_EMS (mak, 6/22/09)
               I_DUM    = IFD + (  IIPAR * ( JFD - 1)  )
     &                  + (  IIPAR * JJPAR * ( LFD - 1 )  )
     &                  + (  IIPAR * JJPAR * LLPAR * ( ICSFD - 1 )  )

               WRITE(6,*) ' GRADNT(FD) = ', GRADNT(I_DUM)

               WRITE(6,*) ' MIN/MAX ICS_SF_ADJ = ', 
     &           MINVAL(ICS_SF_ADJ(:,:,:,:)), 
     &           MAXVAL(ICS_SF_ADJ(:,:,:,:))

            ELSEIF ( LFDTEST .AND. LADJ_EMS ) THEN 

                  WRITE(6,*) ' MIN/MAX EMS_SF_ADJ = ', 
     &              MINVAL(EMS_SF_ADJ(:,:,:,:)), 
     &              MAXVAL(EMS_SF_ADJ(:,:,:,:))

                  ! Strat prod and loss (hml)
                  IF ( LADJ_STRAT ) THEN

                     WRITE(6,*) ' MIN/MAX PROD_SF_ADJ = ',
     &                 MINVAL(PROD_SF_ADJ(:,:,:,:)),
     &                 MAXVAL(PROD_SF_ADJ(:,:,:,:))

                     WRITE(6,*) ' MIN/MAX LOSS_SF_ADJ = ',
     &                 MINVAL(LOSS_SF_ADJ(:,:,:,:)),
     &                 MAXVAL(LOSS_SF_ADJ(:,:,:,:))
                  ENDIF

            ELSEIF ( LICS ) THEN

!               print*, 'gradnt', gradnt(1),
!     &         gradnt(1+iipar*jjpar*llpar*(1)) ,
!     &         gradnt(1+iipar*jjpar*llpar*2),
!     &         gradnt(1+iipar*jjpar*llpar*3)

            ELSEIF ( LADJ_EMS ) THEN

               WRITE(6,*) ' EMS_SF(:,:,1,1) range is ',
     &                 MINVAL( EMS_SF(:,:,1,1) ), ' to ',
     &                 MAXVAL( EMS_SF(:,:,1,1) )

               print*, ' GRADNT range is',
     &               MINVAL( GRADNT ), ' to ',
     &               MAXVAL( GRADNT )

              ! Strat prod and loss (hml)
              IF ( LADJ_STRAT ) THEN
                  WRITE(6,*) ' PROD_SF(:,:,1,1) range is ',
     &                    MINVAL( PROD_SF(:,:,1,1) ), ' to ',
     &                    MAXVAL( PROD_SF(:,:,1,1) )

                  WRITE(6,*) ' LOSS_SF(:,:,1,1) range is ',
     &                    MINVAL( LOSS_SF(:,:,1,1) ), ' to ',
     &                    MAXVAL( LOSS_SF(:,:,1,1) )

                  print*, ' PROD_GRADNT range is',
     &                  MINVAL( GRADNT_P ), ' to ',
     &                  MAXVAL( GRADNT_P )

                  print*, ' LOSS_GRADNT range is',
     &                  MINVAL( GRADNT_L ), ' to ',
     &                  MAXVAL( GRADNT_L )
               ENDIF

            ELSE
              CALL ERROR_STOP( 'ACTIVE VARS not defined!',
     &                   'DISPLAY_STUFF, inverse_mod.f' )
            ENDIF 

            ! For all values of ACTIVE_VARS...
            WRITE(6,*) ' cost function', COST_FUNC 
            IF ( N_CALC > 1 ) THEN 
               WRITE(6,*) ' local change        = ',
     &                  COST_FUNC / COST_FUNC_SAV(N_CALC - 1),
     &                    ' = current / previous '
            ENDIF 
            WRITE(6,*) ' total change so far = ',
     &                  COST_FUNC / COST_FUNC_SAV(1),
     &                  ' = currrent / initial '


        ! Compute an iteration
         CASE( 3 )

            WRITE(6,*) ' COMPUTING NEW VALUES FOR N_CALC = ',
     &                    N_CALC

            IF( LFDTEST .AND. LICS) THEN
               
               WRITE(6,*) ' COMPUTING NEW VALUES FOR N_CALC = ',
     &                    N_CALC
               IF (LICS) THEN
                  WRITE(6,*) ' CURRENT ICS_SF(FD) IS ',
     &                 ICS_SF(IFD,JFD,LFD,ICSFD)
               ENDIF
               IF (LADJ_EMS) THEN
                  WRITE(6,*) ' CURRENT EMS_SF(FD) IS ',
     &              EMS_SF(IFD,JFD,MFD,EMSFD)
               ENDIF


            ELSEIF ( LICS ) THEN
               WRITE(6,*) ' CURRENT ICS_SF(1,1,1,:) IS ',
     &                     ICS_SF(1,1,1,:)

               WRITE(6,*) ' ICS_SF(:,:,1,1) range is ',
     &                    MINVAL(ICS_SF(:,:,1,1) ),
     &                    ' to ', MAXVAL(ICS_SF(:,:,1,1) )

               WRITE(6,*) ' ICS_SF(1,1,:,1) range is ',
     &                    MINVAL(ICS_SF(1,1,:,1) ),
     &                    ' to ', MAXVAL(ICS_SF(1,1,:,1) )

               WRITE(6,*) ' RANGE OF ICS_SF(:,:,:,:) IS ',
     &                    MINVAL(ICS_SF), ' TO ',
     &                    MAXVAL(ICS_SF)

            ELSEIF( LADJ_EMS ) THEN
               
                ! Nothing

            ELSE
              CALL ERROR_STOP( 'ACTIVE VARS not defined!',
     &                   'DISPLAY_STUFF, inverse_mod.f' )
            ENDIF 

         !After the final iteration
         CASE( 4 )

            ! For all values of ACTIVE_VARS...
            WRITE(6,*) 'COST_FUNC = ', COST_FUNC
            IF ( COST_FUNC_SAV(1) > 0d0 )
     &          WRITE(6,*) 'COST_FUNC reduction = ',
     &                      COST_FUNC / COST_FUNC_SAV(1) 

            ! Add gradient diagnostics (dkh, 06/24/09) 
            IF ( LICS ) THEN 
               DO N = 1, N_TRACERS
                  WRITE(6,*) 'MIN ICS_SF_ADJ = ', 
     &               MINVAL(ICS_SF_ADJ(:,:,:,N)), N
                  WRITE(6,*) 'MAX ICS_SF_ADJ = ', 
     &               MAXVAL(ICS_SF_ADJ(:,:,:,N)), N
               ENDDO 
            ENDIF 
            IF ( LADJ_EMS ) THEN 
               DO N = 1, NNEMS
                  WRITE(6,*) 'MIN EMS_SF_ADJ = ', 
     &               MINVAL(EMS_SF_ADJ(:,:,:,N)), N
                  WRITE(6,*) 'MAX EMS_SF_ADJ = ', 
     &               MAXVAL(EMS_SF_ADJ(:,:,:,N)), N
               ENDDO 

               ! strat prod and loss (hml)
               IF ( LADJ_STRAT ) THEN
                  DO N = 1, NSTPL
                     WRITE(6,*) 'MIN PROD_SF_ADJ = ',
     &                  MINVAL(PROD_SF_ADJ(:,:,:,N)), N
                     WRITE(6,*) 'MAX PROD_SF_ADJ = ',
     &                  MAXVAL(PROD_SF_ADJ(:,:,:,N)), N
                     WRITE(6,*) 'MIN LOSS_SF_ADJ = ',
     &                  MINVAL(LOSS_SF_ADJ(:,:,:,N)), N
                     WRITE(6,*) 'MAX LOSS_SF_ADJ = ',
     &                  MAXVAL(LOSS_SF_ADJ(:,:,:,N)), N
                  ENDDO
               ENDIF

            ENDIF 
    
            ! Compile statistics from the finite difference test.
            ! Calculate final gradients after two iterations.
            ! Now only do this for a SPOT test (dkh, 02/21/11) 
            !IF ( LFDTEST .AND. N_CALC == 2 ) THEN 
            IF ( LFD_SPOT .AND. N_CALC == 2 ) THEN 

               IF ( LADJ_EMS ) THEN
                  ! Determine the gradient calculated using the adjoint method
                  ! as an average of the gradient at FD_PERT [ STT_ADJ_FD(1) ]
                  ! and FD_PERT + FD_DIFF [ STT_ADJ_FD(2) ].
                  STT_ADJ_FD(2)  = EMS_SF_ADJ(IFD,JFD,MFD,EMSFD)
                  FINAL_ADJ_GRAD = .5d0 
     &                 * ( STT_ADJ_FD(1) + STT_ADJ_FD(2) )
               ELSEIF ( LICS ) THEN
                  ! Determine the gradient calculated using the adjoint method
                  ! as an average of the gradient at FD_PERT [ STT_ADJ_FD(1) ]
                  ! and FD_PERT + FD_DIFF [ STT_ADJ_FD(2) ].
                  STT_ADJ_FD(2)  = ICS_SF_ADJ(IFD,JFD,LFD,ICSFD)
                  FINAL_ADJ_GRAD = .5d0 
     &                 * ( STT_ADJ_FD(1) + STT_ADJ_FD(2) )
                  
               ENDIF
                   

               ! The finite difference gradient is 
               ! [ J( FD_PERT + FD_DIFF ) - J( FD_PERT ) ] / FD_DIFF
               FINAL_FD_GRAD = ( COST_FUNC - COST_FUNC_SAV(1) )
     &                       / ( FD_DIFF )

               ! Echo results to the screen
               WRITE(6,*) ' ADJOINT gradient = ', FINAL_ADJ_GRAD
               WRITE(6,*) ' FN DIFF gradient = ', FINAL_FD_GRAD
               WRITE(6,*) ' ADJ / FD = ', 
     &                    FINAL_ADJ_GRAD / FINAL_FD_GRAD

            ENDIF 
!
            WRITE(6,*) 'FORCE EXIT AFTER ', N_CALC_STOP,' ITERATIONS.'

         CASE DEFAULT
         ! Nothing
!
      END SELECT


      END SUBROUTINE DISPLAY_STUFF

! needs to be updated 
!!----------------------------------------------------------------------
!!
!!      SUBROUTINE INIT_REGIONAL_EMS
!!
!!********************************************************************************
!! Subroutine INIT_REGIONAL_EMS initializes spatially dependent emissions factors
!! (dkh, 12/04/04)
!!
!! NOTES
!! (1 ) Updated to add random noise. (dkh, 08/27/06)  
!!********************************************************************************
!!
!#     include "CMN_SIZE"   ! Size params
!
!      ! Local variables
!      INTEGER              :: I, J
!      TYPE (XPLEX)     :: RAN
!      
!      !============================================================
!      ! INIT_REGIONAL_EMS begins here!
!      !============================================================
!      WRITE(6,*) ' U S E   S P A T I A L L Y   V A R I A B L E '
!      WRITE(6,*) ' E M I S S I O N S   S C A L I N G S   F O R '
!      WRITE(6,*) ' R E F E R E N C E   C A L C U L A T I O N   '
!      
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, RAN )
!      DO I = 1, IIPAR
!      DO J = 1, JJPAR
!
!         ! Nor Am
!         IF ( I < 28 .AND. J > 28 ) THEN
!             EMS_SF(I,J,1,IDADJEMS_ESOx1) = 0.8D0
!             EMS_SF(I,J,1,IDADJEMS_ESOx2) = 0.8D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx1) = 0.85D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx2) = 0.85D0
!
!         ! Europe
!         ELSEIF ( I > 27 .AND. I < 48 .AND. J > 28 ) THEN
!             EMS_SF(I,J,1,IDADJEMS_ESOx1) = 0.7D0
!             EMS_SF(I,J,1,IDADJEMS_ESOx2) = 0.7D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx1) = 0.95D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx2) = 0.95D0
!      
!         ! Asia / India
!         ELSEIF ( I > 47 .AND. J > 20 ) THEN
!             EMS_SF(I,J,1,IDADJEMS_ESOx1) = 1.3D0
!             EMS_SF(I,J,1,IDADJEMS_ESOx2) = 1.3D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx1) = 1.2D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx2) = 1.2D0
!         
!         ! The rest of the Southern Hemisphere
!         ELSE
!             EMS_SF(I,J,1,IDADJEMS_ESOx1) = 0.75D0
!             EMS_SF(I,J,1,IDADJEMS_ESOx2) = 0.75D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx1) = 0.77D0
!             EMS_SF(I,J,1,IDADJEMS_ENOx2) = 0.77D0
!         
!         ENDIF
!         
!         RAN = DRAN(I+J)
!        
!         ! add a small bit of random variation 
!         EMS_SF(I,J,1,IDADJEMS_ESOx1) = EMS_SF(I,J,1,IDADJEMS_ESOx1)
!     &                                 + RAN / 20
!         EMS_SF(I,J,1,IDADJEMS_ESOx2) = EMS_SF(I,J,1,IDADJEMS_ESOx2)
!     &                                 + RAN / 20
!         EMS_SF(I,J,1,IDADJEMS_ENOx1) = EMS_SF(I,J,1,IDADJEMS_ENOx1) +
!     &                                 + RAN / 20
!         EMS_SF(I,J,1,IDADJEMS_ENOx2) = EMS_SF(I,J,1,IDADJEMS_ENOx2) +
!     &                                 + RAN / 20
!
!      ENDDO
!      ENDDO
!!OMP END PARALLEL DO
!      END SUBROUTINE INIT_REGIONAL_EMS
!!----------------------------------------------------------------------
      SUBROUTINE SET_SF_FORFD
!
!*****************************************************************************
!  Subroutine SET_SF_FORFD is used to initialize ICS_SF during the second 
!   iteration to the orginal value + FD_DIFF. dkh, 02/17/05
!
!  NOTES:
! (1 ) Add support for 2nd order FD calculation
! (2 ) Add support for FD_GLOB option (dkh, 10/11/08) 
! (3 ) Now initialize EMS_SF to FD_BKGRND  (dkh, 10/11/08) 
! (4 ) Change name to SET_SF_FORFD, replace CMN_ADJ, simplify the definition
!       of the FD pert (dkh, ks, mak, cs  06/07/09) 
! (5 ) Now support strat fluxes LADJ_STRAT and add flags to avoid accessing 
!       unallocated arrays (hml, dkh, 02/20/12, adj32_025) 
!*****************************************************************************
!
      ! Reference to f90 modules 
      USE ADJ_ARRAYS_MOD,  ONLY : ICS_SF, ICS_SF0
      USE ADJ_ARRAYS_MOD,  ONLY : EMS_SF, EMS_SF0
      USE ADJ_ARRAYS_MOD,  ONLY : IFD,JFD,LFD,NFD
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_SF, PROD_SF0
      USE ADJ_ARRAYS_MOD,  ONLY : LOSS_SF, LOSS_SF0
      USE ADJ_ARRAYS_MOD,  ONLY : MFD, EMSFD
      USE ADJ_ARRAYS_MOD,  ONLY : FD_DIFF
      USE ADJ_ARRAYS_MOD,  ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,  ONLY : ICSFD
      USE ADJ_ARRAYS_MOD,  ONLY : STRFD
      USE LOGICAL_ADJ_MOD, ONLY : LFD_SPOT
      USE LOGICAL_ADJ_MOD, ONLY : LFD_GLOB
      USE LOGICAL_ADJ_MOD, ONLY : LICS, LADJ_EMS
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_STRAT


#     include "CMN_SIZE"    ! Size params

      !=================================================================
      ! SET_SF_FORFD begins here!
      !=================================================================

      ICS_SF(:,:,:,:)                 = ICS_SF0(:,:,:,:)
      IF ( LADJ_EMS ) EMS_SF(:,:,:,:) = EMS_SF0(:,:,:,:)
      IF ( LADJ_STRAT ) THEN 
         PROD_SF(:,:,:,:)             = PROD_SF0(:,:,:,:)
         LOSS_SF(:,:,:,:)             = LOSS_SF0(:,:,:,:)
      ENDIF 

      ! Nudge the scaling factor value only in the FD cell 
      IF ( LFD_SPOT ) THEN 

         ! for initial conditions :
         IF ( LICS ) THEN 

            IF ( N_CALC == 2 ) THEN
               ICS_SF(IFD,JFD,LFD,ICSFD) = ICS_SF(IFD,JFD,LFD,ICSFD) 
     &                                   + FD_DIFF 
            ELSEIF ( N_CALC == 3 ) THEN
               ICS_SF(IFD,JFD,LFD,ICSFD) = ICS_SF(IFD,JFD,LFD,ICSFD) 
     &                                   - FD_DIFF 
            ENDIF

         ! for boundary conditions :
         ELSEIF ( LADJ_EMS ) THEN 
 
            ! Strat prod and loss (hml)
            IF ( .NOT. LADJ_STRAT )  THEN
               IF ( N_CALC == 2 ) THEN
                  EMS_SF(IFD,JFD,MFD,EMSFD) = EMS_SF(IFD,JFD,MFD,EMSFD)
     &                                      + FD_DIFF
               ELSEIF ( N_CALC == 3 ) THEN
                  EMS_SF(IFD,JFD,MFD,EMSFD) = EMS_SF(IFD,JFD,MFD,EMSFD)
     &                                      - FD_DIFF
               ENDIF

            ELSEIF ( LADJ_STRAT ) THEN
               IF ( N_CALC == 2 ) THEN
                 LOSS_SF(IFD,JFD,MFD,STRFD) = LOSS_SF(IFD,JFD,MFD,STRFD)
     &                                      + FD_DIFF
               ELSEIF ( N_CALC == 3 ) THEN
                 LOSS_SF(IFD,JFD,MFD,STRFD) = LOSS_SF(IFD,JFD,MFD,STRFD)
     &                                      - FD_DIFF
               ENDIF
            ENDIF

         ENDIF 

      ! Perturb thoughout model domain.
      ELSEIF ( LFD_GLOB ) THEN 

         ! for test with no transport:
         print*, 'PERTURB GLOBALLY !!!!'

         IF ( LICS ) THEN 

            IF ( N_CALC == 2 ) THEN
               ICS_SF(:,:,LFD,ICSFD) = ICS_SF(:,:,LFD,ICSFD) + FD_DIFF 
            ELSEIF ( N_CALC == 3 ) THEN
               ICS_SF(:,:,LFD,ICSFD) = ICS_SF(:,:,LFD,ICSFD) - FD_DIFF 
            ENDIF

         ELSEIF ( LADJ_EMS ) THEN 

            ! Strat prod and loss (hml)
            IF ( .NOT. LADJ_STRAT )  THEN
               IF ( N_CALC == 2 ) THEN
                 EMS_SF(:,:,MFD,EMSFD) = EMS_SF(:,:,MFD,EMSFD)
     &                                 + FD_DIFF
               ELSEIF ( N_CALC == 3 ) THEN
                 EMS_SF(:,:,MFD,EMSFD) = EMS_SF(:,:,MFD,EMSFD)
     &                                 - FD_DIFF
               ENDIF

            ELSEIF ( LADJ_STRAT ) THEN
               IF ( N_CALC == 2 ) THEN
                 LOSS_SF(:,:,MFD,STRFD) = LOSS_SF(:,:,MFD,STRFD)
     &                                  + FD_DIFF
               ELSEIF ( N_CALC == 3 ) THEN
                 LOSS_SF(:,:,MFD,STRFD) = LOSS_SF(:,:,MFD,STRFD)
     &                                  - FD_DIFF
               ENDIF
            ENDIF
         ENDIF

      ENDIF 


      ! Return to calling program
      END SUBROUTINE SET_SF_FORFD

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_CFN_FILE( )
!
!******************************************************************************
!  Subroutine MAKE_CFN_FILE creates a cfn.NN file which stores the current 
!   iteration number and cost function value. (dkh, 02/13/06)  
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Module Variable as Output:
!  ============================================================================
!  (1 ) COST_FUNC : Current cost function value
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC, COST_FUNC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR 
      USE FILE_MOD,          ONLY : IOERROR

#     include "CMN_SIZE"  

      ! Local variables 
      CHARACTER(LEN=80)  :: OUTPUT_CFN_FILE 
      CHARACTER(LEN=120) :: REMOVE_CFN_FILE_CMD
      CHARACTER(LEN=80)  :: FILENAME
      INTEGER            :: IOS

      !=================================================================
      ! MAKE_CFN_FILE begins here!
      !=================================================================

      ! Make file name
      OUTPUT_CFN_FILE = 'cfn.NN'

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_CFN_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )

      !=================================================================
      ! Open the cfn file for output 
      !=================================================================

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CFN_FILE: Writing ', a )

      ! Remove any previous cfn files for the current iteration 
      REMOVE_CFN_FILE_CMD = 'rm ' // TRIM (FILENAME) 

      CALL SYSTEM ( TRIM( REMOVE_CFN_FILE_CMD ) )


      ! Open file for input 
      OPEN( 65,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL',
     &      POSITION='APPEND' ) 

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 65, 'write_cost_func:1')

      ! Write iteration number and cost function
      WRITE( 65, *) N_CALC, COST_FUNC


      ! Return to calling program
      END SUBROUTINE MAKE_CFN_FILE
!------------------------------------------------------------------------------

      SUBROUTINE READ_CFN_FILE( )
!
!******************************************************************************
!  Subroutine READ_CFN_FILE reads the value fo the cost function at iteration 
!   NN from the cfn.NN file.  (dkh, 02/13/06)  
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Module variable as Output:
!  ============================================================================
!  (1 ) COST_FUNC : Cost function value
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : COST_FUNC
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR
      USE FILE_MOD,          ONLY : IOERROR
      USE ERROR_MOD,         ONLY : ERROR_STOP

#     include "CMN_SIZE"  

      ! Local variables 
      CHARACTER(LEN=80)  :: OUTPUT_CFN_FILE 
      CHARACTER(LEN=80)  :: FILENAME
      INTEGER            :: N, N_TMP, IOS
      TYPE (XPLEX)             :: CFN_TMP, COST_FUNC_check
      LOGICAL            :: FOUND = .FALSE. 

      !=================================================================
      ! READ_CFN_FILE begins here!
      !=================================================================

      ! Make file name
      OUTPUT_CFN_FILE = 'cfn.NN'

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_CFN_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )

      !=================================================================
      ! Open the cost function file for input
      !=================================================================

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CFN_FILE: Reading ', a )

      ! Open file for input -- readonly
      OPEN( 65,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL',
     &      POSITION='REWIND')

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 65, 'read_cost_func:1')

      ! Read values in file
      READ( 65, *) N_TMP, CFN_TMP

      ! Check to make sure that we're reading the correct file. If so, update
      ! COST_FUNC with the value from the file.    
      IF ( N_TMP == N_CALC) THEN
         COST_FUNC = CFN_TMP
         FOUND     = .TRUE.
      ENDIF 

      ! Error check
      IF ( .NOT. FOUND ) THEN
         CALL ERROR_STOP('Cost function value missing', 'inverse_mod' )
      ENDIF 
        
      ! Return to calling program
      END SUBROUTINE READ_CFN_FILE

!------------------------------------------------------------------------------

      SUBROUTINE SET_OPT_RANGE( )
!
!******************************************************************************
!  Subroutine SET_OPT_RANGE sets the range of the emissions which we 
!  wish to optimize by setting all others to zero. (dkh, 10/17/06)  
! 
!
!  Module variables as Input:
!  ============================================================================
!  (1 ) EMS_SF_ADJ       : All emissions gradients 
!  (2 ) ICS_SF_ADJ       : All tracer gradients 
!  (3 ) OPT_THIS_EMS     : Logial array of emissions to optimize
!  (4 ) OPT_THIS_ICS     : Logial array of initial conditions to optimize
!     
!  Module variables as Output:
!  ============================================================================
!  (1 ) EMS_SF_ADJ       : All emissions gradients 
!  (2 ) ICS_SF_ADJ       : All tracer gradients 
!     
!  NOTES:
! (1 ) Replace CMN_ADJ, update naming, add spatial filter from ks
!       (dkh, ks, mak, cs  06/07/09) 
! (2 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : ICS_SF_ADJ,  EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : NNEMS
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_EMS, OPT_THIS_TRACER
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_PROD
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_LOSS
      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_STRAT
      USE LOGICAL_ADJ_MOD, ONLY : LICS, LADJ_EMS
      USE TRACER_MOD,      ONLY : N_TRACERS

#     include "CMN_SIZE"  ! Size params

      ! Local variables 
      INTEGER I, J, M, N

      !=================================================================
      ! SET_OPT_RANGE begins here!
      !=================================================================

         ! dkh debug
         print*, ' SET_OPT_RANGE: MIN / MAX ICS_SF_ADJ = ',
     &      MINVAL(ICS_SF_ADJ), MAXVAL(ICS_SF_ADJ)

      IF ( LICS ) THEN 
         DO N = 1, N_TRACERS
            IF ( .not.  OPT_THIS_TRACER(N) ) THEN
               ICS_SF_ADJ(:,:,:,N) = 0d0
            ENDIF
         ENDDO
      ENDIF 
  
         ! dkh debug
         print*, ' SET_OPT_RANGE 2 : MIN / MAX ICS_SF_ADJ = ',
     &      MINVAL(ICS_SF_ADJ), MAXVAL(ICS_SF_ADJ)

      ! Zero the gradients of the species we don't want to optimize 
      IF ( LADJ_EMS ) THEN 
         DO N = 1, NNEMS
            IF ( .not. OPT_THIS_EMS(N) ) THEN 
               EMS_SF_ADJ(:,:,:,N) = 0d0
            ENDIF
         ENDDO

         ! Strat prod and loss (hml)
         IF ( LADJ_STRAT ) THEN
            DO N = 1, NSTPL
               IF ( .not. OPT_THIS_PROD(N) ) THEN
                  PROD_SF_ADJ(:,:,:,N) = 0d0
               ENDIF
               IF ( .not. OPT_THIS_LOSS(N) ) THEN
                  LOSS_SF_ADJ(:,:,:,N) = 0d0
               ENDIF
            ENDDO
         ENDIF

      ENDIF 


!      ! Only consider gradients in specific spatial range  
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M, N )    
!      DO N = 1, NNEMS 
!      DO M = 1, MMSCL
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!    
!         ! Zero the gradients which we don't 
!         ! want to optimize 
!!         IF ( (( I < 42 .and. I > 34 ) .and.   ! IN 
!!     &         ( J > 32 .and. J < 39 ))        ! EUROPE 
!!         IF ( (( I > 18 .and. I < 23 ) .and.   ! IN 
!!     &         ( J > 30 .and. J < 35 ))        ! Eastern US
!!         IF ( (( I < 19 .or. I > 22 ) .or.     ! not IN 
!!     &         ( J < 31 .or. J > 34 ))          ! Eastern US
!         IF ( (( I < 12 .or. I > 22 ) .or.     ! not IN 
!     &         ( J < 31 .or. J > 34 ))         ! US
!     &         .or. 
!     &        ( .not. OPT_THIS_EMS(N) )  )   THEN   
!
!            EMS_SF_ADJ(I,J,M,N) = 0d0 
!
!         ENDIF 
!
!      ENDDO 
!      ENDDO 
!      ENDDO 
!      ENDDO 
!!$OMP END PARALLEL DO
!

! old code from ks 
!#if defined ( TES_O3_OBS ) 
!
!      ! Zero the gradients above NLEVS 
!      ICS_FD_ADJ(:,:,NLEVS+1:LLPAR,:)                = 0d0 
! 
!      ! Smoothly drive gradients to zero at poles
!      IF (NLAT_TO_IGNORE > 0) THEN 
! 
!         DO N = 1, N_TRACERS
!         DO L = 1, NLEVS
!         DO J = 1, NLAT_TO_IGNORE
!         DO I = 1, IIPAR
!            TEMP                =  NLAT_TO_IGNORE - J
!            FACTOR              = COS( (TEMP / ( NLAT_TO_IGNORE - 1 ) ) 
!     &                          * ( pi / 2 )**2
!            ICS_FD_ADJ(I,J,L,N) = ICS_FD_ADJ(I,J,L,N) * FACTOR 
!         ENDDO
!         ENDDO
!         ENDDO
!         ENDDO
!
!         DO N = 1, N_TRACERS
!         DO L = 1, NLEVS
!         DO J = JJPAR - NLAT_TO_IGNORE + 1, JJPAR
!         DO I = 1,IIPAR
!            TEMP                =  NLAT_TO_IGNORE - J
!            FACTOR              = COS( (TEMP / ( NLAT_TO_IGNORE - 1 ) ) 
!     &                          * ( pi / 2 )**2
!            ICS_FD_ADJ(I,J,L,N) = ICS_FD_ADJ(I,J,L,N) * FACTOR 
!         ENDDO
!         ENDDO
!         ENDDO
!         ENDDO
!
!      ENDIF 
!#endif 

      ! Return to calling program
      END SUBROUTINE SET_OPT_RANGE

!------------------------------------------------------------------------------
      TYPE (XPLEX) FUNCTION DRAN(K)
C
C     RANDOM NUMBER GENERATOR - BASED ON ALGORITHM 266 BY PIKE AND
C      HILL (MODIFIED BY HANSSON), COMMUNICATIONS OF THE ACM,
C      VOL. 8, NO. 10, OCTOBER 1965.
C
C     THE SINGLE PRECISION VERSION OF THIS SUBPROGRAM IS INTENDED
C     FOR USE ON COMPUTERS WITH FIXED POINT WORDLENGTH OF AT
C     LEAST 29 BITS.  IT IS BEST IF THE FLOATING POINT
C     SIGNIFICAND HAS AT MOST 29 BITS.
C
C     FOLLOWING CODY AND WAITE'S RECOMMENDATION (P .14), WE
C     PRODUCE A PAIR OF RANDOM NUMBERS AND USE RAN1 +
C     2**(-29)*RAN2 IN AN ATTEMPT TO GENERATE ABOUT 58 RANDOM BITS.
C
      INTEGER IY,J,K
      DATA             IY     /100001/
C
      J = K
      IY = IY * 125
      IY = IY - (IY/2796203) * 2796203
      DRAN = XPLX(XPLX(IY)) / 2796203.0D+00
C
      IY = IY * 125
      IY = IY - (IY/2796203) * 2796203
      DRAN = DRAN + (XPLX(XPLX(IY))/2796203.0D+00)/536870912.0D+00
      RETURN
C     ---------- LAST CARD OF DRAN ----------
      END FUNCTION DRAN
! needs to be updated:
!--------------------------------------------------------------------------------
!
!      SUBROUTINE UPDATE_HESSIAN( )
!
!******************************************************************************
!  Subroutine UPDATE_HESSIAN constructs an approximation of the inverse
!  Hessian using the DFP formula (see Muller and Stavrakou, 2005, eqn 18).
!
!  This routine is set up to be used offline so that the Hessian is 
!  only approximated at the end of a convered optimization. To implement, 
!  uncomment code in 3 places in inverse.f
!
!  The initial estimate can be identiy matrix or initial estimate of uncertainty
!
!  It takes too long to consider all possible correlations, so we apply the 
!  following filters:
!   - Only consider corelations between emissions of 
!     - anth SOx (surface and stack)
!     - anth NOx (surface and stack)
!     - anth NH3
!     - natural NH3
!   - Only within the U.S. 
!   - Only in places where ADJ_EMS at first iteration is > 1d-4
!
!   If these filters are changes, the array diminsion HMAX will need to be
!   updated. To determine the size of the MASD parameter, do a dry run,
!   then go back and update. 
!
!  NOTES:
!
!******************************************************************************
!
!
!      print*, ' SUBROUTINE UPDATE_HESSIAN needs to be updated '
!      print*, ' SUBROUTINE UPDATE_HESSIAN needs to be updated '
!      print*, ' SUBROUTINE UPDATE_HESSIAN needs to be updated '
!      print*, ' SUBROUTINE UPDATE_HESSIAN needs to be updated '
!
!      ! Reference to f90 modules
!
!#     include "CMN_SIZE"  
!
!      ! Arguments
!    
!      ! Local variables 
!      INTEGER, PARAMETER :: HMAX = 3675
!
!      INTEGER       :: I, J, M, N, II, JJ, NITR
!
!      TYPE (XPLEX), SAVE  :: USA_MASK(IIPAR,JJPAR)
!
!      INTEGER, SAVE :: IIMAP(IIPAR,JJPAR,MMSCL,NNEMS) = 0d0
!      INTEGER, SAVE :: MAPI(HMAX), MAPJ(HMAX)
!      INTEGER, SAVE :: MAPM(HMAX), MAPN(HMAX)
!
!      TYPE (XPLEX), SAVE  :: EMS_SF_OLD(IIPAR,JJPAR,MMSCL,NNEMS) 
!      TYPE (XPLEX), SAVE  :: ADJ_EMS_OLD(IIPAR,JJPAR,MMSCL,NNEMS)
!      TYPE (XPLEX), SAVE  :: HINV(HMAX,HMAX)
!      LOGICAL, SAVE  :: FIRST = .TRUE. 
!
!      TYPE (XPLEX)        :: S(HMAX),         Y(HMAX),     YTS,  YTHINVY
!      TYPE (XPLEX)        :: YTS_INV,  YTHINVY_INV
!      TYPE (XPLEX)        :: SST(HMAX,HMAX), HINVY(HMAX), YTHINV(HMAX)
!      TYPE (XPLEX)        :: HINVYYTHINV(HMAX,HMAX)
!
!      !=================================================================
!      ! UPDATE_HESSIAN begins here!
!      !=================================================================
!
!      PRINT*, ' UPDATE HESSIAN AT ITERATE ', N_CALC
!
! 
!      IF ( FIRST ) THEN 
!        
!         ! Initialize HINV to the identity matrix (or initial unc. est)
!         HINV(:,:)   = 0d0 
!
!         DO JJ = 1, HMAX 
!         DO II = 1, HMAX 
! 
!           IF ( II == JJ ) HINV(II,II) = 0.3d0 
!   
!         ENDDO
!         ENDDO
!           
!         ! Get USA mask
!         CALL READ_USA_MASK( USA_MASK )
!
!         ! dkh debug
!         print*, ' yea yea eya'
!
!         II = 0 
! 
!         DO N = 1, NNEMS
!         DO M = 1, MMSCL
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Spatial filter
!            ! Only in US:
!!            IF ( USA_MASK(I,J) == 0d0 ) CYCLE 
!            ! Only in places where emissions are nonzero 
!            IF ( ABS(ADJ_EMS(I,J,M,N)) < 1d-4 ) CYCLE 
!
!            IF ( 
!     &        N == IDADJEMS_ESOx1 .or. 
!     &        N == IDADJEMS_ESOx2 .or. 
!     &        N == IDADJEMS_ENOx1 .or. 
!     &        N == IDADJEMS_ENOx2 .or. 
!     &        N == IDADJEMS_ENH3_an  .or. 
!     &        N == IDADJEMS_ENH3_na  
!     &                                 ) THEN 
!   
!
!               ! Update vector index
!               II = II + 1 
!
!               ! Save mapping arrays
!               IIMAP(I,J,M,N) = II          
!               MAPI(II) = I
!               MAPJ(II) = J
!               MAPM(II) = M
!               MAPN(II) = N
!
!            ENDIF 
!
!         ENDDO 
!         ENDDO 
!         ENDDO 
!         ENDDO 
!
!      
!         EMS_ICS_OLD(:,:,:,:) = EMS_ICS(:,:,:,:)
!         ADJ_EMS_OLD(:,:,:,:) = ADJ_EMS(:,:,:,:)
!         print*, ' UPDATE HESSIAN, pts founds = ', II
!         CALL MAKE_HESS_FILE( HINV, USA_MASK , HMAX, IIMAP , 1 )
!         FIRST = .FALSE. 
!
!
!         print*, 'EMS_ICS  = ', EMS_ICS(19,33,1,IDADJEMS_ESOx2)
!         print*, 'EMS_ICS_OLD  = ', EMS_ICS_OLD(19,33,1,IDADJEMS_ESOx2)
!         print*, 'ADJ_EMS  = ', ADJ_EMS(19,33,1,IDADJEMS_ESOx2)
!         print*, 'ADJ_EMS_OLD  = ', ADJ_EMS_OLD(19,33,1,IDADJEMS_ESOx2)
!   
!         RETURN 
!      ENDIF 
!
! 
!      DO II = 1, HMAX 
! 
!         I = MAPI(II)
!         J = MAPJ(II)
!         M = MAPM(II)
!         N = MAPN(II)
!
!         ! find s_k = f_{k+1} - f_{k}
!         S(II) = EMS_ICS(I,J,M,N) - EMS_ICS_OLD(I,J,M,N)
! 
!         ! find y_k = grad_{k+1} - grad_{k}
!         Y(II) = ADJ_EMS(I,J,M,N) - ADJ_EMS_OLD(I,J,M,N)
!
!      ENDDO 
!
!      print*, ' UPDATE HESSIAN, pts founds = ', II
!
!      print*, 'EMS_ICS  = ', EMS_ICS(19,33,1,IDADJEMS_ESOx2)
!      print*, 'EMS_ICS_OLD  = ', EMS_ICS_OLD(19,33,1,IDADJEMS_ESOx2)
!      print*, 'ADJ_EMS  = ', ADJ_EMS(19,33,1,IDADJEMS_ESOx2)
!      print*, 'ADJ_EMS_OLD  = ', ADJ_EMS_OLD(19,33,1,IDADJEMS_ESOx2)
!
!      ! Rotate
!      EMS_ICS_OLD(:,:,:,:) = EMS_ICS(:,:,:,:)
!      ADJ_EMS_OLD(:,:,:,:) = ADJ_EMS(:,:,:,:)
!
!      !----------------------------------------------------------
!      ! Update inverse Hessian 
!      !----------------------------------------------------------
!
!      ! y^T*s 
!      YTS = 0d0 
!      DO II = 1, HMAX
!
!         YTS = YTS + Y(II) * S(II)
!
!      ENDDO 
!   
!      print*, ' YTS = ', YTS , N_CALC
! 
!      ! s * s^T / YTS 
!      DO II = 1, HMAX
!      DO JJ = 1, HMAX
!
!         SST(II,JJ) = S(II) * S(JJ) 
!
!      ENDDO
!      ENDDO
!
!      ! HINV * y
!      DO II = 1, HMAX 
! 
!         HINVY(II) = 0D0 
!        
!         DO JJ = 1, HMAX 
! 
!            HINVY(II) = HINVY(II) + HINV(II,JJ) * Y(JJ)
!
!         ENDDO 
!      ENDDO
!
!      ! y^T * HINV 
!      DO JJ = 1, HMAX 
!
!         YTHINV(JJ) = 0d0 
! 
!         DO II = 1, HMAX
! 
!            YTHINV(JJ) = YTHINV(JJ) + Y(II) * HINV(II,JJ)
!
!         ENDDO
!      ENDDO
!
!
!      ! HINVY * YTHINV 
!      DO JJ = 1, HMAX 
!      DO II = 1, HMAX 
! 
!            HINVYYTHINV(II,JJ) = HINVY(II) * YTHINV(JJ)
!
!      ENDDO 
!      ENDDO 
!
! 
!      ! YT * HINVY
!      YTHINVY = 0d0 
!      DO II = 1, HMAX
!         YTHINVY = YTHINVY + Y(II) * HINVY(II)
!      ENDDO 
!      print*, 'YTHINVY = ', YTHINVY  
!
!      ! HINV = HINV + SST * (1/YTS) - HINVYYTHINV * (1/YTHINVY) 
!      YTS_INV     = 1 / YTS 
!      YTHINVY_INV = 1 / YTHINVY
!      DO JJ = 1, HMAX
!      DO II = 1, HMAX
!         
!         HINV(II,JJ) = HINV(II,JJ) 
!     &               + SST(II,JJ)         * YTS_INV
!     &               - HINVYYTHINV(II,JJ) * YTHINVY_INV
!
!      ENDDO      
!      ENDDO      
!
!      print*, ' MAX HINV = ', MAXVAL(HINV(:,:))
!      print*, ' MIN HINV = ', MINVAL(HINV(:,:))
!
!      NITR = N_CALC 
!
!      CALL MAKE_HESS_FILE( HINV, USA_MASK , HMAX, IIMAP , NITR )
!
!      ! Return to calling program
!      END SUBROUTINE UPDATE_HESSIAN
!!------------------------------------------------------------------------------
! needs to be updated
!
!      SUBROUTINE MAKE_HESS_FILE( HINV, USA_MASK , HMAX, IIMAP, NITR )
!!
!!******************************************************************************
!!  Subroutine MAKE_HESS_FILE creates a binary file of selected elements 
!!  of the approximate inverse hessian. (dkh, 05/15/07) 
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) HINV      : Current estimate of inverse hessian 
!!
!!  Module Variable as Input:
!!  ============================================================================
!!  (1 ) N_CALC    : Current iteration number
!!
!!  NOTES:
!!  (1 ) Just like MAKE_GDT_FILE except
!!        - pass NITR as an argument
!!******************************************************************************
!!
!
!      ! References to F90 modules
!      USE BPCH2_MOD
!      USE ERROR_MOD, ONLY : DEBUG_MSG, ERROR_STOP
!      USE FILE_MOD,  ONLY : IU_RST,      IOERROR
!      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,  ONLY : EXPAND_DATE, GET_TAU
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "CMN_SETUP"  !
!#     include "CMN"        ! LPRT
!#     include "CMN_ADJ"    ! NADJ, OPTDATA_DIR, ACTIVE_VARS
!
! 
!      ! Arguments
!      INTEGER       :: HMAX
!      TYPE (XPLEX)        :: HINV(HMAX,HMAX)
!      TYPE (XPLEX)        :: USA_MASK(IIPAR,JJPAR)
!      INTEGER       :: IIMAP(IIPAR,JJPAR,MMSCL,NNEMS)
!      INTEGER       :: NITR
!
!      ! Local Variables
!      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N, II, JJ
!      INTEGER              :: YYYY, MM, DD,  HH, SS
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: EMS_3D(IIPAR,JJPAR,MMSCL)
!      CHARACTER(LEN=255)   :: FILENAME
!     
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER, PARAMETER   :: HALFPOLAR = 1
!      INTEGER, PARAMETER   :: CENTER180 = 1
!
!      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE
!
!      !=================================================================
!      ! MAKE_HESS_FILE begins here!
!      !=================================================================
!
!      ! Clear intermediate arrays
!      EMS_3D(:,:,:) = 0d0
!
!      ! Hardwire output file for now
!#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
!      OUTPUT_GDT_FILE = 'gctm.invhess.NN'
!#else
!      OUTPUT_GDT_FILE = 'gctm.invhess.NN'
!#endif
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM Adjoint File: ' //
!     &           'Inverse hessian  '
!      UNIT     = 'none'
!      CATEGORY = 'IJ-INVH-'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the adjoint file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output observation file name into a local variable
!      FILENAME = TRIM( OUTPUT_GDT_FILE )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, NITR )
!
!      ! Add the OPTDATA_DIR prefix to the file name
!      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_HESS_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!!      IF ( ACTIVE_VARS == 'TRACERS'.OR. 
!!     &     ACTIVE_VARS == 'FDTEST'       ) THEN
!      IF ( ACTIVE_VARS == 'TRACERS' ) THEN 
!  
!         CALL ERROR_STOP( 'inverse hessian not supported ', 
!     &                    ' MAKE_HESS_FILE, inverse_mod.f')
! 
!      ELSEIF ( ACTIVE_VARS == 'EMISSIONS' .OR.
!     &         ACTIVE_VARS == 'FDTEST'          ) THEN
!
!         ! Reset CATEGORY as labeling in gamap is different 
!         CATEGORY = 'IJ-INVH-'
!
!         !=================================================================
!         ! Write each observed quantity to the observation file
!         !=================================================================
!         DO N = 1, NNEMS
!
!            !Temporarily store quantities in the TRACER array
!            EMS_3D(I,J,M) = 0d0
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M, II )
!            DO M = 1, MMSCL
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!
!
!               II = IIMAP(I,J,M,N)
!               IF ( II == 0 ) CYCLE 
!
!                  IF ( HINV(II,II) > 0 )  THEN 
!                     EMS_3D(I,J,M) = REAL(SQRT(HINV(II,II)))
!                  ELSE 
!                     print*, I, J, M, N, II 
!                     CALL ERROR_STOP('non positive hessian diagonal ', 
!     &                               'inverse_mod.f')
!                  ENDIF 
!               
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &                  J0+1,      1,         EMS_3D )
!
!         ENDDO
!
!         ! Reset CATEGORY as labeling in gamap is different 
!         CATEGORY = 'IJ-COREL'
!
!         !=================================================================
!         ! Write correlation for a given cell
!         !=================================================================
!         DO N = 1, NNEMS
!
!            ! target cell
!            JJ = IIMAP(13,33,1,IDADJEMS_ENH3_an)
!
!            !Temporarily store quantities in the TRACER array
!            EMS_3D(I,J,M) = 0d0
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M, II )
!            DO M = 1, MMSCL
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!
!
!               II = IIMAP(I,J,M,N)
!               !IF ( II == 0 ) CYCLE 
!               IF ( II == 0 ) THEN 
!                  EMS_3D(I,J,M) = 0d0
!               ELSE
!                  EMS_3D(I,J,M) = REAL(HINV(II,JJ)/(SQRT(HINV(II,II))
!     &                          * SQRT(HINV(JJ,JJ))))
!               ENDIF
!
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &                  J0+1,      1,         EMS_3D )
!
!         ENDDO
!      ELSE
!         CALL ERROR_STOP( 'ACTIVE_VARS not defined!',
!     &                    'MAKE_HESS_FILE' )
!      ENDIF
!
!      ! Close file
!      CLOSE( IU_RST )
!
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_HESS_FILE: wrote file' )
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_HESS_FILE

!------------------------------------------------------------------------------
!      SUBROUTINE READ_USA_MASK( USA_MASK )
!!
!!******************************************************************************
!!  Subroutine READ_USA_MASK reads the USA mask from disk.   The USA mask is
!!  the fraction of the grid box (I,J) which lies w/in the continental USA.
!!  (rch, bmy, 11/10/04, 10/3/05)
!!
!!  NOTES:
!!  (1 ) Now can read data for GEOS and GCAP grids (bmy, 8/16/05)
!!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      !USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
!
!#include "CMN_SIZE"
!
!      ! Local variables
!      TYPE (XPLEX)             :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)             :: XTAU
!      TYPE (XPLEX)             :: USA_MASK(IGLOB,JGLOB)
!      CHARACTER(LEN=255) :: FILENAME
!
!      
!      !=================================================================
!      ! READ_USA_MASK begins here!
!      !=================================================================
!      
!      ! File name
!      ! Argg - haven't initialized the forward model yet, so DATA_DIR undefined
!      ! Just put the mask in the home directory
!!      FILENAME = TRIM( DATA_DIR )           //
!!     &           'EPA_NEI_200411/usa_mask.' // GET_NAME_EXT_2D() //
!      FILENAME = 
!     &           'usa_mask.geos' //
!     &           '.'                        // GET_RES_EXT()
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_USA_MASK: Reading ', a )
!      
!      ! Get TAU0 for Jan 1985
!      XTAU  = GET_TAU0( 1, 1, 1985 ) 
!
!      ! USA mask is stored in the bpch file as #2
!      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
!     &                 XTAU,      IGLOB,    JGLOB,
!     &                 1,         ARRAY,    QUIET=.TRUE. )
!
!      ! Cast to TYPE (XPLEX)
!      !CALL TRANSFER_2D( ARRAY(:,:,1), USA_MASK )
!      USA_MASK(:,:) = ARRAY(:,:,1)
!
!      ! Return to calling program
!      END SUBROUTINE READ_USA_MASK
!
!!------------------------------------------------------------------------------
      SUBROUTINE CALC_NOPT

!
!******************************************************************************
!  Subroutine CALC_NOPT calculates the number of paramteres to optimize 
!
!  NOTES:
!  (1 ) Set NOPT for initial conditions to 3D: IIPAR*JJPAR*LLPAR*N_TRACERS to
!       be consistent with other parts of the code (mak, 6/18/09)
!  (2 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!     
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : NOPT 
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS, LICS
      USE TRACER_MOD,        ONLY : N_TRACERS, ITS_A_TAGCO_SIM

#     include "CMN_SIZE"

      !=================================================================
      ! CALC_NOPT begins here!
      !=================================================================

      ! if optimizing both initial emissions and initial conditions
      IF ( LADJ_EMS .AND. LICS ) THEN
         NOPT = IIPAR * JJPAR * MMSCL * NNEMS +
     &          IIPAR * JJPAR *  LLPAR * N_TRACERS

      ! if optimizing emissions only
      ELSEIF ( LADJ_EMS ) THEN

         NOPT = IIPAR * JJPAR * MMSCL * NNEMS

         IF ( ITS_A_TAGCO_SIM() .AND. NNEMS == 2 ) THEN
            NOPT = IIPAR * JJPAR * MMSCL + 1
         ENDIF

         ! Strat prod and loss (hml)
         IF ( LADJ_STRAT ) THEN
            NOPT = NOPT + IIPAR * JJPAR * MMSCL * NSTPL * 2
         ENDIF

      ! if optimizing initial conditions only
      ELSEIF ( LICS ) THEN

         NOPT = IIPAR * JJPAR * LLPAR * N_TRACERS

      ENDIF

      PRINT*, 'Max size of control vector is:', NOPT

      ! Return to calling program
      END SUBROUTINE CALC_NOPT

!------------------------------------------------------------------------------

      SUBROUTINE ITER_CONDITION( IT )
!
!******************************************************************************
!  Subroutine ITER_CONDITION output information which will be used 
!  to determine whether the convergence has been reached (zhe 11/28/10) 
!
!  Variable as Input:
!  ============================================================================
!  (1 ) IT        : Current iteration number
!     
! NOTES:
! (1 ) Place output in DIAGADJ_DIR instead of OPTDATA_DIR (dkh, 02/04/11) 
! (2 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : COST_FUNC_SAV
      USE DIRECTORY_ADJ_MOD,   ONLY : DIAGADJ_DIR
      USE LOGICAL_ADJ_MOD,     ONLY : LADJ_STRAT
      USE LOGICAL_ADJ_MOD,     ONLY : LATF


#     include "CMN_SIZE"            ! Size parameters

      ! Arguments 
      INTEGER                      :: IT


      ! Local variables 
      INTEGER                      :: I
      TYPE (XPLEX)                       :: PG, NG, PS, NS
      CHARACTER(LEN=255)           :: FILENAME
      LOGICAL, SAVE                :: FIRST = .TRUE.

      ! For strat prod and loss (hml)
      TYPE (XPLEX)                       :: PG_P, PG_L, NG_P, NG_L
      TYPE (XPLEX)                       :: PS_P, PS_L, NS_P, NS_L

      !=================================================================
      ! ITER_CONDITION begins here!
      !=================================================================

      PG = 0.0
      NG = 0.0
      PS = 0.0
      NS = 0.0

      ! For strat prod and loss (hml)
      PG_P = 0.0
      NG_p = 0.0
      PS_P = 0.0
      NS_P = 0.0
      PG_L = 0.0
      NG_L = 0.0
      PS_L = 0.0
      NS_L = 0.0

      FILENAME = 'gctm.iteration'
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )

      IF ( FIRST ) THEN
         OPEN (99, FILE = FILENAME, STATUS ='REPLACE')
         WRITE(99, 1001)
         WRITE(99, 1002)
         FIRST = .FALSE.
      ENDIF

      ! For strat prod and loss (hml)
      IF ( LADJ_STRAT ) THEN

         DO I = 1, IIPAR * JJPAR
            IF ( GRADNT_P(I) .GT. 0 .AND. GRADNT_L(I) .GT. 0 ) THEN
               PG_P = PG_P + GRADNT_P(I)
               PG_L = PG_L + GRADNT_L(I)
            ELSE
               NG_P = NG_P + GRADNT_P(I)
               NG_L = NG_L + GRADNT_L(I)
            ENDIF

            IF ( XP(I) .GT. 1 .AND. XL(I) .GT. 1 ) THEN
               PS_P = PS_P + XP(I) - 1
               PS_L = PS_L + XL(I) - 1
            ELSE
               NS_P = NS_P + XP(I) - 1
               NS_L = NS_L + XL(I) - 1
            ENDIF
         ENDDO

         WRITE(99, 1005) IT, LATF, COST_FUNC_SAV(IT),
     &      COST_FUNC_SAV(IT)/COST_FUNC_SAV(1), PG_P, PG_L,
     &      NG_P, NG_L, PS_P, PS_L, NS_P, NS_L

      ELSE

         DO I = 1, IIPAR * JJPAR        
            IF ( GRADNT(I) .GT. 0 ) THEN
               PG = PG + GRADNT(I)
            ELSE
               NG = NG + GRADNT(I)
            ENDIF
         
            IF ( X(I) .GT. 1 ) THEN
               PS = PS + X(I) - 1
            ELSE
               NS = NS + X(I) - 1
            ENDIF
         ENDDO

         WRITE(99, 1003) IT, LATF, COST_FUNC_SAV(IT), 
     &      COST_FUNC_SAV(IT)/COST_FUNC_SAV(1), PG, NG, PS, NS
      
      ENDIF 

 1001 format ('GEOS-CHEM ADJOINT CONVERGNECE CONDITION',/,/,
     + 'IT    = iteration number',/,
     + 'A     = accepted iteration',/,
     + 'F     = cost fun',/,
     + 'FdF0  = cost fun reduction',/,
     + 'PG    = total positive gradient',/,
     + 'NG    = total negative gradient',/,
     + 'PS    = total underestimated scaling factor',/,
     + 'NS    = total overestimated scaling factor',/)
 
 1002 format (/,3x,'IT',2x,'A',7x,'F',10x,'FdF0',9x,'PG',12x,'NG',
     +        10x,'PS',10x,'NS')
 1003 format (3x,i2,2x,L1,2x,E12.6,2x,F8.6,2x,E11.5,2x,
     +        E12.5,2x,F9.2,2x,F10.2)

! Strat prod and loss (hml)
 1004 format (/,3x,'IT',2x,'A',7x,'F',10x,'FdF0',9x,'PG_P',12x,'PG_L',
     +        12x,'NG_P',10x,'NG_L',10x,'PS_P',10x,'PS_L',10x,'NS_P',
     +        10x,'NS_L')
 1005 format (3x,i2,2x,L1,2x,E12.6,2x,F8.6,2x,E11.5,2x,E11.5,2x,
     +        E12.5,2x,E12.5,2x,F9.2,2x,F9.2,2x,F10.2,2x,F10.2)

      ! Return to calling program
      END SUBROUTINE ITER_CONDITION

!--------------------------------------------------------------------------------

      SUBROUTINE MAYBE_DO_GEOS_CHEM_ADJ( )
!
!******************************************************************************
!  Subroutine MAYBE_DO_GEOS_CHEM_ADJ is called for FDTESTS and determines
!  whether or not the adjoint model needs to be run. (dkh, 02/21/11) 
!
!  Module variables as Input:
!  ============================================================================
!  (1 ) LFD_GLOB       (LOGICAL) : Switch to perform global finite diff test
!  (2 ) LFD_SPOT       (LOGICAL) : Switch to perform spot finite diff test
!  (3 ) N_CALC_STOP    (INTEGER) : Current iteration number 
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC_STOP
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_GLOB 
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_SPOT 
      USE GEOS_CHEM_ADJ_MOD, ONLY : DO_GEOS_CHEM_ADJ
      
      !=================================================================
      ! MAYBE_DO_GEOS_CHEM_ADJ begins here!
      !=================================================================

      ! For global finite difference test we compare the average of 
      ! two finite difference sensitivities with an adjoint sensitivity
      ! around the base case. 
      IF ( LFD_GLOB ) THEN

         ! Only calculate the adjoint during the first iteration  
         IF ( N_CALC_STOP == 1 ) THEN

            CALL DO_GEOS_CHEM_ADJ

         ! Don't bother with more than 3 iterations 
         ELSEIF ( N_CALC_STOP > 3 ) THEN

            CALL ERROR_STOP('To many iterations for FD_GLOB',
     &                      'inverse_mod.f'                  )
         ENDIF


      ! For SPOT finite difference test we compare the average of 
      ! two adjoint sensitivities with a finite difference sensitivity 
      ! around the base case + 1/2 FD_DIFF
      ELSEIF ( LFD_SPOT ) THEN 
               
        ! calculate the adjoint during the first and second iteration  
        IF ( N_CALC_STOP == 1 .or. N_CALC_STOP == 2 ) THEN
                   
           CALL DO_GEOS_CHEM_ADJ
                      
        ! Don't bother with more than 2 iteratoins 
         ELSEIF ( N_CALC_STOP > 2 ) THEN
                  
           CALL ERROR_STOP('To many iterations for FD_SPOT',
     &                     'inverse_mod.f'                  )
         ENDIF

      ENDIF

      ! Return to calling program
      END SUBROUTINE MAYBE_DO_GEOS_CHEM_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE INIT_INVERSE
!
!******************************************************************************
!  Subroutine INIT_INVERSE initializes and zeros all allocatable arrays
!  declared in "inverse_mod.f" (dkh, 1/26/05)
!
!  NOTES:
!  (1 ) Now also allocate EMS_ICS_orig (dkh, 03/29/05)
!  (2 ) Now check for incompatible preproc. definitions and ACTIVE_VARS. (dkh, 10/17/06)  
!  (3 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!
!******************************************************************************
!     
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : NOPT
      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL
      USE ADJ_ARRAYS_MOD,  ONLY : MMSCL
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_STRAT
      USE ERROR_MOD, ONLY       : ALLOC_ERR, ERROR_STOP

#     include "CMN_SIZE"        ! Size parameters

      ! Local variables 
      LOGICAL, SAVE            :: IS_INIT = .FALSE.
      INTEGER                  :: AS, I

      !=================================================================
      ! INIT_INVERSE begins here!
      !=================================================================

      ! Return if we have already initialized
      IF ( IS_INIT ) RETURN

      ! Allocate arrays
      ALLOCATE( GRADNT( NOPT ), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRADNT' )
     
      ALLOCATE( X( NOPT ), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'X' )

      IF ( LADJ_STRAT ) THEN 
         ALLOCATE( GRADNT_P( IIPAR*JJPAR*MMSCL*NSTPL ), STAT = AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRADNT_P' )

         ALLOCATE( GRADNT_L( IIPAR*JJPAR*MMSCL*NSTPL ), STAT = AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRADNT_L' )

         ALLOCATE( XP( IIPAR*JJPAR*MMSCL*NSTPL ), STAT = AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'XP' )

         ALLOCATE( XL( IIPAR*JJPAR*MMSCL*NSTPL ), STAT = AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'XL' )

      ENDIF 

      END SUBROUTINE INIT_INVERSE

!------------------------------------------------------------------------------


      ! Return to calling program 
      SUBROUTINE CLEANUP_INVERSE
!
!******************************************************************************
!  Subroutine CLEANUP_INVERE deallocates all previously allocated arrays
!  for inverse_mod -- call at the end of the program (dkh, 1/26/05) 
!
!  NOTES:
!  (1 ) Now also deallocate EMS_ICS_orig (dkh, 03/29/05)
!  (2 ) No longer make EMS_ICS an array in this module (dkh, 06/08/09)
!  (3 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025) 
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_INVERSE begins here!
      !=================================================================
      IF ( ALLOCATED( GRADNT        ) ) DEALLOCATE( GRADNT        )
      IF ( ALLOCATED( X             ) ) DEALLOCATE( X             )

      ! Return to calling program 
      END SUBROUTINE CLEANUP_INVERSE
      
!------------------------------------------------------------------------------

      END MODULE INVERSE_MOD

 
