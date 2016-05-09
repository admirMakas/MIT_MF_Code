! $Id: rpmares_mod.f,v 1.5 2012/03/01 22:00:27 daven Exp $
      MODULE RPMARES_MOD
!
!******************************************************************************
!  Module RPMARES_MOD contains the RPMARES routines, which compute the aerosol
!  thermodynamical equilibrium. (rjp, bdf, bmy, 11/6/02, 6/11/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) HNO3_SAV    (TYPE (XPLEX) ) : Array to save evolving HNO3 concentrations
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_RPMARES       : Bridge between GEOS-CHEM and RPMARES
!  (2 ) GET_HNO3         : Gets evolving HNO3; relaxes to monthly mean every 3h
!  (3 ) SET_HNO3         : Saves HNO3 in an array for the next timestep
!  (4 ) RPMARES          : Driver for RPMARES code
!  (5 ) AWATER           : Computes thermodynamical equilibrium (?)
!  (6 ) POLY4            : Evaluates a 4th order polynomial expression
!  (7 ) POLY6            : Evaluates a 6th order polynomial expression
!  (8 ) CUBIC            : Solver to find cubic roots
!  (9 ) ACTCOF           : Computes activity coefficients
!  (10) INIT_RPMARES     : Initializes and allocates all module arrays
!  (11) CLEANUP_RPMARES  : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by rpmares_mod.f
!  ============================================================================
!  (1 ) dao_mod.f        : Module w/ arrays for DAO met fields
!  (2 ) error_mod.f      : Module w/ NaN and other error check routines
!  (3 ) time_mod.f       : Module w/ routines to compute date & time
!  (4 ) tracer_mod.f     : Module w/ GEOS-CHEM tracer array STT etc.
!  (5 ) tracerid_mod.f   : Module w/ pointers to tracers & emissions
!  (6 ) tropopause_mod.f : Module
!
!  NOTES:
!  (1 ) Added module variables ELAPSED_SEC and HNO3_sav.  Added module routines
!        GET_HNO3, SET_HNO3, INIT_RPMARES, CLEANUP_RPMARES. (bmy, 12/16/02)
!  (2 ) Replace ALOG with F90 intrinsic LOG to facilitate compilation on 
!        COMPAQ/Alpha platform. (bmy, 3/23/03)
!  (3 ) Now references "time_mod.f".  Now removed ELAPSED_SEC module variable
!        since we can get this info from "time_mod.f".  (bmy, 3/24/03)
!  (4 ) Now references "tracer_mod.f" (bmy, 7/20/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now only apply RPMARES to boxes in the troposphere. (bmy, 4/3/08)
!  (7 ) Add bug fix for low ammonia case in RPMARES (phs, bmy, 4/10/08)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "rpmares_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      ! ... except these variables ...
      PRIVATE :: HNO3_sav

      ! ... and these routines 
      ! adj_group: these are needed in rpmares_adj_mod
      !PRIVATE :: POLY4,    POLY6,    CUBIC,   ACTCOF  
      !PRIVATE :: GET_HNO3, SET_HNO3, RPMARES, AWATER 
      PRIVATE :: GET_HNO3, SET_HNO3, RPMARES
      PRIVATE :: INIT_RPMARES

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      TYPE (XPLEX), ALLOCATABLE :: HNO3_sav(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_RPMARES
!
!******************************************************************************
!  Subroutine DO_RPMARES is the interface between the GEOS-CHEM model
!  and the aerosol thermodynamical equilibrium routine in "rpmares.f"
!  (rjp, bdf, bmy, 12/17/01, 4/10/08)
!
!  NOTES
!  (1 ) Bundled into "rpmares_mod.f" (bmy, 11/15/02)
!  (2 ) Now let HNO3 concentration evolve, but relax to monthly mean values
!        every 3 hours, via routines GET_HNO3 and SET_HNO3. (bmy, 12/16/02)
!  (3 ) Now removed NSEC from the arg list -- use GET_ELAPSED_SEC() function
!        from the new "time_mod.f".  Also use GET_MONTH from "time_mod.f".
!        (bmy, 3/24/03)
!  (4 ) Now reference STT, ITS_A_FULLCHEM_SIM, and ITS_AN_AEROSOL_SIM
!        from "tracer_mod.f".  Now references ITS_A_NEW_MONTH from
!        "time_mod.f" (bmy, 7/20/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now limit RPMARES to the tropopause. (bmy, 4/10/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,         ONLY : AIRVOL, RH, T
      USE GLOBAL_HNO3_MOD, ONLY : GET_GLOBAL_HNO3
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC,    GET_MONTH
      USE TIME_MOD,        ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM 
      USE TRACER_MOD,      ONLY : ITS_AN_AEROSOL_SIM, STT
      USE TRACERID_MOD,    ONLY : IDTSO4,             IDTNH3, IDTNH4 
      USE TRACERID_MOD,    ONLY : IDTNIT,             IDTHNO3 
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT
      ! adj_group: (dkh, 09/07/09) 
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD, NFD
      USE CHECKPT_MOD,     ONLY : RP_IN
      USE CHECKPT_MOD,     ONLY : RP_OUT
      USE CHECKPT_MOD,     ONLY : NITR_MAX
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE          :: FIRST     = .TRUE.
      INTEGER, SAVE          :: LASTMONTH = -99
      INTEGER                :: I,    J,     L
      TYPE (XPLEX)           :: ARH,  ATEMP, AVOL, SO4,  ASO4, ANO3    
      TYPE (XPLEX)           :: AH2O, ANH4,  GNH3, GNO3, AHSO4   
      CHARACTER(LEN=255)     :: X 
      ! adj_group: need to track EXIT status (dkh, 09/07/09) 
      INTEGER                :: EXIT

      ! concentration lower limit [ug/m3 ]
      TYPE (XPLEX),  PARAMETER     :: CONMIN = xplex(1.0D-30,0d0)

      !=================================================================
      ! DO_RPMARES begins here!
      !=================================================================

      ! Initialize on first call
      IF ( FIRST ) THEN
         CALL INIT_RPMARES
         FIRST = .FALSE.
      ENDIF
      
      ! Error check tracer ID's
      X = 'DO_RPMARES (rpmares_mod.f)'
      IF ( IDTSO4 == 0 ) CALL ERROR_STOP( 'IDTSO4 is not defined!', X )
      IF ( IDTNH3 == 0 ) CALL ERROR_STOP( 'IDTNH3 is not defined!', X )
      IF ( IDTNH4 == 0 ) CALL ERROR_STOP( 'IDTNH4 is not defined!', X )
      IF ( IDTNIT == 0 ) CALL ERROR_STOP( 'IDTNIT is not defined!', X )
    
      ! Check to see if we have to read in monthly mean HNO3
      IF ( IDTHNO3 == 0 ) THEN

         IF ( ITS_A_FULLCHEM_SIM() ) THEN

            ! Coupled simulation: stop w/ error since we need HNO3
            CALL ERROR_STOP( 'IDTHNO3 is not defined!', X )

         ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! Offline simulation: read monthly mean HNO3
            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_HNO3( GET_MONTH() )
            ENDIF

         ELSE

            ! Otherwise stop w/ error
            CALL ERROR_STOP( 'Invalid simulation type!', X )

         ENDIF
      ENDIF

      ! adj_group: reset iteration counter to 1 (dkh, 9/10/04, 09/07/09) 
      IF ( LADJ ) NITR_MAX(:,:,:) = 1

      !=================================================================
      ! Get equilibrium values of water, ammonium  and nitrate content
      !=================================================================
! adj_group: add EXIT 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,    L,    ATEMP, ARH,  AVOL,  SO4  )
!$OMP+PRIVATE( ANH4, ANO3, GNH3, GNO3,  ASO4, AHSO4, AH2O )
!$OMP+PRIVATE( EXIT ) 
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Skip if we are in the stratosphere (bmy, 4/3/08)
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! Temperature [K], RH [unitless], and volume [m3]
         ATEMP = T(I,J,L)
         ARH   = RH(I,J,L) * 1.d-2
         AVOL  = AIRVOL(I,J,L)

         ! Convert sulfate, ammonimum, gaseous NH3, gaseous HNO3, 
         ! and aerosol NO3  from [kg] to [ug/m3].
         SO4   = MAX( STT(I,J,L,IDTSO4) * 1.d9 / AVOL, CONMIN )
         GNH3  = MAX( STT(I,J,L,IDTNH3) * 1.d9 / AVOL, CONMIN ) 
         ANH4  = MAX( STT(I,J,L,IDTNH4) * 1.d9 / AVOL, CONMIN ) 
         ANO3  = MAX( STT(I,J,L,IDTNIT) * 1.d9 / AVOL, CONMIN )

         ! For coupled simulations, use HNO3 tracer from STT array.
         ! For offline simulations, call GET_HNO3, which lets HNO3
         ! conc's evolve, but relaxes to monthly mean values every 3h.
         IF ( IDTHNO3 > 0 ) THEN
            GNO3 = MAX( STT(I,J,L,IDTHNO3) * 1.d9 / AVOL, CONMIN )
         ELSE
            GNO3 = MAX( GET_HNO3( I, J, L ), CONMIN )
         ENDIF

         !==============================================================
         ! Call the RPMARES code with the following quantities:
         ! 
         ! SO4   : Total sulfate as sulfate                  [ug/m3]
         ! GNO3  : Nitric Acid vapor (actually gaseous HNO3) [ug/m3]
         ! GNH3  : Gas phase ammonia                         [ug/m3]
         ! ARH   : Fractional relative humidity              [unitless]
         ! ATEMP : Temperature                               [K]
         ! ASO4  : Aerosol phase sulfate                     [ug/m3]
         ! AHSO4 : Aerosol phase in bisulfate                [ug/m3]
         ! ANO3  : Aerosol phase nitrate                     [ug/m3]
         ! AH2O  : Aerosol phase water                       [ug/m3]
         ! ANH4  : Aerosol phase ammonium                    [ug/m3]
         !==============================================================
         ! adj_group: call special version for adjoint (dkh, 9/10/04, 09/07/09) 
         ! Now always use RPMARES_FORADJ, and check LADJ therein (dkh, 02/12/12, adj32_001) 
         !IF ( .not. LADJ ) THEN 
         !   CALL RPMARES( SO4,  GNO3,  GNH3, ARH,  ATEMP,
         !                 ASO4, AHSO4, ANO3, AH2O, ANH4 )
         ! ELSE 

         IF ( LADJ ) THEN 
  
            ! adj_group: Checkpoint inputs to RPMARES
            RP_IN(I,J,L,1)  = SO4
            RP_IN(I,J,L,2)  = GNO3
            RP_IN(I,J,L,3)  = GNH3
            RP_IN(I,J,L,4)  = ANO3
            RP_IN(I,J,L,5)  = ANH4
            RP_IN(I,J,L,6)  = ARH
            RP_IN(I,J,L,7)  = ATEMP

         ENDIF 

         ! adj_group: call modified version
         CALL RPMARES_FORADJ( SO4,  GNO3,  GNH3, ARH,  ATEMP,
     &                       ASO4, AHSO4, ANO3, AH2O, ANH4,
     &                          I,    J,     L,    EXIT       )

         IF ( LADJ ) THEN 

            ! adj_group: checkpoint outputs 
            RP_OUT(I,J,L,1) = ASO4
            RP_OUT(I,J,L,2) = AHSO4
            RP_OUT(I,J,L,3) = ANO3
            RP_OUT(I,J,L,4) = AH2O
            RP_OUT(I,J,L,5) = ANH4
            RP_OUT(I,J,L,6) = SO4
            RP_OUT(I,J,L,7) = GNO3
            RP_OUT(I,J,L,8) = GNH3
            RP_OUT(I,J,L,9) = EXIT
       
            IF ( LPRINTFD
     &         .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
               print*, ' After RECOMP_RPMARES ', I, J, L
               print*, ' RP_IN = ', RP_IN(I,J,L,:)
               print*, ' RP_OUT = ', RP_OUT(I,J,L,:)
               print*, ' EXIT =   ', RP_OUT(I,J,L,9), L
            ENDIF
 
         ENDIF 

         ! Convert modified concentrations from [ug/m3] to [kg] 
         ! for ammonium, ammonia, nitric acid (g), and Nitrate         
         ! NOTE: We don't modify the total sulfate mass.
         STT(I,J,L,IDTNH3) = MAX( GNH3 * AVOL * 1.d-9, CONMIN )  
         STT(I,J,L,IDTNH4) = MAX( ANH4 * AVOL * 1.d-9, CONMIN )  
         STT(I,J,L,IDTNIT) = MAX( ANO3 * AVOL * 1.d-9, CONMIN )  

         ! For coupled runs, convert HNO3 [kg] and store in STT.
         ! For offline runs, save evolving HNO3 [ug/m3] for next timestep.
         IF ( IDTHNO3 > 0 ) THEN
            STT(I,J,L,IDTHNO3) = MAX( GNO3 * AVOL * 1.d-9, CONMIN )  
         ELSE
            CALL SET_HNO3( I, J, L, GNO3 )
         ENDIF   
     
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO      
    
      ! Return to calling program
      END SUBROUTINE DO_RPMARES

!------------------------------------------------------------------------------

      SUBROUTINE RECOMP_RPMARES 
!
!******************************************************************************
! Subroutine RECOMP_RPMARES recomputes the resuls of the aerosol thermo
! calculation in order to be used in the adjoint routins, adrpmares. This
! MAY be faster than saving these values to a checkpt file during forward
! integration.  (dkh, 2/08/05)
!
!  NOTES
!  (1 ) Based on DO_RPMARES 
!  (2 ) Updated to GCv8 (dkh, 09/08/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT
      ! adj_group: (dkh, 09/07/09) 
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD, NFD
      USE CHECKPT_MOD,     ONLY : RP_IN
      USE CHECKPT_MOD,     ONLY : RP_OUT
      USE CHECKPT_MOD,     ONLY : NITR_MAX
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE          :: FIRST     = .TRUE.
      INTEGER, SAVE          :: LASTMONTH = -99
      INTEGER                :: I,    J,     L
      TYPE (XPLEX)           :: ARH,  ATEMP, AVOL, SO4,  ASO4, ANO3    
      TYPE (XPLEX)                 :: AH2O, ANH4,  GNH3, GNO3, AHSO4   
      CHARACTER(LEN=255)     :: X 
      ! adj_group: need to track EXIT status (dkh, 09/07/09) 
      INTEGER                :: EXIT

      ! concentration lower limit [ug/m3 ]
      TYPE (XPLEX),  PARAMETER     :: CONMIN = xplex(1.0D-30,0d0)

      !=================================================================
      ! RECOMP_RPMARES begins here!
      !=================================================================

      ! adj_group: reset iteration counter to 1 (dkh, 9/10/04, 09/07/09) 
      IF ( LADJ ) NITR_MAX(:,:,:) = 1

      !=================================================================
      ! Get equilibrium values of water, ammonium  and nitrate content
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,    L,    ATEMP, ARH,  AVOL,  SO4  )
!$OMP+PRIVATE( ANH4, ANO3, GNH3, GNO3,  ASO4, AHSO4, AH2O )
!$OMP+PRIVATE( EXIT )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Skip if we are in the stratosphere (bmy, 4/3/08)
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         SO4   = RP_IN(I,J,L,1)
         GNO3  = RP_IN(I,J,L,2)
         GNH3  = RP_IN(I,J,L,3)
         ANO3  = RP_IN(I,J,L,4)
         ANH4  = RP_IN(I,J,L,5)
         ARH   = RP_IN(I,J,L,6)
         ATEMP = RP_IN(I,J,L,7)

         !==============================================================
         ! Call the RPMARES_FORADJ code with the following quantities:
         !
         ! SO4   : Total sulfate as sulfate                  [ug/m3]
         ! GNO3  : Nitric Acid vapor (actually gaseous HNO3) [ug/m3]
         ! GNH3  : Gas phase ammonia                         [ug/m3]
         ! ARH   : Fractional relative humidity              [unitless]
         ! ATEMP : Temperature                               [K]
         ! ASO4  : Aerosol phase sulfate                     [ug/m3]
         ! AHSO4 : Aerosol phase in bisulfate                [ug/m3]
         ! ANO3  : Aerosol phase nitrate                     [ug/m3]
         ! AH2O  : Aerosol phase water                       [ug/m3]
         ! ANH4  : Aerosol phase ammonium                    [ug/m3]
         !==============================================================
         CALL RPMARES_FORADJ( SO4,  GNO3,  GNH3, ARH,  ATEMP,
     &                       ASO4, AHSO4, ANO3, AH2O, ANH4,
     &                       I,    J,     L,    EXIT )

         RP_OUT(I,J,L,1) = ASO4
         RP_OUT(I,J,L,2) = AHSO4
         RP_OUT(I,J,L,3) = ANO3
         RP_OUT(I,J,L,4) = AH2O
         RP_OUT(I,J,L,5) = ANH4
         RP_OUT(I,J,L,6) = SO4
         RP_OUT(I,J,L,7) = GNO3
         RP_OUT(I,J,L,8) = GNH3
         RP_OUT(I,J,L,9) = EXIT
         
         IF ( LPRINTFD 
     &      .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
            print*, ' After RECOMP_RPMARES ', I, J, L
            print*, ' RP_IN = ', RP_IN(I,J,L,:)
            print*, ' RP_OUT = ', RP_OUT(I,J,L,:)
            print*, ' EXIT =   ', RP_OUT(I,J,L,9), L
         ENDIF
     
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO      
    
      ! Return to calling program
      END SUBROUTINE RECOMP_RPMARES

!------------------------------------------------------------------------------

      FUNCTION GET_HNO3( I, J, L ) RESULT ( HNO3_UGM3 )
!
!******************************************************************************
!  Subroutine GET_HNO3 allows the HNO3 concentrations to evolve with time,
!  but relaxes back to the monthly mean concentrations every 3 hours.
!  (bmy, 12/16/02, 3/24/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box lon, lat, alt indices
!
!  Function Value: 
!  ============================================================================
!  (1  ) HNO3_UGM3 (TYPE (XPLEX) ) : HNO3 concentration in ug/m3
!
!  NOTES:
!  (1 ) Now use function GET_ELAPSED_MIN() from the new "time_mod.f" to
!        get the elapsed minutes since the start of run. (bmy, 3/24/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      TYPE (XPLEX)              :: HNO3_UGM3

      !=================================================================
      ! GET_HNO3 begins here!
      !=================================================================

      ! Relax to monthly mean HNO3 concentrations every 3 hours
      ! Otherwise just return the concentration in HNO3_sav
      IF ( MOD( GET_ELAPSED_MIN(), 180 ) == 0 ) THEN
         HNO3_UGM3 = GET_HNO3_UGM3( I, J, L )
      ELSE
         HNO3_UGM3 = HNO3_sav(I,J,L)
      ENDIF

      ! Return to calling program
      END FUNCTION GET_HNO3

!------------------------------------------------------------------------------

      SUBROUTINE SET_HNO3( I, J, L, HNO3_UGM3 )
!
!******************************************************************************
!  Subroutine SET_HNO3 stores the modified HNO3 value back into the HNO3_sav
!  array for the next timestep. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box lon, lat, alt indices
!  (4  ) HNO3_UGM3 (TYPE (XPLEX) ) : HNO3 concentration in ug/m3
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      TYPE (XPLEX),  INTENT(IN) :: HNO3_UGM3
      
      !=================================================================
      ! SET_HNO3 begins here!
      !=================================================================
      HNO3_sav(I,J,L) = HNO3_UGM3

      ! Return to calling program
      END SUBROUTINE SET_HNO3

!------------------------------------------------------------------------------

      SUBROUTINE RPMARES( SO4,  GNO3,  GNH3, RH,   TEMP,
     &                    ASO4, AHSO4, ANO3, AH2O, ANH4 )
!
!******************************************************************************
!
! Description:
!
!   ARES calculates the chemical composition of a sulfate/nitrate/
!   ammonium/water aerosol based on equilibrium thermodynamics.
!
!   This code considers two regimes depending upon the molar ratio
!   of ammonium to sulfate.
!
!   For values of this ratio less than 2,the code solves a cubic for
!   hydrogen ion molality, H+,  and if enough ammonium and liquid
!   water are present calculates the dissolved nitric acid. For molal
!   ionic strengths greater than 50, nitrate is assumed not to be present.
!
!   For values of the molar ratio of 2 or greater, all sulfate is assumed
!   to be ammonium sulfate and a calculation is made for the presence of
!   ammonium nitrate.
!
!   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!   obtain the activity coefficients. Abandoned -7/30/97 FSB
!
!   The Bromley method of calculating the multicomponent activity coefficients
!    is used in this version 7/30/97 SJR/FSB
!
!   The calculation of liquid water
!   is done in subroutine water. Details for both calculations are given
!   in the respective subroutines.
!
!   Based upon MARS due to
!   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
!   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!
!   and SCAPE due to
!   Kim, Seinfeld, and Saxeena, Aerosol Sience and Technology,
!   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!
! NOTE: All concentrations supplied to this subroutine are TOTAL
!       over gas and aerosol phases
!
! Parameters:
!
!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate 
!  GNO3  : Nitric Acid vapor in MICROGRAMS/M**3 as nitric acid 
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 
!  RH    : Fractional relative humidity 
!  TEMP  : Temperature in Kelvin 
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 
!  AHSO4 : Aerosol phase in bisulfate in MICROGRAMS/M**3 [rjp, 12/12/01]
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 
!
! Revision History:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   11/10/87  Received the first version of the MARS code
!   S.Roselle   12/30/87  Restructured code
!   S.Roselle   2/12/88   Made correction to compute liquid-phase
!                         concentration of H2O2.
!   S.Roselle   5/26/88   Made correction as advised by SAI, for
!                         computing H+ concentration.
!   S.Roselle   3/1/89    Modified to operate with EM2
!   S.Roselle   5/19/89   Changed the maximum ionic strength from
!                         100 to 20, for numerical stability.
!   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!                         using equations for nitrate budget.
!   F.Binkowski 6/18/91   New ammonia poor case which
!                         omits letovicite.
!   F.Binkowski 7/25/91   Rearranged entire code, restructured
!                         ammonia poor case.
!   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!                         as SO4--
!   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
!                         there is only neutralized sulfate (ammonium
!                         sulfate) and sulfuric acid.
!   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement
!                          with the Cohen et al. (1987)  maximum molality
!                          of 36.2 in Table III.( J. Phys Chem (91) page
!                          4569, and Table IV p 4587.)
!   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!                         possibility for denomenator becoming zero;
!                         this involved solving for H+ first.
!                         Note that for a relative humidity
!                          less than 50%, the model assumes that there is no
!                          aerosol nitrate.
!   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)
!                          Redid logic as follows
!                         1. Water algorithm now follows Spann & Richardson
!                         2. Pitzer Multicomponent method used
!                         3. Multicomponent practical osmotic coefficient
!                            use to close iterations.
!                         4. The model now assumes that for a water
!                            mass fraction WFRAC less than 50% there is
!                            no aerosol nitrate.
!   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor
!                         case, and changed the WFRAC criterion to 40%.
!                         For ammonium to sulfate ratio less than 1.0
!                         all ammonium is aerosol and no nitrate aerosol
!                         exists.
!   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!                         allow gas-phase ammonia to exist.
!   F.Binkowski 7/26/95   Changed equilibrium constants to values from
!                         Kim et al. (1993)
!   F.Binkowski 6/27/96   Changed to new water format
!   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
!                         activity coefficients. The binary activity
!                         coefficients
!                         are the same as the previous version
!   F.Binkowski 8/1/97    Changed minimum sulfate from 0.0 to 1.0e-6 i.e.
!                         1 picogram per cubic meter
!   F.Binkowski 2/23/98   Changes to code made by Ingmar Ackermann to
!                         deal with precision problems on workstations 
!                         incorporated in to this version.  Also included
!                         are his improved descriptions of variables. 
!  F. Binkowski 8/28/98   changed logic as follows: 
!                         If iterations fail, initial values of nitrate
!                          are retained. 
!                         Total mass budgets are changed to account for gas
!                         phase returned.
!  F.Binkowski 10/01/98   Removed setting RATIO to 5 for low to 
!                         to zero sulfate sulfate case.
!  F.Binkowski 01/10/2000 reconcile versions
!
!  F.Binkowski 05/17/2000 change to logic for calculating RATIO
!  F.Binkowski 04/09/2001 change for very low values of RATIO,
!                         RATIO < 0.5, no iterative calculations are done
!                         in low ammonia case a MAX(1.0e-10, MSO4) IS
!                         applied, and the iteration count is
!                         reduced to fifty for each iteration loop.
!  R. Yantosca 09/25/2002 Bundled into "rpmares_mod.f".  Declared all REALs
!                          as TYPE (XPLEX)'s.  Cleaned up comments.  Also now force
!                          TYPE (XPLEX) explicitly with "D" exponents.
!  P. Le Sager and        Bug fix for low ammonia case -- prevent floating
!  R. Yantosca 04/10/2008  point underflow and NaN's.
!  P. Le Sager 06/10/2008 Better catch of over/underflow for low ammonia case
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP, IS_SAFE_DIV

      !=================================================================
      ! ARGUMENTS and their descriptions
      !=================================================================
      TYPE (XPLEX) :: SO4              ! Total sulfate in micrograms / m**3
      TYPE (XPLEX) :: GNO3             ! Gas-phase nitric acid in micrograms / m**3
      TYPE (XPLEX) :: GNH3             ! Gas-phase ammonia in micrograms / m**3 
      TYPE (XPLEX) :: RH               ! Fractional relative humidity
      TYPE (XPLEX) :: TEMP             ! Temperature in Kelvin
      TYPE (XPLEX) :: ASO4             ! Aerosol sulfate in micrograms / m**3
      TYPE (XPLEX) :: AHSO4            ! Aerosol bisulfate in micrograms / m**3
      TYPE (XPLEX) :: ANO3             ! Aerosol nitrate in micrograms / m**3
      TYPE (XPLEX) :: AH2O             ! Aerosol liquid water content water in
                                 !   micrograms / m**3
      TYPE (XPLEX) :: ANH4             ! Aerosol ammonium in micrograms / m**3

      !=================================================================
      ! PARAMETERS and their descriptions:
      !=================================================================

      ! Molecular weights
      TYPE (XPLEX), PARAMETER :: MWNACL = xplex(58.44277d0,0d0)               ! NaCl
      TYPE (XPLEX), PARAMETER :: MWNO3  = xplex(62.0049d0,0d0)                ! NO3
      TYPE (XPLEX), PARAMETER :: MWHNO3 = xplex(63.01287d0,0d0)               ! HNO3
      TYPE (XPLEX), PARAMETER :: MWSO4  = xplex(96.0576d0,0d0)                ! SO4
      TYPE (XPLEX), PARAMETER :: MWHSO4 = xplex(MWSO4%r + 1.0080d0,0d0)         ! HSO4
      TYPE (XPLEX), PARAMETER :: MH2SO4 = xplex(98.07354d0,0d0)               ! H2SO4
      TYPE (XPLEX), PARAMETER :: MWNH3  = xplex(17.03061d0,0d0)               ! NH3
      TYPE (XPLEX), PARAMETER :: MWNH4  = xplex(18.03858d0,0d0)               ! NH4
      TYPE (XPLEX), PARAMETER :: MWORG  = xplex(16.0d0,0d0)                   ! Organic Species
      TYPE (XPLEX), PARAMETER :: MWCL   = xplex(35.453d0,0d0)                 ! Chloride
      TYPE (XPLEX), PARAMETER :: MWAIR  = xplex(28.964d0,0d0)                 ! AIR
      TYPE (XPLEX), PARAMETER :: MWLCT  = xplex(3.0d0 * MWNH4%r +          ! Letovicite
     &                              2.0d0 * MWSO4%r + 1.0080d0,0d0)  
      TYPE (XPLEX), PARAMETER :: MWAS=xplex(2.0d0*MWNH4%r+MWSO4%r,0d0)    ! Amm. Sulfate
      TYPE (XPLEX),PARAMETER::MWABS=xplex(MWNH4%r+MWSO4%r+1.0080d0,0d0) ! Amm. Bisulfate

      ! Minimum value of sulfate aerosol concentration
      TYPE (XPLEX), PARAMETER :: MINSO4 = xplex(1.0d-6 / MWSO4%r,0d0) 

      ! Minimum total nitrate cncentration
      TYPE (XPLEX), PARAMETER :: MINNO3 = xplex(1.0d-6 / MWNO3%r,0d0) 

      ! Force a minimum concentration
      TYPE (XPLEX), PARAMETER :: FLOOR  = xplex(1.0d-30 ,0d0)

      ! Tolerances for convergence test.  NOTE: We now have made these
      ! parameters so they don't lose their values (phs, bmy, 4/10/08)
      TYPE (XPLEX), PARAMETER :: TOLER1 = xplex(0.00001d0,0d0)       
      TYPE (XPLEX), PARAMETER :: TOLER2 = xplex(0.001d0,0d0)       

      !=================================================================
      ! SCRATCH LOCAL VARIABLES and their descriptions:
      !=================================================================

      INTEGER :: IRH              ! Index set to percent relative humidity
      INTEGER :: NITR             ! Number of iterations for activity
                                  !   coefficients
      INTEGER :: NNN              ! Loop index for iterations
      INTEGER :: NR               ! Number of roots to cubic equation for
                                  ! H+ ciaprecision
      TYPE (XPLEX)  :: A0               ! Coefficients and roots of
      TYPE (XPLEX)  :: A1               ! Coefficients and roots of
      TYPE (XPLEX)  :: A2               ! Coefficients and roots of
      TYPE (XPLEX)  :: AA               ! Coefficients and discriminant for
                                  ! quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: BAL              ! internal variables ( high ammonia case)
      TYPE (XPLEX)  :: BB               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: BHAT             ! Variables used for ammonia solubility
      TYPE (XPLEX)  :: CC               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: CONVT            ! Factor for conversion of units
      TYPE (XPLEX)  :: DD               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: DISC             ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: EROR             ! Relative error used for convergence test
      TYPE (XPLEX)  :: FNH3             ! "Free ammonia concentration", that
                                  !   which exceeds TWOSO4
      TYPE (XPLEX)  :: GAMAAB           ! Activity Coefficient for (NH4+,
                                  !   HSO4-)GAMS( 2,3 )
      TYPE (XPLEX)  :: GAMAAN           ! Activity coefficient for (NH4+, NO3-)
                                  !   GAMS( 2,2 )
      TYPE (XPLEX)  :: GAMAHAT          ! Variables used for ammonia solubility
      TYPE (XPLEX)  :: GAMANA           ! Activity coefficient for (H+ ,NO3-)
                                  !   GAMS( 1,2 )
      TYPE (XPLEX)  :: GAMAS1           ! Activity coefficient for (2H+, SO4--)
                                  !   GAMS( 1,1 )
      TYPE (XPLEX)  :: GAMAS2           ! Activity coefficient for (H+, HSO4-)
                                  !   GAMS( 1,3 )
      TYPE (XPLEX)  :: GAMOLD           ! used for convergence of iteration
      TYPE (XPLEX)  :: GASQD            ! internal variables ( high ammonia case)
      TYPE (XPLEX)  :: HPLUS            ! Hydrogen ion (low ammonia case) (moles
                                  !   / kg water)
      TYPE (XPLEX)  :: K1A              ! Equilibrium constant for ammonia to
                                  !   ammonium
      TYPE (XPLEX)  :: K2SA             ! Equilibrium constant for
                                  !   sulfate-bisulfate (aqueous)
      TYPE (XPLEX)  :: K3               ! Dissociation constant for ammonium
                                  !   nitrate
      TYPE (XPLEX)  :: KAN              ! Equilibrium constant for ammonium
                                  !   nitrate (aqueous)
      TYPE (XPLEX)  :: KHAT             ! Variables used for ammonia solubility
      TYPE (XPLEX)  :: KNA              ! Equilibrium constant for nitric acid
                                  !   (aqueous)
      TYPE (XPLEX)  :: KPH              ! Henry's Law Constant for ammonia
      TYPE (XPLEX)  :: KW               ! Equilibrium constant for water
                                  !  dissociation
      TYPE (XPLEX)  :: KW2              ! Internal variable using KAN
      TYPE (XPLEX)  :: MAN              ! Nitrate (high ammonia case) (moles /
                                  !   kg water)
      TYPE (XPLEX)  :: MAS              ! Sulfate (high ammonia case) (moles /
                                  !   kg water)
      TYPE (XPLEX)  :: MHSO4            ! Bisulfate (low ammonia case) (moles /
                                  !   kg water)
      TYPE (XPLEX)  :: MNA              ! Nitrate (low ammonia case) (moles / kg
                                  !   water)
      TYPE (XPLEX)  :: MNH4             ! Ammonium (moles / kg water)
      TYPE (XPLEX)  :: MOLNU            ! Total number of moles of all ions
      TYPE (XPLEX)  :: MSO4             ! Sulfate (low ammonia case) (moles / kg
                                  !   water)
      TYPE (XPLEX)  :: PHIBAR           ! Practical osmotic coefficient
      TYPE (XPLEX)  :: PHIOLD           ! Previous value of practical osmotic
                                  !   coefficient used for convergence of
                                  !   iteration
      TYPE (XPLEX)  :: RATIO            ! Molar ratio of ammonium to sulfate
      TYPE (XPLEX)  :: RK2SA            ! Internal variable using K2SA
      TYPE (XPLEX)  :: RKNA             ! Internal variables using KNA
      TYPE (XPLEX)  :: RKNWET           ! Internal variables using KNA
      TYPE (XPLEX)  :: RR1
      TYPE (XPLEX)  :: RR2
      TYPE (XPLEX)  :: STION            ! Ionic strength
      TYPE (XPLEX)  :: T1               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T2               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T21              ! Internal variables of convenience (low
                                  !   ammonia case)
      TYPE (XPLEX)  :: T221             ! Internal variables of convenience (low
                                  !   ammonia case)
      TYPE (XPLEX)  :: T3               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T4               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T6               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: TNH4             ! Total ammonia and ammonium in
                                  !   micromoles / meter ** 3
      TYPE (XPLEX)  :: TNO3             ! Total nitrate in micromoles / meter ** 3
      TYPE (XPLEX)  :: TSO4             ! Total sulfate in micromoles / meter ** 3
      TYPE (XPLEX)  :: TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles
                                  !   / kg water)
      TYPE (XPLEX)  :: WFRAC            ! Water mass fraction
      TYPE (XPLEX)  :: WH2O             ! Aerosol liquid water content (internally)
                                  !   micrograms / meter **3 on output
                                  !   internally it is 10 ** (-6) kg (water)
                                  !   / meter ** 3
                                  !   the conversion factor (1000 g = 1 kg)
                                  !   is applied for AH2O output
      TYPE (XPLEX)  :: WSQD             ! internal variables ( high ammonia case)
      TYPE (XPLEX)  :: XNO3             ! Nitrate aerosol concentration in
                                  ! micromoles / meter ** 3
      TYPE (XPLEX)  :: XXQ              ! Variable used in quadratic solution
      TYPE (XPLEX)  :: YNH4             ! Ammonium aerosol concentration in
                                  !  micromoles / meter** 3
      TYPE (XPLEX)  :: ZH2O             ! Water variable saved in case ionic
                                  !  strength too high.
      TYPE (XPLEX)  :: ZSO4             ! Total sulfate molality - mso4 + mhso4
                                  !  (low ammonia case) (moles / kg water)
      TYPE (XPLEX)  :: CAT( 2 )         ! Array for cations (1, H+); (2, NH4+)
                                  !  (moles / kg water)
      TYPE (XPLEX)  :: AN ( 3 )         ! Array for anions (1, SO4--); (2,
                                  !   NO3-); (3, HSO4-)  (moles / kg water)
      TYPE (XPLEX)  :: CRUTES( 3 )      ! Coefficients and roots of
      TYPE (XPLEX)  :: GAMS( 2, 3 )     ! Array of activity coefficients
      TYPE (XPLEX)  :: TMASSHNO3        ! Total nitrate (vapor and particle) 
      TYPE (XPLEX)  :: GNO3_IN, ANO3_IN                
      
      !=================================================================
      ! RPMARES begins here!
      ! convert into micromoles/m**3
      !=================================================================

      ! For extremely low relative humidity ( less than 1% ) set the 
      ! water content to a minimum and skip the calculation.
      IF ( RH .LT. 0.01 ) THEN
         AH2O = FLOOR
         RETURN
      ENDIF 

      ! total sulfate concentration
      TSO4 = MAX( FLOOR, SO4 / MWSO4  )      
      ASO4 = SO4

      !Cia models3 merge NH3/NH4 , HNO3,NO3 here
      !c *** recommended by Dr. Ingmar Ackermann

      ! total nitrate
      TNO3      = MAX( 0.0d0, ( ANO3 / MWNO3 + GNO3 / MWHNO3 ) )            

      ! total ammonia
      TNH4      = MAX( 0.0d0, ( GNH3 / MWNH3 + ANH4 / MWNH4 )  )

      GNO3_IN   = GNO3
      ANO3_IN   = ANO3
      TMASSHNO3 = MAX( 0.0d0, GNO3 + ANO3 )     

      ! set the  molar ratio of ammonium to sulfate
      RATIO = TNH4 / TSO4

      ! validity check for negative concentration
      IF ( TSO4 < 0.0d0 .OR. TNO3 < 0.0d0 .OR. TNH4 < 0.0d0 ) THEN
          PRINT*, 'TSO4 : ', TSO4
          PRINT*, 'TNO3 : ', TNO3
          PRINT*, 'TNH4 : ', TNH4
          CALL GEOS_CHEM_STOP
      ENDIF   

      ! now set humidity index IRH as a percent
      IRH = NINT( 100.0 * RH )

      ! now set humidity index IRH as a percent
      IRH = MAX(  1, IRH )
      IRH = MIN( 99, IRH )

      !=================================================================
      ! Specify the equilibrium constants at  correct temperature.  
      ! Also change units from ATM to MICROMOLE/M**3 (for KAN, KPH, and K3 )
      ! Values from Kim et al. (1993) except as noted.
      ! Equilibrium constant in Kim et al. (1993)
      !   K = K0 exp[ a(T0/T -1) + b(1+log(T0/T)-T0/T) ], T0 = 298.15 K
      !   K = K0 EXP[ a T3 + b T4 ] in the code here.
      !=================================================================
      CONVT = 1.0d0 / ( 0.082d0 * TEMP )
      T6    = 0.082d-9 *  TEMP
      T1    = 298.0d0 / TEMP
      T2    = LOG( T1 )
      T3    = T1 - 1.0d0
      T4    = 1.0d0 + T2 - T1

      !=================================================================
      ! Equilibrium Relation
      ! 
      ! HSO4-(aq)         = H+(aq)   + SO4--(aq)  ; K2SA
      ! NH3(g)            = NH3(aq)               ; KPH
      ! NH3(aq) + H2O(aq) = NH4+(aq) + OH-(aq)    ; K1A
      ! HNO3(g)           = H+(aq)   + NO3-(aq)   ; KNA
      ! NH3(g) + HNO3(g)  = NH4NO3(s)             ; K3
      ! H2O(aq)           = H+(aq)   + OH-(aq)    ; KW
      !=================================================================
      KNA  = 2.511d+06 *  EXP(  29.17d0 * T3 + 16.83d0 * T4 ) * T6
      K1A  = 1.805d-05 *  EXP(  -1.50d0 * T3 + 26.92d0 * T4 )
      K2SA = 1.015d-02 *  EXP(   8.85d0 * T3 + 25.14d0 * T4 )
      KW   = 1.010d-14 *  EXP( -22.52d0 * T3 + 26.92d0 * T4 )
      KPH  = 57.639d0  *  EXP(  13.79d0 * T3 -  5.39d0 * T4 ) * T6
      !K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
      KHAT =  KPH * K1A / KW
      KAN  =  KNA * KHAT

      ! Compute temperature dependent equilibrium constant for NH4NO3
      ! (from Mozurkewich, 1993)
      K3 = EXP( 118.87d0  - 24084.0d0 / TEMP -  6.025d0  * LOG( TEMP ) )

      ! Convert to (micromoles/m**3) **2
      K3     = K3 * CONVT * CONVT

      WH2O   = 0.0d0
      STION  = 0.0d0
      AH2O   = 0.0d0
      MAS    = 0.0d0
      MAN    = 0.0d0
      HPLUS  = 0.0d0
      NITR   = 0
      NR     = 0
      GAMAAN = 1.0d0
      GAMOLD = 1.0d0

      ! If there is very little sulfate and  nitrate 
      ! set concentrations to a very small value and return.
      IF ( ( TSO4 .LT. MINSO4 ) .AND. ( TNO3 .LT. MINNO3 ) ) THEN
         ASO4  = MAX( FLOOR, ASO4  )
         AHSO4 = MAX( FLOOR, AHSO4 ) ! [rjp, 12/12/01]
         ANO3  = MAX( FLOOR, ANO3  )
         ANH4  = MAX( FLOOR, ANH4  )
         WH2O  = FLOOR
         AH2O  = FLOOR
         GNH3  = MAX( FLOOR, GNH3  )
         GNO3  = MAX( FLOOR, GNO3  )
         
         RETURN
      ENDIF

      !=================================================================
      ! High Ammonia Case
      !=================================================================
      IF ( RATIO .GT. 2.0d0 ) THEN
        
         GAMAAN = 0.1d0

         ! Set up twice the sulfate for future use.
         TWOSO4 = 2.0d0 * TSO4
         XNO3   = 0.0d0
         YNH4   = TWOSO4

         ! Treat different regimes of relative humidity
         !
         ! ZSR relationship is used to set water levels. Units are
         !  10**(-6) kg water/ (cubic meter of air)
         !  start with ammomium sulfate solution without nitrate
       
         CALL AWATER( IRH, TSO4, YNH4, TNO3, AH2O ) !**** note TNO3
         WH2O = 1.0d-3 * AH2O
         ASO4 = TSO4   * MWSO4

         ! In sulfate poor case, Sulfate ion is preferred
         ! Set bisulfate equal to zero [rjp, 12/12/01]
         AHSO4 = 0.0d0
         ANO3  = 0.0d0
         ANH4  = YNH4 * MWNH4
         WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )

        !IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
        IF ( WFRAC .LT. 0.2d0 ) THEN

           ! "dry" ammonium sulfate and ammonium nitrate
           ! compute free ammonia 
           FNH3 = TNH4 - TWOSO4
           CC   = TNO3 * FNH3 - K3

           ! check for not enough to support aerosol
           IF ( CC .LE. 0.0d0 ) THEN
              XNO3 = 0.0d0
           ELSE
              AA   = 1.0d0
              BB   = -( TNO3 + FNH3 )
              DISC = BB * BB - 4.0d0 * CC

              ! Check for complex roots of the quadratic
              ! set retain initial values of nitrate and RETURN 
              ! if complex roots are found
              IF ( DISC .LT. 0.0d0 ) THEN
                 XNO3  = 0.0d0
                 AH2O  = 1000.0d0 * WH2O
                 YNH4  = TWOSO4
                 ASO4  = TSO4 * MWSO4
                 AHSO4 = 0.0d0
                 ANH4  = YNH4 * MWNH4
                 GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
                 GNO3  = GNO3_IN
                 ANO3  = ANO3_IN
                 RETURN
              ENDIF

              ! to get here, BB .lt. 0.0, CC .gt. 0.0 always
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN ( 1.0d0, BB ) * DD )


              ! Since both roots are positive, select smaller root.
              XNO3 = MIN( XXQ / AA, CC / XXQ )

           ENDIF                ! CC .LE. 0.0

           AH2O  = 1000.0d0 * WH2O
           YNH4  = TWOSO4 + XNO3          
           ASO4  = TSO4 * MWSO4
           AHSO4 = FLOOR
           ANO3  = XNO3 * MWNO3
           ANH4  = YNH4 * MWNH4
           GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 )  )      
           GNO3  = MAX( FLOOR, ( TMASSHNO3 - ANO3 ) )    
           RETURN
        ENDIF                  ! WFRAC .LT. 0.2

        ! liquid phase containing completely neutralized sulfate and
        ! some nitrate.  Solve for composition and quantity.
        MAS    = TSO4 / WH2O
        MAN    = 0.0d0
        XNO3   = 0.0d0
        YNH4   = TWOSO4
        PHIOLD = 1.0d0

        !===============================================================
        ! Start loop for iteration
        !
        ! The assumption here is that all sulfate is ammonium sulfate,
        ! and is supersaturated at lower relative humidities.
        !===============================================================
        DO NNN = 1, 50 ! loop count reduced 0409/2001 by FSB

           NITR  = NNN
           GASQD = GAMAAN * GAMAAN
           WSQD  = WH2O * WH2O
           KW2   = KAN * WSQD / GASQD
           AA    = 1.0 - KW2
           BB    = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
           CC    = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

           ! This is a quadratic for XNO3 [MICROMOLES / M**3] 
           ! of nitrate in solution
           DISC = BB * BB - 4.0d0 * AA * CC

           ! Check for complex roots, retain inital values and RETURN
           IF ( DISC .LT. 0.0 ) THEN
              XNO3  = 0.0d0
              AH2O  = 1000.0d0 * WH2O
              YNH4  = TWOSO4
              ASO4  = TSO4 * MWSO4
              AHSO4 = FLOOR     ! [rjp, 12/12/01]
              ANH4  = YNH4 * MWNH4
              GNH3  = MWNH3 * MAX( FLOOR, (TNH4 - YNH4 ) )
              GNO3  = GNO3_IN
              ANO3  = ANO3_IN
              RETURN
           ENDIF
           
           ! Deal with degenerate case (yoj)
           IF ( AA .NE. 0.0d0 ) THEN
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN( 1.0d0, BB ) * DD )
              RR1 = XXQ / AA
              RR2 = CC / XXQ

              ! choose minimum positve root
              IF ( ( RR1 * RR2 ) .LT. 0.0d0 ) THEN
                 XNO3 = MAX( RR1, RR2 )
              ELSE
                 XNO3 = MIN( RR1, RR2 )
              ENDIF
           ELSE
              XNO3 = - CC / BB  ! AA equals zero here.
           ENDIF

           XNO3 = MIN( XNO3, TNO3 )
           
           ! This version assumes no solid sulfate forms (supersaturated )
           ! Now update water
           CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

           ! ZSR relationship is used to set water levels. Units are
           ! 10**(-6) kg water/ (cubic meter of air).  The conversion 
           ! from micromoles to moles is done by the units of WH2O.
           WH2O = 1.0d-3 * AH2O 

           ! Ionic balance determines the ammonium in solution.
           MAN  = XNO3 / WH2O
           MAS  = TSO4 / WH2O
           MNH4 = 2.0d0 * MAS + MAN
           YNH4 = MNH4 * WH2O

           ! MAS, MAN and MNH4 are the aqueous concentrations of sulfate, 
           ! nitrate, and ammonium in molal units (moles/(kg water) ).
           STION    = 3.0d0 * MAS + MAN
           CAT( 1 ) = 0.0d0
           CAT( 2 ) = MNH4
           AN ( 1 ) = MAS
           AN ( 2 ) = MAN
           AN ( 3 ) = 0.0d0
           CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )
           GAMAAN = GAMS( 2, 2 )

           ! Use GAMAAN for convergence control
           EROR   = ABS( GAMOLD - GAMAAN ) / GAMOLD
           GAMOLD = GAMAAN

           ! Check to see if we have a solution
           IF ( EROR .LE. TOLER1 ) THEN
              ASO4  = TSO4 * MWSO4
              AHSO4 = 0.0d0       ! [rjp, 12/12/01]
              ANO3  = XNO3 * MWNO3
              ANH4  = YNH4 * MWNH4
              GNO3  = MAX( FLOOR, ( TMASSHNO3  - ANO3 ) )
              GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
              AH2O  = 1000.0d0 * WH2O
              RETURN
           ENDIF

        ENDDO

        ! If after NITR iterations no solution is found, then:
        ! FSB retain the initial values of nitrate particle and vapor
        ! note whether or not convert all bisulfate to sulfate
        ASO4  = TSO4 * MWSO4
        AHSO4 = FLOOR      
        XNO3  = TNO3 / MWNO3
        YNH4  = TWOSO4
        ANH4  = YNH4 * MWNH4

        CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

        GNO3  = GNO3_IN
        ANO3  = ANO3_IN
        GNH3  = MAX( FLOOR, MWNH3 * (TNH4 - YNH4 ) )
        RETURN

      !================================================================
      ! Low Ammonia Case 
      !
      ! Coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98
      ! modified 04/09/2001
      !       
      ! All cases covered by this logic
      !=================================================================
      ELSE

         WH2O = 0.0d0
         CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
         WH2O = 1.0d-3 * AH2O
         ZH2O = AH2O

         ! convert 10**(-6) kg water/(cubic meter of air) to micrograms 
         ! of water per cubic meter of air (1000 g = 1 kg)
         ! in sulfate rich case, preferred form is HSO4-
         !ASO4 = TSO4 * MWSO4
         ASO4  = FLOOR          ![rjp, 12/12/01]
         AHSO4 = TSO4 * MWSO4   ![rjp, 12/12/01]
         ANH4  = TNH4 * MWNH4
         ANO3  = ANO3_IN
         GNO3  = TMASSHNO3 - ANO3
         GNH3  = FLOOR 
        
         !==============================================================
         ! *** Examine special cases and return if necessary.
         !         
         ! FSB For values of RATIO less than 0.5 do no further 
         ! calculations.  The code will cycle and still predict the 
         ! same amount of ASO4, ANH4, ANO3, AH2O so terminate early 
         ! to swame computation
         !==============================================================
         IF ( RATIO .LT. 0.5d0 ) RETURN ! FSB 04/09/2001 
    
         ! Check for zero water.
         IF ( WH2O .EQ. 0.0d0 ) RETURN
         ZSO4 = TSO4 / WH2O

         ! ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4
         ! do not solve for aerosol nitrate for total sulfate molality
         ! greater than 11.0 because the model parameters break down
         !### IF ( ZSO4 .GT. 11.0 ) THEN
         IF ( ZSO4 .GT. 9.0 ) THEN ! 18 June 97
            RETURN
         ENDIF

         ! *** Calculation may now proceed.
         !
         ! First solve with activity coeffs of 1.0, then iterate.
         PHIOLD = 1.0d0
         GAMANA = 1.0d0
         GAMAS1 = 1.0d0
         GAMAS2 = 1.0d0
         GAMAAB = 1.0d0
         GAMOLD = 1.0d0

         ! All ammonia is considered to be aerosol ammonium.
         MNH4 = TNH4 / WH2O

         ! MNH4 is the molality of ammonium ion.
         YNH4 = TNH4

         ! loop for iteration
         DO NNN = 1, 50    ! loop count reduced 04/09/2001 by FSB
            NITR = NNN

            !------------------------------------------------------------
            ! Add robustness: now check if GAMANA or GAMAS1 is too small
            ! for the division in RKNA and RK2SA. If they are, return w/ 
            ! original values: basically replicate the procedure used 
            ! after the current DO-loop in case of no-convergence
            ! (phs, bmy, rjp, 4/10/08)
            ! Now uses IS_SAFE_DIV to avoid compiler/machine dependency 
            ! and to check for both underlow and overflow. Also 
            ! use REAL4 flag to avoid under/overflow when computing A0 
            ! and A1 from RKNA and RK2SA (phs, 5/28/08)
            !------------------------------------------------------------
            IF ( .NOT. (
     &           IS_SAFE_DIV( GAMAS2, GAMAS1*GAMAS1*GAMAS1, R4=.TRUE. )
     &           .AND. IS_SAFE_DIV( KNA, GAMANA*GAMANA, R4=.TRUE. ) 
     &           ) ) THEN
               
               WRITE(6,*) 'RPMARES: not safe to divide...exit'
               CALL flush(6)
               GOTO 100
            ENDIF
            
            ! set up equilibrium constants including activities
            ! solve the system for hplus first then sulfate & nitrate
            RK2SA  = K2SA * GAMAS2 * GAMAS2 / (GAMAS1 * GAMAS1 * GAMAS1)
            RKNA   = KNA / ( GAMANA * GAMANA )
            RKNWET = RKNA * WH2O
            T21    = ZSO4 - MNH4
            T221   = ZSO4 + T21

            ! set up coefficients for cubic
            A2 = RK2SA + RKNWET - T21
            A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )
     &           - RK2SA * ZSO4 - RKNA * TNO3
            A0 = - (T21 * RK2SA * RKNWET
     &           + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )

            CALL CUBIC ( A2, A1, A0, NR, CRUTES )

            ! Code assumes the smallest positive root is in CRUTES(1)
            ! But, it can be negative (see CUBIC, case of one real root, 
            ! but can also be propagated by over/underflown)... if it is
            ! the case then return with original values (phs, 5/27/08)
            HPLUS = CRUTES( 1 )
            IF (HPLUS <= 0d0) GOTO 100
            BAL   = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0

            ! molality of sulfate ion
            MSO4  = RK2SA * ZSO4 / ( HPLUS + RK2SA ) 

            ! molality of bisulfate ion
            ! MAX added 04/09/2001 by FSB
            MHSO4 = MAX( 1.0d-10, ZSO4 - MSO4 ) 

            ! molality of nitrate ion
            MNA   = RKNA * TNO3 / ( HPLUS + RKNWET ) 
            MNA   = MAX( 0.0d0, MNA )
            MNA   = MIN( MNA, TNO3 / WH2O )
            XNO3  = MNA * WH2O
            ANO3  = MNA * WH2O * MWNO3
            GNO3  = MAX( FLOOR, TMASSHNO3 - ANO3 )
            ASO4  = MSO4 * WH2O * MWSO4 ![rjp, 12/12/01]
            AHSO4 = MHSO4 * WH2O * MWSO4 ![rjp, 12/12/01]

            ! Calculate ionic strength
            STION = 0.5d0 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0d0 * MSO4)

            ! Update water
            CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

            ! Convert 10**(-6) kg water/(cubic meter of air) to micrograms 
            ! of water per cubic meter of air (1000 g = 1 kg)
            WH2O     = 1.0d-3 * AH2O
            CAT( 1 ) = HPLUS
            CAT( 2 ) = MNH4
            AN ( 1 ) = MSO4
            AN ( 2 ) = MNA
            AN ( 3 ) = MHSO4

            CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )

            GAMANA = GAMS( 1, 2 )
            GAMAS1 = GAMS( 1, 1 )
            GAMAS2 = GAMS( 1, 3 )
            GAMAAN = GAMS( 2, 2 )

            ! NOTE: Improved for robustness!
            GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
            BHAT = KHAT * GAMAHAT
            !### EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
            !### PHIOLD = PHIBAR
            EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD
            GAMOLD = GAMAHAT

            ! return with good solution
            IF ( EROR .LE. TOLER2 ) THEN
               RETURN
            ENDIF
            
         ENDDO

         ! after NITR iterations, failure to solve the system
         ! convert all ammonia to aerosol ammonium and return input
         ! values of NO3 and HNO3
 100     ANH4 = TNH4 * MWNH4
         GNH3 = FLOOR
         GNO3 = GNO3_IN
         ANO3 = ANO3_IN
         ASO4 = TSO4 * MWSO4    ! [rjp, 12/17/01]
         AHSO4= FLOOR           ! [rjp, 12/17/01]

         CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
         
         RETURN

      ENDIF                     ! ratio .gt. 2.0

      ! Return to calling program
      END SUBROUTINE RPMARES

!------------------------------------------------------------------------------

      ! adj_group: just like RPMARES_FORADJ but with checkpointing needed 
      ! for the adjoint (dkh, 9/10/04, 09/07/09) 
      SUBROUTINE RPMARES_FORADJ( SO4,  GNO3,  GNH3, RH,   TEMP,
     &                    ASO4, AHSO4, ANO3, AH2O, ANH4, 
     &                    I,    J,     L,    EXIT       )
!
!******************************************************************************
!
! Description:
!
!   ARES calculates the chemical composition of a sulfate/nitrate/
!   ammonium/water aerosol based on equilibrium thermodynamics.
!
!   This code considers two regimes depending upon the molar ratio
!   of ammonium to sulfate.
!
!   For values of this ratio less than 2,the code solves a cubic for
!   hydrogen ion molality, H+,  and if enough ammonium and liquid
!   water are present calculates the dissolved nitric acid. For molal
!   ionic strengths greater than 50, nitrate is assumed not to be present.
!
!   For values of the molar ratio of 2 or greater, all sulfate is assumed
!   to be ammonium sulfate and a calculation is made for the presence of
!   ammonium nitrate.
!
!   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!   obtain the activity coefficients. Abandoned -7/30/97 FSB
!
!   The Bromley method of calculating the multicomponent activity coefficients
!    is used in this version 7/30/97 SJR/FSB
!
!   The calculation of liquid water
!   is done in subroutine water. Details for both calculations are given
!   in the respective subroutines.
!
!   Based upon MARS due to
!   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
!   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!
!   and SCAPE due to
!   Kim, Seinfeld, and Saxeena, Aerosol Sience and Technology,
!   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!
! NOTE: All concentrations supplied to this subroutine are TOTAL
!       over gas and aerosol phases
!
! Parameters:
!
!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate 
!  GNO3  : Nitric Acid vapor in MICROGRAMS/M**3 as nitric acid 
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 
!  RH    : Fractional relative humidity 
!  TEMP  : Temperature in Kelvin 
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 
!  AHSO4 : Aerosol phase in bisulfate in MICROGRAMS/M**3 [rjp, 12/12/01]
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 
!
! Revision History:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   11/10/87  Received the first version of the MARS code
!   S.Roselle   12/30/87  Restructured code
!   S.Roselle   2/12/88   Made correction to compute liquid-phase
!                         concentration of H2O2.
!   S.Roselle   5/26/88   Made correction as advised by SAI, for
!                         computing H+ concentration.
!   S.Roselle   3/1/89    Modified to operate with EM2
!   S.Roselle   5/19/89   Changed the maximum ionic strength from
!                         100 to 20, for numerical stability.
!   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!                         using equations for nitrate budget.
!   F.Binkowski 6/18/91   New ammonia poor case which
!                         omits letovicite.
!   F.Binkowski 7/25/91   Rearranged entire code, restructured
!                         ammonia poor case.
!   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!                         as SO4--
!   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
!                         there is only neutralized sulfate (ammonium
!                         sulfate) and sulfuric acid.
!   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement
!                          with the Cohen et al. (1987)  maximum molality
!                          of 36.2 in Table III.( J. Phys Chem (91) page
!                          4569, and Table IV p 4587.)
!   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!                         possibility for denomenator becoming zero;
!                         this involved solving for H+ first.
!                         Note that for a relative humidity
!                          less than 50%, the model assumes that there is no
!                          aerosol nitrate.
!   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)
!                          Redid logic as follows
!                         1. Water algorithm now follows Spann & Richardson
!                         2. Pitzer Multicomponent method used
!                         3. Multicomponent practical osmotic coefficient
!                            use to close iterations.
!                         4. The model now assumes that for a water
!                            mass fraction WFRAC less than 50% there is
!                            no aerosol nitrate.
!   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor
!                         case, and changed the WFRAC criterion to 40%.
!                         For ammonium to sulfate ratio less than 1.0
!                         all ammonium is aerosol and no nitrate aerosol
!                         exists.
!   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!                         allow gas-phase ammonia to exist.
!   F.Binkowski 7/26/95   Changed equilibrium constants to values from
!                         Kim et al. (1993)
!   F.Binkowski 6/27/96   Changed to new water format
!   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
!                         activity coefficients. The binary activity
!                         coefficients
!                         are the same as the previous version
!   F.Binkowski 8/1/97    Changed minimum sulfate from 0.0 to 1.0e-6 i.e.
!                         1 picogram per cubic meter
!   F.Binkowski 2/23/98   Changes to code made by Ingmar Ackermann to
!                         deal with precision problems on workstations 
!                         incorporated in to this version.  Also included
!                         are his improved descriptions of variables. 
!  F. Binkowski 8/28/98   changed logic as follows: 
!                         If iterations fail, initial values of nitrate
!                          are retained. 
!                         Total mass budgets are changed to account for gas
!                         phase returned.
!  F.Binkowski 10/01/98   Removed setting RATIO to 5 for low to 
!                         to zero sulfate sulfate case.
!  F.Binkowski 01/10/2000 reconcile versions
!
!  F.Binkowski 05/17/2000 change to logic for calculating RATIO
!  F.Binkowski 04/09/2001 change for very low values of RATIO,
!                         RATIO < 0.5, no iterative calculations are done
!                         in low ammonia case a MAX(1.0e-10, MSO4) IS
!                         applied, and the iteration count is
!                         reduced to fifty for each iteration loop.
!  R. Yantosca 09/25/2002 Bundled into "rpmares_mod.f".  Declared all REALs
!                          as TYPE (XPLEX)'s.  Cleaned up comments.  Also now force
!                          TYPE (XPLEX) explicitly with "D" exponents.
!  P. Le Sager and        Bug fix for low ammonia case -- prevent floating
!  R. Yantosca 04/10/2008  point underflow and NaN's.
!  P. Le Sager 06/10/2008 Better catch of over/underflow for low ammonia case
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP, IS_SAFE_DIV
      ! (dkh, 9/10/04, 09/07/09) 
      USE ERROR_MOD, ONLY   : ERROR_STOP
      USE CHECKPT_MOD, ONLY : nitr_max
      USE CHECKPT_MOD, ONLY : gamaan_fwd
      USE CHECKPT_MOD, ONLY : gamold_fwd
      USE CHECKPT_MOD, ONLY : wh2o_fwd
      USE CHECKPT_MOD, ONLY : ynh4_fwd
      USE CHECKPT_MOD, ONLY : eror_fwd
      USE CHECKPT_MOD, ONLY : exit_fwd
      USE CHECKPT_MOD, ONLY : gamana_fwd
      USE CHECKPT_MOD, ONLY : gamas1_fwd
      USE CHECKPT_MOD, ONLY : gamas2_fwd
      USE CHECKPT_MOD, ONLY : NNNMAX

      USE LOGICAL_ADJ_MOD, ONLY : LADJ

      !=================================================================
      ! ARGUMENTS and their descriptions
      !=================================================================
      TYPE (XPLEX) :: SO4              ! Total sulfate in micrograms / m**3
      TYPE (XPLEX) :: GNO3             ! Gas-phase nitric acid in micrograms / m**3
      TYPE (XPLEX) :: GNH3             ! Gas-phase ammonia in micrograms / m**3 
      TYPE (XPLEX) :: RH               ! Fractional relative humidity
      TYPE (XPLEX) :: TEMP             ! Temperature in Kelvin
      TYPE (XPLEX) :: ASO4             ! Aerosol sulfate in micrograms / m**3
      TYPE (XPLEX) :: AHSO4            ! Aerosol bisulfate in micrograms / m**3
      TYPE (XPLEX) :: ANO3             ! Aerosol nitrate in micrograms / m**3
      TYPE (XPLEX) :: AH2O             ! Aerosol liquid water content water in
                                 !   micrograms / m**3
      TYPE (XPLEX) :: ANH4             ! Aerosol ammonium in micrograms / m**3

      ! (dkh, 9/10/04, 09/07/09) 
      INTEGER :: I, J, L         ! Grid cell coords
      INTEGER EXIT

      !=================================================================
      ! PARAMETERS and their descriptions:
      !=================================================================

      ! Molecular weights
      TYPE (XPLEX), PARAMETER :: MWNACL = xplex(58.44277d0,0d0)               ! NaCl
      TYPE (XPLEX), PARAMETER :: MWNO3  = xplex(62.0049d0,0d0)                ! NO3
      TYPE (XPLEX), PARAMETER :: MWHNO3 = xplex(63.01287d0,0d0)               ! HNO3
      TYPE (XPLEX), PARAMETER :: MWSO4  = xplex(96.0576d0,0d0)                ! SO4
      TYPE (XPLEX), PARAMETER :: MWHSO4 = xplex(MWSO4%r+1.0080d0,0d0)         ! HSO4
      TYPE (XPLEX), PARAMETER :: MH2SO4 = xplex(98.07354d0,0d0)               ! H2SO4
      TYPE (XPLEX), PARAMETER :: MWNH3  = xplex(17.03061d0,0d0)               ! NH3
      TYPE (XPLEX), PARAMETER :: MWNH4  = xplex(18.03858d0,0d0)               ! NH4
      TYPE (XPLEX), PARAMETER :: MWORG  = xplex(16.0d0,0d0)                   ! Organic Species
      TYPE (XPLEX), PARAMETER :: MWCL   = xplex(35.453d0,0d0)                 ! Chloride
      TYPE (XPLEX), PARAMETER :: MWAIR  = xplex(28.964d0,0d0)                 ! AIR
      TYPE (XPLEX), PARAMETER :: MWLCT  = xplex(3.0d0 * MWNH4%r +          ! Letovicite
     &                              2.0d0 * MWSO4%r + 1.0080d0,0d0)  
      TYPE (XPLEX), PARAMETER :: MWAS=xplex(2.0d0*MWNH4%r+MWSO4%r,0d0)    ! Amm. Sulfate
      TYPE (XPLEX), PARAMETER::MWABS=xplex(MWNH4%r+MWSO4%r+1.0080d0,0d0) ! Amm. Bisulfate

      ! Minimum value of sulfate aerosol concentration
      TYPE (XPLEX), PARAMETER :: MINSO4 = xplex(1.0d-6 / MWSO4%r,0d0) 

      ! Minimum total nitrate cncentration
      TYPE (XPLEX), PARAMETER :: MINNO3 = xplex(1.0d-6 / MWNO3%r,0d0)  

      ! Force a minimum concentration
      TYPE (XPLEX), PARAMETER :: FLOOR  = xplex(1.0d-30,0d0) 

      ! Tolerances for convergence test.  NOTE: We now have made these
      ! parameters so they don't lose their values (phs, bmy, 4/10/08)
      TYPE (XPLEX), PARAMETER :: TOLER1 = xplex(0.00001d0,0d0)       
      TYPE (XPLEX), PARAMETER :: TOLER2 = xplex(0.001d0,0d0)       

      !=================================================================
      ! SCRATCH LOCAL VARIABLES and their descriptions:
      !=================================================================

      INTEGER :: IRH              ! Index set to percent relative humidity
      INTEGER :: NITR             ! Number of iterations for activity
                                  !   coefficients
      INTEGER :: NNN              ! Loop index for iterations
      INTEGER :: NR               ! Number of roots to cubic equation for
                                  ! H+ ciaprecision
      TYPE (XPLEX)  :: A0               ! Coefficients and roots of
      TYPE (XPLEX)  :: A1               ! Coefficients and roots of
      TYPE (XPLEX)  :: A2               ! Coefficients and roots of
      TYPE (XPLEX)  :: AA               ! Coefficients and discriminant for
                                  ! quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: BAL              ! internal variables ( high ammonia case)
      TYPE (XPLEX)  :: BB               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: BHAT             ! Variables used for ammonia solubility
      TYPE (XPLEX)  :: CC               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: CONVT            ! Factor for conversion of units
      TYPE (XPLEX)  :: DD               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: DISC             ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      TYPE (XPLEX)  :: EROR             ! Relative error used for convergence test
      TYPE (XPLEX)  :: FNH3             ! "Free ammonia concentration", that
                                  !   which exceeds TWOSO4
      TYPE (XPLEX)  :: GAMAAB           ! Activity Coefficient for (NH4+,
                                  !   HSO4-)GAMS( 2,3 )
      TYPE (XPLEX)  :: GAMAAN           ! Activity coefficient for (NH4+, NO3-)
                                  !   GAMS( 2,2 )
      TYPE (XPLEX)  :: GAMAHAT          ! Variables used for ammonia solubility
      TYPE (XPLEX)  :: GAMANA           ! Activity coefficient for (H+ ,NO3-)
                                  !   GAMS( 1,2 )
      TYPE (XPLEX)  :: GAMAS1           ! Activity coefficient for (2H+, SO4--)
                                  !   GAMS( 1,1 )
      TYPE (XPLEX)  :: GAMAS2           ! Activity coefficient for (H+, HSO4-)
                                  !   GAMS( 1,3 )
      TYPE (XPLEX)  :: GAMOLD           ! used for convergence of iteration
      TYPE (XPLEX)  :: GASQD            ! internal variables ( high ammonia case)
      TYPE (XPLEX)  :: HPLUS            ! Hydrogen ion (low ammonia case) (moles
                                  !   / kg water)
      TYPE (XPLEX)  :: K1A              ! Equilibrium constant for ammonia to
                                  !   ammonium
      TYPE (XPLEX)  :: K2SA             ! Equilibrium constant for
                                  !   sulfate-bisulfate (aqueous)
      TYPE (XPLEX)  :: K3               ! Dissociation constant for ammonium
                                  !   nitrate
      TYPE (XPLEX)  :: KAN              ! Equilibrium constant for ammonium
                                  !   nitrate (aqueous)
      TYPE (XPLEX)  :: KHAT             ! Variables used for ammonia solubility
      TYPE (XPLEX)  :: KNA              ! Equilibrium constant for nitric acid
                                  !   (aqueous)
      TYPE (XPLEX)  :: KPH              ! Henry's Law Constant for ammonia
      TYPE (XPLEX)  :: KW               ! Equilibrium constant for water
                                  !  dissociation
      TYPE (XPLEX)  :: KW2              ! Internal variable using KAN
      TYPE (XPLEX)  :: MAN              ! Nitrate (high ammonia case) (moles /
                                  !   kg water)
      TYPE (XPLEX)  :: MAS              ! Sulfate (high ammonia case) (moles /
                                  !   kg water)
      TYPE (XPLEX)  :: MHSO4            ! Bisulfate (low ammonia case) (moles /
                                  !   kg water)
      TYPE (XPLEX)  :: MNA              ! Nitrate (low ammonia case) (moles / kg
                                  !   water)
      TYPE (XPLEX)  :: MNH4             ! Ammonium (moles / kg water)
      TYPE (XPLEX)  :: MOLNU            ! Total number of moles of all ions
      TYPE (XPLEX)  :: MSO4             ! Sulfate (low ammonia case) (moles / kg
                                  !   water)
      TYPE (XPLEX)  :: PHIBAR           ! Practical osmotic coefficient
      TYPE (XPLEX)  :: PHIOLD           ! Previous value of practical osmotic
                                  !   coefficient used for convergence of
                                  !   iteration
      TYPE (XPLEX)  :: RATIO            ! Molar ratio of ammonium to sulfate
      TYPE (XPLEX)  :: RK2SA            ! Internal variable using K2SA
      TYPE (XPLEX)  :: RKNA             ! Internal variables using KNA
      TYPE (XPLEX)  :: RKNWET           ! Internal variables using KNA
      TYPE (XPLEX)  :: RR1
      TYPE (XPLEX)  :: RR2
      TYPE (XPLEX)  :: STION            ! Ionic strength
      TYPE (XPLEX)  :: T1               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T2               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T21              ! Internal variables of convenience (low
                                  !   ammonia case)
      TYPE (XPLEX)  :: T221             ! Internal variables of convenience (low
                                  !   ammonia case)
      TYPE (XPLEX)  :: T3               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T4               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: T6               ! Internal variables for temperature
                                  !   corrections
      TYPE (XPLEX)  :: TNH4             ! Total ammonia and ammonium in
                                  !   micromoles / meter ** 3
      TYPE (XPLEX)  :: TNO3             ! Total nitrate in micromoles / meter ** 3
      TYPE (XPLEX)  :: TSO4             ! Total sulfate in micromoles / meter ** 3
      TYPE (XPLEX)  :: TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles
                                  !   / kg water)
      TYPE (XPLEX)  :: WFRAC            ! Water mass fraction
      TYPE (XPLEX)  :: WH2O             ! Aerosol liquid water content (internally)
                                  !   micrograms / meter **3 on output
                                  !   internally it is 10 ** (-6) kg (water)
                                  !   / meter ** 3
                                  !   the conversion factor (1000 g = 1 kg)
                                  !   is applied for AH2O output
      TYPE (XPLEX)  :: WSQD             ! internal variables ( high ammonia case)
      TYPE (XPLEX)  :: XNO3             ! Nitrate aerosol concentration in
                                  ! micromoles / meter ** 3
      TYPE (XPLEX)  :: XXQ              ! Variable used in quadratic solution
      TYPE (XPLEX)  :: YNH4             ! Ammonium aerosol concentration in
                                  !  micromoles / meter** 3
      TYPE (XPLEX)  :: ZH2O             ! Water variable saved in case ionic
                                  !  strength too high.
      TYPE (XPLEX)  :: ZSO4             ! Total sulfate molality - mso4 + mhso4
                                  !  (low ammonia case) (moles / kg water)
      TYPE (XPLEX)  :: CAT( 2 )         ! Array for cations (1, H+); (2, NH4+)
                                  !  (moles / kg water)
      TYPE (XPLEX)  :: AN ( 3 )         ! Array for anions (1, SO4--); (2,
                                  !   NO3-); (3, HSO4-)  (moles / kg water)
      TYPE (XPLEX)  :: CRUTES( 3 )      ! Coefficients and roots of
      TYPE (XPLEX)  :: GAMS( 2, 3 )     ! Array of activity coefficients
      TYPE (XPLEX)  :: TMASSHNO3        ! Total nitrate (vapor and particle) 
      TYPE (XPLEX)  :: GNO3_IN, ANO3_IN                
      
      !=================================================================
      ! RPMARES_FORADJ begins here!
      ! convert into micromoles/m**3
      !=================================================================

      ! Initialize EXIT to 0 (dkh, 9/10/04, 09/07/09) 
      EXIT = 0

      ! For extremely low relative humidity ( less than 1% ) set the 
      ! water content to a minimum and skip the calculation.
      IF ( RH .LT. 0.01 ) THEN
         AH2O = FLOOR

         ! (dkh, 9/10/04, 09/07/09)  
         EXIT = 1

         RETURN
      ENDIF 

      ! total sulfate concentration
      TSO4 = MAX( FLOOR, SO4 / MWSO4  )      
      ASO4 = SO4

      !Cia models3 merge NH3/NH4 , HNO3,NO3 here
      !c *** recommended by Dr. Ingmar Ackermann

      ! total nitrate
      TNO3      = MAX( 0.0d0, ( ANO3 / MWNO3 + GNO3 / MWHNO3 ) )            

      ! total ammonia
      TNH4      = MAX( 0.0d0, ( GNH3 / MWNH3 + ANH4 / MWNH4 )  )

      GNO3_IN   = GNO3
      ANO3_IN   = ANO3
      TMASSHNO3 = MAX( 0.0d0, GNO3 + ANO3 )     

      ! set the  molar ratio of ammonium to sulfate
      RATIO = TNH4 / TSO4

      ! validity check for negative concentration
      IF ( TSO4 < 0.0d0 .OR. TNO3 < 0.0d0 .OR. TNH4 < 0.0d0 ) THEN
          PRINT*, 'TSO4 : ', TSO4
          PRINT*, 'TNO3 : ', TNO3
          PRINT*, 'TNH4 : ', TNH4
          CALL GEOS_CHEM_STOP
      ENDIF   

      ! now set humidity index IRH as a percent
      IRH = NINT( 100.0 * RH )

      ! now set humidity index IRH as a percent
      IRH = MAX(  1, IRH )
      IRH = MIN( 99, IRH )

      !=================================================================
      ! Specify the equilibrium constants at  correct temperature.  
      ! Also change units from ATM to MICROMOLE/M**3 (for KAN, KPH, and K3 )
      ! Values from Kim et al. (1993) except as noted.
      ! Equilibrium constant in Kim et al. (1993)
      !   K = K0 exp[ a(T0/T -1) + b(1+log(T0/T)-T0/T) ], T0 = 298.15 K
      !   K = K0 EXP[ a T3 + b T4 ] in the code here.
      !=================================================================
      CONVT = 1.0d0 / ( 0.082d0 * TEMP )
      T6    = 0.082d-9 *  TEMP
      T1    = 298.0d0 / TEMP
      T2    = LOG( T1 )
      T3    = T1 - 1.0d0
      T4    = 1.0d0 + T2 - T1

      !=================================================================
      ! Equilibrium Relation
      ! 
      ! HSO4-(aq)         = H+(aq)   + SO4--(aq)  ; K2SA
      ! NH3(g)            = NH3(aq)               ; KPH
      ! NH3(aq) + H2O(aq) = NH4+(aq) + OH-(aq)    ; K1A
      ! HNO3(g)           = H+(aq)   + NO3-(aq)   ; KNA
      ! NH3(g) + HNO3(g)  = NH4NO3(s)             ; K3
      ! H2O(aq)           = H+(aq)   + OH-(aq)    ; KW
      !=================================================================
      KNA  = 2.511d+06 *  EXP(  29.17d0 * T3 + 16.83d0 * T4 ) * T6
      K1A  = 1.805d-05 *  EXP(  -1.50d0 * T3 + 26.92d0 * T4 )
      K2SA = 1.015d-02 *  EXP(   8.85d0 * T3 + 25.14d0 * T4 )
      KW   = 1.010d-14 *  EXP( -22.52d0 * T3 + 26.92d0 * T4 )
      KPH  = 57.639d0  *  EXP(  13.79d0 * T3 -  5.39d0 * T4 ) * T6
      !K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
      KHAT =  KPH * K1A / KW
      KAN  =  KNA * KHAT

      ! Compute temperature dependent equilibrium constant for NH4NO3
      ! (from Mozurkewich, 1993)
      K3 = EXP( 118.87d0  - 24084.0d0 / TEMP -  6.025d0  * LOG( TEMP ) )

      ! Convert to (micromoles/m**3) **2
      K3     = K3 * CONVT * CONVT

      WH2O   = 0.0d0
      STION  = 0.0d0
      AH2O   = 0.0d0
      MAS    = 0.0d0
      MAN    = 0.0d0
      HPLUS  = 0.0d0
      NITR   = 0
      NR     = 0
      GAMAAN = 1.0d0
      GAMOLD = 1.0d0

      ! If there is very little sulfate and  nitrate 
      ! set concentrations to a very small value and return.
      IF ( ( TSO4 .LT. MINSO4 ) .AND. ( TNO3 .LT. MINNO3 ) ) THEN
         ASO4  = MAX( FLOOR, ASO4  )
         AHSO4 = MAX( FLOOR, AHSO4 ) ! [rjp, 12/12/01]
         ANO3  = MAX( FLOOR, ANO3  )
         ANH4  = MAX( FLOOR, ANH4  )
         WH2O  = FLOOR
         AH2O  = FLOOR
         GNH3  = MAX( FLOOR, GNH3  )
         GNO3  = MAX( FLOOR, GNO3  )
         
         ! (dkh, 9/10/04, 09/07/09) 
         EXIT = 2 

         RETURN
      ENDIF

      !=================================================================
      ! High Ammonia Case
      !=================================================================
      IF ( RATIO .GT. 2.0d0 ) THEN
        
         GAMAAN = 0.1d0

         ! Set up twice the sulfate for future use.
         TWOSO4 = 2.0d0 * TSO4
         XNO3   = 0.0d0
         YNH4   = TWOSO4

         ! Treat different regimes of relative humidity
         !
         ! ZSR relationship is used to set water levels. Units are
         !  10**(-6) kg water/ (cubic meter of air)
         !  start with ammomium sulfate solution without nitrate
       
         CALL AWATER( IRH, TSO4, YNH4, TNO3, AH2O ) !**** note TNO3
         WH2O = 1.0d-3 * AH2O
         ASO4 = TSO4   * MWSO4

         ! In sulfate poor case, Sulfate ion is preferred
         ! Set bisulfate equal to zero [rjp, 12/12/01]
         AHSO4 = 0.0d0
         ANO3  = 0.0d0
         ANH4  = YNH4 * MWNH4
         WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )

        ! This statement was commented out long ago, so 
        ! no need to replace the RETURN with an EXIT = 1
        !IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
        IF ( WFRAC .LT. 0.2d0 ) THEN

           ! "dry" ammonium sulfate and ammonium nitrate
           ! compute free ammonia 
           FNH3 = TNH4 - TWOSO4
           CC   = TNO3 * FNH3 - K3

           ! check for not enough to support aerosol
           IF ( CC .LE. 0.0d0 ) THEN
              XNO3 = 0.0d0
           ELSE
              AA   = 1.0d0
              BB   = -( TNO3 + FNH3 )
              DISC = BB * BB - 4.0d0 * CC

              ! Check for complex roots of the quadratic
              ! set retain initial values of nitrate and RETURN 
              ! if complex roots are found
              IF ( DISC .LT. 0.0d0 ) THEN
                 XNO3  = 0.0d0
                 AH2O  = 1000.0d0 * WH2O
                 YNH4  = TWOSO4
                 ASO4  = TSO4 * MWSO4
                 AHSO4 = 0.0d0
                 ANH4  = YNH4 * MWNH4
                 GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
                 GNO3  = GNO3_IN
                 ANO3  = ANO3_IN

                 ! (dkh, 9/10/04, 09/07/09) 
                 EXIT = 3 

                 RETURN
              ENDIF

              ! to get here, BB .lt. 0.0, CC .gt. 0.0 always
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN ( 1.0d0, BB ) * DD )


              ! Since both roots are positive, select smaller root.
              XNO3 = MIN( XXQ / AA, CC / XXQ )

           ENDIF                ! CC .LE. 0.0

           AH2O  = 1000.0d0 * WH2O
           YNH4  = TWOSO4 + XNO3          
           ASO4  = TSO4 * MWSO4
           AHSO4 = FLOOR
           ANO3  = XNO3 * MWNO3
           ANH4  = YNH4 * MWNH4
           GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 )  )      
           GNO3  = MAX( FLOOR, ( TMASSHNO3 - ANO3 ) )    
 
           ! (dkh, 9/10/04, 09/07/09) 
           EXIT = 4 

           RETURN
        ENDIF                  ! WFRAC .LT. 0.2

        ! liquid phase containing completely neutralized sulfate and
        ! some nitrate.  Solve for composition and quantity.
        MAS    = TSO4 / WH2O
        MAN    = 0.0d0
        XNO3   = 0.0d0
        YNH4   = TWOSO4
        PHIOLD = 1.0d0

        !===============================================================
        ! Start loop for iteration
        !
        ! The assumption here is that all sulfate is ammonium sulfate,
        ! and is supersaturated at lower relative humidities.
        !===============================================================
        ! Define max iterations as NNNMAX (dkh, 9/10/04, 09/07/09) 
        !DO NNN = 1, 50 ! loop count reduced 0409/2001 by FSB
        DO NNN = 1, NNNMAX ! loop count reduced 0409/2001 by FSB

           NITR  = NNN
           GASQD = GAMAAN * GAMAAN
           WSQD  = WH2O * WH2O
           KW2   = KAN * WSQD / GASQD
           AA    = 1.0 - KW2
           BB    = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
           CC    = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

           ! This is a quadratic for XNO3 [MICROMOLES / M**3] 
           ! of nitrate in solution
           DISC = BB * BB - 4.0d0 * AA * CC

           ! Check for complex roots, retain inital values and RETURN
           IF ( DISC .LT. 0.0 ) THEN
              XNO3  = 0.0d0
              AH2O  = 1000.0d0 * WH2O
              YNH4  = TWOSO4
              ASO4  = TSO4 * MWSO4
              AHSO4 = FLOOR     ! [rjp, 12/12/01]
              ANH4  = YNH4 * MWNH4
              GNH3  = MWNH3 * MAX( FLOOR, (TNH4 - YNH4 ) )
              GNO3  = GNO3_IN
              ANO3  = ANO3_IN

              ! (dkh, 9/10/04, 09/07/09) 
              EXIT = 5

              RETURN
           ENDIF
           
           ! Deal with degenerate case (yoj)
           IF ( AA .NE. 0.0d0 ) THEN
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN( 1.0d0, BB ) * DD )
              RR1 = XXQ / AA
              RR2 = CC / XXQ

              ! choose minimum positve root
              IF ( ( RR1 * RR2 ) .LT. 0.0d0 ) THEN
                 XNO3 = MAX( RR1, RR2 )
              ELSE
                 XNO3 = MIN( RR1, RR2 )
              ENDIF
           ELSE
              XNO3 = - CC / BB  ! AA equals zero here.
           ENDIF

           XNO3 = MIN( XNO3, TNO3 )
           
           ! This version assumes no solid sulfate forms (supersaturated )
           ! Now update water
           CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

           ! ZSR relationship is used to set water levels. Units are
           ! 10**(-6) kg water/ (cubic meter of air).  The conversion 
           ! from micromoles to moles is done by the units of WH2O.
           WH2O = 1.0d-3 * AH2O 

           ! Ionic balance determines the ammonium in solution.
           MAN  = XNO3 / WH2O
           MAS  = TSO4 / WH2O
           MNH4 = 2.0d0 * MAS + MAN
           YNH4 = MNH4 * WH2O

           ! MAS, MAN and MNH4 are the aqueous concentrations of sulfate, 
           ! nitrate, and ammonium in molal units (moles/(kg water) ).
           STION    = 3.0d0 * MAS + MAN
           CAT( 1 ) = 0.0d0
           CAT( 2 ) = MNH4
           AN ( 1 ) = MAS
           AN ( 2 ) = MAN
           AN ( 3 ) = 0.0d0
           CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )
           GAMAAN = GAMS( 2, 2 )

           ! Use GAMAAN for convergence control
           EROR   = ABS( GAMOLD - GAMAAN ) / GAMOLD
           GAMOLD = GAMAAN

           ! Now check here for adjoint run (dkh, 02/12/12, adj32_001) 
           IF ( LADJ ) THEN 

              !==================================================================
              ! CHECKPOINT (dkh, 9/10/04, 09/07/09) 
              ! Save variables gamaan,wh2o,ynh4,gamold and exit for adjoint 
              ! calculation
              !==================================================================
              gamaan_fwd (I,J,L,NNN) = GAMAAN
              wh2o_fwd   (I,J,L,NNN) = WH2O
              ynh4_fwd   (I,J,L,NNN) = YNH4
              gamold_fwd (I,J,L,NNN) = GAMOLD
              exit_fwd   (I,J,L,NNN) = EXIT

           ENDIF 

           ! Check to see if we have a solution
           IF ( EROR .LE. TOLER1 ) THEN
              ASO4  = TSO4 * MWSO4
              AHSO4 = 0.0d0       ! [rjp, 12/12/01]
              ANO3  = XNO3 * MWNO3
              ANH4  = YNH4 * MWNH4
              GNO3  = MAX( FLOOR, ( TMASSHNO3  - ANO3 ) )
              GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
              AH2O  = 1000.0d0 * WH2O

              ! (dkh, 9/10/04, 09/07/09) 
              EXIT = 6 
              ! check for out of bounds NITR
              IF (NITR > NNNMAX) THEN
                 print*, I, J, L, nitr, nnn
                 CALL ERROR_STOP ( ' NITR > NNNMAX ',
     &                              'RPMARES_FORADJ' )
              ENDIF
 
              ! Now check here for adjoint run (dkh, 02/12/12, adj32_001) 
              IF ( LADJ ) THEN 

                 ! Save more variables for adjoint calc
                 nitr_max (I,J,L)     = NITR
                 exit_fwd (I,J,L,NNN) = exit

              ENDIF 

              RETURN
           ENDIF

        ENDDO

        ! If after NITR iterations no solution is found, then:
        ! FSB retain the initial values of nitrate particle and vapor
        ! note whether or not convert all bisulfate to sulfate
        ASO4  = TSO4 * MWSO4
        AHSO4 = FLOOR      
        XNO3  = TNO3 / MWNO3
        YNH4  = TWOSO4
        ANH4  = YNH4 * MWNH4

        CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

        GNO3  = GNO3_IN
        ANO3  = ANO3_IN
        GNH3  = MAX( FLOOR, MWNH3 * (TNH4 - YNH4 ) )

        ! (dkh, 9/10/04, 09/07/09) 
        EXIT = 7 

        RETURN

      !================================================================
      ! Low Ammonia Case 
      !
      ! Coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98
      ! modified 04/09/2001
      !       
      ! All cases covered by this logic
      !=================================================================
      ELSE

         WH2O = 0.0d0
         CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
         WH2O = 1.0d-3 * AH2O
         ZH2O = AH2O

         ! convert 10**(-6) kg water/(cubic meter of air) to micrograms 
         ! of water per cubic meter of air (1000 g = 1 kg)
         ! in sulfate rich case, preferred form is HSO4-
         !ASO4 = TSO4 * MWSO4
         ASO4  = FLOOR          ![rjp, 12/12/01]
         AHSO4 = TSO4 * MWSO4   ![rjp, 12/12/01]
         ANH4  = TNH4 * MWNH4
         ANO3  = ANO3_IN
         GNO3  = TMASSHNO3 - ANO3
         GNH3  = FLOOR 
        
         !==============================================================
         ! *** Examine special cases and return if necessary.
         !         
         ! FSB For values of RATIO less than 0.5 do no further 
         ! calculations.  The code will cycle and still predict the 
         ! same amount of ASO4, ANH4, ANO3, AH2O so terminate early 
         ! to swame computation
         !==============================================================
         ! (dkh, 9/10/04, 09/07/09) 
         !IF ( RATIO .LT. 0.5d0 ) RETURN ! FSB 04/09/2001 
         IF ( RATIO .LT. 0.5d0 ) THEN 
           
            EXIT = 8 

            RETURN ! FSB 04/09/2001 
         ENDIF 
    
         ! Check for zero water.
         ! (dkh, 9/10/04, 09/07/09) 
         !IF ( WH2O .EQ. 0.0d0 ) RETURN
         IF ( WH2O .EQ. 0.0d0 ) THEN
 
            EXIT = 9 

            RETURN
         ENDIF 

         ZSO4 = TSO4 / WH2O

         ! ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4
         ! do not solve for aerosol nitrate for total sulfate molality
         ! greater than 11.0 because the model parameters break down
         !### IF ( ZSO4 .GT. 11.0 ) THEN
         IF ( ZSO4 .GT. 9.0 ) THEN ! 18 June 97
 
            ! (dkh, 9/10/04, 09/07/09) 
            EXIT = 10 
            RETURN
         ENDIF

         ! *** Calculation may now proceed.
         !
         ! First solve with activity coeffs of 1.0, then iterate.
         PHIOLD = 1.0d0
         GAMANA = 1.0d0
         GAMAS1 = 1.0d0
         GAMAS2 = 1.0d0
         GAMAAB = 1.0d0
         GAMOLD = 1.0d0

         ! All ammonia is considered to be aerosol ammonium.
         MNH4 = TNH4 / WH2O

         ! MNH4 is the molality of ammonium ion.
         YNH4 = TNH4

         ! loop for iteration
         ! (dkh, 9/10/04, 09/07/09) 
         !DO NNN = 1, 50    ! loop count reduced 04/09/2001 by FSB
         DO NNN = 1, NNNMAX    ! loop count reduced 04/09/2001 by FSB
            NITR = NNN

            !------------------------------------------------------------
            ! Add robustness: now check if GAMANA or GAMAS1 is too small
            ! for the division in RKNA and RK2SA. If they are, return w/ 
            ! original values: basically replicate the procedure used 
            ! after the current DO-loop in case of no-convergence
            ! (phs, bmy, rjp, 4/10/08)
            ! Now uses IS_SAFE_DIV to avoid compiler/machine dependency 
            ! and to check for both underlow and overflow. Also 
            ! use REAL4 flag to avoid under/overflow when computing A0 
            ! and A1 from RKNA and RK2SA (phs, 5/28/08)
            !------------------------------------------------------------
            IF ( .NOT. (
     &           IS_SAFE_DIV( GAMAS2, GAMAS1*GAMAS1*GAMAS1, R4=.TRUE. )
     &           .AND. IS_SAFE_DIV( KNA, GAMANA*GAMANA, R4=.TRUE. ) 
     &           ) ) THEN
               
               WRITE(6,*) 'RPMARES: not safe to divide...exit'
               CALL flush(6)
               GOTO 100
            ENDIF
            
            ! set up equilibrium constants including activities
            ! solve the system for hplus first then sulfate & nitrate
            RK2SA  = K2SA * GAMAS2 * GAMAS2 / (GAMAS1 * GAMAS1 * GAMAS1)
            RKNA   = KNA / ( GAMANA * GAMANA )
            RKNWET = RKNA * WH2O
            T21    = ZSO4 - MNH4
            T221   = ZSO4 + T21

            ! set up coefficients for cubic
            A2 = RK2SA + RKNWET - T21
            A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )
     &           - RK2SA * ZSO4 - RKNA * TNO3
            A0 = - (T21 * RK2SA * RKNWET
     &           + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )

            CALL CUBIC ( A2, A1, A0, NR, CRUTES )

            ! Code assumes the smallest positive root is in CRUTES(1)
            ! But, it can be negative (see CUBIC, case of one real root, 
            ! but can also be propagated by over/underflown)... if it is
            ! the case then return with original values (phs, 5/27/08)
            HPLUS = CRUTES( 1 )
            ! (dkh, 9/10/04, 09/07/09) 
            ! reinstate (dkh, 02/12/12, adj32_001) 
            IF (HPLUS <= 0d0) GOTO 100
            BAL   = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0

            ! molality of sulfate ion
            MSO4  = RK2SA * ZSO4 / ( HPLUS + RK2SA ) 

            ! molality of bisulfate ion
            ! MAX added 04/09/2001 by FSB
            MHSO4 = MAX( 1.0d-10, ZSO4 - MSO4 ) 

            ! molality of nitrate ion
            MNA   = RKNA * TNO3 / ( HPLUS + RKNWET ) 
            MNA   = MAX( 0.0d0, MNA )
            MNA   = MIN( MNA, TNO3 / WH2O )
            XNO3  = MNA * WH2O
            ANO3  = MNA * WH2O * MWNO3
            GNO3  = MAX( FLOOR, TMASSHNO3 - ANO3 )
            ASO4  = MSO4 * WH2O * MWSO4 ![rjp, 12/12/01]
            AHSO4 = MHSO4 * WH2O * MWSO4 ![rjp, 12/12/01]

            ! Calculate ionic strength
            STION = 0.5d0 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0d0 * MSO4)

            ! Update water
            CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

            ! Convert 10**(-6) kg water/(cubic meter of air) to micrograms 
            ! of water per cubic meter of air (1000 g = 1 kg)
            WH2O     = 1.0d-3 * AH2O
            CAT( 1 ) = HPLUS
            CAT( 2 ) = MNH4
            AN ( 1 ) = MSO4
            AN ( 2 ) = MNA
            AN ( 3 ) = MHSO4

            CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )

            GAMANA = GAMS( 1, 2 )
            GAMAS1 = GAMS( 1, 1 )
            GAMAS2 = GAMS( 1, 3 )
            GAMAAN = GAMS( 2, 2 )

            ! NOTE: Improved for robustness!
            GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
            BHAT = KHAT * GAMAHAT
            !### EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
            !### PHIOLD = PHIBAR
            EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD
            GAMOLD = GAMAHAT

            ! Now check here for adjoint run (dkh, 02/12/12, adj32_001) 
            IF ( LADJ ) THEN

               !==================================================================
               ! CHECKPOINT (dkh, 9/10/04, 09/07/09) 
               ! Save variables error,gamana,gamas1,gamas2,gamold and wh2o for 
               ! the adjoint calculation
               !================================================================== 
               gamana_fwd(I,J,L,NNN) = GAMANA
               gamas1_fwd(I,J,L,NNN) = GAMAS1
               gamas2_fwd(I,J,L,NNN) = GAMAS2
               gamold_fwd(I,J,L,NNN) = GAMOLD
               wh2o_fwd(I,J,L,NNN)   = WH2O
               eror_fwd(I,J,L,NNN)   = EROR
               exit_fwd(I,J,L,NNN)   = EXIT

            ENDIF 

            ! return with good solution
            IF ( EROR .LE. TOLER2 ) THEN

               ! (dkh, 9/10/04, 09/07/09) 
               EXIT = 11

               ! Now check here for adjoint run (dkh, 02/12/12, adj32_001) 
               IF ( LADJ ) THEN

                  ! Save more variables for adjoint 
                  exit_fwd(I,J,L,NNN) = EXIT

                  nitr_max(I,J,L)     = NITR

               ENDIF 

               ! check for out of bounds NITR
               IF (NITR > NNNMAX) THEN
                  print*, I, J, L, nitr, nnn
                  CALL ERROR_STOP ( ' NITR > NNNMAX ',
     &                              'RPMARES_FORADJ' )                                   
               ENDIF

               RETURN
            ENDIF
            
         ENDDO

         ! after NITR iterations, failure to solve the system
         ! convert all ammonia to aerosol ammonium and return input
         ! values of NO3 and HNO3
 100     ANH4 = TNH4 * MWNH4
         GNH3 = FLOOR
         GNO3 = GNO3_IN
         ANO3 = ANO3_IN
         ASO4 = TSO4 * MWSO4    ! [rjp, 12/17/01]
         AHSO4= FLOOR           ! [rjp, 12/17/01]

         CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
        
         !  (dkh, 9/10/04, 09/07/09) 
         EXIT = 12 

         RETURN

      ENDIF                     ! ratio .gt. 2.0

      ! Return to calling program
      END SUBROUTINE RPMARES_FORADJ

!------------------------------------------------------------------------------

      SUBROUTINE AWATER( IRHX, MSO4, MNH4, MNO3, WH2O )
!
!******************************************************************************
! NOTE!!! wh2o is returned in micrograms / cubic meter
!         mso4,mnh4,mno3 are in microMOLES / cubic meter
!
!  This  version uses polynomials rather than tables, and uses empirical
! polynomials for the mass fraction of solute (mfs) as a function of water
! activity
!   where:
!
!            mfs = ms / ( ms + mw)
!             ms is the mass of solute
!             mw is the mass of water.
!
!  Define y = mw/ ms
!
!  then  mfs = 1 / (1 + y)
!
!    y can then be obtained from the values of mfs as
!
!             y = (1 - mfs) / mfs
!
!
!     the aerosol is assumed to be in a metastable state if the rh is
!     is below the rh of deliquescence, but above the rh of crystallization.
!
!     ZSR interpolation is used for sulfates with x ( the molar ratio of
!     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
!     section 1: 0 <= x < 1
!     section 2: 1 <= x < 1.5
!     section 3: 1.5 <= x < 2.0
!     section 4: 2 <= x
!     In sections 1 through 3, only the sulfates can affect the amount of water
!     on the particles.
!     In section 4, we have fully neutralized sulfate, and extra ammonium which
!     allows more nitrate to be present. Thus, the ammount of water is
!     calculated
!     using ZSR for ammonium sulfate and ammonium nitrate. Crystallization is
!     assumed to occur in sections 2,3,and 4. See detailed discussion below.
!
!
!
! definitions:
!     mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of air)
!      for sulfate, ammonium, and nitrate respectively
!     irhx is the relative humidity (%)
!     wh2o is the returned water amount in micrograms / cubic meter of air
!     x is the molar ratio of ammonium to sulfate
!     y0,y1,y1.5, y2 are the water contents in mass of water/mass of solute
!     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
!     y3 is the value of the mass ratio of water to solute for
!     a pure ammonium nitrate  solution.
!
!
!     coded by Dr. Francis S. Binkowski, 4/8/96.
!
! *** modified 05/30/2000 
!     The use of two values of mfs at an ammonium to sulfate ratio 
!     representative of ammonium sulfate led to an minor inconsistancy 
!     in nitrate behavior as the ratio went from a value less than two
!     to a value greater than two and vice versa with either ammonium 
!     held constant and sulfate changing, or sulfate held constant and 
!     ammonium changing. the value of Chan et al. (1992) is the only value
!     now used. 
!
! *** Modified 09/25/2002
!     Ported into "rpmares_mod.f".  Now declare all variables with TYPE (XPLEX).
!     Also cleaned up comments and made cosmetic changes.  Force double 
!     precision explicitly with "D" exponents. 
!******************************************************************************
!
      ! Arguments
      INTEGER           :: IRHX
      TYPE (XPLEX)            :: MSO4, MNH4, MNO3, WH2O

      ! Local variables
      INTEGER           :: IRH
      TYPE (XPLEX)          :: TSO4,  TNH4,  TNO3,  X,      AW,     AWC
      TYPE (XPLEX)            :: MFS0,  MFS1,  MFS15, Y 
      TYPE (XPLEX)          :: Y0,    Y1,    Y15,   Y2,     Y3,     Y40 
      TYPE (XPLEX)            :: Y140,  Y1540, YC,    MFSSO4, MFSNO3

      ! Molecular weight parameters
      TYPE (XPLEX), PARAMETER :: MWSO4  = xplex(96.0636d0,0d0)
      TYPE (XPLEX), PARAMETER :: MWNH4  = xplex(18.0985d0,0d0)
      TYPE (XPLEX), PARAMETER :: MWNO3  = xplex(62.0649d0,0d0)
      TYPE (XPLEX), PARAMETER :: MW2=xplex(MWSO4%r+2.0d0*MWNH4%r,0d0)
      TYPE (XPLEX), PARAMETER :: MWANO3 = xplex(MWNO3%r + MWNH4%r,0d0) 
      
      !=================================================================
      ! The polynomials use data for aw as a function of mfs from Tang 
      ! and Munkelwitz, JGR 99: 18801-18808, 1994.  The polynomials were 
      ! fit to Tang's values of water activity as a function of mfs.
      ! 
      ! *** coefficients of polynomials fit to Tang and Munkelwitz data
      !     now give mfs as a function of water activity.
      !=================================================================
      TYPE (XPLEX) :: C1(4)  = (/ xplex(0.9995178d0,0d0),
     &                            xplex(  -0.7952896d0,0d0), 
     &                      xplex(0.99683673d0,0d0),
     & xplex( -1.143874d0,0d0) /)

      TYPE (XPLEX) :: C15(4) = (/ xplex(1.697092d0,0d0),
     & xplex( -4.045936d0,0d0), 
     &                      xplex(5.833688d0,0d0),
     & xplex( -3.463783d0,0d0) /)
      
      TYPE (XPLEX) :: C2(4)  = (/ xplex(2.085067d0,0d0),
     & xplex( -6.024139d0,0d0), 
     &                      xplex(8.967967d0,0d0),
     & xplex( -5.002934d0,0d0) /)

      !=================================================================
      ! The following coefficients are a fit to the data in Table 1 of
      !    Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
      !      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
      !
      ! New data fit to data from
      !       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
      !       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
      !       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      !=================================================================
      TYPE (XPLEX)::C0(4)=(/xplex(0.798079d0,0d0),
     &                      xplex(-1.574367d0,0d0), 
     &                      xplex( 2.536686d0,0d0),
     &                      xplex( -1.735297d0,0d0) /)

      !=================================================================
      ! Polynomials for ammonium nitrate and ammonium sulfate are from:
      ! Chan et al.1992, Atmospheric Environment (26A): 1661-1673.
      !=================================================================
      TYPE (XPLEX) :: KNO3(6) = (/  xplex(0.2906d0,0d0),
     &   xplex(6.83665d0,0d0), xplex(-26.9093d0,0d0),
     &                       xplex(46.6983d0,0d0),
     &   xplex(-38.803d0,0d0),    xplex(11.8837d0,0d0) /)
      
      TYPE (XPLEX) :: KSO4(6) = (/xplex(2.27515d0,0d0),
     &  xplex(-11.147d0,0d0),  xplex(36.3369d0,0d0),
     &                       xplex(-64.2134d0,0d0),
     &    xplex(56.8341d0,0d0), xplex(-20.0953d0,0d0) /)

      !=================================================================
      ! AWATER begins here!
      !=================================================================
      
      ! Check range of per cent relative humidity
      IRH  = IRHX
      IRH  = MAX( 1, IRH )
      IRH  = MIN( IRH, 100 )

      ! Water activity = fractional relative humidity
      AW   = DBLE( IRH ) / 100.0d0
      TSO4 = MAX( MSO4 , 0.0d0 )
      TNH4 = MAX( MNH4 , 0.0d0 )
      TNO3 = MAX( MNO3 , 0.0d0 )
      X    = 0.0d0

      ! If there is non-zero sulfate calculate the molar ratio
      ! otherwise check for non-zero nitrate and ammonium
      IF ( TSO4 .GT. 0.0d0 ) THEN
         X = TNH4 / TSO4
      ELSE
         IF ( TNO3 .GT. 0.0d0 .AND. TNH4 .GT. 0.0d0 ) X = 10.0d0
      ENDIF

      ! *** begin screen on x for calculating wh2o
      IF ( X .LT. 1.0d0 ) THEN
         MFS0 = POLY4( C0, AW )
         MFS1 = POLY4( C1, AW )
         Y0   = ( 1.0d0 - MFS0 ) / MFS0
         Y1   = ( 1.0d0 - MFS1 ) / MFS1
         Y    = ( 1.0d0 - X    ) * Y0 + X * Y1

      ELSE IF ( X .LT. 1.5d0 ) THEN

         IF ( IRH .GE. 40 ) THEN
            MFS1  = POLY4( C1,  AW )
            MFS15 = POLY4( C15, AW )
            Y1    = ( 1.0d0 - MFS1  ) / MFS1
            Y15   = ( 1.0d0 - MFS15 ) / MFS15
            Y     = 2.0d0 * ( Y1 * ( 1.5d0 - X ) + Y15 *( X - 1.0d0 ) )

         !==============================================================
         ! Set up for crystalization
         !
         ! Crystallization is done as follows:
         !
         ! For 1.5 <= x, crystallization is assumed to occur 
         ! at rh = 0.4
         !
         ! For x <= 1.0, crystallization is assumed to occur at an 
         ! rh < 0.01, and since the code does not allow ar rh < 0.01, 
         ! crystallization is assumed not to occur in this range.
         !
         ! For 1.0 <= x <= 1.5 the crystallization curve is a straignt 
         ! line from a value of y15 at rh = 0.4 to a value of zero at 
         ! y1. From point B to point A in the diagram.  The algorithm 
         ! does a double interpolation to calculate the amount of
         ! water.
         !
         !        y1(0.40)               y15(0.40)
         !         +                     + Point B
         !
         !
         !
         !
         !         +--------------------+
         !       x=1                   x=1.5
         !      Point A
         !==============================================================
         ELSE

            ! rh along the crystallization curve.
            AWC = 0.80d0 * ( X - 1.0d0 ) 
            Y   = 0.0d0

            ! interpolate using crystalization curve
            IF ( AW .GE. AWC ) THEN 
               MFS1  = POLY4( C1,  xplx(0.40d0) )
               MFS15 = POLY4( C15, xplx(0.40d0) )
               Y140  = ( 1.0d0 - MFS1  ) / MFS1
               Y1540 = ( 1.0d0 - MFS15 ) / MFS15
               Y40   = 2.0d0 * ( Y140  * ( 1.5d0 - X ) + 
     &                           Y1540 * ( X - 1.0d0 ) )

               ! Y along crystallization curve
               YC   = 2.0d0 * Y1540 * ( X - 1.0d0 ) 
               Y    = Y40 - (Y40 - YC) * (0.40d0 - AW) / (0.40d0 - AWC)
            ENDIF              
         ENDIF                 

      ELSE IF ( X .LT. 2.0d0 ) then               ! changed 12/11/2000 by FSB
         Y = 0.0D0

         IF ( IRH .GE. 40 ) THEN 
            MFS15  = POLY4( C15, AW )            
            !MFS2  = POLY4( C2,  AW )
            Y15    = ( 1.0d0 - MFS15 ) / MFS15
            !y2    = ( 1.0d0 - MFS2  ) / MFS2
            MFSSO4 = POLY6( KSO4, AW )             ! Changed 05/30/2000 by FSB
            Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
            Y      = 2.0d0 * (Y15 * (2.0d0 - X) + Y2 * (X - 1.5d0) )
         ENDIF                  

      ELSE                                 ! 2.0 <= x changed 12/11/2000 by FSB

         !==============================================================
         ! Regime where ammonium sulfate and ammonium nitrate are 
         ! in solution.
         ! 
         ! following cf&s for both ammonium sulfate and ammonium nitrate
         ! check for crystallization here. their data indicate a 40% 
         ! value is appropriate.
         !==============================================================
         Y2 = 0.0d0
         Y3 = 0.0d0

         IF ( IRH .GE. 40 ) THEN
            MFSSO4 = POLY6( KSO4, AW )
            MFSNO3 = POLY6( KNO3, AW )
            Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
            Y3     = ( 1.0d0 - MFSNO3 ) / MFSNO3

         ENDIF

      ENDIF                     ! end of checking on x

      !=================================================================
      ! Now set up output of WH2O
      ! WH2O units are micrograms (liquid water) / cubic meter of air
      !=================================================================
      IF ( X .LT. 2.0D0 ) THEN  ! changed 12/11/2000 by FSB

         WH2O =  Y * ( TSO4 * MWSO4 + MWNH4 * TNH4 )

      ELSE

         ! this is the case that all the sulfate is ammonium sulfate
         ! and the excess ammonium forms ammonum nitrate
         WH2O =   Y2 * TSO4 * MW2 + Y3 * TNO3 * MWANO3

      ENDIF

      ! Return to calling program
      END SUBROUTINE AWATER

!------------------------------------------------------------------------------

      FUNCTION POLY4( A, X ) RESULT( Y )
      
      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: A(4), X

      ! Return value
      TYPE (XPLEX)             :: Y

      !=================================================================
      ! POLY4 begins here! 
      !=================================================================
      Y = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) )))

      ! Return to calling program
      END FUNCTION POLY4

!------------------------------------------------------------------------------

      FUNCTION POLY6( A, X ) RESULT( Y )

      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: A(6), X
      
      ! Return value
      TYPE (XPLEX)             :: Y

      !=================================================================
      ! POLY6 begins here! 
      !=================================================================
      Y = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) +
     &           X * ( A(5) + X * ( A(6)  )))))

      ! Return to calling program
      END FUNCTION POLY6

!------------------------------------------------------------------------------

      SUBROUTINE CUBIC( A2, A1, A0, NR, CRUTES )
!
!******************************************************************************
! Subroutine to find the roots of a cubic equation / 3rd order polynomial
! Formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   Dr. Francis S. Binkowski modified the routine on 6/24/91, 8/7/97
! ***
! *** modified 2/23/98 by fsb to incorporate Dr. Ingmar Ackermann's
!     recommendations for setting a0, a1,a2 as TYPE (XPLEX) variables.
!
! Modified by Bob Yantosca (10/15/02) 
! - Now use upper case / white space
! - force TYPE (XPLEX) with "D" exponents
! - updated comments / cosmetic changes
! - now call ERROR_STOP from "error_mod.f" to stop the run safely
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      INTEGER           :: NR
      TYPE (XPLEX)            :: A2, A1, A0
      TYPE (XPLEX)            :: CRUTES(3)

      ! Local variables
      TYPE (XPLEX)            :: QQ,    RR,    A2SQ,  THETA, DUM1, DUM2
      TYPE (XPLEX)            :: PART1, PART2, PART3, RRSQ,  PHI,  YY1
      TYPE (XPLEX)            :: YY2,   YY3,   COSTH, SINTH
      TYPE (XPLEX), PARAMETER :: ONE    = xplex(1.0d0,0d0)
      TYPE (XPLEX), PARAMETER :: SQRT3  = xplex(1.732050808d0,0d0)
      TYPE (XPLEX), PARAMETER :: ONE3RD = xplex(0.333333333d0,0d0)

      !=================================================================
      ! CUBIC begins here!
      !=================================================================
      A2SQ = A2 * A2
      QQ   = ( A2SQ - 3.d0*A1 ) / 9.d0
      RR   = ( A2*( 2.d0*A2SQ - 9.d0*A1 ) + 27.d0*A0 ) / 54.d0

      ! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      DUM1 = QQ * QQ * QQ
      RRSQ = RR * RR
      DUM2 = DUM1 - RRSQ

      IF ( DUM2 .GE. 0.d0 ) THEN

         ! Now we have three real roots
         PHI = SQRT( DUM1 )

         IF ( ABS( PHI ) .LT. 1.d-20 ) THEN
            CRUTES(1) = 0.0d0
            CRUTES(2) = 0.0d0
            CRUTES(3) = 0.0d0
            NR        = 0
            CALL ERROR_STOP( 'PHI < 1d-20', 'CUBIC (rpmares_mod.f)' )
         ENDIF
         
         THETA = ACOS( RR / PHI ) / 3.0d0
         COSTH = COS( THETA )
         SINTH = SIN( THETA )

         ! Use trig identities to simplify the expressions
         ! Binkowski's modification
         PART1     = SQRT( QQ )
         YY1       = PART1 * COSTH
         YY2       = YY1 - A2/3.0d0
         YY3       = SQRT3 * PART1 * SINTH
         CRUTES(3) = -2.0d0*YY1 - A2/3.0d0
         CRUTES(2) = YY2 + YY3
         CRUTES(1) = YY2 - YY3

         ! Set negative roots to a large positive value
         IF ( CRUTES(1) .LT. 0.0d0 ) CRUTES(1) = 1.0d9
         IF ( CRUTES(2) .LT. 0.0d0 ) CRUTES(2) = 1.0d9
         IF ( CRUTES(3) .LT. 0.0d0 ) CRUTES(3) = 1.0d9

         ! Put smallest positive root in crutes(1)
         CRUTES(1) = MIN( CRUTES(1), CRUTES(2), CRUTES(3) )
         NR        = 3

         !  Trap for negative CRUTES (dkh, 01/06/12, adj32_001)) 
         IF ( CRUTES(1) < 0d0 ) THEN 
            IF ( CRUTES(1) .ne. ABS(CRUTES(1) ) ) CRUTES(1)=1d9
            IF ( CRUTES(2) .ne. ABS(CRUTES(2) ) ) CRUTES(2)=1d9
            IF ( CRUTES(3) .ne. ABS(CRUTES(3) ) ) CRUTES(3)=1d9
            CRUTES(1) = MIN( CRUTES(1), CRUTES(2), CRUTES(3) )
            IF ( CRUTES(1) < 0d0 ) THEN 
               print*, ' CRUTES = ', CRUTES         , NR
               CALL ERROR_STOP( ' checking for neg val failed',
     &                       ' rpmares_mod.f ' )
            ENDIF 
         ENDIF 

      ELSE  

         ! Now here we have only one real root
         PART1     = SQRT( RRSQ - DUM1 )
         PART2     = ABS( RR )
         PART3     = ( PART1 + PART2 )**ONE3RD
         CRUTES(1) = -SIGN(ONE,RR) * ( PART3 + (QQ/PART3) ) - A2/3.D0
         CRUTES(2) = 0.D0
         CRUTES(3) = 0.D0
         NR        = 1

      ENDIF
     



      ! Return to calling program
      END SUBROUTINE CUBIC
      
!------------------------------------------------------------------------------

       SUBROUTINE ACTCOF( CAT, AN, GAMA, MOLNU, PHIMULT )
!
!******************************************************************************
!
! DESCRIPTION:
!
!  This subroutine computes the activity coefficients of (2NH4+,SO4--),
!  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
!  multicomponent solution, using Bromley's model and Pitzer's method.
!
! REFERENCES:
!
!   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
!     in aqueous solutions.  AIChE J. 19, 313-320.
!
!   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of
!     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
!
!   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures
!     of strong acids over saline solutions - I HNO3,
!     Atmos. Environ. (22): 91-100
!
!   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
!     and mean activity and osmotic coefficients of 0-100% nitric acid
!     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
!
!   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
!     general equilibrium model for inorganic multicomponent atmospheric
!     aerosols.  Atmos. Environ. 21(11), 2453-2466.
!
!
!
!
! ARGUMENT DESCRIPTION:
!
!     CAT(1) : conc. of H+    (moles/kg)
!     CAT(2) : conc. of NH4+  (moles/kg)
!     AN(1)  : conc. of SO4-- (moles/kg)
!     AN(2)  : conc. of NO3-  (moles/kg)
!     AN(3)  : conc. of HSO4- (moles/kg)
!     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
!     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
!     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
!     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
!     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
!     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
!     MOLNU   : the total number of moles of all ions.
!     PHIMULT : the multicomponent paractical osmotic coefficient.
!
! REVISION HISTORY:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
!                         new routine using a method described by Pilinis
!                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
!   S.Roselle   7/30/97   Modified for use in Models-3
!   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
!   R.Yantosca  9/25/02   Ported into "rpmares_mod.f" for GEOS-CHEM.  Cleaned
!                         up comments, etc.  Also force TYPE (XPLEX) by
!                         declaring REALs as TYPE (XPLEX) and by using "D" exponents.
!******************************************************************************
!

      ! Error codes


     
      !=================================================================
      ! PARAMETERS and their descriptions:
      !=================================================================
      INTEGER, PARAMETER :: NCAT = 2         ! number of cation
      INTEGER, PARAMETER :: NAN  = 3         ! number of anions
      TYPE (XPLEX),  PARAMETER :: XSTAT0 = xplex(0,0d0)       ! Normal, successful completion
      TYPE (XPLEX),  PARAMETER :: XSTAT1 = xplex(1,0d0)       ! File I/O error
      TYPE (XPLEX),  PARAMETER :: XSTAT2 = xplex(2,0d0)       ! Execution error
      TYPE (XPLEX),  PARAMETER :: XSTAT3 = xplex(3,0d0)       ! Special  error

      !=================================================================
      ! ARGUMENTS and their descriptions
      !=================================================================
      TYPE (XPLEX)             :: MOLNU            ! tot # moles of all ions
      TYPE (XPLEX)             :: PHIMULT          ! multicomponent paractical 
                                             !   osmotic coef
      TYPE (XPLEX)             :: CAT(NCAT)        ! cation conc in moles/kg (input)
      TYPE (XPLEX)             :: AN(NAN)          ! anion conc in moles/kg (input)
      TYPE (XPLEX)             :: GAMA(NCAT,NAN)   ! mean molal ionic activity coefs

      !=================================================================
      ! SCRATCH LOCAL VARIABLES and their descriptions:
      !=================================================================
      INTEGER            :: IAN              ! anion indX
      INTEGER            :: ICAT             ! cation indX
      TYPE (XPLEX)             :: FGAMA            !
      TYPE (XPLEX)             :: I                ! ionic strength
      TYPE (XPLEX)             :: R                !
      TYPE (XPLEX)             :: S                !
      TYPE (XPLEX)             :: TA               !
      TYPE (XPLEX)             :: TB               !
      TYPE (XPLEX)             :: TC               !
      TYPE (XPLEX)             :: TEXPV            !
      TYPE (XPLEX)             :: TRM              !
      TYPE (XPLEX)             :: TWOI             ! 2*ionic strength
      TYPE (XPLEX)             :: TWOSRI           ! 2*sqrt of ionic strength
      TYPE (XPLEX)             :: ZBAR             !
      TYPE (XPLEX)             :: ZBAR2            !
      TYPE (XPLEX)             :: ZOT1             !
      TYPE (XPLEX)             :: SRI              ! square root of ionic strength
      TYPE (XPLEX)             :: F2(NCAT)         !
      TYPE (XPLEX)             :: F1(NAN)          !
      TYPE (XPLEX)             :: BGAMA (NCAT,NAN) !
      TYPE (XPLEX)             :: X     (NCAT,NAN) !
      TYPE (XPLEX)             :: M     (NCAT,NAN) ! molality of each electrolyte
      TYPE (XPLEX)             :: LGAMA0(NCAT,NAN) ! binary activity coefficients
      TYPE (XPLEX)             :: Y     (NAN,NCAT) !
      TYPE (XPLEX)             :: BETA0 (NCAT,NAN) ! binary activity coef parameter
      TYPE (XPLEX)             :: BETA1 (NCAT,NAN) ! binary activity coef parameter
      TYPE (XPLEX)             :: CGAMA (NCAT,NAN) ! binary activity coef parameter
      TYPE (XPLEX)             :: V1    (NCAT,NAN) ! # of cations in electrolyte
                                             !   formula
      TYPE (XPLEX)             :: V2    (NCAT,NAN) ! # of anions in electrolyte
                                             !   formula
      ! absolute value of charges of cation
      TYPE (XPLEX)             :: ZP(NCAT)=(/xplex(1.0d0,0d0), 
     &                                       xplex(1.0d0,0d0) /)         

      ! absolute value of charges of anion
      TYPE (XPLEX)             :: ZM(NAN)  = (/ xplex(2.0d0,0d0),
     & xplex( 1.0d0,0d0), xplex(1.0d0,0d0) /)         

      ! Character values.
      CHARACTER(LEN=120)      :: XMSG  = ' '
      CHARACTER(LEN=16), SAVE :: PNAME = ' driver program name'

      !================================================================
      ! *** Sources for the coefficients BETA0, BETA1, CGAMA
      ! (1,1);(1,3)  - Clegg & Brimblecombe (1988)
      ! (2,3)        - Pilinis & Seinfeld (1987), cgama different
      ! (1,2)        - Clegg & Brimblecombe (1990)
      ! (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)
      !================================================================
 
      ! now set the basic constants, BETA0, BETA1, CGAMA
      DATA BETA0(1,1)/xplex(2.98d-2,0d0)/,
     &     BETA1(1,1) /xplex(0.0d0,0d0)/,
     &     CGAMA(1,1)/xplex(4.38d-2,0d0)/                                 ! 2H+SO4-

      DATA BETA0(1,2)/xplex(1.2556d-1,0d0)/,
     &     BETA1(1,2)/xplex(2.8778d-1,0d0)/,
     &     CGAMA(1,2)/xplex(-5.59d-3,0d0)/                               ! HNO3

      DATA BETA0(1,3) / xplex(2.0651d-1,0d0)/,   
     &     BETA1(1,3) / xplex(5.556d-1,0d0)/,
     &     CGAMA(1,3) /xplex(0.0d0,0d0)/                                   ! H+HSO4-

      DATA BETA0(2,1) / xplex(4.6465d-2,0d0)/,   
     &     BETA1(2,1) /xplex(-0.54196d0,0d0)/,
     &     CGAMA(2,1) /xplex(-1.2683d-3,0d0)/                              ! (NH4)2SO4

      DATA BETA0(2,2) /xplex(-7.26224d-3,0d0)/,  
     &     BETA1(2,2) /xplex(-1.168858d0,0d0)/,
     &     CGAMA(2,2) / xplex(3.51217d-5,0d0)/                             ! NH4NO3

      DATA BETA0(2,3) / xplex(4.494d-2,0d0)/,    
     &     BETA1(2,3) / xplex(2.3594d-1,0d0)/,
     &     CGAMA(2,3) /xplex(-2.962d-3,0d0)/                               ! NH4HSO4

      DATA V1(1,1),V2(1,1)/xplex(2.0d0,0d0),xplex(1.0d0,0d0) /     ! 2H+SO4-
      DATA V1(2,1),V2(2,1)/xplex(2.0d0,0d0),xplex(1.0d0,0d0) /     ! (NH4)2SO4
      DATA V1(1,2),V2(1,2)/xplex(1.0d0,0d0),xplex(1.0d0,0d0) /     ! HNO3
      DATA V1(2,2),V2(2,2)/xplex(1.0d0,0d0),xplex(1.0d0,0d0) /     ! NH4NO3
      DATA V1(1,3),V2(1,3)/xplex(1.0d0,0d0),xplex(1.0d0,0d0) /     ! H+HSO4-
      DATA V1(2,3),V2(2,3)/xplex(1.0d0,0d0),xplex(1.0d0,0d0) /     ! NH4HSO4

      !=================================================================
      ! ACTCOF begins here!
      !=================================================================

      ! Compute ionic strength
      I = 0.0d0

      DO ICAT = 1, NCAT
         I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
      ENDDO

      DO IAN = 1, NAN
         I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
      ENDDO

      I = 0.5d0 * I

      ! check for problems in the ionic strength
      IF ( I .EQ. 0.0d0 ) THEN

         DO IAN  = 1, NAN
         DO ICAT = 1, NCAT
            GAMA( ICAT, IAN ) = 0.0d0
         ENDDO
         ENDDO

         XMSG = 'Ionic strength is zero...returning zero activities'
         !CALL M3WARN ( PNAME, 0, 0, XMSG )
         RETURN

      ELSE IF ( I .LT. 0.0d0 ) THEN
         XMSG = 'Ionic strength below zero...negative concentrations'
         ! remove this on prospero (dkh, 11/18/10) 
         write(6,*)xmsg
         call flush(6)
         !CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT1 )
      ENDIF

      ! Compute some essential expressions
      SRI    = SQRT( I )
      TWOSRI = 2.0d0 * SRI
      TWOI   = 2.0d0 * I
      TEXPV  = 1.0d0 - EXP( -TWOSRI ) * ( 1.0d0 + TWOSRI - TWOI )
      R      = 1.0d0 + 0.75d0 * I
      S      = 1.0d0 + 1.5d0  * I
      ZOT1   = 0.511d0 * SRI / ( 1.0d0 + SRI )

      ! Compute binary activity coeffs
      FGAMA = -0.392d0 * ( ( SRI / ( 1.0d0 + 1.2d0 * SRI )
     &      + ( 2.0d0 / 1.2d0 ) * LOG( 1.0d0 + 1.2d0 * SRI ) ) )

      DO ICAT = 1, NCAT
      DO IAN  = 1, NAN

         BGAMA( ICAT, IAN ) = 2.0d0 * BETA0( ICAT, IAN )
     &        + ( 2.0d0 * BETA1( ICAT, IAN ) / ( 4.0d0 * I ) )
     &        * TEXPV

         ! Compute the molality of each electrolyte for given ionic strength
         M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )
     &                   *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0d0
     &                   / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )

         ! Calculate the binary activity coefficients
         LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA
     &        + M( ICAT, IAN )
     &        * ( 2.0d0 * V1( ICAT, IAN ) * V2( ICAT, IAN )
     &        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )
     &        * BGAMA( ICAT, IAN ) )
     &        + M( ICAT, IAN ) * M( ICAT, IAN )
     &        * ( 2.0d0 * ( V1( ICAT, IAN )
     &        * V2( ICAT, IAN ) )**1.5d0
     &        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )
     &        * CGAMA( ICAT, IAN ) ) ) / 2.302585093d0

      ENDDO
      ENDDO

      ! prepare variables for computing the multicomponent activity coeffs
      DO IAN = 1, NAN
      DO ICAT = 1, NCAT
         ZBAR           = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5d0
         ZBAR2          = ZBAR * ZBAR
         Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
         X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
      ENDDO
      ENDDO

      DO IAN = 1, NAN
         F1( IAN ) = 0.0d0
         DO ICAT = 1, NCAT
            F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )
     &                + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
         ENDDO
      ENDDO

      DO ICAT = 1, NCAT
         F2( ICAT ) = 0.0d0
         DO IAN = 1, NAN
            F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0(ICAT, IAN)
     &                 + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
         ENDDO
      ENDDO

      ! now calculate the multicomponent activity coefficients
      DO IAN  = 1, NAN
      DO ICAT = 1, NCAT

         TA  = -ZOT1 * ZP( ICAT ) * ZM( IAN )
         TB  = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
         TC  = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
         TRM = TA + TB * TC
         
         IF ( TRM .GT. 30.0d0 ) THEN
            GAMA( ICAT, IAN ) = 1.0d+30
         ELSE
            GAMA( ICAT, IAN ) = 10.0d0**TRM
         ENDIF

      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ACTCOF

!------------------------------------------------------------------------------
      
      SUBROUTINE INIT_RPMARES
!
!*****************************************************************************
!  Subroutine INIT_RPMARES initializes all module arrays (bmy, 12/16/02)
!
!  NOTES:
!******************************************************************************
!
      ! F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_RPMARES begins here!
      !=================================================================
      ALLOCATE( HNO3_sav( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HNO3_sav' )
      HNO3_sav = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_RPMARES

!-----------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_RPMARES
!
!*****************************************************************************
!  Subroutine CLEANUP_RPMARES deallocates all module arrays (bmy, 12/16/02)
!
!  NOTES:
!******************************************************************************
!
      IF ( ALLOCATED( HNO3_sav ) ) DEALLOCATE( HNO3_sav )

      ! Return to calling program
      END SUBROUTINE CLEANUP_RPMARES

!------------------------------------------------------------------------------

      END MODULE RPMARES_MOD
