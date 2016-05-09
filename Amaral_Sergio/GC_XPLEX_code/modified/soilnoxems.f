! $Id: soilnoxems.f,v 1.4 2010/04/28 21:00:00 daven Exp $
      SUBROUTINE SOILNOXEMS( SUNCOS )
!
!******************************************************************************
!  Subroutine SOILNOXEMS computes the emission of soil and fertilizer NOx
!  for the GEOS-CHEM model. (yhw, gmg, djj, 8/94; bdf, bmy, 10/4/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SUNCOS (TYPE (XPLEX))     : Array for COSINE( solar zenith angle ) [unitless]
!
!  Other Important variables (mostly via common blocks or F90 modules)
!  ============================================================================
!  (1 ) IJLOOP    (INTEGER) : Grid-box number                 
!  (2 ) JLOP      (INTEGER) : Index to IJLOOP; for a given I,J         
!  (3 ) NLAND     (INTEGER) : Total number of land boxes       
!  (4 ) NSRCE     (INTEGER) : Emission timestep [seconds]
!  (5 ) INDEXSOIL (INTEGER) : (I,J) of the grid                         
!  (6 ) IREG      (INTEGER) : Number of landtypes in grid square (I,J) 
!  (7 ) ILAND     (INTEGER) : ID in grid square (I,J) for IREG landtypes 
!  (8 ) IUSE      (INTEGER) : Fraction ((per mil) of box covered by land types
!  (9 ) NCONSOIL  (INTEGER) : Converts from Olson type  -> soil type
!  (10) FRAC      (TYPE (XPLEX) ) : Fraction of grid-box that rained
!  (11) RATE      (TYPE (XPLEX) ) : Rate of total rain fall mm/day             
!  (12) RPULSE    (TYPE (XPLEX) ) : Pulsing rate (computed via "pulsing.f")   
!  (13) SOILTEMP  (TYPE (XPLEX) ) : Temperature factor (external function) 
!  (14) TMMPK     (TYPE (XPLEX) ) : Local air temperature (K), w/ diurnal variation
!  (15) SOILCRF   (TYPE (XPLEX) ) : Soil canopy reduction factor [unitless]
!  (16) WINDSQR   (TYPE (XPLEX) ) : Surface winds squared [m2/s2] (from sfcwindsqr.f)
!  (17) SOILBASE  (TYPE (XPLEX) ) : Emissions         
!  (18) BXHEIGHT  (TYPE (XPLEX) ) : Grid-box height [m]
!  (19) SOILNOX   (TYPE (XPLEX) ) : Output [molec NOx/cm2/s]
!
!  References:
!  ============================================================================
!  (1 ) Yienger and Levy [1995]                                   
!  (2 ) Wang et al [1998],  Global Simulation of tropospheric           
!        O3-NOx-hydrocarbon; JGR Vol 103, pages 10713-10725       
!
!  NOTES:
!  (1 ) Be sure to force TYPE (XPLEX) with the DBLE function and the "D" 
!        exponent, wherever necessary (bmy, 10/6/99)  
!  (2 ) Made JLOP a local variable, so as not to have to reference it from  
!        "comode.h".  "comode.h" should be only for SMVGEAR. (bmy, 10/19/00) 
!  (3 ) Now save soil NOx into GEMISNOX2 array (bdf, bmy, 6/15/01)          
!  (4 ) Replaced IM with IIPAR and JM with JJPAR (bmy, 6/25/02)
!  (5 ) Now reference BXHEIGHT from "dao_mod.f".  Also updated comments and
!        made cosmetic changes. (bmy, 9/18/02)
!  (6 ) Removed NSRCE from the call to "pulsing.f".  Now add I0, J0 as local
!        variables.  Now use functions GET_XOFFSET and GET_YOFFSET from
!        "grid_mod.f". (bmy, 2/11/03)
!  (7 ) Need to pass SUNCOS to SOILCRF when computing FERTDIAG. (bmy, 10/14/03)
!  (8 ) Now references LFUTURE from "logical_mod.f" and GET_FUTURE_SCALE_NOxft
!        from "future_emissions_mod.f".  Now compute future emissions of NOx
!        from soils if necessary. (swu, bmy, 5/30/06)
!  (9 ) Bug fix: future emissions only need to be applied the fertilizer
!        term in the NOx emissions below. (swu, bmy, 10/4/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,              ONLY : BXHEIGHT
      USE DIAG_MOD,             ONLY : AD32_fe,     AD32_so
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxft
      USE GRID_MOD,             ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,          ONLY : LFUTURE
    
      ! adj_group:  add scaling factors (dkh, 11/08/09) 
      USE TIME_MOD,             ONLY : GET_DIRECTION
      USE CHECKPT_MOD,          ONLY : SOILNOX_CHK
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
 
#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! Diagnostic switches & arrays
#     include "CMN_NOX"              ! GEMISNOX2
#     include "CMN_DEP"              ! CANOPYNOX
#     include "commsoil.h"           ! Soil pulsing & wetness variables

      ! Arguments
      TYPE (XPLEX), INTENT(IN)            :: SUNCOS(MAXIJ)

      ! Local variables
      INTEGER                       :: I,    J,       M,      III,    NN  
      INTEGER                       :: K,    IREF,    JREF,   IJLOOP, L
      INTEGER                       :: KBL,  JLOP(IIPAR,JJPAR,1), I0, J0
      TYPE (XPLEX)                 :: TMMP, WINDSQR, RPULSE, ZBL,   DUM
      TYPE (XPLEX)                     :: FERTDIAG(IIPAR,JJPAR), FUT_SCL

      
      ! External functions
      TYPE (XPLEX), EXTERNAL           :: BOXVL,      FERTADD,  PULSING
      TYPE (XPLEX), EXTERNAL           :: SFCWINDSQR, SOILBASE, SOILCRF
      TYPE (XPLEX), EXTERNAL              :: SOILTEMP,   XLTMMP

      !=================================================================
      ! SOILNOXEMS begins here!
      !=================================================================

      ! Get nested-grid offsets (bmy, 2/11/03)
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      ! Initalize
      IJLOOP = 0
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SOILNOX(I,J)  = 0.D0
         FERTDIAG(I,J) = 0.D0     
         IJLOOP        = IJLOOP + 1
         JLOP(I,J,1)   = IJLOOP
      ENDDO
      ENDDO

      ! Call soiltype to determine whether soil is dry or 
      ! wet for all land grid-boxes
      CALL SOILTYPE

      ! Loop over each land grid-box
      DO M = 1, NLAND
         IREF   = INDEXSOIL(1,M)
         JREF   = INDEXSOIL(2,M)
         I      = IREF - I0
         J      = JREF - J0
         IJLOOP = JLOP(I,J,1)

         IF ( (I.GE.1) .AND. (I.LE.IIPAR)  .AND.
     &        (J.GE.1) .AND. (J.LE.JJPAR) ) THEN

            !===========================================================
            ! PULSING FACTOR "FUNCTION PULSING(I,J,M,NSRCE)"
            !
            ! ECO SYSTEM DEPENDENT
            ! TEMPERATURE FACTOR "FUNCTION SOILTEMP(I,J,M,NN)"
            ! BASE EMISSION WITH FERTERLIZATION 
            ! CANOPY REDkUCTION
            ! SOIL NOX EMISSIONS (WATCH OUT FOR TROP. EVERGREEN)
            !===========================================================
            TMMP    = XLTMMP(I,J,IJLOOP)-273.15
            WINDSQR = SFCWINDSQR(I,J)
            RPULSE  = PULSING( I, J, M ) 

            DO K = 1, IREG(IREF,JREF)
               NN = NCONSOIL(ILAND(IREF,JREF,K)+1)

               ! IPCC future emission scenario for NOx from fertilizers
               IF ( LFUTURE ) THEN
                  FUT_SCL = GET_FUTURE_SCALE_NOXft( I, J )
               ELSE
                  FUT_SCL = 1d0
               ENDIF

               ! SOILNOX contains soil NOx emissions in [molec NOx/cm2/s]
               SOILNOX(I,J) = SOILNOX(I,J) +   
     &           ( SOILTEMP(I,J,M,NN,TMMP) * SOILBASE(I,J,M,NN,RPULSE) + 
     &             FERTADD(J,M,NN)         * FUT_SCL ) 
     &          *(1.D0-SOILCRF(I,J,IREF,JREF,IJLOOP,M,NN,K,
     &             WINDSQR,SUNCOS))*XPLX(IUSE(IREF,JREF,K))/1000.D0

               ! Archive fertilizer for ND32 diagnostic (bey)
               FERTDIAG(I,J) = FERTDIAG(I,J) + 
     &             FERTADD(J,M,NN) * FUT_SCL
     &             *(1.D0-SOILCRF(I,J,IREF,JREF,IJLOOP,M,NN,K,
     &              WINDSQR,SUNCOS))*XPLX(IUSE(IREF,JREF,K))/1000.D0
            ENDDO
         ENDIF

        
         ! adj_group:  add scaling and checkpointing of soil nox
         ! (dkh, 02/06/07) 
         IF ( GET_DIRECTION() > 0 .and. LADJ ) THEN


            ! now apply scaling factor in setemis (dkh, 03/30/10) 
            !! Apply scaling factor if active 
            !IF ( IDADJ_ENOxso > 0 ) THEN 
            ! 
            !   IF ( MMSCL /= 1 ) THEN
            !      CALL ERROR_STOP('MMSCL', 'soilnoxems.f')
            !   ENDIF
            !
            !   ! Apply scaling factor 
            !   SOILNOX(I,J) = SOILNOX(I,J)
            !                * EMS_SF(I,J,1,IDADJ_ENOxso)
            ! 
            !ENDIF 

            ! Checkpoint SOILNOX emissions
            SOILNOX_CHK(I,J) = SOILNOX(I,J)

         ELSEIF( LADJ ) THEN 

            ! Overwrite with checkpointed values, thus correctly accounting 
            ! for dependence of PULSING on previous time step. 
            SOILNOX(I,J) = SOILNOX_CHK(I,J)

         ENDIF 


         ! Skip if there are no soil NOx emissions
         IF ( SOILNOX(I,J) .EQ. 0.0 ) GOTO 110

         ! Archive soil NOx and fertilizer NOx emissions [molec NOx/cm2/s]
         IF ( ND32 > 0 ) THEN
            AD32_so(I,J) = AD32_so(I,J) + SOILNOX(I,J)
            AD32_fe(I,J) = AD32_fe(I,J) + FERTDIAG(I,J)
         ENDIF

         ! Spread NOx emission into the boundary layer
         ! NOTE: BXHEIGHT is in m, so BXHEIGHT * 100 is in cm.
         ZBL = 0.D0
         KBL = 1
         DO L = 1, KBL
            ZBL = ZBL + BXHEIGHT(I,J,L)*100.D0
         ENDDO

         ! Store soil NOx in GEMISNOX2, the global NOx emissions array, which
         ! gets passed to SMVGEAR.  GEMISNOX2 has units of [molec NOx/cm3/s].
         GEMISNOX2(I,J) = GEMISNOX2(I,J) + ( SOILNOX(I,J) / ZBL )
           
 110     CONTINUE
      ENDDO ! M

      ! Return to calling program
      END SUBROUTINE SOILNOXEMS
