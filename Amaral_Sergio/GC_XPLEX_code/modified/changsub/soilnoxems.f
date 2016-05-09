! $Id: soilnoxems.f,v 1.1.1.1 2009/06/09 21:51:54 daven Exp $
      SUBROUTINE SOILNOXEMS( SUNCOS )
!
!******************************************************************************
!  Subroutine SOILNOXEMS computes the emission of soil and fertilizer NOx
!  for the GEOS-CHEM model. (yhw, gmg, djj, 8/94; bdf, bmy, 10/4/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SUNCOS (REAL*8)     : Array for COSINE( solar zenith angle ) [unitless]
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
!  (10) FRAC      (REAL*8 ) : Fraction of grid-box that rained
!  (11) RATE      (REAL*8 ) : Rate of total rain fall mm/day             
!  (12) RPULSE    (REAL*8 ) : Pulsing rate (computed via "pulsing.f")   
!  (13) SOILTEMP  (REAL*8 ) : Temperature factor (external function) 
!  (14) TMMPK     (REAL*8 ) : Local air temperature (K), w/ diurnal variation
!  (15) SOILCRF   (REAL*8 ) : Soil canopy reduction factor [unitless]
!  (16) WINDSQR   (REAL*8 ) : Surface winds squared [m2/s2] (from sfcwindsqr.f)
!  (17) SOILBASE  (REAL*8 ) : Emissions         
!  (18) BXHEIGHT  (REAL*8 ) : Grid-box height [m]
!  (19) SOILNOX   (REAL*8 ) : Output [molec NOx/cm2/s]
!
!  References:
!  ============================================================================
!  (1 ) Yienger and Levy [1995]                                   
!  (2 ) Wang et al [1998],  Global Simulation of tropospheric           
!        O3-NOx-hydrocarbon; JGR Vol 103, pages 10713-10725       
!
!  NOTES:
!  (1 ) Be sure to force double precision with the DBLE function and the "D" 
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

      IMPLICIT NONE
 
#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! Diagnostic switches & arrays
#     include "CMN_NOX"              ! GEMISNOX2
#     include "CMN_DEP"              ! CANOPYNOX
#     include "commsoil.h"           ! Soil pulsing & wetness variables

      ! Arguments
      REAL*8, INTENT(IN)            :: SUNCOS(MAXIJ)

      ! Local variables
      INTEGER                       :: I,    J,       M,      III,    NN  
      INTEGER                       :: K,    IREF,    JREF,   IJLOOP, L
      INTEGER                       :: KBL,  JLOP(IIPAR,JJPAR,1), I0, J0
      REAL*8                        :: TMMP, WINDSQR, RPULSE, ZBL,   DUM
      REAL*8                        :: FERTDIAG(IIPAR,JJPAR), FUT_SCL

      
      ! External functions
      REAL*8, EXTERNAL              :: BOXVL,      FERTADD,  PULSING
      REAL*8, EXTERNAL              :: SFCWINDSQR, SOILBASE, SOILCRF
      REAL*8, EXTERNAL              :: SOILTEMP,   XLTMMP

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
     &             WINDSQR,SUNCOS))*DBLE(IUSE(IREF,JREF,K))/1000.D0

               ! Archive fertilizer for ND32 diagnostic (bey)
               FERTDIAG(I,J) = FERTDIAG(I,J) + 
     &             FERTADD(J,M,NN) * FUT_SCL
     &             *(1.D0-SOILCRF(I,J,IREF,JREF,IJLOOP,M,NN,K,
     &              WINDSQR,SUNCOS))*DBLE(IUSE(IREF,JREF,K))/1000.D0
            ENDDO
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
