! $Id: CO_strat_pl_adj.f,v 1.1 2010/05/07 20:39:47 daven Exp $
      SUBROUTINE CO_STRAT_PL_ADJ( COPROD, COLOSS )
!
!******************************************************************************
!  Subroutine CO_STRAT_PL_ADJ computes adjoint of net production of CO above 
!  annual mean tropopause using archived rates for P(CO) and L(CO).
!  (dkh, 05/02/10) 
!
!  Based on forward model  (bnd, qli, bmy, 12/9/99, 10/25/05)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) COPROD : (TYPE (XPLEX)) Zonally averaged P(CO) in [v/v/s]
!  (2 ) COLOSS : (TYPE (XPLEX)) Zonally averaged L(CO) in [1/s]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules 
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE DAO_MOD,        ONLY : AD
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : XNUMOL,   XNUMOLAIR
      USE TRACERID_MOD,   ONLY : IDTCO,           IDTCH2O
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP, GET_MIN_TPAUSE_LEVEL

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
     
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      TYPE (XPLEX),  INTENT(IN) :: COPROD(JJPAR,LLPAR)
      TYPE (XPLEX),  INTENT(IN) :: COLOSS(JJPAR,LLPAR)

      ! Local variables
      INTEGER             :: I, J, L, M, N, LMIN

      TYPE (XPLEX)              :: BAIRDENS, DT, GCO, STTTOGCO

      ! External functions
      TYPE (XPLEX), EXTERNAL    :: BOXVL

      !=================================================================
      ! CO_STRAT_PL_ADJ begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DT = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Loop over all stratospheric grid boxes ( L >= LPAUSE(I,J) ). 
      !
      ! Compute the net CO from the P(CO) and L(CO) rates that are 
      ! stored in the COPROD and COLOSS arrays.
      !
      ! Unit conversion to/from [kg/box] and [molec/cm3] is required.
      ! The conversion factor is STTTOGCO, which is given below.
      !
      !   kg CO       box     |   mole CO   | 6.022e23 molec CO       
      !  ------- * -----------+-------------+-------------------  
      !    box      BOXVL cm3 | 28e-3 kg CO |     mole CO             
      !
      !  =  molec CO
      !     --------
      !       cm3
      !=================================================================

      ! Get the minimum extent of the tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

      DO L = LMIN, LLPAR
      DO J = 1,    JJPAR
      DO I = 1,    IIPAR

         ! Skip tropospheric grid boxes
         IF ( ITS_IN_THE_TROP(I,J,L) ) CYCLE

         ! fwd code:
         !STTTOGCO = 6.022d23 / ( 28d-3 * BOXVL(I,J,L) )
         !GCO = STT(I,J,L,IDTCO) * STTTOGCO
         !GCO = GCO * ( 1d0 - COLOSS(J,L) * DT ) +
         !            ( COPROD(J,L) * DT * BAIRDENS )
         !STT(I,J,L,IDTCO) = GCO / STTTOGCO
         ! adj code (production does not affect adjoint):
         STT_ADJ(I,J,L,IDTCO) = STT_ADJ(I,J,L,IDTCO)
     &                        * ( 1d0 - COLOSS(J,L) * DT )

         ! production does not affect adjoint 
         !STT(I,J,L,IDTCH2O) = STT(I,J,L,IDTCH2O) + 
         !                     XNUMOL(IDTCO) / XNUMOL(IDTCH2O) * 
         !                     COPROD(J,L)   * BAIRDENS        / 
         !                     STTTOGCO


      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CO_STRAT_PL_ADJ
