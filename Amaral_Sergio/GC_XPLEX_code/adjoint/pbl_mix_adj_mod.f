! $Id: pbl_mix_adj_mod.f,v 1.3 2009/11/12 00:45:48 daven Exp $
      MODULE PBL_MIX_ADJ_MOD
!
!******************************************************************************
!  Module PBL_MIX_MOD_ADJ contains adjoint routines and variables used to compute the
!  planetary boundary layer (PBL) height and to mix tracers underneath the 
!  PBL top. (ks, dkh, 07/08/09) 
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_PBL_MIX_ADJ        : Driver routine for PBL mixing
!  (2 ) TURBDAY_ADJ           : Adjoint of TURBDAY
! 
!  NOTES:
!  (1 ) Now modified for GCAP and GEOS-5 met fields (bmy, 5/24/05)
!  (2 ) Remove reference to "CMN" and XTRA2. (bmy, 8/30/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (4 ) Now recalculate IMIX, FPBL rather than checkpoint (dkh, 07/08/09) 
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "pbl_mix_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: DO_PBL_MIX_ADJ

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_PBL_MIX_ADJ( DO_TURBDAY )
!
!******************************************************************************
!  Subroutine DO_PBL_MIX is the driver routine for planetary boundary layer
!  mixing.  The PBL layer height and related quantities are always computed.
!  Complete mixing of tracers underneath the PBL top is toggled by the 
!  DO_TURBDAY switch. (bmy, 2/11/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DO_TURBDAY (LOGICAL) : Switch which turns on PBL mixing of tracers
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE LOGICAL_MOD,    ONLY : LTURB
      USE PBL_MIX_MOD,    ONLY : COMPUTE_PBL_HEIGHT
      USE TRACER_MOD,     ONLY : N_TRACERS, TCVV

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN) :: DO_TURBDAY

      !=================================================================
      ! DO_PBL_MIX_ADJ begins here!
      !=================================================================

      ! Now recompute these rather than checkpoint. (dkh, 07/08/09) 
      ! Compute PBL height and related quantities
      CALL COMPUTE_PBL_HEIGHT

      ! Do complete mixing of tracers in the PBL (if necessary)
      IF ( DO_TURBDAY ) CALL TURBDAY_ADJ( N_TRACERS, TCVV )

      ! Return to calling program
      END SUBROUTINE DO_PBL_MIX_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE TURBDAY_ADJ(NTRC, TCVV)
!
!******************************************************************************
!  Subroutine TURBDAY_ADJ executes the adjoint of the GEOS-CTM dry convection
!  / boundary layer mixing algorithm from TURBDAY. It is a combination of the
!  forward code TURBDAY with TAMC generated adjoint of the loop over N.
!  See notes in turbday.f for info about the original forward code, and
!  below for notes on modifications made for the adjoint version.
!  (dkh, 10/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) NTRC    : Number of tracers used in computation  [1 to NNPAR]
!  (2 ) TCVV    : mol. wt. air / mol. wt. tracer
!
!  Modue variable Input / Output
!  ======================================================================
!  (1 ) STT_ADJ :  Adjoint tracer array
!
!  NOTES:
!  (1 ) Rather than save / write / read the info from the forward run of TURBDAY,
!        we will just recompute most of it, hence most of the original code for
!        TURBDAY is part of ADJ_TURBDAY.  However, some alterations were made to
!        the forward code.
!        Changes to forward code:
!               -  argument list just NTRC and TCVV, which are passed
!                   the values of NADJ and ADJ_TCVV, respectfully.
!               -  add reference to CMN_ADJ for ADJ_STT
!               -  no ND15 diagnostic update, so get rid of USE DAO_MOD
!               -  get rid of XTRA2, LTURB
!               -  get rid of initial print out
!  (2 )  TAMC (and modified TAMC) code is lower case.
!        Changes to TAMC code:
!               -  The varialbes TC_IN and TC_OUT were used to construct the
!                   adjoint, but they are not needed here.  Just replace them
!                   with ADJ_STT
!               -  Replade multiple do loops with ":" operator (so no longer
!                   need integers ip1,ip2,ip3,ip4)
!               -  Force variables explicitly to TYPE (XPLEX) using .d0
!               -  Initialize and update global adjoint variables (adtc, addtc)
!                   before and after the PARALLEL DO loop
!  (3 ) Updated for v8 adjoint (dkh, 07/14/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE DAO_MOD,         ONLY : AD, PBL
      USE GRID_MOD,        ONLY : GET_AREA_M2
      USE TIME_MOD,        ONLY : GET_TS_CONV
      USE PBL_MIX_MOD,     ONLY : GET_IMIX
      USE PBL_MIX_MOD,     ONLY : GET_FPBL
      USE PRESSURE_MOD,    ONLY : GET_PEDGE
 
      ! dkh debug
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD 

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"

      ! Argument variables
      INTEGER,  INTENT(IN)     :: NTRC
      TYPE (XPLEX),   INTENT(IN)     :: TCVV(NTRC)

      ! Local variables
      INTEGER                  :: I, J, L, LTOP, N
      TYPE (XPLEX)                   :: AA,  CC, CC_AA,   BLTOP
      TYPE (XPLEX)                   :: PW,  PS, AREA_M2, DTCONV
      TYPE (XPLEX)                   :: P(0:LLPAR)
      TYPE (XPLEX)                   :: A(IIPAR,JJPAR)
      TYPE (XPLEX)                   :: DTC(IIPAR,JJPAR,LLPAR,NTRC)

      ! Adjoint variables
      TYPE (XPLEX) adcc
      TYPE (XPLEX) adcc_aa
      TYPE (XPLEX) adtc(iipar,jjpar,llpar,ntrc)
      TYPE (XPLEX) adtc_in(iipar,jjpar,llpar,ntrc)
      TYPE (XPLEX) adtc_out(iipar,jjpar,llpar,ntrc)

      !=================================================================
      ! TURBDAY_ADJ begins here!
      !=================================================================

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) '       -- TURBDAY_ADJ'

      ! Don't need DTCONV for adjoint calculation
      ! Convection timestep [s]
      !DTCONV = GET_TS_CONV() * 60d0

      ! We assume full mixing in the boundary layer, so the A
      ! coefficients are 1 everywhere, day & night (bmy, 2/11/03)
      A(:,:) = 1d0

      !----------------------------------------------
      ! SET GLOBAL ADJOINT VARIABLES
      !----------------------------------------------
      adtc(:,:,:,:)      = 0.d0
      adtc_in(:,:,:,:)   = 0.d0
      adtc_out(:,:,:,:)  = STT_ADJ(:,:,:,:)

      ! Loop over Lat/Long grid boxes (I,J)
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, N, AA, CC, CC_AA )
!!$OMP+PRIVATE( adcc, adcc_aa)
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Calculate air mass within PBL at grid box (I,J,L)
         AA = 0.d0
         DO L = 1, GET_IMIX(I,J)-1
            AA = AA + AD(I,J,L)
         ENDDO

         L  = GET_IMIX(I,J)
         AA = AA + AD(I,J,L) * GET_FPBL(I,J)

         ! Loop over tracers
         DO N = 1, NTRC

            !----------------------------------------------
            ! RESET LOCAL ADJOINT VARIABLES
            !----------------------------------------------
            adcc =             0.d0
            adcc_aa =          0.d0

            !----------------------------------------------
            ! ADJOINT ROUTINE BODY
            !----------------------------------------------
            l = GET_IMIX(i,j)
            adtc(i,j,:,n) = adtc(i,j,:,n)+adtc_out(i,j,:,n)
            adtc_out(i,j,:,n) = 0.d0
            adcc_aa = adcc_aa+adtc(i,j,l,n)*a(i,j)*GET_FPBL(i,j)
            adtc(i,j,l,n) = adtc(i,j,l,n)*(1.d0-a(i,j)*GET_FPBL(i,j))
            do l = 1, GET_IMIX(i,j)-1
              adcc_aa = adcc_aa+adtc(i,j,l,n)*a(i,j)
              adtc(i,j,l,n) = adtc(i,j,l,n)*(1.d0-a(i,j))
            end do
            adcc = adcc+adcc_aa/aa
            adcc_aa = 0.d0
            l = GET_IMIX(i,j)
            adtc(i,j,l,n) = adtc(i,j,l,n)+adcc*ad(i,j,l)*GET_FPBL(i,j)
            do l = 1, GET_IMIX(i,j)-1
              adtc(i,j,l,n) = adtc(i,j,l,n)+adcc*ad(i,j,l)
            end do
            adtc_in(i,j,:,n) = adtc_in(i,j,:,n)+adtc(i,j,:,n)
            adtc(i,j,:,n) = 0.d0
          !WRITE(6,*)'AIIIIICIIIII'
         ENDDO    !N
      ENDDO       !I
      ENDDO       !J
!!$OMP END PARALLEL DO
      WRITE(6,*)'AIIIIICIIIII'
      ! Update global adjoint variables
      STT_ADJ(:,:,:,:) = adtc_in(:,:,:,:)
      adtc_in(:,:,:,:) = 0d0
      WRITE(6,*)'AIIIIICIIIII'
      !  Return to calling program
      END SUBROUTINE TURBDAY_ADJ

!------------------------------------------------------------------------------

      ! End of module
      END MODULE PBL_MIX_ADJ_MOD
