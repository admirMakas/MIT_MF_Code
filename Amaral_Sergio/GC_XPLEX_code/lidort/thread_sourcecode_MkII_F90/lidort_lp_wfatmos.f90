!$Id: lidort_lp_wfatmos.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.5 F90                               #
! #  Release Date :   June 2010                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #                                                         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level routines--------------                   #
! #            UPUSER_PROFILEWF                            #
! #            DNUSER_PROFILEWF                            #
! #            MIFLUX_PROFILEWF                            #
! #            LIDORT_LP_CONVERGE                          #
! #                                                        #
! #     Output at quadrature angles ---------              #
! #            QUADPROFILEWF_LEVEL_UP                      #
! #            QUADPROFILEWF_LEVEL_DN                      #
! #            QUADPROFILEWF_OFFGRID_UP                    #
! #            QUADPROFILEWF_OFFGRID_DN                    #
! #                                                        #
! #     Post-processing at user angles --------            #
! #            GET_LP_TOASOURCE                            #
! #            GET_LP_BOASOURCE                            #
! #            LP_WHOLELAYER_STERM_UP                      #
! #            LP_WHOLELAYER_STERM_DN                      #
! #            LP_PARTLAYER_STERM_UP                       #
! #            LP_PARTLAYER_STERM_DN                       #
! #            LP_QUAD_GFUNCMULT                           #
! #                                                        #
! ##########################################################

SUBROUTINE UPUSER_PROFILEWF                                            &
           ( DO_USER_STREAMS,  DO_PLANE_PARALLEL,                      & ! Input
             DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER,   & ! Input
             NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
             UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                  & ! Input
             FOURIER_COMPONENT, IBEAM, VARIATION_INDEX, K_PARAMETERS,  & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,     & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,               & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM,   & ! Input
             AGM, BGP, XPOS, LCON, LCON_XVEC, MCON, MCON_XVEC,         & ! Input
             U_XPOS, U_XNEG, U_WPOS, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
             UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                   & ! Input
             PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,             & ! Input
             UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                   & ! Input
             L_T_DELT_USERM, L_T_UTUP_USERM, HELP_AQ, HELP_BQ,         & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,       & ! Input
             L_T_DELT_EIGEN, L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,         & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                & ! Input
             L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC, & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                        & ! Input
             L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,             & ! Input
             LP_BOA_MSSOURCE, LP_BOA_DBSOURCE,                         & ! Input
             FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,           & ! Output
             PROFILEWF_F, QUADPROFILEWF )                                ! Output

!  Upwelling post-processed Profile Jacobian Fourier component

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL, intent(in)  ::   DO_USER_STREAMS
      LOGICAL, intent(in)  ::   DO_MSMODE_LIDORT
      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  local control flags

      LOGICAL, intent(in)  ::   DO_INCLUDE_MVOUTPUT

!  Control integers

      INTEGER, intent(in)  ::   NSTREAMS
      INTEGER, intent(in)  ::   N_USER_STREAMS
      INTEGER, intent(in)  ::   NLAYERS
      INTEGER, intent(in)  ::   N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL, intent(in)  ::   PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Linearization control

      INTEGER, intent(in)  ::   K_PARAMETERS
      INTEGER, intent(in)  ::   VARIATION_INDEX

!  Input Fourier number and beam index
!  surface factor (2 for m = 0, 1 otherwise). Not required.
!  Flux multiplier = F/4.pi

      INTEGER, intent(in)  ::   FOURIER_COMPONENT, IBEAM
      REAL(kind=8), intent(in)  ::  FLUX_MULTIPLIER

!  Regular Inputs
!  --------------

!  Local flags for the solution saving option

      LOGICAL, intent(in)  ::   DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(kind=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(kind=8), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(kind=8), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Saved Quantitites from the Green's function calculation

      REAL(kind=8), intent(in)  ::  AGM(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=8), intent(in)  ::  XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS_2,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS_2,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=8), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(kind=8), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(kind=8), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(kind=8), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(kind=8), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(kind=8), intent(in)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(kind=8), intent(in)  :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(kind=8), intent(in)  :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(in)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(kind=8), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  :: LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances, solar beam

      REAL(kind=8), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(kind=8), intent(in)  :: LP_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(kind=8), intent(in)  :: LP_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(kind=8), intent(in)  :: LP_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized BOA source terms

      REAL(kind=8), intent(in)  :: LP_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: LP_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Outputs
!  -------

!  Profile weighting functions at quadrature angles

      REAL(kind=8), intent(out) :: QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS,  MAX_USER_LEVELS, &
                                                   MAXSTREAMS,   MAXBEAMS,   MAX_DIRECTIONS )

!  Profile weighting functions at user angles

      REAL(kind=8), intent(out) :: PROFILEWF_F ( MAX_ATMOSWFS,     MAXLAYERS, MAX_USER_LEVELS, &
                                                 MAX_USER_STREAMS, MAXBEAMS,  MAX_DIRECTIONS)

!  Linearized Green's function multipliers for off-grid optical depths

      LOGICAL, intent(out)      :: FLAGS_LP_GMULT(MAX_PARTLAYERS)

!  Linearized Green functions multipliers for off-grid optical depths
!   Will only be calculated as output, if the flag has been set

      REAL(kind=8), intent(out) :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       ::  N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER       ::  UTA, UM, Q, NC, UT, IB, K

      REAL(kind=8)  ::  L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  L_FINAL_SOURCE

!  index

      K  = VARIATION_INDEX
      IB = IBEAM

!  Zero all Fourier component output here (safety)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, K_PARAMETERS
            DO UM = 1, N_USER_STREAMS
              PROFILEWF_F(Q,K,UTA,UM,IB,UPIDX) = ZERO
             ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

!  start the recursion

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_CUMUL_SOURCE(UM,Q) = LP_BOA_MSSOURCE(UM,Q) + LP_BOA_DBSOURCE(UM,Q)
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            CALL LP_WHOLELAYER_STERM_UP                           &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT,                 & ! Input
             NSTREAMS, N_USER_STREAMS, IB, N, K, K_PARAMETERS,    & ! Input
             DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),            & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,       & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,        & ! Input
             HMULT_1,  HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,     & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ, & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,  & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,           & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                   & ! Input
             L_LAYER_SOURCE )                                       ! Output

            IF ( N.EQ.K ) THEN
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)          + &
                         T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,Q)   + &
                       L_T_DELT_USERM(N,UM,Q) *   CUMSOURCE_UP(UM,NC-1)
                ENDDO
              ENDDO
            ELSE
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  + &
                        T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            FLAGS_LP_GMULT(UT) = .TRUE.
            CALL QUADPROFILEWF_OFFGRID_UP                                                  &
          ( DO_PLANE_PARALLEL, NSTREAMS, IB, UTA, UT, N, K, K_PARAMETERS, FLUX_MULTIPLIER, & ! Input
            LAYER_PIS_CUTOFF, T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,  & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,                  & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS, L_ATERM_SAVE, L_BTERM_SAVE,          & ! Input
            UT_GMULT_UP, UT_GMULT_DN, XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,              & ! Input
            L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,                                          & ! Input
            QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                  ! Output
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

            CALL LP_PARTLAYER_STERM_UP                                             &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS, N_USER_STREAMS,        & ! Input
             IB, UT, N, K, K_PARAMETERS, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N), & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,                        & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                         & ! Input
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,      & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,                  & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,                   & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                            & ! Input
             L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,                         & ! Input
             L_LAYER_SOURCE )                                                        ! Output

            IF ( N.EQ.K ) THEN
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)             + &
                        T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q) + &
                      L_T_UTUP_USERM(UT,UM,Q)*   CUMSOURCE_UP(UM,NC)
                  PROFILEWF_F(Q,K,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ELSE
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q) + T_UTUP_USERM(UT,UM) * L_CUMUL_SOURCE(UM,Q)
                  PROFILEWF_F(Q,K,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADPROFILEWF_LEVEL_UP                            &
           ( NLAYERS, NSTREAMS, IB, UTA, NLEVEL, K, K_PARAMETERS,  & ! Input
             FLUX_MULTIPLIER, L_XPOS, L_XNEG, L_WLOWER, L_WUPPER,  & ! Input
             LCON, LCON_XVEC, NCON_XVEC,   T_DELT_EIGEN,           & ! Input
             MCON, MCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN,           & ! Input
             QUADPROFILEWF )                                         ! Output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                PROFILEWF_F(Q,K,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_PROFILEWF

!

SUBROUTINE DNUSER_PROFILEWF                                           &
           ( DO_USER_STREAMS,  DO_PLANE_PARALLEL,                     & ! Input
             DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                   & ! Input
             NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                 & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                 & ! Input
             FOURIER_COMPONENT,                                       & ! Input
             IBEAM, VARIATION_INDEX, K_PARAMETERS, FLUX_MULTIPLIER,   & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,    & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,              & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM,  & ! Input
             AGM, BGP, XPOS, LCON, LCON_XVEC, MCON, MCON_XVEC,        & ! Input
             U_XPOS, U_XNEG, U_WNEG, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
             UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                  & ! Input
             PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,            & ! Input
             UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                  & ! Input
             L_T_DELT_USERM, L_T_UTDN_USERM, HELP_AQ, HELP_BQ,        & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,      & ! Input
             L_T_DELT_EIGEN, L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,        & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,      & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,               & ! Input
             L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,          & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                       & ! Input
             L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,            & ! Input
             LP_TOA_SOURCE,                                           & ! Input
             FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,          & ! Input
             PROFILEWF_F, QUADPROFILEWF )                               ! Output

!  Upwelling post-processed Profile Jacobian Fourier component

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL, intent(in)  ::   DO_USER_STREAMS
      LOGICAL, intent(in)  ::   DO_MSMODE_LIDORT
      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  local control flags

      LOGICAL, intent(in)  ::   DO_INCLUDE_MVOUTPUT

!  Control integers

      INTEGER, intent(in)  ::   NSTREAMS
      INTEGER, intent(in)  ::   N_USER_STREAMS
      INTEGER, intent(in)  ::   N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL, intent(in)  ::   PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Input Fourier number and beam index
!  Flux multiplier = F/4.pi

      INTEGER, intent(in)  ::   FOURIER_COMPONENT, IBEAM
      REAL(kind=8), intent(in)  ::  FLUX_MULTIPLIER

!  Linearization control

      INTEGER, intent(in)  ::   K_PARAMETERS
      INTEGER, intent(in)  ::   VARIATION_INDEX

!  Regular Inputs
!  --------------

!  Local flags for the solution saving option

      LOGICAL, intent(in)  ::   DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)       ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(kind=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initia     I        FLUX_MULTIPLIER,l setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(kind=8), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(kind=8), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(kind=8), intent(in)  ::  AGM(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=8), intent(in)  ::  XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS_2,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS_2,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=8), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(kind=8), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(kind=8), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(kind=8), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(kind=8), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(kind=8), intent(in)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(kind=8), intent(in)  :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(kind=8), intent(in)  :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(kind=8), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  :: LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances, solar beam

      REAL(kind=8), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(kind=8), intent(in)  :: LP_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(kind=8), intent(in)  :: LP_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(kind=8), intent(in)  :: LP_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized TOA source terms

      REAL(kind=8), intent(in)  :: LP_TOA_SOURCE  (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Outputs
!  -------

!  Profile weighting functions at quadrature angles

      REAL(kind=8), intent(out) :: QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS,  MAX_USER_LEVELS, &
                                                   MAXSTREAMS,   MAXBEAMS,   MAX_DIRECTIONS )

!  Profile weighting functions at user angles

      REAL(kind=8), intent(out) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                 MAX_USER_STREAMS, MAXBEAMS,  MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths

      LOGICAL     , intent(inout) :: FLAGS_LP_GMULT(MAX_PARTLAYERS)
      REAL(kind=8), intent(out)   :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out)   :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL       ::  LOCAL_LP_GMULT
      INTEGER       ::  N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER       ::  UTA, UM, Q, NC, UT, IB, K
      REAL(kind=8)  ::  L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  L_FINAL_SOURCE

!  Initialise

      K  = VARIATION_INDEX
      IB = IBEAM

!  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, K_PARAMETERS
            DO UM = 1, 1 
              PROFILEWF_F(Q,K,UTA,UM,IB,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      CALL GET_LP_TOASOURCE ( N_USER_STREAMS, LP_TOA_SOURCE, K_PARAMETERS )

! Initialize

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_CUMUL_SOURCE(UM,Q) = LP_TOA_SOURCE(UM,Q)
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N

            CALL LP_WHOLELAYER_STERM_DN                                         &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS, N_USER_STREAMS,     & ! Input
             IB, N, K, K_PARAMETERS, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),  & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,                     & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                      & ! Input
             HMULT_1,  HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,                   & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,               & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,                & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,                         & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                                 & ! Input
             L_LAYER_SOURCE )                                                     ! Output

            IF ( N.EQ.K ) THEN
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  + &
                       T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) + &
                      L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
                ENDDO
              ENDDO
            ELSE
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
                ENDDO
              ENDDO
            ENDIF

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_LP_GMULT = FLAGS_LP_GMULT(UT)
            CALL QUADPROFILEWF_OFFGRID_DN                                         &
          ( DO_PLANE_PARALLEL, NSTREAMS, IB, UTA, UT, N, K, K_PARAMETERS,         & ! Input
            FLUX_MULTIPLIER, LOCAL_LP_GMULT, LAYER_PIS_CUTOFF,                    & ! Input
            T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,           & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR,  LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,        & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS, L_ATERM_SAVE, L_BTERM_SAVE, & ! Input
            UT_GMULT_UP, UT_GMULT_DN, XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,     & ! Input
            L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,                                 & ! Input
            QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                         ! Output
            FLAGS_LP_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LP_PARTLAYER_STERM_DN                                             &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS, N_USER_STREAMS,        & ! Input
             IB, UT, N, K, K_PARAMETERS, DO_LAYER_SCATTERING(FOURIER_COMPONENT,N), & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,                        & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                         & ! Input
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,      & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,                  & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,                   & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,                            & ! Input
             L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,                         & ! Input
             L_LAYER_SOURCE )                                                        ! Output

            IF ( N.EQ.K ) THEN
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)            + &
                      T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,Q) + &
                    L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,NC)
                PROFILEWF_F(Q,K,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ELSE
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q) + T_UTDN_USERM(UT,UM) * L_CUMUL_SOURCE(UM,Q)
                PROFILEWF_F(Q,K,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADPROFILEWF_LEVEL_DN                  & 
           ( NSTREAMS, IB, UTA, NLEVEL, K, K_PARAMETERS, & ! Input
             FLUX_MULTIPLIER, L_XPOS, L_XNEG, L_WLOWER,  & ! Input
             LCON, MCON, LCON_XVEC,   T_DELT_EIGEN,      & ! Input
             NCON_XVEC,  PCON_XVEC, L_T_DELT_EIGEN,      & ! Input
             QUADPROFILEWF )                               ! Output    
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                PROFILEWF_F(Q,K,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_PROFILEWF

!

SUBROUTINE MIFLUX_PROFILEWF                                        &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,    & ! Input
             THREAD, IB, K, K_PARAMETERS,                          & ! Input
             NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,                 & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,              & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,              & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                           & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,          & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADPROFILEWF,            & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,   & ! Input
             MINT_PROFILEWF, FLUX_PROFILEWF )                        ! Output

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  input arguments
!  ---------------

!  Flags

      LOGICAL, intent(in)  ::   DO_UPWELLING
      LOGICAL, intent(in)  ::   DO_DNWELLING
      LOGICAL, intent(in)  ::   DO_INCLUDE_DIRECTBEAM

!  Thread

      INTEGER, intent(in)  ::   THREAD

!  Index

      INTEGER, intent(in)  ::   IB

!  linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  Control integers

      INTEGER, intent(in)  ::   NSTREAMS
      INTEGER, intent(in)  ::   N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL, intent(in)  ::   PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux factor

      REAL(kind=8), intent(in)  ::  FLUX_FACTOR

!  Quadrature

      REAL(kind=8), intent(in)  ::  QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(kind=8), intent(in)  ::  QUAD_STRMWTS ( MAXSTREAMS )

!  initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)       ::  LAYER_PIS_CUTOFF(MAXBEAMS)

!  local solar zenith angle cosine

      REAL(kind=8), intent(in)  :: LOCAL_CSZA ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(kind=8), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances

      REAL(kind=8), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Quadrature-defined weighting functions

      REAL(kind=8), intent(in)  :: QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                   MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  Output arguments
!  ----------------

!  Mean intensity (actini! flux)

      REAL(kind=8), intent(out) :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXBEAMS,     MAX_DIRECTIONS, MAXTHREADS )

!  Flux

      REAL(kind=8), intent(out) :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXBEAMS,     MAX_DIRECTIONS, MAXTHREADS )

!  local variables
!  ----------------

      INTEGER       ::  I, UTA, UT, Q, N
      REAL(kind=8)  ::  SM, SF, FMU0
      REAL(kind=8)  ::  L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

!  mean intensity and flux
!  -----------------------

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
           SM = SM + QUAD_WEIGHTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX)
           SF = SF + QUAD_STRMWTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX)
          ENDDO
          MINT_PROFILEWF(Q,K,UTA,IB,UPIDX,THREAD) = SM * HALF
          FLUX_PROFILEWF(Q,K,UTA,IB,UPIDX,THREAD) = SF * PI2
         ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN

!  Diffuse term contribution

        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
           SM = SM + QUAD_WEIGHTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX)
           SF = SF + QUAD_STRMWTS(I)*QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX)
          ENDDO
          MINT_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) = SM * HALF
          FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) = SF * PI2
         ENDDO
        ENDDO

!  nothing to do if no solar sources

        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  For the downward direction, add the direct beam contributions

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    L_INITIAL_TRANS is a logarithmi! derivative

            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
              DO Q = 1, K_PARAMETERS
                L_TRANS = LP_T_UTDN_MUBAR(UT,K,IB,Q) + LP_INITIAL_TRANS(N,K,IB,Q) * T_UTDN_MUBAR(UT,IB)
                L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                L_DIRECT_FLUX  = FMU0 * L_TRANS
                MINT_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) = MINT_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) + L_DIRECT_MEANI
                FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) = FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) + L_DIRECT_FLUX
              ENDDO
            ENDIF

!  For the on-grid balues

          ELSE
            N = UTAU_LEVEL_MASK_DN(UTA)
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              IF ( N.GT.0 ) THEN
                FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
                DO Q = 1, K_PARAMETERS
                  L_TRANS = LP_T_DELT_MUBAR(N,K,IB,Q) + LP_INITIAL_TRANS(N,K,IB,Q) * T_DELT_MUBAR(N,IB)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                  L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                  L_DIRECT_FLUX  = FMU0 * L_TRANS
                  MINT_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) = MINT_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) + L_DIRECT_MEANI
                  FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) = FLUX_PROFILEWF(Q,K,UTA,IB,DNIDX,THREAD) + L_DIRECT_FLUX
                ENDDO
              ENDIF
            ENDIF
          ENDIF

!  Close loops

        ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_PROFILEWF

!

SUBROUTINE QUADPROFILEWF_LEVEL_UP                               &
         ( NLAYERS, NSTREAMS, IB, UTA, NL, K, K_PARAMETERS,     & ! Input
           FLUX_MULTIPLIER, L_XPOS, L_XNEG, L_WLOWER, L_WUPPER, & ! Input
           LCON, LCON_XVEC, NCON_XVEC,   T_DELT_EIGEN,          & ! Input
           MCON, MCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN,          & ! Input
           QUADPROFILEWF )                                        ! Output

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments
!  --------------------------

!  Control

      INTEGER, intent(in)  ::   NLAYERS, NSTREAMS

!  Indices

      INTEGER, intent(in)  ::   IB, UTA, NL

!  Flux

      REAL(kind=8), intent(in)  ::  FLUX_MULTIPLIER

!  linearization control

      INTEGER, intent(in)       ::  K, K_PARAMETERS

!  Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  ::  T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  ::  L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(in)  ::  L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  ::  L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!  Quadrature-defined weighting functions

      REAL(kind=8), intent(out) ::  QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER       ::  N, I, I1, AA, Q
      REAL(kind=8)  ::  SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5, FM

!  homogeneous and particular solution contributions SHOM and SPAR

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1
      FM = FLUX_MULTIPLIER

!  For the lowest level

      IF ( NL .EQ. NLAYERS  ) THEN

!  If this is also the layer that is varying, extra contributions

        IF ( K .EQ. NL ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
                HOM2 = LCON_XVEC(I1,AA,NL) * L_T_DELT_EIGEN(AA,NL,Q)
                HOM3 = LCON(AA,NL)*T_DELT_EIGEN(AA,NL)*L_XPOS(I1,AA,NL,Q)
                HOM4 = PCON_XVEC(I1,AA,NL,Q)
                HOM5 = MCON(AA,NL) * L_XNEG(I1,AA,NL,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WLOWER(I1,NL,Q)
              QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  non-varying lowest layer

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
                HOM2 = PCON_XVEC(I1,AA,NL,Q)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = L_WLOWER(I1,NL,Q)
              QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

       ENDIF

!  For other levels in the atmosphere
!  ----------------------------------

      ELSE

!  If this is also the layer that is varying, extra contributions

        IF ( K .EQ. N ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,N,Q) 
                HOM2 = LCON(AA,N) * L_XPOS(I1,AA,N,Q)
                HOM3 = MCON_XVEC(I1,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM4 = MCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XNEG(I1,AA,N,Q)
                HOM5 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WUPPER(I1,N,Q)
              QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  non-varying layer lower than the varying layer

        ELSE IF ( K.LT.N ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,N,Q)
                HOM2 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = L_WUPPER(I1,N,Q)
              QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  non-varying layer higher than the varying layer

         ELSE IF ( K.GT.N ) THEN

           DO I = 1, NSTREAMS
             I1 = I + NSTREAMS
             DO Q = 1, K_PARAMETERS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                 HOM1 = NCON_XVEC(I1,AA,N,Q)
                 HOM2 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                 SHOM = SHOM + HOM1 + HOM2
               ENDDO
               QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FM * SHOM
             ENDDO
           ENDDO

        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADPROFILEWF_LEVEL_UP

!

SUBROUTINE QUADPROFILEWF_LEVEL_DN                     &
         ( NSTREAMS, IB, UTA, NL, K, K_PARAMETERS,    & ! input
           FLUX_MULTIPLIER, L_XPOS, L_XNEG, L_WLOWER, & ! input
           LCON, MCON, LCON_XVEC,   T_DELT_EIGEN,     & ! input
           NCON_XVEC,  PCON_XVEC, L_T_DELT_EIGEN,     & ! input
           QUADPROFILEWF )                              ! Output

!  Downwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments
!  --------------------------

!  Control

      INTEGER, intent(in)  ::   NSTREAMS

!  Indices

      INTEGER, intent(in)  ::   IB, UTA, NL

!  Flux

      REAL(kind=8), intent(in)  ::  FLUX_MULTIPLIER

!  linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  ::  T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  ::  L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(in)  ::  L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  ::  L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!  Quadrature-defined weighting functions

      REAL(kind=8), intent(out) ::  QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER       ::  N, I, AA, Q
      REAL(kind=8)  ::  SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5, FM

!  homogeneous and particular solution contributions SHOM and SPAR

      N = NL
      FM = FLUX_MULTIPLIER

!  Downwelling weighting function at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
             QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = ZERO
          ENDDO
        ENDDO

!  For other levels in the atmosphere
!  ----------------------------------

      ELSE

!  If this is also the layer that is varying, extra contributions

        IF ( K .EQ. N ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                HOM4 = PCON_XVEC(I,AA,N,Q)
                HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WLOWER(I,N,Q)
              QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  non-varying layer lower than the varying layer

        ELSE IF ( K.LT.N ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q) 
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              SPAR = L_WLOWER(I,N,Q)
              QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  non-varying layer higher than the varying layer

        ELSE IF ( K.GT.N ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = PCON_XVEC(I,AA,N,Q) 
                SHOM = SHOM + HOM1 + HOM2
              ENDDO
              QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = FM * SHOM
            ENDDO
          ENDDO

        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADPROFILEWF_LEVEL_DN

!

SUBROUTINE QUADPROFILEWF_OFFGRID_UP                                                       &
          ( DO_PLANE_PARALLEL, NSTREAMS, IB, UTA, UT, N, K, K_PARAMETERS, FLUX_MULTIPLIER, & ! Input
            LAYER_PIS_CUTOFF, T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,  & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,                  & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS, L_ATERM_SAVE, L_BTERM_SAVE,          & ! Input
            UT_GMULT_UP, UT_GMULT_DN, XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,              & ! Input
            L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,                                          & ! Input
            QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                                  ! Output

!  Linearization of Quadrature Jacobians (off-grid only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS

!  Indices

      INTEGER, intent(in)  ::   IB, UTA, UT, N, K, K_PARAMETERS

!  Flux

      REAL(kind=8), intent(in)  ::  FLUX_MULTIPLIER

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(kind=8), intent(in)  ::  T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(kind=8), intent(in)  ::  T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized transmittances

      REAL(kind=8), intent(in)  ::  L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(kind=8), intent(in)  ::  XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(in)  ::  UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  ::  L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!  Quadrature-defined weighting functions

      REAL(kind=8), intent(out) ::  QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  Linearized Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(out) ::  LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) ::  LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       ::  I, I1, AA, Q
      REAL(kind=8)  ::  SPAR, SHOM, HOM1, HOM2, PAR1, PAR2
      REAL(kind=8)  ::  H1, H2, H3, H4, H5, H6, FMULT

!  short hand

      FMULT = FLUX_MULTIPLIER

!  Homogeneous solutions
!  ---------------------

!  solution for N being the varying layer K

      IF ( N.EQ.K ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_XVEC(I1,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
              H2 = NCON_XVEC(I1,AA,N,Q)
              H3 = LCON(AA,N)   * L_XPOS(I1,AA,N,Q)
              H4 = MCON_XVEC(I1,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
              H5 = PCON_XVEC(I1,AA,N,Q)
              H6 = MCON(AA,N)   * L_XNEG(I1,AA,N,Q)
              HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
              HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FMULT * SHOM
          ENDDO
        ENDDO

!  Solution for N not equal to K

      ELSE IF ( N .NE. K ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = T_UTDN_EIGEN(AA,UT) * NCON_XVEC(I1,AA,N,Q)
              HOM2 = T_UTUP_EIGEN(AA,UT) * PCON_XVEC(I1,AA,N,Q)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = FMULT * SHOM
          ENDDO
        ENDDO

      ENDIF

!  get the linearized Green's function multipliers
!  -----------------------------------------------

      CALL LP_QUAD_GFUNCMULT                              &
          ( DO_PLANE_PARALLEL, NSTREAMS,                  & ! Input
            IB, UT, N, K, K_PARAMETERS, LAYER_PIS_CUTOFF, & ! Input
              T_UTUP_EIGEN,     T_UTDN_EIGEN,             & ! Input
            L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,             & ! Input
               T_DELT_MUBAR,      T_UTDN_MUBAR,           & ! Input
            LP_T_DELT_MUBAR,   LP_T_UTDN_MUBAR,           & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS,     & ! Input
            L_ATERM_SAVE, L_BTERM_SAVE,                   & ! Input
               UT_GMULT_UP,    UT_GMULT_DN,               & ! Input
            LP_UT_GMULT_UP, LP_UT_GMULT_DN )                ! Output

!  Sum up the Green's function contributions
!  -----------------------------------------

!  solution for N being the varying layer K

      IF ( N.EQ.K ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = ZERO
            DO AA = 1, NSTREAMS
              PAR1 = L_XPOS(I,AA,N,Q)  *    UT_GMULT_UP(AA,UT) +    &
                       XPOS(I,AA,N)    * LP_UT_GMULT_UP(AA,UT,K,Q)
              PAR2 = L_XPOS(I1,AA,N,Q) *    UT_GMULT_DN(AA,UT) +    &
                       XPOS(I1,AA,N)   * LP_UT_GMULT_DN(AA,UT,K,Q)
              SPAR = SPAR + PAR1 + PAR2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = &
              QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) + FMULT * SPAR
          ENDDO
        ENDDO

!  Solution for N > K (N below the layer that varies)

      ELSE IF ( N.GT.K ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = ZERO
            DO AA = 1, NSTREAMS
              PAR1 = XPOS(I,AA,N)  * LP_UT_GMULT_UP(AA,UT,K,Q)
              PAR2 = XPOS(I1,AA,N) * LP_UT_GMULT_DN(AA,UT,K,Q)
              SPAR = SPAR + PAR1 + PAR2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) = &
              QUADPROFILEWF(Q,K,UTA,I,IB,UPIDX) + FMULT * SPAR
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADPROFILEWF_OFFGRID_UP

!

SUBROUTINE QUADPROFILEWF_OFFGRID_DN                                               &
          ( DO_PLANE_PARALLEL, NSTREAMS, IB, UTA, UT, N, K, K_PARAMETERS,         & ! Input
            FLUX_MULTIPLIER, LOCAL_LP_GMULT, LAYER_PIS_CUTOFF,                    & ! Input
            T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,           & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR,  LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,        & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS, L_ATERM_SAVE, L_BTERM_SAVE, & ! Input
            UT_GMULT_UP, UT_GMULT_DN, XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,     & ! Input
            L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,                                 & ! Input
            QUADPROFILEWF, LP_UT_GMULT_UP, LP_UT_GMULT_DN )                         ! Output

!  Linearization of Quadrature Jacobians (off-grid only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS

!  Indices

      INTEGER, intent(in)  ::   IB, UTA, UT, N, K, K_PARAMETERS

!  Flux

      REAL(kind=8), intent(in)  ::  FLUX_MULTIPLIER

!  Local flag for getting the multipliers

      LOGICAL, intent(in)  ::   LOCAL_LP_GMULT

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(kind=8), intent(in)  ::  T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(kind=8), intent(in)  ::  T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized transmittances

      REAL(kind=8), intent(in)  ::  L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(kind=8), intent(in)  ::  XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(in)  ::  UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  ::  L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!  Quadrature-defined weighting functions

      REAL(kind=8), intent(out) ::  QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS )

!  Linearized Green functions multipliers for off-grid optical depths
!   Will only be calculated as output, if the flag has been set

      REAL(kind=8), intent(inout) ::  LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(inout) ::  LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       ::  I, I1, AA, Q
      REAL(kind=8)  ::  SPAR, SHOM, HOM1, HOM2, PAR1, PAR2
      REAL(kind=8)  ::  H1, H2, H3, H4, H5, H6, FMULT

!  Short hand

      FMULT = FLUX_MULTIPLIER

!  Homogeneous
!  -----------

!  solution for N being the varying layer K

      IF ( N.EQ.K ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_XVEC(I,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
              H2 = NCON_XVEC(I,AA,N,Q)
              H3 = LCON(AA,N)   * L_XPOS(I,AA,N,Q)
              H4 = MCON_XVEC(I,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
              H5 = PCON_XVEC(I,AA,N,Q)
              H6 = MCON(AA,N)   * L_XNEG(I,AA,N,Q)
              HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
              HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = FMULT * SHOM
          ENDDO
        ENDDO

!  Solution for N not equal to K

      ELSE IF ( N.NE.K ) THEN
      
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = T_UTDN_EIGEN(AA,UT) * NCON_XVEC(I,AA,N,Q)
              HOM2 = T_UTUP_EIGEN(AA,UT) * PCON_XVEC(I,AA,N,Q)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = FMULT * SHOM
          ENDDO
        ENDDO

      ENDIF

!  get the linearized Green's function multipliers
!  -----------------------------------------------

!    Will only be done if these are flagged.

      IF ( LOCAL_LP_GMULT ) THEN
        CALL LP_QUAD_GFUNCMULT                            &
          ( DO_PLANE_PARALLEL, NSTREAMS,                  & ! Input
            IB, UT, N, K, K_PARAMETERS, LAYER_PIS_CUTOFF, & ! Input
              T_UTUP_EIGEN,     T_UTDN_EIGEN,             & ! Input
            L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,             & ! Input
               T_DELT_MUBAR,      T_UTDN_MUBAR,           & ! Input
            LP_T_DELT_MUBAR,   LP_T_UTDN_MUBAR,           & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS,     & ! Input
            L_ATERM_SAVE, L_BTERM_SAVE,                   & ! Input
               UT_GMULT_UP,    UT_GMULT_DN,               & ! Input
            LP_UT_GMULT_UP, LP_UT_GMULT_DN )                ! Output
      ENDIF

!  Sum up the Green's function contributions
!  -----------------------------------------

!  solution for N being the varying layer K, or K = 0

      IF ( N.EQ.K ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = ZERO
            DO AA = 1, NSTREAMS
              PAR1 = L_XPOS(I1,AA,N,Q)  *    UT_GMULT_UP(AA,UT) +  &
                       XPOS(I1,AA,N)    * LP_UT_GMULT_UP(AA,UT,K,Q)
              PAR2 = L_XPOS(I,AA,N,Q)   *    UT_GMULT_DN(AA,UT) +  &
                       XPOS(I,AA,N)     * LP_UT_GMULT_DN(AA,UT,K,Q)
              SPAR = SPAR + PAR1 + PAR2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = &
              QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) + FMULT * SPAR
          ENDDO
        ENDDO

!  Solution for N > K (N below the layer that varies)

      ELSE IF ( N.GT.K ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = ZERO
            DO AA = 1, NSTREAMS
              PAR1 = XPOS(I1,AA,N)  * LP_UT_GMULT_UP(AA,UT,K,Q)
              PAR2 = XPOS(I,AA,N)   * LP_UT_GMULT_DN(AA,UT,K,Q)
              SPAR = SPAR + PAR1 + PAR2
            ENDDO
            QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) = &
              QUADPROFILEWF(Q,K,UTA,I,IB,DNIDX) + FMULT * SPAR
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADPROFILEWF_OFFGRID_DN

!

SUBROUTINE GET_LP_TOASOURCE ( N_USER_STREAMS, L_TOA_SOURCE, K_PARAMETERS )

!  Linearized Bottom of the atmosphere source term

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments
!  --------------------

!  control integer

      INTEGER, intent(in)  ::   N_USER_STREAMS

!  Linearization

      INTEGER, intent(in)  ::   K_PARAMETERS

!  output

      REAL(kind=8), intent(out) ::  L_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       ::   UM, Q

!  initialise TOA source function
!  ------------------------------

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_TOA_SOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  Finish

END SUBROUTINE GET_LP_TOASOURCE

!

SUBROUTINE GET_LP_BOASOURCE                                              &
    ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_BRDF_SURFACE,        & ! Input
      DO_USER_STREAMS, NLAYERS, NSTREAMS, N_USER_STREAMS,                & ! Input
      FOURIER_COMPONENT, IBEAM, K, K_PARAMETERS,SURFACE_FACTOR,          & ! Input
      ALBEDO, USER_BRDF_F, QUAD_STRMWTS, USER_DIRECT_BEAM, DELTAU_SLANT, & ! Input
      L_DELTAU_VERT, L_XPOS, L_XNEG, LCON, LCON_XVEC, T_DELT_EIGEN,      & ! Input
      MCON, NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_WLOWER,              & ! Input
      L_BOA_MSSOURCE, L_BOA_DBSOURCE )                                     ! Output

!  Linearized Bottom of the atmosphere source term

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, intent(in)  ::   DO_INCLUDE_SURFACE
      LOGICAL, intent(in)  ::   DO_BRDF_SURFACE
      LOGICAL, intent(in)  ::   DO_INCLUDE_DIRECTBEAM
      LOGICAL, intent(in)  ::   DO_USER_STREAMS

!  control integers

      INTEGER, intent(in)  ::   NLAYERS, NSTREAMS, N_USER_STREAMS

!  surface multiplier, albedo and Fourier/beam indices

      REAL(kind=8), intent(in)  ::  SURFACE_FACTOR, ALBEDO
      INTEGER, intent(in)  ::   FOURIER_COMPONENT
      INTEGER, intent(in)  ::   IBEAM

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected user streams

      REAL(KIND=8), intent(in)  :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Quadrature

      REAL(kind=8), intent(in)  ::  QUAD_STRMWTS ( MAXSTREAMS )

!  linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  ::  T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Direct beam solutions

      REAL(kind=8), intent(in)  ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Linearized deltaQUAD_STRMWTSu

      REAL(kind=8), intent(in)  ::  L_DELTAU_VERT(MAX_ATMOSWFS, MAXLAYERS )

!  Slant optical thickness values

      REAL(kind=8), intent(in)  ::  DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  ::  L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  ::  NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(in)  ::  L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  ::  L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(kind=8), intent(out) ::  L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) ::  L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       ::   M, N, J, I, UM, AA, Q, IB
      REAL(kind=8)  ::  L_DOWN(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  REFLEC, L_BEAM, FAC, KMULT
      REAL(kind=8)  ::  SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

!  Starting section
!  ----------------

!  Fourier number, layer number

      M  = FOURIER_COMPONENT
      N  = NLAYERS
      IB = IBEAM

!  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_MSSOURCE(UM,Q) = ZERO
            L_BOA_DBSOURCE(UM,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Exit if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Linearization of downwelling quadrature field at surface
!    Set reflectance integrand  a(j).x(j).L_DOWN(-j)

!  Two cases:
!  If  K = N, this is also the layer that is varying --> Extras!
!  If  N > K with variations in layer K above N

      IF ( K.EQ.N ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
              HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
              HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
              HOM4 = PCON_XVEC(I,AA,N,Q)
              HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
              SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
            ENDDO
            L_DOWN(I,Q) = ( SHOM + L_WLOWER(I,N,Q) ) * QUAD_STRMWTS(I)
          ENDDO
        ENDDO
      ELSE IF (K.LT.N.AND.K.NE.0) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
              HOM2 = PCON_XVEC(I,AA,N,Q) 
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            L_DOWN(I,Q) = ( SHOM + L_WLOWER(I,N,Q) ) * QUAD_STRMWTS(I)
           ENDDO
        ENDDO
      ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  .. integrate reflectance, same for all user-streams in Lambertian case

      IF ( .not. DO_BRDF_SURFACE ) THEN

        KMULT = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO Q = 1, K_PARAMETERS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + L_DOWN(J,Q)
            ENDDO
            REFLEC = KMULT * REFLEC
            IF ( DO_USER_STREAMS ) THEN
              DO UM = 1, N_USER_STREAMS
                L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      
!  .. integrate reflectance, BRDF case

      ELSE IF ( DO_BRDF_SURFACE ) THEN

        DO Q = 1, K_PARAMETERS
          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,Q) * USER_BRDF_F(M,UM,J)
              ENDDO
              REFLEC = SURFACE_FACTOR * REFLEC
              L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
            ENDDO
          ENDIF
        ENDDO

      ENDIF
        
!  Add direct beam if flagged.
!   This is different from the column WF case

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          FAC = - USER_DIRECT_BEAM(UM,IB) * DELTAU_SLANT(N,K,IB) 
          DO Q = 1, K_PARAMETERS
            L_BEAM = L_DELTAU_VERT(Q,K) * FAC
            L_BOA_DBSOURCE(UM,Q) = L_BEAM
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_LP_BOASOURCE

!

SUBROUTINE LP_WHOLELAYER_STERM_UP                                     &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS,           & ! Input
             N_USER_STREAMS, IB, N, K, K_PARAMETERS, SOURCETERM_FLAG, & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,           & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,            & ! Input
             HMULT_1,  HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,         & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,     & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,      & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,               & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                       & ! Input
             L_LAYERSOURCE )                                            ! Output

!  Linearization of Post-processed multiplier (Whole layers only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  MSMODE flag

      LOGICAL, intent(in)  ::   DO_MSMODE_LIDORT

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB)

      INTEGER, intent(in)  ::   N, IB

!  Source term flag

      LOGICAL, intent(in)  ::   SOURCETERM_FLAG

!  Linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)       ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(kind=8), intent(in)  ::  AGM(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8), intent(in)  ::  HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  forcing term multipliers, whole layer

      REAL(kind=8), intent(in)  ::  EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(kind=8), intent(in)  ::  PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  ::  HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(kind=8), intent(in)  ::  LP_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(kind=8), intent(in)  ::  NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8), intent(in)  ::  L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(kind=8), intent(in)  ::  LP_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(kind=8), intent(out) ::  L_LAYERSOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(kind=8)  ::  LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(kind=8)  ::  MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(kind=8)  ::  NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER       ::  AA, UM, Q
      REAL(kind=8)  ::  SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1
      REAL(kind=8)  ::  L_SD, L_SU, L_WDEL, CONST, WDEL
      REAL(kind=8)  ::  L_UP(MAXSTREAMS), L_DN(MAXSTREAMS)

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  Exit if no source term

      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  These quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Homogeneous solutions
!  =====================

!  Special case when N = K, additional layer contributions
!  Other cases when N not equal to K (only variation of Integ-Cons)

      IF ( N.EQ.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
              H3 = LCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
              H4 = MCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
              H6 = MCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
              SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ELSE IF ( N.NE.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * HMULT_2(AA,UM,N)
              H5 = PCON_UXVEC(UM,AA,Q) * HMULT_1(AA,UM,N)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ENDIF

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Nothing further to do if Varying layer is below N

      IF ( N .LT. K ) RETURN

!  Linearization of Green's function particular solution
!  =====================================================

!  For the case N > K, see later section of code

!  This section for N = K
!  ----------------------

      IF ( N.EQ.K ) THEN

!  Some beam particulars

       WDEL  = T_DELT_MUBAR(N,IB)
       CONST = INITIAL_TRANS(N,IB)

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)
        DO UM = 1, N_USER_STREAMS

!  Linearized Green function multipliers

          DO AA = 1, NSTREAMS
            L_SD = CONST * L_HMULT_2(AA,UM,N,Q)
            L_SD = L_SD -  LP_EMULT_UP(UM,N,K,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LP_GAMMA_M(AA,N,K,Q)
            T1 = T1 * PMULT_UD(AA,UM,N)
            L_DN(AA) = AGM(AA,N) * L_SD + T1 
            L_SU = - CONST * ( L_HMULT_1(AA,UM,N,Q) *   WDEL + &
                                 HMULT_1(AA,UM,N)   * L_WDEL )
            L_SU =  L_SU + LP_EMULT_UP(UM,N,K,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LP_GAMMA_P(AA,N,K,Q)
            T1 = T1 * PMULT_UU(AA,UM,N)
            L_UP(AA) = BGP(AA,N) * L_SU + T1 
          ENDDO

!  Sum the particular solution contributions, add result to total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 =   U_XPOS(UM,AA,N)   * L_DN(AA)
            H2 = L_U_XPOS(UM,AA,N,Q) * PMULT_UD(AA,UM,N)
            H3 =   U_XNEG(UM,AA,N)   * L_UP(AA)
            H4 = L_U_XNEG(UM,AA,N,Q) * PMULT_UU(AA,UM,N)
            SPAR = SPAR + H1 + H2 + H3 + H4
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LP_U_WPOS(UM,N,Q) *    EMULT_UP(UM,N,IB) +   &
                      U_WPOS(UM,N)   * LP_EMULT_UP(UM,N,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
        ENDDO
       ENDIF

!  Finish

       RETURN

!  This section for N > K
!  ----------------------

      ELSE

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        DO UM = 1, N_USER_STREAMS

!  Distinguish between plane-parallel and pseudo-spherical cases

          IF ( DO_PLANE_PARALLEL ) THEN
            T1 = LP_INITIAL_TRANS(N,K,IB,Q)
            DO AA = 1, NSTREAMS
              L_DN(AA) = T1 * PMULT_UD(AA,UM,N)
              L_UP(AA) = T1 * PMULT_UU(AA,UM,N)
            ENDDO
          ELSE
            DO AA = 1, NSTREAMS
              L_SD = HELP_AQ(N,K,IB,Q) * HMULT_2(AA,UM,N)
              L_SD = L_SD -  LP_EMULT_UP(UM,N,K,IB,Q)
              T1 = LP_GAMMA_M(AA,N,K,Q) * PMULT_UD(AA,UM,N)
              L_DN(AA) = AGM(AA,N) * L_SD + T1 
              L_SU = HELP_BQ(N,K,IB,Q)  * HMULT_1(AA,UM,N)
              L_SU =  L_SU + LP_EMULT_UP(UM,N,K,IB,Q)
              T1 = LP_GAMMA_P(AA,N,K,Q) * PMULT_UU(AA,UM,N)
              L_UP(AA) = BGP(AA,N) * L_SU + T1 
            ENDDO
          ENDIF

!  Sum the particular solution contributions, add result to total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 = U_XPOS(UM,AA,N) * L_DN(AA)
            H3 = U_XNEG(UM,AA,N) * L_UP(AA)
            SPAR = SPAR + H1 + H3
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = U_WPOS(UM,N) * LP_EMULT_UP(UM,N,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
         ENDDO
        ENDDO
       ENDIF

!  End clauses for N and K

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_WHOLELAYER_STERM_UP

!

SUBROUTINE LP_WHOLELAYER_STERM_DN                                           &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS, N_USER_STREAMS, & ! Input
             IB, N, K, K_PARAMETERS, SOURCETERM_FLAG,                       & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,                 & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                  & ! Input
             HMULT_1,  HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,               & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,           & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,            & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,                     & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                             & ! Input
             L_LAYERSOURCE )                                                  ! Output

!  Linearization of Post-processed multiplier (Whole layers only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  MSMODE flag

      LOGICAL, intent(in)  ::   DO_MSMODE_LIDORT

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB)

      INTEGER, intent(in)  ::   N, IB

!  Source term flag

      LOGICAL, intent(in)  ::   SOURCETERM_FLAG

!  Linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(kind=8), intent(in)  ::  AGM(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8), intent(in)  ::  HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  forcing term multipliers, whole layer

      REAL(kind=8), intent(in)  ::  EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(kind=8), intent(in)  ::  PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  ::  HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(kind=8), intent(in)  ::  LP_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(kind=8), intent(in)  ::  NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8), intent(in)  ::  L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(kind=8), intent(in)  ::  LP_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(kind=8), intent(out) ::  L_LAYERSOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(kind=8)  ::  LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(kind=8)  ::  MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(kind=8)  ::  NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER       ::  AA, UM, Q
      REAL(kind=8)  ::  SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1
      REAL(kind=8)  ::  L_SD, L_SU, L_WDEL, CONST, WDEL
      REAL(kind=8)  ::  L_UP(MAXSTREAMS), L_DN(MAXSTREAMS)

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  Exit if no source term

      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  Combined quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Homogeneous solutions
!  =====================

!  Special case when N = K, additional layer contributions
!  Other cases when N not equal to K (only variation of Integ-Cons)

      IF ( N.EQ.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
              H3 = LCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
              H4 = MCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
              H6 = MCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
            SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ELSE
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * HMULT_1(AA,UM,N)
              H5 = PCON_UXVEC(UM,AA,Q) * HMULT_2(AA,UM,N)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ENDIF

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Nothing further to do if Varying layer is below N

      IF ( N .LT. K ) RETURN

!  Linearization of Green's function particular solution
!  =====================================================

!  For the case N > K, see later section of code

!  This section for N = K
!  ----------------------

      IF ( N.EQ.K ) THEN

!  Some beam particulars

       WDEL  = T_DELT_MUBAR(N,IB)
       CONST = INITIAL_TRANS(N,IB)

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS

        L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)

!  start local user angle loop and initialize

        DO UM = 1, N_USER_STREAMS

!  Linearized Green;s function multipliers

          DO AA = 1, NSTREAMS
            L_SD = CONST * L_HMULT_1(AA,UM,N,Q)
            L_SD = L_SD - LP_EMULT_DN(UM,N,K,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LP_GAMMA_M(AA,N,K,Q)
            T1 = T1 * PMULT_DD(AA,UM,N)
            L_DN(AA) = AGM(AA,N) * L_SD + T1
            L_SU = - CONST * ( L_HMULT_2(AA,UM,N,Q) *   WDEL + &
                                 HMULT_2(AA,UM,N)   * L_WDEL )
            L_SU =  L_SU + LP_EMULT_DN(UM,N,K,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LP_GAMMA_P(AA,N,K,Q)
            T1 = T1 * PMULT_DU(AA,UM,N)
            L_UP(AA) = BGP(AA,N) * L_SU + T1
          ENDDO

!  Sum the particular solution contributions, add result to total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 =   U_XNEG(UM,AA,N)   * L_DN(AA)
            H2 = L_U_XNEG(UM,AA,N,Q) * PMULT_DD(AA,UM,N)
            H3 =   U_XPOS(UM,AA,N)   * L_UP(AA)
            H4 = L_U_XPOS(UM,AA,N,Q) * PMULT_DU(AA,UM,N)
            SPAR = SPAR + H1 + H2 + H3 + H4
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LP_U_WNEG(UM,N,Q) *    EMULT_DN(UM,N,IB) +    &
                      U_WNEG(UM,N)   * LP_EMULT_DN(UM,N,K,IB,Q)
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
         ENDDO
       ENDIF

!  finish this section

       RETURN

!  This section for N > K
!  ----------------------

      ELSE

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        DO UM = 1, N_USER_STREAMS

!  Distinguish between plane-parallel and pseudo-spherical cases

          IF ( DO_PLANE_PARALLEL ) THEN
            T1 = LP_INITIAL_TRANS(N,K,IB,Q)
            DO AA = 1, NSTREAMS
              L_DN(AA) = T1 * PMULT_DD(AA,UM,N)
              L_UP(AA) = T1 * PMULT_DU(AA,UM,N)
            ENDDO
          ELSE
            DO AA = 1, NSTREAMS
              L_SD = HELP_AQ(N,K,IB,Q) * HMULT_1(AA,UM,N)
              L_SD = L_SD -  LP_EMULT_DN(UM,N,K,IB,Q)
              T1 = LP_GAMMA_M(AA,N,K,Q) * PMULT_DD(AA,UM,N)
              L_DN(AA) = AGM(AA,N) * L_SD + T1 
              L_SU = HELP_BQ(N,K,IB,Q)  * HMULT_2(AA,UM,N)
              L_SU =  L_SU + LP_EMULT_DN(UM,N,K,IB,Q)
              T1 = LP_GAMMA_P(AA,N,K,Q) * PMULT_DU(AA,UM,N)
              L_UP(AA) = BGP(AA,N) * L_SU + T1 
            ENDDO
          ENDIF

!  Sum the particular solution contributions, add result to total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 = U_XNEG(UM,AA,N) * L_DN(AA)
            H3 = U_XPOS(UM,AA,N) * L_UP(AA)
            SPAR = SPAR + H1 + H3
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = U_WNEG(UM,N) * LP_EMULT_DN(UM,N,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
         ENDDO
        ENDDO
       ENDIF

!  End clause for N and K

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_WHOLELAYER_STERM_DN

!

SUBROUTINE LP_PARTLAYER_STERM_UP                                                   &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS, N_USER_STREAMS,        & ! Input
             IB, UT, N, K, K_PARAMETERS, SOURCETERM_FLAG,                          & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,                        & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                         & ! Input
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,      & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,                  & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,                   & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WPOS, NCON, PCON,                            & ! Input
             L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,                         & ! Input
             L_LAYERSOURCE )                                                         ! Output

!  Linearization of Post-processed multiplier (partial layers only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  MSMODE flag

      LOGICAL, intent(in)  ::   DO_MSMODE_LIDORT

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB), offgrid index UT

      INTEGER, intent(in)  ::   N, IB, UT

!  Source term flag

      LOGICAL, intent(in)  ::   SOURCETERM_FLAG

!  Linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(kind=8), intent(in)  ::  AGM(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8), intent(in)  ::  UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  forcing term multipliers, partial layer

      REAL(kind=8), intent(in)  ::  UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (partial layer)

      REAL(kind=8), intent(in)  ::  UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  ::  HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(kind=8), intent(in)  ::  LP_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(kind=8), intent(in)  ::  NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8), intent(in)  ::  L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(kind=8), intent(in)  ::  LP_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  output linearized layer source term

      REAL(kind=8), intent(out) ::  L_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(kind=8)  ::  LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(kind=8)  ::  MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(kind=8)  ::  NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      INTEGER       ::  AA, UM, Q
      REAL(kind=8)  ::  SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1
      REAL(kind=8)  ::  L_SD, L_SU, L_WDEL, CONST, WDEL
      REAL(kind=8)  ::  L_UP(MAXSTREAMS), L_DN(MAXSTREAMS)

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  Exit if no source term

      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  These quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Special case when N = K, additional contributions
!  Other cases when N not equal to K (only variation of Integ-Cons)

      IF ( N.EQ.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_UD(AA,UM,UT,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_UD(AA,UM,UT)
              H3 = LCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
              H3 = H3 * UT_HMULT_UD(AA,UM,UT)
              H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_UU(AA,UM,UT,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_UU(AA,UM,UT)
              H6 = MCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
              H6 = H6 * UT_HMULT_UU(AA,UM,UT)
              SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ELSE IF ( N.NE.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * UT_HMULT_UD(AA,UM,UT)
              H5 = PCON_UXVEC(UM,AA,Q) * UT_HMULT_UU(AA,UM,UT)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ENDIF

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Nothing further to do if Varying layer is below N

      IF ( N .LT. K ) RETURN

!  Linearization of Green's function particular solution
!  =====================================================

!  For the case N > K, see later section of code

!  This section for N = K
!  ----------------------

      IF ( N.EQ.K ) THEN

!  Some beam particulars

       WDEL  = T_DELT_MUBAR(N,IB)
       CONST = INITIAL_TRANS(N,IB)

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)
        DO UM = 1, N_USER_STREAMS

!  Linearized Green's function multipliers

          DO AA = 1, NSTREAMS
            L_SD = CONST * L_UT_HMULT_UD(AA,UM,UT,Q)
            L_SD = L_SD -  LP_UT_EMULT_UP(UM,UT,K,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LP_GAMMA_M(AA,N,K,Q)
            T1 = T1 * UT_PMULT_UD(AA,UM,UT)
            L_DN(AA) = AGM(AA,N) * L_SD + T1 
            L_SU = - CONST * ( L_UT_HMULT_UU(AA,UM,UT,Q) *   WDEL + &
                                 UT_HMULT_UU(AA,UM,UT)   * L_WDEL )
            L_SU =  L_SU + LP_UT_EMULT_UP(UM,UT,K,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LP_GAMMA_P(AA,N,K,Q)
            T1 = T1 * UT_PMULT_UU(AA,UM,UT)
            L_UP(AA) = BGP(AA,N) * L_SU + T1 
          ENDDO

!  Sum the particular solution contributions, Add result to the total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 =   U_XPOS(UM,AA,N)   * L_DN(AA)
            H2 = L_U_XPOS(UM,AA,N,Q) * UT_PMULT_UD(AA,UM,UT)
            H3 =   U_XNEG(UM,AA,N)   * L_UP(AA)
            H4 = L_U_XNEG(UM,AA,N,Q) * UT_PMULT_UU(AA,UM,UT)
            SPAR = SPAR + H1 + H2 + H3 + H4
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LP_U_WPOS(UM,N,Q) *    UT_EMULT_UP(UM,UT,IB) +    &
                      U_WPOS(UM,N)   * LP_UT_EMULT_UP(UM,UT,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
        ENDDO
       ENDIF

!  Finish this section and return

       RETURN

!  This section for N > K
!  ----------------------

      ELSE

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        DO UM = 1, N_USER_STREAMS

!  Linearized Green's function multipliers
!    Distinguish between plane-parallel and pseudo-spherical cases

          IF ( .NOT.DO_PLANE_PARALLEL ) THEN
            DO AA = 1, NSTREAMS
              L_SD = HELP_AQ(N,K,IB,Q) * UT_HMULT_UD(AA,UM,UT) - LP_UT_EMULT_UP(UM,UT,K,IB,Q)
              T1 = LP_GAMMA_M(AA,N,K,Q) * UT_PMULT_UD(AA,UM,UT)
              L_DN(AA) = AGM(AA,N) * L_SD + T1
              L_SU = HELP_BQ(N,K,IB,Q) * UT_HMULT_UU(AA,UM,UT) + LP_UT_EMULT_UP(UM,UT,K,IB,Q)
              T1 = LP_GAMMA_P(AA,N,K,Q) * UT_PMULT_UU(AA,UM,UT)
              L_UP(AA) = BGP(AA,N) * L_SU + T1
            ENDDO
          ELSE
            DO AA = 1, NSTREAMS
              T1 = LP_INITIAL_TRANS(N,K,IB,Q)
              L_DN(AA) = T1 * UT_PMULT_UD(AA,UM,UT)
              L_UP(AA) = T1 * UT_PMULT_UU(AA,UM,UT)
            ENDDO
          ENDIF

!  Sum the particular solution contributions, Add result to the total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 = U_XPOS(UM,AA,N) * L_DN(AA)
            H3 = U_XNEG(UM,AA,N) * L_UP(AA)
            SPAR = SPAR + H1 + H3
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = U_WPOS(UM,N)  * LP_UT_EMULT_UP(UM,UT,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
         ENDDO
        ENDDO
       ENDIF

!  End clause for N and K

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_PARTLAYER_STERM_UP

!

SUBROUTINE LP_PARTLAYER_STERM_DN                                                   &
           ( DO_PLANE_PARALLEL, DO_MSMODE_LIDORT, NSTREAMS, N_USER_STREAMS,        & ! Input
             IB, UT, N, K, K_PARAMETERS, SOURCETERM_FLAG,                          & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,                        & ! Input
             AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                         & ! Input
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,      & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, HELP_AQ, HELP_BQ,                  & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,                   & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WNEG, NCON, PCON,                            & ! Input
             L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,                         & ! Input
             L_LAYERSOURCE )                                                         ! Output

!  Linearization of Post-processed multiplier (partial layers only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  MSMODE flag

      LOGICAL, intent(in)  ::   DO_MSMODE_LIDORT

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB), offgrid index UT

      INTEGER, intent(in)  ::   N, IB, UT

!  Source term flag

      LOGICAL, intent(in)  ::   SOURCETERM_FLAG

!  Linearization control

      INTEGER, intent(in)  ::   K, K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(kind=8), intent(in)  ::  AGM(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(kind=8), intent(in)  ::  U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(kind=8), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8), intent(in)  ::   UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::   UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  forcing term multipliers, partial layer

      REAL(kind=8), intent(in)  ::  UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (partial layer)

      REAL(kind=8), intent(in)  ::  UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  ::  HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8), intent(in)  ::  L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(kind=8), intent(in)  ::  LP_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(kind=8), intent(in)  ::  NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8), intent(in)  ::  L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(kind=8), intent(in)  ::  LP_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  output linearized layer source term

      REAL(kind=8), intent(out) ::  L_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(kind=8)  ::  LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(kind=8)  ::  MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(kind=8)  ::  NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  ::  PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      INTEGER       ::  AA, UM, Q
      REAL(kind=8)  ::  SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1
      REAL(kind=8)  ::  L_SD, L_SU, L_WDEL, CONST, WDEL
      REAL(kind=8)  ::  L_UP(MAXSTREAMS), L_DN(MAXSTREAMS)

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  Exit if no source term

      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  Combined quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Special case when N = K, additional contributions
!  Other cases when N not equal to K (only variation of Integ-Cons)

      IF ( N.EQ.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_DD(AA,UM,UT,Q)
              H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_DD(AA,UM,UT)
              H3 = LCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
              H3 = H3 * UT_HMULT_DD(AA,UM,UT)
              H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_DU(AA,UM,UT,Q)
              H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_DU(AA,UM,UT)
              H6 = MCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
              H6 = H6 * UT_HMULT_DU(AA,UM,UT)
              SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ELSE IF ( N.NE.K ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              H2 = NCON_UXVEC(UM,AA,Q) * UT_HMULT_DD(AA,UM,UT)
              H5 = PCON_UXVEC(UM,AA,Q) * UT_HMULT_DU(AA,UM,UT)
              SHOM = SHOM + H2 + H5
            ENDDO
            L_LAYERSOURCE(UM,Q) = SHOM
          ENDDO
        ENDDO
      ENDIF

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Nothing further to do if Varying layer is below N

      IF ( N .LT. K ) RETURN

!  Linearization of Green's function particular solution
!  =====================================================

!  For the case N > K, see later section of code

!  This section for N = K
!  ----------------------

      IF ( N.EQ.K ) THEN

!  Some beam particulars

       WDEL  = T_DELT_MUBAR(N,IB)
       CONST = INITIAL_TRANS(N,IB)

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)
        DO UM = 1, N_USER_STREAMS

!  Linearized Green's function multipliers

          DO AA = 1, NSTREAMS
            L_SD = CONST * L_UT_HMULT_DD(AA,UM,UT,Q)
            L_SD = L_SD -  LP_UT_EMULT_DN(UM,UT,K,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LP_GAMMA_M(AA,N,K,Q)
            T1 = T1 * UT_PMULT_DD(AA,UM,UT)
            L_DN(AA) = AGM(AA,N) * L_SD + T1 
            L_SU = - CONST * ( L_UT_HMULT_DU(AA,UM,UT,Q) *   WDEL + &
                                 UT_HMULT_DU(AA,UM,UT)   * L_WDEL )
            L_SU =  L_SU + LP_UT_EMULT_DN(UM,UT,K,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LP_GAMMA_P(AA,N,K,Q)
            T1 = T1 * UT_PMULT_DU(AA,UM,UT)
            L_UP(AA) = BGP(AA,N) * L_SU + T1 
          ENDDO

!  Sum the particular solution contributions, Add result to the total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 =   U_XNEG(UM,AA,N)   * L_DN(AA)
            H2 = L_U_XNEG(UM,AA,N,Q) * UT_PMULT_DD(AA,UM,UT)
            H3 =   U_XPOS(UM,AA,N)   * L_UP(AA)
            H4 = L_U_XPOS(UM,AA,N,Q) * UT_PMULT_DU(AA,UM,UT)
            SPAR = SPAR + H1 + H2 + H3 + H4
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LP_U_WNEG(UM,N,Q) *    UT_EMULT_DN(UM,UT,IB) +   &
                      U_WNEG(UM,N)   * LP_UT_EMULT_DN(UM,UT,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
        ENDDO
       ENDIF

!  Finish this section and return

       RETURN

!  This section for N > K
!  ----------------------

      ELSE

!  Start parameter and local user angle loops

       DO Q = 1, K_PARAMETERS
        DO UM = 1, N_USER_STREAMS

!  Linearized Green's function multipliers
!    Distinguish between plane-parallel and pseudo-spherical cases

          IF ( .NOT.DO_PLANE_PARALLEL ) THEN
            DO AA = 1, NSTREAMS
              L_SD = HELP_AQ(N,K,IB,Q) * UT_HMULT_DD(AA,UM,UT) - LP_UT_EMULT_DN(UM,UT,K,IB,Q)
              T1 = LP_GAMMA_M(AA,N,K,Q) * UT_PMULT_DD(AA,UM,UT)
              L_DN(AA) = AGM(AA,N) * L_SD + T1
              L_SU = HELP_BQ(N,K,IB,Q) * UT_HMULT_DU(AA,UM,UT) + LP_UT_EMULT_DN(UM,UT,K,IB,Q)
              T1 = LP_GAMMA_P(AA,N,K,Q) * UT_PMULT_DU(AA,UM,UT)
              L_UP(AA) = BGP(AA,N) * L_SU + T1
            ENDDO
          ELSE
            DO AA = 1, NSTREAMS
              T1 = LP_INITIAL_TRANS(N,K,IB,Q)
              L_DN(AA) = T1 * UT_PMULT_DD(AA,UM,UT)
              L_UP(AA) = T1 * UT_PMULT_DU(AA,UM,UT)
            ENDDO
          ENDIF

!  Sum the particular solution contributions, Add result to the total

          SPAR = ZERO
          DO AA = 1, NSTREAMS
            H1 = U_XNEG(UM,AA,N) * L_DN(AA)
            H3 = U_XPOS(UM,AA,N) * L_UP(AA)
            SPAR = SPAR + H1 + H3
          ENDDO
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
       ENDDO

!  Add single scatter term if flagged

       IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = U_WNEG(UM,N)  * LP_UT_EMULT_DN(UM,UT,K,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
         ENDDO
        ENDDO
       ENDIF

!  End clause for N and K

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_PARTLAYER_STERM_DN

!

SUBROUTINE LP_QUAD_GFUNCMULT                                &
          ( DO_PLANE_PARALLEL, NSTREAMS,                    & ! Input
            IB, UT, N, K, K_PARAMETERS, LAYER_PIS_CUTOFF,   & ! Input
              T_UTUP_EIGEN,     T_UTDN_EIGEN,               & ! Input
            L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,               & ! Input
               T_DELT_MUBAR,      T_UTDN_MUBAR,             & ! Input
            LP_T_DELT_MUBAR,   LP_T_UTDN_MUBAR,             & ! Input
            LP_GAMMA_M, LP_GAMMA_P, LP_INITIAL_TRANS,       & ! Input
            L_ATERM_SAVE, L_BTERM_SAVE,                     & ! Input
               UT_GMULT_UP,    UT_GMULT_DN,                 & ! Input
            LP_UT_GMULT_UP, LP_UT_GMULT_DN )                  ! Output

!  Linearization of Quadrature Green's function multiplier (off-grid only)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL, intent(in)  ::   DO_PLANE_PARALLEL

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS

!  Beam index, offgrid indices, linearization control

      INTEGER, intent(in)  ::   IB, N, UT, K, K_PARAMETERS

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  ::   LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(kind=8), intent(in)  ::  T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  ::  T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(kind=8), intent(in)  ::  T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized transmittances

      REAL(kind=8), intent(in)  ::  L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized initial transittance factors for solar beams.

      REAL(kind=8), intent(in)  ::  LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  ::  LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  ::  LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  ::  LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  ::  L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  ::  L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(in)  ::  UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  ::  UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  output arguments
!  ----------------

!  Linearized Green functions multipliers for off-grid optical depths

      REAL(kind=8), intent(out) ::  LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) ::  LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       ::   Q, AA
      REAL(kind=8)  ::  SD, SU, TD, TU, T0
      REAL(kind=8)  ::  ZX_DN, ZX_UP, ZW, WX, WDEL, L_TI
      REAL(kind=8)  ::  L_ZX_DN, L_ZX_UP, L_ZW, L_WX, L_WDEL

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            LP_UT_GMULT_UP(AA,UT,K,Q) = ZERO
            LP_UT_GMULT_DN(AA,UT,K,Q) = ZERO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Layer constant terms

      WX    = T_UTDN_MUBAR(UT,IB)
      WDEL  = T_DELT_MUBAR(N,IB)

!  For the Pseudo-spherical (average secant) multipliers
!  =====================================================

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  If the layer containing UT is the same as the varying layer

        IF ( N .EQ. K ) THEN

          DO Q = 1, K_PARAMETERS
            L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)
            L_WX   = LP_T_UTDN_MUBAR(UT,K,IB,Q)
            DO AA = 1, NSTREAMS
              ZX_DN = T_UTDN_EIGEN(AA,UT)
              ZX_UP = T_UTUP_EIGEN(AA,UT)
              ZW    = WDEL * ZX_UP
              L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
              L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
              L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP            
              SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
              TD = L_ATERM_SAVE(AA,N,Q) + LP_GAMMA_M(AA,N,K,Q) + SD
              LP_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)
              SU  = ( L_WX    - L_ZW ) / ( WX    - ZW )
              TU = L_BTERM_SAVE(AA,N,Q) + LP_GAMMA_P(AA,N,K,Q) + SU
              LP_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)
            ENDDO
          ENDDO

!  If the varying layer is above layer N

        ELSE IF ( N.GT.K ) THEN

          DO Q = 1, K_PARAMETERS
            L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)
            L_WX   = LP_T_UTDN_MUBAR(UT,K,IB,Q)
            L_TI   = LP_INITIAL_TRANS(N,K,IB,Q)
            DO AA = 1, NSTREAMS
              ZX_DN = T_UTDN_EIGEN(AA,UT)
              ZX_UP = T_UTUP_EIGEN(AA,UT)
              ZW    = WDEL * ZX_UP
              L_ZW  = L_WDEL * ZX_UP
              SD  = - L_WX / ( ZX_DN - WX )
              TD = LP_GAMMA_M(AA,N,K,Q) + SD + L_TI
              LP_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)
              SU  = ( L_WX  - L_ZW ) / ( WX - ZW )
              TU = LP_GAMMA_P(AA,N,K,Q) + SU + L_TI
              LP_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)
            ENDDO
          ENDDO

        ENDIF

!  Plane parallel case
!  ===================

      ELSE IF ( DO_PLANE_PARALLEL ) THEN

!  If the layer containing UT is the same as the varying layer

        IF ( N.EQ.K ) THEN

          DO Q = 1, K_PARAMETERS
            L_WDEL = LP_T_DELT_MUBAR(N,K,IB,Q)
            L_WX   = LP_T_UTDN_MUBAR(UT,K,IB,Q)
            DO AA = 1, NSTREAMS
              ZX_DN = T_UTDN_EIGEN(AA,UT)
              ZX_UP = T_UTUP_EIGEN(AA,UT)
              ZW    = WDEL * ZX_UP 
              L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
              L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
              L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP  
              SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
              TD = L_ATERM_SAVE(AA,N,Q) + LP_GAMMA_M(AA,N,K,Q) + SD
              LP_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)
              SU  = ( L_WX - L_ZW ) / ( WX - ZW )
              TU = L_BTERM_SAVE(AA,N,Q) + LP_GAMMA_P(AA,N,K,Q) + SU
              LP_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)
            ENDDO
          ENDDO

!  If the varying layer is above layer N, where N > K

        ELSE IF ( N.GT.K ) THEN

          DO Q = 1, K_PARAMETERS
            T0 = LP_INITIAL_TRANS(N,K,IB,Q)
            DO AA = 1, NSTREAMS            
              LP_UT_GMULT_DN(AA,UT,K,Q) = T0 * UT_GMULT_DN(AA,UT)
              LP_UT_GMULT_UP(AA,UT,K,Q) = T0 * UT_GMULT_UP(AA,UT)
            ENDDO
          ENDDO

        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_QUAD_GFUNCMULT

!

SUBROUTINE LIDORT_LP_CONVERGE                         &
           ( DO_UPWELLING, DO_SSFULL,                 & ! Input       
             DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,     & ! Input
             DO_NO_AZIMUTH, AZMFAC, LOCAL_N_USERAZM,  & ! Input
             N_USER_STREAMS, N_USER_LEVELS, NLAYERS,  & ! Input
             IBEAM, THREAD, FOURIER_COMPONENT,        & ! Input
             UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,   & ! Input
             LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Input
             PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input
             PROFILEWF )                                ! Output

!  Just upgrades the weighting function Fourier cosine series

      IMPLICIT NONE

!  Include file

      INCLUDE 'LIDORT.PARS_F90'

!  input variables
!  ---------------

!  Local flags

      LOGICAL, intent(in)  ::   DO_NO_AZIMUTH
      LOGICAL, intent(in)  ::   DO_UPWELLING
      LOGICAL, intent(in)  ::   DO_SSFULL
      LOGICAL, intent(in)  ::   DO_SSCORR_NADIR
      LOGICAL, intent(in)  ::   DO_SSCORR_OUTGOING

!  FOurier component and thread, and beam

      INTEGER, intent(in)  ::   FOURIER_COMPONENT, THREAD, IBEAM

!  Control integers

      INTEGER, intent(in)  ::   NLAYERS
      INTEGER, intent(in)  ::   N_USER_LEVELS
      INTEGER, intent(in)  ::   N_USER_STREAMS

!  Directional control

      INTEGER, intent(in)  ::   N_DIRECTIONS
      INTEGER, intent(in)  ::   WHICH_DIRECTIONS(2)

!  Bookkeeping: Offsets for geometry indexing

      INTEGER, intent(in)  ::   UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local number of azimuths and azimuth factors

      INTEGER, intent(in)       ::   LOCAL_N_USERAZM
      REAL(kind=8), intent(in)  ::  AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Linearization control

      LOGICAL, intent(in)  ::   LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, intent(in)  ::   LAYER_VARY_NUMBER ( MAXLAYERS )

!  Fourier-component Profile weighting functions at user angles

      REAL(kind=8), intent(in)  ::  PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                  MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter Profile weighting functions at user angles

      REAL(kind=8), intent(in)  ::  PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                   MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct-bounce Profile weighting functions at user angles

      REAL(kind=8), intent(in)  ::  PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                   MAX_GEOMETRIES )

!  output
!  ------

      REAL(kind=8), intent(out) ::  PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,      MAX_USER_LEVELS, &
                                                MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

!  local variables

      INTEGER       ::   I, IDIR, UT, UA, Q, W, V, N

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version 2.3)

!  Bulk/column atmospheri! weighting functions (Version 3.3)
!  ---------------------------------------------------------

!  This section newly written, 26 September 2006, installed 14 May 2007.

!  Diffuse field at all output angles

        IF ( .NOT. DO_SSFULL ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W,THREAD) = PROFILEWF_F(Q,N,UT,I,IBEAM,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W,THREAD) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    PROFILEWF(Q,N,UT,V,W,THREAD) = PROFILEWF(Q,N,UT,V,W,THREAD) + PROFILEWF_SS(Q,N,UT,V,W)
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling

        IF ( DO_UPWELLING ) THEN
         IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
          DO N = 1, NLAYERS
           IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                   PROFILEWF(Q,N,UT,V,UPIDX,THREAD) = PROFILEWF(Q,N,UT,V,UPIDX,THREAD) + PROFILEWF_DB(Q,N,UT,V)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF
        ENDIF

!  If no_azimuth, then exit

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output

        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               V = UMOFF(IBEAM,I) + UA
               PROFILEWF(Q,N,UT,V,W,THREAD) = PROFILEWF(Q,N,UT,V,W,THREAD) + &
                     PROFILEWF_F(Q,N,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

!  end Fourier clause

      ENDIF
 
      
!  Finish

      RETURN
END SUBROUTINE LIDORT_LP_CONVERGE
