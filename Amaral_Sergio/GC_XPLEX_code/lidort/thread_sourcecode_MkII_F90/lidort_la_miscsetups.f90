!$Id: lidort_la_miscsetups.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LA_DELTAMSCALE                            #
! #            LIDORT_LA_PREPTRANS                              #
! #                                                             #
! ###############################################################

SUBROUTINE LIDORT_LA_DELTAMSCALE                               &
       ( DO_DELTAM_SCALING, NLAYERS, NMOMENTS, NBEAMS, THREAD, & ! Input
         CHAPMAN_FACTORS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,  & ! Input
         DELTAU_VERT_INPUT,    L_DELTAU_VERT_INPUT,            & ! Input
         OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,            & ! Input
         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,         & ! Input
         OMEGA_MOMS, FAC1, TRUNC_FACTOR,                       & ! Input
         L_DELTAU_VERT, L_DELTAU_SLANT,                        & ! Output
         L_OMEGA_TOTAL, L_OMEGA_MOMS, L_PHASMOMS_TOTAL,        & ! Output
         L_TRUNC_FACTOR, DO_PHASFUNC_VARIATION )                 ! Output

!  Linearization of the deltam scaling

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Inputs
!  ======

!  deltam scaling flag

      LOGICAL, intent(in)  ::  DO_DELTAM_SCALING

!  Number of layers, moments and beams

      INTEGER, intent(in)  ::  NLAYERS, NBEAMS, NMOMENTS

!  Thread number

      INTEGER, intent(in)  ::  THREAD

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=8), intent(in)  :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Linearization control

      LOGICAL, intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical dephs before scaling (Threaded input)

      REAL(kind=8), intent(in)  :: OMEGA_TOTAL_INPUT ( MAXLAYERS, MAXTHREADS )
      REAL(kind=8), intent(in)  :: DELTAU_VERT_INPUT ( MAXLAYERS, MAXTHREADS )
      REAL(kind=8), intent(in)  :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Linearized Input optical dephs before scaling (Threaded input)

      REAL(kind=8), intent(in)  :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS, MAXTHREADS )
      REAL(kind=8), intent(in)  :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS, MAXTHREADS )
      REAL(kind=8), intent(in)  :: L_PHASMOMS_TOTAL_INPUT &
            ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Input optical property after delta-M scaling

      REAL(kind=8), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(kind=8), intent(in)  :: TRUNC_FACTOR(MAXLAYERS)
      REAL(kind=8), intent(in)  :: FAC1(MAXLAYERS)

!  Subroutine output arguments
!  ===========================

!  Linearized optical properties after delta-M scaling 

      LOGICAL     , intent(out)  :: DO_PHASFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=8), intent(out)  :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=8), intent(out)  :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )

      REAL(kind=8), intent(out)  :: L_PHASMOMS_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS )
      REAL(kind=8), intent(out)  :: L_OMEGA_MOMS     ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )
      REAL(kind=8), intent(out)  :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)
      
!  Linearized truncation factor

      REAL(kind=8), intent(out)  :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  local variables
!  ---------------

      REAL(kind=8)  :: BLD, DL, OF1, F, F1, UQ, EQ
      REAL(kind=8)  :: FZM, T1, L_2, L_1, VAR_L
      REAL(kind=8)  :: ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
      INTEGER       ::   N, Q, NM1, L, K, IB
      LOGICAL       :: LOOP

!  Number of moments for truncation

      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
      ELSE
        NM1 = NMOMENTS
      ENDIF

!  Set phase function linearization flag
!    Dimensioning bug, 21 March 2007. Care with using NM1

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LOOP = .TRUE.
            L = 0
            DO WHILE (LOOP.AND.L.LT.NM1)
              L = L + 1
              LOOP = ( DABS(L_PHASMOMS_TOTAL_INPUT(Q,L,N,THREAD)) &
                         .LT.1000.0*SMALLNUM )
            ENDDO
            DO_PHASFUNC_VARIATION(Q,N) = .NOT.LOOP
          ENDDO
        ENDIF
      ENDDO

!  DELTAM SCALING
!  ==============

      IF ( DO_DELTAM_SCALING ) THEN

!  start layer loop

        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            OF1 = ( ONE - FAC1(N) ) / FAC1(N)
            F   = TRUNC_FACTOR(N)
            F1  = ONE - F

            DO Q = 1, LAYER_VARY_NUMBER(N)

!  scale phase function linearization additionally

              IF ( DO_PHASFUNC_VARIATION(Q,N) ) THEN

                UQ  = L_OMEGA_TOTAL_INPUT(Q,N,THREAD)
                EQ  = L_DELTAU_VERT_INPUT(Q,N,THREAD)
                ZMQ = L_PHASMOMS_TOTAL_INPUT(Q,NM1,N,THREAD)
                FZM = F * ZMQ
                L_TRUNC_FACTOR(Q,N) = FZM
                UZQ_SCALE = ( UQ + ZMQ ) * OF1
                ZQ_SCALE = FZM / F1
                L_OMEGA_TOTAL(Q,N)      = UQ + UZQ_SCALE - ZQ_SCALE
                L_DELTAU_VERT(Q,N)      = EQ - UZQ_SCALE
                L_PHASMOMS_TOTAL(Q,0,N) = ZERO
                DO L = 1, NMOMENTS
                  ZLQ = L_PHASMOMS_TOTAL_INPUT(Q,L,N,THREAD)
                  DL  = DBLE(2*L+1)
                  BLD = PHASMOMS_TOTAL_INPUT(L,N,THREAD) / DL
                  T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                  L_PHASMOMS_TOTAL(Q,L,N) = T1 + ZQ_SCALE
                ENDDO

!  No phase function linearization
!   Zero all linearized phase function quantities now;

              ELSE

                UQ = L_OMEGA_TOTAL_INPUT(Q,N,THREAD)
                EQ = L_DELTAU_VERT_INPUT(Q,N,THREAD)
                L_TRUNC_FACTOR(Q,N) = ZERO
                UQ_SCALE = UQ * OF1
                L_OMEGA_TOTAL(Q,N) = UQ + UQ_SCALE
                L_DELTAU_VERT(Q,N) = EQ - UQ_SCALE
                DO L = 0, NMOMENTS
                  L_PHASMOMS_TOTAL(Q,L,N) = ZERO
                ENDDO

              ENDIF

!  End parameter loop

            ENDDO

!  End layer loop

          ENDIF
        ENDDO

!  Scale layer path thickness values (linearized)

        DO N = 1, NLAYERS
         DO K = 1, N
          IF ( LAYER_VARY_FLAG(K) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(K)
            L_1 = L_TRUNC_FACTOR(Q,K) *   OMEGA_TOTAL_INPUT(K,THREAD)   &
                  + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K,THREAD)
            L_2 =  - L_1 *   DELTAU_VERT_INPUT(K,THREAD) +  &
                 FAC1(K) * L_DELTAU_VERT_INPUT(Q,K,THREAD)
            DO IB = 1, NBEAMS
             L_DELTAU_SLANT(Q,N,K,IB) = L_2 * CHAPMAN_FACTORS(N,K,IB)
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDDO

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

!  Optical thickness

        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_DELTAU_VERT(Q,N)  = L_DELTAU_VERT_INPUT(Q,N,THREAD)
            ENDDO
          ENDIF
        ENDDO

!  Scattering variables just copied

        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_TRUNC_FACTOR(Q,N) = ZERO
              L_OMEGA_TOTAL(Q,N)  = L_OMEGA_TOTAL_INPUT(Q,N,THREAD)
              IF ( DO_PHASFUNC_VARIATION(Q,N) ) THEN
                DO L = 0, MAXMOMENTS
                  L_PHASMOMS_TOTAL(Q,L,N) = L_PHASMOMS_TOTAL_INPUT(Q,L,N,THREAD)
                ENDDO
              ELSE
                DO L = 0, MAXMOMENTS
                  L_PHASMOMS_TOTAL(Q,L,N) = ZERO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO

!  Scale layer path thickness values (linearized)

        DO N = 1, NLAYERS
         DO K = 1, N
          IF ( LAYER_VARY_FLAG(K) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(K)
            L_2 = L_DELTAU_VERT_INPUT(Q,K,THREAD)
            DO IB = 1, NBEAMS
             L_DELTAU_SLANT(Q,N,K,IB) = L_2 * CHAPMAN_FACTORS(N,K,IB)
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDDO

!  End delta-m clause

      ENDIF

!  phase moment-weighted OMEGA and linearizations
!  Including phase function linearization

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
            DO L = 0, NMOMENTS
              VAR_L = L_OMEGA_TOTAL   (Q,N) + L_PHASMOMS_TOTAL(Q,L,N)
              L_OMEGA_MOMS(Q,N,L) = OMEGA_MOMS(N,L) * VAR_L
             ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Finish module

      RETURN
END SUBROUTINE LIDORT_LA_DELTAMSCALE

!

SUBROUTINE LIDORT_LA_PREPTRANS                                   &
        ( DO_SOLUTION_SAVING, DO_USER_STREAMS,                   & ! Input
          NLAYERS, NSTREAMS, N_USER_STREAMS,                     & ! Input
          N_PARTLAYERS, PARTLAYERS_LAYERIDX,                     & ! Input
          QUAD_STREAMS, USER_SECANTS,                            & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                    & ! Input
          DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,               & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,        & ! Input
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,          & ! Input
          L_T_DELT_DISORDS, L_T_DISORDS_UTUP,  L_T_DISORDS_UTDN, & ! Output
          L_T_DELT_USERM,   L_T_UTUP_USERM,    L_T_UTDN_USERM )    ! Output 

!  General linearization of the Transmittances

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Inputs
!  ------

!  Flags

      LOGICAL, intent(in)  :: DO_SOLUTION_SAVING
      LOGICAL, intent(in)  :: DO_USER_STREAMS

!  Control numbers

      INTEGER, intent(in)  ::   NSTREAMS, N_USER_STREAMS

!  Layer control

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS

!  output optical depth indices

      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_USER_LEVELS)

!  User stream cosines

      REAL(kind=8), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Quadrature

      REAL(kind=8), intent(in)  :: QUAD_STREAMS( MAXSTREAMS )

!  Linearization control

      LOGICAL, intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical depths after delta-M scaling

      REAL(kind=8), intent(in)  :: DELTAU_VERT    ( MAXLAYERS )
      REAL(kind=8), intent(in)  :: PARTAU_VERT    ( MAX_PARTLAYERS )

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(kind=8), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(kind=8), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(kind=8), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(kind=8), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Linearized Optical depths

      REAL(kind=8), intent(in) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Outputs
!  -------

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(kind=8), intent(out) :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(out) :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       :: N, Q, UT, UM, I
      REAL(kind=8)  :: VD, VU, TRANS, UX, TRANS_D, TRANS_U
      REAL(kind=8)  :: L_TAU, L_TDEL, L_TD, L_TU, LDN, LUP, XT

!  Linearization of discrete ordinate transmittances
!  Only required for the solution saving option
!    (automati! for BVP telescoping)
!  Completed by R. Spurr, RTSOLUTIONS Inc., 30 August 2005

      IF ( DO_SOLUTION_SAVING ) THEN

!  whole layers

        DO N = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TAU = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)
            DO I = 1, NSTREAMS 
              L_TDEL = - L_TAU / QUAD_STREAMS(I)
              L_T_DELT_DISORDS(I,N,Q) = T_DELT_DISORDS(I,N) * L_TDEL
            ENDDO
          ENDDO
        ENDDO

!  Partial layers

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = PARTAU_VERT(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TD = - L_DELTAU_VERT(Q,N) * XT
            L_TU = - L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - XT )
            DO I = 1, NSTREAMS
              LDN = L_TD / QUAD_STREAMS(I)
              LUP = L_TU / QUAD_STREAMS(I)
              L_T_DISORDS_UTDN(I,UT,Q) = T_DISORDS_UTDN(I,UT) * LDN
              L_T_DISORDS_UTUP(I,UT,Q) = T_DISORDS_UTUP(I,UT) * LUP
            ENDDO
          ENDDO
        ENDDO

!  end solution saving option

      ENDIF

!  Linearization of Transmittance factors for User Streams
!  =======================================================

!  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

!  Whole Layer transmittance factors
!  ---------------------------------

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS = T_DELT_USERM(N,UM) * USER_SECANTS(UM) * DELTAU_VERT(N)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Partial Layer transmittance factors for off-grid optical depths
!  ---------------------------------------------------------------

      DO UT = 1, N_PARTLAYERS
        N  = PARTLAYERS_LAYERIDX(UT)
        UX = PARTAU_VERT(UT)
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS_D = T_UTDN_USERM(UT,UM) * USER_SECANTS(UM)
            TRANS_U = T_UTUP_USERM(UT,UM) * USER_SECANTS(UM)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              VD = L_DELTAU_VERT(Q,N) * UX
              VU = L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - UX )
              L_T_UTDN_USERM(UT,UM,Q) = - TRANS_D * VD
              L_T_UTUP_USERM(UT,UM,Q) = - TRANS_U * VU
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  debug

!  Finish

      RETURN
END SUBROUTINE LIDORT_LA_PREPTRANS

