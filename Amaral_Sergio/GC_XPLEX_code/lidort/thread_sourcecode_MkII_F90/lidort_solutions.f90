!$Id: lidort_solutions.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

!module lidort_solutions

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #     Homogeneous solution                                    #
! #                                                             #
! #              LIDORT_HOM_SOLUTION                            #
! #              LIDORT_HOM_EIGENTRANS                          #
! #              LIDORT_HOM_NORMS                               #
! #              LIDORT_HOM_USERSOLUTION                        #
! #              HMULT_MASTER (master)                          #
! #                                                             #
! #     Green's function solution                               #
! #                                                             #
! #              LIDORT_GBEAM_SOLUTION                          #
! #              LIDORT_GBEAM_USERSOLUTION                      #
! #              QUAD_GFUNCMULT                                 #
! #                                                             #
! ###############################################################

!contains

SUBROUTINE LIDORT_HOM_SOLUTION                                  &
         ( DO_SOLUTION_SAVING, NSTREAMS, NMOMENTS,              & ! Inputs
           GIVEN_LAYER, FOURIER, DO_LAYER_SCATTERING,           & ! Inputs
           OMEGA_MOMS, QUAD_STREAMS, QUAD_WEIGHTS,              & ! Inputs 
           PLMI_PLMJ_P, PLMI_PLMJ_M,                            & ! Inputs
           SAB, DAB, EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE, & ! output
           KEIGEN, XPOS, XNEG,                                  & ! output
           STATUS, MESSAGE, TRACE )                               ! output

!  Numerical solution of Eigenproblem.

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Solution saving flag

      LOGICAL, intent(in) :: DO_SOLUTION_SAVING

!  Number of streams and moments

      INTEGER, intent(in)  :: NSTREAMS
      INTEGER, intent(in)  :: NMOMENTS

!  Quadrature

      REAL(KIND=8), intent(in)  :: QUAD_STREAMS( MAXSTREAMS )
      REAL(KIND=8), intent(in)  :: QUAD_WEIGHTS( MAXSTREAMS )

!  Given layer index and Fourier number (inputs)

      INTEGER, intent(in)  :: GIVEN_LAYER
      INTEGER, intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Saved array involving product of OMEGA and phase function moments

      REAL(KIND=8), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  Legendre polynomial products

      REAL(KIND=8), intent(in)  :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8), intent(in)  :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(KIND=8), intent(out) :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8), intent(out) :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8), intent(out) :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8), intent(out) :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8), intent(out) :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!  (Positive) Eigenvalues

      REAL(KIND=8), intent(out) :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(KIND=8), intent(out) :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  local matrices for eigenvalue computation

      REAL(KIND=8)    :: EIGENMAT(MAXSTREAMS,MAXSTREAMS)

!  (output from Eigenpackage module ASYMTX)

      REAL(KIND=8)    :: KSQ(MAXSTREAMS), WK(MAXSTREAMS_2)
      REAL(KIND=8)    :: EVEC(MAXSTREAMS,MAXSTREAMS)
      INTEGER         :: IER
      LOGICAL         :: ASYMTX_FAILURE

!  Miscellaneous local variables

      INTEGER         :: I, J, I1, L, N, M, AA, K
      REAL(KIND=8)    :: DP, DM, SUM, FAC, KVAL, NORM, XINV
      CHARACTER*3     :: CN, CI

!  initialise exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then the eigenvectors are just the
!  discrete ordinates, the eigenvectors are unit vectors in the
!  discrete ordinate directions.

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            XINV = ONE / QUAD_STREAMS(I)
            DO AA = 1, NSTREAMS
              XPOS(I,AA,N)  = ZERO
              XPOS(I1,AA,N) = ZERO
              XNEG(I1,AA,N) = ZERO
              XNEG(I,AA,N)  = ZERO
            ENDDO
            KEIGEN(I,N)  = XINV
            XPOS(I,I,N)  = ONE
            XNEG(I1,I,N) = ONE
          ENDDO
          RETURN
        ENDIF
      ENDIF

!  Scattering solutions
!  ====================

!  Construct Eigenmatrix
!  ---------------------

!  zero the Eigenmatrix

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          EIGENMAT(I,J) = ZERO
        ENDDO
      ENDDO

!  Develop Sum and Difference matrices

      DO I = 1, NSTREAMS
        XINV = ONE/QUAD_STREAMS(I)
        DO J = 1, NSTREAMS
          FAC = XINV * HALF * QUAD_WEIGHTS(J)
          DP = ZERO
          DM = ZERO
          DO L = M, NMOMENTS
            DP = DP + PLMI_PLMJ_P(I,J,L) * OMEGA_MOMS(N,L)
            DM = DM + PLMI_PLMJ_M(I,J,L) * OMEGA_MOMS(N,L)
          ENDDO
          SAB(I,J) = FAC * ( DP + DM ) 
          DAB(I,J) = FAC * ( DP - DM )
        ENDDO
        SAB(I,I) = SAB(I,I) - XINV
        DAB(I,I) = DAB(I,I) - XINV
      ENDDO

!  Compute Eigenmatrix

      DO I = 1, NSTREAMS 
        DO J = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM + DAB(I,K) * SAB(K,J)
          ENDDO
          EIGENMAT(I,J) = SUM
        ENDDO
      ENDDO

!  save Eigenmatrix (original is destroyed by ASMTYX)

      DO I = 1, NSTREAMS 
        DO J = 1, NSTREAMS
          EIGENMAT_SAVE(I,J) = EIGENMAT(I,J)
        ENDDO
      ENDDO

!  Eigensolution package
!  ---------------------

!  Let's see how we get on with the DISORT package
!  tested 19 May 1999. Gives same results as Analytical approach.

!  second test using DGEEV (LAPACK module) - Worked against ASYMTX.
!  However, DGEEV is 1.7 - 2.0 times slower because it must look for
!  complex roots. [ASYMTX was specially written to avoid this search].
!  Also first call to DGEEV is 20 times slower. Tested June 24, 1999.
!  Conclusion stick with ASYMTX for now.

      CALL  ASYMTX ( EIGENMAT, NSTREAMS, MAXSTREAMS, MAXSTREAMS, &
                     EVEC, KSQ, IER, WK, MESSAGE, ASYMTX_FAILURE )

!  Exception handling 1

      IF ( ASYMTX_FAILURE  ) THEN
        WRITE(CN,'(I3)')N
        TRACE   ='ASYMTX error in LIDORT_HOM_SOLUTION, Layer='//CN
        STATUS  = LIDORT_SERIOUS
        RETURN
      ENDIF

! Exception handling 2

      IF ( IER.GT.0 ) THEN
        WRITE(CI,'(I3)')IER
        WRITE(CN,'(I3)')N
        MESSAGE = 'eigenvalue '//CI//' has not converged'
        TRACE   = 'ASYMTX error in LIDORT_HOM_SOLUTION, Layer='//CN
        STATUS = LIDORT_SERIOUS
        RETURN
      ENDIF

!  second test using DGEEV (LAPACK module). Here is what to use.
!        CALL DGEEV   ( 'N', 'V', NSTREAMS, EIGENMAT,
!     O        MAXSTREAMS, REAL_KSQ, IMAG_KSQ,
!     O        LEFT_EVEC, MAXSTREAMS, EVEC, MAXSTREAMS,
!     W        WORK, LWORK, LAPACK_INFO )

!  Find solution eigenvectors XPOS, XNEG for all eigenvalues
!  ---------------------------------------------------------

      DO AA = 1, NSTREAMS

!  Store positive values in output array for each layer

        KVAL = DSQRT(KSQ(AA))
        KEIGEN(AA,N) = KVAL

!  Normalize eigenvectors to 1

        NORM = ZERO
        DO I = 1, NSTREAMS
          NORM = NORM + EVEC(I,AA)*EVEC(I,AA)
        ENDDO
        NORM = DSQRT(NORM)

!  Find normalized eigenvector EIGENVEC_SAVE

        DO I = 1, NSTREAMS
          EIGENVEC_SAVE(I,AA) = EVEC(I,AA)/NORM
        ENDDO

!  Find difference eigenvector DIFVE! (Siewert's notation)

        DO I = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM - SAB(I,K) * EIGENVEC_SAVE(K,AA)
          ENDDO
          DIFVEC_SAVE(I,AA) = SUM / KVAL
          ENDDO

!  assign original evectors; first N are "DOWN", last N are "UP" (streams)

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XPOS(I,AA,N)  = HALF * ( EIGENVEC_SAVE(I,AA) + DIFVEC_SAVE(I,AA) )
          XPOS(I1,AA,N) = HALF * ( EIGENVEC_SAVE(I,AA) - DIFVEC_SAVE(I,AA) )
        ENDDO

!  debug
!   --linearization check
!        if ( do_debug_write ) then
!          write(88,'(3i4,1p20e15.7)')M,N,AA,(XPOS(I,AA,N),I=1,NSTREAMS)
!        endif

!  Use symmetry properties to set -ve eigenvectors

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XNEG(I1,AA,N) = XPOS(I,AA,N)
          XNEG(I,AA,N)  = XPOS(I1,AA,N)
        ENDDO

!  End eigenstream loop

      ENDDO

!  debug

!       k = 97
!       if (do_fdtest)k=98
!        if ( m.lt.3.and.n.gt.0) then
!        write(k,'(2i5)')M,N
!        DO AA = 1, NSTREAMS
!          WRITE(k,'(A3,I5,1p12e20.8)')'***',AA,
!     *          (XPOS(I,AA,N),I=1,NSTREAMS/2)
!        ENDDO
!        endif

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_SOLUTION

!

SUBROUTINE LIDORT_HOM_EIGENTRANS                                      &
    ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,           & ! input
      PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! input
      FOURIER, DO_LAYER_SCATTERING, DELTAU_VERT, PARTAU_VERT, KEIGEN, & ! input
      T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                 & ! Output
      T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN )                    ! Output

!  Eigenproblem, transmittance matrices

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Solution saving flag

      LOGICAL, intent(in)  :: DO_SOLUTION_SAVING

!  Number of streams and layers

      INTEGER, intent(in)  :: NSTREAMS
      INTEGER, intent(in)  :: NLAYERS

!  Number of output levels

      INTEGER, intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth mask

      LOGICAL, intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  Fourier number (input)

      INTEGER, intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Optical depths

      REAL(KIND=8), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(KIND=8), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )

!  (Positive) Eigenvalues

      REAL(KIND=8), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(KIND=8), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(KIND=8), intent(out) :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(out) :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER         :: N, M, AA, UT, UTA
      REAL(KIND=8)    :: HELP, TAU_UP, TAU_DN
    
!  Fourier

      M = FOURIER

!  Layer transmittances for the eigenvalues
!  ========================================

!  start the layer loop

      DO N = 1, NLAYERS

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND. &
              .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          DO AA = 1, NSTREAMS
            T_DELT_EIGEN(AA,N) = T_DELT_DISORDS(AA,N)
          ENDDO

!  Otherwise, get the full set of Eigenstream transmittance factors

        ELSE
          DO AA = 1, NSTREAMS
            HELP = KEIGEN(AA,N)*DELTAU_VERT(N)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DELT_EIGEN(AA,N) = ZERO
            ELSE
              T_DELT_EIGEN(AA,N) = DEXP(-HELP)
            ENDIF
          ENDDO
        ENDIF

!  end layer loop

      ENDDO

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

!  Code completed by R. Spurr, RT SOLUTIONS Inc. 30 August 2005

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        N  = PARTLAYERS_LAYERIDX(UT)

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND. &
              .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          DO AA = 1, NSTREAMS
            T_UTDN_EIGEN(AA,UT) = T_DISORDS_UTDN(AA,UT)
            T_UTUP_EIGEN(AA,UT) = T_DISORDS_UTUP(AA,UT)
          ENDDO

!  Otherwise, Compute the Eigenstream transmittance factors

        ELSE

          TAU_DN = PARTAU_VERT(UT)
          TAU_UP = DELTAU_VERT(N) - TAU_DN
          DO AA = 1, NSTREAMS
            HELP = KEIGEN(AA,N) * TAU_DN
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTDN_EIGEN(AA,UT) = ZERO
            ELSE
              T_UTDN_EIGEN(AA,UT) = DEXP(-HELP)
            ENDIF
            HELP = KEIGEN(AA,N) * TAU_UP
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTUP_EIGEN(AA,UT) = ZERO
            ELSE
              T_UTUP_EIGEN(AA,UT) = DEXP(-HELP)
            ENDIF
          ENDDO

        ENDIF

!  end loop over off-grid optical depths

       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_EIGENTRANS

!

SUBROUTINE LIDORT_HOM_NORMS                       &
         ( NSTREAMS, NLAYERS, QUAD_STRMWTS, XPOS, & ! Input
           NORM_SAVED )                             ! Output

!  Eigenproblem, solution norms for Green's function

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Number of streams and layers

      INTEGER, intent(in)  :: NSTREAMS
      INTEGER, intent(in)  :: NLAYERS

!  Quadrature

      REAL(KIND=8), intent(in)  :: QUAD_STRMWTS( MAXSTREAMS )

!  Eigenvector solutions

      REAL(KIND=8), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  outputs
!  -------

!  Saved quantities for the Green function solution

      REAL(KIND=8), intent(out) :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER         :: N, AA, J, J1
      REAL(KIND=8)    :: T1, T2, NORM
    
!  For all layers, save the norms

      DO N = 1, NLAYERS
        DO AA = 1, NSTREAMS
          NORM = ZERO
          DO J = 1, NSTREAMS
            J1 = J + NSTREAMS
            T1 = XPOS(J,AA,N)  * XPOS(J,AA,N)
            T2 = XPOS(J1,AA,N) * XPOS(J1,AA,N)
            NORM = NORM + QUAD_STRMWTS(J)*(T1-T2)
          ENDDO
          NORM_SAVED(N,AA) = NORM
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_NORMS

!

SUBROUTINE HMULT_MASTER                                                &
       ( DO_UPWELLING, DO_DNWELLING,                                   & ! Input
         NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,             & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input
         T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,                 & ! Input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                 & ! Input
         ZETA_M, ZETA_P, HMULT_1, HMULT_2,                             & ! Output
         UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD )            ! Output

!  Homogeneous solution multipliers

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine arguments
!  --------------------

!  Direction flags

      LOGICAL, intent(in)  :: DO_UPWELLING
      LOGICAL, intent(in)  :: DO_DNWELLING

!  Number of streams and layers

      INTEGER, intent(in)  :: NSTREAMS
      INTEGER, intent(in)  :: N_USER_STREAMS
      INTEGER, intent(in)  :: NLAYERS

!  Number of output levels

      INTEGER, intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth mask

      LOGICAL, intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  User stream cosines

      REAL(KIND=8), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL, intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(KIND=8), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(KIND=8), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(KIND=8), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(KIND=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  coefficient functions for user-defined angles

      REAL(KIND=8), intent(in)  :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      
!  Output = Global multipliers
!  ===========================

!  Integrated homogeneous solution multipliers, whole layer

      REAL(KIND=8), intent(out) :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(KIND=8), intent(out) :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(out) :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(out) :: UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(out) :: UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

      INTEGER         :: UM, K, N, UT, UTA
      REAL(KIND=8)    :: UDEL, SM
      REAL(KIND=8)    :: ZDEL, ZUDEL, THETA_1, THETA_2
      REAL(KIND=8)    :: UX_UP, UX_DN, ZX_UP, ZX_DN
      REAL(KIND=8)    :: THETA_DN, THETA_UP
      LOGICAL         :: STERM_EXIST(MAXLAYERS)

!  Existence flags

      DO N = 1, NLAYERS
        STERM_EXIST(N) = ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) )
      ENDDO

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      DO N = 1, NLAYERS
        IF ( STERM_EXIST(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM    = USER_SECANTS(UM)
            DO K = 1, NSTREAMS
              ZDEL    = T_DELT_EIGEN(K,N)
              ZUDEL   = ZDEL * UDEL
              THETA_2 = ONE    - ZUDEL
              THETA_1 = ZDEL - UDEL
              HMULT_1(K,UM,N) = SM * THETA_1 * ZETA_M(K,UM,N)
              HMULT_2(K,UM,N) = SM * THETA_2 * ZETA_P(K,UM,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  partial layer multipliers
!  -------------------------

      DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          IF ( STERM_EXIST(N) ) THEN

!  Upwelling

            IF ( DO_UPWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                UX_UP = T_UTUP_USERM(UT,UM)
                SM    = USER_SECANTS(UM)
                DO K = 1, NSTREAMS
                  ZDEL     = T_DELT_EIGEN(K,N)
                  ZX_UP    = T_UTUP_EIGEN(K,UT)
                  ZX_DN    = T_UTDN_EIGEN(K,UT)
                  THETA_DN = ZX_DN - ZDEL * UX_UP
                  THETA_UP = ZX_UP - UX_UP
                  UT_HMULT_UD(K,UM,UT) = SM * THETA_DN * ZETA_P(K,UM,N)
                  UT_HMULT_UU(K,UM,UT) = SM * THETA_UP * ZETA_M(K,UM,N)
                ENDDO
              ENDDO
            ENDIF

!  Downwelling

            IF ( DO_DNWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                UX_DN = T_UTDN_USERM(UT,UM)
                SM    = USER_SECANTS(UM)
                DO K = 1, NSTREAMS
                  ZDEL     = T_DELT_EIGEN(K,N)
                  ZX_UP    = T_UTUP_EIGEN(K,UT)
                  ZX_DN    = T_UTDN_EIGEN(K,UT)
                  THETA_DN = ZX_DN - UX_DN
                  THETA_UP = ZX_UP - ZDEL * UX_DN
                  UT_HMULT_DD(K,UM,UT) = SM * THETA_DN * ZETA_M(K,UM,N)
                  UT_HMULT_DU(K,UM,UT) = SM * THETA_UP * ZETA_P(K,UM,N)
                ENDDO
              ENDDO
            ENDIF

!  end loop over partial layers

          ENDIF
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE HMULT_MASTER

!

SUBROUTINE LIDORT_GBEAM_SOLUTION                                      &                      
         ( NSTREAMS, NSTREAMS_2, NMOMENTS, GIVEN_LAYER, FOURIER,      & ! input
           FLUX_FACTOR, IBEAM, DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF, & ! input
           INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,               & ! input
           QUAD_WTS, OMEGA_MOMS, KEIGEN, XPOS, T_DELT_EIGEN,          & ! input
           NORM_SAVED, PLMI_X0_P, PLMI_X0_M,                          & ! input
           GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE,        & ! Output
           CFUNC, DFUNC, AGM, BGP, GFUNC_UP, GFUNC_DN,                & ! Output
           WUPPER, WLOWER )                                             ! Output

!  Green's function beam particular integral, one layer only.
!  Uses coefficient expansion of attenuation.

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine input arguments
!  ==========================

!  Number of streams and moments

      INTEGER, intent(in)  :: NSTREAMS, NSTREAMS_2
      INTEGER, intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER, intent(in)  :: GIVEN_LAYER
      INTEGER, intent(in)  :: FOURIER

!  Solar Flux

      REAL(KIND=8), intent(in)  :: FLUX_FACTOR

!  Beam index

      INTEGER, intent(in)  :: IBEAM

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(KIND=8), intent(in)  :: QUAD_WTS( MAXSTREAMS )

!  Saved array involving product of OMEGA and phase function moments

      REAL(KIND=8), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  initial tramsittance factors for solar beams.

      REAL(KIND=8), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(KIND=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Average-secant for solar beams

      REAL(KIND=8), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  transmittance factors for +/- eigenvalues

      REAL(KIND=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Legendre polynomial products

      REAL(KIND=8), intent(in)  :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(KIND=8), intent(in)  :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  (Positive) Eigenvalues

      REAL(KIND=8), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(KIND=8), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(KIND=8), intent(in)  :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)

!  subroutine output arguments
!  ===========================

!  Green function solution
!  ***********************

!  Saved quantities for the Green function solution

      REAL(KIND=8), intent(out) :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)
      REAL(KIND=8), intent(out) :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: BGP(MAXSTREAMS,MAXLAYERS)

!  Layer C and D functions

      REAL(KIND=8), intent(out) :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: DFUNC(MAXSTREAMS,MAXLAYERS)

!  Green function Multipliers for solution
!         ( GFUNC_DN = CFUN! * ATERM_SAVE )
!         ( GFUNC_UP = DFUN! * BTERM_SAVE )
      
      REAL(KIND=8), intent(out) :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(KIND=8), intent(out) :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the boundaries

      REAL(KIND=8), intent(out) :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(KIND=8), intent(out) :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  local variables
!  ===============

      INTEGER         :: AA, L, I, I1, M, N
      REAL(KIND=8)    :: TPA, TMA, SUM_LA, SUM_LB
      REAL(KIND=8)    :: S_P_U, S_P_L, S_M_U, S_M_L

      REAL(KIND=8)    :: RHO_M, RHO_P, F1
      REAL(KIND=8)    :: CONST, SECBAR, WDEL, ZDEL, ZWDEL

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
              .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, NSTREAMS_2
          WUPPER(I,N) = ZERO
          WLOWER(I,N) = ZERO
        ENDDO
        RETURN
      ENDIF

!  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS(N,IBEAM)

!  Gamma constants

      DO AA = 1, NSTREAMS
        RHO_M = SECBAR - KEIGEN(AA,N)
        RHO_P = SECBAR + KEIGEN(AA,N)
        GAMMA_P(AA,N) = ONE / RHO_P
        GAMMA_M(AA,N) = ONE / RHO_M
      ENDDO

!  3. Optical depth integrations for the discrete ordinate solution
!     =============================================================

      WDEL    = T_DELT_MUBAR(N,IBEAM)
      DO AA = 1, NSTREAMS
        ZDEL  = T_DELT_EIGEN(AA,N)
        ZWDEL = ZDEL * WDEL
        CFUNC(AA,N)  = ( ZDEL - WDEL ) * GAMMA_M(AA,N)
        DFUNC(AA,N)  = ( ONE - ZWDEL ) * GAMMA_P(AA,N)
      ENDDO

!  4. Form quantities independent of optical depth
!     ============================================

!  set up help arrays (independent of eigenvector)

      F1 = FLUX_FACTOR / PI4
      DO I = 1, NSTREAMS
        TPA = ZERO
        TMA = ZERO
        DO L = M, NMOMENTS
          TPA = TPA + PLMI_X0_P(I,L,N,IBEAM) * OMEGA_MOMS(N,L)
          TMA = TMA + PLMI_X0_M(I,L,N,IBEAM) * OMEGA_MOMS(N,L)
        ENDDO
        DPI(I) = TPA * F1
        DMI(I) = TMA * F1
      ENDDO

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE

      DO AA = 1, NSTREAMS
        SUM_LA = ZERO
        SUM_LB = ZERO
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TPA = QUAD_WTS(I)*(DPI(I)*XPOS(I,AA,N)+DMI(I)*XPOS(I1,AA,N))
          TMA = QUAD_WTS(I)*(DMI(I)*XPOS(I,AA,N)+DPI(I)*XPOS(I1,AA,N))
          SUM_LA  = SUM_LA + TPA
          SUM_LB  = SUM_LB + TMA
        ENDDO
        ATERM_SAVE(AA,N) = SUM_LA / NORM_SAVED(N,AA)
        BTERM_SAVE(AA,N) = SUM_LB / NORM_SAVED(N,AA)
        AGM(AA,N)  = ATERM_SAVE(AA,N) * GAMMA_M(AA,N)
        BGP(AA,N)  = BTERM_SAVE(AA,N) * GAMMA_P(AA,N)
      ENDDO

!  5. Green function multipliers
!     ==========================

!  For each eigenstream

      DO AA = 1, NSTREAMS
        GFUNC_DN(AA,N) = CFUNC(AA,N) * ATERM_SAVE(AA,N) * CONST
        GFUNC_UP(AA,N) = DFUNC(AA,N) * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  6. Set particular integral from Green function expansion
!     =====================================================

!  particular integrals at lower and upper boundaries

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        S_P_U = ZERO
        S_P_L = ZERO
        S_M_U = ZERO
        S_M_L = ZERO
        DO AA = 1, NSTREAMS
          S_P_U = S_P_U + GFUNC_UP(AA,N)*XPOS(I1,AA,N)
          S_M_U = S_M_U + GFUNC_UP(AA,N)*XPOS(I,AA,N)
          S_P_L = S_P_L + GFUNC_DN(AA,N)*XPOS(I,AA,N)
          S_M_L = S_M_L + GFUNC_DN(AA,N)*XPOS(I1,AA,N)
        ENDDO
        WUPPER(I,N)  = S_P_U
        WUPPER(I1,N) = S_M_U
        WLOWER(I1,N) = S_M_L
        WLOWER(I,N)  = S_P_L
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_GBEAM_SOLUTION

!

SUBROUTINE LIDORT_HOM_USERSOLUTION                              &
         ( NSTREAMS, N_USER_STREAMS, NMOMENTS,                  & ! input
           GIVEN_LAYER, FOURIER, DO_LAYER_SCATTERING,           & ! input
           USER_SECANTS, KEIGEN, XPOS, XNEG,                    & ! input
           OMEGA_MOMS, WT_LEGP, WT_LEGM, U_LEG_P,               & ! input
           ZETA_M, ZETA_P, U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )   ! output

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine input arguments
!  --------------------------

!  Number of streams and moments

      INTEGER, intent(in)  :: NSTREAMS
      INTEGER, intent(in)  :: N_USER_STREAMS
      INTEGER, intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER, intent(in)  :: GIVEN_LAYER
      INTEGER, intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  User stream cosines

      REAL(KIND=8), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Saved array involving product of OMEGA and phase function moments

      REAL(KIND=8), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  (Positive) Eigenvalues

      REAL(KIND=8), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(KIND=8), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Polynomial-weight Legendre products

      REAL(KIND=8), intent(in)  :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8), intent(in)  :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(KIND=8), intent(in)  :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Subroutine output arguments
!  ---------------------------

!  coefficient functions for user-defined angles

      REAL(KIND=8), intent(out) :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Saved help variables

      REAL(KIND=8), intent(out) :: U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8), intent(out) :: U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(KIND=8), intent(out) :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER         :: UM, J, J1, L, N, M, AA
      REAL(KIND=8)    :: SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2
      REAL(KIND=8)    :: RHO_P, RHO_M, SECMUI, ULP

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Zeta constants (always required)

      DO UM = 1, N_USER_STREAMS
        SECMUI = USER_SECANTS(UM)
        DO AA = 1, NSTREAMS
          RHO_P = SECMUI + KEIGEN(AA,N)
          RHO_M = SECMUI - KEIGEN(AA,N)
          ZETA_P(AA,UM,N) = ONE / RHO_P
          ZETA_M(AA,UM,N) = ONE / RHO_M
        ENDDO
      ENDDO

!  If there is no scattering, zero the user solutions and exit

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = 1, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            U_XPOS(UM,AA,N) = ZERO
            U_XNEG(UM,AA,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each eigenvector

      DO AA = 1, NSTREAMS

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

        DO L = M, NMOMENTS
          SUM_POS = ZERO
          SUM_NEG = ZERO
          DO  J = 1, NSTREAMS
            J1 = J + NSTREAMS
            POS1 = XPOS(J1,AA,N) * WT_LEGP(J,L)
            POS2 = XPOS(J,AA,N)  * WT_LEGM(J,L)
            NEG1 = XNEG(J1,AA,N) * WT_LEGP(J,L)
            NEG2 = XNEG(J,AA,N)  * WT_LEGM(J,L)
            SUM_POS = SUM_POS + POS1 + POS2
            SUM_NEG = SUM_NEG + NEG1 + NEG2
          ENDDO
          U_HELP_P(AA,L) = SUM_POS
          U_HELP_M(AA,L) = SUM_NEG
        ENDDO

!  Now sum over all harmonic contributions at each user-defined stream

        DO UM = 1, N_USER_STREAMS
          SUM_POS = ZERO
          SUM_NEG = ZERO
          DO L = M, NMOMENTS
            ULP = U_LEG_P(UM,L) * OMEGA_MOMS(N,L)
            SUM_POS = SUM_POS + U_HELP_P(AA,L) * ULP
            SUM_NEG = SUM_NEG + U_HELP_M(AA,L) * ULP
          ENDDO
          U_XPOS(UM,AA,N) = SUM_POS
          U_XNEG(UM,AA,N) = SUM_NEG
        ENDDO

!  end eigenvector loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_USERSOLUTION

!

SUBROUTINE LIDORT_GBEAM_USERSOLUTION                             &
         ( DO_UPWELLING, DO_DNWELLING, N_USER_STREAMS, NMOMENTS, & ! input
           GIVEN_LAYER, FOURIER, IBEAM, FLUX_FACTOR,             & ! input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,                & ! input
           OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,                 & ! input
           U_WPOS1, U_WNEG1, W_HELP )                              ! Output

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  subroutine input arguments
!  --------------------------

!  Direction flags

      LOGICAL, intent(in)  :: DO_UPWELLING
      LOGICAL, intent(in)  :: DO_DNWELLING

!  Number of streams

      INTEGER, intent(in)  :: N_USER_STREAMS

!  Number of moments

      INTEGER, intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER, intent(in)  :: GIVEN_LAYER
      INTEGER, intent(in)  :: FOURIER
      INTEGER, intent(in)  :: IBEAM

!  Solar Flux

      REAL(KIND=8), intent(in)  :: FLUX_FACTOR

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Saved array involving product of OMEGA and phase function moments

      REAL(KIND=8), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  Legendre functions on User defined polar angles

      REAL(KIND=8), intent(in)  :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(KIND=8), intent(in)  :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Legendre polynomials, LEG0_M holds stored quantities.

      REAL(KIND=8), intent(in)  :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Subroutine output arguments
!  ---------------------------

!  Saved help variables

      REAL(KIND=8), intent(out) :: W_HELP(0:MAXMOMENTS)

!  Particular beam solutions at user-defined stream angles

      REAL(KIND=8), intent(out) :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8), intent(out) :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER         :: UM, L, N, M
      REAL(KIND=8)    :: POS1, POS2, F1

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the user solutions and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
              .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = 1, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            U_WPOS1(UM,N) = ZERO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            U_WNEG1(UM,N) = ZERO
          ENDIF
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  ====================

!  For each moment do inner sum over computational angles

      F1 = FLUX_FACTOR / PI4
      DO L = M, NMOMENTS
        W_HELP(L) = LEG0_M(L,N,IBEAM)*OMEGA_MOMS(N,L)*F1
      ENDDO

!  Now sum over all harmonic contributions at each user-defined stream

      DO UM = 1, N_USER_STREAMS
        POS1 = ZERO
        POS2 = ZERO
        DO L = M, NMOMENTS
          POS1 = POS1 + W_HELP(L)*U_LEG_P(UM,L)
          POS2 = POS2 + W_HELP(L)*U_LEG_M(UM,L)
        ENDDO
        U_WPOS1(UM,N) = POS1
        U_WNEG1(UM,N) = POS2
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_GBEAM_USERSOLUTION

!

SUBROUTINE QUAD_GFUNCMULT                                           &
          ( IB, UT, N, NSTREAMS, INITIAL_TRANS, LAYER_PIS_CUTOFF,   & ! Input
            T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR, T_UTDN_MUBAR, & ! Input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input
            UT_GMULT_UP, UT_GMULT_DN )                                ! Output

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input arguments
!  ===============

!  Number of streams

      INTEGER, intent(in)  :: NSTREAMS

!  layer and beam indices

      INTEGER, intent(in)  :: N, UT, IB

!  initial tramsittance factors for solar beams.

      REAL(KIND=8), intent(in) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(KIND=8), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(KIND=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,MAXBEAMS )
      REAL(KIND=8), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Holding arrays for Multiplier coefficients

      REAL(KIND=8), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(KIND=8), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  output arguments
!  ----------------

!  Green functions multipliers for off-grid optical depths

      REAL(KIND=8), intent(out) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(out) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

      INTEGER         :: AA
      REAL(KIND=8)    :: SD, SU, WDEL 
      REAL(KIND=8)    :: ZX_DN, ZX_UP, ZW, WX, CONST

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          UT_GMULT_DN(AA,UT) = ZERO
          UT_GMULT_UP(AA,UT) = ZERO
        ENDDO
        RETURN
      ENDIF
        
!  Layer constant terms

      WX    = T_UTDN_MUBAR(UT,IB)
      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

!  Tau integration without coefficients (average secant approximation)

      DO AA = 1, NSTREAMS
        ZX_DN = T_UTDN_EIGEN(AA,UT)
        ZX_UP = T_UTUP_EIGEN(AA,UT)
        ZW    = WDEL * ZX_UP
        SD =  ( ZX_DN - WX ) * GAMMA_M(AA,N)
        SU =  ( WX    - ZW ) * GAMMA_P(AA,N)
        UT_GMULT_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
        UT_GMULT_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUAD_GFUNCMULT

!  End Module

!end module lidort_solutions

