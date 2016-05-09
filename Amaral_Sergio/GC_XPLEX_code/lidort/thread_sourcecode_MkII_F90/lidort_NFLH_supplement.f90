!$Id: lidort_NFLH_supplement.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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
! #  This Version :   3.5                                   #
! #  Release Date :   June 2010                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.1)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED code             (3.5)         #
! #                                                         #
! #       NFLH Supplement to V3.5 added 07 April 2010       #
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
! #            LC_BVP_SOLUTION_MASTER_NFLH                      #
! #            LC_BVP_COLUMN_SETUP_NFLH                         #
! #                                                             #
! #            GET_LC_BOASOURCE_NFLH                            #
! #                                                             #
! #            LIDORT_LC_SSCORR_OUTGOING_NFLH (master)          #
! #                 LCH_og_integration_up_NF                    #
! #                 LCH_og_integration_dn_NF                    #
! #                 LCH_og_attenuations_NF                      #
! #                                                             #
! #            CHAPMAN_FUNCTION_NF                              #
! #            CHAPMAN_FUNCTION_NF_PLUS                         #
! #                                                             #
! #            sphergeom_NF                                     #
! #            LH_sphergeom_NF                                  #
! #                                                             #
! #            LIDORT_LC_PREPTRANS_NFLH                         #
! #                                                             #
! ###############################################################

SUBROUTINE LC_BVP_SOLUTION_MASTER_NFLH                                 &
            ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,               & ! Inputs
              DO_BRDF_SURFACE, DO_PLANE_PARALLEL,                      & ! Inputs
              NLAYERS, NSTREAMS, NSTREAMS_2,                           & ! Inputs
              NTOTAL, N_SUBDIAG, N_SUPDIAG, N_WEIGHTFUNCS,             & ! Inputs
              FOURIER_COMPONENT, IBEAM, QUAD_STRMWTS, SURFACE_FACTOR,  & ! Inputs
              ALBEDO, BRDF_F, LC_DIRECT_BEAM,                          & ! Inputs
              DO_LAYER_SCATTERING, INITIAL_TRANS, LAYER_PIS_CUTOFF,    & ! Inputs
              T_DELT_EIGEN, T_DELT_MUBAR, XPOS, XNEG,                  & ! Inputs
              GAMMA_M, GAMMA_P, GFUNC_UP, GFUNC_DN,                    & ! Inputs
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                        & ! Inputs
              LCON, MCON, LCONMASK, MCONMASK,                          & ! Inputs
              L_ATERM_SAVE, L_BTERM_SAVE,                              & ! Inputs
              L_KEIGEN, L_XPOS, L_XNEG, L_T_DELT_EIGEN,                & ! Inputs
              LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,    & ! Inputs
              LC_GAMMA_M, LC_GAMMA_P, L_WUPPER, L_WLOWER,              & ! Outputs
              NCON, PCON, NCON_XVEC, PCON_XVEC,                        & ! Output
              STATUS, MESSAGE, TRACE )                                   ! Output

!  Linearization of the Boundary Problem Solution

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  inclusion flags

      LOGICAL, intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, intent(in)  :: DO_INCLUDE_SURFACE

!  BRDF flag

      LOGICAL, intent(in)  :: DO_BRDF_SURFACE

!  Flag

      LOGICAL, intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER, intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Number of layers

      INTEGER, intent(in)  :: NLAYERS

!  BVProblem Band matrix control

      INTEGER, intent(in)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Fourier component and beam number

      INTEGER, intent(in)  :: FOURIER_COMPONENT, IBEAM

!  Quadrature input

      REAL(kind=8), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface factor = 1+delta(m,0). Albedo

      REAL(kind=8), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, multiplied by strmwts
!    ( New code, 23 March 2010 )
!    incident quadrature directions, reflected quadrature streams

      REAL(kind=8), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Linearization control

      INTEGER, intent(in)  :: N_WEIGHTFUNCS

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  linearized Direct beam solutions

      REAL(kind=8), intent(in)  :: LC_DIRECT_BEAM ( MAX_ATMOSWFS, MAXSTREAMS, MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=8), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Green function Multipliers for solution
      
      REAL(kind=8), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Matrix, Band-matrix

      REAL(kind=8), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(kind=8), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER, intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER, intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER, intent(in)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER, intent(in)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(kind=8), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearized (Positive) Eigenvalues

      REAL(kind=8), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(kind=8), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  linearization of pseudo-spherical approximation

      REAL(kind=8), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of solar beam layer transmittances

      REAL(kind=8), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(out) :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(out) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(kind=8)  :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(kind=8)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  error tracing variables

      INTEGER       :: INFO
      CHARACTER*3   :: CI

!  Other local help variables 

      INTEGER       :: I, I1, Q, N, AA

!  initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '
      
!  Linearization of the regular BVP case
!  =====================================

!  Set up the column vectors for Column linearizations
!  ---------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'

      CALL LC_BVP_COLUMN_SETUP_NFLH                               &
       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,               & ! Inputs
         DO_BRDF_SURFACE, DO_PLANE_PARALLEL,                      & ! Inputs
         NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL, N_WEIGHTFUNCS,    & ! Inputs
         FOURIER_COMPONENT, IBEAM, QUAD_STRMWTS,                  & ! Inputs
         SURFACE_FACTOR, ALBEDO, BRDF_F, LC_DIRECT_BEAM,          & ! Inputs
         DO_LAYER_SCATTERING, INITIAL_TRANS, LAYER_PIS_CUTOFF,    & ! Inputs
         XPOS, XNEG, T_DELT_EIGEN, T_DELT_MUBAR,                  & ! Inputs
         GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, LCON, MCON,        & ! Inputs
         L_ATERM_SAVE, L_BTERM_SAVE,                              & ! Inputs
         L_KEIGEN, L_XPOS, L_XNEG, L_T_DELT_EIGEN,                & ! Inputs
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,    & ! Inputs
         LC_GAMMA_M, LC_GAMMA_P, L_WUPPER, L_WLOWER,              & ! Outputs
         COL2_WF, SCOL2_WF )                                        ! Outputs

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_WEIGHTFUNCS, &
             BANDMAT2, MAXBANDTOTAL, IPIVOT,  COL2_WF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for DGBTRS call in LC_BVP_SOLUTION_MASTER_NFLH'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(I,N,Q) = COL2_WF(LCONMASK(I,N),Q)
              PCON(I,N,Q) = COL2_WF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for 1-layer: DGETRS call in LC_BVP_SOLUTION_MASTER_NFLH'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_WEIGHTFUNCS
            NCON(I,N,Q) = SCOL2_WF(I,Q)
            PCON(I,N,Q) = SCOL2_WF(I1,Q)
          ENDDO
        ENDDO

      ENDIF

!  linearized BVP results
!  ======================

!  Associated quantities

      DO N = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          DO AA = 1, NSTREAMS
             DO Q = 1, N_WEIGHTFUNCS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  finish

      RETURN
END SUBROUTINE LC_BVP_SOLUTION_MASTER_NFLH

!

SUBROUTINE LC_BVP_COLUMN_SETUP_NFLH                               &
       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,               & ! Inputs
         DO_BRDF_SURFACE, DO_PLANE_PARALLEL,                      & ! Inputs
         NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL, N_WEIGHTFUNCS,    & ! Inputs
         FOURIER_COMPONENT, IPARTIC, QUAD_STRMWTS,SURFACE_FACTOR, & ! Inputs
         ALBEDO, BRDF_F, LC_DIRECT_BEAM,                          & ! Inputs
         DO_LAYER_SCATTERING, INITIAL_TRANS, LAYER_PIS_CUTOFF,    & ! Inputs
         XPOS, XNEG, T_DELT_EIGEN, T_DELT_MUBAR,                  & ! Inputs
         GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, LCON, MCON,        & ! Inputs
         L_ATERM_SAVE, L_BTERM_SAVE,                              & ! Inputs
         L_KEIGEN, L_XPOS, L_XNEG, L_T_DELT_EIGEN,                & ! Inputs
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,    & ! Inputs
         LC_GAMMA_M, LC_GAMMA_P, L_WUPPER, L_WLOWER,              & ! Outputs
         COL2_WF, SCOL2_WF )                                        ! Outputs

!  Linearized column vector setup

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  input arguments
!  ---------------

!  Flag

      LOGICAL, intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER, intent(in)  :: NSTREAMS, NSTREAMS_2

!  Number of layers

      INTEGER, intent(in)  :: NLAYERS, NTOTAL

!  surface and direct beam flags

      LOGICAL, intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL, intent(in)  :: DO_BRDF_SURFACE
      LOGICAL, intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Fourier component

      INTEGER, intent(in)  :: FOURIER_COMPONENT

!  surface factors

      REAL(kind=8), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, multiplied by strmwts
!    ( New code, 23 March 2010 )
!    incident quadrature directions, reflected quadrature streams

      REAL(kind=8), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Linearization control

      INTEGER, intent(in)  :: N_WEIGHTFUNCS

!  Beam solution number

      INTEGER, intent(in)  :: IPARTIC

!  Quadrature input

      REAL(kind=8), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Local flags for the solution saving option

      LOGICAL, intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Linearized Direct beam solutions

      REAL(kind=8), intent(in)  :: LC_DIRECT_BEAM ( MAX_ATMOSWFS, MAXSTREAMS, MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=8), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Transmittance factors for average secant stream

      REAL(kind=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Green function Multipliers for solution
      
      REAL(kind=8), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(kind=8), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearized (Positive) Eigenvalues

      REAL(kind=8), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(kind=8), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  linearization of pseudo-spherical approximation

      REAL(kind=8), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of solar beam layer transmittances

      REAL(kind=8), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8), intent(out) :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(kind=8), intent(out) :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  help variables

      INTEGER       :: Q, AA, N, N1, I, I1, CM, C0
      REAL(kind=8)  :: CPOS, CNEG, L_HOM, L_PARTIC
      REAL(kind=8)  :: L_HOMD, L_HOMU
      LOGICAL       :: MODIFIED_BOUNDARY

!  Local linearized reflectance arrays

      REAL(kind=8)  :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  initialise
!  ----------

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  Top of first layer (TOA), UPPER boundary condition
!  --------------------------------------------------

      N = 1

!  Get the linearized beam solution for the first layer

      CALL LC_BEAMSOLUTION_NEQK                                 &
        ( DO_PLANE_PARALLEL, NSTREAMS, NSTREAMS_2,              & ! Inputs
          FOURIER_COMPONENT, IPARTIC, N, N_WEIGHTFUNCS,         & ! Inputs
          INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING, & ! Inputs
          T_DELT_EIGEN,  T_DELT_MUBAR, XPOS,                    & ! Inputs
          GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Inputs
          LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_AVERAGE_SECANT, & ! Inputs
          L_KEIGEN, L_T_DELT_EIGEN, L_XPOS,                     & ! Inputs
          L_ATERM_SAVE, L_BTERM_SAVE,                           & ! Inputs
          LC_GAMMA_M, LC_GAMMA_P, L_WUPPER, L_WLOWER )            ! Outputs

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO I = 1, NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          L_PARTIC = - L_WUPPER(I,N,Q)
          L_HOM    = ZERO
          DO AA = 1, NSTREAMS
            CPOS = L_XPOS(I,AA,N,Q)
            CNEG =   T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) +  &
                  L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
            L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
          ENDDO
          COL2_WF(I,Q) = L_PARTIC - L_HOM
        ENDDO
      ENDDO
 
!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer

        CALL LC_BEAMSOLUTION_NEQK                               &
        ( DO_PLANE_PARALLEL, NSTREAMS, NSTREAMS_2,              & ! Inputs
          FOURIER_COMPONENT, IPARTIC, N1, N_WEIGHTFUNCS,        & ! Inputs
          INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING, & ! Inputs
          T_DELT_EIGEN,  T_DELT_MUBAR, XPOS,                    & ! Inputs
          GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Inputs
          LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_AVERAGE_SECANT, & ! Inputs
          L_KEIGEN, L_T_DELT_EIGEN, L_XPOS,                     & ! Inputs
          L_ATERM_SAVE, L_BTERM_SAVE,                           & ! Inputs
          LC_GAMMA_M, LC_GAMMA_P, L_WUPPER, L_WLOWER )            ! Outputs

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_WEIGHTFUNCS
            L_HOMD = ZERO
            L_HOMU = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N1,Q)
              CNEG =   T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + &
                     L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
              L_HOMU = L_HOMU + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + &
                     L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOMD = L_HOMD + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            L_PARTIC      = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
            L_HOM         = L_HOMU - L_HOMD
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
          ENDDO
        ENDDO

!  End layer

      ENDDO

!  LOWER layer 
!  -----------

      N = NLAYERS
      MODIFIED_BOUNDARY = .TRUE.

!  get the linearized downward-reflected term

      CALL LC_BVP_SURFACE_SETUP                                     &
      ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, NLAYERS,     & ! Inputs
        MODIFIED_BOUNDARY, FOURIER_COMPONENT, SURFACE_FACTOR,       & ! Inputs
        ALBEDO, BRDF_F, N_WEIGHTFUNCS, QUAD_STRMWTS, T_DELT_EIGEN,  & ! Inputs
        L_T_DELT_EIGEN, XPOS, L_XPOS, L_XNEG, L_WLOWER,             & ! Inputs
        R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )                           ! Outputs

!  Compute the solution

      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      DO I = 1, NSTREAMS
        CM = C0 + I
        I1 = I + NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)
          L_HOM    = ZERO
          DO AA = 1, NSTREAMS
            CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
            CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
            CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
            L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          ENDDO
          COL2_WF(CM,Q) = - L_PARTIC - L_HOM
        ENDDO
      ENDDO

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          DO Q = 1, N_WEIGHTFUNCS
           COL2_WF(CM,Q) = COL2_WF(CM,Q) + LC_DIRECT_BEAM(Q,I,IPARTIC)
          ENDDO
        ENDDO
      ENDIF

!  Single layer

      IF ( NLAYERS .eq. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_WEIGHTFUNCS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  finish

      RETURN
END SUBROUTINE LC_BVP_COLUMN_SETUP_NFLH

!

SUBROUTINE GET_LC_BOASOURCE_NFLH                                  &
         ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,             & ! Inputs
           DO_BRDF_SURFACE, DO_USER_STREAMS,                      & ! Inputs
           NLAYERS, NSTREAMS, N_USER_STREAMS,                     & ! Inputs
           FOURIER_COMPONENT, IBEAM, K_PARAMETERS,                & ! Inputs
           SURFACE_FACTOR, ALBEDO, USER_BRDF_F, QUAD_STRMWTS,     & ! Inputs
           LC_USER_DIRECT_BEAM,  L_XPOS, L_XNEG,                  & ! Inputs
           LCON, LCON_XVEC, T_DELT_EIGEN,                         & ! Inputs
           MCON, NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_WLOWER,  & ! Inputs
           L_BOA_MSSOURCE, L_BOA_DBSOURCE )                         ! Outputs

!  Linearized Bottom of the atmosphere source term

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL, intent(in)  :: DO_BRDF_SURFACE
      LOGICAL, intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, intent(in)  :: DO_USER_STREAMS

!  control integers

      INTEGER, intent(in)  :: NLAYERS, NSTREAMS, N_USER_STREAMS

!  surface multiplier, albedo and Fourier/beam indices

      REAL(kind=8), intent(in)  :: SURFACE_FACTOR, ALBEDO
      INTEGER, intent(in)  :: FOURIER_COMPONENT
      INTEGER, intent(in)  :: IBEAM

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected user streams

      REAL(kind=8), intent(in)  :: USER_BRDF_F ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Quadrature

      REAL(kind=8), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  linearization control

      INTEGER, intent(in)  :: K_PARAMETERS

!  Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(kind=8), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Direct beam solutions

      REAL(kind=8), intent(in)  :: LC_USER_DIRECT_BEAM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized Eigensolutions

      REAL(kind=8), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(kind=8), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(kind=8), intent(out) :: L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(out) :: L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       :: M, N, J, I, UM, AA, Q, IB
      REAL(kind=8)  :: L_DOWN(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  :: REFLEC, KMULT
      REAL(kind=8)  :: SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

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
          L_DOWN(I,Q) = SHOM + L_WLOWER(I,N,Q)
          L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_STRMWTS(I)
        ENDDO
      ENDDO

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

!  Add direct beam if flagged

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_DBSOURCE(UM,Q) = LC_USER_DIRECT_BEAM(Q,UM,IB)
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_LC_BOASOURCE_NFLH

!

SUBROUTINE LIDORT_LC_SSCORR_OUTGOING_NFLH                              &
       ( DO_UPWELLING, DO_DNWELLING, DO_COLUMN_LINEARIZATION,          & ! Inputs
         DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING,                      & ! Inputs
         NLAYERS, NMOMENTS_INPUT, NMOMENTS, NBEAMS,                    & ! Inputs
         N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                & ! Inputs
         N_TOTALCOLUMN_WFS, LHMASK_COLUMN_WFS, THREAD, SSFLUX,         & ! Inputs
         BEAM_SZAS_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,    & ! Inputs
         N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  & ! Inputs
         UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, DO_PHFUNC_VARIATION,  & ! Inputs
         HEIGHT_GRID, LH_HEIGHT_GRID, EARTH_RADIUS,                    & ! Inputs
         TRUNC_FACTOR, DELTAU_VERT, L_DELTAU_VERT,                     & ! Inputs
         OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,                    & ! Inputs
         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                 & ! Input
         UP_LOSTRANS, BOA_ATTN, L_UP_LOSTRANS, L_BOA_ATTN,             & ! Output
         INTENSITY_SS, COLUMNWF_SS,                                    & ! Output
          STATUS, MESSAGE, TRACE)                                        ! Output

!  Single scatter exact calculation for the outgoing LOS
!   This one with optional linearizations
!         - NEW for version 3.2

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft, January 23rd 2007
!   Validated against TOMRAD, 29 March 2007.
!   Partial layer output added September 2007

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  directional control

      LOGICAL, intent(in)  :: DO_UPWELLING
      LOGICAL, intent(in)  :: DO_DNWELLING

!  Logical control for overall linearizations

      LOGICAL, intent(in)  :: DO_COLUMN_LINEARIZATION

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations

      LOGICAL, intent(in)  :: DO_SSCORR_TRUNCATION

!  Deltam scaling flag

      LOGICAL, intent(in)  :: DO_DELTAM_SCALING

!  Thread number

      INTEGER, intent(in)  :: THREAD

!  FLux

      REAL(kind=8), intent(in)  :: SSFLUX

!  number of computational layers

      INTEGER, intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER, intent(in)  :: NBEAMS

!  number of Legendre phase function expansion moments

      INTEGER, intent(in)  :: NMOMENTS_INPUT
      INTEGER, intent(in)  :: NMOMENTS

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER, intent(in)  :: N_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER, intent(in)  :: N_USER_STREAMS

!  Number of user levels

      INTEGER, intent(in)  :: N_USER_LEVELS

!  Control for atmospheric linearizations

      INTEGER, intent(in)  :: N_TOTALCOLUMN_WFS

!  Column linearization control inputs

      logical, intent(in)  :: LHmask_column_wfs ( MAX_ATMOSWFS )

!  Adjusted geometries

      REAL(kind=8), intent(in)  :: BEAM_SZAS_ADJUST    (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      REAL(kind=8), intent(in)  :: USER_RELAZMS_ADJUST (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      REAL(kind=8), intent(in)  :: USER_ANGLES_ADJUST  (MAX_USER_STREAMS)

!   Offsets for geometry indexing

      INTEGER, intent(in)  :: N_GEOMETRIES
      INTEGER, intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Layer masks for doing integrated source terms

      LOGICAL, intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  output optical depth masks and indices
!  off-grid optical depths (values, masks, indices)

      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  multilayer optical property (bulk) inputs

      REAL(kind=8), intent(in)  :: OMEGA_TOTAL_INPUT ( MAXLAYERS, MAXTHREADS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(kind=8), intent(in)  :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(kind=8), intent(in)  :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(kind=8), intent(in)  :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS,MAXTHREADS)

!  Saved array for truncation factor 

      REAL(kind=8), intent(in)  :: TRUNC_FACTOR ( MAXLAYERS )

!  Linearized input optical properties after delta-M scaling

      LOGICAL     , intent(in)  :: DO_PHFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )

!  multilayer optical depth inputs

      REAL(kind=8), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(kind=8), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Earth radius and height grid

      REAL(kind=8), intent(in)  :: EARTH_RADIUS
      REAL(kind=8), intent(in)  :: HEIGHT_GRID    ( 0:MAXLAYERS )
      REAL(kind=8), intent(in)  :: LH_HEIGHT_GRID ( MAX_ATMOSWFS, 0:MAXLAYERS )

!  Output arguments
!  ----------------

!  single scatter results

      REAL(kind=8), intent(out) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Transmittances. Output required for later on.

      REAL(kind=8), intent(out) :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(kind=8), intent(out) :: BOA_ATTN(MAX_GEOMETRIES)

!  COLUMN WF single scatter results

      REAL(kind=8), intent(out) :: COLUMNWF_SS ( MAX_ATMOSWFS, MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Linearized output. Output required for later on.

      REAL(kind=8), intent(out) :: L_UP_LOSTRANS  (MAXLAYERS, MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(kind=8), intent(out) :: L_BOA_ATTN(MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Exception handling. Updated, 18 May 2010

      integer      , intent(out) :: status
      character*(*), intent(out) :: message, trace

!  Outgoing sphericity stuff
!  -------------------------

!  Downwelling tranmsittances + linearizations

      REAL(kind=8) :: DN_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(kind=8) :: L_DN_LOSTRANS(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Multipliers

      REAL(kind=8) :: UP_MULTIPLIERS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(kind=8) :: DN_MULTIPLIERS(MAXLAYERS,MAX_GEOMETRIES)

!  All linearized multipliers

      REAL(kind=8) :: L_UP_MULTIPLIERS (MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)
      REAL(kind=8) :: L_DN_MULTIPLIERS (MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Geometry routine inputs and outputs
!  -----------------------------------

!  local angles

      real(kind=8) :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      integer      :: ntraverse   (0:maxlayers)
      real(kind=8) :: sunpaths    (0:maxlayers,maxlayers)
      real(kind=8) :: LH_sunpaths (0:maxlayers,maxlayers)
      real(kind=8) :: radii       (0:maxlayers)
      real(kind=8) :: alpha_all   (0:maxlayers)
      real(kind=8) :: LH_alpha_all(0:maxlayers)

!  Other (incidental) geometrical output

      real(kind=8) :: cosscat_up (0:maxlayers)
      real(kind=8) :: cosscat_dn (0:maxlayers)

!  Extinction and var tms

      real(kind=8) :: extinction   (maxlayers)
      real(kind=8) :: L_extinction (maxlayers,max_atmoswfs)
      real(kind=8) :: var_tms      (maxlayers,max_atmoswfs)

!  local variables
!  ---------------

      real(kind=8) :: attn   ( 0:maxlayers )
      real(kind=8) :: L_attn ( 0:maxlayers, max_atmoswfs )

!  Saved Legendre polynomials

      REAL(kind=8) :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      REAL(kind=8) :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(kind=8) :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      REAL(kind=8) :: SSFDEL ( MAXLAYERS )
      REAL(kind=8) :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  Exact Phase function calculations

      REAL(kind=8) :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      REAL(kind=8) :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Exact Phase function calculations

      REAL(kind=8) :: L_EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8) :: L_EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS,MAX_ATMOSWFS)

!  Cumulative single scatter source terms

      REAL(kind=8) :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=8) :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  Linearized Cumulative single scatter source terms

      REAL(kind=8) :: L_SS_CUMSOURCE(MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Local array of height derivatives

      real(kind=8) :: local_lh_height_grid(0:maxlayers)

!  Indices

      INTEGER      :: N, NUT, NSTART, NUT_PREV, NLEVEL, L
      INTEGER      :: UTA, UM, IA, NC, IB, V, NM1
      integer      :: k, q, T

!  help variables (double precision)

      REAL(kind=8) :: FINAL_SOURCE, SS_LAYERSOURCE
      REAL(kind=8) :: SSCORRECTION, COSSCAT, LEGPOLY, VAR1
      REAL(kind=8) :: DF1(0:MAXMOMENTS_INPUT)
      REAL(kind=8) :: DF2(0:MAXMOMENTS_INPUT)
      REAL(kind=8) :: UVAR, AVAR, FT1, FT2, LSS, L_PHASMOM, L_HELP
      REAL(kind=8) :: L_FINAL_SOURCE, L_SSCORRECTION, DNM1, EXT
      REAL(kind=8) :: FACT, DNL1, FDNL1, GK11, HELP, HELP1, HELP2

!  Help variables for exception handling

      character*3    :: cv
      logical        :: fail

!  Set up operations
!  -----------------

!  set exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Thread number

      T = THREAD

!  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N,T)
          TMS(N) = OMEGA_TOTAL_INPUT(N,T) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N,T)
        ENDIF
      ENDDO

!  linearizations of the TMS factors
!  ---------------------------------

!  Create Linearized TMS factor for each layer
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
!    Distinguish between Deltam case or not.
!  Only if required

      IF ( DO_COLUMN_LINEARIZATION ) THEN

!  layer loop

       DO K = 1, NLAYERS

!  Deltam_scaling. Extra contributions if phase function moments are varying.
    
         IF ( DO_DELTAM_SCALING ) THEN
          NM1 = NMOMENTS+1
          FT1 = TRUNC_FACTOR(K) * TMS(K)
          FT2 = ONE + FT1
          DO Q = 1, N_TOTALCOLUMN_WFS
           IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
            AVAR = UVAR + L_PHASMOMS_TOTAL_INPUT(Q,NM1,K,T)
            VAR_TMS(K,Q) = UVAR + FT1 * AVAR
           ELSE
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
            VAR_TMS(K,Q) = UVAR * FT2
           ENDIF
          ENDDO
         ENDIF

!  No delta-M scaling, just copy

         IF ( .NOT. DO_DELTAM_SCALING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
           UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
           VAR_TMS(K,Q) = UVAR
          ENDDO
         ENDIF

!  End layer loop

       ENDDO

!  End linearization clause

      ENDIF

!  Additional Delta-M scaling
!  --------------------------

!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN

!  Basic scaling

        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N,T) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO

!  Linearization, Only change if SCATMAT variation is set

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
            DO Q = 1, N_TOTALCOLUMN_WFS
              IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                FT1 = ONE / ( ONE- SSFDEL(K) )
                L_SSFDEL(K,Q) = PHASMOMS_TOTAL_INPUT(NM1,K,T) * &
                        L_PHASMOMS_TOTAL_INPUT(Q,NM1,K,T) / DNM1
                VAR_TMS(K,Q) = VAR_TMS(K,Q) - FT1 * L_SSFDEL(K,Q)
              ENDIF
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Create extinctions
!  ------------------

!  Local height variable array

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO q = 1, n_totalcolumn_wfs
          if ( LHMASK_COLUMN_WFS(Q) ) THEN
            do k = 0, nlayers
              local_lh_height_grid(k) = LH_HEIGHT_GRID(Q,K)
            enddo
          endif
        enddo
      ENDIF

!  Use basic definitions

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

!  Linearized extinctions

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO K = 1, NLAYERS
          HELP = 1.0d0 / (HEIGHT_GRID(K-1) - HEIGHT_GRID(K))
          EXT  = EXTINCTION(K)
          DO Q = 1, N_TOTALCOLUMN_WFS
            IF ( LHMASK_COLUMN_WFS(Q) ) THEN
              L_HELP = (LH_HEIGHT_GRID(Q,K-1)-LH_HEIGHT_GRID(Q,K))*HELP
              L_EXTINCTION(K,Q) = EXT * ( L_DELTAU_VERT(Q,K)- L_HELP )
            ELSE
              L_EXTINCTION(K,Q) = L_DELTAU_VERT(Q,K) * EXT
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Source function contributions
!  =============================

!  Start the main loop over all solar and viewing geometries

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = BEAM_SZAS_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = UMOFF(IB,UM) + IA
    
!  Call to geometry

!  Call to geometry routine, path distances etc....

            call LH_sphergeom_NF                           &
       ( maxlayers, nlayers, earth_radius, v,              & ! Inputs
         height_grid, local_lh_height_grid,                & ! Inputs
         alpha_boa, theta_boa, phi_boa,                    & ! Inputs
         ntraverse,  cosscat_up, cosscat_dn, radii,        & ! Outputs
         sunpaths, alpha_all,  LH_sunpaths,  LH_alpha_all, & ! Outputs
         fail, message )                                     ! Outputs

!  Excpetion handling. Updated 18 May 2010

            if ( fail ) then
              write(cv,'(I3)')v
              trace =  'Error from LH_sphergeom_NF, geometry # '//cv
              status = lidort_serious
              return
            endif

!  debug geometry
!            do n = 0, nlayers
!              write(14,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
!            enddo
!            pause

!  Get the attenuations

            call LCH_og_attenuations_NF                 & 
           ( maxlayers, max_atmoswfs, nlayers,          & ! Inputs
             n_totalcolumn_wfs, LHmask_column_wfs,      & ! Inputs
             extinction, L_extinction,                  & ! Inputs
             sunpaths, LH_sunpaths, ntraverse,          & ! Inputs 
             attn, L_attn )                               ! outputs

!  debug
!              if ( v.eq.10) then
!               do n = 0, nlayers
!                if (n_totalcolumn_wfs .ne.0 ) then
!                 write(88,*)n,attn(n),L_attn(n,3)
!                else
!                 write(88,*)n,attn(n)
!                endif
!               enddo
!             endif

!  Set BOA attenuation output

            boa_attn(v) = attn(nlayers)
            do q = 1, n_totalcolumn_wfs
              l_boa_attn(q,v) = L_attn(nlayers,q)
            enddo
      
!  Upwelling calculation
!  ---------------------

            IF ( DO_UPWELLING ) THEN

!  Multipliers and transmittances + linearizations

              call LCH_og_integration_up_NF                          &
            ( maxlayers, max_atmoswfs,                               & ! Inputs
              nlayers, n_totalcolumn_wfs, LHmask_column_wfs,         & ! Inputs
              extinction, l_extinction, radii, local_lh_height_grid, & ! Inputs
              alpha_all, LH_alpha_all, attn, l_attn,                 & ! Inputs
              up_multipliers(1,v), up_lostrans(1,v),                 & ! Outputs
              l_up_multipliers(1,1,v), l_up_lostrans(1,1,v) )          ! Outputs

!  DEBUG
!              if ( v.eq.10) then
!               do n = 1, nlayers
!                if (n_totalcolumn_wfs .ne.0 ) then
!                 write(8,*)n,up_lostrans(n,v),up_multipliers(n,v)
!               write(8,*)n,l_up_lostrans(n,3,v),l_up_multipliers(n,3,v)
!                else
!                   write(8,*)n,up_lostrans(n,v),up_multipliers(n,v)
!                endif
!               enddo
!             endif

!  DEbug
!              write(35,*)v
!              do n = 1, nlayers
!                write(35,'(1p2e18.10)')  up_lostrans(n,v),up_multipliers(n,v)
!              enddo
!              if ( v.eq.8)pause

!              if ( do_column_linearization ) then
!              do n = 1, nlayers
!          write(36,'(2i4,1p20e18.10)')v,n,  up_lostrans(n,v),up_multipliers(n,v), &
!                  l_up_lostrans(n,3,v),l_up_multipliers(n,3,v)
!              enddo
!              else
!                do n = 1, nlayers
!                  write(35,'(2i4,1p2e18.10)')v,n, up_lostrans(n,v),up_multipliers(n,v)
!                 enddo
!              endif
!         pause'gronkdon'

!  legendre polynomials

              COSSCAT = COSSCAT_UP(NLAYERS)
              SS_PLEG_UP(V,1,0) = ONE
              SS_PLEG_UP(V,1,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_UP(V,1,L) = DF1(L) * SS_PLEG_UP(V,1,L-1) * COSSCAT - DF2(L) * SS_PLEG_UP(V,1,L-2)
              ENDDO

!  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(N) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(N) * DNL1
                      FACT  = ONE - SSFDEL(N)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,N,T)-FDNL1) / FACT
                    ELSE
                      GK11 = PHASMOMS_TOTAL_INPUT(L,N,T)
                    ENDIF 
                    LEGPOLY = SS_PLEG_UP(V,1,L)
                    HELP = HELP + GK11 * LEGPOLY
                  ENDDO
                  EXACTSCAT_UP(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

!  Linearized phase functions
!   --extra terms for Phase function moment variations
!   -- must add TMS correction factor linearization

              IF ( DO_COLUMN_LINEARIZATION ) THEN
               DO K = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(K)) THEN
                  DO Q = 1, N_TOTALCOLUMN_WFS
                   VAR1 = EXACTSCAT_UP(V,K) * VAR_TMS(K,Q)
                   IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                    HELP = ZERO
                    DO L = 0, NMOMENTS_INPUT
                     IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(K) * DNL1
                      FACT  = ONE - SSFDEL(K)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,K,T)-FDNL1)/FACT
                      HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                      HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                      L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                     ELSE
                      L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                     ENDIF 
                     LEGPOLY = SS_PLEG_UP(V,1,L)
                     HELP = HELP + L_PHASMOM*LEGPOLY
                    ENDDO
                    L_EXACTSCAT_UP(V,K,Q) = HELP * TMS(K) + VAR1
                   ELSE
                    L_EXACTSCAT_UP(V,K,Q) = VAR1
                   ENDIF
                  ENDDO
                ENDIF
               ENDDO
              ENDIF

!  End upwelling clause

            ENDIF

!  Downwelling calculation
!  -----------------------

            IF ( DO_DNWELLING ) THEN

!  Multipliers, transmittances + linearizations

              call LCH_og_integration_dn_NF                          &
            ( maxlayers, max_atmoswfs,                               & ! Inputs
              nlayers, n_totalcolumn_wfs, LHmask_column_wfs,         & ! Inputs
              extinction, l_extinction, radii, local_lh_height_grid, & ! Inputs
              alpha_all, LH_alpha_all, attn, l_attn,                 & ! Inputs
              dn_multipliers(1,v), dn_lostrans(1,v),                 & ! outputs
              l_dn_multipliers(1,1,v), l_dn_lostrans(1,1,v) )          ! outputs

!  Debug
!              do n = 1, nlayers
!                write(89,'(i4,1p2e20.10)')n,l_dn_multipliers(n,1,1),L_dn_multipliers(n,2,1)
!              enddo
!              pause

!  Legendre polynomials

              COSSCAT = COSSCAT_DN(NLAYERS)
              SS_PLEG_DN(V,1,0) = ONE
              SS_PLEG_DN(V,1,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_DN(V,1,L) = DF1(L) * SS_PLEG_DN(V,1,L-1) * COSSCAT - DF2(L) * SS_PLEG_DN(V,1,L-2)
              ENDDO

!  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(N) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(N) * DNL1
                      FACT  = ONE - SSFDEL(N)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,N,T)-FDNL1) / FACT
                    ELSE
                      GK11 = PHASMOMS_TOTAL_INPUT(L,N,T)
                    ENDIF 
                    LEGPOLY = SS_PLEG_DN(V,1,L)
                    HELP = HELP + GK11 * LEGPOLY
                  ENDDO
                  EXACTSCAT_DN(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

!  Linearized phase functions
!   --extra terms for Phase function moment variations
!   -- must add TMS correction factor linearization

              IF ( DO_COLUMN_LINEARIZATION ) THEN
               DO K = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(K)) THEN
                  DO Q = 1, N_TOTALCOLUMN_WFS
                   VAR1 = EXACTSCAT_DN(V,K) * VAR_TMS(K,Q)
                   IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                    HELP = ZERO
                    DO L = 0, NMOMENTS_INPUT
                     IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(K) * DNL1
                      FACT  = ONE - SSFDEL(K)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,K,T)-FDNL1)/FACT
                      HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                      HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                      L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                     ELSE
                      L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                     ENDIF 
                     LEGPOLY = SS_PLEG_DN(V,1,L)
                     HELP = HELP + L_PHASMOM*LEGPOLY
                    ENDDO
                    L_EXACTSCAT_DN(V,K,Q) = HELP * TMS(K) + VAR1
                   ELSE
                    L_EXACTSCAT_DN(V,K,Q) = VAR1
                   ENDIF
                  ENDDO
                ENDIF 
               ENDDO
              ENDIF

!  End Downwelling clause

            ENDIF

!   Finish geometry loops

          ENDDO
        ENDDO
      ENDDO

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

      IF ( DO_UPWELLING ) THEN

!  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_UP(V,NC) = ZERO
        ENDDO

!  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!  Multiplier using new integration scheme

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_UP(V,N)
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS(N,V)
              SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE + UP_LOSTRANS(N,V) * SS_CUMSOURCE_UP(V,NC-1)
            ENDDO
          ENDDO

!  Offgrid output-------
!    Add additional partial layer source term = Exact Phase Func * Multiplier
!    Get final cumulative source and set the Single scatter results
!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

!          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!            UT = PARTLAYERS_OUTINDEX(UTA)
!            N  = PARTLAYERS_LAYERIDX(UT)
!            DO V = 1, N_GEOMETRIES
!              HELP           = EXACTSCAT_UP(V,N)
!              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
!              SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,NC)
!              TRANS          = UP_LOSTRANS_UT(UT,V)
!              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
!              SSCORRECTION   = SSFLUX * FINAL_SOURCE
!              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
!            ENDDO
!          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_UP(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
!          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

!  Recurrence relation for the DOWNWELLING intensity
!  =================================================

      IF ( DO_DNWELLING ) THEN

!  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_DN(V,NC) = ZERO
        ENDDO

!  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!      Multiplier by new integration method

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              HELP =  EXACTSCAT_DN(V,N)
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS(N,V)
              SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE + DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,NC-1)
            ENDDO
          ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multiplier
!    Set final cumulative source and Correct the intensity
!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

!          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!            UT = PARTLAYERS_OUTINDEX(UTA)
!            N  = PARTLAYERS_LAYERIDX(UT)
!            DO V = 1, N_GEOMETRIES
!              HELP = EXACTSCAT_DN(V,N)
!              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS_UT(UT,V)
!              SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,NC)
!              TRANS          = DN_LOSTRANS_UT(UT,V)
!              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
!              SSCORRECTION   = SSFLUX * FINAL_SOURCE
!              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
!            ENDDO
!          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_DN(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
!          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

!  end optical depth loop and Downwelling intensity clause

        ENDDO
      ENDIF

!  Recurrence relation for the UPWELLING Jacobians
!  ===============================================

      IF ( DO_UPWELLING .AND. DO_COLUMN_LINEARIZATION ) THEN

!  Start the main layer variation loop
!  -----------------------------------
   
!  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_SS_CUMSOURCE(Q,V) = ZERO
         ENDDO
        ENDDO

!  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!  ----------------------------------------

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT = NLEVEL + 1

!  Cumulative single scatter source terms to layer NUT
!  ---------------------------------------------------

!    1. Get layer source terms = Exact scattering * Multiplier
!    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                LSS = EXACTSCAT_UP(V,N)   * L_UP_MULTIPLIERS(N,Q,V)     &
                  + L_EXACTSCAT_UP(V,N,Q) *   UP_MULTIPLIERS(N,V)
                L_SS_CUMSOURCE(Q,V) =  LSS                              &
                  +  UP_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V)          &
                  + L_UP_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_UP(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

!  Offgrid output----------
!  Set final cumulative source and Single scatter Weighting function

!          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!            UT = PARTLAYERS_OUTINDEX(UTA)
!            N  = PARTLAYERS_LAYERIDX(UT)
!            DO V = 1, N_GEOMETRIES
!              DO Q = 1, N_TOTALCOLUMN_WFS
!                LSS = EXACTSCAT_UP(V,N)  * L_UP_MULTIPLIERS_UT(UT,Q,V) &
!               +   L_EXACTSCAT_UP(V,N,Q) *   UP_MULTIPLIERS_UT(UT,V)
!                L_FINAL_SOURCE =  LSS &
!               +   UP_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V) &
!               + L_UP_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_UP(V,NC)
!                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
!                 COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
!              ENDDO
!            ENDDO

!  Ongrid output---------
!  just set to the cumulative source term 

!          ELSE
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
            ENDDO
!          ENDIF

!  Check for updating the recursion

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end loop over output optical depths

        ENDDO

!  end Upwelling Jacobian clause

      ENDIF

!  Recurrence relation for the DOWNWELLING Jacobians
!  =================================================

      IF ( DO_DNWELLING .AND. DO_COLUMN_LINEARIZATION ) THEN

!  Start the main layer variation loop
!  -----------------------------------

!  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1,N_TOTALCOLUMN_WFS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

!  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

!  Cumulative single scatter source terms to layer NUT
!  ---------------------------------------------------

!    1. Get layer source terms = Exact scattering * Multiplier
!    2. Loop over layers working downwards to NUT

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                LSS = EXACTSCAT_DN(V,N)   * L_DN_MULTIPLIERS(N,Q,V)     &
                  + L_EXACTSCAT_DN(V,N,Q) *   DN_MULTIPLIERS(N,V)
                L_SS_CUMSOURCE(Q,V) =  LSS                              &
                  +  DN_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V)          &
                  + L_DN_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_DN(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

!  Offgrid output----------
!  Set final cumulative source and Single scatter Weighting function

!          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!            UT = PARTLAYERS_OUTINDEX(UTA)
!            N  = PARTLAYERS_LAYERIDX(UT)
!            DO V = 1, N_GEOMETRIES
!              DO Q = 1, N_TOTALCOLUMN_WFS
!                LSS = EXACTSCAT_DN(V,N)  * L_DN_MULTIPLIERS_UT(UT,Q,V) &
!                +  L_EXACTSCAT_DN(V,N,Q) *   DN_MULTIPLIERS_UT(UT,V)
!                L_FINAL_SOURCE =  LSS &
!                +  DN_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V) &
!                + L_DN_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_DN(V,NC)
!                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
!                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
!              ENDDO
!            ENDDO

!  Ongrid output---------
!  just set to the cumulative source term 

!          ELSE
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
            ENDDO
!          ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over output optical depths

        ENDDO

!  end Downwelling Jacobian clause

      ENDIF

!  DEBUG

!        do v = 1, 36
!          write(99,1234) INTENSITY_SS(1,v,1),COLUMNWF_SS(3,1,v,1)
!        enddo
!      pause
! 1234 format(1p3e20.10)

!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_SSCORR_OUTGOING_NFLH

!

subroutine LCH_og_integration_up_NF                      &
       ( maxlayers, maxvars,                             & ! Inputs
         nlayers, n_totalcolumn_wfs, LHmask_column_wfs,  & ! Inputs
         extinction, LC_extinction, radii, LH_radii,     & ! Inputs
         alpha_all, LH_alpha_all, attn, LC_attn,         & ! Inputs
         emultipliers, lostrans,                         & ! Inputs     
         LC_emultipliers, LC_lostrans )                    ! Outputs

      implicit none

!  Does the optical depth integration over layers.

!    No Partial layer integration
!    No fine layering here
!    Column weighting functions only.
!    LH-mask array --> which parameter has additional LH variation

!  inputs
!  ------

!  dimensioning

      integer, intent(in)  :: maxlayers
      integer, intent(in)  :: maxvars

!  control

      integer, intent(in)  :: nlayers

!  Column linearization control inputs

      logical, intent(in)  :: LHmask_column_wfs ( maxvars)
      integer, intent(in)  :: n_totalcolumn_wfs

!  radii

      real(kind=8), intent(in)  :: radii      (0:maxlayers)
      real(kind=8), intent(in)  :: LH_radii   (0:maxlayers)

!  line of sight angles

      real(kind=8), intent(in)  :: alpha_all     (0:maxlayers)
      real(kind=8), intent(in)  :: LH_alpha_all  (0:maxlayers)

!  Extinction inputs

      real(kind=8), intent(in)  :: extinction    (maxlayers)
      real(kind=8), intent(in)  :: LC_extinction (maxlayers, maxvars)

!  Attenuations

      real(kind=8), intent(in)  :: attn      ( 0:maxlayers )

!  Linearized attenuations

      real(kind=8), intent(in)  :: LC_attn ( 0:maxlayers, maxvars )

!  outputs
!  -------

!  Regular quantities

      real(kind=8), intent(out) :: emultipliers   (maxlayers)
      real(kind=8), intent(out) :: lostrans       (maxlayers)

!  Linearized quantities

      real(kind=8), intent(out) :: LC_emultipliers   ( maxlayers, maxvars )
      real(kind=8), intent(out) :: LC_lostrans       ( maxlayers, maxvars )

!  help variables
!  --------------

      integer      :: n, q
      real(kind=8) :: l_func_1, l_func_2, l_tran_1, l_csqt_1
      real(kind=8) :: esum, emult, l_emult, kn, l_kn
      real(kind=8) :: salpha, calpha, raycon
      real(kind=8) :: argm_1, csq_1, csqt_1, csq_2, cot_1, cot_2
      real(kind=8) :: func_1, func_2, tran_1, S

      real(kind=8) :: LH_raycon, LH_kn, LH_argm_1, LH_csq_1
      real(kind=8) :: LH_csq_2, LH_cot_1, LH_cot_2, LHS, l_sum

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(kind=8), parameter  :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)    = 0.0d0
        emultipliers(n) = 0.0d0
      enddo

!  Column weighting functions

      do n = 1, nlayers
        do q = 1, n_totalcolumn_wfs
          LC_emultipliers(n,q) = 0.0d0
          LC_lostrans(n,q)     = 0.0d0
        enddo
      enddo

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      calpha = dcos(alpha_all(0))
      raycon    = radii(0) * salpha
      LH_raycon = LH_radii(0) * salpha + radii(0) * calpha * LH_alpha_all(0)

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1  = 1.0d0 / salpha / salpha
      cot_1  = calpha / salpha

      LH_csq_1 = -2.0d0 * csq_1 * cot_1 * LH_alpha_all(n)
      LH_cot_1 = - csq_1 * LH_alpha_all(n)

!  Start layer loop
!  ================

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

        S   = alpha_all(n) - alpha_all(n-1)
        LHS =  LH_alpha_all(n) - LH_alpha_all(n-1)

        LH_csq_2 = -2.0d0 * csq_2 * cot_2 * LH_alpha_all(n-1)
        LH_cot_2 = - csq_2 * LH_alpha_all(n-1)

!  set up

        kn     =    raycon * extinction(n)
        LH_kn  = LH_raycon * extinction(n)

        argm_1 = cot_2 - cot_1
        tran_1 = dexp ( - kn * argm_1 )
        csqt_1 = csq_1 * tran_1

        LH_argm_1 = LH_cot_2 - LH_cot_1

!  Intensity multipliers + transmittance
!  -------------------------------------

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic scattering multiplier

        func_2 = attn(n-1) *  csq_2 
        func_1 = attn(n)   * csqt_1
        esum  = 0.5d0 * ( func_1 + func_2 )
        emult = esum * kn
        emultipliers(n) = emult * S

!  Column Linearizations
!  ---------------------

!  Elastic

        do q = 1, n_totalcolumn_wfs
          l_kn   = raycon * LC_extinction(n,q)
          if ( LHmask_column_wfs(q) ) then
            l_kn = l_kn + LH_kn
            l_tran_1 = - ( l_kn * argm_1 + kn * LH_argm_1 ) * tran_1
            l_csqt_1 = csq_1 * l_tran_1 + LH_csq_1 * tran_1 
            l_func_2 = LC_attn(n-1,q) * csq_2  + attn(n-1) * LH_csq_2
            l_func_1 = LC_attn(n,q)   * csqt_1 + attn(n)   * l_csqt_1
            l_sum    = 0.5d0 * ( l_func_1 + l_func_2 )
            l_emult  = l_kn * esum + kn * l_sum
            LC_lostrans(n,q)     = l_tran_1
            LC_emultipliers(n,q) = S * l_emult + LHS * emult
          else
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_csqt_1 = csq_1 * l_tran_1
            l_func_2 = LC_attn(n-1,q) * csq_2
            l_func_1 = LC_attn(n,q) *   csqt_1  + attn(n)   * l_csqt_1
            l_sum    = 0.5d0 * ( l_func_1 + l_func_2 )
            l_emult  = l_kn * esum + kn * l_sum
            LC_lostrans(n,q)     = l_tran_1
            LC_emultipliers(n,q) = S * l_emult
          endif
        enddo

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2
        LH_csq_1 = LH_csq_2
        LH_cot_1 = LH_cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
END subroutine LCH_og_integration_up_NF

!

subroutine LCH_og_integration_dn_NF                      &
       ( maxlayers, maxvars,                             & ! Inputs
         nlayers, n_totalcolumn_wfs, LHmask_column_wfs,  & ! Inputs
         extinction, LC_extinction, radii, LH_radii,     & ! Inputs
         alpha_all, LH_alpha_all, attn, LC_attn,         & ! Inputs
         emultipliers, lostrans,                         & ! Inputs     
         LC_emultipliers, LC_lostrans )                    ! Outputs

      implicit none

!  Does the optical depth integration over layers.

!    No Partial layer integration
!    No fine layering here
!    Column weighting functions only.
!    LH-mask array --> which parameter has additional LH variation

!  inputs
!  ------

!  dimensioning

      integer, intent(in)  :: maxlayers
      integer, intent(in)  :: maxvars

!  control

      integer, intent(in)  :: nlayers

!  Column linearization control inputs

      logical, intent(in)  :: LHmask_column_wfs ( maxvars)
      integer, intent(in)  :: n_totalcolumn_wfs

!  radii

      real(kind=8), intent(in)  :: radii      (0:maxlayers)
      real(kind=8), intent(in)  :: LH_radii   (0:maxlayers)

!  line of sight angles

      real(kind=8), intent(in)  :: alpha_all     (0:maxlayers)
      real(kind=8), intent(in)  :: LH_alpha_all  (0:maxlayers)

!  Extinction inputs

      real(kind=8), intent(in)  :: extinction    (maxlayers)
      real(kind=8), intent(in)  :: LC_extinction (maxlayers, maxvars)

!  Attenuations

      real(kind=8), intent(in)  :: attn      ( 0:maxlayers )

!  Linearized attenuations

      real(kind=8), intent(in)  :: LC_attn ( 0:maxlayers, maxvars )

!  outputs
!  -------

!  Regular quantities

      real(kind=8), intent(out) :: emultipliers   (maxlayers)
      real(kind=8), intent(out) :: lostrans       (maxlayers)

!  Linearized quantities

      real(kind=8), intent(out) :: LC_emultipliers   ( maxlayers, maxvars )
      real(kind=8), intent(out) :: LC_lostrans       ( maxlayers, maxvars )

!  help variables
!  --------------

      integer      :: n, q
      real(kind=8) :: l_func_1, l_func_2, l_tran_1, l_csqt_1
      real(kind=8) :: esum, emult, l_emult, kn, l_kn
      real(kind=8) :: salpha, calpha, raycon
      real(kind=8) :: argm_1, csq_1, csqt_1, csq_2, cot_1, cot_2
      real(kind=8) :: func_1, func_2, tran_1, S

      real(kind=8) :: LH_raycon, LH_kn, LH_argm_1, LH_csq_1
      real(kind=8) :: LH_csq_2, LH_cot_1, LH_cot_2, LHS, l_sum

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(kind=8), parameter  :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)    = 0.0d0
        emultipliers(n) = 0.0d0
      enddo

      do n = 1, nlayers
        do q = 1, n_totalcolumn_wfs
          LC_emultipliers(n,q) = 0.0d0
          LC_lostrans(n,q)     = 0.0d0
        enddo
      enddo

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      calpha = dcos(alpha_all(0))
      raycon  = radii(0) * salpha
      LH_raycon = LH_radii(0) * salpha + radii(0) * calpha * LH_alpha_all(0)

!  Work down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1  = 1.0d0 / salpha / salpha
      cot_1  = calpha / salpha

      LH_csq_1 = -2.0d0 * csq_1 * cot_1 * LH_alpha_all(n)
      LH_cot_1 = - csq_1 * LH_alpha_all(n)

!  Start layer loop
!  ================

      do n = 1, nlayers

!  Save some quantities

        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

        S    = alpha_all(n) - alpha_all(n-1)
        LHS  =  LH_alpha_all(n) - LH_alpha_all(n-1)

        LH_csq_2 = -2.0d0 * csq_2 * cot_2 * LH_alpha_all(n)
        LH_cot_2 = - csq_2 * LH_alpha_all(n)

!  set up

        kn     =    raycon * extinction(n)
        LH_kn  = LH_raycon * extinction(n)
        argm_1 = cot_1 - cot_2
        tran_1 = dexp ( - kn * argm_1 )
        csqt_1 = csq_1 * tran_1
        LH_argm_1 = LH_cot_1 - LH_cot_2

!  Intensity multipliers + transmittance
!  -------------------------------------

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic scattering multiplier

        func_1 = attn(n-1) * csqt_1
        func_2 = attn(n)   * csq_2
        esum   = 0.5d0 * ( func_1 + func_2 )
        emult  = esum * kn
        emultipliers(n) = emult * S

!  Column Linearizations
!  ---------------------

!  Elastic

        do q = 1, n_totalcolumn_wfs
          l_kn   = raycon * LC_extinction(n,q)
          if ( LHmask_column_wfs(q) ) then
            l_kn = l_kn + LH_kn
            l_tran_1 = - ( l_kn * argm_1 + kn * LH_argm_1 ) * tran_1
            l_csqt_1 = csq_1 * l_tran_1 + LH_csq_1 * tran_1 
            l_func_2 = LC_attn(n,q)   * csq_2   + attn(n)   * LH_csq_2
            l_func_1 = LC_attn(n-1,q) * csqt_1  + attn(n-1) * l_csqt_1
            l_sum    = 0.5d0 * ( l_func_1 + l_func_2 )
            l_emult  = l_kn * esum + kn * l_sum
            LC_lostrans(n,q)     = l_tran_1
            LC_emultipliers(n,q) = S * l_emult +  LHS * emult
          else
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_csqt_1 = csq_1 * l_tran_1
            l_func_2 = LC_attn(n,q) * csq_2
            l_func_1 = LC_attn(n-1,q) *   csqt_1 + attn(n-1)   * l_csqt_1
            l_sum    = 0.5d0 * ( l_func_1 + l_func_2 )
            l_emult  = l_kn * esum + kn * l_sum
            LC_lostrans(n,q)     = l_tran_1
            LC_emultipliers(n,q) = S * l_emult
          endif
        enddo

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2
        LH_csq_1 = LH_csq_2
        LH_cot_1 = LH_cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
END subroutine LCH_og_integration_dn_NF

!
             
subroutine LCH_og_attenuations_NF              &
       ( maxlayers, maxvars, nlayers,          & ! Inputs
         n_totalcolumn_wfs, LHmask_column_wfs, & ! Inputs
         extinction, LC_extinction,            & ! Inputs
         sunpaths, LH_sunpaths, ntraverse,     & ! Inputs
         attn, LC_attn )                         ! Outputs

      implicit none

!  Does attenuations
!    LH-mask array --> which parameter has additional LH variation

!  inputs
!  ------

!  dimensioning

      integer, intent(in)  :: maxlayers
      integer, intent(in)  :: maxvars

!  control

      integer, intent(in)  :: nlayers

!  Column linearization control inputs

      logical, intent(in)  :: LHmask_column_wfs ( maxvars)
      integer, intent(in)  :: n_totalcolumn_wfs

!  Whole layers

      integer     , intent(in)  :: ntraverse(0:maxlayers)
      real(kind=8), intent(in)  :: sunpaths(0:maxlayers,maxlayers)
      real(kind=8), intent(in)  :: LH_sunpaths(0:maxlayers,maxlayers)

!  Extinction inputs 

      real(kind=8), intent(in)  :: extinction   (maxlayers)
      real(kind=8), intent(in)  :: LC_extinction (maxlayers,maxvars )

!  outputs
!  -------

      real(kind=8), intent(out) :: attn    ( 0:maxlayers )
      real(kind=8), intent(out) :: LC_attn ( 0:maxlayers,maxvars )

!  help variables
!  --------------

      integer      :: n, k, q
      real(kind=8) :: tau, l_iop, sum
      
!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(kind=8), parameter  :: LOCAL_CUTOFF = 32.0D0

!  Zero output

      do n = 0, nlayers
        attn(n) = 0.0d0
        do q = 1, n_totalcolumn_wfs
          LC_attn(n,q) = 0.0d0
        enddo
      enddo

!  attenuation functions, whole layers

      do n = 0, nlayers
        tau     = 0.0d0
        do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo

!  Linearized attenuations
!  -----------------------

!  Linearized attenuation factors. Column linearization

      do q = 1, n_totalcolumn_wfs
        if ( LHmask_column_wfs(q) ) then
          do n = 0, nlayers
            sum = 0.0d0
            do k = 1, ntraverse(n)
              l_iop = LC_extinction(k,q) * sunpaths(n,k) + &
                         extinction(k)   * LH_sunpaths(n,k)
              sum = sum + l_iop
            enddo
            LC_attn(n,q) = - sum * attn(n) 
          enddo
        else
          do n = 0, nlayers
            sum = 0.0d0
            do k = 1, ntraverse(n)
              l_iop = LC_extinction(k,q) * sunpaths(n,k)
              sum = sum + l_iop
            enddo
            LC_attn(n,q) = - sum * attn(n) 
          enddo
        endif
      enddo

!  Finish

      return
end subroutine LCH_og_attenuations_NF

!

SUBROUTINE CHAPMAN_FUNCTION_NF                      &
          ( MAXLAYERS, DO_PLANE_PARALLEL, NLAYERS,  & ! Inputs
            COS_SZA, EARTH_RADIUS, HEIGHT_GRID,     & ! Inputs
            CHAPMAN_FACTORS )                         ! Output

!  This is the Chapman function calculation of the slant path
!  optical thickness array TAUTHICK_INPUT.

!  This module calculates TAUTHICK_INPUT internally inside LRRS
!  saving you the job of doing it yourself.

!  You must specify the Earth_radius and the height grid in order
!  to make this work.

!  This is a straightforward geometrical calculation, and is only
!  valid for a NON-REFRACTIVE atmosphere.

!  The non-refractive condition will be checked before this module
!  is called

      IMPLICIT NONE

!  Input arguemnts
!  ---------------

      LOGICAL     , intent(in)  :: DO_PLANE_PARALLEL
      INTEGER     , intent(in)  :: MAXLAYERS, NLAYERS
      REAL(kind=8), intent(in)  :: COS_SZA
      REAL(kind=8), intent(in)  :: EARTH_RADIUS
      REAL(kind=8), intent(in)  :: HEIGHT_GRID(0:MAXLAYERS)

!  output arguemnts
!  ----------------

      REAL(kind=8), intent(out) :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: N, M
      REAL(kind=8)  :: GM_TOA, HELP1, HELP2
      REAL(kind=8)  :: H(0:MAXLAYERS), DELZ(MAXLAYERS) 
      REAL(kind=8)  :: STH, CTH, DELS, S1, S0

!  get spherical optical depths
!  ----------------------------

!  Prepare spherical attenuation (shell geometry)

      IF ( .NOT.DO_PLANE_PARALLEL ) THEN

        GM_TOA = DSQRT ( 1.0D0 - COS_SZA * COS_SZA )
        DO N = 0, NLAYERS
          H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
        ENDDO
        DO N = 1, NLAYERS
          DELZ(N) = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        ENDDO

        DO N = 1, NLAYERS
          STH = GM_TOA * H(N)/H(0)
          CTH = DSQRT ( 1.0d0 - STH * STH )
          S0 = 0.0d0
          HELP1 = H(0)*CTH
          HELP2 = -H(0)*H(0)*STH*STH
          DO M = 1, N
            S1 = HELP1 - DSQRT(HELP2 + H(M)*H(M))
            DELS = S1 - S0
            CHAPMAN_FACTORS(N,M) = DELS / DELZ(M)
            S0 = S1
          ENDDO
        ENDDO

!  Plane parallel

      ELSE

        DO N = 1, NLAYERS
          DO M = 1, N
            CHAPMAN_FACTORS(N,M) = 1.0d0 / COS_SZA
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE CHAPMAN_FUNCTION_NF

!

SUBROUTINE CHAPMAN_FUNCTION_NF_PLUS                 &
          ( MAXLAYERS, MAXPARS, DO_PLANE_PARALLEL,  & ! Inputs
            NLAYERS, NPARS, COS_SZA, EARTH_RADIUS,  & ! Inputs
            HEIGHT_GRID, L_HEIGHT_GRID,             & ! Inputs
            CHAPMAN_FACTORS, L_CHAPMAN_FACTORS )      ! Outputs

!  This is the Chapman function calculation of the slant path
!  optical thickness array TAUTHICK_INPUT.

!  This module calculates TAUTHICK_INPUT internally inside LRRS
!  saving you the job of doing it yourself.

!  You must specify the Earth_radius and the height grid in order
!  to make this work.

!  This is a straightforward geometrical calculation, and is only
!  valid for a NON-REFRACTIVE atmosphere.

!  The non-refractive condition will be checked before this module
!  is called

      IMPLICIT NONE

!  Input arguments
!  ---------------

      LOGICAL     , intent(in)  :: DO_PLANE_PARALLEL
      INTEGER     , intent(in)  :: MAXLAYERS, MAXPARS, NLAYERS, NPARS
      REAL(kind=8), intent(in)  :: COS_SZA
      REAL(kind=8), intent(in)  :: EARTH_RADIUS
      REAL(kind=8), intent(in)  :: HEIGHT_GRID  (0:MAXLAYERS)
      REAL(kind=8), intent(in)  :: L_HEIGHT_GRID(MAXPARS,0:MAXLAYERS)

!  output arguemnts
!  ----------------

      REAL(kind=8), intent(out) :: CHAPMAN_FACTORS  (MAXLAYERS,MAXLAYERS)
      REAL(kind=8), intent(out) :: L_CHAPMAN_FACTORS(MAXPARS,MAXLAYERS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: N, M, Q
      REAL(kind=8)  :: GM_TOA, HELP1, HELP2, HELP3
      REAL(kind=8)  :: H(0:MAXLAYERS), DELZ(MAXLAYERS) 
      REAL(kind=8)  :: STH, CTH, DELS, S1, S0

      REAL(kind=8)  :: L_HELP1(MAXPARS), L_HELP2(MAXPARS), L_HELP3
      REAL(kind=8)  :: L_H   (MAXPARS,0:MAXLAYERS)
      REAL(kind=8)  :: L_DELZ(MAXPARS,MAXLAYERS) 
      REAL(kind=8)  :: L_STH(MAXPARS), L_CTH(MAXPARS), L_DELS
      REAL(kind=8)  :: L_S1(MAXPARS), L_S0(MAXPARS)

!  Initialize

      DO N = 1, NLAYERS
        DO M = 1, N
          CHAPMAN_FACTORS(N,M) = 0.0d0
          DO Q = 1, NPARS
            L_CHAPMAN_FACTORS(Q,N,M) = 0.0d0
          ENDDO
        ENDDO
      ENDDO

!  Prepare spherical attenuation (shell geometry)

      IF ( .NOT.DO_PLANE_PARALLEL ) THEN

        GM_TOA = DSQRT ( 1.0D0 - COS_SZA * COS_SZA )
        DO N = 0, NLAYERS
          H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
          DO Q = 1, NPARS
            L_H(Q,N) = L_HEIGHT_GRID(Q,N)
          ENDDO
        ENDDO
        DO N = 1, NLAYERS
          DELZ(N) = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
          DO Q = 1, NPARS
            L_DELZ(Q,N) = L_HEIGHT_GRID(Q,N-1) - L_HEIGHT_GRID(Q,N)
          ENDDO
        ENDDO

        DO N = 1, NLAYERS
          STH = GM_TOA * H(N)/H(0)
          CTH = DSQRT ( 1.0d0 - STH * STH )
          DO Q = 1, NPARS
             L_STH(Q) = STH * ( (L_H(Q,N)/H(N)) - (L_H(Q,0)/H(0)) )
             L_CTH(Q) = - STH * L_STH(Q) / CTH
          ENDDO
          S0 = 0.0d0
          HELP1 = H(0)*CTH
          HELP2 = -H(0)*H(0)*STH*STH
          DO Q = 1, NPARS
            L_S0(Q) = 0.0d0
            L_HELP1(Q) = L_H(Q,0)*CTH + H(0) * L_CTH(Q)
            L_HELP2(Q) = 2.0d0*HELP2*((L_H(Q,0)/H(0))+(L_STH(Q)/STH))
          ENDDO
          DO M = 1, N
            HELP3 = DSQRT(HELP2 + H(M)*H(M))
            S1 = HELP1 - HELP3
            DELS = S1 - S0
            CHAPMAN_FACTORS(N,M) = DELS / DELZ(M)
            DO Q = 1, NPARS
              L_HELP3 = 0.5d0*(L_HELP2(Q)+2.0d0*H(M)*L_H(Q,M))/ HELP3
              L_S1(Q) = L_HELP1(Q) - L_HELP3
              L_DELS  = L_S1(Q) - L_S0(Q)
              L_CHAPMAN_FACTORS(Q,N,M) = CHAPMAN_FACTORS(N,M) * &
                          ( (L_DELS/DELS) - (L_DELZ(Q,M)/DELZ(M)) )
            ENDDO
            S0 = S1
            DO Q = 1, NPARS
              L_S0(Q) = L_S1(Q)
!              L_CHAPMAN_FACTORS(Q,N,M) = 0.0d0         ! TEMPO
            ENDDO
          ENDDO
        ENDDO

!  Plane parallel

      ELSE

        DO N = 1, NLAYERS
          DO M = 1, N
            CHAPMAN_FACTORS(N,M) = 1.0d0 / COS_SZA
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE CHAPMAN_FUNCTION_NF_PLUS

!

subroutine LH_sphergeom_NF                               &
       ( maxlayers, nlayers, eradius, v,                 & ! inputs
         heights, L_hgt, alpha_boa, theta_boa, phi_boa,  & ! outputs
         ntraverse,  cosscat_up, cosscat_dn, radii,      & ! outputs
         sunpaths, alpha_all,  L_sunpaths,  L_alpha_all, & ! outputs
         fail, message )                                   ! outputs

!  Completely stand-alone geometry routine for the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.1. September 2007, Partial layer geometries added
!  Version 2.2. July 2009, This version without partials

      implicit none

!  inputs

      integer     , intent(in)    :: maxlayers, v
      integer     , intent(in)    :: nlayers
      real(kind=8), intent(in)    :: heights (0:maxlayers), L_hgt (0:maxlayers)
      real(kind=8), intent(in)    :: eradius, alpha_boa, theta_boa
      real(kind=8), intent(inout) :: phi_boa

!  main outputs (geometry)

      integer     , intent(out) :: ntraverse   (0:maxlayers)
      real(kind=8), intent(out) :: sunpaths   (0:maxlayers,maxlayers)
      real(kind=8), intent(out) :: alpha_all  (0:maxlayers)
      real(kind=8), intent(out) :: L_sunpaths (0:maxlayers,maxlayers)
      real(kind=8), intent(out) :: L_alpha_all(0:maxlayers)
      real(kind=8), intent(out) :: cosscat_up (0:maxlayers)
      real(kind=8), intent(out) :: cosscat_dn (0:maxlayers)
      real(kind=8), intent(out) :: radii      (0:maxlayers)

!  Status output

      logical      , intent(out) :: fail
      character*(*), intent(out) :: message

!  Other (incidental) geometrical output

      real(kind=8) :: lospaths(maxlayers)
      real(kind=8) :: theta_all  (0:maxlayers)
      real(kind=8) :: phi_all    (0:maxlayers)

      real(kind=8) :: L_theta_all  (0:maxlayers)
      real(kind=8) :: RR   (maxlayers)
      real(kind=8) :: L_RR (maxlayers)

!  Local

      logical          direct_sun
      integer          n, k, krad, n1
      real(kind=8) :: deg_to_rad, ex, ey, ez, px, pz, b
      real(kind=8) :: salpha_boa, calpha_boa, sphi_boa
      real(kind=8) :: stheta_boa, ctheta_boa, cphi_boa
      real(kind=8) :: ksi, cksi, sksi, tangr, fac, fac1
      real(kind=8) :: ctheta, stheta, calpha, salpha, cphi
      real(kind=8) :: sth0, th0, sks1, ks1, sth1, th1
      real(kind=8) :: xicum, xicum0, cxicum, sxicum

      real(kind=8) :: L_px, L_pz, L_b
      real(kind=8) :: L_ksi, L_cksi, L_sksi, L_tangr
      real(kind=8) :: L_ctheta, L_stheta, L_calpha, L_salpha
      real(kind=8) :: L_sth0, L_th0, L_ks1
      real(kind=8) :: L_xicum, L_cxicum, L_sxicum, L_sth1, L_th1

!  Initialise output

      fail = .false.
      message = ' '

!  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif

!  zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        alpha_all(n)   = 0.0d0
        L_alpha_all(n) = 0.0d0
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k)   = 0.0d0
         L_sunpaths(n,k) = 0.0d0
        enddo
      enddo

!  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))
      cosscat_up (nlayers) = - calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa 
      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa 

!  Radii
!  -----

!  layer levels

      do n = 0, nlayers
        radii(n)   = eradius + heights(n)
      enddo

!  Radii ratios

      do n = 1, nlayers
        RR(n)     = radii(n) / radii(n-1)
        L_RR(n)   = (L_hgt(n)-RR(n)*L_hgt(n-1))/radii(n-1)
      enddo

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up(n-1) = cosscat_up(n)
          cosscat_dn(n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)

!  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k)   = radii(k-1)-radii(k)
              L_sunpaths(n,k) = L_hgt(k-1)-L_hgt(k)
            enddo
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
!  ----start
            sth0 = stheta_boa
            th0  = theta_all(n)
            L_sth0 = 0.0d0
            L_th0  = 0.0d0
!  Work upwards------------------------------
            do k = n, 1, -1
!  ----sunpaths
              sth1 = sth0 * RR(k)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sks1 = dsin(ks1)
              sunpaths(n,k) = sks1*radii(k)/sth1
!  ----Lsunpaths
              L_sth1 = L_sth0 * RR(k) + sth0 * L_RR(k)
              L_th1  = L_sth1 / dcos(th1)
              L_ks1  = L_th0 - L_th1
              L_sunpaths(n,k) = dcos(ks1)* L_ks1*radii(k)/sth1 + &
                                sks1     * L_hgt(k)  /sth1 -     &
                                L_sth1 * sunpaths(n,k)   /sth1
!  ----replace
              sth0 = sth1
              th0  = th1
              L_sth0 = L_sth0
              L_th0  = L_th1
            enddo
          endif

!  End main layer loop

        enddo

!  Return, as everything now done

        return

!  end regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing sphericity geometry
!  ============================

!  define Unit solar vector at BOA

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k)   = radii(k-1)-radii(k)
          L_sunpaths(nlayers,k) = L_hgt(k-1)-L_hgt(k)
        enddo
      else
!  ----start
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        L_sth0 = 0.0d0
        L_th0  = 0.0d0
!  Work upwards------------------------------
        do k = nlayers, 1, -1
!  ----sunpaths
          sth1 = sth0 * RR(k)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sks1 = dsin(ks1)
          sunpaths(nlayers,k) = sks1*radii(k)/sth1
!  ----Lsunpaths
          L_sth1 = L_sth0 * RR(k) + sth0 * L_RR(k)
          L_th1  = L_sth1 / dcos(th1)
          L_ks1  = L_th0 - L_th1
          L_sunpaths(nlayers,k) = dcos(ks1)* L_ks1*radii(k)/sth1 + &
                               sks1     * L_hgt(k)  /sth1 -        &
                               L_sth1 * sunpaths(nlayers,k)/sth1
!  ----replace
          sth0 = sth1
          th0  = th1
          L_sth0 = L_sth1
          L_th0  = L_th1

!          write(55,'(3i4,1p2e20.11)')v,nlayers, &
!              k,sunpaths(nlayers,k),L_sunpaths(nlayers,k) 

        enddo
      endif

!  initialise los cumulative angle

      xicum    = 0.0d0
      L_xicum  = 0.0d0

!  set TOA direct illumination flag

      direct_sun = .true.

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1
  
!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

        L_salpha = (salpha_boa*L_hgt(nlayers)-salpha*L_hgt(n))/radii(n)
        L_alpha_all(n)  = L_salpha / calpha
        L_calpha        = - salpha * L_alpha_all(n)

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

        L_ksi = L_alpha_all(n1) - L_alpha_all(n)
        L_sksi =  cksi * L_ksi
        L_cksi = -sksi * L_ksi
        L_xicum =  L_xicum + L_ksi

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         L_theta_all(n) = L_xicum
         L_ctheta = - stheta * L_theta_all(n)
         L_stheta =   ctheta * L_theta_all(n)
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         sxicum = dsin(xicum)
         cxicum = dcos(xicum)
         px = - radii(n) * sxicum
         pz =   radii(n) * cxicum
         b   = ex*px + ez*pz

         L_sxicum =   cxicum * L_xicum
         L_cxicum = - sxicum * L_xicum
         L_px = - L_hgt(n) * sxicum - radii(n) * L_sxicum
         L_pz =   L_hgt(n) * cxicum + radii(n) * L_cxicum
         L_b = ex * L_px + ez * L_pz

         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)


         L_ctheta = - ( L_b + ctheta * L_hgt(n) ) / radii(n)
         L_stheta =  - ctheta * L_ctheta / stheta
         L_theta_all(n) = - L_ctheta / stheta

        endif

!  Fix phi by using constancy of scatter angle
!  Only for the scattering up directions..................

        cosscat_up(n) = cosscat_up(n+1)
        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)
         phi_all(n)     = dacos(cphi)
        endif

!        if ( v.eq.2)write(*,*)v,n,theta_all(n),L_theta_all(n)
!        if ( v.eq.2)write(*,*)v,n,alpha_all(n),L_alpha_all(n)

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!   Work up from level n to TOA
!    Layer top calculation gets left out at TOA

        if ( direct_sun ) then
         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          L_sth0 = L_stheta
          L_th0  = L_theta_all(n)
          do k = n, 1, -1

           sth1 = sth0 * RR(k)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sks1 = dsin(ks1)
           sunpaths(n,k) = sks1*radii(k)/sth1

           L_sth1 = L_sth0 * RR(k) + sth0 * L_RR(k)
           L_th1  = L_sth1 / dcos(th1)
           L_ks1  = L_th0 - L_th1
           L_sunpaths(n,k) = dcos(ks1)* L_ks1 * radii(k) /sth1 + &
                                sks1  * L_hgt(k)  /sth1 -        &
                               L_sth1 * sunpaths(n,k) /sth1

!          write(57,'(3i4,1p2e20.11)')v,n, &
!              k,sunpaths(n,k),L_sunpaths(n,k) 

           sth0 = sth1
           th0  = th1
           L_sth0 = L_sth1
           L_th0  = L_th1
          enddo
         endif
        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr   = stheta*radii(n)
         L_tangr = L_stheta*radii(n) + stheta*L_hgt(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)
         L_sth0 = ( L_tangr - sth0 * L_hgt(0) ) / radii(0)
         L_th0 = L_sth0 / dcos(th0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1

          L_sth1 = L_sth0 * RR(k) + sth0 * L_RR(k)
          L_th1  = L_sth1 / dcos(th1)
          L_ks1  = L_th0 - L_th1
          L_sunpaths(n,k) = dcos(ks1)* L_ks1*radii(k) /sth1 + &
                               sks1  * L_hgt(k)  /sth1 -      &
                              L_sth1 * sunpaths(n,k) /sth1

          sth0 = sth1
          th0  = th1
          L_sth0 = L_sth1
          L_th0  = L_th1
         enddo
        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr   = stheta*radii(n)
         L_tangr = L_stheta*radii(n) + stheta*L_hgt(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)
         L_sth0 = ( L_tangr - sth0 * L_hgt(0) ) / radii(0)
         L_th0 = L_sth0 / dcos(th0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = sth0 / RR(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          fac1 = fac / sth1
          sunpaths(n,k) = fac1 * dsin(ks1)*radii(k-1)

          L_sth1 = ( L_sth0 - sth1 * L_RR(k) ) / RR(k)
          L_th1  = L_sth1 / dcos(th1)
          L_ks1  = L_th0 - L_th1
          L_sunpaths(n,k) = dcos(ks1)* L_ks1 * radii(k-1) /fac1 + &
                                      sks1  * L_hgt(k-1) /fac1 -  &
                                  L_sth1 * sunpaths(n,k) /fac1

          sth0 = sth1
          th0  = th1
          L_sth0 = L_sth1
          L_th0  = L_th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)   =   2.0d0 * radii(krad) * dcos(th1)   
         L_sunpaths(n,krad+1) = - 2.0d0 * radii(krad) * sth1 * L_th1  &
                                + 2.0d0 * L_hgt(krad) * dcos(th1)

!  Complete tangent point calculation

        endif

!  End layer loop

      enddo

!  Finish

      return
end subroutine LH_sphergeom_NF

!

subroutine sphergeom_NF                               &
       ( maxlayers, nlayers, eradius,                 & ! inputs
         heights, alpha_boa, theta_boa, phi_boa,      & ! inputs
         ntraverse,  cosscat_up, cosscat_dn, radii,   & ! outputs
         sunpaths, alpha_all,                         & ! outputs
         fail, message )                                ! outputs

!  Completely stand-alone geometry routine for the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.1. September 2007, Partial layer geometries added
!  Version 2.2. July 2009, This version without partials

      implicit none

!  inputs

      integer     , intent(in)    :: maxlayers
      integer     , intent(in)    :: nlayers
      real(kind=8), intent(in)    :: eradius, heights (0:maxlayers)
      real(kind=8), intent(in)    :: alpha_boa, theta_boa
      real(kind=8), intent(inout) :: phi_boa

!  main outputs (geometry)

      integer     , intent(out) :: ntraverse   (0:maxlayers)
      real(kind=8), intent(out) :: sunpaths   (0:maxlayers,maxlayers)
      real(kind=8), intent(out) :: alpha_all  (0:maxlayers)
      real(kind=8), intent(out) :: cosscat_up (0:maxlayers)
      real(kind=8), intent(out) :: cosscat_dn (0:maxlayers)
      real(kind=8), intent(out) :: radii      (0:maxlayers)

!  Status output

      logical      , intent(out) :: fail
      character*(*), intent(out) :: message

!  Other (incidental) geometrical output

      real(kind=8) :: lospaths(maxlayers)
      real(kind=8) :: theta_all  (0:maxlayers)
      real(kind=8) :: phi_all    (0:maxlayers)

!  Local

      logical      :: direct_sun
      integer      :: n, k, krad, n1
      real(kind=8) :: deg_to_rad, ex, ey, ez, px, py, pz
      real(kind=8) :: salpha_boa, calpha_boa, sphi_boa
      real(kind=8) :: stheta_boa, ctheta_boa, cphi_boa
      real(kind=8) :: ksi, cksi, sksi, xicum, xicum0, tangr, fac
      real(kind=8) :: ctheta, stheta, calpha, salpha, cphi
      real(kind=8) :: b, sth0, th0, ks1, sth1, th1

!  Initialise output

      fail = .false.
      message = ' '

!  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif

!  zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = 0.0d0
        enddo
      enddo


!  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))
      cosscat_up (nlayers) = - calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa 
      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa 

!  Radii
!  -----

!  layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up(n-1) = cosscat_up(n)
          cosscat_dn(n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)

!  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
          endif

!  End main layer loop

        enddo

!  Return, as everything now done

        return

!  end regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing spehricity geometry
!  ============================

!  define Unit solar vector at BOA

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  initialise los cumulative angle

      xicum  = 0.0d0

!  set TOA direct illumination flag

      direct_sun = .true.

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1
  
!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
        endif

!  Fix phi by using constancy of scatter angle
!  Only for the scattering up directions..................

        cosscat_up(n) = cosscat_up(n+1)
        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)
         phi_all(n)     = dacos(cphi)
        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!   Work up from level n to TOA
!    Layer top calculation gets left out at TOA

        if ( direct_sun ) then
         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif
        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)   

!  Complete tangent point calculation

        endif

!  End layer loop

      enddo

!  Finish

      return
end subroutine sphergeom_NF

!

SUBROUTINE LIDORT_LC_PREPTRANS_NFLH                         &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS, & ! Inputs
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,           & ! Inputs
          DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,          & ! Inputs
          CHAPMAN_FACTORS, L_CHAPMAN_FACTORS,               & ! Inputs
          AVERAGE_SECANT, LAYER_PIS_CUTOFF,                 & ! Inputs
          T_DELT_MUBAR,   T_UTDN_MUBAR,                     & ! Inputs
          LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                & ! Outputs
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT )               ! Outputs
  
!  Profile linearization of transmittances

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Inputs
!  ------

!  Plane-parallel flag

      LOGICAL, intent(in)  :: DO_PLANE_PARALLEL

!  Control integers

      INTEGER, intent(in)  :: NLAYERS, NBEAMS
      INTEGER, intent(in)  :: N_PARTLAYERS

!  output optical depth masks and indices

      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX (MAX_USER_LEVELS)

!  Linearization control

      INTEGER, intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical depths after delta-M scaling and Chapman function

      REAL(kind=8), intent(in)  :: DELTAU_VERT    ( MAXLAYERS )
      REAL(kind=8), intent(in)  :: PARTAU_VERT    ( MAX_PARTLAYERS )

!  Input Chapman Factors

      REAL(kind=8), intent(in)  :: CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, NBEAMS )

!  Last layer to include Particular integral solution

      INTEGER     , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant for solar beams.

      REAL(kind=8), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(kind=8), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Optical depths

      REAL(kind=8), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized Input Chapman Factors

      REAL(kind=8), intent(in)  :: L_CHAPMAN_FACTORS ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Outputs
!  -------

!  Linearized transmittances, solar beam

      REAL(kind=8), intent(out) :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(out) :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=8), intent(out) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(out) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER       :: N, Q, UT, K, IB
      REAL(kind=8)  :: WDEL, VAR, RHO, FAC, DELT, LAMDA, SUM, SUM2

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        N = 1
        DO Q = 1, LAYER_VARY_NUMBER(N)
          LC_INITIAL_TRANS(N,IB,Q) = ZERO
        ENDDO
        DO N = 2, NLAYERS
         IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           SUM = ZERO  
           DO K = 1, N-1
            FAC =   CHAPMAN_FACTORS(N-1,K,IB) * L_DELTAU_VERT(Q,K) + &
                  L_CHAPMAN_FACTORS(Q,N-1,K,IB)
            SUM = SUM + FAC * DELTAU_VERT(K) 
           ENDDO
           LC_INITIAL_TRANS(N,IB,Q) = - SUM
          ENDDO
         ELSE
          DO Q = 1, LAYER_VARY_NUMBER(N)
           LC_INITIAL_TRANS(N,IB,Q) = ZERO
          ENDDO
         ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_AVERAGE_SECANT(N,IB,Q) = L_CHAPMAN_FACTORS(Q,N,N,IB)
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DELT  = DELTAU_VERT(N)
              LAMDA = AVERAGE_SECANT(N,IB)
              FAC   =  - LAMDA
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
              ENDDO
              DO Q = 1, LAYER_VARY_NUMBER(K)
                SUM = ZERO
                DO K = 1, N
                  FAC = DELTAU_VERT(K) * ( L_CHAPMAN_FACTORS(Q,N,K,IB) &
                       + CHAPMAN_FACTORS(N,K,IB) * L_DELTAU_VERT(Q,K))
                  SUM = SUM + FAC
                ENDDO
                SUM2 = ZERO
                DO K = 1, N-1
                  FAC = DELTAU_VERT(K) * ( L_CHAPMAN_FACTORS(Q,N-1,K,IB) &
                       + CHAPMAN_FACTORS(N-1,K,IB) * L_DELTAU_VERT(Q,K))
                  SUM2 = SUM2 + FAC
                ENDDO
                LC_AVERAGE_SECANT(N,IB,Q) = LC_AVERAGE_SECANT(N,IB,Q) + ( SUM - SUM2 ) / DELT
              ENDDO
            ELSE
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!      DO N = 1, NLAYERS
!       DO Q = 1, LAYER_VARY_NUMBER(N)
!         write(*,*)n,q,LC_AVERAGE_SECANT(N,1,Q),AVERAGE_SECANT(N,1)
!       enddo
!      enddo

!  Linearization of Whole layer Transmittance factors
!  ==================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         LAMDA = AVERAGE_SECANT(N,IB)
         FAC   = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_DELT_MUBAR(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS

!  zero it

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_T_UTDN_MUBAR(UT,IB,Q) = ZERO
          ENDDO
        ENDDO

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_UTDN_MUBAR(UT,IB,Q) =  L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_PREPTRANS_NFLH
