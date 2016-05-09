!$Id: lidort_brdf_ls_supplement.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED code             (3.5)         #
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
! #            LIDORT_BRDF_INPUTS_PLUS                          #
! #                                                             #
! #            LIDORT_BRDF_LS_MASTER (master), calling          #
! #                                                             #
! #              LIDORT_BRDF_MAKER_PLUS                         #
! #              LIDORT_BRDF_LS_FOURIER                         #
! #                                                             #
! #                                                             #
! ###############################################################

SUBROUTINE LIDORT_BRDF_LS_MASTER                                           &
        ( DO_USER_STREAMS, DO_SHADOW_EFFECT, DO_GLITTER_DBMS,              & ! Inputs
          DO_SURFACE_EMISSION, N_BRDF_KERNELS, WHICH_BRDF,                 & ! Inputs
          LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS,             & ! Inputs
          N_BRDF_PARAMETERS, BRDF_PARAMETERS,                              & ! Inputs
          DO_KERNEL_PARAMS_WFS, DO_KERNEL_FACTOR_WFS, DO_KPARAMS_DERIVS,   & ! Inputs
          N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS, N_SURFACE_WFS,         & ! Inputs
          NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,                & ! Inputs
          BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,                      & ! Inputs
          BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,                    & ! Outputs
          LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,        & ! Outputs
          EXACTDB_BRDFUNC,  LS_EXACTDB_BRDFUNC,                            & ! Outputs
          EMISSIVITY, USER_EMISSIVITY, LS_EMISSIVITY, LS_USER_EMISSIVITY )   ! Outputs

!  Prepares (linearizations of) the bidirectional reflectance functions
!  necessary for LIDORT.

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input arguments
!  ===============

!  User stream Control

      LOGICAL, intent(in) :: DO_USER_STREAMS

!  Surface emission

      LOGICAL, intent(in) :: DO_SURFACE_EMISSION

!   Number and index-list of bidirectional functions

      INTEGER, intent(in) :: N_BRDF_KERNELS
      INTEGER, intent(in) :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER     , intent(in) :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(kind=8), intent(in) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL, intent(in) :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(kind=8), intent(in) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER, intent(in) :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, intent(in) :: DO_SHADOW_EFFECT

!  Multiple reflectance correction for direct beam flag
!              (only for GLITTER type kernels)

      LOGICAL, intent(in) :: DO_GLITTER_DBMS

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL, intent(in) :: DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL, intent(in) :: DO_KERNEL_PARAMS_WFS ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, intent(in) :: DO_KPARAMS_DERIVS  ( MAX_BRDF_KERNELS )

!  number of surfaceweighting functions

      INTEGER, intent(inout) :: N_SURFACE_WFS
      INTEGER, intent(in)    :: N_KERNEL_FACTOR_WFS
      INTEGER, intent(in)    :: N_KERNEL_PARAMS_WFS

!  Local angle control

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: N_USER_RELAZMS

!  Angles

      REAL(kind=8), intent(in) ::  BEAM_SZAS         (MAXBEAMS)
      REAL(kind=8), intent(in) ::  USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Output arguments
!  ================

!  Exact (direct bounce) BRDF (same all threads)
    
      REAL(kind=8), intent(out) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(kind=8), intent(out) :: BRDF_F_0 ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: BRDF_F   ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Exact (direct bounce) BRDF (same all threads)
    
      REAL(kind=8), intent(out) :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized Fourier components of BRDF, (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(kind=8), intent(out) :: LS_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: LS_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Fourier components of emissivity

      REAL(kind=8), intent(out) :: USER_EMISSIVITY(0:MAXMOMENTS,MAX_USER_STREAMS)
      REAL(kind=8), intent(out) :: EMISSIVITY     (0:MAXMOMENTS,MAXSTREAMS      )

!  Fourier components of emissivity

      REAL(kind=8), intent(out) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS,0:MAXMOMENTS,MAX_USER_STREAMS)
      REAL(kind=8), intent(out) :: LS_EMISSIVITY      ( MAX_SURFACEWFS,0:MAXMOMENTS,MAXSTREAMS      )

!  BRDF functions
!  --------------

!  ordinary BRDF without derivatives

      EXTERNAL       LAMBERTIAN_FUNCTION
      EXTERNAL       ROSSTHIN_FUNCTION
      EXTERNAL       ROSSTHICK_FUNCTION
      EXTERNAL       LISPARSE_FUNCTION
      EXTERNAL       LIDENSE_FUNCTION
      EXTERNAL       HAPKE_FUNCTION
      EXTERNAL       ROUJEAN_FUNCTION
      EXTERNAL       RAHMAN_FUNCTION
      EXTERNAL       COXMUNK_FUNCTION
      EXTERNAL       COXMUNK_FUNCTION_DB

!  new for Version 3.4R

      EXTERNAL       BREONVEG_FUNCTION
      EXTERNAL       BREONSOIL_FUNCTION

!  Hapke old uses exact DISORT code
!      EXTERNAL       HAPKE_FUNCTION_OLD

!  BRDFs with derivatives

      EXTERNAL       LISPARSE_FUNCTION_PLUS
      EXTERNAL       LIDENSE_FUNCTION_PLUS
      EXTERNAL       HAPKE_FUNCTION_PLUS
      EXTERNAL       RAHMAN_FUNCTION_PLUS
      EXTERNAL       COXMUNK_FUNCTION_PLUS
      EXTERNAL       COXMUNK_FUNCTION_PLUS_DB

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8)  :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  DB Kernel values

      REAL(kind=8)  :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8)  :: EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8)  :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Linearized Exact DB values

      REAL(kind=8)  :: D_DBKERNEL_BRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      REAL(kind=8)  :: PHIANG(MAX_USER_RELAZMS)
      REAL(kind=8)  :: COSPHI(MAX_USER_RELAZMS)
      REAL(kind=8)  :: SINPHI(MAX_USER_RELAZMS)

!  SZAs

      REAL(kind=8)  :: SZASURCOS(MAXBEAMS)
      REAL(kind=8)  :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates

      REAL(kind=8)  :: QUAD_STREAMS(MAXSTREAMS)
      REAL(kind=8)  :: QUAD_WEIGHTS(MAXSTREAMS)
      REAL(kind=8)  :: QUAD_SINES  (MAXSTREAMS)
      REAL(kind=8)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Viewing zenith streams

      REAL(kind=8)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(kind=8)  :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER       :: NBRDF_HALF
      REAL(kind=8)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8)  :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(kind=8)  :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8)  :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8)  :: SXE_BRDF ( MAXSTHALF_BRDF )

!  Azimuth factors

      REAL(kind=8)  :: BRDF_AZMFAC(MAXSTREAMS_BRDF)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: LOCAL_BRDF_F   ( MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8)  :: LOCAL_BRDF_F_0 ( MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8)  :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8)  :: LOCAL_USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8)  :: LOCAL_EMISSIVITY      ( MAXSTREAMS       )
      REAL(kind=8)  :: LOCAL_USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Local Derivative-kernel Fourier components
!  ==========================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8)  :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8)  :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8)  :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8)  :: D_LOCAL_EMISSIVITY      ( MAX_BRDF_PARAMETERS, MAXSTREAMS       )
      REAL(kind=8)  :: D_LOCAL_USER_EMISSIVITY ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS )

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=8)       :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  help

      INTEGER            :: K, Q, P, I, J, IB, UI, UM, IA, M
      INTEGER            :: QOFFSET ( MAX_BRDF_KERNELS)
      INTEGER            :: LOCAL_BRDF_NPARS, NMOMENTS
      REAL(kind=8)       :: LOCAL_BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL            :: LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )
      REAL(kind=8)       :: MUX, DELFAC, HELP_A, SUM, FF
      LOGICAL            :: ADD_FOURIER
      LOGICAL, parameter :: DO_BRDFQUAD_GAUSSIAN = .true.

!  Set up Quadrature streams

      CALL GAULEG(0.0d0,1.0d0, QUAD_STREAMS, QUAD_WEIGHTS, NSTREAMS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = DSQRT(1.0d0-QUAD_STREAMS(I)*QUAD_STREAMS(I))
        QUAD_STRMWTS(I) = QUAD_STREAMS(I) * QUAD_WEIGHTS(I)
      enddo

!  Number of Fourier components to calculate

      NMOMENTS = 2 * NSTREAMS - 1

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams
!    Warning, this shoudl be the BOA angle. OK for the non-refractive case.
!        
      DO IB = 1, NBEAMS
        MUX =  DCOS(BEAM_SZAS(IB)*DEG_TO_RAD)
        SZASURCOS(IB) = MUX
        SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
      ENDDO

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = DCOS(USER_ANGLES_INPUT(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = DSQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
        COSPHI(IA) = DCOS(PHIANG(IA))
        SINPHI(IA) = DCOS(PHIANG(IA))
      ENDDO

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian                                          &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )       ! Outputs
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid                                         &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )       ! Outputs
      ENDIF

!  Number of weighting functions, and offset

      Q = 0
      QOFFSET(1) = 0
      DO K = 1, N_BRDF_KERNELS
        IF ( DO_KERNEL_FACTOR_WFS(K) ) Q = Q + 1
        DO P = 1, N_BRDF_PARAMETERS(K)
          IF ( DO_KERNEL_PARAMS_WFS(K,P) ) Q = Q + 1
        ENDDO
        IF ( K.LT.N_BRDF_KERNELS ) QOFFSET(K+1) = Q
      ENDDO

      N_SURFACE_WFS = N_KERNEL_FACTOR_WFS + N_KERNEL_PARAMS_WFS
      IF ( Q .ne. N_SURFACE_WFS ) stop'bookkeeping wrong'

!  Initialise BRDF function
!  ------------------------

      DO M = 0, NMOMENTS
        DO I = 1, NSTREAMS
          DO IB = 1, NBEAMS
            BRDF_F_0(M,I,IB) = ZERO
          ENDDO
          DO J = 1, NSTREAMS
            BRDF_F(M,I,J) = ZERO
          ENDDO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            DO IB = 1, NBEAMS
              USER_BRDF_F_0(M,UM,IB) = ZERO
            ENDDO
            DO J = 1, NSTREAMS
                USER_BRDF_F(M,UM,J) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Compute Exact Direct Beam BRDF

      DO IA = 1, N_USER_RELAZMS
        DO IB = 1, NBEAMS
          DO UM = 1, N_USER_STREAMS
            EXACTDB_BRDFUNC(UM,IA,IB) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Initialize surface emissivity
!  -----------------------------

      IF ( DO_SURFACE_EMISSION ) THEN
        DO M = 0, NMOMENTS
          DO I = 1, NSTREAMS
            EMISSIVITY(M,I) = ONE
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_EMISSIVITY(M,UI) = ONE
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Copy parameter variables into local quantities

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO P = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(P) = BRDF_PARAMETERS(K,P)
        ENDDO
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO P = 1, MAX_BRDF_PARAMETERS
            LOCAL_BRDF_DERIVS(P) = DO_KERNEL_PARAMS_WFS(K,P)
          ENDDO
        ENDIF

!  Lambertian kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. LAMBERTIAN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER &
           ( LAMBERTIAN_FUNCTION, LAMBERTIAN_FUNCTION,             & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER &
           ( ROSSTHIN_FUNCTION, ROSSTHIN_FUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL LIDORT_BRDF_MAKER &
           ( ROSSTHICK_FUNCTION, ROSSTHICK_FUNCTION,               & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS &
           ( LISPARSE_FUNCTION_PLUS, LISPARSE_FUNCTION_PLUS,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL LIDORT_BRDF_MAKER &
           ( LISPARSE_FUNCTION, LISPARSE_FUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS &
           ( LIDENSE_FUNCTION_PLUS, LIDENSE_FUNCTION_PLUS,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL LIDORT_BRDF_MAKER &
           ( LIDENSE_FUNCTION, LIDENSE_FUNCTION,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS &
           ( HAPKE_FUNCTION_PLUS, HAPKE_FUNCTION_PLUS,                       & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL LIDORT_BRDF_MAKER &
           ( HAPKE_FUNCTION, HAPKE_FUNCTION,                       & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL LIDORT_BRDF_MAKER &
           ( ROUJEAN_FUNCTION, ROUJEAN_FUNCTION,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Rahman kernel: (3 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS &
           ( RAHMAN_FUNCTION_PLUS, RAHMAN_FUNCTION_PLUS,                     & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL LIDORT_BRDF_MAKER &
           ( RAHMAN_FUNCTION, RAHMAN_FUNCTION,                     & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Cox-Munk kernel: (2 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
         IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0d0
         IF (  DO_GLITTER_DBMS ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS &
           ( COXMUNK_FUNCTION_PLUS, COXMUNK_FUNCTION_PLUS_DB,                & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL LIDORT_BRDF_MAKER &
           ( COXMUNK_FUNCTION, COXMUNK_FUNCTION_DB,                & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
         ELSE
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL LIDORT_BRDF_MAKER_PLUS &
           ( COXMUNK_FUNCTION_PLUS, COXMUNK_FUNCTION_PLUS,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL LIDORT_BRDF_MAKER &
           ( COXMUNK_FUNCTION, COXMUNK_FUNCTION,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
         ENDIF
        ENDIF

!  Breon Vegetation kernel (0 free parameter)

        IF ( WHICH_BRDF(K) .EQ. BREONVEG_IDX ) THEN
          CALL LIDORT_BRDF_MAKER &
           ( BREONVEG_FUNCTION, BREONVEG_FUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Breon Soil kernel (0 free parameter)

        IF ( WHICH_BRDF(K) .EQ. BREONSOIL_IDX ) THEN
          CALL LIDORT_BRDF_MAKER &
           ( BREONSOIL_FUNCTION, BREONSOIL_FUNCTION,               & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Exact BRDFUNC
!  ------------

!  Compute Exact Direct Beam BRDF

        DO IA = 1, N_USER_RELAZMS
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              EXACTDB_BRDFUNC(UM,IA,IB) = EXACTDB_BRDFUNC(UM,IA,IB) &
                + BRDF_FACTORS(K) * DBKERNEL_BRDFUNC(UM,IA,IB)
            ENDDO
          ENDDO
        ENDDO

!  Linearization w.r.t Kernel Factor

        FF = BRDF_FACTORS(K)
        Q  = QOFFSET(K)
        IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
          Q = Q + 1
          DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                LS_EXACTDB_BRDFUNC(Q,UM,IA,IB) = DBKERNEL_BRDFUNC(UM,IA,IB)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Linearization w.r.t Kernel parameters

        DO P = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(P) ) THEN
            Q = Q + 1
            DO IA = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                LS_EXACTDB_BRDFUNC(Q,UM,IA,IB) = BRDF_FACTORS(K) * D_DBKERNEL_BRDFUNC(P,UM,IA,IB)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  Fourier Work now
!  ================

        DO M = 0, NMOMENTS

!  Fourier addition flag

          ADD_FOURIER = ( .not.LAMBERTIAN_KERNEL_FLAG(K) .or. &
                          (LAMBERTIAN_KERNEL_FLAG(K).AND.M.EQ.0) )

!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = ONE
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I)
            ENDDO
          ELSE
            DELFAC   = TWO
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I) * DCOS ( M * X_BRDF(I) )
            ENDDO
          ENDIF

!  Call

          CALL LIDORT_BRDF_FOURIER                                      &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LAMBERTIAN_KERNEL_FLAG(K), BRDF_FACTORS(K), M, DELFAC,       & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,           & ! Inputs
           EBRDFUNC, USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,      & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,                      & ! Outputs
           LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                      ! Outputs

!  Linear call

          IF ( LOCAL_BRDF_NPARS .gt. 0 ) then
            CALL LIDORT_BRDF_LS_FOURIER                                 &                
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,                         & ! Inputs
           LAMBERTIAN_KERNEL_FLAG(K), BRDF_FACTORS(K), M, DELFAC,       & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           D_BRDFUNC,  D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0,   & ! Inputs
           D_EBRDFUNC, D_USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,  & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,                       & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0,                  & ! Outputs
           D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )                 ! Outputs
          ENDIF

!  Spherical albedo (debug only)

          IF ( M .EQ. 0 ) THEN
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A = ZERO
              DO I = 1, NSTREAMS
               SUM = ZERO
               DO J = 1, NSTREAMS
                SUM = SUM + LOCAL_BRDF_F(I,J) * QUAD_STRMWTS(J)
               ENDDO
               HELP_A = HELP_A + SUM * QUAD_STRMWTS(I)
              ENDDO
              SPHERICAL_ALBEDO(K) = HELP_A*FOUR
             ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature reflectance)
!  ------------------------------------------------

!  Basic Kernel sum

            DO I = 1, NSTREAMS
              DO IB = 1, NBEAMS
                BRDF_F_0(M,I,IB) = BRDF_F_0(M,I,IB) &
                     + BRDF_FACTORS(K) * LOCAL_BRDF_F_0(I,IB)
              ENDDO
              DO J = 1, NSTREAMS
                BRDF_F(M,I,J) = BRDF_F(M,I,J) &
                     + BRDF_FACTORS(K) * LOCAL_BRDF_F(I,J)
              ENDDO
            ENDDO

!  Linearization w.r.t Kernel Factor

            FF = BRDF_FACTORS(K)
            Q  = QOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              Q = Q + 1
              DO I = 1, NSTREAMS
                DO IB = 1, NBEAMS
                  LS_BRDF_F_0(Q,M,I,IB) = LOCAL_BRDF_F_0(I,IB)
                ENDDO
                DO J = 1, NSTREAMS
                  LS_BRDF_F(Q,M,I,J) =  LOCAL_BRDF_F(I,J)
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                Q = Q + 1
                DO I = 1, NSTREAMS
                  DO IB = 1, NBEAMS
                    LS_BRDF_F_0(Q,M,I,IB) = FF*D_LOCAL_BRDF_F_0(P,I,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    LS_BRDF_F(Q,M,I,J) = FF*D_LOCAL_BRDF_F(P,I,J)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Kernel combinations (for user-stream reflectance)
!  -------------------------------------------------

!  Basci kernel summation

            IF ( DO_USER_STREAMS ) THEN
              DO UM = 1, N_USER_STREAMS
                DO IB = 1, NBEAMS
                  USER_BRDF_F_0(M,UM,IB) = USER_BRDF_F_0(M,UM,IB) &
                     + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F_0(UM,IB)
                ENDDO
                DO J = 1, NSTREAMS
                  USER_BRDF_F(M,UM,J) = USER_BRDF_F(M,UM,J) &
                       + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F(UM,J)
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel Factor

            Q  = QOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
             Q = Q + 1
             DO UM = 1, N_USER_STREAMS
              DO IB = 1, NBEAMS
               LS_USER_BRDF_F_0(Q,M,UM,IB) = LOCAL_USER_BRDF_F_0(UM,IB)
              ENDDO
              DO J = 1, NSTREAMS
               LS_USER_BRDF_F(Q,M,UM,J) =  LOCAL_USER_BRDF_F(UM,J)
              ENDDO
             ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                Q = Q + 1
                DO UM = 1, N_USER_STREAMS
                  DO IB = 1, NBEAMS
                    LS_USER_BRDF_F_0(Q,M,UM,IB) = FF * D_LOCAL_USER_BRDF_F_0(P,UM,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    LS_USER_BRDF_F(Q,M,UM,J) = FF * D_LOCAL_USER_BRDF_F(P,UM,J)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Total emissivities
!  ------------------

!  only if flagged

            IF ( DO_SURFACE_EMISSION ) THEN

!  Basci kernel contributions

              DO I = 1, NSTREAMS
               EMISSIVITY(M,I) = EMISSIVITY(M,I) - LOCAL_EMISSIVITY(I)
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  USER_EMISSIVITY(M,UI) = USER_EMISSIVITY(M,UI) - LOCAL_USER_EMISSIVITY(UI)
                ENDDO
              ENDIF

!  Linearization w.r.t Kernel Factor

              Q  = QOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                DO I = 1, NSTREAMS
                  LS_EMISSIVITY(Q,M,I) = - LOCAL_EMISSIVITY(I) / FF
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                  DO UI = 1, N_USER_STREAMS
                    LS_USER_EMISSIVITY(Q,M,UI) = - LOCAL_USER_EMISSIVITY(UI) / FF
                  ENDDO
                ENDIF
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  DO I = 1, NSTREAMS
                    LS_EMISSIVITY(Q,M,I) = - D_LOCAL_EMISSIVITY(P,I)
                  ENDDO
                  IF ( DO_USER_STREAMS ) THEN
                    DO UI = 1, N_USER_STREAMS
                      LS_USER_EMISSIVITY(Q,M,UI) = - D_LOCAL_USER_EMISSIVITY(P,UI)
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO

!  End emissivity clause

            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_LS_MASTER

!

SUBROUTINE LIDORT_BRDF_MAKER_PLUS                                            &                                   
           ( BRDF_FUNCTION_PLUS, BRDF_FUNCTION_PLUS_DB,                      & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                            & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                         & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                   & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,           & ! Inputs
             EXACTDB_BRDFUNC, BRDFUNC, USER_BRDFUNC,                         & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_EXACTDB_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                   & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs

      implicit none

!  Prepares the bidirectional reflectance scatter matrices

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input arguments
!  ===============

!  BRDF functions (external calls)

      EXTERNAL         BRDF_FUNCTION_PLUS
      EXTERNAL         BRDF_FUNCTION_PLUS_DB

!  Local flags

      LOGICAL, intent(in) :: DO_USER_STREAMS
      LOGICAL, intent(in) :: DO_SURFACE_EMISSION

!  Local angle control

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: N_USER_RELAZMS

!  Local angles

      REAL(kind=8), intent(in) ::  PHIANG(MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  COSPHI(MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  SINPHI(MAX_USER_RELAZMS)

      REAL(kind=8), intent(in) ::  SZASURCOS(MAXBEAMS)
      REAL(kind=8), intent(in) ::  SZASURSIN(MAXBEAMS)

      REAL(kind=8), intent(in) ::  QUAD_STREAMS(MAXSTREAMS)
      REAL(kind=8), intent(in) ::  QUAD_SINES  (MAXSTREAMS)

      REAL(kind=8), intent(in) ::  USER_STREAMS(MAX_USER_STREAMS)
      REAL(kind=8), intent(in) ::  USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      INTEGER     , intent(in) ::  NSTREAMS_BRDF
      INTEGER     , intent(in) ::  NBRDF_HALF
      REAL(kind=8), intent(in) ::  X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8), intent(in) ::  SXE_BRDF ( MAXSTHALF_BRDF )

!  Local number of parameters and local parameter array

      INTEGER     , intent(in) ::  LOCAL_BRDF_NPARS
      REAL(kind=8), intent(in) ::  LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )
      LOGICAL     , intent(in) ::  LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(kind=8), intent(out) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8), intent(out) :: EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(kind=8), intent(out) :: D_EXACTDB_BRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8), intent(out) :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER          :: Q, I, UI, J, K, KE, IB
      REAL(kind=8)     :: DFUNC ( MAX_BRDF_PARAMETERS )

!  Exact DB calculation
!  --------------------

      DO K = 1, N_USER_RELAZMS
         DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
               CALL BRDF_FUNCTION_PLUS_DB &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  & ! Inputs
                 LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         & ! Inputs
                 USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),        & ! Inputs
                 EXACTDB_BRDFUNC(UI,K,IB), DFUNC )                         ! Output
              DO Q = 1, LOCAL_BRDF_NPARS
                D_EXACTDB_BRDFUNC(Q,UI,K,IB)  = DFUNC(Q)
              ENDDO
            ENDDO
         ENDDO
      ENDDO

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  &
                 LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 &
                 SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),          &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC_0(I,IB,K), DFUNC )
            DO Q = 1, LOCAL_BRDF_NPARS
              D_BRDFUNC_0(Q,I,IB,K)  = DFUNC(Q)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  &
                 LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 &
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),        &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC(I,J,K), DFUNC )
            DO Q = 1, LOCAL_BRDF_NPARS
              D_BRDFUNC(Q,I,J,K)  = DFUNC(Q)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  &
                 LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 &
                 CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),            &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 EBRDFUNC(I,KE,K), DFUNC )
              DO Q = 1, LOCAL_BRDF_NPARS
                D_EBRDFUNC(Q,I,KE,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam

        DO IB = 1, NBEAMS
          DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  &
                 LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 USER_BRDFUNC_0(UI,IB,K), DFUNC )
              DO Q = 1, LOCAL_BRDF_NPARS
                D_USER_BRDFUNC_0(Q,UI,IB,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  &
                   LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 &
                   QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),       &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_BRDFUNC(UI,J,K), DFUNC )
              DO Q = 1, LOCAL_BRDF_NPARS
                D_USER_BRDFUNC(Q,UI,J,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS,                  &
                 LOCAL_BRDF_PARS,     LOCAL_BRDF_DERIVS,                 &
                   CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),         &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),    &
                   USER_EBRDFUNC(UI,KE,K), DFUNC )
                DO Q = 1, LOCAL_BRDF_NPARS
                  D_USER_EBRDFUNC(Q,UI,KE,K)  = DFUNC(Q)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_MAKER_PLUS

!

SUBROUTINE LIDORT_BRDF_LS_FOURIER                                       &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,                         & ! Inputs
           LAMBERTIAN_FLAG, FACTOR, M, DELFAC,                          & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0,    & ! Inputs
           D_EBRDFUNC, D_USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,  & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,                       & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0,                  & ! Outputs
           D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )                 ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input arguments
!  ===============

!  Control

      LOGICAL, intent(in)      :: LAMBERTIAN_FLAG
      LOGICAL, intent(in)      :: DO_USER_STREAMS
      LOGICAL, intent(in)      :: DO_SURFACE_EMISSION
      REAL(kind=8), intent(in) :: DELFAC, FACTOR
      INTEGER, intent(in)      :: M, LOCAL_BRDF_NPARS
      LOGICAL, intent(in)      :: LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Local numbers

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(kind=8), intent(in) :: BRDF_AZMFAC ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  A_BRDF      ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) :: BAX_BRDF    ( MAXSTHALF_BRDF  )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(in) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF ) 
      REAL(kind=8), intent(in) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(in) :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS ,MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=8), intent(in) :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output: Derivative-kernel Fourier components
!  ============================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8), intent(out) :: D_LOCAL_EMISSIVITY      ( MAX_BRDF_PARAMETERS, MAXSTREAMS       )
      REAL(kind=8), intent(out) :: D_LOCAL_USER_EMISSIVITY ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS )

!  local variables
!  ===============

      INTEGER         :: I, UI, J, K, KPHI, IB, Q
      REAL(kind=8)    :: SUM, REFL, HELP

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
            DO IB = 1, NBEAMS
              DO I = 1, NSTREAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + D_BRDFUNC_0(Q,I,IB,K)*BRDF_AZMFAC(K)
                ENDDO
                D_LOCAL_BRDF_F_0(Q,I,IB) = SUM * HELP
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + D_BRDFUNC(Q,I,J,K) * BRDF_AZMFAC(K)
                ENDDO
                D_LOCAL_BRDF_F(Q,I,J) = SUM * HELP
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                    SUM = SUM+D_USER_BRDFUNC_0(Q,UI,IB,K)*BRDF_AZMFAC(K)
                  ENDDO
                  D_LOCAL_USER_BRDF_F_0(Q,UI,IB) = SUM * HELP
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO UI = 1, N_USER_STREAMS
                DO J = 1, NSTREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                    SUM = SUM + D_USER_BRDFUNC(Q,UI,J,K)*BRDF_AZMFAC(K)
                  ENDDO
                  D_LOCAL_USER_BRDF_F(Q,UI,J) = SUM * HELP
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!  Emissivity
!  ----------

!  Assumed to exist only for the total intensity
!        (first element of Stokes Vector) - is this right ??????

      IF ( DO_SURFACE_EMISSION ) THEN

!  Lambertian case

        IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO I = 1, NSTREAMS
                D_LOCAL_EMISSIVITY(Q,I) = ZERO
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  D_LOCAL_USER_EMISSIVITY(Q,UI) = ZERO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!  bidirectional reflectance

        IF ( .not. LAMBERTIAN_FLAG ) THEN

!  Quadrature polar directions

          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO I = 1, NSTREAMS
                REFL = ZERO
                DO KPHI= 1, NSTREAMS_BRDF
                  SUM = ZERO
                  DO K = 1, NBRDF_HALF
                    SUM = SUM + D_EBRDFUNC(Q,I,K,KPHI) * BAX_BRDF(K)
                  ENDDO
                  REFL = REFL + A_BRDF(KPHI) * SUM
                ENDDO
                D_LOCAL_EMISSIVITY(Q,I) = REFL * FACTOR
              ENDDO
            ENDIF
          ENDDO

!   user-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO Q = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
                DO UI = 1, N_USER_STREAMS
                  REFL = ZERO
                  DO KPHI= 1, NSTREAMS_BRDF
                    SUM = ZERO
                    DO K = 1, NBRDF_HALF
                      SUM = SUM+D_USER_EBRDFUNC(Q,UI,K,KPHI)*BAX_BRDF(K)
                    ENDDO
                    REFL = REFL + A_BRDF(KPHI) * SUM
                  ENDDO
                  D_LOCAL_USER_EMISSIVITY(Q,UI) = REFL * FACTOR
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Not lambertian

        ENDIF

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_LS_FOURIER

!

SUBROUTINE LIDORT_BRDF_INPUTS_PLUS ( FILNAM,                              & ! Inputs
           DO_USER_STREAMS, DO_BRDF_SURFACE, DO_SURFACE_EMISSION,         & ! Outputs
           N_BRDF_KERNELS, WHICH_BRDF, BRDF_NAMES,                        & ! Outputs
           LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS,           & ! Outputs
           N_BRDF_PARAMETERS, BRDF_PARAMETERS,                            & ! Outputs
           DO_SHADOW_EFFECT, DO_GLITTER_DBMS,                             & ! Outputs
           DO_KERNEL_PARAMS_WFS, DO_KERNEL_FACTOR_WFS, DO_KPARAMS_DERIVS, & ! Outputs
           N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS, N_SURFACE_WFS,       & ! Outputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,              & ! Outputs
           BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,                    & ! Outputs
           STATUS, NMESSAGES, MESSAGES, ACTIONS )                           ! Outputs

!  Input routine for BRDF program

      implicit none

!  Include file of Dimensions

      INCLUDE 'LIDORT.PARS_F90'

!  Module arguments (input filename)

      CHARACTER*(*), intent(in)  :: FILNAM

!  stream angle flag

      LOGICAL, intent(out) :: DO_USER_STREAMS

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL, intent(out) :: DO_BRDF_SURFACE

!  Surface emission

      LOGICAL, intent(out) :: DO_SURFACE_EMISSION

!   Number and index-list and names of bidirectional functions

      INTEGER, intent(out)      :: N_BRDF_KERNELS
      INTEGER, intent(out)      :: WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER*10, intent(out) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER     , intent(out) :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(kind=8), intent(out) ::  BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL, intent(out) :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(kind=8), intent(out) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER, intent(out) :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, intent(out) :: DO_SHADOW_EFFECT

!  Multiple reflectance correction for direct beam flag
!              (only for GLITTER type kernels)

      LOGICAL, intent(out) :: DO_GLITTER_DBMS

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL, intent(out) :: DO_KERNEL_FACTOR_WFS  ( MAX_BRDF_KERNELS )
      LOGICAL, intent(out) :: DO_KERNEL_PARAMS_WFS  ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, intent(out) :: DO_KPARAMS_DERIVS  ( MAX_BRDF_KERNELS )

!  number of surface weighting functions

      INTEGER, intent(out) :: N_SURFACE_WFS
      INTEGER, intent(out) :: N_KERNEL_FACTOR_WFS
      INTEGER, intent(out) :: N_KERNEL_PARAMS_WFS

!  Number of discrete ordinate streams

      INTEGER, intent(out) :: NSTREAMS

!  Local angle control

      INTEGER, intent(out) :: NBEAMS
      INTEGER, intent(out) :: N_USER_STREAMS
      INTEGER, intent(out) :: N_USER_RELAZMS

!  Angles

      REAL(kind=8), intent(out) :: BEAM_SZAS         (MAXBEAMS)
      REAL(kind=8), intent(out) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(kind=8), intent(out) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Exception handling. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER, intent(out)        :: STATUS
      INTEGER, intent(out)        :: NMESSAGES
      CHARACTER*(*), intent(out)  :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(out)  :: ACTIONS (0:MAX_MESSAGES)

!  local variables
!  ===============

      CHARACTER(Len=8), parameter ::  PREFIX = 'LIDORT -' 

      INTEGER           :: DUM_INDEX, DUM_NPARS
      CHARACTER(Len=10) :: DUM_NAME
      LOGICAL           :: ERROR
      CHARACTER(Len=80) :: PAR_STR
      LOGICAL           :: GFINDPAR
      INTEGER           :: I, J, K, L, FILUNIT, LEN_STRING, NM

      EXTERNAL             GFINDPAR
      EXTERNAL             LEN_STRING

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0

      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO
      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES_INPUT(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0

      DO_SHADOW_EFFECT    = .FALSE.
      DO_GLITTER_DBMS     = .FALSE.
      DO_SURFACE_EMISSION = .FALSE.

      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        BRDF_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

      N_SURFACE_WFS  = 0
      N_KERNEL_FACTOR_WFS = 0
      N_KERNEL_PARAMS_WFS = 0
      DO K = 1, MAX_BRDF_KERNELS
        DO_KPARAMS_DERIVS(K) = .false.
        DO_KERNEL_FACTOR_WFS(K) = .FALSE.
        DO L = 1, MAX_BRDF_PARAMETERS
          DO_KERNEL_PARAMS_WFS(K,L) = .FALSE.
        ENDDO
      ENDDO

!  Read Angle stuff
!  ================

!  user-defined Stream angle

      PAR_STR = 'User-defined stream angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Solar beams
!  ===========

!  number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'TOA solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Azimuth angles
!  ==============

!  Number of angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  =  'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  Angles

      PAR_STR = 'User-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  User defined stream angles (should be positive)
!  ==========================

      IF ( DO_USER_STREAMS ) THEN

!  Number of angles

        PAR_STR = 'Number of user-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  Angles

        PAR_STR = 'User-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  BRDF inputs

      IF ( DO_BRDF_SURFACE ) THEN

!  number of kernels, check this value

        PAR_STR = 'Number of bidirectional reflectance kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of BRDF Kernels > maximum dimension (=3)'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of bidirectional reflectance streams'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          NM = NM + 1
          MESSAGES(NM) =  'Number of  BRDF streams > maximum dimension'
          ACTIONS(NM)  =  'Re-set input value or increase MAXSTREAMS_BRDF dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  Main kernel input

        PAR_STR = 'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
               BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I), &
              N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,3)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )

!  Set the Lambertian kernel flags

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
            LAMBERTIAN_KERNEL_FLAG(I) = .true.
          ENDIF
        ENDDO

!  Shadowing input (for Cox-Munk type)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do shadow effect for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Multiple reflectance DB correction (for Cox-Munk type)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do multiple reflectance for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_GLITTER_DBMS
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Linearized input
!  ----------------

        PAR_STR = 'Kernels, indices, # pars, Jacobian flags'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,57,ERR=998) &
            DUM_NAME, DUM_INDEX,DUM_NPARS,DO_KERNEL_FACTOR_WFS(I), &
                (DO_KERNEL_PARAMS_WFS(I,J),J=1,3)

            IF ( DUM_NAME .NE. BRDF_NAMES(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Kernel name not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel name'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

            IF ( DUM_INDEX .NE. WHICH_BRDF(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Index name not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel Index'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

            IF ( DUM_NPARS .NE. N_BRDF_PARAMETERS(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input Number of BRDF parameters not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of N_BRDF_PARAMETERS'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

!  Compute total number of pars

            IF ( DO_KERNEL_FACTOR_WFS(I) ) THEN
              N_KERNEL_FACTOR_WFS = N_KERNEL_FACTOR_WFS  + 1
            ENDIF
            DO J = 1, N_BRDF_PARAMETERS(I)
              IF ( DO_KERNEL_PARAMS_WFS(I,J) ) THEN
                N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1
              ENDIF
            ENDDO
            DO_KPARAMS_DERIVS(I) = (N_KERNEL_PARAMS_WFS.GT.0)

          ENDDO
          N_SURFACE_WFS = N_KERNEL_FACTOR_WFS+N_KERNEL_PARAMS_WFS
        ENDIF

        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 57     FORMAT( A10, I3, I2, 1X, L2, 2X, 3L2 ) 

!  Check total number of BRDF weighting functions is not out of bounds

        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Surface WFs > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SURFACEWFS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  End BRDF clause
   
      ENDIF

!  Successful finish

      CLOSE(FILUNIT)
      RETURN

!  Open file error

300   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      RETURN

!  line read error - abort immediately

998   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_INPUTS_PLUS

