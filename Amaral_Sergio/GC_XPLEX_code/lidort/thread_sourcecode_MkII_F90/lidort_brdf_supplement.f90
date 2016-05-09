!$Id: lidort_brdf_supplement.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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
! #            LIDORT_BRDF_INPUTS                               #
! #                                                             #
! #            LIDORT_BRDF_MASTER (master), calling             #
! #                                                             #
! #              LIDORT_BRDF_MAKER                              #
! #              BRDF_QUADRATURE_Gaussian                       #
! #              BRDF_QUADRATURE_Trapezoid (not used)           #
! #              LIDORT_BRDF_FOURIER                            #
! #                                                             #
! ###############################################################

SUBROUTINE LIDORT_BRDF_MASTER                                  &
        ( DO_USER_STREAMS, DO_SHADOW_EFFECT, DO_GLITTER_DBMS,  & ! Inputs
          DO_SURFACE_EMISSION, N_BRDF_KERNELS, WHICH_BRDF,     & ! Inputs
          LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS, & ! Inputs
          N_BRDF_PARAMETERS, BRDF_PARAMETERS,                  & ! Inputs
          NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,    & ! Inputs
          BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,          & ! Inputs
          BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,        & ! Outputs
          EXACTDB_BRDFUNC, EMISSIVITY, USER_EMISSIVITY )         ! Outputs

!  Prepares the bidirectional reflectance functions
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

!  Fourier components of emissivity

      REAL(kind=8), intent(out) :: USER_EMISSIVITY(0:MAXMOMENTS,MAX_USER_STREAMS)
      REAL(kind=8), intent(out) :: EMISSIVITY     (0:MAXMOMENTS,MAXSTREAMS      )

!  BRDF External functions
!  =======================

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

      REAL(kind=8)       :: LOCAL_BRDF_F   ( MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8)       :: LOCAL_BRDF_F_0 ( MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8)       :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8)       :: LOCAL_USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8)       :: LOCAL_EMISSIVITY      ( MAXSTREAMS       )
      REAL(kind=8)       :: LOCAL_USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=8)       :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  help

      INTEGER            :: K, B, I, J, IB, UI, UM, IA, M
      INTEGER            :: LOCAL_BRDF_NPARS, NMOMENTS
      REAL(kind=8)       :: LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )
      REAL(kind=8)       :: MUX, DELFAC, HELP_A, SUM
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

!  Local variables

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(B) = BRDF_PARAMETERS(K,B)
        ENDDO

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

!  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
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

!  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
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

!  Rahman kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
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

!  Cox-Munk kernel: (2 free parameters).
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0d0
          IF ( DO_GLITTER_DBMS ) THEN
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

!  Compute Exact Direct Beam BRDF

        DO IA = 1, N_USER_RELAZMS
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              EXACTDB_BRDFUNC(UM,IA,IB) = EXACTDB_BRDFUNC(UM,IA,IB) + BRDF_FACTORS(K) * DBKERNEL_BRDFUNC(UM,IA,IB)
            ENDDO
          ENDDO
        ENDDO

!  Fourier Work now
!  ================

        DO M = 0, NMOMENTS

!  Fourier addition flag

          ADD_FOURIER = ( .not.LAMBERTIAN_KERNEL_FLAG(K) .or. (LAMBERTIAN_KERNEL_FLAG(K).AND.M.EQ.0) )

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

!  Spherical albedo

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

            DO I = 1, NSTREAMS
              DO IB = 1, NBEAMS
                BRDF_F_0(M,I,IB) = BRDF_F_0(M,I,IB) + BRDF_FACTORS(K) * LOCAL_BRDF_F_0(I,IB)
              ENDDO
              DO J = 1, NSTREAMS
                BRDF_F(M,I,J) = BRDF_F(M,I,J) + BRDF_FACTORS(K) * LOCAL_BRDF_F(I,J)
              ENDDO
            ENDDO

!  Kernel combinations (for user-stream reflectance)

            IF ( DO_USER_STREAMS ) THEN
              DO UM = 1, N_USER_STREAMS
                DO IB = 1, NBEAMS
                  USER_BRDF_F_0(M,UM,IB) = USER_BRDF_F_0(M,UM,IB) + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F_0(UM,IB)
                ENDDO
                DO J = 1, NSTREAMS
                  USER_BRDF_F(M,UM,J) = USER_BRDF_F(M,UM,J) + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F(UM,J)
                ENDDO
              ENDDO
            ENDIF

!  Total emissivities

            IF ( DO_SURFACE_EMISSION ) THEN
              DO I = 1, NSTREAMS
               EMISSIVITY(M,I) = EMISSIVITY(M,I) - LOCAL_EMISSIVITY(I)
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  USER_EMISSIVITY(M,UI) = USER_EMISSIVITY(M,UI) - LOCAL_USER_EMISSIVITY(UI)
                ENDDO
              ENDIF
            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_MASTER

SUBROUTINE BRDF_QUADRATURE_GAUSSIAN                                          &
        ( DO_BRDF_SURFEMISSION, NSTREAMS_BRDF, NBRDF_HALF,                   & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )     ! Outputs

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input
!  =====

!  Emission flag

      LOGICAL, intent(in) :: DO_BRDF_SURFEMISSION

!  Number of streams

      INTEGER, intent(in) :: NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      REAL(kind=8), intent(out) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: A_BRDF  ( MAXSTREAMS_BRDF )

!  For emission calculations

      REAL(kind=8), intent(out) :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8), intent(out) :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8), intent(out) :: SXE_BRDF ( MAXSTHALF_BRDF )

!  local variables
!  ---------------

      INTEGER          :: I, I1, K

!  BRDF quadrature (Gauss-Legendre)
!  ---------------

!  Save these quantities for efficient coding

      CALL GAULEG ( ZERO, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
        DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
          X_BRDF(I1) = - X_BRDF(I)
          A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = DSQRT(ONE-X_BRDF(I)*X_BRDF(I))
        ENDDO
        DO I = 1, NSTREAMS_BRDF
       X_BRDF(I) = PIE * X_BRDF(I)
       CX_BRDF(I) = DCOS ( X_BRDF(I) )
       SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_BRDF_SURFEMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BRDF_QUADRATURE_GAUSSIAN

!

SUBROUTINE BRDF_QUADRATURE_TRAPEZOID                                         &
        ( DO_BRDF_SURFEMISSION, NSTREAMS_BRDF, NBRDF_HALF,                   & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )     ! Outputs

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input
!  =====

!  Emission flag

      LOGICAL, intent(in) :: DO_BRDF_SURFEMISSION

!  Number of streams

      INTEGER, intent(in) :: NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      REAL(kind=8), intent(out) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: A_BRDF  ( MAXSTREAMS_BRDF )

!  For emission calculations

      REAL(kind=8), intent(out) :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8), intent(out) :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8), intent(out) :: SXE_BRDF ( MAXSTHALF_BRDF )

!  local variables
!  ---------------

      INTEGER       :: I, I1, K
      REAL(kind=8)  :: DF1, DEL

!  BRDF quadrature (Trapezium)
!  ---------------

!  Save these quantities for efficient coding

      DF1 = DBLE(NSTREAMS_BRDF - 1 )
      DEL = TWO * PIE / DF1
        DO I = 1, NSTREAMS_BRDF
        I1 = I - 1
        X_BRDF(I) = DBLE(I1) * DEL - PIE
        X_BRDF(I) = DBLE(I1) * DEL
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
        CXE_BRDF(I) = CX_BRDF(I)
        SXE_BRDF(I) = SX_BRDF(I)
      ENDDO
        DO I = 2, NSTREAMS_BRDF - 1
        A_BRDF(I)  = DEL / PIE
      ENDDO
      A_BRDF(1)              = DEL * HALF / PIE
      A_BRDF(NSTREAMS_BRDF)  = DEL * HALF / PIE

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_BRDF_SURFEMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BRDF_QUADRATURE_TRAPEZOID

!

SUBROUTINE LIDORT_BRDF_MAKER                                       &                                   
           ( BRDF_FUNCTION, BRDF_FUNCTION_DB,                      & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                  & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,         & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                    & ! Inputs
             EXACTDB_BRDFUNC, BRDFUNC, USER_BRDFUNC,               & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs

      implicit none

!  Prepares the bidirectional reflectance scatter matrices

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Input arguments
!  ===============

!  BRDF functions (external calls)

      EXTERNAL         BRDF_FUNCTION
      EXTERNAL         BRDF_FUNCTION_DB

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

!  local variables
!  ---------------

      INTEGER        :: I, UI, J, K, KE, IB

!  Exact DB calculation
!  --------------------

      DO K = 1, N_USER_RELAZMS
         DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
               CALL BRDF_FUNCTION_DB &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         & ! Inputs
                 USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),        & ! Inputs
                 EXACTDB_BRDFUNC(UI,K,IB) )                                ! Output
            ENDDO
         ENDDO
      ENDDO

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
       DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
         CALL BRDF_FUNCTION &
             ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, &
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),          &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
               BRDFUNC_0(I,IB,K) )
        ENDDO
       ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, &
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),        &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC(I,J,K) )
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
                 ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, &
                   CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),            &
                   QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                   EBRDFUNC(I,KE,K) )
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
            CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 USER_BRDFUNC_0(UI,IB,K) )
          ENDDO
         ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
                 ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, &
                   QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),       &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_BRDFUNC(UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
                 ( MAX_BRDF_PARAMETERS, LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, &
                   CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),           &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_EBRDFUNC(UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_MAKER

!

SUBROUTINE LIDORT_BRDF_FOURIER                                          &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LAMBERTIAN_FLAG, FACTOR, M, DELFAC,                          & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,           & ! Inputs
           EBRDFUNC, USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,      & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,                      & ! Outputs
           LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                      ! Outputs

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
      INTEGER, intent(in)      :: M

!  Local numbers

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(kind=8), intent(in) ::  BRDF_AZMFAC ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  A_BRDF      ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  BAX_BRDF    ( MAXSTHALF_BRDF  )

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(in) ::  BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(in) ::  USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=8), intent(in) ::  EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: LOCAL_BRDF_F   ( MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: LOCAL_BRDF_F_0 ( MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: LOCAL_USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8), intent(out) :: LOCAL_EMISSIVITY      ( MAXSTREAMS       )
      REAL(kind=8), intent(out) :: LOCAL_USER_EMISSIVITY ( MAX_USER_STREAMS )

!  local variables
!  ===============

      INTEGER       :: I, UI, J, K, KPHI, IB
      REAL(kind=8)  :: SUM, REFL, HELP

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO IB = 1, NBEAMS
          DO I = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC_0(I,IB,K)*BRDF_AZMFAC(K)
            ENDDO
            LOCAL_BRDF_F_0(I,IB) = SUM * HELP
          ENDDO
        ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
        DO IB = 1, NBEAMS
          DO I = 1, NSTREAMS
            LOCAL_BRDF_F_0(I,IB) = ONE
          ENDDO
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC(I,J,K) * BRDF_AZMFAC(K)
            ENDDO
            LOCAL_BRDF_F(I,J) = SUM * HELP
          ENDDO
        ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            LOCAL_BRDF_F(I,J) = ONE
          ENDDO
        ENDDO
      ENDIF

!  debug information

!      IF ( DO_DEBUG_WRITE ) THEN
!        WRITE(555,'(A)')'BRDF_1 Fourier 0 quad values'
!        IF ( FOURIER .EQ. 0 ) THEN
!          DO I = 1, NSTREAMS
!          WRITE(555,'(1PE12.5,3x,1P10E12.5)') BIREFLEC_0(1,I,1),(BIREFLEC(1,I,J),J=1,NSTREAMS)
!         ENDDO
!        ENDIF
!      ENDIF

!  albedo check, always calculate the spherical albedo.
!   (Plane albedo calculations are commented out)


!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + USER_BRDFUNC_0(UI,IB,K) * BRDF_AZMFAC(K)
              ENDDO
              LOCAL_USER_BRDF_F_0(UI,IB) = SUM * HELP
            ENDDO
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              LOCAL_USER_BRDF_F_0(UI,IB) = ONE
            ENDDO
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + USER_BRDFUNC(UI,J,K) * BRDF_AZMFAC(K)
              ENDDO
              LOCAL_USER_BRDF_F(UI,J) = SUM * HELP
            ENDDO
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              LOCAL_USER_BRDF_F(UI,J) = ONE
            ENDDO
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
          DO I = 1, NSTREAMS
            LOCAL_EMISSIVITY(I) = FACTOR
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              LOCAL_USER_EMISSIVITY(UI) = FACTOR
            ENDDO
          ENDIF
        ENDIF

!  bidirectional reflectance

        IF ( .not. LAMBERTIAN_FLAG ) THEN

!  Quadrature polar directions

          DO I = 1, NSTREAMS
            REFL = ZERO
            DO KPHI= 1, NSTREAMS_BRDF
              SUM = ZERO
              DO K = 1, NBRDF_HALF
                SUM = SUM + EBRDFUNC(I,K,KPHI) * BAX_BRDF(K)
              ENDDO
              REFL = REFL + A_BRDF(KPHI) * SUM
            ENDDO
            LOCAL_EMISSIVITY(I) = REFL * FACTOR
          ENDDO

!   user-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              REFL = ZERO
              DO KPHI= 1, NSTREAMS_BRDF
                SUM = ZERO
                DO K = 1, NBRDF_HALF
                  SUM = SUM + USER_EBRDFUNC(UI,K,KPHI)*BAX_BRDF(K)
                ENDDO
                REFL = REFL + A_BRDF(KPHI) * SUM
              ENDDO
              LOCAL_USER_EMISSIVITY(UI) = REFL * FACTOR
            ENDDO
          ENDIF

        ENDIF

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_BRDF_FOURIER

!

SUBROUTINE LIDORT_BRDF_INPUTS ( FILNAM,                             & ! Inputs
           DO_USER_STREAMS, DO_BRDF_SURFACE, DO_SURFACE_EMISSION,   & ! Outputs
           N_BRDF_KERNELS, WHICH_BRDF, BRDF_NAMES,                  & ! Outputs
           LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS,     & ! Outputs
           N_BRDF_PARAMETERS, BRDF_PARAMETERS,                      & ! Outputs
           DO_SHADOW_EFFECT, DO_GLITTER_DBMS,                       & ! Outputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,        & ! Outputs
           BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,              & ! Outputs
           STATUS, NMESSAGES, MESSAGES, ACTIONS )                     ! Outputs

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

      LOGICAL           :: ERROR
      CHARACTER(Len=80) :: PAR_STR
      LOGICAL           :: GFINDPAR
      INTEGER           :: I, K, L, FILUNIT, LEN_STRING, NM

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
END SUBROUTINE LIDORT_BRDF_INPUTS

