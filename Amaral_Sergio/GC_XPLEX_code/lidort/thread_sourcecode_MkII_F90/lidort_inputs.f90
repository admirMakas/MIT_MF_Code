!$Id: lidort_inputs.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #        --           -            -        -        -    #
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_INPUTS_BASIC_MASTER (master), calling:    #
! #             LIDORT_INIT_INPUTS_BASIC                        #
! #             LIDORT_READ_INPUTS_BASIC                        #
! #                                                             #
! #            LIDORT_CHECK_INPUTS_BASIC                        #
! #            LIDORT_CHECK_INPUTS_THREAD                       #
! #            LIDORT_DERIVE_INPUTS_BASIC                       #
! #                                                             #
! ###############################################################

SUBROUTINE LIDORT_INPUTS_BASIC_MASTER                              &
          ( FILNAM,                                                & ! INPUT
            DO_SOLAR_SOURCES,                                      & ! output
            DO_UPWELLING,            DO_DNWELLING,                 & ! output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,              & ! output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,           & ! output
            DO_SSFULL,               DO_SSCORR_TRUNCATION,         & ! output
            DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,              & ! output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                & ! output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,       & ! output 
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,           & ! output  
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,            & ! output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,           & ! output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                & ! output
            NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,        & ! output
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! output
            FLUX_FACTOR, LIDORT_ACCURACY,                          & ! output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,  & ! output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_LEVELS,     & ! output
            STATUS, NMESSAGES, MESSAGES, ACTIONS )                   ! Output

!  Read all control inputs for LIDORT

      IMPLICIT NONE

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Input argument (input filename)
!  -------------------------------

      CHARACTER*(*), intent(in)   :: FILNAM

!  OUTPUTS
!  -------

!  1. CONTROL FLAGS
!  ================

!  Basic top-level control
!   -- Thermal to be reintroduced later 2010

      LOGICAL, intent(out) ::    DO_SOLAR_SOURCES
!      LOGICAL, intent(out) ::    DO_THERMAL_EMISSION

!  directional control

      LOGICAL, intent(out) ::    DO_UPWELLING
      LOGICAL, intent(out) ::    DO_DNWELLING

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL, intent(out) ::    DO_FULLRAD_MODE

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL, intent(out) ::    DO_USER_STREAMS


!  single scatter corrections. (Includes direct beam reflection)
!    - NADIR    : Plane-parallel line-of-sight
!    - OUTGOING : Line-of-sight in curved atmosphere

      LOGICAL, intent(out) ::    DO_SSCORR_NADIR
      LOGICAL, intent(out) ::    DO_SSCORR_OUTGOING

!  Flag for Full-up single scatter calculation. (No MS field)
!    One of the above two SSCORR flags must be set

      LOGICAL, intent(out) ::    DO_SSFULL

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations. **** Use with CAUTION.

      LOGICAL, intent(out) ::    DO_SSCORR_TRUNCATION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL, intent(out) ::    DO_PLANE_PARALLEL

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL, intent(out) ::    DO_BRDF_SURFACE

!  mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL, intent(out) ::    DO_ADDITIONAL_MVOUT

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL, intent(out) ::    DO_MVOUT_ONLY

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL, intent(out) ::    DO_CHAPMAN_FUNCTION

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL, intent(out) ::    DO_REFRACTIVE_GEOMETRY

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL, intent(out) ::    DO_DELTAM_SCALING

!  double convergence test flag

      LOGICAL, intent(out) ::    DO_DOUBLE_CONVTEST

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL, intent(out) ::    DO_SOLUTION_SAVING
      LOGICAL, intent(out) ::    DO_BVP_TELESCOPING

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL, intent(out) ::    DO_RAYLEIGH_ONLY
      LOGICAL, intent(out) ::    DO_ISOTROPIC_ONLY

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL, intent(out) ::    DO_NO_AZIMUTH
      LOGICAL, intent(out) ::    DO_ALL_FOURIER


!  2. CONTROL INTEGERS
!  ===================

!  Number of discrete ordinate streams

      INTEGER, intent(out) :: NSTREAMS

!  number of computational layers

      INTEGER, intent(out) :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER, intent(out) :: NFINELAYERS

!  number of Legendre phase function expansion moments

      INTEGER, intent(out) :: NMOMENTS_INPUT

!  number of solar beams to be processed

      INTEGER, intent(out) :: NBEAMS

!  Number of user-defined relative azimuths

      INTEGER, intent(out) :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER, intent(out) :: N_USER_STREAMS

!  Number of User-defined vertical levels for  output

      INTEGER, intent(out) :: N_USER_LEVELS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(KIND=8), intent(out) :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      REAL(KIND=8), intent(out) :: LIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      REAL(KIND=8), intent(out) :: ZENITH_TOLERANCE

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(KIND=8), intent(out) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(KIND=8), intent(out) :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(KIND=8), intent(out) :: GEOMETRY_SPECHEIGHT

!  BOA solar zenith angles (degrees)

      REAL(KIND=8), intent(out) :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(KIND=8), intent(out) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined viewing zenith angles input (degrees) 

      REAL(KIND=8), intent(out) :: USER_ANGLES  (MAX_USER_STREAMS)

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(KIND=8), intent(out) :: USER_LEVELS   (MAX_USER_LEVELS)

!  4. OUTPUT EXCEPTIONS
!  ====================

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out) :: STATUS
      INTEGER      , intent(out) :: NMESSAGES
      CHARACTER*(*), intent(out) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(out) :: ACTIONS (0:MAX_MESSAGES)

!  ============================================================
!  ============================================================

!  SPECIAL NOTE on variable GEOMETRY_SPECHEIGHT

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (wiht or without an exact single-scatter correction) uses the same
!    set of angles all the up the nadir from the bottom of atmosphere.

!     This height is normally set equal to the height at the lowest
!     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
!     In this case, no adjustment to the geometrical inputs is needed
!     for the outgoing sphericity correction.

!     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
!     then an adjustment to the geometrical inputs is needed for the
!     outgoing sphericity correction. This adjustment is internal and
!     the model output will still be given at the geometrical angles
!     as specified by the user, even though these angles may not be the
!     ones at which the calculations were done. This situation will occur
!     when we are given a BOA geometry but we want to make a calculation
!     for a reflecting surface (such as a cloud-top) which is above the
!     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
!     height HEIGHT_GRID(NLAYERS) = cloud-top.

!     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
!     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
!     the calculation will go through without the adjustment. A warning
!     about this incorrect input choice will be sent to LOGFILE.

!  ============================================================
!  ============================================================

!  local variables
!  ---------------

      INTEGER            :: STATUS_SUB, FILUNIT, LEN_STRING
      EXTERNAL           LEN_STRING

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  initialize variables

      CALL LIDORT_INIT_INPUTS_BASIC                                &
          ( DO_SOLAR_SOURCES,                                      & ! output
            DO_UPWELLING,            DO_DNWELLING,                 & ! output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,              & ! output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,           & ! output
            DO_SSFULL,               DO_SSCORR_TRUNCATION,         & ! output
            DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,              & ! output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                & ! output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,       & ! output 
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,           & ! output  
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,            & ! output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,           & ! output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                & ! output
            NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,        & ! output
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! output
            FLUX_FACTOR, LIDORT_ACCURACY,                          & ! output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,  & ! output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_LEVELS )      ! output

!  Open file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  read standard inputs

      CALL LIDORT_READ_INPUTS_BASIC                                &
          ( DO_SOLAR_SOURCES,                                      & ! input/output
            DO_UPWELLING,            DO_DNWELLING,                 & ! input/output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,              & ! input/output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,           & ! input/output
            DO_SSFULL,               DO_SSCORR_TRUNCATION,         & ! input/output
            DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,              & ! input/output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                & ! input/output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,       & ! input/output 
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,           & ! input/output  
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,            & ! input/output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,           & ! input/output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                & ! input/output
            NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,        & ! input/output
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! input/output
            FLUX_FACTOR, LIDORT_ACCURACY,                          & ! input/output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,  & ! input/output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_LEVELS,     & ! input/output
            STATUS_SUB, NMESSAGES, MESSAGES, ACTIONS )               ! Output 

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        CLOSE(FILUNIT)
        RETURN
      ENDIF

!  normal return

      CLOSE(FILUNIT)
      RETURN
      
!  Open file error

300   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right File!!'
      CLOSE(FILUNIT)
      RETURN

!  Finish

END SUBROUTINE LIDORT_INPUTS_BASIC_MASTER  


SUBROUTINE LIDORT_INIT_INPUTS_BASIC                                &
          ( DO_SOLAR_SOURCES,                                      & ! output
            DO_UPWELLING,            DO_DNWELLING,                 & ! output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,              & ! output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,           & ! output
            DO_SSFULL,               DO_SSCORR_TRUNCATION,         & ! output
            DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,              & ! output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                & ! output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,       & ! output 
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,           & ! output  
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,            & ! output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,           & ! output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                & ! output
            NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,        & ! output
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! output
            FLUX_FACTOR, LIDORT_ACCURACY,                          & ! output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,  & ! output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_LEVELS )      ! output

!  Initialises all control inputs for LIDORT
!  ------------------------------------------

      IMPLICIT NONE

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  output Arguments
!  ================

      LOGICAL, intent(out)  ::    DO_SOLAR_SOURCES
      LOGICAL, intent(out)  ::    DO_UPWELLING
      LOGICAL, intent(out)  ::    DO_DNWELLING

      LOGICAL, intent(out)  ::    DO_FULLRAD_MODE
      LOGICAL, intent(out)  ::    DO_USER_STREAMS

      LOGICAL, intent(out)  ::    DO_SSCORR_NADIR
      LOGICAL, intent(out)  ::    DO_SSCORR_OUTGOING
      LOGICAL, intent(out)  ::    DO_SSFULL
      LOGICAL, intent(out)  ::    DO_SSCORR_TRUNCATION

      LOGICAL, intent(out)  ::    DO_PLANE_PARALLEL
      LOGICAL, intent(out)  ::    DO_BRDF_SURFACE

      LOGICAL, intent(out)  ::    DO_ADDITIONAL_MVOUT
      LOGICAL, intent(out)  ::    DO_MVOUT_ONLY

      LOGICAL, intent(out)  ::    DO_CHAPMAN_FUNCTION
      LOGICAL, intent(out)  ::    DO_REFRACTIVE_GEOMETRY

      LOGICAL, intent(out)  ::    DO_DELTAM_SCALING
      LOGICAL, intent(out)  ::    DO_DOUBLE_CONVTEST

      LOGICAL, intent(out)  ::    DO_SOLUTION_SAVING
      LOGICAL, intent(out)  ::    DO_BVP_TELESCOPING

      LOGICAL, intent(out)  ::    DO_RAYLEIGH_ONLY
      LOGICAL, intent(out)  ::    DO_ISOTROPIC_ONLY

      LOGICAL, intent(out)  ::    DO_NO_AZIMUTH
      LOGICAL, intent(out)  ::    DO_ALL_FOURIER

      INTEGER, intent(out)  :: NSTREAMS
      INTEGER, intent(out)  :: NLAYERS
      INTEGER, intent(out)  :: NFINELAYERS
      INTEGER, intent(out)  :: NMOMENTS_INPUT

      INTEGER, intent(out)  :: NBEAMS
      INTEGER, intent(out)  :: N_USER_RELAZMS
      INTEGER, intent(out)  :: N_USER_STREAMS
      INTEGER, intent(out)  :: N_USER_LEVELS

      REAL(KIND=8), intent(out)  :: FLUX_FACTOR
      REAL(KIND=8), intent(out)  :: LIDORT_ACCURACY
      REAL(KIND=8), intent(out)  :: EARTH_RADIUS
      REAL(KIND=8), intent(out)  :: RFINDEX_PARAMETER
      REAL(KIND=8), intent(out)  :: GEOMETRY_SPECHEIGHT

      REAL(KIND=8), intent(out)  :: BEAM_SZAS     (MAXBEAMS)
      REAL(KIND=8), intent(out)  :: USER_RELAZMS  (MAX_USER_RELAZMS)
      REAL(KIND=8), intent(out)  :: USER_ANGLES   (MAX_USER_STREAMS)
      REAL(KIND=8), intent(out)  :: USER_LEVELS   (MAX_USER_LEVELS)

!  Local variables

      INTEGER         I

      DO_SOLAR_SOURCES  = .FALSE.

      DO_UPWELLING      = .FALSE.
      DO_DNWELLING      = .FALSE.

      DO_FULLRAD_MODE   = .FALSE.
      DO_USER_STREAMS   = .FALSE.

      DO_SSCORR_NADIR      = .FALSE.
      DO_SSCORR_OUTGOING   = .FALSE.
      DO_SSCORR_TRUNCATION = .FALSE.
      DO_SSFULL            = .FALSE.

      DO_PLANE_PARALLEL      = .FALSE.
      DO_BRDF_SURFACE        = .FALSE.

      DO_ADDITIONAL_MVOUT = .FALSE.
      DO_MVOUT_ONLY       = .FALSE.

      DO_CHAPMAN_FUNCTION    = .FALSE.
      DO_REFRACTIVE_GEOMETRY = .FALSE.

      DO_DELTAM_SCALING  = .FALSE.
      DO_DOUBLE_CONVTEST = .FALSE.

      DO_RAYLEIGH_ONLY  = .FALSE.
      DO_ISOTROPIC_ONLY = .FALSE.

      DO_SOLUTION_SAVING = .FALSE.
      DO_BVP_TELESCOPING = .FALSE.

      DO_NO_AZIMUTH      = .FALSE.
      DO_ALL_FOURIER     = .FALSE.

!  basic integer inputs

      NSTREAMS = 0
      NLAYERS  = 0
      NFINELAYERS    = 0
      NMOMENTS_INPUT = 0

      NBEAMS   = 0
      N_USER_STREAMS = 0
      N_USER_RELAZMS = 0
      N_USER_LEVELS  = 0

!  basic numbers

      FLUX_FACTOR      = ZERO
      LIDORT_ACCURACY  = ZERO

      EARTH_RADIUS      = ZERO
      RFINDEX_PARAMETER = ZERO
      GEOMETRY_SPECHEIGHT = ZERO

      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO

      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES(I) = ZERO
      ENDDO

      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

      DO I = 1, MAX_USER_LEVELS
        USER_LEVELS(I) = ZERO
      ENDDO

! finish

      RETURN
END SUBROUTINE LIDORT_INIT_INPUTS_BASIC

!

SUBROUTINE LIDORT_READ_INPUTS_BASIC                                &
          ( DO_SOLAR_SOURCES,                                      & ! input/output
            DO_UPWELLING,            DO_DNWELLING,                 & ! input/output
            DO_FULLRAD_MODE,         DO_USER_STREAMS,              & ! input/output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,           & ! input/output
            DO_SSFULL,               DO_SSCORR_TRUNCATION,         & ! input/output
            DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,              & ! input/output
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                & ! input/output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,       & ! input/output 
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,           & ! input/output  
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,            & ! input/output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,           & ! input/output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                & ! input/output
            NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,        & ! input/output
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! input/output
            FLUX_FACTOR, LIDORT_ACCURACY,                          & ! input/output
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,  & ! input/output
            BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_LEVELS,     & ! input/output
            STATUS, NMESSAGES, MESSAGES, ACTIONS )                   ! Output 

!  Read all control inputs for LIDORT

      IMPLICIT NONE

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  OUTPUTS
!  -------

!  1. CONTROL FLAGS
!  ================

!  Basic top-level control
!   -- Thermal to be reintroduced later 2010

      LOGICAL, intent(inout) ::    DO_SOLAR_SOURCES
!      LOGICAL, intent(inout) ::    DO_THERMAL_EMISSION

!  directional control

      LOGICAL, intent(inout) ::    DO_UPWELLING
      LOGICAL, intent(inout) ::    DO_DNWELLING

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL, intent(inout) ::    DO_FULLRAD_MODE

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL, intent(inout) ::    DO_USER_STREAMS


!  single scatter corrections. (Includes direct beam reflection)
!    - NADIR    : Plane-parallel line-of-sight
!    - OUTGOING : Line-of-sight in curved atmosphere

      LOGICAL, intent(inout) ::    DO_SSCORR_NADIR
      LOGICAL, intent(inout) ::    DO_SSCORR_OUTGOING

!  Flag for Full-up single scatter calculation. (No MS field)
!    One of the above two SSCORR flags must be set

      LOGICAL, intent(inout) ::    DO_SSFULL

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations. **** Use with CAUTION.

      LOGICAL, intent(inout) ::    DO_SSCORR_TRUNCATION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL, intent(inout) ::    DO_PLANE_PARALLEL

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL, intent(inout) ::    DO_BRDF_SURFACE

!  mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL, intent(inout) ::    DO_ADDITIONAL_MVOUT

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL, intent(inout) ::    DO_MVOUT_ONLY

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL, intent(inout) ::    DO_CHAPMAN_FUNCTION

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL, intent(inout) ::    DO_REFRACTIVE_GEOMETRY

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL, intent(inout) ::    DO_DELTAM_SCALING

!  double convergence test flag

      LOGICAL, intent(inout) ::    DO_DOUBLE_CONVTEST

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL, intent(inout) ::    DO_SOLUTION_SAVING
      LOGICAL, intent(inout) ::    DO_BVP_TELESCOPING

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL, intent(inout) ::    DO_RAYLEIGH_ONLY
      LOGICAL, intent(inout) ::    DO_ISOTROPIC_ONLY

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL, intent(inout) ::    DO_NO_AZIMUTH
      LOGICAL, intent(inout) ::    DO_ALL_FOURIER

!  2. CONTROL INTEGERS
!  ===================

!  Number of discrete ordinate streams

      INTEGER, intent(inout) :: NSTREAMS

!  number of computational layers

      INTEGER, intent(inout) :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER, intent(inout) :: NFINELAYERS

!  number of Legendre phase function expansion moments

      INTEGER, intent(inout) :: NMOMENTS_INPUT

!  number of solar beams to be processed

      INTEGER, intent(inout) :: NBEAMS

!  Number of user-defined relative azimuths

      INTEGER, intent(inout) :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER, intent(inout) :: N_USER_STREAMS

!  Number of User-defined vertical levels for  output

      INTEGER, intent(inout) :: N_USER_LEVELS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(KIND=8), intent(inout) :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      REAL(KIND=8), intent(inout) :: LIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      REAL(KIND=8), intent(inout) :: ZENITH_TOLERANCE

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(KIND=8), intent(inout) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(KIND=8), intent(inout) :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(KIND=8), intent(inout) :: GEOMETRY_SPECHEIGHT

!  BOA solar zenith angles (degrees)

      REAL(KIND=8), intent(inout) :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(KIND=8), intent(inout) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined viewing zenith angles input (degrees) 

      REAL(KIND=8), intent(inout) :: USER_ANGLES  (MAX_USER_STREAMS)

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(KIND=8), intent(inout) :: USER_LEVELS   (MAX_USER_LEVELS)

!  4. OUTPUT EXCEPTIONS
!  ====================

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out) :: STATUS
      INTEGER      , intent(out) :: NMESSAGES
      CHARACTER*(*), intent(out) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(out) :: ACTIONS (0:MAX_MESSAGES)

!  ============================================================
!  ============================================================

!  SPECIAL NOTE on variable GEOMETRY_SPECHEIGHT

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (wiht or without an exact single-scatter correction) uses the same
!    set of angles all the up the nadir from the bottom of atmosphere.

!     This height is normally set equal to the height at the lowest
!     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
!     In this case, no adjustment to the geometrical inputs is needed
!     for the outgoing sphericity correction.

!     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
!     then an adjustment to the geometrical inputs is needed for the
!     outgoing sphericity correction. This adjustment is internal and
!     the model output will still be given at the geometrical angles
!     as specified by the user, even though these angles may not be the
!     ones at which the calculations were done. This situation will occur
!     when we are given a BOA geometry but we want to make a calculation
!     for a reflecting surface (such as a cloud-top) which is above the
!     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
!     height HEIGHT_GRID(NLAYERS) = cloud-top.

!     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
!     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
!     the calculation will go through without the adjustment. A warning
!     about this incorrect input choice will be sent to LOGFILE.

!  ============================================================
!  ============================================================

!  local variables

      CHARACTER*8, parameter     ::  PREFIX = 'LIDORT -'

      LOGICAL         :: ERROR
      CHARACTER * 80  :: PAR_STR
      LOGICAL         :: GFINDPAR
      INTEGER         :: I, FILUNIT, LEN_STRING, NM

!  External functions

      EXTERNAL        GFINDPAR
      EXTERNAL        LEN_STRING

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

!  These are already initialized in calling routine
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
!      ACTIONS(0)      = 'No Action required for this Task'

      ERROR  = .FALSE.
      NM     = 0

!  file unit

      FILUNIT = LIDORT_INUNIT

!  1. READ ALL THE CONTROL FlAGS (BOOLEAN INPUTS)
!  ==============================================

!  operation modes
!  ---------------

!  Basic control

      PAR_STR = 'Include solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  directional output control

      PAR_STR = 'Upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  full Radiance calculation

      PAR_STR = 'Do full Radiance calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  user-defined viewing zenith angle

      PAR_STR = 'User-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
         READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  SS Corrections
!  --------------

!  SSCORR_NADIR

      PAR_STR = 'Do nadir single scatter correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_NADIR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  SSCORR_OUTGOING

      PAR_STR = 'Do outgoing single scatter correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_OUTGOING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Full-up single scatter calculation

      PAR_STR = 'Do Full-up single scatter calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSFULL
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Additional deltam scaling for either single-scatter corrections

      IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING.OR.DO_SSFULL )  THEN
        PAR_STR = 'Additional Delta-M scaling for SS correction?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_TRUNCATION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  other major control
!  -------------------

!  Pseudo-spherical control

      IF ( DO_SOLAR_SOURCES ) THEN
        PAR_STR = 'Plane-parallel treatment of direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  BRDF surface control

      PAR_STR = 'Bi-directional surface reflection?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
        READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )


!  THIS TO BE REINTRODUCED *****************************************
!  multiple scatter source function output control
!      PAR_STR = 'Output multiple scatter layer source functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) SAVE_LAYER_MSST
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  mean-value output control
!  -------------------------

      PAR_STR = 'Generate mean value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Generate only mean value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solar beam control
!  ------------------

!  Other control options for the solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  internal Chapman function calculation

       PAR_STR = 'Perform internal Chapman function calculation?'
       IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_CHAPMAN_FUNCTION
       CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  refractive atmosphere

       PAR_STR = 'Beam path with refractive atmosphere?'
       IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
          READ (FILUNIT,*,ERR=998) DO_REFRACTIVE_GEOMETRY
       CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End control

      ENDIF

!  performance control
!  -------------------

!  Delta-M scaling
!    Should only be set for solar beam sources

      IF ( DO_SOLAR_SOURCES ) THEN
       PAR_STR = 'Include Delta-M scaling?'
       IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
       CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  double convergence test

      PAR_STR = 'Perform double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solution saving mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Solution saving mode?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_SOLUTION_SAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Boundary value problem (BVP) telescope mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Boundary value problem telescoping mode?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_BVP_TELESCOPING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  numerical control (azimuth series)
!  ----------------------------------

!  scatterers and phase function control

      PAR_STR='Rayleigh atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_RAYLEIGH_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR='Isotropic atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  No azimuth dependence (TRUE means Fourier m = 0 only )

      PAR_STR = 'No azimuth dependence in the solution?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All possible Fourier components (2N-1). Debug only

      PAR_STR = 'Compute all Fourier components?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Write control Commented out, 02 June 2010
!  -----------------------------------------

!  output write flags

!      PAR_STR = 'Input control write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,*,ERR=998) DO_WRITE_INPUT
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Input scenario write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,*,ERR=998) DO_WRITE_SCENARIO
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Fourier component output write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,*,ERR=998) DO_WRITE_FOURIER
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Results write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  output filenames

!      PAR_STR = 'filename for input write'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,'(a)',ERR=998) INPUT_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'filename for scenario write'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,'(a)',ERR=998) SCENARIO_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Fourier output filename'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,'(a)',ERR=998) FOURIER_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'filename for main output'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
!           READ (FILUNIT,'(a)',ERR=998) RESULTS_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2. READ ALL THE CONTROL NUMBERS (INTEGER INPUTS)
!  ================================================

!  streams/layers/finelayers/moments (INTEGER input)
!  -------------------------------------------------

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Number of atmospheric layers'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) NLAYERS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_SSCORR_OUTGOING  .or. DO_SSFULL ) THEN
        PAR_STR = 'Number of fine layers (outgoing sphericity correction only)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
             READ (FILUNIT,*,ERR=998) NFINELAYERS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

      PAR_STR = 'Number of input Legendre phase function coefficients'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) NMOMENTS_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions
!  ------------------------------------------------------

!     New Exception handling code, 18 May 2010

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of half-space streams" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension in LIDORT.PARS'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of atmospheric layers" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension in LIDORT.PARS'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

      IF ( DO_SSCORR_OUTGOING .or. DO_SSFULL ) THEN
        IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of fine layers..." > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in LIDORT.PARS'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF
      ENDIF

      IF ( NMOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of input Legendre phase function coefficients" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension in LIDORT.PARS'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Geometry inputs
!  ---------------

!  1. number of Solar zenith angles
!     ---- check not exceeding dimensioned number

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of solar zenith angles" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension in LIDORT.PARS'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  2. always need some azimuth angles if Azimuth flag set
!      --- Set number of azimuths locally to 1 if flag not set
!     ---- check not exceeding dimensioned number

      IF ( .NOT. DO_NO_AZIMUTH ) THEN
        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of ...relative azimuths" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension in LIDORT.PARS'
          STATUS       = LIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF
      ELSE
        N_USER_RELAZMS = 1
      ENDIF

!  3. User defined viewing zenith angles (should be positive)
!     ---- check not exceeding dimensioned number

      IF ( DO_USER_STREAMS ) THEN
        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of ...viewing zenith angles" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension in LIDORT.PARS'
          STATUS       = LIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF
      ENDIF

!  4. Number of output levels
!     ---- check not exceeding dimensioned number

      PAR_STR = 'Number of user-defined vertical output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
              READ (FILUNIT,*,ERR=998) N_USER_LEVELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of ...output levels" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension in LIDORT.PARS'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  3. READ ALL THE FLOATING-POINT INPUTS
!  =====================================

!  Flux constant

      PAR_STR = 'Solar flux constant'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) FLUX_FACTOR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  accuracy criterion

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) LIDORT_ACCURACY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

! C  Earth radius and RF parameter
!  -- only for Chapman function calculation

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        PAR_STR = 'Earth radius (km)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
             READ (FILUNIT,*,ERR=998) EARTH_RADIUS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          PAR_STR = 'Refractive index parameter'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
               READ (FILUNIT,*,ERR=998) RFINDEX_PARAMETER
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF
      ENDIF

!  Input Geometry specfication height

      IF ( DO_SSCORR_OUTGOING ) THEN
        PAR_STR = 'Input geometry specification height [km]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) GEOMETRY_SPECHEIGHT
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Geometry inputs
!  ---------------

!  BOA solar zenith angle inputs

      PAR_STR = 'BOA solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  always need some azimuth angles if Azimuth flag set

      IF ( .NOT. DO_NO_AZIMUTH ) THEN
        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  User defined stream (viewing) angles (should be positive)

      IF ( DO_USER_STREAMS ) THEN
        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  vertical output levels

      PAR_STR = 'User-defined vertical output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_LEVELS
          READ (FILUNIT,*,ERR=998) USER_LEVELS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  normal return

      NMESSAGES = NM
      RETURN

!  line read error - abort immediately

998   CONTINUE
      NM = NM + 1
      STATUS       = LIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: '//PAR_STR(1:LEN_STRING(PAR_STR))
      ACTIONS(NM)  = 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_READ_INPUTS_BASIC

!

SUBROUTINE LIDORT_CHECK_INPUT_BASIC                                                 &
      ( DO_SOLAR_SOURCES, DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY, DO_SSFULL,           & ! Input
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                              & ! Input
        DO_CHAPMAN_FUNCTION, DO_REFRACTIVE_GEOMETRY, DO_NO_AZIMUTH,                 & ! Input/Output
        DO_USER_STREAMS, DO_DELTAM_SCALING, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,    & ! Input/Output
        DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, & ! Input/Output
        NLAYERS, NFINELAYERS, NSTREAMS,  NMOMENTS_INPUT,                            & ! Input         
        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                      & ! Input
        BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_LEVELS,                          & ! Input
        EARTH_RADIUS, HEIGHT_GRID,                                                  & ! Input
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                                        ! output

!  Check the non-threaded inputs

      IMPLICIT NONE

!  include file

      INCLUDE 'LIDORT.PARS_F90'

!  Module inputs
!  -------------

!  single scatter corrections

      LOGICAL, intent(inout)  ::    DO_SSCORR_NADIR       ! May be re-set with a Warning
      LOGICAL, intent(inout)  ::    DO_SSCORR_OUTGOING    ! May be re-set with a Warning

!  Flag for Full-up single scatter calculation

      LOGICAL, intent(in)     ::    DO_SSFULL

!  Basic top-level control

      LOGICAL, intent(in)     ::    DO_SOLAR_SOURCES

!  particular solution control
!    Removed in stripped down version, March 2008
!      LOGICAL, intent(in)  ::    DO_CLASSICAL_SOLUTION

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(in)     ::    DO_PLANE_PARALLEL
      LOGICAL, intent(inout)  ::    DO_REFRACTIVE_GEOMETRY    ! May be re-set with a Warning
      LOGICAL, intent(inout)  ::    DO_CHAPMAN_FUNCTION       ! May be re-set with a Warning

!  scatterers and phase function control

      LOGICAL, intent(inout)   ::    DO_RAYLEIGH_ONLY    ! May be re-set with a Warning
      LOGICAL, intent(inout)   ::    DO_ISOTROPIC_ONLY   ! May be re-set with a Warning
      LOGICAL, intent(inout)   ::    DO_NO_AZIMUTH       ! May be re-set with a Warning

!  Performance control

      LOGICAL, intent(inout)  ::    DO_DELTAM_SCALING    ! May be re-set with a Warning
      LOGICAL, intent(inout)  ::    DO_SOLUTION_SAVING   ! May be re-set with a Warning
      LOGICAL, intent(inout)  ::    DO_BVP_TELESCOPING   ! May be re-set with a Warning

!  directional control

      LOGICAL, intent(in)     ::    DO_UPWELLING
      LOGICAL, intent(in)     ::    DO_DNWELLING

!  stream angle flag

      LOGICAL, intent(inout)  ::    DO_USER_STREAMS       ! May be re-set with a Warning

!  mean value control

      LOGICAL, intent(inout)  ::    DO_ADDITIONAL_MVOUT    ! May be re-set with a Warning
      LOGICAL, intent(inout)  ::    DO_MVOUT_ONLY          ! May be re-set with a Warning

!  Number of discrete ordinate streams

      INTEGER, intent(in)  :: NSTREAMS

!  number of computational layers

      INTEGER, intent(in)  :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing spherical correction algorithm)

      INTEGER, intent(in)  :: NFINELAYERS

!  number of Legendre phase function expansion moments

      INTEGER, intent(inout)  :: NMOMENTS_INPUT               ! May be re-set with a Warning

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER, intent(in)       :: N_USER_RELAZMS
      REAL(KIND=8), intent(in)  :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined zenith angle input 

      INTEGER, intent(in)       :: N_USER_STREAMS
      REAL(KIND=8), intent(in)  :: USER_ANGLES  (MAX_USER_STREAMS)

!  number of solar beams to be processed

      INTEGER, intent(in)       :: NBEAMS

!  TOA solar zenith angles

      REAL(KIND=8), intent(in)  :: BEAM_SZAS ( MAXBEAMS )

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(KIND=8), intent(in)  :: EARTH_RADIUS

!  height grid

      REAL(KIND=8), intent(in)  :: HEIGHT_GRID  ( 0:MAXLAYERS )

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      INTEGER, intent(in)       :: N_USER_LEVELS
      REAL(KIND=8), intent(in)  :: USER_LEVELS   (MAX_USER_LEVELS)

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  local variables

      REAL(KIND=8)   :: XT
      INTEGER        :: I, N, UTA, NSTART, NALLSTREAMS, NM
      CHARACTER*2    :: C2
      LOGICAL        :: LOOP

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of LIDORT Basic Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM     = NMESSAGES

!  Check top level options, set warnings
!  =====================================

!  Check Solar sources present

      IF ( .NOT.DO_SOLAR_SOURCES ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input:  DO_SOLAR_SOURCES not set'
        ACTIONS(NM)  = 'Abort: must set source flag DO_SOLAR_SOURCES'
        STATUS = LIDORT_SERIOUS
      ENDIF

!  Dimensioning  Basic layer/streams/moments

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams NSTREAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

      IF ( DO_SSCORR_OUTGOING .or. DO_SSFULL ) THEN
        IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of fine layers NFINELAYERS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension'
          STATUS       = LIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF
      ENDIF

      IF ( NMOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Number of Legendre moments NMOMENTS_INPUT > maximum dimension'
        ACTIONS(NM) = 'Re-set input value or increase MAXFINELAYERS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  3. Basic geometry

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles NBEAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

      IF ( DO_USER_STREAMS ) THEN
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_STREAMS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS       = LIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  4. Vertical levels output

      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of user vertical output levels N_USER_LEVELS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  check inputs (both file-read and derived)
!  -----------------------------------------

!  Check Chapman function options

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( .NOT. DO_CHAPMAN_FUNCTION ) THEN
        IF ( DO_PLANE_PARALLEL ) THEN
          NM = NM + 1
          MESSAGES(NM) =  'Chapman Function not set in plane parallel mode'
          ACTIONS(NM)  = 'Warning: Chapman function set internally'
          STATUS       = LIDORT_WARNING
          DO_CHAPMAN_FUNCTION = .TRUE.
        ELSE
          NM = NM + 1
          MESSAGES(NM) = 'Chapman Function not set, pseudo-spherical'
          ACTIONS(NM)  = 'You need to set DELTAU_SLANT_INPUT values directly'
          STATUS       = LIDORT_SERIOUS
        ENDIF
       ENDIF
      ENDIF

!  Check sphericity corrections....Cannot both be turned on
!    --------- New code 31 January 2007

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_SSCORR_NADIR .and. DO_SSCORR_OUTGOING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Cannot have both single scatter corrections on'
          ACTIONS(NM)  = 'Turn off DO_SSCORR_NADIR and/or DO_SSCORR_OUTGOING'
          STATUS       = LIDORT_SERIOUS
        ENDIF     
      ENDIF

!  Check beam mode operation. Warning

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_PLANE_PARALLEL ) THEN
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         NM = NM + 1
         MESSAGES(NM)  = 'Bad input: plane-parallel and refractive flags both set'
         ACTIONS(NM)   = 'Warning: turn off Refraction internally'
         STATUS        = LIDORT_WARNING
         DO_REFRACTIVE_GEOMETRY = .FALSE.
        ENDIF
       ENDIF
      ENDIF

!  check consistency of mean value input control
!  ---------------------------------------------

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: Cannot have both mean-value flags set'
          ACTIONS(NM)   = 'Warning: disable DO_MVOUT_ONLY flag internally'
          STATUS        = LIDORT_WARNING
          DO_MVOUT_ONLY = .FALSE.
        ENDIF
      ENDIF

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( .NOT. DO_NO_AZIMUTH ) THEN
            NM = NM + 1
            MESSAGES(NM)  = 'Bad input: Mean-value option requires NO azimuth'
            ACTIONS(NM)   = 'Warning: DO_NO_AZIMUTH flag was set internally'
            STATUS        = LIDORT_WARNING
            DO_NO_AZIMUTH = .TRUE.
          ENDIF
          IF ( DO_USER_STREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM)    = 'Bad input: Mean-value option needs quadratures only'
            ACTIONS(NM)     = 'Warning: DO_USER_STREAMS flag disabled internally'
            STATUS          = LIDORT_WARNING
            DO_USER_STREAMS = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  Check consistency of BVP_TELESCOPING and SOLUTION_SAVING flags
!  ---Warning. Set solution-saving internally

      IF (DO_BVP_TELESCOPING.AND..NOT.DO_SOLUTION_SAVING) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Bad input: BVP telescoping -> solution saving must be set'
        ACTIONS(NM)  = 'Warning:  Solution saveing was set internally'
        STATUS       = LIDORT_WARNING
        DO_SOLUTION_SAVING = .TRUE.
      ENDIF

!  Check consistency of Rayleigh-only and Isotropic-only cases

      IF ( DO_RAYLEIGH_ONLY .AND. DO_ISOTROPIC_ONLY ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Isotropic_only & Rayleigh-only flags both set'
        ACTIONS(NM)  = 'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
        STATUS = LIDORT_SERIOUS
      ENDIF

!  -----------Note the following in the scalar code -------------------

!  no Delta-M scaling with Rayleigh only
!   ---Warning. Turn off delta-M scaling.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: No delta-M scaling with Rayleigh-only'
          ACTIONS(NM)   = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS        = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  no Delta-M scaling with Isotropic only
!   ---Warning. Turn off delta-M scaling.

      IF ( DO_ISOTROPIC_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: No delta-M scaling with Isotropic-only'
          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS       = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Bad input: no directional input is set'
        ACTIONS(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS       = LIDORT_SERIOUS
      ENDIF

!  check number of input Legendre moments
!  ======================================

!  (general scattering case)

      IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY &
             .AND..NOT. DO_SSFULL ) THEN

        IF ( DO_DELTAM_SCALING ) THEN
          IF ( NMOMENTS_INPUT.LT.2*NSTREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM)   = 'Bad input: Fewer than 2N moments with delta-M'
            ACTIONS(NM)    = 'Warning: Re-set NMOMENTS_INPUT to 2N internally'
            STATUS         = LIDORT_WARNING
            NMOMENTS_INPUT = 2*NSTREAMS
          ENDIF
        ELSE
          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
            IF ( NMOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
              NM = NM + 1
              MESSAGES(NM)   = 'Bad input: Fewer than 2N-1 moments without delta-M'
              ACTIONS(NM)    = 'Warning: Re-set NMOMENTS_INPUT to 2N-1 internally'
              STATUS         = LIDORT_WARNING
              NMOMENTS_INPUT = 2*NSTREAMS - 1
            ENDIF
          ENDIF
        ENDIF

      ELSE

!  Checks for Rayleigh only option
!   All warnings.

        IF ( DO_RAYLEIGH_ONLY ) THEN

          IF ( NMOMENTS_INPUT.NE.2 ) THEN
            NM = NM + 1
            MESSAGES(NM)   = 'Bad input: Rayleigh-only phase momemts not = 2'
            ACTIONS(NM)    = 'Warning: Set NMOMENTS_INPUT = 2 internally'
            STATUS         = LIDORT_WARNING
            NMOMENTS_INPUT = 2
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            NM = NM + 1
            MESSAGES(NM)  = 'Bad input: Bvp telescoping not possible, Rayleigh only'
            ACTIONS(NM)   = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS        = LIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Solution saving not possible, Rayleigh only'
            ACTIONS(NM)  = 'Warning: Turn off SOLUTION_SAVING internally'
            STATUS       = LIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
          ENDIF

        ENDIF

!  Checks for Isotropic only option
!   All warning messages

        IF ( DO_ISOTROPIC_ONLY ) THEN

          IF ( NMOMENTS_INPUT.NE.0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: No phase function moments for isotropic-only'
            ACTIONS(NM)  = 'Warning: NMOMENTS_INPUT = 0 was set internally'
            STATUS       = LIDORT_WARNING
            NMOMENTS_INPUT = 0
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Bvp telescoping not possible, Isotropic only'
            ACTIONS(NM)  = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS       = LIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Solution saving not possible, Isotropic only'
            ACTIONS(NM)  = 'Warning" Turn off SOLUTION_SAVING internally'
            STATUS       = LIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
          ENDIF

        ENDIF

      ENDIF

!  reset solution saving and bvp telescoping flags

      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND. &
          ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND. &
          ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

!  BVP telescoping doesn't work with non-Lambertian surfaces
!   Reason, not yet coded. Theory already worked out.
!    ---WARNING. BVP telescoping Flag turned off
!     ---Not required in the basic setup

!      IF (  DO_BVP_TELESCOPING ) THEN
!        IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
!          NM = NM + 1
!          MESSAGES(NM) =  'BVP telescoping must be disabled, non-Lambertian'
!          ACTIONS(NM) = 'Warning: DO_BVP_TELESCOPING turned off internally'
!          STATUS = LIDORT_WARNING
!          DO_BVP_TELESCOPING = .FALSE.
!          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
!        ENDIF
!      ENDIF
   
!  Check azimuth-only conditions
!  =============================

!  Check no-Azimuth flag
!    ---WARNING. Do-no-Azimuth Flag turned on

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_USER_STREAMS .AND. N_USER_STREAMS.EQ.1 ) THEN
          IF ( USER_ANGLES(1) .EQ. ZERO ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: zenith-sky output requires no azimuth'
            ACTIONS(NM)  = 'Warning: DO_NO_AZIMUTH flag set true internally'
            STATUS       = LIDORT_WARNING
          ENDIF
        ENDIF
      ENDIF

!  Checks for Isotropic only option
!    ---WARNING. Do-no-Azimuth Flag turned on

      IF ( DO_ISOTROPIC_ONLY ) THEN
        IF ( .NOT.DO_NO_AZIMUTH ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: no azimuth dependence for isotropic_only'
          ACTIONS(NM)   = 'Warning: DO_NO_AZIMUTH turned on internally'
          STATUS        = LIDORT_WARNING
          DO_NO_AZIMUTH = .TRUE.
        ENDIF
      ENDIF

!  Check single scattering correction and Do no Azimuth
!    ---WARNING. Do-no-Azimuth Flag turned off

      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        IF ( DO_NO_AZIMUTH ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: need azimuth dependence for SS corrections'
          ACTIONS(NM)   = 'Warning: DO_NO_AZIMUTH turned off internally'
          STATUS        = LIDORT_WARNING
          DO_NO_AZIMUTH = .FALSE.
        ENDIF
      ENDIF

!  Check: OLD single scattering correction and Do Rayleigh
!    ---WARNING. SS Flag turned off
!  Check only required for the diffuse field calculations (Version 2.3)

      IF ( DO_SSCORR_NADIR ) THEN
        IF ( DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
          NM = NM + 1
          MESSAGES(NM)    = 'Bad input: No SS correction for Rayleigh only'
          ACTIONS(NM)     = 'Warning: DO_SSCORR_NADIR turned off internally'
          STATUS          = LIDORT_WARNING
          DO_SSCORR_NADIR = .FALSE.
        ENDIF
      ENDIF

!  Full-up single scatter, enabled 25 September 2007.
!   Single scatter corrections must be turned on

      IF ( DO_SSFULL ) THEN
        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN  
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Full SS, must have one SSCORR flag set'
          ACTIONS(NM)  = 'Full SS: default to use outgoing SS correction'
          STATUS       = LIDORT_WARNING
          DO_SSCORR_NADIR    = .FALSE.
          DO_SSCORR_OUTGOING = .TRUE.
        ENDIF
      ENDIF

!  Full-up single scatter, enabled 25 September 2007.
!   Diffuse-field Delta-M scaling must be turned off

      IF ( DO_SSFULL ) THEN
        IF ( DO_DELTAM_SCALING ) THEN        
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Full SS, diffuse-field delta-M on'
          ACTIONS(NM)  = 'Full SS: default to deltam_scaling = false'
          STATUS       = LIDORT_WARNING
          DO_DELTAM_SCALING   = .FALSE.
        ENDIF
      ENDIF

!  check viewing geometry input
!  ============================

!  Check earth radius (Chapman function only)
!    ---WARNING. Default value of 6371.0 will be set

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          IF ( EARTH_RADIUS.LT.6320.0D0 .OR.  &
               EARTH_RADIUS.GT.6420.0D0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Earth radius outside of [6320-6420]'
            ACTIONS(NM)  = 'Re-set value'
            STATUS       = LIDORT_SERIOUS
          ENDIF
        ENDIF
      ENDIF

!  Check dimensioning on Legendre numbers (refractive geometry only)

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        NALLSTREAMS = NBEAMS*NLAYERS + NSTREAMS + N_USER_STREAMS
        IF ( NALLSTREAMS .GT. MAX_ALLSTRMS_P1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Dimensioning error for refractive beam angles'
          ACTIONS(NM)  = 'Increase dimension MAX_ALLSTRMS_P1'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Check solar zenith angle input

      DO I = 1, NBEAMS
        IF ( BEAM_SZAS(I) .LT. ZERO .OR.             &
             BEAM_SZAS(I) .GE. 90.0D0 ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) =  'Bad input: out-of-range beam angle, no. '//C2
          ACTIONS(NM) = 'Look at BEAM_SZAS input, should be < 90 & > 0'
          STATUS = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check zenith tolerance input
!    ---WARNING. Default of 0.001 will be set

!      IF ( ZENITH_TOLERANCE.LE.ZERO .OR.
!           ZENITH_TOLERANCE.GT.0.001 ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Bad input: Zenith tolerance level out of bounds'
!        ACTIONS(NM)  = 'Warning: ZENITH_TOLERANCE set to 0.001 internally'
!        STATUS       = LIDORT_WARNING
!        ZENITH_TOLERANCE = 0.001D0
!      ENDIF

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0D0   .OR. &
             USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTIONS(NM)  = 'Look at azimuth angle input, should be in [0,360]'
          LOOP         = .FALSE.
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_STREAMS ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
          I = I + 1
          IF ( USER_ANGLES(I) .GT. 90.0   .OR. &
               USER_ANGLES(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range user stream, no. '//C2
            ACTIONS(NM)  = 'Look at user-defined angle input'
            LOOP         = .FALSE.
            STATUS       = LIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Check height grid input (Chapman function only)

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.NLAYERS)
          I = I + 1
          IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Height-grid not monotonically decreasing; Layer '//C2
            ACTIONS(NM)  = 'Look at Height-grid input'
            LOOP         = .FALSE.
            STATUS       = LIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Check vertical outputs
!  ----------------------

!  check vertical output levels (should always be within atmosphere!)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_LEVELS)
        I = I + 1
        IF ( USER_LEVELS(I) .GT. DBLE(NLAYERS) .OR. &
            USER_LEVELS(I) .LT. ZERO )  THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Out of range for level choice # '//C2
          ACTIONS(NM)  = 'Re-set level output '
          LOOP         = .FALSE.
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  check repetition of vertical output choices

      UTA = 0
      LOOP = .TRUE.
      DO WHILE ( LOOP .AND. UTA .LT. N_USER_LEVELS )
        UTA = UTA + 1
        XT = USER_LEVELS(UTA)
        NSTART = 0
        DO N = 1, N_USER_LEVELS
          IF ( XT .EQ. USER_LEVELS(N)) NSTART = NSTART + 1
        ENDDO
        IF ( NSTART .NE. 1 ) THEN
          NM = NM + 1
          LOOP         = .FALSE.
          MESSAGES(NM) = 'Bad input: repetition of vertical output choice'
          ACTIONS(NM)  = 'Re-set level output '
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_CHECK_INPUT_BASIC

!

SUBROUTINE LIDORT_CHECK_INPUT_THREAD                                &
          ( NLAYERS, THREAD, LAMBERTIAN_ALBEDO, DELTAU_VERT_INPUT,  & ! Input
           OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                 & ! Input
           STATUS, NMESSAGES, MESSAGES, ACTIONS )                     ! output

!  Check the threaded optical property inputs

      IMPLICIT NONE

!  include file

      INCLUDE 'LIDORT.PARS_F90'

!  Module input
!  ------------

      INTEGER, intent(in)  :: NLAYERS, THREAD

!  multilayer optical property (bulk) inputs

      REAL(KIND=8), intent(in)  :: OMEGA_TOTAL_INPUT  ( MAXLAYERS, MAXTHREADS )
      REAL(KIND=8), intent(in)  :: DELTAU_VERT_INPUT  ( MAXLAYERS, MAXTHREADS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(KIND=8), intent(in)  ::  PHASMOMS_TOTAL_INPUT &
                  ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Lambertian Surface control

      REAL(KIND=8), intent(in)  :: LAMBERTIAN_ALBEDO (MAXTHREADS)

!  Module output
!  -------------

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  local variables

      INTEGER          :: L, NM
      CHARACTER*3      :: C3, WTHREAD

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of LIDORT Threaded Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM     = NMESSAGES

!  thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

!  check Thread-dependent optical property inputs
!  ----------------------------------------------

!  make sure the Lambertian surface is in range

      IF ( LAMBERTIAN_ALBEDO(THREAD) .LT. ZERO .OR. &
          LAMBERTIAN_ALBEDO(THREAD) .GT. ONE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Lambertian albedo not in range [0,1]'
        ACTIONS(NM)  = 'Check albedo input, thread # '//wthread
        STATUS       = LIDORT_SERIOUS
      ENDIF

!  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT_INPUT(L,THREAD).LE.ZERO ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTIONS(NM)  = 'Check optical thickness input, thread # '//wthread
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  check single scatter albedos

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L,THREAD).GT.ONE-OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input, thread # '//wthread
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  solar beam, cannot be too small

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L,THREAD).LT.OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input, thread # '//wthread
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO
      
!  Check first phase function moments

      DO L = 1, NLAYERS
        IF ( PHASMOMS_TOTAL_INPUT(0,L,THREAD).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'First phase moment not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First phase moment input, thread # '//wthread
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_CHECK_INPUT_THREAD

!

      SUBROUTINE LIDORT_DERIVE_INPUT_BASIC                           &
      ( DO_FULLRAD_MODE, DO_USER_STREAMS,                            & ! Input
        DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                         & ! Input
        DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,              & ! Input
        DO_UPWELLING, DO_DNWELLING, DO_REFRACTIVE_GEOMETRY,          & ! Input
        DO_DOUBLE_CONVTEST, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,  & ! Input or Output
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,                          & ! Input or Output
        NLAYERS, NSTREAMS, NBEAMS, NMOMENTS_INPUT,                   & ! Input
        N_USER_STREAMS, N_USER_RELAZMS, USER_ANGLES_INPUT,           & ! Input
        N_USER_LEVELS,  USER_LEVELS, BEAM_SZAS, SZA_LOCAL_INPUT,     & ! Input
        DO_MSMODE_LIDORT, NMOMENTS, BEAM_COSINES, SUNLAYER_COSINES,  & ! Output
        N_CONVTESTS, N_DIRECTIONS, WHICH_DIRECTIONS,                 & ! Output
        NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_PARTLAYERS,      & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                    & ! Output
        USER_ANGLES, USER_STREAMS, USER_SECANTS,                     & ! Output
        PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX,PARTLAYERS_LAYERIDX,  & ! Output
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_VALUES,   & ! Output
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN )                       ! Output


!  This is the basic bookkeeping routine

!  Continuation lines :
!     I  is INTENT(in)
!     M  is INTENT(inout)
!     O  is INTENT(out)

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  List of Input arguments
!  -----------------------

!  Full Radiance  calculation

      LOGICAL, intent(in)  ::    DO_FULLRAD_MODE

!  single scatter and direct beam corrections

      LOGICAL, intent(in)  ::    DO_SSCORR_NADIR
      LOGICAL, intent(in)  ::    DO_SSCORR_OUTGOING

!  Flag for Full-up single scatter calculation

      LOGICAL, intent(in)  ::    DO_SSFULL

!  stream angle flag

      LOGICAL, intent(in)  ::    DO_USER_STREAMS

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(in)  ::    DO_REFRACTIVE_GEOMETRY

!  scatterers and phase function control

      LOGICAL, intent(in)  ::    DO_RAYLEIGH_ONLY
      LOGICAL, intent(in)  ::    DO_ISOTROPIC_ONLY

!  directional control

      LOGICAL, intent(in)  ::    DO_UPWELLING
      LOGICAL, intent(in)  ::    DO_DNWELLING

!  Number of discrete ordinate streams

      INTEGER, intent(in)  :: NSTREAMS

!  number of computational layers

      INTEGER, intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER, intent(in)  :: NBEAMS

!  number of Legendre phase function expansion moments

      INTEGER, intent(in)  :: NMOMENTS_INPUT

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      REAL(KIND=8), intent(in)  :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAXBEAMS )

!  Input solar zenith angles

      REAL(KIND=8), intent(in)  :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER     , intent(in)  :: N_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER     , intent(in)  :: N_USER_STREAMS
      REAL(KIND=8), intent(in)  :: USER_ANGLES_INPUT  (MAX_USER_STREAMS)

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      INTEGER     , intent(in)  :: N_USER_LEVELS
      REAL(KIND=8), intent(in)  :: USER_LEVELS   (MAX_USER_LEVELS)

!  Input modified arguments
!  ------------------------

!               Intent(inout) to this routine

!  double convergence test flag

      LOGICAL, intent(inout)  ::    DO_DOUBLE_CONVTEST

!  Performance enhancement

      LOGICAL, intent(inout)  ::    DO_SOLUTION_SAVING
      LOGICAL, intent(inout)  ::    DO_BVP_TELESCOPING

!  mean value control

      LOGICAL, intent(inout)  ::    DO_ADDITIONAL_MVOUT
      LOGICAL, intent(inout)  ::    DO_MVOUT_ONLY

!  Pure output arguments
!  ---------------------

!  Mode of operation

      LOGICAL, intent(out)  :: DO_MSMODE_LIDORT

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER, intent(out)  :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER, intent(out)  :: NSTREAMS_2
      INTEGER, intent(out)  :: NTOTAL
      INTEGER, intent(out)  :: N_SUBDIAG, N_SUPDIAG

!  Number of directions (1 or 2) and directional array

      INTEGER, intent(out)  :: N_DIRECTIONS
      INTEGER, intent(out)  :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Local solar zenith angles Cosines 

      REAL(KIND=8), intent(out)  :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)
      REAL(KIND=8), intent(out)  :: BEAM_COSINES(MAXBEAMS)

!  Quadrature weights and abscissae, and product

      REAL(KIND=8), intent(out)  :: QUAD_STREAMS (MAXSTREAMS)
      REAL(KIND=8), intent(out)  :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(KIND=8), intent(out)  :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(KIND=8), intent(out)  :: USER_ANGLES  (MAX_USER_STREAMS)
      REAL(KIND=8), intent(out)  :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(KIND=8), intent(out)  :: USER_SECANTS  (MAX_USER_STREAMS)

!  output optical depth masks and indices

      LOGICAL, intent(out)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(out)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(out)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(out)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  off-grid optical depths (values, masks, indices)

      INTEGER, intent(out)      :: N_PARTLAYERS
      INTEGER, intent(out)      :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)
      REAL(KIND=8), intent(out) :: PARTLAYERS_VALUES       (MAX_PARTLAYERS)

!  Number of convergences

      INTEGER, intent(out)  :: N_CONVTESTS

!  Layer masks for doing integrated source terms

      LOGICAL, intent(out)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(out)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  local variables
!  ---------------

      INTEGER               :: INDEX_ANGLES ( MAX_USER_STREAMS )
      REAL(KIND=8)          :: ALL_ANGLES ( MAX_USER_STREAMS ), DT, RT
      REAL(KIND=8)          :: MU1, MU2
      INTEGER               :: I, UT, N, UTA, NSTART
      INTEGER               :: N_ALLLAYERS_UP, N_ALLLAYERS_DN

!  set additional numbers (derived input)
!  ======================

!  Mode of operation

      DO_MSMODE_LIDORT = .FALSE.
      IF ( DO_FULLRAD_MODE ) THEN
        IF ( .NOT. DO_SSFULL ) THEN 
          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
            DO_MSMODE_LIDORT = .TRUE.
          ENDIF
        ENDIF
      ELSE
        IF ( .NOT. DO_SSFULL ) THEN
          DO_MSMODE_LIDORT = .TRUE.
        ENDIF
      ENDIF

!  SSFULL flag cancels some other flags

      IF ( DO_SSFULL ) THEN
        DO_DOUBLE_CONVTEST = .FALSE.
        DO_SOLUTION_SAVING = .FALSE.
        DO_BVP_TELESCOPING = .FALSE.
        DO_ADDITIONAL_MVOUT = .FALSE.
        DO_MVOUT_ONLY       = .FALSE.
      ENDIF

!  Directional indices

      IF ( DO_UPWELLING .AND. DO_DNWELLING ) THEN
        N_DIRECTIONS = 2
        WHICH_DIRECTIONS(1) = UPIDX
        WHICH_DIRECTIONS(2) = DNIDX
      ELSE
        N_DIRECTIONS = 1
        WHICH_DIRECTIONS(2) = 0
        IF ( DO_UPWELLING ) THEN
          WHICH_DIRECTIONS(1) = UPIDX
        ELSE IF ( DO_DNWELLING) THEN
          WHICH_DIRECTIONS(1) = DNIDX
        ENDIF
      ENDIF

!  Number of moments

      IF ( DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = 2
      ENDIF
      IF ( DO_ISOTROPIC_ONLY ) THEN
        NMOMENTS = 0
      ENDIF
      IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY ) THEN
        NMOMENTS = MIN ( 2 * NSTREAMS - 1, NMOMENTS_INPUT )
      ENDIF

!  total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS

!  Set Quadrature abscissae and weights

      CALL GAULEG(ZERO,ONE,QUAD_STREAMS,QUAD_WEIGHTS,NSTREAMS)

!  set auxiliary quantities

      DO I = 1, NSTREAMS
        QUAD_STRMWTS(I) = QUAD_STREAMS(I)*QUAD_WEIGHTS(I)
      ENDDO

!  size of boundary value problem matrices and vectors

      NTOTAL = NLAYERS*2*NSTREAMS

!  number of sub and super diagonals in band matrix (boundary value problem)

      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTREAMS - 1
        N_SUPDIAG = 2*NSTREAMS - 1
      ELSE
        N_SUBDIAG = 3*NSTREAMS - 1
        N_SUPDIAG = 3*NSTREAMS - 1
      ENDIF

!  Set average cosines in the refractive geometry case

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO I = 1, NBEAMS
          MU1 = DCOS(SZA_LOCAL_INPUT(0,I)*DEG_TO_RAD)
          DO N = 1, NLAYERS
            MU2 = DCOS(SZA_LOCAL_INPUT(N,I)*DEG_TO_RAD)
            SUNLAYER_COSINES(N,I) = HALF * ( MU1 + MU2 )
            MU1 = MU2
          ENDDO
        ENDDO
      ENDIF

!  Set cosines in the non-refractive geometry case (same all layers)

      IF ( .not. DO_REFRACTIVE_GEOMETRY ) THEN
        DO I = 1, NBEAMS
          MU1 = DCOS(BEAM_SZAS(I)*DEG_TO_RAD)
          BEAM_COSINES(I) = DCOS(BEAM_SZAS(I)*DEG_TO_RAD)
        ENDDO
      ENDIF

!  Set the angle masks
!  ===================

!  initialize
!  ----------

!  Rank the output angles

      IF ( DO_USER_STREAMS ) THEN
        IF ( N_USER_STREAMS .NE. 1 ) THEN
          DO I = 1, N_USER_STREAMS
            ALL_ANGLES(I) = USER_ANGLES_INPUT(I)
          ENDDO
          CALL INDEXX ( N_USER_STREAMS, ALL_ANGLES, INDEX_ANGLES )
          DO I = 1, N_USER_STREAMS
            USER_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
          ENDDO
        ELSE
          USER_ANGLES(1) = USER_ANGLES_INPUT(1)
        ENDIF
      ENDIF

!  User stream cosines and secants

      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_ANGLES(I))
          USER_SECANTS(I) = ONE / USER_STREAMS(I)
        ENDDO
      ENDIF

!  number of tests to be applied for convergence

      N_CONVTESTS = N_USER_RELAZMS * N_USER_STREAMS * N_DIRECTIONS
      N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS 

!  Sort out User vertical level outputs
!  ------------------------------------

!  Sort in ascending order

      IF ( N_USER_LEVELS .GT. 1 ) THEN
        CALL HPSORT(N_USER_LEVELS,USER_LEVELS)
      ENDIF

!  mark all output levels not equal to layer boundary values

      NSTART = 0
      UT = 0
      DO UTA = 1, N_USER_LEVELS
        DT = USER_LEVELS(UTA)
        RT = DT - DBLE(INT(DT))
        N = INT(DT) + 1
        IF ( RT.GT.ZERO) THEN
          UT = UT + 1
          PARTLAYERS_OUTFLAG(UTA)  = .TRUE.
          PARTLAYERS_OUTINDEX(UTA) = UT
          PARTLAYERS_LAYERIDX(UT)  = N
          UTAU_LEVEL_MASK_UP(UTA) = N
          UTAU_LEVEL_MASK_DN(UTA) = N - 1
          PARTLAYERS_VALUES(UT)    = RT
        ELSE
          PARTLAYERS_OUTFLAG(UTA)  = .FALSE.
          PARTLAYERS_OUTINDEX(UTA) =   0
          UTAU_LEVEL_MASK_UP(UTA) = N - 1
          UTAU_LEVEL_MASK_DN(UTA) = N - 1
        ENDIF
      ENDDO
      N_PARTLAYERS = UT

!  Set masking and number of layer source terms
!  --------------------------------------------

!   .. for upwelling

      IF ( DO_UPWELLING ) THEN
        DO N = 1, NLAYERS
          STERM_LAYERMASK_UP(N) = .FALSE.
        ENDDO
        UTA = 1
        UT  = 1
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_ALLLAYERS_UP = UTAU_LEVEL_MASK_UP(UTA) + 1
        ELSE
          N_ALLLAYERS_UP = PARTLAYERS_LAYERIDX(UT)
        ENDIF
        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
      ENDIF

!   .. for downwelling

      IF ( DO_DNWELLING ) THEN
        DO N = 1, NLAYERS
          STERM_LAYERMASK_DN(N) = .FALSE.
        ENDDO
        UTA = N_USER_LEVELS
        UT  = N_PARTLAYERS
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_ALLLAYERS_DN = UTAU_LEVEL_MASK_DN(UTA)
        ELSE
          N_ALLLAYERS_DN = PARTLAYERS_LAYERIDX(UT)
        ENDIF
        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_DERIVE_INPUT_BASIC

