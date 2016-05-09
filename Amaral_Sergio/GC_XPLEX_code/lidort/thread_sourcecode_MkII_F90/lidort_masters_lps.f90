!$Id: lidort_masters_lps.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LPS_MASTER (top-level master)             #
! #            LIDORT_LPS_FOURIER_MASTER                        #
! #                                                             #
! ###############################################################

SUBROUTINE LIDORT_LPS_MASTER                                            & 
          ( THREAD,                  DO_SOLAR_SOURCES,                  & ! Input
            DO_UPWELLING,            DO_DNWELLING,                      & ! Input
            DO_FULLRAD_MODE,         DO_USER_STREAMS,                   & ! Input/Output
            DO_SSCORR_NADIR,         DO_SSCORR_OUTGOING,                & ! Input/Output
            DO_SSFULL,               DO_SSCORR_TRUNCATION,              & ! Input/Output
            DO_PLANE_PARALLEL,       DO_BRDF_SURFACE,                   & ! Input
            DO_ADDITIONAL_MVOUT,     DO_MVOUT_ONLY,                     & ! Input/Output
            DO_CHAPMAN_FUNCTION,     DO_REFRACTIVE_GEOMETRY,            & ! Input/Output
            DO_DELTAM_SCALING,       DO_DOUBLE_CONVTEST,                & ! Input/Output  
            DO_RAYLEIGH_ONLY,        DO_ISOTROPIC_ONLY,                 & ! Input/Output
            DO_SOLUTION_SAVING,      DO_BVP_TELESCOPING,                & ! Input/Output
            DO_ALL_FOURIER,          DO_NO_AZIMUTH,                     & ! Input/Output
            DO_PROFILE_LINEARIZATION,DO_SURFACE_LINEARIZATION,          & ! Input
            NSTREAMS, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,             & ! Input (Nmoment_input re-set)
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,      & ! Input
            LAYER_VARY_FLAG,  LAYER_VARY_NUMBER,  N_SURFACE_WFS,        & ! Input
            FLUX_FACTOR, LIDORT_ACCURACY,                               & ! Input
            EARTH_RADIUS, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,       & ! Input
            BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_LEVELS,    & ! Input
            FINEGRID, HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,     & ! Input
            DELTAU_VERT_INPUT,    L_DELTAU_VERT_INPUT,                  & ! Input
            OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,                  & ! Input
            PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,               & ! Input
            LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC,     & ! Input
               BRDF_F_0,    BRDF_F,    USER_BRDF_F_0,    USER_BRDF_F,   & ! Input
            LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,   & ! Input
            INTENSITY, MEAN_INTENSITY, FLUX_INTEGRAL,                   & ! Output
            N_GEOMETRIES, FOURIER_SAVED,                                & ! Output
            PROFILEWF, MINT_PROFILEWF, FLUX_PROFILEWF,                  & ! Output
            SURFACEWF, MINT_SURFACEWF, FLUX_SURFACEWF,                  & ! Output
            STATUS_INPUTCHECK,  NCHECKMESSAGES, CHECKMESSAGES, ACTIONS, & ! Output
            STATUS_CALCULATION, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )      ! Output

!  master program for Radiances + Profile Jacobians Calculation

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine input arguments (Intent: (in) or (inout) )
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Thread number

      INTEGER, intent(in)  ::          THREAD

!  1. CONTROL FLAGS
!  ================

!  Basic top-level control
!   -- Thermal to be reintroduced later 2010

      LOGICAL, intent(in) ::    DO_SOLAR_SOURCES
!      LOGICAL, intent(in) ::    DO_THERMAL_EMISSION

!  directional control

      LOGICAL, intent(in) ::    DO_UPWELLING
      LOGICAL, intent(in) ::    DO_DNWELLING

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL, intent(in) ::    DO_FULLRAD_MODE

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL, intent(inout) ::    DO_USER_STREAMS

!  single scatter corrections. (Includes direct beam reflection)
!    - NADIR    : Plane-parallel line-of-sight
!    - OUTGOING : Line-of-sight in curved atmosphere

      LOGICAL, intent(inout) ::    DO_SSCORR_NADIR         ! May be re-set after Checking
      LOGICAL, intent(inout) ::    DO_SSCORR_OUTGOING      ! May be re-set after Checking

!  Flag for Full-up single scatter calculation. (No MS field)
!    One of the above two SSCORR flags must be set

      LOGICAL, intent(in)    ::    DO_SSFULL

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations. **** Use with CAUTION.

      LOGICAL, intent(inout) ::    DO_SSCORR_TRUNCATION         ! May be re-set after Checking

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL, intent(in) ::    DO_PLANE_PARALLEL

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL, intent(in) ::    DO_BRDF_SURFACE

!  mean value control (1). If set --> Flux output AS WELL AS Intensities
   
      LOGICAL, intent(inout) ::    DO_ADDITIONAL_MVOUT     ! May be re-set after Checking

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL, intent(inout) ::    DO_MVOUT_ONLY           ! May be re-set after Checking

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL, intent(inout) ::    DO_CHAPMAN_FUNCTION     ! May be re-set after Checking

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL, intent(inout) ::    DO_REFRACTIVE_GEOMETRY  ! May be re-set after Checking

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL, intent(inout) ::    DO_DELTAM_SCALING       ! May be re-set after Checking

!  double convergence test flag

      LOGICAL, intent(inout) ::    DO_DOUBLE_CONVTEST      ! May be re-set after Checking

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL, intent(inout) ::    DO_SOLUTION_SAVING      ! May be re-set after Checking
      LOGICAL, intent(inout) ::    DO_BVP_TELESCOPING      ! May be re-set after Checking

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL, intent(inout) ::    DO_RAYLEIGH_ONLY        ! May be re-set after Checking
      LOGICAL, intent(inout) ::    DO_ISOTROPIC_ONLY       ! May be re-set after Checking

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL, intent(inout) ::    DO_NO_AZIMUTH           ! May be re-set after Checking
      LOGICAL, intent(inout) ::    DO_ALL_FOURIER

!  Control linearization
 
      LOGICAL, intent(in)  ::  DO_PROFILE_LINEARIZATION
      LOGICAL, intent(in)  ::  DO_SURFACE_LINEARIZATION

!  2. CONTROL INTEGERS
!  ===================

!  Number of discrete ordinate streams

      INTEGER, intent(in) :: NSTREAMS

!  number of computational layers

      INTEGER, intent(in) :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER, intent(in) :: NFINELAYERS

!  number of Legendre phase function expansion moments

      INTEGER, intent(inout) :: NMOMENTS_INPUT          ! May be re-set after Checking

!  number of solar beams to be processed

      INTEGER, intent(in) :: NBEAMS

!  Number of user-defined relative azimuths

      INTEGER, intent(in) :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER, intent(in) :: N_USER_STREAMS

!  Number of User-defined vertical levels for  output

      INTEGER, intent(in) :: N_USER_LEVELS

!  Linearization control
!  ---------------------

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, intent(in)  ::  LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER, intent(in)  ::  LAYER_VARY_NUMBER(MAXLAYERS)

!  Total number of Surface Jacobians

      INTEGER, intent(in)  ::  N_SURFACE_WFS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(KIND=8), intent(in) :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      REAL(KIND=8), intent(in) :: LIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      REAL(KIND=8), intent(in) :: ZENITH_TOLERANCE

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(KIND=8), intent(in) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(KIND=8), intent(in) :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(KIND=8), intent(in) :: GEOMETRY_SPECHEIGHT

!  BOA solar zenith angles (degrees)

      REAL(KIND=8), intent(in) :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(KIND=8), intent(in) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined viewing zenith angles input (degrees) 

      REAL(KIND=8), intent(in) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(KIND=8), intent(in) :: USER_LEVELS   (MAX_USER_LEVELS)

!  4. ATMOSPHERIC INPUTS
!  =====================

!  PTH inputs
!  ----------

!  multilayer Height inputs
!   (only required for the Chapman function calculations)

      REAL(KIND=8), intent(in) :: HEIGHT_GRID ( 0:MAXLAYERS )

!  multilayer atmospheric inputs (P andT)
!   (only required for the Chapman function calculations, refractive geometry)

      REAL(KIND=8), intent(in) :: PRESSURE_GRID    ( 0:MAXLAYERS )
      REAL(KIND=8), intent(in) :: TEMPERATURE_GRID ( 0:MAXLAYERS )

!  Number of fine layer gradations 
!    (only for Chapman function calculations with refractive geometry)

      INTEGER, intent(in)      :: FINEGRID ( MAXLAYERS )

!  optical property inputs
!  -----------------------

!  multilayer optical property (bulk) inputs

      REAL(KIND=8), intent(in) :: OMEGA_TOTAL_INPUT  ( MAXLAYERS, MAXTHREADS )
      REAL(KIND=8), intent(in) :: DELTAU_VERT_INPUT  ( MAXLAYERS, MAXTHREADS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(KIND=8), intent(in) :: PHASMOMS_TOTAL_INPUT &
           ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(KIND=8), intent(in) :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(KIND=8), intent(in) :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(KIND=8), intent(in) :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS,MAXTHREADS)

!  5. Surface inputs (New BRDF arrays, 23 March 2010)
!  ==================================================

!  Lambertian Surface control

      REAL(KIND=8), intent(in) :: LAMBERTIAN_ALBEDO (MAXTHREADS)

!  Exact (direct bounce) BRDF (same all threads)
     
      REAL(KIND=8), intent(in) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized Exact (direct bounce) BRDF (same all threads)
     
      REAL(KIND=8), intent(in) :: LS_EXACTDB_BRDFUNC &
            ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(KIND=8), intent(in) :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: BRDF_F        ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(KIND=8), intent(in) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Fourier components of BRDF, in the following order (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(KIND=8), intent(in) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(KIND=8), intent(in) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Subroutine output arguments (Intent(out))
!  +++++++++++++++++++++++++++++++++++++++++

!  Main results
!  ------------

!  Intensity Results at all angles and optical depths

      REAL(KIND=8), intent(out) :: INTENSITY ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

!  Results for mean-value output

      REAL(KIND=8), intent(out) :: MEAN_INTENSITY (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
      REAL(KIND=8), intent(out) :: FLUX_INTEGRAL  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)

!  Profile weighting functions

      REAL(kind=8), intent(out) :: PROFILEWF ( MAX_ATMOSWFS,   MAXLAYERS,      MAX_USER_LEVELS, &
                                               MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

!  mean intensity weighting functions

      REAL(kind=8), intent(out) :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
 
!  flux weighting functions

      REAL(kind=8), intent(out) :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

!  Surface weighting functions at user angles

      REAL(kind=8), intent(out) :: SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS,  &
                                               MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

!  Flux and mean intensity surface weighting functions

      REAL(kind=8), intent(out) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
      REAL(kind=8), intent(out) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

!  Bookkeeping output
!  ------------------

!  Fourier numbers used

      INTEGER, intent(out)  ::  FOURIER_SAVED ( MAXBEAMS, MAXTHREADS ) 

!  Number of geometries (bookkeeping output)

      INTEGER, intent(out)  ::  N_GEOMETRIES

!  Exception handling 
!  ------------------

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER, intent(out)         :: STATUS_INPUTCHECK
      INTEGER, intent(out)         :: NCHECKMESSAGES
      CHARACTER*(*), intent(out)   :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(out)   :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER, intent(out)         :: STATUS_CALCULATION
      CHARACTER*(*), intent(out)   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

! ######################################################################

!                       Local Arguments
!                       +++++++++++++++

! ######################################################################

!  Local bookkeeping
!  ----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode of operation

      LOGICAL         :: DO_MSMODE_LIDORT

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      REAL(KIND=8)    :: QUAD_STREAMS (MAXSTREAMS)
      REAL(KIND=8)    :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(KIND=8)    :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(KIND=8)    :: USER_ANGLES  (MAX_USER_STREAMS)
      REAL(KIND=8)    :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(KIND=8)    :: USER_SECANTS  (MAX_USER_STREAMS)

!  output optical depth masks and indices

      LOGICAL         :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER         :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  off-grid optical depths (values, masks, indices)

      INTEGER         :: N_PARTLAYERS
      INTEGER         :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)
      REAL(KIND=8)    :: PARTLAYERS_VALUES       (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL         :: STERM_LAYERMASK_DN(MAXLAYERS)

!  indexing numbers

      INTEGER         :: N_VIEWING

!  Offsets for geometry indexing

      INTEGER         :: IBOFF(MAXBEAMS)
      INTEGER         :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(KIND=8)    :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)

!  Local solar zenith angles Cosines (regular case)

      REAL(KIND=8)    :: BEAM_COSINES(MAXBEAMS)

!  solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS
      INTEGER         :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests

      INTEGER         :: N_CONVTESTS

!  Adjusted geometries. New, 2007.
!  -------------------------------

!               Intent(Out) from the Adjust-geometry routine
!               Intent(In)  to   the Correction      routine

      REAL(KIND=8)    :: USER_ANGLES_ADJUST (MAX_USER_STREAMS)
      REAL(KIND=8)    :: BEAM_SZAS_ADJUST   (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      REAL(KIND=8)    :: USER_RELAZMS_ADJUST(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option

      INTEGER         :: LAYER_MAXMOMENTS (MAXLAYERS)

!  Initial transmittances * (secants)

      REAL(KIND=8)    :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(KIND=8)    :: TRUNC_FACTOR(MAXLAYERS)
      REAL(KIND=8)    :: FAC1(MAXLAYERS)

!  Derived Slant optical thickness inputs

      REAL(KIND=8)    :: TAUSLANT    ( 0:MAXLAYERS, MAXBEAMS )
      REAL(KIND=8)    :: TAUGRID     ( 0:MAXLAYERS )
      REAL(KIND=8)    :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(KIND=8)    :: DELTAU_SLANT_UNSCALED( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      REAL(KIND=8)    :: OMEGA_TOTAL    ( MAXLAYERS )
      REAL(KIND=8)    :: PHASMOMS_TOTAL ( 0:MAXMOMENTS, MAXLAYERS )

!  Linearized Input optical depths after delta-M scaling

      REAL(kind=8)    :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=8)    :: L_PHASMOMS_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS )

!  Derived Slant optical thickness inputs

      REAL(kind=8)    :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)
      
!  Linearized truncation factor

      REAL(kind=8)    :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  L'Hopital's rule logical variables

      LOGICAL         :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  coefficient functions for user-defined angles

      REAL(KIND=8)    :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(KIND=8)    :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Arrays for SS/DB corrections, 
!  =============================

!               Intent(in) to the Converge and LP_converge routines

!  Single scatter intensity

      REAL(KIND=8)    :: INTENSITY_SS &
          (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Direct bounce intensity

      REAL(KIND=8)    :: INTENSITY_DB (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  Single scatter Profile weighting functions at user angles

      REAL(kind=8)    :: PROFILEWF_SS  ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                         MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Direct Bounce Profile weighting functions at user angles

      REAL(kind=8)    :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Direct Bounce surface weighting functions at user angles

      REAL(kind=8)    :: SURFACEWF_DB  ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Fourier component output
!  ========================

!               Intent(Out) from the Fourier routine
!               Intent(in)  to   the Converge and LP_converge routines

!  Fourier comonents User-defined solutions

      REAL(KIND=8)    :: INTENSITY_F &
        (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Fourier-component Profile weighting functions at user angles

      REAL(kind=8)    :: PROFILEWF_F ( MAX_ATMOSWFS,    MAXLAYERS, MAX_USER_LEVELS, &
                                       MAX_USER_STREAMS,MAXBEAMS,  MAX_DIRECTIONS)

!  Fourier-component surface weighting functions at user angles

      REAL(kind=8)    :: SURFACEWF_F ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
                                       MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Help arrays from the SS/DB correction routines
!  ==============================================

!  Saved Legendre polynomials

      REAL(KIND=8)    :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      REAL(KIND=8)    :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(KIND=8)    :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      REAL(KIND=8)    :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(KIND=8)    :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      REAL(KIND=8)    :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)
 
!  Cumulative single scatter source terms

      REAL(KIND=8)    :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(KIND=8)    :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  Atmospheric attenuation before reflection

      REAL(KIND=8)    :: ATTN_DB_SAVE(MAX_GEOMETRIES)

!  Cumulative direct bounce source terms

      REAL(KIND=8)    :: DB_CUMSOURCE(MAX_GEOMETRIES,0:MAXLAYERS)

!  Outgoing sphericity stuff
!  Whole and part-layer LOS transmittance factors

      REAL(KIND=8)    :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(KIND=8)    :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(KIND=8)    :: BOA_ATTN(MAX_GEOMETRIES)

!  Outgoing sphericity stuff
!    Linearized output. Output required for later on.

      REAL(KIND=8)    :: L_UP_LOSTRANS   (MAXLAYERS,     MAX_ATMOSWFS,MAX_GEOMETRIES)
      REAL(KIND=8)    :: L_UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(kind=8)    :: L_BOA_ATTN(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      REAL(KIND=8)    :: DELTAU_VERT    ( MAXLAYERS )
      REAL(KIND=8)    :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(KIND=8)    :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )

!  Linearized input optical properties after delta-M scaling

      LOGICAL         :: DO_PHFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=8)    :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=8)    :: L_OMEGA_MOMS   ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      REAL(KIND=8)    :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAXBEAMS )

!  Chapman factors (from pseudo-spherical geometry)

      REAL(KIND=8)    :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER         :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(KIND=8)    :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(KIND=8)    :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(KIND=8)    :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=8)    :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(kind=8)    :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!     Solar beam attenuation

      REAL(KIND=8)    :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

!  Thread_dependent derive inputs
!  ------------------------------

!               Intent(InOut) To the Fourier routine

!  reflectance flags

      LOGICAL         :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL         :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL         :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      INTEGER         :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER         :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      LOGICAL         :: DO_BVTEL_INITIAL
      INTEGER         :: NLAYERS_TEL
      INTEGER         :: ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER         :: N_BVTELMATRIX_SIZE

!  set up for band matrix compression

      INTEGER         :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      INTEGER         :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!               Intent(In) To the Fourier routine

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(KIND=8)    :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage

      REAL(KIND=8)    :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(KIND=8)    :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      REAL(KIND=8)    :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage 

      REAL(KIND=8)    :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(KIND=8)    :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(KIND=8)    :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(KIND=8)    :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(KIND=8)    :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(KIND=8)    :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(KIND=8)    :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(KIND=8)    :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(KIND=8)    :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam

      REAL(kind=8)    :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8)    :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Multiplier arrays
!  -----------------

!               Intent(In) To the Fourier routine

!  forcing term multipliers (saved for whole atmosphere)

      REAL(KIND=8)    :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(KIND=8)    :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      REAL(KIND=8)    :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(KIND=8)    :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Linearized whole layer multipliers

      REAL(kind=8)    :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8)    :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(kind=8)    :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8)    :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(kind=8)    :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8)    :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Help variables
!  ==============

!  Local flags

      LOGICAL         :: LOCAL_DO_NO_AZIMUTH
      LOGICAL         :: SAVE_DO_NO_AZIMUTH
      LOGICAL         :: LOCAL_ITERATION
      LOGICAL         :: ADJUST_SURFACE
      LOGICAL         :: DO_PROCESS_FOURIER

!  Local fourier index

      INTEGER         :: FOURIER_COMPONENT
      INTEGER         :: N_FOURIER_COMPONENTS

!  local integers

      INTEGER         :: UA, UM, IB, TESTCONV, L
      INTEGER         :: LOCAL_N_USERAZM, STATUS_SUB
!      INTEGER         :: IUNIT, SUNIT, FUNIT, RUNIT

!  Flux multiplier

      REAL(KIND=8)    :: SS_FLUX_MULTIPLIER

!  Modified eradius

      REAL(KIND=8)    :: modified_eradius

!  Fourier cosine arguments

      REAL(KIND=8)    :: AZM_ARGUMENT, DFC
      REAL(KIND=8)    :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Local convergence control

      INTEGER         :: IBEAM_COUNT, IBEAM
      LOGICAL         :: BEAM_ITERATION ( MAXBEAMS )
      INTEGER         :: BEAM_TESTCONV  ( MAXBEAMS )

!  Local linearization variables

      INTEGER         :: LAYER_TO_VARY, N_PARAMETERS

!  Local error handling

      LOGICAL         :: FAIL
      CHARACTER*3     :: WTHREAD

!  ######################################################################
!  ######################################################################
!  ######################################################################

!  initialize output status
!  ------------------------

!  Main flags

      STATUS_CALCULATION = LIDORT_SUCCESS
      STATUS_INPUTCHECK  = LIDORT_SUCCESS

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

!  Single scatter correction: flux multiplier
!    Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

! dkh debug
! >>>
!!  Check input Basic. This could be put outside the thread loop.
!
!      CALL LIDORT_CHECK_INPUT_BASIC                                                 &
!      ( DO_SOLAR_SOURCES, DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY, DO_SSFULL,           & ! Input
!        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                              & ! Input
!        DO_CHAPMAN_FUNCTION, DO_REFRACTIVE_GEOMETRY, DO_NO_AZIMUTH,                 & ! Input/Output
!        DO_USER_STREAMS, DO_DELTAM_SCALING, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,    & ! Input/Output
!        DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, & ! Input/Output
!        NLAYERS, NFINELAYERS, NSTREAMS,  NMOMENTS_INPUT,                            & ! Input
!        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                      & ! Input
!        BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_LEVELS,                    & ! Input
!        EARTH_RADIUS, HEIGHT_GRID,                                                  & ! Input
!        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )                          ! output
!
!!  Exception handling
!
!      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
!        STATUS_INPUTCHECK = LIDORT_SERIOUS
!        RETURN
!      ELSE IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
!        STATUS_INPUTCHECK = LIDORT_WARNING
!      ENDIF
!
!!  Check input threaded values (IOPs, albedo)
!
!      CALL LIDORT_CHECK_INPUT_THREAD                               &
!         ( NLAYERS, THREAD, LAMBERTIAN_ALBEDO, DELTAU_VERT_INPUT,  & ! Input
!           OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                & ! Input
!           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS)       ! output
!
!!  Exception handling
!
!      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
!        STATUS_INPUTCHECK = LIDORT_SERIOUS
!        RETURN
!      ENDIF
! <<<

!  if there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
      ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  save some offsets for indexing geometries
!     These results are now stored in LIDORT_BOOKKEEP.VARS

      N_VIEWING    = N_USER_STREAMS * LOCAL_N_USERAZM
      N_GEOMETRIES = NBEAMS * N_VIEWING

      DO IBEAM = 1, NBEAMS
        IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        DO UM = 1, N_USER_STREAMS
          UMOFF(IBEAM,UM) = IBOFF(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
        END DO
      END DO

!  Geometry adjustment
!  -------------------

!  Adjust surface condition

      ADJUST_SURFACE = .FALSE.
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
         ADJUST_SURFACE = .TRUE.
        ENDIF
      ENDIF

!  Perform adjustment

      modified_eradius = earth_radius + GEOMETRY_SPECHEIGHT
      CALL multi_outgoing_adjustgeom                                   & 
        ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS,                & ! Input
          N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                  & ! Input
          height_grid(nlayers), modified_eradius, adjust_surface,      & ! Input
          user_angles_input,  beam_szas, user_relazms,                 & ! Input
          user_angles_adjust, beam_szas_adjust, user_relazms_adjust,   & ! Output
          fail, message, trace_1 )                                       ! Output

!  Exception handling

      if ( fail ) then
        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
        TRACE_3 = ' ** LIDORT_LPS_MASTER, THREAD # '//wthread
        STATUS_CALCULATION = LIDORT_SERIOUS
        return
      endif

!  Chapman function calculation
!  ----------------------------

!  Could be done outside the thread loop

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN

          CALL LIDORT_CHAPMAN                             &
          ( DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,    & ! Input
            NLAYERS, NBEAMS, FINEGRID, BEAM_SZAS,         & ! Input
            EARTH_RADIUS, RFINDEX_PARAMETER,              & ! Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, & ! Input
            CHAPMAN_FACTORS, SZA_LOCAL_INPUT,             & ! output
            FAIL, MESSAGE, TRACE_1 )                        ! output

          IF (FAIL) THEN
            TRACE_2 = 'Direct call in LIDORT_LPS_MASTER'
            TRACE_3 = ' ** LIDORT_LPS_MASTER, THREAD # '//wthread
            STATUS_CALCULATION = LIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF
      ENDIF

!  write input variables
!  ---------------------

!  open file, call standard input write, close file

!      IF ( DO_WRITE_INPUT ) THEN
!        IUNIT = LIDORT_INUNIT
!        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')
!        CALL LIDORT_WRITEINPUT ( THREAD, IUNIT )
!        CLOSE(IUNIT)
!      ENDIF

!  Get derived inputs
!  ==================

!  Miscellaneous and layer input.
!    This could be put outside the thread loop
!   ( 6/21/10) Important point: Must use ADJUSTED VZangles as input here.

      CALL LIDORT_DERIVE_INPUT_BASIC                                 &
      ( DO_FULLRAD_MODE, DO_USER_STREAMS,                            & ! Input
        DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                         & ! Input
        DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,              & ! Input
        DO_UPWELLING, DO_DNWELLING, DO_REFRACTIVE_GEOMETRY,          & ! Input
        DO_DOUBLE_CONVTEST, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,  & ! Input or Output
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,                          & ! Input or Output
        NLAYERS, NSTREAMS, NBEAMS, NMOMENTS_INPUT,                   & ! Input
        N_USER_STREAMS, N_USER_RELAZMS, USER_ANGLES_ADJUST,          & ! Input
        N_USER_LEVELS,  USER_LEVELS, BEAM_SZAS, SZA_LOCAL_INPUT,     & ! Input
        DO_MSMODE_LIDORT, NMOMENTS, BEAM_COSINES, SUNLAYER_COSINES,  & ! Output
        N_CONVTESTS, N_DIRECTIONS, WHICH_DIRECTIONS,                 & ! Output
        NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_PARTLAYERS,      & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                    & ! Output
        USER_ANGLES, USER_STREAMS, USER_SECANTS,                     & ! Output
        PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX,PARTLAYERS_LAYERIDX,  & ! Output
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_VALUES,   & ! Output
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN )                       ! Output

!  #################
!  Set up operations
!  #################

!   MISCSETUPS (5 subroutines)  :
!   -----------------------------

!       Performance Setup,
!       Delta-M scaling,
!       average-secant formulation,
!       transmittance setup
!       Beam solution multipliers

!  Performance set-up

      CALL LIDORT_PERFORMANCE_SETUP                                        &
        ( DO_SOLAR_SOURCES, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,           & ! input
          DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,                          & ! input
          DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                             & ! input
          THREAD, NLAYERS, NMOMENTS, NMOMENTS_INPUT, PHASMOMS_TOTAL_INPUT, & ! input
          LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG,         & ! Output
          STATUS_SUB, MESSAGE, TRACE_1)                                      ! Output

!  Exception handling on this module

      IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        TRACE_2 = 'Direct call in LIDORT_LPS_MASTER'
        TRACE_3 = ' ** LIDORT_LPS_MASTER, THREAD # '//wthread
        STATUS_CALCULATION = LIDORT_WARNING
      ENDIF 

!  Delta-m scaling of input quantities

      CALL LIDORT_DELTAMSCALE                                        &
       ( DO_SOLAR_SOURCES, DO_DELTAM_SCALING,                        & ! input
         NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS, N_USER_LEVELS,     & ! input
         PARTLAYERS_OUTFLAG, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & ! input
         THREAD, CHAPMAN_FACTORS, DELTAU_VERT_INPUT,                 & ! input
         OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                    & ! input
         DELTAU_VERT, PARTAU_VERT, TAUGRID,                          & ! Output
         OMEGA_TOTAL, OMEGA_MOMS, PHASMOMS_TOTAL,                    & ! Output
         FAC1, TRUNC_FACTOR,                                         & ! Output
         DELTAU_SLANT, DELTAU_SLANT_UNSCALED )                         ! Output

!  Prepare quasi-spherical attenuation

      CALL LIDORT_QSPREP                                                &
       ( DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,   & ! input
         NLAYERS, NBEAMS, BEAM_COSINES, SUNLAYER_COSINES,               & ! input
         TAUGRID, DELTAU_VERT, DELTAU_SLANT,                            & ! input
         INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,   & ! Output
         TAUSLANT, SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM )            ! Output

!  Transmittances and Transmittance factors

      CALL LIDORT_PREPTRANS                                            &
        ( DO_SOLAR_SOURCES, DO_SOLUTION_SAVING, DO_USER_STREAMS,       & ! input
          NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS, N_PARTLAYERS,     & ! input
          QUAD_STREAMS, DELTAU_VERT, PARTLAYERS_LAYERIDX, PARTAU_VERT, & ! input
          USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,        & ! input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,             & ! input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output
          T_DELT_MUBAR,   T_UTUP_MUBAR,   T_UTDN_MUBAR,                & ! Output
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                & ! Output
          ITRANS_USERM  )                                                ! Output

!  Linearization of Delta-M scaled inputs

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_LA_DELTAMSCALE                                &
       ( DO_DELTAM_SCALING, NLAYERS, NMOMENTS, NBEAMS, THREAD,    & ! Input
         CHAPMAN_FACTORS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     & ! Input
         DELTAU_VERT_INPUT,    L_DELTAU_VERT_INPUT,               & ! Input
         OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,               & ! Input
         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,            & ! Input
         OMEGA_MOMS, FAC1, TRUNC_FACTOR,                          & ! Input
         L_DELTAU_VERT, L_DELTAU_SLANT,                           & ! Output
         L_OMEGA_TOTAL, L_OMEGA_MOMS, L_PHASMOMS_TOTAL,           & ! Output
         L_TRUNC_FACTOR, DO_PHFUNC_VARIATION )                      ! Output
      ENDIF

!  Linearization of transmittances 

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_LA_PREPTRANS                                 &
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
      ENDIF

!  Linearization of pseudo-spherical setup

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_LP_PREPTRANS                                  &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,       & ! Input
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                 & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,  & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,        & ! Input
          T_DELT_MUBAR, T_UTDN_MUBAR,                             & ! Input
          LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR, LP_INITIAL_TRANS,    & ! Output
          LP_AVERAGE_SECANT, HELP_AQ, HELP_BQ )                     ! Output
      ENDIF

!   EMULT_MASTER  : Beam source function multipliers. Not required for the
!                  Full SS calculation in outgoing mode

      IF ( DO_SOLAR_SOURCES  ) THEN
        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

          CALL EMULT_MASTER                                       & 
       ( DO_UPWELLING, DO_DNWELLING, N_USER_STREAMS, NBEAMS,      & ! input
         NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,              & ! input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! input
         DELTAU_VERT, PARTAU_VERT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,            & ! input
         ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF,          & ! input
         SIGMA_M, SIGMA_P, EMULT_HOPRULE,                         & ! Output
         EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )             ! Output


          IF ( DO_PROFILE_LINEARIZATION ) THEN
            CALL LP_EMULT_MASTER                                        &
         ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,               & ! Input
           NLAYERS, N_PARTLAYERS, N_USER_STREAMS, NBEAMS,               & ! Input
           STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, PARTLAYERS_LAYERIDX, & ! Input
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                          & ! Input
           LAYER_PIS_CUTOFF, ITRANS_USERM,                              & ! Input
           T_DELT_USERM, T_UTUP_USERM, T_DELT_MUBAR, USER_SECANTS,      & ! Input
           DELTAU_VERT, PARTAU_VERT,  L_DELTAU_VERT, EMULT_HOPRULE,     & ! Input
           LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                         & ! Input
           LP_T_DELT_MUBAR, L_T_DELT_USERM,                             & ! Input
           LP_T_UTDN_MUBAR, L_T_UTDN_USERM,                             & ! Input
           SIGMA_P, EMULT_UP, UT_EMULT_UP,                              & ! Input
           SIGMA_M, EMULT_DN, UT_EMULT_DN,                              & ! Input
           LP_EMULT_UP, LP_UT_EMULT_UP,                                 & ! Output
           LP_EMULT_DN, LP_UT_EMULT_DN )                                  ! Output
          ENDIF

        ENDIF
      ENDIF

!  #####################
!  Correction operations 
!  #####################

!  Single scatter corrections (pre-calculation)
!  ============================================

!     Must be done after MISCSETUPS and EMULT_MASTER, as we need
!      multipliers and transmittance factors for SUN and LOS paths.
!      Code added 6 May 2005. Replaces call in Master routine.
!      Version 3.1. Added call to the new outgoing sphericity correction

      IF ( DO_USER_STREAMS ) THEN

!  regular nadir-scattering SS correction
!  --------------------------------------

        IF ( DO_SSCORR_NADIR ) THEN

!  Intensity

          CALL LIDORT_SSCORR_NADIR                                       &
        ( DO_UPWELLING, DO_DNWELLING, DO_SSCORR_TRUNCATION,              & ! Input
          DO_DELTAM_SCALING, DO_REFRACTIVE_GEOMETRY, THREAD,             & ! Input
          SS_FLUX_MULTIPLIER, NLAYERS, NMOMENTS_INPUT, NBEAMS,           & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                 & ! Input
          BEAM_SZAS, SUNLAYER_COSINES, USER_STREAMS, USER_RELAZMS,       & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,   & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                        & ! Input
          TRUNC_FACTOR, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,         & ! Input
          EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                  & ! Input
          T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                      & ! Input
          SS_PLEG_UP, SS_PLEG_DN, TMS, SSFDEL,                           & ! Output  
          EXACTSCAT_UP, EXACTSCAT_DN,                                    & ! Output  
          SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, INTENSITY_SS )                 ! Output  

!  profile linearization

          IF ( DO_PROFILE_LINEARIZATION ) THEN

            DO LAYER_TO_VARY = 1, NLAYERS

              IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN

                N_PARAMETERS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

                CALL LIDORT_LP_SSCORR_NADIR                                &
        ( DO_UPWELLING, DO_DNWELLING, DO_SSCORR_TRUNCATION,                & ! Input
          DO_DELTAM_SCALING, NLAYERS, NMOMENTS_INPUT, NMOMENTS,            & ! Input
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,           & ! Input
          LAYER_TO_VARY, LAYER_VARY_NUMBER, THREAD, SS_FLUX_MULTIPLIER,    & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,     & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,    & ! Input
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, DO_PHFUNC_VARIATION,     & ! Input
          TRUNC_FACTOR, PHASMOMS_TOTAL_INPUT,                              & ! Input
          L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                     & ! Input
          SS_PLEG_UP, SS_PLEG_DN, TMS, SSFDEL,                             & ! Input
          EXACTSCAT_UP, EXACTSCAT_DN, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN,    & ! Input
          EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                    & ! Input
          LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,        & ! Input
          T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                        & ! Input
          L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                  & ! Input
          PROFILEWF_SS )                                                     ! Output

              ENDIF

            ENDDO

          ENDIF

        ENDIF

!  New outgoing sphericity correction
!  ----------------------------------

        IF ( DO_SSCORR_OUTGOING ) THEN

          CALL LIDORT_LP_SSCORR_OUTGOING                                  & 
        ( DO_UPWELLING, DO_DNWELLING, DO_PROFILE_LINEARIZATION,           & ! Input
          DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING,                        & ! Input
          NLAYERS, NFINELAYERS, NMOMENTS_INPUT, NMOMENTS,                 & ! Input
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,          & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER, THREAD, SS_FLUX_MULTIPLIER, & ! Input 
          BEAM_SZAS_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,      & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input 
          N_PARTLAYERS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,           & ! Input
          HEIGHT_GRID, EARTH_RADIUS, DO_PHFUNC_VARIATION,                 & ! Input
          TRUNC_FACTOR, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,          & ! Input
          OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,                      & ! Input
          PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                   & ! Input
            UP_LOSTRANS,   UP_LOSTRANS_UT,   BOA_ATTN,                    & ! Output
          L_UP_LOSTRANS, L_UP_LOSTRANS_UT, L_BOA_ATTN,                    & ! Output
          INTENSITY_SS, PROFILEWF_SS,                                     & ! Output
          STATUS_SUB, MESSAGE, TRACE_1  )                                  ! Output  

          IF ( STATUS_SUB .ne. LIDORT_SUCCESS ) THEN
            TRACE_2 = 'Error from LP_LIDORT_SSCORR_OUTGOING'
            TRACE_3 = ' ** LIDORT_LPS_MASTER, THREAD # '//wthread
            STATUS_CALCULATION = LIDORT_SERIOUS
            RETURN
          ENDIF 

        ENDIF

!  Exact (Direct bounce) surface term if single scatter is being done
!     renamed after introduction of BRDF option, 23 March 2010

        IF ( DO_SSCORR_NADIR .or. DO_SSCORR_OUTGOING ) THEN

! Intensity term

          CALL LIDORT_DBCORRECTION                                       &
        ( DO_SSCORR_OUTGOING, DO_BRDF_SURFACE,                           & ! Input
          DO_REFRACTIVE_GEOMETRY, DO_UPWELLING,                          & ! Input
          DO_REFLECTED_DIRECTBEAM, SS_FLUX_MULTIPLIER, NLAYERS, NBEAMS,  & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                 & ! Input
          BEAM_SZAS, BEAM_SZAS_ADJUST, SZA_LOCAL_INPUT,                  & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP,                       & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          LAMBERTIAN_ALBEDO(THREAD), EXACTDB_BRDFUNC,                    & ! Input
          SOLAR_BEAM_OPDEP, BOA_ATTN,                                    & ! Input
          T_DELT_USERM, UP_LOSTRANS, T_UTUP_USERM, UP_LOSTRANS_UT,       & ! Input
          INTENSITY_DB, ATTN_DB_SAVE, DB_CUMSOURCE )                       ! Output

!  Profile LInearization

          IF ( DO_PROFILE_LINEARIZATION ) THEN

            CALL LIDORT_LAP_DBCORRECTION                                    &
        ( DO_SSCORR_OUTGOING, DO_BRDF_SURFACE,                              & ! Input
          DO_REFRACTIVE_GEOMETRY, DO_UPWELLING,                             & ! Input
          DO_REFLECTED_DIRECTBEAM, SS_FLUX_MULTIPLIER, NLAYERS, NBEAMS,     & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                    & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                               & ! Input
          BEAM_SZAS, BEAM_SZAS_ADJUST, SZA_LOCAL_INPUT,                     & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP,                          & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input
          SOLAR_BEAM_OPDEP, DELTAU_SLANT, L_DELTAU_VERT, L_BOA_ATTN,        & ! Input
          T_DELT_USERM, T_UTUP_USERM, L_T_DELT_USERM, L_T_UTUP_USERM,       & ! Input
          UP_LOSTRANS, UP_LOSTRANS_UT, L_UP_LOSTRANS, L_UP_LOSTRANS_UT,     & ! Input
          LAMBERTIAN_ALBEDO(THREAD), EXACTDB_BRDFUNC, DB_CUMSOURCE,         & ! Input 
          PROFILEWF_DB )                                                      ! Output

          ENDIF

!  Surface  linearization

          IF ( DO_SURFACE_LINEARIZATION ) THEN

            CALL LIDORT_LS_DBCORRECTION                                 &
        ( DO_SSCORR_OUTGOING, DO_UPWELLING, DO_BRDF_SURFACE,            & ! Input
          DO_REFLECTED_DIRECTBEAM, NLAYERS, NBEAMS, N_SURFACE_WFS,      & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP, SS_FLUX_MULTIPLIER,  & ! Input
          LS_EXACTDB_BRDFUNC, ATTN_DB_SAVE,                             & ! Input  
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          T_DELT_USERM, UP_LOSTRANS, T_UTUP_USERM, UP_LOSTRANS_UT,      & ! Input
          SURFACEWF_DB )                                                  ! Output

          ENDIF

        ENDIF

!  End User streams clause

      ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  Local flags

      IF ( .NOT. DO_SOLAR_SOURCES  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_ISOTROPIC_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_MVOUT_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_SSFULL  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
        DO_PROCESS_FOURIER  = .false.
      ELSE
        DO_PROCESS_FOURIER  = .true.
      ENDIF

!  set Fourier number (2 for Rayleigh only, otherwise 2*Nstreams))
!    Local do no azimuth skips convergence

      IF ( LOCAL_DO_NO_AZIMUTH ) THEN
        N_FOURIER_COMPONENTS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIER_COMPONENTS = 2
        ELSE
          N_FOURIER_COMPONENTS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  Fourier loop
!  ============

!  initialize

      LOCAL_ITERATION   = .TRUE.
      FOURIER_COMPONENT = -1
      TESTCONV          = 0

!  set up solar beam flags

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

!  start loop
!  ----------

      DO WHILE ( LOCAL_ITERATION .AND. FOURIER_COMPONENT.LT.N_FOURIER_COMPONENTS )

!  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

!  Local start of user-defined streams
!  Now set = 1. Fudging in earlier versions caused Havoc.
!  Here is the older version fudge
!        IF ( FOURIER_COMPONENT. EQ. 0 ) THEN
!          LOCAL_UM_START = 1
!        ELSE
!          LOCAL_UM_START = N_OUT_STREAMS - N_CONV_STREAMS + 1
!        ENDIF
!        LOCAL_UM_START = 1

!  azimuth cosine factor, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          DO UA = 1, LOCAL_N_USERAZM
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,IB,UA) * DFC
                AZMFAC(UM,IB,UA)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Main call to Lidort Fourier module
!  ----------------------------------

!  Only if flagged

        IF ( DO_PROCESS_FOURIER ) THEN

!        write(*,*)' ..calculating fourier component',FOURIER_COMPONENT

          CALL LIDORT_LPS_FOURIER_MASTER                                 &
       ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,                   & ! Input
         DO_BRDF_SURFACE, DO_PLANE_PARALLEL,                             & ! input
         DO_USER_STREAMS, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,        & ! input/Output
         DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_REFRACTIVE_GEOMETRY,     & ! input/Output
         NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, NBEAMS,       & ! input
         DO_MSMODE_LIDORT, DO_MULTIBEAM, THREAD, FOURIER_COMPONENT,      & ! input
         DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,             & ! input 
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                             & ! input 
         FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES,                    & ! input
         NMOMENTS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! input
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                       & ! input
         USER_STREAMS, USER_SECANTS, LAMBERTIAN_ALBEDO,                  & ! input
         N_SURFACE_WFS, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,    & ! input
         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,       & ! input
         PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX,PARTLAYERS_LAYERIDX,     & ! input
         UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                         & ! input     
         STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                         & ! input    
         DELTAU_VERT, PARTAU_VERT, OMEGA_MOMS, DELTAU_SLANT,             & ! input 
         SOLAR_BEAM_OPDEP, L_DELTAU_VERT, L_OMEGA_MOMS,                  & ! input 
         INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,    & ! input
         DO_REFLECTED_DIRECTBEAM,DO_LAYER_SCATTERING,BVP_REGULAR_FLAG,   & ! input/Output
         LCONMASK, MCONMASK, BMAT_ROWMASK, BTELMAT_ROWMASK,              & ! input/Output    
         DO_BVTEL_INITIAL,NLAYERS_TEL,ACTIVE_LAYERS,N_BVTELMATRIX_SIZE,  & ! input/Output 
         T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN, T_DELT_MUBAR,   & ! input 
         T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,         & ! input
         EMULT_UP,    EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                & ! input
         LP_INITIAL_TRANS, LP_AVERAGE_SECANT, HELP_AQ, HELP_BQ,          & ! input
         L_T_DELT_DISORDS,  L_T_DISORDS_UTUP,  L_T_DISORDS_UTDN,         & ! input
         LP_T_DELT_MUBAR,   LP_T_UTDN_MUBAR,                             & ! input
         L_T_DELT_USERM,    L_T_UTUP_USERM,    L_T_UTDN_USERM,           & ! input
         LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,       & ! input
         INTENSITY_F, MEAN_INTENSITY, FLUX_INTEGRAL,                     & ! Output 
         PROFILEWF_F, MINT_PROFILEWF, FLUX_PROFILEWF,                    & ! Output
         SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                    & ! Output
         STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                           ! Output

!  error handling

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            TRACE_3 = ' ** LIDORT_LPS_MASTER, THREAD # '//wthread
            STATUS_CALCULATION = LIDORT_SERIOUS
            RETURN
          ENDIF

!  End fourier processing

        ENDIF

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!   -- only done for beams which are still not converged
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

        IBEAM_COUNT = 0
        DO IBEAM = 1, NBEAMS
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

!  Convergence and radiance summation

            CALL LIDORT_CONVERGE                                      &
           ( DO_UPWELLING, DO_NO_AZIMUTH,                             & ! Input
             DO_RAYLEIGH_ONLY, DO_ALL_FOURIER,                        & ! Input
             DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,          & ! Input
             DO_DOUBLE_CONVTEST, N_CONVTESTS, LIDORT_ACCURACY,        & ! Input
             N_USER_STREAMS, N_USER_LEVELS, N_USER_RELAZMS,           & ! Input
             NSTREAMS, IBEAM, THREAD, FOURIER_COMPONENT,              & ! Input
             UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,  & ! Input
             AZMFAC, INTENSITY_F, INTENSITY_SS, INTENSITY_DB,         & ! Input
             INTENSITY, FOURIER_SAVED,                                & ! output
             BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )              ! output

!  Fourier summation of Profile Jacobians

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              CALL LIDORT_LP_CONVERGE                 &
           ( DO_UPWELLING, DO_SSFULL,                 & ! Input       
             DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,     & ! Input
             DO_NO_AZIMUTH, AZMFAC, LOCAL_N_USERAZM,  & ! Input
             N_USER_STREAMS, N_USER_LEVELS, NLAYERS,  & ! Input
             IBEAM, THREAD, FOURIER_COMPONENT,        & ! Input
             UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,   & ! Input
             LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Input
             PROFILEWF_F, PROFILEWF_SS, PROFILEWF_DB, & ! Input
             PROFILEWF )                                ! Output
            ENDIF

!  Fourier summation of Surface Jacobians

           IF ( DO_SURFACE_LINEARIZATION ) THEN
              CALL LIDORT_LS_CONVERGE                                      &
              ( DO_UPWELLING, DO_BRDF_SURFACE,                             & ! Input
                DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,            & ! Input
                N_USER_STREAMS, N_USER_LEVELS, N_SURFACE_WFS,              & ! Input
                LOCAL_N_USERAZM, AZMFAC, IBEAM, THREAD, FOURIER_COMPONENT, & ! Input
                UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,                     & ! Input
                SURFACEWF_F,  SURFACEWF_DB,                                & ! Input
                SURFACEWF )                                                  ! Output
           ENDIF

!  Check number of beams already converged

            IF ( BEAM_ITERATION(IBEAM) ) THEN
              IBEAM_COUNT = IBEAM_COUNT + 1
            ELSE
              DO L = FOURIER_COMPONENT+1,MAXFOURIER
                DO_MULTIBEAM (IBEAM,L) = .FALSE.
              ENDDO
            ENDIF

!  end beam count loop

          ENDIF
        END DO

!  If all beams have converged, stop iteration

        IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

!  Fourier output
!  --------------

!  Open file if Fourier = 0
!  Write Standard Fourier output
!  Close file if iteration has finished

!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field

!        IF ( DO_WRITE_FOURIER ) THEN
!          FUNIT = LIDORT_FUNIT
!          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
!            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
!          ENDIF
!          CALL LIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
!          IF ( DO_PROFILE_LINEARIZATION ) THEN
!            CALL LIDORT_LP_WRITEFOURIER ( DO_INCLUDE_SURFACE, FUNIT, FOURIER_COMPONENT )
!          ENDIF
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  end iteration loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Major result output
!  -------------------

!      IF ( DO_WRITE_RESULTS ) THEN
!        RUNIT = LIDORT_RESUNIT
!        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')
!        CALL LIDORT_WRITERESULTS ( RUNIT )
!        IF ( DO_PROFILE_LINEARIZATION ) THEN
!          CALL LIDORT_LP_WRITERESULTS ( RUNIT, LOCAL_N_USERAZM )
!        ENDIF
!        CLOSE(RUNIT)
!      ENDIF

!  Geophysical input (scenario) write
!  ----------------------------------

!      IF ( DO_WRITE_SCENARIO ) THEN
!        SUNIT = LIDORT_SCENUNIT
!        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')
!        CALL LIDORT_WRITESCEN ( SUNIT )
!        IF ( DO_PROFILE_LINEARIZATION ) THEN
!          CALL LIDORT_LP_WRITESCEN ( SUNIT )
!        ENDIF
!        CLOSE(SUNIT)
!      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LPS_MASTER

!

SUBROUTINE LIDORT_LPS_FOURIER_MASTER                                     &
       ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,                   & ! Input
         DO_BRDF_SURFACE, DO_PLANE_PARALLEL,                             & ! input
         DO_USER_STREAMS, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,        & ! input/Output
         DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_REFRACTIVE_GEOMETRY,     & ! input/Output
         NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, NBEAMS,       & ! input
         DO_MSMODE_LIDORT, DO_MULTIBEAM, THREAD, FOURIER_COMPONENT,      & ! input
         DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,             & ! input 
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                             & ! input 
         FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES,                    & ! input
         NMOMENTS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! input
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                       & ! input
         USER_STREAMS, USER_SECANTS, LAMBERTIAN_ALBEDO,                  & ! input
         N_SURFACE_WFS, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,    & ! input
         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,       & ! input
         PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX,PARTLAYERS_LAYERIDX,     & ! input
         UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                         & ! input     
         STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                         & ! input    
         DELTAU_VERT, PARTAU_VERT, OMEGA_MOMS, DELTAU_SLANT,             & ! input 
         SOLAR_BEAM_OPDEP, L_DELTAU_VERT, L_OMEGA_MOMS,                  & ! input 
         INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,    & ! input
         DO_REFLECTED_DIRECTBEAM,DO_LAYER_SCATTERING,BVP_REGULAR_FLAG,   & ! input/Output
         LCONMASK, MCONMASK, BMAT_ROWMASK, BTELMAT_ROWMASK,              & ! input/Output    
         DO_BVTEL_INITIAL,NLAYERS_TEL,ACTIVE_LAYERS,N_BVTELMATRIX_SIZE,  & ! input/Output 
         T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN, T_DELT_MUBAR,   & ! input 
         T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,         & ! input
         EMULT_UP,    EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                & ! input
         LP_INITIAL_TRANS, LP_AVERAGE_SECANT, HELP_AQ, HELP_BQ,          & ! input
         L_T_DELT_DISORDS,  L_T_DISORDS_UTUP,  L_T_DISORDS_UTDN,         & ! input
         LP_T_DELT_MUBAR,   LP_T_UTDN_MUBAR,                             & ! input
         L_T_DELT_USERM,    L_T_UTUP_USERM,    L_T_UTDN_USERM,           & ! input
         LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,       & ! input
         INTENSITY_F, MEAN_INTENSITY, FLUX_INTEGRAL,                     & ! Output 
         PROFILEWF_F, MINT_PROFILEWF, FLUX_PROFILEWF,                    & ! Output
         SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                    & ! Output
         STATUS, MESSAGE, TRACE_1, TRACE_2 )                               ! Output

!  Complete Fourier component calculation for the Profile-Jacobian Code.

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE 'LIDORT.PARS_F90'

!  input
!  -----

!  Basic top-level control

      LOGICAL, intent(in)  :: DO_SOLAR_SOURCES

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(in)  :: DO_PLANE_PARALLEL

!  stream angle flag

      LOGICAL, intent(in)  :: DO_USER_STREAMS

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Surface control (New, 23 March 2010)

      LOGICAL, intent(in)  :: DO_BRDF_SURFACE

!  directional control

      LOGICAL, intent(in)  :: DO_UPWELLING
      LOGICAL, intent(in)  :: DO_DNWELLING

!  Number of discrete ordinate streams

      INTEGER, intent(in)  :: NSTREAMS

!  number of computational layers

      INTEGER, intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER, intent(in)  :: NBEAMS

!  Local input solar zenith angles Cosines

      REAL(KIND=8), intent(in)  :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)
      REAL(KIND=8), intent(in)  :: BEAM_COSINES(MAXLAYERS,MAXBEAMS)

!  User-defined zenith angle input 

      INTEGER, intent(in)  :: N_USER_STREAMS

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      INTEGER, intent(in)  :: N_USER_LEVELS

!  mean value control

      LOGICAL, intent(in)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, intent(in)  :: DO_MVOUT_ONLY

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, intent(in)  :: LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER, intent(in)  :: LAYER_VARY_NUMBER(MAXLAYERS)

!  Total number of Surface Jacobians

      INTEGER, intent(in)  ::  N_SURFACE_WFS

!  Control linearization
 
      LOGICAL, intent(in)  ::  DO_PROFILE_LINEARIZATION
      LOGICAL, intent(in)  ::  DO_SURFACE_LINEARIZATION

!  Input modified arguments
!  ------------------------

!  Performance enhancement

      LOGICAL, intent(inout)  :: DO_SOLUTION_SAVING
      LOGICAL, intent(inout)  :: DO_BVP_TELESCOPING

!  Bookkeeping arguments
!  ---------------------

!  Mode of operation

      LOGICAL, intent(in)  :: DO_MSMODE_LIDORT

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER, intent(in)  :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER, intent(in)  :: NSTREAMS_2
      INTEGER, intent(in)  :: NTOTAL
      INTEGER, intent(in)  :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      REAL(KIND=8), intent(in)   :: QUAD_STREAMS (MAXSTREAMS)
      REAL(KIND=8), intent(in)   :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(KIND=8), intent(in)   :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(KIND=8), intent(in)   :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(KIND=8), intent(in)   :: USER_SECANTS  (MAX_USER_STREAMS)

!  Offgrid output optical depth masks and indices

      LOGICAL, intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL, intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  solar beam flags (always internal)

      LOGICAL, intent(in)  :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Thread number

      INTEGER, intent(in)  :: THREAD

!  Input Fourier component number

      INTEGER, intent(in)  :: FOURIER_COMPONENT

!  Absoute flux factor

      REAL(KIND=8), intent(in)  :: FLUX_FACTOR

!  Lambertian Surface control

      REAL(KIND=8), intent(in)  :: LAMBERTIAN_ALBEDO (MAXTHREADS)

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(KIND=8), intent(in)  :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in)  :: BRDF_F        ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(KIND=8), intent(in)  :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in)  :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Fourier components of BRDF, in the following order (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(KIND=8), intent(in) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

      REAL(KIND=8), intent(in) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Outputs
!  =======

!  Fourier component solutions

      REAL(KIND=8), intent(Out) :: INTENSITY_F (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Results for mean-value output

      REAL(KIND=8), intent(Out) :: MEAN_INTENSITY (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
      REAL(KIND=8), intent(Out) :: FLUX_INTEGRAL  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)

!  Fourier-component Profile weighting functions at user angles

      REAL(kind=8), intent(out) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                 MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  mean intensity weighting functions

      REAL(kind=8), intent(out) :: MINT_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
 
!  flux weighting functions

      REAL(kind=8), intent(out) :: FLUX_PROFILEWF ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                    MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

!  Surface weighting functions at user angles

      REAL(kind=8), intent(out) :: SURFACEWF_F ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
                                                 MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Flux and mean intensity Surface weighting functions

      REAL(kind=8), intent(out) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                    MAXBEAMS,       MAX_DIRECTIONS, MAXTHREADS)
      REAL(kind=8), intent(out) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                    MAXBEAMS,       MAX_DIRECTIONS, MAXTHREADS)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER, intent(out)         :: STATUS
      CHARACTER*(*), intent(out)   :: MESSAGE, TRACE_1, TRACE_2

!  Arrays required at the Top level
!  ================================

!                    Intent(In) to the Fourier routine

!  Input optical depths after delta-M scaling

      REAL(KIND=8), intent(in) :: DELTAU_VERT    ( MAXLAYERS )
      REAL(KIND=8), intent(in) :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(KIND=8), intent(in) :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )
      REAL(kind=8), intent(in) :: DELTAU_SLANT   ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Linearized input optical properties after delta-M scaling

      REAL(kind=8), intent(in) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=8), intent(in) :: L_OMEGA_MOMS  ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  Last layer to include Particular integral solution

      INTEGER, intent(in)      :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(KIND=8), intent(in) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(KIND=8), intent(in) :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=8), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(kind=8), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized help arrays for Green's function

      REAL(kind=8), intent(in)  :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Solar beam attenuation

      REAL(KIND=8), intent(in)   :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

!  Thread_dependent derive inputs
!  ------------------------------

!                    Intent(InOut) to the Fourier routine

!  reflectance flags

      LOGICAL, intent(InOut)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL, intent(InOut)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL, intent(InOut)  :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      INTEGER, intent(inout)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER, intent(inout)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      LOGICAL, intent(inout)  :: DO_BVTEL_INITIAL
      INTEGER, intent(inout)  :: NLAYERS_TEL
      INTEGER, intent(inout)  :: ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, intent(inout)  :: N_BVTELMATRIX_SIZE

!  set up for band matrix compression

      INTEGER, intent(inout)  :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      INTEGER, intent(inout)  :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!                    Intent(In) to the Fourier routine

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(KIND=8), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(KIND=8), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(KIND=8), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(KIND=8), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(KIND=8), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(KIND=8), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!                 Intent(In) to Fourier Routine
 
!  forcing term multipliers (saved for whole atmosphere)

      REAL(KIND=8), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(KIND=8), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      REAL(KIND=8), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(KIND=8), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Linearized Transmittance Setups
!  -------------------------------

!                    Intent(In) to the Fourier routine

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(kind=8), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=8), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized transmittances, solar beam

      REAL(kind=8), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Beam multipliers
!  ---------------------------

!                    Intent(In) to the Fourier routine

!  Linearized whole layer multipliers

      REAL(kind=8), intent(in)  :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(kind=8), intent(in)  :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=8), intent(in)  :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local Arrays for argument passing
!  =================================

!  Atmospheric attenuation

      REAL(KIND=8)    :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions

      REAL(KIND=8)    :: DIRECT_BEAM      ( MAXSTREAMS,       MAXBEAMS )
      REAL(KIND=8)    :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(KIND=8)    :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8)    :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(KIND=8)    :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: GAMMA_P(MAXSTREAMS,MAXLAYERS)
 
!  Integrated homogeneous solution multipliers, whole layer

      REAL(KIND=8)    :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8)    :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(KIND=8)    :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Legendre Setups
!  ---------------

!  At quadrature angles

      REAL(KIND=8)    :: LEG_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8)    :: LEG_M(MAXSTREAMS,0:MAXMOMENTS)

!  At beam angles. LEG0_M holds stored quantities.

      REAL(KIND=8)    :: LEG0_P(0:MAXMOMENTS)
      REAL(KIND=8)    :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Legendre polynomial products

      REAL(KIND=8)    :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8)    :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8)    :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(KIND=8)    :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Polynomial-weight products

      REAL(KIND=8)    :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8)    :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(KIND=8)    :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(KIND=8)    :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(KIND=8)    :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8)    :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8)    :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8)    :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8)    :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!  (Positive) Eigenvalues

      REAL(KIND=8)    :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(KIND=8)    :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(KIND=8)    :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(KIND=8)    :: U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(KIND=8)    :: U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(KIND=8)    :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Help arrays for reflected solutions

      REAL(KIND=8)    :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(KIND=8)    :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  Boundary Value Problem
!  ----------------------

!  Single Matrix, Band-matrices

      REAL(KIND=8)    :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(KIND=8)    :: BANDMAT2   (MAXBANDTOTAL,MAXTOTAL)
      REAL(KIND=8)    :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER         :: IPIVOT     (MAXTOTAL)
      INTEGER         :: SIPIVOT    (MAXSTREAMS_2)
      INTEGER         :: IPIVOTTEL  (MAXTOTAL)

!  particular integrals
!  --------------------

!  General beam solutions at the boundaries

      REAL(KIND=8)    :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(KIND=8)    :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Help array for reflected solutions

      REAL(KIND=8)    :: H_WLOWER(MAXSTREAMS)

!  Green's function particular integral arrays

      REAL(KIND=8)    :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)
      REAL(KIND=8)    :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)
      REAL(KIND=8)    :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: BGP(MAXSTREAMS,MAXLAYERS)

!  Layer ! and D functions
!  Green function Multipliers for solution

      REAL(KIND=8)    :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(KIND=8)    :: W_HELP(0:MAXMOMENTS)

!  Particular beam solutions at user-defined stream angles

      REAL(KIND=8)    :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8)    :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(KIND=8)    :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8)    :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8)    :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(KIND=8)    :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(KIND=8)    :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      LOGICAL         :: FLAGS_GMULT(MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(KIND=8)    :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  output from Boundary Value Problem
!  ----------------------------------

!  Column vectors for solving BCs

      REAL(KIND=8)    :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(KIND=8)    :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(KIND=8)    :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Solution constants of integration, and related quantities

      REAL(KIND=8)    :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(KIND=8)    :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Post-processing variables
!  -------------------------

!  BOA source terms

      REAL(KIND=8)    :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(KIND=8)    :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(KIND=8)    :: IDOWNSURF(MAXSTREAMS)

!  TOA source term

      REAL(KIND=8)    :: TOA_SOURCE(MAX_USER_STREAMS)

!  Quadrature-defined solutions

      REAL(KIND=8)    :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(KIND=8)    :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(KIND=8)    :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Local LInearized Arrays for argument passing, each Fourier or M = 0
!  ===================================================================

!  Linearized (Positive) Eigenvalues

      REAL(kind=8)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvector solutions

      REAL(kind=8)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=8)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(kind=8)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized norms for the Green function solution

      REAL(kind=8)  :: L_NORM_SAVED(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(kind=8)  :: L_HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(kind=8)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=8)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Particular beam solutions at user-defined stream angles

      REAL(kind=8)  :: LP_U_WPOS1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: LP_U_WNEG1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(kind=8)  :: LP_GAMMA_M (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: LP_GAMMA_P (MAXSTREAMS,MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=8)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=8)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=8)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized TOA and BOA source terms

      REAL(kind=8)  :: LP_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  :: LP_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(kind=8)  :: LP_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Profile weighting functions at quadrature angles

      REAL(kind=8)  :: QUADPROFILEWF ( MAX_ATMOSWFS, MAXLAYERS,  MAX_USER_LEVELS, &
                                       MAXSTREAMS,   MAXBEAMS,   MAX_DIRECTIONS )

!  Linearized Green's function multipliers for off-grid optical depths

      LOGICAL       :: FLAGS_LP_GMULT(MAX_PARTLAYERS)
      REAL(kind=8)  :: LP_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=8)  :: LP_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAXLAYERS,MAX_ATMOSWFS)

!  Local help variables
!  --------------------

!  Local direct beam reflectance

      LOGICAL         :: DO_LOCAL_DIRECTBEAM ( MAXBEAMS )

!  Indices

      INTEGER         :: LAYER, IBEAM

!  local inclusion flags

      LOGICAL         :: DO_INCLUDE_MVOUTPUT
      LOGICAL         :: DO_INCLUDE_DIRECTBEAM
      LOGICAL         :: DO_INCLUDE_SURFACE
      LOGICAL         :: DO_INCLUDE_SURFACEWF

!  Flux multiplier and Fourier component numbers

      REAL(KIND=8)    :: FLUX_MULTIPLIER
      REAL(KIND=8)    :: DELTA_FACTOR
      REAL(KIND=8)    :: SURFACE_FACTOR, ALBEDO

!  Local linearization control

      LOGICAL         :: DOVARY
      LOGICAL         :: DO_SIMULATION_ONLY
      INTEGER         :: N_PARAMETERS, LAYER_TO_VARY, N_LAYER_WFS

!  error tracing

      INTEGER         :: STATUS_SUB

!  progress

      logical, parameter :: do_write_screen = .false.
      character*2        :: CF

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Set local flags
!  ---------------

!  Local simulation only

      DO_SIMULATION_ONLY = ( .NOT. DO_PROFILE_LINEARIZATION .AND. .NOT. DO_SURFACE_LINEARIZATION )

!  Albedo (Lambertian case). Dark surface = Lambertian with albedo 0

      ALBEDO = LAMBERTIAN_ALBEDO(THREAD)

!  Surface flag (for inclusion of some kind of reflecting boundary)
!    shoudl be true for every component if the BRDF case 

      DO_INCLUDE_SURFACEWF = .TRUE.
      DO_INCLUDE_SURFACE   = .TRUE.
      IF ( .not. DO_BRDF_SURFACE ) THEN
        IF ( FOURIER_COMPONENT .NE. 0 ) THEN
          DO_INCLUDE_SURFACE   = .FALSE.
          DO_INCLUDE_SURFACEWF = .FALSE.
        ELSE
          IF ( ALBEDO .EQ. ZERO ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)

      DO IBEAM = 1, NBEAMS
        DO_LOCAL_DIRECTBEAM(IBEAM) = .FALSE.
      ENDDO
      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_LOCAL_DIRECTBEAM(IBEAM) = DO_REFLECTED_DIRECTBEAM(IBEAM)
        ENDDO
      ENDIF

!  surface reflectance factors

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Flux multipliers
!   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUX_MULTIPLIER   = DELTA_FACTOR
      
!  inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

!  Initialise BVP telescoping. Important.

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        DO_BVTEL_INITIAL = DO_BVP_TELESCOPING
      ENDIF

!  Reflected Direct beam attenuation

      CALL LIDORT_DIRECTBEAM                                &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,            & ! input
            DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS,        & ! input
            NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS,      & ! input
            FOURIER_COMPONENT, DELTA_FACTOR, FLUX_FACTOR,   & ! input
            BEAM_COSINES, SUNLAYER_COSINES,                 & ! input
            ALBEDO, BRDF_F_0, USER_BRDF_F_0,                & ! input
            SOLAR_BEAM_OPDEP, DO_LOCAL_DIRECTBEAM,          & ! input
            ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM )       ! Output

!  Get Legendre polynomials for this Fourier component

      CALL LIDORT_LEGENDRE_SETUP                                       &
         ( DO_REFRACTIVE_GEOMETRY, FOURIER_COMPONENT,                  & ! Input
           NSTREAMS, NBEAMS, NMOMENTS, NLAYERS,                        & ! Input
           BEAM_COSINES, SUNLAYER_COSINES, QUAD_STREAMS, QUAD_WEIGHTS, & ! Input
           PLMI_PLMJ_P, PLMI_PLMJ_M, PLMI_X0_P, PLMI_X0_M,             & ! Output
           LEG_P, LEG_M, LEG0_P, LEG0_M, WT_LEGP, WT_LEGM )              ! Output

      IF ( DO_USER_STREAMS ) THEN
        CALL LIDORT_USERLEGENDRE_SETUP         &
           ( N_USER_STREAMS, NMOMENTS,         & ! Input
             USER_STREAMS, FOURIER_COMPONENT,  & ! Input
             U_LEG_P, U_LEG_M )                  ! Output
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

!  Start layer loop

      DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer

        CALL LIDORT_HOM_SOLUTION                                  &
          ( DO_SOLUTION_SAVING, NSTREAMS, NMOMENTS,               & ! Input
            LAYER, FOURIER_COMPONENT, DO_LAYER_SCATTERING,        & ! Input
            OMEGA_MOMS, QUAD_STREAMS, QUAD_WEIGHTS,               & ! Input
            PLMI_PLMJ_P, PLMI_PLMJ_M,                             & ! Input
            SAB, DAB, EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,  & ! Output
            KEIGEN, XPOS, XNEG,                                   & ! Output
            STATUS_SUB, MESSAGE, TRACE_1 )                          ! Output                 

!  .. error tracing

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 =  'Called in LIDORT_LPS_FOURIER_MASTER, Fourier component '//CF
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Get Post-processing ("user") solutions for this layer

        IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
          IF ( DO_USER_STREAMS ) THEN
            CALL LIDORT_HOM_USERSOLUTION                        &
          ( NSTREAMS, N_USER_STREAMS, NMOMENTS,                 & ! Input
            LAYER, FOURIER_COMPONENT, DO_LAYER_SCATTERING,      & ! Input
            USER_SECANTS, KEIGEN, XPOS, XNEG,                   & ! Input
            OMEGA_MOMS, WT_LEGP, WT_LEGM, U_LEG_P,              & ! Input
            ZETA_M, ZETA_P, U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )  ! Output
          ENDIF
        ENDIF

!  Linearized solutions
!  --------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Parameter control

          DOVARY        = LAYER_VARY_FLAG(LAYER)
          N_PARAMETERS  = LAYER_VARY_NUMBER(LAYER)

!  Linearized discrete ordinate solutions

          CALL LIDORT_L_HOM_SOLUTION                      &
         ( DO_SOLUTION_SAVING, NSTREAMS, NMOMENTS,        & ! Input
           LAYER, FOURIER_COMPONENT, DO_LAYER_SCATTERING, & ! Input
           DOVARY, N_PARAMETERS, L_OMEGA_MOMS,            & ! Input
           QUAD_STREAMS, QUAD_WEIGHTS,                    & ! Input
           PLMI_PLMJ_P, PLMI_PLMJ_M, KEIGEN, SAB, DAB,    & ! Input
           EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,     & ! Input
           L_KEIGEN, L_XPOS, L_XNEG,                      & ! Input
           STATUS_SUB, MESSAGE, TRACE_1 )                   ! Output                 

!  .. error tracing

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER_COMPONENT
            TRACE_2 =  'Called in LIDORT_LPS_FOURIER_MASTER, Fourier component '//CF
            STATUS  = LIDORT_SERIOUS
            RETURN
          ENDIF

!  Get Linearizations of ("user") solutions for this layer

          CALL LIDORT_L_HOM_USERSOLUTION                    & 
         ( NSTREAMS, N_USER_STREAMS, NMOMENTS, LAYER,       & ! Input
           FOURIER_COMPONENT, DOVARY, N_PARAMETERS,         & ! Input
           STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
           DO_LAYER_SCATTERING, OMEGA_MOMS, L_OMEGA_MOMS,   & ! Input
           L_XPOS, L_XNEG, WT_LEGP, WT_LEGM,                & ! Input
           U_LEG_P, U_HELP_P, U_HELP_M,                     & ! Input
           L_U_XPOS, L_U_XNEG  )                              ! Output

!  End linearization clause

        ENDIF

!  end layer loop

      ENDDO

!  prepare eigenstream tranmsittances, and linearizations

      CALL LIDORT_HOM_EIGENTRANS                                         &
        ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,          & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          FOURIER_COMPONENT, DO_LAYER_SCATTERING,                        & ! Input
          DELTAU_VERT, PARTAU_VERT, KEIGEN,                              & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                & ! Input
          T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN )                   ! Output

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_L_HOM_EIGENTRANS                                                 &
       ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,                       & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,               & ! Input
         FOURIER_COMPONENT, DO_LAYER_SCATTERING, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! Input
         DELTAU_VERT,   PARTAU_VERT,  KEIGEN, L_DELTAU_VERT, L_KEIGEN,               & ! Input
         T_DELT_EIGEN,     T_UTUP_EIGEN,     T_UTDN_EIGEN,                           & ! Input
         L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,                       & ! Input
         L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN )                          ! Output

      ENDIF

!  Prepare solution norms, and linearizations

      CALL LIDORT_HOM_NORMS                          &
           ( NSTREAMS, NLAYERS, QUAD_STRMWTS, XPOS,  & ! Input
             NORM_SAVED )                              ! Output

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL LIDORT_L_HOM_NORMS                               &
          ( NSTREAMS, NLAYERS, QUAD_STRMWTS,                  & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, XPOS, L_XPOS, & ! Input
            L_NORM_SAVED )                                      ! Output
      ENDIF

!  Prepare homogeneous solution multipliers, and linearizations

      CALL HMULT_MASTER                                                 &
        ( DO_UPWELLING, DO_DNWELLING,                                   & ! Input
          NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,             & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input
          T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN,                 & ! Input
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                 & ! Input
          ZETA_M, ZETA_P,                                               & ! Input
          HMULT_1,        HMULT_2,                                      & ! Output
          UT_HMULT_UU, UT_HMULT_UD,                                     & ! Output
          UT_HMULT_DU, UT_HMULT_DD )                                      ! Output

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL L_HMULT_MASTER                                             &
       ( DO_UPWELLING, DO_DNWELLING,                                    & ! Input
         NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,              & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input 
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, ZETA_M, ZETA_P, L_KEIGEN,  & ! Input
         T_DELT_EIGEN, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Input
         HMULT_1, UT_HMULT_UU, UT_HMULT_UD,                             & ! Input
         HMULT_2, UT_HMULT_DU, UT_HMULT_DD,                             & ! Input
         L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,            & ! Input
         L_T_DELT_USERM,   L_T_UTUP_USERM,   L_T_UTDN_USERM,            & ! Input
         L_HMULT_1, L_UT_HMULT_UU, L_UT_HMULT_UD,                       & ! Output
         L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD )                        ! Output
      ENDIF

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!      write(*,*)'bvpmatrix',BVP_REGULAR_FLAG(FOURIER_COMPONENT)

!  standard case using compression of band matrices, etc..

      IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

        CALL BVP_MATRIXSETUP_MASTER                                            &
    ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER_COMPONENT,                  & ! Inputs
      NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! Inputs
      QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F, XPOS, XNEG, T_DELT_EIGEN,  & ! Inputs
      H_XPOS, H_XNEG, LCONMASK, MCONMASK,                                      & ! Output
      BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                          & ! Output
      STATUS_SUB, MESSAGE, TRACE_1 )                                             ! Output

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error from BVP_MATRIXSETUP_MASTER, '//   &
               'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
          STATUS = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Telescoped case

      ELSE

        CALL BVPTEL_MATRIXSETUP_MASTER                                &
           ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS,                   & ! Input
             NSTREAMS_2, N_SUPDIAG, N_SUBDIAG, FOURIER_COMPONENT,     & ! Input
             DO_LAYER_SCATTERING, XPOS, XNEG, T_DELT_EIGEN,           & ! Input
             DO_BVTEL_INITIAL, NLAYERS_TEL,                           & ! Output
             ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,                       & ! Output
             BTELMAT_ROWMASK, BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, & ! Output
             STATUS_SUB, MESSAGE, TRACE_1 )                             ! Output

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '// &
           'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
          STATUS = LIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Only calculate if still not converged

        IF ( DO_MULTIBEAM(IBEAM,FOURIER_COMPONENT) ) THEN

!  Solar beam Particular solutions (Green's function)
!  --------------------------------------------------

!  start layer loop

          DO LAYER = 1, NLAYERS

!  Control
            
            IF ( DO_PROFILE_LINEARIZATION ) THEN
              DOVARY        = LAYER_VARY_FLAG(LAYER)
              N_PARAMETERS  = LAYER_VARY_NUMBER(LAYER)
            ENDIF

!  Green's function solution

            CALL LIDORT_GBEAM_SOLUTION                                  &
          ( NSTREAMS, NSTREAMS_2, NMOMENTS, LAYER, FOURIER_COMPONENT,   & ! Input
            FLUX_FACTOR, IBEAM, DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,  & ! Input
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                & ! Input
            QUAD_WEIGHTS, OMEGA_MOMS, KEIGEN, XPOS, T_DELT_EIGEN,       & ! Input
            NORM_SAVED, PLMI_X0_P, PLMI_X0_M,                           & ! Input
            GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE,         & ! Output
            CFUNC, DFUNC, AGM, BGP, GFUNC_UP, GFUNC_DN,                 & ! Output
            WUPPER, WLOWER )                                              ! Output

!  Linearized Green's function solution

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              CALL LIDORT_L_GBEAM_SOLUTION                   &
         ( NSTREAMS, NMOMENTS, LAYER, FOURIER_COMPONENT,     & ! Input
           FLUX_FACTOR, IBEAM, DOVARY, N_PARAMETERS,         & ! Input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,            & ! Input
           QUAD_WEIGHTS, L_OMEGA_MOMS, PLMI_X0_P, PLMI_X0_M, & ! Input
           NORM_SAVED, XPOS, L_NORM_SAVED, L_XPOS,           & ! Input
           DMI, DPI, ATERM_SAVE, BTERM_SAVE,                 & ! Input
           L_ATERM_SAVE, L_BTERM_SAVE )                        ! Output
            ENDIF

!  Post-processing solutions

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN

              IF ( DO_USER_STREAMS ) THEN

                CALL LIDORT_GBEAM_USERSOLUTION                    &
          ( DO_UPWELLING, DO_DNWELLING, N_USER_STREAMS, NMOMENTS, & ! Input
            LAYER, FOURIER_COMPONENT, IBEAM, FLUX_FACTOR,         & ! Input
            DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,                & ! Input
            OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,                 & ! Input
            U_WPOS1, U_WNEG1, W_HELP )                              ! Output

                IF ( DO_PROFILE_LINEARIZATION ) THEN
                  CALL LIDORT_L_GBEAM_USERSOLUTION             &
         ( DO_UPWELLING, DO_DNWELLING,                         & ! Input
           N_USER_STREAMS, NMOMENTS, LAYER, FOURIER_COMPONENT, & ! Input
           IBEAM, FLUX_FACTOR, DOVARY, N_PARAMETERS,           & ! Input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,              & ! Input
           L_OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,             & ! Input
           LP_U_WPOS1, LP_U_WNEG1 )                              ! Output
                ENDIF

              END IF

            END IF

!  end layer loop

          END DO

!  Solve boundary value problem
!  ----------------------------

!       if ( do_write_screen) write(*,*)'bvp solution',ibeam

!  Get the Regular BVP solution

          IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

           CALL BVP_SOLUTION_MASTER                                                 &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_LOCAL_DIRECTBEAM(IBEAM), & ! Input
            NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,     & ! Input
            FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR, ALBEDO, BRDF_F,        & ! Input
            QUAD_STRMWTS, XPOS, XNEG, WUPPER, WLOWER, DIRECT_BEAM,           & ! Input
            BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                & ! Input
            H_WLOWER, COL2, SCOL2, LCON, MCON,                               & ! Output
            LCON_XVEC, MCON_XVEC, STATUS_SUB, MESSAGE, TRACE_1 )               ! Output

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             write(CF,'(I2)')FOURIER_COMPONENT
             TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
                'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
             STATUS = LIDORT_SERIOUS
             RETURN
           ENDIF 

!  Get the telescoped boundary value result

          ELSE

           CALL BVPTEL_SOLUTION_MASTER                            &
             ( IBEAM, NSTREAMS, NLAYERS,                          & ! Input
               NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, WUPPER, WLOWER,  & ! Input
               XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
               NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
               BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,            & ! Input
               COLTEL2, SCOL2, LCON, MCON,                        & ! Output
               LCON_XVEC,  MCON_XVEC,                             & ! Output
               STATUS_SUB, MESSAGE, TRACE_1 )                       ! Output

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             write(CF,'(I2)')FOURIER_COMPONENT
             TRACE_2 = 'Error return from BVPTEL_SOLUTION_MASTER, '// &
                'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
             STATUS = LIDORT_SERIOUS
             RETURN
           ENDIF

          ENDIF

!  Radiance Field Post Processing
!  ------------------------------

!  Direct beam inclusion flag:
!   This now has the DBCORRECTION option: if the DBCORRECTION Flag
!   is set, then we will be doing exact calculations of the reflected
!   directbeam, so we do not need to include it in the Post-processing.
!   However, the direct beam will need to be included in the basi! RT
!   solution (the BVP), and this is controlled separately by the
!   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
!     R. Spurr, RT Solutions, Inc., 19 August 2005.

!          DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
!             (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))

          DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. DO_LOCAL_DIRECTBEAM(IBEAM) ) &
                              .and. .not. DO_MSMODE_LIDORT

!  upwelling

          IF ( DO_UPWELLING ) THEN

            CALL GET_BOASOURCE                                &             
          ( DO_USER_STREAMS, DO_INCLUDE_SURFACE,              & ! Input
            DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! Input
            NSTREAMS, NLAYERS, N_USER_STREAMS, IBEAM,         & ! Input
            FOURIER_COMPONENT, SURFACE_FACTOR, QUAD_STRMWTS,  & ! Input
            LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER,       & ! Input
            ALBEDO, USER_BRDF_F, USER_DIRECT_BEAM,            & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE, IDOWNSURF )          ! Output

            CALL UPUSER_INTENSITY                                     &
            ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,    & ! Input
              DO_INCLUDE_MVOUTPUT, FOURIER_COMPONENT,                 & ! Input
              NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,       & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                & ! Input 
              UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                & ! Input
              IBEAM, FLUX_MULTIPLIER,                                 & ! Input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,   & ! Input
              T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,             & ! Input
              T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM, & ! Input
              GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input   
              XPOS, WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC, & ! Input          
              U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,    & ! Input               
              UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                 & ! Input      
              BOA_SOURCE, DIRECT_BOA_SOURCE,                          & ! Input                      
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,           & ! Output                
              FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                  & ! Output
              INTENSITY_F, QUADINTENS, CUMSOURCE_UP )                   ! Output

          ENDIF

!  Downwelling

          IF ( DO_DNWELLING ) THEN

            CALL GET_TOASOURCE ( N_USER_STREAMS, TOA_SOURCE )

            CALL DNUSER_INTENSITY                                      &
           ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,      & ! input
             DO_INCLUDE_MVOUTPUT, FOURIER_COMPONENT,                   & ! input
             NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                  & ! input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                  & ! input
             IBEAM, FLUX_MULTIPLIER,                                   & ! input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,     & ! input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,               & ! input
             T_DELT_MUBAR,  T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM,  & ! input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                 & ! input
             XPOS, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,           & ! input
             U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,      & ! input
             UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN, TOA_SOURCE,       & ! input
             PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,             & ! Output
             FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                    & ! Output
             INTENSITY_F, QUADINTENS, CUMSOURCE_DN )                     ! Output

          ENDIF

!  mean value output

          IF ( DO_INCLUDE_MVOUTPUT ) THEN

            CALL MIFLUX_INTENSITY                                     &
           ( DO_UPWELLING, DO_DNWELLING, DO_LOCAL_DIRECTBEAM(IBEAM),  & ! input
             THREAD, IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,     & ! input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                 & ! input
             QUAD_WEIGHTS, QUAD_STRMWTS,                              & ! input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,             & ! input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,                  & ! input
             MEAN_INTENSITY, FLUX_INTEGRAL )                            ! Output

          ENDIF

!  Finish this Beam solution, if only intensity is required

          IF ( .not. DO_SIMULATION_ONLY ) THEN

!  Profile weighting functions
!  ===========================

          IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

           DO LAYER_TO_VARY = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
             N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

!  4A. Solve the linearized BVP
!  ----------------------------

!  Regular BVP linearization

             IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

!  Get the linearized BVP solution for this beam component

              CALL LP_BVP_SOLUTION_MASTER                               & 
            ( DO_INCLUDE_SURFACE, DO_LOCAL_DIRECTBEAM(IBEAM),           & ! Input
              DO_BRDF_SURFACE,DO_PLANE_PARALLEL, NLAYERS, NSTREAMS,     & ! Input
              NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,                 & ! Input
              LAYER_TO_VARY, N_LAYER_WFS, FOURIER_COMPONENT,            & ! Input
              IBEAM, QUAD_STRMWTS, SURFACE_FACTOR,                      & ! Input
              ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,                & ! Input
              DO_LAYER_SCATTERING, INITIAL_TRANS, LAYER_PIS_CUTOFF,     & ! Input
              T_DELT_EIGEN, T_DELT_MUBAR, XPOS, XNEG,                   & ! Input
              GAMMA_M, GAMMA_P, GFUNC_UP, GFUNC_DN,                     & ! Input
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                         & ! Input
              LCON, MCON, LCONMASK, MCONMASK,                           & ! Input
              L_DELTAU_VERT, L_ATERM_SAVE, L_BTERM_SAVE,                & ! Input
              L_KEIGEN, L_XPOS, L_XNEG, L_T_DELT_EIGEN,                 & ! Input
              LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,     & ! Input
              LP_GAMMA_M, LP_GAMMA_P, L_WUPPER, L_WLOWER,               & ! Output
              NCON, PCON, NCON_XVEC, PCON_XVEC,                         & ! Output
              STATUS_SUB, MESSAGE, TRACE_1 )                              ! Output

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from LP_BVP_SOLUTION_MASTER, '// &
                   'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
                STATUS = LIDORT_SERIOUS
                RETURN
              ENDIF

!  telescoped BVP linearization

             ELSE

              CALL LP_BVPTEL_SOLUTION_MASTER                      &
          ( DO_PLANE_PARALLEL, NSTREAMS, NSTREAMS_2,              & ! Input
            NLAYERS, N_SUPDIAG, N_SUBDIAG,                        & ! Input
            LAYER_TO_VARY, N_LAYER_WFS, FOURIER_COMPONENT, IBEAM, & ! Input
            ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,       & ! Input
            DO_LAYER_SCATTERING, INITIAL_TRANS, LAYER_PIS_CUTOFF, & ! Input
            T_DELT_DISORDS,  T_DELT_EIGEN, T_DELT_MUBAR,          & ! Input
            XPOS, XNEG, GAMMA_M, GAMMA_P, GFUNC_UP, GFUNC_DN,     & ! Input
            BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,               & ! Input
            LCON, MCON, LCON_XVEC, MCON_XVEC,                     & ! Input
            L_T_DELT_DISORDS, L_ATERM_SAVE, L_BTERM_SAVE,         & ! Input
            L_KEIGEN, L_XPOS, L_XNEG, L_T_DELT_EIGEN,             & ! Input
            LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input
            LP_GAMMA_M, LP_GAMMA_P, L_WUPPER, L_WLOWER,           & ! Output
            NCON, PCON, NCON_XVEC, PCON_XVEC,                     & ! Output
            STATUS_SUB, MESSAGE, TRACE_1 )                          ! Output

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from LP_BVPTEL_SOLUTION_MASTER, '// &
                   'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
                STATUS = LIDORT_SERIOUS
                RETURN
              ENDIF

             ENDIF

!  Post-processing for Column weighting functions (upwelling)
!  ----------------------------------------------------------

             IF ( DO_UPWELLING ) THEN

!  Linearized BOA source term

              CALL GET_LP_BOASOURCE                                         &
    ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_BRDF_SURFACE,           & ! Input
      DO_USER_STREAMS, NLAYERS, NSTREAMS, N_USER_STREAMS,                   & ! Input
      FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS, SURFACE_FACTOR, & ! Input
      ALBEDO, USER_BRDF_F, QUAD_STRMWTS, USER_DIRECT_BEAM, DELTAU_SLANT,    & ! Input
      L_DELTAU_VERT, L_XPOS, L_XNEG, LCON, LCON_XVEC, T_DELT_EIGEN,         & ! Input
      MCON, NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_WLOWER,                 & ! Input
      LP_BOA_MSSOURCE, LP_BOA_DBSOURCE )                                      ! Output

!  Linearized upwelling weighting functions

               CALL UPUSER_PROFILEWF                                   &
           ( DO_USER_STREAMS,  DO_PLANE_PARALLEL,                      & ! Input
             DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER,   & ! Input
             NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
             UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                  & ! Input
             FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS,     & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,     & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,               & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM,   & ! Input
             AGM, BGP, XPOS, LCON, LCON_XVEC, MCON, MCON_XVEC,         & ! Input
             U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,      & ! Input
             UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                   & ! Input
             PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,             & ! Input
             UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                   & ! Input
             L_T_DELT_USERM, L_T_UTUP_USERM, HELP_AQ, HELP_BQ,         & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,       & ! Input
             L_T_DELT_EIGEN, L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,         & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WPOS1, NCON, PCON,               & ! Input
             L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC, & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                        & ! Input
             L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,             & ! Input
             LP_BOA_MSSOURCE, LP_BOA_DBSOURCE,                         & ! Input
             FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,           & ! Output
             PROFILEWF_F, QUADPROFILEWF )                                ! Output

             ENDIF

!  Post-processing for Column weighting functions (Downwelling)
!  ------------------------------------------------------------

             IF ( DO_DNWELLING ) THEN

!  TOA source term linearized

              CALL GET_LP_TOASOURCE ( N_USER_STREAMS, LP_TOA_SOURCE, N_LAYER_WFS )

!  Downwelling weighting functions

              CALL DNUSER_PROFILEWF                                   &
           ( DO_USER_STREAMS,  DO_PLANE_PARALLEL,                     & ! Input
             DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT,                   & ! Input
             NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                 & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                 & ! Input
             FOURIER_COMPONENT,                                       & ! Input
             IBEAM, LAYER_TO_VARY, N_LAYER_WFS, FLUX_MULTIPLIER,      & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,    & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,              & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM,  & ! Input
             AGM, BGP, XPOS, LCON, LCON_XVEC, MCON, MCON_XVEC,        & ! Input
             U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,     & ! Input
             UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                  & ! Input
             PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,            & ! Input
             UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                  & ! Input
             L_T_DELT_USERM, L_T_UTDN_USERM, HELP_AQ, HELP_BQ,        & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,      & ! Input
             L_T_DELT_EIGEN, L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,        & ! Input
             LP_GAMMA_M, LP_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,      & ! Input
             L_U_XPOS, L_U_XNEG, LP_U_WNEG1, NCON, PCON,              & ! Input
             L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,          & ! Input
             L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                       & ! Input
             L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN,            & ! Input
             LP_TOA_SOURCE,                                           & ! Input
             FLAGS_LP_GMULT, LP_UT_GMULT_UP, LP_UT_GMULT_DN,          & ! Input
             PROFILEWF_F, QUADPROFILEWF )                               ! Output

             ENDIF

!  Post-processing for Profile weighting functions (Mean/Flux)
!  -----------------------------------------------------------

             IF ( DO_INCLUDE_MVOUTPUT ) THEN
              CALL MIFLUX_PROFILEWF                                   &
           ( DO_UPWELLING, DO_DNWELLING, DO_LOCAL_DIRECTBEAM(IBEAM),  & ! Input
             THREAD, IBEAM, LAYER_TO_VARY, N_LAYER_WFS,               & ! Input
             NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,                    & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                 & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                              & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,             & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADPROFILEWF,               & ! Input
             LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,      & ! Input
             MINT_PROFILEWF, FLUX_PROFILEWF )                           ! Output

             ENDIF

!  Finish loop over layers with variation

            ENDIF
           ENDDO

!  End atmospheric profile weighting functions

          ENDIF

!   Surface Property Jacobians (will be done regardless of Dark surface condition)
!  ===========================

          IF ( DO_SURFACE_LINEARIZATION ) THEN
           IF ( DO_INCLUDE_SURFACEWF ) THEN

!  Surface WF master

            CALL SURFACEWF_MASTER                                               &
          ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_STREAMS,       & ! Input
            DO_LOCAL_DIRECTBEAM(IBEAM), DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT,  & ! Input
            NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_USER_LEVELS,               & ! Input
            NSTREAMS, NSTREAMS_2, N_USER_STREAMS, N_SURFACE_WFS,                & ! Input
            PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                             & ! Input
            FOURIER_COMPONENT, IBEAM, THREAD, FLUX_MULTIPLIER,                  & ! Input
            SURFACE_FACTOR, ALBEDO, LS_BRDF_F, LS_BRDF_F_0,                     & ! Input
            USER_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0,                      & ! Input
            QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF,                  & ! Input
            LCONMASK, MCONMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,               & ! Input
            LCON, MCON, XPOS, XNEG, H_XPOS, H_XNEG, H_WLOWER,                   & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                           & ! Input
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                           & ! Input
            U_XPOS, U_XNEG, HMULT_1, HMULT_2,                                   & ! Input
            UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,                 & ! Input
            SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                        & ! Output
            STATUS_SUB, MESSAGE, TRACE_1 )                                        ! Output

!  Error handling

             IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER_COMPONENT
               TRACE_2 = 'Error return from SURFACEWF_MASTER, '// &
                  'Called in LIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS
               RETURN
             ENDIF

!  end of Surface weighting functions

           ENDIF
          ENDIF

!  End complete weighting function clause (NON simulation-only!)

          ENDIF

!  End loop over beam solutions

        END IF
      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE LIDORT_LPS_FOURIER_MASTER

