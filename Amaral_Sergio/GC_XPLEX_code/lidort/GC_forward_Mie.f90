!$Id: GC_forward_Mie.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
subroutine GC_forward_Mie &
     ( NPM_MAX, NCONC, DIAMETER, VARIANCE, REFINDEX, WAVELENGTH, & ! input
       BEXT, SSA, PMOMS, FAIL, MESSAGES )                          ! Output

!  Mie modules. Constants only.

  use Mie_precision
  use Mie_constants

!  Implicit none

  implicit none

!  Input arguments
!  ===============

   integer        , intent(in) :: NPM_MAX       !  Maximum number of phase expansion coeffs

   real(kind=8)   , intent(in) :: NCONC         !  Number concentration in mol/cm3
   real(kind=8)   , intent(in) :: DIAMETER      !  Mode diameter in Microns
   real(kind=8)   , intent(in) :: VARIANCE      !  Mode  variance

   complex(kind=8), intent(in) :: REFINDEX      !  Refractive index
   real(kind=8)   , intent(in) :: WAVELENGTH    !  wavelength in Microns

!  output
!  ======

!  optical properties

   real(kind=8)   , intent(out) :: BEXT             ! Extinction coefficient in mm^-1
   real(kind=8)   , intent(out) :: SSA              ! Single scattering albedo
   real(kind=8)   , intent(out) :: PMOMS(0:NPM_MAX) ! phase function expansion coeffs.

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: messages(4)

!  Mie code variables
!  ==================

!  Dimensioning input for Mie code

   INTEGER, PARAMETER :: max_Mie_angles     = 700
   INTEGER, PARAMETER :: max_Mie_sizes      = 20
   INTEGER, PARAMETER :: max_Mie_points     = 350
   INTEGER, PARAMETER :: max_Mie_distpoints = 20

!  This is an extreme case. Use with Caution - can cause Segmentation errors
!     Needs lots of memory

!  INTEGER, PARAMETER :: max_Mie_angles     = 5441
!  INTEGER, PARAMETER :: max_Mie_sizes      = 20
!  INTEGER, PARAMETER :: max_Mie_points     = 2720
!  INTEGER, PARAMETER :: max_Mie_distpoints = 20

!  Distribution parameters

   INTEGER            :: mie_idis
   REAL    (KIND=dp)  :: mie_pars(3)

!  Numerical control of PSD

   INTEGER            :: mie_nblocks
   INTEGER            :: mie_nweights
   REAL    (KIND=dp)  :: mie_cutoff
   REAL    (KIND=dp)  :: xparticle_limit
   COMPLEX (KIND=dp)  :: m_complex

!  general control variables

   LOGICAL           :: do_use_cutoff
   LOGICAL           :: do_angular_variation, do_bulk_only
   LOGICAL           :: do_external_angles, do_coeffct_angles

!  Angles

   INTEGER           :: n_coeffct_angles
   REAL    (KIND=dp) :: coeff_cosines(max_Mie_angles)
   REAL    (KIND=dp) :: coeff_weights(max_Mie_angles)
   INTEGER           :: n_external_angles
   REAL    (KIND=dp) :: external_angle_cosines(max_Mie_angles)
   REAL    (KIND=dp) :: mie_angle_cosines(max_Mie_angles)

!  Wavelength

   REAL    (KIND=dp)  :: mie_wavelength

!  Maximum and Minimum Radius

   REAL    (KIND=dp)  :: mie_rmax, mie_rmin

!  Full Mie output, wavelength-saved

   REAL    (KIND=dp)  :: Mie_bulk(4)
   REAL    (KIND=dp)  :: Mie_expcoeffs(6,0:max_Mie_angles)

!  5 Distribution parameters

   REAL    (KIND=dp)  :: MIE_DIST(5)

!  4 F-matrix values at user-defined or coefficient angles

   REAL    (KIND=dp)  :: MIE_FMAT(4,max_Mie_angles)

!  Other local variables
!  =====================

   integer      :: L
   logical      :: startup
   real(kind=8) :: norm

!  initialize output
!  =================

   BEXT  = 0.0d0
   SSA   = 0.0d0
   PMOMS = 0.0d0

   FAIL          = .false.
   MESSAGES(1:4) = ' '

!  Set the Mie program control inputs
!  ==================================

!  PSD quadrature

   mie_nblocks  = 5
   mie_nweights = 20

!  particle limit

   xparticle_limit = 2000.0d0
   ! dkh 
   !xparticle_limit = 10000.0d0

!  Cutoff control

   ! dkh 
   !mie_cutoff    = 1.0d-8
   mie_cutoff    = 1.0d-5
   do_use_cutoff = .true.
   do_use_cutoff = .false.
   mie_rmin = 0.001d0
   mie_rmax = 2.0d0

!  Flags

   do_external_angles   = .false.
   do_coeffct_angles    = .true.
   do_angular_variation = .true.
   do_bulk_only         = .false.

!  initialize angle inputs

   n_external_angles = 0
   mie_angle_cosines = 0.0d0

!  log normal choices

   mie_idis    = 4
   mie_pars(1) = 0.5d0 * DIAMETER
!   mie_pars(2) = DEXP(VARIANCE)
   mie_pars(2) = VARIANCE
   mie_pars(3) = 0.0d0

!  Wavelength and refractive index

   m_complex      = REFINDEX
   mie_wavelength = WAVELENGTH

!  Set start-up flag

   startup = .true.

!  MIE Operation
!  =============

!  Call to the Mie code

   call Mie_main                                                               & !-------MIE CALL
     (   max_Mie_angles, max_Mie_sizes, max_Mie_points, max_Mie_distpoints,    & ! Dimensioning
         do_external_angles, do_coeffct_angles,  do_use_cutoff,                & ! I
         mie_idis, mie_pars, startup, mie_nblocks, mie_nweights, mie_cutoff,   & ! I
         n_external_angles, external_angle_cosines,                            & ! I
         n_coeffct_angles, coeff_cosines, coeff_weights,                       & ! I
         m_complex, xparticle_limit, mie_wavelength, mie_rmax, mie_rmin,       & ! I
         Mie_bulk, MIE_dist, MIE_fmat,                                         & ! O
         messages(1), messages(2), messages(3), fail )                           ! O

!  Exception handling

   if ( fail ) then
      messages(4) = 'Failure from the Mie call'
      return
   endif

!  Develop the phase function coefficients

   CALL develop                                                             & !--------DEVELOP CALL
      ( max_Mie_angles, n_coeffct_angles, n_coeffct_angles,                 & ! I
        coeff_cosines, coeff_weights, MIE_fmat, Mie_expcoeffs )              ! I/O

!  Output Interpretation
!  =====================

!  Extinction coefficient --
!    Multiply the extinction efficiency by the geometric cross-section
!    Multiply by the number concentration normalization

   NORM = NCONC / MIE_DIST(1)
   BEXT = NORM * MIE_BULK(1) * MIE_DIST(2) 

!  Singel scattering albedo

   SSA = MIE_BULK(4)

   !write(*,*)'dist 1 = ', MIE_DIST(1)
   !write(*,*)'dist 2 = ', MIE_DIST(2)

!  Preserve only the first 20 expansion coefficients

   PMOMS(0) = 1.0d0
   do L = 1, NPM_MAX
      PMOMS(L) = Mie_expcoeffs(1,L)
   enddo

!  Finish

   RETURN
end subroutine GC_forward_Mie
