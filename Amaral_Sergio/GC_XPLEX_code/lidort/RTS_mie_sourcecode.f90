!$Id: RTS_mie_sourcecode.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
SUBROUTINE Mie_main & 
     (   max_Mie_angles, max_Mie_sizes,                             & ! D
         max_Mie_points, max_Mie_distpoints,                        & ! D
         do_external_angles, do_coeffct_angles, do_use_cutoff,      & ! I     
         idis, nr_parameters, startup,                              & ! I
         nblocks, nweights, cutoff,                                 & ! I
         n_external_angles, external_angle_cosines,                 & ! I
         n_coeffct_angles, coeff_cosines, coeff_weights,            & ! I
         m_complex, xparticle_limit, wavelength, rmax, rmin,        & ! I
         mie_bulk, dist, fmat,                                      & ! O
         message, trace, action, failmie )                            ! O

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three, c_i

!  implicit none statement

  IMPLICIT NONE

!  Dimensioning input

  INTEGER, INTENT (IN) :: max_Mie_angles, max_Mie_sizes
  INTEGER, INTENT (IN) :: max_Mie_points, max_Mie_distpoints

!  input

  LOGICAL          , INTENT (IN) :: do_external_angles
  LOGICAL          , INTENT (IN) :: do_coeffct_angles
  LOGICAL          , INTENT (IN) :: do_use_cutoff

  INTEGER          , INTENT (IN) :: idis
  REAL    (KIND=dp), INTENT (IN) :: nr_parameters(3)

  INTEGER          , INTENT (IN) :: nblocks
  INTEGER          , INTENT (IN) :: nweights
  REAL    (KIND=dp), INTENT (IN) :: cutoff

  INTEGER          , INTENT (IN) :: n_external_angles
  REAL    (KIND=dp), INTENT (IN) :: external_angle_cosines(max_Mie_angles)

  LOGICAL          , INTENT (INOUT) :: startup
  INTEGER          , INTENT (INOUT) :: n_coeffct_angles
  REAL    (KIND=dp), INTENT (INOUT) :: coeff_cosines(max_Mie_angles)
  REAL    (KIND=dp), INTENT (INOUT) :: coeff_weights(max_Mie_angles)

  COMPLEX (KIND=dp), INTENT (IN) :: m_complex
  REAL    (KIND=dp), INTENT (IN) :: xparticle_limit
  REAL    (KIND=dp), INTENT (IN) :: wavelength

!  output

  REAL    (KIND=dp), INTENT (OUT)   :: MIE_BULK(4)
  REAL    (KIND=dp), INTENT (OUT)   :: DIST(5)
  REAL    (KIND=dp), INTENT (OUT)   :: FMAT(4,max_Mie_angles)

  LOGICAL          , INTENT (OUT)   :: failmie
  CHARACTER*(*)    , INTENT (OUT)   :: message, trace, action
  REAL    (KIND=dp), INTENT (INOUT) :: rmin, rmax

!  Local Mie output 

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes)                  :: q_ext, q_sca, asym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles, max_Mie_sizes)  :: splus,sminus

!  local variables for Mie code

  CHARACTER*5 :: char5
  LOGICAL     :: do_angular_variation
  LOGICAL     :: failmm, faild
  INTEGER     :: i, angle, n_angles
  INTEGER     :: iblock, n_sizes, kf

  REAL    (KIND=dp)  :: factor_0, factor_1, d_pi
  REAL    (KIND=dp)  :: rstart, rfinis, help, rblock
  REAL    (KIND=dp)  :: quad, quadr2, quadr3, quadr4
  REAL    (KIND=dp)  :: ndens, gxsec, reff, volume, veff
  REAL    (KIND=dp)  :: Qext, Qsca, Qasy, ssalbedo
  REAL    (KIND=dp)  :: f(4)
  REAL    (KIND=dp)  :: angle_cosines(max_Mie_angles)

!  redundant variables
!  REAL    (KIND=dp)  :: xeff, xeff_d(3)
!  LOGICAL            :: fail
 
  COMPLEX (KIND=dp)  :: sp, sm, csp, csm
  COMPLEX (KIND=dp)  :: c_mi

  REAL    (KIND=dp)  :: xpart_root3, xparticle
  COMPLEX (KIND=dp)  :: y_argument
  INTEGER            :: limmax, limstop, limsize

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes)   :: particle_sizes
  REAL    (KIND=dp), DIMENSION (max_Mie_sizes)   :: rquad, weights, nr

!  start up
!  --------

  c_mi = - c_i

  n_sizes = nweights

  d_pi = 4.0_dp * ATAN(d_one)
  factor_0 = d_two * d_pi / wavelength
  factor_1 = wavelength / factor_0

!  Zeroing
!  -------

  trace   = ' '
  message = ' '
  action  = ' '
  failmie = .FALSE.
  failmm = .false.

  Qext = d_zero
  Qsca = d_zero
  Qasy = d_zero

  ndens = d_zero
  gxsec = d_zero
  reff  = d_zero
  veff  = d_zero

!  limiting radii calculation

  if ( do_use_cutoff ) then
    CALL rminmax ( idis, nr_parameters, cutoff, rmin, rmax, message, failmm )
    IF ( failmm ) THEN
       failmie = .TRUE.
       trace   = 'Trace   : First Check in Mie Main. Failed to find radii extrema'
       action  = 'Action  : Consult with R. Spurr' 
       RETURN
    END IF
  endif

!  Check limiting radii if set externally

  if ( .not. do_use_cutoff ) then
    if ( rmin.lt.0.0d0 ) then
      failmm = .true.
      message = 'External Rmin value < 0, out of bounds'
    else if ( rmax .le. 0.0d0 ) then
      failmm = .true.
      message = 'External Rmax value =< 0, out of bounds'
    else if ( rmin .ge. rmax ) then
      failmm = .true.
      message = 'External Rmin >= Rmax, Cannot be possible!'
    endif
    if ( failmm ) then
       trace   = 'Trace   : First Check in Mie Main. User Rmin/Rmax wrong'
       action  = 'Action  : Change input values of Rmin and Rmax' 
       RETURN
    END IF
  endif    

!  number of blocks

  rblock = ( rmax - rmin ) / DBLE(nblocks)

!  limiting number of terms for coefficient computation

  xparticle  = factor_0  * rmax
  y_argument = xparticle * m_complex
  limstop    = 2
  IF ( xparticle > 0.02) THEN
     xpart_root3 = xparticle ** ( d_one / d_three )
     IF ( xparticle <= 8.0_dp ) THEN
        limstop = xparticle + 4.0_dp * xpart_root3 + d_two
     ELSE IF ( xparticle < 4200.0_dp ) THEN
        limstop = xparticle + 4.05_dp * xpart_root3 + d_two
     ELSE
        limstop = xparticle + 4.0_dp * xpart_root3 + d_two
     END IF
  END IF
  limmax = nint(max(DBLE(limstop),ABS(y_argument)) + 15.0_dp)

!  Dimensioning and exception handling checks
!  ------------------------------------------

!  return if size limit exceeded

  IF ( xparticle > xparticle_limit ) THEN
     failmie = .TRUE.
     limsize = int(xparticle) + 1
     write(char5,'(i5)')limsize
     message = 'Message : error size parameter overflow'
     trace   = 'Trace   : Second check in Mie_main'
     action  = 'Action  : In configuration file, increase cutoff or '// &
                          'increase xparticle_limit to at least '//char5
     RETURN
  END IF

!  return if maximum number of terms too great 

  IF ( limstop > max_Mie_points )  THEN
     failmie = .TRUE.
     write(char5,'(i5)')limstop
     message = 'Message : Insufficient dimensioning for maximum number of terms'
     trace   = 'Trace   : Third check in Mie_main'
     action  = 'Action  : Increase max_Mie_points in calling program to at least '//char5
     RETURN
  END IF

!  And again, Dave recurrence

  IF ( limmax > max_Mie_points )  THEN
     failmie = .TRUE.
     write(char5,'(i5)')limmax
     message = 'Message : Insufficient dimensioning for maximum number of terms (Dave recurrence)'
     trace   = 'Trace   : Fourth check in Mie_main'
     action  = 'Action  : Increase max_Mie_points in calling program to at least '//char5
     RETURN
  END IF

!  Compute the number of angles required for coefficient computation

  IF ( do_coeffct_angles .AND. startup ) THEN
     n_coeffct_angles = 2*limmax + 2
     if ( n_coeffct_angles > max_Mie_angles ) then
       failmie = .true.
       write(char5,'(i5)')n_coeffct_angles
       message = 'Message : Dimensioning error for number of terms for coefficient computation'
       trace   = 'Trace   : Fifth check in Mie Main'
       action  = 'Action  : Increase value of max_Mie_angles in calling program to at least '//char5
       return
     endif
  ENDIF

!  Compute the angles required for coefficient computation
!  -------------------------------------------------------

!  Quadrature (only need to do it once)

  IF ( do_coeffct_angles .AND. startup ) THEN
     n_coeffct_angles = 2*limmax + 2
     CALL mie_gauleg ( max_Mie_angles, n_coeffct_angles, -1.0_dp, 1.0_dp, & ! Input
                         coeff_cosines, coeff_weights )                     ! Output
     startup = .FALSE.
  END IF

! Overall cosines

  IF ( do_coeffct_angles ) THEN
     n_angles = n_coeffct_angles
     DO angle = 1, n_angles
        angle_cosines(angle) = coeff_cosines(angle)
     END DO
     do_angular_variation = .TRUE.
  ELSE IF ( do_external_angles ) THEN
     n_angles = n_external_angles
     DO angle = 1, n_angles
        angle_cosines(angle) = external_angle_cosines(angle)
     END DO
     do_angular_variation = .TRUE.
  ELSE
     n_angles = 0
     do_angular_variation = .FALSE.
  END IF

!  zero the angular input, if flagged

  IF ( do_angular_variation ) THEN
    DO angle = 1, n_angles
      DO kf = 1, 4
        fmat(kf,angle) = d_zero
      END DO
    END DO
  END IF

!  start integration
!  -----------------

  DO iblock = 1, nblocks

    rstart = rmin + ( iblock-1) * rblock
    rfinis = rstart + rblock

    CALL mie_gauleg ( max_Mie_sizes, n_sizes, rstart, rfinis, rquad, weights )

!  prepare particle sizes
    
    DO i = 1, n_sizes
      particle_sizes(i) = factor_0 * rquad(i)
    ENDDO

!  Call to coefficients
!    WARNING - easy o get segmentation fault before this call
!              Check dimensioning first, Use lots of memory    !!!!!!!

    CALL mie_coeffs                                              &
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation,                                   & ! Input
         n_angles, n_sizes, m_complex,                           & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext, q_sca, asym, splus, sminus )                       ! Output
 
!  size distribution and derivatives

    CALL sizedis &
    ( max_Mie_sizes, idis, nr_parameters, rquad, n_sizes, &
            nr,  message, faild )

    IF ( faild ) THEN
      failmie = faild
      write(char5,'(i5)')iblock
      trace  = 'Trace   : Sixth check in Mie_main. Subroutine sizedis failed for block number '//char5
      action = 'Action  : Consult with R. Spurr' 
      RETURN
    END IF

!  Integration over particle sizes within block
!  --------------------------------------------

    DO i = 1, n_sizes

!  Number density, geometric cross-section, 3rd and 4th powers

      quad   = nr(i) * weights(i)
      quadr2 = quad   * rquad(i) * rquad(i)
      quadr3 = quadr2 * rquad(i)
      quadr4 = quadr3 * rquad(i)
      ndens  = ndens + quad
      gxsec  = gxsec + quadr2
      reff   = reff  + quadr3
      veff   = veff  + quadr4

!  Basic coefficients

      Qext  = Qext  + quad * q_ext(i)
      Qsca  = Qsca  + quad * q_sca(i)
      Qasy  = Qasy  + quad * asym(i)

!  angular variation loop

      IF ( do_angular_variation ) THEN
        DO angle = 1, n_angles
          sp  = splus(angle,i)
          sm  = sminus(angle,i)
          csp = CONJG(sp)
          csm = CONJG(sm)
          f(1) =   REAL   ( sp * csp + sm * csm )
          f(2) = - REAL   ( sm * csp + sp * csm )
          f(3) =   REAL   ( sp * csp - sm * csm )
          f(4) =   REAL ( ( sm * csp - sp * csm ) * c_mi )
          DO kf = 1, 4
            FMAT(kf,angle) = FMAT(kf,angle) + quad*f(kf) 
          ENDDO
        END DO
      END IF

!  Finish integration loops

    END DO
  END DO

!  Final Assignations
!  ------------------

!  F matrix stuff

  IF ( do_angular_variation ) THEN
    DO angle = 1, n_angles
      DO kf = 1, 4
        FMAT(kf,angle) = d_half * FMAT(kf,angle) / Qsca
      END DO
    END DO
  END IF

!  geometric cross-section

  gxsec   = d_pi * gxsec

!  asymmetry parameter

  Qasy = d_two * Qasy / Qsca

!  basic coefficients

  Qsca  = Qsca * factor_1
  Qext  = Qext * factor_1
  Qsca  = Qsca/gxsec
  Qext  = Qext/gxsec

!  single scattering albedo

  ssalbedo = Qsca/Qext

!  geometrical quantities

  volume= (4.0_dp/3.0_dp) * d_pi * reff
  reff  = d_pi * reff / gxsec

!  Variance output

  help  = d_pi / gxsec / reff / reff 
  veff  = help * veff
  veff  = veff - d_one

!  Particle size parameter output

!  xeff  = factor_0 * reff

!  Final assignation

  MIE_BULK(1) = Qext
  MIE_BULK(2) = Qsca
  MIE_BULK(3) = Qasy
  MIE_BULK(4) = ssalbedo

  DIST(1) = ndens
  DIST(2) = gxsec
  DIST(3) = volume
  DIST(4) = reff
!  DIST(5) = xeff
  DIST(5) = veff

  RETURN
END SUBROUTINE Mie_main


SUBROUTINE mie_coeffs                                            & !
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation, n_angles, n_sizes, m_complex,     & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext,  q_sca,  asym,  splus,  sminus )                   ! Output

! name:
!       mie_coeffs

! purpose:
!       calculates the scattering parameters of a series of particles
!       using the mie scattering theory. FOR USE WITH POLYDISPERSE COD

! inputs:
!       particle_sizes:       array of particle size parameters
!       angle_cosines:        array of angle cosines
!       m_complex:            the complex refractive index of the particles
!       n_angles, n_sizes:    number of scattering angles, number of particle sizes
!       do_angular_variation: flag for S+/S- output


! outputs (1):
!       q_ext:       the extinction efficiency
!       q_sca:       the scattering efficiency
!       asym:        the asymmetry parameter
!       splus:       the first amplitude function
!       sminus:      the second amplitude function

! modification history
!       g. thomas IDL Mie code       (February  2004). Basic Monodisperse derivatives.
!       r. spurr  F90 Mie code       ( October  2004). Extension all Polydisperse derivatives.
!       r. spurr  Exception handling (September 2008). Exception handling removed.

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three, &
                               c_zero, c_one, c_i

!  implicit none statement

  IMPLICIT NONE

!  Dimensioning input

  INTEGER, INTENT (IN) :: max_Mie_angles, max_Mie_sizes
  INTEGER, INTENT (IN) :: max_Mie_points

!  input

  LOGICAL          , INTENT (IN) :: do_angular_variation
  INTEGER          , INTENT (IN) :: n_angles, n_sizes 
  COMPLEX (KIND=dp), INTENT (IN) :: m_complex

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes),  INTENT (IN) :: particle_sizes
  REAL    (KIND=dp), DIMENSION (max_Mie_angles), INTENT (IN) :: angle_cosines

!  output (1)

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes),                 INTENT (OUT) :: q_ext, q_sca, asym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles, max_Mie_sizes), INTENT (OUT) :: splus, sminus

!  local variables for Mie code

  INTEGER            :: size, angle, n, nm1, nstop(max_Mie_sizes), nmax, maxstop
  REAL    (KIND=dp)  :: xparticle, xpart_root3
  REAL    (KIND=dp)  :: xinv, xinvsq, two_d_xsq, xinv_dx
  REAL    (KIND=dp)  :: dn, dnp1, dnm1, dnsq, dnnp1, tnp1, tnm1, hnp1, hnm1
  REAL    (KIND=dp)  :: cos_x, sin_x, psi0, psi1, chi0, chi1, psi, chi
  REAL    (KIND=dp)  :: s, t, tau_n, factor, forward, bckward

  COMPLEX (KIND=dp)  :: inverse_m, y_argument, yinv, yinvsq, a1, zeta, zeta1
  COMPLEX (KIND=dp)  :: an, bn, an_star, bn_star, anm1, bnm1, bnm1_star
  COMPLEX (KIND=dp)  :: biga_divs_m, biga_mult_m, noverx, aterm, bterm
  COMPLEX (KIND=dp)  :: facplus, facminus, c_mi
  COMPLEX (KIND=dp)  :: an_denom, bn_denom

!  redundant variables
!  REAL    (KIND=dp), INTENT (IN) :: xparticle_limit           ! subroutine argument
!  REAL    (KIND=dp)  :: four_d_xsq
!  CHARACTER*4        :: char4
!  COMPLEX (KIND=dp)  :: common, an_denom_dm, s1, s2
!  INTEGER            :: nmax_end

!  local arrays

  REAL    (KIND=dp), DIMENSION (max_Mie_angles)                :: pi_n, pi_nm1
  COMPLEX (KIND=dp), DIMENSION (max_Mie_points)                :: biga
  REAL    (KIND=dp), DIMENSION (max_Mie_angles,max_Mie_points) :: polyplus, polyminus

!  Initial section
!  ---------------

  maxstop = 0

  c_mi = -c_i

  DO size = 1, n_sizes

!  particle size

    xparticle  = particle_sizes (size)

!  assign number of terms and maximum

    IF ( xparticle < 0.02) THEN
      nstop(size) = 2
    ELSE
      xpart_root3 = xparticle ** ( d_one / d_three )
      IF ( xparticle <= 8.0_dp ) THEN
        nstop(size) = xparticle + 4.0_dp * xpart_root3 + d_two
      ELSE IF ( xparticle < 4200.0_dp ) THEN
        nstop(size) = xparticle + 4.05_dp * xpart_root3 + d_two
      ELSE
      nstop(size) = xparticle + 4.0_dp * xpart_root3 + d_two
      END IF
    END IF
    maxstop = max(nstop(size),maxstop)

  END DO

!  phase function expansion polynomials
!    ---> initialise phase function Legendre polynomials
!    ---> Recurrence phase function Legendre polynomials

  IF ( do_angular_variation ) THEN

    DO angle = 1, n_angles
      pi_nm1(angle) = d_zero
      pi_n(angle)   = d_one
    END DO

    DO n = 1, maxstop
      nm1 = n - 1
      dn    = dble(n)
      dnp1  = dn   + d_one
      forward = dnp1 / dn
      DO angle = 1, n_angles
        s = angle_cosines(angle) * pi_n(angle)
        t = s - pi_nm1(angle)
        tau_n = dn*t - pi_nm1(angle)
        polyplus(angle,n)  = pi_n(angle) + tau_n
        polyminus(angle,n) = pi_n(angle) - tau_n
        pi_nm1(angle) = pi_n(angle)
        pi_n(angle)   = s + t*forward
      END DO
    END DO

  END IF

!  start loop over particle sizes
!  ------------------------------

  DO size = 1, n_sizes

!  initialize output

    asym(size)  = d_zero
    q_ext(size) = d_zero
    q_sca(size) = d_zero

!  some auxiliary quantities

    xparticle  = particle_sizes (size)
    xinv       = d_one / xparticle
    xinvsq     = xinv * xinv
    two_d_xsq  = d_two * xinvsq
    xinv_dx    = - d_two * xinv
 
    inverse_m  =     c_one / m_complex
    y_argument = xparticle * m_complex
    yinv       = d_one / y_argument 
    yinvsq     = yinv * yinv

!  Biga = ratio derivative, recurrence due to J. Dave

    nmax = nint(max(dble(nstop(size)),abs(y_argument)) + 15.0_dp)
    biga(nmax) = c_zero
    DO n = nmax-1, 1,-1
       a1      = dble(n+1) / y_argument
       biga(n) = a1 - c_one / (a1+biga(n+1))
    END DO

! initialize Riccati-Bessel functions

    tnp1 = d_one
    cos_x = COS(xparticle)
    sin_x = SIN(xparticle)
    psi0 = cos_x
    psi1 = sin_x
    chi1 =-cos_x
    chi0 = sin_x
    zeta1 = CMPLX(psi1,chi1,kind=dp)

!  initialise sp and sm

    IF ( do_angular_variation ) THEN
      DO angle = 1, n_angles
         splus(angle,size)  = c_zero
         sminus(angle,size) = c_zero
      END DO
    END IF

!  main loop

    DO n = 1, nstop(size)

!  various factors

      dn    = dble(n)
      dnp1  = dn   + d_one
      dnm1  = dn   - d_one
      tnp1  = tnp1 + d_two
      tnm1  = tnp1 - d_two

      dnsq  = dn * dn
      dnnp1 = dnsq + dn
      factor  = tnp1 / dnnp1
      bckward = dnm1 / dn

!  Ricatti - Bessel recurrence

      psi = tnm1 * psi1/xparticle - psi0
      chi = tnm1 * chi1/xparticle - chi0
      zeta = CMPLX(psi,chi,kind=dp)

!  a(n) and b(n)

      biga_divs_m = biga(n) * inverse_m
      biga_mult_m = biga(n) * m_complex
      noverx = CMPLX(dn/xparticle,d_zero,kind=dp)
      aterm = biga_divs_m + noverx
      bterm = biga_mult_m + noverx

      an_denom = (aterm * zeta - zeta1)
      bn_denom = (bterm * zeta - zeta1)
      an =   ( aterm*psi-psi1 ) / an_denom
      bn =   ( bterm*psi-psi1 ) / bn_denom
      an_star = CONJG(an)
      bn_star = CONJG(bn)

!  basic coefficients
!  ------------------

!  Q coefficients

      q_ext(size) = q_ext(size) + tnp1 * REAL ( an + bn )   
      q_sca(size) = q_sca(size) + tnp1 * REAL ( an*CONJG(an) + bn*CONJG(bn) )

!  asymmetry parameter

      IF ( n > 1 ) THEN
        hnp1 = bckward * dnp1
        hnm1 = tnm1  / (dnsq - dn)
        asym(size) = asym(size) &
           +  hnp1 * REAL ( anm1*an_star + bnm1*bn_star) &
           +  hnm1 * REAL ( anm1*bnm1_star) 
      END IF

!  Upgrades
!  --------

!  upgrade an/bn recurrences (only for asymmetry parameter)

      anm1      = an
      bnm1      = bn
      bnm1_star = bn_star

!  upgrade Ricatti-Bessel recurrences

      psi0 = psi1
      psi1 = psi
      chi0 = chi1
      chi1 = chi
      zeta1 = CMPLX(psi1,chi1,kind=dp)

!  S+/S- function stuff
!  --------------------

      IF ( do_angular_variation ) THEN
        facplus  = factor * ( an + bn )
        facminus = factor * ( an - bn )
        DO angle = 1, n_angles
          splus(angle,size)  = splus(angle,size)  + facplus  * polyplus(angle,n)
          sminus(angle,size) = sminus(angle,size) + facminus * polyminus(angle,n)
        END DO
      END IF

!  end sum loop

    END DO

!  End loop and finish
!  -------------------

!  end loop over particle sizes

  END DO

!  finish

  RETURN
END SUBROUTINE mie_coeffs


!  Contains the following modules

!     sizedist
!     sizedist_nod
!     gammafunction
!     gauleg
!     rminmax

SUBROUTINE sizedis                                      &
    ( max_Mie_distpoints, idis, par, radius, numradius, &
      nwithr, message, faild )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************
!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!* subroutine arguments

  INTEGER          , INTENT (IN)  :: max_Mie_distpoints
  REAL    (KIND=dp), INTENT (IN)  :: par(3)
  INTEGER          , INTENT (IN)  :: idis, numradius
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: faild

  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (IN)  :: radius
  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (OUT) :: nwithr

!* local variables

  INTEGER         :: i
  REAL  (KIND=dp) :: pi,r,logr,root2p
  REAL  (KIND=dp) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3
  REAL  (KIND=dp) :: b1,b2,logb1,logb2,rc
  REAL  (KIND=dp) :: logrg,logsi,logsi_inv,gamma,gamma1,rg
  REAL  (KIND=dp) :: rmin,rmax,fac1,fac2,aperg, alpha2
  REAL  (KIND=dp) :: n1, n2, C, logC, logC1, logC2

  REAL    (KIND=dp) :: gammln, dummy
  CHARACTER*70      :: message_gamma
  LOGICAL           :: fail
  character*1       :: cdis

!  check

      faild = .FALSE.

      if (idis == 0 ) RETURN
      IF ( IDIS > 8 ) THEN
        faild = .TRUE.
        message = 'illegal index in sizedis'
        RETURN
      END IF

!  setup

      pi     = dacos(-1.d0)
      root2p = dsqrt(pi+pi)

!  IDIS = 1 : TWO-PARAMETER GAMMA with alpha and b given

  IF ( idis == 1 ) THEN

      alpha  = par(1)
      b      = par(2)
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
         r  = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given

  ELSE IF ( idis == 2 ) THEN

      alpha  = d_one/par(2) - d_three
      b      = d_one/(par(1)*par(2))
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 3 : BIMODAL GAMMA with equal mode weights

  ELSE IF ( idis == 3 ) THEN

      alpha  = d_one/par(3) - d_three
      b1     = d_one/(par(1)*par(3))
      b2     = d_one/(par(2)*par(3))
      alpha1 = alpha + d_one
      CALL gammafunction ( alpha1, .false., gammln, dummy, fail, message_gamma )
      logb1 = LOG(b1)
      logb2 = LOG(b2)
      logC1 = alpha1*logb1 - gammln
      logC2 = alpha1*logb2 - gammln

      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC1 + alpha*logr
         arg2 = logC2 + alpha*logr
         n1   = EXP(arg1 - b1*r)
         n2   = EXP(arg2 - b2*r)
         nwithr(i) = d_half * ( n1 + n2 )
      END DO
  
!  4  LOG-NORMAL with rg and sigma given

  ELSE IF ( idis == 4 ) THEN

      logrg = dlog(par(1))
      logsi = dabs(dlog(par(2)))
      logsi_inv = d_one / logsi
      C      = logsi_inv / root2p
      DO i = 1, numradius
         r     = radius(i)
         logr  = LOG(r)
         arg   = ( logr - logrg ) / logsi
         argsq = arg * arg
         nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  5 LOG-NORMAL with reff and veff given                               *

  ELSE IF ( idis == 5 ) THEN

      alpha1 = d_one + par(2)
      alpha2 = dlog(alpha1)
      rg     = par(1)/(d_one+par(2))**2.5_dp
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      DO i = 1, numradius
         r     = radius(i)
         logr  = LOG(r)
         arg   = ( logr - logrg ) / logsi
         argsq = arg * arg
         nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  6 POWER LAW                               *

  ELSE IF ( idis == 6 ) THEN

      alpha = par(1)
      rmin  = par(2)
      rmax  = par(3)
      alpha1 = alpha - d_one
      fac1 = (rmax/rmin)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * rmax**alpha1 * fac2
      DO i = 1, numradius
         r     = radius(i)
         if ( (r < rmax) .and. (r > rmin) ) then
            nwithr(i)    = C*r**(-alpha)
         else
            nwithr(i)    = d_zero
         endif
      END DO

!  7 MODIFIED GAMMA with alpha, rc and gamma given

  ELSE IF ( idis == 7 ) THEN

      alpha = par(1)
      rc    = par(2)
      gamma = par(3)
      b     = alpha / (gamma*rc**gamma)
      logb  = LOG(b)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         r3   = b*r ** gamma
         nwithr(i) = EXP ( arg1 - r3 )
      END DO

!  8 MODIFIED GAMMA with alpha, b and gamma given

  ELSE IF ( idis == 8 ) THEN

      alpha = par(1)
      b     = par(2)
      gamma = par(3)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      logb = LOG(b)
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         r3   = r ** gamma
         nwithr(i) = EXP ( arg1 - b*r3 )
      END DO

  END IF

!  normal return

  RETURN

!  special return

240  CONTINUE
  faild = .TRUE.
  write(cdis,'(I1)')idis
  message = message_gamma(1:LEN(message_gamma))//', distribution : '//cdis
  RETURN

END SUBROUTINE sizedis



SUBROUTINE sizedis_nod &
    ( idis, par, numradius, radius, nwithr, message, failnod )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************
!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!* subroutine arguments

  REAL    (KIND=dp), INTENT (IN)  :: par(3)
  INTEGER          , INTENT (IN)  :: idis, numradius
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: failnod

  REAL    (KIND=dp), DIMENSION (numradius), INTENT (IN)  :: radius
  REAL    (KIND=dp), DIMENSION (numradius), INTENT (OUT) :: nwithr

!* local variables

  LOGICAL         :: deriv
  INTEGER         :: i
  REAL  (KIND=dp) :: pi,r,logr,root2p
  REAL  (KIND=dp) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3

  REAL  (KIND=dp) :: b1,b2,logb1,logb2,rc, gammln, dgammln
  REAL  (KIND=dp) :: logrg,logsi,logsi_inv,gamma,gamma1,rg
  REAL  (KIND=dp) :: rmin,rmax,fac1,fac2,aperg,alpha2
  REAL  (KIND=dp) :: n1, n2, C, logC, logC1, logC2

  CHARACTER*70      :: message_gamma
  LOGICAL           :: fail
  character*1       :: cdis

!  check

  failnod = .FALSE.
  if (idis == 0 ) RETURN
  IF ( IDIS > 8 ) THEN
     failnod = .TRUE.
     message = 'illegal index in sizedis'
     RETURN
  END IF

!  setup

  pi     = dacos(-1.d0)
  root2p = dsqrt(pi+pi)
  deriv  = .FALSE.

!  IDIS = 1 : TWO-PARAMETER GAMMA with alpha and b given

  IF ( idis == 1 ) THEN

      alpha  = par(1)
      b      = par(2)
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
        r  = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given

  ELSE IF ( idis == 2 ) THEN

      alpha  = d_one/par(2) - d_three
      b      = d_one/(par(1)*par(2))
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 3 : BIMODAL GAMMA with equal mode weights

  ELSE IF ( idis == 3 ) THEN

      alpha  = d_one/par(3) - d_three
      b1     = d_one/(par(1)*par(3))
      b2     = d_one/(par(2)*par(3))
      alpha1 = alpha + d_one
      CALL gammafunction ( alpha1, deriv, gammln, dgammln, fail, message_gamma )
      logb1 = LOG(b1)
      logb2 = LOG(b2)
      logC1 = alpha1*logb1 - gammln
      logC2 = alpha1*logb2 - gammln

      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC1 + alpha*logr
        arg2 = logC2 + alpha*logr
        n1   = EXP(arg1 - b1*r)
        n2   = EXP(arg2 - b2*r)
        nwithr(i) = d_half * ( n1 + n2 )
      END DO
  
!  4  LOG-NORMAL with rg and sigma given

  ELSE IF ( idis == 4 ) THEN

      logrg = dlog(par(1))
      logsi = dabs(dlog(par(2)))
      logsi_inv = d_one / logsi
      C      = logsi_inv / root2p
      DO i = 1, numradius
        r     = radius(i)
        logr  = LOG(r)
        arg   = ( logr - logrg ) / logsi
        argsq = arg * arg
        nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  5 LOG-NORMAL with reff and veff given                               *

  ELSE IF ( idis == 5 ) THEN

      alpha1 = d_one + par(2)
      alpha2 = dlog(alpha1)
      rg     = par(1)/(d_one+par(2))**2.5_dp
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      DO i = 1, numradius
        r     = radius(i)
        logr  = LOG(r)
        arg   = ( logr - logrg ) / logsi
        argsq = arg * arg
        nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  6 POWER LAW                               *

  ELSE IF ( idis == 6 ) THEN

      alpha = par(1)
      rmin  = par(2)
      rmax  = par(3)
      alpha1 = alpha - d_one
      fac1 = (rmax/rmin)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * rmax**alpha1 * fac2
      DO i = 1, numradius
        r     = radius(i)
        if ( (r < rmax) .and. (r > rmin) ) then
          nwithr(i)    = C*r**(-alpha)
        else
          nwithr(i)    = d_zero
        endif
      END DO

!  7 MODIFIED GAMMA with alpha, rc and gamma given

  ELSE IF ( idis == 7 ) THEN

      alpha = par(1)
      rc    = par(2)
      gamma = par(3)
      b     = alpha / (gamma*rc**gamma)
      logb  = LOG(b)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        r3   = b*r ** gamma
        nwithr(i) = EXP ( arg1 - r3 )
      END DO

!  8 MODIFIED GAMMA with alpha, b and gamma given

  ELSE IF ( idis == 8 ) THEN

      alpha = par(1)
      b     = par(2)
      gamma = par(3)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      logb = LOG(b)
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        r3   = r ** gamma
        nwithr(i) = EXP ( arg1 - b*r3 )
      END DO

  END IF

!  normal return

  RETURN

!  special return

240  CONTINUE
     failnod = .TRUE.
  write(cdis,'(I1)')idis
  message = message_gamma(1:LEN(message_gamma))//', distribution : '//cdis
  RETURN

END SUBROUTINE sizedis_nod


SUBROUTINE rminmax( idis, par, cutoff, rmin, rmax, message, fail )

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!  subroutine arguments

  REAL    (KIND=dp), INTENT (IN)  :: par(3), cutoff
  INTEGER          , INTENT (IN)  :: idis
  CHARACTER*(*)    , INTENT (OUT) :: message
  REAL    (KIND=dp), INTENT (OUT) :: rmin, rmax
  LOGICAL          , INTENT (OUT) :: fail

!  local variables

  REAL    (KIND=dp)            :: r(1), nwithr(1),ref,rref,sef,r0,r1,eps
  LOGICAL                      :: failnod

!************************************************************************
!*  Find the integration bounds rmin and rmax for the integration over  *
!*  a size distribution. These bounds are chosen such that the size     *
!*  distribution falls below the user specified cutoff. It is essential *
!*  that the size distribution is normalized such that the integral     *
!*  over all r is equal to one !                                        *
!*  This is programmed rather clumsy and will in the future be changed  *
!************************************************************************

  fail = .FALSE.
  eps = 1.0E-10_dp

  IF (idis == 0) THEN
     rmin= par(1)
     rmax= par(1)
     RETURN
  ELSE IF ( idis == 1)  THEN
     sef  = d_one/SQRT(par(2)+d_three)
     ref  = d_one/(sef*sef*par(2))
     rref = ref
  ELSE IF ( idis == 2)  THEN
     ref = par(1)
     sef = SQRT(par(2))
     rref= ref
  ELSE IF ( idis == 3)  THEN
     sef = SQRT(par(3))
     ref = MAX(par(1),par(2))+sef
     rref= MAX(par(1),par(2))
  ELSE IF ( idis == 4)  THEN
     sef = SQRT(EXP(LOG(par(2))**d_two)-d_one)
     ref = par(1)*(d_one+sef*sef)**0.4_dp
     rref= ref
  ELSE IF ( idis == 5)  THEN
     ref = par(1)
     sef = SQRT(ref)
     rref= ref
  ELSE IF ( idis == 6)  THEN
     rmin= par(2)
     rmax= par(3)
     RETURN
  ELSE IF ( idis == 7)  THEN
     ref = par(2)
     sef = d_two*ref
     rref= d_half*ref
  ELSE IF ( idis == 8)  THEN
     ref = (par(1)/(par(2)*par(3)))**par(3)
     sef = d_two*ref
     rref= d_half*ref
  END IF

!************************************************************************
!*  search for a value of r such that the size distribution
!*  is less than the cutoff. Start the search at ref+sef which          *
!*  guarantees that such a value will be found on the TAIL of the       *
!*  distribution.                                                       *
!************************************************************************

  r(1) = ref+sef
  r0   = ref
200 CONTINUE
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r0   = r(1)
     r(1) = d_two*r(1)
     goto 200
  END IF
  r1 = r(1)

!************************************************************************
!*  Now the size distribution assumes the cutoff value somewhere        *
!*  between r0 and r1  Use bisection to find the corresponding r        *
!************************************************************************

300 CONTINUE
  r(1) = d_half*(r0+r1)
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r0   = r(1)
  ELSE
     r1 = r(1)
  END IF
  IF ((r1-r0) > eps) GOTO 300
  rmax = d_half*(r0+r1)

!************************************************************************
!*  Search for a value of r on the low end of the size distribution     *
!*  such that the distribution falls below the cutoff. There is no      *
!*  guarantee that such a value exists, so use an extra test to see if  *
!*  the search comes very near to r = 0                                 *
!************************************************************************

  r1 = rref
  r0 = d_zero
400 CONTINUE
  r(1) = d_half*r1
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r1 = r(1)
     IF (r1 > eps) GOTO 400
  ELSE
     r0 = r(1)
  END IF

!************************************************************************
!*  Possibly the size distribution goes through cutoff between r0       *
!*  and r1 try to find the exact value of r where this happens by       *
!*  bisection.                                                          *
!*  In case there is no solution, the algorithm will terminate soon.    *
!************************************************************************

500 CONTINUE
  r(1) = d_half*(r0+r1)
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r1 = r(1)
  ELSE
     r0 = r(1)
  END IF
  IF ( (r1-r0) > eps) GOTO 500
  IF (r1 <= eps) THEN
     rmin = d_zero
  ELSE
     rmin = d_half*(r0+r1)
  END IF

! normal return

  RETURN

!  error return

899  CONTINUE
  fail = .TRUE.
  RETURN
END SUBROUTINE rminmax

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE gammafunction &
     ( xarg, do_derivative, gammln, dgammln, gammafail, message )

!************************************************************************
!*  Return the value of the natural logarithm of the gamma function.    *
!*  and its derivative (the Digamma function)                           *
!*  The argument xarg must be real and positive.                        *
!*  This function is documented in :                                    *
!*                                                                      *
!*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
!*  page 157 (ISBN 0-521-30811)                                         *
!*                                                                      *
!*  When the argument xarg is between zero and one, the relation (6.1.4)*
!*  on page 156 of the book by Press is used.                           *
!*                                         V.L. Dolman April 18 1989    *
!************************************************************************

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two

  IMPLICIT NONE

!* subroutine arguments

  REAL    (KIND=dp), INTENT (IN)  :: xarg
  LOGICAL          , INTENT (IN)  :: do_derivative
  REAL    (KIND=dp), INTENT (OUT) :: gammln, dgammln
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: gammafail

!* local parameters and data

  REAL    (KIND=dp), PARAMETER :: gammaf_eps = 1.d-10
  REAL    (KIND=dp), PARAMETER :: gammaf_fpf = 5.5D0
  REAL    (KIND=dp) :: cof(6),stp, c0
  DATA cof /    76.18009172947146D0,    &
               -86.50532032941677D0,    &
                24.01409824083091D0,    &
                -1.231739572450155D0,   &
                 0.1208650973866179D-2, &
                -0.5395239384953D-5      /
  DATA c0  / 1.000000000190015D0  /
  DATA stp / 2.5066282746310005D0 /

!* local variables

  INTEGER         :: j
  REAL  (KIND=dp) :: x,xx,xxx,tmp,x1,x2,logtmp,pi,ser,dser,gtmp,dgtmp,pix

!*  initialize output

  message   = ' '
  gammln    = d_zero
  dgammln   = d_zero
  gammafail = .FALSE.

!* check for bad input

  IF (xarg <= d_zero) THEN
     message = ' gammafunction: called with negative argument xarg'
     gammafail = .TRUE.
     RETURN
  END IF
  IF (ABS(xarg-d_one) < gammaf_eps) THEN
     message = ' gammafunction: argument too close to one '
     gammafail = .TRUE.
     RETURN
  END IF

!*  set up

  pi = 4.0_dp * ATAN(d_one)
  IF (xarg .ge. d_one) THEN
     xxx = xarg
  ELSE
     xxx = xarg + d_two
  END IF

!* Numerical Recipes stuff

  xx = xxx - d_one
  x1 = xx  + gammaf_fpf
  x2 = xx  + d_half

  logtmp = LOG(x1)
  tmp = x2*logtmp-x1

  ser = c0
  x = xx
  DO j =1, 6
     x = x + d_one
     ser = ser+cof(j)/x
  END DO
  gtmp = tmp + LOG(stp*ser)

!* derivative of gammln

  IF ( do_derivative ) THEN
     dser = d_zero
     x = xx
     DO j = 1, 6
        x = x+d_one
        dser = dser+cof(j)/x/x
     END DO
     dgtmp = logtmp - (5.0_dp/x1) - (dser/ser)
  END IF

!* assign output

  IF ( do_derivative ) THEN
     IF  (xarg > d_one) THEN
        gammln  = gtmp
        dgammln = dgtmp
     ELSE
        pix = pi*(d_one-xarg)
        gammln  = LOG(pix/SIN(pix))-gtmp
        dgammln = - dgtmp - (pi*COS(pix)/SIN(pix)) - (d_one/(d_one-xarg))
     END IF
  ELSE
     IF  (xarg > d_one) THEN
        gammln  = gtmp
     ELSE
        pix = pi*(d_one-xarg)
        gammln  = LOG(pix/SIN(pix))-gtmp
     END IF
  END IF

  RETURN
end SUBROUTINE gammafunction

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE mie_gauleg(maxn,n,x1,x2,x,w)

  USE Mie_precision

  IMPLICIT NONE

!* subroutine arguments

  INTEGER          , INTENT (IN)  :: maxn, n
  REAL    (KIND=dp), INTENT (IN)  :: x1,x2
  REAL    (KIND=dp), INTENT (OUT) :: x(maxn),w(maxn)

  INTEGER            :: i, m, j
  REAL    (KIND=dp)  ::  xm,xl,p1,p2,p3,pp,z,z1,eps

  eps=3.0e-14_dp
  m=(n+1)/2
  xm=0.5_dp*(x2+x1)
  xl=0.5_dp*(x2-x1)

  DO i=1,m
     z=COS(3.1415926540_dp*(DBLE(i)-0.250_dp)/(n+0.50_dp))
1  CONTINUE
     p1=1.0_dp
     p2=0.0_dp

     DO j=1,n
        p3=p2
        p2=p1
        p1=((2.0_dp*DBLE(j)-1.0_dp)*z*p2-(DBLE(j)-1.0_dp)*p3)/DBLE(j)
     END DO

     pp=n*(z*p1-p2)/(z*z-1.0_dp)
     z1=z
     z=z1-p1/pp
     if(dabs(z-z1) > eps) GO TO 1
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
     w(n+1-i)=w(i)
  END DO

  RETURN
END SUBROUTINE mie_gauleg


SUBROUTINE develop ( max_Mie_angles, ncoeffs, nangles, &
                     cosines, weights, FMAT, expcoeffs )

!  Based on the Meerhoff Mie code

!************************************************************************
!*  Calculate the expansion coefficients of the scattering matrix in    *
!*  generalized spherical functions by numerical integration over the   *
!*  scattering angle.                                   *
!************************************************************************

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three, d_four

!  implicit none statement

  IMPLICIT NONE

!  input

  INTEGER          , INTENT (IN) :: max_Mie_angles
  INTEGER          , INTENT (IN) :: ncoeffs, nangles
 
  REAL    (KIND=dp), INTENT (IN) :: cosines(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: weights(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: FMAT(4,max_Mie_angles)

!  output

  REAL    (KIND=dp), INTENT (OUT) :: expcoeffs(6,0:max_Mie_angles)

!  local variables

  REAL    (KIND=dp) :: P00(max_Mie_angles,2)
  REAL    (KIND=dp) :: P02(max_Mie_angles,2)
  REAL    (KIND=dp) :: P22(max_Mie_angles,2)
  REAL    (KIND=dp) :: P2m2(max_Mie_angles,2)
  REAL    (KIND=dp) :: fmatw(4,max_Mie_angles)

  INTEGER           :: i, j, l, lnew, lold, itmp
  INTEGER           :: index_11, index_12, index_22, index_33, index_34, index_44 
  REAL    (KIND=dp) :: dl, dl2, qroot6, fac1, fac2, fac3, fl,&
                       sql4, sql41, twol1, tmp1, tmp2, denom, &
                       alfap, alfam

!  Initialization

  qroot6 = -0.25_dp*SQRT(6.0_dp)
  index_11 = 1
  index_12 = 2
  index_22 = 3
  index_33 = 4
  index_34 = 5
  index_44 = 6

  DO j = 0, ncoeffs
    DO i = 1, 6
      expcoeffs(i,j) = d_zero
    END DO
  END DO

!  Multiply the scattering matrix F with the weights w for all angles  *
!  We do this here because otherwise it should be done for each l      *

  DO i = 1, 4
    DO j = 1, nangles
      fmatw(i,j) = weights(j)*FMAT(i,j)
    END DO
  END DO

!  Start loop over the coefficient index l                             *
!  first update generalized spherical functions, then calculate coefs. *
!  lold and lnew are pointer-like indices used in recurrence           *

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

    IF (l == 0) THEN

      dl   = d_zero
      DO  i=1, nangles
        P00(i,lold) = d_one
        P00(i,lnew) = d_zero
        P02(i,lold) = d_zero
        P22(i,lold) = d_zero
        P2m2(i,lold)= d_zero
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO

    ELSE

      dl   = DBLE(l)
      dl2  = dl * dl
      fac1 = (d_two*dl-d_one)/dl
      fac2 = (dl-d_one)/dl
      DO  i=1, nangles
        P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
      END DO

    ENDIF

    IF (l == 2) THEN

      DO  i=1, nangles
        P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
        P22(i,lold) = 0.25_dp*(d_one+cosines(i))*(d_one+cosines(i))
        P2m2(i,lold)= 0.25_dp*(d_one-cosines(i))*(d_one-cosines(i))
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO
      sql41 = d_zero

    ELSE IF (l > 2) THEN

      sql4  = sql41
      sql41 = dsqrt(dl2-d_four)
      twol1 = 2.D0*dl - d_one
      tmp1  = twol1/sql41
      tmp2  = sql4/sql41
      denom = (dl-d_one)*(dl2-d_four)
      fac1  = twol1*(dl-d_one)*dble(l)/denom
      fac2  = 4.D0*twol1/denom
      fac3  = dl*((dl-d_one)*(dl-d_one)-d_four)/denom
      DO i=1, nangles
        P02(i,lold) = tmp1*cosines(i)*P02(i,lnew)         - tmp2*P02(i,lold)
        P22(i,lold) = (fac1*cosines(i)-fac2)*P22(i,lnew)  - fac3*P22(i,lold)
        P2m2(i,lold)= (fac1*cosines(i)+fac2)*P2m2(i,lnew) - fac3*P2m2(i,lold)
      END DO

    END IF

    itmp = lnew
    lnew = lold
    lold = itmp
    alfap = d_zero
    alfam = d_zero

    fl = dl+d_half
    do i=1, nangles
      expcoeffs(index_11,l) = expcoeffs(index_11,l) + P00(i,lnew)*fmatw(1,i)
      alfap = alfap + P22(i,lnew)  * (fmatw(1,i)+fmatw(3,i))
      alfam = alfam + P2m2(i,lnew) * (fmatw(1,i)-fmatw(3,i))
      expcoeffs(index_44,l) = expcoeffs(index_44,l) + P00(i,lnew)*fmatw(3,i)
      expcoeffs(index_12,l) = expcoeffs(index_12,l) + P02(i,lnew)*fmatw(2,i)
      expcoeffs(index_34,l) = expcoeffs(index_34,l) + P02(i,lnew)*fmatw(4,i)
    END DO
    expcoeffs(index_11,l) =  fl*expcoeffs(index_11,l)
    expcoeffs(index_22,l) =  fl*d_half*(alfap+alfam)
    expcoeffs(index_33,l) =  fl*d_half*(alfap-alfam)
    expcoeffs(index_44,l) =  fl*expcoeffs(index_44,l)
    expcoeffs(index_12,l) =  fl*expcoeffs(index_12,l)
    expcoeffs(index_34,l) =  fl*expcoeffs(index_34,l)
  END DO

  RETURN
END SUBROUTINE develop



SUBROUTINE expand ( max_Mie_angles, ncoeffs, nangles, cosines, expcoeffs, FMAT )

!  Based on the Meerhoff Mie code

!  Use the expansion coefficients of the scattering matrix in 
!  generalized spherical functions to expand F matrix

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_one, d_two, d_four

!  implicit none statement

  IMPLICIT NONE

!  input

  INTEGER          , INTENT (IN) :: max_Mie_angles
  INTEGER          , INTENT (IN) :: ncoeffs, nangles
  REAL    (KIND=dp), INTENT (IN) :: cosines(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: expcoeffs(6,0:max_Mie_angles)

!  output

  REAL    (KIND=dp), INTENT (OUT) :: FMAT(4,max_Mie_angles)

!  local variables

  REAL    (KIND=dp) :: P00(max_Mie_angles,2)
  REAL    (KIND=dp) :: P02(max_Mie_angles,2)

  INTEGER           :: i, j, l, lnew, lold, itmp
  INTEGER           :: index_11, index_12, index_34, index_44 
  REAL    (KIND=dp) :: dl, qroot6, fac1, fac2, sql4, sql41, tmp1, tmp2

!  Initialization

  qroot6 = -0.25_dp*SQRT(6.0_dp)
  index_11 = 1
  index_12 = 2
  index_34 = 5
  index_44 = 6

!  Set scattering matrix F to zero

  DO j = 1, 4
    DO i = 1, nangles
      FMAT(j,i) = d_zero
    END DO
  END DO

!  Start loop over the coefficient index l
!  first update generalized spherical functions, then calculate coefs.
!  lold and lnew are pointer-like indices used in recurrence 

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

     IF ( l == 0) THEN

!  Adding paper Eqs. (76) and (77) with m=0

        DO i=1, nangles
           P00(i,lold) = d_one
           P00(i,lnew) = d_zero
           P02(i,lold) = d_zero
           P02(i,lnew) = d_zero
        END DO

     ELSE

        dl   = DBLE(l)
        fac1 = (d_two*dl-d_one)/dl
        fac2 = (dl-d_one)/dl

! Adding paper Eq. (81) with m=0

        DO i=1, nangles
           P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
        END DO

     END IF

     IF ( l == 2) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

        DO i=1, nangles
           P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
           P02(i,lnew) = d_zero
        END DO
        sql41 = d_zero

     ELSE IF ( l > 2) THEN

! Adding paper Eq. (82) with m=0

        sql4  = sql41
        sql41 = dsqrt(dl*dl-d_four)
        tmp1  = (d_two*dl-d_one)/sql41
        tmp2  = sql4/sql41

        DO i=1, nangles
           P02(i,lold) = tmp1*cosines(i)*P02(i,lnew) - tmp2*P02(i,lold)
        END DO

     END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

     itmp = lnew
     lnew = lold
     lold = itmp

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44

     DO i=1, nangles
        FMAT(1,i) = FMAT(1,i) + expcoeffs(index_11,l)*P00(i,lnew)
        FMAT(2,i) = FMAT(2,i) + expcoeffs(index_12,l)*P02(i,lnew)
        FMAT(3,i) = FMAT(3,i) + expcoeffs(index_44,l)*P00(i,lnew)
        FMAT(4,i) = FMAT(4,i) + expcoeffs(index_34,l)*P02(i,lnew)
     END DO

  END DO

  RETURN
END SUBROUTINE expand

