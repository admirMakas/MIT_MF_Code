!$Id: RTS_mie_sourcecode_plus.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
SUBROUTINE Mie_main_plus & 
     (   max_Mie_angles, max_Mie_sizes,                             & ! D
         max_Mie_points, max_Mie_distpoints,                        & ! D
         do_external_angles, do_coeffct_angles,                     & ! I 
         do_use_cutoff, do_m_derivatives,                           & ! I     
         idis, nr_parameters, do_nr_derivatives, startup,           & ! I
         nblocks, nweights, cutoff,                                 & ! I
         n_external_angles, external_angle_cosines,                 & ! I
         n_coeffct_angles, coeff_cosines, coeff_weights,            & ! I
         m_complex, xparticle_limit, wavelength, rmax, rmin,        & ! I
         mie_bulk, mie_bulk_d, dist, dist_d, fmat, fmat_d,          & ! O
         message, trace, action, failmie )                            ! O

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three

!  implicit none statement

  IMPLICIT NONE

!  Dimensioning input

  INTEGER, INTENT (IN) :: max_Mie_angles, max_Mie_sizes
  INTEGER, INTENT (IN) :: max_Mie_points, max_Mie_distpoints

!  input

  LOGICAL          , INTENT (IN) :: do_external_angles
  LOGICAL          , INTENT (IN) :: do_coeffct_angles
  LOGICAL          , INTENT (IN) :: do_m_derivatives
  LOGICAL          , INTENT (IN) :: do_use_cutoff

  INTEGER          , INTENT (IN) :: idis
  REAL    (KIND=dp), INTENT (IN) :: nr_parameters(3)
  LOGICAL          , INTENT (IN) :: do_nr_derivatives(3)

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

  REAL    (KIND=dp), INTENT (OUT)   :: MIE_BULK(4), MIE_BULK_D(4,5)
  REAL    (KIND=dp), INTENT (OUT)   :: DIST(5), DIST_D(5,3)
  REAL    (KIND=dp), INTENT (OUT)   :: FMAT(4,max_Mie_angles), FMAT_D(4,5,max_Mie_angles)

  LOGICAL          , INTENT (OUT)   :: failmie
  CHARACTER*(*)    , INTENT (OUT)   :: message, trace, action
  REAL    (KIND=dp), INTENT (INOUT) :: rmin, rmax

!  Local Mie output 

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes)                  :: q_ext, q_sca, asym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles, max_Mie_sizes)  :: splus,sminus
  REAL    (KIND=dp), DIMENSION (max_Mie_sizes,2)                :: dq_ext, dq_sca, dasym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles,max_Mie_sizes,2) :: dsplus, dsminus

!  local variables for Mie code

  CHARACTER*5 :: char5
  LOGICAL     :: do_Mie_linearization, do_angular_variation
  LOGICAL     :: failmm, faild
  INTEGER     :: i, kd, md, mdoff, angle, n_angles
  INTEGER     :: iblock, n_sizes, kf, nderivs, nkderivs

  REAL    (KIND=dp)  :: factor_0, factor_1, d_pi
  REAL    (KIND=dp)  :: rstart, rfinis, help, rblock
  REAL    (KIND=dp)  :: quad, quadr2, quadr3, quadr4, quad_d, quadr2_d, quadr3_d, quadr4_d
  REAL    (KIND=dp)  :: ndens, gxsec, reff, volume, veff
  REAL    (KIND=dp)  :: ndens_d(3), gxsec_d(3), reff_d(3), volume_d(3), veff_d(3)
  REAL    (KIND=dp)  :: Qext, Qsca, Qasy, ssalbedo
  REAL    (KIND=dp)  :: Qext_d(5), Qsca_d(5), Qasy_d(5)
  REAL    (KIND=dp)  :: f(4), df(4), ssalbedo_d(5)
  REAL    (KIND=dp)  :: angle_cosines(max_Mie_angles)

!  redundant variables
!  REAL    (KIND=dp)  :: xeff, xeff_d(3)
!  LOGICAL            :: fail
 
  COMPLEX (KIND=dp)  :: sp, sm, csp, csm, sp_dmr, sm_dmr, csp_dmr, csm_dmr
  COMPLEX (KIND=dp)  :: c_i, c_mi

  REAL    (KIND=dp)  :: xpart_root3, xparticle
  COMPLEX (KIND=dp)  :: y_argument
  INTEGER            :: limmax, limstop, limsize

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes)   :: particle_sizes
  REAL    (KIND=dp), DIMENSION (max_Mie_sizes)   :: rquad, weights, nr
  REAL    (KIND=dp), DIMENSION (max_Mie_sizes,3) :: nr_derivs

!  start up
!  --------

  c_i  = ( 0.0_dp, 1.0_dp )
  c_mi = - c_i

  do_Mie_linearization = do_m_derivatives
  DO kd = 1, 3
    do_Mie_linearization = do_Mie_linearization .or. do_nr_derivatives(kd)
  END DO
  mdoff = 0
  IF ( do_m_derivatives ) mdoff = 2  
  nkderivs = 0
  DO kd = 1, 3
    IF ( do_nr_derivatives(kd) ) nkderivs = nkderivs + 1
  END DO
  nderivs = mdoff + nkderivs
  n_sizes = nweights

  d_pi = 4.0_dp * ATAN(d_one)
  factor_0 = d_two * d_pi / wavelength
  factor_1 = wavelength / factor_0

!  Zeroing
!  -------

  message = ' '
  trace   = ' '
  action  = ' '
  failmie = .FALSE.
  failmm  = .FALSE.

  Qext = d_zero
  Qsca = d_zero
  Qasy = d_zero

  ndens = d_zero
  gxsec = d_zero
  reff  = d_zero
  veff  = d_zero

  IF ( do_Mie_linearization ) THEN
    DO md = 1, nderivs
      Qext_d(md)  = d_zero
      Qsca_d(md)  = d_zero
      Qasy_d(md)  = d_zero
    END DO
  END IF

  IF ( do_Mie_linearization ) THEN
    DO kd = 1, nkderivs
      ndens_d(kd) = d_zero
      gxsec_d(kd) = d_zero
      reff_d(kd)  = d_zero
      veff_d(kd)  = d_zero
    END DO
  END IF

!  limiting radii

  IF ( do_use_cutoff ) THEN
     CALL rminmax ( idis, nr_parameters, cutoff, rmin, rmax, message, failmm )
     IF ( failmm ) THEN
       failmie = .TRUE.
       trace   = 'Trace   : First Check in Mie Main +. Failed to find radii extrema'
       action  = 'Action  : Consult with R. Spurr' 
       RETURN
     END IF
  END  IF

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
       trace   = 'Trace   : First Check in Mie Main +. User Rmin/Rmax wrong'
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

!  debug
!  write(*,*)xparticle, rmax, rmin, wavelength
!  write(*,*)limmax, limstop, max_Mie_points

!  Dimensioning and exception handling checks
!  ------------------------------------------

!  return if size limit exceeded

  IF ( xparticle > xparticle_limit ) THEN
     failmie = .TRUE.
     limsize = int(xparticle) + 1
     write(char5,'(i5)')limsize
     message = 'Message : error size parameter overflow'
     trace   = 'Trace   : Second check in Mie_main_plus'
     action  = 'Action  : In configuration file, increase cutoff or '// &
                          'increase xparticle_limit to at least '//char5
     RETURN
  END IF

!  return if maximum number of terms too great 

  IF ( limstop > max_Mie_points )  THEN
     failmie = .TRUE.
     write(char5,'(i5)')limstop
     message = 'Message : Insufficient dimensioning for maximum number of terms'
     trace   = 'Trace   : Third check in Mie_main_plus'
     action  = 'Action  : Increase max_Mie_points in calling program to at least '//char5
     RETURN
  END IF

!  And again, Dave recurrence

  IF ( limmax > max_Mie_points )  THEN
     failmie = .TRUE.
     write(char5,'(i5)')limmax
     message = 'Message : Insufficient dimensioning for maximum number of terms (Dave recurrence)'
     trace   = 'Trace   : Fourth check in Mie_main_plus'
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
       trace   = 'Trace   : Fifth check in Mie Main_plus'
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
    IF ( do_Mie_linearization ) THEN
      DO angle = 1, n_angles
        DO kf = 1, 4
          DO md = 1, nderivs
            fmat_d(kf,md,angle) = d_zero
          END DO
        END DO
      END DO
    END IF
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

!  Get the basic Mie output

    CALL mie_coeffs_d &
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation, do_m_derivatives,                 & ! Input
         n_angles, n_sizes, m_complex,                           & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext,  q_sca,  asym,  splus,   sminus,                 & ! Output
         dq_ext, dq_sca, dasym, dsplus,  dsminus )                 ! Output

!  Non-linearized part
!      CALL mie_coeffs &
!       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
!         do_angular_variation,                                   & ! Input
!         n_angles, n_sizes, m_complex,                           & ! Input
!         particle_sizes, angle_cosines,                          & ! Input
!         q_ext, q_sca, asym, splus, sminus )                       ! Output
 
!  size distribution and derivatives

    CALL sizedis_plus &
    ( max_Mie_sizes, idis, nr_parameters, do_nr_derivatives, rquad, n_sizes, &
            nr, nr_derivs, message, faild )

    IF ( faild ) THEN
      failmie = faild
      write(char5,'(i5)')iblock
      trace  = 'Trace   : Sixth check in Mie_main_plus. Subroutine sizedis+ failed, block number '//char5
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

!  Basic coefficient derivatives w.r.t refractive index parameters

      IF ( do_m_derivatives ) THEN
        DO md = 1, 2
          Qext_d(md)  = Qext_d(md) + quad * dq_ext(i,md)
          Qsca_d(md)  = Qsca_d(md) + quad * dq_sca(i,md)
          Qasy_d(md)  = Qasy_d(md) + quad * dasym(i,md)
        END DO
      END IF

!  Basic coefficient derivatives w.r.t distribution parameters
!    distribution derivatives there as a check

      DO kd = 1, nkderivs
        IF ( do_nr_derivatives(kd) ) THEN
          md = kd + mdoff
          quad_d   = nr_derivs(i,kd) * weights(i)
          quadr2_d = quad_d   * rquad(i) * rquad(i)
          quadr3_d = quadr2_d * rquad(i)
          quadr4_d = quadr3_d * rquad(i)
          ndens_d(kd) = ndens_d(kd) + quad_d
          gxsec_d(kd) = gxsec_d(kd) + quadr2_d
          reff_d(kd)  = reff_d(kd)  + quadr3_d
          veff_d(kd)  = veff_d(kd)  + quadr4_d
          Qext_d(md)  = Qext_d(md)  + quad_d * q_ext(i)
          Qsca_d(md)  = Qsca_d(md)  + quad_d * q_sca(i)
          Qasy_d(md)  = Qasy_d(md)  + quad_d * asym(i)
        END IF
      END DO

!  angular variation loop

      IF ( do_angular_variation ) THEN

        DO angle = 1, n_angles

!  basic F-matrix

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

!  F-matrix derivatiaves w.r.t. refractive index parameters

          IF ( do_m_derivatives ) THEN
            DO md = 1, 2
              sp_dmr = dsplus(angle,i,md)
              sm_dmr = dsminus(angle,i,md)
              csp_dmr = CONJG(sp_dmr)
              csm_dmr = CONJG(sm_dmr)
              df(1) =   REAL ( sp * csp_dmr + sm * csm_dmr   &
                                 +  sp_dmr * csp + sm_dmr * csm )
              df(2) = - REAL ( sm * csp_dmr + sp * csm_dmr   &
                                +  sm_dmr * csp + sp_dmr * csm )
              df(3) =   REAL ( sp * csp_dmr - sm * csm_dmr   &
                                +  sp_dmr * csp - sm_dmr * csm )
              df(4) =   REAL ( ( sm * csp_dmr - sp * csm_dmr &
                                  +  sm_dmr * csp - sp_dmr * csm ) * c_mi )
              DO kf = 1, 4
                FMAT_D(kf,md,angle) = FMAT_D(kf,md,angle) + quad * df(kf) 
              ENDDO
            ENDDO
          END IF

!  F-matrix derivatives w.r.t distribution parameters

          DO kd = 1, nkderivs
            IF ( do_nr_derivatives(kd) ) THEN
              quad_d   = nr_derivs(i,kd) * weights(i)
              md = mdoff + kd
              DO kf = 1, 4
                FMAT_D(kf,md,angle) = FMAT_D(kf,md,angle) + quad_d * f(kf) 
              ENDDO
            END IF
          END DO

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
    IF ( do_Mie_linearization ) THEN
      DO angle = 1, n_angles
        DO kf = 1, 4
          DO md = 1, nderivs
            help = d_half * FMAT_D(kf,md,angle) - Qsca_d(md) * FMAT(kf,angle)
            FMAT_D(kf,md,angle) = help / Qsca
          END DO
        END DO
      END DO
    END IF
  END IF

!  geometric cross-section and derivatives

  gxsec   = d_pi * gxsec
  IF ( do_Mie_linearization ) THEN
    DO kd = 1, nkderivs
      gxsec_d(kd) = d_pi * gxsec_d(kd)
    END DO
  END IF

!  asymmetry parameter = derivatives

  Qasy = d_two * Qasy / Qsca
  IF ( do_Mie_linearization ) THEN
    DO md = 1, nderivs
      Qasy_d(md) = ( d_two * Qasy_d(md) - Qasy * Qsca_d(md) ) / Qsca
    END DO
  END IF

!  basic coefficients

  Qsca  = Qsca * factor_1
  Qext  = Qext * factor_1
  IF ( do_Mie_linearization ) THEN
    DO md = 1, nderivs
      Qext_d(md) = factor_1 * Qext_d(md)
      Qsca_d(md) = factor_1 * Qsca_d(md)
    END DO
  END IF

  Qsca  = Qsca/gxsec
  Qext  = Qext/gxsec
  IF ( do_Mie_linearization ) THEN
    IF ( do_m_derivatives ) THEN
      DO md = 1, 2
        Qext_d(md) = Qext_d(md) / gxsec
        Qsca_d(md) = Qsca_d(md) / gxsec
      END DO
    END IF
    DO kd = 1, nkderivs
      md = mdoff + kd
      IF ( do_nr_derivatives(kd) ) THEN
        Qext_d(md) = ( Qext_d(md) -  Qext * gxsec_d(kd) ) / gxsec
        Qsca_d(md) = ( Qsca_d(md) -  Qsca * gxsec_d(kd) ) / gxsec
      END IF
    END DO         
  END IF

!  single scattering albedo

  ssalbedo = Qsca/Qext
  IF ( do_Mie_linearization ) THEN
    DO md = 1, nderivs
      ssalbedo_d(md) = ( Qsca_d(md)-ssalbedo*Qext_d(md) ) / Qext
    END DO
  END IF

!  geometrical quantities

  volume= (4.0_dp/3.0_dp) * d_pi * reff
  IF ( do_Mie_linearization ) THEN
    DO kd = 1, nkderivs
      IF ( do_nr_derivatives(kd) ) THEN
        volume_d(kd) = volume * reff_d(kd) / reff
      END IF
    END DO
  END IF

  reff  = d_pi * reff / gxsec
  IF ( do_Mie_linearization ) THEN
    DO kd = 1, nkderivs
      IF ( do_nr_derivatives(kd) ) THEN
        reff_d(kd) = ( d_pi * reff_d(kd) - reff * gxsec_d(kd) ) / gxsec
      END IF
    END DO
  END IF

!  Variance output

  help  = d_pi / gxsec / reff / reff 
  veff  = help * veff
  IF ( do_Mie_linearization ) THEN
    DO kd = 1, nkderivs
      IF ( do_nr_derivatives(kd) ) THEN
        veff_d(kd) = help * veff_d(kd) - veff * ((gxsec_d(kd)/gxsec)+(d_two*reff_d(kd)/reff))
      END IF
    END DO
  END IF
  veff  = veff - d_one

!  Particle size parameter output

!  xeff  = factor_0 * reff
!  IF ( do_Mie_linearization ) THEN
!    DO kd = 1, nkderivs
!      IF ( do_nr_derivatives(kd) ) THEN
!        xeff_d(kd) = factor_0 * reff_d(kd) 
!      END IF
!    END DO
!  END IF

!  Final assignation

  MIE_BULK(1) = Qext
  MIE_BULK(2) = Qsca
  MIE_BULK(3) = Qasy
  MIE_BULK(4) = ssalbedo
  IF ( do_Mie_linearization ) THEN
    DO md = 1, nderivs
      MIE_BULK_D(1,md) = Qext_d(md)
      MIE_BULK_D(2,md) = Qsca_d(md)
      MIE_BULK_D(3,md) = Qasy_d(md)
      MIE_BULK_D(4,md) = ssalbedo_d(md)
    END DO
  END IF

  DIST(1) = ndens
  DIST(2) = gxsec
  DIST(3) = volume
  DIST(4) = reff
!  DIST(5) = xeff
  DIST(5) = veff

  IF ( do_Mie_linearization ) THEN
    DO kd = 1, nkderivs
      IF ( do_nr_derivatives(kd) ) THEN
        DIST_D(1,kd) = ndens_d(kd)   
        DIST_D(2,kd) = gxsec_d(kd)   
        DIST_D(3,kd) = volume_d(kd)   
        DIST_D(4,kd) = reff_d(kd)   
!        DIST_D(5,kd) = xeff_d(kd)
        DIST_D(5,kd) = veff_d(kd)
      END IF
    END DO
  END IF

  RETURN
END SUBROUTINE Mie_main_plus


SUBROUTINE mie_coeffs_d                                          &
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation, do_m_derivatives,                 & ! Input
         n_angles, n_sizes, m_complex,                           & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext,  q_sca,  asym,  splus,  sminus,                  & ! Output
         dq_ext, dq_sca, dasym, dsplus, dsminus )                  ! Output

! name:
!       mie_coeffs_d

! purpose:
!       calculates the scattering parameters of a series of particles
!       using the mie scattering theory. FOR USE WITH POLYDISPERSE CODE

! inputs:
!       particle_sizes:       array of particle size parameters
!       angle_cosines:        array of angle cosines
!       m_complex:            the complex refractive index of the particles
!       n_angles, n_sizes:    number of scattering angles, number of particle sizes
!       do_angular_variation: flag for S+/S- output
!       do_m_derivatives:     differentiation flag  w.r.t refractive index

! outputs (1):
!       q_ext:       the extinction efficiency
!       q_sca:       the scattering efficiency
!       asym:        the asymmetry parameter
!       splus:       the first amplitude function
!       sminus:      the second amplitude function

! outputs (2):
!       dq_ext:       derivatives of the extinction efficiency
!       dq_sca:       derivatives of the scattering efficiency
!       dasym:        derivatives of the asymmetry parameter
!       dsplus:       derivatives of the first amplitude function
!       dsminus:      derivatives of the second amplitude function

! modification history
!       g. thomas IDL Mie code       (February  2004). Basic Monodisperse derivatives.
!       r. spurr  F90 Mie code       ( October  2004). Extension all Polydisperse derivatives.
!       r. spurr  Exception handling (September 2008). Exception handling removed.

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three

!  implicit none statement

  IMPLICIT NONE

!  Dimensioning input

  INTEGER, INTENT (IN) :: max_Mie_angles, max_Mie_sizes
  INTEGER, INTENT (IN) :: max_Mie_points

!  input

  LOGICAL          , INTENT (IN) :: do_angular_variation, do_m_derivatives
  INTEGER          , INTENT (IN) :: n_angles, n_sizes 
  COMPLEX (KIND=dp), INTENT (IN) :: m_complex

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes),  INTENT (IN) :: particle_sizes
  REAL    (KIND=dp), DIMENSION (max_Mie_angles), INTENT (IN) :: angle_cosines

!  output (1)

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes),                 INTENT (OUT) :: q_ext, q_sca, asym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles, max_Mie_sizes), INTENT (OUT) :: splus, sminus

!  output (2)

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes,2),                INTENT (OUT) :: dq_ext, dq_sca, dasym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles,max_Mie_sizes,2), INTENT (OUT) :: dsplus, dsminus

!  CHARACTER*(*) , INTENT (OUT) :: message, action
!  LOGICAL       , INTENT (OUT) :: fail

!  local variables for Mie code

  INTEGER            :: size, angle, n, nm1, nstop(max_Mie_sizes), nmax, maxstop
  REAL    (KIND=dp)  :: xparticle, xpart_root3
  REAL    (KIND=dp)  :: xinv, xinvsq, two_d_xsq, xinv_dx
  REAL    (KIND=dp)  :: dn, dnp1, dnm1, dnsq, dnnp1, tnp1, tnm1, hnp1, hnm1
  REAL    (KIND=dp)  :: cos_x, sin_x, psi0, psi1, chi0, chi1, psi, chi
  REAL    (KIND=dp)  :: s, t, tau_n, factor, forward, bckward

  COMPLEX (KIND=dp)  :: inverse_m, y_argument, a1, zeta, zeta1
  COMPLEX (KIND=dp)  :: an, bn, an_star, bn_star, anm1, bnm1, bnm1_star, yinv, yinvsq
  COMPLEX (KIND=dp)  :: biga_divs_m, biga_mult_m, noverx, aterm, bterm
  COMPLEX (KIND=dp)  :: facplus, facminus, dfacplus, dfacminus, c_zero, c_one, c_i, c_mi
  COMPLEX (KIND=dp)  :: an_denom, bn_denom, an_denom_dm, common, a_num_dm, b_num_dm
  COMPLEX (KIND=dp)  :: an_dmr, an_dmi, bn_dmr, bn_dmi
  COMPLEX (KIND=dp)  :: an_star_dmr, an_star_dmi, bn_star_dmr, bn_star_dmi
  COMPLEX (KIND=dp)  :: anm1_dmr, anm1_dmi, bnm1_dmr, bnm1_dmi, bnm1_star_dmr, bnm1_star_dmi

!  redundant variables
!  REAL    (KIND=dp), INTENT (IN) :: xparticle_limit           ! subroutine argument
!  REAL    (KIND=dp)  :: four_d_xsq
!  CHARACTER*4        :: char4
!  COMPLEX (KIND=dp)  :: an_dm, bn_dm, s1, s2, ds1, ds2, sbck
!  INTEGER            :: nmax_end

!  local arrays

  REAL    (KIND=dp), DIMENSION (max_Mie_angles)                :: pi_n, pi_nm1
  COMPLEX (KIND=dp), DIMENSION (max_Mie_points)                :: biga, biga_sq, biga_sq1
  REAL    (KIND=dp), DIMENSION (max_Mie_angles,max_Mie_points) :: polyplus, polyminus

!  Initial section
!  ---------------

  maxstop = 0

  c_zero = ( 0.0_dp, 0.0_dp )
  c_one  = ( 1.0_dp, 0.0_dp )
  c_i    = ( 0.0_dp, 1.0_dp )
  c_mi   = - c_i

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

    IF ( do_m_derivatives ) THEN
      dq_ext(size,1) = d_zero
      dq_ext(size,2) = d_zero
      dq_sca(size,1) = d_zero
      dq_sca(size,2) = d_zero
      dasym(size,1)  = d_zero
      dasym(size,2)  = d_zero
    END IF

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
    IF ( do_m_derivatives ) THEN
       DO n = 1, nmax
          biga_sq(n) = biga(n) * biga(n)
          biga_sq1(n) = c_one + biga_sq(n)
       END DO
    END IF

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
      IF ( do_m_derivatives ) THEN
        DO angle = 1, n_angles
           dsplus(angle,size,1)  = c_zero
           dsminus(angle,size,1) = c_zero
           dsplus(angle,size,2)  = c_zero
           dsminus(angle,size,2) = c_zero
        END DO
      END IF
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

!  derivatives of an and bn w.r.t refractive index

      IF ( do_m_derivatives ) THEN
        an_denom_dm = an_denom * m_complex  
        common = dnnp1 * yinv - y_argument*biga_sq1(n) 
        a_num_dm =  common - biga(n)
        b_num_dm =  common + biga(n)
        an_dmi = a_num_dm / an_denom_dm / an_denom_dm
        bn_dmi = b_num_dm / bn_denom    / bn_denom
        an_dmr = c_mi * an_dmi
        bn_dmr = c_mi * bn_dmi
        an_star_dmr = CONJG(an_dmr)
        bn_star_dmr = CONJG(bn_dmr)
        an_star_dmi = CONJG(an_dmi)
        bn_star_dmi = CONJG(bn_dmi)
      END IF

!  basic coefficients
!  ------------------

!  Q coefficients

      q_ext(size) = q_ext(size) + tnp1 * REAL ( an + bn )   
      q_sca(size) = q_sca(size) + tnp1 * REAL ( an*CONJG(an) + bn*CONJG(bn) )

!  derivatives of Q coefficients w.r.t. complex refractive index

      IF ( do_m_derivatives ) THEN
        dq_ext(size,1) = dq_ext(size,1) + tnp1 * REAL ( an_dmr + bn_dmr )   
        dq_ext(size,2) = dq_ext(size,2) + tnp1 * REAL ( an_dmi + bn_dmi )   
        dq_sca(size,1) = dq_sca(size,1) + d_two * tnp1 * REAL ( an*an_star_dmr + bn_dmr*bn_star )
        dq_sca(size,2) = dq_sca(size,2) + d_two * tnp1 * REAL ( an*an_star_dmi + bn_dmi*bn_star )
      END IF

!  asymmetry parameter

      IF ( n > 1 ) THEN
        hnp1 = bckward * dnp1
        hnm1 = tnm1  / (dnsq - dn)
        asym(size) = asym(size) &
           +  hnp1 * REAL ( anm1*an_star + bnm1*bn_star) &
           +  hnm1 * REAL ( anm1*bnm1_star) 
        IF ( do_m_derivatives ) THEN
          dasym(size,1) = dasym(size,1) &
           +  hnp1 * REAL ( anm1 * an_star_dmr + anm1_dmr * an_star + &
                            bnm1 * bn_star_dmr + bnm1_dmr * bn_star ) &
           +  hnm1 * REAL ( anm1 * bnm1_star_dmr + anm1_dmr * bnm1_star)
          dasym(size,2) = dasym(size,2) &
           +  hnp1 * REAL ( anm1 * an_star_dmi + anm1_dmi * an_star + &
                            bnm1 * bn_star_dmi + bnm1_dmi * bn_star ) &
           +  hnm1 * REAL ( anm1 * bnm1_star_dmi + anm1_dmi * bnm1_star)
        END IF
      END IF

!  Upgrades
!  --------

!  upgrade an/bn recurrences (only for asymmetry parameter)

      anm1      = an
      bnm1      = bn
      bnm1_star = bn_star
      IF ( do_m_derivatives ) THEN
        anm1_dmr      = an_dmr
        bnm1_dmr      = bn_dmr
        bnm1_star_dmr = bn_star_dmr
        anm1_dmi      = an_dmi
        bnm1_dmi      = bn_dmi
        bnm1_star_dmi = bn_star_dmi
      END IF

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
        IF ( do_m_derivatives ) THEN
          dfacplus  = factor * ( an_dmr + bn_dmr )
          dfacminus = factor * ( an_dmr - bn_dmr )
          DO angle = 1, n_angles
            dsplus(angle,size,1)  = dsplus(angle,size,1)  + dfacplus  * polyplus(angle,n)
            dsminus(angle,size,1) = dsminus(angle,size,1) + dfacminus * polyminus(angle,n)
          END DO
          dfacplus  = factor * ( an_dmi + bn_dmi )
          dfacminus = factor * ( an_dmi - bn_dmi )
          DO angle = 1, n_angles
            dsplus(angle,size,2)  = dsplus(angle,size,2)  + dfacplus  * polyplus(angle,n)
            dsminus(angle,size,2) = dsminus(angle,size,2) + dfacminus * polyminus(angle,n)
          END DO
        END IF
      END IF

!  end sum loop

    END DO

!  End loop and finish
!  -------------------

!  end loop over particle sizes

  END DO

!  finish

  RETURN
END SUBROUTINE mie_coeffs_d

!  Contains the following modules

!     sizedist_plus
!     gammafunction
!     gauleg
!     rminmax


SUBROUTINE sizedis_plus                                        &
    ( max_Mie_distpoints, idis, par, deriv, radius, numradius, &
      nwithr, nwithr_d, message, faild )

!  Contains the following modules

!     sizedist_nod
!     sizedist
!     gammafunction
!     gauleg
!     rminmax

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
  LOGICAL          , INTENT (IN)  :: deriv(3)
  INTEGER          , INTENT (IN)  :: idis, numradius
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: faild

  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (IN)  :: radius
  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (OUT) :: nwithr
  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints,3), INTENT (OUT) :: nwithr_d

!* local variables

  INTEGER         :: i
  REAL  (KIND=dp) :: pi,r,logr,root2p
  REAL  (KIND=dp) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3
  REAL  (KIND=dp) :: b1,b2,b11,b13,b22,b23,logb1,logb2,rc
  REAL  (KIND=dp) :: logrg,logsi,logsi_inv,fac_d1,gamma,gamma1,rg
  REAL  (KIND=dp) :: rmin,rmax,fac1,fac2,aperg
  REAL  (KIND=dp) :: alpha2, fac_d2a
  REAL  (KIND=dp) :: n1, n2, n1_d1, n1_d3, n2_d2, n2_d3

!  redundant variables
!  REAL  (KIND=dp) :: sigfac, logC1_d2

  REAL  (KIND=dp) :: C, logC, logC_d1, logC_d2, logC_d3
  REAL  (KIND=dp) :: logC1, logC2, logC1_d1, logC1_d3, logC2_d2, logC2_d3

  REAL    (KIND=dp) :: gammln, dgammln
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
      CALL gammafunction ( alpha1, deriv(1), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240
      logC  = alpha1*logb - gammln

      IF ( deriv(1) .and. deriv(2)  ) then
        logC_d1 = logb - dgammln
        logC_d2 = alpha1/b
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
          nwithr_d(i,2) = ( logC_d2 - r )    * nwithr(i)
          nwithr_d(i,1) = ( logC_d1 + logr ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r  = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
        END DO
      END IF

!  IDIS = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given

  ELSE IF ( idis == 2 ) THEN

      alpha  = d_one/par(2) - d_three
      b      = d_one/(par(1)*par(2))
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv(2), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = alpha1*logb - gammln

      IF ( deriv(1) .and. deriv(2) ) then
        b1      = b / par(1)
        b2      = b / par(2) 
        logC_d1 = - alpha1 / par(1)
        logC_d2 = ( dgammln - logb - alpha1*par(2) ) / par(2) / par(2)
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
          nwithr_d(i,1) = ( logC_d1 + b1*r  )  * nwithr(i)
          nwithr_d(i,2) = ( logC_d2 - logr/par(2)/par(2) + b2*r )  * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
        END DO
      END IF

!  IDIS = 3 : BIMODAL GAMMA with equal mode weights

  ELSE IF ( idis == 3 ) THEN

      alpha  = d_one/par(3) - d_three
      b1     = d_one/(par(1)*par(3))
      b2     = d_one/(par(2)*par(3))
      alpha1 = alpha + d_one
      CALL gammafunction ( alpha1, deriv(3), gammln, dgammln, fail, message_gamma )
      logb1 = LOG(b1)
      logb2 = LOG(b2)
      logC1 = alpha1*logb1 - gammln
      logC2 = alpha1*logb2 - gammln

      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        b11      = b1 / par(1)
        b13      = b1 / par(3)
        b22      = b2 / par(2)
        b23      = b2 / par(3)
        logC1_d1 = - alpha1 / par(1)
        logC1_d3 = ( dgammln - logb1 - alpha1*par(3) ) / par(3) / par(3)
        logC2_d2 = - alpha1 / par(2)
        logC2_d3 = ( dgammln - logb2 - alpha1*par(3) ) / par(3) / par(3)
        DO i = 1, numradius
          r  = radius(i)
          logr = LOG(r)
          arg1 = logC1 + alpha*logr
          arg2 = logC2 + alpha*logr
          n1   = EXP(arg1 - b1*r)
          n2   = EXP(arg2 - b2*r)
          nwithr(i) = d_half * ( n1 + n2 )
          n1_d1 = ( logC1_d1 + b11*r  )                 * n1
          n1_d3 = ( logC1_d3 - logr/par(3)/par(3) + b13*r ) * n1
          n2_d2 = ( logC2_d2 + b22*r  )                 * n2
          n2_d3 = ( logC2_d3 - logr/par(3)/par(3) + b23*r ) * n2
          nwithr_d(i,1) = d_half * n1_d1
          nwithr_d(i,2) = d_half * n2_d2
          nwithr_d(i,3) = d_half * ( n1_d3 + n2_d3 ) 
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC1 + alpha*logr
          arg2 = logC2 + alpha*logr
          n1   = EXP(arg1 - b1*r)
          n2   = EXP(arg2 - b2*r)
          nwithr(i) = d_half * ( n1 + n2 )
        END DO
      END IF
  
!  4  LOG-NORMAL with rg and sigma given

  ELSE IF ( idis == 4 ) THEN

      logrg = dlog(par(1))
      logsi = dabs(dlog(par(2)))
      logsi_inv = d_one / logsi
      C      = logsi_inv / root2p
      IF ( deriv(1) .and. deriv(2) ) then
        logC_d2 = - logsi_inv / par(2)
        fac_d1  =   logsi_inv / par(1)
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
          nwithr_d(i,1) = arg * fac_d1 * nwithr(i)
          nwithr_d(i,2) = logC_d2 * ( d_one - argsq ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
        END DO
      END IF

!  5 LOG-NORMAL with reff and veff given                               *

  ELSE IF ( idis == 5 ) THEN

      alpha1 = d_one + par(2)
      alpha2 = dlog(alpha1)
      rg     = par(1)/(d_one+par(2))**2.5_dp
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      IF ( deriv(1) .and. deriv(2) ) then
        logC_d2 = - d_half / alpha2 / alpha1
        fac_d1  =   logsi_inv / par(1)
        fac_d2a  =  - 2.5_dp * logsi_inv / alpha1
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
          nwithr_d(i,1) = arg * fac_d1 * nwithr(i)
          nwithr_d(i,2) = ( arg * fac_d2a + logC_d2*(d_one-argsq) ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
        END DO
      END IF

!  6 POWER LAW                               *

  ELSE IF ( idis == 6 ) THEN

      alpha = par(1)
      rmin  = par(2)
      rmax  = par(3)
      alpha1 = alpha - d_one
      fac1 = (rmax/rmin)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * rmax**alpha1 * fac2
      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        logC_d1 = (d_one/alpha1) + LOG(par(3)) - fac1 * fac2 * LOG(par(3)/par(2))
        DO i = 1, numradius
          r     = radius(i)
          if ( (r < rmax) .and. (r > rmin) ) then
            nwithr(i)    = C*r**(-alpha)
            nwithr_d(i,1) = ( logC_d1 - log(r) ) * nwithr(i)
          else
            nwithr(i)    = d_zero
            nwithr_d(i,1) = d_zero
          endif
        END DO
      ELSE
        DO i = 1, numradius
          r     = radius(i)
          if ( (r < rmax) .and. (r > rmin) ) then
            nwithr(i)    = C*r**(-alpha)
          else
            nwithr(i)    = d_zero
          endif
        END DO
      END IF

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
      CALL gammafunction ( aperg, deriv(1), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        logC_d1 = ( logb - dgammln ) * gamma1 + aperg/par(1)
        logC_d2 = - aperg * gamma / par(2)
        logC_d3 = gamma1 - aperg * ( logb - dgammln ) * gamma1 - aperg * (gamma1 + LOG(par(2)) ) 
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = b * r ** gamma
          nwithr(i) = EXP ( arg1 - r3 )
          nwithr_d(i,1) = ( logC_d1 + logr - r3/par(1) ) * nwithr(i)
          nwithr_d(i,2) = ( logC_d2 + r3*gamma/par(2) )  * nwithr(i)
          nwithr_d(i,3) = ( logC_d3 + r3*(gamma1+LOG(par(2))-logr) ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = b*r ** gamma
          nwithr(i) = EXP ( arg1 - r3 )
        END DO
      END IF

!  8 MODIFIED GAMMA with alpha, b and gamma given

  ELSE IF ( idis == 8 ) THEN

      alpha = par(1)
      b     = par(2)
      gamma = par(3)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      logb = LOG(b)
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv(1), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        b1      = b / par(1)
        b2      = b / par(2)
        logC_d1 = ( logb - dgammln ) * gamma1
        logC_d2 = aperg / b
        logC_d3 = gamma1 - aperg * logC_d1
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = r ** gamma
          nwithr(i) = EXP ( arg1 - b*r3 )
          nwithr_d(i,1) = ( logC_d1 + logr )       * nwithr(i)
          nwithr_d(i,2) = ( logC_d2 - r3 )         * nwithr(i)
          nwithr_d(i,3) = ( logC_d3 - b*logr*r3 )  * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = r ** gamma
          nwithr(i) = EXP ( arg1 - b*r3 )
        END DO
      END IF

  END IF

!  normal return

  RETURN

!  special return

240  CONTINUE
  faild = .TRUE.
  write(cdis,'(I1)')idis
  message = message_gamma(1:LEN(message_gamma))//', distribution : '//cdis
  RETURN

END SUBROUTINE sizedis_plus


SUBROUTINE develop_d ( max_Mie_angles, ncoeffs, nangles, nderivs, do_Mie_linearization, &
                       cosines, weights, FMAT, FMAT_D,                                  &
                       expcoeffs, expcoeffs_d )

!  Based on the Meerhoff Mie code
!   Linearization additions by R. Spurr, October 2004

!************************************************************************
!*  Calculate the expansion coefficients of the scattering matrix in    *
!*  generalized spherical functions by numerical integration over the   *
!*  scattering angle.  AND derivatives.                                 *
!************************************************************************

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_half, d_one, d_two, d_three, d_four

!  implicit none statement

  IMPLICIT NONE

!  input

  INTEGER          , INTENT (IN) :: max_Mie_angles
  LOGICAL          , INTENT (IN) :: do_Mie_linearization
  INTEGER          , INTENT (IN) :: ncoeffs, nangles, nderivs
 
  REAL    (KIND=dp), INTENT (IN) :: cosines(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: weights(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: FMAT(4,max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: FMAT_D(4,5,max_Mie_angles)

!  output

  REAL    (KIND=dp), INTENT (OUT) :: expcoeffs(6,0:max_Mie_angles)
  REAL    (KIND=dp), INTENT (OUT) :: expcoeffs_d(6,5,0:max_Mie_angles)

!  local variables

  REAL    (KIND=dp) :: P00(max_Mie_angles,2)
  REAL    (KIND=dp) :: P02(max_Mie_angles,2)
  REAL    (KIND=dp) :: P22(max_Mie_angles,2)
  REAL    (KIND=dp) :: P2m2(max_Mie_angles,2)
  REAL    (KIND=dp) :: fmatw(4,max_Mie_angles)
  REAL    (KIND=dp) :: fmatwd(4,5,max_Mie_angles)

  INTEGER           :: i, k, j, l, lnew, lold, itmp
  INTEGER           :: index_11, index_12, index_22, index_33, index_34, index_44 
  REAL    (KIND=dp) :: dl, dl2, qroot6, fac1, fac2, fac3, fl,&
                       sql4, sql41, twol1, tmp1, tmp2, denom, &
                       alfap, alfam, alfapd(5), alfamd(5)

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
  IF ( do_Mie_linearization ) THEN
    DO j = 0, ncoeffs
      DO k = 1, nderivs
        DO i = 1, 6
          expcoeffs_d(i,k,j) = d_zero
        END DO
      END DO
    END DO
  END IF

!  Multiply the scattering matrix F with the weights w for all angles  *
!  We do this here because otherwise it should be done for each l      *

  DO i = 1, 4
    DO j = 1, nangles
      fmatw(i,j) = weights(j)*FMAT(i,j)
    END DO
  END DO
  IF ( do_Mie_linearization ) THEN
    DO i = 1, 4
      DO k = 1, nderivs
        DO j = 1, nangles
          fmatwd(i,k,j) = weights(j)*FMAT_D(i,k,j)
        END DO
      END DO
    END DO
  END IF

!  Start loop over the coefficient index l                             *
!  first update generalized spherical functions, then calculate coefs. *
!  lold and lnew are pointer-like indices used in recurrence           *

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

    dl   = DBLE(l)

    IF (l == 0) THEN

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
    IF ( do_Mie_linearization ) THEN
      DO k = 1, nderivs
        alfapd(k) = d_zero
        alfamd(k) = d_zero
      END DO
    END IF

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

    IF ( do_Mie_linearization ) THEN
      DO k = 1, nderivs
        DO i=1, nangles
          expcoeffs_d(index_11,k,l) = expcoeffs_d(index_11,k,l) + P00(i,lnew)*fmatwd(1,k,i)
          alfapd(k) = alfapd(k) + P22(i,lnew)  * (fmatwd(1,k,i)+fmatwd(3,k,i))
          alfamd(k) = alfamd(k) + P2m2(i,lnew) * (fmatwd(1,k,i)-fmatwd(3,k,i))
          expcoeffs_d(index_44,k,l) = expcoeffs_d(index_44,k,l) + P00(i,lnew)*fmatwd(3,k,i)
          expcoeffs_d(index_12,k,l) = expcoeffs_d(index_12,k,l) + P02(i,lnew)*fmatwd(2,k,i)
          expcoeffs_d(index_34,k,l) = expcoeffs_d(index_34,k,l) + P02(i,lnew)*fmatwd(4,k,i)
        END DO
        expcoeffs_d(index_11,k,l) =  fl*expcoeffs_d(index_11,k,l)
        expcoeffs_d(index_22,k,l) =  fl*d_half*(alfapd(k)+alfamd(k))
        expcoeffs_d(index_33,k,l) =  fl*d_half*(alfapd(k)-alfamd(k))
        expcoeffs_d(index_44,k,l) =  fl*expcoeffs_d(index_44,k,l)
        expcoeffs_d(index_12,k,l) =  fl*expcoeffs_d(index_12,k,l)
        expcoeffs_d(index_34,k,l) =  fl*expcoeffs_d(index_34,k,l)
      END DO
    END IF

  END DO

  RETURN
END SUBROUTINE develop_d




SUBROUTINE expand_d ( max_Mie_angles, ncoeffs, nangles, nderivs, do_Mie_linearization, &
                      cosines, expcoeffs, expcoeffs_d, FMAT, FMAT_D )

!  Based on the Meerhoff Mie code
!   Linearization additions by R. Spurr, November 2004

!  Use the expansion coefficients of the scattering matrix in
!  generalized spherical functions to exapnd F matrix and derivative

!  modules

  USE Mie_precision
  USE MIE_constants,    ONLY : d_zero, d_one, d_two, d_four

!  implicit none statement

  IMPLICIT NONE

!  input

  INTEGER          , INTENT (IN) :: max_Mie_angles
  LOGICAL          , INTENT (IN) :: do_Mie_linearization
  INTEGER          , INTENT (IN) :: ncoeffs, nangles, nderivs
 
  REAL    (KIND=dp), INTENT (IN) :: cosines(max_Mie_angles)

  REAL    (KIND=dp), INTENT (IN) :: expcoeffs(6,0:max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: expcoeffs_d(6,5,0:max_Mie_angles)

!  output

  REAL    (KIND=dp), INTENT (OUT) :: FMAT(4,max_Mie_angles)
  REAL    (KIND=dp), INTENT (OUT) :: FMAT_D(4,5,max_Mie_angles)

!  local variables

  REAL    (KIND=dp) :: P00(max_Mie_angles,2)
  REAL    (KIND=dp) :: P02(max_Mie_angles,2)

  INTEGER           :: i, k, j, l, lnew, lold, itmp
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
  IF ( do_Mie_linearization ) THEN
     DO k = 1, nderivs
        DO i = 1, nangles
           DO j = 1, 4
              FMAT_D(j,k,i) = d_zero
           END DO
        END DO
     END DO
  END IF

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

     IF ( do_Mie_linearization ) THEN
        DO k = 1, nderivs
           DO i = 1, nangles
              FMAT_D(1,k,i) = FMAT_D(1,k,i) + expcoeffs_d(index_11,k,l)*P00(i,lnew)
              FMAT_D(2,k,i) = FMAT_D(2,k,i) + expcoeffs_d(index_12,k,l)*P02(i,lnew)
              FMAT_D(3,k,i) = FMAT_D(3,k,i) + expcoeffs_d(index_44,k,l)*P00(i,lnew)
              FMAT_D(4,k,i) = FMAT_D(4,k,i) + expcoeffs_d(index_34,k,l)*P02(i,lnew)
           END DO
        END DO
     ENDIF

  END DO

  RETURN
END SUBROUTINE expand_d

