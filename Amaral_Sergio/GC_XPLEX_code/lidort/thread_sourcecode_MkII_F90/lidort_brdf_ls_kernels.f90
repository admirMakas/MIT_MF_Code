!$Id: lidort_brdf_ls_kernels.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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
! #            LISPARSE_FUNCTION_PLUS                           #
! #            LIDENSE_FUNCTION_PLUS                            #
! #            HAPKE_FUNCTION_PLUS                              #
! #            RAHMAN_FUNCTION_PLUS                             #
! #            COXMUNK_FUNCTION_PLUS                            #
! #            COXMUNK_FUNCTION_PLUS_DB                         #
! #                                                             #
! ###############################################################

SUBROUTINE HAPKE_FUNCTION_PLUS                       &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
              HAPKE_KERNEL, HAPKE_DERIVATIVES )

      implicit none

!  include file of constants amd dimensions

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: HAPKE_KERNEL
      REAL(kind=8), intent(out) :: HAPKE_DERIVATIVES ( MAXPARS )

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's model
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  Local variables

      INTEGER       :: J
      REAL(kind=8)  :: CTHETA, THETA, PHASE
      REAL(kind=8)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(kind=8)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(kind=8)  :: HELP_J, GHELP_J, TERM_J
      REAL(kind=8)  :: HELP_I, GHELP_I, TERM_I
      REAL(kind=8)  :: TI_TJ, DT1, DT2
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL       = ZERO
      DO J = 1, NPARS
        HAPKE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION

!  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_KERNEL / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
        HAPKE_DERIVATIVES(1) = DT1 + DT2 * REFLEC
      ENDIF

!  Hotspot  derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_DERIVATIVES(2) = DT1 * REFLEC * PHASE
      ENDIF

!  empirical factor derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        DT1 = B_HOT / B0_EMPIR 
        HAPKE_DERIVATIVES(3) = DT1 * REFLEC * PHASE
      ENDIF

!  Finish

      RETURN
END SUBROUTINE HAPKE_FUNCTION_PLUS

!

SUBROUTINE LISPARSE_FUNCTION_PLUS                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
              LISPARSE_KERNEL, LISPARSE_DERIVATIVES )

      implicit none

!  include file of constants amd dimensions

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: LISPARSE_KERNEL
      REAL(kind=8), intent(out) :: LISPARSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER       :: J
      REAL(kind=8)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=8)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=8)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=8)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(kind=8)  :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL       = ZERO
      DO J = 1, NPARS
        LISPARSE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_DERIVATIVES(1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_DERIVATIVES(2) = HALF * DX_P - DX_QR
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LISPARSE_FUNCTION_PLUS

!

SUBROUTINE LIDENSE_FUNCTION_PLUS                     &
             ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
               XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
               LIDENSE_KERNEL, LIDENSE_DERIVATIVES )

      implicit none

!  include file of constants amd dimensions

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: LIDENSE_KERNEL
      REAL(kind=8), intent(out) :: LIDENSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER       :: J
      REAL(kind=8)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=8)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=8)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=8)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(kind=8)  :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL       = ZERO
      DO J = 1, NPARS
        LIDENSE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LIDENSE_DERIVATIVES(1) = - P_QR * DY_Q / Q 
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_P_QR = ( DX_P / P ) - ( DX_R / R ) - ( DX_Q / Q )
        LIDENSE_DERIVATIVES(2) = P_QR * DX_P_QR
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDENSE_FUNCTION_PLUS

!

SUBROUTINE RAHMAN_FUNCTION_PLUS                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
              RAHMAN_KERNEL, RAHMAN_DERIVATIVES )

      implicit none

!  include file of constants amd dimensions

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: RAHMAN_KERNEL
      REAL(kind=8), intent(out) :: RAHMAN_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER       :: J
      REAL(kind=8)  :: T_INC, T_REF, DT1, DT2, D_FACT, D_HELPM
      REAL(kind=8)  :: CXI, DELTA, K1_SQ, FACT, D_K0, D_K1, D_K2
      REAL(kind=8)  :: GEOM, PHASE, RFAC, RFAC1, K0, K1, K2, HELPR
      REAL(kind=8)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT, HELPG, HELPM
      REAL(kind=8), PARAMETER  :: SMALL = 1.0d-04

!  Initialise

      RAHMAN_KERNEL = ZERO
      DO J = 1, NPARS
        RAHMAN_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0d0 - K0 )
      FACT = FACT * ( 1.0d0 - K1 ) / ( 1.0d0 + K1 ) / ( 1.0d0 + K1 )
      GEOM = ( 2.0d0 * XJ * XJ * XJ ) ** ( K2 - 1.0d0 )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are 
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(PHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      HELPM = ONE - K1_SQ 
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVES(1) = RAHMAN_KERNEL * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_DERIVATIVES(2) = RAHMAN_KERNEL * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVES(3) = RAHMAN_KERNEL * D_K2
      ENDIF

!  Finish

      RETURN
END SUBROUTINE RAHMAN_FUNCTION_PLUS

!

SUBROUTINE COXMUNK_FUNCTION_PLUS                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
              COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: COXMUNK_KERNEL
      REAL(kind=8), intent(out) :: COXMUNK_DERIVATIVES ( MAXPARS )


!  Critical exponent taken out

      REAL(kind=8), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      INTEGER       :: J
      REAL(kind=8)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(kind=8)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(kind=8)  :: XPHI, CKPHI, T1_R, T2_R, DCOT_R, T1_I, T2_I, DCOT_I
      REAL(kind=8)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(kind=8)  :: SHADOWI, SHADOWR, SHADOW
      REAL(kind=8)  :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      REAL(kind=8)  :: D_S1, D_S2, D_T1, D_T2
      REAL(kind=8)  :: D_SHADOWI, D_SHADOWR, D_SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL     = ZERO
      DO J = 1, NPARS
        COXMUNK_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2                 
      B = ( XI + XJ ) / A                  
      IF ( B .GT. ONE ) B = ONE           
      A  = PIO2 - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA            
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_DERIVATIVES(1) = - COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  square refractive index derivative
!  --This section of code was formerly at the end of routine
!    -- Now moved here before the shadowing option
!    -- otherwise derivative will not get done
!         Bug found by V. Natraj in VLIDORT. 02 February 2007.

      IF ( DO_DERIV_PARS(2) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_DERIVATIVES(2) = COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code
!  -----------

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
      T2_I   = DERFC(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
      T2_R   = DERFC(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!  Update Scalar derivatives
!  -------------------------

!  add the shadow derivative to inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I 
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )
        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
        COXMUNK_DERIVATIVES(1) = COXMUNK_DERIVATIVES(1) * SHADOW + &
                                 COXMUNK_KERNEL * D_SHADOW / SHADOW
      ENDIF

!  Refractive index derivative, update

      IF ( DO_DERIV_PARS(2) ) THEN
        COXMUNK_DERIVATIVES(2) = COXMUNK_DERIVATIVES(2) * SHADOW 
      ENDIF

!  Finish

      RETURN
END SUBROUTINE COXMUNK_FUNCTION_PLUS

!

SUBROUTINE COXMUNK_FUNCTION_PLUS_DB                 &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
              COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: COXMUNK_KERNEL
      REAL(kind=8), intent(out) :: COXMUNK_DERIVATIVES ( MAXPARS )

!  local variables
!  ---------------

!  help variables

      integer       :: n, k, i, i1, q, N_phiquad_HALF
      REAL(kind=8)  :: XM, SXM, sum_pr, pr, dfunc(2)
      REAL(kind=8)  :: sumr, w_p, reflec_0, reflec_1
      REAL(KIND=8)  :: phi_sub1, cphi_sub1, sphi_sub1
      REAL(KIND=8)  :: phi_sub2, cphi_sub2, sphi_sub2
      REAL(KIND=8)  :: d_reflec_0(3), d_reflec_1(3)

!  Local quadrature stuff

      integer, parameter :: max_msrs_muquad = 40
      integer, parameter :: max_msrs_phiquad = 50
      integer            :: n_muquad, n_phiquad

!  arrays

      REAL(kind=8)  :: X_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: W_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: SX_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: WXX_MUQUAD(max_msrs_muquad)

      REAL(kind=8)  :: X_PHIQUAD (max_msrs_phiquad)
      REAL(kind=8)  :: W_PHIQUAD (max_msrs_phiquad)

      REAL(kind=8)  :: R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: D_R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad,2)
      REAL(kind=8)  :: D_R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad,2)

!  Safety first zeroing

      REFLEC_0 = ZERO
      REFLEC_1 = ZERO
      COXMUNK_KERNEL = ZERO

      DO Q = 1, NPARS
        IF ( DO_DERIV_PARS(Q) ) THEN
          COXMUNK_DERIVATIVES(Q) = ZERO
          D_REFLEC_0(Q) = ZERO
          D_REFLEC_1(Q) = ZERO
        ENDIF
      ENDDO

!  Air to water, Polar quadrature

      n_muquad = 40
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 100
      N_phiquad_HALF = N_PHIQUAD / 2
      CALL GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
      DO I = 1, N_PHIQUAD_HALF
        I1 = I + N_PHIQUAD_HALF
        X_PHIQUAD(I1) = - X_PHIQUAD(I)
        W_PHIQUAD(I1) =   W_PHIQUAD(I)
      ENDDO
      DO I = 1, N_PHIQUAD
        X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
      ENDDO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION_PLUS                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
              REFLEC_0, D_REFLEC_0 )

!  Quadrature output for first order R/T calculations 

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          SPHI_SUB1 = DSIN(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          SPHI_SUB2 = DSIN(PHI_SUB2)
          CALL COXMUNK_FUNCTION_PLUS                              &
             ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,               &
               XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,  &
               R0_OUT_QUAD(K,N), DFUNC )
          D_R0_OUT_QUAD(K,N,1) = DFUNC(1)
          D_R0_OUT_QUAD(K,N,2) = DFUNC(2)
          CALL COXMUNK_FUNCTION_PLUS                              &
             ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,               &
               XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,  &
               R0_QUAD_IN(K,N),  DFUNC )
          D_R0_QUAD_IN(K,N,1) = DFUNC(1)
          D_R0_QUAD_IN(K,N,2) = DFUNC(2)
        ENDDO
      ENDDO

!  compute the next order

      SUMR = ZERO
      DO K = 1, n_muquad
        SUM_PR = ZERO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          SUM_PR = SUM_PR + W_P * R0_QUAD_IN(K,N) * R0_OUT_QUAD(K,N)
        ENDDO
        SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
      ENDDO
      REFLEC_1 = SUMR

!  Compute total

      COXMUNK_KERNEL = REFLEC_0 + REFLEC_1

!  Derivatives

      DO Q = 1, 2
       IF ( DO_DERIV_PARS(Q) ) THEN
        SUMR = ZERO
        DO K = 1, n_muquad
         PR = ZERO
         DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          PR = PR + W_P *   R0_QUAD_IN(K,N)   * D_R0_OUT_QUAD(K,N,Q) &
                  + W_P * D_R0_QUAD_IN(K,N,Q) *   R0_OUT_QUAD(K,N)
         ENDDO
         SUMR = SUMR + PR * WXX_MUQUAD(K)
        ENDDO
        D_REFLEC_1(Q) = SUMR
        COXMUNK_DERIVATIVES(Q) = D_REFLEC_0(Q) + D_REFLEC_1(Q)
       ENDIF
      ENDDO

!      write(34,'(1p6e14.5)')reflec_0, reflec_1, pars(1), &
!         dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
END SUBROUTINE COXMUNK_FUNCTION_PLUS_DB
