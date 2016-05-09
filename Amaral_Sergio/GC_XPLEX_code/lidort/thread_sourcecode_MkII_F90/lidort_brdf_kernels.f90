!$Id: lidort_brdf_kernels.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
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
! #            LAMBERTIAN_FUNCTION                              #
! #            ROSSTHIN_FUNCTION                                #
! #            ROSSTHICK_FUNCTION                               #
! #            LISPARSE_FUNCTION                                #
! #            LIDENSE_FUNCTION                                 #
! #            ROUJEAN_FUNCTION                                 #
! #            HAPKE_FUNCTION                                   #
! #            RAHMAN_FUNCTION                                  #
! #            COXMUNK_FUNCTION                                 #
! #            COXMUNK_FUNCTION_DB                              #
! #                                                             #
! #  These two kernels introduced for Version 3.4R              #
! #                                                             #
! #            BREONVEG_FUNCTION                                #
! #            BREONSOIL_FUNCTION                               #
! #                                                             #
! ###############################################################

SUBROUTINE ROSSTHIN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHIN_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: ROSSTHIN_KERNEL

!  Local variables

      REAL(kind=8)  :: DS1, DS2, CKSI, SKSI, KSI, FUNC
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      ROSSTHIN_KERNEL = ZERO
      XPHI = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_KERNEL = FUNC - PIO2

!  Finish

      RETURN
END SUBROUTINE ROSSTHIN_FUNCTION

!

SUBROUTINE ROSSTHICK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHICK_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: ROSSTHICK_KERNEL

!  Local variables

      REAL(kind=8)  :: DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      ROSSTHICK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_KERNEL = FUNC - PIO4

!  Finish

      RETURN
END SUBROUTINE ROSSTHICK_FUNCTION

!

SUBROUTINE ROUJEAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROUJEAN_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: ROUJEAN_KERNEL

!  Local variables

      REAL(kind=8)  :: DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      REAL(kind=8)  :: XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      ROUJEAN_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      XPHI_C = XPHI
      IF ( XPHI .GT. PIE )  XPHI_C = TWO*PIE - XPHI
      IF ( XPHI .LT. ZERO ) XPHI_C = - XPHI

      IF ( SXI .LT. ZERO ) THEN
        XPHI_R  = ( PIE - XPHI_C )
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
        TXI =  - ( SXI / XI )
      ELSE
        TXI =   ( SXI / XI )
        XPHI_R  = XPHI_C
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
      ENDIF

      TXJ =  ( SXJ / XJ )
      DS1 = TWO * TXJ * TXI
      DS2 = TXJ + TXI
      DS3 = TXJ*TXJ  + TXI*TXI
      PHIFAC = ( ( PIE - XPHI_R ) * CXPHI_R + SXPHI_R ) / PI4
      S1 = PHIFAC * DS1
      S2 = ( DS2 + DSQRT ( DS3 - DS1 * CXPHI_R ) ) / PIE
      ROUJEAN_KERNEL = S1 - S2

!  Finish

      RETURN
END SUBROUTINE ROUJEAN_FUNCTION

!

SUBROUTINE LISPARSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LISPARSE_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: LISPARSE_KERNEL

!  local variables

      REAL(kind=8)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=8)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=8)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=8)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL = ZERO
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

!  contributions P and R

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

!  Finish

      RETURN
END SUBROUTINE LISPARSE_FUNCTION

!

SUBROUTINE LIDENSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LIDENSE_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: LIDENSE_KERNEL

!  local variables

      REAL(kind=8)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=8)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=8)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=8)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL = ZERO
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

!  contributions P and R

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

!  Finish

      RETURN
END SUBROUTINE LIDENSE_FUNCTION

!

SUBROUTINE HAPKE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: HAPKE_KERNEL

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

!  local variables

      REAL(kind=8)  :: CTHETA, THETA, PHASE
      REAL(kind=8)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(kind=8)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(kind=8)  :: HELP_J, TERM_J, HELP_I, TERM_I
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. ONE ) CTHETA = ONE
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
      TERM_J   = ( ONE + HELP_J ) / ( ONE + HELP_J * GAMMA )
      HELP_I   = TWO * XI
      TERM_I   = ( ONE + HELP_I ) / ( ONE + HELP_I * GAMMA )

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TERM_J * TERM_I - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION
 
!  Finish

      RETURN
END SUBROUTINE HAPKE_FUNCTION

!

SUBROUTINE RAHMAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, RAHMAN_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: RAHMAN_KERNEL

!  Revision. 24 October 2007.
!  --------------------------

!    * Limiting cases and hotspot evaluation.
!    * Revision based on the DISORT_2 code
!    *  In Disort, this kernel is known as the RPV^ BRDF.

!     The RPV reference is:
!       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
!       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
!       With NOAA Advanced Very High Resolution Radiometer Data,
!       J. Geophys. Res., 98, 20791-20801.

!  The hotspot should occur when XI = XJ and PHI = 180.

!  local variables

      REAL(kind=8)  :: T_INC, T_REF, DT1, DT2 
      REAL(kind=8)  :: CXI, DELTA, K1_SQ, FACT
      REAL(kind=8)  :: GEOM, PHASE, RFAC, K0, K1, K2
      REAL(kind=8)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT
      REAL(kind=8), PARAMETER  :: SMALL = 1.0d-04

!  Initial section
!  ---------------

!  Initialise output

      RAHMAN_KERNEL = ZERO

!  Limiting case, formerly
!      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  Azimuth convettion

      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( TWO - K0 )
      FACT = FACT * ( ONE - K1 ) / ( ONE + K1 ) / ( ONE + K1 )
      GEOM = ( TWO * XJ * XJ * XJ ) ** ( K2 - ONE )
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
      FACT  = ( ONE + K1_SQ + TWO * K1 * CXI ) ** ONEP5
      PHASE = ( ONE - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      RFAC = ( ONE - K0 ) / ( ONE + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * ( ONE + RFAC ) * GEOM

!  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

!  Finish

      RETURN
END SUBROUTINE RAHMAN_FUNCTION

!

SUBROUTINE HAPKE_FUNCTION_OLD  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: HAPKE_KERNEL

!  Hapke Kernel function.
!   From DISORT code; used as validation.

!  local variables

      REAL(kind=8)  :: XPHI, CKPHI
      REAL             MU, MUP, DUMMY, DPHI, HAPKEKER
      EXTERNAL         HAPKEKER

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      DUMMY = ZERO
      HAPKE_KERNEL = HAPKEKER( DUMMY, DUMMY, MU, MUP, DPHI )

!  Finish

      RETURN
END SUBROUTINE HAPKE_FUNCTION_OLD

!

REAL FUNCTION hapkeker( WVNMLO, WVNMHI, MU, MUP, DPHI )

      implicit none

!      Supplies surface bi-directional reflectivity.
!
!      NOTE 1: Bidirectional reflectivity in DISORT is defined
!              by Eq. 39 in STWL.
!      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
!              angles) are positive.
!
!  INPUT:
!
!    WVNMLO : Lower wavenumber (inv cm) of spectral interval
!
!    WVNMHI : Upper wavenumber (inv cm) of spectral interval
!
!    MU     : Cosine of angle of reflection (positive)
!
!    MUP    : Cosine of angle of incidence (positive)
!
!    DPHI   : Difference of azimuth angles of incidence and reflection
!                (radians)
!
!  LOCAL VARIABLES:
!
!    IREF   : bidirectional reflectance options
!             1 - Hapke's BDR model
!
!    B0     : empirical factor to account for the finite size of
!             particles in Hapke's BDR model
!
!    B      : term that accounts for the opposition effect
!             (retroreflectance, hot spot) in Hapke's BDR model
!
!    CTHETA : cosine of phase angle in Hapke's BDR model
!
!    GAMMA  : albedo factor in Hapke's BDR model
!
!    H0     : H( mu0 ) in Hapke's BDR model
!
!    H      : H( mu ) in Hapke's BDR model
!
!    HH     : angular width parameter of opposition effect in Hapke's
!             BDR model
!
!    P      : scattering phase function in Hapke's BDR model
!
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model
!
!    W      : single scattering albedo in Hapke's BDR model
!
!
!   Called by- DREF, SURFAC
! +-------------------------------------------------------------------+
!     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
!     ..
!     .. Local Scalars ..

      INTEGER   IREF
      REAL      B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
!     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

!                              ** Hapke's BRDF model (times Pi/Mu0)
!                              ** (Hapke, B., Theory of reflectance
!                              ** and emittance spectroscopy, Cambridge
!                              ** University Press, 1993, Eq. 8.89 on
!                              ** page 233. Parameters are from
!                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5 * COS( DPHI )
         THETA = ACOS( CTHETA )

         P    = 1. + 0.5 * CTHETA

         HH   = 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )

         W = 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

         hapkeker = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END IF

      RETURN
END FUNCTION hapkeker

!

SUBROUTINE COXMUNK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, COXMUNK_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: COXMUNK_KERNEL

!  Critical exponent taken out

      REAL(kind=8), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(kind=8)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(kind=8)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(kind=8)  :: XPHI, CKPHI
      REAL(kind=8)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(kind=8)  :: SHADOWI, SHADOWR, SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL = ZERO
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
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A = PIO2 - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB  / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/DSQRT(ONE-XXI)
      T1   = DEXP(-DCOT*DCOT*S2)
      T2   = DERFC(DCOT*S3)
      SHADOWI = HALF*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/DSQRT(ONE-XXJ)
      T1   = DEXP(-DCOT*DCOT*S2)
      T2   = DERFC(DCOT*S3)
      SHADOWR = HALF*(S1*T1/DCOT-T2)

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW
 
!     Finish

      RETURN
END SUBROUTINE COXMUNK_FUNCTION

!

SUBROUTINE COXMUNK_FUNCTION_DB  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, COXMUNK_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: COXMUNK_KERNEL

!  local variables
!  ---------------

!  help variables

      integer          n, k, i, i1, N_phiquad_HALF
      REAL(kind=8)  :: XM, SXM, sum_pr
      REAL(kind=8)  :: sumr, w_p, reflec_0, reflec_1
      REAL(KIND=8)  :: phi_sub1, cphi_sub1, sphi_sub1
      REAL(KIND=8)  :: phi_sub2, cphi_sub2, sphi_sub2

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

!  Safety first zeroing

      REFLEC_0 = ZERO
      REFLEC_1 = ZERO
      COXMUNK_KERNEL = ZERO

!  Air to water, Polar quadrature

      n_muquad = 20
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 40
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

      CALL COXMUNK_FUNCTION &
         ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, REFLEC_0 )

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
          CALL COXMUNK_FUNCTION &
   ( MAXPARS, NPARS, PARS, XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2, R0_OUT_QUAD(K,N) )
          CALL COXMUNK_FUNCTION &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, R0_QUAD_IN(K,N) )
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

!      write(34,'(1p6e14.5)')reflec_0, reflec_1, coxmunk_kernel, &
!            dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
END SUBROUTINE COXMUNK_FUNCTION_DB   

!

SUBROUTINE BREONVEG_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BREONVEG_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: BREONVEG_KERNEL

!  Local variables

      REAL(kind=8)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, FP
      REAL(kind=8)  :: XPHI, CKPHI, ATTEN, PROJECTIONS
      REAL(kind=8)  :: sgamma, cgamma, calpha, calpha_sq, salpha
      REAL(kind=8)  :: PLEAF, GS, GV, FP0

!  Data coefficients

      REAL(kind=8)  :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098, 0.011187479, 0.043329567, 0.19262991/
   
!  F-.M. Breon vegetation model (2009)

!  Initialise

      BREONVEG_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

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
!   PARS(1) = refractive index squared

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(1) * Z2
      H2 = DSQRT ( PARS(1) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = HALF * ( RP*RP + RL*RL )

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)
      
!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
 
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = dsqrt(one - calpha_sq)

! Projection of leaf surface to outgoing direction

      gv = PLAGIOPHILE_COEFFS(1) + xi * &
          (PLAGIOPHILE_COEFFS(2) + xi * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xj * &
          (PLAGIOPHILE_COEFFS(2) + xj * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))
      
! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0d0 * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25d0 * PLEAF * FP / xi / xj / PROJECTIONS

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BREONVEG_KERNEL = Fp0 * atten

!     Finish

      RETURN
END SUBROUTINE BREONVEG_FUNCTION

!

SUBROUTINE BREONSOIL_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BREONSOIL_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: BREONSOIL_KERNEL

!  Local variables

      REAL(kind=8)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, FP
      REAL(kind=8)  :: XPHI, CKPHI, ATTEN, FP1
      REAL(kind=8)  :: sgamma, cgamma, calpha, calpha_sq, salpha
   
!  F-.M. Breon Soil model (2009)

!  Initialise

      BREONSOIL_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

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
!   PARS(1) = refractive index squared

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(1) * Z2
      H2 = DSQRT ( PARS(1) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = HALF * ( RP*RP + RL*RL )

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)
      
!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
 
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = dsqrt(one - calpha_sq)

! Polarization model for soil

      Fp1 = 0.25d0 * FP / xi / xj

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BREONSOIL_KERNEL = Fp1 * atten

!     Finish

      RETURN
END SUBROUTINE BREONSOIL_FUNCTION

!

SUBROUTINE LAMBERTIAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LAMBERTIAN_KERNEL )

      implicit none

!  include file of constants

      INCLUDE 'LIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=8), intent(out) :: LAMBERTIAN_KERNEL
   
!  Lambertian kernel

      LAMBERTIAN_KERNEL = ONE
  
!  Finish

      RETURN
END SUBROUTINE LAMBERTIAN_FUNCTION
