!$Id: RTS_mie_modules.f90,v 1.1 2010/07/30 23:47:04 daven Exp $
!  Contains the following modules
!       MODULE Mie_precision
!       MODULE Mie_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Mie_precision

  IMPLICIT NONE

! Define KIND variables for single and double precision

  INTEGER, PARAMETER :: sprec  = KIND( 1.0 )
  INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

END MODULE Mie_precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Mie_constants

  use Mie_precision, ONLY: dp

  IMPLICIT NONE

! Numbers such as 1, 2, pie, etc.......

  REAL (KIND=dp), PARAMETER :: d_zero  = 0.0_dp
  REAL (KIND=dp), PARAMETER :: d_one   = 1.0_dp
  REAL (KIND=dp), PARAMETER :: d_two   = 2.0_dp
  REAL (KIND=dp), PARAMETER :: d_three = 3.0_dp
  REAL (KIND=dp), PARAMETER :: d_four  = 4.0_dp
  REAL (KIND=dp), PARAMETER :: d_half  = 0.5_dp

  COMPLEX (KIND=dp), PARAMETER :: c_zero = ( d_zero, d_zero )
  COMPLEX (KIND=dp), PARAMETER :: c_one =  ( d_one,  d_zero )
  COMPLEX (KIND=dp), PARAMETER :: c_i =    ( d_zero,  d_one )

!  Mie version

  CHARACTER (LEN=10), PARAMETER :: Mie_f90_version = 'F90_Mie_V1'

END MODULE Mie_constants
