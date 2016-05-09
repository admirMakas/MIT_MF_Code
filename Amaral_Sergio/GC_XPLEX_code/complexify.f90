!******************************************************************************
!	Written for 'complexify.py 1.3'
!	J.R.R.A.Martins 1999
!       21-Apr-00  Fixed tan, sinh, cosh
!                  sign now returns complex
!                  added log10 and nint
!                  changed ==, /= and >= -- see comments below
!       20-May-00  added cosd, sind, and epsilon
!       11-Jul-00  took away cosd, sind (they are reserved, but not
!                  intrinsic functions in F90)
!       21-Jul-00  converted all trig functions to the value/derivative
!                  formulas -- not general complex number formulas
!       15-Aug-00  Fixed bug in atan formula and added the rest of the
!                  _ci and _ic cominations to the relational operators.
!                  P. Sturdza
!                  
!******************************************************************************
!
! Assume all code is compiled with double precision (-r8 compiler flag)
!

!TODO:
!     more typ combinations: cc, cr, rc, ic ?
!     check all fcns
!

module complexify
  USE MYTYPE
  implicit none
! ABS
  interface abs
     module procedure abs_c
     module procedure abs_xxxx
  end interface
! ACOS
  interface acos
     module procedure acos_c
  end interface
! ASIN
  interface asin
     module procedure asin_c
  end interface

! ATAN
  interface atan
     module procedure atan_c
  end interface

! ATAN2
  interface atan2
     module procedure atan2_cc
  end interface

! COSH
  interface cosh
     module procedure cosh_c
  end interface

! MAX (limited to 2-4 complex args, 2 mixed args)
  interface max
     module procedure max_cc
     module procedure max_cr
     module procedure max_rc
     module procedure max_ccc     ! added because of DFLUX.f
     module procedure max_cccc     ! added because of DFLUX.f
  end interface

! MIN (limited to 2-4 complex args, 2 mixed args)
  interface min
     module procedure min_cc
     module procedure min_cr
     module procedure min_rc
     module procedure min_ccc
     module procedure min_cccc
  end interface

! SIGN
  interface sign
     module procedure sign_cc
     module procedure sign_cr
     module procedure sign_rc
  end interface

! DIM
  interface dim
     module procedure dim_cc
     module procedure dim_cr
     module procedure dim_rc
  end interface

! SINH
  interface sinh
     module procedure sinh_c
  end interface
  
! TAN
  interface tan
     module procedure tan_c
  end interface
  
! TANH
  interface tanh
     module procedure tanh_c
  end interface

! LOG10
  interface log10
     module procedure log10_c
  end interface

! NINT
  interface nint
     module procedure nint_c
  end interface
! ISNAN
  interface ISNAN
     module procedure isnan_c
  end interface
!MINVAL
  interface minval
     module procedure minval_c
     module procedure minval_2c
     module procedure minval_3c
     module procedure minval_4c
  end interface

!MAXVAL
  interface maxval
     module procedure maxval_c
     module procedure maxval_2c
     module procedure maxval_3c
     module procedure maxval_4c
  end interface

! TINY
  interface tiny
     module procedure tiny_c
  end interface

! MOD 
  interface mod
     module procedure mod_c
     module procedure mod_c2
  end interface

! <
  interface operator (<)
     module procedure lt_cc
     module procedure lt_cr
     module procedure lt_rc
     module procedure lt_ci
     module procedure lt_ic
  end interface

! <=
  interface operator (<=)
     module procedure le_cc
     module procedure le_cr
     module procedure le_rc
     module procedure le_ci
     module procedure le_ic
     module procedure le_xc
     module procedure le_cx
  end interface

! >
  interface operator (>)
     module procedure gt_cc
     module procedure gt_cr
     module procedure gt_rc
     module procedure gt_ci
     module procedure gt_ic
  end interface
! >= 
  interface operator (>=)
     module procedure ge_cc
     module procedure ge_cr
     module procedure ge_rc
     module procedure ge_ci
     module procedure ge_ic
  end interface

! ==
  interface operator (==)
     module procedure eq_cc
     module procedure eq_cr
     module procedure eq_rc
     module procedure eq_ci
  end interface
  interface operator (.ceq.)
     module procedure eq_cc
     module procedure eq_rr
     module procedure eq_ii
     module procedure eq_aa
     module procedure eq_cr
     module procedure eq_rc
     module procedure eq_ci
     module procedure eq_ic
     module procedure eq_ir
     module procedure eq_ri
  end interface
  interface operator (.ne.)
     module procedure ne_cc
     module procedure ne_cr
     module procedure ne_rc
     module procedure ne_ci
     module procedure ne_ic
  end interface
  interface assignment (=)
     module procedure eq_xr4
     module procedure eq_x1r41
     module procedure eq_x2r42
     module procedure eq_x3r43
     module procedure eq_x4r44
     module procedure eq_x5r45
     module procedure eq_x6r46
     module procedure eq_x7r47
     module procedure eq_x1r81
     module procedure eq_x1r8
     module procedure eq_x1i
     module procedure eq_x2i
     module procedure eq_x3i
     module procedure eq_i2x2
     module procedure eq_i1x1
     module procedure eq_x2r8
     module procedure eq_x3r8
     module procedure eq_x4r8
     module procedure eq_x5r8
     module procedure eq_x2r82
     module procedure eq_x3r83
     module procedure eq_x4r84
     module procedure eq_x5r85
     module procedure eq_x6r86
     module procedure eq_x7r87
     module procedure eq_xr8
     module procedure eq_xi
     module procedure eq_x1i1
     module procedure eq_xx
     module procedure eq_x1x1
     module procedure eq_x2x2
     module procedure eq_x3x3
     module procedure eq_x4x4
     module procedure eq_ix
  end interface
  interface operator (-)
    module procedure sub_xx
    module procedure sub_x1x1
    module procedure sub_x2x2
    module procedure sub_x2_xscal
    module procedure sub_x3_xscal
    module procedure sub_x1_iscal
    module procedure sub_x1_rscal
    module procedure sub_xr
    module procedure sub_rx 
    module procedure sub_x
    module procedure sub_x3
    module procedure sub_xi
    module procedure sub_ix
    module procedure sub_xxi
    module procedure sub_x2i
  end interface 
  interface operator (/)
    module procedure div_xx
    module procedure div_x1x1
    module procedure div_x2x2
    module procedure div_x3x3
    module procedure div_x3i3
    module procedure div_xr
    module procedure div_xi
    module procedure div_rx
    module procedure div_x_xscal
    module procedure div_x_rscal
    module procedure div_xx_scal
    module procedure div_xx_rscal
    module procedure div_xx_iscal
    module procedure div_xxx_scal
    module procedure div_xxx_iscal
    module procedure div_xxxx_iscal
    module procedure div_xxx_rscal
    module procedure div_xr_scal
    module procedure div_ix
  end interface
  interface operator (**)
    module procedure pow_xi
    module procedure pow_x2_i
    module procedure pow_rx
    module procedure pow_xx
    module procedure pow_xr
    module procedure pow_x1_r
  end interface
  interface operator (*)
    module procedure mult_xx
    module procedure mult_x1x1
    module procedure mult_x2x2
    module procedure mult_x3x3
    module procedure mult_x4x4
    module procedure mult_xr
    module procedure mult_rx
    module procedure mult_xi
    module procedure mult_ix
    module procedure mult_x_xscal
    module procedure mult_x_rscal
    module procedure mult_xx_xscal
    module procedure mult_xx_rscal
    module procedure mult_xxx_rscal
    module procedure mult_iscal_xxx
    module procedure mult_iscal_xx
    module procedure mult_xxx_xscal
    module procedure mult_xscal_x
    module procedure mult_rscal_x
    module procedure mult_xscal_xx
    module procedure mult_xscal_xxx
    module procedure mult_rscal_xx
  end interface
  interface operator (+)
    module procedure add_x
    module procedure add_xx
    module procedure add_x1x1
    module procedure add_x2x2
    module procedure add_xscal_x2
    module procedure add_x2_xscal
    module procedure add_x3_xscal
    module procedure add_x2i2
    module procedure add_x3x3
    module procedure add_x4x4
    module procedure add_x2r
    module procedure add_x1r
    module procedure add_xr
    module procedure add_xi
    module procedure add_rx
    module procedure add_ix
  end interface 
  interface exp
    module procedure exp_x
    module procedure exp_x1
  end interface
  interface log
    module procedure log_x
    module procedure log_x4
    module procedure log_x1
  end interface
  interface int
    module procedure int_x
  end interface
  interface sin
    module procedure sin_x
  end interface
  interface cos
    module procedure cos_x
  end interface
  interface sum
    module procedure sum_x
    module procedure sum_xx
    module procedure sum_xxx
    module procedure sum_xxxx
  end interface 
  interface sqrt
    module procedure sqrt_x
  end interface
  interface ceiling
    module procedure ceiling_x
  end interface 
  interface floor
    module procedure floor_x
  end interface
  ! XPLX
      interface XPLX
        module procedure xplx_i
        module procedure xplx_i2
        module procedure xplx_r
        module procedure xplx_r2r2
        module procedure xplx_r1r1
        module procedure xplx_r4r4
        module procedure xplx_x
      end interface

contains

!******************************************************************************
!
!   Function definitions
!
!******************************************************************************

! FLOOR
  integer function floor_x(z)
   type (xplex) :: z    
   floor_x = floor(z%r)
  end function floor_x  
! CEILING
  integer function ceiling_x(z)
   type (xplex) :: z
   ceiling_x = ceiling(z%r)
  end function ceiling_x
! SQRT
  type (xplex) function sqrt_x(z)
    type (xplex) :: z
    real*8 ::zr,zi
    zr=z%r
    zi=z%i
    sqrt_x%r = sqrt(zr)
    sqrt_x%i = 0d0!zi/2 * ( zr**(-0.5)  )
    return
  end function sqrt_x
! SUM
  type (xplex) function sum_x(z)
    type (xplex),dimension(:)::z
    sum_x%r = sum(z%r)
    sum_x%i = 0d0!sum(z%i)
    return
  end function sum_x  
  type (xplex) function sum_xx(z)
    type (xplex),dimension(:,:)::z
    integer ::i,j
    sum_xx%r = sum(z%r)
    sum_xx%i = 0d0!sum(z%i)
    return
  end function sum_xx
  type (xplex) function sum_xxx(z)
    type (xplex),dimension(:,:,:)::z
    integer ::i,j,k
    sum_xxx%r = sum(z%r)
    sum_xxx%i = 0d0!sum(z%i)
    return
  end function sum_xxx
  type (xplex) function sum_xxxx(z)
    type (xplex),dimension(:,:,:,:)::z
    integer ::i,j,k,l
    sum_xxxx%r = sum(z%r)
    sum_xxxx%i = 0d0!sum(z%i)
    return
  end function sum_xxxx
! SIN
  type (xplex) function sin_x(z)
    type (xplex) :: z,iz
    real*8 :: zr,zi
    zr = z%r
    zi = z%i
    sin_x%r =sin(zr)
    sin_x%i =0d0!cos(zr)*zi
    return
  end function sin_x
! COS
  type (xplex) function cos_x(z)
    type (xplex) :: z
    real*8 :: zr,zi
    zr=z%r
    zi=z%i
    cos_x%r =cos(zr)
    cos_x%i =0d0!-zi*sin(zr)
    return
  end function cos_x
! INT
  integer function int_x(z)
    type (xplex),intent(in) ::z
    int_x = int(z%r)
    return
  end function int_x
!===============================
  subroutine eq_ix(int1,z)
    type (xplex),intent(in)::z
    integer,intent(inout)::int1
    int1 = int(z%r)
    return
  end subroutine eq_ix
  subroutine eq_xi(z,int1)
    type (xplex),intent(inout)::z
    integer,intent(in)::int1
    z%r = dble(int1)
    z%i = 0d0
    return
  end subroutine eq_xi
  subroutine eq_xr4(z,r)
    type (xplex),intent(inout)::z
    real*4,intent(in)::r
    z%r = dble(r)
    z%i = 0d0
    return
  end subroutine eq_xr4
  subroutine eq_x1r41(z,r)
    type (xplex),intent(inout)::z(:)
    real*4,intent(in)::r(:)
    z(:)%r = dble(r(:))
    z(:)%i = 0d0
    return
  end subroutine eq_x1r41
  subroutine eq_x1r81(z,r)
    type (xplex),intent(inout)::z(:)
    real*8,intent(in)::r(:)
    z(:)%r = dble(r(:))
    z(:)%i = 0d0
    return
  end subroutine eq_x1r81
  subroutine eq_x1i1(z,r)
    type (xplex),intent(inout)::z(:)
    integer,intent(in)::r(:)
    z(:)%r = dble(r(:))
    z(:)%i = 0d0
    return
  end subroutine eq_x1i1
  subroutine eq_x1r8(z,r)
    type (xplex),intent(inout)::z(:)
    real*8,intent(in)::r
    z(:)%r = r
    z(:)%i = 0d0
  end subroutine eq_x1r8
  subroutine eq_x1i(z,r)
    type (xplex),intent(inout)::z(:)
    integer,intent(in)::r
    z(:)%r = r
    z(:)%i = 0d0
    return
  end subroutine eq_x1i
  subroutine eq_x2r42(z,r)
    type (xplex),intent(inout)::z(:,:)
    real*4,intent(in)::r(:,:)
    z(:,:)%r = dble(r(:,:))
    z(:,:)%i = 0d0
    return
  end subroutine eq_x2r42
  subroutine eq_x2r82(z,r)
    type (xplex),intent(inout)::z(:,:)
    real*8,intent(in)::r(:,:)
    z(:,:)%r = dble(r(:,:))
    z(:,:)%i = 0d0
    return
  end subroutine eq_x2r82
  subroutine eq_x2i(z,r)
    type (xplex),intent(inout)::z(:,:)
    integer,intent(in)::r
    z(:,:)%r = r
    z(:,:)%i = 0d0
    return
  end subroutine eq_x2i
  subroutine eq_i2x2(int1,z)
    integer,intent(inout)::int1(:,:)
    type (xplex),intent(in)::z(:,:)
    int1(:,:) = z(:,:)%r
    return
  end subroutine eq_i2x2
  subroutine eq_i1x1(int1,z)
    integer,intent(inout)::int1(:)
    type (xplex),intent(in)::z(:)
    int1(:) = z(:)%r
    return
  end subroutine eq_i1x1
  subroutine eq_x2r8(z,r)
    type (xplex),intent(inout)::z(:,:)
    real*8,intent(in)::r
    z(:,:)%r = r
    z(:,:)%i = 0d0
    return
  end subroutine eq_x2r8
  subroutine eq_x3r43(z,r)
    type (xplex),intent(inout)::z(:,:,:)
    real*4,intent(in)::r(:,:,:)
    z(:,:,:)%r = dble(r(:,:,:))
    z(:,:,:)%i = 0d0 
    return
  end subroutine eq_x3r43
  subroutine eq_x3r8(z,r)
    type (xplex),intent(inout)::z(:,:,:)
    real*8,intent(in)::r
    z(:,:,:)%r = r
    z(:,:,:)%i = 0d0
    return
  end subroutine eq_x3r8
  subroutine eq_x3i(z,r)
    type (xplex),intent(inout)::z(:,:,:)
    integer,intent(in)::r
    z(:,:,:)%r = r
    z(:,:,:)%i = 0d0
    return
  end subroutine eq_x3i
  subroutine eq_x3r83(z,r)
    type (xplex),intent(inout)::z(:,:,:)
    real*8,intent(in)::r(:,:,:)
    z(:,:,:)%r = dble(r(:,:,:))
    z(:,:,:)%i = 0d0
    return
  end subroutine eq_x3r83
  subroutine eq_x4r44(z,r)
    type (xplex),intent(inout)::z(:,:,:,:)
    real*4,intent(in)::r(:,:,:,:)
    z(:,:,:,:)%r = dble(r(:,:,:,:))
    z(:,:,:,:)%i = 0d0
    return
  end subroutine eq_x4r44
  subroutine eq_x4r84(z,r)
    type (xplex),intent(inout)::z(:,:,:,:)
    real*8,intent(in)::r(:,:,:,:)
    z(:,:,:,:)%r = dble(r(:,:,:,:))
    z(:,:,:,:)%i = 0d0
    return
  end subroutine eq_x4r84
  subroutine eq_x4r8(z,r)
    type (xplex),intent(inout)::z(:,:,:,:)
    real*8,intent(in)::r
    integer ::i,j,k,l
    z(:,:,:,:)%r = r
    z(:,:,:,:)%i = 0d0
    return
  end subroutine eq_x4r8
  subroutine eq_x5r8(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:)
    real*8,intent(in)::r
    z(:,:,:,:,:)%r = r
    z(:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x5r8
  subroutine eq_x5r45(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:)
    real*4,intent(in)::r(:,:,:,:,:)
    z(:,:,:,:,:)%r = dble(r(:,:,:,:,:))
    z(:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x5r45
  subroutine eq_x5r85(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:)
    real*8,intent(in)::r(:,:,:,:,:)
    z(:,:,:,:,:)%r = dble(r(:,:,:,:,:))
    z(:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x5r85
  subroutine eq_x6r46(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:,:)
    real*4,intent(in)::r(:,:,:,:,:,:)
    z(:,:,:,:,:,:)%r = dble(r(:,:,:,:,:,:))
    z(:,:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x6r46
  subroutine eq_x6r86(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:,:)
    real*8,intent(in)::r(:,:,:,:,:,:)
    z(:,:,:,:,:,:)%r = dble(r(:,:,:,:,:,:))
    z(:,:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x6r86
  subroutine eq_x7r47(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:,:,:)
    real*4,intent(in)::r(:,:,:,:,:,:,:)
    z(:,:,:,:,:,:,:)%r = dble(r(:,:,:,:,:,:,:))
    z(:,:,:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x7r47
  subroutine eq_x7r87(z,r)
    type (xplex),intent(inout)::z(:,:,:,:,:,:,:)
    real*8,intent(in)::r(:,:,:,:,:,:,:)
    z(:,:,:,:,:,:,:)%r = dble(r(:,:,:,:,:,:,:))
    z(:,:,:,:,:,:,:)%i = 0d0
    return
  end subroutine eq_x7r87
  subroutine eq_xr8(z,r)
    type (xplex),intent(inout)::z
    real*8,intent(in)::r
    z%r = r
    z%i = 0d0
  end subroutine eq_xr8
  subroutine eq_xx(z1,z2)
    type (xplex),intent(inout)::z1
    type (xplex),intent(in)::z2
    z1%r = z2%r
    z1%i = z2%i
  end subroutine eq_xx
  subroutine eq_x1x1(z1,z2)
    type (xplex),intent(inout)::z1(:)
    type (xplex),intent(in)::z2(:)
    z1(:)%r = z2(:)%r
    z1(:)%i = z2(:)%i
    return
  end subroutine eq_x1x1
  subroutine eq_x2x2(z1,z2)
    type (xplex),intent(inout)::z1(:,:)
    type (xplex),intent(in)::z2(:,:)
    z1(:,:)%r = z2(:,:)%r
    z1(:,:)%i = z2(:,:)%i
    return
  end subroutine eq_x2x2
  subroutine eq_x3x3(z1,z2)
    type (xplex),intent(inout)::z1(:,:,:)
    type (xplex),intent(in)::z2(:,:,:)
    z1(:,:,:)%r = z2(:,:,:)%r
    z1(:,:,:)%i = z2(:,:,:)%i
    return
  end subroutine eq_x3x3
  subroutine eq_x4x4(z1,z2)
    type (xplex),intent(inout)::z1(:,:,:,:)
    type (xplex),intent(in)::z2(:,:,:,:)
    z1(:,:,:,:)%r = z2(:,:,:,:)%r
    z1(:,:,:,:)%i = z2(:,:,:,:)%i
    return
  end subroutine eq_x4x4
! LOG
  type (xplex) function log_x(z)
    type (xplex) :: z
    real*8 :: zr,zi
    zr = z%r
    zi = z%i
    log_x%r = log(zr)
    log_x%i =0d0! sign(max(abs(zi/zr),0d0),zi*zr)
    return
  end function log_x
  function log_x4(z) result(C)
    type (xplex),intent(in) :: z(:,:,:,:)
    real*8 :: zr(size(z,1),size(z,2),size(z,3),size(z,4))
    real*8 :: zi(size(z,1),size(z,2),size(z,3),size(z,4)) 
    real*8 :: theta1(size(z,1),size(z,2),size(z,3),size(z,4))
    type(xplex)::C(size(z,1),size(z,2),size(z,3),size(z,4))
    zr(:,:,:,:) = z(:,:,:,:)%r
    zi(:,:,:,:) = z(:,:,:,:)%i
    theta1(:,:,:,:) = zi(:,:,:,:)* (zr(:,:,:,:)**-1.)
    C(:,:,:,:)%r = log(zr(:,:,:,:))
    C(:,:,:,:)%i =0d0! sign(max(abs(zi/zr),0d0),zi*zr)
    return
  end function log_x4
  function log_x1(z) result(C)
    type (xplex):: z(:)
    real*8 :: zr(size(z,1)),zi(size(z,1)),theta1(size(z,1))
    type(xplex)::C(size(z,1))
    zr(:) = z(:)%r
    zi(:) = z(:)%i
    theta1(:) = zi(:) * ( zr(:)**-1.)
    C(:)%r = log(zr(:))
    C(:)%i = 0d0!sign(max(abs(zi/zr),0d0),zi*zr)
    return
  end function log_x1
! EXP
  type (xplex) function exp_x(z)
    type (xplex)::z
    real*8 :: temp1,temp2,temp3,temp4,temp5,zr,zi
    zr = z%r
    zi = z%i
    temp1 = exp(zr)
    temp2 = 1d0
    temp3 = zi
    temp4 = temp1*temp2
    temp5 = max(temp1*temp3,0d0)
    exp_x = xplex(temp4,0d0)!temp5)
    return
  end function exp_x
  function exp_x1(z) result(C)
    type (xplex)::z(:)
    type(xplex):: C(size(z,1))
    real*8 :: temp1(size(z,1)),temp2(size(z,1)),temp3(size(z,1))
    real*8 :: temp4(size(z,1)),temp5(size(z,1)),zr(size(z,1))
    real*8 :: zi(size(z,1))
    zr(:) = z(:)%r
    zi(:) = z(:)%i
    temp1(:) = exp(zr(:))
    temp2(:) = 1d0
    temp3(:) = zi(:)
    temp4(:) = temp1(:)*temp2(:)
    temp5(:) = max(temp1*temp3,0d0)
    C(:)%r = temp4(:)
    C(:)%i = 0d0!temp5(:)
    return
  end function exp_x1
! ** ** ** ** ** ** ** ** ** ** ** ** **
  type (xplex) function pow_xi(z,int1)
    type (xplex),intent(in) :: z
    integer ,intent(in)::int1
    integer::i
    pow_xi%r = z%r**int1
    pow_xi%i = 0d0!z%r**(int1-1) * int1 * z%i
    return
   end function pow_xi
   function pow_x2_i(z,int1) result(C)
    type (xplex),intent(in) :: z(:,:)
    integer,intent(in) ::int1
    type (xplex):: C(size(z,1),size(z,2))
    real*8 :: zr(size(z,1),size(z,2)),zi(size(z,1),size(z,2))
     C(:,:)%r=z(:,:)%r**int1
     C(:,:)%i=0d0!z(:,:)%i**(int1-1) * int1 * z(:,:)%i
    return
   end function pow_x2_i
   type (xplex) function pow_rx(r,z)
    type (xplex), intent(in) :: z
    real*8,intent(in) ::r
    real*8 :: z1r,z1i,z2r,z2i
    z1r = r
    z2r = 0d0
    z2r = z%r
    z2i = z%i
    pow_rx%r = z1r**z2r !* cos(z2i*log(z1r))
    pow_rx%i = 0d0!z1r**z2r * sin(z2i*log(z1r))
    return
   end function pow_rx
   type (xplex) function pow_xx(z1,z2)
    type (xplex), intent(in) :: z1,z2
    real*8 :: z1i,z1r,z2i,z2r
    z1r = z1%r
    z1i = z1%i
    z2r = z2%r
    z2i = z2%i 
    pow_xx%r = z1r**z2r
    pow_xx%i = 0d0!z1r**(z2r-1d0)*z1i
    return
   end function pow_xx
  type (xplex) function pow_xr(z1,z2)
    type (xplex), intent(in) :: z1
    real*8, intent(in)::z2
    real*8 :: z1i,z1r,z2i,z2r
    z1r = z1%r
    z1i = z1%i
    z2r = z2
    z2i = 0d0
    pow_xr%r = z1r**z2r
    pow_xr%i = 0d0!z1r**(z2r-1) * z2r * z1i
    return
   end function pow_xr
   function pow_x1_r(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:)
    type (xplex):: C(size(z1,1))
    integer::i
    real*8, intent(in)::z2
    real*8 :: z1i(size(z1,1)),z1r(size(z1,1)),z2i,z2r
    z1r(:) = z1(:)%r
    z1i(:) = z1(:)%i
    z2r    = z2
    z2i    = 0d0
    C(:)%r = z1r(:)**z2r
    C(:)%i = 0d0!z1r(:)**(z2r-1) * z2r * z1i(:)
    return
   end function pow_x1_r 
! -------------------------------------------------------
  type (xplex) function sub_x(z1)
    type (xplex),intent(in)::z1
    sub_x%r = -z1%r
    sub_x%i = 0d0!-z1%i
   end function sub_x
  function sub_x3(z1) result(C)
    type (xplex),intent(in)::z1(:,:,:)
    type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
    C(:,:,:)%r = -z1(:,:,:)%r
    C(:,:,:)%i = 0d0!-z1(:,:,:)%i
    return
   end function sub_x3
  type (xplex) function sub_xx(z1,z2)
    type (xplex), intent(in) :: z1,z2
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = z2%r
    z2i = z2%i
    temp1 = z1r-z2r
    temp2 = z1i-z2i
    sub_xx%r =temp1
    sub_xx%i =0d0!temp2
    return
  end function sub_xx
  function sub_x2x2(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:,:),z2(:,:)
    integer::i,j
    type (xplex):: C(size(z1,1),size(z1,2))
    C(:,:)%r =z1(:,:)%r-z2(:,:)%r
    C(:,:)%i =0d0!z1(:,:)%i-z2(:,:)%i
    return
  end function sub_x2x2
  function sub_x1x1(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:),z2(:)
    type (xplex):: C(size(z1,1))
    C(:)%r =z1(:)%r-z2(:)%r
    C(:)%i =0d0!z1(:)%i-z2(:)%i
    return
  end function sub_x1x1
  function sub_x2_xscal(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:,:)
    type (xplex), intent(in) ::z2
    type (xplex):: C(size(z1,1),size(z1,2))
    C(:,:)%r =z1(:,:)%r-z2%r
    C(:,:)%i =0d0!z1(:,:)%i-z2%i
    return
  end function sub_x2_xscal
  function sub_x1_iscal(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:)
    integer, intent(in) ::z2
    type (xplex):: C(size(z1,1))
    C(:)%r =z1(:)%r - dble(z2)
    C(:)%i =0d0!z1(:)%i - 0d0
    return
  end function sub_x1_iscal
  function sub_x1_rscal(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:)
    real*8, intent(in) ::z2
    type (xplex):: C(size(z1,1))
    C(:)%r =z1(:)%r-z2
    C(:)%i =0d0!z1(:)%i-0d0
    return
  end function sub_x1_rscal

  function sub_x3_xscal(z1,z2) result(C)
    type (xplex), intent(in) :: z1(:,:,:)
    type (xplex), intent(in) ::z2
    type (xplex):: C(size(z1,1),size(z1,2),size(z1,3))
    C(:,:,:)%r =z1(:,:,:)%r-z2%r
    C(:,:,:)%i =0d0!z1(:,:,:)%i-z2%i
    return
  end function sub_x3_xscal
  function sub_x2i(z1,int1) result(C)
    type (xplex),intent(in) :: z1(:,:)
    integer, intent(in):: int1
    integer :: C(size(z1,1),size(z1,2))
    C(:,:) = z1(:,:)%r-dble(int1)
    return
  end function sub_x2i
  type (xplex) function sub_xr(z1,r)
    type (xplex), intent(in) :: z1
    real*8, intent(in) ::r
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = r
    z2i = 0d0
    temp1 = z1r-z2r
    temp2 = z1i-z2i
    sub_xr%r =temp1
    sub_xr%i =0d0!temp2
    return
  end function sub_xr
  type (xplex) function sub_xi(z1,r)
    type (xplex), intent(in) :: z1
    integer, intent(in) ::r
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = r
    z2i = 0d0
    temp1 = z1r-z2r
    temp2 = z1i-z2i
    sub_xi%r =temp1
    sub_xi%i =0d0!temp2
    return
  end function sub_xi
  type (xplex) function sub_ix(int1,z2)
    integer, intent(in) ::int1
    type (xplex), intent(in) :: z2
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = int1
    z1i = 0d0
    z2r = z2%r
    z2i = z2%i
    temp1 = z1r-z2r
    temp2 = z1i-z2i
    sub_ix%r =temp1
    sub_ix%i =0d0!temp2
    return
  end function sub_ix
  integer function sub_xxi(z1,r)
    type (xplex), intent(in) :: z1
    integer, intent(in) ::r
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = r
    z2i = 0d0
    temp1 = int(z1r-z2r)
    temp2 = int(z1i-z2i)
    sub_xxi =temp1
    return
  end function sub_xxi

  type (xplex) function sub_rx(r,z2)
    type (xplex), intent(in) :: z2
    real*8,intent(in) :: r
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = r
    z1i = 0d0
    z2r = z2%r
    z2i = z2%i
    temp1 = z1r-z2r
    temp2 = z1i-z2i
    sub_rx%r =temp1
    sub_rx%i =0d0!temp2
    return
  end function sub_rx

! /,intrinsic///////////////////////
  type (xplex) function div_xx(n,d)
    type (xplex),intent(in) :: n,d
    real*8 :: nr,ni,dr,di,r1,r2,theta1,theta2,RN,RD
    nr = n%r
    ni = n%i
    dr = d%r
    di = d%i
    div_xx%r = nr/dr
    div_xx%i = 0d0!ni!/dr! - (nr*di/dr)*(1d0/dr)
    return
  end function div_xx
  function div_x2x2(n,d) result(C)
    type (xplex), intent(in) :: n(:,:),d(:,:)
    real*8 :: nr(size(n,1),size(n,2)),ni(size(n,1),size(n,2))
    real*8 :: dr(size(n,1),size(n,2)),di(size(n,1),size(n,2))
    TYPE (XPLEX) :: C(size(n,1),size(n,2))
    nr(:,:) = n(:,:)%r
    ni(:,:) = n(:,:)%i
    dr(:,:) = d(:,:)%r
    di(:,:) = d(:,:)%i
    C(:,:)%r = nr(:,:)/dr(:,:)
    C(:,:)%i = 0d0!ni(:,:)!/dr(:,:)!-(nr(:,:)*di(:,:)/dr(:,:))*(1d0/dr(:,:))
    return
  end function div_x2x2
  function div_x1x1(n,d) result(C)
    type (xplex), intent(in) :: n(:),d(:)
    real*8 :: nr(size(n,1)),ni(size(n,1)),dr(size(n,1)),di(size(n,1))
    TYPE (XPLEX) :: C(size(n,1))
    nr(:) = n(:)%r
    ni(:) = n(:)%i
    dr(:) = d(:)%r
    di(:) = d(:)%i
    C(:)%r = nr(:)/dr(:)
    C(:)%i = 0d0!ni(:)!/dr(:)! - (nr(:)*di(:)/dr(:))*(1d0/dr(:))
    return
  end function div_x1x1
  function div_x3x3(n,d) result(C)
    type (xplex), intent(in) :: n(:,:,:),d(:,:,:)
    real*8 :: nr(size(n,1),size(n,2),size(n,3))
    real*8 :: ni(size(n,1),size(n,2),size(n,3))
    real*8 :: dr(size(n,1),size(n,2),size(n,3))
    real*8 :: di(size(n,1),size(n,2),size(n,3))
    TYPE (XPLEX) :: C(size(n,1),size(n,2),size(n,3))
    nr(:,:,:) = n(:,:,:)%r
    ni(:,:,:) = n(:,:,:)%i
    dr(:,:,:) = d(:,:,:)%r
    di(:,:,:) = d(:,:,:)%i
    C(:,:,:)%r = nr(:,:,:)/dr(:,:,:)
    C(:,:,:)%i = 0d0!ni(:,:,:)!/dr(:,:,:)!-(nr(:,:,:)*di(:,:,:)/dr(:,:,:)) &
                                    !*(1d0/dr(:,:,:))
    return
  end function div_x3x3
  function div_x3i3(n,d) result(C)
    type (xplex), intent(in) :: n(:,:,:)
    integer,intent(in)::d(:,:,:)
    real*8 :: nr(size(n,1),size(n,2),size(n,3))
    real*8 :: ni(size(n,1),size(n,2),size(n,3))
    real*8 :: dr(size(n,1),size(n,2),size(n,3))
    real*8 :: di(size(n,1),size(n,2),size(n,3))
    TYPE (XPLEX) :: C(size(n,1),size(n,2),size(n,3))
    nr(:,:,:) = n(:,:,:)%r
    ni(:,:,:) = n(:,:,:)%i
    dr(:,:,:) = d(:,:,:)
    di(:,:,:) = 0d0
    C(:,:,:)%r = nr(:,:,:)/dr(:,:,:)
    C(:,:,:)%i = 0d0!ni(:,:,:)!/dr(:,:,:)
    return
  end function div_x3i3
  type (xplex) function div_ix(n,d)
    type (xplex), intent(in) :: d
    integer , intent(in) ::n
    real*8 :: nr,ni,dr,di
    nr = n
    ni = 0d0
    dr = d%r
    di = d%i
    div_ix%r = nr/dr
    div_ix%i =0d0!-(nr*di/dr)*(1d0/dr)
    return
  end function div_ix
  function div_xx_scal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:)
    type (xplex), intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2))
    real*8 :: nr(size(n,1),size(n,2)),ni(size(n,1),size(n,2)),dr,di
    nr(:,:) = n(:,:)%r
    ni(:,:) = n(:,:)%i
    dr = d%r
    di = d%i
    C(:,:)%r = nr(:,:)/dr
    C(:,:)%i = 0d0!ni(:,:)!/dr!-(nr(:,:)*di/dr)*(1d0/dr)
    return
  end function div_xx_scal
  function div_x_xscal(n,d) result(C)
    type (xplex), intent(in) :: n(:)
    type (xplex), intent(in) :: d
    type (xplex) :: C(size(n,1))
    real*8 :: nr(size(n,1)),ni(size(n,1)),dr,di
    nr(:) = n(:)%r
    ni(:) = n(:)%i
    dr = d%r
    di = d%i
    C(:)%r = nr(:)/dr
    C(:)%i = 0d0!ni(:)!/dr!-(nr(:)*di/dr)*(1d0/dr)
    return
  end function div_x_xscal
  function div_x_rscal(n,d) result(C)
    type (xplex), intent(in) :: n(:)
    real*8, intent(in) :: d
    type (xplex) :: C(size(n,1))
    real*8 :: nr(size(n,1)),ni(size(n,1)),dr,di
    nr(:) = n(:)%r
    ni(:) = n(:)%i
    dr = d
    di = 0d0
    C(:)%r = nr(:)/dr
    C(:)%i = 0d0!ni(:)!/dr
    return
  end function div_x_rscal
  function div_xx_rscal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:)
    real*8, intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2))
    real*8 :: nr(size(n,1),size(n,2)),ni(size(n,1),size(n,2)),dr,di
    nr(:,:) = n(:,:)%r
    ni(:,:) = n(:,:)%i
    dr = d
    di = 0d0
    C(:,:)%r = nr(:,:)/dr
    C(:,:)%i = 0d0!ni(:,:)!/dr
    return
  end function div_xx_rscal
  function div_xx_iscal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:)
    integer, intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2))
    real*8 :: nr(size(n,1),size(n,2)),ni(size(n,1),size(n,2)),dr,di
    nr(:,:) = n(:,:)%r
    ni(:,:) = n(:,:)%i
    dr = d
    di = 0d0
    C(:,:)%r = nr(:,:)/dr
    C(:,:)%i = 0d0!ni(:,:)!/dr
    return
  end function div_xx_iscal
  function div_xxx_scal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:,:)
    type (xplex), intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2),size(n,3))
    real*8 :: nr(size(n,1),size(n,2),size(n,3))
    real*8 :: ni(size(n,1),size(n,2),size(n,3)),dr,di
    nr(:,:,:) = n(:,:,:)%r
    ni(:,:,:) = n(:,:,:)%i
    dr = d%r
    di = d%i
    C(:,:,:)%r = nr(:,:,:)/dr
    C(:,:,:)%i = 0d0!ni(:,:,:)!/di!-(nr(:,:,:)*di/dr)*(1d0/dr)
    return
  end function div_xxx_scal
  function div_xxx_iscal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:,:)
    integer, intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2),size(n,3))
    real*8 :: nr(size(n,1),size(n,2),size(n,3))
    real*8 :: ni(size(n,1),size(n,2),size(n,3)),dr,di
    nr(:,:,:) = n(:,:,:)%r
    ni(:,:,:) = n(:,:,:)%i
    dr = d
    di = 0d0
    C(:,:,:)%r = nr(:,:,:)/dr
    C(:,:,:)%i = 0d0!ni(:,:,:)!/dr
    return
  end function div_xxx_iscal
  function div_xxxx_iscal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:,:,:)
    integer, intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2),size(n,3),size(n,4))
    real*8 :: nr(size(n,1),size(n,2),size(n,3),size(n,4))
    real*8 :: ni(size(n,1),size(n,2),size(n,3),size(n,4)),dr,di
    nr(:,:,:,:) = n(:,:,:,:)%r
    ni(:,:,:,:) = n(:,:,:,:)%i
    dr = d
    di = 0d0
    C(:,:,:,:)%r = nr(:,:,:,:)/dr
    C(:,:,:,:)%i = 0d0!ni(:,:,:,:)!/dr
    return
  end function div_xxxx_iscal  
  function div_xxx_rscal(n,d) result(C)
    type (xplex), intent(in) :: n(:,:,:)
    real*8, intent(in) :: d
    type (xplex) :: C(size(n,1),size(n,2),size(n,3))
    real*8 :: nr(size(n,1),size(n,2),size(n,3))
    real*8 :: ni(size(n,1),size(n,2),size(n,3)),dr,di
    nr(:,:,:) = n(:,:,:)%r
    ni(:,:,:) = n(:,:,:)%i
    dr = d
    di = 0d0
    C(:,:,:)%r = nr(:,:,:)/dr
    C(:,:,:)%i = 0d0!ni(:,:,:)!/dr
    return
  end function div_xxx_rscal
  function div_xr_scal(n,d) result(C)
    type (xplex), intent(in) :: n(:)
    real*8, intent(in) :: d
    type (xplex):: C(size(n,1))
    real*8 :: nr(size(n,1)),ni(size(n,1)),dr,di
    nr(:) = n(:)%r
    ni(:) = n(:)%i
    dr = d
    di = 0d0
    C(:)%r = nr(:)/dr
    C(:)%i = 0d0!ni(:)!/dr
    return
  end function div_xr_scal
  type (xplex) function div_xi_scal(n,d)
    type (xplex), intent(in) :: n
    integer, intent(in) :: d
    integer :: i,j
    real*8 :: nr,ni,dr,di
    nr = n%r
    ni = n%i
    dr = dble(d)
    di = 0d0
    div_xi_scal%r = nr/dr
    div_xi_scal%i =0d0! ni!/dr
    return
  end function div_xi_scal
  type (xplex) function div_rx(n,d)
    type (xplex), intent(in) :: d
    real*8, intent(in)::n
    real*8 :: nr,ni,dr,di
    nr = n
    ni = 0d0
    dr = d%r
    di = d%i
    div_rx%r = nr/dr
    div_rx%i = 0d0!-(nr*di/dr)*(1d0/dr)
    return
  end function div_rx
  type (xplex) function div_xr(n,d)
    type (xplex), intent(in) :: n
    real*8,intent(in)::d
    real*8 :: nr,ni,dr,di
    nr = n%r
    ni = n%i
    dr = dble(d)
    di = 0d0
    div_xr%r = nr/dr
    div_xr%i = 0d0!ni!/dr
    return
  end function div_xr
  type (xplex) function div_xi(n,d)
    type (xplex), intent(in) :: n
    integer,intent(in)::d
    real*8 :: nr,ni,dr,di
    nr = n%r
    ni = n%i
    dr = dble(d)
    di = 0d0
    div_xi%r = nr/dr
    div_xi%i = 0d0!ni!/dr
    return
  end function div_xi

! ++++++++++++++++++++++++++++++++++
  type (XPLEX) function add_x(z)
    type (xplex), intent(in)::z
    add_x = z
    return
  end function add_x
  type (XPLEX) function add_xx(z1,z2)
    type (XPLEX), intent(in) :: z1,z2
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = z2%r
    z2i = z2%i
    temp1 = z1r+z2r
    temp2 = z1i+z2i
    add_xx%r = temp1
    add_xx%i = 0d0!temp2
    return  
    end function add_xx
  function add_x1x1(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:),z2(:)
    type (XPLEX):: C(size(z1,1))
    C(:)%r = z1(:)%r+z2(:)%r
    C(:)%i = 0d0!z1(:)%i+z2(:)%i
    return
    end function add_x1x1
  function add_x2x2(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:,:),z2(:,:)
    type (XPLEX):: C(size(z1,1),size(z1,2))
    C(:,:)%r = z1(:,:)%r+z2(:,:)%r
    C(:,:)%i = 0d0!z1(:,:)%i+z2(:,:)%i
    return
    end function add_x2x2
  function add_xscal_x2(z1,z2) result (C)
    type (XPLEX), intent(in) :: z2(:,:)
    type (XPLEX),intent(in) :: z1
    type (XPLEX):: C(size(z2,1),size(z2,2))
    C(:,:)%r = z1%r+z2(:,:)%r
    C(:,:)%i = 0d0!z1%i+z2(:,:)%i
    return
    end function add_xscal_x2
  function add_x2_xscal(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:,:)
    type (XPLEX),intent(in) :: z2
    type (XPLEX):: C(size(z1,1),size(z1,2))
    C(:,:)%r = z1(:,:)%r+z2%r
    C(:,:)%i = 0d0!z1(:,:)%i+z2%i
    return
    end function add_x2_xscal
  function add_x3_xscal(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:,:,:)
    type (XPLEX),intent(in) :: z2
    type (XPLEX):: C(size(z1,1),size(z1,2),size(z1,3))
    C(:,:,:)%r = z1(:,:,:)%r+z2%r
    C(:,:,:)%i = 0d0!z1(:,:,:)%i+z2%i
    return
    end function add_x3_xscal
  function add_x2i2(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:,:)
    integer,intent(in)::z2(:,:)
    type (XPLEX):: C(size(z1,1),size(z1,2))
    C(:,:)%r = z1(:,:)%r+z2(:,:)
    C(:,:)%i = 0d0!z1(:,:)%i+0d0
    return
    end function add_x2i2
  function add_x3x3(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:,:,:),z2(:,:,:)
    type (XPLEX)::C(size(z1,1),size(z1,2),size(z1,3))
    C(:,:,:)%r = z1(:,:,:)%r+z2(:,:,:)%r
    C(:,:,:)%i = 0d0!z1(:,:,:)%i+z2(:,:,:)%i
    return
    end function add_x3x3
  function add_x4x4(z1,z2) result (C)
    type (XPLEX), intent(in) :: z1(:,:,:,:),z2(:,:,:,:)
    type (XPLEX):: C(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
    C(:,:,:,:)%r = z1(:,:,:,:)%r+z2(:,:,:,:)%r
    C(:,:,:,:)%i = 0d0!z1(:,:,:,:)%i+z2(:,:,:,:)%i
    return
    end function add_x4x4
  function add_x2r(z1,r) result (C)
    type (XPLEX), intent(in) :: z1(:,:)
    real*8,intent(in)::r
    type (XPLEX):: C(size(z1,1),size(z1,2))
    C(:,:)%r = z1(:,:)%r+r
    C(:,:)%i = 0d0!z1(:,:)%i+0d0
    return
    end function add_x2r
  function add_x1r(z1,r) result (C)
    type (XPLEX), intent(in) :: z1(:)
    real*8,intent(in)::r
    type (XPLEX):: C(size(z1,1))
    C(:)%r = z1(:)%r+r
    C(:)%i = 0d0!z1(:)%i+0d0
    return
    end function add_x1r
  type (XPLEX) function add_xr(z1,z2)
    type (XPLEX), intent(in) :: z1
    real*8,intent(in) :: z2
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = z2
    z2i = 0d0
    temp1 = z1r+z2r
    temp2 = z1i+z2i
    add_xr%r = temp1
    add_xr%i = 0d0!temp2
    return
    end function add_xr
  type (XPLEX) function add_rx(z1,z2)
    type (XPLEX), intent(in) :: z2
    real*8,intent(in) :: z1
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1
    z1i = 0d0
    z2r = z2%r
    z2i = z2%i
    temp1 = z1r+z2r
    temp2 = z1i+z2i
    add_rx%r = temp1
    add_rx%i = 0d0!temp2
    return
    end function add_rx
  type (XPLEX) function add_xi(z1,z2)
    type (XPLEX), intent(in) :: z1
    integer,intent(in) :: z2
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1%r
    z1i = z1%i
    z2r = z2
    z2i = 0d0
    temp1 = z1r+z2r
    temp2 = z1i+z2i
    add_xi%r = temp1
    add_xi%i = 0d0!temp2
    return
    end function add_xi
  type (XPLEX) function add_ix(z1,z2)
    integer,intent(in) :: z1
    type (xplex),intent(in) :: z2
    real*8 :: z1r,z1i,z2r,z2i,temp1,temp2
    z1r = z1
    z1i = 0d0
    z2r = z2%r
    z2i = z2%i
    temp1 = z1r+z2r
    temp2 = z1i+z2i
    add_ix%r = temp1
    add_ix%i = 0d0!temp2
    return
    end function add_ix
! * * * * * * * * *
  function mult_xx_rscal(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:)
    real*8, intent(in):: z2
    type (XPLEX) ::C(size(z1,1),size(z1,2))
    real*8 :: z1r(size(z1,1),size(z1,2)),z1i(size(z1,1),size(z1,2))
    real*8 :: z2r,z2i
    z1r(:,:) = z1(:,:)%r
    z1i(:,:) = z1(:,:)%i
    z2r = z2
    z2i = 0d0
    C(:,:)%r= z1r(:,:)*z2r
    C(:,:)%i= 0d0!z1i(:,:)*z2r+z2i*z1r(:,:)
    return
    end function mult_xx_rscal
  function mult_x_rscal(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:)
    real*8, intent(in):: z2
    type (XPLEX) ::C(size(z1,1))
    real*8 :: z1r(size(z1,1)),z1i(size(z1,1)),z2r,z2i
    z1r(:) = z1(:)%r
    z1i(:) = z1(:)%i
    z2r = z2
    z2i = 0d0
    C(:)%r= z1r(:)*z2r
    C(:)%i= 0d0!z1i(:)*z2r+z2i*z1r(:)
    return
    end function mult_x_rscal
  function mult_x_xscal(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:)
    type (xplex), intent(in):: z2
    type (XPLEX) ::C(size(z1,1))
    real*8 :: z1r(size(z1,1)),z1i(size(z1,1)),z2r,z2i
    z1r(:) = z1(:)%r
    z1i(:) = z1(:)%i
    z2r = z2%r
    z2i = z2%i
    C(:)%r= z1r(:)*z2r
    C(:)%i= 0d0!z1i(:)*z2r+z2i*z1r(:)
    return
    end function mult_x_xscal
  function mult_xxx_rscal(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:,:)
    real*8, intent(in):: z2
    type (XPLEX) ::C(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z1r(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z1i(size(z1,1),size(z1,2),size(z1,3)),z2r,z2i
    z1r(:,:,:) = z1(:,:,:)%r
    z1i(:,:,:) = z1(:,:,:)%i
    z2r = z2
    z2i = 0d0
    C(:,:,:)%r= z1r(:,:,:)*z2r
    C(:,:,:)%i= 0d0!z1i(:,:,:)*z2r+z2i*z1r(:,:,:)
    return
    end function mult_xxx_rscal
  function mult_iscal_xxx(z1,z2) result(C)
    integer, intent(in) :: z1
    type (xplex), intent(in):: z2(:,:,:)
    type (XPLEX) ::C(size(z2,1),size(z2,2),size(z2,3))
    real*8 :: z1r,z1i,z2r(size(z2,1),size(z2,2),size(z2,3))
    real*8 :: z2i(size(z2,1),size(z2,2),size(z2,3))
    z1r = z1
    z1i = 0d0
    z2r(:,:,:) = z2(:,:,:)%r
    z2i(:,:,:) = z2(:,:,:)%i
    C(:,:,:)%r= z1r*z2r(:,:,:)
    C(:,:,:)%i= 0d0!z1i*z2r(:,:,:)+z2i(:,:,:)*z1r
    return
    end function mult_iscal_xxx
  function mult_iscal_xx(z1,z2) result(C)
    integer, intent(in) :: z1
    type (xplex), intent(in):: z2(:,:)
    type (XPLEX) ::C(size(z2,1),size(z2,2))
    real*8 :: z1r,z1i,z2r(size(z2,1),size(z2,2))
    real*8 :: z2i(size(z2,1),size(z2,2))
    z1r = z1
    z1i = 0d0
    z2r(:,:) = z2(:,:)%r
    z2i(:,:) = z2(:,:)%i
    C(:,:)%r= z1r*z2r(:,:)
    C(:,:)%i= 0d0!z1i*z2r(:,:)+z2i(:,:)*z1r
    return
    end function mult_iscal_xx
  function mult_xxx_xscal(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:,:)
    type (XPLEX), intent(in):: z2
    type (XPLEX) ::C(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z1r(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z1i(size(z1,1),size(z1,2),size(z1,3)),z2r,z2i
    z1r(:,:,:) = z1(:,:,:)%r
    z1i(:,:,:) = z1(:,:,:)%i
    z2r = z2%r
    z2i = z2%i
    C(:,:,:)%r= z1r(:,:,:)*z2r
    C(:,:,:)%i= 0d0!z1i(:,:,:)*z2r+z2i*z1r(:,:,:)
    return
    end function mult_xxx_xscal
  function mult_xscal_xxx(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1
    type (XPLEX), intent(in):: z2(:,:,:)
    type (XPLEX) ::C(size(z2,1),size(z2,2),size(z2,3))
    real*8 :: z1r,z1i,z2r(size(z2,1),size(z2,2),size(z2,3))
    real*8 :: z2i(size(z2,1),size(z2,2),size(z2,3))
    z1r = z1%r
    z1i = z1%i
    z2r(:,:,:) = z2(:,:,:)%r
    z2i(:,:,:) = z2(:,:,:)%i
    C(:,:,:)%r= z1r*z2r(:,:,:)
    C(:,:,:)%i= 0d0!z1i*z2r(:,:,:)+z2i(:,:,:)*z1r
    return
    end function mult_xscal_xxx
  function mult_xscal_xx(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1
    type (XPLEX), intent(in):: z2(:,:)
    type (XPLEX) ::C(size(z2,1),size(z2,2))
    real*8 :: z1r,z1i,z2r(size(z2,1),size(z2,2))
    real*8 :: z2i(size(z2,1),size(z2,2))
    z1r = z1%r
    z1i = z1%i
    z2r(:,:) = z2(:,:)%r
    z2i(:,:) = z2(:,:)%i
    C(:,:)%r= z1r*z2r(:,:)
    C(:,:)%i= 0d0!z1i*z2r(:,:)+z2i(:,:)*z1r
    return
    end function mult_xscal_xx
  function mult_xscal_x(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1
    type (XPLEX), intent(in):: z2(:)
    type (XPLEX) ::C(size(z2,1))
    real*8 :: z1r,z1i,z2r(size(z2,1)),z2i(size(z2,1))
    z1r = z1%r
    z1i = z1%i
    z2r(:) = z2(:)%r
    z2i(:) = z2(:)%i
    C(:)%r= z1r*z2r
    C(:)%i= 0d0!z1i*z2r(:)+z2i(:)*z1r
    return
    end function mult_xscal_x
  function mult_rscal_x(z1,z2) result(C)
    real*8, intent(in) :: z1
    type (XPLEX), intent(in):: z2(:)
    type (XPLEX) ::C(size(z2,1))
    real*8 :: z1r,z1i,z2r(size(z2,1)),z2i(size(z2,1))
    z1r = z1
    z1i = 0d0
    z2r(:) = z2(:)%r
    z2i(:) = z2(:)%i
    C(:)%r= z1r*z2r(:)
    C(:)%i= 0d0!z1i*z2r(:)+z2i(:)*z1r
    return
    end function mult_rscal_x
  function mult_rscal_xx(z1,z2) result(C)
    real*8,intent(in)::z1
    type (xplex), intent(in):: z2(:,:)
    type (XPLEX) ::C(size(z2,1),size(z2,2))
    real*8 :: z1r,z1i,z2r(size(z2,1),size(z2,2))
    real*8 :: z2i(size(z2,1),size(z2,2))
    z1r = z1
    z1i = 0d0
    z2r(:,:) = z2(:,:)%r
    z2i(:,:) = z2(:,:)%i
    C(:,:)%r= z1r*z2r(:,:)
    C(:,:)%i= 0d0!z1i*z2r(:,:)+z2i(:,:)*z1r
    return
    end function mult_rscal_xx
  function mult_xx_xscal(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:)
    type (xplex), intent(in):: z2
    type (XPLEX) ::C(size(z1,1),size(z1,2))
    real*8 :: z1r(size(z1,1),size(z1,2))
    real*8 :: z1i(size(z1,1),size(z1,2)),z2r,z2i
    z1r(:,:) = z1(:,:)%r
    z1i(:,:) = z1(:,:)%i
    z2r = z2%r
    z2i = z2%i
    C(:,:)%r= z1r(:,:)*z2r
    C(:,:)%i= 0d0!z1i(:,:)*z2r+z2i*z1r(:,:)
    return
    end function mult_xx_xscal
  type (XPLEX) function mult_xx(z1,z2)
    type (XPLEX), intent(in) :: z1,z2
    real*8 :: z1r,z1i,z2r,z2i,r1,r2,theta1,&
              theta2,RZ1,RZ2
    z1r = z1%r
    z1i = z1%i
    z2r = z2%r
    z2i = z2%i
    mult_xx%r= z1r*z2r
    mult_xx%i= 0d0!z1i*z2r+z2i*z1r
    return
    end function mult_xx
  function mult_x2x2(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:),z2(:,:)
    type (XPLEX):: C(size(z1,1),size(z1,2))
    real*8 :: z1r(size(z1,1),size(z1,2)),z1i(size(z1,1),size(z1,2))
    real*8 :: z2r(size(z1,1),size(z1,2)),z2i(size(z1,1),size(z1,2))
    z1r(:,:) = z1(:,:)%r
    z1i(:,:) = z1(:,:)%i
    z2r(:,:) = z2(:,:)%r
    z2i(:,:) = z2(:,:)%i
    C(:,:)%r= z1r(:,:)*z2r(:,:)
    C(:,:)%i= 0d0!z1i(:,:)*z2r(:,:)+z2i(:,:)*z1r(:,:)
    return
    end function mult_x2x2
  function mult_x1x1(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:),z2(:)
    type (XPLEX):: C(size(z1,1))
    real*8 :: z1r(size(z1,1)),z1i(size(z1,1))
    real*8 :: z2r(size(z1,1)),z2i(size(z1,1))
    z1r(:) = z1(:)%r
    z1i(:) = z1(:)%i
    z2r(:) = z2(:)%r
    z2i(:) = z2(:)%i
    C(:)%r= z1r(:)*z2r(:)
    C(:)%i= 0d0!z1i(:)*z2r(:)+z2i(:)*z1r(:)
    return
    end function mult_x1x1
  function mult_x3x3(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:,:),z2(:,:,:)
    type (XPLEX):: C(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z1r(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z1i(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z2r(size(z1,1),size(z1,2),size(z1,3))
    real*8 :: z2i(size(z1,1),size(z1,2),size(z1,3))
    z1r(:,:,:) = z1(:,:,:)%r
    z1i(:,:,:) = z1(:,:,:)%i
    z2r(:,:,:) = z2(:,:,:)%r
    z2i(:,:,:) = z2(:,:,:)%i
    C(:,:,:)%r= z1r(:,:,:)*z2r(:,:,:)
    C(:,:,:)%i= 0d0!z1i(:,:,:)*z2r(:,:,:)+z2i(:,:,:)*z1r(:,:,:)
    return
    end function mult_x3x3
  function mult_x4x4(z1,z2) result(C)
    type (XPLEX), intent(in) :: z1(:,:,:,:),z2(:,:,:,:)
    type (XPLEX):: C(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
    real*8 :: z1r(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
    real*8 :: z1i(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
    real*8 :: z2r(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
    real*8 :: z2i(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
    z1r(:,:,:,:) = z1(:,:,:,:)%r
    z1i(:,:,:,:) = z1(:,:,:,:)%i
    z2r(:,:,:,:) = z2(:,:,:,:)%r
    z2i(:,:,:,:) = z2(:,:,:,:)%i
    C(:,:,:,:)%r= z1r(:,:,:,:)*z2r(:,:,:,:)
    C(:,:,:,:)%i= 0d0!z1i(:,:,:,:)*z2r(:,:,:,:)+z2i(:,:,:,:)*z1r(:,:,:,:)
    return
    end function mult_x4x4
  type (XPLEX) function mult_xi(z1,z2)
    type (XPLEX), intent(in) :: z1
    integer, intent(in)::z2
    real*8 :: z1r,z1i,z2r,z2i
    z1r = z1%r
    z1i = z1%i
    z2r = z2
    z2i = 0d0
    mult_xi%r= z1r*z2r
    mult_xi%i= 0d0!z1i*z2r+z2i*z1r
    end function mult_xi
  type (XPLEX) function mult_ix(z1,z2)
    type (XPLEX), intent(in) :: z2
    integer, intent(in)::z1
    real*8 :: z1r,z1i,z2r,z2i
    z2r = z2%r
    z2i = z2%i
    z1r = z1
    z1i = 0d0
    mult_ix%r= z1r*z2r
    mult_ix%i= 0d0!z1i*z2r+z2i*z1r
    end function mult_ix
  type (XPLEX) function mult_xr(z1,z2)
    type (XPLEX), intent(in) :: z1
    real*8,intent(in)::z2
    real*8 :: z1r,z1i,z2r,z2i
    z1r = z1%r
    z1i = z1%i
    z2r = z2
    z2i = 0d0
    mult_xr%r=z1r*z2r
    mult_xr%i=0d0!z2r*z1i+z2i*z1r
    end function mult_xr
  type (XPLEX) function mult_rx(z1,z2)
    type (XPLEX), intent(in) :: z2
    real*8,intent(in) ::z1
    real*8 :: z1r,z1i,z2r,z2i
    z1r = z1
    z1i = 0d0
    z2r = z2%r
    z2i = z2%i
    mult_rx%r= z1r*z2r
    mult_rx%i= 0d0!z1i*z2r+z2i*z1r
    end function mult_rx

! ABS, intrinsic
  type (xplex) function abs_c(val)
    type(xplex), intent(in) :: val
    abs_c = val
    if (val%r < 0) then
    abs_c%r = -val%r
    abs_c%i = 0d0!-val%i
    endif
    return
  end function abs_c
  function abs_xxxx(val) result(C)
    type(xplex), intent(in) :: val(:,:,:,:)
    type(xplex)::C(size(val,1),size(val,3),size(val,3),size(val,4))
    integer::i,j,k,l
!!$OMP PARALLEL DO private(i,j,k,l)
    do l=1,size(val,4)
    do k=1,size(val,3)
    do j=1,size(val,2)
    do i=1,size(val,1)
    C(i,j,k,l)%r = val(i,j,k,l)%r
    C(i,j,k,l)%i = 0d0!val(i,j,k,l)%i
    if (val(i,j,k,l)%r < 0d0) then
    C(i,j,k,l)%r = -val(i,j,k,l)%r
    C(i,j,k,l)%i = 0d0!-val(i,j,k,l)%i
    endif
    enddo
    enddo
    enddo
    enddo
!!$OMP END PARALLEL DO
    return
  end function abs_xxxx

! ACOS
  type (xplex) function acos_c(z)
    type (xplex), intent(in) :: z
    real*8 :: zr,zi
    zr=z%r
    zi=z%i
    acos_c%r = acos(zr)
    acos_c%i = 0d0!-zi/sqrt(1d0-zr**2))
    return
  end function acos_c

! ASIN
  type (xplex) function asin_c(z)
    type (xplex), intent(in) :: z
    real*8 :: zr,zi
    zr=z%r
    zi=z%i
    asin_c%r =asin(zr)
    asin_c%i =0d0! zi/sqrt(1d0-zr**2))  
    return
  end function asin_c

! ATAN
  type (xplex) function atan_c(z)
    type (xplex), intent(in) :: z
    atan_c%r = atan(z%r)
    atan_c%i =  0d0!z%i/(1d0+z%r**2)) 
   return
  end function atan_c
  
! ATAN2
  type (xplex) function atan2_cc(csn, ccs)
    type (xplex), intent(in) :: csn, ccs
    atan2_cc%r = atan2(csn%r,ccs%r)
    atan2_cc%i = 0d0
    return
  end function atan2_cc

! COSH
  type (xplex) function cosh_c(z)
    type (xplex), intent(in) :: z
    cosh_c%r=cosh(z%r)
    cosh_c%i=0d0!(exp(z)+exp(-z))/2d0
    return
  end function cosh_c

! SINH
  type (xplex) function sinh_c(z)
    type (xplex), intent(in) :: z
    sinh_c%r=sinh(z%r)
    sinh_c%i=0d0!(exp(z)-exp(-z))/2d0
    return
  end function sinh_c

! TAN
  type (xplex) function tan_c(z)
    type (xplex), intent(in) :: z
    tan_c%r=tan(z%r)
    tan_c%i=0d0!z%i/cos(z%r)**2)
    return
  end function tan_c
  
! TANH
  type (xplex) function tanh_c(a)
    type (xplex), intent(in) :: a
    tanh_c%r=tanh(a%r)
    tanh_c%i=0d0
    return
  end function tanh_c

! MAX, intrinsic
  type (xplex) function max_cc(val1, val2)
    type (xplex), intent(in) :: val1, val2
    if (val1%r > val2%r) then
      max_cc = val1
    else
      max_cc = val2
    endif
    return
  end function max_cc
  type (xplex) function max_cr(val1, val2)
    type (xplex), intent(in) :: val1    
    real*8, intent(in) :: val2    
    if (val1%r > val2) then
      max_cr = val1
    else
      max_cr = xplex(val2, 0.d0)
    endif
    return
  end function max_cr
  type (xplex) function max_rc(val1, val2)
    real*8, intent(in) :: val1
    type (XPLEX), intent(in) :: val2
    if (val1 > val2%r) then
      max_rc = xplex(val1, 0.d0)
    else
      max_rc = val2
    endif
    return
  end function max_rc
  type (xplex) function max_ccc(val1, val2, val3)
    type (xplex), intent(in) :: val1, val2, val3
    if (val1%r > val2%r) then
      max_ccc = val1
    else
      max_ccc = val2
    endif
    if (val3%r > max_ccc%r) then
      max_ccc = val3
    endif
    return
  end function max_ccc
  function max_cccc(val1, val2, val3, val4)
    type (xplex), intent(in) :: val1, val2, val3, val4
    type (xplex) max_cccc
    type (xplex) max_cccc2
    if (val1%r > val2%r) then
      max_cccc = val1
    else
      max_cccc = val2
    endif
    if (val3%r > val4%r) then
      max_cccc2 = val3
    else
      max_cccc2 = val4
    endif
    if (max_cccc2%r > max_cccc%r) then
      max_cccc = max_cccc2
    endif
    return
  end function max_cccc

! MIN, intrinsic
  type (xplex) function min_cc(val1, val2)
    type (xplex), intent(in) :: val1, val2
    if (val1%r < val2%r) then
      min_cc = val1
    else
      min_cc = val2
    endif
    return
  end function min_cc
  type (xplex) function min_cr(val1, val2)
    type (xplex), intent(in) :: val1    
    real*8, intent(in) :: val2    
    if (val1%r < val2) then
      min_cr = val1
    else
      min_cr = xplex(val2, 0.d0)
    endif
    return
  end function min_cr
  type (XPLEX) function min_rc(val1, val2)
    real*8, intent(in) :: val1
    type (XPLEX), intent(in) :: val2
    if (val1 < val2%r) then
      min_rc = xplex(val1, 0.d0)
    else
      min_rc = val2
    endif
    return
  end function min_rc
  type (xplex) function min_ccc(val1, val2, val3)
    type (xplex), intent(in) :: val1, val2, val3
    if (val1%r < val2%r) then
      min_ccc = val1
    else
      min_ccc = val2
    endif
    if (val3%r < min_ccc%r) then
      min_ccc = val3
    endif
    return
  end function min_ccc
  function min_cccc(val1, val2, val3, val4)
    type (xplex), intent(in) :: val1, val2, val3, val4
    type (xplex) min_cccc
    type (xplex) min_cccc2
    if (val1%r < val2%r) then
      min_cccc = val1
    else
      min_cccc = val2
    endif
    if (val3%r < val4%r) then
      min_cccc2 = val3
    else
      min_cccc2 = val4
    endif
    if (min_cccc2%r < min_cccc%r) then
      min_cccc = min_cccc2
    endif
    return
  end function min_cccc

  
! SIGN, intrinsic, assume that val1 is always a complex*16
!                  in reality could be int
  type (xplex) function sign_cc(val1, val2)
    type (xplex), intent(in) :: val1, val2
    real*8  sign
    if (val2%r < 0.d0) then
      sign = -1.d0
    else
      sign = 1.d0
    endif
    sign_cc = sign * val1
    return
  end function sign_cc
  type (xplex) function sign_cr(val1, val2)
    type (xplex), intent(in) :: val1
    real*8, intent(in) :: val2
    real*8 sign
    if (val2 < 0.d0) then
      sign = -1.d0
    else
      sign = 1.d0
    endif
    sign_cr = sign * val1
    return
  end function sign_cr
  type (xplex) function sign_rc(val1, val2)
    real*8, intent(in) :: val1
    type (xplex), intent(in) :: val2
    real*8 sign
    if (val2 < 0.d0) then
      sign = -1.d0
    else
      sign = 1.d0
    endif
    sign_rc = sign * val1
    return
  end function sign_rc

! DIM, intrinsic
  type (xplex) function dim_cc(val1, val2)
    type (xplex), intent(in) :: val1, val2
    if (val1%r > val2%r) then
      dim_cc = val1 - val2
    else
      dim_cc = xplex(0.d0, 0.d0)
    endif
    return
  end function dim_cc
  type (xplex) function dim_cr(val1, val2)
    type (xplex), intent(in) :: val1
    real*8, intent(in) :: val2
    if (val1%r > val2) then
      dim_cr = val1 - xplex(val2, 0.d0)
    else
      dim_cr = xplex(0.d0, 0.d0)
    endif
    return
  end function dim_cr
  type (xplex) function dim_rc(val1, val2)
    real*8, intent(in) :: val1
    type (xplex), intent(in) :: val2
    if (val1 > val2%r) then
      dim_rc = xplex(val1, 0.d0) - val2
    else
      dim_rc = xplex(0.d0, 0.d0)
    endif
    return
  end function dim_rc
  
! LOG10
  type (xplex) function log10_c(z)
    type (xplex), intent(in) :: z
    log10_c%r=log(z%r)/log(10.d0)
    log10_c%i=0d0
  end function log10_c

! NINT
  integer function nint_c(z)
    type (xplex), intent(in) :: z
    nint_c=nint(z%r)
  end function nint_c

! EPSILON !! bad news ulness compiled with -r8
  type (xplex) function epsilon_c(z)
    type (xplex), intent(in) :: z
    epsilon_c%r=epsilon(z%r)
    epsilon_c%i=0d0
  end function epsilon_c

! <, .lt.
  logical function lt_cc(lhs, rhs)
    type (xplex), intent(in) :: lhs, rhs
    lt_cc = lhs%r < rhs%r
  end function lt_cc
  logical function lt_cr(lhs, rhs)
    type (xplex), intent(in) :: lhs
    real*8, intent(in) :: rhs
    lt_cr = lhs%r < rhs
  end function lt_cr
  logical function lt_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    lt_rc = lhs < rhs%r
  end function lt_rc
  logical function lt_ci(lhs, rhs)
    type (xplex), intent(in) :: lhs
    integer, intent(in) :: rhs
    lt_ci = lhs%r < rhs
  end function lt_ci
  logical function lt_ic(lhs, rhs)
    integer, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    lt_ic = lhs < rhs%r
  end function lt_ic

! <=, .le.
  logical function le_cc(lhs, rhs)
    type (xplex), intent(in) :: lhs, rhs
    le_cc = lhs%r <= rhs%r
  end function le_cc
  logical function le_cr(lhs, rhs)
    type (xplex), intent(in) :: lhs
    real*8, intent(in) :: rhs
    le_cr = lhs%r <= rhs
  end function le_cr
  logical function le_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    le_rc = lhs <= rhs%r
  end function le_rc
  logical function le_ci(lhs, rhs)
    type (xplex), intent(in) :: lhs
    integer, intent(in) :: rhs
    le_ci = lhs%r <= rhs
  end function le_ci
  logical function le_ic(lhs, rhs)
    integer, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    le_ic = lhs <= rhs%r
  end function le_ic
  logical function le_xc(lhs, rhs)
    type (xplex), intent(in) :: lhs
    complex*16, intent(in) :: rhs
    le_xc = lhs%r <= dble(rhs)
  end function le_xc
  logical function le_cx(lhs, rhs)
    complex*16, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    le_cx = dble(lhs) <= rhs%r
  end function le_cx

! >, .gt.
  logical function gt_cc(lhs, rhs)
    type (xplex), intent(in) :: lhs, rhs
    gt_cc = lhs%r > rhs%r
  end function gt_cc
  logical function gt_cr(lhs, rhs)
    type (xplex), intent(in) :: lhs
    real*8, intent(in) :: rhs
    gt_cr = lhs%r > rhs
  end function gt_cr
  logical function gt_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    gt_rc = lhs > rhs%r
  end function gt_rc
  logical function gt_ci(lhs, rhs)
    type (xplex), intent(in) :: lhs
    integer, intent(in) :: rhs
    gt_ci = lhs%r > rhs
  end function gt_ci
  logical function gt_ic(lhs, rhs)
    integer, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    gt_ic = lhs > rhs%r
  end function gt_ic

!! here are the redefined ones:
! >=, .ge.
  logical function ge_cc(lhs, rhs)
    type (xplex), intent(in) :: lhs, rhs
    ge_cc = lhs%r >= rhs%r
  end function ge_cc
  logical function ge_rr(lhs, rhs)
    real*8, intent(in) :: lhs, rhs
    ge_rr = lhs >= rhs
  end function ge_rr
  logical function ge_ii(lhs, rhs)
    integer, intent(in) :: lhs, rhs
    ge_ii = lhs >= rhs
  end function ge_ii
  logical function ge_aa(lhs, rhs)
    character(len=*), intent(in) :: lhs, rhs
    ge_aa = lhs >= rhs
  end function ge_aa
  logical function ge_cr(lhs, rhs)
    type (xplex), intent(in) :: lhs
    real*8, intent(in) :: rhs
    ge_cr = lhs%r >= rhs
  end function ge_cr
  logical function ge_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    ge_rc = lhs >= rhs%r
  end function ge_rc
  logical function ge_ci(lhs, rhs)
    type (xplex), intent(in) :: lhs
    integer, intent(in) :: rhs
    ge_ci = lhs%r >= rhs
  end function ge_ci
  logical function ge_ic(lhs, rhs)
    integer, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    ge_ic = lhs >= rhs%r
  end function ge_ic
  logical function ge_ir(lhs, rhs)
    integer, intent(in) :: lhs
    real*8, intent(in) :: rhs
    ge_ir = lhs >= rhs
  end function ge_ir
  logical function ge_ri(lhs, rhs)
    real*8, intent(in) :: lhs
    integer, intent(in) :: rhs
    ge_ri = lhs >= rhs
  end function ge_ri

! ==, .eq.
  logical function eq_cc(lhs, rhs)
    type (xplex), intent(in) :: lhs, rhs
    eq_cc = lhs%r == rhs%r
  end function eq_cc
  logical function eq_rr(lhs, rhs)
    real*8, intent(in) :: lhs, rhs
    eq_rr = lhs == rhs
  end function eq_rr
  logical function eq_ii(lhs, rhs)
    integer, intent(in) :: lhs, rhs
    eq_ii = lhs == rhs
  end function eq_ii
  logical function eq_aa(lhs, rhs)
    character(len=*), intent(in) :: lhs, rhs
    eq_aa = lhs == rhs
  end function eq_aa
  logical function eq_cr(lhs, rhs)
    type (xplex), intent(in) :: lhs
    real*8, intent(in) :: rhs
    eq_cr = lhs%r == rhs
  end function eq_cr
  logical function eq_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    eq_rc = lhs == rhs%r
  end function eq_rc
  logical function eq_ci(lhs, rhs)
    type (xplex), intent(in) :: lhs
    integer, intent(in) :: rhs
    eq_ci = lhs%r == rhs
  end function eq_ci
  logical function eq_ic(lhs, rhs)
    integer, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    eq_ic = lhs == rhs%r
  end function eq_ic
  logical function eq_ir(lhs, rhs)
    integer, intent(in) :: lhs
    real*8, intent(in) :: rhs
    eq_ir = lhs == rhs
  end function eq_ir
  logical function eq_ri(lhs, rhs)
    real*8, intent(in) :: lhs
    integer, intent(in) :: rhs
    eq_ri = lhs == rhs
  end function eq_ri

! /=, .ne.
  logical function ne_cc(lhs, rhs)
    type (xplex), intent(in) :: lhs, rhs
    ne_cc = lhs%r /= rhs%r
  end function ne_cc
  logical function ne_rr(lhs, rhs)
    real*8, intent(in) :: lhs, rhs
    ne_rr = lhs /= rhs
  end function ne_rr
  logical function ne_ii(lhs, rhs)
    integer, intent(in) :: lhs, rhs
    ne_ii = lhs /= rhs
  end function ne_ii
  logical function ne_aa(lhs, rhs)
    character(len=*), intent(in) :: lhs, rhs
    ne_aa = lhs /= rhs
  end function ne_aa
  logical function ne_cr(lhs, rhs)
    type (xplex), intent(in) :: lhs
    real*8, intent(in) :: rhs
    ne_cr = lhs%r /= rhs
  end function ne_cr
  logical function ne_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    ne_rc = lhs /= rhs%r
  end function ne_rc
  logical function ne_ci(lhs, rhs)
    type (xplex), intent(in) :: lhs
    integer, intent(in) :: rhs
    ne_ci = lhs%r /= rhs
  end function ne_ci
  logical function ne_ic(lhs, rhs)
    integer, intent(in) :: lhs
    type (xplex), intent(in) :: rhs
    ne_ic = lhs /= rhs%r
  end function ne_ic
  logical function ne_ir(lhs, rhs)
    integer, intent(in) :: lhs
    real*8, intent(in) :: rhs
    ne_ir = lhs /= rhs
  end function ne_ir
  logical function ne_ri(lhs, rhs)
    real*8, intent(in) :: lhs
    integer, intent(in) :: rhs
    ne_ri = lhs /= rhs
  end function ne_ri

! ISNAN
  logical elemental function isnan_c(var)
    type (xplex), intent(in) :: var
    isnan_c = isnan(var%r) .OR. isnan(var%i)
  end function isnan_c
  

! FP_CLASS
!  integer function fp_class_c(var)
!    complex*16, intent(in) ::var
!    fp_class_c = fp_class(dble(var))
!  end function fp_class_c

! ABS
!  integer function exponent_c(var)
!    complex*16, intent(in) ::var
!    exponent_c = exponent(dble(var))
!  end function exponent_c

! MAXABS
!  integer function maxexponent_c(var)
!    complex*16, intent(in) ::var
!    maxexponent_c = maxexponent(dble(var))
!  end function maxexponent_c

! MINABS
!  integer function minexponent_c(var)
!    complex*16, intent(in) ::var
!    minexponent_c = minexponent(dble(var))
!  end function minexponent_c

! MINVAL
  type (xplex) function minval_c(array)
    type (xplex), intent(in) ::array(:)
    integer :: ii
    real*8 :: min_real, min_comp
    min_real = minval(array(:)%r)
    do ii=1,size(array,1)
       if (array(ii)%r.eq.min_real) then
           min_comp = array(ii)%i
       endif
    enddo
    minval_c = xplex(min_real,min_comp)
  end function minval_c
  type (xplex) function minval_2c(array)
    type (xplex), intent(in) ::array(:,:)
    integer :: ii,jj
    real*8 :: min_real, min_comp
    min_real = minval(array(:,:)%r)
    do ii=1,size(array,1)
       do jj=1,size(array,2)
          if (array(ii,jj)%r.eq.min_real) then
              min_comp = array(ii,jj)%i
          endif
       enddo
    enddo
    minval_2c = xplex(min_real,min_comp)
  end function minval_2c
  type (xplex) function minval_3c(array)
    type (xplex), intent(in) ::array(:,:,:)
    integer :: ii,jj,kk
    real*8 :: min_real, min_comp
    min_real = minval(array(:,:,:)%r)
    do ii=1,size(array,1)
       do jj=1,size(array,2)
          do kk=1,size(array,3)
             if (array(ii,jj,kk)%r.eq.min_real) then
                 min_comp = array(ii,jj,kk)%i
             endif
          enddo
       enddo
    enddo
    minval_3c = xplex(min_real,min_comp)
  end function minval_3c
  type (xplex) function minval_4c(array)
    type (xplex), intent(in) ::array(:,:,:,:)
    integer :: ii,jj,kk,ll
    real*8 :: min_real, min_comp
    min_real = minval(array(:,:,:,:)%r)
    do ii=1,size(array,1)
       do jj=1,size(array,2)
          do kk=1,size(array,3)
             do ll=1,size(array,4)
       if (array(ii,jj,kk,ll)%r.eq.min_real) then
           min_comp = array(ii,jj,kk,ll)%i
       endif
             enddo
          enddo
       enddo
    enddo
    minval_4c = xplex(min_real,min_comp)
  end function minval_4c

! MAXVAL
  type (xplex) function maxval_c(array)
    type (xplex), intent(in) ::array(:)
    integer :: ii
    real*8 :: max_real, max_comp
    max_real = maxval(array(:)%r)
    do ii=1,size(array)
       if (array(ii)%r.eq.max_real) then
           max_comp = array(ii)%i
       endif
    enddo
    maxval_c = xplex(max_real,max_comp)
  end function maxval_c
  type (xplex) function maxval_2c(array)
    type (xplex), intent(in) ::array(:,:)
    integer :: ii,jj
    real*8 :: max_real, max_comp
    max_real = maxval(array(:,:)%r)
    do ii=1,size(array,1)
       do jj=1,size(array,2)
       if (array(ii,jj)%r.eq.max_real) then
           max_comp = array(ii,jj)%i
       endif
       enddo
    enddo
    maxval_2c = xplex(max_real,max_comp)
  end function maxval_2c
  type (xplex) function maxval_3c(array)
    type (xplex), intent(in) ::array(:,:,:)
    integer :: ii,jj,kk
    real*8 :: max_real, max_comp
    max_real = maxval(array(:,:,:)%r)
    do ii=1,size(array,1)
       do jj=1,size(array,2)
          do kk=1,size(array,3)
       if (array(ii,jj,kk)%r.eq.max_real) then
           max_comp = array(ii,jj,kk)%i
       endif
          enddo
       enddo
    enddo
    maxval_3c = xplex(max_real,max_comp)
  end function maxval_3c
  type (xplex) function maxval_4c(array)
    type (xplex), intent(in) ::array(:,:,:,:)
    integer :: ii,jj,kk,ll
    real*8 :: max_real, max_comp
    max_real = maxval(array(:,:,:,:)%r)
    do ii=1,size(array,1)
       do jj=1,size(array,2)
          do kk=1,size(array,3)
             do ll=1,size(array,4)
       if (array(ii,jj,kk,ll)%r.eq.max_real) then
           max_comp = array(ii,jj,kk,ll)%i
       endif
             enddo
          enddo
       enddo
    enddo
    maxval_4c = xplex(max_real,max_comp)
  end function maxval_4c

! TINY   
  type (xplex) function tiny_c(var)
    type (xplex), intent(in):: var
    tiny_c%r = tiny(var%r)
    tiny_c%i = 0d0
  end function tiny_c

! MOD
  real*8 function mod_c(var1,var2)
    type (xplex), intent(in) :: var1
    real*8, intent(in) :: var2
    mod_c = mod(var1%r,var2)
  end function mod_c
  real*8 function mod_c2(var1,var2)
    type (xplex), intent(in) :: var1
    type (xplex), intent(in) :: var2
    mod_c2 = mod(var1%r,var2%r)
  end function mod_c2
! XPLX
      type (xplex) function xplx_i(i)
        integer,intent(in)::i
        xplx_i = xplex(dble(i),0d0)
      end function xplx_i
      function xplx_i2(i) result(C)
        integer,intent(in)::i(:,:)
        type (xplex)::C(size(i,1),size(i,2))
        C%r = dble(i)
        C%i = 0d0
      end function xplx_i2
      function xplx_r2r2(r1,r2) result(C)
        real*8,intent(in)::r1(:,:),r2(:,:)
        type (xplex)::C(size(r1,1),size(r1,2))
        C%r = r1
        C%i = r2
      end function xplx_r2r2
      function xplx_r4r4(r1,r2) result(C)
        real*8,intent(in)::r1(:,:,:,:),r2(:,:,:,:)
        type (xplex)::C(size(r1,1),size(r1,2),size(r1,3),size(r1,4))
        C%r = r1
        C%i = r2
      end function xplx_r4r4
      function xplx_r1r1(r1,r2) result(C)
        real*8,intent(in)::r1(:),r2(:)
        type (xplex)::C(size(r1,1))
        C%r = r1
        C%i = r2
      end function xplx_r1r1
      type (xplex) function xplx_r(r)
        real*8,intent(in)::r
        xplx_r = xplex(r,0d0)
      end function xplx_r
      type (xplex) function xplx_x(z)
        type (xplex),intent(in)::z
        xplx_x = xplex(z%r,z%i)
      end function xplx_x

end module complexify
