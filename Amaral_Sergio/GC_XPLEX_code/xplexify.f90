      module complexify
      
      use MYTYPE
      IMPLICIT NONE
   
      ! - - - - - - - -
      interface operator (-)
        module procedure minus_x
        module procedure minus_x3
        module procedure minus_xx
        module procedure minus_x2x2
        module procedure minus_x1x1
        module procedure minus_x1x
        module procedure minus_x1i
        module procedure minus_x2x
        module procedure minus_x3x
        module procedure minus_xr
        module procedure minus_rx
        module procedure minus_xi
        module procedure minus_x2i
        module procedure minus_ix
      end interface
      ! / / / / / / / / 
      interface operator (/)
        module procedure div_xx
        module procedure div_x2x
        module procedure div_x2r
        module procedure div_x3r
        module procedure div_x3i
        module procedure div_x2x2
        module procedure div_x1x1
        module procedure div_x1x
        module procedure div_x3x3
        module procedure div_x3x
        module procedure div_x4i
        module procedure div_xr
        module procedure div_x1r
        module procedure div_xi
        module procedure div_rx
        module procedure div_ix
      end interface
      ! MIN
      interface min
        module procedure min_rx
        module procedure min_xx
        module procedure min_xr
        module procedure min_xxx
      end interface
      ! MAX
      interface max
        module procedure max_rx
        module procedure max_xr
        module procedure max_xx
        module procedure max_xxx
      end interface
      ! ** ** ** ** ** **
      interface operator (**)
        module procedure pow_xi
        module procedure pow_x2i
        module procedure pow_rx
        module procedure pow_xx
        module procedure pow_xr
        module procedure pow_x1r
      end interface
      ! .LT.
      interface operator (.lt.)
        module procedure lt_xx
        module procedure lt_rx
        module procedure lt_xr
        module procedure lt_xi
      end interface
      ! .GT.
      interface operator (.gt.)
        module procedure gt_xr
        module procedure gt_xx
        module procedure gt_xi
        module procedure gt_ix
      end interface
      ! .GE.
      interface operator (.ge.)
        module procedure ge_xr
        module procedure ge_xx
        module procedure ge_xi
      end interface
      ! .LE.
      interface operator (.le.)
        module procedure le_xi
        module procedure le_xx
        module procedure le_xr
        module procedure le_rx
        module procedure le_ix
      end interface
      ! .NE.
      interface operator (.ne.)
        module procedure ne_xx
        module procedure ne_xr
        module procedure ne_xi
      end interface
      ! ==
      interface operator (==)
        module procedure eq_xi
        module procedure eq_xx
        module procedure eq_rx
        module procedure eq_xr
      end interface
      ! * * * * * * * * *
      interface operator (*)
        module procedure mult_xx
        module procedure mult_x1x1
        module procedure mult_x1x
        module procedure mult_x1r
        module procedure mult_xx1
        module procedure mult_xx2
        module procedure mult_xx3
        module procedure mult_x2x
        module procedure mult_x2x2
        module procedure mult_x3x3
        module procedure mult_x4x4
        module procedure mult_x3x
        module procedure mult_x2r
        module procedure mult_x3r
        module procedure mult_rx
        module procedure mult_rx2
        module procedure mult_rx1
        module procedure mult_xr
        module procedure mult_ix
        module procedure mult_ix2
        module procedure mult_xi
      end interface
      ! + + + + + + + + +
      interface operator (+)
        module procedure add_xx
        module procedure add_x1x1
        module procedure add_xx1
        module procedure add_x1x
        module procedure add_x1i
        module procedure add_x1r
        module procedure add_xx2
        module procedure add_x2x2
        module procedure add_x2x
        module procedure add_x2r
        module procedure add_x3x3
        module procedure add_x4x4
        module procedure add_x3x
        module procedure add_xr
        module procedure add_xi
        module procedure add_x2i2
        module procedure add_rx
        module procedure add_ix
      end interface
      ! ISNAN
      interface ISNAN
        module procedure isnan_x
      end interface
      ! = 
      interface assignment (=)
        module procedure asgn_xx
        module procedure asgn_x1x1
        module procedure asgn_x2x2
        module procedure asgn_x3x3
        module procedure asgn_x4x4
        module procedure asgn_xr
        module procedure asgn_x1r
        module procedure asgn_x2r
        module procedure asgn_x3r
        module procedure asgn_x4r
        module procedure asgn_x5r
        module procedure asgn_xr4
        module procedure asgn_xi
        module procedure asgn_x1i1
        module procedure asgn_x2i2
        module procedure asgn_x3i3
        module procedure asgn_x4i4
        module procedure asgn_x1i
        module procedure asgn_x2i
        module procedure asgn_x3i
        module procedure asgn_x4i
        module procedure asgn_ix
        module procedure asgn_i1x1
        module procedure asgn_i2x2
        module procedure asgn_i3x3
        module procedure asgn_i4x4
        module procedure asgn_x1r41
        module procedure asgn_x2r42
        module procedure asgn_x3r43
        module procedure asgn_x4r44
        module procedure asgn_x5r45
        module procedure asgn_x6r46
        module procedure asgn_x7r47
        module procedure asgn_x1r81
        module procedure asgn_x2r82
        module procedure asgn_x3r83
        module procedure asgn_x4r84
        module procedure asgn_x5r85
        module procedure asgn_x6r86
        module procedure asgn_x7r87
      end interface
      ! EXP
      interface EXP
        module procedure exp_x
        module procedure exp_x1
      end interface
      ! LOG
      interface LOG
        module procedure log_x
        module procedure log_x1
        module procedure log_x4
      end interface
      ! LOG10
      interface LOG10
        module procedure log10_x
      end interface
      ! ABS
      interface ABS
        module procedure abs_x
        module procedure abs_x4
      end interface
      ! INT
      interface INT
        module procedure int_x
      end interface
      ! NINT
      interface NINT
        module procedure nint_x
      end interface
      ! SIN
      interface SIN
        module procedure sin_x
      end interface
      ! COS
      interface COS
        module procedure cos_x
      end interface
      ! SUM
      interface SUM
        module procedure sum_x
        module procedure sum_xx
        module procedure sum_xxx
        module procedure sum_xxxx
      end interface
      ! MINVAL
      interface MINVAL
        module procedure minval_x
        module procedure minval_xx
        module procedure minval_xxx
        module procedure minval_xxxx
      end interface  
      ! MAXVAL
      interface MAXVAL
        module procedure maxval_xx
        module procedure maxval_xxxx
        module procedure maxval_x
        module procedure maxval_xxx
      end interface
      ! TINY
      interface TINY
        module procedure tiny_x
      end interface
      ! SQRT
      interface SQRT
        module procedure sqrt_x
      end interface
      ! MOD
      interface MOD
        module procedure mod_xr
        module procedure mod_xx
      end interface
      ! SIGN
      interface SIGN
        module procedure sign_xx
        module procedure sign_rx
      end interface
      ! XPLX
      interface XPLX
        module procedure xplx_i
        module procedure xplx_i2
        module procedure xplx_r
        module procedure xplx_r2r2
        module procedure xplx_r1r1
        module procedure xplx_x
      end interface
      ! ACOS
      interface ACOS
        module procedure acos_x
      end interface
      ! CEILING
      interface CEILING
        module procedure ceiling_x
      end interface
      ! ATAN
      interface ATAN
        module procedure atan_x
      end interface
      ! FLOOR
      interface FLOOR
        module procedure floor_x
      end interface
      CONTAINS
! FLOOR
      integer function floor_x(z)
        type(xplex),intent(in)::z
        floor_x = floor(z%r)
      end function floor_x
! ATAN
      type (xplex) function atan_x(z)
        type(xplex),intent(in)::z
        atan_x%r = atan(z%r)
        atan_x%i = z%i/sqrt(1+z%r**2)
      end function atan_x
! CEILING
      integer function ceiling_x(z)
        type (xplex),intent(in)::z
        ceiling_x = ceiling(z%r)
      end function ceiling_x
! ACOS
      type (xplex) function acos_x(z)
        type (xplex), intent(in):: z
        acos_x%r = acos(z%r)
        acos_x%i = -z%i/sqrt(1-z%r**2)
      end function acos_x
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
! SIGN
      type (xplex) function sign_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        real*8::z1r,z1i
        z1r = z1%r
        z1i = z1%i
        sign_xx%r = sign(z1%r,z2%r)
        if (sign_xx%r==z1r) then
            sign_xx%i = z1i
        else
            sign_xx%i = -z1i
        endif
      end function sign_xx
      type (xplex) function sign_rx(r,z)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        sign_rx%r = sign(r,z%r)
        sign_rx%i = 0d0
      end function sign_rx
! MOD
      type (xplex) function mod_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        mod_xr%r = mod(z%r,r)
        mod_xr%i = 0d0
      end function mod_xr
      type (xplex) function mod_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        mod_xx%r = mod(z1%r,z2%r)
        mod_xx%i = 0d0
      end function mod_xx

! SQRT
      type (xplex) function sqrt_x(z)
        type (xplex),intent(in)::z
        sqrt_x%r = dble(sqrt(dcmplx(z%r,z%i)))
        sqrt_x%i = dimag(sqrt(dcmplx(z%r,z%i)))
        !if (z%r==0d0) then
        !   sqrt_x%i = 0d0
        !else
        !   sqrt_x%i = z%i/(2*sqrt(z%r))
        !endif
      end function sqrt_x

! TINY
      type (xplex) function tiny_x(z)
        type (xplex),intent(in)::z
        tiny_x%r = tiny(z%r)
        tiny_x%i = 0d0
      end function tiny_x

! COS
      type (xplex) function cos_x(z)
        type (xplex),intent(in)::z
        cos_x%r = dble(cos(dcmplx(z%r,z%i)))
        cos_x%i = dimag(cos(dcmplx(z%r,z%i)))
        !cos_x%r = cos(z%r)
        !cos_x%i =-sin(z%r)*z%i
      end function cos_x
! MINVAL
      type (xplex) function minval_x(z)
        type (xplex),intent(in),dimension(:)::z
        integer::i
        minval_x%r=minval(z%r)
        do i=1,size(z,1)
           if (minval_x%r==z(i)%r) minval_x%i = z(i)%i
        enddo
      end function minval_x
      type (xplex) function minval_xx(z)
        type (xplex),intent(in),dimension(:,:)::z
        integer::i,j
        minval_xx%r=minval(z%r)
        do j=1,size(z,2)
        do i=1,size(z,1)
           if (minval_xx%r==z(i,j)%r) minval_xx%i = z(i,j)%i
        enddo
        enddo
      end function minval_xx
      type (xplex) function minval_xxx(z)
        type (xplex),intent(in),dimension(:,:,:)::z
        integer::i,j,l
        minval_xxx%r=minval(z%r)
        do l=1,size(z,3)
        do j=1,size(z,2)
        do i=1,size(z,1)
           if (minval_xxx%r==z(i,j,l)%r) minval_xxx%i = z(i,j,l)%i
        enddo
        enddo
        enddo
      end function minval_xxx
      type (xplex) function minval_xxxx(z)
        type (xplex),intent(in),dimension(:,:,:,:)::z
        integer::i,j,k,l
        minval_xxxx%r=minval(z%r)
        do l=1,size(z,4)
        do k=1,size(z,3)
        do j=1,size(z,2)
        do i=1,size(z,1)
           if (minval_xxxx%r==z(i,j,k,l)%r) minval_xxxx%i = z(i,j,k,l)%i
        enddo
        enddo
        enddo
        enddo
      end function minval_xxxx
! MAXVAL
      type (xplex) function maxval_x(z)
        type (xplex),intent(in),dimension(:)::z
        integer::i
        maxval_x%r=maxval(z%r)
        do i=1,size(z,1)
           if (maxval_x%r==z(i)%r) maxval_x%i = z(i)%i
        enddo
      end function maxval_x
      type (xplex) function maxval_xx(z)
        type (xplex),intent(in),dimension(:,:)::z
        integer::i,j
        maxval_xx%r=maxval(z%r)
        do j=1,size(z,2)
        do i=1,size(z,1)
           if (maxval_xx%r==z(i,j)%r) maxval_xx%i = z(i,j)%i
        enddo
        enddo
      end function maxval_xx
      type (xplex) function maxval_xxx(z)
        type (xplex),intent(in),dimension(:,:,:)::z
        integer::i,j,l
        maxval_xxx%r=maxval(z%r)
        do l=1,size(z,3)
        do j=1,size(z,2)
        do i=1,size(z,1)
           if (maxval_xxx%r==z(i,j,l)%r) maxval_xxx%i = z(i,j,l)%i
        enddo
        enddo
        enddo
      end function maxval_xxx
      type (xplex) function maxval_xxxx(z)
        type (xplex),intent(in),dimension(:,:,:,:)::z
        integer::i,j,k,l
        maxval_xxxx%r=maxval(z%r)
        do l=1,size(z,4)
        do k=1,size(z,3)
        do j=1,size(z,2)
        do i=1,size(z,1)
           if (maxval_xxxx%r==z(i,j,k,l)%r) maxval_xxxx%i = z(i,j,k,l)%i
        enddo
        enddo
        enddo
        enddo
      end function maxval_xxxx
! SUM
      type (xplex) function sum_x(z)
        type (xplex), intent(in),dimension(:)::z    
        sum_x%r = sum(z%r)
        sum_x%i = sum(z%i)
      end function sum_x
      type (xplex) function sum_xx(z)
        type (xplex), intent(in),dimension(:,:)::z
        sum_xx%r = sum(z%r)
        sum_xx%i = sum(z%i)
      end function sum_xx
      type (xplex) function sum_xxx(z)
        type (xplex), intent(in),dimension(:,:,:)::z
        sum_xxx%r = sum(z%r)
        sum_xxx%i = sum(z%i)
      end function sum_xxx
      type (xplex) function sum_xxxx(z)
        type (xplex), intent(in),dimension(:,:,:,:)::z
        sum_xxxx%r = sum(z%r)
        sum_xxxx%i = sum(z%i)
      end function sum_xxxx
! SIN
      type (xplex) function sin_x(z)
        type (xplex),intent(in)::z
        sin_x%r = dble(sin(dcmplx(z%r,z%i)))
        sin_x%i = dimag(sin(dcmplx(z%r,z%i)))
        !sin_x%r = sin(z%r)
        !sin_x%i = cos(z%r)*z%i
      end function sin_x
! NINT
      integer function nint_x(z)
        type (xplex),intent(in)::z
        nint_x = nint(z%r)
      end function nint_x
! INT 
      integer function int_x(z)
        type (xplex),intent(in)::z
        int_x = int(z%r)
      end function int_x
! ABS
      type (xplex) function abs_x(z)
        type (xplex),intent(in)::z
        abs_x%r = z%r
        abs_x%i = z%i
        if (abs_x%r<0d0) abs_x=-abs_x
      end function abs_x 
      function abs_x4(z) result(C)
        type (xplex),intent(in)::z(:,:,:,:)
        type (xplex)::C(size(z,1),size(z,2),size(z,3),size(z,4))
        C%r = z%r
        C%i = z%i
        where (C%r<0d0) 
           C%r=-C%r
           C%i=-C%i
        endwhere
      end function abs_x4
! LOG
      type (xplex) function log_x(z)
        type (xplex), intent(in)::z
        log_x%r = dble(log(dcmplx(z%r,z%i)))
        log_x%i = dimag(log(dcmplx(z%r,z%i)))
        !log_x%r = log(z%r)
        !log_x%i = 0d0!z%i/z%r
      end function log_x
      function log_x4(z) result(C)
        type (xplex), intent(in)::z(:,:,:,:)
        type (xplex)::C(size(z,1),size(z,2),size(z,3),size(z,4))
        C%r = dble(log(dcmplx(z%r,z%i)))
        C%i = dimag(log(dcmplx(z%r,z%i)))
        !C%r = log(z%r)
        !C%i = z%i/z%r
      end function log_x4
      function log_x1(z) result(C)
        type (xplex), intent(in)::z(:)
        type (xplex)::C(size(z,1))
        C%r = dble(log(dcmplx(z%r,z%i)))
        C%i = dimag(log(dcmplx(z%r,z%i)))
        !C%r = log(z%r)
        !C%i = z%i/z%r
      end function log_x1
! LOG10
      type (xplex) function log10_x(z)
        type (xplex), intent(in)::z
        log10_x%r = dble( log(dcmplx(z%r,z%i))/log(10.0d0))
        log10_x%i = dimag( log(dcmplx(z%r,z%i))/log(10.0d0))
        !log10_x%r = log(z%r)/log(10.0d0)
        !log10_x%i = (z%i/z%r)/log(10.0d0)
      end function log10_x 
! EXP
      type (xplex) function exp_x(z)
        type (xplex),intent(in)::z
        exp_x%r = dble(exp(dcmplx(z%r,z%i)))
        exp_x%i = dimag(exp(dcmplx(z%r,z%i)))
        !exp_x%r = exp(z%r)
        !exp_x%i = exp(z%r)*z%i
      end function exp_x
      function exp_x1(z) result(C)
        type (xplex),intent(in)::z(:)
        type (xplex)::C(size(z,1))
        C%r = dble(exp(dcmplx(z%r,z%i)))
        C%i = dble(exp(dcmplx(z%r,z%i)))
        !C%r = exp(z%r)
        !C%i = exp(z%r)*z%i
      end function exp_x1

! = 
      subroutine asgn_xx(z1,z2)
        type(xplex),intent(inout)::z1
        type(xplex),intent(in)::z2
        z1%r = z2%r
        z1%i = z2%i
      end subroutine asgn_xx
      subroutine asgn_x1x1(z1,z2)
        type(xplex),intent(inout)::z1(:)
        type(xplex),intent(in)::z2(:)
        z1%r = z2%r
        z1%i = z2%i
      end subroutine asgn_x1x1
      subroutine asgn_x2x2(z1,z2)
        type(xplex),intent(inout)::z1(:,:)
        type(xplex),intent(in)::z2(:,:)
        z1%r = z2%r
        z1%i = z2%i
      end subroutine asgn_x2x2
      subroutine asgn_x3x3(z1,z2)
        type(xplex),intent(inout)::z1(:,:,:)
        type(xplex),intent(in)::z2(:,:,:)
        z1%r = z2%r
        z1%i = z2%i
      end subroutine asgn_x3x3
      subroutine asgn_x4x4(z1,z2)
        type(xplex),intent(inout)::z1(:,:,:,:)
        type(xplex),intent(in)::z2(:,:,:,:)
        z1%r = z2%r
        z1%i = z2%i
      end subroutine asgn_x4x4


      subroutine asgn_xr(z,r)
        type(xplex),intent(inout)::z 
        real*8,intent(in)::r
        z%r = r
        z%i = 0d0
      end subroutine asgn_xr
      subroutine asgn_x1r(z,r)
        type(xplex),intent(inout)::z(:)
        real*8,intent(in)::r
        z%r = r
        z%i = 0d0
      end subroutine asgn_x1r
      subroutine asgn_x2r(z,r)
        type(xplex),intent(inout)::z(:,:)
        real*8,intent(in)::r
        z%r = r
        z%i = 0d0
      end subroutine asgn_x2r
      subroutine asgn_x3r(z,r)
        type(xplex),intent(inout)::z(:,:,:)
        real*8,intent(in)::r
        z%r = r
        z%i = 0d0
      end subroutine asgn_x3r
      subroutine asgn_x4r(z,r)
        type(xplex),intent(inout)::z(:,:,:,:)
        real*8,intent(in)::r
        z%r = r
        z%i = 0d0
      end subroutine asgn_x4r
      subroutine asgn_x5r(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:)
        real*8,intent(in)::r
        z%r = r
        z%i = 0d0
      end subroutine asgn_x5r
      subroutine asgn_xr4(z,r4)
        type(xplex),intent(inout)::z
        real*4,intent(in)::r4
        z%r = dble(r4)
        z%i = 0d0
      end subroutine asgn_xr4
      subroutine asgn_xi(z,i)
        type(xplex),intent(inout)::z
        integer,intent(in)::i
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_xi
      subroutine asgn_x1i1(z,i)
        type(xplex),intent(inout)::z(:)
        integer,intent(in)::i(:)
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x1i1
      subroutine asgn_x2i2(z,i)
        type(xplex),intent(inout)::z(:,:)
        integer,intent(in)::i(:,:)
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x2i2
      subroutine asgn_x3i3(z,i)
        type(xplex),intent(inout)::z(:,:,:)
        integer,intent(in)::i(:,:,:)
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x3i3
      subroutine asgn_x4i4(z,i)
        type(xplex),intent(inout)::z(:,:,:,:)
        integer,intent(in)::i(:,:,:,:)
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x4i4
      subroutine asgn_x1i(z,i)
        type(xplex),intent(inout)::z(:)
        integer,intent(in)::i
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x1i
      subroutine asgn_x2i(z,i)
        type(xplex),intent(inout)::z(:,:)
        integer,intent(in)::i
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x2i
      subroutine asgn_x3i(z,i)
        type(xplex),intent(inout)::z(:,:,:)
        integer,intent(in)::i
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x3i
      subroutine asgn_x4i(z,i)
        type(xplex),intent(inout)::z(:,:,:,:)
        integer,intent(in)::i
        z%r = dble(i)
        z%i = 0d0
      end subroutine asgn_x4i
      subroutine asgn_ix(i,z)
        integer,intent(inout)::i
        type(xplex),intent(in)::z
        i = int(z%r)
      end subroutine asgn_ix
      subroutine asgn_i1x1(i,z)
        integer,intent(inout)::i(:)
        type(xplex),intent(in)::z(:)
        i = int(z%r)
      end subroutine asgn_i1x1
      subroutine asgn_i2x2(i,z)
        integer,intent(inout)::i(:,:)
        type(xplex),intent(in)::z(:,:)
        i = int(z%r)
      end subroutine asgn_i2x2
      subroutine asgn_i3x3(i,z)
        integer,intent(inout)::i(:,:,:)
        type(xplex),intent(in)::z(:,:,:)
        i = int(z%r)
      end subroutine asgn_i3x3
      subroutine asgn_i4x4(i,z)
        integer,intent(inout)::i(:,:,:,:)
        type(xplex),intent(in)::z(:,:,:,:)
        i = int(z%r)
      end subroutine asgn_i4x4
      subroutine asgn_x1r81(z,r)
        type(xplex),intent(inout)::z(:)
        real*8,intent(in)::r(:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x1r81
      subroutine asgn_x2r82(z,r)
        type(xplex),intent(inout)::z(:,:)
        real*8,intent(in)::r(:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x2r82
      subroutine asgn_x3r83(z,r)
        type(xplex),intent(inout)::z(:,:,:)
        real*8,intent(in)::r(:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x3r83
      subroutine asgn_x4r84(z,r)
        type(xplex),intent(inout)::z(:,:,:,:)
        real*8,intent(in)::r(:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x4r84
      subroutine asgn_x5r85(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:)
        real*8,intent(in)::r(:,:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x5r85
      subroutine asgn_x6r86(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:,:)
        real*8,intent(in)::r(:,:,:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x6r86
      subroutine asgn_x7r87(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:,:,:)
        real*8,intent(in)::r(:,:,:,:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x7r87
      subroutine asgn_x1r41(z,r)
        type(xplex),intent(inout)::z(:)
        real*4,intent(in)::r(:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x1r41
      subroutine asgn_x2r42(z,r)
        type(xplex),intent(inout)::z(:,:)
        real*4,intent(in)::r(:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x2r42
      subroutine asgn_x3r43(z,r)
        type(xplex),intent(inout)::z(:,:,:)
        real*4,intent(in)::r(:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x3r43
      subroutine asgn_x4r44(z,r)
        type(xplex),intent(inout)::z(:,:,:,:)
        real*4,intent(in)::r(:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x4r44
      subroutine asgn_x5r45(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:)
        real*4,intent(in)::r(:,:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x5r45
      subroutine asgn_x6r46(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:,:)
        real*4,intent(in)::r(:,:,:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x6r46
      subroutine asgn_x7r47(z,r)
        type(xplex),intent(inout)::z(:,:,:,:,:,:,:)
        real*4,intent(in)::r(:,:,:,:,:,:,:)
        z%r = r
        z%i = 0d0
      end subroutine asgn_x7r47
! ==
      logical function eq_xi(z,i)
        type(xplex),intent(in)::z
        integer,intent(in)::i
        eq_xi = z%r==i
      end function eq_xi
      logical function eq_xx(z1,z2)
        type(xplex),intent(in)::z1,z2
        eq_xx = z1%r==z2%r
      end function eq_xx
      logical function eq_rx(r,z)
        type(xplex),intent(in)::z
        real*8,intent(in)::r
        eq_rx = r==z%r
      end function eq_rx
      logical function eq_xr(z,r)
        type(xplex),intent(in)::z
        real*8,intent(in)::r
        eq_xr = z%r==r
      end function eq_xr

! ISNAN
      logical function isnan_x(z)
        type (xplex),intent(in)::z
        isnan_x = isnan(z%r).or.isnan(z%i)
      end function isnan_x

! - - - - - - - - - - - - - - - - - - - - -
      type (xplex) function minus_x(z)
        type (xplex), intent(in)::z
        minus_x%r = -z%r
        minus_x%i = -z%i
      end function minus_x
      function minus_x3(z) result(C)
        type (xplex), intent(in)::z(:,:,:)
        type(xplex)::C(size(z,1),size(z,2),size(z,3))
        C%r = -z%r
        C%i = -z%i
      end function minus_x3
      type (xplex) function minus_xx(z1,z2)
        type (xplex), intent(in)::z1,z2
        minus_xx%r = z1%r-z2%r
        minus_xx%i = z1%i-z2%i
      end function minus_xx
      function minus_x2x2(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:),z2(:,:)
        type (xplex)::C(size(z1,1),size(z2,2))
        C%r = z1%r-z2%r
        C%i = z1%i-z2%i
      end function minus_x2x2
      function minus_x1x1(z1,z2) result(C)
        type (xplex), intent(in)::z1(:),z2(:)
        type (xplex)::C(size(z1,1))
        C%r = z1%r-z2%r
        C%i = z1%i-z2%i
      end function minus_x1x1
      function minus_x1x(z1,z2) result(C)
        type (xplex), intent(in)::z1(:),z2
        type (xplex)::C(size(z1,1))
        C%r = z1%r-z2%r
        C%i = z1%i-z2%i
      end function minus_x1x
      function minus_x1i(z1,i) result(C)
        type (xplex), intent(in)::z1(:)
        integer,intent(in)::i
        type (xplex)::C(size(z1,1))
        C%r = z1%r-dble(i)
        C%i = z1%i
      end function minus_x1i
      function minus_x2x(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:),z2
        type (xplex)::C(size(z1,1),size(z1,2))
        C%r = z1%r-z2%r
        C%i = z1%i-z2%i
      end function minus_x2x
      function minus_x3x(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:,:),z2
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
        C%r = z1%r-z2%r
        C%i = z1%i-z2%i
      end function minus_x3x
      type (xplex) function minus_xr(z,r)
        type (xplex), intent(in)::z
        real*8, intent(in)::r
        minus_xr%r = z%r-r
        minus_xr%i = z%i
      end function minus_xr
      type (xplex) function minus_rx(r,z)
        type (xplex),intent(in)::z
        real*8 , intent(in)::r
        minus_rx%r = r-z%r
        minus_rx%i =  -z%i
      end function minus_rx
      type (xplex) function minus_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        minus_xi%r = z%r - dble(i)
        minus_xi%i = z%i
      end function minus_xi
      function minus_x2i(z,i) result(C)
        type (xplex),intent(in)::z(:,:)
        type (xplex)::C(size(z,1),size(z,2))
        integer,intent(in)::i
        C(:,:)%r = z(:,:)%r - dble(i)
        C(:,:)%i = z(:,:)%i
      end function minus_x2i
      type (xplex) function minus_ix(i,z)
        type (xplex),intent(in)::z
        integer , intent(in)::i
        minus_ix%r = dble(i)-z%r
        minus_ix%i =  -z%i
      end function minus_ix
! / / / / / / / / / / / / / / / / / / / / / 
      type (xplex) function div_xx(z1,z2)
        type (xplex), intent(in)::z1,z2
        div_xx%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        div_xx%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !div_xx%r = z1%r/z2%r
        !div_xx%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        if (z1%r==1D-99.or.z1%r==1D-30.or.z1%r==1D-32.or.z1%r==1D-35    &
            .or. z1%r==1D-20) div_xx%i = 0d0 
      end function div_xx
      function div_x2x(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:),z2
        type (xplex)::C(size(z1,1),size(z1,2))
        C%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !C%r = z1%r/z2%r
        !C%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        where (z1%r==1D-99) C%i = 0d0
      end function div_x2x
      function div_x2r(z1,r) result(C)
        type (xplex), intent(in)::z1(:,:)
        real*8,intent(in)::r
        type (xplex)::C(size(z1,1),size(z1,2))
        C%r = dble(dcmplx(z1%r,z1%i)/r)
        C%i = dimag(dcmplx(z1%r,z1%i)/r)
        !C%r = z1%r/r
        !C%i = z1%i/r
      end function div_x2r
      function div_x3x(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:,:),z2
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
        C%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !C%r = z1%r/z2%r
        !C%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        where (z1%r==1D-99) C%i = 0d0
      end function div_x3x
      function div_x4i(z1,i) result(C)
        type (xplex), intent(in)::z1(:,:,:,:)
        integer,intent(in)::i
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
        C%r = dble(dcmplx(z1%r,z1%i)/i)
        C%i = dimag(dcmplx(z1%r,z1%i)/i)
        !C%r = z1%r/i
        !C%i = z1%i/i
      end function div_x4i
      function div_x2x2(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:),z2(:,:)
        type (xplex)::C(size(z1,1),size(z1,2))
        C%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !C%r = z1%r/z2%r
        !C%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        where (z1%r==1D-99) C%i = 0d0
      end function div_x2x2
      function div_x1x(z1,z2) result(C)
        type (xplex), intent(in)::z1(:),z2
        type (xplex)::C(size(z1,1))
        C%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !C%r = z1%r/z2%r
        !C%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        where (z1%r==1D-99) C%i = 0d0
      end function div_x1x
      function div_x1x1(z1,z2) result(C)
        type (xplex), intent(in)::z1(:),z2(:)
        type (xplex)::C(size(z1,1))
        C%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !C%r = z1%r/z2%r
        !C%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        where (z1%r==1D-99) C%i = 0d0
      end function div_x1x1
      function div_x3x3(z1,z2) result(C)
        type (xplex), intent(in)::z1(:,:,:),z2(:,:,:)
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
        C%r = dble(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)/dcmplx(z2%r,z2%i))
        !C%r = z1%r/z2%r
        !C%i = z1%i/z2%r-(z1%r*z2%i/z2%r)/z2%r!(z1%r/z2%r)*(z2%i/z2%r)
        where (z1%r==1D-99) C%i = 0d0
      end function div_x3x3
      function div_x3i(z1,i) result(C)
        type (xplex), intent(in)::z1(:,:,:)
        integer,intent(in)::i
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
        C%r = dble(dcmplx(z1%r,z1%i)/i)
        C%i = dimag(dcmplx(z1%r,z1%i)/i)
        !C%r = z1%r/i
        !C%i = z1%i/i
      end function div_x3i
      type (xplex) function div_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        div_xr%r = dble(dcmplx(z%r,z%i)/r)
        div_xr%i = dimag(dcmplx(z%r,z%i)/r)
        !div_xr%r = z%r/r
        !div_xr%i = z%i/r
      end function div_xr
      function div_x1r(z,r) result(C)
        type (xplex),intent(in)::z(:)
        real*8,intent(in)::r
        type (xplex)::C(size(z,1))
        C%r = dble(dcmplx(z%r,z%i)/r)
        C%i = dimag(dcmplx(z%r,z%i)/r)
        !C%r = z%r/r
        !C%i = z%i/r
      end function div_x1r
      function div_x3r(z,r) result(C)
        type (xplex),intent(in)::z(:,:,:)
        real*8,intent(in)::r
        type (xplex)::C(size(z,1),size(z,2),size(z,3))
        C%r = dble(dcmplx(z%r,z%i)/r)
        C%i = dimag(dcmplx(z%r,z%i)/r)
        !C%r = z%r/r
        !C%i = z%i/r
      end function div_x3r
      type (xplex) function div_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        div_xi%r = dble(dcmplx(z%r,z%i)/i)
        div_xi%i = dimag(dcmplx(z%r,z%i)/i)
        !div_xi%r = z%r/i
        !div_xi%i = z%i/i
      end function div_xi
      type (xplex) function div_rx(r,z)
        real*8, intent(in)::r
        type (xplex),intent(in)::z
        div_rx%r = dble(r/dcmplx(z%r,z%i))
        div_rx%i = dimag(r/dcmplx(z%r,z%i))
        !div_rx%r = r/z%r
        !div_rx%i = 0d0-(r*z%i/z%r)/z%r!(r/z%r)*(z%i/z%r)
        if (r==1D-99) div_rx%i = 0d0
      end function div_rx
      type (xplex) function div_ix(i,z)
        integer, intent(in)::i
        type (xplex),intent(in)::z
        div_ix%r = dble(i/dcmplx(z%r,z%i))
        div_ix%i = dimag(i/dcmplx(z%r,z%i))
        !div_ix%r = i/z%r
        !div_ix%i = 0d0-(i*z%i/z%r)/z%r!(i/z%r)*(z%i/z%r)
        if (i==1D-99) div_ix%i = 0d0
      end function div_ix
 
! MIN
      type (xplex) function min_rx(r,z)
        real*8,intent(in) :: r
        type (xplex),intent(in)::z
        min_rx%r = min(r,z%r)
        if (min_rx%r==z%r) then 
           min_rx%i = z%i
        else
           min_rx%i = 0d0
        endif
      end function min_rx
      type (xplex) function min_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        min_xx%r = min(z1%r,z2%r)
        if (min_xx%r==z1%r) then
           min_xx%i = z1%i
        else
           min_xx%i = z2%i
        endif
      end function min_xx
      type (xplex) function min_xxx(z1,z2,z3)
        type (xplex),intent(in)::z1,z2,z3
        type (xplex)::min_xx
        min_xx%r = min(z1%r,z2%r)
        if (min_xx%r==z1%r) then
           min_xx%i = z1%i
        else
           min_xx%i = z2%i
        endif
        min_xxx%r = min(min_xx%r,z3%r)
        if (min_xxx%r==min_xx%r) then
           min_xxx%i = min_xx%i
        else
           min_xxx%i = z3%i
        endif
      end function min_xxx
      type (xplex) function min_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        min_xr%r = min(z%r,r)
        if (min_xr%r==z%r) then
           min_xr%i = z%i
        else
           min_xr%i = 0d0
        endif
      end function min_xr

! MAX
      type (xplex) function max_rx(r,z)
        real*8,intent(in) :: r
        type (xplex),intent(in)::z
        max_rx%r = max(r,z%r)
        if (max_rx%r==z%r) then
           max_rx%i = z%i
        else
           max_rx%i = 0d0
        endif
      end function max_rx
      type (xplex) function max_xr(z,r)
        real*8,intent(in) :: r
        type (xplex),intent(in)::z
        max_xr%r = max(z%r,r)
        if (max_xr%r==z%r) then
           max_xr%i = z%i
        else
           max_xr%i = 0d0
        endif
      end function max_xr
      type (xplex) function max_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        max_xx%r = max(z1%r,z2%r)
        if (max_xx%r==z1%r) then
           max_xx%i = z1%i
        else
           max_xx%i = z2%i
        endif
      end function max_xx
      type (xplex) function max_xxx(z1,z2,z3)
        type (xplex),intent(in)::z1,z2,z3
        type (xplex)::max_xx
        max_xx%r = max(z1%r,z2%r)
        if (max_xx%r==z1%r) then
           max_xx%i = z1%i
        else
           max_xx%i = z2%i
        endif
        max_xxx%r = max(max_xx%r,z3%r)
        if (max_xxx%r==max_xx%r) then
           max_xxx%i = max_xx%i
        else
           max_xxx%i = z3%i
        endif
      end function max_xxx
! ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
      type (xplex) function pow_xi(z,i)
        type (xplex),intent(in) :: z
        integer,intent(in):: i
        pow_xi%r = dble(dcmplx(z%r,z%i)**i)
        pow_xi%i = dimag(dcmplx(z%r,z%i)**i)
        !pow_xi%r = z%r**i
        !pow_xi%i = z%r**(i-1)*i*z%i
      end function pow_xi
      function pow_x2i(z,i) result(C)
        type (xplex),intent(in) :: z(:,:)
        integer,intent(in):: i
        type (xplex)::C(size(z,1),size(z,2))
        C%r = dble(dcmplx(z%r,z%i)**i)
        C%i = dimag(dcmplx(z%r,z%i)**i)
        !C%r = z%r**i
        !C%i = z%r**(i-1)*i*z%i
      end function pow_x2i
      type (xplex) function pow_rx(r,z)
        type (xplex),intent(in) :: z
        real*8,intent(in):: r
        pow_rx%r = dble(r**dcmplx(z%r,z%i))
        pow_rx%i = dimag(r**dcmplx(z%r,z%i))
        !pow_rx%r = r**z%r
        !pow_rx%i = r**z%r*sin(z%i*log(r))!dimag(dcmplx(r,0d0)**dcmplx(z%r,z%i))
        !if (r==0d0) pow_rx%i=0d0
      end function pow_rx
      type (xplex) function pow_xx(z1,z2)
        type (xplex),intent(in) :: z1,z2
        !pow_xx%r = dble(dcmplx(z1%r,z1%i)**dcmplx(z2%r,z2%i))
        !pow_xx%i = dimag(dcmplx(z1%r,z1%i)**dcmplx(z2%r,z2%i))
        pow_xx%r = z1%r**z2%r
        pow_xx%i = z1%r**z2%r*sin(z2%i*log(z1%r))*sin(z2%r*z1%i/z1%r)!z1%r**z2%r*sin(z2%r*z1%i/z1%r+z2%i*log(z1%r))!   dimag(dcmplx(z1%r,z1%i)**dcmplx(z2%r,z2%i))
        if(z1%r==0d0.or.z2%r==0d0 ) pow_xx%i=0d0
      end function pow_xx
      type (xplex) function pow_xr(z,r)
        type (xplex),intent(in) :: z
        real*8,intent(in)::r
        pow_xr%r = dble(dcmplx(z%r,z%i)**r)
        pow_xr%i = dimag(dcmplx(z%r,z%i)**r)
        !pow_xr%r = z%r**r
        !pow_xr%i = z%r**r*sin(r*z%i/z%r)
        !if (z%r==0d0) pow_xr%i=0d0
      end function pow_xr
      function pow_x1r(z,r) result(C)
        type (xplex),intent(in) :: z(:)
        real*8,intent(in)::r
        type (xplex)::C(size(z,1))
        C%r = dble(dcmplx(z%r,z%i)**r)
        C%i = dimag(dcmplx(z%r,z%i)**r)
        !C%r = z%r**r
        !C%i = z%r**r*sin(r*z%i/z%r)
        !where (z%r==0d0) C%i=0d0
      end function pow_x1r
! .LT.
      logical function lt_xx(z1,z2)
        type (xplex),intent(in) :: z1,z2
        lt_xx = z1%r.lt.z2%r
      end function lt_xx
      logical function lt_rx(r,z)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        lt_rx = r.lt.z%r
      end function lt_rx
      logical function lt_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        lt_xr = z%r.lt.r
      end function lt_xr
      logical function lt_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        lt_xi = z%r.lt.i
      end function lt_xi
! .GT.
      logical function gt_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        gt_xr = z%r.gt.r
      end function gt_xr
      logical function gt_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        gt_xx = z1%r.gt.z2%r
      end function gt_xx
      logical function gt_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        gt_xi = z%r.gt.i
      end function gt_xi
      logical  function gt_ix(i,z)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        gt_ix = i.gt.z%r
      end function gt_ix
! .GE. 
      logical function ge_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        ge_xr = z%r.ge.r
      end function ge_xr
      logical function ge_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        ge_xi = z%r.ge.i
      end function ge_xi
      logical function ge_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        ge_xx = z1%r.ge.z2%r
      end function ge_xx
! .LE.
      logical function le_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        le_xi = z%r.le.i
      end function le_xi
      logical function le_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        le_xx = z1%r.le.z2%r
      end function le_xx
      logical function le_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        le_xr = z%r.le.r
      end function le_xr
      logical function le_rx(r,z)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        le_rx = r.le.z%r
      end function le_rx
      logical function le_ix(i,z)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        le_ix = i.le.z%r
      end function le_ix
! .NE. 
      logical function ne_xx(z1,z2)
        type (xplex),intent(in)::z1,z2
        ne_xx = z1%r.ne.z2%r
      end function ne_xx
      logical function ne_xr(z,r)
        type (xplex),intent(in)::z
        real*8,intent(in)::r
        ne_xr = z%r.ne.r
      end function ne_xr
      logical function ne_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        ne_xi = z%r.ne.i
      end function ne_xi
! * * * * * * * * * * * * * * * * * * * * * * * 
      type (xplex) function mult_xx(z1,z2)
        type (xplex),intent(in) :: z1,z2
        mult_xx%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        mult_xx%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !mult_xx%r = z1%r*z2%r
        !mult_xx%i = z1%i*z2%r+z1%r*z2%i
      end function mult_xx
      function mult_x2x2(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:),z2(:,:)
        type (xplex)::C(size(z1,1),size(z1,2))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x2x2
      function mult_x1x1(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:),z2(:)
        type (xplex)::C(size(z1,1))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x1x1
      function mult_x1x(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:),z2
        type (xplex)::C(size(z1,1))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x1x
      function mult_x1r(z1,r) result(C)
        type (xplex),intent(in) :: z1(:)
        real*8,intent(in)::r
        type (xplex)::C(size(z1,1))
        C%r = dble(dcmplx(z1%r,z1%i)*r)
        C%i = dimag(dcmplx(z1%r,z1%i)*r)
        !C%r = z1%r*r
        !C%i = z1%i*r
      end function mult_x1r
      function mult_xx1(z1,z2) result(C)
        type (xplex),intent(in) :: z1,z2(:)
        type (xplex)::C(size(z2,1))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_xx1
      function mult_xx2(z1,z2) result(C)
        type (xplex),intent(in) :: z1,z2(:,:)
        type (xplex)::C(size(z2,1),size(z2,2))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_xx2
      function mult_xx3(z1,z2) result(C)
        type (xplex),intent(in) :: z1,z2(:,:,:)
        type (xplex)::C(size(z2,1),size(z2,2),size(z2,3))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_xx3
      function mult_x2x(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:),z2
        type (xplex)::C(size(z1,1),size(z1,2))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x2x
      function mult_x3x(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:,:),z2
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x3x
      function mult_x3x3(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:,:),z2(:,:,:)
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x3x3
      function mult_x4x4(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:,:,:),z2(:,:,:,:)
        type (xplex)::C(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
        C%r = dble(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        C%i = dimag(dcmplx(z1%r,z1%i)*dcmplx(z2%r,z2%i))
        !C%r = z1%r*z2%r
        !C%i = z1%i*z2%r+z1%r*z2%i
      end function mult_x4x4
      type (xplex) function mult_rx(r,z)
        real*8,intent(in)::r
        type(xplex),intent(in)::z
        mult_rx%r = dble(r*dcmplx(z%r,z%i))
        mult_rx%i = dimag(r*dcmplx(z%r,z%i))
        !mult_rx%r = r*z%r
        !mult_rx%i = r*z%i
      end function mult_rx
      function mult_rx2(r,z) result(C)
        real*8,intent(in)::r
        type(xplex),intent(in)::z(:,:)
        type(xplex)::C(size(z,1),size(z,2))
        C%r = dble(r*dcmplx(z%r,z%i))
        C%i = dimag(r*dcmplx(z%r,z%i))
        !C%r = r*z%r
        !C%i = r*z%i
      end function mult_rx2
      function mult_rx1(r,z) result(C)
        real*8,intent(in)::r
        type(xplex),intent(in)::z(:)
        type(xplex)::C(size(z,1))
        C%r = dble(r*dcmplx(z%r,z%i))
        C%i = dimag(r*dcmplx(z%r,z%i))
        !C%r = r*z%r
        !C%i = r*z%i
      end function mult_rx1
      type (xplex) function mult_xr(z,r)
        type(xplex),intent(in)::z
        real*8,intent(in)::r
        mult_xr%r = dble(dcmplx(z%r,z%i)*r)
        mult_xr%i = dimag(dcmplx(z%r,z%i)*r)
        !mult_xr%r = z%r*r
        !mult_xr%i = z%i*r
      end function mult_xr
      function mult_x2r(z,r) result(C)
        type(xplex),intent(in)::z(:,:)
        real*8,intent(in)::r
        type(xplex)::C(size(z,1),size(z,2))
        C%r = dble(dcmplx(z%r,z%i)*r)
        C%i = dimag(dcmplx(z%r,z%i)*r)
        !C%r = z%r*r
        !C%i = z%i*r
      end function mult_x2r
      function mult_x3r(z,r) result(C)
        type(xplex),intent(in)::z(:,:,:)
        real*8,intent(in)::r
        type(xplex)::C(size(z,1),size(z,2),size(z,3))
        C%r = dble(dcmplx(z%r,z%i)*r)
        C%i = dimag(dcmplx(z%r,z%i)*r)
        !C%r = z%r*r
        !C%i = z%i*r
      end function mult_x3r
      type (xplex) function mult_ix(i,z)
        integer,intent(in)::i
        type(xplex),intent(in)::z
        mult_ix%r = dble(dble(i)*dcmplx(z%r,z%i))
        mult_ix%i = dimag(dble(i)*dcmplx(z%r,z%i))
        !mult_ix%r = dble(i)*z%r
        !mult_ix%i = dble(i)*z%i
      end function mult_ix
      function mult_ix2(i,z) result(C)
        integer,intent(in)::i
        type(xplex),intent(in)::z(:,:)
        type(xplex)::C(size(z,1),size(z,2))
        C%r = dble(dble(i)*dcmplx(z%r,z%i))
        C%i = dimag(dble(i)*dcmplx(z%r,z%i))
        !C%r = dble(i)*z%r
        !C%i = dble(i)*z%i
      end function mult_ix2
      type (xplex) function mult_xi(z,i)
        type(xplex),intent(in)::z
        integer,intent(in)::i
        mult_xi%r = dble(dcmplx(z%r,z%i)*dble(i))
        mult_xi%i = dimag(dcmplx(z%r,z%i)*dble(i))
        !mult_xi%r = z%r*dble(i)
        !mult_xi%i = z%i*dble(i)
      end function mult_xi
! + + + + + + + + + + + + + + + + + + + + + + +
      type (xplex) function add_xx(z1,z2)
        type (xplex),intent(in) :: z1,z2
        add_xx%r = z1%r+z2%r
        add_xx%i = z1%i+z2%i
      end function add_xx
      function add_x1x1(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:),z2(:)
        type (xplex):: C(size(z1,1))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x1x1
      function add_x1x(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:),z2
        type (xplex):: C(size(z1,1))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x1x
      function add_x1i(z1,i) result(C)
        type (xplex),intent(in) :: z1(:)
        integer,intent(in)::i
        type (xplex):: C(size(z1,1))
        C%r = z1%r+dble(i)
        C%i = z1%i
      end function add_x1i
      function add_x1r(z1,r) result(C)
        type (xplex),intent(in) :: z1(:)
        real*8,intent(in)::r
        type (xplex):: C(size(z1,1))
        C%r = z1%r+r
        C%i = z1%i
      end function add_x1r
      function add_xx1(z1,z2) result(C)
        type (xplex),intent(in) :: z1,z2(:)
        type (xplex):: C(size(z2,1))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_xx1
      function add_xx2(z1,z2) result(C)
        type (xplex),intent(in) :: z1,z2(:,:)
        type (xplex):: C(size(z2,1),size(z2,2))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_xx2
      function add_x2x2(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:),z2(:,:)
        type (xplex):: C(size(z1,1),size(z1,2))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x2x2
      function add_x2x(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:),z2
        type (xplex):: C(size(z1,1),size(z1,2))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x2x
      function add_x2r(z1,r) result(C)
        type (xplex),intent(in) :: z1(:,:)
        real*8,intent(in)::r
        type (xplex):: C(size(z1,1),size(z1,2))
        C%r = z1%r+r
        C%i = z1%i
      end function add_x2r
      function add_x3x3(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:,:),z2(:,:,:)
        type (xplex):: C(size(z1,1),size(z1,2),size(z1,3))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x3x3
      function add_x4x4(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:,:,:),z2(:,:,:,:)
        type (xplex):: C(size(z1,1),size(z1,2),size(z1,3),size(z1,4))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x4x4
      function add_x3x(z1,z2) result(C)
        type (xplex),intent(in) :: z1(:,:,:),z2
        type (xplex):: C(size(z1,1),size(z1,2),size(z1,3))
        C%r = z1%r+z2%r
        C%i = z1%i+z2%i
      end function add_x3x
      type (xplex) function add_xr(z,r)
        type (xplex),intent(in)::z
        real*8, intent(in)::r
        add_xr%r = z%r+r
        add_xr%i = z%i
      end function add_xr
      type (xplex) function add_xi(z,i)
        type (xplex),intent(in)::z
        integer,intent(in)::i
        add_xi%r = z%r+dble(i)
        add_xi%i = z%i
      end function add_xi
      function add_x2i2(z,i) result(C)
        type (xplex),intent(in)::z(:,:)
        integer,intent(in)::i(:,:)
        type (xplex)::C(size(z,1),size(z,2))
        C%r = z%r+dble(i)
        C%i = z%i
      end function add_x2i2
      type (xplex) function add_rx(r,z)
        type (xplex),intent(in)::z
        real*8, intent(in)::r
        add_rx%r = z%r+r
        add_rx%i = z%i
      end function add_rx
      type (xplex) function add_ix(i,z)
        type (xplex),intent(in)::z
        integer, intent(in)::i
        add_ix%r = z%r+dble(i)
        add_ix%i = z%i
      end function add_ix

      end module complexify
