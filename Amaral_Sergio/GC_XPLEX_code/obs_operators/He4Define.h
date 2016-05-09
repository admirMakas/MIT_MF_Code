! $Id: He4Define.h,v 1.1 2009/07/03 01:55:32 daven Exp $
!
! "He4Define.h" -- sets Cpp flags for HD4 and HDF-EOS4 code (bmy, 3/21/08)

!-----------------------------------------------------
! Pick the compiler/architecture that you are using
!-----------------------------------------------------
!#define INTEL32 'INTEL32'
#define INTEL64  'INTEL64'
!#define SGI32   'SGI32'
!#define SGI64   'SGI64'
!#define SUN32   'SUN32'
!#define SUN64   'SUN64'

!-----------------------------------------------------
! For HDF5-EOS you may have to pass 64-bit integer 
! (INTEGER*8) dimensions to some routines.  If so,
! then you can #define NEED_INT_64.
!
! However, for HDF4-EOS, it seems like 32-bit 
! integers (INTEGER*4) are used for dimensioning
! variables.  
!-----------------------------------------------------

#if defined( INTEL32 ) || defined( INTEL64 )

! IFORT compiler 
#define NEED_INT_32 'NEED_INT_32'

#elif defined( SGI32 ) || defined( SGI64 )

! 32-bit SGI 
#define NEED_INT_32 'NEED_INT_32'

#elif defined( SUN32 ) || defined( SUN64 )

! SPARC or SunStudio compiolers
#define NEED_INT_32 'NEED_INT_32'

#endif
