! $Id: He4ErrorModule.f90,v 1.1 2009/06/18 19:53:07 daven Exp $
MODULE He4ErrorModule

  !========================================================================
  ! Module He4ErrorModule contains error check routines for the
  ! Fortran code that reads HDF4 and HDF-EOS4 data from disk.
  ! (bmy, 7/26/00, 3/21/08)
  ! 
  ! Module Methods:
  ! -----------------------------------------------------------------------
  ! (1 ) He4AllocErr   : Prints an error message for allocating arrays
  ! (2 ) He4ErrMsg     : Prints an error message and halts execution
  ! (3 ) He4Msg        : Prints a message and flushes buffer
  ! (4 ) He4CheckValue : Checks a value for NaN or Infinity condition
  ! (5 ) ItIsNan       : Checks for NaN
  ! (6 ) ItIsFinite    : Checks for Infinity
  !
  ! NOTES:
  ! (1 ) Now use intrinsic functions ISNAN and FP_CLASS to test for
  !       NaN and Infinity on IFORT compiler.  These functions were not
  !       available in the older EFC compiler. (bmy, 8/14/07)
  ! (2 ) Now uses updated flags from He4Define.h (bmy, 3/21/08)
  !========================================================================
  
  USE MYTYPE
  USE COMPLEXIFY
  IMPLICIT NONE
  
  !-------------------------------
  ! PRIVATE / PUBLIC DECLARATIONS
  !-------------------------------

  ! Make everything PRIVATE ...
  PRIVATE

  ! ... except these routines
  PUBLIC :: He4AllocErr
  PUBLIC :: He4ErrMsg
  PUBLIC :: He4Msg
  PUBLIC :: He4CheckValue
  PUBLIC :: ItIsNan
  PUBLIC :: ItIsFinite

  !-------------------------------
  ! MODULE ROUTINES
  !-------------------------------

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE He4AllocErr( arrayName )

    !======================================================================
    ! Subroutine He4AllocErr stops program execution upon an error
    ! allocating arrays. (bmy, 1/17/06)
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) arrayName (CHARACTER) : Name of array
    !
    ! NOETS:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: arrayName

    !----------------------------
    ! He4AllocErr begins here!
    !----------------------------

    ! Write info
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 100   ) TRIM( arrayName )
    WRITE( 6, 110   )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    CALL FLUSH( 6 )

    ! Exit
    CALL EXIT( 1 )
    
    ! FORMAT strings
100 FORMAT( 'Allocation error for array ', a )
110 FORMAT( 'STOP in allocErr ("Hdf4ErrorModule.f90")' )
      
  END SUBROUTINE He4AllocErr

!------------------------------------------------------------------------------

  SUBROUTINE He4ErrMsg( msg, loc )

    !======================================================================
    ! Subroutine He4ErrMsg halts displays an error message and halts
    ! program execution. (bmy, 1/17/06)
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) msg (CHARACTER) : Error message to display
    ! (2 ) loc (CHARACTER) : Location where the error occurred
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: msg
    CHARACTER(LEN=*), INTENT(IN) :: loc

    !--------------------------
    ! He4ErrMsg begins here!
    !--------------------------

    ! Print error message
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, '(a)' ) TRIM( msg )
    WRITE( 6, 100   ) TRIM( loc )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    CALL FLUSH( 6 )
    
    ! Exit simulation
    CALL EXIT( 1 )
    
    ! FORMAT string
100 FORMAT( 'STOP in ', a )

  END SUBROUTINE He4ErrMsg

!------------------------------------------------------------------------------

  SUBROUTINE He4Msg( str )

    !======================================================================
    ! Subroutine He4Msg prints a string and flushes the output buffer.
    ! (bmy, 1/17/06)
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) str (CHARACTER) : Message to display
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: str

    !---------------------
    ! He4Msg begins here!
    !---------------------

    ! Print message
    WRITE( 6, '(a)' ) TRIM( str )
    CALL flush( 6 )
      
  END SUBROUTINE He4Msg

!-----------------------------------------------------------------------------

  SUBROUTINE He4CheckValue( value, name, loc )

    !======================================================================
    ! Subroutine He4CheckValue tests a value for IEEE NaN or Infinity.
    ! (bmy, 1/17/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)   ) : value to be tested
    ! (2 ) name  (CHARACTER) : name of the variable 
    ! (3 ) loc   (INTEGER  ) : Grid box location (/i,j,l,t/)
    !======================================================================

    ! Arguments
    TYPE (XPLEX),           INTENT(IN) :: value
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER,          INTENT(IN) :: loc(4)

    ! If VALUE is NaN, stop w/ error message
    IF ( itIsNaN( value ) ) THEN
!!$OMP CRITICAL
       WRITE( 6, 100 ) TRIM( name ), loc
 100   FORMAT( a, ' is NaN at grid box: ', 4i4, '!' )
       STOP
!!$OMP END CRITICAL
    ENDIF

    ! If VALUE is +/- Infinity, stop w/ error message
    IF ( .not. itIsFinite( value ) ) THEN
!!$OMP CRITICAL
       WRITE( 6, 110 ) TRIM( name ), loc
 110   FORMAT( a, ' is +/- Infinity at grid box: ', 4i4, '!' )
       STOP
!!$OMP END CRITICAL
    ENDIF

  END SUBROUTINE He4CheckValue

!-----------------------------------------------------------------------------

  FUNCTION ItIsNan( value ) RESULT( itIsANaN )

    !===================================================================
    ! Subroutine itIsNaN tests a value for IEEE NaN on SGI, Altix, 
    ! Linux, or Sun platforms. (bmy, 1/17/06, 8/14/07)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)) : value to be tested
    !
    ! NOTES:
    ! (1 ) Add error checking for Sun/SPARC compiler (bmy, 2/15/07)    
    ! (2 ) Now use FP_CLASS function for IFORT compiler (bmy, 8/14/07)
    !===================================================================

#include "He4Define.h"

    ! Argument
    TYPE (XPLEX), INTENT(IN) :: value
    LOGICAL            :: itIsANaN

    !----------------------
    ! ItIsNan begins here!
    !----------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    itIsANaN = IEEE_IS_NAN( value )   

#elif defined( INTEL32 ) || defined( INTEL64 )

    ! Use Intel/IFORT intrinsic function ISNAN
    itIsANan = ISNAN( value )

#elif defined( SUN32 ) || defined( SUN64 )

    ! Declare Sun intrinsic IR_ISNAN as an external function
    INTEGER, EXTERNAL :: IR_ISNAN

    ! Test if VALUE is a NaN
    ItIsANan = ( IR_ISNAN( value ) /= 0 )

#endif

  END FUNCTION ItIsNan

!-----------------------------------------------------------------------------

  FUNCTION ItIsFinite( value ) RESULT( itIsAFinite )

    !===================================================================
    ! Subroutine itIsFinite tests a value for IEEE Finite on SGI,
    ! Altix, Linux, or Sun platforms. (bmy, 1/17/06, 8/14/07)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)) : value to be tested
    !
    ! NOTES:
    ! (1 ) Add error checking for Sun/SPARC compiler (bmy, 2/15/07)
    ! (2 ) Now use FP_CLASS function for IFORT compiler (bmy, 8/14/07)
    !===================================================================

#   include "He4Define.h"

#if defined( INTEL32 ) || defined( INTEL64 ) 
#   include "fordef.for"
    INTEGER :: fpc
#endif

    ! Arguments
    TYPE (XPLEX), INTENT(IN) :: value
    LOGICAL            :: itIsAFinite

    !-------------------------
    ! itisFinite begins here!
    !-------------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    itIsAFinite = IEEE_FINITE( value )   

#elif defined( INTEL32 ) || defined( INTEL64 )  

    ! Get the floating point type class for VALUE
    fpc         = FP_CLASS( value )

    ! VALUE is infinite if it is either +Inf or -Inf
    ! Also flag an error if VALUE is a signaling or quiet NaN
    itIsAFinite = ( fpc /= FOR_K_FP_POS_INF .and. &
                    fpc /= FOR_K_FP_NEG_INF .and. &
                    fpc /= FOR_K_FP_SNAN    .and. &
                    fpc /= FOR_K_FP_QNAN         )

#elif defined( SUN32 ) || defined( SUN64 )  

    ! Declare Sun intrinsic IR_FINITE as an external function
    INTEGER, EXTERNAL :: IR_FINITE

    ! Test if VALUE is a finite number
    ItIsAFinite = ( IR_FINITE( value ) /= 0 )

#endif

  END FUNCTION ItIsFinite

!-----------------------------------------------------------------------------

END MODULE He4ErrorModule
