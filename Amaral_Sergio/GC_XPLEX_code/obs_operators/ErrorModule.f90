! $Id: ErrorModule.f90,v 1.2 2009/06/23 06:47:07 daven Exp $
MODULE ErrorModule

  !========================================================================
  ! Module ErrorModule contains error check routines for the Fortran code 
  ! that reads netCDF data from disk. (bmy, 2/15/07)
  ! 
  ! Module Methods:
  ! -----------------------------------------------------------------------
  ! (1 ) AllocErr           : Prints an error message for allocating arrays
  ! (2 ) ErrMsg             : Prints an error message and halts execution
  ! (3 ) Msg                : Prints a message and flushes buffer
  ! (4 ) CheckValue         : Checks a value for NaN or Infinity condition
  ! (5 ) ReplaceNanAndInfR4 : Replaces a TYPE (XPLEX) Nan/Inf value w/ other data
  ! (6 ) ReplaceNanAndInfR8 : Checks for NaN
  ! (7 ) ItIsFiniteR4       : Checks a TYPE (XPLEX) value for Infinity
  ! (8 ) ItIsFiniteR8       : Checks a TYPE (XPLEX) value for Infinity
  ! (9 ) ItIsNanR4          : Checks a TYPE (XPLEX) value for NaN
  ! (10) ItIsNanR8          : Checks a TYPE (XPLEX) value for NaN
  !
  ! Module Interfaces:
  ! -----------------------------------------------------------------------
  ! (1 ) ReplaceNanAndInf   : Overloads ReplaceNanAndInfR4, ReplaceNanAndInfR8
  ! (2 ) ItIsFinite         : Overloads ItIsFiniteR4, ItIsFiniteR8
  ! (3 ) ItIsNan            : Overloads ItIsNanR4, ItIsNanR8
  !
  ! NOTES:
  ! (1 ) Adapted from He4ErrorModule.f90 (bmy, 2/15/07)
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
  PUBLIC :: AllocErr
  PUBLIC :: ErrMsg
  PUBLIC :: Msg
  PUBLIC :: CheckValue
  PUBLIC :: ItIsNan
  PUBLIC :: ItIsFinite
  PUBLIC :: ReplaceNanAndInf

  !-------------------------------
  ! MODULE INTERFACES
  !-------------------------------
  
  INTERFACE ReplaceNanAndInf
!     MODULE PROCEDURE ReplaceNanAndInfR4
     MODULE PROCEDURE ReplaceNanAndInfR8
  END INTERFACE

  INTERFACE ItIsNan
!     MODULE PROCEDURE ItIsNanR4
     MODULE PROCEDURE ItIsNanR8
  END INTERFACE

  INTERFACE ItIsFinite
!     MODULE PROCEDURE ItIsFiniteR4
     MODULE PROCEDURE ItIsFiniteR8
  END INTERFACE

  !-------------------------------
  ! MODULE ROUTINES
  !-------------------------------

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE AllocErr( arrayName )

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
      
  END SUBROUTINE AllocErr

!------------------------------------------------------------------------------

  SUBROUTINE ErrMsg( msg, loc )

    !======================================================================
    ! Subroutine ErrMsg halts displays an error message and halts
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
    ! ErrMsg begins here!
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

  END SUBROUTINE ErrMsg

!------------------------------------------------------------------------------

  SUBROUTINE Msg( str )

    !======================================================================
    ! Subroutine Msg prints a string and flushes the output buffer.
    ! (bmy, 1/17/06)
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) str (CHARACTER) : Message to display
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: str

    !---------------------
    ! Msg begins here!
    !---------------------

    ! Print message
    WRITE( 6, '(a)' ) TRIM( str )
    CALL flush( 6 )
      
  END SUBROUTINE Msg

!-----------------------------------------------------------------------------

  SUBROUTINE CheckValue( value, name, loc )

    !======================================================================
    ! Subroutine CheckValue tests a value for IEEE NaN or Infinity.
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

  END SUBROUTINE CheckValue

!-----------------------------------------------------------------------------

  SUBROUTINE ReplaceNanAndInfR4( value, replacement )

    !======================================================================
    ! Subroutine ReplaceNaNandInfR4 replaces a NaN or infinity TYPE (XPLEX) 
    ! value with a replacement value.  You can use this to assign missing
    ! data flags such as -9999 to NaN or infinity values. (bmy, 2/15/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) value       (TYPE (XPLEX)) : Value to be tested
    ! (2 ) replacement (TYPE (XPLEX)) : Replacement value
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (1 ) value       (TYPE (XPLEX)) : Value is overwritten and returned
    !======================================================================

    ! Arguments
    TYPE (XPLEX), INTENT(INOUT) :: value
    TYPE (XPLEX), INTENT(IN)    :: replacement

    !----------------------------------
    ! ReplaceNanAndInfR4 begins here!
    !----------------------------------
    
    IF ( ItIsNan( value ) ) THEN
       value = replacement
    ELSE IF ( .not. ItIsFinite( value ) ) THEN
       value = replacement
    ENDIF

  END SUBROUTINE ReplaceNanAndInfR4

!-----------------------------------------------------------------------------

  SUBROUTINE ReplaceNanAndInfR8( value, replacement )

    !======================================================================
    ! Subroutine ReplaceNaNandInfR8 replaces a NaN or infinity TYPE (XPLEX)
    ! value with a replacement value.  You can use this to assign missing
    ! data flags such as -9999 to NaN or infinity values. (bmy, 2/15/07)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) value       (TYPE (XPLEX)) : Value to be tested
    ! (2 ) replacement (TYPE (XPLEX)) : Replacement value
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (1 ) value       (TYPE (XPLEX)) : Value is overwritten and returned
    !======================================================================

    ! Arguments
    TYPE (XPLEX), INTENT(INOUT) :: value
    TYPE (XPLEX), INTENT(IN)    :: replacement

    !----------------------------------
    ! ReplaceNanAndInfR8 begins here!
    !----------------------------------
    
    IF ( ItIsNan( value ) ) THEN
       value = replacement
    ELSE IF ( .not. ItIsFinite( value ) ) THEN
       value = replacement
    ENDIF

  END SUBROUTINE ReplaceNanAndInfR8

!-----------------------------------------------------------------------------

  FUNCTION ItIsNanR4( value ) RESULT( itIsANaN )

    !===================================================================
    ! Subroutine ItIsNanR4 tests a TYPE (XPLEX) value for IEEE NaN on SGI, 
    ! Altix, Linux, or Sun platforms.  (bmy, 2/15/07)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)) : value to be tested
    !===================================================================

#   include "define.h"

    ! Argument
    TYPE (XPLEX), INTENT(IN) :: value
    LOGICAL            :: itIsANaN

    !-------------------------
    ! ItIsNanR4 begins here!
    !-------------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    ItIsANan = IEEE_IS_NAN( value )   

#elif defined( ALTIX ) || defined( PC )

    ! Declare IS_NAN as an external function
    INTEGER, EXTERNAL  :: IS_NAN

    ! For LINUX or IFORT compilers, use C routine "is_nan" to test for NaN
    ! VALUE must be cast to DBLE since "is_nan" only takes doubles.
    ItIsANan = ( IS_NAN( DCMPLX( value ) ) /= 0 )

#elif defined( SPARC )

    ! Declare Sun intrinsic IR_ISNAN as an external function
    INTEGER, EXTERNAL :: IR_ISNAN

    ! Test if VALUE is a NaN
    ItIsANan = ( IR_ISNAN( value ) /= 0 )

#endif

  END FUNCTION ItIsNanR4

!-----------------------------------------------------------------------------

  FUNCTION ItIsNanR8( value ) RESULT( ItIsANan )

    !===================================================================
    ! Subroutine ItIsNanR8 tests a TYPE (XPLEX) value for IEEE NaN on SGI, 
    ! Altix, Linux, or Sun platforms. (bmy, 2/15/07)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)) : value to be tested
    !===================================================================

    ! Argument
    TYPE (XPLEX), INTENT(IN) :: value
    LOGICAL            :: ItIsANan

    !-------------------------
    ! ItisNanR8 begins here!
    !-------------------------

#if   defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    ItIsANan = IEEE_IS_NAN( value )   

#elif defined( ALTIX ) || defined( PC )

    ! Declare IS_NAN as an external function
    INTEGER, EXTERNAL  :: IS_NAN

    ! For LINUX or IFORT compilers, use C routine "is_nan" to test for NaN
    ! VALUE must be cast to DBLE since "is_nan" only takes doubles.
    ItIsANan = ( is_nan( value ) /= 0 )

#elif defined( SPARC )

    ! Declare ID_ISNAN as an external function
    INTEGER, EXTERNAL :: ID_ISNAN

    ! Test if VALUE is a NaN
    ItIsANan = ( ID_ISNAN( value ) /= 0 )

#endif

  END FUNCTION ItIsNanR8

!-----------------------------------------------------------------------------

  FUNCTION ItIsFiniteR4( value ) RESULT( itIsAFinite )

    !===================================================================
    ! Subroutine ItIsFiniteR4 tests a TYPE (XPLEX) value for IEEE Finite on 
    ! SGI, Altix, Linux, or Sun platforms. (bmy, 2/15/07)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)) : value to be tested
    !===================================================================

    ! Arguments
    TYPE (XPLEX), INTENT(IN) :: value
    LOGICAL            :: ItIsAFinite

    !----------------------------
    ! ItIsFiniteR4 begins here!
    !----------------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    ItIsAFinite = IEEE_FINITE( value )   

#elif defined( ALTIX ) || defined( PC )  

    ! Declare IS_FINITE as an external function
    INTEGER, EXTERNAL  :: IS_FINITE

    ! For LINUX or INTEL_FC compilers, use C routine "is_finite" to test if 
    ! VALUE is finite.   VALUE must be cast to DBLE since "is_finite" only
    ! takes doubles.
    ItIsAFinite = ( IS_FINITE( DCMPLX( value ) ) /= 0 )

#elif defined( SPARC )

    ! Declare Sun intrinsic IR_FINITE as an external function
    INTEGER, EXTERNAL :: IR_FINITE

    ! Test if VALUE is a finite number
    ItIsAFinite = ( IR_FINITE( VALUE ) /= 0 )

#endif

  END FUNCTION ItIsFiniteR4

!-----------------------------------------------------------------------------

  FUNCTION ItIsFiniteR8( value ) RESULT( itIsAFinite )

    !===================================================================
    ! Subroutine ItIsFiniteR8 tests a TYPE (XPLEX) value for IEEE Finite on 
    ! SGI, Altix, Linux, or Sun platforms. (bmy, 2/15/07)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (TYPE (XPLEX)) : value to be tested
    !===================================================================

    ! Arguments
    TYPE (XPLEX), INTENT(IN) :: value
    LOGICAL            :: ItIsAFinite

    !----------------------------
    ! ItIsFiniteR4 begins here!
    !----------------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    ItIsAFinite = IEEE_FINITE( value )   

#elif defined( ALTIX ) || defined( PC )  

    ! Declare IS_FINITE as an external function
    INTEGER, EXTERNAL  :: IS_FINITE

    ! For Altix or Linux compilers, use C routine 
    ! "is_finite" to test if VALUE is finite.   
    ItIsAFinite = ( IS_FINITE( value ) /= 0 )

#elif defined( SPARC )

    ! Declare Sun intrinsic ID_FINITE as an external function
    INTEGER, EXTERNAL :: ID_FINITE

    ! Test if VALUE is a finite number
    ItIsAFinite = ( ID_FINITE( VALUE ) /= 0 )

#endif

  END FUNCTION ItIsFiniteR8

!-----------------------------------------------------------------------------

END MODULE ErrorModule
