! $Id: cleanup_adj.f,v 1.3 2012/03/01 22:00:26 daven Exp $
      SUBROUTINE CLEANUP_ADJ
!
!******************************************************************************
!  Subroutine CLEANUP_ADJ deallocates the memory assigned to dynamic allocatable 
!  arrays in adjoint model routines (dkh, 06/12/09) 
!
!  NOTES:
!  (1 ) Based on CLEANUP
!  (2 ) Add support for CH4 (kjw, dkh, 02/12/12, adj32_023) 
!******************************************************************************
!
      ! References to F90 modules 
      USE ADJ_ARRAYS_MOD,          ONLY : CLEANUP_ADJ_ARRAYS
      USE GLOBAL_CH4_ADJ_MOD,      ONLY : CLEANUP_GLOBAL_CH4_ADJ
      USE POPULATION_MOD,          ONLY : CLEANUP_POPULATION_MOD


      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! CLEANUP_ADJ begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) 
 100  FORMAT( '     - CLEANUP_ADJ: deallocating arrays now...' )

      ! Call cleanup routines from individual F90 modules
      CALL CLEANUP_ADJ_ARRAYS
      CALL CLEANUP_GLOBAL_CH4_ADJ
      CALL CLEANUP_POPULATION_MOD
      

      ! Return to calling program
      END SUBROUTINE CLEANUP_ADJ
