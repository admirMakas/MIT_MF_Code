!$Id: directory_adj_mod.f,v 1.2 2009/06/12 01:44:48 daven Exp $
      MODULE DIRECTORY_ADJ_MOD
!
!******************************************************************************
!  Module DIRECTORY_ADJ_MOD contains the directory path variables used by 
!  GEOS-CHEM adjoint code. (adj_group, 6/07/09)
!     
!  Module Variables:
!  ============================================================================
!  (1 ) OPTDATA_DIR      (CHAR*255) : Directory with necessary adj output
!  (2 ) ADJTMP_DIR       (CHAR*255) : Directory with temp/intermediate adj output
!  (3 ) DIAGADJ_DIR      (CHAR*255) : Directory with adj diagnostics 
!
!  NOTES:
!*****************************************************************************

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      CHARACTER(LEN=255) :: OPTDATA_DIR
      CHARACTER(LEN=255) :: ADJTMP_DIR
      CHARACTER(LEN=255) :: DIAGADJ_DIR

      ! End of module
      END MODULE DIRECTORY_ADJ_MOD
