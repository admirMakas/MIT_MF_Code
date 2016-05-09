! $Id: rdland.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      SUBROUTINE RDLAND
!
!******************************************************************************
!  Subroutine RDLAND reads the land types and fractions (times 1000)
!  from the "vegtype.global" file.  (yhw, gmg, djj, 1994; bmy, 12/20/04)
!
!  Common-block variables from header file "CMN_DEP":
!  ============================================================================
!  (1 ) FRCLND(I,J)        : Land fraction (0.0 - 1.0)
!  (2 ) IREG(I,J)          : Number of landtypes in each grid box
!  (3 ) ILAND(I,J,LDT)     : Land type ID for element LDT =1, IREG(I,J)
!  (4 ) IUSE(I,J,LDT)      : Fraction (per mil) of gridbox area occupied by 
!                             land type element LDT
!
!  Common-block variables from header file "CMN_VEL":
!  ============================================================================
!  (1 ) IJREG(IJLOOP)      : 2-D (I*J, LDT) version of IJREG  (for DEPVEL)
!  (2 ) IJLAND(IJLOOP,LDT) : 2-D (I*J, LDT) version of IJLAND (for DEPVEL)
!  (3 ) IJUSE(IJLOOP,LDT)  : 2-D (I*J, LDT) version of IJUSE  (for DEPVEL)
!
!  NOTES:
!  (1 ) Now read the "vegtype.global" file from the leaf_area_index_200412 
!        subdirectory of DATA_DIR.  This is the same Olson land map as was
!        used previously.  Also updated comments and added standard GEOS-CHEM 
!        program documentation header. (tmf, bmy, 12/6/04)
!  (2 ) Now read the "vegtype.global" file from the leaf_area_index_200412
!        subdirectory if LAVHRRLAI=T.  Also updated comments and added 
!        standard GEOS-CHEM program documentation header. (bmy, 12/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GRID_MOD,      ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,   ONLY : LAVHRRLAI

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DEP"   ! FRCLND, IREG, ILAND, IUSE
#     include "CMN_VEL"   ! IJREG, IJLAND, IJUSE

      ! Local variables
      INTEGER :: I, J, K, IJLOOP, IREF, JREF
      INTEGER :: I0, J0

      ! For filename
      CHARACTER(LEN=255) :: FILENAME
 
      !=================================================================
      ! RDLAND begins here!
      !=================================================================

      ! Get nested-grid offsets (bmy, 2/11/03)
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      ! Read the "vegtype.global" from the proper directory
      ! depending on the setting of the LAVHRRLAI flag (bmy, 12/20/04)
      IF ( LAVHRRLAI ) THEN
         FILENAME = TRIM( DATA_DIR ) // 
     &              'leaf_area_index_200412/vegtype.global'
      ELSE
         FILENAME = TRIM( DATA_DIR ) // 
     &              'leaf_area_index_200202/vegtype.global'         
      ENDIF


      WRITE( 6, 50 ) TRIM( FILENAME )
 50   FORMAT( '     - RDLAND: Reading ', a )

      ! Open the file 
      OPEN( 65, FILE=TRIM( FILENAME ), STATUS='OLD',
     &          FORM='FORMATTED',      ERR=700 )

      ! Read data
 100  READ(65,101,end=110,ERR=800) I,J,IREG(I,J),
     &     (ILAND(I,J,K),K=1,IREG(I,J)),
     &     (IUSE(I,J,K),K=1,IREG(I,J))
 101  FORMAT(20I4)
      GO TO 100

      ! Process data into arrays
 110  CONTINUE
      CLOSE (65)
      IJLOOP = 0
      DO 500 J = 1, JJPAR
         JREF = J + J0
         DO 400 I = 1, IIPAR
            FRCLND(I,J) = 1000.d0
            IREF = I + I0
            IJLOOP = IJLOOP + 1
            IJREG(IJLOOP) = IREG(IREF,JREF)
            DO 300 K=1,IJREG(IJLOOP)
               IJLAND(IJLOOP,K) = ILAND(IREF,JREF,K)
               IJUSE(IJLOOP,K)  = IUSE(IREF,JREF,K)
               IF (IJLAND(IJLOOP,K) .EQ. 0 )
     &              FRCLND(I,J) = FRCLND(I,J) - IJUSE(IJLOOP,K)
 300        CONTINUE
            FRCLND(I,J) = FRCLND(I,J) / 1000.d0
 400     CONTINUE
 500  CONTINUE
      
      ! Return
      RETURN

      ! Trap File open error
 700  CONTINUE
      CALL ERROR_STOP( 'Error opening "vegtype.global"', 'rdland.f' )
      
      ! Trap file read error
 800  CONTINUE
      CALL ERROR_STOP( 'Error reading "vegtype.global"', 'rdland.f' )
      print*,'FRCLND',FRCLND
      ! Return to calling program
      END SUBROUTINE RDLAND
