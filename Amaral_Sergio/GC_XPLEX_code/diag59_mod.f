! $Id: diag59_mod.f,v 1.1 2011/02/23 00:08:47 daven Exp $
      MODULE DIAG59_MOD
!
!******************************************************************************
!  Module DIAG59_MOD contains arrays and routines for archiving the ND59
!  diagnostic -- concentrations of NH3 [ug/m3]. (lz,10/07/10)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD59 (TYPE (XPLEX))  : Array for NH3 concentrations [ug/m3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) DIAG59         : Archives quantities for diagnostic
!  (2 ) ZERO_DIAG59    : Sets all module arrays to zero
!  (3 ) WRITE_DIAG59   : Writes data in module arrays to bpch file
!  (4 ) INIT_DIAG59    : Allocates all module arrays
!  (5 ) CLEANUP_DIAG59 : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag03_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f    : Module w/ NaN and other error check routines
!  (3 ) file_mod.f     : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f     : Module w/ horizontal grid information
!  (5 ) pressure_mod.f : Module w/ routines to compute P(I,J,L)
!  (6 ) time_mod.f     : Module w/ routines to compute date & time
!
!  NOTES:
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag59_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND59, LD59

      ! Parameters
      INTEGER, PARAMETER   :: PD59 = 6

      ! Arrays
      TYPE (XPLEX),  ALLOCATABLE :: AD59(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE DIAG59
!
!******************************************************************************
!  Subroutine DIAG59 archives NH3 concentrations [ug/m3] for the ND59
!  diagnostic. (lz,10/07/10)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AIRVOL, T
      !USE DIAG_MOD,     ONLY : LTOTH
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTNH3
      USE TRACERID_MOD, ONLY : IDTNH4
      USE TRACERID_MOD, ONLY : IDTNIT
      USE TRACERID_MOD, ONLY : IDTSO4
      USE TRACERID_MOD, ONLY : IDTBCPI,IDTBCPO

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! NDxx flags

      ! Local variables
      INTEGER               :: I,      J,    L
      TYPE (XPLEX)                :: FACTOR, PRES

      ! Factor for computing standard volume
      TYPE (XPLEX),PARAMETER::STD_VOL_FAC=xplex(1013.25d0/273.15d0,0d0)
     
      !================================================================= 
      ! DIAG59 begins here! 
      !================================================================= 

      ! Error check
       IF ( IDTNH3 == 0 ) RETURN
       IF ( IDTNH4 == 0 ) RETURN
       IF ( IDTNIT == 0 ) RETURN
       IF ( IDTSO4 == 0 ) RETURN
       IF ( IDTBCPI == 0 ) RETURN
       IF ( IDTBCPO == 0 ) RETURN

      ! Loop over grid boxes     
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FACTOR, PRES )  
      DO L = 1, LD59  
      DO J = 1, JJPAR 
      DO I = 1, IIPAR

         ! Conversion factor from [kg] --> [ug/m3]
         ! (LTOTH=1 if between OTH_HR1 and OTH_HR2, LTOTH=0 otherwise)
         !FACTOR        = 1d9 / AIRVOL(I,J,L) * LTOTH(I,J) 

         ! Conversion factor from [kg] --> [ug/m3]
         FACTOR        = 1d9 / AIRVOL(I,J,L)

         ! NH3 [ug/m3]
         AD59(I,J,L,1) = AD59(I,J,L,1) + 
     &                   ( STT(I,J,L,IDTNH3) * FACTOR)

         ! NH4 [ug/m3]
         AD59(I,J,L,2) = AD59(I,J,L,2) +
     &                   ( STT(I,J,L,IDTNH4) * FACTOR)

         ! NIT [ug/m3]
         AD59(I,J,L,3) = AD59(I,J,L,3) +
     &                   ( STT(I,J,L,IDTNIT) * FACTOR)

         ! SO4 [ug/m3]
         AD59(I,J,L,4) = AD59(I,J,L,4) +
     &                   ( STT(I,J,L,IDTSO4) * FACTOR)

         ! BCPI [ug/m3]
         AD59(I,J,L,5) = AD59(I,J,L,5) +
     &                   ( STT(I,J,L,IDTBCPI) * FACTOR)

         ! BCPO [ug/m3]
         AD59(I,J,L,6) = AD59(I,J,L,6) +
     &                   ( STT(I,J,L,IDTBCPO) * FACTOR)


      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO 

      ! Return to calling program
      END SUBROUTINE DIAG59

!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG59
!
!******************************************************************************
!  Subroutine ZERO_DIAG59 zeroes the ND03 diagnostic arrays. 
!  (dkh, bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! ZERO_DIAG59 begins here!
      !=================================================================

      ! Exit if ND59 is turned off
      IF ( ND59 == 0 ) RETURN

      ! Zero arrays
      AD59(:,:,:,:) = 0d0

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG59

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG59
!
!******************************************************************************
!  Subroutine WRITE_DIAG03 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 5/22/06, 9/5/06)
!
!   # : Field   : Description                 : Units    : Scale factor
!  -----------------------------------------------------------------------
!  (1 ) IJ-ugm3 : NH3                         : ug/m3    : SCALE_OTH
!  (2 ) IJ-ugm3 : NH4                         : ug/m3    : SCALE_OTH
!  (3 ) IJ-ugm3 : NIT                         : ug/m3    : SCALE_OTH
!  (4 ) IJ-ugm3 : SO4                         : ug/m3    : SCALE_OTH
!  (5 ) IJ-ugm3 : BCPI                        : ug/m3    : SCALE_OTH
!  (6 ) IJ-ugm3 : BCPO                        : ug/m3    : SCALE_OTH
!
!  NOTES:
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems  on SUN 4100 platform
!        (bmy, 9/5/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      !USE DIAG_MOD,     ONLY : CTOTH
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_DYN,  GET_DIAGb,  GET_DIAGe

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! TINDEX

      ! Local variables
      INTEGER               :: CENTER180, HALFPOLAR
      INTEGER               :: L,         M,         N
      INTEGER               :: IFIRST,    JFIRST,    LFIRST        
      TYPE (XPLEX)                :: LONRES,    LATRES
      TYPE (XPLEX)                :: ARRAY(IIPAR,JJPAR,LLPAR)
      !TYPE (XPLEX)                :: SCALE(IIPAR,JJPAR)
      TYPE (XPLEX)                :: SCALE
      TYPE (XPLEX)                :: DIAGb,     DIAGe
      CHARACTER(LEN=20)     :: MODELNAME 
      CHARACTER(LEN=40)     :: CATEGORY
      CHARACTER(LEN=40)     :: RESERVED
      CHARACTER(LEN=40)     :: UNIT

      !=================================================================
      ! WRITE_DIAG59 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND59 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      !SCALE     = FLOAT( CTOTH ) + TINY( 1d0 )
      !SCALE     = DBLE( GET_CT_DYN() ) + TINY( 1d0 )
      SCALE     = DBLE( GET_CT_DYN() ) + TINY( 1d0 )

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! debug
      !print*, ' LD59 = ', LD59
      !print*, ' some values of AD59 = ', AD59(20,20,:,:)

      ! Loop over ND03 diagnostic tracers
      DO M = 1, TMAX(44)

         ! Define quantities
         N        = TINDEX(59,M)
         CATEGORY = 'IJ-ugm3'
         UNIT ='ug/m3'

         IF ( N == 0 ) CYCLE 

         ! Apply scale factor
         DO L = 1, LD59
            ARRAY(:,:,L) = AD59(:,:,L,N) / SCALE
         ENDDO

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD59,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD59) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG59

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG59
!
!******************************************************************************
!  Subroutine INIT_DIAG59 allocates all module arrays (bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LSOA

#     include "CMN_SIZE"    ! Size parameters 

      ! Local variables
      INTEGER              :: AS
      
      !=================================================================
      ! INIT_DIAG42 begins here!
      !=================================================================

      ! Turn off ND59 if NH3 tracers are not used
!      IF ( .not. LNH3 ) THEN
!         ND59 = 0
!         RETURN
!      ENDIF

!      ! debug
!      print*, ' check ND59 = ', ND59 

      ! Exit if ND59 is turned off
      IF ( ND59 == 0 ) RETURN

      ! Number of levels to save for this diagnostic
      LD59 = MIN( ND59, LLPAR )

!     print*, ' check LD59 = ', LD59 

      ! 2-D array ("LFLASH-$")
      ALLOCATE( AD59( IIPAR, JJPAR, LD59, PD59 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59' )

      ! Zero arrays
      CALL ZERO_DIAG59

      ! Return to calling program
      END SUBROUTINE INIT_DIAG59

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG59
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG59 deallocates all module arrays (bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG59 begins here!
      !=================================================================
      IF ( ALLOCATED( AD59 ) ) DEALLOCATE( AD59 ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG59

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG59_MOD
