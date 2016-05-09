! $Id: diag42_mod.f,v 1.1.1.1 2009/06/09 21:51:53 daven Exp $
      MODULE DIAG42_MOD
!
!******************************************************************************
!  Module DIAG42_MOD contains arrays and routines for archiving the ND42
!  diagnostic -- secondary organic aerosols [ug/m3]. (dkh,bmy,5/22/06,3/29/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD42 (TYPE (XPLEX))  : Array for SOA concentrations [ug/m3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) DIAG42         : Archives quantities for diagnostic
!  (2 ) ZERO_DIAG42    : Sets all module arrays to zero
!  (3 ) WRITE_DIAG42   : Writes data in module arrays to bpch file
!  (4 ) INIT_DIAG42    : Allocates all module arrays
!  (5 ) CLEANUP_DIAG42 : Deallocates all module arrays
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
!  (2 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  (3 ) Add diagnostics for SOAG and SOAM (tmf, 1/7/09)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag42_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND42, LD42

      ! Parameters
      INTEGER, PARAMETER   :: PD42 = 14

      ! Arrays
      TYPE (XPLEX),  ALLOCATABLE :: AD42(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE DIAG42
!
!******************************************************************************
!  Subroutine DIAG42 archives SOA concentrations [ug/m3] for the ND42
!  diagnostic. (dkh, bmy, 5/22/06, 3/29/07)
!
!  NOTES:
!  (1 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AIRVOL, T
      !USE DIAG_MOD,     ONLY : LTOTH
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSOA1, IDTSOA2, IDTSOA3, IDTSOA4
      USE TRACERID_MOD, ONLY : IDTOCPI, IDTOCPO
      USE TRACERID_MOD, ONLY : IDTSOAG, IDTSOAM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! NDxx flags

      ! Local variables
      INTEGER               :: I,      J,    L
      TYPE (XPLEX)                :: FACTOR, PRES

      ! Factor for computing standard volume
      TYPE (XPLEX),PARAMETER::STD_VOL_FAC=xplex(1013.25d0/273.15d0,0d0)
     
      !================================================================= 
      ! DIAG42 begins here! 
      !================================================================= 

      ! Error check
      IF ( IDTSOA1 == 0 ) RETURN
      IF ( IDTSOA2 == 0 ) RETURN
      IF ( IDTSOA3 == 0 ) RETURN
      IF ( IDTSOA4 == 0 ) RETURN
      IF ( IDTOCPO == 0 ) RETURN
      IF ( IDTOCPI == 0 ) RETURN

      ! Loop over grid boxes     
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FACTOR, PRES )  
      DO L = 1, LD42  
      DO J = 1, JJPAR 
      DO I = 1, IIPAR

         ! Conversion factor from [kg] --> [ug/m3]
         ! (LTOTH=1 if between OTH_HR1 and OTH_HR2, LTOTH=0 otherwise)
         !FACTOR        = 1d9 / AIRVOL(I,J,L) * LTOTH(I,J) 

         ! Conversion factor from [kg] --> [ug/m3]
         FACTOR        = 1d9 / AIRVOL(I,J,L)

         ! SOA1 [ug/m3]
         AD42(I,J,L,1) = AD42(I,J,L,1)        + 
     &                   ( STT(I,J,L,IDTSOA1) * FACTOR )
 
         ! SOA2 [ug/m3]
         AD42(I,J,L,2) = AD42(I,J,L,2)        + 
     &                   ( STT(I,J,L,IDTSOA2) * FACTOR )

         ! SOA3 [ug/m3]
         AD42(I,J,L,3) = AD42(I,J,L,3)        + 
     &                   ( STT(I,J,L,IDTSOA3) * FACTOR )

         ! SOA4 [ug/m3]
         AD42(I,J,L,4) = AD42(I,J,L,4)        + 
     &                   ( STT(I,J,L,IDTSOA4) * FACTOR )

         ! Sum of original 3 SOA types [ug/m3]
         AD42(I,J,L,5) = AD42(I,J,L,5)        + 
     &                   ( STT(I,J,L,IDTSOA1) + 
     &                     STT(I,J,L,IDTSOA2) +  
     &                     STT(I,J,L,IDTSOA3) ) * FACTOR

         ! Sum of SOA1 to SOA4
         AD42(I,J,L,6) = AD42(I,J,L,6)        + 
     &                   ( STT(I,J,L,IDTSOA1) + 
     &                     STT(I,J,L,IDTSOA2) + 
     &                     STT(I,J,L,IDTSOA3) + 
     &                     STT(I,J,L,IDTSOA4) ) * FACTOR

         ! Sum of primary OC + SOA1 to SOA4 [ug C/m3] 
         ! Use higher ratio (2.1) of molecular weight of
         ! organic mass per carbon mass accounting for non-carbon
         ! components attached to OC [Turpin and Lim, 2001] 
         AD42(I,J,L,7) = AD42(I,J,L,7)          +
     &                   ( ( STT(I,J,L,IDTSOA1) + 
     &                       STT(I,J,L,IDTSOA2) + 
     &                       STT(I,J,L,IDTSOA3) + 
     &                       STT(I,J,L,IDTSOA4) )   / 2.1d0
     &                   + ( STT(I,J,L,IDTOCPO) + 
     &                       STT(I,J,L,IDTOCPI) ) ) * FACTOR

         ! Sum of PRIMARY OC + SOA1 to SOA4 [ug C/m3] at STP
         PRES          = GET_PCENTER( I, J, L )
         AD42(I,J,L,8) = AD42(I,J,L,7) * STD_VOL_FAC * T(I,J,L) / PRES

!--------------------------------------------------------
! Additional diagnostics for SOAG, SOAM (tmf, 12/8/07) 
! Assume SOAG mass = GLYX mass, SOAM mass = MGLY mass
! Test if SOAG and SOAM are simulated (ccc, 12/18/08)
!--------------------------------------------------------
         IF ( IDTSOAG /= 0 .AND. IDTSOAM /=0 ) THEN
            ! SOAG [ug total mass /m3]
            AD42(I,J,L,9) = AD42(I,J,L,9)        + 
     &                      ( STT(I,J,L,IDTSOAG) * 1.d0 * FACTOR )

            ! SOAM [ug total mass /m3]
            AD42(I,J,L,10) = AD42(I,J,L,10)        + 
     &                      ( STT(I,J,L,IDTSOAM) * 1.d0 * FACTOR )


            ! Sum of SOA1 to SOA4, SOAG, SOAM (tmf, 1/31/07)
            AD42(I,J,L,11) = AD42(I,J,L,11)        + 
     &                      ( STT(I,J,L,IDTSOA1) + 
     &                        STT(I,J,L,IDTSOA2) + 
     &                        STT(I,J,L,IDTSOA3) + 
     &                        STT(I,J,L,IDTSOA4)  +
     &                      ( STT(I,J,L,IDTSOAG) * 1.d0 ) +
     &                      ( STT(I,J,L,IDTSOAM) * 1.d0 )) * FACTOR 

            ! Sum of SOA1 to SOA4, SOAG, SOAM in carbon (tmf, 1/31/07) 
            ! Except SOAG is 0.41 carbon, SOAM is 0.5 carbon
            AD42(I,J,L,12) = AD42(I,J,L,12)          +
     &                      ( ( STT(I,J,L,IDTSOA1) + 
     &                          STT(I,J,L,IDTSOA2) + 
     &                          STT(I,J,L,IDTSOA3) + 
     &                          STT(I,J,L,IDTSOA4) )   / 2.1d0 +
     &                         ( STT(I,J,L,IDTSOAG) * 0.41D0 ) +
     &                         ( STT(I,J,L,IDTSOAM) * 0.50D0 ) +
     &                        ( STT(I,J,L,IDTOCPO) + 
     &                          STT(I,J,L,IDTOCPI) ) ) * FACTOR

            ! Sum of SOA1 to SOA4, SOAG, SOAM at STP [ug/sm3 STP] (tmf, 1/31/07)  
            PRES          = GET_PCENTER( I, J, L )
            AD42(I,J,L,13) = AD42(I,J,L,11) * STD_VOL_FAC * T(I,J,L) 
     &                       / PRES

            ! Sum of all OC [ug C/sm3] at STP (including SOAG, SOAM)
            AD42(I,J,L,14) = AD42(I,J,L,12) * STD_VOL_FAC * T(I,J,L) 
     &                       / PRES
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO 

      ! Return to calling program
      END SUBROUTINE DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG42
!
!******************************************************************************
!  Subroutine ZERO_DIAG42 zeroes the ND03 diagnostic arrays. 
!  (dkh, bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! ZERO_DIAG42 begins here!
      !=================================================================

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Zero arrays
      AD42(:,:,:,:) = 0d0

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG42
!
!******************************************************************************
!  Subroutine WRITE_DIAG03 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 5/22/06, 9/5/06)
!
!   # : Field    : Description                 : Units    : Scale factor
!  -----------------------------------------------------------------------
!  (1 ) IJ-SOA-$ : SOA1                        : ug/m3    : SCALE_OTH
!  (2 ) IJ-SOA-$ : SOA2                        : ug/m3    : SCALE_OTH
!  (3 ) IJ-SOA-$ : SOA3                        : ug/m3    : SCALE_OTH
!  (4 ) IJ-SOA-$ : SOA4                        : ug/m3    : SCALE_OTH
!  (5 ) IJ-SOA-$ : SOA1 + SOA2 + SOA3          : ug/m3    : SCALE_OTH
!  (6 ) IJ-SOA-$ : SOA1 + SOA2 + SOA3 + SOA4   : ug/m3    : SCALE_OTH
!  (7 ) IJ-SOA-$ : Sum of all Org Carbon       : ug C/m3  : SCALE_OTH
!  (8 ) IJ-SOA-$ : Sum of all Org Carbon @ STP : ug C/sm3 : SCALE_OTH
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
      ! WRITE_DIAG42 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND42 == 0 ) RETURN

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

      ! Loop over ND03 diagnostic tracers
      DO M = 1, TMAX(42)

         ! Define quantities
         N        = TINDEX(42,M)
         CATEGORY = 'IJ-SOA-$'

         ! Pick proper unit
         SELECT CASE ( N )
            CASE( 7 )
               UNIT = 'ug C/m3'
            CASE( 8 )
               UNIT = 'ug C/sm3'
            CASE( 12 )
               UNIT = 'ug C/m3'
            CASE( 13 )
               UNIT = 'ug/sm3'
            CASE( 14 )
               UNIT = 'ug C/sm3'
            CASE DEFAULT
               UNIT = 'ug/m3'
         END SELECT

         ! Apply scale factor
         DO L = 1, LD42
            !ARRAY(:,:,L) = AD42(:,:,L,N) / SCALE(:,:)
            ARRAY(:,:,L) = AD42(:,:,L,N) / SCALE
         ENDDO

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD42,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD42) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG42
!
!******************************************************************************
!  Subroutine INIT_DIAG42 allocates all module arrays (bmy, 5/22/06)
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

      ! Turn off ND42 if SOA tracers are not used
      IF ( .not. LSOA ) THEN
         ND42 = 0
         RETURN
      ENDIF

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Number of levels to save for this diagnostic
      LD42 = MIN( ND42, LLPAR )

      ! 2-D array ("LFLASH-$")
      ALLOCATE( AD42( IIPAR, JJPAR, LD42, PD42 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD42' )

      ! Zero arrays
      CALL ZERO_DIAG42

      ! Return to calling program
      END SUBROUTINE INIT_DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG42
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG42 deallocates all module arrays (bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG42 begins here!
      !=================================================================
      IF ( ALLOCATED( AD42 ) ) DEALLOCATE( AD42 ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG42

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG42_MOD
