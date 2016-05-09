! $Id: comode_mod.f,v 1.8 2011/02/23 00:08:48 daven Exp $
      MODULE COMODE_MOD
!
!******************************************************************************
!  Module COMODE_MOD contains allocatable arrays for SMVGEAR that were
!  previously contained in common blocks in header file "comode.h".
!  (bmy, 8/31/00, 9/28/04)
!  
!  In case you were wondering, "comode" stands for:
!     "COMmon blocks: Ordinary Differential Equations"
!  
!  Module Variables:
!  ============================================================================
!  (1 ) ABSHUM     : array for absolute humidity [H2O molec/cm3]
!  (2 ) AIRDENS    : array for air density [molec/cm3]
!  (3 ) CSPEC      : array of chemical species concentration [molec/cm3]
!  (3a) CSPEC_FULL : array of chemical species for full potential troposphere
!  (4 ) CSUMA    : array for time of sunrise/sunset, measured from midnight [s]
!  (5 ) CSUMC    : array for temporary storage 
!  (6 ) ERADIUS  : array for aerosol or dust radii [cm]
!  (7 ) ERRMX2   : array for storing stiffness values 
!  (8 ) IXSAVE   : array of grid box longitude indices
!  (9 ) IYSAVE   : array of grid box latitude indices
!  (10) IZSAVE   : array of grid box altitude indices
!  (11) JLOP     : array of 1-D grid box indices
!  (12) PRESS3   : array for grid box pressure [mb]
!  (13) REMIS    : array for emissions from GEOS-CHEM [molec/cm3] 
!  (14) T3       : array for grid box temperature [K]
!  (15) TAREA    : array for surface area of aerosol or dust [cm2/cm3]
!  (16) VOLUME   : array for grid box volume [cm3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) INIT_COMODE    : allocates memory for arrays
!  (2 ) CLEANUP_COMODE : deallocates memory for arrays
!
!  GEOS-CHEM modules referenced by comode_mod.f
!  ============================================================================
!  (1 ) error_mod.f    : Module containing NaN and other error check routines
!
!  NOTES:
!  (1 ) Now zero CSPEC after allocating memory (bmy, 9/8/00)
!  (2 ) Now declare more SMVGEAR arrays allocatable (bmy, 10/19/00)
!  (3 ) Updated comments (bmy, 9/4/01)
!  (4 ) Now make ERADIUS, TAREA 2-D arrays, for het chem (bmy, 11/15/01)
!  (5 ) DARSFCA is now obsolete, remove it.  Now allocate ERADIUS and
!        TAREA arrays to be of size (ITLOOP,NDUST+NAER).  (rvm, bmy, 2/27/02)
!  (5 ) Removed obsolete code from 2/02 (bmy, 4/15/02)
!  (6 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (7 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (8 ) Now add CSUMA, CSUMC, ERRMX2 arrays for SMVGEAR II (bmy, 7/18/03)
!  (9 ) Now also references "tracer_mod.f" (bmy, 9/28/04)
!  (10) Add WTAREA and WERADIUS variables. 
!       For SOA production from reactive uptake of dicarbonyls, 
!       archived WTAREA and WERADIUS should include dusts, 
!       but excludes BCPO and OCPO (tmf, ccc, 1/7/09)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      TYPE (XPLEX),  ALLOCATABLE :: ABSHUM(:) 
      TYPE (XPLEX),  ALLOCATABLE :: AIRDENS(:) 
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_FULL(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSUMA(:) 
      TYPE (XPLEX),  ALLOCATABLE :: CSUMC(:) 
      TYPE (XPLEX),  ALLOCATABLE :: ERADIUS(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: ERRMX2(:) 
      INTEGER, ALLOCATABLE :: IXSAVE(:)
      INTEGER, ALLOCATABLE :: IYSAVE(:)
      INTEGER, ALLOCATABLE :: IZSAVE(:)
      INTEGER, ALLOCATABLE :: JLOP(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: PRESS3(:)      
      TYPE (XPLEX),  ALLOCATABLE :: REMIS(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: T3(:)      
      TYPE (XPLEX),  ALLOCATABLE :: TAREA(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: VOLUME(:)      
      TYPE (XPLEX),  ALLOCATABLE :: WTAREA(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: WERADIUS(:,:)

      !/---------------------------------------------\!
      ! ADJ_GROUP: Adding more variables to be used   !
      ! for checkpointing and adjoint calculations    !
      !\---------------------------------------------/|
      TYPE (XPLEX),  ALLOCATABLE :: R_KPP(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_ADJ(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_FOR_KPP(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_ORIG(:,:)
      !TYPE (XPLEX),  ALLOCATABLE :: CSPEC_ADJ_FOR_KPP(:,:)   
      TYPE (XPLEX),  ALLOCATABLE :: HSAVE(:,:,:)
      ! Add CSPEC_PRIOR, CHK_CSPEC, CSPEC_FOR_KPP_ADJ (dkh, 06/11/09) 
      ! and O3_AFTER_CHEM (dkh, 06/12/09) 
      ! and NO2_AFTER_CHEM (dkh, 06/14/09) 
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_PRIOR(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CHK_CSPEC(:,:)
      !TYPE (XPLEX),  ALLOCATABLE :: CSPEC_FOR_KPP_ADJ(:,:)
 
      ! Replace these with CSPEC_AFTER_CHEM and 
      ! CSPEC_AFTER_CHEM_ADJ (dkh, 02/09/11) 
      !!TYPE (XPLEX),  ALLOCATABLE :: O3_AFTER_CHEM(:)
      !!TYPE (XPLEX),  ALLOCATABLE :: NO2_AFTER_CHEM(:)
      !!TYPE (XPLEX),  ALLOCATABLE :: NO2_AFTER_CHEM_ADJ(:)
      !!TYPE (XPLEX),  ALLOCATABLE :: CSPEC_ADJ_FORCE(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_AFTER_CHEM(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_AFTER_CHEM_ADJ(:,:)


      ! LVARTROP support for adj (dkh, 01/26/11)
      TYPE (XPLEX),  ALLOCATABLE :: CSPEC_FULL_PRIOR(:,:,:,:)
      INTEGER, ALLOCATABLE :: ISAVE_PRIOR(:,:)
      INTEGER              :: NTLOOP_PRIOR


      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
      
!------------------------------------------------------------------------------
      
      SUBROUTINE INIT_COMODE
!
!******************************************************************************
!  Subroutine INIT_COMODE allocates memory for allocatable arrays that were 
!  previously contained in common blocks in "comode.h". (bmy, 8/31/00, 9/28/04)
!
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Cosmetic chagnes (bmy, 2/27/03)
!  (3 ) Now allocate CSUMA, CSUMC, ERRMX2; cosmetic changes (bmy, 7/18/03)
!  (4 ) Now allocate certain arrays for offline aerosol sim (bmy, 9/28/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE TRACER_MOD, ONLY : ITS_AN_AEROSOL_SIM, ITS_A_FULLCHEM_SIM

#     include "CMN_SIZE"
#     include "comode.h" 

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_COMODE begins here!
      !=================================================================
      WRITE( 6, 100 )
 100  FORMAT( '     - INIT_COMODE: Allocating arrays for SMVGEAR...' )

      !----------------------------------
      ! FULL CHEMISTRY SIMULATION
      !----------------------------------
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
      
         ALLOCATE( ABSHUM( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ABSHUM' )
         ABSHUM = 0d0
      
         ALLOCATE( AIRDENS( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDENS' )
         AIRDENS = 0d0      

         ALLOCATE( CSPEC( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC' )
         CSPEC = 0d0

         ALLOCATE( CSPEC_FULL( ILONG, ILAT, IPVERT, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_FULL' )
         CSPEC_FULL = 0d0

         ALLOCATE( CSUMA( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSUMA' )
         CSUMA = 0d0
      
         ALLOCATE( CSUMC( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSUMC' )
         CSUMC = 0d0

         ALLOCATE( ERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERADIUS' )
         ERADIUS = 0d0      

         ALLOCATE( ERRMX2( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERRMX2' )
         ERRMX2 = 0d0
           
         ALLOCATE( IXSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IXSAVE' )
         IXSAVE = 0
      
         ALLOCATE( IYSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IYSAVE' )
         IYSAVE = 0
      
         ALLOCATE( IZSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IZSAVE' )
         IZSAVE = 0
      
         ALLOCATE( JLOP( ILONG, ILAT, IPVERT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP' )
         JLOP = 0
      
         ALLOCATE( PRESS3( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRESS3' )
         PRESS3 = 0d0
      
         ALLOCATE( REMIS( ITLOOP, MAXGL3 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'REMIS' )
         REMIS = 0d0
      
         ALLOCATE( T3( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'T3' )
         T3 = 0d0

         ALLOCATE( TAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAREA' )
         TAREA = 0d0      
      
         ALLOCATE( VOLUME( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'VOLUME' )
         VOLUME = 0d0

         ALLOCATE( WTAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'WTAREA' )
         WTAREA = 0d0      

         ALLOCATE( WERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'WERADIUS' )
         WERADIUS = 0d0    

         !/---------------------------------------------\!
         ! ADJ_GROUP: Adding more variables to be used   !
         ! for checkpointing and adjoint calculations    !
         !\---------------------------------------------/|
         ALLOCATE( R_KPP( ITLOOP, NMTRATE ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'R_KPP' )
         R_KPP = 0d0

         ALLOCATE( CSPEC_ADJ( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_ADJ' )
         CSPEC_ADJ = 0d0

         !ALLOCATE( CSPEC_ADJ_FOR_KPP( ITLOOP, IGAS ), STAT=AS )
         !IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_ADJ_FOR_KPP' )
         !CSPEC_ADJ_FOR_KPP = 0d0

         ALLOCATE( CSPEC_FOR_KPP( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_FOR_KPP' )
         CSPEC_FOR_KPP = 0d0

         ALLOCATE( CSPEC_ORIG( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_ORIG' )
         CSPEC_ORIG = 0d0

         ALLOCATE( HSAVE( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'HSAVE' )
         HSAVE = 0.d0

        ! LVARTROP support for adj (dkh, 01/26/11)
         ALLOCATE( CSPEC_FULL_PRIOR( ILONG, ILAT, IPVERT, IGAS ),
     &      STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_FULL' )
         CSPEC_FULL_PRIOR = 0d0

         ALLOCATE( ISAVE_PRIOR( ITLOOP, 3 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ISAVE_PRIOR' )
         ISAVE_PRIOR = 0


         ! Add CSPEC_PRIOR (dkh, 06/11/09) 
         ALLOCATE( CSPEC_PRIOR( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_PRIOR' )
         CSPEC_PRIOR = 0d0

         ! Add CSPEC_PRIOR (dkh, 06/11/09) 
         ALLOCATE( CHK_CSPEC( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_CSPEC' )
         CHK_CSPEC = 0d0

         ! Add CSPEC_FOR_KPP_ADJ (dkh, 06/11/09) 
         !ALLOCATE( CSPEC_FOR_KPP_ADJ( ITLOOP, IGAS ), STAT=AS )
         !IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_FOR_KPP_ADJ' )
         !CSPEC_FOR_KPP_ADJ = 0d0

         ! Now use CSPEC_AFTER_CHEM and CSPEC_AFTER_CHEM_ADJ
         ! (dkh, 02/09/11) 
         !!! Add O3_AFTER_CHEM (dkh, 06/12/09) 
         !!ALLOCATE( O3_AFTER_CHEM( ITLOOP ), STAT=AS )
         !!IF ( AS /= 0 ) CALL ALLOC_ERR( 'O3_AFTER_CHEM' )
         !!O3_AFTER_CHEM = 0d0
         !!
         !!! Add NO2_AFTER_CHEM (dkh, 06/14/09) 
         !!ALLOCATE( NO2_AFTER_CHEM( ITLOOP ), STAT=AS )
         !!IF ( AS /= 0 ) CALL ALLOC_ERR( 'NO2_AFTER_CHEM' )
         !!NO2_AFTER_CHEM = 0d0
         !! Add NO2_AFTER_CHEM (dkh, 06/14/09) 
         !!
         !! Add NO2_AFTER_CHEM_ADJ  (dkh, 07/31/09) 
         !!ALLOCATE( NO2_AFTER_CHEM_ADJ( ITLOOP ), STAT=AS )
         !!IF ( AS /= 0 ) CALL ALLOC_ERR( 'NO2_AFTER_CHEM_ADJ' )
         !!NO2_AFTER_CHEM_ADJ = 0d0 
         !!
         !!! Add CSPEC_ADJ_FORCE (dkh, 07/31/09) 
         !!ALLOCATE( CSPEC_ADJ_FORCE( ITLOOP, IGAS ), STAT=AS )
         !!IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC_ADJ_FORCE' )
         !!CSPEC_ADJ_FORCE = 0d0


      ENDIF

      !----------------------------------
      ! OFFLINE AEROSOL SIMULATION
      !----------------------------------
      IF ( ITS_AN_AEROSOL_SIM() ) THEN

         ALLOCATE( ABSHUM( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ABSHUM' )
         ABSHUM = 0d0
      
         ALLOCATE( AIRDENS( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDENS' )
         AIRDENS = 0d0      

         ALLOCATE( ERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERADIUS' )
         ERADIUS = 0d0      

         ALLOCATE( IXSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IXSAVE' )
         IXSAVE = 0
      
         ALLOCATE( IYSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IYSAVE' )
         IYSAVE = 0
      
         ALLOCATE( IZSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IZSAVE' )
         IZSAVE = 0
      
         ALLOCATE( JLOP( ILONG, ILAT, IPVERT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP' )
         JLOP = 0

         ALLOCATE( TAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAREA' )
         TAREA = 0d0      
         
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_COMODE

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_COMODE
!
!******************************************************************************
!  Subroutine CLEANUP_COMODE deallocates memory from allocatable arrays 
!  that were previously contained in common blocks in "comode.h" 
!  (bmy, 8/31/00, 7/18/03)
!
!  NOTES:
!  (1 ) Now deallocate CSPEC, CSUMA, ERRMX2; cosmetic changes (bmy, 7/18/03)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_COMODE begins here!
      !=================================================================
      IF ( ALLOCATED( ABSHUM     ) ) DEALLOCATE( ABSHUM  )
      IF ( ALLOCATED( AIRDENS    ) ) DEALLOCATE( AIRDENS )
      IF ( ALLOCATED( CSPEC      ) ) DEALLOCATE( CSPEC   )
      IF ( ALLOCATED( CSPEC_FULL ) ) DEALLOCATE( CSPEC_FULL)
      IF ( ALLOCATED( CSUMA      ) ) DEALLOCATE( CSUMA   )
      IF ( ALLOCATED( CSUMC      ) ) DEALLOCATE( CSUMC   )
      IF ( ALLOCATED( ERADIUS    ) ) DEALLOCATE( ERADIUS )
      IF ( ALLOCATED( ERRMX2     ) ) DEALLOCATE( ERRMX2  )
      IF ( ALLOCATED( IXSAVE     ) ) DEALLOCATE( IXSAVE  )
      IF ( ALLOCATED( IYSAVE     ) ) DEALLOCATE( IYSAVE  )
      IF ( ALLOCATED( IZSAVE     ) ) DEALLOCATE( IZSAVE  )
      IF ( ALLOCATED( JLOP       ) ) DEALLOCATE( JLOP    )
      IF ( ALLOCATED( PRESS3     ) ) DEALLOCATE( PRESS3  )     
      IF ( ALLOCATED( REMIS      ) ) DEALLOCATE( REMIS   )
      IF ( ALLOCATED( T3         ) ) DEALLOCATE( T3      )     
      IF ( ALLOCATED( TAREA      ) ) DEALLOCATE( TAREA   )
      IF ( ALLOCATED( VOLUME     ) ) DEALLOCATE( VOLUME  )  
      IF ( ALLOCATED( WTAREA     ) ) DEALLOCATE( WTAREA  )
      IF ( ALLOCATED( WERADIUS   ) ) DEALLOCATE( WERADIUS )

      !/---------------------------------------------\!
      ! ADJ_GROUP: Adding more variables to be used   !
      ! for checkpointing and adjoint calculations    !
      !\---------------------------------------------/|
      IF ( ALLOCATED( CSPEC_ADJ  ) )    DEALLOCATE( CSPEC_ADJ         )
      IF ( ALLOCATED( CSPEC_FOR_KPP ) ) DEALLOCATE( CSPEC_FOR_KPP     )
      IF ( ALLOCATED( CSPEC_ORIG ) )    DEALLOCATE( CSPEC_ORIG        )
      IF ( ALLOCATED( R_KPP       ) )   DEALLOCATE( R_KPP             ) 
      !IF ( ALLOCATED( CSPEC_ADJ_FOR_KPP ) ) 
      !&                                  DEALLOCATE( CSPEC_ADJ_FOR_KPP )  
      IF ( ALLOCATED( HSAVE          ) ) DEALLOCATE( HSAVE       )   
      ! add CSPEC_PRIOR, CHK_CSPEC, CSPEC_FOR_KPP_ADJ (dkh, 06/11/09) 
      ! and O3_AFTER_CHEM (dkh, 06/12/09) 
      ! and NO2_AFTER_CHEM 
      !IF ( ALLOCATED( CSPEC_FOR_KPP_ADJ ) ) 
      !&    DEALLOCATE( CSPEC_FOR_KPP_ADJ )
      IF ( ALLOCATED( CSPEC_PRIOR       ) ) DEALLOCATE( CSPEC_PRIOR    )
      IF ( ALLOCATED( CHK_CSPEC         ) ) DEALLOCATE( CHK_CSPEC      )
      !IF ( ALLOCATED( O3_AFTER_CHEM     ) ) DEALLOCATE( O3_AFTER_CHEM  )
      !IF ( ALLOCATED( NO2_AFTER_CHEM    ) ) DEALLOCATE( NO2_AFTER_CHEM )
      !IF ( ALLOCATED( NO2_AFTER_CHEM_ADJ) ) 
      !&    DEALLOCATE( NO2_AFTER_CHEM_ADJ )
      !IF ( ALLOCATED( CSPEC_ADJ_FORCE   ) ) DEALLOCATE( CSPEC_ADJ_FORCE)
      IF ( ALLOCATED( CSPEC_AFTER_CHEM ) ) DEALLOCATE( CSPEC_AFTER_CHEM)
      IF ( ALLOCATED( CSPEC_AFTER_CHEM_ADJ ) ) 
     &    DEALLOCATE( CSPEC_AFTER_CHEM_ADJ )

      ! LVARTROP support for adj (dkh, 01/26/11)
      IF ( ALLOCATED( CSPEC_FULL_PRIOR  ) ) DEALLOCATE(CSPEC_FULL_PRIOR)
      IF ( ALLOCATED( ISAVE_PRIOR  ) )      DEALLOCATE(ISAVE_PRIOR)

      

      ! Return to calling program
      END SUBROUTINE CLEANUP_COMODE

!------------------------------------------------------------------------------

      END MODULE COMODE_MOD

