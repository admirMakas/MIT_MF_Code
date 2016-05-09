      MODULE ADJ_EMISARRAY_MOD

!******************************************************************************
! 
!******************************************************************************
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "carbon_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE              :: REMIS_ADJ,  DEPSAV_ADJ
      PRIVATE              :: EMISRN_ADJ, EMISRRN_ADJ

      ! PRIVATE module routines 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      ! Allocatable variables
      REAL*8, ALLOCATABLE  :: RRATE_ADJ(:,:)
      REAL*8, ALLOCATABLE  :: TAREA_ADJ(:,:)
      REAL*8, ALLOCATABLE  :: ERADIUS_ADJ(:,:)
      REAL*8, ALLOCATABLE  :: REMIS_ADJ(:,:)
      REAL*8, ALLOCATABLE  :: DEPSAV_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRRN_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRN_ADJ(:,:,:)
      !REAL*8, ALLOCATABLE  :: GEMISNOX_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: EMIS_LI_NOx_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: GEMISNOX2_ADJ(:,:)
      REAL*8, ALLOCATABLE  :: BURNEMIS_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: BIOFUEL_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRR_ADJ(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRRB_ADJ(:,:,:)
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS


      SUBROUTINE INIT_ADJ_ANTHROEMS
!     
!******************************************************************************
!  Subroutine INIT_ADJ_ANTHROEMS initializes all module arrays (dkh, 06/01/06)  
!
!  NOTES:
!  (1 ) Add ADJ_GEMISNOX.  (dkh, 02/10/07) 
!  (2 ) Add ADJ_GEMISNOX2. (dkh, 03/20/08) 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR, ERROR_STOP
      USE DRYDEP_MOD, ONLY : NUMDEP
      USE BIOMASS_MOD, ONLY : NBIOTRCE 
      USE BIOFUEL_MOD, ONLY : NBFTRACE
      
#     include "CMN_SIZE" ! Size parameters
#     include "comode.h" ! ITLOOP


      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_ADJ_ANTHROEMS begins here!
      !=================================================================
      
      ! Return if we already allocated arrays
      IF ( IS_INIT ) RETURN
     
      ! Check to make sure that NOXLEVEL and NOXEXTENT are both 2
      IF ( ( NOXLEVELS /= 2 ) .or. ( NOXEXTENT /= 2 ) ) 
     &     CALL ERROR_STOP( 'Invalid NOXLEVELS', 
     &   ' INIT_ADJ_ANTHROEMS in adj_anthroems_mod.f' )                                    
 
      ALLOCATE( _DEPSAV_ADJ( IIPAR, JJPAR, NUMDEP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DEPSAV_ADJ' )
      DEPSAV_ADJ = 0d0

      ALLOCATE( RRATE_ADJ( ITLOOP, NMTRATE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RRATE_ADJ' )
      RRATE_ADJ = 0d0

      ALLOCATE( TAREA_ADJ( ITLOOP, NDUST + NAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAREA_ADJ' )
      TAREA_ADJ = 0d0

      ALLOCATE( ERADIUS_ADJ( ITLOOP, NDUST + NAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERADIUS_ADJ' )
      ERADIUS_ADJ = 0d0

      ALLOCATE( REMIS_ADJ( ITLOOP, NEMIS(NCSURBAN) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REMIS_ADJ' )
      REMIS_ADJ = 0d0

      ALLOCATE( EMISRRN_ADJ( IIPAR, JJPAR, NOXEXTENT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISRRN_ADJ' )
      EMISRRN_ADJ = 0d0

      ALLOCATE( EMISRN_ADJ( IGLOB, JGLOB, NOXEXTENT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISRN_ADJ' )
      EMISRN_ADJ = 0d0

      !ALLOCATE( GEMISNOX_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEMISNOX_ADJ' )
      !GEMISNOX_ADJ = 0d0

      ALLOCATE( EMIS_LI_NOx_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_LI_NOx_ADJ' )
      EMIS_LI_NOx_ADJ = 0d0
      
      ALLOCATE( GEMISNOX2_ADJ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEMISNOX2_ADJ' )
      GEMISNOX2_ADJ = 0d0

      ALLOCATE( BIOFUEL_ADJ( IIPAR, JJPAR, NBFTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOFUEL_ADJ' )
      BIOFUEL_ADJ = 0d0

      ALLOCATE( BURNEMIS_ADJ( IIPAR, JJPAR, NBIOTRCE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BURNEMIS_ADJ' )
      BURNEMIS_ADJ = 0d0

      ALLOCATE( EMISRR_ADJ( IIPAR, JJPAR, 2:NEMPARA+NEMPARB ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISRR_ADJ' )
      EMISRR_ADJ = 0d0

      ALLOCATE( EMISRRB_ADJ( IIPAR, JJPAR, 2:NEMPARA+NEMPARB ), STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISRRB_ADJ' )
      EMISRRB_ADJ = 0d0

      ! Reset IS_INIT
      IS_INIT = .TRUE. 

      ! Return to calling progam 
      END SUBROUTINE INIT_ADJ_ANTHROEMS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_ADJ_ANTHROEMS
!     
!******************************************************************************
!  Subroutine CLEANUP_ADJ_ANTHROEMS deallocates all module arrays  
!   (dkh, 06/01/06)  
!
!  NOTES:
!  (1 ) Add ADJ_GEMISNOX.  (dkh, 02/10/07) 
!  (2 ) Add ADJ_GEMISNOX2. (dkh, 03/20/08) 
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_ADJ_ANTHROEMS begins here!
      !=================================================================
      IF ( ALLOCATED( RRATE_ADJ     ) ) DEALLOCATE( RRATE_ADJ     )
      IF ( ALLOCATED( TAREA_ADJ     ) ) DEALLOCATE( TAREA_ADJ     )
      IF ( ALLOCATED( ERADIUS_ADJ   ) ) DEALLOCATE( ERADIUS_ADJ   )
      IF ( ALLOCATED( REMIS_ADJ     ) ) DEALLOCATE( REMIS_ADJ     )
      IF ( ALLOCATED( DEPSAV_ADJ    ) ) DEALLOCATE( DEPSAV_ADJ    )
      IF ( ALLOCATED( EMISRRN_ADJ   ) ) DEALLOCATE( EMISRRN_ADJ   )
      IF ( ALLOCATED( EMISRN_ADJ    ) ) DEALLOCATE( EMISRN_ADJ    )
      !IF ( ALLOCATED( GEMISNOX_ADJ  ) ) DEALLOCATE( GEMISNOX_ADJ  )
      IF ( ALLOCATED( EMIS_LI_NOx_ADJ  ) )
     &                    DEALLOCATE( EMIS_LI_NOX_ADJ  )
      IF ( ALLOCATED( GEMISNOX2_ADJ ) ) DEALLOCATE( GEMISNOX2_ADJ )
      IF ( ALLOCATED( BIOFUEL_ADJ   ) ) DEALLOCATE( BIOFUEL_ADJ   )
      IF ( ALLOCATED( BURNEMIS_ADJ  ) ) DEALLOCATE( BURNEMIS_ADJ  )
      IF ( ALLOCATED( EMISRR_ADJ    ) ) DEALLOCATE( EMISRR_ADJ    )
      IF ( ALLOCATED( EMISRRB_ADJ   ) ) DEALLOCATE( EMISRRB_ADJ   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_ADJ_ANTHROEMS

!------------------------------------------------------------------------------
      END MODULE ADJ_EMISARRAY_MOD
