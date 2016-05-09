!$Id: population_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
       MODULE POPULATION_MOD
!
!******************************************************************************
!  Module POPULATION_MOD contains code for incorporating  population weighting 
!  into cost functions / exposure metrics. Population data taken from:
!
!    Center for International Earth Science Information Network (CIESIN), 
!    Columbia University; and Centro Internacional de Agricultura Tropical 
!    (CIAT). 2005. Gridded Population of the World, Version 3 (GPWv3): 
!    Population Count Grid. Palisades, NY: Socioeconomic Data and Applications 
!    Center (SEDAC), Columbia University. 
!    Available at http://sedac.ciesin.columbia.edu/gpw. 2/11/2012.
!
!  Steven Vogel, jk, dkh, 02/04/2012, adj32_024
!     
!  Module Variables:
!  ============================================================================
!  (1 ) POP_REDUCED   (TYPE (XPLEX)) : Array of census population
!
!  Module Routines:
!  ===========================================================================
!  (1 ) POP_WEIGHT_COST        : Computes population weighted cost function
!  (2 ) READ_IN_POPULATION     : Reads in gridded population data file
!  (3 ) INIT_POPULATOIN_MOD    : Allocates & zeroes module arrays
!  (4 ) CLEANUP_POPULATION_MOD : Deallocates module arrays
!
!  NOTES:
!
!*****************************************************************************
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      PUBLIC

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      TYPE (XPLEX), ALLOCATABLE :: POP_REDUCED(:,:)
      LOGICAL             :: WGT_MORT 

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE POP_WEIGHT_COST
!
!******************************************************************************
! This subroutine based on CALC_ADJ_FORCE_FOR_SENSE in geos_chem_adj_mod.f
! Calculates population weighted cost function when called in 
! geos_chem_adj_mod.f
!
!******************************************************************************
!
      ! References to F90 modules 
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC, COST_FUNC
      USE ADJ_ARRAYS_MOD,       ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : GET_WEIGHT
      USE ADJ_ARRAYS_MOD,       ONLY : NSPAN
      USE ADJ_ARRAYS_MOD,       ONLY : OBS_THIS_TRACER
      USE CHECKPT_MOD,          ONLY : CHK_STT
      USE DAO_MOD,              ONLY : AIRVOL, AD
      USE LOGICAL_MOD,          ONLY : LPRT
      USE TRACER_MOD,           ONLY : N_TRACERS

      ! Header files
#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      TYPE (XPLEX)          :: ADJ_FORCE(IIPAR,JJPAR,LLPAR,N_TRACERS)
      INTEGER             :: I, J, L, N
      LOGICAL, SAVE       :: FIRST = .TRUE.

      TYPE (XPLEX), DIMENSION(IIPAR,JJPAR,LLPAR,N_TRACERS) :: COST_NUMM
      TYPE (XPLEX), DIMENSION(IIPAR,JJPAR)                 :: DENOMM_POP
      TYPE (XPLEX), DIMENSION(IIPAR,JJPAR,LLPAR)           :: DENOMM_VOL
      TYPE (XPLEX)              :: NEW_COST_SCALAR
      TYPE (XPLEX)              :: POP_TOT
      TYPE (XPLEX)              :: VOL_TOT
      TYPE (XPLEX)              :: FACTORR

      !=================================================================
      ! POP_WEIGHT_COST begins here!
      !=================================================================

      ! Get population data     
      IF ( FIRST ) THEN
         CALL INIT_POPULATION_MOD 
         CALL READ_IN_POPULATION
         FIRST = .FALSE.
      ENDIF

      IF ( LPRT ) THEN  
         print*, 'SEV DEBUG = ', maxval(POP_REDUCED)
         print*, 'SEV DEBUG = ', minval(POP_REDUCED)
      ENDIF 

      ! Initialze cost fnc variables
      NEW_COST_SCALAR = 0d0
      COST_NUMM       = 0d0
      DENOMM_POP      = 0d0
      DENOMM_VOL      = 0d0


!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,    L,   N)
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Determine the contribution to the cost function in each grid cell
         ! from each species
         COST_NUMM(I,J,L,N)   = GET_WEIGHT(I,J,L,N)
     &                        * CHK_STT(I,J,L,N)
     &                        * POP_REDUCED(I,J)

         ! Set denominator population and volume 
         IF ( GET_WEIGHT(I,J,L,N) > 0 ) THEN
            DENOMM_POP(I,J)  =
     &                        POP_REDUCED(I,J)

            DENOMM_VOL(I,J,L)  =
     &                        AIRVOL(I,J,L)
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      POP_TOT = SUM(DENOMM_POP)
      VOL_TOT = SUM(DENOMM_VOL)

      IF ( LPRT ) THEN  
         print*, 'SEV DEBUG TOTAL VOLUME = ', VOL_TOT
         print*, 'SEV DEBUG TOTAL POP = ', POP_TOT
#if   defined ( GRID2x25 )
         print*, 'SEV DEBUG DENOM POP= ', DENOMM_POP(117,64)
#endif 
      ENDIF 

      FACTORR  = 1d9 / ( POP_TOT * VOL_TOT * NSPAN )

      ! calculate total average population weighted ug/m3
      ! accounting for the contribution at this hour to the 
      ! final average 
      DO N = 1, N_TRACERS

         IF ( OBS_THIS_TRACER(N) ) THEN
            NEW_COST_SCALAR = NEW_COST_SCALAR
     &                      + SUM(COST_NUMM(:,:,:,N)) * FACTORR

         ENDIF

      ENDDO

      ! Update cost function
      COST_FUNC = COST_FUNC + NEW_COST_SCALAR

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,    L,   N)
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR


         ! Force the adjoint variables x with dJ/dx=1
         ADJ_FORCE(I,J,L,N)   = GET_WEIGHT(I,J,L,N)
     &                        * POP_REDUCED(I,J)
     &                        * FACTORR


         STT_ADJ(I,J,L,N)   = STT_ADJ(I,J,L,N) + ADJ_FORCE(I,J,L,N)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
       END SUBROUTINE POP_WEIGHT_COST

!------------------------------------------------------------------------------
 
       SUBROUTINE READ_IN_POPULATION
!        
!******************************************************************************
!  Subroutine READ_IN_POPULATION reads in gridded population data. 
!  by Steven Vogel, based on code from Jamin Koo (dkh, 02/13/12, adj32_024) 
!
!  NOTES:
!
!******************************************************************************
!

      ! References to F90 modules 
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LPRT
#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      CHARACTER(LEN=255)     :: FNAME 
      INTEGER                :: IOS, IOS2

      !=================================================================
      ! READ_IN_POPULATION begins here!
      !=================================================================

      ! Generate population data filename 

      IF ( WGT_MORT ) THEN 
         FNAME = TRIM( DATA_DIR ) // 'population_201202/' // 
     &             'world_mort_premult.' // GET_RES_EXT()
      ELSE 
         FNAME = TRIM( DATA_DIR ) // 'population_201202/' // 
     &             'world_population.' // GET_RES_EXT()

      ENDIF

      ! Read the population from ascii file.
      WRITE( 6, '(a)' ) '  Reading in population from ', FNAME

      OPEN( UNIT=11, FILE=FNAME, STATUS='OLD', IOSTAT=IOS)

      IF ( IOS /= 0 ) THEN
         PRINT *, 'ERROR opening weight'
      ELSE
         READ( UNIT=11, FMT=*, IOSTAT=IOS2 ) POP_REDUCED
         IF ( IOS2 < 0 ) THEN 
            WRITE( 6, '(a)' ) '  Unexpected End of File encountered  '
         ELSE IF ( IOS > 0 ) THEN
            WRITE( 6, '(a)' ) '  Error occurred reading pop data!  '
         ENDIF
      ENDIF

      CLOSE( UNIT=11 )
      PRINT *, 'POP ', POP_REDUCED(52,25)
      PRINT *, 'POP ', POP_REDUCED(35,27)
      IF ( LPRT ) THEN 
         PRINT *, 'sum of population',         sum(POP_REDUCED)
         PRINT *, 'Population Grid Test Max',  maxval(POP_REDUCED)
         PRINT *, 'Population Grid Test Min',  minval(POP_REDUCED)
         PRINT *, 'Population Grid Test Size', size(POP_REDUCED)
      ENDIF 

      ! Return to calling program
      END SUBROUTINE READ_IN_POPULATION      
!------------------------------------------------------------------------------

      SUBROUTINE INIT_POPULATION_MOD
!
!******************************************************************************
!  Subroutine INIT_POPULATION_MOD initializes and zeros all allocatable arrays
!  declared in "population_mod.f"
!     
!  NOTES:
!     
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ALLOC_ERR

#     include "CMN_SIZE"      ! Size parameters

      ! local variables 
      INTEGER                :: AS
      
      !=================================================================
      ! INIT_POPULATION_MOD
      !=================================================================

      ALLOCATE( POP_REDUCED(IIPAR,JJPAR) , STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'POP_REDUCED' ) 
      POP_REDUCED = 0d0
               
      ! Return to calling program
      END SUBROUTINE INIT_POPULATION_MOD
                  
!-----------------------------------------------------------------------------
      SUBROUTINE CLEANUP_POPULATION_MOD
!              
!******************************************************************************
!  Subroutine CLEANUP_POPULATION_MOD deallocates all previously allocated arrays 
!              
!  NOTES:         
!
!******************************************************************************
!              
      !=================================================================
      ! CLEANUP_POPULATION_MOD begins here!
      !=================================================================
      IF ( ALLOCATED( POP_REDUCED ) ) DEALLOCATE( POP_REDUCED )
               
      ! Return to calling program
      END SUBROUTINE CLEANUP_POPULATION_MOD
                  
!------------------------------------------------------------------------------
      
      END MODULE POPULATION_MOD
