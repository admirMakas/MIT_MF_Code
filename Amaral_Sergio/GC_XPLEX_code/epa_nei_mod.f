!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: epa_nei_mod
!
! !DESCRIPTION: !  Module EPA_NEI_MOD contains variables and routines to read 
!  the weekday/weekend emissions from the EPA/NEI99 emissions inventory.
!\\
!\\
! !INTERFACE: 
!
      MODULE EPA_NEI_MOD
!
! !USES:
! 
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_EPA_NEI
      PUBLIC  :: EMISS_EPA_NEI
      PUBLIC  :: GET_USA_MASK
      PUBLIC  :: GET_EPA_ANTHRO
      PUBLIC  :: GET_EPA_BIOFUEL 
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: EMISS_EPA_NEI_AN
      PRIVATE :: EMISS_EPA_NEI_BF
      PRIVATE :: INIT_EPA_NEI
      PRIVATE :: READ_USA_MASK
      PRIVATE :: READ_USA_MASK_05x0666
      PRIVATE :: READ_EPA
      PRIVATE :: READ_EPA_05x0666
      PRIVATE :: TOTAL_ANTHRO_Tg
      PRIVATE :: TOTAL_BIOFUEL_Tg
!
! !REMARKS:
!  When you use EPA/NEI99 emissions (LNEI99=T), then the routines in this
!  module will read and process in the anthropogenic and biofuel emissions
!  from the EPA/NEI99 inventory.
!                                                                             .
!  When you use the EPA/NEI05 emissions (LNEI99=F and LNEI05=T), then the
!  routines in this module will read and process only the biofuel emissions
!  from the EPA/NEI99 inventory.  The EPA/NEI05 inventory only contains
!  anthropogenic emissions, so we are constrained to use the EPA/NEI99 
!  biofuel data.
! 
! !REVISION HISTORY:
!  (1 ) Prevent out of bounds errors in routines TOTAL_ANTHRO_TG and 
!        TOTAL_BIOFUEL_TG (bmy, 1/26/05)
!  (2 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now replace FMOL with TRACER_MW_KG (bmy, 10/25/05) 
!  (5 ) Now modified for IPCC future emissions (swu, bmy, 5/30/06)
!  (6 ) Now use int'annual scalars (amv, 08/24/07)
!  (7 ) Now can use ICARTT-based modifications from Hudman (phs, 1/26/09)
!        (a) Hudman et al., 2007, J. Geophys. Res., 112, D12S05, 
!             doi:10.1029/2006JD007912
!        (b) Hudman et al., 2008, Geophys. Res. Lett., 35, L04801, 
!             doi:10.1029/2007GL032393
!  (8 ) Now can read 0.5 x 0.667 nested grid emissions (amv, bmy, 12/18/09)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! USA Mask
      TYPE (XPLEX),  ALLOCATABLE :: USA_MASK(:,:)

      ! Fossil fuel arrays -- weekday
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_NOX(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_CO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_ALK4(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_ACET(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_MEK(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_PRPE(:,:)      
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_C2H6(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_C3H8(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_CH2O(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_NH3(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_SO2(:,:)      
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_AN_SO4(:,:)      
      
      ! Fossil fuel arrays -- weekend
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_NOX(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_CO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_ALK4(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_ACET(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_MEK(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_PRPE(:,:)      
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_C2H6(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_C3H8(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_CH2O(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_NH3(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_SO2(:,:) 
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_AN_SO4(:,:) 

      ! Biofuel arrays -- weekday
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_NOX(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_CO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_ALK4(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_ACET(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_MEK(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_PRPE(:,:)      
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_C2H6(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_C3H8(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_CH2O(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_NH3(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_SO2(:,:)  
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WD_BF_SO4(:,:)  

      ! Biofuel arrays -- weekend
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_NOX(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_CO(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_ALK4(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_ACET(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_MEK(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_PRPE(:,:)      
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_C2H6(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_C3H8(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_CH2O(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_NH3(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_SO2(:,:)  
      TYPE (XPLEX),  ALLOCATABLE :: EPA_WE_BF_SO4(:,:)  

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
! Prior to 2/7/11:
! Preserve original code here (bmy, 2/7/11)
!      
!      SUBROUTINE EMISS_EPA_NEI
!!
!!******************************************************************************
!!  Subroutine EMISS_EPA_NEI reads all EPA emissions from disk at the start
!!  of a new month. (rch, bmy, 11/10/04, 12/18/09)
!!
!!  NOTES:
!!  (1 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!!  (2 ) Modified for IPCC future emissions (swu, bmy, 5/30/06)
!!  (3 ) Now can read 0.5 x 0.667 nested grid emissions (amv, bmy, 12/18/09)
!!******************************************************************************
!!                               
!      ! References to F90 modules
!      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_ALK4ff  
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C2H6ff
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C3H8ff
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_CObf  
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff  
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3bf
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxbf
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff 
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_PRPEff
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2bf
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_TONEff  
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCbf
!      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCff
!      USE LOGICAL_MOD,          ONLY : LFUTURE, LICARTT
!      USE TIME_MOD,             ONLY : EXPAND_DATE,     GET_MONTH
!
!#     include "CMN_SIZE"             ! Size parameters
!
!      ! Local variables
!      LOGICAL, SAVE                 :: FIRST = .TRUE.
!      INTEGER                       :: I, J, THISMONTH, YYYYMMDD
!      TYPE (XPLEX)                        :: ALK4ff, C2H6ff, C3H8ff, COff
!      TYPE (XPLEX)                        :: NH3an,  NOxff,  PRPEff, SO2ff
!      TYPE (XPLEX)                        :: TONEff, VOCff,  CObf,   NH3bf
!      TYPE (XPLEX)                        :: NOxbf,  SO2bf,  VOCbf 
!      CHARACTER(LEN=255)            :: FILENAME
!
!      !=================================================================
!      ! EMISS_EPA_NEI begins here!
!      !=================================================================
!      
!      ! First-time initialization
!      IF ( FIRST ) THEN
!
!         ! Allocate arrays
!         CALL INIT_EPA_NEI
!
!         ! Read mask over the USA
!#if   defined(GRID05_0666)
!         CALL READ_USA_MASK_05x0666
!#else
!         CALL READ_USA_MASK
!#endif
!
!         ! Reset first-time flag
!         FIRST = .FALSE.
!      ENDIF
!
!      ! Get the current month
!      THISMONTH = GET_MONTH()
!
!      ! Get date for 1999 emissions
!      YYYYMMDD = 19990000 + ( THISMONTH * 100 ) + 01
!      
!      ! Echo info
!      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!      WRITE( 6, 100   ) 
! 100  FORMAT( 'E P A  /  N E I   E M I S S I O N S',
!     &       '  -- Baseline Year: 1999', / )
!    
!      !=================================================================
!      ! Read EPA weekday average anthropogenic emissions
!      !=================================================================
!
!      ! Weekday anthro file name
!      IF ( LICARTT ) THEN
!         
!        FILENAME = TRIM( DATA_DIR )                         // 
!     &             'EPA_NEI_200806/wkday_avg_an.YYYYMM.'    //
!     &             GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
!      ELSE
!         FILENAME = TRIM( DATA_DIR )                         // 
!     &              'EPA_NEI_200708/wkday_avg_an.YYYYMM.'    //
!     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
!      ENDIF
!         
!      ! Replace date in filename
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )
!
!      ! Read weekday data 
!#if defined(GRID05x0666)
!      CALL READ_EPA_05x0666( FILENAME,
!     &               EPA_WD_AN_NOX,  EPA_WD_AN_CO,   EPA_WD_AN_ALK4,
!     &               EPA_WD_AN_ACET, EPA_WD_AN_MEK,  EPA_WD_AN_PRPE,
!     &               EPA_WD_AN_C3H8, EPA_WD_AN_CH2O, EPA_WD_AN_C2H6,
!     &               EPA_WD_AN_SO2,  EPA_WD_AN_SO4,  EPA_WD_AN_NH3, 0)
!#else
!      CALL READ_EPA( FILENAME,       
!     &               EPA_WD_AN_NOX,  EPA_WD_AN_CO,   EPA_WD_AN_ALK4, 
!     &               EPA_WD_AN_ACET, EPA_WD_AN_MEK,  EPA_WD_AN_PRPE, 
!     &               EPA_WD_AN_C3H8, EPA_WD_AN_CH2O, EPA_WD_AN_C2H6, 
!     &               EPA_WD_AN_SO2,  EPA_WD_AN_SO4,  EPA_WD_AN_NH3 )
!#endif
!
!      !=================================================================
!      ! Read EPA weekend average anthropogenic emissions
!      !=================================================================
!
!      ! Weekday anthro file name
!      IF ( LICARTT ) THEN
!         
!        FILENAME = TRIM( DATA_DIR )                         // 
!     &             'EPA_NEI_200806/wkend_avg_an.YYYYMM.'    //
!     &             GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
!      ELSE
!         FILENAME = TRIM( DATA_DIR )                         // 
!     &              'EPA_NEI_200708/wkend_avg_an.YYYYMM.'    //
!     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
!      ENDIF
!
!      ! Replace date in filename
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )
!
!      ! Read weekend data 
!#if defined(GRID05x0666)
!      CALL READ_EPA_05x0666( FILENAME,
!     &               EPA_WE_AN_NOX,  EPA_WE_AN_CO,   EPA_WE_AN_ALK4,
!     &               EPA_WE_AN_ACET, EPA_WE_AN_MEK,  EPA_WE_AN_PRPE,
!     &               EPA_WE_AN_C3H8, EPA_WE_AN_CH2O, EPA_WE_AN_C2H6,
!     &               EPA_WE_AN_SO2,  EPA_WE_AN_SO4,  EPA_WE_AN_NH3, 0)
!#else
!      CALL READ_EPA( FILENAME,       
!     &               EPA_WE_AN_NOX,  EPA_WE_AN_CO,   EPA_WE_AN_ALK4, 
!     &               EPA_WE_AN_ACET, EPA_WE_AN_MEK,  EPA_WE_AN_PRPE, 
!     &               EPA_WE_AN_C3H8, EPA_WE_AN_CH2O, EPA_WE_AN_C2H6, 
!     &               EPA_WE_AN_SO2,  EPA_WE_AN_SO4,  EPA_WE_AN_NH3 )
!#endif
!
!      !=================================================================
!      ! Read EPA weekday average biofuel emissions
!      !=================================================================
!
!      ! Weekday biofuel file name
!      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'EPA_NEI_200411/wkday_avg_bf.YYYYMM.'    //
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
!
!      ! Replace date in filename
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )
!
!      ! Read weekday data 
!#if defined(GRID05x0666)
!      CALL READ_EPA_05x0666( FILENAME,
!     &               EPA_WD_BF_NOX,  EPA_WD_BF_CO,   EPA_WD_BF_ALK4,
!     &               EPA_WD_BF_ACET, EPA_WD_BF_MEK,  EPA_WD_BF_PRPE,
!     &               EPA_WD_BF_C3H8, EPA_WD_BF_CH2O, EPA_WD_BF_C2H6,
!     &               EPA_WD_BF_SO2,  EPA_WD_BF_SO4,  EPA_WD_BF_NH3, 1)
!#else
!      CALL READ_EPA( FILENAME,       
!     &               EPA_WD_BF_NOX,  EPA_WD_BF_CO,   EPA_WD_BF_ALK4, 
!     &               EPA_WD_BF_ACET, EPA_WD_BF_MEK,  EPA_WD_BF_PRPE, 
!     &               EPA_WD_BF_C3H8, EPA_WD_BF_CH2O, EPA_WD_BF_C2H6, 
!     &               EPA_WD_BF_SO2,  EPA_WD_BF_SO4,  EPA_WD_BF_NH3 )
!#endif
!
!      !=================================================================
!      ! Read EPA weekend average biofuel emissions
!      !=================================================================
!
!      ! Weekend biofuel file name
!      FILENAME = TRIM( DATA_DIR )                         // 
!     &           'EPA_NEI_200411/wkend_avg_bf.YYYYMM.'    //
!     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()  
!
!      ! Replace date in filename
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )
!
!      ! Read weekend data 
!#if defined(GRID05x0666)
!      CALL READ_EPA_05x0666( FILENAME,
!     &               EPA_WE_BF_NOX,  EPA_WE_BF_CO,   EPA_WE_BF_ALK4,
!     &               EPA_WE_BF_ACET, EPA_WE_BF_MEK,  EPA_WE_BF_PRPE,
!     &               EPA_WE_BF_C3H8, EPA_WE_BF_CH2O, EPA_WE_BF_C2H6,
!     &               EPA_WE_BF_SO2,  EPA_WE_BF_SO4,  EPA_WE_BF_NH3, 1)
!#else
!      CALL READ_EPA( FILENAME,       
!     &               EPA_WE_BF_NOX,  EPA_WE_BF_CO,   EPA_WE_BF_ALK4, 
!     &               EPA_WE_BF_ACET, EPA_WE_BF_MEK,  EPA_WE_BF_PRPE, 
!     &               EPA_WE_BF_C3H8, EPA_WE_BF_CH2O, EPA_WE_BF_C2H6, 
!     &               EPA_WE_BF_SO2,  EPA_WE_BF_SO4,  EPA_WE_BF_NH3 )
!#endif
!
!      !=================================================================
!      ! Apply USA Mask (keep emissions over US, zero elsewhere)
!      !=================================================================
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,     J,     ALK4ff, C2H6ff, C3H8ff, COff  )
!!$OMP+PRIVATE( NH3an, NOxff, PRPEff, SO2ff,  TONEff, VOCff )         
!!$OMP+PRIVATE( CObf,  NH3bf, NOxbf,  SO2bf,  VOCbf         )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Weekday avg anthro
!         EPA_WD_AN_NOX (I,J) = EPA_WD_AN_NOX (I,J) * USA_MASK(I,J)
!         EPA_WD_AN_CO  (I,J) = EPA_WD_AN_CO  (I,J) * USA_MASK(I,J) 
!         EPA_WD_AN_ALK4(I,J) = EPA_WD_AN_ALK4(I,J) * USA_MASK(I,J) 
!         EPA_WD_AN_ACET(I,J) = EPA_WD_AN_ACET(I,J) * USA_MASK(I,J) 
!         EPA_WD_AN_MEK (I,J) = EPA_WD_AN_MEK (I,J) * USA_MASK(I,J) 
!         EPA_WD_AN_PRPE(I,J) = EPA_WD_AN_PRPE(I,J) * USA_MASK(I,J)      
!         EPA_WD_AN_C3H8(I,J) = EPA_WD_AN_C3H8(I,J) * USA_MASK(I,J) 
!         EPA_WD_AN_CH2O(I,J) = EPA_WD_AN_CH2O(I,J) * USA_MASK(I,J) 
!         EPA_WD_AN_C2H6(I,J) = EPA_WD_AN_C2H6(I,J) * USA_MASK(I,J)
!         EPA_WD_AN_SO2 (I,J) = EPA_WD_AN_SO2 (I,J) * USA_MASK(I,J)   
!         EPA_WD_AN_SO4 (I,J) = EPA_WD_AN_SO4 (I,J) * USA_MASK(I,J)  
!         EPA_WD_AN_NH3 (I,J) = EPA_WD_AN_NH3 (I,J) * USA_MASK(I,J) 
!                             
!         ! Weekend avg anthro     
!         EPA_WE_AN_NOX (I,J) = EPA_WE_AN_NOX (I,J) * USA_MASK(I,J)
!         EPA_WE_AN_CO  (I,J) = EPA_WE_AN_CO  (I,J) * USA_MASK(I,J)
!         EPA_WE_AN_ALK4(I,J) = EPA_WE_AN_ALK4(I,J) * USA_MASK(I,J)
!         EPA_WE_AN_ACET(I,J) = EPA_WE_AN_ACET(I,J) * USA_MASK(I,J)
!         EPA_WE_AN_MEK (I,J) = EPA_WE_AN_MEK (I,J) * USA_MASK(I,J)
!         EPA_WE_AN_PRPE(I,J) = EPA_WE_AN_PRPE(I,J) * USA_MASK(I,J)   
!         EPA_WE_AN_C3H8(I,J) = EPA_WE_AN_C3H8(I,J) * USA_MASK(I,J)
!         EPA_WE_AN_CH2O(I,J) = EPA_WE_AN_CH2O(I,J) * USA_MASK(I,J)
!         EPA_WE_AN_C2H6(I,J) = EPA_WE_AN_C2H6(I,J) * USA_MASK(I,J)
!         EPA_WE_AN_SO2 (I,J) = EPA_WE_AN_SO2 (I,J) * USA_MASK(I,J)
!         EPA_WE_AN_SO4 (I,J) = EPA_WE_AN_SO4 (I,J) * USA_MASK(I,J)
!         EPA_WE_AN_NH3 (I,J) = EPA_WE_AN_NH3 (I,J) * USA_MASK(I,J)
!                              
!         ! Weekday avg biofuel    
!         EPA_WD_BF_NOX (I,J) = EPA_WD_BF_NOX (I,J) * USA_MASK(I,J)
!         EPA_WD_BF_CO  (I,J) = EPA_WD_BF_CO  (I,J) * USA_MASK(I,J)
!         EPA_WD_BF_ALK4(I,J) = EPA_WD_BF_ALK4(I,J) * USA_MASK(I,J)
!         EPA_WD_BF_ACET(I,J) = EPA_WD_BF_ACET(I,J) * USA_MASK(I,J)
!         EPA_WD_BF_MEK (I,J) = EPA_WD_BF_MEK (I,J) * USA_MASK(I,J)
!         EPA_WD_BF_PRPE(I,J) = EPA_WD_BF_PRPE(I,J) * USA_MASK(I,J)
!         EPA_WD_BF_C3H8(I,J) = EPA_WD_BF_C3H8(I,J) * USA_MASK(I,J)
!         EPA_WD_BF_CH2O(I,J) = EPA_WD_BF_CH2O(I,J) * USA_MASK(I,J)    
!         EPA_WD_BF_C2H6(I,J) = EPA_WD_BF_C2H6(I,J) * USA_MASK(I,J)
!         EPA_WD_BF_SO2 (I,J) = EPA_WD_BF_SO2 (I,J) * USA_MASK(I,J)
!         EPA_WD_BF_SO4 (I,J) = EPA_WD_BF_SO4 (I,J) * USA_MASK(I,J)
!         EPA_WD_BF_NH3 (I,J) = EPA_WD_BF_NH3 (I,J) * USA_MASK(I,J)
!
!         ! Weekend avg biofuel    
!         EPA_WE_BF_NOX (I,J) = EPA_WE_BF_NOX (I,J) * USA_MASK(I,J)
!         EPA_WE_BF_CO  (I,J) = EPA_WE_BF_CO  (I,J) * USA_MASK(I,J)
!         EPA_WE_BF_ALK4(I,J) = EPA_WE_BF_ALK4(I,J) * USA_MASK(I,J)
!         EPA_WE_BF_ACET(I,J) = EPA_WE_BF_ACET(I,J) * USA_MASK(I,J)
!         EPA_WE_BF_MEK (I,J) = EPA_WE_BF_MEK (I,J) * USA_MASK(I,J)
!         EPA_WE_BF_PRPE(I,J) = EPA_WE_BF_PRPE(I,J) * USA_MASK(I,J)   
!         EPA_WE_BF_C3H8(I,J) = EPA_WE_BF_C3H8(I,J) * USA_MASK(I,J)
!         EPA_WE_BF_CH2O(I,J) = EPA_WE_BF_CH2O(I,J) * USA_MASK(I,J)
!         EPA_WE_BF_C2H6(I,J) = EPA_WE_BF_C2H6(I,J) * USA_MASK(I,J)
!         EPA_WE_BF_SO2 (I,J) = EPA_WE_BF_SO2 (I,J) * USA_MASK(I,J)
!         EPA_WE_BF_SO4 (I,J) = EPA_WE_BF_SO4 (I,J) * USA_MASK(I,J)
!         EPA_WE_BF_NH3 (I,J) = EPA_WE_BF_NH3 (I,J) * USA_MASK(I,J)
!
!         !----------------------------------------------
!         ! Compute IPCC future emissions (if necessary)
!         !----------------------------------------------
!         IF ( LFUTURE .and. USA_MASK(I,J) > 0d0 ) THEN
!
!            ! Future anthro scale factors
!            ALK4ff = GET_FUTURE_SCALE_ALK4ff( I, J )
!            C2H6ff = GET_FUTURE_SCALE_C2H6ff( I, J )
!            C3H8ff = GET_FUTURE_SCALE_C3H8ff( I, J )
!            COff   = GET_FUTURE_SCALE_COff(   I, J )
!            NH3an  = GET_FUTURE_SCALE_NH3an(  I, J )
!            NOxff  = GET_FUTURE_SCALE_NOxff(  I, J )
!            PRPEff = GET_FUTURE_SCALE_PRPEff( I, J )
!            TONEff = GET_FUTURE_SCALE_TONEff( I, J )
!            SO2ff  = GET_FUTURE_SCALE_SO2ff(  I, J )
!            VOCff  = GET_FUTURE_SCALE_VOCff(  I, J )
!
!            ! Future biofuel scale factors            
!            CObf   = GET_FUTURE_SCALE_CObf(   I, J )
!            NH3bf  = GET_FUTURE_SCALE_NH3bf(  I, J )
!            NOxbf  = GET_FUTURE_SCALE_NOXbf(  I, J )
!            SO2bf  = GET_FUTURE_SCALE_SO2bf(  I, J )
!            VOCbf  = GET_FUTURE_SCALE_VOCbf(  I, J )
!
!            ! Future weekday avg anthro
!            EPA_WD_AN_NOX (I,J) = EPA_WD_AN_NOX (I,J) * NOxff
!            EPA_WD_AN_CO  (I,J) = EPA_WD_AN_CO  (I,J) * COff
!            EPA_WD_AN_ALK4(I,J) = EPA_WD_AN_ALK4(I,J) * ALK4ff
!            EPA_WD_AN_ACET(I,J) = EPA_WD_AN_ACET(I,J) * TONEff
!            EPA_WD_AN_MEK (I,J) = EPA_WD_AN_MEK (I,J) * TONEff
!            EPA_WD_AN_PRPE(I,J) = EPA_WD_AN_PRPE(I,J) * PRPEff
!            EPA_WD_AN_C3H8(I,J) = EPA_WD_AN_C3H8(I,J) * C3H8ff 
!            EPA_WD_AN_CH2O(I,J) = EPA_WD_AN_CH2O(I,J) * VOCff
!            EPA_WD_AN_C2H6(I,J) = EPA_WD_AN_C2H6(I,J) * C2H6ff
!            EPA_WD_AN_SO2 (I,J) = EPA_WD_AN_SO2 (I,J) * SO2ff
!            EPA_WD_AN_SO4 (I,J) = EPA_WD_AN_SO4 (I,J) * SO2ff
!            EPA_WD_AN_NH3 (I,J) = EPA_WD_AN_NH3 (I,J) * NH3an 
!                             
!            ! Weekend avg anthro     
!            EPA_WE_AN_NOX (I,J) = EPA_WE_AN_NOX (I,J) * NOxff
!            EPA_WE_AN_CO  (I,J) = EPA_WE_AN_CO  (I,J) * COff
!            EPA_WE_AN_ALK4(I,J) = EPA_WE_AN_ALK4(I,J) * ALK4ff
!            EPA_WE_AN_ACET(I,J) = EPA_WE_AN_ACET(I,J) * TONEff
!            EPA_WE_AN_MEK (I,J) = EPA_WE_AN_MEK (I,J) * TONEff
!            EPA_WE_AN_PRPE(I,J) = EPA_WE_AN_PRPE(I,J) * PRPEff
!            EPA_WE_AN_C3H8(I,J) = EPA_WE_AN_C3H8(I,J) * C3H8ff
!            EPA_WE_AN_CH2O(I,J) = EPA_WE_AN_CH2O(I,J) * VOCff
!            EPA_WE_AN_C2H6(I,J) = EPA_WE_AN_C2H6(I,J) * C2H6ff 
!            EPA_WE_AN_SO2 (I,J) = EPA_WE_AN_SO2 (I,J) * SO2ff
!            EPA_WE_AN_SO4 (I,J) = EPA_WE_AN_SO4 (I,J) * SO2ff
!            EPA_WE_AN_NH3 (I,J) = EPA_WE_AN_NH3 (I,J) * NH3an 
!                              
!            ! Weekday avg biofuel    
!            EPA_WD_BF_NOX (I,J) = EPA_WD_BF_NOX (I,J) * NOxbf
!            EPA_WD_BF_CO  (I,J) = EPA_WD_BF_CO  (I,J) * CObf
!            EPA_WD_BF_ALK4(I,J) = EPA_WD_BF_ALK4(I,J) * VOCbf
!            EPA_WD_BF_ACET(I,J) = EPA_WD_BF_ACET(I,J) * VOCbf
!            EPA_WD_BF_MEK (I,J) = EPA_WD_BF_MEK (I,J) * VOCbf
!            EPA_WD_BF_PRPE(I,J) = EPA_WD_BF_PRPE(I,J) * VOCbf 
!            EPA_WD_BF_C3H8(I,J) = EPA_WD_BF_C3H8(I,J) * VOCbf 
!            EPA_WD_BF_CH2O(I,J) = EPA_WD_BF_CH2O(I,J) * VOCbf     
!            EPA_WD_BF_C2H6(I,J) = EPA_WD_BF_C2H6(I,J) * VOCbf 
!            EPA_WD_BF_SO2 (I,J) = EPA_WD_BF_SO2 (I,J) * SO2bf 
!            EPA_WD_BF_SO4 (I,J) = EPA_WD_BF_SO4 (I,J) * SO2bf
!            EPA_WD_BF_NH3 (I,J) = EPA_WD_BF_NH3 (I,J) * NH3bf
!            
!            ! Weekend avg biofuel    
!            EPA_WE_BF_NOX (I,J) = EPA_WE_BF_NOX (I,J) * NOxbf
!            EPA_WE_BF_CO  (I,J) = EPA_WE_BF_CO  (I,J) * CObf
!            EPA_WE_BF_ALK4(I,J) = EPA_WE_BF_ALK4(I,J) * VOCbf 
!            EPA_WE_BF_ACET(I,J) = EPA_WE_BF_ACET(I,J) * VOCbf 
!            EPA_WE_BF_MEK (I,J) = EPA_WE_BF_MEK (I,J) * VOCbf 
!            EPA_WE_BF_PRPE(I,J) = EPA_WE_BF_PRPE(I,J) * VOCbf   
!            EPA_WE_BF_C3H8(I,J) = EPA_WE_BF_C3H8(I,J) * VOCbf 
!            EPA_WE_BF_CH2O(I,J) = EPA_WE_BF_CH2O(I,J) * VOCbf 
!            EPA_WE_BF_C2H6(I,J) = EPA_WE_BF_C2H6(I,J) * VOCbf 
!            EPA_WE_BF_SO2 (I,J) = EPA_WE_BF_SO2 (I,J) * SO2bf
!            EPA_WE_BF_SO4 (I,J) = EPA_WE_BF_SO4 (I,J) * SO2bf
!            EPA_WE_BF_NH3 (I,J) = EPA_WE_BF_NH3 (I,J) * NH3bf
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Print totals to the log file
!      CALL TOTAL_ANTHRO_TG( THISMONTH )
!      CALL TOTAL_BIOFUEL_TG( THISMONTH )
!
!      ! Fancy output
!      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!
!      ! Return to calling program
!      END SUBROUTINE EMISS_EPA_NEI
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_epa_nei
!
! !DESCRIPTION: Subroutine EMISS\_EPA\_NEI reads all EPA emissions from disk 
!  at the start of a new month.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_EPA_NEI
!
! !USES:
!
      USE LOGICAL_MOD, ONLY : LNEI05 
      USE LOGICAL_MOD, ONLY : LNEI99
!
! !REMARKS:
!  Read EPA/NEI99 anthro + EPA/NEI99 biofuel if LNEI99=T, or
!  Read EPA/NEI05 anthro + EPA/NEI99 biofuel if LNEI99=F but LNEI05=T
! 
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  07 Feb 2011 - R. Yantosca - Now reads anthro & biofuel emissions separately
!   7 Feb 2011 - R. Yantosca - Now reads biofuel only when LNEI05=T
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      
      !=================================================================
      ! First-time initialization
      !=================================================================
      IF ( FIRST ) THEN

         ! Allocate arrays
         CALL INIT_EPA_NEI

#if   defined( GRID05_0666 )
         ! Read mask over the USA, nested-grid resolution
         CALL READ_USA_MASK_05x0666
#else
         ! Read mask over the USA, global resolution
         CALL READ_USA_MASK
#endif

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read data
      !=================================================================

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) 
 100  FORMAT( 'E P A  /  N E I 9 9   E M I S S I O N S', / )

      IF ( LNEI99 ) THEN

         ! If we are using EPA/NEI99, then we will read both the anthro
         ! and biofuel emissions data from disk for each month.
         WRITE( 6, 120 )
         CALL EMISS_EPA_NEI_AN
         CALL EMISS_EPA_NEI_BF

      ELSE 

         ! If we are not using EPA/NEI99, then check if we are using 
         ! EPA/NEI05.  The EPA/NEI05 inventory only has anthro emissions.  
         ! In this case, we must read the EPA/NEI99 biofuel emissions, since 
         ! these data are better defined over the USA than the "default"
         ! biomass data from Yevich & Logan 2003.
         IF ( LNEI05 ) THEN
            WRITE( 6, 130 )
            CALL EMISS_EPA_NEI_BF
         ENDIF

      ENDIF

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! FORMATs
 120  FORMAT( '%%% Using EPA/NEI99 anthro + EPA/NEI99 biofuel %%%', / )
 130  FORMAT( '%%% Using EPA/NEI05 anthro + EPA/NEI99 biofuel %%%', / )

      END SUBROUTINE EMISS_EPA_NEI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_epa_nei_an
!
! !DESCRIPTION: Subroutine EMISS\_EPA\_NEI\_AN reads only the EPA/NEI99
!  anthropogenic emissions from disk at the start of a new month.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_EPA_NEI_AN
!
! !USES:
!
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,            ONLY : GET_RES_EXT
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_ALK4ff  
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C2H6ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C3H8ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff  
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_PRPEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_TONEff  
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCff
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE LOGICAL_MOD,          ONLY : LICARTT
      USE TIME_MOD,             ONLY : EXPAND_DATE
      USE TIME_MOD,             ONLY : GET_MONTH

#     include "CMN_SIZE"             ! Size parameters
!
! !REMARKS:
!  Split off from EMISS_EPA_NEI.
! 
! !REVISION HISTORY: 
!  (1 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (2 ) Modified for IPCC future emissions (swu, bmy, 5/30/06)
!  (3 ) Now can read 0.5 x 0.667 nested grid emissions (amv, bmy, 12/18/09)
!  07 Feb 2011 - R. Yantosca - Now only read anthro emissions.
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I,      J,      THISMONTH, YYYYMMDD
      TYPE (XPLEX)             :: ALK4ff, C2H6ff, C3H8ff,    COff
      TYPE (XPLEX)             :: NH3an,  NOxff,  PRPEff,    SO2ff
      TYPE (XPLEX)             :: TONEff, VOCff 
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! EMISS_EPA_NEI_AN begins here!
      !=================================================================
      
      ! Get the current month
      THISMONTH = GET_MONTH()

      ! Get date for 1999 emissions
      YYYYMMDD = 19990000 + ( THISMONTH * 100 ) + 01
      
      !=================================================================
      ! Read EPA weekday average anthropogenic emissions
      !=================================================================

      ! Weekday anthro file name
      IF ( LICARTT ) THEN
         FILENAME = TRIM( DATA_DIR )                         // 
     &              'EPA_NEI_200806/wkday_avg_an.YYYYMM.'    //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
      ELSE
         FILENAME = TRIM( DATA_DIR )                         // 
     &              'EPA_NEI_200708/wkday_avg_an.YYYYMM.'    //
     &               GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
      ENDIF
         
      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

#if   defined( GRID05x0666 )

      ! Read weekday data, nested grids
      CALL READ_EPA_05x0666( FILENAME,
     &               EPA_WD_AN_NOX,  EPA_WD_AN_CO,   EPA_WD_AN_ALK4,
     &               EPA_WD_AN_ACET, EPA_WD_AN_MEK,  EPA_WD_AN_PRPE,
     &               EPA_WD_AN_C3H8, EPA_WD_AN_CH2O, EPA_WD_AN_C2H6,
     &               EPA_WD_AN_SO2,  EPA_WD_AN_SO4,  EPA_WD_AN_NH3, 0)

#else

      ! Read weekday data, global grids
      CALL READ_EPA( FILENAME,       
     &               EPA_WD_AN_NOX,  EPA_WD_AN_CO,   EPA_WD_AN_ALK4, 
     &               EPA_WD_AN_ACET, EPA_WD_AN_MEK,  EPA_WD_AN_PRPE, 
     &               EPA_WD_AN_C3H8, EPA_WD_AN_CH2O, EPA_WD_AN_C2H6, 
     &               EPA_WD_AN_SO2,  EPA_WD_AN_SO4,  EPA_WD_AN_NH3 )

#endif

      !=================================================================
      ! Read EPA weekend average anthropogenic emissions
      !=================================================================

      ! Weekday anthro file name
      IF ( LICARTT ) THEN
         FILENAME = TRIM( DATA_DIR )                         // 
     &              'EPA_NEI_200806/wkend_avg_an.YYYYMM.'    //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
      ELSE
         FILENAME = TRIM( DATA_DIR )                         // 
     &              'EPA_NEI_200708/wkend_avg_an.YYYYMM.'    //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     
      ENDIF

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

#if   defined( GRID05x0666 )

      ! Read weekend data, nested grid
      CALL READ_EPA_05x0666( FILENAME,
     &               EPA_WE_AN_NOX,  EPA_WE_AN_CO,   EPA_WE_AN_ALK4,
     &               EPA_WE_AN_ACET, EPA_WE_AN_MEK,  EPA_WE_AN_PRPE,
     &               EPA_WE_AN_C3H8, EPA_WE_AN_CH2O, EPA_WE_AN_C2H6,
     &               EPA_WE_AN_SO2,  EPA_WE_AN_SO4,  EPA_WE_AN_NH3, 0)

#else

      ! Read weekend data, global grids
      CALL READ_EPA( FILENAME,       
     &               EPA_WE_AN_NOX,  EPA_WE_AN_CO,   EPA_WE_AN_ALK4, 
     &               EPA_WE_AN_ACET, EPA_WE_AN_MEK,  EPA_WE_AN_PRPE, 
     &               EPA_WE_AN_C3H8, EPA_WE_AN_CH2O, EPA_WE_AN_C2H6, 
     &               EPA_WE_AN_SO2,  EPA_WE_AN_SO4,  EPA_WE_AN_NH3 )

#endif

      !=================================================================
      ! Apply USA Mask (keep emissions over US, zero elsewhere)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,     ALK4ff, C2H6ff, C3H8ff, COff  )
!$OMP+PRIVATE( NH3an, NOxff, PRPEff, SO2ff,  TONEff, VOCff )         
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Weekday avg anthro
         EPA_WD_AN_NOX (I,J) = EPA_WD_AN_NOX (I,J) * USA_MASK(I,J)
         EPA_WD_AN_CO  (I,J) = EPA_WD_AN_CO  (I,J) * USA_MASK(I,J) 
         EPA_WD_AN_ALK4(I,J) = EPA_WD_AN_ALK4(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_ACET(I,J) = EPA_WD_AN_ACET(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_MEK (I,J) = EPA_WD_AN_MEK (I,J) * USA_MASK(I,J) 
         EPA_WD_AN_PRPE(I,J) = EPA_WD_AN_PRPE(I,J) * USA_MASK(I,J)      
         EPA_WD_AN_C3H8(I,J) = EPA_WD_AN_C3H8(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_CH2O(I,J) = EPA_WD_AN_CH2O(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_C2H6(I,J) = EPA_WD_AN_C2H6(I,J) * USA_MASK(I,J)
         EPA_WD_AN_SO2 (I,J) = EPA_WD_AN_SO2 (I,J) * USA_MASK(I,J)   
         EPA_WD_AN_SO4 (I,J) = EPA_WD_AN_SO4 (I,J) * USA_MASK(I,J)  
         EPA_WD_AN_NH3 (I,J) = EPA_WD_AN_NH3 (I,J) * USA_MASK(I,J) 
                             
         ! Weekend avg anthro     
         EPA_WE_AN_NOX (I,J) = EPA_WE_AN_NOX (I,J) * USA_MASK(I,J)
         EPA_WE_AN_CO  (I,J) = EPA_WE_AN_CO  (I,J) * USA_MASK(I,J)
         EPA_WE_AN_ALK4(I,J) = EPA_WE_AN_ALK4(I,J) * USA_MASK(I,J)
         EPA_WE_AN_ACET(I,J) = EPA_WE_AN_ACET(I,J) * USA_MASK(I,J)
         EPA_WE_AN_MEK (I,J) = EPA_WE_AN_MEK (I,J) * USA_MASK(I,J)
         EPA_WE_AN_PRPE(I,J) = EPA_WE_AN_PRPE(I,J) * USA_MASK(I,J)   
         EPA_WE_AN_C3H8(I,J) = EPA_WE_AN_C3H8(I,J) * USA_MASK(I,J)
         EPA_WE_AN_CH2O(I,J) = EPA_WE_AN_CH2O(I,J) * USA_MASK(I,J)
         EPA_WE_AN_C2H6(I,J) = EPA_WE_AN_C2H6(I,J) * USA_MASK(I,J)
         EPA_WE_AN_SO2 (I,J) = EPA_WE_AN_SO2 (I,J) * USA_MASK(I,J)
         EPA_WE_AN_SO4 (I,J) = EPA_WE_AN_SO4 (I,J) * USA_MASK(I,J)
         EPA_WE_AN_NH3 (I,J) = EPA_WE_AN_NH3 (I,J) * USA_MASK(I,J)
                              
         !----------------------------------------------
         ! Compute IPCC future emissions (if necessary)
         !----------------------------------------------
         IF ( LFUTURE .and. USA_MASK(I,J) > 0d0 ) THEN

            ! Future anthro scale factors
            ALK4ff = GET_FUTURE_SCALE_ALK4ff( I, J )
            C2H6ff = GET_FUTURE_SCALE_C2H6ff( I, J )
            C3H8ff = GET_FUTURE_SCALE_C3H8ff( I, J )
            COff   = GET_FUTURE_SCALE_COff(   I, J )
            NH3an  = GET_FUTURE_SCALE_NH3an(  I, J )
            NOxff  = GET_FUTURE_SCALE_NOxff(  I, J )
            PRPEff = GET_FUTURE_SCALE_PRPEff( I, J )
            TONEff = GET_FUTURE_SCALE_TONEff( I, J )
            SO2ff  = GET_FUTURE_SCALE_SO2ff(  I, J )
            VOCff  = GET_FUTURE_SCALE_VOCff(  I, J )

            ! Future weekday avg anthro
            EPA_WD_AN_NOX (I,J) = EPA_WD_AN_NOX (I,J) * NOxff
            EPA_WD_AN_CO  (I,J) = EPA_WD_AN_CO  (I,J) * COff
            EPA_WD_AN_ALK4(I,J) = EPA_WD_AN_ALK4(I,J) * ALK4ff
            EPA_WD_AN_ACET(I,J) = EPA_WD_AN_ACET(I,J) * TONEff
            EPA_WD_AN_MEK (I,J) = EPA_WD_AN_MEK (I,J) * TONEff
            EPA_WD_AN_PRPE(I,J) = EPA_WD_AN_PRPE(I,J) * PRPEff
            EPA_WD_AN_C3H8(I,J) = EPA_WD_AN_C3H8(I,J) * C3H8ff 
            EPA_WD_AN_CH2O(I,J) = EPA_WD_AN_CH2O(I,J) * VOCff
            EPA_WD_AN_C2H6(I,J) = EPA_WD_AN_C2H6(I,J) * C2H6ff
            EPA_WD_AN_SO2 (I,J) = EPA_WD_AN_SO2 (I,J) * SO2ff
            EPA_WD_AN_SO4 (I,J) = EPA_WD_AN_SO4 (I,J) * SO2ff
            EPA_WD_AN_NH3 (I,J) = EPA_WD_AN_NH3 (I,J) * NH3an 
                             
            ! Weekend avg anthro     
            EPA_WE_AN_NOX (I,J) = EPA_WE_AN_NOX (I,J) * NOxff
            EPA_WE_AN_CO  (I,J) = EPA_WE_AN_CO  (I,J) * COff
            EPA_WE_AN_ALK4(I,J) = EPA_WE_AN_ALK4(I,J) * ALK4ff
            EPA_WE_AN_ACET(I,J) = EPA_WE_AN_ACET(I,J) * TONEff
            EPA_WE_AN_MEK (I,J) = EPA_WE_AN_MEK (I,J) * TONEff
            EPA_WE_AN_PRPE(I,J) = EPA_WE_AN_PRPE(I,J) * PRPEff
            EPA_WE_AN_C3H8(I,J) = EPA_WE_AN_C3H8(I,J) * C3H8ff
            EPA_WE_AN_CH2O(I,J) = EPA_WE_AN_CH2O(I,J) * VOCff
            EPA_WE_AN_C2H6(I,J) = EPA_WE_AN_C2H6(I,J) * C2H6ff 
            EPA_WE_AN_SO2 (I,J) = EPA_WE_AN_SO2 (I,J) * SO2ff
            EPA_WE_AN_SO4 (I,J) = EPA_WE_AN_SO4 (I,J) * SO2ff
            EPA_WE_AN_NH3 (I,J) = EPA_WE_AN_NH3 (I,J) * NH3an 
                              
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Print totals to the log file
      CALL TOTAL_ANTHRO_TG( THISMONTH )

      END SUBROUTINE EMISS_EPA_NEI_AN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_epa_nei_bf
!
! !DESCRIPTION: Subroutine EMISS\_EPA\_NEI\_BF reads only the EPA/NEI99
!  biofuel emissions from disk at the start of a new month.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_EPA_NEI_BF
!
! !USES:
!
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,            ONLY : GET_RES_EXT
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_CObf  
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3bf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxbf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2bf
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCbf
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE LOGICAL_MOD,          ONLY : LICARTT
      USE TIME_MOD,             ONLY : EXPAND_DATE
      USE TIME_MOD,             ONLY : GET_MONTH

#     include "CMN_SIZE"             ! Size parameters

!
! !REMARKS:
!  Split off from EMISS_EPA_NEI.
! 
! !REVISION HISTORY: 
!  (1 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (2 ) Modified for IPCC future emissions (swu, bmy, 5/30/06)
!  (3 ) Now can read 0.5 x 0.667 nested grid emissions (amv, bmy, 12/18/09)
!  07 Feb 2011 - R. Yantosca - Now only read biofuel emissions
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I,     J,     THISMONTH, YYYYMMDD
      TYPE (XPLEX)          :: CObf,  NH3bf, NOxbf,     SO2bf,   VOCbf 
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! EMISS_EPA_NEI_BF begins here!
      !=================================================================
      
      ! Get the current month
      THISMONTH = GET_MONTH()

      ! Get date for 1999 emissions
      YYYYMMDD = 19990000 + ( THISMONTH * 100 ) + 01
      
      !=================================================================
      ! Read EPA weekday average biofuel emissions
      !=================================================================

      ! Weekday biofuel file name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'EPA_NEI_200411/wkday_avg_bf.YYYYMM.'    //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

#if   defined( GRID05x0666 )

      ! Read weekday data, nested grids
      CALL READ_EPA_05x0666( FILENAME,
     &               EPA_WD_BF_NOX,  EPA_WD_BF_CO,   EPA_WD_BF_ALK4,
     &               EPA_WD_BF_ACET, EPA_WD_BF_MEK,  EPA_WD_BF_PRPE,
     &               EPA_WD_BF_C3H8, EPA_WD_BF_CH2O, EPA_WD_BF_C2H6,
     &               EPA_WD_BF_SO2,  EPA_WD_BF_SO4,  EPA_WD_BF_NH3, 1)

#else

      ! Read weekday data, global grids
      CALL READ_EPA( FILENAME,       
     &               EPA_WD_BF_NOX,  EPA_WD_BF_CO,   EPA_WD_BF_ALK4, 
     &               EPA_WD_BF_ACET, EPA_WD_BF_MEK,  EPA_WD_BF_PRPE, 
     &               EPA_WD_BF_C3H8, EPA_WD_BF_CH2O, EPA_WD_BF_C2H6, 
     &               EPA_WD_BF_SO2,  EPA_WD_BF_SO4,  EPA_WD_BF_NH3 )

#endif

      !=================================================================
      ! Read EPA weekend average biofuel emissions
      !=================================================================

      ! Weekend biofuel file name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'EPA_NEI_200411/wkend_avg_bf.YYYYMM.'    //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()  

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

#if   defined( GRID05x0666 )

      ! Read weekend data, nested-grids
      CALL READ_EPA_05x0666( FILENAME,
     &               EPA_WE_BF_NOX,  EPA_WE_BF_CO,   EPA_WE_BF_ALK4,
     &               EPA_WE_BF_ACET, EPA_WE_BF_MEK,  EPA_WE_BF_PRPE,
     &               EPA_WE_BF_C3H8, EPA_WE_BF_CH2O, EPA_WE_BF_C2H6,
     &               EPA_WE_BF_SO2,  EPA_WE_BF_SO4,  EPA_WE_BF_NH3, 1 )

#else

      ! Read weekend data, global grids
      CALL READ_EPA( FILENAME,       
     &               EPA_WE_BF_NOX,  EPA_WE_BF_CO,   EPA_WE_BF_ALK4, 
     &               EPA_WE_BF_ACET, EPA_WE_BF_MEK,  EPA_WE_BF_PRPE, 
     &               EPA_WE_BF_C3H8, EPA_WE_BF_CH2O, EPA_WE_BF_C2H6, 
     &               EPA_WE_BF_SO2,  EPA_WE_BF_SO4,  EPA_WE_BF_NH3 )

#endif

      !=================================================================
      ! Apply USA Mask (keep emissions over US, zero elsewhere)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CObf, NH3bf, NOxbf, SO2bf, VOCbf )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                              
         ! Weekday avg biofuel    
         EPA_WD_BF_NOX (I,J) = EPA_WD_BF_NOX (I,J) * USA_MASK(I,J)
         EPA_WD_BF_CO  (I,J) = EPA_WD_BF_CO  (I,J) * USA_MASK(I,J)
         EPA_WD_BF_ALK4(I,J) = EPA_WD_BF_ALK4(I,J) * USA_MASK(I,J)
         EPA_WD_BF_ACET(I,J) = EPA_WD_BF_ACET(I,J) * USA_MASK(I,J)
         EPA_WD_BF_MEK (I,J) = EPA_WD_BF_MEK (I,J) * USA_MASK(I,J)
         EPA_WD_BF_PRPE(I,J) = EPA_WD_BF_PRPE(I,J) * USA_MASK(I,J)
         EPA_WD_BF_C3H8(I,J) = EPA_WD_BF_C3H8(I,J) * USA_MASK(I,J)
         EPA_WD_BF_CH2O(I,J) = EPA_WD_BF_CH2O(I,J) * USA_MASK(I,J)    
         EPA_WD_BF_C2H6(I,J) = EPA_WD_BF_C2H6(I,J) * USA_MASK(I,J)
         EPA_WD_BF_SO2 (I,J) = EPA_WD_BF_SO2 (I,J) * USA_MASK(I,J)
         EPA_WD_BF_SO4 (I,J) = EPA_WD_BF_SO4 (I,J) * USA_MASK(I,J)
         EPA_WD_BF_NH3 (I,J) = EPA_WD_BF_NH3 (I,J) * USA_MASK(I,J)

         ! Weekend avg biofuel    
         EPA_WE_BF_NOX (I,J) = EPA_WE_BF_NOX (I,J) * USA_MASK(I,J)
         EPA_WE_BF_CO  (I,J) = EPA_WE_BF_CO  (I,J) * USA_MASK(I,J)
         EPA_WE_BF_ALK4(I,J) = EPA_WE_BF_ALK4(I,J) * USA_MASK(I,J)
         EPA_WE_BF_ACET(I,J) = EPA_WE_BF_ACET(I,J) * USA_MASK(I,J)
         EPA_WE_BF_MEK (I,J) = EPA_WE_BF_MEK (I,J) * USA_MASK(I,J)
         EPA_WE_BF_PRPE(I,J) = EPA_WE_BF_PRPE(I,J) * USA_MASK(I,J)   
         EPA_WE_BF_C3H8(I,J) = EPA_WE_BF_C3H8(I,J) * USA_MASK(I,J)
         EPA_WE_BF_CH2O(I,J) = EPA_WE_BF_CH2O(I,J) * USA_MASK(I,J)
         EPA_WE_BF_C2H6(I,J) = EPA_WE_BF_C2H6(I,J) * USA_MASK(I,J)
         EPA_WE_BF_SO2 (I,J) = EPA_WE_BF_SO2 (I,J) * USA_MASK(I,J)
         EPA_WE_BF_SO4 (I,J) = EPA_WE_BF_SO4 (I,J) * USA_MASK(I,J)
         EPA_WE_BF_NH3 (I,J) = EPA_WE_BF_NH3 (I,J) * USA_MASK(I,J)

         !----------------------------------------------
         ! Compute IPCC future emissions (if necessary)
         !----------------------------------------------
         IF ( LFUTURE .and. USA_MASK(I,J) > 0d0 ) THEN

            ! Future biofuel scale factors            
            CObf   = GET_FUTURE_SCALE_CObf(   I, J )
            NH3bf  = GET_FUTURE_SCALE_NH3bf(  I, J )
            NOxbf  = GET_FUTURE_SCALE_NOXbf(  I, J )
            SO2bf  = GET_FUTURE_SCALE_SO2bf(  I, J )
            VOCbf  = GET_FUTURE_SCALE_VOCbf(  I, J )
                              
            ! Weekday avg biofuel    
            EPA_WD_BF_NOX (I,J) = EPA_WD_BF_NOX (I,J) * NOxbf
            EPA_WD_BF_CO  (I,J) = EPA_WD_BF_CO  (I,J) * CObf
            EPA_WD_BF_ALK4(I,J) = EPA_WD_BF_ALK4(I,J) * VOCbf
            EPA_WD_BF_ACET(I,J) = EPA_WD_BF_ACET(I,J) * VOCbf
            EPA_WD_BF_MEK (I,J) = EPA_WD_BF_MEK (I,J) * VOCbf
            EPA_WD_BF_PRPE(I,J) = EPA_WD_BF_PRPE(I,J) * VOCbf 
            EPA_WD_BF_C3H8(I,J) = EPA_WD_BF_C3H8(I,J) * VOCbf 
            EPA_WD_BF_CH2O(I,J) = EPA_WD_BF_CH2O(I,J) * VOCbf     
            EPA_WD_BF_C2H6(I,J) = EPA_WD_BF_C2H6(I,J) * VOCbf 
            EPA_WD_BF_SO2 (I,J) = EPA_WD_BF_SO2 (I,J) * SO2bf 
            EPA_WD_BF_SO4 (I,J) = EPA_WD_BF_SO4 (I,J) * SO2bf
            EPA_WD_BF_NH3 (I,J) = EPA_WD_BF_NH3 (I,J) * NH3bf
            
            ! Weekend avg biofuel    
            EPA_WE_BF_NOX (I,J) = EPA_WE_BF_NOX (I,J) * NOxbf
            EPA_WE_BF_CO  (I,J) = EPA_WE_BF_CO  (I,J) * CObf
            EPA_WE_BF_ALK4(I,J) = EPA_WE_BF_ALK4(I,J) * VOCbf 
            EPA_WE_BF_ACET(I,J) = EPA_WE_BF_ACET(I,J) * VOCbf 
            EPA_WE_BF_MEK (I,J) = EPA_WE_BF_MEK (I,J) * VOCbf 
            EPA_WE_BF_PRPE(I,J) = EPA_WE_BF_PRPE(I,J) * VOCbf   
            EPA_WE_BF_C3H8(I,J) = EPA_WE_BF_C3H8(I,J) * VOCbf 
            EPA_WE_BF_CH2O(I,J) = EPA_WE_BF_CH2O(I,J) * VOCbf 
            EPA_WE_BF_C2H6(I,J) = EPA_WE_BF_C2H6(I,J) * VOCbf 
            EPA_WE_BF_SO2 (I,J) = EPA_WE_BF_SO2 (I,J) * SO2bf
            EPA_WE_BF_SO4 (I,J) = EPA_WE_BF_SO4 (I,J) * SO2bf
            EPA_WE_BF_NH3 (I,J) = EPA_WE_BF_NH3 (I,J) * NH3bf
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Print totals to the log file
      CALL TOTAL_BIOFUEL_TG( THISMONTH )

      END SUBROUTINE EMISS_EPA_NEI_BF
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_epa
!
! !DESCRIPTION: Subroutine READ\_EPA reads an EPA data file (biomass or anthro)
!  from disk.  The entire file is read through on one pass for better I/O 
!  optimization.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EPA( FILENAME, NOX,  CO,   ALK4, ACET, MEK,       
     &                     PRPE,     C3H8, CH2O, C2H6, SO2,  SO4,  NH3 )
!
! !USES:
!
      USE BPCH2_MOD,        ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,         ONLY : IU_FILE, IOERROR
      USE LOGICAL_MOD,      ONLY : LICARTT
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      USE TIME_MOD,         ONLY : GET_YEAR
   
#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_O3"               ! FSCLYR
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME            ! File to read
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      TYPE (XPLEX),           INTENT(INOUT) :: NOX(IIPAR,JJPAR)    ! NOx  data
      TYPE (XPLEX),           INTENT(INOUT) :: CO(IIPAR,JJPAR)     ! CO   data
      TYPE (XPLEX),           INTENT(INOUT) :: ALK4(IIPAR,JJPAR)   ! ALK4 data
      TYPE (XPLEX),           INTENT(INOUT) :: ACET(IIPAR,JJPAR)   ! ACET data
      TYPE (XPLEX),           INTENT(INOUT) :: MEK(IIPAR,JJPAR)    ! MEK  data
      TYPE (XPLEX),           INTENT(INOUT) :: PRPE(IIPAR,JJPAR)   ! PRPE data  
      TYPE (XPLEX),           INTENT(INOUT) :: C3H8(IIPAR,JJPAR)   ! C3H8 data
      TYPE (XPLEX),           INTENT(INOUT) :: CH2O(IIPAR,JJPAR)   ! CH2O data
      TYPE (XPLEX),           INTENT(INOUT) :: C2H6(IIPAR,JJPAR)   ! C2H6 data
      TYPE (XPLEX),           INTENT(INOUT) :: SO2(IIPAR,JJPAR)    ! SO2  data
      TYPE (XPLEX),           INTENT(INOUT) :: SO4(IIPAR,JJPAR)    ! SO4  data
      TYPE (XPLEX),           INTENT(INOUT) :: NH3(IIPAR,JJPAR)    ! NH3  data
!
! !REVISION HISTORY: 
!  01 Jul 2004 - R. Hudman   - Initial version
!  (1 ) now apply yearly scale factor (amv, phs, 3/10/08)
!  (2 ) Now accounts for FSCLYR (phs, 3/17/08)
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,  J,  L,  N,  IOS
      INTEGER           :: NTRACER,   NSKIP
      INTEGER           :: HALFPOLAR, CENTER180
      INTEGER           :: NI,        NJ,        NL
      INTEGER           :: IFIRST,    JFIRST,    LFIRST
      INTEGER           :: SCALEYEAR, BASEYEAR
      TYPE (XPLEX)            :: ARRAY(IGLOB,JGLOB,1)
      REAL*4            :: D_ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)            :: LONRES,    LATRES
      REAL*4            :: D_LONRES, D_LATRES
      TYPE (XPLEX)            :: ZTAU0,     ZTAU1
      REAL*8             :: D_ZTAU0, D_ZTAU1
      TYPE (XPLEX)            :: SC(IIPAR,JJPAR)
      CHARACTER(LEN=20) :: MODELNAME
      CHARACTER(LEN=40) :: CATEGORY
      CHARACTER(LEN=40) :: UNIT     
      CHARACTER(LEN=40) :: RESERVED

      !=================================================================
      ! READ_EPA begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_EPA: Reading ', a )

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      ! Read the entire file in one pass (for I/O optimization)
      DO 

         ! Read 1st data block header line
         READ( IU_FILE, IOSTAT=IOS ) 
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LATRES = (D_LATRES)
         LONRES = (D_LONRES)
         ! Check for EOF or errors
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:2' )

         ! Read 2nd data block header line
         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0 = (D_ZTAU0)
         ZTAU1 = (D_ZTAU1)
         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:3' )

         ! Read data
         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( D_ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         ARRAY(:,:,:) = (D_ARRAY(:,:,:))
         ! Error check 
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:4' )

         !==============================================================
         ! Save into tracer arrays
         !==============================================================
         SELECT CASE ( NTRACER )
            CASE( 1  )
               CALL TRANSFER_2D( ARRAY(:,:,1), NOx  )
            CASE( 4  )
               CALL TRANSFER_2D( ARRAY(:,:,1), CO   )
            CASE( 5  )
               CALL TRANSFER_2D( ARRAY(:,:,1), ALK4 )
            CASE( 9  )
               CALL TRANSFER_2D( ARRAY(:,:,1), ACET )
            CASE( 10 )
               CALL TRANSFER_2D( ARRAY(:,:,1), MEK  )
            CASE( 18 )
               CALL TRANSFER_2D( ARRAY(:,:,1), PRPE )
            CASE( 19 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C3H8 )
            CASE( 20 )
               CALL TRANSFER_2D( ARRAY(:,:,1), CH2O )
            CASE( 21 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C2H6 )
            CASE( 26 )
               CALL TRANSFER_2D( ARRAY(:,:,1), SO2  )
            CASE( 27 )
               CALL TRANSFER_2D( ARRAY(:,:,1), SO4  )
            CASE( 29 )
               CALL TRANSFER_2D( ARRAY(:,:,1), NH3  )
            CASE DEFAULT
               ! Nothing
         END SELECT

      ENDDO

      ! Close file
      CLOSE( IU_FILE )


      ! Get/Apply annual scalar factor (amv, 08/24/07)
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      BASEYEAR = 1999
      IF ( LICARTT ) BASEYEAR = 2004
      
      CALL GET_ANNUAL_SCALAR( 71, BASEYEAR, SCALEYEAR, SC)
      NOx(:,:) = NOx(:,:) * SC(:,:)

      CALL GET_ANNUAL_SCALAR( 72, BASEYEAR, SCALEYEAR, SC)
      CO(:,:) = CO(:,:) * SC(:,:)

      CALL GET_ANNUAL_SCALAR( 73, 1999, SCALEYEAR, SC)
      SO2(:,:) = SO2(:,:) * SC(:,:)

      CALL GET_ANNUAL_SCALAR( 73, 1999, SCALEYEAR, SC)
      SO4(:,:) = SO4(:,:) * SC(:,:)

      END SUBROUTINE READ_EPA
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_epa_05x0666
!
! !DESCRIPTION: Subroutine READ\_EPA reads an EPA data file (biomass or anthro,
!  0.5 x 0.666 resolution) from disk.  The entire file is read through on one 
!  pass for better I/O  optimization.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EPA_05x0666( FILENAME, NOX,  CO,   ALK4, ACET, 
     &                             MEK,      PRPE, C3H8, CH2O, C2H6, 
     &                             SO2,      SO4,  NH3,  BF )

!
! !USES:
!
      USE BPCH2_MOD,        ONLY : READ_BPCH2, GET_TAU0
      USE LOGICAL_MOD,      ONLY : LICARTT
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      USE TIME_MOD,         ONLY : GET_YEAR, GET_MONTH

#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_O3"               ! FSCLYR
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME            ! File to read
      INTEGER,          INTENT(IN)    :: BF                  ! Read biofuels
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      TYPE (XPLEX),           INTENT(INOUT) :: NOX(IIPAR,JJPAR)    ! NOx  data 
      TYPE (XPLEX),           INTENT(INOUT) :: CO(IIPAR,JJPAR)     ! CO   data
      TYPE (XPLEX),           INTENT(INOUT) :: ALK4(IIPAR,JJPAR)   ! ALK4 data
      TYPE (XPLEX),           INTENT(INOUT) :: ACET(IIPAR,JJPAR)   ! ACET data
      TYPE (XPLEX),           INTENT(INOUT) :: MEK(IIPAR,JJPAR)    ! MEK  data
      TYPE (XPLEX),           INTENT(INOUT) :: PRPE(IIPAR,JJPAR)   ! PRPE data
      TYPE (XPLEX),           INTENT(INOUT) :: C3H8(IIPAR,JJPAR)   ! C3H8 data
      TYPE (XPLEX),           INTENT(INOUT) :: CH2O(IIPAR,JJPAR)   ! CH2O data
      TYPE (XPLEX),           INTENT(INOUT) :: C2H6(IIPAR,JJPAR)   ! C2H6 data
      TYPE (XPLEX),           INTENT(INOUT) :: SO2(IIPAR,JJPAR)    ! SO2  data
      TYPE (XPLEX),           INTENT(INOUT) :: SO4(IIPAR,JJPAR)    ! SO4  data
      TYPE (XPLEX),           INTENT(INOUT) :: NH3(IIPAR,JJPAR)    ! NH3  data
!
! !REMARKS:
!  Modified for nested grids (A. van Donkelaar)
! 
! !REVISION HISTORY: 
!  (1 ) now apply yearly scale factor (amv, phs, 3/10/08)
!  (2 ) Now accounts for FSCLYR (phs, 3/17/08)
!  07 Feb 2011 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,  J,  L,  N,  IOS
      INTEGER           :: NTRACER,   NSKIP
      INTEGER           :: HALFPOLAR, CENTER180
      INTEGER           :: SCALEYEAR, BASEYEAR
      TYPE (XPLEX)            :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)            :: LONRES,    LATRES
      TYPE (XPLEX)            :: ZTAU0,     ZTAU1, XTAU
      TYPE (XPLEX)            :: SC(IIPAR,JJPAR)
      CHARACTER(LEN=20) :: MODELNAME
      CHARACTER(LEN=40) :: CATEGORY
      CHARACTER(LEN=40) :: UNIT
      CHARACTER(LEN=40) :: RESERVED
      CHARACTER(LEN=8)  :: CAT

      !=================================================================
      ! READ_EPA begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_EPA: Reading ', a )

      XTAU = GET_TAU0(GET_MONTH(), 1, 1999)

      CAT = 'ANTHSRCE'
      IF (BF .EQ. 1) CAT = 'BIOFSRCE'

      CALL READ_BPCH2( FILENAME, CAT, 1, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      NOx(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT, 4, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      CO(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT, 5, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      ALK4(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT, 9, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      ACET(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,10, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      MEK(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,18, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      PRPE(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,19, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      C3H8(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,21, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      C2H6(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,26, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      SO2(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,27, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      SO4(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, CAT,29, XTAU,
     &    IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      NH3(:,:) = ARRAY(:,:,1)

      ! Get/Apply annual scalar factor (amv, 08/24/07)
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      BASEYEAR = 1999
      IF ( LICARTT ) BASEYEAR = 2004

      CALL GET_ANNUAL_SCALAR( 71, BASEYEAR, SCALEYEAR, SC)
      NOx(:,:) = NOx(:,:) * SC(:,:)

      CALL GET_ANNUAL_SCALAR( 72, BASEYEAR, SCALEYEAR, SC)
      CO(:,:) = CO(:,:) * SC(:,:)

      CALL GET_ANNUAL_SCALAR( 73, 1999, SCALEYEAR, SC)
      SO2(:,:) = SO2(:,:) * SC(:,:)

      CALL GET_ANNUAL_SCALAR( 73, 1999, SCALEYEAR, SC)
      SO4(:,:) = SO4(:,:) * SC(:,:)

      END SUBROUTINE READ_EPA_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_usa_mask
!
! !DESCRIPTION: Subroutine READ\_USA\_MASK reads the USA mask from disk.   The USA mask is
!  the fraction of the grid box (I,J) which lies w/in the continental USA.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_USA_MASK

!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE LOGICAL_MOD,   ONLY : LCAC,            LBRAVO
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
!
! !REMARKS:
!  From README file:
!    "The mask files have been updated to deal correctly with the boxes on
!     the border. They are simply set to include boxes where EPA emissions
!     are non-zero. We assume that BRAVO and CAC are also used when EPA is
!     used. If you use neither BRAVO nor CAC, you should use the older masks
!     to avoid missing emissions along the borders. Masks for the
!     case that either BRAVO or CAC (but not both) is used w/ EPA have not
!     been produced" (phs, 12/23/08)
! 
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman - Initial version
!  (1 ) Now can read data for GEOS and GCAP grids (bmy, 8/16/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Read larger mask for correct overlap with BRAVO and CAC, if these
!       regional inventories are used. 
!  (4 ) Temporary fix (until larger masks at 1x1 and 0.5x0.667, and cut to
!       the NA window, are available) nested NA runs.
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)             :: ARRAY(IGLOB,JGLOB,1)
      TYPE (XPLEX)             :: XTAU
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=14 ) :: MASK_DIR

      !=================================================================
      ! READ_USA_MASK begins here!
      !=================================================================

      ! File name
      IF ( LCAC .AND. LBRAVO ) THEN
         MASK_DIR = 'EPA_NEI_200806'
      ELSE IF ( LCAC .OR. LBRAVO ) THEN
         WRITE( 6, * ) '!! WARNING !! : WITH EPA, EITHER NONE OR ' //
     $        ' BOTH BRAVO AND CAC INVENTORIES SHOULD BE USED !!!'
         MASK_DIR = 'EPA_NEI_200806'
      ELSE
         MASK_DIR = 'EPA_NEI_200411'
      ENDIF

#if   defined( NESTED_NA )
      MASK_DIR = 'EPA_NEI_200411'
#endif

      ! Create filename
      FILENAME = TRIM( DATA_DIR )         //
     &           MASK_DIR // '/usa_mask.' // GET_NAME_EXT_2D() //
     &           '.'                      // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_USA_MASK: Reading ', a )

      ! Get TAU0 for Jan 1985
      XTAU  = GET_TAU0( 1, 1, 1985 )

      ! USA mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 XTAU,      IGLOB,    JGLOB,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to TYPE (XPLEX)
      CALL TRANSFER_2D( ARRAY(:,:,1), USA_MASK )

      END SUBROUTINE READ_USA_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_usa_mask_05x0666
!
! !DESCRIPTION: Subroutine READ\_USA\_MASK\_05x0666 reads the USA mask from 
!  disk, at 0.5 x 0.666 resolution.  The USA mask is the fraction of the grid 
!  box (I,J) which lies w/in the continental USA.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_USA_MASK_05x0666
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE LOGICAL_MOD,   ONLY : LCAC,            LBRAVO
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters
!
! !REMARKS:
!  From README file:
!    "The mask files have been updated to deal correctly with the boxes on
!     the border. They are simply set to include boxes where EPA emissions
!     are non-zero. We assume that BRAVO and CAC are also used when EPA is
!     used. If you use neither BRAVO nor CAC, you should use the older masks
!     to avoid missing emissions along the borders. Masks for the
!     case that either BRAVO or CAC (but not both) is used w/ EPA have not
!     been produced" (phs, 12/23/08) 
! 
! !REVISION HISTORY: 
!  (1 ) Now can read data for GEOS and GCAP grids (bmy, 8/16/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Read larger mask for correct overlap with BRAVO and CAC, if these
!       regional inventories are used.
!  (4 ) Temporary fix (until larger masks at 1x1 and 0.5x0.667, and cut to
!       the NA window, are available) nested NA runs.
!  07 Feb 2011 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE (XPLEX)             :: ARRAY(IIPAR,JJPAR,1)
      TYPE (XPLEX)             :: XTAU
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=14 ) :: MASK_DIR

      !=================================================================
      ! READ_USA_MASK begins here!
      !=================================================================

      ! File name
      IF ( LCAC .AND. LBRAVO ) THEN
         MASK_DIR = 'EPA_NEI_200806'
      ELSE IF ( LCAC .OR. LBRAVO ) THEN
         WRITE( 6, * ) '!! WARNING !! : WITH EPA, EITHER NONE OR ' //
     $        ' BOTH BRAVO AND CAC INVENTORIES SHOULD BE USED !!!'
         MASK_DIR = 'EPA_NEI_200806'
      ELSE
         MASK_DIR = 'EPA_NEI_200411'
      ENDIF

#if   defined( NESTED_NA )
      MASK_DIR = 'EPA_NEI_200411'
#endif
      
      ! Create filename
      FILENAME = TRIM( DATA_DIR )         //
     &           MASK_DIR // '/usa_mask.' // GET_NAME_EXT_2D() //
     &           '.'                      // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_USA_MASK: Reading ', a )

      ! Get TAU0 for Jan 1985
      XTAU  = GET_TAU0( 1, 1, 1985 )

      ! USA mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to TYPE (XPLEX)
      CALL TRANSFER_2D( ARRAY(:,:,1), USA_MASK )

      END SUBROUTINE READ_USA_MASK_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the amount of EPA/NEI99
!  anthropogenic emissions that are emitted each month in Tg or Tg C.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_Tg( THISMONTH )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"                 ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISMONTH   ! Current month (1-12)
! 
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  (1) Scale factors were determined by Jennifer Logan (jal@io.harvard.edu),
!      Bryan Duncan (bnd@io.harvard.edu), and Daniel Jacob (djj@io.harvard.edu)
!  (2) Now replace DXYP(J)*1d4 with routine GET_AREA_CM2 from "grid_mod.f".
!       (bmy, 2/4/03)
!  (3) Prevent out of bounds error when tracers are undefined (bmy, 1/25/05)
!  (4) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5) Now replace FMOL with TRACER_MW_KG (bmy, 10/25/05) 

!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J
      TYPE (XPLEX)           :: WD_NOX,  WD_CO,   WD_ALK4, WD_ACET
      TYPE (XPLEX)           :: WD_MEK,  WD_PRPE, WD_C2H6, WD_C3H8
      TYPE (XPLEX)           :: WD_CH2O, WD_NH3,  WD_SO2,  WD_SO4,  A
      TYPE (XPLEX)           :: WE_NOX,  WE_CO,   WE_ALK4, WE_ACET
      TYPE (XPLEX)           :: WE_MEK,  WE_PRPE, WE_C2H6, WE_C3H8
      TYPE (XPLEX)           :: WE_CH2O, WE_NH3,  WE_SO2,  WE_SO4
      TYPE (XPLEX)           :: F_NOX,   F_CO,    F_ALK4,  F_ACET
      TYPE (XPLEX)           :: F_MEK,   F_PRPE,  F_C2H6,  F_C3H8
      TYPE (XPLEX)           :: F_CH2O,  F_SO2,   F_SO4,   F_NH3
      CHARACTER(LEN=6) :: UNIT

      ! Days per month
      INTEGER          :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                               31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Summing variables for weekday avg anthro
      WD_NOX  = 0d0 
      WD_CO   = 0d0 
      WD_ALK4 = 0d0 
      WD_ACET = 0d0 
      WD_MEK  = 0d0 
      WD_PRPE = 0d0 
      WD_C2H6 = 0d0 
      WD_C3H8 = 0d0 
      WD_CH2O = 0d0 
      WD_NH3  = 0d0 
      WD_SO2  = 0d0 
      WD_SO4  = 0d0 

      ! Summing variables for weekend avg anthro
      WE_NOX  = 0d0 
      WE_CO   = 0d0 
      WE_ALK4 = 0d0 
      WE_ACET = 0d0 
      WE_MEK  = 0d0 
      WE_PRPE = 0d0 
      WE_C2H6 = 0d0 
      WE_C3H8 = 0d0 
      WE_CH2O = 0d0 
      WE_NH3  = 0d0 
      WE_SO2  = 0d0 
      WE_SO4  = 0d0 

      ! Molecular weights
      F_NOX   = 0d0   
      F_CO    = 0d0 
      F_ALK4  = 0d0 
      F_ACET  = 0d0 
      F_MEK   = 0d0 
      F_PRPE  = 0d0 
      F_C2H6  = 0d0 
      F_C3H8  = 0d0 
      F_CH2O  = 0d0 
      F_SO2   = 0d0 
      F_SO4   = 0d0 
      F_NH3   = 0d0 

      ! Prevent array out of bounds error for undefined tracers
      IF ( IDTNOX  > 0 ) F_NOX  = TRACER_MW_KG(IDTNOX )
      IF ( IDTCO   > 0 ) F_CO   = TRACER_MW_KG(IDTCO  )
      IF ( IDTALK4 > 0 ) F_ALK4 = TRACER_MW_KG(IDTALK4)
      IF ( IDTACET > 0 ) F_ACET = TRACER_MW_KG(IDTACET)
      IF ( IDTMEK  > 0 ) F_MEK  = TRACER_MW_KG(IDTMEK )
      IF ( IDTPRPE > 0 ) F_PRPE = TRACER_MW_KG(IDTPRPE)
      IF ( IDTC2H6 > 0 ) F_C2H6 = TRACER_MW_KG(IDTC2H6)
      IF ( IDTC3H8 > 0 ) F_C3H8 = TRACER_MW_KG(IDTC3H8)
      IF ( IDTCH2O > 0 ) F_CH2O = TRACER_MW_KG(IDTCH2O)
      IF ( IDTSO2  > 0 ) F_SO2  = TRACER_MW_KG(IDTSO2 )
      IF ( IDTSO4  > 0 ) F_SO4  = TRACER_MW_KG(IDTSO4 )
      IF ( IDTNH3  > 0 ) F_NH3  = TRACER_MW_KG(IDTNH3 )

      !=================================================================
      ! Sum anthropogenic emissions
      !=================================================================

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Weekday avg emissions
            WD_NOX  = WD_NOX  + EPA_WD_AN_NOX (I,J) * A * F_NOX 
            WD_CO   = WD_CO   + EPA_WD_AN_CO  (I,J) * A * F_CO  
            WD_ALK4 = WD_ALK4 + EPA_WD_AN_ALK4(I,J) * A * F_ALK4
            WD_ACET = WD_ACET + EPA_WD_AN_ACET(I,J) * A * F_ACET
            WD_MEK  = WD_MEK  + EPA_WD_AN_MEK (I,J) * A * F_MEK 
            WD_PRPE = WD_PRPE + EPA_WD_AN_PRPE(I,J) * A * F_PRPE
            WD_C2H6 = WD_C2H6 + EPA_WD_AN_C2H6(I,J) * A * F_C2H6
            WD_C3H8 = WD_C3H8 + EPA_WD_AN_C3H8(I,J) * A * F_C3H8
            WD_CH2O = WD_CH2O + EPA_WD_AN_CH2O(I,J) * A * F_CH2O
            WD_SO2  = WD_SO2  + EPA_WD_AN_SO2 (I,J) * A * F_SO2 
            WD_SO4  = WD_SO4  + EPA_WD_AN_SO4 (I,J) * A * F_SO4 
            WD_NH3  = WD_NH3  + EPA_WD_AN_NH3 (I,J) * A * F_NH3 

            ! Weekend avg emissions
            WE_NOX  = WE_NOX  + EPA_WE_AN_NOX (I,J) * A * F_NOX 
            WE_CO   = WE_CO   + EPA_WE_AN_CO  (I,J) * A * F_CO  
            WE_ALK4 = WE_ALK4 + EPA_WE_AN_ALK4(I,J) * A * F_ALK4
            WE_ACET = WE_ACET + EPA_WE_AN_ACET(I,J) * A * F_ACET
            WE_MEK  = WE_MEK  + EPA_WE_AN_MEK (I,J) * A * F_MEK 
            WE_PRPE = WE_PRPE + EPA_WE_AN_PRPE(I,J) * A * F_PRPE
            WE_C2H6 = WE_C2H6 + EPA_WE_AN_C2H6(I,J) * A * F_C2H6
            WE_C3H8 = WE_C3H8 + EPA_WE_AN_C3H8(I,J) * A * F_C3H8
            WE_CH2O = WE_CH2O + EPA_WE_AN_CH2O(I,J) * A * F_CH2O
            WE_SO2  = WE_SO2  + EPA_WE_AN_SO2 (I,J) * A * F_SO2 
            WE_SO4  = WE_SO4  + EPA_WE_AN_SO4 (I,J) * A * F_SO4 
            WE_NH3  = WE_NH3  + EPA_WE_AN_NH3 (I,J) * A * F_NH3 

         ENDDO
      ENDDO
 
      !=================================================================
      ! Print info
      !=================================================================
      
      ! Weekday avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 100   ) 'NOx ', THISMONTH, WD_NOX,  '  '
      WRITE( 6, 100   ) 'CO  ', THISMONTH, WD_CO,   '  '
      WRITE( 6, 100   ) 'ALK4', THISMONTH, WD_ALK4, ' C'
      WRITE( 6, 100   ) 'ACET', THISMONTH, WD_ACET, ' C'
      WRITE( 6, 100   ) 'MEK ', THISMONTH, WD_MEK,  ' C'
      WRITE( 6, 100   ) 'PRPE', THISMONTH, WD_PRPE, ' C'
      WRITE( 6, 100   ) 'C3H8', THISMONTH, WD_C3H8, ' C'
      WRITE( 6, 100   ) 'CH2O', THISMONTH, WD_CH2O, '  '
      WRITE( 6, 100   ) 'C2H6', THISMONTH, WD_C2H6, ' C'
      WRITE( 6, 100   ) 'SO2 ', THISMONTH, WD_SO2,  '  '
      WRITE( 6, 100   ) 'SO4 ', THISMONTH, WD_SO4,  '  '
      WRITE( 6, 100   ) 'NH3 ', THISMONTH, WD_NH3,  '  '
 100  FORMAT( 'Total weekday avg anthro ', a4, ' for 1999/', 
     &         i2.2, ': ', 2f13.6, ' Tg', a2 )

      ! Weekend avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NOx ', THISMONTH, WE_NOX,  '  '
      WRITE( 6, 110   ) 'CO  ', THISMONTH, WE_CO,   '  '
      WRITE( 6, 110   ) 'ALK4', THISMONTH, WE_ALK4, ' C'
      WRITE( 6, 110   ) 'ACET', THISMONTH, WE_ACET, ' C'
      WRITE( 6, 110   ) 'MEK ', THISMONTH, WE_MEK,  ' C'
      WRITE( 6, 110   ) 'PRPE', THISMONTH, WE_PRPE, ' C'
      WRITE( 6, 110   ) 'C3H8', THISMONTH, WE_C3H8, ' C'
      WRITE( 6, 110   ) 'CH2O', THISMONTH, WE_CH2O, '  '
      WRITE( 6, 110   ) 'C2H6', THISMONTH, WE_C2H6, ' C'
      WRITE( 6, 110   ) 'SO2 ', THISMONTH, WE_SO2,  '  '
      WRITE( 6, 110   ) 'SO4 ', THISMONTH, WE_SO4,  '  '
      WRITE( 6, 110   ) 'NH3 ', THISMONTH, WE_NH3,  '  '
 110  FORMAT( 'Total weekend avg anthro ', a4, ' for 1999/', 
     &         i2.2, ': ', 2f13.6, ' Tg', a2 )

      END SUBROUTINE TOTAL_ANTHRO_TG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_biofuel_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_BIOFUEL\_Tg prints the amount of EPA/NEI99
!  biofuel emissions that are emitted each month in Tg or Tg C.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_BIOFUEL_Tg( THISMONTH )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"  ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISMONTH   ! Current month (1-12)
! 
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  (1 ) Prevent out of bounds error when tracers are undefined (bmy, 1/25/05)
!  (2 ) Now replace FMOL with TRACER_MW_KG
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J
      TYPE (XPLEX)           :: WD_NOX,  WD_CO,   WD_ALK4, WD_ACET
      TYPE (XPLEX)           :: WD_MEK,  WD_PRPE, WD_C2H6, WD_C3H8
      TYPE (XPLEX)           :: WD_CH2O, WD_NH3,  WD_SO2,  WD_SO4,  A
      TYPE (XPLEX)           :: WE_NOX,  WE_CO,   WE_ALK4, WE_ACET
      TYPE (XPLEX)           :: WE_MEK,  WE_PRPE, WE_C2H6, WE_C3H8
      TYPE (XPLEX)           :: WE_CH2O, WE_NH3,  WE_SO2,  WE_SO4
      TYPE (XPLEX)           :: F_NOX,   F_CO,    F_ALK4,  F_ACET
      TYPE (XPLEX)           :: F_MEK,   F_PRPE,  F_C2H6,  F_C3H8
      TYPE (XPLEX)           :: F_CH2O,  F_SO2,   F_SO4,   F_NH3
      CHARACTER(LEN=6) :: UNIT

      ! Days per month
      INTEGER          :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                               31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_BIOFUEL_TG begins here!
      !=================================================================

      ! Summing variables for weekday avg anthro
      WD_NOX  = 0d0 
      WD_CO   = 0d0 
      WD_ALK4 = 0d0 
      WD_ACET = 0d0 
      WD_MEK  = 0d0 
      WD_PRPE = 0d0 
      WD_C2H6 = 0d0 
      WD_C3H8 = 0d0 
      WD_CH2O = 0d0 
      WD_NH3  = 0d0 
      WD_SO2  = 0d0 
      WD_SO4  = 0d0 

      ! Summing variables for weekend avg anthro
      WE_NOX  = 0d0 
      WE_CO   = 0d0 
      WE_ALK4 = 0d0 
      WE_ACET = 0d0 
      WE_MEK  = 0d0 
      WE_PRPE = 0d0 
      WE_C2H6 = 0d0 
      WE_C3H8 = 0d0 
      WE_CH2O = 0d0 
      WE_NH3  = 0d0 
      WE_SO2  = 0d0 
      WE_SO4  = 0d0 

      ! Molecular weights
      F_NOX   = 0d0   
      F_CO    = 0d0 
      F_ALK4  = 0d0 
      F_ACET  = 0d0 
      F_MEK   = 0d0 
      F_PRPE  = 0d0 
      F_C2H6  = 0d0 
      F_C3H8  = 0d0 
      F_CH2O  = 0d0 
      F_SO2   = 0d0 
      F_SO4   = 0d0 
      F_NH3   = 0d0 

      ! Prevent array out of bounds error for undefined tracers
      IF ( IDTNOX  > 0 ) F_NOX  = TRACER_MW_KG(IDTNOX )
      IF ( IDTCO   > 0 ) F_CO   = TRACER_MW_KG(IDTCO  )
      IF ( IDTALK4 > 0 ) F_ALK4 = TRACER_MW_KG(IDTALK4)
      IF ( IDTACET > 0 ) F_ACET = TRACER_MW_KG(IDTACET)
      IF ( IDTMEK  > 0 ) F_MEK  = TRACER_MW_KG(IDTMEK )
      IF ( IDTPRPE > 0 ) F_PRPE = TRACER_MW_KG(IDTPRPE)
      IF ( IDTC2H6 > 0 ) F_C2H6 = TRACER_MW_KG(IDTC2H6)
      IF ( IDTC3H8 > 0 ) F_C3H8 = TRACER_MW_KG(IDTC3H8)
      IF ( IDTCH2O > 0 ) F_CH2O = TRACER_MW_KG(IDTCH2O)
      IF ( IDTSO2  > 0 ) F_SO2  = TRACER_MW_KG(IDTSO2 )
      IF ( IDTSO4  > 0 ) F_SO4  = TRACER_MW_KG(IDTSO4 )
      IF ( IDTNH3  > 0 ) F_NH3  = TRACER_MW_KG(IDTNH3 )

      !=================================================================
      ! Sum biofuel emissions
      !=================================================================

      ! Loop over surface boxes
      DO J = 1, JJPAR

         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23

         DO I = 1, IIPAR

            ! Weekday avg emissions
            WD_NOX  = WD_NOX  + EPA_WD_BF_NOX (I,J) * A * F_NOX
            WD_CO   = WD_CO   + EPA_WD_BF_CO  (I,J) * A * F_CO  
            WD_ALK4 = WD_ALK4 + EPA_WD_BF_ALK4(I,J) * A * F_ALK4
            WD_ACET = WD_ACET + EPA_WD_BF_ACET(I,J) * A * F_ACET
            WD_MEK  = WD_MEK  + EPA_WD_BF_MEK (I,J) * A * F_MEK 
            WD_PRPE = WD_PRPE + EPA_WD_BF_PRPE(I,J) * A * F_PRPE
            WD_C2H6 = WD_C2H6 + EPA_WD_BF_C2H6(I,J) * A * F_C2H6
            WD_C3H8 = WD_C3H8 + EPA_WD_BF_C3H8(I,J) * A * F_C3H8
            WD_CH2O = WD_CH2O + EPA_WD_BF_CH2O(I,J) * A * F_CH2O
            WD_SO2  = WD_SO2  + EPA_WD_BF_SO2 (I,J) * A * F_SO2 
            WD_SO4  = WD_SO4  + EPA_WD_BF_SO4 (I,J) * A * F_SO4 
            WD_NH3  = WD_NH3  + EPA_WD_BF_NH3 (I,J) * A * F_NH3 
            
            ! Weekend avg emissions
            WE_NOX  = WE_NOX  + EPA_WE_BF_NOX (I,J) * A * F_NOX 
            WE_CO   = WE_CO   + EPA_WE_BF_CO  (I,J) * A * F_CO  
            WE_ALK4 = WE_ALK4 + EPA_WE_BF_ALK4(I,J) * A * F_ALK4
            WE_ACET = WE_ACET + EPA_WE_BF_ACET(I,J) * A * F_ACET
            WE_MEK  = WE_MEK  + EPA_WE_BF_MEK (I,J) * A * F_MEK 
            WE_PRPE = WE_PRPE + EPA_WE_BF_PRPE(I,J) * A * F_PRPE
            WE_C2H6 = WE_C2H6 + EPA_WE_BF_C2H6(I,J) * A * F_C2H6
            WE_C3H8 = WE_C3H8 + EPA_WE_BF_C3H8(I,J) * A * F_C3H8
            WE_CH2O = WE_CH2O + EPA_WE_BF_CH2O(I,J) * A * F_CH2O
            WE_SO2  = WE_SO2  + EPA_WE_BF_SO2 (I,J) * A * F_SO2 
            WE_SO4  = WE_SO4  + EPA_WE_BF_SO4 (I,J) * A * F_SO4 
            WE_NH3  = WE_NH3  + EPA_WE_BF_NH3 (I,J) * A * F_NH3 
            
         ENDDO
      ENDDO
 
      !=================================================================
      ! Print info
      !=================================================================
      
      ! Weekday avg biofuel
      WRITE( 6, '(a)' )
      WRITE( 6, 100   ) 'NOx ', THISMONTH, WD_NOX,  '  '
      WRITE( 6, 100   ) 'CO  ', THISMONTH, WD_CO,   '  '
      WRITE( 6, 100   ) 'ALK4', THISMONTH, WD_ALK4, ' C'
      WRITE( 6, 100   ) 'ACET', THISMONTH, WD_ACET, ' C'
      WRITE( 6, 100   ) 'MEK ', THISMONTH, WD_MEK,  ' C'
      WRITE( 6, 100   ) 'PRPE', THISMONTH, WD_PRPE, ' C'
      WRITE( 6, 100   ) 'C3H8', THISMONTH, WD_C3H8, ' C'
      WRITE( 6, 100   ) 'CH2O', THISMONTH, WD_CH2O, '  '
      WRITE( 6, 100   ) 'C2H6', THISMONTH, WD_C2H6, ' C'
      WRITE( 6, 100   ) 'SO2 ', THISMONTH, WD_SO2,  '  '
      WRITE( 6, 100   ) 'SO4 ', THISMONTH, WD_SO4,  '  '
      WRITE( 6, 100   ) 'NH3 ', THISMONTH, WD_NH3,  '  '
 100  FORMAT( 'Total weekday avg biofuel ', a4, ' for 1999/', 
     &         i2.2, ': ', 2f13.6, ' Tg', a2 )

      ! Weekend avg biofuel
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NOx ', THISMONTH, WE_NOX,  '  '
      WRITE( 6, 110   ) 'CO  ', THISMONTH, WE_CO,   '  '
      WRITE( 6, 110   ) 'ALK4', THISMONTH, WE_ALK4, ' C'
      WRITE( 6, 110   ) 'ACET', THISMONTH, WE_ACET, ' C'
      WRITE( 6, 110   ) 'MEK ', THISMONTH, WE_MEK,  ' C'
      WRITE( 6, 110   ) 'PRPE', THISMONTH, WE_PRPE, ' C'
      WRITE( 6, 110   ) 'C3H8', THISMONTH, WE_C3H8, ' C'
      WRITE( 6, 110   ) 'CH2O', THISMONTH, WE_CH2O, '  '
      WRITE( 6, 110   ) 'C2H6', THISMONTH, WE_C2H6, ' C'
      WRITE( 6, 110   ) 'SO2 ', THISMONTH, WE_SO2,  '  '
      WRITE( 6, 110   ) 'SO4 ', THISMONTH, WE_SO4,  '  '
      WRITE( 6, 110   ) 'NH3 ', THISMONTH, WE_NH3,  '  '
 110  FORMAT( 'Total weekend avg biofuel ', a4, ' for 1999/', 
     &         i2.2, ': ', 2f13.6, ' Tg', a2 )

      END SUBROUTINE TOTAL_BIOFUEL_TG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_usa_mask
!
! !DESCRIPTION: Function GET\_USA\_MASK returns the value of the USA mask 
!  (i.e. the fraction of a grid box which lies w/in the continental USA) 
!  at a given (I,J) location. (rch, bmy, 11/10/04)

!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_USA_MASK( I, J ) RESULT( USA )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I     ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J     ! GEOS-Chem latitude index
!
! !RETURN VALUE:
!
      TYPE (XPLEX)              :: USA   ! Fraction of box in continental USA
!
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_USA_MASK begins here!
      !=================================================================
      USA = USA_MASK(I,J)
      
      END FUNCTION GET_USA_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_epa_anthro
!
! !DESCRIPTION: Function GET\_EPA\_ANTHRO returns the EPA/NEI99 weekday avg or 
!  weekend avg anthropogenic emissions at a (I,J) location.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EPA_ANTHRO( I, J, N, WEEKDAY ) RESULT( EPA_NEI )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"               ! Size parameters
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: WEEKDAY   ! =T if it's a weekday, =F if weekend
      INTEGER, INTENT(IN) :: I         ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J         ! GEOS-Chem latitude index
      INTEGER, INTENT(IN) :: N         ! GEOS-Chem tracer number
!
! !RETURN VALUE:
!
      TYPE (XPLEX)              :: EPA_NEI   ! Anthropogenic Emissions [molec/cm2/s] 
                                       !  or [atoms C/cm2/s] for some HC's
! 
! !REVISION HISTORY: 
!  07 Feb 2011 - R. Yantosca - Initial version
!  10 Nov 2004 - R. Hudman   - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY.  Also remove reference
!        to BPCH2_MOD and TRACERID_MOD, they're not needed.  (bmy, 10/3/05)
!  (2 ) Default value changed to -1 to identify tracers without EPA/NEI 
!        emissions. (hotp, ccc, 5/29/09)

!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_EPA_ANTHRO begins here!
      !=================================================================

      ! Return either weekday or weekend avg emissions
      IF ( WEEKDAY ) THEN 

         !--------------------
         ! Weekday avg anthro
         !--------------------
         IF ( N == IDTNOX ) THEN
            EPA_NEI = EPA_WD_AN_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WD_AN_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WD_AN_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WD_AN_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WD_AN_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WD_AN_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WD_AN_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WD_AN_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WD_AN_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WD_AN_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WD_AN_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WD_AN_NH3(I,J)

         ELSE
            ! Some anthropogenic species don't have EPA/NEI emissions. 
            ! We need to keep background emissions. (hotp, ccc, 5/29/09)
            EPA_NEI = -1d0
            
         ENDIF

      ELSE

         !--------------------
         ! Weekend avg anthro
         !--------------------
         IF (  N == IDTNOX  ) THEN
            EPA_NEI = EPA_WE_AN_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WE_AN_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WE_AN_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WE_AN_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WE_AN_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WE_AN_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WE_AN_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WE_AN_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WE_AN_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WE_AN_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WE_AN_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WE_AN_NH3(I,J)

         ELSE
            ! Some anthropogenic species don't have EPA/NEI emissions. 
            ! We need to keep background emissions. (hotp, ccc, 5/29/09)
            EPA_NEI = -1d0
            
         ENDIF

      ENDIF

      END FUNCTION GET_EPA_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_epa_biofuel
!
! !DESCRIPTION: Function GET\_EPA\_BIOFUEL returns the EPA/NEI99 weekday avg 
!  or weekend avg biofuel emissions at a (I,J) location. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EPA_BIOFUEL( I, J, N, WEEKDAY ) RESULT( EPA_NEI )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"               ! Size parameters
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: WEEKDAY   ! =T if it's a weekday, =F if weekend
      INTEGER, INTENT(IN) :: I         ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J         ! GEOS-Chem latitude index
      INTEGER, INTENT(IN) :: N         ! GEOS-Chem tracer number
!
! !RETURN VALUE:
!
      TYPE (XPLEX)              :: EPA_NEI
!
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY.  Also remove reference
!        to BPCH2_MOD and TRACERID_MOD, they're not needed.  (bmy, 10/3/05)
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_EPA_BIOFUEL begins here!
      !=================================================================

      ! Return either weekday or weekend avg emissions
      IF ( WEEKDAY ) THEN 

         !---------------------
         ! Weekday avg biofuel
         !---------------------
         IF ( N == IDTNOX ) THEN
            EPA_NEI = EPA_WD_BF_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WD_BF_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WD_BF_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WD_BF_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WD_BF_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WD_BF_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WD_BF_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WD_BF_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WD_BF_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WD_BF_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WD_BF_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WD_BF_NH3(I,J)

         ELSE
            EPA_NEI = 0d0
            
         ENDIF

      ELSE

         !---------------------
         ! Weekend avg biofuel
         !---------------------
         IF (  N == IDTNOX  ) THEN
            EPA_NEI = EPA_WE_BF_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WE_BF_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WE_BF_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WE_BF_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WE_BF_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WE_BF_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WE_BF_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WE_BF_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WE_BF_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WE_BF_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WE_BF_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WE_BF_NH3(I,J)

         ELSE
            EPA_NEI = 0d0
            
         ENDIF

      ENDIF

      END FUNCTION GET_EPA_BIOFUEL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_epa_nei
!
! !DESCRIPTION: Subroutine INIT\_EPA\_NEI allocates and zeroes all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_EPA_NEI
!
! !USES:
!
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LNEI99
      USE LOGICAL_MOD, ONLY : LNEI05

#     include "CMN_SIZE"  ! Size parameters
!
! !REMARKS:
!  Read EPA/NEI99 anthro + EPA/NEI99 biofuel if LNEI99=T, or
!  Read EPA/NEI05 anthro + EPA/NEI99 biofuel if LNEI99=F but LNEI05=T
! 
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      
      !=================================================================
      ! INIT_EPA_NEI begins here!
      !=================================================================

      !------------------------------------
      ! Prior to 2/7/11:
      ! Return if we LNEI99 = .FALSE.
      !IF ( .not. LNEI99 ) RETURN
      !------------------------------------

      ! Return if we have not selected LNEI99 (or LNEI05, since we
      ! must read the EPA/NEI99 biofuels when using EPA/NEI05 anthro)
      ! (bmy, 2/7/11)
      IF ( .not. LNEI99 ) THEN
         IF  ( .not. LNEI05 ) RETURN
      ENDIF

      ! USA Mask
      ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
      USA_MASK = 0d0

      !-----------------------
      ! Anthro - weekday avg
      !-----------------------
      IF ( .not. LNEI05 ) THEN 
   
         ALLOCATE( EPA_WD_AN_NOX( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_NOX' )
         EPA_WD_AN_NOX = 0d0
   
         ALLOCATE( EPA_WD_AN_CO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_CO' )
         EPA_WD_AN_CO = 0d0
   
         ALLOCATE( EPA_WD_AN_ALK4( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_ALK4' )
         EPA_WD_AN_ALK4 = 0d0
   
         ALLOCATE( EPA_WD_AN_ACET( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_ACET' )
         EPA_WD_AN_ACET = 0d0
   
         ALLOCATE( EPA_WD_AN_MEK( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_MEK' )
         EPA_WD_AN_MEK = 0d0
   
         ALLOCATE( EPA_WD_AN_PRPE( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_PRPE' )
         EPA_WD_AN_PRPE = 0d0
   
         ALLOCATE( EPA_WD_AN_C2H6( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_C2H6' )
         EPA_WD_AN_C2H6 = 0d0
   
         ALLOCATE( EPA_WD_AN_C3H8( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_C3H8' )
         EPA_WD_AN_C3H8 = 0d0
   
         ALLOCATE( EPA_WD_AN_CH2O( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_CH2O' )
         EPA_WD_AN_CH2O = 0d0
   
         ALLOCATE( EPA_WD_AN_NH3( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_NH3' )
         EPA_WD_AN_NH3 = 0d0
   
         ALLOCATE( EPA_WD_AN_SO2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_SO2' )
         EPA_WD_AN_SO2 = 0d0
   
         ALLOCATE( EPA_WD_AN_SO4( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_SO4' )
         EPA_WD_AN_SO4 = 0d0

      ENDIF

      !-----------------------
      ! Anthro - weekend avg
      !-----------------------
      IF ( .not. LNEI05 ) THEN

         ALLOCATE( EPA_WE_AN_NOX( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_NOX' )
         EPA_WE_AN_NOX = 0d0
   
         ALLOCATE( EPA_WE_AN_CO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_CO' )
         EPA_WE_AN_CO = 0d0
   
         ALLOCATE( EPA_WE_AN_ALK4( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_ALK4' )
         EPA_WE_AN_ALK4 = 0d0
   
         ALLOCATE( EPA_WE_AN_ACET( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_ACET' )
         EPA_WE_AN_ACET = 0d0
   
         ALLOCATE( EPA_WE_AN_MEK( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_MEK' )
         EPA_WE_AN_MEK = 0d0
   
         ALLOCATE( EPA_WE_AN_PRPE( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_PRPE' )
         EPA_WE_AN_PRPE = 0d0
   
         ALLOCATE( EPA_WE_AN_C2H6( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_C2H6' )
         EPA_WE_AN_C2H6 = 0d0
   
         ALLOCATE( EPA_WE_AN_C3H8( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_C3H8' )
         EPA_WE_AN_C3H8 = 0d0
   
         ALLOCATE( EPA_WE_AN_CH2O( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_CH2O' )
         EPA_WE_AN_CH2O = 0d0
   
         ALLOCATE( EPA_WE_AN_NH3( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_NH3' )
         EPA_WE_AN_NH3 = 0d0
   
         ALLOCATE( EPA_WE_AN_SO2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_SO2' )
         EPA_WE_AN_SO2 = 0d0
   
         ALLOCATE( EPA_WE_AN_SO4( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_SO4' )
         EPA_WE_AN_SO4 = 0d0

      ENDIF

      !-----------------------
      ! Biofuel - weekday avg
      !-----------------------
      ALLOCATE( EPA_WD_BF_NOX(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_NOX' )
      EPA_WD_BF_NOX = 0d0

      ALLOCATE( EPA_WD_BF_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_CO' )
      EPA_WD_BF_CO = 0d0

      ALLOCATE( EPA_WD_BF_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_ALK4' )
      EPA_WD_BF_ALK4 = 0d0

      ALLOCATE( EPA_WD_BF_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_ACET' )
      EPA_WD_BF_ACET = 0d0

      ALLOCATE( EPA_WD_BF_MEK(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_MEK' )
      EPA_WD_BF_MEK = 0d0

      ALLOCATE( EPA_WD_BF_PRPE( IIPAR, JJPAR ), STAT=AS )     
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_PRPE' )
      EPA_WD_BF_PRPE = 0d0

      ALLOCATE( EPA_WD_BF_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_C2H6' )
      EPA_WD_BF_C2H6 = 0d0

      ALLOCATE( EPA_WD_BF_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_C3H8' )
      EPA_WD_BF_C3H8 = 0d0

      ALLOCATE( EPA_WD_BF_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_CH2O' )
      EPA_WD_BF_CH2O = 0d0

      ALLOCATE( EPA_WD_BF_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_NH3' )
      EPA_WD_BF_NH3 = 0d0

      ALLOCATE( EPA_WD_BF_SO2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_SO2' )
      EPA_WD_BF_SO2 = 0d0

      ALLOCATE( EPA_WD_BF_SO4( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_SO4' )
      EPA_WD_BF_SO4 = 0d0

      !-----------------------
      ! Biofuel - weekend avg
      !-----------------------
      ALLOCATE( EPA_WE_BF_NOX(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_NOX' )
      EPA_WE_BF_NOX = 0d0

      ALLOCATE( EPA_WE_BF_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_CO' )
      EPA_WE_BF_CO = 0d0

      ALLOCATE( EPA_WE_BF_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_ALK4' )
      EPA_WE_BF_ALK4 = 0d0

      ALLOCATE( EPA_WE_BF_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_ACET' )
      EPA_WE_BF_ACET = 0d0

      ALLOCATE( EPA_WE_BF_MEK(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_MEK' )
      EPA_WE_BF_MEK = 0d0

      ALLOCATE( EPA_WE_BF_PRPE( IIPAR, JJPAR ), STAT=AS )     
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_PRPE' )
      EPA_WE_BF_PRPE = 0d0

      ALLOCATE( EPA_WE_BF_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_C2H6' )
      EPA_WE_BF_C2H6 = 0d0

      ALLOCATE( EPA_WE_BF_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_C3H8' )
      EPA_WE_BF_C3H8 = 0d0

      ALLOCATE( EPA_WE_BF_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_CH2O' )
      EPA_WE_BF_CH2O = 0d0

      ALLOCATE( EPA_WE_BF_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_NH3' )
      EPA_WE_BF_NH3 = 0d0

      ALLOCATE( EPA_WE_BF_SO2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_SO2' )
      EPA_WE_BF_SO2 = 0d0

      ALLOCATE( EPA_WE_BF_SO4( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_SO4' )
      EPA_WE_BF_SO4 = 0d0

      END SUBROUTINE INIT_EPA_NEI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_epa_nei
!
! !DESCRIPTION: Subroutine CLEANUP\_EPA\_NEI deallocates all module arrays 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_EPA_NEI
! 
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman   - Initial version
!  07 Feb 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_EPA_NEI begins here!
      !=================================================================

      ! USA mask
      IF ( ALLOCATED( USA_MASK       ) ) DEALLOCATE( USA_MASK       )

      ! Fossil fuel -- weekday
      IF ( ALLOCATED( EPA_WD_AN_NOX  ) ) DEALLOCATE( EPA_WD_AN_NOX  )
      IF ( ALLOCATED( EPA_WD_AN_CO   ) ) DEALLOCATE( EPA_WD_AN_CO   )
      IF ( ALLOCATED( EPA_WD_AN_ALK4 ) ) DEALLOCATE( EPA_WD_AN_ALK4 )
      IF ( ALLOCATED( EPA_WD_AN_ACET ) ) DEALLOCATE( EPA_WD_AN_ACET )
      IF ( ALLOCATED( EPA_WD_AN_MEK  ) ) DEALLOCATE( EPA_WD_AN_MEK  )
      IF ( ALLOCATED( EPA_WD_AN_PRPE ) ) DEALLOCATE( EPA_WD_AN_PRPE )
      IF ( ALLOCATED( EPA_WD_AN_C2H6 ) ) DEALLOCATE( EPA_WD_AN_C2H6 )
      IF ( ALLOCATED( EPA_WD_AN_C3H8 ) ) DEALLOCATE( EPA_WD_AN_C3H8 )
      IF ( ALLOCATED( EPA_WD_AN_C2H6 ) ) DEALLOCATE( EPA_WD_AN_C2H6 )
      IF ( ALLOCATED( EPA_WD_AN_NH3  ) ) DEALLOCATE( EPA_WD_AN_NH3  )
      IF ( ALLOCATED( EPA_WD_AN_SO2  ) ) DEALLOCATE( EPA_WD_AN_SO2  )
      IF ( ALLOCATED( EPA_WD_AN_SO4  ) ) DEALLOCATE( EPA_WD_AN_SO4  )

      ! Fossil fuel -- weekend
      IF ( ALLOCATED( EPA_WE_AN_NOX  ) ) DEALLOCATE( EPA_WE_AN_NOX  )
      IF ( ALLOCATED( EPA_WE_AN_CO   ) ) DEALLOCATE( EPA_WE_AN_CO   )
      IF ( ALLOCATED( EPA_WE_AN_ALK4 ) ) DEALLOCATE( EPA_WE_AN_ALK4 )
      IF ( ALLOCATED( EPA_WE_AN_ACET ) ) DEALLOCATE( EPA_WE_AN_ACET )
      IF ( ALLOCATED( EPA_WE_AN_MEK  ) ) DEALLOCATE( EPA_WE_AN_MEK  )
      IF ( ALLOCATED( EPA_WE_AN_PRPE ) ) DEALLOCATE( EPA_WE_AN_PRPE )
      IF ( ALLOCATED( EPA_WE_AN_C2H6 ) ) DEALLOCATE( EPA_WE_AN_C2H6 )
      IF ( ALLOCATED( EPA_WE_AN_C3H8 ) ) DEALLOCATE( EPA_WE_AN_C3H8 )
      IF ( ALLOCATED( EPA_WE_AN_C2H6 ) ) DEALLOCATE( EPA_WE_AN_C2H6 )
      IF ( ALLOCATED( EPA_WE_AN_NH3  ) ) DEALLOCATE( EPA_WE_AN_NH3  )
      IF ( ALLOCATED( EPA_WE_AN_SO2  ) ) DEALLOCATE( EPA_WE_AN_SO2  )
      IF ( ALLOCATED( EPA_WE_AN_SO4  ) ) DEALLOCATE( EPA_WE_AN_SO4  )

      ! Biofuel -- weekday
      IF ( ALLOCATED( EPA_WD_BF_NOX  ) ) DEALLOCATE( EPA_WD_BF_NOX  )
      IF ( ALLOCATED( EPA_WD_BF_CO   ) ) DEALLOCATE( EPA_WD_BF_CO   )
      IF ( ALLOCATED( EPA_WD_BF_ALK4 ) ) DEALLOCATE( EPA_WD_BF_ALK4 )
      IF ( ALLOCATED( EPA_WD_BF_ACET ) ) DEALLOCATE( EPA_WD_BF_ACET )
      IF ( ALLOCATED( EPA_WD_BF_MEK  ) ) DEALLOCATE( EPA_WD_BF_MEK  )
      IF ( ALLOCATED( EPA_WD_BF_PRPE ) ) DEALLOCATE( EPA_WD_BF_PRPE )
      IF ( ALLOCATED( EPA_WD_BF_C2H6 ) ) DEALLOCATE( EPA_WD_BF_C2H6 )
      IF ( ALLOCATED( EPA_WD_BF_C3H8 ) ) DEALLOCATE( EPA_WD_BF_C3H8 )
      IF ( ALLOCATED( EPA_WD_BF_C2H6 ) ) DEALLOCATE( EPA_WD_BF_C2H6 )
      IF ( ALLOCATED( EPA_WD_BF_NH3  ) ) DEALLOCATE( EPA_WD_BF_NH3  )
      IF ( ALLOCATED( EPA_WD_BF_SO2  ) ) DEALLOCATE( EPA_WD_BF_SO2  )
      IF ( ALLOCATED( EPA_WD_BF_SO4  ) ) DEALLOCATE( EPA_WD_BF_SO4  )

      ! Biofuel -- weekend
      IF ( ALLOCATED( EPA_WE_BF_NOX  ) ) DEALLOCATE( EPA_WE_BF_NOX  )
      IF ( ALLOCATED( EPA_WE_BF_CO   ) ) DEALLOCATE( EPA_WE_BF_CO   )
      IF ( ALLOCATED( EPA_WE_BF_ALK4 ) ) DEALLOCATE( EPA_WE_BF_ALK4 )
      IF ( ALLOCATED( EPA_WE_BF_ACET ) ) DEALLOCATE( EPA_WE_BF_ACET )
      IF ( ALLOCATED( EPA_WE_BF_MEK  ) ) DEALLOCATE( EPA_WE_BF_MEK  )
      IF ( ALLOCATED( EPA_WE_BF_PRPE ) ) DEALLOCATE( EPA_WE_BF_PRPE )
      IF ( ALLOCATED( EPA_WE_BF_C2H6 ) ) DEALLOCATE( EPA_WE_BF_C2H6 )
      IF ( ALLOCATED( EPA_WE_BF_C3H8 ) ) DEALLOCATE( EPA_WE_BF_C3H8 )
      IF ( ALLOCATED( EPA_WE_BF_C2H6 ) ) DEALLOCATE( EPA_WE_BF_C2H6 )
      IF ( ALLOCATED( EPA_WE_BF_NH3  ) ) DEALLOCATE( EPA_WE_BF_NH3  )
      IF ( ALLOCATED( EPA_WE_BF_SO2  ) ) DEALLOCATE( EPA_WE_BF_SO2  )
      IF ( ALLOCATED( EPA_WE_BF_SO4  ) ) DEALLOCATE( EPA_WE_BF_SO4  )

      END SUBROUTINE CLEANUP_EPA_NEI
!EOC
      END MODULE EPA_NEI_MOD
