!$Id: tes_nh3_mod.f,v 1.7 2011/02/23 00:08:48 daven Exp $(?)
      MODULE TES_NH3_MOD
  
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
 
      ! Parameters
      INTEGER, PARAMETER           :: LLTES  = 15
      INTEGER, PARAMETER           :: NFOR   = 30
      INTEGER, PARAMETER           :: MAXTES = 3991
      LOGICAL, PARAMETER           :: LTES_NIGHT = .FALSE. 


      ! Record to store data from each TES obs
      TYPE TES_NH3_OBS 
         INTEGER                                 :: NYMD
         INTEGER                                 :: NHMS
         INTEGER                                 :: LLNT
         INTEGER                                 :: QFLAG
         INTEGER                                 :: DFLAG
         TYPE (XPLEX)                                  :: LAT
         TYPE (XPLEX)                                  :: LON
         TYPE (XPLEX),  DIMENSION(LLTES)               :: NH3
         TYPE (XPLEX),  DIMENSION(LLTES,LLTES)         :: AVG_KERNEL
         TYPE (XPLEX),  DIMENSION(LLTES)               :: PRES
         TYPE (XPLEX),  DIMENSION(LLTES,LLTES)         :: OER_INV
         TYPE (XPLEX),  DIMENSION(LLTES)               :: PRIOR
         TYPE (XPLEX),  DIMENSION(LLTES,LLTES)         :: AVK_VMR
         TYPE (XPLEX),  DIMENSION(LLTES,LLTES)         :: OEI_VMR
         CHARACTER(LEN=255)                      :: FILENAME
         TYPE (XPLEX)                                  :: BLVMR
         TYPE (XPLEX),  DIMENSION(LLTES)               :: BLVMR_WGT
      ENDTYPE TES_NH3_OBS  

      TYPE(TES_NH3_OBS)                          :: TES(MAXTES)



      ! Allocatable variables
      TYPE (XPLEX), ALLOCATABLE                        :: NH3_SAVE(:,:)


      

      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_TES_NH3_OBS( NH3,     AVG_KERNEL, PRES, 
     &                             OER_INV, FILENAME,   LLNT,   
     &                             PRIOR,   QFLAG,      DFLAG,
     &                             LAT,     LON,        BLVMR,
     &                             BLVMR_WGT                   )
!
!******************************************************************************
!  Subroutine READ_TES_NH3_OBS reads the file and passes back info contained
!  therein. (dkh, 02/19/09) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME   (TYPE (XPLEX)) : TES observation filename to read
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NH3        (TYPE (XPLEX)) : Retrieval NH3                        [ppmv]
!  (2 ) AVG_KERNEL (TYPE (XPLEX)) : Retrieval averaging kernel           [lnvmr/lnvmr]
!  (3 ) PRES       (TYPE (XPLEX)) : Retrieval pressure levels            [mbar]
!  (4 ) OER_INV    (TYPE (XPLEX)) : TES retrieval error matrix inverse   [((lnvmr))^-2]
!  (5 ) LLNT       (TYPE (XPLEX)) : TES retrieval number of levels                  
!  (6 ) PRIOR      (TYPE (XPLEX)) : TES retrieval prior                  [ppmv]     
!  (7 ) QFLAG     (INTEGER) : TES retrieval quality flag
!  (8 ) DFLAG     (INTEGER) : TES retrieval diagnostic flag
!  (9 ) LAT        (TYPE (XPLEX)) : TES retrieval latitude
!  (10) LON        (TYPE (XPLEX)) : TES retrieval longitude
!  (11) BLVMR      (TYPE (XPLEX)) : TES retrieval BLVMR                  [ppmv]
!  (12) BLVMR_WGT  (TYPE (XPLEX)) : TES retrieval BLVMR mapping          [lnblvr/lnppmv]
!     
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!  (2 ) Add QFLAG and DFLAG (dkh, 11/06/09) 
!  (3 ) Add LAT and LON (dkh, 11/30/09) 
!  (4 ) Add BLVMR and BLVMR_WGT (dkh, 11/02/10) 
!  (5 ) Use averageing kernel in lnvmr (AVG_KERNEL) instead of vmr (AVK_VMR)
!  (6 ) Use observation error cov in lnvmr (OER_INV) instead of lnvrm (OEI_VMR)
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE LOGICAL_ADJ_MOD,        ONLY : LTES_BLVMR 
      USE NETCDF 

#     include "CMN_SIZE" 

      ! Arguments
      TYPE (XPLEX),             INTENT(OUT) :: NH3(LLTES)
      TYPE (XPLEX),             INTENT(OUT) :: AVG_KERNEL(LLTES,LLTES)
      !TYPE (XPLEX)              INTENT(OUT) :: AVK_VMR(LLTES,LLTES)
      TYPE (XPLEX),             INTENT(OUT) :: PRES(LLTES)
      !TYPE (XPLEX),             INTENT(OUT) :: OEI_VMR(LLTES,LLTES)
      TYPE (XPLEX),             INTENT(OUT) :: OER_INV(LLTES,LLTES)
      TYPE (XPLEX),             INTENT(OUT) :: PRIOR(LLTES)
      INTEGER,            INTENT(OUT) :: LLNT
      INTEGER,            INTENT(OUT) :: QFLAG
      INTEGER,            INTENT(OUT) :: DFLAG
      TYPE (XPLEX),             INTENT(OUT) :: LAT
      TYPE (XPLEX),             INTENT(OUT) :: LON
      TYPE (XPLEX),             INTENT(OUT) :: BLVMR
      TYPE (XPLEX),             INTENT(OUT) :: BLVMR_WGT(LLTES)
      CHARACTER(LEN=255), INTENT(IN)  :: FILENAME
    
      ! Local variables
      CHARACTER(LEN=255)              :: READ_FILENAME
      CHARACTER(LEN=5)                :: TMP
      TYPE (XPLEX)                          :: BLVMR_TMP(3)
      TYPE (XPLEX)                          :: BLVMR_WGT_TMP(3,LLTES)

      ! netCDF id's 
      INTEGER                         :: FID, VARID, DIMID, ATTID
      INTEGER                         :: NBL
    
      ! Loop indexes, and error handling.
      INTEGER                         :: L,   LL

      ! dkh debug
      INTEGER :: dim_ids(10), n_dims

      !=================================================================
      ! READ_TES_NH3_OBS begins here!
      !=================================================================

      ! Initialize
      NH3        = 0d0
      !AVK_VMR    = 0d0 
      AVG_KERNEL = 0d0 
      PRES       = 0d0 
      !OEI_VMR    = 0d0 
      OER_INV    = 0d0 
      LLNT       = 0
      PRIOR      = 0d0
      BLVMR_WGT  = 0d0
      BLVMR      = -999d0

      ! Construct complete filename 
      READ_FILENAME = TRIM( DATA_DIR ) // 
     &                TRIM( '../TES_NH3/' )   //
     &                TRIM( 'tes_nh3_gs_July_2009_for_paper/' ) //
     &                TRIM( FILENAME ) 

      WRITE(6,*) '    - READ_TES_NH3_OBS: reading file: ', READ_FILENAME


      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )
   
      ! READ number of retrievals, LLTES 
      CALL CHECK( NF90_INQ_DIMID ( FID, "nretv", DIMID        ), 101 )
      CALL CHECK( NF90_INQUIRE_DIMENSION( FID,DIMID, TMP, LLNT), 102 )

      ! READ NH3 column, NH3
      CALL CHECK( NF90_INQ_VARID( FID, "xretv", VARID       ), 1 )
      CALL CHECK( NF90_GET_VAR  ( FID, VARID,   NH3(1:LLNT) ), 2 )
    
      ! READ averaging kernal, AVG_KERNEL
      CALL CHECK( NF90_INQ_VARID( FID, "avg_kernel", VARID  ), 3 )
      CALL CHECK( NF90_GET_VAR  ( FID, VARID,        
     &                            AVG_KERNEL(1:LLNT,1:LLNT) ), 4 )

      ! READ pressure levels, PRES
      CALL CHECK( NF90_INQ_VARID( FID, "pressure", VARID        ), 5 )
      CALL CHECK( NF90_GET_VAR  ( FID, VARID,      PRES(1:LLNT) ), 6 )

      ! READ inverse observational error, OER_INV
      ! (note: a priori error already subtracted from this term)
      CALL CHECK( NF90_INQ_VARID( FID, "noise_error", VARID   ), 7 )
      CALL CHECK( NF90_GET_VAR  ( FID, VARID,         
     &                            OER_INV(1:LLNT,1:LLNT)      ), 8 )

      ! READ apriori NH3 column, PRIOR
      CALL CHECK( NF90_INQ_VARID( FID, "xa", VARID            ), 9  )
      CALL CHECK( NF90_GET_VAR  ( FID, VARID,   PRIOR(1:LLNT) ), 10 )

      ! READ quality flag QFLAG
      ! 1 = successful
      ! 0 = failed.  For reason why, see DFLAG. 
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL,  "Quality_Flag", 
     &                          QFLAG  ), 12 )

      ! READ diagnostic flag DFLAG
      !  0 = converged, but DOFS < 0.5, though thermal contrast ok 
      ! -1 = converged, but DOFS < 0.5 & thermal contrast poor
      ! -2 = did not converge
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL, "Diagnostic_Flag",  
     &                          DFLAG  ), 14 )

      ! READ latitude LAT 
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL, "Latitude",
     &                          LAT    ), 15 )

      ! READ longitude LON
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL, "Longitude",
     &                          LON    ), 16 )  

      !! READ averaging kernal for vmr retv, AVK_VMR        
      !CALL CHECK( NF90_INQ_VARID( FID, "ak_vmr",  VARID  ), 17 )
      !CALL CHECK( NF90_GET_VAR  ( FID, VARID,        
      !&                            AVK_VMR(1:LLNT,1:LLNT) ), 18 )

!      ! READ inverse observational error, OEI_VMR, based on vmr retv
!      CALL CHECK( NF90_INQ_VARID( FID, "h_vmr",  VARID   ), 19 )
!      CALL CHECK( NF90_GET_VAR  ( FID, VARID,         
!     &                            OEI_VMR(1:LLNT,1:LLNT) ), 20 )

      ! READ BLVMR only for good retrievals, fill othewise
      IF ( QFLAG == 1 .and. DFLAG == 1 .and. LTES_BLVMR ) THEN

         ! dkh debug
         print*, ' quality retv; look for BLVMR '

         ! get NBL
         CALL CHECK( NF90_INQ_DIMID ( FID, "nbl", DIMID      ), 201 )
         CALL CHECK( NF90_INQUIRE_DIMENSION( FID,DIMID, TMP, NBL), 212 )
 
         ! Get BLVMR_TMP
         CALL CHECK( NF90_INQ_VARID( FID, "blvmr",  VARID   ), 21 )
         CALL CHECK( NF90_GET_VAR  ( FID, VARID, BLVMR_TMP(1:NBL) ),
     &      222 )

         ! READ BLVMR_WGT only if BLVMR not = -999
         IF ( BLVMR_TMP(1) > 0 ) THEN  
 
            ! Get BLVMR_WGT_TMP
            CALL CHECK( NF90_INQ_VARID( FID, "blvmr_wgt",  VARID ), 
     &                                                             232 )

! dkh debug
!         CALL CHECK( NF90_INQUIRE_VARIABLE( FID,VARID, ndims = n_dims,
!     &                                      dimids = dim_ids ), 212 )
         
            CALL CHECK( NF90_GET_VAR  ( FID, VARID, 
     &         BLVMR_WGT_TMP(1:NBL,1:LLNT) ), 333 )
           
            ! It is possible to have multiple BLVMRs when DOF > 1.2.
    
            IF ( NBL == 1 ) THEN 

               ! Keep only the first one. 
               BLVMR        = BLVMR_TMP(1) 
               BLVMR_WGT(:) = BLVMR_WGT_TMP(1,:)

            ELSEIF ( NBL > 1 ) THEN 
               ! Keep the last one (dkh, 01/21/11) 
               BLVMR        = BLVMR_TMP(NBL) 
               BLVMR_WGT(:) = BLVMR_WGT_TMP(NBL,:)
 
            ENDIF 
 
         ELSE 

            ! Fill value 
            BLVMR = -999D0 
            
         ENDIF 
      ELSE 

         ! Fill value 
         BLVMR     = -999D0 
         BLVMR_WGT = 0d0

      ENDIF 

      ! Check the data.
      DO L = 1, LLNT
         print*, ' NH3 = ', NH3(L), L
      ENDDO
    
      DO L  = 1, LLNT
         !print*, ' diag(AVG_KERNEL) = ', AVK_VMR(L,L)
         print*, ' diag(AVG_KERNEL) = ', AVG_KERNEL(L,L)
      ENDDO

      IF ( LTES_BLVMR ) print*, ' BLVMR = ', BLVMR

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )
    
      ! Return to calling program
      END SUBROUTINE READ_TES_NH3_OBS
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_TES_NH3_OBS( FILENAME, BLVMR_GC, 
     &                             NH3_GC,   NH3_HAT,  LLNT  )
     
!
!******************************************************************************
!  Subroutine MAKE_TES_NH3_OBS adds BLVMR_GC to the TES data file. (dkh, 11/02/10) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME   (TYPE (XPLEX)) : TES observation filename to read
!  (2 ) BLVMR_GC   (TYPE (XPLEX)) : Boundary layer volume mixing ration in GC [ppmv]
!  (3 ) NH3_GC     (TYPE (XPLEX)) : GEOS-Chem NH3 concentrations on TES grid  [ppmv]
!  (4 ) NH3_HAT    (TYPE (XPLEX)) : GEOS-Chem NH3 after applying TES obs op   [ppmv]
!  (5 ) LLNT      (INTEGER) : Number of TES levels for current retrieval 
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE NETCDF 

#     include "CMN_SIZE" 

      ! Arguments
      CHARACTER(LEN=255), INTENT(IN)  :: FILENAME
      TYPE (XPLEX),             INTENT(IN)  :: BLVMR_GC
      TYPE (XPLEX),             INTENT(IN)  :: NH3_GC(LLTES)
      TYPE (XPLEX),             INTENT(IN)  :: NH3_HAT(LLTES)
      INTEGER,            INTENT(IN)  :: LLNT
    
      ! Local variables
      CHARACTER(LEN=255)              :: READ_FILENAME
      CHARACTER(LEN=5)                :: TMP

      ! netCDF id's 
      INTEGER                         :: FID, VARID, DIMID, ATTID
      INTEGER                         :: VARID2, VARID3, DIMID2
    
      ! Loop indexes, and error handling.
      INTEGER                         :: L,   LL

      !=================================================================
      ! MAKE_TES_NH3_OBS begins here!
      !=================================================================

      ! Construct complete filename 
      READ_FILENAME = TRIM( DATA_DIR ) // 
     &                TRIM( '../TES_NH3/' )   //
     &                TRIM( 'tes_nh3_gs_July_2009_for_paper/' ) //
     &                TRIM( FILENAME ) 

      WRITE(6,*) '    - MAKE_TES_NH3_OBS: reading file: ', READ_FILENAME


      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_WRITE, FID ), 0 )

      ! READ nbl
      CALL CHECK( NF90_INQ_DIMID ( FID, "nbl", DIMID        ), 101 )

      ! READ nretv
      CALL CHECK( NF90_INQ_DIMID ( FID, "nretv", DIMID2     ), 102 )

      ! Place file into define mode
      CALL CHECK( NF90_REDEF( FID ), 1 )

      ! define new BLVMR_GC variable 
      CALL CHECK( NF90_DEF_VAR( FID, "blvmr_gc", NF90_FLOAT, DIMID,
     &     VARID ), 2 )
     
      ! also define new NH3_GC variable 
      CALL CHECK( NF90_DEF_VAR( FID, "nh3_gc", NF90_FLOAT, DIMID2,
     &     VARID2), 3 )
     
      ! also define new NH3_HAT variable 
      CALL CHECK( NF90_DEF_VAR( FID, "nh3_hat", NF90_FLOAT, DIMID2,
     &     VARID3), 4 )
     
      ! end define mode
      CALL CHECK( NF90_ENDDEF( FID ), 5 )

      ! put values in BLVMR_GC
      CALL CHECK( NF90_PUT_VAR( FID, VARID, BLVMR_GC ), 6 )

      ! put values in NH3_GC
      CALL CHECK( NF90_PUT_VAR( FID, VARID2, NH3_GC(1:LLNT) ), 7 )

      ! put values in NH3_HAT
      CALL CHECK( NF90_PUT_VAR( FID, VARID3, NH3_HAT(1:LLNT) ), 8 )

      ! close the file 
      CALL CHECK( NF90_CLOSE( FID ), 9999 )
    
      ! Return to calling program
      END SUBROUTINE MAKE_TES_NH3_OBS
!------------------------------------------------------------------------------

      SUBROUTINE CHECK( STATUS, LOCATION )
!
!******************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries routines
!  (dkh, 02/15/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call    
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was made   
!     
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules 
      USE ERROR_MOD,    ONLY  : ERROR_STOP
      USE NETCDF 
 
      ! Arguments
      INTEGER, INTENT(IN)    :: STATUS 
      INTEGER, INTENT(IN)    :: LOCATION
    
      !=================================================================
      ! CHECK begins here!
      !=================================================================

      IF ( STATUS /= NF90_NOERR ) THEN 
        WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
        WRITE(6,*) 'At location = ', LOCATION 
        CALL ERROR_STOP('netCDF error', 'tes_nh3_mod')
      ENDIF 

      ! Return to calling program
      END SUBROUTINE CHECK

!------------------------------------------------------------------------------

      SUBROUTINE CALC_TES_NH3_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_TES_NH3_FORCE calculates the adjoint forcing from the TES
!  NH3 observations and updates the cost function. (dkh, 02/15/09)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (TYPE (XPLEX)) : Cost funciton                        [unitless]
!     
!     
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!  (1 ) Add more diagnostics.  Now read and write doubled NH3 (dkh, 11/08/09) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AD
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE LOGICAL_ADJ_MOD,    ONLY : LTES_BLVMR
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
      USE TRACERID_MOD,       ONLY : IDTNH3
      USE TRACER_MOD,         ONLY : STT
      USE TRACER_MOD,         ONLY : TCVV

#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER                     :: NTSTART, NTSTOP, NT 
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL,     LLNT 
      TYPE (XPLEX)                      :: GC_PRES(LLPAR)
      TYPE (XPLEX)                      :: GC_NH3_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: GC_NH3(LLTES)
      TYPE (XPLEX)                      :: GC_PSURF
      TYPE (XPLEX)                      :: MAP(LLPAR,LLTES)
      TYPE (XPLEX)                      :: NH3_HAT(LLTES)
      TYPE (XPLEX)                      :: NH3_PERT(LLTES)
      TYPE (XPLEX)                      :: FORCE(LLTES)
      TYPE (XPLEX)                      :: DIFF(LLTES)
      TYPE (XPLEX)                      :: NEW_COST(MAXTES)
      TYPE (XPLEX)                      :: OLD_COST

      TYPE (XPLEX)                      :: GC_NH3_NATIVE_DBL(LLPAR)
      TYPE (XPLEX)                      :: GC_NH3_DBL(LLTES)
      TYPE (XPLEX)                      :: NH3_HAT_DBL(LLTES)
      TYPE (XPLEX)                      :: NH3_PERT_DBL(LLTES)

      TYPE (XPLEX)                      :: ADJ_GC_NH3_NATIVE(LLPAR)
      TYPE (XPLEX)                      :: ADJ_NH3_HAT(LLTES)
      TYPE (XPLEX)                      :: ADJ_NH3_PERT(LLTES)
      TYPE (XPLEX)                      :: ADJ_GC_NH3(LLTES)
      TYPE (XPLEX)                      :: ADJ_DIFF(LLTES)
      TYPE (XPLEX)                      :: BLVMR_GC
      TYPE (XPLEX)                      :: BLVMR_TC
      TYPE (XPLEX)                      :: TMP1
   
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME



      !=================================================================
      ! CALC_TES_NH3_FORCE begins here!
      !=================================================================

      print*, '     - CALC_TES_NH3_FORCE '

      NEW_COST = 0d0 

      !IF ( FIRST ) THEN 
      !   CALL READ_NH3_FILE( )
      !   FIRST = .FALSE. 
      !ENDIF 

      ! Open files for output
      IF ( FIRST ) THEN
         FILENAME = 'pres.NN.m' 
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'gc_nh3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 102,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'tes_nh3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 103,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'apriori.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 104,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'diff.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 105,      FILE=TRIM( FILENAME    ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'force.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 106,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'nt_ll.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 107,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'adj_nh3_pert.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 108,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'adj_gc_nh3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 109,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'exp_nh3_hat.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 110,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FILENAME = 'exp_nh3_hat_dbl.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC ) 
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 111,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         FIRST = .FALSE.

      ENDIF


      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Get range of TES observations the occur during the current hour
      CALL GET_NT_RANGE( GET_NYMD(), GET_NHMS(), NTSTART, NTSTOP ) 

      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 

         print*, ' No matching TES NH3 obs for this hour'
        RETURN
      ENDIF 

      ! Loop over TES observations (first do a sequential loop for 
      ! reading the data)
      DO NT  = NTSTART, NTSTOP, -1

         print*, '     - CALC_TES_NH3_FORCE: reading record ', NT 

         CALL READ_TES_NH3_OBS( TES(NT)%NH3(:), TES(NT)%AVG_KERNEL(:,:), 
     &                          TES(NT)%PRES(:),  TES(NT)%OER_INV(:,:), 
     &                          TES(NT)%FILENAME, TES(NT)%LLNT,
     &                          TES(NT)%PRIOR(:), TES(NT)%QFLAG, 
     &                          TES(NT)%DFLAG,                         
     &                          TES(NT)%LAT,      TES(NT)%LON,
     &                          TES(NT)%BLVMR,    TES(NT)%BLVMR_WGT(:) )

     
         print*, ' retrieved data for lat / lon : ', TES(NT)%LAT, 
     &      TES(NT)%LON


      ENDDO

! need to update this in order to do i/o with this loop parallel 
!!      ! Now do a parallel loop for analyzing data 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, MAP, LLNT, IIJJ,  I, J,  L,   LL    )
!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_NH3,  DIFF  )
!!$OMP+PRIVATE( GC_NH3_NATIVE, NH3_PERT, NH3_HAT, FORCE )
!!$OMP+PRIVATE( ADJ_GC_NH3_NATIVE,       ADJ_GC_NH3     )
!!$OMP+PRIVATE( ADJ_NH3_PERT,            ADJ_NH3_HAT    )
!!$OMP+PRIVATE( ADJ_DIFF                                )
      DO NT  = NTSTART, NTSTOP, -1

         print*, '     - CALC_TES_NH3_FORCE: analyzing record ', NT 

         ! Skip spiky retrievals that look cloud contaminated 
         !IF ( NT == 415 .or. NT == 668 ) THEN
         !   print*, ' SKIPPING record owing to cloud ', NT
         !   CYCLE
         !ENDIF
         IF ( .not. LTES_NIGHT .and. TES(NT)%NHMS < 150000 ) THEN
            print*, ' SKIPPING night time retrievals ', NT
            CYCLE
         ENDIF


         ! Check quality of retrieval 
         IF ( TES(NT)%QFLAG .ne. 1 ) THEN

            IF ( TES(NT)%DFLAG .ne. 0 ) THEN
               print*, ' SKIPPING record ', NT
               print*, ' QFLAG, DFLAG = ', TES(NT)%QFLAG, TES(NT)%DFLAG
               CYCLE
            ELSE

               print*, ' QFLAG = 0 but DFLAG = 0 '

            ENDIF

         ENDIF

         ! For safety, initialize these up to LLTES 
         GC_NH3(:)       = 0d0 
         MAP(:,:)        = 0d0 
         ADJ_NH3_HAT(:)  = 0d0 
         FORCE(:)        = 0d0 


         ! Copy LLNT to make coding a bit cleaner
         LLNT = TES(NT)%LLNT

         ! Get grid box of current record
         IIJJ = GET_IJ( DCMPLX(TES(NT)%LON),DCMPLX(TES(NT)%LAT) )
         I    = IIJJ(1)
         J    = IIJJ(2)

         ! Get GC pressure levels (mbar) 
         DO L = 1, LLPAR
            GC_PRES(L) = GET_PCENTER(I,J,L)
         ENDDO

         ! Get GC surface pressure (mbar) 
         GC_PSURF = GET_PEDGE(I,J,1) 

         
         ! Calculate the interpolation weight matrix 
         MAP(1:LLPAR,1:LLNT) 
     &      = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
     &                    LLNT,  TES(NT)%PRES(1:LLNT), GC_PSURF  )


         ! Get NH3 values at native model resolution
         GC_NH3_NATIVE(:)  = CHK_STT(I,J,:,IDTNH3) 
        
         print*, 'I,J = ', I, J

         ! Convert from kg/box to ppm
         GC_NH3_NATIVE(:) = GC_NH3_NATIVE(:) * TCVV(IDTNH3)
     &                    / AD(I,J,:) * 1d6

         NH3_SAVE(:,NT) = GC_NH3_NATIVE(:)

! skip for real data  
!         ! Get NH3 values from doubled emissions run [ppmv]
!         GC_NH3_NATIVE_DBL(:) 
!     &      = GET_DOUBLED_NH3( GET_NYMD(),          GET_NHMS(), 
!     &                         REAL(TES(NT)%LON,4), REAL(TES(NT)%LAT,4))

         ! Interpolate GC NH3 column to TES grid 
         DO LL = 1, LLNT
            GC_NH3(LL) = 0d0 
            DO L = 1, LLPAR 
               GC_NH3(LL) = GC_NH3(LL) 
     &                    + MAP(L,LL) * GC_NH3_NATIVE(L) 
            ENDDO
         ENDDO

         !print*, ' gc_nh3 =', gc_nh3(:)

! skip for real data
!         ! Interpolate doubled GC NH3 column to TES grid 
!         DO LL = 1, LLNT
!            GC_NH3_DBL(LL) = 0d0 
!            DO L = 1, LLPAR 
!               GC_NH3_DBL(LL) = GC_NH3_DBL(LL) 
!     &                    + MAP(L,LL) * GC_NH3_NATIVE_DBL(L) 
!            ENDDO
!         ENDDO

         ! dkh debug: compare profiles:
         print*, ' GC_PRES, GC_native_NH3 '
         WRITE(6,100) (GC_PRES(L), GC_NH3_NATIVE(L), 
     &                   L = LLPAR, 1, -1 )
         print*, ' TES_PRES, GC_NH3  ' 
         WRITE(6,100) (TES(NT)%PRES(LL), GC_NH3(LL), LL = LLNT, 1, -1 ) 
 100  FORMAT(1X,F16.8,1X,F16.8)


         !--------------------------------------------------------------
         ! Apply TES observation operator
         !
         !   x_hat = ln(x_a) + A_k ( ln(x_m) - ln(x_a) ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by TES [lnvmr]
         !    x_a   = TES apriori column               [ppmv]
         !    x_m   = GC modeled column                [ppmv]
         !    A_k   = TES averaging kernel             [lnvmr/lnvmr]
         !--------------------------------------------------------------
 
         ! fill
         DO L = 1, LLNT 
           GC_NH3(L) = MAX(GC_NH3(L),0.00000001)
         ENDDO
   
         ! x_m - x_a
         DO L = 1, LLNT 
           NH3_PERT(L) = LOG(GC_NH3(L)) - LOG(TES(NT)%PRIOR(L))
         ENDDO
     
         ! x_a + A_k * ( x_m - x_a )  
         DO L = 1, LLNT
            NH3_HAT(L)    = 0d0 
            DO LL = 1, LLNT
               NH3_HAT(L) = NH3_HAT(L) 
     &                    + TES(NT)%AVG_KERNEL(LL,L) * NH3_PERT(LL) 
!     &                    + TES(NT)%AVK_VMR(L,LL) * NH3_PERT(LL) 
            ENDDO
            NH3_HAT(L)    = NH3_HAT(L) + LOG(TES(NT)%PRIOR(L))
         ENDDO

! skip for real data
!         ! x_m - x_a for doubled 
!         DO L = 1, LLNT 
!           GC_NH3_DBL(L)   = MAX(GC_NH3_DBL(L), 1d-10)
!           NH3_PERT_DBL(L) = LOG(GC_NH3_DBL(L)) - LOG(TES(NT)%PRIOR(L))
!         ENDDO
!     
!         ! x_a + A_k * ( x_m - x_a )  
!         DO L = 1, LLNT
!            NH3_HAT_DBL(L)    = 0d0 
!            DO LL = 1, LLNT
!               NH3_HAT_DBL(L) = NH3_HAT_DBL(L) 
!!     &                        + TES(NT)%AVG_KERNEL(L,LL) 
!     &                        + TES(NT)%AVG_KERNEL(LL,L) 
!     &                        * NH3_PERT_DBL(LL) 
!            ENDDO
!            NH3_HAT_DBL(L)    = NH3_HAT_DBL(L) + LOG(TES(NT)%PRIOR(L))
!         ENDDO

         !--------------------------------------------------------------
         ! Calculate cost function, given S is error on ln(vmr)
         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ ln(model - obs ]
         !--------------------------------------------------------------

         ! Calculate difference between modeled and observed profile
         DO L = 1, LLNT 
            DIFF(L) = NH3_HAT(L) -  LOG( TES(NT)%NH3(L) )
         ENDDO 
          
         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF 
         DO L = 1, LLNT
            FORCE(L)     = 0d0 
            ! try using just the diagonal of OEI_VMR
            !DO LL = 1, LLNT
            !   FORCE(L)  = FORCE(L) + TES(NT)%OEI_VMR(L,LL) * DIFF(LL)
            !ENDDO
            !FORCE(L)  = TES(NT)%OEI_VMR(L,L) * DIFF(L)
            ! try using diag and a minimum error of 0.01 ppb --> 1/0.00001^2 --> 1d10
            !IF ( TES(NT)%OEI_VMR(L,L)  < 1d10 ) THEN 
            !   FORCE(L)  = TES(NT)%OEI_VMR(L,L) * DIFF(L)
            !ELSE 
            !   FORCE(L)  = 1d10 * DIFF(L)
            !ENDIF 

            ! put a cap on the error, but use the entire matrix
            DO LL = 1, LLNT

               ! Cap values

               IF ( MAXVAL(TES(NT)%NH3(:)) < 0.001 ) THEN
                 IF ( L == LL ) THEN
                  TES(NT)%OER_INV(L,LL)
     &               = MIN(0.5d0 , ABS( TES(NT)%OER_INV(L,LL) ) )
                 ELSE
                  TMP1
     &               = MIN(0.25d0, ABS( TES(NT)%OER_INV(L,LL) ) )

                  ! retain sign 
                  TES(NT)%OER_INV(L,LL)
     &               = SIGN( TMP1, TES(NT)%OER_INV(L,LL)  )
                 ENDIF
               ENDIF

               FORCE(L)  = FORCE(L) + TES(NT)%OER_INV(L,LL) * DIFF(LL)

            ENDDO

            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
         ENDDO

         ! dkh debug: compare profiles:
         print*, ' TES_PRIOR, NH3_HAT, NH3_TES in vmr'
         WRITE(6,101) (TES(NT)%PRIOR(L), EXP(NH3_HAT(L)),TES(NT)%NH3(L),
     &             L,    L = LLNT, 1, -1 )
         print*, ' TES_PRIOR, NH3_HAT, NH3_TES in lnvmr'
         WRITE(6,101) (LOG(TES(NT)%PRIOR(L)), NH3_HAT(L), 
     &     LOG(TES(NT)%NH3(L)), L,  L = LLNT, 1, -1 )
 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1x,i3)

         !--------------------------------------------------------------
         ! Calculate BLVMR_GC and save to the data file 
         !--------------------------------------------------------------
         
         ! Only do this for records which have a BLVMR
         IF( LTES_BLVMR .and.  TES(NT)%BLVMR > 0 ) THEN 

            BLVMR_GC = 0d0 
            BLVMR_TC = 0d0 

            DO L = 1, LLNT
 
               BLVMR_GC = BLVMR_GC 
     &                  + TES(NT)%BLVMR_WGT(L)
     &                  * NH3_HAT(L)

               BLVMR_TC = BLVMR_TC 
     &                  + TES(NT)%BLVMR_WGT(L)
     &                  * TES(NT)%NH3(L) 
            ENDDO 
 
            BLVMR_GC = EXP( BLVMR_GC )
            BLVMR_TC = EXP( BLVMR_TC )
             
            ! Append to data file
            CALL MAKE_TES_NH3_OBS( TES(NT)%FILENAME, BLVMR_GC, 
     &                             GC_NH3,           EXP(NH3_HAT(:)), 
     &                             LLNT                               )

            print*, ' BLVMRs = ', TES(NT)%BLVMR, BLVMR_GC, BLVMR_TC, NT

         ENDIF 


         !--------------------------------------------------------------
         ! Begin adjoint calculations 
         !--------------------------------------------------------------

         ! dkh debug
         print*, 'DIFF , FORCE, diag OEI ' 
         WRITE(6,102) (DIFF(L), FORCE(L), TES(NT)%OER_INV(L,L),
     &       L = LLNT, 1, -1 )
 102  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)

         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         ADJ_DIFF(:) = FORCE(:) 

         ! Adjoint of difference
         DO L = 1, LLNT 
            ADJ_NH3_HAT(L) =  ADJ_DIFF(L)
         ENDDO 

         ! adjoint of TES operator
         DO L  = 1, LLNT
            ADJ_NH3_PERT(L)    = 0d0
            DO LL = 1, LLNT
               ADJ_NH3_PERT(L) = ADJ_NH3_PERT(L) 
     &                         + TES(NT)%AVG_KERNEL(L,LL)
     &                         * ADJ_NH3_HAT(LL)

           ENDDO
         ENDDO
 
         ! Adjoint of x_m - x_a
         DO L = 1, LLNT 
           ! fwd code:
           !NH3_PERT(L) = LOG(GC_NH3(L)) - LOG(TES(NT)%PRIOR(L)))
           ! adj code:
           ADJ_GC_NH3(L) = ADJ_NH3_PERT(L) / GC_NH3(L)
         ENDDO

         ! dkh debug
         print*, 'ADJ_NH3_HAT, ADJ_NH3_PERT, ADJ_GC_NH3 '
         WRITE(6,103) (ADJ_NH3_HAT(L), ADJ_NH3_PERT(L), ADJ_GC_NH3(L), 
     &       L = LLNT, 1, -1 )
 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)

         ! adjoint of interpolation 
         DO L  = 1, LLPAR
            ADJ_GC_NH3_NATIVE(L) = 0d0 
            DO LL = 1, LLNT
               ADJ_GC_NH3_NATIVE(L) = ADJ_GC_NH3_NATIVE(L)
     &                              + MAP(L,LL) * ADJ_GC_NH3(LL)
            ENDDO
         ENDDO
         
         ! Adjoint of unit conversion 
         ADJ_GC_NH3_NATIVE(:) = ADJ_GC_NH3_NATIVE(:) * TCVV(IDTNH3)
     &                        / AD(I,J,:) * 1d6

         ! Pass adjoint back to adjoint tracer array
         STT_ADJ(I,J,:,IDTNH3)  = STT_ADJ(I,J,:,IDTNH3)
     &                          + ADJ_GC_NH3_NATIVE(:)

         ! dkh debug
         print*, ' adj_stt force = ', ADJ_GC_NH3_NATIVE(:)

         ! dkh debug
         print*, 'ADJ_GC_NH3_NATIVE, ADJ_GC_NH3_NATIVE conv, ADJ_STT ' 
         WRITE(6,104) (ADJ_GC_NH3_NATIVE(L) * AD(I,J,L) 
     &                                      / 1d6 / TCVV(IDTNH3), 
     &       ADJ_GC_NH3_NATIVE(L), STT_ADJ(I,J,L,IDTNH3),
     &       L = LLPAR, 1, -1 )
 104  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)


         WRITE(101,112) ( TES(NT)%PRES(LL),  LL=LLNT,1,-1)
         WRITE(102,112) ( GC_NH3(LL),        LL=LLNT,1,-1)
         WRITE(103,112) ( TES(NT)%NH3(LL),   LL=LLNT,1,-1)
         WRITE(104,112) ( TES(NT)%PRIOR(LL), LL=LLNT,1,-1)
         WRITE(105,112) ( DIFF(LL),          LL=LLNT,1,-1)
         WRITE(106,112) ( FORCE(LL),         LL=LLNT,1,-1)
         WRITE(107,111) NT, LLNT
         WRITE(108,112) ( ADJ_NH3_PERT(LL),  LL=LLNT,1,-1)
         WRITE(109,112) ( ADJ_GC_NH3(LL),    LL=LLNT,1,-1)
         WRITE(110,112) ( EXP(NH3_HAT(LL)),  LL=LLNT,1,-1)
         ! skip for real data
         !WRITE(111,110) ( EXP(NH3_HAT_DBL(LL)),  LL=LLNT,1,-1)
 110     FORMAT(F16.8,1X)
 111     FORMAT(i4,1X,i4,1x)
 112     FORMAT(D14.6,1X)



      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function 
      COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))

      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, ' TES contribution           = ', COST_FUNC - OLD_COST  



      ! dkh debug
      !IF ( NTSTOP == 1 ) THEN 
      !   CALL MAKE_NH3_FILE( )
      !ENDIF 

      ! Return to calling program
      END SUBROUTINE CALC_TES_NH3_FORCE

!!------------------------------------------------------------------------------
!
!      SUBROUTINE CALC_TES_NH3_FORCE_FD( COST_FUNC, PERT, ADJ )
!!
!!******************************************************************************
!!  Subroutine to test CALC_TES_NH3_FORCE (dkh, 09/30/10).
!! 
!!  It has a few changes from the original, noted with 'for FD test'
!!
!!  Call from CALC_FORCE_FOR_OBS:
!!      USE TES_NH3_MOD,          ONLY : CALC_TES_NH3_FORCE_FD
!!...
!!      ! for FD test
!!      TYPE (XPLEX)              :: COST_FUNC_0
!!      TYPE (XPLEX)              :: COST_FUNC_1
!!      TYPE (XPLEX)              :: COST_FUNC_2
!!      TYPE (XPLEX)              :: PERT(LLPAR)
!!      TYPE (XPLEX)              :: ADJ(LLPAR)
!!      TYPE (XPLEX)              :: FD(LLPAR)
!!      TYPE (XPLEX)              :: ADJ_SAVE(LLPAR)
!!...
!!      PERT(:) = 1D0
!!      CALL CALC_TES_NH3_FORCE_FD( COST_FUNC_0, PERT, ADJ )
!!      ADJ_SAVE(:) = ADJ(:)
!!      print*, 'do3:  COST_FUNC_0 = ', COST_FUNC_0
!!      DO L = 1, LLPAR
!!         PERT(:) = 1D0
!!         PERT(L) = 1.1
!!         COST_FUNC = 0D0
!!         CALL CALC_TES_NH3_FORCE_FD( COST_FUNC_1, PERT, ADJ )
!!         PERT(L) = 0.9
!!         COST_FUNC = 0D0
!!         CALL CALC_TES_NH3_FORCE_FD( COST_FUNC_2, PERT, ADJ )
!!         FD(L)       = ( COST_FUNC_1 - COST_FUNC_2 ) / 0.2d0
!!         print*, 'do3:  FD  = ', FD(L), L
!!         print*, 'do3:  ADJ = ', ADJ_SAVE(L), L
!!         print*, 'do3:  COST = ', COST_FUNC_2, COST_FUNC_1, L
!!         print*, 'do3:  FD / ADJ ', FD(L) / ADJ_SAVE(L) , L
!!      ENDDO
!!
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
!      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
!      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
!      USE CHECKPT_MOD,        ONLY : CHK_STT
!      USE DAO_MOD,            ONLY : AD
!      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
!      USE GRID_MOD,           ONLY : GET_IJ
!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
!      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
!      USE TRACERID_MOD,       ONLY : IDTNH3
!      USE TRACER_MOD,         ONLY : STT
!      USE TRACER_MOD,         ONLY : TCVV
!
!#     include      "CMN_SIZE"      ! Size params
!
!      ! Arguments
!      TYPE (XPLEX), INTENT(INOUT)       :: COST_FUNC
!   
!      ! for FD test
!      TYPE (XPLEX), INTENT(IN)          :: PERT(LLPAR) 
!      TYPE (XPLEX), INTENT(OUT)         :: ADJ(LLPAR)
!
!      ! Local variables 
!      INTEGER                     :: NTSTART, NTSTOP, NT 
!      INTEGER                     :: IIJJ(2), I,      J
!      INTEGER                     :: L,       LL,     LLNT 
!      TYPE (XPLEX)                      :: GC_PRES(LLPAR)
!      TYPE (XPLEX)                      :: GC_NH3_NATIVE(LLPAR)
!      TYPE (XPLEX)                      :: GC_NH3(LLTES)
!      TYPE (XPLEX)                      :: GC_PSURF
!      TYPE (XPLEX)                      :: MAP(LLPAR,LLTES)
!      TYPE (XPLEX)                      :: NH3_HAT(LLTES)
!      TYPE (XPLEX)                      :: NH3_PERT(LLTES)
!      TYPE (XPLEX)                      :: FORCE(LLTES)
!      TYPE (XPLEX)                      :: DIFF(LLTES)
!      TYPE (XPLEX)                      :: NEW_COST(MAXTES)
!      TYPE (XPLEX)                      :: OLD_COST
!
!      TYPE (XPLEX)                      :: GC_NH3_NATIVE_DBL(LLPAR)
!      TYPE (XPLEX)                      :: GC_NH3_DBL(LLTES)
!      TYPE (XPLEX)                      :: NH3_HAT_DBL(LLTES)
!      TYPE (XPLEX)                      :: NH3_PERT_DBL(LLTES)
!
!      TYPE (XPLEX)                      :: ADJ_GC_NH3_NATIVE(LLPAR)
!      TYPE (XPLEX)                      :: ADJ_NH3_HAT(LLTES)
!      TYPE (XPLEX)                      :: ADJ_NH3_PERT(LLTES)
!      TYPE (XPLEX)                      :: ADJ_GC_NH3(LLTES)
!      TYPE (XPLEX)                      :: ADJ_DIFF(LLTES)
!   
!      LOGICAL, SAVE               :: FIRST = .TRUE. 
!      INTEGER                     :: IOS
!      CHARACTER(LEN=255)          :: FILENAME
!
!
!
!      !=================================================================
!      ! CALC_TES_NH3_FORCE begins here!
!      !=================================================================
!
!      print*, '     - CALC_TES_NH3_FORCE '
!
!      NEW_COST = 0d0 
!
!      !IF ( FIRST ) THEN 
!      !   CALL READ_NH3_FILE( )
!      !   FIRST = .FALSE. 
!      !ENDIF 
!
!      ! Open files for output
!      IF ( FIRST ) THEN
!         FILENAME = 'pres.NN.m' 
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 101,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'gc_nh3.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 102,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'tes_nh3.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 103,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'apriori.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 104,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'diff.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 105,      FILE=TRIM( FILENAME    ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'force.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 106,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'nt_ll.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 107,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'adj_nh3_pert.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 108,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'adj_gc_nh3.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 109,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'exp_nh3_hat.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 110,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         FILENAME = 'exp_nh3_hat_dbl.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC ) 
!         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 111,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         !for FD test: move this to later 
!!         FIRST = .FALSE.
!
!
!      ENDIF
!
!
!      ! Save a value of the cost function first
!      OLD_COST = COST_FUNC
!
!      ! Get range of TES observations the occur during the current hour
!      ! for FD test 
!!      CALL GET_NT_RANGE( GET_NYMD(), GET_NHMS(), NTSTART, NTSTOP ) 
!!     
!!
!!      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN 
!!
!!         print*, ' No matching TES NH3 obs for this hour'
!!        RETURN
!!      ENDIF 
!!
!!      ! Loop over TES observations (first do a sequential loop for 
!!      ! reading the data)
!!      DO NT  = NTSTART, NTSTOP, -1
!       NTSTART = 1
!       NTSTOP  = 1
!       NT      = 1
!       IF ( FIRST ) THEN 
!       
!         print*, '     - CALC_TES_NH3_FORCE: reading record ', NT 
!
!         CALL READ_TES_NH3_OBS( TES(NT)%NH3(:), TES(NT)%AVK_VMR(:,:), 
!     &                          TES(NT)%PRES(:),  TES(NT)%OEI_VMR(:,:),     
!     &                          TES(NT)%FILENAME, TES(NT)%LLNT,
!     &                          TES(NT)%PRIOR(:), TES(NT)%QFLAG, 
!     &                          TES(NT)%DFLAG,                         
!     &                          TES(NT)%LAT,      TES(NT)%LON          )
!
!     
!         print*, ' retrieved data for lat / lon : ', TES(NT)%LAT, 
!     &      TES(NT)%LON
!
!
!       ! for FD test 
!!      ENDDO
!         FIRST = .FALSE.
!       ENDIF 
!
!! need to update this in order to do i/o with this loop parallel 
!!!      ! Now do a parallel loop for analyzing data 
!!!$OMP PARALLEL DO
!!!$OMP+DEFAULT( SHARED )
!!!$OMP+PRIVATE( NT, MAP, LLNT, IIJJ,  I, J,  L,   LL    )
!!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_NH3,  DIFF  )
!!!$OMP+PRIVATE( GC_NH3_NATIVE, NH3_PERT, NH3_HAT, FORCE )
!!!$OMP+PRIVATE( ADJ_GC_NH3_NATIVE,       ADJ_GC_NH3     )
!!!$OMP+PRIVATE( ADJ_NH3_PERT,            ADJ_NH3_HAT    )
!!!$OMP+PRIVATE( ADJ_DIFF                                )
!      DO NT  = NTSTART, NTSTOP, -1
!
!         print*, '     - CALC_TES_NH3_FORCE: analyzing record ', NT 
!
!         ! Skip spiky retrievals that look cloud contaminated 
!         IF ( NT == 415 .or. NT == 668 ) THEN
!            print*, ' SKIPPING record owing to cloud ', NT
!            CYCLE
!         ENDIF
!
!         ! Check quality of retrieval 
!         IF ( TES(NT)%QFLAG .ne. 1 ) THEN
!
!            IF ( TES(NT)%DFLAG .ne. 0 ) THEN
!               print*, ' SKIPPING record ', NT
!               print*, ' QFLAG, DFLAG = ', TES(NT)%QFLAG, TES(NT)%DFLAG
!               CYCLE
!            ELSE
!
!               print*, ' QFLAG = 0 but DFLAG = 0 '
!
!            ENDIF
!
!         ENDIF
!
!         ! For safety, initialize these up to LLTES 
!         GC_NH3(:)       = 0d0 
!         MAP(:,:)        = 0d0 
!         ADJ_NH3_HAT(:)  = 0d0 
!         FORCE(:)        = 0d0 
!
!
!         ! Copy LLNT to make coding a bit cleaner
!         LLNT = TES(NT)%LLNT
!
!         ! Get grid box of current record
!         IIJJ = GET_IJ( REAL(TES(NT)%LON,4), REAL(TES(NT)%LAT,4) )
!         I    = IIJJ(1)
!         J    = IIJJ(2)
!
!         ! Get GC pressure levels (mbar) 
!         DO L = 1, LLPAR
!            GC_PRES(L) = GET_PCENTER(I,J,L)
!         ENDDO
!
!         ! Get GC surface pressure (mbar) 
!         GC_PSURF = GET_PEDGE(I,J,1) 
!
!         
!         ! Calculate the interpolation weight matrix 
!         MAP(1:LLPAR,1:LLNT) 
!     &      = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF, 
!     &                    LLNT,  TES(NT)%PRES(1:LLNT), GC_PSURF  )
!
!
!         ! Get NH3 values at native model resolution
!         ! for FD test 
!         !GC_NH3_NATIVE(:)  = CHK_STT(I,J,:,IDTNH3) 
!         GC_NH3_NATIVE(:)  = CHK_STT(I,J,:,IDTNH3) * PERT(:)
!        
!         print*, 'I,J = ', I, J
!
!         ! Convert from kg/box to ppm
!         GC_NH3_NATIVE(:) = GC_NH3_NATIVE(:) * TCVV(IDTNH3)
!     &                    / AD(I,J,:) * 1d6
!
!         NH3_SAVE(:,NT) = GC_NH3_NATIVE(:)
!
!! skip for real data  
!!         ! Get NH3 values from doubled emissions run [ppmv]
!!         GC_NH3_NATIVE_DBL(:) 
!!     &      = GET_DOUBLED_NH3( GET_NYMD(),          GET_NHMS(), 
!!     &                         REAL(TES(NT)%LON,4), REAL(TES(NT)%LAT,4))
!
!         ! Interpolate GC NH3 column to TES grid 
!         DO LL = 1, LLNT
!            GC_NH3(LL) = 0d0 
!            DO L = 1, LLPAR 
!               GC_NH3(LL) = GC_NH3(LL) 
!     &                    + MAP(L,LL) * GC_NH3_NATIVE(L) 
!            ENDDO
!         ENDDO
!
!         print*, ' gc_nh3 =', gc_nh3(:)
!
!! skip for real data
!!         ! Interpolate doubled GC NH3 column to TES grid 
!!         DO LL = 1, LLNT
!!            GC_NH3_DBL(LL) = 0d0 
!!            DO L = 1, LLPAR 
!!               GC_NH3_DBL(LL) = GC_NH3_DBL(LL) 
!!     &                    + MAP(L,LL) * GC_NH3_NATIVE_DBL(L) 
!!            ENDDO
!!         ENDDO
!
!         ! dkh debug: compare profiles:
!         print*, ' GC_PRES, GC_native_NH3 '
!         WRITE(6,100) (GC_PRES(L), GC_NH3_NATIVE(L), 
!     &                   L = LLPAR, 1, -1 )
!         print*, ' TES_PRES, GC_NH3  ' 
!         WRITE(6,100) (TES(NT)%PRES(LL), GC_NH3(LL), LL = LLNT, 1, -1 ) 
! 100  FORMAT(1X,F16.8,1X,F16.8)
!
!
!         !--------------------------------------------------------------
!         ! Apply TES observation operator
!         !
!         !   x_hat = x_a + A_k ( x_m - x_a ) 
!         !  
!         !  where  
!         !    x_hat = GC modeled column as seen by TES [ppmv]
!         !    x_a   = TES apriori column               [ppmv]
!         !    x_m   = GC modeled column                [ppmv]
!         !    A_k   = TES averaging kernel 
!         !--------------------------------------------------------------
! 
!         ! x_m - x_a
!         DO L = 1, LLNT 
!           NH3_PERT(L) = GC_NH3(L) - TES(NT)%PRIOR(L)
!         ENDDO
!     
!         ! x_a + A_k * ( x_m - x_a )  
!         DO L = 1, LLNT
!            NH3_HAT(L)    = 0d0 
!            DO LL = 1, LLNT
!               NH3_HAT(L) = NH3_HAT(L) 
!     &                    + TES(NT)%AVK_VMR(LL,L) * NH3_PERT(LL) 
!!     &                    + TES(NT)%AVK_VMR(L,LL) * NH3_PERT(LL) 
!            ENDDO
!            NH3_HAT(L)    = NH3_HAT(L) + TES(NT)%PRIOR(L)
!         ENDDO
!
!! skip for real data
!!         ! x_m - x_a for doubled 
!!         DO L = 1, LLNT 
!!           GC_NH3_DBL(L)   = MAX(GC_NH3_DBL(L), 1d-10)
!!           NH3_PERT_DBL(L) = LOG(GC_NH3_DBL(L)) - LOG(TES(NT)%PRIOR(L))
!!         ENDDO
!!     
!!         ! x_a + A_k * ( x_m - x_a )  
!!         DO L = 1, LLNT
!!            NH3_HAT_DBL(L)    = 0d0 
!!            DO LL = 1, LLNT
!!               NH3_HAT_DBL(L) = NH3_HAT_DBL(L) 
!!!     &                        + TES(NT)%AVG_KERNEL(L,LL) 
!!     &                        + TES(NT)%AVG_KERNEL(LL,L) 
!!     &                        * NH3_PERT_DBL(LL) 
!!            ENDDO
!!            NH3_HAT_DBL(L)    = NH3_HAT_DBL(L) + LOG(TES(NT)%PRIOR(L))
!!         ENDDO
!
!         !--------------------------------------------------------------
!         ! Calculate cost function, given S is error on ln(vmr)
!         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ ln(model - obs ]
!         !--------------------------------------------------------------
!
!         ! Calculate difference between modeled and observed profile
!         DO L = 1, LLNT 
!            DIFF(L) = NH3_HAT(L) -  TES(NT)%NH3(L) 
!         ENDDO 
!          
!         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF 
!         DO L = 1, LLNT
!            FORCE(L)     = 0d0 
!            ! try using just the diagonal of OEI_VMR
!            !DO LL = 1, LLNT
!            !   FORCE(L)  = FORCE(L) + TES(NT)%OEI_VMR(L,LL) * DIFF(LL)
!            !ENDDO
!            !FORCE(L)  = TES(NT)%OEI_VMR(L,L) * DIFF(L)
!            ! try using diag and a minimum error of 0.01 ppb --> 1/0.00001^2 --> 1d10
!            IF ( TES(NT)%OEI_VMR(L,L)  < 1d10 ) THEN 
!               FORCE(L)  = TES(NT)%OEI_VMR(L,L) * DIFF(L)
!            ELSE 
!               FORCE(L)  = 1d10 * DIFF(L)
!            ENDIF 
!            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
!         ENDDO
!
!         ! dkh debug: compare profiles:
!         print*, ' TES_PRIOR, NH3_HAT, NH3_TES in vmr'
!         WRITE(6,101) (TES(NT)%PRIOR(L), NH3_HAT(L),TES(NT)%NH3(L),
!     &             L,    L = LLNT, 1, -1 )
!         print*, ' TES_PRIOR, NH3_HAT, NH3_TES in lnvmr'
!!         WRITE(6,101) (LOG(TES(NT)%PRIOR(L)), NH3_HAT(L), 
!!     &     LOG(TES(NT)%NH3(L)), L,  L = LLNT, 1, -1 )
! 101  FORMAT(1X,F16.8,1X,F16.8,1X,F16.8,1x,i3)
!
!         !--------------------------------------------------------------
!         ! Begin adjoint calculations 
!         !--------------------------------------------------------------
!
!         ! dkh debug
!         print*, 'DIFF , FORCE, diag OEI ' 
!         WRITE(6,102) (DIFF(L), FORCE(L), TES(NT)%OEI_VMR(L,L),
!     &       L = LLNT, 1, -1 )
! 102  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)
!
!         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
!         ADJ_DIFF(:) = FORCE(:) 
!
!         ! Adjoint of difference
!         DO L = 1, LLNT 
!            ADJ_NH3_HAT(L) =  ADJ_DIFF(L)
!         ENDDO 
!
!         ! adjoint of TES operator
!         DO L  = 1, LLNT
!            ADJ_NH3_PERT(L)    = 0d0
!            DO LL = 1, LLNT
!               ADJ_NH3_PERT(L) = ADJ_NH3_PERT(L) 
!     &                         + TES(NT)%AVK_VMR(L,LL) 
!!     &                         + TES(NT)%AVK_VMR(LL,L) 
!     &                         * ADJ_NH3_HAT(LL)
!           ENDDO
!         ENDDO
! 
!         ! Adjoint of x_m - x_a
!         DO L = 1, LLNT 
!           ! fwd code:
!           !NH3_PERT(L) = GC_NH3(L) - TES(NT)%PRIOR(L))
!           ! adj code:
!           ADJ_GC_NH3(L) = ADJ_NH3_PERT(L)
!         ENDDO
!
!         ! dkh debug
!         print*, 'ADJ_NH3_HAT, ADJ_NH3_PERT, ADJ_GC_NH3 '
!         WRITE(6,103) (ADJ_NH3_HAT(L), ADJ_NH3_PERT(L), ADJ_GC_NH3(L), 
!     &       L = LLNT, 1, -1 )
! 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)
!
!         ! adjoint of interpolation 
!         DO L  = 1, LLPAR
!            ADJ_GC_NH3_NATIVE(L) = 0d0 
!            DO LL = 1, LLNT
!               ADJ_GC_NH3_NATIVE(L) = ADJ_GC_NH3_NATIVE(L)
!     &                              + MAP(L,LL) * ADJ_GC_NH3(LL)
!            ENDDO
!         ENDDO
!         
!         ! Adjoint of unit conversion 
!         ADJ_GC_NH3_NATIVE(:) = ADJ_GC_NH3_NATIVE(:) * TCVV(IDTNH3)
!     &                        / AD(I,J,:) * 1d6
!
!         ! Pass adjoint back to adjoint tracer array
!         STT_ADJ(I,J,:,IDTNH3)  = STT_ADJ(I,J,:,IDTNH3)
!     &                          + ADJ_GC_NH3_NATIVE(:)
!
!         ! for FD test
!         ADJ(:) = ADJ_GC_NH3_NATIVE(:) * CHK_STT(I,J,:,IDTNH3)
!
!         ! dkh debug
!         print*, ' adj_stt force = ', ADJ_GC_NH3_NATIVE(:)
!
!         ! dkh debug
!         print*, 'ADJ_GC_NH3_NATIVE, ADJ_GC_NH3_NATIVE conv, ADJ_STT ' 
!         WRITE(6,104) (ADJ_GC_NH3_NATIVE(L) * AD(I,J,L) 
!     &                                      / 1d6 / TCVV(IDTNH3), 
!     &       ADJ_GC_NH3_NATIVE(L), STT_ADJ(I,J,L,IDTNH3),
!     &       L = LLPAR, 1, -1 )
! 104  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)
!
!
!         WRITE(101,112) ( TES(NT)%PRES(LL),  LL=LLNT,1,-1)
!         WRITE(102,112) ( GC_NH3(LL),        LL=LLNT,1,-1)
!         WRITE(103,112) ( TES(NT)%NH3(LL),   LL=LLNT,1,-1)
!         WRITE(104,112) ( TES(NT)%PRIOR(LL), LL=LLNT,1,-1)
!         WRITE(105,112) ( DIFF(LL),          LL=LLNT,1,-1)
!         WRITE(106,112) ( FORCE(LL),         LL=LLNT,1,-1)
!         WRITE(107,111) NT, LLNT
!         WRITE(108,112) ( ADJ_NH3_PERT(LL),  LL=LLNT,1,-1)
!         WRITE(109,112) ( ADJ_GC_NH3(LL),    LL=LLNT,1,-1)
!         WRITE(110,112) ( NH3_HAT(LL),       LL=LLNT,1,-1)
!         ! skip for real data
!         !WRITE(111,110) ( EXP(NH3_HAT_DBL(LL)),  LL=LLNT,1,-1)
! 110     FORMAT(F16.8,1X)
! 111     FORMAT(i4,1X,i4,1x)
! 112     FORMAT(D14.6,1X)
!
!
!
!      ENDDO  ! NT
!!!$OMP END PARALLEL DO
!
!      ! Update cost function 
!      ! for FD test
!      !COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))
!      COST_FUNC = SUM(NEW_COST(NTSTOP:NTSTART))
!
!      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
!      print*, ' TES contribution           = ', COST_FUNC - OLD_COST  
!
!      ! dkh debug
!      !IF ( NTSTOP == 1 ) THEN 
!      !   CALL MAKE_NH3_FILE( )
!      !ENDIF 
!
!      ! Return to calling program
!      END SUBROUTINE CALC_TES_NH3_FORCE_FD
!
!!------------------------------------------------------------------------------

      SUBROUTINE GET_NT_RANGE( GCNYMD, GCNHMS, NTSTART, NTSTOP )
!
!******************************************************************************
!  Subroutine GET_NT_RANGE 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) GCNYMD (INTEGER) : Current model YYYYMMDD 
!  (2 ) GCNHMS (INTEGER) : Current model HHMMSS  
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) NTSTART (INTEGER) : TES record number at which to start
!  (1 ) NTSTOP  (INTEGER) : TES record number at which to stop
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN)   :: GCNYMD 
      INTEGER, INTENT(IN)   :: GCNHMS
      INTEGER, INTENT(OUT)  :: NTSTART
      INTEGER, INTENT(OUT)  :: NTSTOP
    
      ! Local variables 
      INTEGER, SAVE         :: NTSAVE = MAXTES 
      !INTEGER, SAVE         :: NTSAVE = 148
      LOGICAL               :: FOUND_ALL_RECORDS 
      INTEGER               :: NTEST

      !=================================================================
      ! GET_NT_RANGE begins here!
      !=================================================================


      ! Initialize 
      FOUND_ALL_RECORDS  = .FALSE. 
      NTSTART            = 0
      NTSTOP             = 0

      print*, ' GET_NT_RANGE ', GCNYMD, GCNHMS
      print*, ' NTSAVE ', NTSAVE

      ! All records have been read already 
      IF ( NTSAVE == 0 ) THEN 

         print*, 'All records have been read already '
         RETURN 

      ! No records reached yet
      ELSEIF ( TES(NTSAVE)%NYMD < GCNYMD ) THEN 
           
      
         print*, 'No records reached yet'
         RETURN

      ! Model day matches day of existing records
      ELSEIF ( TES(NTSAVE)%NYMD == GCNYMD ) THEN 

         print*,' Model day matches day of existing records'
         IF ( TES(NTSAVE)%NHMS + 7000 < GCNHMS ) THEN 

            ! No records reached yet
            print*, 'But no records reached yet'
            RETURN

         ! Model hour (+/- 30 min) matches existing records
         ELSEIF ( TES(NTSAVE)%NHMS + 7000 >=  GCNHMS ) THEN 
      
            ! Starting record found
            NTSTART = NTSAVE   

            print*, ' Starting : TES%NYMD, TES%NHMS ', 
     &               TES(NTSTART)%NYMD, TES(NTSTART)%NHMS, NTSTART
 
            ! Now search forward to find stopping record
            NTEST = NTSTART

            DO WHILE ( FOUND_ALL_RECORDS == .FALSE. ) 
              
               ! Advance to the next record
               NTEST = NTEST - 1  
           
               ! Stop if we reach the earliest available record 
               IF ( NTEST == 0 ) THEN 
           
                  NTSTOP            = NTEST + 1
                  FOUND_ALL_RECORDS = .TRUE.

                  print*, ' Records found '
                  print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP

                  ! Reset NTSAVE 
                  NTSAVE = NTEST

               ! When the combined test date rounded up to the nearest
               ! half hour is smaller than the current model date, the 
               ! stopping record has been passed. 
               ELSEIF (  TES(NTEST)%NHMS + 7000 <  GCNHMS  ) THEN 
          
                  print*, ' Testing : TES%NYMD, TES%NHMS ', 
     &                  TES(NTEST)%NYMD, TES(NTEST)%NHMS, NTEST
 
                  NTSTOP            = NTEST + 1 
                  FOUND_ALL_RECORDS = .TRUE. 

                  print*, ' Records found '
                  print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP
                  
                  ! Reset NTSAVE 
                  NTSAVE = NTEST 

               ! When test day is earlier than the model day and the 
               ! hour is less than 23:30, the  stopping record has been passed. 
               ELSEIF (  TES(NTEST)%NYMD  <  GCNYMD .and. 
     &                   TES(NTEST)%NHMS  <  233000       ) THEN 
          
                  NTSTOP            = NTEST + 1 
                  FOUND_ALL_RECORDS = .TRUE. 

                  print*, ' Records found '
                  print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP
                  
                  ! Reset NTSAVE 
                  NTSAVE = NTEST 
               ELSE 
                   print*, ' still looking ', NTEST 
                  
               ENDIF 
                 
            ENDDO 
 
         ENDIF 

      ELSE

         CALL ERROR_STOP('problem', 'GET_NT_RANGE' ) 

      ENDIF 

      ! Return to calling program
      END SUBROUTINE GET_NT_RANGE

!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP,
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )
     *         RESULT      ( HINTERPZ )
!
!******************************************************************************
!  Function GET_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!     
!  NOTES:
!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod. 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

      ! Arguments
      INTEGER            :: LGC_TOP, LTM_TOP
      TYPE (XPLEX)             :: GC_PRESC(LGC_TOP)
      TYPE (XPLEX)             :: TM_PRESC(LTM_TOP) 
      TYPE (XPLEX)             :: GC_SURFP
      TYPE (XPLEX)             :: TM_SURFP
 
      ! Return value 
      TYPE (XPLEX)             :: HINTERPZ(LGC_TOP, LTM_TOP)

      ! Local variables 
      INTEGER  :: LGC, LTM
      TYPE (XPLEX)   :: DIFF, DELTA_SURFP
      TYPE (XPLEX)   :: LOW, HI

      !=================================================================
      ! GET_HINTERPZ_2 begins here!
      !=================================================================

      HINTERPZ(:,:) = 0D0 
  
!      ! Rescale GC grid according to TM surface pressure
!!         p1_A =     (a1 + b1 (ps_A - PTOP))
!!         p2_A =     (a2 + b2 (ps_A - PTOP))
!!         p1_B =     (a + b (ps_B - PTOP))
!!         p2_B =    *(a + b (ps_B - PTOP))
!!         pc_A = 0.5(a1+a2 +(b1+b2)*(ps_A - PTOP))
!!         pc_B = 0.5(a1+a2 +(b1+b2)*(ps_B - PTOP))
!!         pc_B - pc_A = 0.5(b1_b2)(ps_B-ps_A)
!!         pc_B = 0.5(b1_b2)(ps_B-ps_A) + pc_A
!      DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
!
!      DO LGC = 1, LGC_TOP
!         GC_PRESC(LGC) = ( GET_BP(LGC) + GET_BP(LGC+1))
!     &               * DELTA_SURFP + GC_PRESC(LGC)
!         IF (GC_PRESC(LGC) < 0) THEN 
!            CALL ERROR_STOP( 'highly unlikey', 
!     &                       'read_sciano2_mod.f')
!         ENDIF 
!
!      ENDDO 
      

      ! Loop over each pressure level of TM grid
      DO LTM = 1, LTM_TOP
 
         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, LGC_TOP - 1

            LOW = GC_PRESC(LGC+1)
            HI  = GC_PRESC(LGC)
            IF (LGC == 0) HI = TM_SURFP

            ! Linearly interpolate value on the LTM grid 
            IF ( TM_PRESC(LTM) <= HI .and. 
     &           TM_PRESC(LTM)  > LOW) THEN 

               DIFF                = HI - LOW  
               HINTERPZ(LGC+1,LTM) = ( HI - TM_PRESC(LTM)  ) / DIFF
               HINTERPZ(LGC  ,LTM) = ( TM_PRESC(LTM) - LOW ) / DIFF


            ENDIF 
 
            ! dkh debug
            !print*, 'LGC,LTM,HINT', LGC, LTM, HINTERPZ(LGC,LTM)

          ENDDO
       ENDDO

       ! Correct for case where the lowest TES level pressure is lower than the
       ! lowest GC layer pressure.  In this case, just 1:1 map. 
 
       ! Bug fix:  a more general version allows for multiples TES pressure
       ! levels to exist below the lowest GC pressure.  (dm, dkh, 09/30/10) 
       ! OLD code:
       !IF ( TM_PRESC(1) > GC_PRESC(1) ) THEN
       !   HINTERPZ(1,1)         = 1D0 
       !   HINTERPZ(2:LGC_TOP,1) = 0D0 
       !ENDIF
       ! New code:
       ! Loop over each pressure level of TM grid
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             HINTERPZ(1,LTM)         = 1D0
             HINTERPZ(2:LGC_TOP,LTM) = 0D0
          ENDIF
       ENDDO 

      ! Return to calling program
      END FUNCTION GET_INTMAP

!------------------------------------------------------------------------------
      SUBROUTINE MAKE_NH3_FILE(  )
!
!******************************************************************************
!  Subroutine MAKE_NH3_FILE saves NH3 profiles that correspond to time and
!  place of TES NH3 obs. (dkh, 03/01/09) 
!
!  Module variables as Input:
!  ============================================================================
!  (1 ) NH3_SAVE (TYPE (XPLEX)) : NH3 profiles                             [ppmv]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE ERROR_MOD,        ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : EXPAND_DATE

#     include "CMN_SIZE"    ! Size params
      
      ! Local variables    
      INTEGER              :: I, J, I0, J0, L, NT
      CHARACTER(LEN=120)   :: FILENAME
      TYPE (XPLEX)               :: DAT(1,LLPAR,MAXTES)
      INTEGER, PARAMETER   :: IUN = 88 
      
      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      !=================================================================
      ! MAKE_NH3_FILE begins here!
      !=================================================================
      
      FILENAME = TRIM( 'nh3.bpch' )
      
      ! Append data directory prefix
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )
      
      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'NH3 profile '
      CATEGORY = 'IJ-AVE-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UNIT     = 'ppmv'

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================


      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_NH3_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IUN, FILENAME, TITLE )

      ! Temporarily store data in DAT as REAL4
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( NT ) 
      DO NT = 1, MAXTES

         DAT(1,:,NT) = DCMPLX(NH3_SAVE(:,NT))

      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IUN,       MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UNIT,      1d0,       1d0,       RESERVED,
     &            1,         LLPAR,     MAXTES,     I0+1,
     &            J0+1,      1,         DAT )

      ! Close file
      CLOSE( IUN )        

      print*, ' NH3_SAVE sum write = ', SUM(NH3_SAVE(:,:))

      ! Return to calling program
      END SUBROUTINE MAKE_NH3_FILE

!------------------------------------------------------------------------------
      SUBROUTINE READ_NH3_FILE(  )
!
!******************************************************************************
!  Subroutine READ_NH3_FILE reads the GC modeled NH3 profiles that correspond
!  to the TES NH3 times and locations. (dkh, 03/01/09) 
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD,         ONLY : READ_BPCH2
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR


#     include "CMN_SIZE"          ! Size parameters

      ! Local variables
      TYPE (XPLEX)                     :: DAT(1,LLPAR,MAXTES)
      CHARACTER(LEN=255)         :: FILENAME

      !=================================================================
      ! READ_USA_MASK begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( ADJTMP_DIR )           //
     &           'nh3.bpch'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_NH3_FILE: Reading ', a )
      
      
      ! USA mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'IJ-AVE-$', 1,
     &                 1d0,            1,        LLPAR, 
     &                 MAXTES,    DAT,      QUIET=.TRUE. )
      
      ! Cast to TYPE (XPLEX)
      NH3_SAVE(:,:) = DAT(1,:,:)
      
      print*, ' NH3_SAVE sum read = ', SUM(NH3_SAVE(:,:))

      ! Return to calling program
      END SUBROUTINE READ_NH3_FILE

!-----------------------------------------------------------------------------
      FUNCTION GET_DOUBLED_NH3( NYMD, NHMS, LON, LAT ) RESULT( NH3_DBL )
!
!******************************************************************************
!  Subroutine GET_DOUBLED_NH3 reads and returns the nh3 profiles from 
!  model run with doubled emissions. (dkh, 11/08/09) 
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD,         ONLY : READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR
      USE GRID_MOD,          ONLY : GET_IJ
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_TAU


#     include "CMN_SIZE"          ! Size parameters

      ! Arguments   
      INTEGER                    :: NYMD, NHMS
      TYPE (XPLEX)                     :: LON,  LAT
      
      ! Function arg 
      TYPE (XPLEX)                     :: NH3_DBL(LLPAR)

      ! Local variables
      TYPE (XPLEX)                     :: DAT(144,91,20)
      CHARACTER(LEN=255)         :: FILENAME
      INTEGER                    :: IIJJ(2)

      !=================================================================
      ! GET_DOUBLED_NH3 begins here!
      !=================================================================

      ! filename
      FILENAME = 'nh3.YYYYMMDD.hhmm'

      ! Expand filename
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      ! Full path to file
      FILENAME = TRIM( DATA_DIR )           //
     &           'doubled_nh3/'             // 
     &           TRIM( FILENAME )           //
     &           TRIM( '00'     )           

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_DOUBLED_NH3: Reading ', a )
      
      ! dkh debug
      print*, ' GET_TAU() = ', GET_TAU()
      
      ! Get data
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 29, 
     &                 GET_TAU(), 144,      91, 
     &                 20,        DAT,      QUIET=.FALSE. )
      
      ! Now use GET_IJ in grid_mod.f (dkh, 02/16/11) 
      !IIJJ = GET_IJ_2x25( LON, LAT )
      IIJJ = GET_IJ( LON, LAT )

 
      print*, ' found doubled in I/J = ', IIJJ

      ! just the column for the present location, and convert ppb to ppm
      NH3_DBL(1:20)     = DCMPLX(DAT(IIJJ(1),IIJJ(2),:)) / 1000d0 
      NH3_DBL(21:LLPAR) = 0d0 
     
      print*, ' NH3_DBL = ', NH3_DBL
 
      ! Return to calling program
      END FUNCTION GET_DOUBLED_NH3

! Now we use GET_IJ in grid_mod.f (dkh, 02/16/11) 
!!------------------------------------------------------------------------------
!      FUNCTION GET_IJ_2x25( LON, LAT ) RESULT ( IIJJ )
!
!!
!!******************************************************************************
!!  Subroutine GET_IJ_2x25 returns I and J index from the 2 x 2.5 grid for a 
!!  LON, LAT coord. (dkh, 11/08/09) 
!! 
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) LON (TYPE (XPLEX)) : Longitude                          [degrees]
!!  (2 ) LAT (TYPE (XPLEX)) : Latitude                           [degrees]
!!     
!!  Function result
!!  ============================================================================
!!  (1 ) IIJJ(1) (INTEGER) : Long index                    [none]
!!  (2 ) IIJJ(2) (INTEGER) : Lati index                    [none]
!!     
!!  NOTES:
!!
!!******************************************************************************
!!     
!      ! Reference to f90 modules
!      USE ERROR_MOD,    ONLY : ERROR_STOP
!
!      ! Arguments
!      TYPE (XPLEX)    :: LAT, LON
!      
!      ! Return
!      INTEGER :: I, J, IIJJ(2)
!      
!      ! Local variables 
!      TYPE (XPLEX)              :: TLON, TLAT, DLON, DLAT
!      TYPE (XPLEX),  PARAMETER  :: DISIZE = 2.5d0
!      TYPE (XPLEX),  PARAMETER  :: DJSIZE = 2.0d0
!      INTEGER, PARAMETER  :: IIMAX  = 144
!      INTEGER, PARAMETER  :: JJMAX  = 91
!      
!      
!      !=================================================================
!      ! GET_IJ_2x25 begins here!
!      !=================================================================
!
!      TLON = 180d0 + LON + DISIZE
!      TLAT =  90d0 + LAT + DJSIZE
!      
!      I = TLON / DISIZE
!      J = TLAT / DJSIZE
!
!      
!      IF ( TLON / DISIZE - REAL(I)  >= 0.5d0 ) THEN
!         I = I + 1
!      ENDIF
!      
!      IF ( TLAT / DJSIZE - REAL(J)  >= 0.5d0 ) THEN
!         J = J + 1
!      ENDIF
!
!      
!      ! Longitude wraps around
!      !IF ( I == 73 ) I = 1 
!      IF ( I == ( IIMAX + 1 ) ) I = 1
!      
!      ! Check for impossible values 
!      IF ( I > IIMAX .or. J > JJMAX .or. 
!     &     I < 1     .or. J < 1          ) THEN
!         CALL ERROR_STOP('Error finding grid box', 'GET_IJ_2x25')
!      ENDIF
!      
!      IIJJ(1) = I
!      IIJJ(2) = J
!      
!      ! Return to calling program
!      END FUNCTION GET_IJ_2x25
!
!!-----------------------------------------------------------------------------

      SUBROUTINE READ_TES_BLVMR( FILENAME, LAT, LON,  
     &                           BLVMR,    BLVMR_GC          )
!
!******************************************************************************
!  Subroutine READ_TES_BLVMR reads the BLVMR values in a TES netcdf file. 
!  (dkh, 01/20/11) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME   (TYPE (XPLEX)) : TES observation filename to read
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) LAT        (TYPE (XPLEX)) : TES retrieval latitude
!  (2 ) LON        (TYPE (XPLEX)) : TES retrieval longitude
!  (3 ) BLVMR      (TYPE (XPLEX)) : TES retrieval BLVMR                  [ppmv]
!  (4 ) BLVMR_GC   (TYPE (XPLEX)) : GEOS-Chem BLVMR                      [ppmv]
!     
!  NOTES:
! 
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE NETCDF 

#     include "CMN_SIZE" 

      ! Arguments
      CHARACTER(LEN=255), INTENT(IN)  :: FILENAME
      TYPE (XPLEX),             INTENT(OUT) :: LAT
      TYPE (XPLEX),             INTENT(OUT) :: LON
      TYPE (XPLEX),             INTENT(OUT) :: BLVMR
      TYPE (XPLEX),             INTENT(OUT) :: BLVMR_GC
    
      ! Local variables
      CHARACTER(LEN=255)              :: READ_FILENAME
      CHARACTER(LEN=5)                :: TMP
      TYPE (XPLEX)                          :: BLVMR_TMP(3)
      TYPE (XPLEX)                          :: BLVMR_GC_TMP(3)
      INTEGER                         :: QFLAG
      INTEGER                         :: DFLAG

      ! netCDF id's 
      INTEGER                         :: FID, VARID, DIMID, ATTID
      INTEGER                         :: NBL
    
      ! Loop indexes, and error handling.
      INTEGER                         :: L,   LL

      ! dkh debug
      INTEGER :: dim_ids(10), n_dims

      !=================================================================
      ! READ_TES_BLVMR begins here!
      !=================================================================

      ! Initialize
      BLVMR_GC  = -999d0
      BLVMR     = -999d0

      ! Construct complete filename 
      READ_FILENAME = TRIM( DATA_DIR ) // 
     &                TRIM( '../TES_NH3/' )   //
     &                TRIM( 'tes_nh3_gs_January_2006_for_paper/' ) //
     &                TRIM( FILENAME ) 

      WRITE(6,*) '    - READ_TES_NH3_OBS: reading file: ', READ_FILENAME


      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )
   
      ! READ quality flag QFLAG
      ! 1 = successful
      ! 0 = failed.  For reason why, see DFLAG. 
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL,  "Quality_Flag", 
     &                          QFLAG  ), 12 )

      ! READ diagnostic flag DFLAG
      !  0 = converged, but DOFS < 0.5, though thermal contrast ok 
      ! -1 = converged, but DOFS < 0.5 & thermal contrast poor
      ! -2 = did not converge
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL, "Diagnostic_Flag",  
     &                          DFLAG  ), 14 )

      ! READ latitude LAT 
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL, "Latitude",
     &                          LAT    ), 15 )

      ! READ longitude LON
      CALL CHECK( NF90_GET_ATT( FID, NF90_GLOBAL, "Longitude",
     &                          LON    ), 16 )  

      ! READ BLVMR only for good retrievals, fill othewise
      IF ( QFLAG == 1 .and. DFLAG == 1 ) THEN

         ! dkh debug
         print*, ' quality retv; look for BLVMR '

         ! get NBL
         CALL CHECK( NF90_INQ_DIMID ( FID, "nbl", DIMID      ), 201 )
         CALL CHECK( NF90_INQUIRE_DIMENSION( FID,DIMID, TMP, NBL), 212 )
 
         ! Get BLVMR_TMP
         CALL CHECK( NF90_INQ_VARID( FID, "blvmr",  VARID   ), 21 )
         CALL CHECK( NF90_GET_VAR  ( FID, VARID, BLVMR_TMP(1:NBL) ),
     &      222 )

         ! READ BLVMR_GC only if BLVMR not = -999
         IF ( BLVMR_TMP(1) > 0 ) THEN  
 
            ! Get BLVMR_GC
            CALL CHECK( NF90_INQ_VARID( FID, "blvmr_gc",  VARID ), 
     &                                                             232 )

            CALL CHECK( NF90_GET_VAR  ( FID, VARID, 
     &         BLVMR_GC_TMP(1:NBL) ), 333 )
           
            ! It is possible to have two BLVMRs when DOF > 1.2.
            IF ( NBL == 1 ) THEN
            
               ! Keep only the first one. 
               BLVMR        = BLVMR_TMP(1)
               
            ELSEIF ( NBL > 1 ) THEN
               ! Keep the last one (dkh, 01/21/11) 
               BLVMR        = BLVMR_TMP(NBL)

            ENDIF

  
         ELSE 

            ! Fill value 
            BLVMR    = -999D0 
            BLVMR_GC = -999D0 
            
         ENDIF 
      ELSE 

         ! Fill value 
         BLVMR    = -999D0 
         BLVMR_GC = -999D0 

      ENDIF 

      print*, ' BLVMR = ', BLVMR, BLVMR_GC, LAT, LON, FILENAME 

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )
    
      ! Return to calling program
      END SUBROUTINE READ_TES_BLVMR
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_TES_BLVMR( )
!
!******************************************************************************
!  Subroutine MAKE_TES_BLVMR makes a binary punch file out of BLVMR values
!  (dkh, 01/20/11) 
! 
!  NOTES:
!  
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_IJ, GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : GET_TAU

#     include "CMN_SIZE"

      ! Arguments
    
      ! Local variables 
      INTEGER :: I, J, I0, J0, IIJJ(2), NTES 
      CHARACTER(LEN=120) :: WRITE_FILENAME
      CHARACTER(LEN=255) :: FILENAME
      TYPE (XPLEX)  :: BLVMR_BAR(IIPAR,JJPAR),  BLVMR_GC_BAR(IIPAR,JJPAR)
      TYPE (XPLEX)  :: DELTA_BLVMR,           DELTA_BLVMR_GC
      INTEGER :: N(IIPAR,JJPAR)

      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1 
      INTEGER              :: IU_BPCH

      ! Variables read from file
      TYPE (XPLEX)   :: LAT, LON
      TYPE (XPLEX)   :: BLVMR
      TYPE (XPLEX)   :: BLVMR_GC

      !=================================================================
      ! MAKE_TES_BLVMR begins here!
      !=================================================================
      BLVMR_BAR(:,:)     = 0d0
      BLVMR_GC_BAR(:,:)  = 0d0
      N(:,:)             = 0d0

 
      DO NTES = 1, MAXTES 

         FILENAME = TES(NTES)%FILENAME

         WRITE(6,*) ' MAKE_TES_BLVMR: reading file: ', FILENAME, NTES 
 
         CALL READ_TES_BLVMR( FILENAME, LAT, LON, BLVMR, BLVMR_GC ) 
         
         ! dkh debug
         print*, ' LAT      = ', LAT      
         print*, ' LON      = ', LON      
         print*, ' BLVMR    = ', BLVMR
         print*, ' BLVMR_GC = ', BLVMR_GC

         ! Quality check
         IF (
     &         BLVMR    < 0 .or. 
     &         BLVMR_GC < 0       ) CYCLE
 
         ! Get grid box 
         IIJJ = GET_IJ(DCMPLX(LON),DCMPLX(LAT))
         I    = IIJJ(1)
         J    = IIJJ(2)
   
         ! Update local count 
         N(I,J) = N(I,J) + 1
      
         ! Update mean
         DELTA_BLVMR    = BLVMR    - BLVMR_BAR(I,J)
         DELTA_BLVMR_GC = BLVMR_GC - BLVMR_GC_BAR(I,J)
   
         BLVMR_BAR(I,J)     = BLVMR_BAR(I,J)    + DELTA_BLVMR    /N(I,J)
         BLVMR_GC_BAR(I,J)  = BLVMR_GC_BAR(I,J) + DELTA_BLVMR_GC /N(I,J)
      
      ENDDO

      WRITE(6,*) 'Done reading TES data files '
      WRITE(6,*) 'Number of good data points found: ', SUM(N(:,:))

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'BLVMR   data file ' 
      CATEGORY = 'BLVMR'    
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UNIT     = 'ppm'

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the file for output -- binary punch format
      !=================================================================

      ! Add ADJ_DIR prefix to filename
      WRITE_FILENAME = TRIM( 'blvmr.bpch' ) 

      WRITE( 6, 100 ) TRIM( WRITE_FILENAME )
 100  FORMAT( '     - MAKE_TES_BLVMR: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, WRITE_FILENAME, TITLE )

      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,     I0+1,
     &            J0+1,      1,      DCMPLX(BLVMR_BAR) )

      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  2,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,     I0+1,
     &            J0+1,      1,      DCMPLX(BLVMR_GC_BAR) )


      ! Close file
      CLOSE( IU_BPCH )

      ! Return to calling program
      END SUBROUTINE MAKE_TES_BLVMR

!------------------------------------------------------------------------------
      SUBROUTINE INIT_TES_NH3
!
!*****************************************************************************
!  Subroutine INIT_TES_NH3 deallocates all module arrays.  (dkh, 02/15/09) 
!        
!  NOTES:
!
!******************************************************************************
!     
      USE ERROR_MOD,  ONLY : ALLOC_ERR

#     include "CMN_SIZE"   ! IIPAR, JJPAR      
       
      ! Local variables
      INTEGER :: AS

      !=================================================================      
      ! INIT_TES_NH3 begins here
      !================================================================= 

      ! dkh debug
      print*, ' INIT_TES_NH3'

      ALLOCATE( NH3_SAVE( LLPAR, MAXTES ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3_SAVE' ) 
      NH3_SAVE = 0d0

      ! tes_nh3_gs_July_2009_for_paper
      TES(1)%NYMD = 20090701
      TES(2)%NYMD = 20090701
      TES(3)%NYMD = 20090701
      TES(4)%NYMD = 20090701
      TES(5)%NYMD = 20090701
      TES(6)%NYMD = 20090701
      TES(7)%NYMD = 20090701
      TES(8)%NYMD = 20090701
      TES(9)%NYMD = 20090701
      TES(10)%NYMD = 20090701
      TES(11)%NYMD = 20090701
      TES(12)%NYMD = 20090701
      TES(13)%NYMD = 20090701
      TES(14)%NYMD = 20090701
      TES(15)%NYMD = 20090701
      TES(16)%NYMD = 20090701
      TES(17)%NYMD = 20090701
      TES(18)%NYMD = 20090701
      TES(19)%NYMD = 20090701
      TES(20)%NYMD = 20090701
      TES(21)%NYMD = 20090701
      TES(22)%NYMD = 20090701
      TES(23)%NYMD = 20090701
      TES(24)%NYMD = 20090701
      TES(25)%NYMD = 20090701
      TES(26)%NYMD = 20090701
      TES(27)%NYMD = 20090701
      TES(28)%NYMD = 20090701
      TES(29)%NYMD = 20090701
      TES(30)%NYMD = 20090701
      TES(31)%NYMD = 20090701
      TES(32)%NYMD = 20090701
      TES(33)%NYMD = 20090701
      TES(34)%NYMD = 20090701
      TES(35)%NYMD = 20090701
      TES(36)%NYMD = 20090701
      TES(37)%NYMD = 20090701
      TES(38)%NYMD = 20090701
      TES(39)%NYMD = 20090701
      TES(40)%NYMD = 20090701
      TES(41)%NYMD = 20090701
      TES(42)%NYMD = 20090701
      TES(43)%NYMD = 20090701
      TES(44)%NYMD = 20090701
      TES(45)%NYMD = 20090701
      TES(46)%NYMD = 20090701
      TES(47)%NYMD = 20090701
      TES(48)%NYMD = 20090701
      TES(49)%NYMD = 20090701
      TES(50)%NYMD = 20090701
      TES(51)%NYMD = 20090701
      TES(52)%NYMD = 20090701
      TES(53)%NYMD = 20090701
      TES(54)%NYMD = 20090701
      TES(55)%NYMD = 20090701
      TES(56)%NYMD = 20090701
      TES(57)%NYMD = 20090701
      TES(58)%NYMD = 20090701
      TES(59)%NYMD = 20090701
      TES(60)%NYMD = 20090701
      TES(61)%NYMD = 20090701
      TES(62)%NYMD = 20090701
      TES(63)%NYMD = 20090701
      TES(64)%NYMD = 20090701
      TES(65)%NYMD = 20090701
      TES(66)%NYMD = 20090701
      TES(67)%NYMD = 20090701
      TES(68)%NYMD = 20090701
      TES(69)%NYMD = 20090701
      TES(70)%NYMD = 20090701
      TES(71)%NYMD = 20090701
      TES(72)%NYMD = 20090701
      TES(73)%NYMD = 20090701
      TES(74)%NYMD = 20090701
      TES(75)%NYMD = 20090701
      TES(76)%NYMD = 20090701
      TES(77)%NYMD = 20090701
      TES(78)%NYMD = 20090701
      TES(79)%NYMD = 20090701
      TES(80)%NYMD = 20090701
      TES(81)%NYMD = 20090701
      TES(82)%NYMD = 20090701
      TES(83)%NYMD = 20090701
      TES(84)%NYMD = 20090701
      TES(85)%NYMD = 20090701
      TES(86)%NYMD = 20090701
      TES(87)%NYMD = 20090701
      TES(88)%NYMD = 20090701
      TES(89)%NYMD = 20090701
      TES(90)%NYMD = 20090701
      TES(91)%NYMD = 20090701
      TES(92)%NYMD = 20090701
      TES(93)%NYMD = 20090701
      TES(94)%NYMD = 20090701
      TES(95)%NYMD = 20090701
      TES(96)%NYMD = 20090701
      TES(97)%NYMD = 20090701
      TES(98)%NYMD = 20090701
      TES(99)%NYMD = 20090701
      TES(100)%NYMD = 20090701
      TES(101)%NYMD = 20090701
      TES(102)%NYMD = 20090701
      TES(103)%NYMD = 20090701
      TES(104)%NYMD = 20090701
      TES(105)%NYMD = 20090701
      TES(106)%NYMD = 20090701
      TES(107)%NYMD = 20090701
      TES(108)%NYMD = 20090701
      TES(109)%NYMD = 20090701
      TES(110)%NYMD = 20090701
      TES(111)%NYMD = 20090701
      TES(112)%NYMD = 20090701
      TES(113)%NYMD = 20090701
      TES(114)%NYMD = 20090701
      TES(115)%NYMD = 20090701
      TES(116)%NYMD = 20090701
      TES(117)%NYMD = 20090701
      TES(118)%NYMD = 20090701
      TES(119)%NYMD = 20090701
      TES(120)%NYMD = 20090701
      TES(121)%NYMD = 20090701
      TES(122)%NYMD = 20090701
      TES(123)%NYMD = 20090701
      TES(124)%NYMD = 20090701
      TES(125)%NYMD = 20090701
      TES(126)%NYMD = 20090701
      TES(127)%NYMD = 20090701
      TES(128)%NYMD = 20090701
      TES(129)%NYMD = 20090701
      TES(130)%NYMD = 20090701
      TES(131)%NYMD = 20090701
      TES(132)%NYMD = 20090701
      TES(133)%NYMD = 20090701
      TES(134)%NYMD = 20090701
      TES(135)%NYMD = 20090701
      TES(136)%NYMD = 20090701
      TES(137)%NYMD = 20090701
      TES(138)%NYMD = 20090701
      TES(139)%NYMD = 20090701
      TES(140)%NYMD = 20090701
      TES(141)%NYMD = 20090701
      TES(142)%NYMD = 20090701
      TES(143)%NYMD = 20090701
      TES(144)%NYMD = 20090701
      TES(145)%NYMD = 20090701
      TES(146)%NYMD = 20090701
      TES(147)%NYMD = 20090701
      TES(148)%NYMD = 20090701
      TES(149)%NYMD = 20090702
      TES(150)%NYMD = 20090702
      TES(151)%NYMD = 20090702
      TES(152)%NYMD = 20090702
      TES(153)%NYMD = 20090702
      TES(154)%NYMD = 20090702
      TES(155)%NYMD = 20090702
      TES(156)%NYMD = 20090702
      TES(157)%NYMD = 20090702
      TES(158)%NYMD = 20090702
      TES(159)%NYMD = 20090702
      TES(160)%NYMD = 20090702
      TES(161)%NYMD = 20090702
      TES(162)%NYMD = 20090702
      TES(163)%NYMD = 20090702
      TES(164)%NYMD = 20090702
      TES(165)%NYMD = 20090702
      TES(166)%NYMD = 20090702
      TES(167)%NYMD = 20090702
      TES(168)%NYMD = 20090702
      TES(169)%NYMD = 20090702
      TES(170)%NYMD = 20090702
      TES(171)%NYMD = 20090702
      TES(172)%NYMD = 20090702
      TES(173)%NYMD = 20090702
      TES(174)%NYMD = 20090702
      TES(175)%NYMD = 20090702
      TES(176)%NYMD = 20090702
      TES(177)%NYMD = 20090702
      TES(178)%NYMD = 20090702
      TES(179)%NYMD = 20090702
      TES(180)%NYMD = 20090702
      TES(181)%NYMD = 20090702
      TES(182)%NYMD = 20090702
      TES(183)%NYMD = 20090702
      TES(184)%NYMD = 20090702
      TES(185)%NYMD = 20090702
      TES(186)%NYMD = 20090702
      TES(187)%NYMD = 20090702
      TES(188)%NYMD = 20090702
      TES(189)%NYMD = 20090702
      TES(190)%NYMD = 20090702
      TES(191)%NYMD = 20090702
      TES(192)%NYMD = 20090702
      TES(193)%NYMD = 20090702
      TES(194)%NYMD = 20090702
      TES(195)%NYMD = 20090702
      TES(196)%NYMD = 20090702
      TES(197)%NYMD = 20090702
      TES(198)%NYMD = 20090702
      TES(199)%NYMD = 20090702
      TES(200)%NYMD = 20090702
      TES(201)%NYMD = 20090702
      TES(202)%NYMD = 20090702
      TES(203)%NYMD = 20090702
      TES(204)%NYMD = 20090702
      TES(205)%NYMD = 20090702
      TES(206)%NYMD = 20090702
      TES(207)%NYMD = 20090702
      TES(208)%NYMD = 20090702
      TES(209)%NYMD = 20090702
      TES(210)%NYMD = 20090702
      TES(211)%NYMD = 20090702
      TES(212)%NYMD = 20090702
      TES(213)%NYMD = 20090702
      TES(214)%NYMD = 20090702
      TES(215)%NYMD = 20090702
      TES(216)%NYMD = 20090702
      TES(217)%NYMD = 20090702
      TES(218)%NYMD = 20090702
      TES(219)%NYMD = 20090702
      TES(220)%NYMD = 20090702
      TES(221)%NYMD = 20090702
      TES(222)%NYMD = 20090702
      TES(223)%NYMD = 20090702
      TES(224)%NYMD = 20090702
      TES(225)%NYMD = 20090702
      TES(226)%NYMD = 20090702
      TES(227)%NYMD = 20090702
      TES(228)%NYMD = 20090702
      TES(229)%NYMD = 20090702
      TES(230)%NYMD = 20090702
      TES(231)%NYMD = 20090702
      TES(232)%NYMD = 20090702
      TES(233)%NYMD = 20090702
      TES(234)%NYMD = 20090702
      TES(235)%NYMD = 20090702
      TES(236)%NYMD = 20090702
      TES(237)%NYMD = 20090702
      TES(238)%NYMD = 20090702
      TES(239)%NYMD = 20090702
      TES(240)%NYMD = 20090702
      TES(241)%NYMD = 20090702
      TES(242)%NYMD = 20090702
      TES(243)%NYMD = 20090702
      TES(244)%NYMD = 20090702
      TES(245)%NYMD = 20090702
      TES(246)%NYMD = 20090702
      TES(247)%NYMD = 20090702
      TES(248)%NYMD = 20090702
      TES(249)%NYMD = 20090702
      TES(250)%NYMD = 20090702
      TES(251)%NYMD = 20090702
      TES(252)%NYMD = 20090702
      TES(253)%NYMD = 20090702
      TES(254)%NYMD = 20090702
      TES(255)%NYMD = 20090702
      TES(256)%NYMD = 20090702
      TES(257)%NYMD = 20090702
      TES(258)%NYMD = 20090702
      TES(259)%NYMD = 20090702
      TES(260)%NYMD = 20090702
      TES(261)%NYMD = 20090702
      TES(262)%NYMD = 20090702
      TES(263)%NYMD = 20090702
      TES(264)%NYMD = 20090702
      TES(265)%NYMD = 20090702
      TES(266)%NYMD = 20090702
      TES(267)%NYMD = 20090702
      TES(268)%NYMD = 20090702
      TES(269)%NYMD = 20090702
      TES(270)%NYMD = 20090702
      TES(271)%NYMD = 20090702
      TES(272)%NYMD = 20090702
      TES(273)%NYMD = 20090702
      TES(274)%NYMD = 20090703
      TES(275)%NYMD = 20090703
      TES(276)%NYMD = 20090703
      TES(277)%NYMD = 20090703
      TES(278)%NYMD = 20090703
      TES(279)%NYMD = 20090703
      TES(280)%NYMD = 20090703
      TES(281)%NYMD = 20090703
      TES(282)%NYMD = 20090703
      TES(283)%NYMD = 20090703
      TES(284)%NYMD = 20090703
      TES(285)%NYMD = 20090703
      TES(286)%NYMD = 20090703
      TES(287)%NYMD = 20090703
      TES(288)%NYMD = 20090703
      TES(289)%NYMD = 20090703
      TES(290)%NYMD = 20090703
      TES(291)%NYMD = 20090703
      TES(292)%NYMD = 20090703
      TES(293)%NYMD = 20090703
      TES(294)%NYMD = 20090703
      TES(295)%NYMD = 20090703
      TES(296)%NYMD = 20090703
      TES(297)%NYMD = 20090703
      TES(298)%NYMD = 20090703
      TES(299)%NYMD = 20090703
      TES(300)%NYMD = 20090703
      TES(301)%NYMD = 20090703
      TES(302)%NYMD = 20090703
      TES(303)%NYMD = 20090703
      TES(304)%NYMD = 20090703
      TES(305)%NYMD = 20090703
      TES(306)%NYMD = 20090703
      TES(307)%NYMD = 20090703
      TES(308)%NYMD = 20090703
      TES(309)%NYMD = 20090703
      TES(310)%NYMD = 20090703
      TES(311)%NYMD = 20090703
      TES(312)%NYMD = 20090703
      TES(313)%NYMD = 20090703
      TES(314)%NYMD = 20090703
      TES(315)%NYMD = 20090703
      TES(316)%NYMD = 20090703
      TES(317)%NYMD = 20090703
      TES(318)%NYMD = 20090703
      TES(319)%NYMD = 20090703
      TES(320)%NYMD = 20090703
      TES(321)%NYMD = 20090703
      TES(322)%NYMD = 20090703
      TES(323)%NYMD = 20090703
      TES(324)%NYMD = 20090703
      TES(325)%NYMD = 20090703
      TES(326)%NYMD = 20090703
      TES(327)%NYMD = 20090703
      TES(328)%NYMD = 20090703
      TES(329)%NYMD = 20090703
      TES(330)%NYMD = 20090703
      TES(331)%NYMD = 20090703
      TES(332)%NYMD = 20090703
      TES(333)%NYMD = 20090703
      TES(334)%NYMD = 20090703
      TES(335)%NYMD = 20090703
      TES(336)%NYMD = 20090703
      TES(337)%NYMD = 20090703
      TES(338)%NYMD = 20090703
      TES(339)%NYMD = 20090703
      TES(340)%NYMD = 20090703
      TES(341)%NYMD = 20090703
      TES(342)%NYMD = 20090703
      TES(343)%NYMD = 20090703
      TES(344)%NYMD = 20090703
      TES(345)%NYMD = 20090703
      TES(346)%NYMD = 20090703
      TES(347)%NYMD = 20090703
      TES(348)%NYMD = 20090703
      TES(349)%NYMD = 20090703
      TES(350)%NYMD = 20090703
      TES(351)%NYMD = 20090703
      TES(352)%NYMD = 20090703
      TES(353)%NYMD = 20090703
      TES(354)%NYMD = 20090703
      TES(355)%NYMD = 20090703
      TES(356)%NYMD = 20090703
      TES(357)%NYMD = 20090703
      TES(358)%NYMD = 20090703
      TES(359)%NYMD = 20090703
      TES(360)%NYMD = 20090703
      TES(361)%NYMD = 20090703
      TES(362)%NYMD = 20090703
      TES(363)%NYMD = 20090703
      TES(364)%NYMD = 20090703
      TES(365)%NYMD = 20090703
      TES(366)%NYMD = 20090703
      TES(367)%NYMD = 20090703
      TES(368)%NYMD = 20090703
      TES(369)%NYMD = 20090703
      TES(370)%NYMD = 20090703
      TES(371)%NYMD = 20090703
      TES(372)%NYMD = 20090703
      TES(373)%NYMD = 20090703
      TES(374)%NYMD = 20090703
      TES(375)%NYMD = 20090703
      TES(376)%NYMD = 20090703
      TES(377)%NYMD = 20090703
      TES(378)%NYMD = 20090703
      TES(379)%NYMD = 20090703
      TES(380)%NYMD = 20090703
      TES(381)%NYMD = 20090703
      TES(382)%NYMD = 20090703
      TES(383)%NYMD = 20090703
      TES(384)%NYMD = 20090703
      TES(385)%NYMD = 20090703
      TES(386)%NYMD = 20090703
      TES(387)%NYMD = 20090703
      TES(388)%NYMD = 20090703
      TES(389)%NYMD = 20090703
      TES(390)%NYMD = 20090703
      TES(391)%NYMD = 20090703
      TES(392)%NYMD = 20090703
      TES(393)%NYMD = 20090703
      TES(394)%NYMD = 20090703
      TES(395)%NYMD = 20090703
      TES(396)%NYMD = 20090703
      TES(397)%NYMD = 20090703
      TES(398)%NYMD = 20090703
      TES(399)%NYMD = 20090703
      TES(400)%NYMD = 20090703
      TES(401)%NYMD = 20090703
      TES(402)%NYMD = 20090703
      TES(403)%NYMD = 20090703
      TES(404)%NYMD = 20090703
      TES(405)%NYMD = 20090703
      TES(406)%NYMD = 20090703
      TES(407)%NYMD = 20090703
      TES(408)%NYMD = 20090703
      TES(409)%NYMD = 20090704
      TES(410)%NYMD = 20090704
      TES(411)%NYMD = 20090704
      TES(412)%NYMD = 20090704
      TES(413)%NYMD = 20090704
      TES(414)%NYMD = 20090704
      TES(415)%NYMD = 20090704
      TES(416)%NYMD = 20090704
      TES(417)%NYMD = 20090704
      TES(418)%NYMD = 20090704
      TES(419)%NYMD = 20090704
      TES(420)%NYMD = 20090704
      TES(421)%NYMD = 20090704
      TES(422)%NYMD = 20090704
      TES(423)%NYMD = 20090704
      TES(424)%NYMD = 20090704
      TES(425)%NYMD = 20090704
      TES(426)%NYMD = 20090704
      TES(427)%NYMD = 20090704
      TES(428)%NYMD = 20090704
      TES(429)%NYMD = 20090704
      TES(430)%NYMD = 20090704
      TES(431)%NYMD = 20090704
      TES(432)%NYMD = 20090704
      TES(433)%NYMD = 20090704
      TES(434)%NYMD = 20090704
      TES(435)%NYMD = 20090704
      TES(436)%NYMD = 20090704
      TES(437)%NYMD = 20090704
      TES(438)%NYMD = 20090704
      TES(439)%NYMD = 20090704
      TES(440)%NYMD = 20090704
      TES(441)%NYMD = 20090704
      TES(442)%NYMD = 20090704
      TES(443)%NYMD = 20090704
      TES(444)%NYMD = 20090704
      TES(445)%NYMD = 20090704
      TES(446)%NYMD = 20090704
      TES(447)%NYMD = 20090704
      TES(448)%NYMD = 20090704
      TES(449)%NYMD = 20090704
      TES(450)%NYMD = 20090704
      TES(451)%NYMD = 20090704
      TES(452)%NYMD = 20090704
      TES(453)%NYMD = 20090704
      TES(454)%NYMD = 20090704
      TES(455)%NYMD = 20090704
      TES(456)%NYMD = 20090704
      TES(457)%NYMD = 20090704
      TES(458)%NYMD = 20090704
      TES(459)%NYMD = 20090704
      TES(460)%NYMD = 20090704
      TES(461)%NYMD = 20090704
      TES(462)%NYMD = 20090704
      TES(463)%NYMD = 20090704
      TES(464)%NYMD = 20090704
      TES(465)%NYMD = 20090704
      TES(466)%NYMD = 20090704
      TES(467)%NYMD = 20090704
      TES(468)%NYMD = 20090704
      TES(469)%NYMD = 20090704
      TES(470)%NYMD = 20090704
      TES(471)%NYMD = 20090704
      TES(472)%NYMD = 20090704
      TES(473)%NYMD = 20090704
      TES(474)%NYMD = 20090704
      TES(475)%NYMD = 20090704
      TES(476)%NYMD = 20090704
      TES(477)%NYMD = 20090704
      TES(478)%NYMD = 20090704
      TES(479)%NYMD = 20090704
      TES(480)%NYMD = 20090704
      TES(481)%NYMD = 20090704
      TES(482)%NYMD = 20090704
      TES(483)%NYMD = 20090704
      TES(484)%NYMD = 20090704
      TES(485)%NYMD = 20090704
      TES(486)%NYMD = 20090704
      TES(487)%NYMD = 20090704
      TES(488)%NYMD = 20090704
      TES(489)%NYMD = 20090704
      TES(490)%NYMD = 20090704
      TES(491)%NYMD = 20090704
      TES(492)%NYMD = 20090704
      TES(493)%NYMD = 20090704
      TES(494)%NYMD = 20090704
      TES(495)%NYMD = 20090704
      TES(496)%NYMD = 20090704
      TES(497)%NYMD = 20090704
      TES(498)%NYMD = 20090704
      TES(499)%NYMD = 20090704
      TES(500)%NYMD = 20090704
      TES(501)%NYMD = 20090704
      TES(502)%NYMD = 20090704
      TES(503)%NYMD = 20090704
      TES(504)%NYMD = 20090704
      TES(505)%NYMD = 20090704
      TES(506)%NYMD = 20090704
      TES(507)%NYMD = 20090704
      TES(508)%NYMD = 20090704
      TES(509)%NYMD = 20090704
      TES(510)%NYMD = 20090704
      TES(511)%NYMD = 20090704
      TES(512)%NYMD = 20090704
      TES(513)%NYMD = 20090704
      TES(514)%NYMD = 20090704
      TES(515)%NYMD = 20090704
      TES(516)%NYMD = 20090704
      TES(517)%NYMD = 20090704
      TES(518)%NYMD = 20090705
      TES(519)%NYMD = 20090705
      TES(520)%NYMD = 20090705
      TES(521)%NYMD = 20090705
      TES(522)%NYMD = 20090705
      TES(523)%NYMD = 20090705
      TES(524)%NYMD = 20090705
      TES(525)%NYMD = 20090705
      TES(526)%NYMD = 20090705
      TES(527)%NYMD = 20090705
      TES(528)%NYMD = 20090705
      TES(529)%NYMD = 20090705
      TES(530)%NYMD = 20090705
      TES(531)%NYMD = 20090705
      TES(532)%NYMD = 20090705
      TES(533)%NYMD = 20090705
      TES(534)%NYMD = 20090705
      TES(535)%NYMD = 20090705
      TES(536)%NYMD = 20090705
      TES(537)%NYMD = 20090705
      TES(538)%NYMD = 20090705
      TES(539)%NYMD = 20090705
      TES(540)%NYMD = 20090705
      TES(541)%NYMD = 20090705
      TES(542)%NYMD = 20090705
      TES(543)%NYMD = 20090705
      TES(544)%NYMD = 20090705
      TES(545)%NYMD = 20090705
      TES(546)%NYMD = 20090705
      TES(547)%NYMD = 20090705
      TES(548)%NYMD = 20090705
      TES(549)%NYMD = 20090705
      TES(550)%NYMD = 20090705
      TES(551)%NYMD = 20090705
      TES(552)%NYMD = 20090705
      TES(553)%NYMD = 20090705
      TES(554)%NYMD = 20090705
      TES(555)%NYMD = 20090705
      TES(556)%NYMD = 20090705
      TES(557)%NYMD = 20090705
      TES(558)%NYMD = 20090705
      TES(559)%NYMD = 20090705
      TES(560)%NYMD = 20090705
      TES(561)%NYMD = 20090705
      TES(562)%NYMD = 20090705
      TES(563)%NYMD = 20090705
      TES(564)%NYMD = 20090705
      TES(565)%NYMD = 20090705
      TES(566)%NYMD = 20090705
      TES(567)%NYMD = 20090705
      TES(568)%NYMD = 20090705
      TES(569)%NYMD = 20090705
      TES(570)%NYMD = 20090705
      TES(571)%NYMD = 20090705
      TES(572)%NYMD = 20090705
      TES(573)%NYMD = 20090705
      TES(574)%NYMD = 20090705
      TES(575)%NYMD = 20090705
      TES(576)%NYMD = 20090705
      TES(577)%NYMD = 20090705
      TES(578)%NYMD = 20090705
      TES(579)%NYMD = 20090705
      TES(580)%NYMD = 20090705
      TES(581)%NYMD = 20090705
      TES(582)%NYMD = 20090705
      TES(583)%NYMD = 20090705
      TES(584)%NYMD = 20090705
      TES(585)%NYMD = 20090705
      TES(586)%NYMD = 20090705
      TES(587)%NYMD = 20090705
      TES(588)%NYMD = 20090705
      TES(589)%NYMD = 20090705
      TES(590)%NYMD = 20090705
      TES(591)%NYMD = 20090705
      TES(592)%NYMD = 20090705
      TES(593)%NYMD = 20090705
      TES(594)%NYMD = 20090705
      TES(595)%NYMD = 20090705
      TES(596)%NYMD = 20090705
      TES(597)%NYMD = 20090705
      TES(598)%NYMD = 20090705
      TES(599)%NYMD = 20090705
      TES(600)%NYMD = 20090705
      TES(601)%NYMD = 20090705
      TES(602)%NYMD = 20090705
      TES(603)%NYMD = 20090705
      TES(604)%NYMD = 20090705
      TES(605)%NYMD = 20090705
      TES(606)%NYMD = 20090705
      TES(607)%NYMD = 20090705
      TES(608)%NYMD = 20090705
      TES(609)%NYMD = 20090705
      TES(610)%NYMD = 20090705
      TES(611)%NYMD = 20090705
      TES(612)%NYMD = 20090705
      TES(613)%NYMD = 20090705
      TES(614)%NYMD = 20090705
      TES(615)%NYMD = 20090705
      TES(616)%NYMD = 20090705
      TES(617)%NYMD = 20090705
      TES(618)%NYMD = 20090705
      TES(619)%NYMD = 20090705
      TES(620)%NYMD = 20090705
      TES(621)%NYMD = 20090705
      TES(622)%NYMD = 20090705
      TES(623)%NYMD = 20090705
      TES(624)%NYMD = 20090705
      TES(625)%NYMD = 20090705
      TES(626)%NYMD = 20090705
      TES(627)%NYMD = 20090705
      TES(628)%NYMD = 20090705
      TES(629)%NYMD = 20090705
      TES(630)%NYMD = 20090705
      TES(631)%NYMD = 20090705
      TES(632)%NYMD = 20090705
      TES(633)%NYMD = 20090705
      TES(634)%NYMD = 20090705
      TES(635)%NYMD = 20090705
      TES(636)%NYMD = 20090705
      TES(637)%NYMD = 20090705
      TES(638)%NYMD = 20090705
      TES(639)%NYMD = 20090705
      TES(640)%NYMD = 20090705
      TES(641)%NYMD = 20090705
      TES(642)%NYMD = 20090705
      TES(643)%NYMD = 20090705
      TES(644)%NYMD = 20090705
      TES(645)%NYMD = 20090705
      TES(646)%NYMD = 20090705
      TES(647)%NYMD = 20090705
      TES(648)%NYMD = 20090705
      TES(649)%NYMD = 20090705
      TES(650)%NYMD = 20090705
      TES(651)%NYMD = 20090705
      TES(652)%NYMD = 20090705
      TES(653)%NYMD = 20090705
      TES(654)%NYMD = 20090705
      TES(655)%NYMD = 20090705
      TES(656)%NYMD = 20090705
      TES(657)%NYMD = 20090705
      TES(658)%NYMD = 20090705
      TES(659)%NYMD = 20090705
      TES(660)%NYMD = 20090705
      TES(661)%NYMD = 20090705
      TES(662)%NYMD = 20090705
      TES(663)%NYMD = 20090705
      TES(664)%NYMD = 20090705
      TES(665)%NYMD = 20090705
      TES(666)%NYMD = 20090705
      TES(667)%NYMD = 20090705
      TES(668)%NYMD = 20090705
      TES(669)%NYMD = 20090705
      TES(670)%NYMD = 20090705
      TES(671)%NYMD = 20090705
      TES(672)%NYMD = 20090705
      TES(673)%NYMD = 20090705
      TES(674)%NYMD = 20090705
      TES(675)%NYMD = 20090705
      TES(676)%NYMD = 20090705
      TES(677)%NYMD = 20090706
      TES(678)%NYMD = 20090706
      TES(679)%NYMD = 20090706
      TES(680)%NYMD = 20090706
      TES(681)%NYMD = 20090706
      TES(682)%NYMD = 20090706
      TES(683)%NYMD = 20090706
      TES(684)%NYMD = 20090706
      TES(685)%NYMD = 20090706
      TES(686)%NYMD = 20090706
      TES(687)%NYMD = 20090706
      TES(688)%NYMD = 20090706
      TES(689)%NYMD = 20090706
      TES(690)%NYMD = 20090706
      TES(691)%NYMD = 20090706
      TES(692)%NYMD = 20090706
      TES(693)%NYMD = 20090706
      TES(694)%NYMD = 20090706
      TES(695)%NYMD = 20090706
      TES(696)%NYMD = 20090706
      TES(697)%NYMD = 20090706
      TES(698)%NYMD = 20090706
      TES(699)%NYMD = 20090706
      TES(700)%NYMD = 20090706
      TES(701)%NYMD = 20090706
      TES(702)%NYMD = 20090706
      TES(703)%NYMD = 20090706
      TES(704)%NYMD = 20090706
      TES(705)%NYMD = 20090706
      TES(706)%NYMD = 20090706
      TES(707)%NYMD = 20090706
      TES(708)%NYMD = 20090706
      TES(709)%NYMD = 20090706
      TES(710)%NYMD = 20090706
      TES(711)%NYMD = 20090706
      TES(712)%NYMD = 20090706
      TES(713)%NYMD = 20090706
      TES(714)%NYMD = 20090706
      TES(715)%NYMD = 20090706
      TES(716)%NYMD = 20090706
      TES(717)%NYMD = 20090706
      TES(718)%NYMD = 20090706
      TES(719)%NYMD = 20090706
      TES(720)%NYMD = 20090706
      TES(721)%NYMD = 20090706
      TES(722)%NYMD = 20090706
      TES(723)%NYMD = 20090706
      TES(724)%NYMD = 20090706
      TES(725)%NYMD = 20090706
      TES(726)%NYMD = 20090706
      TES(727)%NYMD = 20090706
      TES(728)%NYMD = 20090706
      TES(729)%NYMD = 20090706
      TES(730)%NYMD = 20090706
      TES(731)%NYMD = 20090706
      TES(732)%NYMD = 20090706
      TES(733)%NYMD = 20090706
      TES(734)%NYMD = 20090706
      TES(735)%NYMD = 20090706
      TES(736)%NYMD = 20090706
      TES(737)%NYMD = 20090706
      TES(738)%NYMD = 20090706
      TES(739)%NYMD = 20090706
      TES(740)%NYMD = 20090706
      TES(741)%NYMD = 20090706
      TES(742)%NYMD = 20090706
      TES(743)%NYMD = 20090706
      TES(744)%NYMD = 20090706
      TES(745)%NYMD = 20090706
      TES(746)%NYMD = 20090706
      TES(747)%NYMD = 20090706
      TES(748)%NYMD = 20090706
      TES(749)%NYMD = 20090706
      TES(750)%NYMD = 20090706
      TES(751)%NYMD = 20090706
      TES(752)%NYMD = 20090706
      TES(753)%NYMD = 20090706
      TES(754)%NYMD = 20090706
      TES(755)%NYMD = 20090706
      TES(756)%NYMD = 20090706
      TES(757)%NYMD = 20090706
      TES(758)%NYMD = 20090706
      TES(759)%NYMD = 20090706
      TES(760)%NYMD = 20090706
      TES(761)%NYMD = 20090706
      TES(762)%NYMD = 20090706
      TES(763)%NYMD = 20090706
      TES(764)%NYMD = 20090706
      TES(765)%NYMD = 20090706
      TES(766)%NYMD = 20090706
      TES(767)%NYMD = 20090706
      TES(768)%NYMD = 20090706
      TES(769)%NYMD = 20090706
      TES(770)%NYMD = 20090706
      TES(771)%NYMD = 20090706
      TES(772)%NYMD = 20090706
      TES(773)%NYMD = 20090706
      TES(774)%NYMD = 20090706
      TES(775)%NYMD = 20090706
      TES(776)%NYMD = 20090706
      TES(777)%NYMD = 20090706
      TES(778)%NYMD = 20090706
      TES(779)%NYMD = 20090706
      TES(780)%NYMD = 20090706
      TES(781)%NYMD = 20090706
      TES(782)%NYMD = 20090706
      TES(783)%NYMD = 20090706
      TES(784)%NYMD = 20090706
      TES(785)%NYMD = 20090706
      TES(786)%NYMD = 20090706
      TES(787)%NYMD = 20090706
      TES(788)%NYMD = 20090706
      TES(789)%NYMD = 20090706
      TES(790)%NYMD = 20090706
      TES(791)%NYMD = 20090706
      TES(792)%NYMD = 20090706
      TES(793)%NYMD = 20090706
      TES(794)%NYMD = 20090706
      TES(795)%NYMD = 20090706
      TES(796)%NYMD = 20090706
      TES(797)%NYMD = 20090706
      TES(798)%NYMD = 20090706
      TES(799)%NYMD = 20090706
      TES(800)%NYMD = 20090706
      TES(801)%NYMD = 20090706
      TES(802)%NYMD = 20090706
      TES(803)%NYMD = 20090706
      TES(804)%NYMD = 20090709
      TES(805)%NYMD = 20090709
      TES(806)%NYMD = 20090709
      TES(807)%NYMD = 20090709
      TES(808)%NYMD = 20090709
      TES(809)%NYMD = 20090709
      TES(810)%NYMD = 20090709
      TES(811)%NYMD = 20090709
      TES(812)%NYMD = 20090709
      TES(813)%NYMD = 20090709
      TES(814)%NYMD = 20090709
      TES(815)%NYMD = 20090709
      TES(816)%NYMD = 20090709
      TES(817)%NYMD = 20090709
      TES(818)%NYMD = 20090709
      TES(819)%NYMD = 20090709
      TES(820)%NYMD = 20090709
      TES(821)%NYMD = 20090709
      TES(822)%NYMD = 20090709
      TES(823)%NYMD = 20090709
      TES(824)%NYMD = 20090709
      TES(825)%NYMD = 20090709
      TES(826)%NYMD = 20090709
      TES(827)%NYMD = 20090709
      TES(828)%NYMD = 20090709
      TES(829)%NYMD = 20090709
      TES(830)%NYMD = 20090709
      TES(831)%NYMD = 20090709
      TES(832)%NYMD = 20090709
      TES(833)%NYMD = 20090709
      TES(834)%NYMD = 20090709
      TES(835)%NYMD = 20090709
      TES(836)%NYMD = 20090709
      TES(837)%NYMD = 20090709
      TES(838)%NYMD = 20090709
      TES(839)%NYMD = 20090709
      TES(840)%NYMD = 20090709
      TES(841)%NYMD = 20090709
      TES(842)%NYMD = 20090709
      TES(843)%NYMD = 20090709
      TES(844)%NYMD = 20090709
      TES(845)%NYMD = 20090709
      TES(846)%NYMD = 20090709
      TES(847)%NYMD = 20090709
      TES(848)%NYMD = 20090709
      TES(849)%NYMD = 20090709
      TES(850)%NYMD = 20090709
      TES(851)%NYMD = 20090709
      TES(852)%NYMD = 20090709
      TES(853)%NYMD = 20090709
      TES(854)%NYMD = 20090709
      TES(855)%NYMD = 20090709
      TES(856)%NYMD = 20090709
      TES(857)%NYMD = 20090709
      TES(858)%NYMD = 20090709
      TES(859)%NYMD = 20090709
      TES(860)%NYMD = 20090709
      TES(861)%NYMD = 20090709
      TES(862)%NYMD = 20090709
      TES(863)%NYMD = 20090709
      TES(864)%NYMD = 20090709
      TES(865)%NYMD = 20090709
      TES(866)%NYMD = 20090709
      TES(867)%NYMD = 20090709
      TES(868)%NYMD = 20090709
      TES(869)%NYMD = 20090709
      TES(870)%NYMD = 20090709
      TES(871)%NYMD = 20090709
      TES(872)%NYMD = 20090709
      TES(873)%NYMD = 20090709
      TES(874)%NYMD = 20090709
      TES(875)%NYMD = 20090709
      TES(876)%NYMD = 20090709
      TES(877)%NYMD = 20090709
      TES(878)%NYMD = 20090709
      TES(879)%NYMD = 20090709
      TES(880)%NYMD = 20090709
      TES(881)%NYMD = 20090709
      TES(882)%NYMD = 20090709
      TES(883)%NYMD = 20090709
      TES(884)%NYMD = 20090709
      TES(885)%NYMD = 20090709
      TES(886)%NYMD = 20090709
      TES(887)%NYMD = 20090709
      TES(888)%NYMD = 20090709
      TES(889)%NYMD = 20090709
      TES(890)%NYMD = 20090709
      TES(891)%NYMD = 20090709
      TES(892)%NYMD = 20090709
      TES(893)%NYMD = 20090709
      TES(894)%NYMD = 20090709
      TES(895)%NYMD = 20090709
      TES(896)%NYMD = 20090709
      TES(897)%NYMD = 20090709
      TES(898)%NYMD = 20090709
      TES(899)%NYMD = 20090709
      TES(900)%NYMD = 20090709
      TES(901)%NYMD = 20090709
      TES(902)%NYMD = 20090709
      TES(903)%NYMD = 20090709
      TES(904)%NYMD = 20090709
      TES(905)%NYMD = 20090709
      TES(906)%NYMD = 20090709
      TES(907)%NYMD = 20090709
      TES(908)%NYMD = 20090709
      TES(909)%NYMD = 20090709
      TES(910)%NYMD = 20090709
      TES(911)%NYMD = 20090709
      TES(912)%NYMD = 20090709
      TES(913)%NYMD = 20090709
      TES(914)%NYMD = 20090709
      TES(915)%NYMD = 20090709
      TES(916)%NYMD = 20090709
      TES(917)%NYMD = 20090709
      TES(918)%NYMD = 20090709
      TES(919)%NYMD = 20090709
      TES(920)%NYMD = 20090709
      TES(921)%NYMD = 20090709
      TES(922)%NYMD = 20090709
      TES(923)%NYMD = 20090709
      TES(924)%NYMD = 20090709
      TES(925)%NYMD = 20090709
      TES(926)%NYMD = 20090709
      TES(927)%NYMD = 20090709
      TES(928)%NYMD = 20090709
      TES(929)%NYMD = 20090709
      TES(930)%NYMD = 20090709
      TES(931)%NYMD = 20090709
      TES(932)%NYMD = 20090709
      TES(933)%NYMD = 20090709
      TES(934)%NYMD = 20090709
      TES(935)%NYMD = 20090709
      TES(936)%NYMD = 20090709
      TES(937)%NYMD = 20090709
      TES(938)%NYMD = 20090709
      TES(939)%NYMD = 20090709
      TES(940)%NYMD = 20090709
      TES(941)%NYMD = 20090709
      TES(942)%NYMD = 20090709
      TES(943)%NYMD = 20090709
      TES(944)%NYMD = 20090709
      TES(945)%NYMD = 20090709
      TES(946)%NYMD = 20090709
      TES(947)%NYMD = 20090709
      TES(948)%NYMD = 20090709
      TES(949)%NYMD = 20090709
      TES(950)%NYMD = 20090709
      TES(951)%NYMD = 20090709
      TES(952)%NYMD = 20090710
      TES(953)%NYMD = 20090710
      TES(954)%NYMD = 20090710
      TES(955)%NYMD = 20090710
      TES(956)%NYMD = 20090710
      TES(957)%NYMD = 20090710
      TES(958)%NYMD = 20090710
      TES(959)%NYMD = 20090710
      TES(960)%NYMD = 20090710
      TES(961)%NYMD = 20090710
      TES(962)%NYMD = 20090710
      TES(963)%NYMD = 20090710
      TES(964)%NYMD = 20090710
      TES(965)%NYMD = 20090710
      TES(966)%NYMD = 20090710
      TES(967)%NYMD = 20090710
      TES(968)%NYMD = 20090710
      TES(969)%NYMD = 20090710
      TES(970)%NYMD = 20090710
      TES(971)%NYMD = 20090710
      TES(972)%NYMD = 20090710
      TES(973)%NYMD = 20090710
      TES(974)%NYMD = 20090710
      TES(975)%NYMD = 20090710
      TES(976)%NYMD = 20090710
      TES(977)%NYMD = 20090710
      TES(978)%NYMD = 20090710
      TES(979)%NYMD = 20090710
      TES(980)%NYMD = 20090710
      TES(981)%NYMD = 20090710
      TES(982)%NYMD = 20090710
      TES(983)%NYMD = 20090710
      TES(984)%NYMD = 20090710
      TES(985)%NYMD = 20090710
      TES(986)%NYMD = 20090710
      TES(987)%NYMD = 20090710
      TES(988)%NYMD = 20090710
      TES(989)%NYMD = 20090710
      TES(990)%NYMD = 20090710
      TES(991)%NYMD = 20090710
      TES(992)%NYMD = 20090710
      TES(993)%NYMD = 20090710
      TES(994)%NYMD = 20090710
      TES(995)%NYMD = 20090710
      TES(996)%NYMD = 20090710
      TES(997)%NYMD = 20090710
      TES(998)%NYMD = 20090710
      TES(999)%NYMD = 20090710
      TES(1000)%NYMD = 20090710
      TES(1001)%NYMD = 20090710
      TES(1002)%NYMD = 20090710
      TES(1003)%NYMD = 20090710
      TES(1004)%NYMD = 20090710
      TES(1005)%NYMD = 20090710
      TES(1006)%NYMD = 20090710
      TES(1007)%NYMD = 20090710
      TES(1008)%NYMD = 20090710
      TES(1009)%NYMD = 20090710
      TES(1010)%NYMD = 20090710
      TES(1011)%NYMD = 20090710
      TES(1012)%NYMD = 20090710
      TES(1013)%NYMD = 20090710
      TES(1014)%NYMD = 20090710
      TES(1015)%NYMD = 20090710
      TES(1016)%NYMD = 20090710
      TES(1017)%NYMD = 20090710
      TES(1018)%NYMD = 20090710
      TES(1019)%NYMD = 20090710
      TES(1020)%NYMD = 20090710
      TES(1021)%NYMD = 20090710
      TES(1022)%NYMD = 20090710
      TES(1023)%NYMD = 20090710
      TES(1024)%NYMD = 20090710
      TES(1025)%NYMD = 20090710
      TES(1026)%NYMD = 20090710
      TES(1027)%NYMD = 20090710
      TES(1028)%NYMD = 20090710
      TES(1029)%NYMD = 20090710
      TES(1030)%NYMD = 20090710
      TES(1031)%NYMD = 20090710
      TES(1032)%NYMD = 20090710
      TES(1033)%NYMD = 20090710
      TES(1034)%NYMD = 20090710
      TES(1035)%NYMD = 20090710
      TES(1036)%NYMD = 20090710
      TES(1037)%NYMD = 20090710
      TES(1038)%NYMD = 20090710
      TES(1039)%NYMD = 20090710
      TES(1040)%NYMD = 20090710
      TES(1041)%NYMD = 20090710
      TES(1042)%NYMD = 20090710
      TES(1043)%NYMD = 20090710
      TES(1044)%NYMD = 20090710
      TES(1045)%NYMD = 20090710
      TES(1046)%NYMD = 20090710
      TES(1047)%NYMD = 20090710
      TES(1048)%NYMD = 20090710
      TES(1049)%NYMD = 20090710
      TES(1050)%NYMD = 20090710
      TES(1051)%NYMD = 20090710
      TES(1052)%NYMD = 20090710
      TES(1053)%NYMD = 20090710
      TES(1054)%NYMD = 20090710
      TES(1055)%NYMD = 20090710
      TES(1056)%NYMD = 20090710
      TES(1057)%NYMD = 20090710
      TES(1058)%NYMD = 20090710
      TES(1059)%NYMD = 20090710
      TES(1060)%NYMD = 20090710
      TES(1061)%NYMD = 20090710
      TES(1062)%NYMD = 20090710
      TES(1063)%NYMD = 20090710
      TES(1064)%NYMD = 20090710
      TES(1065)%NYMD = 20090710
      TES(1066)%NYMD = 20090710
      TES(1067)%NYMD = 20090710
      TES(1068)%NYMD = 20090710
      TES(1069)%NYMD = 20090710
      TES(1070)%NYMD = 20090710
      TES(1071)%NYMD = 20090710
      TES(1072)%NYMD = 20090710
      TES(1073)%NYMD = 20090710
      TES(1074)%NYMD = 20090710
      TES(1075)%NYMD = 20090710
      TES(1076)%NYMD = 20090710
      TES(1077)%NYMD = 20090710
      TES(1078)%NYMD = 20090710
      TES(1079)%NYMD = 20090710
      TES(1080)%NYMD = 20090710
      TES(1081)%NYMD = 20090710
      TES(1082)%NYMD = 20090710
      TES(1083)%NYMD = 20090710
      TES(1084)%NYMD = 20090710
      TES(1085)%NYMD = 20090710
      TES(1086)%NYMD = 20090710
      TES(1087)%NYMD = 20090710
      TES(1088)%NYMD = 20090710
      TES(1089)%NYMD = 20090710
      TES(1090)%NYMD = 20090710
      TES(1091)%NYMD = 20090711
      TES(1092)%NYMD = 20090711
      TES(1093)%NYMD = 20090711
      TES(1094)%NYMD = 20090711
      TES(1095)%NYMD = 20090711
      TES(1096)%NYMD = 20090711
      TES(1097)%NYMD = 20090711
      TES(1098)%NYMD = 20090711
      TES(1099)%NYMD = 20090711
      TES(1100)%NYMD = 20090711
      TES(1101)%NYMD = 20090711
      TES(1102)%NYMD = 20090711
      TES(1103)%NYMD = 20090711
      TES(1104)%NYMD = 20090711
      TES(1105)%NYMD = 20090711
      TES(1106)%NYMD = 20090711
      TES(1107)%NYMD = 20090711
      TES(1108)%NYMD = 20090711
      TES(1109)%NYMD = 20090711
      TES(1110)%NYMD = 20090711
      TES(1111)%NYMD = 20090711
      TES(1112)%NYMD = 20090711
      TES(1113)%NYMD = 20090711
      TES(1114)%NYMD = 20090711
      TES(1115)%NYMD = 20090711
      TES(1116)%NYMD = 20090711
      TES(1117)%NYMD = 20090711
      TES(1118)%NYMD = 20090711
      TES(1119)%NYMD = 20090711
      TES(1120)%NYMD = 20090711
      TES(1121)%NYMD = 20090711
      TES(1122)%NYMD = 20090711
      TES(1123)%NYMD = 20090711
      TES(1124)%NYMD = 20090711
      TES(1125)%NYMD = 20090711
      TES(1126)%NYMD = 20090711
      TES(1127)%NYMD = 20090711
      TES(1128)%NYMD = 20090711
      TES(1129)%NYMD = 20090711
      TES(1130)%NYMD = 20090711
      TES(1131)%NYMD = 20090711
      TES(1132)%NYMD = 20090711
      TES(1133)%NYMD = 20090711
      TES(1134)%NYMD = 20090711
      TES(1135)%NYMD = 20090711
      TES(1136)%NYMD = 20090711
      TES(1137)%NYMD = 20090711
      TES(1138)%NYMD = 20090711
      TES(1139)%NYMD = 20090711
      TES(1140)%NYMD = 20090711
      TES(1141)%NYMD = 20090711
      TES(1142)%NYMD = 20090711
      TES(1143)%NYMD = 20090711
      TES(1144)%NYMD = 20090711
      TES(1145)%NYMD = 20090711
      TES(1146)%NYMD = 20090711
      TES(1147)%NYMD = 20090711
      TES(1148)%NYMD = 20090711
      TES(1149)%NYMD = 20090711
      TES(1150)%NYMD = 20090711
      TES(1151)%NYMD = 20090711
      TES(1152)%NYMD = 20090711
      TES(1153)%NYMD = 20090711
      TES(1154)%NYMD = 20090711
      TES(1155)%NYMD = 20090711
      TES(1156)%NYMD = 20090711
      TES(1157)%NYMD = 20090711
      TES(1158)%NYMD = 20090711
      TES(1159)%NYMD = 20090711
      TES(1160)%NYMD = 20090711
      TES(1161)%NYMD = 20090711
      TES(1162)%NYMD = 20090711
      TES(1163)%NYMD = 20090711
      TES(1164)%NYMD = 20090711
      TES(1165)%NYMD = 20090711
      TES(1166)%NYMD = 20090711
      TES(1167)%NYMD = 20090711
      TES(1168)%NYMD = 20090711
      TES(1169)%NYMD = 20090711
      TES(1170)%NYMD = 20090711
      TES(1171)%NYMD = 20090711
      TES(1172)%NYMD = 20090711
      TES(1173)%NYMD = 20090711
      TES(1174)%NYMD = 20090711
      TES(1175)%NYMD = 20090711
      TES(1176)%NYMD = 20090711
      TES(1177)%NYMD = 20090711
      TES(1178)%NYMD = 20090711
      TES(1179)%NYMD = 20090711
      TES(1180)%NYMD = 20090711
      TES(1181)%NYMD = 20090711
      TES(1182)%NYMD = 20090711
      TES(1183)%NYMD = 20090711
      TES(1184)%NYMD = 20090711
      TES(1185)%NYMD = 20090711
      TES(1186)%NYMD = 20090711
      TES(1187)%NYMD = 20090711
      TES(1188)%NYMD = 20090711
      TES(1189)%NYMD = 20090711
      TES(1190)%NYMD = 20090711
      TES(1191)%NYMD = 20090711
      TES(1192)%NYMD = 20090711
      TES(1193)%NYMD = 20090711
      TES(1194)%NYMD = 20090711
      TES(1195)%NYMD = 20090711
      TES(1196)%NYMD = 20090711
      TES(1197)%NYMD = 20090711
      TES(1198)%NYMD = 20090711
      TES(1199)%NYMD = 20090711
      TES(1200)%NYMD = 20090711
      TES(1201)%NYMD = 20090711
      TES(1202)%NYMD = 20090711
      TES(1203)%NYMD = 20090711
      TES(1204)%NYMD = 20090711
      TES(1205)%NYMD = 20090711
      TES(1206)%NYMD = 20090711
      TES(1207)%NYMD = 20090711
      TES(1208)%NYMD = 20090711
      TES(1209)%NYMD = 20090711
      TES(1210)%NYMD = 20090711
      TES(1211)%NYMD = 20090711
      TES(1212)%NYMD = 20090711
      TES(1213)%NYMD = 20090711
      TES(1214)%NYMD = 20090711
      TES(1215)%NYMD = 20090711
      TES(1216)%NYMD = 20090711
      TES(1217)%NYMD = 20090711
      TES(1218)%NYMD = 20090711
      TES(1219)%NYMD = 20090711
      TES(1220)%NYMD = 20090711
      TES(1221)%NYMD = 20090711
      TES(1222)%NYMD = 20090711
      TES(1223)%NYMD = 20090711
      TES(1224)%NYMD = 20090711
      TES(1225)%NYMD = 20090711
      TES(1226)%NYMD = 20090711
      TES(1227)%NYMD = 20090711
      TES(1228)%NYMD = 20090711
      TES(1229)%NYMD = 20090711
      TES(1230)%NYMD = 20090711
      TES(1231)%NYMD = 20090711
      TES(1232)%NYMD = 20090711
      TES(1233)%NYMD = 20090711
      TES(1234)%NYMD = 20090711
      TES(1235)%NYMD = 20090711
      TES(1236)%NYMD = 20090711
      TES(1237)%NYMD = 20090711
      TES(1238)%NYMD = 20090711
      TES(1239)%NYMD = 20090711
      TES(1240)%NYMD = 20090711
      TES(1241)%NYMD = 20090711
      TES(1242)%NYMD = 20090711
      TES(1243)%NYMD = 20090711
      TES(1244)%NYMD = 20090711
      TES(1245)%NYMD = 20090711
      TES(1246)%NYMD = 20090711
      TES(1247)%NYMD = 20090712
      TES(1248)%NYMD = 20090712
      TES(1249)%NYMD = 20090712
      TES(1250)%NYMD = 20090712
      TES(1251)%NYMD = 20090712
      TES(1252)%NYMD = 20090712
      TES(1253)%NYMD = 20090712
      TES(1254)%NYMD = 20090712
      TES(1255)%NYMD = 20090712
      TES(1256)%NYMD = 20090712
      TES(1257)%NYMD = 20090712
      TES(1258)%NYMD = 20090712
      TES(1259)%NYMD = 20090712
      TES(1260)%NYMD = 20090712
      TES(1261)%NYMD = 20090712
      TES(1262)%NYMD = 20090712
      TES(1263)%NYMD = 20090712
      TES(1264)%NYMD = 20090712
      TES(1265)%NYMD = 20090712
      TES(1266)%NYMD = 20090712
      TES(1267)%NYMD = 20090712
      TES(1268)%NYMD = 20090712
      TES(1269)%NYMD = 20090712
      TES(1270)%NYMD = 20090712
      TES(1271)%NYMD = 20090712
      TES(1272)%NYMD = 20090712
      TES(1273)%NYMD = 20090712
      TES(1274)%NYMD = 20090712
      TES(1275)%NYMD = 20090712
      TES(1276)%NYMD = 20090712
      TES(1277)%NYMD = 20090712
      TES(1278)%NYMD = 20090712
      TES(1279)%NYMD = 20090712
      TES(1280)%NYMD = 20090712
      TES(1281)%NYMD = 20090712
      TES(1282)%NYMD = 20090712
      TES(1283)%NYMD = 20090712
      TES(1284)%NYMD = 20090712
      TES(1285)%NYMD = 20090712
      TES(1286)%NYMD = 20090712
      TES(1287)%NYMD = 20090712
      TES(1288)%NYMD = 20090712
      TES(1289)%NYMD = 20090712
      TES(1290)%NYMD = 20090712
      TES(1291)%NYMD = 20090712
      TES(1292)%NYMD = 20090712
      TES(1293)%NYMD = 20090712
      TES(1294)%NYMD = 20090712
      TES(1295)%NYMD = 20090712
      TES(1296)%NYMD = 20090712
      TES(1297)%NYMD = 20090712
      TES(1298)%NYMD = 20090712
      TES(1299)%NYMD = 20090712
      TES(1300)%NYMD = 20090712
      TES(1301)%NYMD = 20090712
      TES(1302)%NYMD = 20090712
      TES(1303)%NYMD = 20090712
      TES(1304)%NYMD = 20090712
      TES(1305)%NYMD = 20090712
      TES(1306)%NYMD = 20090712
      TES(1307)%NYMD = 20090712
      TES(1308)%NYMD = 20090712
      TES(1309)%NYMD = 20090712
      TES(1310)%NYMD = 20090712
      TES(1311)%NYMD = 20090712
      TES(1312)%NYMD = 20090712
      TES(1313)%NYMD = 20090712
      TES(1314)%NYMD = 20090712
      TES(1315)%NYMD = 20090712
      TES(1316)%NYMD = 20090712
      TES(1317)%NYMD = 20090712
      TES(1318)%NYMD = 20090712
      TES(1319)%NYMD = 20090712
      TES(1320)%NYMD = 20090712
      TES(1321)%NYMD = 20090712
      TES(1322)%NYMD = 20090712
      TES(1323)%NYMD = 20090712
      TES(1324)%NYMD = 20090712
      TES(1325)%NYMD = 20090712
      TES(1326)%NYMD = 20090712
      TES(1327)%NYMD = 20090712
      TES(1328)%NYMD = 20090712
      TES(1329)%NYMD = 20090712
      TES(1330)%NYMD = 20090712
      TES(1331)%NYMD = 20090712
      TES(1332)%NYMD = 20090712
      TES(1333)%NYMD = 20090712
      TES(1334)%NYMD = 20090712
      TES(1335)%NYMD = 20090712
      TES(1336)%NYMD = 20090712
      TES(1337)%NYMD = 20090712
      TES(1338)%NYMD = 20090712
      TES(1339)%NYMD = 20090712
      TES(1340)%NYMD = 20090712
      TES(1341)%NYMD = 20090712
      TES(1342)%NYMD = 20090712
      TES(1343)%NYMD = 20090712
      TES(1344)%NYMD = 20090712
      TES(1345)%NYMD = 20090712
      TES(1346)%NYMD = 20090712
      TES(1347)%NYMD = 20090712
      TES(1348)%NYMD = 20090712
      TES(1349)%NYMD = 20090712
      TES(1350)%NYMD = 20090712
      TES(1351)%NYMD = 20090712
      TES(1352)%NYMD = 20090712
      TES(1353)%NYMD = 20090712
      TES(1354)%NYMD = 20090712
      TES(1355)%NYMD = 20090712
      TES(1356)%NYMD = 20090712
      TES(1357)%NYMD = 20090712
      TES(1358)%NYMD = 20090712
      TES(1359)%NYMD = 20090712
      TES(1360)%NYMD = 20090712
      TES(1361)%NYMD = 20090712
      TES(1362)%NYMD = 20090712
      TES(1363)%NYMD = 20090712
      TES(1364)%NYMD = 20090712
      TES(1365)%NYMD = 20090712
      TES(1366)%NYMD = 20090712
      TES(1367)%NYMD = 20090712
      TES(1368)%NYMD = 20090712
      TES(1369)%NYMD = 20090712
      TES(1370)%NYMD = 20090712
      TES(1371)%NYMD = 20090712
      TES(1372)%NYMD = 20090712
      TES(1373)%NYMD = 20090713
      TES(1374)%NYMD = 20090713
      TES(1375)%NYMD = 20090713
      TES(1376)%NYMD = 20090713
      TES(1377)%NYMD = 20090713
      TES(1378)%NYMD = 20090713
      TES(1379)%NYMD = 20090713
      TES(1380)%NYMD = 20090713
      TES(1381)%NYMD = 20090713
      TES(1382)%NYMD = 20090713
      TES(1383)%NYMD = 20090713
      TES(1384)%NYMD = 20090713
      TES(1385)%NYMD = 20090713
      TES(1386)%NYMD = 20090713
      TES(1387)%NYMD = 20090713
      TES(1388)%NYMD = 20090713
      TES(1389)%NYMD = 20090713
      TES(1390)%NYMD = 20090713
      TES(1391)%NYMD = 20090713
      TES(1392)%NYMD = 20090713
      TES(1393)%NYMD = 20090713
      TES(1394)%NYMD = 20090713
      TES(1395)%NYMD = 20090713
      TES(1396)%NYMD = 20090713
      TES(1397)%NYMD = 20090713
      TES(1398)%NYMD = 20090713
      TES(1399)%NYMD = 20090713
      TES(1400)%NYMD = 20090713
      TES(1401)%NYMD = 20090713
      TES(1402)%NYMD = 20090713
      TES(1403)%NYMD = 20090713
      TES(1404)%NYMD = 20090713
      TES(1405)%NYMD = 20090713
      TES(1406)%NYMD = 20090713
      TES(1407)%NYMD = 20090713
      TES(1408)%NYMD = 20090713
      TES(1409)%NYMD = 20090713
      TES(1410)%NYMD = 20090713
      TES(1411)%NYMD = 20090713
      TES(1412)%NYMD = 20090713
      TES(1413)%NYMD = 20090713
      TES(1414)%NYMD = 20090713
      TES(1415)%NYMD = 20090713
      TES(1416)%NYMD = 20090713
      TES(1417)%NYMD = 20090713
      TES(1418)%NYMD = 20090713
      TES(1419)%NYMD = 20090713
      TES(1420)%NYMD = 20090713
      TES(1421)%NYMD = 20090713
      TES(1422)%NYMD = 20090713
      TES(1423)%NYMD = 20090713
      TES(1424)%NYMD = 20090713
      TES(1425)%NYMD = 20090713
      TES(1426)%NYMD = 20090713
      TES(1427)%NYMD = 20090713
      TES(1428)%NYMD = 20090713
      TES(1429)%NYMD = 20090713
      TES(1430)%NYMD = 20090713
      TES(1431)%NYMD = 20090713
      TES(1432)%NYMD = 20090713
      TES(1433)%NYMD = 20090713
      TES(1434)%NYMD = 20090713
      TES(1435)%NYMD = 20090713
      TES(1436)%NYMD = 20090713
      TES(1437)%NYMD = 20090713
      TES(1438)%NYMD = 20090713
      TES(1439)%NYMD = 20090713
      TES(1440)%NYMD = 20090713
      TES(1441)%NYMD = 20090713
      TES(1442)%NYMD = 20090713
      TES(1443)%NYMD = 20090713
      TES(1444)%NYMD = 20090713
      TES(1445)%NYMD = 20090713
      TES(1446)%NYMD = 20090713
      TES(1447)%NYMD = 20090713
      TES(1448)%NYMD = 20090713
      TES(1449)%NYMD = 20090713
      TES(1450)%NYMD = 20090713
      TES(1451)%NYMD = 20090713
      TES(1452)%NYMD = 20090713
      TES(1453)%NYMD = 20090713
      TES(1454)%NYMD = 20090713
      TES(1455)%NYMD = 20090713
      TES(1456)%NYMD = 20090713
      TES(1457)%NYMD = 20090713
      TES(1458)%NYMD = 20090713
      TES(1459)%NYMD = 20090713
      TES(1460)%NYMD = 20090713
      TES(1461)%NYMD = 20090713
      TES(1462)%NYMD = 20090713
      TES(1463)%NYMD = 20090713
      TES(1464)%NYMD = 20090713
      TES(1465)%NYMD = 20090713
      TES(1466)%NYMD = 20090713
      TES(1467)%NYMD = 20090713
      TES(1468)%NYMD = 20090713
      TES(1469)%NYMD = 20090713
      TES(1470)%NYMD = 20090713
      TES(1471)%NYMD = 20090713
      TES(1472)%NYMD = 20090713
      TES(1473)%NYMD = 20090713
      TES(1474)%NYMD = 20090713
      TES(1475)%NYMD = 20090713
      TES(1476)%NYMD = 20090713
      TES(1477)%NYMD = 20090713
      TES(1478)%NYMD = 20090713
      TES(1479)%NYMD = 20090713
      TES(1480)%NYMD = 20090713
      TES(1481)%NYMD = 20090713
      TES(1482)%NYMD = 20090713
      TES(1483)%NYMD = 20090713
      TES(1484)%NYMD = 20090713
      TES(1485)%NYMD = 20090713
      TES(1486)%NYMD = 20090713
      TES(1487)%NYMD = 20090713
      TES(1488)%NYMD = 20090713
      TES(1489)%NYMD = 20090713
      TES(1490)%NYMD = 20090713
      TES(1491)%NYMD = 20090713
      TES(1492)%NYMD = 20090713
      TES(1493)%NYMD = 20090713
      TES(1494)%NYMD = 20090713
      TES(1495)%NYMD = 20090713
      TES(1496)%NYMD = 20090713
      TES(1497)%NYMD = 20090713
      TES(1498)%NYMD = 20090713
      TES(1499)%NYMD = 20090713
      TES(1500)%NYMD = 20090713
      TES(1501)%NYMD = 20090713
      TES(1502)%NYMD = 20090713
      TES(1503)%NYMD = 20090713
      TES(1504)%NYMD = 20090714
      TES(1505)%NYMD = 20090714
      TES(1506)%NYMD = 20090714
      TES(1507)%NYMD = 20090714
      TES(1508)%NYMD = 20090714
      TES(1509)%NYMD = 20090714
      TES(1510)%NYMD = 20090714
      TES(1511)%NYMD = 20090714
      TES(1512)%NYMD = 20090714
      TES(1513)%NYMD = 20090714
      TES(1514)%NYMD = 20090714
      TES(1515)%NYMD = 20090714
      TES(1516)%NYMD = 20090714
      TES(1517)%NYMD = 20090714
      TES(1518)%NYMD = 20090714
      TES(1519)%NYMD = 20090714
      TES(1520)%NYMD = 20090714
      TES(1521)%NYMD = 20090714
      TES(1522)%NYMD = 20090714
      TES(1523)%NYMD = 20090714
      TES(1524)%NYMD = 20090714
      TES(1525)%NYMD = 20090714
      TES(1526)%NYMD = 20090714
      TES(1527)%NYMD = 20090714
      TES(1528)%NYMD = 20090714
      TES(1529)%NYMD = 20090714
      TES(1530)%NYMD = 20090714
      TES(1531)%NYMD = 20090714
      TES(1532)%NYMD = 20090714
      TES(1533)%NYMD = 20090714
      TES(1534)%NYMD = 20090714
      TES(1535)%NYMD = 20090714
      TES(1536)%NYMD = 20090714
      TES(1537)%NYMD = 20090714
      TES(1538)%NYMD = 20090714
      TES(1539)%NYMD = 20090714
      TES(1540)%NYMD = 20090714
      TES(1541)%NYMD = 20090714
      TES(1542)%NYMD = 20090714
      TES(1543)%NYMD = 20090714
      TES(1544)%NYMD = 20090714
      TES(1545)%NYMD = 20090714
      TES(1546)%NYMD = 20090714
      TES(1547)%NYMD = 20090714
      TES(1548)%NYMD = 20090714
      TES(1549)%NYMD = 20090714
      TES(1550)%NYMD = 20090714
      TES(1551)%NYMD = 20090714
      TES(1552)%NYMD = 20090714
      TES(1553)%NYMD = 20090714
      TES(1554)%NYMD = 20090714
      TES(1555)%NYMD = 20090714
      TES(1556)%NYMD = 20090714
      TES(1557)%NYMD = 20090714
      TES(1558)%NYMD = 20090714
      TES(1559)%NYMD = 20090714
      TES(1560)%NYMD = 20090714
      TES(1561)%NYMD = 20090714
      TES(1562)%NYMD = 20090714
      TES(1563)%NYMD = 20090714
      TES(1564)%NYMD = 20090714
      TES(1565)%NYMD = 20090714
      TES(1566)%NYMD = 20090714
      TES(1567)%NYMD = 20090714
      TES(1568)%NYMD = 20090714
      TES(1569)%NYMD = 20090714
      TES(1570)%NYMD = 20090714
      TES(1571)%NYMD = 20090714
      TES(1572)%NYMD = 20090714
      TES(1573)%NYMD = 20090714
      TES(1574)%NYMD = 20090714
      TES(1575)%NYMD = 20090714
      TES(1576)%NYMD = 20090714
      TES(1577)%NYMD = 20090714
      TES(1578)%NYMD = 20090714
      TES(1579)%NYMD = 20090714
      TES(1580)%NYMD = 20090714
      TES(1581)%NYMD = 20090714
      TES(1582)%NYMD = 20090714
      TES(1583)%NYMD = 20090714
      TES(1584)%NYMD = 20090714
      TES(1585)%NYMD = 20090714
      TES(1586)%NYMD = 20090714
      TES(1587)%NYMD = 20090714
      TES(1588)%NYMD = 20090714
      TES(1589)%NYMD = 20090714
      TES(1590)%NYMD = 20090714
      TES(1591)%NYMD = 20090714
      TES(1592)%NYMD = 20090714
      TES(1593)%NYMD = 20090714
      TES(1594)%NYMD = 20090714
      TES(1595)%NYMD = 20090714
      TES(1596)%NYMD = 20090714
      TES(1597)%NYMD = 20090714
      TES(1598)%NYMD = 20090714
      TES(1599)%NYMD = 20090714
      TES(1600)%NYMD = 20090714
      TES(1601)%NYMD = 20090714
      TES(1602)%NYMD = 20090714
      TES(1603)%NYMD = 20090714
      TES(1604)%NYMD = 20090714
      TES(1605)%NYMD = 20090714
      TES(1606)%NYMD = 20090714
      TES(1607)%NYMD = 20090714
      TES(1608)%NYMD = 20090714
      TES(1609)%NYMD = 20090714
      TES(1610)%NYMD = 20090714
      TES(1611)%NYMD = 20090714
      TES(1612)%NYMD = 20090714
      TES(1613)%NYMD = 20090714
      TES(1614)%NYMD = 20090714
      TES(1615)%NYMD = 20090714
      TES(1616)%NYMD = 20090714
      TES(1617)%NYMD = 20090714
      TES(1618)%NYMD = 20090714
      TES(1619)%NYMD = 20090714
      TES(1620)%NYMD = 20090714
      TES(1621)%NYMD = 20090714
      TES(1622)%NYMD = 20090714
      TES(1623)%NYMD = 20090714
      TES(1624)%NYMD = 20090714
      TES(1625)%NYMD = 20090714
      TES(1626)%NYMD = 20090714
      TES(1627)%NYMD = 20090714
      TES(1628)%NYMD = 20090714
      TES(1629)%NYMD = 20090714
      TES(1630)%NYMD = 20090714
      TES(1631)%NYMD = 20090714
      TES(1632)%NYMD = 20090714
      TES(1633)%NYMD = 20090714
      TES(1634)%NYMD = 20090714
      TES(1635)%NYMD = 20090714
      TES(1636)%NYMD = 20090714
      TES(1637)%NYMD = 20090714
      TES(1638)%NYMD = 20090714
      TES(1639)%NYMD = 20090714
      TES(1640)%NYMD = 20090714
      TES(1641)%NYMD = 20090714
      TES(1642)%NYMD = 20090714
      TES(1643)%NYMD = 20090715
      TES(1644)%NYMD = 20090715
      TES(1645)%NYMD = 20090715
      TES(1646)%NYMD = 20090715
      TES(1647)%NYMD = 20090715
      TES(1648)%NYMD = 20090715
      TES(1649)%NYMD = 20090715
      TES(1650)%NYMD = 20090715
      TES(1651)%NYMD = 20090715
      TES(1652)%NYMD = 20090715
      TES(1653)%NYMD = 20090715
      TES(1654)%NYMD = 20090715
      TES(1655)%NYMD = 20090715
      TES(1656)%NYMD = 20090715
      TES(1657)%NYMD = 20090715
      TES(1658)%NYMD = 20090715
      TES(1659)%NYMD = 20090715
      TES(1660)%NYMD = 20090715
      TES(1661)%NYMD = 20090715
      TES(1662)%NYMD = 20090715
      TES(1663)%NYMD = 20090715
      TES(1664)%NYMD = 20090715
      TES(1665)%NYMD = 20090715
      TES(1666)%NYMD = 20090715
      TES(1667)%NYMD = 20090715
      TES(1668)%NYMD = 20090715
      TES(1669)%NYMD = 20090715
      TES(1670)%NYMD = 20090715
      TES(1671)%NYMD = 20090715
      TES(1672)%NYMD = 20090715
      TES(1673)%NYMD = 20090715
      TES(1674)%NYMD = 20090715
      TES(1675)%NYMD = 20090715
      TES(1676)%NYMD = 20090715
      TES(1677)%NYMD = 20090715
      TES(1678)%NYMD = 20090715
      TES(1679)%NYMD = 20090715
      TES(1680)%NYMD = 20090715
      TES(1681)%NYMD = 20090715
      TES(1682)%NYMD = 20090715
      TES(1683)%NYMD = 20090715
      TES(1684)%NYMD = 20090715
      TES(1685)%NYMD = 20090715
      TES(1686)%NYMD = 20090715
      TES(1687)%NYMD = 20090715
      TES(1688)%NYMD = 20090715
      TES(1689)%NYMD = 20090715
      TES(1690)%NYMD = 20090715
      TES(1691)%NYMD = 20090715
      TES(1692)%NYMD = 20090715
      TES(1693)%NYMD = 20090715
      TES(1694)%NYMD = 20090715
      TES(1695)%NYMD = 20090715
      TES(1696)%NYMD = 20090715
      TES(1697)%NYMD = 20090715
      TES(1698)%NYMD = 20090715
      TES(1699)%NYMD = 20090715
      TES(1700)%NYMD = 20090715
      TES(1701)%NYMD = 20090715
      TES(1702)%NYMD = 20090715
      TES(1703)%NYMD = 20090715
      TES(1704)%NYMD = 20090715
      TES(1705)%NYMD = 20090715
      TES(1706)%NYMD = 20090715
      TES(1707)%NYMD = 20090715
      TES(1708)%NYMD = 20090715
      TES(1709)%NYMD = 20090715
      TES(1710)%NYMD = 20090715
      TES(1711)%NYMD = 20090715
      TES(1712)%NYMD = 20090715
      TES(1713)%NYMD = 20090715
      TES(1714)%NYMD = 20090715
      TES(1715)%NYMD = 20090715
      TES(1716)%NYMD = 20090715
      TES(1717)%NYMD = 20090715
      TES(1718)%NYMD = 20090715
      TES(1719)%NYMD = 20090715
      TES(1720)%NYMD = 20090715
      TES(1721)%NYMD = 20090715
      TES(1722)%NYMD = 20090715
      TES(1723)%NYMD = 20090715
      TES(1724)%NYMD = 20090715
      TES(1725)%NYMD = 20090715
      TES(1726)%NYMD = 20090715
      TES(1727)%NYMD = 20090715
      TES(1728)%NYMD = 20090715
      TES(1729)%NYMD = 20090715
      TES(1730)%NYMD = 20090715
      TES(1731)%NYMD = 20090715
      TES(1732)%NYMD = 20090715
      TES(1733)%NYMD = 20090715
      TES(1734)%NYMD = 20090715
      TES(1735)%NYMD = 20090715
      TES(1736)%NYMD = 20090715
      TES(1737)%NYMD = 20090715
      TES(1738)%NYMD = 20090715
      TES(1739)%NYMD = 20090715
      TES(1740)%NYMD = 20090715
      TES(1741)%NYMD = 20090715
      TES(1742)%NYMD = 20090715
      TES(1743)%NYMD = 20090715
      TES(1744)%NYMD = 20090715
      TES(1745)%NYMD = 20090715
      TES(1746)%NYMD = 20090715
      TES(1747)%NYMD = 20090715
      TES(1748)%NYMD = 20090715
      TES(1749)%NYMD = 20090715
      TES(1750)%NYMD = 20090715
      TES(1751)%NYMD = 20090715
      TES(1752)%NYMD = 20090715
      TES(1753)%NYMD = 20090715
      TES(1754)%NYMD = 20090715
      TES(1755)%NYMD = 20090715
      TES(1756)%NYMD = 20090715
      TES(1757)%NYMD = 20090715
      TES(1758)%NYMD = 20090715
      TES(1759)%NYMD = 20090715
      TES(1760)%NYMD = 20090715
      TES(1761)%NYMD = 20090715
      TES(1762)%NYMD = 20090715
      TES(1763)%NYMD = 20090715
      TES(1764)%NYMD = 20090715
      TES(1765)%NYMD = 20090715
      TES(1766)%NYMD = 20090715
      TES(1767)%NYMD = 20090715
      TES(1768)%NYMD = 20090715
      TES(1769)%NYMD = 20090715
      TES(1770)%NYMD = 20090715
      TES(1771)%NYMD = 20090715
      TES(1772)%NYMD = 20090715
      TES(1773)%NYMD = 20090715
      TES(1774)%NYMD = 20090715
      TES(1775)%NYMD = 20090715
      TES(1776)%NYMD = 20090715
      TES(1777)%NYMD = 20090715
      TES(1778)%NYMD = 20090715
      TES(1779)%NYMD = 20090715
      TES(1780)%NYMD = 20090715
      TES(1781)%NYMD = 20090715
      TES(1782)%NYMD = 20090715
      TES(1783)%NYMD = 20090715
      TES(1784)%NYMD = 20090715
      TES(1785)%NYMD = 20090716
      TES(1786)%NYMD = 20090716
      TES(1787)%NYMD = 20090716
      TES(1788)%NYMD = 20090716
      TES(1789)%NYMD = 20090716
      TES(1790)%NYMD = 20090716
      TES(1791)%NYMD = 20090716
      TES(1792)%NYMD = 20090716
      TES(1793)%NYMD = 20090716
      TES(1794)%NYMD = 20090716
      TES(1795)%NYMD = 20090716
      TES(1796)%NYMD = 20090716
      TES(1797)%NYMD = 20090716
      TES(1798)%NYMD = 20090716
      TES(1799)%NYMD = 20090716
      TES(1800)%NYMD = 20090716
      TES(1801)%NYMD = 20090716
      TES(1802)%NYMD = 20090716
      TES(1803)%NYMD = 20090716
      TES(1804)%NYMD = 20090716
      TES(1805)%NYMD = 20090716
      TES(1806)%NYMD = 20090716
      TES(1807)%NYMD = 20090716
      TES(1808)%NYMD = 20090716
      TES(1809)%NYMD = 20090716
      TES(1810)%NYMD = 20090716
      TES(1811)%NYMD = 20090716
      TES(1812)%NYMD = 20090716
      TES(1813)%NYMD = 20090716
      TES(1814)%NYMD = 20090716
      TES(1815)%NYMD = 20090716
      TES(1816)%NYMD = 20090716
      TES(1817)%NYMD = 20090716
      TES(1818)%NYMD = 20090716
      TES(1819)%NYMD = 20090716
      TES(1820)%NYMD = 20090716
      TES(1821)%NYMD = 20090716
      TES(1822)%NYMD = 20090716
      TES(1823)%NYMD = 20090716
      TES(1824)%NYMD = 20090716
      TES(1825)%NYMD = 20090716
      TES(1826)%NYMD = 20090716
      TES(1827)%NYMD = 20090716
      TES(1828)%NYMD = 20090716
      TES(1829)%NYMD = 20090716
      TES(1830)%NYMD = 20090716
      TES(1831)%NYMD = 20090716
      TES(1832)%NYMD = 20090716
      TES(1833)%NYMD = 20090716
      TES(1834)%NYMD = 20090716
      TES(1835)%NYMD = 20090716
      TES(1836)%NYMD = 20090716
      TES(1837)%NYMD = 20090716
      TES(1838)%NYMD = 20090716
      TES(1839)%NYMD = 20090716
      TES(1840)%NYMD = 20090716
      TES(1841)%NYMD = 20090716
      TES(1842)%NYMD = 20090716
      TES(1843)%NYMD = 20090716
      TES(1844)%NYMD = 20090716
      TES(1845)%NYMD = 20090716
      TES(1846)%NYMD = 20090716
      TES(1847)%NYMD = 20090716
      TES(1848)%NYMD = 20090716
      TES(1849)%NYMD = 20090716
      TES(1850)%NYMD = 20090716
      TES(1851)%NYMD = 20090716
      TES(1852)%NYMD = 20090716
      TES(1853)%NYMD = 20090716
      TES(1854)%NYMD = 20090716
      TES(1855)%NYMD = 20090716
      TES(1856)%NYMD = 20090716
      TES(1857)%NYMD = 20090716
      TES(1858)%NYMD = 20090716
      TES(1859)%NYMD = 20090716
      TES(1860)%NYMD = 20090716
      TES(1861)%NYMD = 20090716
      TES(1862)%NYMD = 20090716
      TES(1863)%NYMD = 20090716
      TES(1864)%NYMD = 20090716
      TES(1865)%NYMD = 20090716
      TES(1866)%NYMD = 20090716
      TES(1867)%NYMD = 20090716
      TES(1868)%NYMD = 20090716
      TES(1869)%NYMD = 20090716
      TES(1870)%NYMD = 20090716
      TES(1871)%NYMD = 20090716
      TES(1872)%NYMD = 20090716
      TES(1873)%NYMD = 20090716
      TES(1874)%NYMD = 20090716
      TES(1875)%NYMD = 20090716
      TES(1876)%NYMD = 20090716
      TES(1877)%NYMD = 20090716
      TES(1878)%NYMD = 20090716
      TES(1879)%NYMD = 20090716
      TES(1880)%NYMD = 20090716
      TES(1881)%NYMD = 20090716
      TES(1882)%NYMD = 20090716
      TES(1883)%NYMD = 20090716
      TES(1884)%NYMD = 20090716
      TES(1885)%NYMD = 20090716
      TES(1886)%NYMD = 20090716
      TES(1887)%NYMD = 20090716
      TES(1888)%NYMD = 20090716
      TES(1889)%NYMD = 20090716
      TES(1890)%NYMD = 20090716
      TES(1891)%NYMD = 20090716
      TES(1892)%NYMD = 20090716
      TES(1893)%NYMD = 20090716
      TES(1894)%NYMD = 20090716
      TES(1895)%NYMD = 20090716
      TES(1896)%NYMD = 20090716
      TES(1897)%NYMD = 20090716
      TES(1898)%NYMD = 20090716
      TES(1899)%NYMD = 20090716
      TES(1900)%NYMD = 20090716
      TES(1901)%NYMD = 20090716
      TES(1902)%NYMD = 20090716
      TES(1903)%NYMD = 20090716
      TES(1904)%NYMD = 20090716
      TES(1905)%NYMD = 20090716
      TES(1906)%NYMD = 20090716
      TES(1907)%NYMD = 20090716
      TES(1908)%NYMD = 20090716
      TES(1909)%NYMD = 20090716
      TES(1910)%NYMD = 20090716
      TES(1911)%NYMD = 20090716
      TES(1912)%NYMD = 20090716
      TES(1913)%NYMD = 20090716
      TES(1914)%NYMD = 20090716
      TES(1915)%NYMD = 20090716
      TES(1916)%NYMD = 20090716
      TES(1917)%NYMD = 20090716
      TES(1918)%NYMD = 20090716
      TES(1919)%NYMD = 20090716
      TES(1920)%NYMD = 20090716
      TES(1921)%NYMD = 20090716
      TES(1922)%NYMD = 20090716
      TES(1923)%NYMD = 20090716
      TES(1924)%NYMD = 20090716
      TES(1925)%NYMD = 20090716
      TES(1926)%NYMD = 20090717
      TES(1927)%NYMD = 20090717
      TES(1928)%NYMD = 20090717
      TES(1929)%NYMD = 20090717
      TES(1930)%NYMD = 20090717
      TES(1931)%NYMD = 20090717
      TES(1932)%NYMD = 20090717
      TES(1933)%NYMD = 20090717
      TES(1934)%NYMD = 20090717
      TES(1935)%NYMD = 20090717
      TES(1936)%NYMD = 20090717
      TES(1937)%NYMD = 20090717
      TES(1938)%NYMD = 20090717
      TES(1939)%NYMD = 20090717
      TES(1940)%NYMD = 20090717
      TES(1941)%NYMD = 20090717
      TES(1942)%NYMD = 20090717
      TES(1943)%NYMD = 20090717
      TES(1944)%NYMD = 20090717
      TES(1945)%NYMD = 20090717
      TES(1946)%NYMD = 20090717
      TES(1947)%NYMD = 20090717
      TES(1948)%NYMD = 20090717
      TES(1949)%NYMD = 20090717
      TES(1950)%NYMD = 20090717
      TES(1951)%NYMD = 20090717
      TES(1952)%NYMD = 20090717
      TES(1953)%NYMD = 20090717
      TES(1954)%NYMD = 20090717
      TES(1955)%NYMD = 20090717
      TES(1956)%NYMD = 20090717
      TES(1957)%NYMD = 20090717
      TES(1958)%NYMD = 20090717
      TES(1959)%NYMD = 20090717
      TES(1960)%NYMD = 20090717
      TES(1961)%NYMD = 20090717
      TES(1962)%NYMD = 20090717
      TES(1963)%NYMD = 20090717
      TES(1964)%NYMD = 20090717
      TES(1965)%NYMD = 20090717
      TES(1966)%NYMD = 20090717
      TES(1967)%NYMD = 20090717
      TES(1968)%NYMD = 20090717
      TES(1969)%NYMD = 20090717
      TES(1970)%NYMD = 20090717
      TES(1971)%NYMD = 20090717
      TES(1972)%NYMD = 20090717
      TES(1973)%NYMD = 20090717
      TES(1974)%NYMD = 20090717
      TES(1975)%NYMD = 20090717
      TES(1976)%NYMD = 20090717
      TES(1977)%NYMD = 20090717
      TES(1978)%NYMD = 20090717
      TES(1979)%NYMD = 20090717
      TES(1980)%NYMD = 20090717
      TES(1981)%NYMD = 20090717
      TES(1982)%NYMD = 20090717
      TES(1983)%NYMD = 20090717
      TES(1984)%NYMD = 20090717
      TES(1985)%NYMD = 20090717
      TES(1986)%NYMD = 20090717
      TES(1987)%NYMD = 20090717
      TES(1988)%NYMD = 20090717
      TES(1989)%NYMD = 20090717
      TES(1990)%NYMD = 20090717
      TES(1991)%NYMD = 20090717
      TES(1992)%NYMD = 20090717
      TES(1993)%NYMD = 20090717
      TES(1994)%NYMD = 20090717
      TES(1995)%NYMD = 20090717
      TES(1996)%NYMD = 20090717
      TES(1997)%NYMD = 20090717
      TES(1998)%NYMD = 20090717
      TES(1999)%NYMD = 20090717
      TES(2000)%NYMD = 20090717
      TES(2001)%NYMD = 20090717
      TES(2002)%NYMD = 20090717
      TES(2003)%NYMD = 20090717
      TES(2004)%NYMD = 20090717
      TES(2005)%NYMD = 20090717
      TES(2006)%NYMD = 20090717
      TES(2007)%NYMD = 20090717
      TES(2008)%NYMD = 20090717
      TES(2009)%NYMD = 20090717
      TES(2010)%NYMD = 20090717
      TES(2011)%NYMD = 20090717
      TES(2012)%NYMD = 20090717
      TES(2013)%NYMD = 20090717
      TES(2014)%NYMD = 20090717
      TES(2015)%NYMD = 20090717
      TES(2016)%NYMD = 20090717
      TES(2017)%NYMD = 20090717
      TES(2018)%NYMD = 20090717
      TES(2019)%NYMD = 20090717
      TES(2020)%NYMD = 20090717
      TES(2021)%NYMD = 20090717
      TES(2022)%NYMD = 20090717
      TES(2023)%NYMD = 20090717
      TES(2024)%NYMD = 20090717
      TES(2025)%NYMD = 20090717
      TES(2026)%NYMD = 20090717
      TES(2027)%NYMD = 20090717
      TES(2028)%NYMD = 20090717
      TES(2029)%NYMD = 20090717
      TES(2030)%NYMD = 20090717
      TES(2031)%NYMD = 20090717
      TES(2032)%NYMD = 20090717
      TES(2033)%NYMD = 20090717
      TES(2034)%NYMD = 20090717
      TES(2035)%NYMD = 20090717
      TES(2036)%NYMD = 20090717
      TES(2037)%NYMD = 20090717
      TES(2038)%NYMD = 20090717
      TES(2039)%NYMD = 20090717
      TES(2040)%NYMD = 20090717
      TES(2041)%NYMD = 20090717
      TES(2042)%NYMD = 20090717
      TES(2043)%NYMD = 20090717
      TES(2044)%NYMD = 20090717
      TES(2045)%NYMD = 20090717
      TES(2046)%NYMD = 20090717
      TES(2047)%NYMD = 20090717
      TES(2048)%NYMD = 20090717
      TES(2049)%NYMD = 20090717
      TES(2050)%NYMD = 20090717
      TES(2051)%NYMD = 20090717
      TES(2052)%NYMD = 20090717
      TES(2053)%NYMD = 20090717
      TES(2054)%NYMD = 20090717
      TES(2055)%NYMD = 20090717
      TES(2056)%NYMD = 20090717
      TES(2057)%NYMD = 20090717
      TES(2058)%NYMD = 20090717
      TES(2059)%NYMD = 20090717
      TES(2060)%NYMD = 20090717
      TES(2061)%NYMD = 20090717
      TES(2062)%NYMD = 20090717
      TES(2063)%NYMD = 20090717
      TES(2064)%NYMD = 20090717
      TES(2065)%NYMD = 20090717
      TES(2066)%NYMD = 20090717
      TES(2067)%NYMD = 20090717
      TES(2068)%NYMD = 20090717
      TES(2069)%NYMD = 20090717
      TES(2070)%NYMD = 20090717
      TES(2071)%NYMD = 20090717
      TES(2072)%NYMD = 20090717
      TES(2073)%NYMD = 20090718
      TES(2074)%NYMD = 20090718
      TES(2075)%NYMD = 20090718
      TES(2076)%NYMD = 20090718
      TES(2077)%NYMD = 20090718
      TES(2078)%NYMD = 20090718
      TES(2079)%NYMD = 20090718
      TES(2080)%NYMD = 20090718
      TES(2081)%NYMD = 20090718
      TES(2082)%NYMD = 20090718
      TES(2083)%NYMD = 20090718
      TES(2084)%NYMD = 20090718
      TES(2085)%NYMD = 20090718
      TES(2086)%NYMD = 20090718
      TES(2087)%NYMD = 20090718
      TES(2088)%NYMD = 20090718
      TES(2089)%NYMD = 20090718
      TES(2090)%NYMD = 20090718
      TES(2091)%NYMD = 20090718
      TES(2092)%NYMD = 20090718
      TES(2093)%NYMD = 20090718
      TES(2094)%NYMD = 20090718
      TES(2095)%NYMD = 20090718
      TES(2096)%NYMD = 20090718
      TES(2097)%NYMD = 20090718
      TES(2098)%NYMD = 20090718
      TES(2099)%NYMD = 20090718
      TES(2100)%NYMD = 20090718
      TES(2101)%NYMD = 20090718
      TES(2102)%NYMD = 20090718
      TES(2103)%NYMD = 20090718
      TES(2104)%NYMD = 20090718
      TES(2105)%NYMD = 20090718
      TES(2106)%NYMD = 20090718
      TES(2107)%NYMD = 20090718
      TES(2108)%NYMD = 20090718
      TES(2109)%NYMD = 20090718
      TES(2110)%NYMD = 20090718
      TES(2111)%NYMD = 20090718
      TES(2112)%NYMD = 20090718
      TES(2113)%NYMD = 20090718
      TES(2114)%NYMD = 20090718
      TES(2115)%NYMD = 20090718
      TES(2116)%NYMD = 20090718
      TES(2117)%NYMD = 20090718
      TES(2118)%NYMD = 20090718
      TES(2119)%NYMD = 20090718
      TES(2120)%NYMD = 20090718
      TES(2121)%NYMD = 20090718
      TES(2122)%NYMD = 20090718
      TES(2123)%NYMD = 20090718
      TES(2124)%NYMD = 20090718
      TES(2125)%NYMD = 20090718
      TES(2126)%NYMD = 20090718
      TES(2127)%NYMD = 20090718
      TES(2128)%NYMD = 20090718
      TES(2129)%NYMD = 20090718
      TES(2130)%NYMD = 20090718
      TES(2131)%NYMD = 20090718
      TES(2132)%NYMD = 20090718
      TES(2133)%NYMD = 20090718
      TES(2134)%NYMD = 20090718
      TES(2135)%NYMD = 20090718
      TES(2136)%NYMD = 20090718
      TES(2137)%NYMD = 20090718
      TES(2138)%NYMD = 20090718
      TES(2139)%NYMD = 20090718
      TES(2140)%NYMD = 20090718
      TES(2141)%NYMD = 20090718
      TES(2142)%NYMD = 20090718
      TES(2143)%NYMD = 20090718
      TES(2144)%NYMD = 20090718
      TES(2145)%NYMD = 20090718
      TES(2146)%NYMD = 20090718
      TES(2147)%NYMD = 20090718
      TES(2148)%NYMD = 20090718
      TES(2149)%NYMD = 20090718
      TES(2150)%NYMD = 20090718
      TES(2151)%NYMD = 20090718
      TES(2152)%NYMD = 20090718
      TES(2153)%NYMD = 20090718
      TES(2154)%NYMD = 20090718
      TES(2155)%NYMD = 20090718
      TES(2156)%NYMD = 20090718
      TES(2157)%NYMD = 20090718
      TES(2158)%NYMD = 20090718
      TES(2159)%NYMD = 20090718
      TES(2160)%NYMD = 20090718
      TES(2161)%NYMD = 20090718
      TES(2162)%NYMD = 20090718
      TES(2163)%NYMD = 20090718
      TES(2164)%NYMD = 20090718
      TES(2165)%NYMD = 20090718
      TES(2166)%NYMD = 20090718
      TES(2167)%NYMD = 20090718
      TES(2168)%NYMD = 20090718
      TES(2169)%NYMD = 20090718
      TES(2170)%NYMD = 20090718
      TES(2171)%NYMD = 20090718
      TES(2172)%NYMD = 20090718
      TES(2173)%NYMD = 20090718
      TES(2174)%NYMD = 20090718
      TES(2175)%NYMD = 20090718
      TES(2176)%NYMD = 20090718
      TES(2177)%NYMD = 20090718
      TES(2178)%NYMD = 20090718
      TES(2179)%NYMD = 20090718
      TES(2180)%NYMD = 20090718
      TES(2181)%NYMD = 20090718
      TES(2182)%NYMD = 20090718
      TES(2183)%NYMD = 20090718
      TES(2184)%NYMD = 20090718
      TES(2185)%NYMD = 20090718
      TES(2186)%NYMD = 20090718
      TES(2187)%NYMD = 20090718
      TES(2188)%NYMD = 20090718
      TES(2189)%NYMD = 20090719
      TES(2190)%NYMD = 20090719
      TES(2191)%NYMD = 20090719
      TES(2192)%NYMD = 20090719
      TES(2193)%NYMD = 20090719
      TES(2194)%NYMD = 20090719
      TES(2195)%NYMD = 20090719
      TES(2196)%NYMD = 20090719
      TES(2197)%NYMD = 20090719
      TES(2198)%NYMD = 20090719
      TES(2199)%NYMD = 20090719
      TES(2200)%NYMD = 20090719
      TES(2201)%NYMD = 20090719
      TES(2202)%NYMD = 20090719
      TES(2203)%NYMD = 20090719
      TES(2204)%NYMD = 20090719
      TES(2205)%NYMD = 20090719
      TES(2206)%NYMD = 20090719
      TES(2207)%NYMD = 20090719
      TES(2208)%NYMD = 20090719
      TES(2209)%NYMD = 20090719
      TES(2210)%NYMD = 20090719
      TES(2211)%NYMD = 20090719
      TES(2212)%NYMD = 20090719
      TES(2213)%NYMD = 20090719
      TES(2214)%NYMD = 20090719
      TES(2215)%NYMD = 20090719
      TES(2216)%NYMD = 20090719
      TES(2217)%NYMD = 20090719
      TES(2218)%NYMD = 20090719
      TES(2219)%NYMD = 20090719
      TES(2220)%NYMD = 20090719
      TES(2221)%NYMD = 20090719
      TES(2222)%NYMD = 20090719
      TES(2223)%NYMD = 20090719
      TES(2224)%NYMD = 20090719
      TES(2225)%NYMD = 20090719
      TES(2226)%NYMD = 20090719
      TES(2227)%NYMD = 20090719
      TES(2228)%NYMD = 20090719
      TES(2229)%NYMD = 20090719
      TES(2230)%NYMD = 20090719
      TES(2231)%NYMD = 20090719
      TES(2232)%NYMD = 20090719
      TES(2233)%NYMD = 20090719
      TES(2234)%NYMD = 20090719
      TES(2235)%NYMD = 20090719
      TES(2236)%NYMD = 20090719
      TES(2237)%NYMD = 20090719
      TES(2238)%NYMD = 20090719
      TES(2239)%NYMD = 20090719
      TES(2240)%NYMD = 20090719
      TES(2241)%NYMD = 20090719
      TES(2242)%NYMD = 20090719
      TES(2243)%NYMD = 20090719
      TES(2244)%NYMD = 20090719
      TES(2245)%NYMD = 20090719
      TES(2246)%NYMD = 20090719
      TES(2247)%NYMD = 20090719
      TES(2248)%NYMD = 20090719
      TES(2249)%NYMD = 20090719
      TES(2250)%NYMD = 20090719
      TES(2251)%NYMD = 20090719
      TES(2252)%NYMD = 20090719
      TES(2253)%NYMD = 20090719
      TES(2254)%NYMD = 20090719
      TES(2255)%NYMD = 20090719
      TES(2256)%NYMD = 20090719
      TES(2257)%NYMD = 20090719
      TES(2258)%NYMD = 20090719
      TES(2259)%NYMD = 20090719
      TES(2260)%NYMD = 20090719
      TES(2261)%NYMD = 20090719
      TES(2262)%NYMD = 20090719
      TES(2263)%NYMD = 20090719
      TES(2264)%NYMD = 20090719
      TES(2265)%NYMD = 20090719
      TES(2266)%NYMD = 20090719
      TES(2267)%NYMD = 20090719
      TES(2268)%NYMD = 20090719
      TES(2269)%NYMD = 20090719
      TES(2270)%NYMD = 20090719
      TES(2271)%NYMD = 20090719
      TES(2272)%NYMD = 20090719
      TES(2273)%NYMD = 20090719
      TES(2274)%NYMD = 20090719
      TES(2275)%NYMD = 20090719
      TES(2276)%NYMD = 20090719
      TES(2277)%NYMD = 20090719
      TES(2278)%NYMD = 20090719
      TES(2279)%NYMD = 20090719
      TES(2280)%NYMD = 20090719
      TES(2281)%NYMD = 20090719
      TES(2282)%NYMD = 20090719
      TES(2283)%NYMD = 20090719
      TES(2284)%NYMD = 20090719
      TES(2285)%NYMD = 20090719
      TES(2286)%NYMD = 20090719
      TES(2287)%NYMD = 20090719
      TES(2288)%NYMD = 20090719
      TES(2289)%NYMD = 20090719
      TES(2290)%NYMD = 20090719
      TES(2291)%NYMD = 20090719
      TES(2292)%NYMD = 20090719
      TES(2293)%NYMD = 20090719
      TES(2294)%NYMD = 20090719
      TES(2295)%NYMD = 20090719
      TES(2296)%NYMD = 20090719
      TES(2297)%NYMD = 20090719
      TES(2298)%NYMD = 20090719
      TES(2299)%NYMD = 20090719
      TES(2300)%NYMD = 20090719
      TES(2301)%NYMD = 20090719
      TES(2302)%NYMD = 20090719
      TES(2303)%NYMD = 20090719
      TES(2304)%NYMD = 20090719
      TES(2305)%NYMD = 20090719
      TES(2306)%NYMD = 20090719
      TES(2307)%NYMD = 20090719
      TES(2308)%NYMD = 20090719
      TES(2309)%NYMD = 20090719
      TES(2310)%NYMD = 20090719
      TES(2311)%NYMD = 20090719
      TES(2312)%NYMD = 20090719
      TES(2313)%NYMD = 20090719
      TES(2314)%NYMD = 20090719
      TES(2315)%NYMD = 20090719
      TES(2316)%NYMD = 20090719
      TES(2317)%NYMD = 20090719
      TES(2318)%NYMD = 20090719
      TES(2319)%NYMD = 20090719
      TES(2320)%NYMD = 20090719
      TES(2321)%NYMD = 20090719
      TES(2322)%NYMD = 20090719
      TES(2323)%NYMD = 20090719
      TES(2324)%NYMD = 20090720
      TES(2325)%NYMD = 20090720
      TES(2326)%NYMD = 20090720
      TES(2327)%NYMD = 20090720
      TES(2328)%NYMD = 20090720
      TES(2329)%NYMD = 20090720
      TES(2330)%NYMD = 20090720
      TES(2331)%NYMD = 20090720
      TES(2332)%NYMD = 20090720
      TES(2333)%NYMD = 20090720
      TES(2334)%NYMD = 20090720
      TES(2335)%NYMD = 20090720
      TES(2336)%NYMD = 20090720
      TES(2337)%NYMD = 20090720
      TES(2338)%NYMD = 20090720
      TES(2339)%NYMD = 20090720
      TES(2340)%NYMD = 20090720
      TES(2341)%NYMD = 20090720
      TES(2342)%NYMD = 20090720
      TES(2343)%NYMD = 20090720
      TES(2344)%NYMD = 20090720
      TES(2345)%NYMD = 20090720
      TES(2346)%NYMD = 20090720
      TES(2347)%NYMD = 20090720
      TES(2348)%NYMD = 20090720
      TES(2349)%NYMD = 20090720
      TES(2350)%NYMD = 20090720
      TES(2351)%NYMD = 20090720
      TES(2352)%NYMD = 20090720
      TES(2353)%NYMD = 20090720
      TES(2354)%NYMD = 20090720
      TES(2355)%NYMD = 20090720
      TES(2356)%NYMD = 20090720
      TES(2357)%NYMD = 20090720
      TES(2358)%NYMD = 20090720
      TES(2359)%NYMD = 20090720
      TES(2360)%NYMD = 20090720
      TES(2361)%NYMD = 20090720
      TES(2362)%NYMD = 20090720
      TES(2363)%NYMD = 20090720
      TES(2364)%NYMD = 20090720
      TES(2365)%NYMD = 20090720
      TES(2366)%NYMD = 20090720
      TES(2367)%NYMD = 20090720
      TES(2368)%NYMD = 20090720
      TES(2369)%NYMD = 20090720
      TES(2370)%NYMD = 20090720
      TES(2371)%NYMD = 20090720
      TES(2372)%NYMD = 20090720
      TES(2373)%NYMD = 20090720
      TES(2374)%NYMD = 20090720
      TES(2375)%NYMD = 20090720
      TES(2376)%NYMD = 20090720
      TES(2377)%NYMD = 20090720
      TES(2378)%NYMD = 20090720
      TES(2379)%NYMD = 20090720
      TES(2380)%NYMD = 20090720
      TES(2381)%NYMD = 20090720
      TES(2382)%NYMD = 20090720
      TES(2383)%NYMD = 20090720
      TES(2384)%NYMD = 20090720
      TES(2385)%NYMD = 20090720
      TES(2386)%NYMD = 20090720
      TES(2387)%NYMD = 20090720
      TES(2388)%NYMD = 20090720
      TES(2389)%NYMD = 20090720
      TES(2390)%NYMD = 20090720
      TES(2391)%NYMD = 20090720
      TES(2392)%NYMD = 20090720
      TES(2393)%NYMD = 20090720
      TES(2394)%NYMD = 20090720
      TES(2395)%NYMD = 20090720
      TES(2396)%NYMD = 20090720
      TES(2397)%NYMD = 20090720
      TES(2398)%NYMD = 20090720
      TES(2399)%NYMD = 20090720
      TES(2400)%NYMD = 20090720
      TES(2401)%NYMD = 20090720
      TES(2402)%NYMD = 20090720
      TES(2403)%NYMD = 20090720
      TES(2404)%NYMD = 20090720
      TES(2405)%NYMD = 20090720
      TES(2406)%NYMD = 20090720
      TES(2407)%NYMD = 20090720
      TES(2408)%NYMD = 20090720
      TES(2409)%NYMD = 20090720
      TES(2410)%NYMD = 20090720
      TES(2411)%NYMD = 20090720
      TES(2412)%NYMD = 20090720
      TES(2413)%NYMD = 20090720
      TES(2414)%NYMD = 20090720
      TES(2415)%NYMD = 20090720
      TES(2416)%NYMD = 20090720
      TES(2417)%NYMD = 20090720
      TES(2418)%NYMD = 20090720
      TES(2419)%NYMD = 20090720
      TES(2420)%NYMD = 20090720
      TES(2421)%NYMD = 20090720
      TES(2422)%NYMD = 20090720
      TES(2423)%NYMD = 20090720
      TES(2424)%NYMD = 20090720
      TES(2425)%NYMD = 20090720
      TES(2426)%NYMD = 20090720
      TES(2427)%NYMD = 20090720
      TES(2428)%NYMD = 20090720
      TES(2429)%NYMD = 20090720
      TES(2430)%NYMD = 20090720
      TES(2431)%NYMD = 20090720
      TES(2432)%NYMD = 20090721
      TES(2433)%NYMD = 20090721
      TES(2434)%NYMD = 20090721
      TES(2435)%NYMD = 20090721
      TES(2436)%NYMD = 20090721
      TES(2437)%NYMD = 20090721
      TES(2438)%NYMD = 20090721
      TES(2439)%NYMD = 20090721
      TES(2440)%NYMD = 20090721
      TES(2441)%NYMD = 20090721
      TES(2442)%NYMD = 20090721
      TES(2443)%NYMD = 20090721
      TES(2444)%NYMD = 20090721
      TES(2445)%NYMD = 20090721
      TES(2446)%NYMD = 20090721
      TES(2447)%NYMD = 20090721
      TES(2448)%NYMD = 20090721
      TES(2449)%NYMD = 20090721
      TES(2450)%NYMD = 20090721
      TES(2451)%NYMD = 20090721
      TES(2452)%NYMD = 20090721
      TES(2453)%NYMD = 20090721
      TES(2454)%NYMD = 20090721
      TES(2455)%NYMD = 20090721
      TES(2456)%NYMD = 20090721
      TES(2457)%NYMD = 20090721
      TES(2458)%NYMD = 20090721
      TES(2459)%NYMD = 20090721
      TES(2460)%NYMD = 20090721
      TES(2461)%NYMD = 20090721
      TES(2462)%NYMD = 20090721
      TES(2463)%NYMD = 20090721
      TES(2464)%NYMD = 20090721
      TES(2465)%NYMD = 20090721
      TES(2466)%NYMD = 20090721
      TES(2467)%NYMD = 20090721
      TES(2468)%NYMD = 20090721
      TES(2469)%NYMD = 20090721
      TES(2470)%NYMD = 20090721
      TES(2471)%NYMD = 20090721
      TES(2472)%NYMD = 20090721
      TES(2473)%NYMD = 20090721
      TES(2474)%NYMD = 20090721
      TES(2475)%NYMD = 20090721
      TES(2476)%NYMD = 20090721
      TES(2477)%NYMD = 20090721
      TES(2478)%NYMD = 20090721
      TES(2479)%NYMD = 20090721
      TES(2480)%NYMD = 20090721
      TES(2481)%NYMD = 20090721
      TES(2482)%NYMD = 20090721
      TES(2483)%NYMD = 20090721
      TES(2484)%NYMD = 20090721
      TES(2485)%NYMD = 20090721
      TES(2486)%NYMD = 20090721
      TES(2487)%NYMD = 20090721
      TES(2488)%NYMD = 20090721
      TES(2489)%NYMD = 20090721
      TES(2490)%NYMD = 20090721
      TES(2491)%NYMD = 20090721
      TES(2492)%NYMD = 20090721
      TES(2493)%NYMD = 20090721
      TES(2494)%NYMD = 20090721
      TES(2495)%NYMD = 20090721
      TES(2496)%NYMD = 20090721
      TES(2497)%NYMD = 20090721
      TES(2498)%NYMD = 20090721
      TES(2499)%NYMD = 20090721
      TES(2500)%NYMD = 20090721
      TES(2501)%NYMD = 20090721
      TES(2502)%NYMD = 20090721
      TES(2503)%NYMD = 20090721
      TES(2504)%NYMD = 20090721
      TES(2505)%NYMD = 20090721
      TES(2506)%NYMD = 20090721
      TES(2507)%NYMD = 20090721
      TES(2508)%NYMD = 20090721
      TES(2509)%NYMD = 20090721
      TES(2510)%NYMD = 20090721
      TES(2511)%NYMD = 20090721
      TES(2512)%NYMD = 20090721
      TES(2513)%NYMD = 20090721
      TES(2514)%NYMD = 20090721
      TES(2515)%NYMD = 20090721
      TES(2516)%NYMD = 20090721
      TES(2517)%NYMD = 20090721
      TES(2518)%NYMD = 20090721
      TES(2519)%NYMD = 20090721
      TES(2520)%NYMD = 20090721
      TES(2521)%NYMD = 20090721
      TES(2522)%NYMD = 20090721
      TES(2523)%NYMD = 20090721
      TES(2524)%NYMD = 20090721
      TES(2525)%NYMD = 20090721
      TES(2526)%NYMD = 20090721
      TES(2527)%NYMD = 20090721
      TES(2528)%NYMD = 20090721
      TES(2529)%NYMD = 20090721
      TES(2530)%NYMD = 20090721
      TES(2531)%NYMD = 20090721
      TES(2532)%NYMD = 20090721
      TES(2533)%NYMD = 20090721
      TES(2534)%NYMD = 20090721
      TES(2535)%NYMD = 20090721
      TES(2536)%NYMD = 20090721
      TES(2537)%NYMD = 20090721
      TES(2538)%NYMD = 20090721
      TES(2539)%NYMD = 20090721
      TES(2540)%NYMD = 20090721
      TES(2541)%NYMD = 20090721
      TES(2542)%NYMD = 20090721
      TES(2543)%NYMD = 20090721
      TES(2544)%NYMD = 20090721
      TES(2545)%NYMD = 20090721
      TES(2546)%NYMD = 20090721
      TES(2547)%NYMD = 20090721
      TES(2548)%NYMD = 20090721
      TES(2549)%NYMD = 20090721
      TES(2550)%NYMD = 20090721
      TES(2551)%NYMD = 20090721
      TES(2552)%NYMD = 20090721
      TES(2553)%NYMD = 20090721
      TES(2554)%NYMD = 20090721
      TES(2555)%NYMD = 20090721
      TES(2556)%NYMD = 20090721
      TES(2557)%NYMD = 20090721
      TES(2558)%NYMD = 20090721
      TES(2559)%NYMD = 20090721
      TES(2560)%NYMD = 20090721
      TES(2561)%NYMD = 20090721
      TES(2562)%NYMD = 20090721
      TES(2563)%NYMD = 20090721
      TES(2564)%NYMD = 20090721
      TES(2565)%NYMD = 20090721
      TES(2566)%NYMD = 20090721
      TES(2567)%NYMD = 20090721
      TES(2568)%NYMD = 20090721
      TES(2569)%NYMD = 20090721
      TES(2570)%NYMD = 20090721
      TES(2571)%NYMD = 20090721
      TES(2572)%NYMD = 20090721
      TES(2573)%NYMD = 20090721
      TES(2574)%NYMD = 20090721
      TES(2575)%NYMD = 20090721
      TES(2576)%NYMD = 20090721
      TES(2577)%NYMD = 20090721
      TES(2578)%NYMD = 20090721
      TES(2579)%NYMD = 20090721
      TES(2580)%NYMD = 20090721
      TES(2581)%NYMD = 20090721
      TES(2582)%NYMD = 20090721
      TES(2583)%NYMD = 20090721
      TES(2584)%NYMD = 20090721
      TES(2585)%NYMD = 20090721
      TES(2586)%NYMD = 20090721
      TES(2587)%NYMD = 20090721
      TES(2588)%NYMD = 20090721
      TES(2589)%NYMD = 20090721
      TES(2590)%NYMD = 20090721
      TES(2591)%NYMD = 20090721
      TES(2592)%NYMD = 20090721
      TES(2593)%NYMD = 20090721
      TES(2594)%NYMD = 20090721
      TES(2595)%NYMD = 20090721
      TES(2596)%NYMD = 20090721
      TES(2597)%NYMD = 20090721
      TES(2598)%NYMD = 20090721
      TES(2599)%NYMD = 20090721
      TES(2600)%NYMD = 20090721
      TES(2601)%NYMD = 20090721
      TES(2602)%NYMD = 20090721
      TES(2603)%NYMD = 20090721
      TES(2604)%NYMD = 20090722
      TES(2605)%NYMD = 20090722
      TES(2606)%NYMD = 20090722
      TES(2607)%NYMD = 20090722
      TES(2608)%NYMD = 20090722
      TES(2609)%NYMD = 20090722
      TES(2610)%NYMD = 20090722
      TES(2611)%NYMD = 20090722
      TES(2612)%NYMD = 20090722
      TES(2613)%NYMD = 20090722
      TES(2614)%NYMD = 20090722
      TES(2615)%NYMD = 20090722
      TES(2616)%NYMD = 20090722
      TES(2617)%NYMD = 20090722
      TES(2618)%NYMD = 20090722
      TES(2619)%NYMD = 20090722
      TES(2620)%NYMD = 20090722
      TES(2621)%NYMD = 20090722
      TES(2622)%NYMD = 20090722
      TES(2623)%NYMD = 20090722
      TES(2624)%NYMD = 20090722
      TES(2625)%NYMD = 20090722
      TES(2626)%NYMD = 20090722
      TES(2627)%NYMD = 20090722
      TES(2628)%NYMD = 20090722
      TES(2629)%NYMD = 20090722
      TES(2630)%NYMD = 20090722
      TES(2631)%NYMD = 20090722
      TES(2632)%NYMD = 20090722
      TES(2633)%NYMD = 20090722
      TES(2634)%NYMD = 20090722
      TES(2635)%NYMD = 20090722
      TES(2636)%NYMD = 20090722
      TES(2637)%NYMD = 20090722
      TES(2638)%NYMD = 20090722
      TES(2639)%NYMD = 20090722
      TES(2640)%NYMD = 20090722
      TES(2641)%NYMD = 20090722
      TES(2642)%NYMD = 20090722
      TES(2643)%NYMD = 20090722
      TES(2644)%NYMD = 20090722
      TES(2645)%NYMD = 20090722
      TES(2646)%NYMD = 20090722
      TES(2647)%NYMD = 20090722
      TES(2648)%NYMD = 20090722
      TES(2649)%NYMD = 20090722
      TES(2650)%NYMD = 20090722
      TES(2651)%NYMD = 20090722
      TES(2652)%NYMD = 20090722
      TES(2653)%NYMD = 20090722
      TES(2654)%NYMD = 20090722
      TES(2655)%NYMD = 20090722
      TES(2656)%NYMD = 20090722
      TES(2657)%NYMD = 20090722
      TES(2658)%NYMD = 20090722
      TES(2659)%NYMD = 20090722
      TES(2660)%NYMD = 20090722
      TES(2661)%NYMD = 20090722
      TES(2662)%NYMD = 20090722
      TES(2663)%NYMD = 20090722
      TES(2664)%NYMD = 20090722
      TES(2665)%NYMD = 20090722
      TES(2666)%NYMD = 20090722
      TES(2667)%NYMD = 20090722
      TES(2668)%NYMD = 20090722
      TES(2669)%NYMD = 20090722
      TES(2670)%NYMD = 20090722
      TES(2671)%NYMD = 20090722
      TES(2672)%NYMD = 20090722
      TES(2673)%NYMD = 20090722
      TES(2674)%NYMD = 20090722
      TES(2675)%NYMD = 20090722
      TES(2676)%NYMD = 20090722
      TES(2677)%NYMD = 20090722
      TES(2678)%NYMD = 20090722
      TES(2679)%NYMD = 20090722
      TES(2680)%NYMD = 20090722
      TES(2681)%NYMD = 20090722
      TES(2682)%NYMD = 20090722
      TES(2683)%NYMD = 20090722
      TES(2684)%NYMD = 20090722
      TES(2685)%NYMD = 20090722
      TES(2686)%NYMD = 20090722
      TES(2687)%NYMD = 20090722
      TES(2688)%NYMD = 20090722
      TES(2689)%NYMD = 20090722
      TES(2690)%NYMD = 20090722
      TES(2691)%NYMD = 20090722
      TES(2692)%NYMD = 20090722
      TES(2693)%NYMD = 20090722
      TES(2694)%NYMD = 20090722
      TES(2695)%NYMD = 20090722
      TES(2696)%NYMD = 20090722
      TES(2697)%NYMD = 20090722
      TES(2698)%NYMD = 20090722
      TES(2699)%NYMD = 20090722
      TES(2700)%NYMD = 20090722
      TES(2701)%NYMD = 20090722
      TES(2702)%NYMD = 20090722
      TES(2703)%NYMD = 20090722
      TES(2704)%NYMD = 20090722
      TES(2705)%NYMD = 20090722
      TES(2706)%NYMD = 20090722
      TES(2707)%NYMD = 20090722
      TES(2708)%NYMD = 20090722
      TES(2709)%NYMD = 20090722
      TES(2710)%NYMD = 20090722
      TES(2711)%NYMD = 20090722
      TES(2712)%NYMD = 20090722
      TES(2713)%NYMD = 20090722
      TES(2714)%NYMD = 20090722
      TES(2715)%NYMD = 20090722
      TES(2716)%NYMD = 20090722
      TES(2717)%NYMD = 20090722
      TES(2718)%NYMD = 20090722
      TES(2719)%NYMD = 20090722
      TES(2720)%NYMD = 20090722
      TES(2721)%NYMD = 20090722
      TES(2722)%NYMD = 20090722
      TES(2723)%NYMD = 20090722
      TES(2724)%NYMD = 20090723
      TES(2725)%NYMD = 20090723
      TES(2726)%NYMD = 20090723
      TES(2727)%NYMD = 20090723
      TES(2728)%NYMD = 20090723
      TES(2729)%NYMD = 20090723
      TES(2730)%NYMD = 20090723
      TES(2731)%NYMD = 20090723
      TES(2732)%NYMD = 20090723
      TES(2733)%NYMD = 20090723
      TES(2734)%NYMD = 20090723
      TES(2735)%NYMD = 20090723
      TES(2736)%NYMD = 20090723
      TES(2737)%NYMD = 20090723
      TES(2738)%NYMD = 20090723
      TES(2739)%NYMD = 20090723
      TES(2740)%NYMD = 20090723
      TES(2741)%NYMD = 20090723
      TES(2742)%NYMD = 20090723
      TES(2743)%NYMD = 20090723
      TES(2744)%NYMD = 20090723
      TES(2745)%NYMD = 20090723
      TES(2746)%NYMD = 20090723
      TES(2747)%NYMD = 20090723
      TES(2748)%NYMD = 20090723
      TES(2749)%NYMD = 20090723
      TES(2750)%NYMD = 20090723
      TES(2751)%NYMD = 20090723
      TES(2752)%NYMD = 20090723
      TES(2753)%NYMD = 20090723
      TES(2754)%NYMD = 20090723
      TES(2755)%NYMD = 20090723
      TES(2756)%NYMD = 20090723
      TES(2757)%NYMD = 20090723
      TES(2758)%NYMD = 20090723
      TES(2759)%NYMD = 20090723
      TES(2760)%NYMD = 20090723
      TES(2761)%NYMD = 20090723
      TES(2762)%NYMD = 20090723
      TES(2763)%NYMD = 20090723
      TES(2764)%NYMD = 20090723
      TES(2765)%NYMD = 20090723
      TES(2766)%NYMD = 20090723
      TES(2767)%NYMD = 20090723
      TES(2768)%NYMD = 20090723
      TES(2769)%NYMD = 20090723
      TES(2770)%NYMD = 20090723
      TES(2771)%NYMD = 20090723
      TES(2772)%NYMD = 20090723
      TES(2773)%NYMD = 20090723
      TES(2774)%NYMD = 20090723
      TES(2775)%NYMD = 20090723
      TES(2776)%NYMD = 20090723
      TES(2777)%NYMD = 20090723
      TES(2778)%NYMD = 20090723
      TES(2779)%NYMD = 20090723
      TES(2780)%NYMD = 20090723
      TES(2781)%NYMD = 20090723
      TES(2782)%NYMD = 20090723
      TES(2783)%NYMD = 20090723
      TES(2784)%NYMD = 20090723
      TES(2785)%NYMD = 20090723
      TES(2786)%NYMD = 20090723
      TES(2787)%NYMD = 20090723
      TES(2788)%NYMD = 20090723
      TES(2789)%NYMD = 20090723
      TES(2790)%NYMD = 20090723
      TES(2791)%NYMD = 20090723
      TES(2792)%NYMD = 20090723
      TES(2793)%NYMD = 20090723
      TES(2794)%NYMD = 20090723
      TES(2795)%NYMD = 20090723
      TES(2796)%NYMD = 20090723
      TES(2797)%NYMD = 20090723
      TES(2798)%NYMD = 20090723
      TES(2799)%NYMD = 20090723
      TES(2800)%NYMD = 20090723
      TES(2801)%NYMD = 20090723
      TES(2802)%NYMD = 20090723
      TES(2803)%NYMD = 20090723
      TES(2804)%NYMD = 20090723
      TES(2805)%NYMD = 20090723
      TES(2806)%NYMD = 20090723
      TES(2807)%NYMD = 20090723
      TES(2808)%NYMD = 20090723
      TES(2809)%NYMD = 20090723
      TES(2810)%NYMD = 20090723
      TES(2811)%NYMD = 20090723
      TES(2812)%NYMD = 20090723
      TES(2813)%NYMD = 20090723
      TES(2814)%NYMD = 20090723
      TES(2815)%NYMD = 20090723
      TES(2816)%NYMD = 20090723
      TES(2817)%NYMD = 20090723
      TES(2818)%NYMD = 20090723
      TES(2819)%NYMD = 20090723
      TES(2820)%NYMD = 20090723
      TES(2821)%NYMD = 20090723
      TES(2822)%NYMD = 20090723
      TES(2823)%NYMD = 20090723
      TES(2824)%NYMD = 20090723
      TES(2825)%NYMD = 20090723
      TES(2826)%NYMD = 20090723
      TES(2827)%NYMD = 20090723
      TES(2828)%NYMD = 20090723
      TES(2829)%NYMD = 20090723
      TES(2830)%NYMD = 20090723
      TES(2831)%NYMD = 20090723
      TES(2832)%NYMD = 20090723
      TES(2833)%NYMD = 20090723
      TES(2834)%NYMD = 20090723
      TES(2835)%NYMD = 20090723
      TES(2836)%NYMD = 20090723
      TES(2837)%NYMD = 20090723
      TES(2838)%NYMD = 20090723
      TES(2839)%NYMD = 20090723
      TES(2840)%NYMD = 20090723
      TES(2841)%NYMD = 20090723
      TES(2842)%NYMD = 20090723
      TES(2843)%NYMD = 20090723
      TES(2844)%NYMD = 20090723
      TES(2845)%NYMD = 20090723
      TES(2846)%NYMD = 20090723
      TES(2847)%NYMD = 20090723
      TES(2848)%NYMD = 20090723
      TES(2849)%NYMD = 20090723
      TES(2850)%NYMD = 20090723
      TES(2851)%NYMD = 20090723
      TES(2852)%NYMD = 20090723
      TES(2853)%NYMD = 20090723
      TES(2854)%NYMD = 20090723
      TES(2855)%NYMD = 20090723
      TES(2856)%NYMD = 20090723
      TES(2857)%NYMD = 20090723
      TES(2858)%NYMD = 20090723
      TES(2859)%NYMD = 20090723
      TES(2860)%NYMD = 20090723
      TES(2861)%NYMD = 20090723
      TES(2862)%NYMD = 20090723
      TES(2863)%NYMD = 20090723
      TES(2864)%NYMD = 20090723
      TES(2865)%NYMD = 20090723
      TES(2866)%NYMD = 20090723
      TES(2867)%NYMD = 20090723
      TES(2868)%NYMD = 20090723
      TES(2869)%NYMD = 20090723
      TES(2870)%NYMD = 20090723
      TES(2871)%NYMD = 20090723
      TES(2872)%NYMD = 20090723
      TES(2873)%NYMD = 20090723
      TES(2874)%NYMD = 20090723
      TES(2875)%NYMD = 20090723
      TES(2876)%NYMD = 20090723
      TES(2877)%NYMD = 20090723
      TES(2878)%NYMD = 20090723
      TES(2879)%NYMD = 20090724
      TES(2880)%NYMD = 20090724
      TES(2881)%NYMD = 20090724
      TES(2882)%NYMD = 20090724
      TES(2883)%NYMD = 20090724
      TES(2884)%NYMD = 20090724
      TES(2885)%NYMD = 20090724
      TES(2886)%NYMD = 20090724
      TES(2887)%NYMD = 20090724
      TES(2888)%NYMD = 20090724
      TES(2889)%NYMD = 20090724
      TES(2890)%NYMD = 20090724
      TES(2891)%NYMD = 20090724
      TES(2892)%NYMD = 20090724
      TES(2893)%NYMD = 20090724
      TES(2894)%NYMD = 20090724
      TES(2895)%NYMD = 20090724
      TES(2896)%NYMD = 20090724
      TES(2897)%NYMD = 20090724
      TES(2898)%NYMD = 20090724
      TES(2899)%NYMD = 20090724
      TES(2900)%NYMD = 20090724
      TES(2901)%NYMD = 20090724
      TES(2902)%NYMD = 20090724
      TES(2903)%NYMD = 20090724
      TES(2904)%NYMD = 20090724
      TES(2905)%NYMD = 20090724
      TES(2906)%NYMD = 20090724
      TES(2907)%NYMD = 20090724
      TES(2908)%NYMD = 20090724
      TES(2909)%NYMD = 20090724
      TES(2910)%NYMD = 20090724
      TES(2911)%NYMD = 20090724
      TES(2912)%NYMD = 20090724
      TES(2913)%NYMD = 20090724
      TES(2914)%NYMD = 20090724
      TES(2915)%NYMD = 20090724
      TES(2916)%NYMD = 20090724
      TES(2917)%NYMD = 20090724
      TES(2918)%NYMD = 20090724
      TES(2919)%NYMD = 20090724
      TES(2920)%NYMD = 20090724
      TES(2921)%NYMD = 20090724
      TES(2922)%NYMD = 20090724
      TES(2923)%NYMD = 20090724
      TES(2924)%NYMD = 20090724
      TES(2925)%NYMD = 20090724
      TES(2926)%NYMD = 20090724
      TES(2927)%NYMD = 20090724
      TES(2928)%NYMD = 20090724
      TES(2929)%NYMD = 20090724
      TES(2930)%NYMD = 20090724
      TES(2931)%NYMD = 20090724
      TES(2932)%NYMD = 20090724
      TES(2933)%NYMD = 20090724
      TES(2934)%NYMD = 20090724
      TES(2935)%NYMD = 20090724
      TES(2936)%NYMD = 20090724
      TES(2937)%NYMD = 20090724
      TES(2938)%NYMD = 20090724
      TES(2939)%NYMD = 20090724
      TES(2940)%NYMD = 20090724
      TES(2941)%NYMD = 20090724
      TES(2942)%NYMD = 20090724
      TES(2943)%NYMD = 20090724
      TES(2944)%NYMD = 20090724
      TES(2945)%NYMD = 20090724
      TES(2946)%NYMD = 20090724
      TES(2947)%NYMD = 20090724
      TES(2948)%NYMD = 20090724
      TES(2949)%NYMD = 20090724
      TES(2950)%NYMD = 20090724
      TES(2951)%NYMD = 20090724
      TES(2952)%NYMD = 20090724
      TES(2953)%NYMD = 20090724
      TES(2954)%NYMD = 20090724
      TES(2955)%NYMD = 20090724
      TES(2956)%NYMD = 20090724
      TES(2957)%NYMD = 20090724
      TES(2958)%NYMD = 20090724
      TES(2959)%NYMD = 20090724
      TES(2960)%NYMD = 20090724
      TES(2961)%NYMD = 20090724
      TES(2962)%NYMD = 20090724
      TES(2963)%NYMD = 20090724
      TES(2964)%NYMD = 20090724
      TES(2965)%NYMD = 20090724
      TES(2966)%NYMD = 20090724
      TES(2967)%NYMD = 20090724
      TES(2968)%NYMD = 20090724
      TES(2969)%NYMD = 20090724
      TES(2970)%NYMD = 20090724
      TES(2971)%NYMD = 20090724
      TES(2972)%NYMD = 20090724
      TES(2973)%NYMD = 20090724
      TES(2974)%NYMD = 20090724
      TES(2975)%NYMD = 20090724
      TES(2976)%NYMD = 20090724
      TES(2977)%NYMD = 20090724
      TES(2978)%NYMD = 20090724
      TES(2979)%NYMD = 20090724
      TES(2980)%NYMD = 20090724
      TES(2981)%NYMD = 20090724
      TES(2982)%NYMD = 20090724
      TES(2983)%NYMD = 20090724
      TES(2984)%NYMD = 20090724
      TES(2985)%NYMD = 20090724
      TES(2986)%NYMD = 20090724
      TES(2987)%NYMD = 20090724
      TES(2988)%NYMD = 20090724
      TES(2989)%NYMD = 20090724
      TES(2990)%NYMD = 20090724
      TES(2991)%NYMD = 20090724
      TES(2992)%NYMD = 20090724
      TES(2993)%NYMD = 20090724
      TES(2994)%NYMD = 20090724
      TES(2995)%NYMD = 20090724
      TES(2996)%NYMD = 20090724
      TES(2997)%NYMD = 20090724
      TES(2998)%NYMD = 20090724
      TES(2999)%NYMD = 20090724
      TES(3000)%NYMD = 20090724
      TES(3001)%NYMD = 20090724
      TES(3002)%NYMD = 20090724
      TES(3003)%NYMD = 20090724
      TES(3004)%NYMD = 20090724
      TES(3005)%NYMD = 20090724
      TES(3006)%NYMD = 20090724
      TES(3007)%NYMD = 20090725
      TES(3008)%NYMD = 20090725
      TES(3009)%NYMD = 20090725
      TES(3010)%NYMD = 20090725
      TES(3011)%NYMD = 20090725
      TES(3012)%NYMD = 20090725
      TES(3013)%NYMD = 20090725
      TES(3014)%NYMD = 20090725
      TES(3015)%NYMD = 20090725
      TES(3016)%NYMD = 20090725
      TES(3017)%NYMD = 20090725
      TES(3018)%NYMD = 20090725
      TES(3019)%NYMD = 20090725
      TES(3020)%NYMD = 20090725
      TES(3021)%NYMD = 20090725
      TES(3022)%NYMD = 20090725
      TES(3023)%NYMD = 20090725
      TES(3024)%NYMD = 20090725
      TES(3025)%NYMD = 20090725
      TES(3026)%NYMD = 20090725
      TES(3027)%NYMD = 20090725
      TES(3028)%NYMD = 20090725
      TES(3029)%NYMD = 20090725
      TES(3030)%NYMD = 20090725
      TES(3031)%NYMD = 20090725
      TES(3032)%NYMD = 20090725
      TES(3033)%NYMD = 20090725
      TES(3034)%NYMD = 20090725
      TES(3035)%NYMD = 20090725
      TES(3036)%NYMD = 20090725
      TES(3037)%NYMD = 20090725
      TES(3038)%NYMD = 20090725
      TES(3039)%NYMD = 20090725
      TES(3040)%NYMD = 20090725
      TES(3041)%NYMD = 20090725
      TES(3042)%NYMD = 20090725
      TES(3043)%NYMD = 20090725
      TES(3044)%NYMD = 20090725
      TES(3045)%NYMD = 20090725
      TES(3046)%NYMD = 20090725
      TES(3047)%NYMD = 20090725
      TES(3048)%NYMD = 20090725
      TES(3049)%NYMD = 20090725
      TES(3050)%NYMD = 20090725
      TES(3051)%NYMD = 20090725
      TES(3052)%NYMD = 20090725
      TES(3053)%NYMD = 20090725
      TES(3054)%NYMD = 20090725
      TES(3055)%NYMD = 20090725
      TES(3056)%NYMD = 20090725
      TES(3057)%NYMD = 20090725
      TES(3058)%NYMD = 20090725
      TES(3059)%NYMD = 20090725
      TES(3060)%NYMD = 20090725
      TES(3061)%NYMD = 20090725
      TES(3062)%NYMD = 20090725
      TES(3063)%NYMD = 20090725
      TES(3064)%NYMD = 20090725
      TES(3065)%NYMD = 20090725
      TES(3066)%NYMD = 20090725
      TES(3067)%NYMD = 20090725
      TES(3068)%NYMD = 20090725
      TES(3069)%NYMD = 20090725
      TES(3070)%NYMD = 20090725
      TES(3071)%NYMD = 20090725
      TES(3072)%NYMD = 20090725
      TES(3073)%NYMD = 20090725
      TES(3074)%NYMD = 20090725
      TES(3075)%NYMD = 20090725
      TES(3076)%NYMD = 20090725
      TES(3077)%NYMD = 20090725
      TES(3078)%NYMD = 20090725
      TES(3079)%NYMD = 20090725
      TES(3080)%NYMD = 20090725
      TES(3081)%NYMD = 20090725
      TES(3082)%NYMD = 20090725
      TES(3083)%NYMD = 20090725
      TES(3084)%NYMD = 20090725
      TES(3085)%NYMD = 20090725
      TES(3086)%NYMD = 20090725
      TES(3087)%NYMD = 20090725
      TES(3088)%NYMD = 20090725
      TES(3089)%NYMD = 20090725
      TES(3090)%NYMD = 20090725
      TES(3091)%NYMD = 20090725
      TES(3092)%NYMD = 20090725
      TES(3093)%NYMD = 20090725
      TES(3094)%NYMD = 20090725
      TES(3095)%NYMD = 20090725
      TES(3096)%NYMD = 20090725
      TES(3097)%NYMD = 20090725
      TES(3098)%NYMD = 20090725
      TES(3099)%NYMD = 20090725
      TES(3100)%NYMD = 20090725
      TES(3101)%NYMD = 20090725
      TES(3102)%NYMD = 20090725
      TES(3103)%NYMD = 20090725
      TES(3104)%NYMD = 20090725
      TES(3105)%NYMD = 20090725
      TES(3106)%NYMD = 20090725
      TES(3107)%NYMD = 20090725
      TES(3108)%NYMD = 20090725
      TES(3109)%NYMD = 20090725
      TES(3110)%NYMD = 20090725
      TES(3111)%NYMD = 20090725
      TES(3112)%NYMD = 20090725
      TES(3113)%NYMD = 20090725
      TES(3114)%NYMD = 20090725
      TES(3115)%NYMD = 20090725
      TES(3116)%NYMD = 20090725
      TES(3117)%NYMD = 20090725
      TES(3118)%NYMD = 20090725
      TES(3119)%NYMD = 20090725
      TES(3120)%NYMD = 20090725
      TES(3121)%NYMD = 20090725
      TES(3122)%NYMD = 20090725
      TES(3123)%NYMD = 20090725
      TES(3124)%NYMD = 20090725
      TES(3125)%NYMD = 20090725
      TES(3126)%NYMD = 20090725
      TES(3127)%NYMD = 20090725
      TES(3128)%NYMD = 20090725
      TES(3129)%NYMD = 20090725
      TES(3130)%NYMD = 20090725
      TES(3131)%NYMD = 20090725
      TES(3132)%NYMD = 20090725
      TES(3133)%NYMD = 20090725
      TES(3134)%NYMD = 20090725
      TES(3135)%NYMD = 20090725
      TES(3136)%NYMD = 20090725
      TES(3137)%NYMD = 20090725
      TES(3138)%NYMD = 20090725
      TES(3139)%NYMD = 20090725
      TES(3140)%NYMD = 20090725
      TES(3141)%NYMD = 20090725
      TES(3142)%NYMD = 20090725
      TES(3143)%NYMD = 20090725
      TES(3144)%NYMD = 20090725
      TES(3145)%NYMD = 20090725
      TES(3146)%NYMD = 20090725
      TES(3147)%NYMD = 20090725
      TES(3148)%NYMD = 20090725
      TES(3149)%NYMD = 20090725
      TES(3150)%NYMD = 20090725
      TES(3151)%NYMD = 20090725
      TES(3152)%NYMD = 20090725
      TES(3153)%NYMD = 20090725
      TES(3154)%NYMD = 20090725
      TES(3155)%NYMD = 20090725
      TES(3156)%NYMD = 20090725
      TES(3157)%NYMD = 20090725
      TES(3158)%NYMD = 20090725
      TES(3159)%NYMD = 20090725
      TES(3160)%NYMD = 20090725
      TES(3161)%NYMD = 20090725
      TES(3162)%NYMD = 20090725
      TES(3163)%NYMD = 20090726
      TES(3164)%NYMD = 20090726
      TES(3165)%NYMD = 20090726
      TES(3166)%NYMD = 20090726
      TES(3167)%NYMD = 20090726
      TES(3168)%NYMD = 20090726
      TES(3169)%NYMD = 20090726
      TES(3170)%NYMD = 20090726
      TES(3171)%NYMD = 20090726
      TES(3172)%NYMD = 20090726
      TES(3173)%NYMD = 20090726
      TES(3174)%NYMD = 20090726
      TES(3175)%NYMD = 20090726
      TES(3176)%NYMD = 20090726
      TES(3177)%NYMD = 20090726
      TES(3178)%NYMD = 20090726
      TES(3179)%NYMD = 20090726
      TES(3180)%NYMD = 20090726
      TES(3181)%NYMD = 20090726
      TES(3182)%NYMD = 20090726
      TES(3183)%NYMD = 20090726
      TES(3184)%NYMD = 20090726
      TES(3185)%NYMD = 20090726
      TES(3186)%NYMD = 20090726
      TES(3187)%NYMD = 20090726
      TES(3188)%NYMD = 20090726
      TES(3189)%NYMD = 20090726
      TES(3190)%NYMD = 20090726
      TES(3191)%NYMD = 20090726
      TES(3192)%NYMD = 20090726
      TES(3193)%NYMD = 20090726
      TES(3194)%NYMD = 20090726
      TES(3195)%NYMD = 20090726
      TES(3196)%NYMD = 20090726
      TES(3197)%NYMD = 20090726
      TES(3198)%NYMD = 20090726
      TES(3199)%NYMD = 20090726
      TES(3200)%NYMD = 20090726
      TES(3201)%NYMD = 20090726
      TES(3202)%NYMD = 20090726
      TES(3203)%NYMD = 20090726
      TES(3204)%NYMD = 20090726
      TES(3205)%NYMD = 20090726
      TES(3206)%NYMD = 20090726
      TES(3207)%NYMD = 20090726
      TES(3208)%NYMD = 20090726
      TES(3209)%NYMD = 20090726
      TES(3210)%NYMD = 20090726
      TES(3211)%NYMD = 20090726
      TES(3212)%NYMD = 20090726
      TES(3213)%NYMD = 20090726
      TES(3214)%NYMD = 20090726
      TES(3215)%NYMD = 20090726
      TES(3216)%NYMD = 20090726
      TES(3217)%NYMD = 20090726
      TES(3218)%NYMD = 20090726
      TES(3219)%NYMD = 20090726
      TES(3220)%NYMD = 20090726
      TES(3221)%NYMD = 20090726
      TES(3222)%NYMD = 20090726
      TES(3223)%NYMD = 20090726
      TES(3224)%NYMD = 20090726
      TES(3225)%NYMD = 20090726
      TES(3226)%NYMD = 20090726
      TES(3227)%NYMD = 20090726
      TES(3228)%NYMD = 20090726
      TES(3229)%NYMD = 20090726
      TES(3230)%NYMD = 20090726
      TES(3231)%NYMD = 20090726
      TES(3232)%NYMD = 20090726
      TES(3233)%NYMD = 20090726
      TES(3234)%NYMD = 20090726
      TES(3235)%NYMD = 20090726
      TES(3236)%NYMD = 20090726
      TES(3237)%NYMD = 20090726
      TES(3238)%NYMD = 20090726
      TES(3239)%NYMD = 20090726
      TES(3240)%NYMD = 20090726
      TES(3241)%NYMD = 20090726
      TES(3242)%NYMD = 20090726
      TES(3243)%NYMD = 20090726
      TES(3244)%NYMD = 20090726
      TES(3245)%NYMD = 20090726
      TES(3246)%NYMD = 20090726
      TES(3247)%NYMD = 20090726
      TES(3248)%NYMD = 20090726
      TES(3249)%NYMD = 20090726
      TES(3250)%NYMD = 20090726
      TES(3251)%NYMD = 20090726
      TES(3252)%NYMD = 20090726
      TES(3253)%NYMD = 20090726
      TES(3254)%NYMD = 20090726
      TES(3255)%NYMD = 20090726
      TES(3256)%NYMD = 20090726
      TES(3257)%NYMD = 20090726
      TES(3258)%NYMD = 20090726
      TES(3259)%NYMD = 20090726
      TES(3260)%NYMD = 20090726
      TES(3261)%NYMD = 20090726
      TES(3262)%NYMD = 20090726
      TES(3263)%NYMD = 20090726
      TES(3264)%NYMD = 20090726
      TES(3265)%NYMD = 20090726
      TES(3266)%NYMD = 20090726
      TES(3267)%NYMD = 20090726
      TES(3268)%NYMD = 20090726
      TES(3269)%NYMD = 20090726
      TES(3270)%NYMD = 20090726
      TES(3271)%NYMD = 20090726
      TES(3272)%NYMD = 20090726
      TES(3273)%NYMD = 20090726
      TES(3274)%NYMD = 20090726
      TES(3275)%NYMD = 20090726
      TES(3276)%NYMD = 20090726
      TES(3277)%NYMD = 20090726
      TES(3278)%NYMD = 20090726
      TES(3279)%NYMD = 20090726
      TES(3280)%NYMD = 20090726
      TES(3281)%NYMD = 20090726
      TES(3282)%NYMD = 20090726
      TES(3283)%NYMD = 20090726
      TES(3284)%NYMD = 20090726
      TES(3285)%NYMD = 20090726
      TES(3286)%NYMD = 20090726
      TES(3287)%NYMD = 20090726
      TES(3288)%NYMD = 20090726
      TES(3289)%NYMD = 20090726
      TES(3290)%NYMD = 20090727
      TES(3291)%NYMD = 20090727
      TES(3292)%NYMD = 20090727
      TES(3293)%NYMD = 20090727
      TES(3294)%NYMD = 20090727
      TES(3295)%NYMD = 20090727
      TES(3296)%NYMD = 20090727
      TES(3297)%NYMD = 20090727
      TES(3298)%NYMD = 20090727
      TES(3299)%NYMD = 20090727
      TES(3300)%NYMD = 20090727
      TES(3301)%NYMD = 20090727
      TES(3302)%NYMD = 20090727
      TES(3303)%NYMD = 20090727
      TES(3304)%NYMD = 20090727
      TES(3305)%NYMD = 20090727
      TES(3306)%NYMD = 20090727
      TES(3307)%NYMD = 20090727
      TES(3308)%NYMD = 20090727
      TES(3309)%NYMD = 20090727
      TES(3310)%NYMD = 20090727
      TES(3311)%NYMD = 20090727
      TES(3312)%NYMD = 20090727
      TES(3313)%NYMD = 20090727
      TES(3314)%NYMD = 20090727
      TES(3315)%NYMD = 20090727
      TES(3316)%NYMD = 20090727
      TES(3317)%NYMD = 20090727
      TES(3318)%NYMD = 20090727
      TES(3319)%NYMD = 20090727
      TES(3320)%NYMD = 20090727
      TES(3321)%NYMD = 20090727
      TES(3322)%NYMD = 20090727
      TES(3323)%NYMD = 20090727
      TES(3324)%NYMD = 20090727
      TES(3325)%NYMD = 20090727
      TES(3326)%NYMD = 20090727
      TES(3327)%NYMD = 20090727
      TES(3328)%NYMD = 20090727
      TES(3329)%NYMD = 20090727
      TES(3330)%NYMD = 20090727
      TES(3331)%NYMD = 20090727
      TES(3332)%NYMD = 20090727
      TES(3333)%NYMD = 20090727
      TES(3334)%NYMD = 20090727
      TES(3335)%NYMD = 20090727
      TES(3336)%NYMD = 20090727
      TES(3337)%NYMD = 20090727
      TES(3338)%NYMD = 20090727
      TES(3339)%NYMD = 20090727
      TES(3340)%NYMD = 20090727
      TES(3341)%NYMD = 20090727
      TES(3342)%NYMD = 20090727
      TES(3343)%NYMD = 20090727
      TES(3344)%NYMD = 20090727
      TES(3345)%NYMD = 20090727
      TES(3346)%NYMD = 20090727
      TES(3347)%NYMD = 20090727
      TES(3348)%NYMD = 20090727
      TES(3349)%NYMD = 20090727
      TES(3350)%NYMD = 20090727
      TES(3351)%NYMD = 20090727
      TES(3352)%NYMD = 20090727
      TES(3353)%NYMD = 20090727
      TES(3354)%NYMD = 20090727
      TES(3355)%NYMD = 20090727
      TES(3356)%NYMD = 20090727
      TES(3357)%NYMD = 20090727
      TES(3358)%NYMD = 20090727
      TES(3359)%NYMD = 20090727
      TES(3360)%NYMD = 20090727
      TES(3361)%NYMD = 20090727
      TES(3362)%NYMD = 20090727
      TES(3363)%NYMD = 20090727
      TES(3364)%NYMD = 20090727
      TES(3365)%NYMD = 20090727
      TES(3366)%NYMD = 20090727
      TES(3367)%NYMD = 20090727
      TES(3368)%NYMD = 20090727
      TES(3369)%NYMD = 20090727
      TES(3370)%NYMD = 20090727
      TES(3371)%NYMD = 20090727
      TES(3372)%NYMD = 20090727
      TES(3373)%NYMD = 20090727
      TES(3374)%NYMD = 20090727
      TES(3375)%NYMD = 20090727
      TES(3376)%NYMD = 20090727
      TES(3377)%NYMD = 20090727
      TES(3378)%NYMD = 20090727
      TES(3379)%NYMD = 20090727
      TES(3380)%NYMD = 20090727
      TES(3381)%NYMD = 20090727
      TES(3382)%NYMD = 20090727
      TES(3383)%NYMD = 20090727
      TES(3384)%NYMD = 20090727
      TES(3385)%NYMD = 20090727
      TES(3386)%NYMD = 20090727
      TES(3387)%NYMD = 20090727
      TES(3388)%NYMD = 20090727
      TES(3389)%NYMD = 20090727
      TES(3390)%NYMD = 20090727
      TES(3391)%NYMD = 20090727
      TES(3392)%NYMD = 20090727
      TES(3393)%NYMD = 20090727
      TES(3394)%NYMD = 20090727
      TES(3395)%NYMD = 20090727
      TES(3396)%NYMD = 20090727
      TES(3397)%NYMD = 20090727
      TES(3398)%NYMD = 20090727
      TES(3399)%NYMD = 20090727
      TES(3400)%NYMD = 20090727
      TES(3401)%NYMD = 20090727
      TES(3402)%NYMD = 20090727
      TES(3403)%NYMD = 20090727
      TES(3404)%NYMD = 20090727
      TES(3405)%NYMD = 20090727
      TES(3406)%NYMD = 20090727
      TES(3407)%NYMD = 20090727
      TES(3408)%NYMD = 20090727
      TES(3409)%NYMD = 20090727
      TES(3410)%NYMD = 20090727
      TES(3411)%NYMD = 20090727
      TES(3412)%NYMD = 20090727
      TES(3413)%NYMD = 20090727
      TES(3414)%NYMD = 20090727
      TES(3415)%NYMD = 20090727
      TES(3416)%NYMD = 20090727
      TES(3417)%NYMD = 20090727
      TES(3418)%NYMD = 20090727
      TES(3419)%NYMD = 20090727
      TES(3420)%NYMD = 20090727
      TES(3421)%NYMD = 20090727
      TES(3422)%NYMD = 20090727
      TES(3423)%NYMD = 20090727
      TES(3424)%NYMD = 20090727
      TES(3425)%NYMD = 20090727
      TES(3426)%NYMD = 20090727
      TES(3427)%NYMD = 20090727
      TES(3428)%NYMD = 20090727
      TES(3429)%NYMD = 20090727
      TES(3430)%NYMD = 20090727
      TES(3431)%NYMD = 20090727
      TES(3432)%NYMD = 20090727
      TES(3433)%NYMD = 20090727
      TES(3434)%NYMD = 20090727
      TES(3435)%NYMD = 20090727
      TES(3436)%NYMD = 20090727
      TES(3437)%NYMD = 20090727
      TES(3438)%NYMD = 20090727
      TES(3439)%NYMD = 20090727
      TES(3440)%NYMD = 20090727
      TES(3441)%NYMD = 20090727
      TES(3442)%NYMD = 20090727
      TES(3443)%NYMD = 20090728
      TES(3444)%NYMD = 20090728
      TES(3445)%NYMD = 20090728
      TES(3446)%NYMD = 20090728
      TES(3447)%NYMD = 20090728
      TES(3448)%NYMD = 20090728
      TES(3449)%NYMD = 20090728
      TES(3450)%NYMD = 20090728
      TES(3451)%NYMD = 20090728
      TES(3452)%NYMD = 20090728
      TES(3453)%NYMD = 20090728
      TES(3454)%NYMD = 20090728
      TES(3455)%NYMD = 20090728
      TES(3456)%NYMD = 20090728
      TES(3457)%NYMD = 20090728
      TES(3458)%NYMD = 20090728
      TES(3459)%NYMD = 20090728
      TES(3460)%NYMD = 20090728
      TES(3461)%NYMD = 20090728
      TES(3462)%NYMD = 20090728
      TES(3463)%NYMD = 20090728
      TES(3464)%NYMD = 20090728
      TES(3465)%NYMD = 20090728
      TES(3466)%NYMD = 20090728
      TES(3467)%NYMD = 20090728
      TES(3468)%NYMD = 20090728
      TES(3469)%NYMD = 20090728
      TES(3470)%NYMD = 20090728
      TES(3471)%NYMD = 20090728
      TES(3472)%NYMD = 20090728
      TES(3473)%NYMD = 20090728
      TES(3474)%NYMD = 20090728
      TES(3475)%NYMD = 20090728
      TES(3476)%NYMD = 20090728
      TES(3477)%NYMD = 20090728
      TES(3478)%NYMD = 20090728
      TES(3479)%NYMD = 20090728
      TES(3480)%NYMD = 20090728
      TES(3481)%NYMD = 20090728
      TES(3482)%NYMD = 20090728
      TES(3483)%NYMD = 20090728
      TES(3484)%NYMD = 20090728
      TES(3485)%NYMD = 20090728
      TES(3486)%NYMD = 20090728
      TES(3487)%NYMD = 20090728
      TES(3488)%NYMD = 20090728
      TES(3489)%NYMD = 20090728
      TES(3490)%NYMD = 20090728
      TES(3491)%NYMD = 20090728
      TES(3492)%NYMD = 20090728
      TES(3493)%NYMD = 20090728
      TES(3494)%NYMD = 20090728
      TES(3495)%NYMD = 20090728
      TES(3496)%NYMD = 20090728
      TES(3497)%NYMD = 20090728
      TES(3498)%NYMD = 20090728
      TES(3499)%NYMD = 20090728
      TES(3500)%NYMD = 20090728
      TES(3501)%NYMD = 20090728
      TES(3502)%NYMD = 20090728
      TES(3503)%NYMD = 20090728
      TES(3504)%NYMD = 20090728
      TES(3505)%NYMD = 20090728
      TES(3506)%NYMD = 20090728
      TES(3507)%NYMD = 20090728
      TES(3508)%NYMD = 20090728
      TES(3509)%NYMD = 20090728
      TES(3510)%NYMD = 20090728
      TES(3511)%NYMD = 20090728
      TES(3512)%NYMD = 20090728
      TES(3513)%NYMD = 20090728
      TES(3514)%NYMD = 20090728
      TES(3515)%NYMD = 20090728
      TES(3516)%NYMD = 20090728
      TES(3517)%NYMD = 20090728
      TES(3518)%NYMD = 20090728
      TES(3519)%NYMD = 20090728
      TES(3520)%NYMD = 20090728
      TES(3521)%NYMD = 20090728
      TES(3522)%NYMD = 20090728
      TES(3523)%NYMD = 20090728
      TES(3524)%NYMD = 20090728
      TES(3525)%NYMD = 20090728
      TES(3526)%NYMD = 20090728
      TES(3527)%NYMD = 20090728
      TES(3528)%NYMD = 20090728
      TES(3529)%NYMD = 20090728
      TES(3530)%NYMD = 20090728
      TES(3531)%NYMD = 20090728
      TES(3532)%NYMD = 20090728
      TES(3533)%NYMD = 20090728
      TES(3534)%NYMD = 20090728
      TES(3535)%NYMD = 20090728
      TES(3536)%NYMD = 20090728
      TES(3537)%NYMD = 20090728
      TES(3538)%NYMD = 20090728
      TES(3539)%NYMD = 20090728
      TES(3540)%NYMD = 20090728
      TES(3541)%NYMD = 20090728
      TES(3542)%NYMD = 20090728
      TES(3543)%NYMD = 20090728
      TES(3544)%NYMD = 20090728
      TES(3545)%NYMD = 20090728
      TES(3546)%NYMD = 20090728
      TES(3547)%NYMD = 20090728
      TES(3548)%NYMD = 20090728
      TES(3549)%NYMD = 20090728
      TES(3550)%NYMD = 20090728
      TES(3551)%NYMD = 20090728
      TES(3552)%NYMD = 20090728
      TES(3553)%NYMD = 20090728
      TES(3554)%NYMD = 20090728
      TES(3555)%NYMD = 20090728
      TES(3556)%NYMD = 20090728
      TES(3557)%NYMD = 20090728
      TES(3558)%NYMD = 20090728
      TES(3559)%NYMD = 20090728
      TES(3560)%NYMD = 20090728
      TES(3561)%NYMD = 20090728
      TES(3562)%NYMD = 20090728
      TES(3563)%NYMD = 20090728
      TES(3564)%NYMD = 20090728
      TES(3565)%NYMD = 20090728
      TES(3566)%NYMD = 20090728
      TES(3567)%NYMD = 20090728
      TES(3568)%NYMD = 20090728
      TES(3569)%NYMD = 20090728
      TES(3570)%NYMD = 20090728
      TES(3571)%NYMD = 20090728
      TES(3572)%NYMD = 20090728
      TES(3573)%NYMD = 20090728
      TES(3574)%NYMD = 20090728
      TES(3575)%NYMD = 20090729
      TES(3576)%NYMD = 20090729
      TES(3577)%NYMD = 20090729
      TES(3578)%NYMD = 20090729
      TES(3579)%NYMD = 20090729
      TES(3580)%NYMD = 20090729
      TES(3581)%NYMD = 20090729
      TES(3582)%NYMD = 20090729
      TES(3583)%NYMD = 20090729
      TES(3584)%NYMD = 20090729
      TES(3585)%NYMD = 20090729
      TES(3586)%NYMD = 20090729
      TES(3587)%NYMD = 20090729
      TES(3588)%NYMD = 20090729
      TES(3589)%NYMD = 20090729
      TES(3590)%NYMD = 20090729
      TES(3591)%NYMD = 20090729
      TES(3592)%NYMD = 20090729
      TES(3593)%NYMD = 20090729
      TES(3594)%NYMD = 20090729
      TES(3595)%NYMD = 20090729
      TES(3596)%NYMD = 20090729
      TES(3597)%NYMD = 20090729
      TES(3598)%NYMD = 20090729
      TES(3599)%NYMD = 20090729
      TES(3600)%NYMD = 20090729
      TES(3601)%NYMD = 20090729
      TES(3602)%NYMD = 20090729
      TES(3603)%NYMD = 20090729
      TES(3604)%NYMD = 20090729
      TES(3605)%NYMD = 20090729
      TES(3606)%NYMD = 20090729
      TES(3607)%NYMD = 20090729
      TES(3608)%NYMD = 20090729
      TES(3609)%NYMD = 20090729
      TES(3610)%NYMD = 20090729
      TES(3611)%NYMD = 20090729
      TES(3612)%NYMD = 20090729
      TES(3613)%NYMD = 20090729
      TES(3614)%NYMD = 20090729
      TES(3615)%NYMD = 20090729
      TES(3616)%NYMD = 20090729
      TES(3617)%NYMD = 20090729
      TES(3618)%NYMD = 20090729
      TES(3619)%NYMD = 20090729
      TES(3620)%NYMD = 20090729
      TES(3621)%NYMD = 20090729
      TES(3622)%NYMD = 20090729
      TES(3623)%NYMD = 20090729
      TES(3624)%NYMD = 20090729
      TES(3625)%NYMD = 20090729
      TES(3626)%NYMD = 20090729
      TES(3627)%NYMD = 20090729
      TES(3628)%NYMD = 20090729
      TES(3629)%NYMD = 20090729
      TES(3630)%NYMD = 20090729
      TES(3631)%NYMD = 20090729
      TES(3632)%NYMD = 20090729
      TES(3633)%NYMD = 20090729
      TES(3634)%NYMD = 20090729
      TES(3635)%NYMD = 20090729
      TES(3636)%NYMD = 20090729
      TES(3637)%NYMD = 20090729
      TES(3638)%NYMD = 20090729
      TES(3639)%NYMD = 20090729
      TES(3640)%NYMD = 20090729
      TES(3641)%NYMD = 20090729
      TES(3642)%NYMD = 20090729
      TES(3643)%NYMD = 20090729
      TES(3644)%NYMD = 20090729
      TES(3645)%NYMD = 20090729
      TES(3646)%NYMD = 20090729
      TES(3647)%NYMD = 20090729
      TES(3648)%NYMD = 20090729
      TES(3649)%NYMD = 20090729
      TES(3650)%NYMD = 20090729
      TES(3651)%NYMD = 20090729
      TES(3652)%NYMD = 20090729
      TES(3653)%NYMD = 20090729
      TES(3654)%NYMD = 20090729
      TES(3655)%NYMD = 20090729
      TES(3656)%NYMD = 20090729
      TES(3657)%NYMD = 20090729
      TES(3658)%NYMD = 20090729
      TES(3659)%NYMD = 20090729
      TES(3660)%NYMD = 20090729
      TES(3661)%NYMD = 20090729
      TES(3662)%NYMD = 20090729
      TES(3663)%NYMD = 20090729
      TES(3664)%NYMD = 20090729
      TES(3665)%NYMD = 20090729
      TES(3666)%NYMD = 20090729
      TES(3667)%NYMD = 20090729
      TES(3668)%NYMD = 20090729
      TES(3669)%NYMD = 20090729
      TES(3670)%NYMD = 20090729
      TES(3671)%NYMD = 20090729
      TES(3672)%NYMD = 20090729
      TES(3673)%NYMD = 20090729
      TES(3674)%NYMD = 20090729
      TES(3675)%NYMD = 20090729
      TES(3676)%NYMD = 20090729
      TES(3677)%NYMD = 20090729
      TES(3678)%NYMD = 20090729
      TES(3679)%NYMD = 20090729
      TES(3680)%NYMD = 20090729
      TES(3681)%NYMD = 20090729
      TES(3682)%NYMD = 20090729
      TES(3683)%NYMD = 20090729
      TES(3684)%NYMD = 20090729
      TES(3685)%NYMD = 20090729
      TES(3686)%NYMD = 20090729
      TES(3687)%NYMD = 20090729
      TES(3688)%NYMD = 20090729
      TES(3689)%NYMD = 20090729
      TES(3690)%NYMD = 20090729
      TES(3691)%NYMD = 20090729
      TES(3692)%NYMD = 20090729
      TES(3693)%NYMD = 20090729
      TES(3694)%NYMD = 20090729
      TES(3695)%NYMD = 20090729
      TES(3696)%NYMD = 20090729
      TES(3697)%NYMD = 20090729
      TES(3698)%NYMD = 20090729
      TES(3699)%NYMD = 20090729
      TES(3700)%NYMD = 20090729
      TES(3701)%NYMD = 20090729
      TES(3702)%NYMD = 20090729
      TES(3703)%NYMD = 20090729
      TES(3704)%NYMD = 20090729
      TES(3705)%NYMD = 20090730
      TES(3706)%NYMD = 20090730
      TES(3707)%NYMD = 20090730
      TES(3708)%NYMD = 20090730
      TES(3709)%NYMD = 20090730
      TES(3710)%NYMD = 20090730
      TES(3711)%NYMD = 20090730
      TES(3712)%NYMD = 20090730
      TES(3713)%NYMD = 20090730
      TES(3714)%NYMD = 20090730
      TES(3715)%NYMD = 20090730
      TES(3716)%NYMD = 20090730
      TES(3717)%NYMD = 20090730
      TES(3718)%NYMD = 20090730
      TES(3719)%NYMD = 20090730
      TES(3720)%NYMD = 20090730
      TES(3721)%NYMD = 20090730
      TES(3722)%NYMD = 20090730
      TES(3723)%NYMD = 20090730
      TES(3724)%NYMD = 20090730
      TES(3725)%NYMD = 20090730
      TES(3726)%NYMD = 20090730
      TES(3727)%NYMD = 20090730
      TES(3728)%NYMD = 20090730
      TES(3729)%NYMD = 20090730
      TES(3730)%NYMD = 20090730
      TES(3731)%NYMD = 20090730
      TES(3732)%NYMD = 20090730
      TES(3733)%NYMD = 20090730
      TES(3734)%NYMD = 20090730
      TES(3735)%NYMD = 20090730
      TES(3736)%NYMD = 20090730
      TES(3737)%NYMD = 20090730
      TES(3738)%NYMD = 20090730
      TES(3739)%NYMD = 20090730
      TES(3740)%NYMD = 20090730
      TES(3741)%NYMD = 20090730
      TES(3742)%NYMD = 20090730
      TES(3743)%NYMD = 20090730
      TES(3744)%NYMD = 20090730
      TES(3745)%NYMD = 20090730
      TES(3746)%NYMD = 20090730
      TES(3747)%NYMD = 20090730
      TES(3748)%NYMD = 20090730
      TES(3749)%NYMD = 20090730
      TES(3750)%NYMD = 20090730
      TES(3751)%NYMD = 20090730
      TES(3752)%NYMD = 20090730
      TES(3753)%NYMD = 20090730
      TES(3754)%NYMD = 20090730
      TES(3755)%NYMD = 20090730
      TES(3756)%NYMD = 20090730
      TES(3757)%NYMD = 20090730
      TES(3758)%NYMD = 20090730
      TES(3759)%NYMD = 20090730
      TES(3760)%NYMD = 20090730
      TES(3761)%NYMD = 20090730
      TES(3762)%NYMD = 20090730
      TES(3763)%NYMD = 20090730
      TES(3764)%NYMD = 20090730
      TES(3765)%NYMD = 20090730
      TES(3766)%NYMD = 20090730
      TES(3767)%NYMD = 20090730
      TES(3768)%NYMD = 20090730
      TES(3769)%NYMD = 20090730
      TES(3770)%NYMD = 20090730
      TES(3771)%NYMD = 20090730
      TES(3772)%NYMD = 20090730
      TES(3773)%NYMD = 20090730
      TES(3774)%NYMD = 20090730
      TES(3775)%NYMD = 20090730
      TES(3776)%NYMD = 20090730
      TES(3777)%NYMD = 20090730
      TES(3778)%NYMD = 20090730
      TES(3779)%NYMD = 20090730
      TES(3780)%NYMD = 20090730
      TES(3781)%NYMD = 20090730
      TES(3782)%NYMD = 20090730
      TES(3783)%NYMD = 20090730
      TES(3784)%NYMD = 20090730
      TES(3785)%NYMD = 20090730
      TES(3786)%NYMD = 20090730
      TES(3787)%NYMD = 20090730
      TES(3788)%NYMD = 20090730
      TES(3789)%NYMD = 20090730
      TES(3790)%NYMD = 20090730
      TES(3791)%NYMD = 20090730
      TES(3792)%NYMD = 20090730
      TES(3793)%NYMD = 20090730
      TES(3794)%NYMD = 20090730
      TES(3795)%NYMD = 20090730
      TES(3796)%NYMD = 20090730
      TES(3797)%NYMD = 20090730
      TES(3798)%NYMD = 20090730
      TES(3799)%NYMD = 20090730
      TES(3800)%NYMD = 20090730
      TES(3801)%NYMD = 20090730
      TES(3802)%NYMD = 20090730
      TES(3803)%NYMD = 20090730
      TES(3804)%NYMD = 20090730
      TES(3805)%NYMD = 20090730
      TES(3806)%NYMD = 20090730
      TES(3807)%NYMD = 20090730
      TES(3808)%NYMD = 20090730
      TES(3809)%NYMD = 20090730
      TES(3810)%NYMD = 20090730
      TES(3811)%NYMD = 20090730
      TES(3812)%NYMD = 20090730
      TES(3813)%NYMD = 20090730
      TES(3814)%NYMD = 20090730
      TES(3815)%NYMD = 20090730
      TES(3816)%NYMD = 20090730
      TES(3817)%NYMD = 20090730
      TES(3818)%NYMD = 20090730
      TES(3819)%NYMD = 20090730
      TES(3820)%NYMD = 20090730
      TES(3821)%NYMD = 20090730
      TES(3822)%NYMD = 20090730
      TES(3823)%NYMD = 20090730
      TES(3824)%NYMD = 20090730
      TES(3825)%NYMD = 20090730
      TES(3826)%NYMD = 20090730
      TES(3827)%NYMD = 20090730
      TES(3828)%NYMD = 20090730
      TES(3829)%NYMD = 20090730
      TES(3830)%NYMD = 20090730
      TES(3831)%NYMD = 20090730
      TES(3832)%NYMD = 20090730
      TES(3833)%NYMD = 20090730
      TES(3834)%NYMD = 20090730
      TES(3835)%NYMD = 20090730
      TES(3836)%NYMD = 20090730
      TES(3837)%NYMD = 20090730
      TES(3838)%NYMD = 20090730
      TES(3839)%NYMD = 20090730
      TES(3840)%NYMD = 20090730
      TES(3841)%NYMD = 20090731
      TES(3842)%NYMD = 20090731
      TES(3843)%NYMD = 20090731
      TES(3844)%NYMD = 20090731
      TES(3845)%NYMD = 20090731
      TES(3846)%NYMD = 20090731
      TES(3847)%NYMD = 20090731
      TES(3848)%NYMD = 20090731
      TES(3849)%NYMD = 20090731
      TES(3850)%NYMD = 20090731
      TES(3851)%NYMD = 20090731
      TES(3852)%NYMD = 20090731
      TES(3853)%NYMD = 20090731
      TES(3854)%NYMD = 20090731
      TES(3855)%NYMD = 20090731
      TES(3856)%NYMD = 20090731
      TES(3857)%NYMD = 20090731
      TES(3858)%NYMD = 20090731
      TES(3859)%NYMD = 20090731
      TES(3860)%NYMD = 20090731
      TES(3861)%NYMD = 20090731
      TES(3862)%NYMD = 20090731
      TES(3863)%NYMD = 20090731
      TES(3864)%NYMD = 20090731
      TES(3865)%NYMD = 20090731
      TES(3866)%NYMD = 20090731
      TES(3867)%NYMD = 20090731
      TES(3868)%NYMD = 20090731
      TES(3869)%NYMD = 20090731
      TES(3870)%NYMD = 20090731
      TES(3871)%NYMD = 20090731
      TES(3872)%NYMD = 20090731
      TES(3873)%NYMD = 20090731
      TES(3874)%NYMD = 20090731
      TES(3875)%NYMD = 20090731
      TES(3876)%NYMD = 20090731
      TES(3877)%NYMD = 20090731
      TES(3878)%NYMD = 20090731
      TES(3879)%NYMD = 20090731
      TES(3880)%NYMD = 20090731
      TES(3881)%NYMD = 20090731
      TES(3882)%NYMD = 20090731
      TES(3883)%NYMD = 20090731
      TES(3884)%NYMD = 20090731
      TES(3885)%NYMD = 20090731
      TES(3886)%NYMD = 20090731
      TES(3887)%NYMD = 20090731
      TES(3888)%NYMD = 20090731
      TES(3889)%NYMD = 20090731
      TES(3890)%NYMD = 20090731
      TES(3891)%NYMD = 20090731
      TES(3892)%NYMD = 20090731
      TES(3893)%NYMD = 20090731
      TES(3894)%NYMD = 20090731
      TES(3895)%NYMD = 20090731
      TES(3896)%NYMD = 20090731
      TES(3897)%NYMD = 20090731
      TES(3898)%NYMD = 20090731
      TES(3899)%NYMD = 20090731
      TES(3900)%NYMD = 20090731
      TES(3901)%NYMD = 20090731
      TES(3902)%NYMD = 20090731
      TES(3903)%NYMD = 20090731
      TES(3904)%NYMD = 20090731
      TES(3905)%NYMD = 20090731
      TES(3906)%NYMD = 20090731
      TES(3907)%NYMD = 20090731
      TES(3908)%NYMD = 20090731
      TES(3909)%NYMD = 20090731
      TES(3910)%NYMD = 20090731
      TES(3911)%NYMD = 20090731
      TES(3912)%NYMD = 20090731
      TES(3913)%NYMD = 20090731
      TES(3914)%NYMD = 20090731
      TES(3915)%NYMD = 20090731
      TES(3916)%NYMD = 20090731
      TES(3917)%NYMD = 20090731
      TES(3918)%NYMD = 20090731
      TES(3919)%NYMD = 20090731
      TES(3920)%NYMD = 20090731
      TES(3921)%NYMD = 20090731
      TES(3922)%NYMD = 20090731
      TES(3923)%NYMD = 20090731
      TES(3924)%NYMD = 20090731
      TES(3925)%NYMD = 20090731
      TES(3926)%NYMD = 20090731
      TES(3927)%NYMD = 20090731
      TES(3928)%NYMD = 20090731
      TES(3929)%NYMD = 20090731
      TES(3930)%NYMD = 20090731
      TES(3931)%NYMD = 20090731
      TES(3932)%NYMD = 20090731
      TES(3933)%NYMD = 20090731
      TES(3934)%NYMD = 20090731
      TES(3935)%NYMD = 20090731
      TES(3936)%NYMD = 20090731
      TES(3937)%NYMD = 20090731
      TES(3938)%NYMD = 20090731
      TES(3939)%NYMD = 20090731
      TES(3940)%NYMD = 20090731
      TES(3941)%NYMD = 20090731
      TES(3942)%NYMD = 20090731
      TES(3943)%NYMD = 20090731
      TES(3944)%NYMD = 20090731
      TES(3945)%NYMD = 20090731
      TES(3946)%NYMD = 20090731
      TES(3947)%NYMD = 20090731
      TES(3948)%NYMD = 20090731
      TES(3949)%NYMD = 20090731
      TES(3950)%NYMD = 20090731
      TES(3951)%NYMD = 20090731
      TES(3952)%NYMD = 20090731
      TES(3953)%NYMD = 20090731
      TES(3954)%NYMD = 20090731
      TES(3955)%NYMD = 20090731
      TES(3956)%NYMD = 20090731
      TES(3957)%NYMD = 20090731
      TES(3958)%NYMD = 20090731
      TES(3959)%NYMD = 20090731
      TES(3960)%NYMD = 20090731
      TES(3961)%NYMD = 20090731
      TES(3962)%NYMD = 20090731
      TES(3963)%NYMD = 20090731
      TES(3964)%NYMD = 20090731
      TES(3965)%NYMD = 20090731
      TES(3966)%NYMD = 20090731
      TES(3967)%NYMD = 20090731
      TES(3968)%NYMD = 20090731
      TES(3969)%NYMD = 20090731
      TES(3970)%NYMD = 20090731
      TES(3971)%NYMD = 20090731
      TES(3972)%NYMD = 20090731
      TES(3973)%NYMD = 20090731
      TES(3974)%NYMD = 20090731
      TES(3975)%NYMD = 20090731
      TES(3976)%NYMD = 20090731
      TES(3977)%NYMD = 20090731
      TES(3978)%NYMD = 20090731
      TES(3979)%NYMD = 20090731
      TES(3980)%NYMD = 20090731
      TES(3981)%NYMD = 20090731
      TES(3982)%NYMD = 20090731
      TES(3983)%NYMD = 20090731
      TES(3984)%NYMD = 20090731
      TES(3985)%NYMD = 20090731
      TES(3986)%NYMD = 20090731
      TES(3987)%NYMD = 20090731
      TES(3988)%NYMD = 20090731
      TES(3989)%NYMD = 20090731
      TES(3990)%NYMD = 20090731
      TES(3991)%NYMD = 20090731
      
      
      TES(1)%NHMS = 101700
      TES(2)%NHMS = 101700
      TES(3)%NHMS = 101800
      TES(4)%NHMS = 101900
      TES(5)%NHMS = 102000
      TES(6)%NHMS = 102100
      TES(7)%NHMS = 102100
      TES(8)%NHMS = 102200
      TES(9)%NHMS = 102200
      TES(10)%NHMS = 102300
      TES(11)%NHMS = 102300
      TES(12)%NHMS = 104100
      TES(13)%NHMS = 104100
      TES(14)%NHMS = 104200
      TES(15)%NHMS = 114000
      TES(16)%NHMS = 114300
      TES(17)%NHMS = 114300
      TES(18)%NHMS = 114400
      TES(19)%NHMS = 114700
      TES(20)%NHMS = 114800
      TES(21)%NHMS = 115100
      TES(22)%NHMS = 115300
      TES(23)%NHMS = 115500
      TES(24)%NHMS = 115600
      TES(25)%NHMS = 115700
      TES(26)%NHMS = 115700
      TES(27)%NHMS = 115800
      TES(28)%NHMS = 115900
      TES(29)%NHMS = 115900
      TES(30)%NHMS = 120000
      TES(31)%NHMS = 120100
      TES(32)%NHMS = 132600
      TES(33)%NHMS = 132700
      TES(34)%NHMS = 132700
      TES(35)%NHMS = 155500
      TES(36)%NHMS = 155800
      TES(37)%NHMS = 155800
      TES(38)%NHMS = 160200
      TES(39)%NHMS = 160200
      TES(40)%NHMS = 163800
      TES(41)%NHMS = 163800
      TES(42)%NHMS = 163800
      TES(43)%NHMS = 163900
      TES(44)%NHMS = 163900
      TES(45)%NHMS = 164000
      TES(46)%NHMS = 164000
      TES(47)%NHMS = 165500
      TES(48)%NHMS = 165500
      TES(49)%NHMS = 165600
      TES(50)%NHMS = 165600
      TES(51)%NHMS = 165700
      TES(52)%NHMS = 171600
      TES(53)%NHMS = 171700
      TES(54)%NHMS = 171700
      TES(55)%NHMS = 172200
      TES(56)%NHMS = 172900
      TES(57)%NHMS = 172900
      TES(58)%NHMS = 172900
      TES(59)%NHMS = 173700
      TES(60)%NHMS = 174100
      TES(61)%NHMS = 181000
      TES(62)%NHMS = 181100
      TES(63)%NHMS = 181200
      TES(64)%NHMS = 181200
      TES(65)%NHMS = 181400
      TES(66)%NHMS = 181400
      TES(67)%NHMS = 181500
      TES(68)%NHMS = 181500
      TES(69)%NHMS = 181600
      TES(70)%NHMS = 181700
      TES(71)%NHMS = 181700
      TES(72)%NHMS = 181800
      TES(73)%NHMS = 181800
      TES(74)%NHMS = 181900
      TES(75)%NHMS = 181900
      TES(76)%NHMS = 182000
      TES(77)%NHMS = 182000
      TES(78)%NHMS = 182000
      TES(79)%NHMS = 182100
      TES(80)%NHMS = 182300
      TES(81)%NHMS = 182400
      TES(82)%NHMS = 182600
      TES(83)%NHMS = 183000
      TES(84)%NHMS = 183100
      TES(85)%NHMS = 183100
      TES(86)%NHMS = 183200
      TES(87)%NHMS = 183300
      TES(88)%NHMS = 183500
      TES(89)%NHMS = 183500
      TES(90)%NHMS = 183600
      TES(91)%NHMS = 185400
      TES(92)%NHMS = 185400
      TES(93)%NHMS = 185500
      TES(94)%NHMS = 185500
      TES(95)%NHMS = 185500
      TES(96)%NHMS = 185600
      TES(97)%NHMS = 185600
      TES(98)%NHMS = 185700
      TES(99)%NHMS = 190400
      TES(100)%NHMS = 190500
      TES(101)%NHMS = 190500
      TES(102)%NHMS = 190500
      TES(103)%NHMS = 190800
      TES(104)%NHMS = 190900
      TES(105)%NHMS = 200400
      TES(106)%NHMS = 200500
      TES(107)%NHMS = 200600
      TES(108)%NHMS = 200700
      TES(109)%NHMS = 200800
      TES(110)%NHMS = 200900
      TES(111)%NHMS = 200900
      TES(112)%NHMS = 201000
      TES(113)%NHMS = 201200
      TES(114)%NHMS = 201300
      TES(115)%NHMS = 201400
      TES(116)%NHMS = 201500
      TES(117)%NHMS = 201600
      TES(118)%NHMS = 203200
      TES(119)%NHMS = 203300
      TES(120)%NHMS = 203300
      TES(121)%NHMS = 203500
      TES(122)%NHMS = 203600
      TES(123)%NHMS = 203600
      TES(124)%NHMS = 203600
      TES(125)%NHMS = 203700
      TES(126)%NHMS = 203700
      TES(127)%NHMS = 204100
      TES(128)%NHMS = 204200
      TES(129)%NHMS = 204500
      TES(130)%NHMS = 215300
      TES(131)%NHMS = 221300
      TES(132)%NHMS = 221300
      TES(133)%NHMS = 221800
      TES(134)%NHMS = 223300
      TES(135)%NHMS = 223300
      TES(136)%NHMS = 223300
      TES(137)%NHMS = 223400
      TES(138)%NHMS = 233300
      TES(139)%NHMS = 233400
      TES(140)%NHMS = 235000
      TES(141)%NHMS = 235100
      TES(142)%NHMS = 235100
      TES(143)%NHMS = 235200
      TES(144)%NHMS = 235300
      TES(145)%NHMS = 235400
      TES(146)%NHMS = 235500
      TES(147)%NHMS = 235600
      TES(148)%NHMS = 235600
      TES(149)%NHMS = 000100
      TES(150)%NHMS = 000100
      TES(151)%NHMS = 000200
      TES(152)%NHMS = 000300
      TES(153)%NHMS = 000500
      TES(154)%NHMS = 000500
      TES(155)%NHMS = 000600
      TES(156)%NHMS = 000700
      TES(157)%NHMS = 000800
      TES(158)%NHMS = 000800
      TES(159)%NHMS = 000900
      TES(160)%NHMS = 000900
      TES(161)%NHMS = 001000
      TES(162)%NHMS = 001000
      TES(163)%NHMS = 001100
      TES(164)%NHMS = 012900
      TES(165)%NHMS = 013200
      TES(166)%NHMS = 013300
      TES(167)%NHMS = 013400
      TES(168)%NHMS = 013600
      TES(169)%NHMS = 014200
      TES(170)%NHMS = 014300
      TES(171)%NHMS = 014300
      TES(172)%NHMS = 014400
      TES(173)%NHMS = 030800
      TES(174)%NHMS = 031000
      TES(175)%NHMS = 031700
      TES(176)%NHMS = 040500
      TES(177)%NHMS = 040500
      TES(178)%NHMS = 040600
      TES(179)%NHMS = 040700
      TES(180)%NHMS = 040700
      TES(181)%NHMS = 040800
      TES(182)%NHMS = 042600
      TES(183)%NHMS = 042700
      TES(184)%NHMS = 042800
      TES(185)%NHMS = 042800
      TES(186)%NHMS = 042900
      TES(187)%NHMS = 050400
      TES(188)%NHMS = 050400
      TES(189)%NHMS = 050400
      TES(190)%NHMS = 050500
      TES(191)%NHMS = 050500
      TES(192)%NHMS = 050600
      TES(193)%NHMS = 050700
      TES(194)%NHMS = 050700
      TES(195)%NHMS = 050800
      TES(196)%NHMS = 051000
      TES(197)%NHMS = 051100
      TES(198)%NHMS = 051200
      TES(199)%NHMS = 051200
      TES(200)%NHMS = 054800
      TES(201)%NHMS = 055200
      TES(202)%NHMS = 055300
      TES(203)%NHMS = 055400
      TES(204)%NHMS = 055500
      TES(205)%NHMS = 060600
      TES(206)%NHMS = 060600
      TES(207)%NHMS = 060900
      TES(208)%NHMS = 060900
      TES(209)%NHMS = 061000
      TES(210)%NHMS = 062700
      TES(211)%NHMS = 062800
      TES(212)%NHMS = 062800
      TES(213)%NHMS = 063000
      TES(214)%NHMS = 064100
      TES(215)%NHMS = 064300
      TES(216)%NHMS = 064400
      TES(217)%NHMS = 064400
      TES(218)%NHMS = 064500
      TES(219)%NHMS = 064500
      TES(220)%NHMS = 073900
      TES(221)%NHMS = 074000
      TES(222)%NHMS = 074100
      TES(223)%NHMS = 074300
      TES(224)%NHMS = 074400
      TES(225)%NHMS = 074500
      TES(226)%NHMS = 074600
      TES(227)%NHMS = 074600
      TES(228)%NHMS = 074800
      TES(229)%NHMS = 080700
      TES(230)%NHMS = 080700
      TES(231)%NHMS = 081000
      TES(232)%NHMS = 081000
      TES(233)%NHMS = 081100
      TES(234)%NHMS = 081100
      TES(235)%NHMS = 081200
      TES(236)%NHMS = 081200
      TES(237)%NHMS = 081700
      TES(238)%NHMS = 092000
      TES(239)%NHMS = 092100
      TES(240)%NHMS = 092300
      TES(241)%NHMS = 092300
      TES(242)%NHMS = 092400
      TES(243)%NHMS = 092400
      TES(244)%NHMS = 092500
      TES(245)%NHMS = 092500
      TES(246)%NHMS = 092600
      TES(247)%NHMS = 092700
      TES(248)%NHMS = 092700
      TES(249)%NHMS = 094400
      TES(250)%NHMS = 094500
      TES(251)%NHMS = 094600
      TES(252)%NHMS = 094600
      TES(253)%NHMS = 094700
      TES(254)%NHMS = 094700
      TES(255)%NHMS = 094800
      TES(256)%NHMS = 094800
      TES(257)%NHMS = 094800
      TES(258)%NHMS = 094900
      TES(259)%NHMS = 095000
      TES(260)%NHMS = 095100
      TES(261)%NHMS = 095100
      TES(262)%NHMS = 104500
      TES(263)%NHMS = 105000
      TES(264)%NHMS = 105100
      TES(265)%NHMS = 105200
      TES(266)%NHMS = 105900
      TES(267)%NHMS = 110000
      TES(268)%NHMS = 110100
      TES(269)%NHMS = 110200
      TES(270)%NHMS = 110600
      TES(271)%NHMS = 112300
      TES(272)%NHMS = 112300
      TES(273)%NHMS = 114400
      TES(274)%NHMS = 100400
      TES(275)%NHMS = 100500
      TES(276)%NHMS = 100600
      TES(277)%NHMS = 100800
      TES(278)%NHMS = 100900
      TES(279)%NHMS = 100900
      TES(280)%NHMS = 101000
      TES(281)%NHMS = 101100
      TES(282)%NHMS = 102700
      TES(283)%NHMS = 102800
      TES(284)%NHMS = 102900
      TES(285)%NHMS = 103000
      TES(286)%NHMS = 113300
      TES(287)%NHMS = 113300
      TES(288)%NHMS = 113300
      TES(289)%NHMS = 113800
      TES(290)%NHMS = 113900
      TES(291)%NHMS = 114000
      TES(292)%NHMS = 114400
      TES(293)%NHMS = 114400
      TES(294)%NHMS = 114500
      TES(295)%NHMS = 114600
      TES(296)%NHMS = 114600
      TES(297)%NHMS = 114700
      TES(298)%NHMS = 114800
      TES(299)%NHMS = 131500
      TES(300)%NHMS = 132400
      TES(301)%NHMS = 132700
      TES(302)%NHMS = 134600
      TES(303)%NHMS = 154300
      TES(304)%NHMS = 154600
      TES(305)%NHMS = 154600
      TES(306)%NHMS = 154700
      TES(307)%NHMS = 154900
      TES(308)%NHMS = 155000
      TES(309)%NHMS = 162600
      TES(310)%NHMS = 162700
      TES(311)%NHMS = 162700
      TES(312)%NHMS = 162800
      TES(313)%NHMS = 162800
      TES(314)%NHMS = 164300
      TES(315)%NHMS = 170500
      TES(316)%NHMS = 171900
      TES(317)%NHMS = 172000
      TES(318)%NHMS = 172300
      TES(319)%NHMS = 172800
      TES(320)%NHMS = 172800
      TES(321)%NHMS = 175900
      TES(322)%NHMS = 175900
      TES(323)%NHMS = 180000
      TES(324)%NHMS = 180000
      TES(325)%NHMS = 180100
      TES(326)%NHMS = 180200
      TES(327)%NHMS = 180200
      TES(328)%NHMS = 180400
      TES(329)%NHMS = 180700
      TES(330)%NHMS = 180800
      TES(331)%NHMS = 180800
      TES(332)%NHMS = 180900
      TES(333)%NHMS = 181000
      TES(334)%NHMS = 181100
      TES(335)%NHMS = 181900
      TES(336)%NHMS = 182000
      TES(337)%NHMS = 182100
      TES(338)%NHMS = 182300
      TES(339)%NHMS = 184300
      TES(340)%NHMS = 184300
      TES(341)%NHMS = 184800
      TES(342)%NHMS = 184900
      TES(343)%NHMS = 185000
      TES(344)%NHMS = 185400
      TES(345)%NHMS = 185500
      TES(346)%NHMS = 185800
      TES(347)%NHMS = 185900
      TES(348)%NHMS = 195200
      TES(349)%NHMS = 195500
      TES(350)%NHMS = 195600
      TES(351)%NHMS = 195600
      TES(352)%NHMS = 195700
      TES(353)%NHMS = 195800
      TES(354)%NHMS = 200100
      TES(355)%NHMS = 200300
      TES(356)%NHMS = 202000
      TES(357)%NHMS = 202000
      TES(358)%NHMS = 202100
      TES(359)%NHMS = 202100
      TES(360)%NHMS = 202200
      TES(361)%NHMS = 202200
      TES(362)%NHMS = 202400
      TES(363)%NHMS = 202500
      TES(364)%NHMS = 202600
      TES(365)%NHMS = 203000
      TES(366)%NHMS = 203000
      TES(367)%NHMS = 203100
      TES(368)%NHMS = 214100
      TES(369)%NHMS = 214100
      TES(370)%NHMS = 214200
      TES(371)%NHMS = 214200
      TES(372)%NHMS = 214300
      TES(373)%NHMS = 220000
      TES(374)%NHMS = 220100
      TES(375)%NHMS = 220300
      TES(376)%NHMS = 220600
      TES(377)%NHMS = 220600
      TES(378)%NHMS = 220700
      TES(379)%NHMS = 222000
      TES(380)%NHMS = 233900
      TES(381)%NHMS = 234000
      TES(382)%NHMS = 234000
      TES(383)%NHMS = 234100
      TES(384)%NHMS = 234100
      TES(385)%NHMS = 234200
      TES(386)%NHMS = 234300
      TES(387)%NHMS = 234300
      TES(388)%NHMS = 234400
      TES(389)%NHMS = 234400
      TES(390)%NHMS = 234500
      TES(391)%NHMS = 234500
      TES(392)%NHMS = 234700
      TES(393)%NHMS = 234900
      TES(394)%NHMS = 234900
      TES(395)%NHMS = 235100
      TES(396)%NHMS = 235200
      TES(397)%NHMS = 235200
      TES(398)%NHMS = 235300
      TES(399)%NHMS = 235300
      TES(400)%NHMS = 235400
      TES(401)%NHMS = 235400
      TES(402)%NHMS = 235500
      TES(403)%NHMS = 235500
      TES(404)%NHMS = 235600
      TES(405)%NHMS = 235700
      TES(406)%NHMS = 235700
      TES(407)%NHMS = 235800
      TES(408)%NHMS = 235800
      TES(409)%NHMS = 000400
      TES(410)%NHMS = 000400
      TES(411)%NHMS = 012000
      TES(412)%NHMS = 012000
      TES(413)%NHMS = 012200
      TES(414)%NHMS = 012400
      TES(415)%NHMS = 012500
      TES(416)%NHMS = 012800
      TES(417)%NHMS = 013000
      TES(418)%NHMS = 022500
      TES(419)%NHMS = 023700
      TES(420)%NHMS = 023800
      TES(421)%NHMS = 025600
      TES(422)%NHMS = 025600
      TES(423)%NHMS = 035300
      TES(424)%NHMS = 035400
      TES(425)%NHMS = 035400
      TES(426)%NHMS = 035500
      TES(427)%NHMS = 035900
      TES(428)%NHMS = 041500
      TES(429)%NHMS = 041500
      TES(430)%NHMS = 041600
      TES(431)%NHMS = 041600
      TES(432)%NHMS = 041700
      TES(433)%NHMS = 041700
      TES(434)%NHMS = 041700
      TES(435)%NHMS = 041800
      TES(436)%NHMS = 045200
      TES(437)%NHMS = 045200
      TES(438)%NHMS = 045300
      TES(439)%NHMS = 045400
      TES(440)%NHMS = 045600
      TES(441)%NHMS = 045700
      TES(442)%NHMS = 045800
      TES(443)%NHMS = 045800
      TES(444)%NHMS = 045900
      TES(445)%NHMS = 045900
      TES(446)%NHMS = 050000
      TES(447)%NHMS = 053900
      TES(448)%NHMS = 054400
      TES(449)%NHMS = 054800
      TES(450)%NHMS = 054900
      TES(451)%NHMS = 055000
      TES(452)%NHMS = 055200
      TES(453)%NHMS = 055300
      TES(454)%NHMS = 055400
      TES(455)%NHMS = 055500
      TES(456)%NHMS = 055500
      TES(457)%NHMS = 055600
      TES(458)%NHMS = 055600
      TES(459)%NHMS = 055700
      TES(460)%NHMS = 055700
      TES(461)%NHMS = 062500
      TES(462)%NHMS = 062700
      TES(463)%NHMS = 063100
      TES(464)%NHMS = 063200
      TES(465)%NHMS = 063200
      TES(466)%NHMS = 063400
      TES(467)%NHMS = 063400
      TES(468)%NHMS = 072900
      TES(469)%NHMS = 073000
      TES(470)%NHMS = 073400
      TES(471)%NHMS = 073400
      TES(472)%NHMS = 073500
      TES(473)%NHMS = 073500
      TES(474)%NHMS = 073500
      TES(475)%NHMS = 073600
      TES(476)%NHMS = 075200
      TES(477)%NHMS = 075300
      TES(478)%NHMS = 075500
      TES(479)%NHMS = 075500
      TES(480)%NHMS = 075600
      TES(481)%NHMS = 075700
      TES(482)%NHMS = 075700
      TES(483)%NHMS = 075900
      TES(484)%NHMS = 080000
      TES(485)%NHMS = 080400
      TES(486)%NHMS = 080500
      TES(487)%NHMS = 090700
      TES(488)%NHMS = 091200
      TES(489)%NHMS = 091400
      TES(490)%NHMS = 091400
      TES(491)%NHMS = 091500
      TES(492)%NHMS = 093300
      TES(493)%NHMS = 093300
      TES(494)%NHMS = 093400
      TES(495)%NHMS = 093500
      TES(496)%NHMS = 093600
      TES(497)%NHMS = 093600
      TES(498)%NHMS = 093700
      TES(499)%NHMS = 093800
      TES(500)%NHMS = 093800
      TES(501)%NHMS = 094000
      TES(502)%NHMS = 104100
      TES(503)%NHMS = 104200
      TES(504)%NHMS = 104800
      TES(505)%NHMS = 104800
      TES(506)%NHMS = 105000
      TES(507)%NHMS = 105100
      TES(508)%NHMS = 105100
      TES(509)%NHMS = 105200
      TES(510)%NHMS = 105200
      TES(511)%NHMS = 105300
      TES(512)%NHMS = 105300
      TES(513)%NHMS = 105300
      TES(514)%NHMS = 105400
      TES(515)%NHMS = 111100
      TES(516)%NHMS = 111100
      TES(517)%NHMS = 111200
      TES(518)%NHMS = 095600
      TES(519)%NHMS = 095700
      TES(520)%NHMS = 101500
      TES(521)%NHMS = 101500
      TES(522)%NHMS = 101600
      TES(523)%NHMS = 101600
      TES(524)%NHMS = 101700
      TES(525)%NHMS = 101700
      TES(526)%NHMS = 101800
      TES(527)%NHMS = 101900
      TES(528)%NHMS = 102000
      TES(529)%NHMS = 103700
      TES(530)%NHMS = 111900
      TES(531)%NHMS = 112100
      TES(532)%NHMS = 112100
      TES(533)%NHMS = 112200
      TES(534)%NHMS = 112300
      TES(535)%NHMS = 112900
      TES(536)%NHMS = 113100
      TES(537)%NHMS = 113100
      TES(538)%NHMS = 113200
      TES(539)%NHMS = 113200
      TES(540)%NHMS = 113300
      TES(541)%NHMS = 113500
      TES(542)%NHMS = 113500
      TES(543)%NHMS = 113600
      TES(544)%NHMS = 115300
      TES(545)%NHMS = 125000
      TES(546)%NHMS = 125300
      TES(547)%NHMS = 125300
      TES(548)%NHMS = 125400
      TES(549)%NHMS = 125600
      TES(550)%NHMS = 125700
      TES(551)%NHMS = 125900
      TES(552)%NHMS = 125900
      TES(553)%NHMS = 130000
      TES(554)%NHMS = 130000
      TES(555)%NHMS = 130100
      TES(556)%NHMS = 130200
      TES(557)%NHMS = 130900
      TES(558)%NHMS = 130900
      TES(559)%NHMS = 131000
      TES(560)%NHMS = 131200
      TES(561)%NHMS = 131300
      TES(562)%NHMS = 131300
      TES(563)%NHMS = 131400
      TES(564)%NHMS = 131500
      TES(565)%NHMS = 131500
      TES(566)%NHMS = 135200
      TES(567)%NHMS = 140100
      TES(568)%NHMS = 140200
      TES(569)%NHMS = 153400
      TES(570)%NHMS = 153400
      TES(571)%NHMS = 153500
      TES(572)%NHMS = 153600
      TES(573)%NHMS = 153600
      TES(574)%NHMS = 153900
      TES(575)%NHMS = 161500
      TES(576)%NHMS = 165300
      TES(577)%NHMS = 165500
      TES(578)%NHMS = 165600
      TES(579)%NHMS = 165600
      TES(580)%NHMS = 170600
      TES(581)%NHMS = 170700
      TES(582)%NHMS = 171100
      TES(583)%NHMS = 171300
      TES(584)%NHMS = 171600
      TES(585)%NHMS = 174700
      TES(586)%NHMS = 174700
      TES(587)%NHMS = 174700
      TES(588)%NHMS = 174900
      TES(589)%NHMS = 175000
      TES(590)%NHMS = 175000
      TES(591)%NHMS = 175100
      TES(592)%NHMS = 175100
      TES(593)%NHMS = 175200
      TES(594)%NHMS = 175200
      TES(595)%NHMS = 175300
      TES(596)%NHMS = 175400
      TES(597)%NHMS = 175400
      TES(598)%NHMS = 175500
      TES(599)%NHMS = 175700
      TES(600)%NHMS = 175800
      TES(601)%NHMS = 175800
      TES(602)%NHMS = 180800
      TES(603)%NHMS = 180900
      TES(604)%NHMS = 181000
      TES(605)%NHMS = 181100
      TES(606)%NHMS = 182900
      TES(607)%NHMS = 183200
      TES(608)%NHMS = 183400
      TES(609)%NHMS = 183500
      TES(610)%NHMS = 183700
      TES(611)%NHMS = 183800
      TES(612)%NHMS = 183800
      TES(613)%NHMS = 183900
      TES(614)%NHMS = 184200
      TES(615)%NHMS = 184700
      TES(616)%NHMS = 193900
      TES(617)%NHMS = 193900
      TES(618)%NHMS = 194300
      TES(619)%NHMS = 194600
      TES(620)%NHMS = 194600
      TES(621)%NHMS = 194700
      TES(622)%NHMS = 194700
      TES(623)%NHMS = 194800
      TES(624)%NHMS = 194900
      TES(625)%NHMS = 195000
      TES(626)%NHMS = 200800
      TES(627)%NHMS = 200800
      TES(628)%NHMS = 200900
      TES(629)%NHMS = 200900
      TES(630)%NHMS = 201000
      TES(631)%NHMS = 201200
      TES(632)%NHMS = 201400
      TES(633)%NHMS = 201400
      TES(634)%NHMS = 201500
      TES(635)%NHMS = 201700
      TES(636)%NHMS = 201700
      TES(637)%NHMS = 201800
      TES(638)%NHMS = 201800
      TES(639)%NHMS = 212500
      TES(640)%NHMS = 212500
      TES(641)%NHMS = 212700
      TES(642)%NHMS = 212800
      TES(643)%NHMS = 212900
      TES(644)%NHMS = 212900
      TES(645)%NHMS = 214700
      TES(646)%NHMS = 214800
      TES(647)%NHMS = 214900
      TES(648)%NHMS = 214900
      TES(649)%NHMS = 215000
      TES(650)%NHMS = 215400
      TES(651)%NHMS = 215500
      TES(652)%NHMS = 230900
      TES(653)%NHMS = 232600
      TES(654)%NHMS = 232900
      TES(655)%NHMS = 232900
      TES(656)%NHMS = 233000
      TES(657)%NHMS = 233000
      TES(658)%NHMS = 233100
      TES(659)%NHMS = 233100
      TES(660)%NHMS = 233200
      TES(661)%NHMS = 233900
      TES(662)%NHMS = 234100
      TES(663)%NHMS = 234100
      TES(664)%NHMS = 234200
      TES(665)%NHMS = 234200
      TES(666)%NHMS = 234300
      TES(667)%NHMS = 234300
      TES(668)%NHMS = 234400
      TES(669)%NHMS = 234400
      TES(670)%NHMS = 234500
      TES(671)%NHMS = 234500
      TES(672)%NHMS = 234600
      TES(673)%NHMS = 234600
      TES(674)%NHMS = 234700
      TES(675)%NHMS = 234900
      TES(676)%NHMS = 235100
      TES(677)%NHMS = 010400
      TES(678)%NHMS = 010500
      TES(679)%NHMS = 010600
      TES(680)%NHMS = 010600
      TES(681)%NHMS = 010700
      TES(682)%NHMS = 010800
      TES(683)%NHMS = 010800
      TES(684)%NHMS = 011000
      TES(685)%NHMS = 011000
      TES(686)%NHMS = 011400
      TES(687)%NHMS = 011500
      TES(688)%NHMS = 011600
      TES(689)%NHMS = 011800
      TES(690)%NHMS = 011800
      TES(691)%NHMS = 011900
      TES(692)%NHMS = 011900
      TES(693)%NHMS = 022500
      TES(694)%NHMS = 022600
      TES(695)%NHMS = 024800
      TES(696)%NHMS = 024800
      TES(697)%NHMS = 024900
      TES(698)%NHMS = 024900
      TES(699)%NHMS = 025300
      TES(700)%NHMS = 025400
      TES(701)%NHMS = 025600
      TES(702)%NHMS = 034800
      TES(703)%NHMS = 040000
      TES(704)%NHMS = 040000
      TES(705)%NHMS = 040200
      TES(706)%NHMS = 040300
      TES(707)%NHMS = 040300
      TES(708)%NHMS = 040400
      TES(709)%NHMS = 040400
      TES(710)%NHMS = 040500
      TES(711)%NHMS = 040600
      TES(712)%NHMS = 044000
      TES(713)%NHMS = 044100
      TES(714)%NHMS = 044100
      TES(715)%NHMS = 044100
      TES(716)%NHMS = 044300
      TES(717)%NHMS = 044400
      TES(718)%NHMS = 044400
      TES(719)%NHMS = 044500
      TES(720)%NHMS = 044600
      TES(721)%NHMS = 052200
      TES(722)%NHMS = 052300
      TES(723)%NHMS = 052400
      TES(724)%NHMS = 053000
      TES(725)%NHMS = 053100
      TES(726)%NHMS = 053100
      TES(727)%NHMS = 053200
      TES(728)%NHMS = 053300
      TES(729)%NHMS = 053300
      TES(730)%NHMS = 053500
      TES(731)%NHMS = 053600
      TES(732)%NHMS = 053700
      TES(733)%NHMS = 053800
      TES(734)%NHMS = 053800
      TES(735)%NHMS = 053800
      TES(736)%NHMS = 053900
      TES(737)%NHMS = 054100
      TES(738)%NHMS = 054100
      TES(739)%NHMS = 054200
      TES(740)%NHMS = 054200
      TES(741)%NHMS = 054300
      TES(742)%NHMS = 054300
      TES(743)%NHMS = 054400
      TES(744)%NHMS = 054500
      TES(745)%NHMS = 060400
      TES(746)%NHMS = 061500
      TES(747)%NHMS = 061800
      TES(748)%NHMS = 061900
      TES(749)%NHMS = 062000
      TES(750)%NHMS = 062000
      TES(751)%NHMS = 062100
      TES(752)%NHMS = 062200
      TES(753)%NHMS = 062200
      TES(754)%NHMS = 070700
      TES(755)%NHMS = 070800
      TES(756)%NHMS = 071400
      TES(757)%NHMS = 072000
      TES(758)%NHMS = 072100
      TES(759)%NHMS = 072100
      TES(760)%NHMS = 072200
      TES(761)%NHMS = 072200
      TES(762)%NHMS = 072300
      TES(763)%NHMS = 072300
      TES(764)%NHMS = 072400
      TES(765)%NHMS = 074000
      TES(766)%NHMS = 074000
      TES(767)%NHMS = 074100
      TES(768)%NHMS = 074100
      TES(769)%NHMS = 074200
      TES(770)%NHMS = 074200
      TES(771)%NHMS = 074300
      TES(772)%NHMS = 074400
      TES(773)%NHMS = 074400
      TES(774)%NHMS = 074400
      TES(775)%NHMS = 074500
      TES(776)%NHMS = 074600
      TES(777)%NHMS = 074700
      TES(778)%NHMS = 075100
      TES(779)%NHMS = 075200
      TES(780)%NHMS = 085500
      TES(781)%NHMS = 090000
      TES(782)%NHMS = 090100
      TES(783)%NHMS = 091900
      TES(784)%NHMS = 092000
      TES(785)%NHMS = 092000
      TES(786)%NHMS = 092000
      TES(787)%NHMS = 092100
      TES(788)%NHMS = 092100
      TES(789)%NHMS = 092200
      TES(790)%NHMS = 092400
      TES(791)%NHMS = 092400
      TES(792)%NHMS = 092500
      TES(793)%NHMS = 092600
      TES(794)%NHMS = 092600
      TES(795)%NHMS = 102700
      TES(796)%NHMS = 103000
      TES(797)%NHMS = 103600
      TES(798)%NHMS = 103700
      TES(799)%NHMS = 103800
      TES(800)%NHMS = 103900
      TES(801)%NHMS = 103900
      TES(802)%NHMS = 104100
      TES(803)%NHMS = 105900
      TES(804)%NHMS = 092700
      TES(805)%NHMS = 093100
      TES(806)%NHMS = 093200
      TES(807)%NHMS = 093200
      TES(808)%NHMS = 093300
      TES(809)%NHMS = 093300
      TES(810)%NHMS = 093300
      TES(811)%NHMS = 095000
      TES(812)%NHMS = 095200
      TES(813)%NHMS = 095300
      TES(814)%NHMS = 095300
      TES(815)%NHMS = 095300
      TES(816)%NHMS = 095400
      TES(817)%NHMS = 095400
      TES(818)%NHMS = 095500
      TES(819)%NHMS = 095500
      TES(820)%NHMS = 104900
      TES(821)%NHMS = 105000
      TES(822)%NHMS = 105100
      TES(823)%NHMS = 105100
      TES(824)%NHMS = 105700
      TES(825)%NHMS = 105700
      TES(826)%NHMS = 105800
      TES(827)%NHMS = 110600
      TES(828)%NHMS = 110700
      TES(829)%NHMS = 111000
      TES(830)%NHMS = 111000
      TES(831)%NHMS = 111100
      TES(832)%NHMS = 111100
      TES(833)%NHMS = 111200
      TES(834)%NHMS = 112900
      TES(835)%NHMS = 112900
      TES(836)%NHMS = 112900
      TES(837)%NHMS = 122600
      TES(838)%NHMS = 122600
      TES(839)%NHMS = 122700
      TES(840)%NHMS = 122900
      TES(841)%NHMS = 123000
      TES(842)%NHMS = 123000
      TES(843)%NHMS = 123100
      TES(844)%NHMS = 123200
      TES(845)%NHMS = 123300
      TES(846)%NHMS = 123300
      TES(847)%NHMS = 123400
      TES(848)%NHMS = 123400
      TES(849)%NHMS = 123900
      TES(850)%NHMS = 124000
      TES(851)%NHMS = 124500
      TES(852)%NHMS = 124700
      TES(853)%NHMS = 130800
      TES(854)%NHMS = 133500
      TES(855)%NHMS = 133500
      TES(856)%NHMS = 141700
      TES(857)%NHMS = 141700
      TES(858)%NHMS = 164500
      TES(859)%NHMS = 172600
      TES(860)%NHMS = 172700
      TES(861)%NHMS = 172800
      TES(862)%NHMS = 172800
      TES(863)%NHMS = 173000
      TES(864)%NHMS = 173200
      TES(865)%NHMS = 173200
      TES(866)%NHMS = 173300
      TES(867)%NHMS = 174400
      TES(868)%NHMS = 174400
      TES(869)%NHMS = 174500
      TES(870)%NHMS = 174500
      TES(871)%NHMS = 174500
      TES(872)%NHMS = 174600
      TES(873)%NHMS = 174600
      TES(874)%NHMS = 174700
      TES(875)%NHMS = 174800
      TES(876)%NHMS = 180400
      TES(877)%NHMS = 180500
      TES(878)%NHMS = 180600
      TES(879)%NHMS = 180700
      TES(880)%NHMS = 180700
      TES(881)%NHMS = 180800
      TES(882)%NHMS = 181000
      TES(883)%NHMS = 181300
      TES(884)%NHMS = 181300
      TES(885)%NHMS = 181400
      TES(886)%NHMS = 182000
      TES(887)%NHMS = 182000
      TES(888)%NHMS = 182100
      TES(889)%NHMS = 182100
      TES(890)%NHMS = 185600
      TES(891)%NHMS = 185700
      TES(892)%NHMS = 185700
      TES(893)%NHMS = 185900
      TES(894)%NHMS = 185900
      TES(895)%NHMS = 191200
      TES(896)%NHMS = 191400
      TES(897)%NHMS = 191900
      TES(898)%NHMS = 191900
      TES(899)%NHMS = 192100
      TES(900)%NHMS = 192100
      TES(901)%NHMS = 192100
      TES(902)%NHMS = 192200
      TES(903)%NHMS = 192300
      TES(904)%NHMS = 192300
      TES(905)%NHMS = 194300
      TES(906)%NHMS = 194400
      TES(907)%NHMS = 194400
      TES(908)%NHMS = 194500
      TES(909)%NHMS = 194700
      TES(910)%NHMS = 194900
      TES(911)%NHMS = 195100
      TES(912)%NHMS = 195100
      TES(913)%NHMS = 195300
      TES(914)%NHMS = 210100
      TES(915)%NHMS = 210200
      TES(916)%NHMS = 210200
      TES(917)%NHMS = 210200
      TES(918)%NHMS = 210300
      TES(919)%NHMS = 212200
      TES(920)%NHMS = 212200
      TES(921)%NHMS = 212300
      TES(922)%NHMS = 212400
      TES(923)%NHMS = 212500
      TES(924)%NHMS = 212500
      TES(925)%NHMS = 212600
      TES(926)%NHMS = 212600
      TES(927)%NHMS = 212700
      TES(928)%NHMS = 212700
      TES(929)%NHMS = 212800
      TES(930)%NHMS = 212900
      TES(931)%NHMS = 230100
      TES(932)%NHMS = 230100
      TES(933)%NHMS = 230200
      TES(934)%NHMS = 230300
      TES(935)%NHMS = 230400
      TES(936)%NHMS = 230600
      TES(937)%NHMS = 230700
      TES(938)%NHMS = 230800
      TES(939)%NHMS = 231500
      TES(940)%NHMS = 231500
      TES(941)%NHMS = 231500
      TES(942)%NHMS = 231600
      TES(943)%NHMS = 231600
      TES(944)%NHMS = 231700
      TES(945)%NHMS = 231900
      TES(946)%NHMS = 232000
      TES(947)%NHMS = 232000
      TES(948)%NHMS = 232100
      TES(949)%NHMS = 232200
      TES(950)%NHMS = 232300
      TES(951)%NHMS = 232500
      TES(952)%NHMS = 001300
      TES(953)%NHMS = 004000
      TES(954)%NHMS = 004000
      TES(955)%NHMS = 004100
      TES(956)%NHMS = 004100
      TES(957)%NHMS = 004100
      TES(958)%NHMS = 004200
      TES(959)%NHMS = 004200
      TES(960)%NHMS = 004300
      TES(961)%NHMS = 004300
      TES(962)%NHMS = 004400
      TES(963)%NHMS = 004400
      TES(964)%NHMS = 004500
      TES(965)%NHMS = 004500
      TES(966)%NHMS = 004600
      TES(967)%NHMS = 004800
      TES(968)%NHMS = 004900
      TES(969)%NHMS = 004900
      TES(970)%NHMS = 005000
      TES(971)%NHMS = 005100
      TES(972)%NHMS = 005200
      TES(973)%NHMS = 005600
      TES(974)%NHMS = 005600
      TES(975)%NHMS = 005600
      TES(976)%NHMS = 005700
      TES(977)%NHMS = 005700
      TES(978)%NHMS = 014000
      TES(979)%NHMS = 021800
      TES(980)%NHMS = 022200
      TES(981)%NHMS = 022200
      TES(982)%NHMS = 022300
      TES(983)%NHMS = 022400
      TES(984)%NHMS = 022400
      TES(985)%NHMS = 022500
      TES(986)%NHMS = 023000
      TES(987)%NHMS = 023100
      TES(988)%NHMS = 023200
      TES(989)%NHMS = 023200
      TES(990)%NHMS = 033800
      TES(991)%NHMS = 033800
      TES(992)%NHMS = 033900
      TES(993)%NHMS = 034100
      TES(994)%NHMS = 041500
      TES(995)%NHMS = 041600
      TES(996)%NHMS = 041700
      TES(997)%NHMS = 041700
      TES(998)%NHMS = 041800
      TES(999)%NHMS = 041800
      TES(1000)%NHMS = 041900
      TES(1001)%NHMS = 041900
      TES(1002)%NHMS = 045500
      TES(1003)%NHMS = 045800
      TES(1004)%NHMS = 045900
      TES(1005)%NHMS = 051400
      TES(1006)%NHMS = 051700
      TES(1007)%NHMS = 051800
      TES(1008)%NHMS = 051900
      TES(1009)%NHMS = 052000
      TES(1010)%NHMS = 055000
      TES(1011)%NHMS = 055100
      TES(1012)%NHMS = 055400
      TES(1013)%NHMS = 055400
      TES(1014)%NHMS = 055500
      TES(1015)%NHMS = 055600
      TES(1016)%NHMS = 055700
      TES(1017)%NHMS = 055800
      TES(1018)%NHMS = 055900
      TES(1019)%NHMS = 060100
      TES(1020)%NHMS = 060200
      TES(1021)%NHMS = 060200
      TES(1022)%NHMS = 060300
      TES(1023)%NHMS = 060300
      TES(1024)%NHMS = 060400
      TES(1025)%NHMS = 060400
      TES(1026)%NHMS = 060500
      TES(1027)%NHMS = 060500
      TES(1028)%NHMS = 064100
      TES(1029)%NHMS = 064100
      TES(1030)%NHMS = 064500
      TES(1031)%NHMS = 064600
      TES(1032)%NHMS = 064800
      TES(1033)%NHMS = 064900
      TES(1034)%NHMS = 065100
      TES(1035)%NHMS = 065300
      TES(1036)%NHMS = 065300
      TES(1037)%NHMS = 065500
      TES(1038)%NHMS = 065600
      TES(1039)%NHMS = 065700
      TES(1040)%NHMS = 065700
      TES(1041)%NHMS = 065800
      TES(1042)%NHMS = 065800
      TES(1043)%NHMS = 065900
      TES(1044)%NHMS = 071700
      TES(1045)%NHMS = 071700
      TES(1046)%NHMS = 071700
      TES(1047)%NHMS = 071800
      TES(1048)%NHMS = 071900
      TES(1049)%NHMS = 071900
      TES(1050)%NHMS = 072000
      TES(1051)%NHMS = 072000
      TES(1052)%NHMS = 072100
      TES(1053)%NHMS = 072200
      TES(1054)%NHMS = 072200
      TES(1055)%NHMS = 072500
      TES(1056)%NHMS = 082300
      TES(1057)%NHMS = 082800
      TES(1058)%NHMS = 082800
      TES(1059)%NHMS = 082800
      TES(1060)%NHMS = 082900
      TES(1061)%NHMS = 083200
      TES(1062)%NHMS = 083500
      TES(1063)%NHMS = 083500
      TES(1064)%NHMS = 083700
      TES(1065)%NHMS = 085700
      TES(1066)%NHMS = 085800
      TES(1067)%NHMS = 090000
      TES(1068)%NHMS = 090100
      TES(1069)%NHMS = 090200
      TES(1070)%NHMS = 090300
      TES(1071)%NHMS = 090300
      TES(1072)%NHMS = 090300
      TES(1073)%NHMS = 100600
      TES(1074)%NHMS = 100700
      TES(1075)%NHMS = 100900
      TES(1076)%NHMS = 101000
      TES(1077)%NHMS = 101100
      TES(1078)%NHMS = 101100
      TES(1079)%NHMS = 101200
      TES(1080)%NHMS = 101400
      TES(1081)%NHMS = 101400
      TES(1082)%NHMS = 101500
      TES(1083)%NHMS = 101500
      TES(1084)%NHMS = 103300
      TES(1085)%NHMS = 103400
      TES(1086)%NHMS = 103400
      TES(1087)%NHMS = 103500
      TES(1088)%NHMS = 103600
      TES(1089)%NHMS = 103600
      TES(1090)%NHMS = 103700
      TES(1091)%NHMS = 091700
      TES(1092)%NHMS = 091800
      TES(1093)%NHMS = 091900
      TES(1094)%NHMS = 092100
      TES(1095)%NHMS = 092100
      TES(1096)%NHMS = 093700
      TES(1097)%NHMS = 093800
      TES(1098)%NHMS = 093900
      TES(1099)%NHMS = 094000
      TES(1100)%NHMS = 094200
      TES(1101)%NHMS = 094200
      TES(1102)%NHMS = 094300
      TES(1103)%NHMS = 094300
      TES(1104)%NHMS = 094400
      TES(1105)%NHMS = 094400
      TES(1106)%NHMS = 094600
      TES(1107)%NHMS = 104400
      TES(1108)%NHMS = 104500
      TES(1109)%NHMS = 104600
      TES(1110)%NHMS = 105100
      TES(1111)%NHMS = 105300
      TES(1112)%NHMS = 105400
      TES(1113)%NHMS = 105600
      TES(1114)%NHMS = 105700
      TES(1115)%NHMS = 105700
      TES(1116)%NHMS = 105800
      TES(1117)%NHMS = 105800
      TES(1118)%NHMS = 105900
      TES(1119)%NHMS = 111600
      TES(1120)%NHMS = 111700
      TES(1121)%NHMS = 111800
      TES(1122)%NHMS = 111800
      TES(1123)%NHMS = 121400
      TES(1124)%NHMS = 121400
      TES(1125)%NHMS = 121600
      TES(1126)%NHMS = 121800
      TES(1127)%NHMS = 121800
      TES(1128)%NHMS = 121900
      TES(1129)%NHMS = 121900
      TES(1130)%NHMS = 122000
      TES(1131)%NHMS = 122100
      TES(1132)%NHMS = 122100
      TES(1133)%NHMS = 122400
      TES(1134)%NHMS = 122400
      TES(1135)%NHMS = 122500
      TES(1136)%NHMS = 122500
      TES(1137)%NHMS = 123000
      TES(1138)%NHMS = 123300
      TES(1139)%NHMS = 123400
      TES(1140)%NHMS = 123500
      TES(1141)%NHMS = 123600
      TES(1142)%NHMS = 163200
      TES(1143)%NHMS = 163400
      TES(1144)%NHMS = 163600
      TES(1145)%NHMS = 163800
      TES(1146)%NHMS = 171300
      TES(1147)%NHMS = 171400
      TES(1148)%NHMS = 171400
      TES(1149)%NHMS = 171500
      TES(1150)%NHMS = 171500
      TES(1151)%NHMS = 171600
      TES(1152)%NHMS = 171600
      TES(1153)%NHMS = 171700
      TES(1154)%NHMS = 171700
      TES(1155)%NHMS = 171800
      TES(1156)%NHMS = 171800
      TES(1157)%NHMS = 171800
      TES(1158)%NHMS = 171900
      TES(1159)%NHMS = 171900
      TES(1160)%NHMS = 172000
      TES(1161)%NHMS = 173200
      TES(1162)%NHMS = 173200
      TES(1163)%NHMS = 173500
      TES(1164)%NHMS = 175300
      TES(1165)%NHMS = 175300
      TES(1166)%NHMS = 175400
      TES(1167)%NHMS = 175500
      TES(1168)%NHMS = 175500
      TES(1169)%NHMS = 175600
      TES(1170)%NHMS = 175700
      TES(1171)%NHMS = 180000
      TES(1172)%NHMS = 180700
      TES(1173)%NHMS = 180800
      TES(1174)%NHMS = 180800
      TES(1175)%NHMS = 180900
      TES(1176)%NHMS = 181000
      TES(1177)%NHMS = 184700
      TES(1178)%NHMS = 184700
      TES(1179)%NHMS = 184800
      TES(1180)%NHMS = 184900
      TES(1181)%NHMS = 185400
      TES(1182)%NHMS = 185500
      TES(1183)%NHMS = 185600
      TES(1184)%NHMS = 185600
      TES(1185)%NHMS = 185700
      TES(1186)%NHMS = 185700
      TES(1187)%NHMS = 185800
      TES(1188)%NHMS = 185800
      TES(1189)%NHMS = 190400
      TES(1190)%NHMS = 190700
      TES(1191)%NHMS = 190900
      TES(1192)%NHMS = 191000
      TES(1193)%NHMS = 191100
      TES(1194)%NHMS = 191200
      TES(1195)%NHMS = 193000
      TES(1196)%NHMS = 193100
      TES(1197)%NHMS = 193200
      TES(1198)%NHMS = 193300
      TES(1199)%NHMS = 193400
      TES(1200)%NHMS = 194000
      TES(1201)%NHMS = 194000
      TES(1202)%NHMS = 194100
      TES(1203)%NHMS = 194400
      TES(1204)%NHMS = 194400
      TES(1205)%NHMS = 204500
      TES(1206)%NHMS = 204500
      TES(1207)%NHMS = 204600
      TES(1208)%NHMS = 204700
      TES(1209)%NHMS = 204800
      TES(1210)%NHMS = 205000
      TES(1211)%NHMS = 205100
      TES(1212)%NHMS = 205100
      TES(1213)%NHMS = 205300
      TES(1214)%NHMS = 205300
      TES(1215)%NHMS = 211000
      TES(1216)%NHMS = 211100
      TES(1217)%NHMS = 211200
      TES(1218)%NHMS = 211300
      TES(1219)%NHMS = 211400
      TES(1220)%NHMS = 211600
      TES(1221)%NHMS = 211700
      TES(1222)%NHMS = 211700
      TES(1223)%NHMS = 212000
      TES(1224)%NHMS = 224800
      TES(1225)%NHMS = 224900
      TES(1226)%NHMS = 224900
      TES(1227)%NHMS = 225000
      TES(1228)%NHMS = 225000
      TES(1229)%NHMS = 225100
      TES(1230)%NHMS = 225800
      TES(1231)%NHMS = 225900
      TES(1232)%NHMS = 225900
      TES(1233)%NHMS = 230000
      TES(1234)%NHMS = 230300
      TES(1235)%NHMS = 230300
      TES(1236)%NHMS = 230400
      TES(1237)%NHMS = 230500
      TES(1238)%NHMS = 230600
      TES(1239)%NHMS = 230700
      TES(1240)%NHMS = 230700
      TES(1241)%NHMS = 230800
      TES(1242)%NHMS = 230800
      TES(1243)%NHMS = 230900
      TES(1244)%NHMS = 230900
      TES(1245)%NHMS = 231000
      TES(1246)%NHMS = 231100
      TES(1247)%NHMS = 002700
      TES(1248)%NHMS = 002800
      TES(1249)%NHMS = 002800
      TES(1250)%NHMS = 002900
      TES(1251)%NHMS = 003000
      TES(1252)%NHMS = 003000
      TES(1253)%NHMS = 003200
      TES(1254)%NHMS = 003200
      TES(1255)%NHMS = 003600
      TES(1256)%NHMS = 003600
      TES(1257)%NHMS = 003700
      TES(1258)%NHMS = 003700
      TES(1259)%NHMS = 003800
      TES(1260)%NHMS = 004100
      TES(1261)%NHMS = 004200
      TES(1262)%NHMS = 004200
      TES(1263)%NHMS = 004500
      TES(1264)%NHMS = 004600
      TES(1265)%NHMS = 004700
      TES(1266)%NHMS = 004800
      TES(1267)%NHMS = 004800
      TES(1268)%NHMS = 020600
      TES(1269)%NHMS = 021000
      TES(1270)%NHMS = 021000
      TES(1271)%NHMS = 021300
      TES(1272)%NHMS = 021400
      TES(1273)%NHMS = 022100
      TES(1274)%NHMS = 030700
      TES(1275)%NHMS = 040400
      TES(1276)%NHMS = 044000
      TES(1277)%NHMS = 044500
      TES(1278)%NHMS = 044500
      TES(1279)%NHMS = 044600
      TES(1280)%NHMS = 044600
      TES(1281)%NHMS = 044800
      TES(1282)%NHMS = 050400
      TES(1283)%NHMS = 050800
      TES(1284)%NHMS = 053900
      TES(1285)%NHMS = 053900
      TES(1286)%NHMS = 054000
      TES(1287)%NHMS = 054000
      TES(1288)%NHMS = 054100
      TES(1289)%NHMS = 054100
      TES(1290)%NHMS = 054200
      TES(1291)%NHMS = 054200
      TES(1292)%NHMS = 054300
      TES(1293)%NHMS = 054300
      TES(1294)%NHMS = 054400
      TES(1295)%NHMS = 054500
      TES(1296)%NHMS = 054800
      TES(1297)%NHMS = 054900
      TES(1298)%NHMS = 054900
      TES(1299)%NHMS = 054900
      TES(1300)%NHMS = 055000
      TES(1301)%NHMS = 055000
      TES(1302)%NHMS = 055100
      TES(1303)%NHMS = 055200
      TES(1304)%NHMS = 055200
      TES(1305)%NHMS = 055300
      TES(1306)%NHMS = 055300
      TES(1307)%NHMS = 062000
      TES(1308)%NHMS = 062100
      TES(1309)%NHMS = 062100
      TES(1310)%NHMS = 062200
      TES(1311)%NHMS = 062200
      TES(1312)%NHMS = 062300
      TES(1313)%NHMS = 062400
      TES(1314)%NHMS = 063000
      TES(1315)%NHMS = 063000
      TES(1316)%NHMS = 063400
      TES(1317)%NHMS = 063600
      TES(1318)%NHMS = 063700
      TES(1319)%NHMS = 063700
      TES(1320)%NHMS = 063800
      TES(1321)%NHMS = 064300
      TES(1322)%NHMS = 064400
      TES(1323)%NHMS = 064400
      TES(1324)%NHMS = 064500
      TES(1325)%NHMS = 064500
      TES(1326)%NHMS = 064500
      TES(1327)%NHMS = 064600
      TES(1328)%NHMS = 064600
      TES(1329)%NHMS = 070800
      TES(1330)%NHMS = 071200
      TES(1331)%NHMS = 071300
      TES(1332)%NHMS = 081500
      TES(1333)%NHMS = 081600
      TES(1334)%NHMS = 081600
      TES(1335)%NHMS = 081700
      TES(1336)%NHMS = 082000
      TES(1337)%NHMS = 082100
      TES(1338)%NHMS = 082100
      TES(1339)%NHMS = 082200
      TES(1340)%NHMS = 082300
      TES(1341)%NHMS = 082500
      TES(1342)%NHMS = 082500
      TES(1343)%NHMS = 084300
      TES(1344)%NHMS = 084400
      TES(1345)%NHMS = 084400
      TES(1346)%NHMS = 084600
      TES(1347)%NHMS = 084700
      TES(1348)%NHMS = 084800
      TES(1349)%NHMS = 084900
      TES(1350)%NHMS = 084900
      TES(1351)%NHMS = 085000
      TES(1352)%NHMS = 085100
      TES(1353)%NHMS = 085100
      TES(1354)%NHMS = 085200
      TES(1355)%NHMS = 085200
      TES(1356)%NHMS = 100200
      TES(1357)%NHMS = 100200
      TES(1358)%NHMS = 100300
      TES(1359)%NHMS = 100300
      TES(1360)%NHMS = 100400
      TES(1361)%NHMS = 102100
      TES(1362)%NHMS = 102100
      TES(1363)%NHMS = 102200
      TES(1364)%NHMS = 102200
      TES(1365)%NHMS = 102300
      TES(1366)%NHMS = 102300
      TES(1367)%NHMS = 102400
      TES(1368)%NHMS = 102400
      TES(1369)%NHMS = 102500
      TES(1370)%NHMS = 102500
      TES(1371)%NHMS = 102600
      TES(1372)%NHMS = 102600
      TES(1373)%NHMS = 103300
      TES(1374)%NHMS = 103400
      TES(1375)%NHMS = 104400
      TES(1376)%NHMS = 104500
      TES(1377)%NHMS = 104600
      TES(1378)%NHMS = 104700
      TES(1379)%NHMS = 110500
      TES(1380)%NHMS = 110500
      TES(1381)%NHMS = 120200
      TES(1382)%NHMS = 120200
      TES(1383)%NHMS = 120300
      TES(1384)%NHMS = 120300
      TES(1385)%NHMS = 120400
      TES(1386)%NHMS = 120400
      TES(1387)%NHMS = 120500
      TES(1388)%NHMS = 120500
      TES(1389)%NHMS = 120600
      TES(1390)%NHMS = 120600
      TES(1391)%NHMS = 120800
      TES(1392)%NHMS = 121000
      TES(1393)%NHMS = 121100
      TES(1394)%NHMS = 121200
      TES(1395)%NHMS = 122100
      TES(1396)%NHMS = 122200
      TES(1397)%NHMS = 122300
      TES(1398)%NHMS = 122400
      TES(1399)%NHMS = 122400
      TES(1400)%NHMS = 122500
      TES(1401)%NHMS = 135200
      TES(1402)%NHMS = 135300
      TES(1403)%NHMS = 135800
      TES(1404)%NHMS = 143600
      TES(1405)%NHMS = 160100
      TES(1406)%NHMS = 160200
      TES(1407)%NHMS = 160200
      TES(1408)%NHMS = 160200
      TES(1409)%NHMS = 161900
      TES(1410)%NHMS = 161900
      TES(1411)%NHMS = 162200
      TES(1412)%NHMS = 162200
      TES(1413)%NHMS = 162200
      TES(1414)%NHMS = 162300
      TES(1415)%NHMS = 162600
      TES(1416)%NHMS = 170000
      TES(1417)%NHMS = 170100
      TES(1418)%NHMS = 170100
      TES(1419)%NHMS = 170200
      TES(1420)%NHMS = 170200
      TES(1421)%NHMS = 170400
      TES(1422)%NHMS = 170500
      TES(1423)%NHMS = 170500
      TES(1424)%NHMS = 170600
      TES(1425)%NHMS = 170600
      TES(1426)%NHMS = 171900
      TES(1427)%NHMS = 172100
      TES(1428)%NHMS = 173900
      TES(1429)%NHMS = 174000
      TES(1430)%NHMS = 174300
      TES(1431)%NHMS = 174400
      TES(1432)%NHMS = 175500
      TES(1433)%NHMS = 175800
      TES(1434)%NHMS = 183400
      TES(1435)%NHMS = 183500
      TES(1436)%NHMS = 183500
      TES(1437)%NHMS = 183500
      TES(1438)%NHMS = 183600
      TES(1439)%NHMS = 183600
      TES(1440)%NHMS = 183700
      TES(1441)%NHMS = 184100
      TES(1442)%NHMS = 184100
      TES(1443)%NHMS = 184300
      TES(1444)%NHMS = 184400
      TES(1445)%NHMS = 184500
      TES(1446)%NHMS = 184600
      TES(1447)%NHMS = 184600
      TES(1448)%NHMS = 184800
      TES(1449)%NHMS = 185300
      TES(1450)%NHMS = 185300
      TES(1451)%NHMS = 185400
      TES(1452)%NHMS = 185600
      TES(1453)%NHMS = 185600
      TES(1454)%NHMS = 185700
      TES(1455)%NHMS = 185800
      TES(1456)%NHMS = 185800
      TES(1457)%NHMS = 185900
      TES(1458)%NHMS = 190000
      TES(1459)%NHMS = 190000
      TES(1460)%NHMS = 190000
      TES(1461)%NHMS = 190100
      TES(1462)%NHMS = 191900
      TES(1463)%NHMS = 192000
      TES(1464)%NHMS = 192000
      TES(1465)%NHMS = 192100
      TES(1466)%NHMS = 192100
      TES(1467)%NHMS = 192500
      TES(1468)%NHMS = 192600
      TES(1469)%NHMS = 192600
      TES(1470)%NHMS = 192700
      TES(1471)%NHMS = 192800
      TES(1472)%NHMS = 192900
      TES(1473)%NHMS = 192900
      TES(1474)%NHMS = 203000
      TES(1475)%NHMS = 203000
      TES(1476)%NHMS = 203100
      TES(1477)%NHMS = 203200
      TES(1478)%NHMS = 203400
      TES(1479)%NHMS = 203500
      TES(1480)%NHMS = 203600
      TES(1481)%NHMS = 203900
      TES(1482)%NHMS = 203900
      TES(1483)%NHMS = 205800
      TES(1484)%NHMS = 205900
      TES(1485)%NHMS = 210000
      TES(1486)%NHMS = 210200
      TES(1487)%NHMS = 210500
      TES(1488)%NHMS = 210700
      TES(1489)%NHMS = 223600
      TES(1490)%NHMS = 223700
      TES(1491)%NHMS = 223700
      TES(1492)%NHMS = 223800
      TES(1493)%NHMS = 223800
      TES(1494)%NHMS = 223900
      TES(1495)%NHMS = 224100
      TES(1496)%NHMS = 224400
      TES(1497)%NHMS = 224600
      TES(1498)%NHMS = 224600
      TES(1499)%NHMS = 224800
      TES(1500)%NHMS = 225200
      TES(1501)%NHMS = 225200
      TES(1502)%NHMS = 225700
      TES(1503)%NHMS = 235700
      TES(1504)%NHMS = 001700
      TES(1505)%NHMS = 002000
      TES(1506)%NHMS = 002100
      TES(1507)%NHMS = 002400
      TES(1508)%NHMS = 002500
      TES(1509)%NHMS = 002500
      TES(1510)%NHMS = 002600
      TES(1511)%NHMS = 002800
      TES(1512)%NHMS = 002900
      TES(1513)%NHMS = 002900
      TES(1514)%NHMS = 003000
      TES(1515)%NHMS = 003100
      TES(1516)%NHMS = 003200
      TES(1517)%NHMS = 003200
      TES(1518)%NHMS = 003300
      TES(1519)%NHMS = 003400
      TES(1520)%NHMS = 003400
      TES(1521)%NHMS = 003600
      TES(1522)%NHMS = 003700
      TES(1523)%NHMS = 015400
      TES(1524)%NHMS = 015500
      TES(1525)%NHMS = 015600
      TES(1526)%NHMS = 015800
      TES(1527)%NHMS = 020400
      TES(1528)%NHMS = 020400
      TES(1529)%NHMS = 020500
      TES(1530)%NHMS = 020700
      TES(1531)%NHMS = 020800
      TES(1532)%NHMS = 020800
      TES(1533)%NHMS = 024800
      TES(1534)%NHMS = 031200
      TES(1535)%NHMS = 042700
      TES(1536)%NHMS = 042700
      TES(1537)%NHMS = 042900
      TES(1538)%NHMS = 042900
      TES(1539)%NHMS = 043000
      TES(1540)%NHMS = 043000
      TES(1541)%NHMS = 043200
      TES(1542)%NHMS = 043300
      TES(1543)%NHMS = 043300
      TES(1544)%NHMS = 043700
      TES(1545)%NHMS = 044800
      TES(1546)%NHMS = 044800
      TES(1547)%NHMS = 045400
      TES(1548)%NHMS = 052700
      TES(1549)%NHMS = 052800
      TES(1550)%NHMS = 052800
      TES(1551)%NHMS = 052900
      TES(1552)%NHMS = 053000
      TES(1553)%NHMS = 053000
      TES(1554)%NHMS = 053100
      TES(1555)%NHMS = 053300
      TES(1556)%NHMS = 053400
      TES(1557)%NHMS = 053500
      TES(1558)%NHMS = 053500
      TES(1559)%NHMS = 053600
      TES(1560)%NHMS = 053700
      TES(1561)%NHMS = 053700
      TES(1562)%NHMS = 053800
      TES(1563)%NHMS = 053800
      TES(1564)%NHMS = 053900
      TES(1565)%NHMS = 053900
      TES(1566)%NHMS = 053900
      TES(1567)%NHMS = 054000
      TES(1568)%NHMS = 060800
      TES(1569)%NHMS = 061000
      TES(1570)%NHMS = 061100
      TES(1571)%NHMS = 061200
      TES(1572)%NHMS = 061700
      TES(1573)%NHMS = 061700
      TES(1574)%NHMS = 061800
      TES(1575)%NHMS = 062400
      TES(1576)%NHMS = 062500
      TES(1577)%NHMS = 062500
      TES(1578)%NHMS = 062500
      TES(1579)%NHMS = 062800
      TES(1580)%NHMS = 063000
      TES(1581)%NHMS = 063200
      TES(1582)%NHMS = 063300
      TES(1583)%NHMS = 063300
      TES(1584)%NHMS = 063400
      TES(1585)%NHMS = 063400
      TES(1586)%NHMS = 065400
      TES(1587)%NHMS = 065500
      TES(1588)%NHMS = 070100
      TES(1589)%NHMS = 070500
      TES(1590)%NHMS = 070700
      TES(1591)%NHMS = 070700
      TES(1592)%NHMS = 070800
      TES(1593)%NHMS = 080300
      TES(1594)%NHMS = 080300
      TES(1595)%NHMS = 080400
      TES(1596)%NHMS = 080600
      TES(1597)%NHMS = 080900
      TES(1598)%NHMS = 080900
      TES(1599)%NHMS = 081000
      TES(1600)%NHMS = 081000
      TES(1601)%NHMS = 081100
      TES(1602)%NHMS = 081200
      TES(1603)%NHMS = 081300
      TES(1604)%NHMS = 083000
      TES(1605)%NHMS = 083100
      TES(1606)%NHMS = 083100
      TES(1607)%NHMS = 083200
      TES(1608)%NHMS = 083200
      TES(1609)%NHMS = 083600
      TES(1610)%NHMS = 083700
      TES(1611)%NHMS = 083700
      TES(1612)%NHMS = 083800
      TES(1613)%NHMS = 083800
      TES(1614)%NHMS = 083900
      TES(1615)%NHMS = 084000
      TES(1616)%NHMS = 094500
      TES(1617)%NHMS = 095100
      TES(1618)%NHMS = 095200
      TES(1619)%NHMS = 100800
      TES(1620)%NHMS = 100800
      TES(1621)%NHMS = 101100
      TES(1622)%NHMS = 101200
      TES(1623)%NHMS = 101300
      TES(1624)%NHMS = 101300
      TES(1625)%NHMS = 101400
      TES(1626)%NHMS = 101500
      TES(1627)%NHMS = 111400
      TES(1628)%NHMS = 111500
      TES(1629)%NHMS = 111500
      TES(1630)%NHMS = 111700
      TES(1631)%NHMS = 111800
      TES(1632)%NHMS = 111800
      TES(1633)%NHMS = 112300
      TES(1634)%NHMS = 112400
      TES(1635)%NHMS = 112700
      TES(1636)%NHMS = 112700
      TES(1637)%NHMS = 112800
      TES(1638)%NHMS = 112800
      TES(1639)%NHMS = 112900
      TES(1640)%NHMS = 112900
      TES(1641)%NHMS = 113000
      TES(1642)%NHMS = 114700
      TES(1643)%NHMS = 103100
      TES(1644)%NHMS = 103200
      TES(1645)%NHMS = 103200
      TES(1646)%NHMS = 103200
      TES(1647)%NHMS = 103300
      TES(1648)%NHMS = 103500
      TES(1649)%NHMS = 103500
      TES(1650)%NHMS = 105200
      TES(1651)%NHMS = 105200
      TES(1652)%NHMS = 105200
      TES(1653)%NHMS = 105300
      TES(1654)%NHMS = 105300
      TES(1655)%NHMS = 105400
      TES(1656)%NHMS = 115200
      TES(1657)%NHMS = 115300
      TES(1658)%NHMS = 115300
      TES(1659)%NHMS = 115400
      TES(1660)%NHMS = 115500
      TES(1661)%NHMS = 115500
      TES(1662)%NHMS = 115600
      TES(1663)%NHMS = 115800
      TES(1664)%NHMS = 115900
      TES(1665)%NHMS = 120000
      TES(1666)%NHMS = 120100
      TES(1667)%NHMS = 120300
      TES(1668)%NHMS = 120400
      TES(1669)%NHMS = 120500
      TES(1670)%NHMS = 120900
      TES(1671)%NHMS = 121000
      TES(1672)%NHMS = 121000
      TES(1673)%NHMS = 121100
      TES(1674)%NHMS = 121100
      TES(1675)%NHMS = 121200
      TES(1676)%NHMS = 133800
      TES(1677)%NHMS = 134000
      TES(1678)%NHMS = 134100
      TES(1679)%NHMS = 134600
      TES(1680)%NHMS = 134800
      TES(1681)%NHMS = 154900
      TES(1682)%NHMS = 154900
      TES(1683)%NHMS = 160800
      TES(1684)%NHMS = 160900
      TES(1685)%NHMS = 161300
      TES(1686)%NHMS = 161400
      TES(1687)%NHMS = 164900
      TES(1688)%NHMS = 164900
      TES(1689)%NHMS = 165000
      TES(1690)%NHMS = 165000
      TES(1691)%NHMS = 165100
      TES(1692)%NHMS = 165100
      TES(1693)%NHMS = 165200
      TES(1694)%NHMS = 165300
      TES(1695)%NHMS = 170800
      TES(1696)%NHMS = 170900
      TES(1697)%NHMS = 172700
      TES(1698)%NHMS = 172700
      TES(1699)%NHMS = 173000
      TES(1700)%NHMS = 173100
      TES(1701)%NHMS = 173100
      TES(1702)%NHMS = 173200
      TES(1703)%NHMS = 173600
      TES(1704)%NHMS = 174400
      TES(1705)%NHMS = 174400
      TES(1706)%NHMS = 174500
      TES(1707)%NHMS = 174900
      TES(1708)%NHMS = 175000
      TES(1709)%NHMS = 182200
      TES(1710)%NHMS = 182300
      TES(1711)%NHMS = 182300
      TES(1712)%NHMS = 182400
      TES(1713)%NHMS = 182500
      TES(1714)%NHMS = 182600
      TES(1715)%NHMS = 182900
      TES(1716)%NHMS = 183000
      TES(1717)%NHMS = 183000
      TES(1718)%NHMS = 183100
      TES(1719)%NHMS = 183100
      TES(1720)%NHMS = 183200
      TES(1721)%NHMS = 183200
      TES(1722)%NHMS = 183300
      TES(1723)%NHMS = 183500
      TES(1724)%NHMS = 183600
      TES(1725)%NHMS = 183800
      TES(1726)%NHMS = 184000
      TES(1727)%NHMS = 184200
      TES(1728)%NHMS = 184300
      TES(1729)%NHMS = 184400
      TES(1730)%NHMS = 184600
      TES(1731)%NHMS = 184700
      TES(1732)%NHMS = 184900
      TES(1733)%NHMS = 190600
      TES(1734)%NHMS = 190600
      TES(1735)%NHMS = 190700
      TES(1736)%NHMS = 190700
      TES(1737)%NHMS = 190800
      TES(1738)%NHMS = 190800
      TES(1739)%NHMS = 190900
      TES(1740)%NHMS = 190900
      TES(1741)%NHMS = 191200
      TES(1742)%NHMS = 191400
      TES(1743)%NHMS = 201700
      TES(1744)%NHMS = 201700
      TES(1745)%NHMS = 201800
      TES(1746)%NHMS = 201900
      TES(1747)%NHMS = 202200
      TES(1748)%NHMS = 202300
      TES(1749)%NHMS = 202400
      TES(1750)%NHMS = 202500
      TES(1751)%NHMS = 202500
      TES(1752)%NHMS = 202500
      TES(1753)%NHMS = 202600
      TES(1754)%NHMS = 202700
      TES(1755)%NHMS = 202700
      TES(1756)%NHMS = 202800
      TES(1757)%NHMS = 202800
      TES(1758)%NHMS = 204500
      TES(1759)%NHMS = 204500
      TES(1760)%NHMS = 204500
      TES(1761)%NHMS = 204600
      TES(1762)%NHMS = 204700
      TES(1763)%NHMS = 204700
      TES(1764)%NHMS = 204800
      TES(1765)%NHMS = 204900
      TES(1766)%NHMS = 204900
      TES(1767)%NHMS = 205000
      TES(1768)%NHMS = 205200
      TES(1769)%NHMS = 205400
      TES(1770)%NHMS = 205400
      TES(1771)%NHMS = 222400
      TES(1772)%NHMS = 222400
      TES(1773)%NHMS = 222500
      TES(1774)%NHMS = 222600
      TES(1775)%NHMS = 222600
      TES(1776)%NHMS = 222600
      TES(1777)%NHMS = 222700
      TES(1778)%NHMS = 222700
      TES(1779)%NHMS = 222800
      TES(1780)%NHMS = 222800
      TES(1781)%NHMS = 222900
      TES(1782)%NHMS = 223200
      TES(1783)%NHMS = 223800
      TES(1784)%NHMS = 223900
      TES(1785)%NHMS = 000300
      TES(1786)%NHMS = 000300
      TES(1787)%NHMS = 000400
      TES(1788)%NHMS = 000500
      TES(1789)%NHMS = 000500
      TES(1790)%NHMS = 000600
      TES(1791)%NHMS = 000600
      TES(1792)%NHMS = 000800
      TES(1793)%NHMS = 001100
      TES(1794)%NHMS = 001100
      TES(1795)%NHMS = 001300
      TES(1796)%NHMS = 001500
      TES(1797)%NHMS = 001600
      TES(1798)%NHMS = 001700
      TES(1799)%NHMS = 001800
      TES(1800)%NHMS = 001900
      TES(1801)%NHMS = 001900
      TES(1802)%NHMS = 001900
      TES(1803)%NHMS = 002000
      TES(1804)%NHMS = 002100
      TES(1805)%NHMS = 002200
      TES(1806)%NHMS = 002200
      TES(1807)%NHMS = 002500
      TES(1808)%NHMS = 002500
      TES(1809)%NHMS = 002600
      TES(1810)%NHMS = 002700
      TES(1811)%NHMS = 014100
      TES(1812)%NHMS = 014200
      TES(1813)%NHMS = 014300
      TES(1814)%NHMS = 014300
      TES(1815)%NHMS = 014400
      TES(1816)%NHMS = 014400
      TES(1817)%NHMS = 014800
      TES(1818)%NHMS = 015300
      TES(1819)%NHMS = 015500
      TES(1820)%NHMS = 015600
      TES(1821)%NHMS = 041700
      TES(1822)%NHMS = 041700
      TES(1823)%NHMS = 041700
      TES(1824)%NHMS = 041800
      TES(1825)%NHMS = 041900
      TES(1826)%NHMS = 042000
      TES(1827)%NHMS = 042100
      TES(1828)%NHMS = 042100
      TES(1829)%NHMS = 042200
      TES(1830)%NHMS = 042200
      TES(1831)%NHMS = 042200
      TES(1832)%NHMS = 042300
      TES(1833)%NHMS = 042500
      TES(1834)%NHMS = 042500
      TES(1835)%NHMS = 043800
      TES(1836)%NHMS = 043800
      TES(1837)%NHMS = 044200
      TES(1838)%NHMS = 044200
      TES(1839)%NHMS = 051600
      TES(1840)%NHMS = 051700
      TES(1841)%NHMS = 051800
      TES(1842)%NHMS = 051800
      TES(1843)%NHMS = 051900
      TES(1844)%NHMS = 051900
      TES(1845)%NHMS = 052000
      TES(1846)%NHMS = 052000
      TES(1847)%NHMS = 052100
      TES(1848)%NHMS = 052100
      TES(1849)%NHMS = 052200
      TES(1850)%NHMS = 052300
      TES(1851)%NHMS = 052400
      TES(1852)%NHMS = 055700
      TES(1853)%NHMS = 060300
      TES(1854)%NHMS = 060400
      TES(1855)%NHMS = 060500
      TES(1856)%NHMS = 060600
      TES(1857)%NHMS = 060600
      TES(1858)%NHMS = 061200
      TES(1859)%NHMS = 061300
      TES(1860)%NHMS = 061400
      TES(1861)%NHMS = 061900
      TES(1862)%NHMS = 061900
      TES(1863)%NHMS = 062000
      TES(1864)%NHMS = 062000
      TES(1865)%NHMS = 062100
      TES(1866)%NHMS = 064000
      TES(1867)%NHMS = 064000
      TES(1868)%NHMS = 064200
      TES(1869)%NHMS = 064900
      TES(1870)%NHMS = 065200
      TES(1871)%NHMS = 065200
      TES(1872)%NHMS = 065300
      TES(1873)%NHMS = 065300
      TES(1874)%NHMS = 065400
      TES(1875)%NHMS = 075000
      TES(1876)%NHMS = 075100
      TES(1877)%NHMS = 075500
      TES(1878)%NHMS = 075700
      TES(1879)%NHMS = 075800
      TES(1880)%NHMS = 075900
      TES(1881)%NHMS = 080000
      TES(1882)%NHMS = 080000
      TES(1883)%NHMS = 080000
      TES(1884)%NHMS = 082100
      TES(1885)%NHMS = 082100
      TES(1886)%NHMS = 082200
      TES(1887)%NHMS = 082200
      TES(1888)%NHMS = 082400
      TES(1889)%NHMS = 082400
      TES(1890)%NHMS = 082500
      TES(1891)%NHMS = 082500
      TES(1892)%NHMS = 082800
      TES(1893)%NHMS = 082800
      TES(1894)%NHMS = 082900
      TES(1895)%NHMS = 093300
      TES(1896)%NHMS = 093500
      TES(1897)%NHMS = 093700
      TES(1898)%NHMS = 093900
      TES(1899)%NHMS = 093900
      TES(1900)%NHMS = 095600
      TES(1901)%NHMS = 095600
      TES(1902)%NHMS = 095700
      TES(1903)%NHMS = 095700
      TES(1904)%NHMS = 095700
      TES(1905)%NHMS = 095800
      TES(1906)%NHMS = 095800
      TES(1907)%NHMS = 095900
      TES(1908)%NHMS = 095900
      TES(1909)%NHMS = 100000
      TES(1910)%NHMS = 100000
      TES(1911)%NHMS = 100100
      TES(1912)%NHMS = 100100
      TES(1913)%NHMS = 100200
      TES(1914)%NHMS = 100300
      TES(1915)%NHMS = 105600
      TES(1916)%NHMS = 110200
      TES(1917)%NHMS = 110300
      TES(1918)%NHMS = 110300
      TES(1919)%NHMS = 110600
      TES(1920)%NHMS = 110600
      TES(1921)%NHMS = 111500
      TES(1922)%NHMS = 111600
      TES(1923)%NHMS = 111700
      TES(1924)%NHMS = 111800
      TES(1925)%NHMS = 111800
      TES(1926)%NHMS = 101700
      TES(1927)%NHMS = 101800
      TES(1928)%NHMS = 102000
      TES(1929)%NHMS = 102100
      TES(1930)%NHMS = 102100
      TES(1931)%NHMS = 102200
      TES(1932)%NHMS = 102200
      TES(1933)%NHMS = 102300
      TES(1934)%NHMS = 104000
      TES(1935)%NHMS = 104000
      TES(1936)%NHMS = 104100
      TES(1937)%NHMS = 104100
      TES(1938)%NHMS = 104100
      TES(1939)%NHMS = 104200
      TES(1940)%NHMS = 113900
      TES(1941)%NHMS = 113900
      TES(1942)%NHMS = 114000
      TES(1943)%NHMS = 114100
      TES(1944)%NHMS = 114100
      TES(1945)%NHMS = 114200
      TES(1946)%NHMS = 114400
      TES(1947)%NHMS = 114500
      TES(1948)%NHMS = 114500
      TES(1949)%NHMS = 114700
      TES(1950)%NHMS = 114800
      TES(1951)%NHMS = 115000
      TES(1952)%NHMS = 115500
      TES(1953)%NHMS = 115600
      TES(1954)%NHMS = 115600
      TES(1955)%NHMS = 115700
      TES(1956)%NHMS = 115700
      TES(1957)%NHMS = 115800
      TES(1958)%NHMS = 115900
      TES(1959)%NHMS = 120000
      TES(1960)%NHMS = 120100
      TES(1961)%NHMS = 120100
      TES(1962)%NHMS = 120200
      TES(1963)%NHMS = 132600
      TES(1964)%NHMS = 132700
      TES(1965)%NHMS = 132800
      TES(1966)%NHMS = 133100
      TES(1967)%NHMS = 133500
      TES(1968)%NHMS = 133600
      TES(1969)%NHMS = 133800
      TES(1970)%NHMS = 133900
      TES(1971)%NHMS = 155800
      TES(1972)%NHMS = 160100
      TES(1973)%NHMS = 160200
      TES(1974)%NHMS = 165600
      TES(1975)%NHMS = 171700
      TES(1976)%NHMS = 171800
      TES(1977)%NHMS = 172200
      TES(1978)%NHMS = 173700
      TES(1979)%NHMS = 173900
      TES(1980)%NHMS = 173900
      TES(1981)%NHMS = 174000
      TES(1982)%NHMS = 174100
      TES(1983)%NHMS = 181200
      TES(1984)%NHMS = 181200
      TES(1985)%NHMS = 181400
      TES(1986)%NHMS = 181500
      TES(1987)%NHMS = 181600
      TES(1988)%NHMS = 181600
      TES(1989)%NHMS = 181700
      TES(1990)%NHMS = 181800
      TES(1991)%NHMS = 181800
      TES(1992)%NHMS = 181900
      TES(1993)%NHMS = 181900
      TES(1994)%NHMS = 182000
      TES(1995)%NHMS = 182000
      TES(1996)%NHMS = 182100
      TES(1997)%NHMS = 182200
      TES(1998)%NHMS = 182600
      TES(1999)%NHMS = 183500
      TES(2000)%NHMS = 183500
      TES(2001)%NHMS = 183600
      TES(2002)%NHMS = 185300
      TES(2003)%NHMS = 185400
      TES(2004)%NHMS = 185400
      TES(2005)%NHMS = 185500
      TES(2006)%NHMS = 185500
      TES(2007)%NHMS = 185600
      TES(2008)%NHMS = 185800
      TES(2009)%NHMS = 190200
      TES(2010)%NHMS = 190300
      TES(2011)%NHMS = 190300
      TES(2012)%NHMS = 190300
      TES(2013)%NHMS = 190900
      TES(2014)%NHMS = 191000
      TES(2015)%NHMS = 191000
      TES(2016)%NHMS = 200400
      TES(2017)%NHMS = 200400
      TES(2018)%NHMS = 200500
      TES(2019)%NHMS = 200700
      TES(2020)%NHMS = 200800
      TES(2021)%NHMS = 200900
      TES(2022)%NHMS = 201000
      TES(2023)%NHMS = 201000
      TES(2024)%NHMS = 201200
      TES(2025)%NHMS = 201300
      TES(2026)%NHMS = 201400
      TES(2027)%NHMS = 201400
      TES(2028)%NHMS = 201500
      TES(2029)%NHMS = 201500
      TES(2030)%NHMS = 201600
      TES(2031)%NHMS = 203200
      TES(2032)%NHMS = 203300
      TES(2033)%NHMS = 203400
      TES(2034)%NHMS = 203400
      TES(2035)%NHMS = 203600
      TES(2036)%NHMS = 203600
      TES(2037)%NHMS = 203700
      TES(2038)%NHMS = 203700
      TES(2039)%NHMS = 203900
      TES(2040)%NHMS = 204100
      TES(2041)%NHMS = 204100
      TES(2042)%NHMS = 204200
      TES(2043)%NHMS = 204300
      TES(2044)%NHMS = 204400
      TES(2045)%NHMS = 215400
      TES(2046)%NHMS = 221100
      TES(2047)%NHMS = 221200
      TES(2048)%NHMS = 221200
      TES(2049)%NHMS = 221300
      TES(2050)%NHMS = 221300
      TES(2051)%NHMS = 221400
      TES(2052)%NHMS = 221400
      TES(2053)%NHMS = 221500
      TES(2054)%NHMS = 221700
      TES(2055)%NHMS = 222000
      TES(2056)%NHMS = 222100
      TES(2057)%NHMS = 222100
      TES(2058)%NHMS = 223200
      TES(2059)%NHMS = 223300
      TES(2060)%NHMS = 223300
      TES(2061)%NHMS = 223400
      TES(2062)%NHMS = 223400
      TES(2063)%NHMS = 223500
      TES(2064)%NHMS = 233300
      TES(2065)%NHMS = 235000
      TES(2066)%NHMS = 235100
      TES(2067)%NHMS = 235200
      TES(2068)%NHMS = 235200
      TES(2069)%NHMS = 235600
      TES(2070)%NHMS = 235600
      TES(2071)%NHMS = 235800
      TES(2072)%NHMS = 235900
      TES(2073)%NHMS = 000000
      TES(2074)%NHMS = 000100
      TES(2075)%NHMS = 000100
      TES(2076)%NHMS = 000200
      TES(2077)%NHMS = 000200
      TES(2078)%NHMS = 000500
      TES(2079)%NHMS = 000600
      TES(2080)%NHMS = 000700
      TES(2081)%NHMS = 000700
      TES(2082)%NHMS = 000700
      TES(2083)%NHMS = 000800
      TES(2084)%NHMS = 000900
      TES(2085)%NHMS = 001000
      TES(2086)%NHMS = 001200
      TES(2087)%NHMS = 001400
      TES(2088)%NHMS = 012900
      TES(2089)%NHMS = 013200
      TES(2090)%NHMS = 013300
      TES(2091)%NHMS = 013500
      TES(2092)%NHMS = 013600
      TES(2093)%NHMS = 014000
      TES(2094)%NHMS = 014200
      TES(2095)%NHMS = 023800
      TES(2096)%NHMS = 040500
      TES(2097)%NHMS = 040500
      TES(2098)%NHMS = 040500
      TES(2099)%NHMS = 040600
      TES(2100)%NHMS = 040800
      TES(2101)%NHMS = 042700
      TES(2102)%NHMS = 042700
      TES(2103)%NHMS = 050300
      TES(2104)%NHMS = 050500
      TES(2105)%NHMS = 050600
      TES(2106)%NHMS = 050600
      TES(2107)%NHMS = 050600
      TES(2108)%NHMS = 050700
      TES(2109)%NHMS = 050800
      TES(2110)%NHMS = 050900
      TES(2111)%NHMS = 051000
      TES(2112)%NHMS = 051100
      TES(2113)%NHMS = 051100
      TES(2114)%NHMS = 051200
      TES(2115)%NHMS = 054400
      TES(2116)%NHMS = 054400
      TES(2117)%NHMS = 054800
      TES(2118)%NHMS = 055100
      TES(2119)%NHMS = 055200
      TES(2120)%NHMS = 055300
      TES(2121)%NHMS = 055400
      TES(2122)%NHMS = 055400
      TES(2123)%NHMS = 055900
      TES(2124)%NHMS = 060000
      TES(2125)%NHMS = 060300
      TES(2126)%NHMS = 060300
      TES(2127)%NHMS = 060600
      TES(2128)%NHMS = 060600
      TES(2129)%NHMS = 064100
      TES(2130)%NHMS = 064200
      TES(2131)%NHMS = 064300
      TES(2132)%NHMS = 064300
      TES(2133)%NHMS = 064400
      TES(2134)%NHMS = 064500
      TES(2135)%NHMS = 064500
      TES(2136)%NHMS = 073800
      TES(2137)%NHMS = 073800
      TES(2138)%NHMS = 073900
      TES(2139)%NHMS = 074400
      TES(2140)%NHMS = 074400
      TES(2141)%NHMS = 074500
      TES(2142)%NHMS = 074500
      TES(2143)%NHMS = 074600
      TES(2144)%NHMS = 074700
      TES(2145)%NHMS = 074800
      TES(2146)%NHMS = 074800
      TES(2147)%NHMS = 080400
      TES(2148)%NHMS = 080600
      TES(2149)%NHMS = 080600
      TES(2150)%NHMS = 080900
      TES(2151)%NHMS = 081100
      TES(2152)%NHMS = 081200
      TES(2153)%NHMS = 091900
      TES(2154)%NHMS = 092200
      TES(2155)%NHMS = 092300
      TES(2156)%NHMS = 092400
      TES(2157)%NHMS = 092500
      TES(2158)%NHMS = 092500
      TES(2159)%NHMS = 092600
      TES(2160)%NHMS = 092600
      TES(2161)%NHMS = 092700
      TES(2162)%NHMS = 094300
      TES(2163)%NHMS = 094400
      TES(2164)%NHMS = 094400
      TES(2165)%NHMS = 094500
      TES(2166)%NHMS = 094500
      TES(2167)%NHMS = 094600
      TES(2168)%NHMS = 094700
      TES(2169)%NHMS = 094700
      TES(2170)%NHMS = 094800
      TES(2171)%NHMS = 094900
      TES(2172)%NHMS = 095000
      TES(2173)%NHMS = 104400
      TES(2174)%NHMS = 104500
      TES(2175)%NHMS = 105000
      TES(2176)%NHMS = 105100
      TES(2177)%NHMS = 105900
      TES(2178)%NHMS = 110000
      TES(2179)%NHMS = 110000
      TES(2180)%NHMS = 110200
      TES(2181)%NHMS = 110300
      TES(2182)%NHMS = 110300
      TES(2183)%NHMS = 110400
      TES(2184)%NHMS = 110500
      TES(2185)%NHMS = 110500
      TES(2186)%NHMS = 110600
      TES(2187)%NHMS = 112300
      TES(2188)%NHMS = 114300
      TES(2189)%NHMS = 100400
      TES(2190)%NHMS = 100500
      TES(2191)%NHMS = 100800
      TES(2192)%NHMS = 100800
      TES(2193)%NHMS = 100900
      TES(2194)%NHMS = 100900
      TES(2195)%NHMS = 100900
      TES(2196)%NHMS = 101000
      TES(2197)%NHMS = 101000
      TES(2198)%NHMS = 102700
      TES(2199)%NHMS = 102700
      TES(2200)%NHMS = 102900
      TES(2201)%NHMS = 102900
      TES(2202)%NHMS = 102900
      TES(2203)%NHMS = 103000
      TES(2204)%NHMS = 103200
      TES(2205)%NHMS = 112800
      TES(2206)%NHMS = 112900
      TES(2207)%NHMS = 113000
      TES(2208)%NHMS = 113100
      TES(2209)%NHMS = 113100
      TES(2210)%NHMS = 113200
      TES(2211)%NHMS = 113200
      TES(2212)%NHMS = 113300
      TES(2213)%NHMS = 113400
      TES(2214)%NHMS = 113400
      TES(2215)%NHMS = 113600
      TES(2216)%NHMS = 113700
      TES(2217)%NHMS = 113900
      TES(2218)%NHMS = 114000
      TES(2219)%NHMS = 114400
      TES(2220)%NHMS = 114400
      TES(2221)%NHMS = 114500
      TES(2222)%NHMS = 114500
      TES(2223)%NHMS = 131100
      TES(2224)%NHMS = 131400
      TES(2225)%NHMS = 132300
      TES(2226)%NHMS = 154700
      TES(2227)%NHMS = 154700
      TES(2228)%NHMS = 162600
      TES(2229)%NHMS = 162800
      TES(2230)%NHMS = 170500
      TES(2231)%NHMS = 171900
      TES(2232)%NHMS = 172300
      TES(2233)%NHMS = 172600
      TES(2234)%NHMS = 172700
      TES(2235)%NHMS = 172700
      TES(2236)%NHMS = 172800
      TES(2237)%NHMS = 180000
      TES(2238)%NHMS = 180000
      TES(2239)%NHMS = 180200
      TES(2240)%NHMS = 180300
      TES(2241)%NHMS = 180300
      TES(2242)%NHMS = 180400
      TES(2243)%NHMS = 180400
      TES(2244)%NHMS = 180500
      TES(2245)%NHMS = 180500
      TES(2246)%NHMS = 180600
      TES(2247)%NHMS = 180600
      TES(2248)%NHMS = 180600
      TES(2249)%NHMS = 180700
      TES(2250)%NHMS = 180700
      TES(2251)%NHMS = 180800
      TES(2252)%NHMS = 180800
      TES(2253)%NHMS = 181000
      TES(2254)%NHMS = 181900
      TES(2255)%NHMS = 182000
      TES(2256)%NHMS = 182000
      TES(2257)%NHMS = 182100
      TES(2258)%NHMS = 182100
      TES(2259)%NHMS = 182100
      TES(2260)%NHMS = 182200
      TES(2261)%NHMS = 184100
      TES(2262)%NHMS = 184200
      TES(2263)%NHMS = 184200
      TES(2264)%NHMS = 184400
      TES(2265)%NHMS = 184400
      TES(2266)%NHMS = 184600
      TES(2267)%NHMS = 184800
      TES(2268)%NHMS = 184900
      TES(2269)%NHMS = 184900
      TES(2270)%NHMS = 185000
      TES(2271)%NHMS = 185400
      TES(2272)%NHMS = 185400
      TES(2273)%NHMS = 185500
      TES(2274)%NHMS = 185700
      TES(2275)%NHMS = 185800
      TES(2276)%NHMS = 185800
      TES(2277)%NHMS = 195600
      TES(2278)%NHMS = 195700
      TES(2279)%NHMS = 200000
      TES(2280)%NHMS = 200200
      TES(2281)%NHMS = 200200
      TES(2282)%NHMS = 200200
      TES(2283)%NHMS = 202000
      TES(2284)%NHMS = 202100
      TES(2285)%NHMS = 202100
      TES(2286)%NHMS = 202200
      TES(2287)%NHMS = 202200
      TES(2288)%NHMS = 202300
      TES(2289)%NHMS = 202700
      TES(2290)%NHMS = 214100
      TES(2291)%NHMS = 215900
      TES(2292)%NHMS = 220000
      TES(2293)%NHMS = 220000
      TES(2294)%NHMS = 220100
      TES(2295)%NHMS = 220200
      TES(2296)%NHMS = 220300
      TES(2297)%NHMS = 220500
      TES(2298)%NHMS = 221000
      TES(2299)%NHMS = 222000
      TES(2300)%NHMS = 233800
      TES(2301)%NHMS = 233800
      TES(2302)%NHMS = 234000
      TES(2303)%NHMS = 234100
      TES(2304)%NHMS = 234100
      TES(2305)%NHMS = 234200
      TES(2306)%NHMS = 234300
      TES(2307)%NHMS = 234300
      TES(2308)%NHMS = 234400
      TES(2309)%NHMS = 234400
      TES(2310)%NHMS = 234500
      TES(2311)%NHMS = 235000
      TES(2312)%NHMS = 235100
      TES(2313)%NHMS = 235200
      TES(2314)%NHMS = 235300
      TES(2315)%NHMS = 235300
      TES(2316)%NHMS = 235500
      TES(2317)%NHMS = 235500
      TES(2318)%NHMS = 235600
      TES(2319)%NHMS = 235600
      TES(2320)%NHMS = 235700
      TES(2321)%NHMS = 235700
      TES(2322)%NHMS = 235800
      TES(2323)%NHMS = 235800
      TES(2324)%NHMS = 000300
      TES(2325)%NHMS = 011900
      TES(2326)%NHMS = 011900
      TES(2327)%NHMS = 012000
      TES(2328)%NHMS = 012000
      TES(2329)%NHMS = 012100
      TES(2330)%NHMS = 012400
      TES(2331)%NHMS = 012500
      TES(2332)%NHMS = 013000
      TES(2333)%NHMS = 023700
      TES(2334)%NHMS = 023800
      TES(2335)%NHMS = 025600
      TES(2336)%NHMS = 030600
      TES(2337)%NHMS = 035400
      TES(2338)%NHMS = 035400
      TES(2339)%NHMS = 035400
      TES(2340)%NHMS = 041100
      TES(2341)%NHMS = 041400
      TES(2342)%NHMS = 041400
      TES(2343)%NHMS = 041400
      TES(2344)%NHMS = 041500
      TES(2345)%NHMS = 041800
      TES(2346)%NHMS = 045100
      TES(2347)%NHMS = 045300
      TES(2348)%NHMS = 045400
      TES(2349)%NHMS = 045400
      TES(2350)%NHMS = 045500
      TES(2351)%NHMS = 045500
      TES(2352)%NHMS = 045500
      TES(2353)%NHMS = 045600
      TES(2354)%NHMS = 045600
      TES(2355)%NHMS = 045700
      TES(2356)%NHMS = 045800
      TES(2357)%NHMS = 045900
      TES(2358)%NHMS = 050000
      TES(2359)%NHMS = 053200
      TES(2360)%NHMS = 053900
      TES(2361)%NHMS = 054700
      TES(2362)%NHMS = 054800
      TES(2363)%NHMS = 054900
      TES(2364)%NHMS = 055000
      TES(2365)%NHMS = 055100
      TES(2366)%NHMS = 055100
      TES(2367)%NHMS = 055100
      TES(2368)%NHMS = 055200
      TES(2369)%NHMS = 055300
      TES(2370)%NHMS = 055500
      TES(2371)%NHMS = 055600
      TES(2372)%NHMS = 055700
      TES(2373)%NHMS = 061500
      TES(2374)%NHMS = 061600
      TES(2375)%NHMS = 062500
      TES(2376)%NHMS = 062700
      TES(2377)%NHMS = 062700
      TES(2378)%NHMS = 062800
      TES(2379)%NHMS = 062900
      TES(2380)%NHMS = 062900
      TES(2381)%NHMS = 063000
      TES(2382)%NHMS = 063100
      TES(2383)%NHMS = 063200
      TES(2384)%NHMS = 063200
      TES(2385)%NHMS = 063200
      TES(2386)%NHMS = 063300
      TES(2387)%NHMS = 072500
      TES(2388)%NHMS = 072600
      TES(2389)%NHMS = 072600
      TES(2390)%NHMS = 073000
      TES(2391)%NHMS = 073100
      TES(2392)%NHMS = 073300
      TES(2393)%NHMS = 073300
      TES(2394)%NHMS = 073300
      TES(2395)%NHMS = 073400
      TES(2396)%NHMS = 073500
      TES(2397)%NHMS = 073500
      TES(2398)%NHMS = 073600
      TES(2399)%NHMS = 075300
      TES(2400)%NHMS = 075400
      TES(2401)%NHMS = 075400
      TES(2402)%NHMS = 075500
      TES(2403)%NHMS = 075600
      TES(2404)%NHMS = 075700
      TES(2405)%NHMS = 075800
      TES(2406)%NHMS = 075800
      TES(2407)%NHMS = 075900
      TES(2408)%NHMS = 075900
      TES(2409)%NHMS = 080000
      TES(2410)%NHMS = 080300
      TES(2411)%NHMS = 080400
      TES(2412)%NHMS = 080400
      TES(2413)%NHMS = 091300
      TES(2414)%NHMS = 091300
      TES(2415)%NHMS = 093100
      TES(2416)%NHMS = 093200
      TES(2417)%NHMS = 093300
      TES(2418)%NHMS = 093300
      TES(2419)%NHMS = 093400
      TES(2420)%NHMS = 093400
      TES(2421)%NHMS = 093500
      TES(2422)%NHMS = 093600
      TES(2423)%NHMS = 093700
      TES(2424)%NHMS = 093800
      TES(2425)%NHMS = 093800
      TES(2426)%NHMS = 093900
      TES(2427)%NHMS = 103800
      TES(2428)%NHMS = 104800
      TES(2429)%NHMS = 105200
      TES(2430)%NHMS = 111000
      TES(2431)%NHMS = 111000
      TES(2432)%NHMS = 095600
      TES(2433)%NHMS = 101400
      TES(2434)%NHMS = 101500
      TES(2435)%NHMS = 101500
      TES(2436)%NHMS = 101700
      TES(2437)%NHMS = 101700
      TES(2438)%NHMS = 101800
      TES(2439)%NHMS = 101900
      TES(2440)%NHMS = 111600
      TES(2441)%NHMS = 111700
      TES(2442)%NHMS = 111800
      TES(2443)%NHMS = 111900
      TES(2444)%NHMS = 111900
      TES(2445)%NHMS = 111900
      TES(2446)%NHMS = 112300
      TES(2447)%NHMS = 112400
      TES(2448)%NHMS = 112400
      TES(2449)%NHMS = 112600
      TES(2450)%NHMS = 112900
      TES(2451)%NHMS = 113000
      TES(2452)%NHMS = 113100
      TES(2453)%NHMS = 113100
      TES(2454)%NHMS = 113200
      TES(2455)%NHMS = 113200
      TES(2456)%NHMS = 113300
      TES(2457)%NHMS = 113300
      TES(2458)%NHMS = 113400
      TES(2459)%NHMS = 113400
      TES(2460)%NHMS = 113400
      TES(2461)%NHMS = 113700
      TES(2462)%NHMS = 125000
      TES(2463)%NHMS = 125000
      TES(2464)%NHMS = 125200
      TES(2465)%NHMS = 125300
      TES(2466)%NHMS = 125300
      TES(2467)%NHMS = 125400
      TES(2468)%NHMS = 125500
      TES(2469)%NHMS = 125500
      TES(2470)%NHMS = 125600
      TES(2471)%NHMS = 125700
      TES(2472)%NHMS = 125900
      TES(2473)%NHMS = 130000
      TES(2474)%NHMS = 130000
      TES(2475)%NHMS = 130000
      TES(2476)%NHMS = 130100
      TES(2477)%NHMS = 130200
      TES(2478)%NHMS = 131000
      TES(2479)%NHMS = 131100
      TES(2480)%NHMS = 131100
      TES(2481)%NHMS = 131200
      TES(2482)%NHMS = 135200
      TES(2483)%NHMS = 140100
      TES(2484)%NHMS = 153400
      TES(2485)%NHMS = 153400
      TES(2486)%NHMS = 153500
      TES(2487)%NHMS = 153700
      TES(2488)%NHMS = 153900
      TES(2489)%NHMS = 170600
      TES(2490)%NHMS = 171000
      TES(2491)%NHMS = 171200
      TES(2492)%NHMS = 171500
      TES(2493)%NHMS = 171500
      TES(2494)%NHMS = 174900
      TES(2495)%NHMS = 174900
      TES(2496)%NHMS = 175100
      TES(2497)%NHMS = 175100
      TES(2498)%NHMS = 175200
      TES(2499)%NHMS = 175300
      TES(2500)%NHMS = 175400
      TES(2501)%NHMS = 175500
      TES(2502)%NHMS = 175500
      TES(2503)%NHMS = 175600
      TES(2504)%NHMS = 175600
      TES(2505)%NHMS = 175700
      TES(2506)%NHMS = 175700
      TES(2507)%NHMS = 175800
      TES(2508)%NHMS = 175800
      TES(2509)%NHMS = 175900
      TES(2510)%NHMS = 180800
      TES(2511)%NHMS = 180900
      TES(2512)%NHMS = 180900
      TES(2513)%NHMS = 181000
      TES(2514)%NHMS = 181000
      TES(2515)%NHMS = 181000
      TES(2516)%NHMS = 181100
      TES(2517)%NHMS = 182900
      TES(2518)%NHMS = 183000
      TES(2519)%NHMS = 183000
      TES(2520)%NHMS = 183000
      TES(2521)%NHMS = 183100
      TES(2522)%NHMS = 183100
      TES(2523)%NHMS = 183200
      TES(2524)%NHMS = 183400
      TES(2525)%NHMS = 183400
      TES(2526)%NHMS = 183500
      TES(2527)%NHMS = 183500
      TES(2528)%NHMS = 183700
      TES(2529)%NHMS = 183800
      TES(2530)%NHMS = 183800
      TES(2531)%NHMS = 183900
      TES(2532)%NHMS = 183900
      TES(2533)%NHMS = 184000
      TES(2534)%NHMS = 184000
      TES(2535)%NHMS = 184100
      TES(2536)%NHMS = 192100
      TES(2537)%NHMS = 192100
      TES(2538)%NHMS = 193900
      TES(2539)%NHMS = 193900
      TES(2540)%NHMS = 194000
      TES(2541)%NHMS = 194400
      TES(2542)%NHMS = 194500
      TES(2543)%NHMS = 194700
      TES(2544)%NHMS = 194700
      TES(2545)%NHMS = 194800
      TES(2546)%NHMS = 194900
      TES(2547)%NHMS = 200800
      TES(2548)%NHMS = 200800
      TES(2549)%NHMS = 200900
      TES(2550)%NHMS = 200900
      TES(2551)%NHMS = 201000
      TES(2552)%NHMS = 201200
      TES(2553)%NHMS = 201300
      TES(2554)%NHMS = 201300
      TES(2555)%NHMS = 201400
      TES(2556)%NHMS = 201500
      TES(2557)%NHMS = 201500
      TES(2558)%NHMS = 201600
      TES(2559)%NHMS = 212300
      TES(2560)%NHMS = 212300
      TES(2561)%NHMS = 212700
      TES(2562)%NHMS = 212700
      TES(2563)%NHMS = 212800
      TES(2564)%NHMS = 212900
      TES(2565)%NHMS = 213000
      TES(2566)%NHMS = 214700
      TES(2567)%NHMS = 214700
      TES(2568)%NHMS = 214800
      TES(2569)%NHMS = 214800
      TES(2570)%NHMS = 214900
      TES(2571)%NHMS = 214900
      TES(2572)%NHMS = 215000
      TES(2573)%NHMS = 215100
      TES(2574)%NHMS = 215200
      TES(2575)%NHMS = 215200
      TES(2576)%NHMS = 215200
      TES(2577)%NHMS = 215500
      TES(2578)%NHMS = 215600
      TES(2579)%NHMS = 232800
      TES(2580)%NHMS = 232900
      TES(2581)%NHMS = 232900
      TES(2582)%NHMS = 233000
      TES(2583)%NHMS = 233000
      TES(2584)%NHMS = 233100
      TES(2585)%NHMS = 233100
      TES(2586)%NHMS = 233800
      TES(2587)%NHMS = 233900
      TES(2588)%NHMS = 233900
      TES(2589)%NHMS = 234000
      TES(2590)%NHMS = 234000
      TES(2591)%NHMS = 234000
      TES(2592)%NHMS = 234100
      TES(2593)%NHMS = 234200
      TES(2594)%NHMS = 234300
      TES(2595)%NHMS = 234400
      TES(2596)%NHMS = 234500
      TES(2597)%NHMS = 234600
      TES(2598)%NHMS = 234700
      TES(2599)%NHMS = 234800
      TES(2600)%NHMS = 234900
      TES(2601)%NHMS = 235000
      TES(2602)%NHMS = 235000
      TES(2603)%NHMS = 235100
      TES(2604)%NHMS = 010500
      TES(2605)%NHMS = 010600
      TES(2606)%NHMS = 010700
      TES(2607)%NHMS = 010700
      TES(2608)%NHMS = 010800
      TES(2609)%NHMS = 011000
      TES(2610)%NHMS = 011000
      TES(2611)%NHMS = 011400
      TES(2612)%NHMS = 011500
      TES(2613)%NHMS = 011500
      TES(2614)%NHMS = 011600
      TES(2615)%NHMS = 011600
      TES(2616)%NHMS = 011800
      TES(2617)%NHMS = 011800
      TES(2618)%NHMS = 011900
      TES(2619)%NHMS = 012100
      TES(2620)%NHMS = 022500
      TES(2621)%NHMS = 022600
      TES(2622)%NHMS = 022600
      TES(2623)%NHMS = 024500
      TES(2624)%NHMS = 024800
      TES(2625)%NHMS = 024800
      TES(2626)%NHMS = 024900
      TES(2627)%NHMS = 025400
      TES(2628)%NHMS = 025500
      TES(2629)%NHMS = 025500
      TES(2630)%NHMS = 025600
      TES(2631)%NHMS = 035300
      TES(2632)%NHMS = 040300
      TES(2633)%NHMS = 040300
      TES(2634)%NHMS = 040300
      TES(2635)%NHMS = 040400
      TES(2636)%NHMS = 040400
      TES(2637)%NHMS = 040500
      TES(2638)%NHMS = 043900
      TES(2639)%NHMS = 044000
      TES(2640)%NHMS = 044100
      TES(2641)%NHMS = 044200
      TES(2642)%NHMS = 044300
      TES(2643)%NHMS = 044300
      TES(2644)%NHMS = 044300
      TES(2645)%NHMS = 044400
      TES(2646)%NHMS = 044500
      TES(2647)%NHMS = 044500
      TES(2648)%NHMS = 052200
      TES(2649)%NHMS = 052300
      TES(2650)%NHMS = 052400
      TES(2651)%NHMS = 052800
      TES(2652)%NHMS = 053000
      TES(2653)%NHMS = 053100
      TES(2654)%NHMS = 053100
      TES(2655)%NHMS = 053200
      TES(2656)%NHMS = 053600
      TES(2657)%NHMS = 053800
      TES(2658)%NHMS = 053900
      TES(2659)%NHMS = 054000
      TES(2660)%NHMS = 054200
      TES(2661)%NHMS = 054400
      TES(2662)%NHMS = 054400
      TES(2663)%NHMS = 054500
      TES(2664)%NHMS = 060400
      TES(2665)%NHMS = 061600
      TES(2666)%NHMS = 061600
      TES(2667)%NHMS = 061700
      TES(2668)%NHMS = 061900
      TES(2669)%NHMS = 062000
      TES(2670)%NHMS = 062100
      TES(2671)%NHMS = 062100
      TES(2672)%NHMS = 062200
      TES(2673)%NHMS = 070700
      TES(2674)%NHMS = 071300
      TES(2675)%NHMS = 071400
      TES(2676)%NHMS = 071500
      TES(2677)%NHMS = 072100
      TES(2678)%NHMS = 072100
      TES(2679)%NHMS = 072100
      TES(2680)%NHMS = 074000
      TES(2681)%NHMS = 074000
      TES(2682)%NHMS = 074100
      TES(2683)%NHMS = 074100
      TES(2684)%NHMS = 074200
      TES(2685)%NHMS = 074200
      TES(2686)%NHMS = 074300
      TES(2687)%NHMS = 074300
      TES(2688)%NHMS = 074400
      TES(2689)%NHMS = 074400
      TES(2690)%NHMS = 074600
      TES(2691)%NHMS = 074700
      TES(2692)%NHMS = 074700
      TES(2693)%NHMS = 074800
      TES(2694)%NHMS = 075100
      TES(2695)%NHMS = 075200
      TES(2696)%NHMS = 075200
      TES(2697)%NHMS = 090100
      TES(2698)%NHMS = 090200
      TES(2699)%NHMS = 090200
      TES(2700)%NHMS = 092100
      TES(2701)%NHMS = 092100
      TES(2702)%NHMS = 092200
      TES(2703)%NHMS = 092200
      TES(2704)%NHMS = 092300
      TES(2705)%NHMS = 092300
      TES(2706)%NHMS = 092400
      TES(2707)%NHMS = 092400
      TES(2708)%NHMS = 092500
      TES(2709)%NHMS = 092600
      TES(2710)%NHMS = 092600
      TES(2711)%NHMS = 092700
      TES(2712)%NHMS = 092800
      TES(2713)%NHMS = 102700
      TES(2714)%NHMS = 103600
      TES(2715)%NHMS = 103700
      TES(2716)%NHMS = 103800
      TES(2717)%NHMS = 103800
      TES(2718)%NHMS = 103900
      TES(2719)%NHMS = 105800
      TES(2720)%NHMS = 105900
      TES(2721)%NHMS = 105900
      TES(2722)%NHMS = 105900
      TES(2723)%NHMS = 110000
      TES(2724)%NHMS = 093700
      TES(2725)%NHMS = 094300
      TES(2726)%NHMS = 094300
      TES(2727)%NHMS = 094500
      TES(2728)%NHMS = 100200
      TES(2729)%NHMS = 100300
      TES(2730)%NHMS = 100400
      TES(2731)%NHMS = 100400
      TES(2732)%NHMS = 100500
      TES(2733)%NHMS = 100500
      TES(2734)%NHMS = 100600
      TES(2735)%NHMS = 100700
      TES(2736)%NHMS = 100700
      TES(2737)%NHMS = 100700
      TES(2738)%NHMS = 110200
      TES(2739)%NHMS = 111000
      TES(2740)%NHMS = 111000
      TES(2741)%NHMS = 111100
      TES(2742)%NHMS = 112000
      TES(2743)%NHMS = 112100
      TES(2744)%NHMS = 112100
      TES(2745)%NHMS = 112200
      TES(2746)%NHMS = 112200
      TES(2747)%NHMS = 112300
      TES(2748)%NHMS = 112300
      TES(2749)%NHMS = 123800
      TES(2750)%NHMS = 123800
      TES(2751)%NHMS = 124100
      TES(2752)%NHMS = 124200
      TES(2753)%NHMS = 124300
      TES(2754)%NHMS = 124400
      TES(2755)%NHMS = 124400
      TES(2756)%NHMS = 124400
      TES(2757)%NHMS = 124500
      TES(2758)%NHMS = 124600
      TES(2759)%NHMS = 124600
      TES(2760)%NHMS = 124900
      TES(2761)%NHMS = 124900
      TES(2762)%NHMS = 125100
      TES(2763)%NHMS = 125200
      TES(2764)%NHMS = 130000
      TES(2765)%NHMS = 130100
      TES(2766)%NHMS = 130200
      TES(2767)%NHMS = 130300
      TES(2768)%NHMS = 132100
      TES(2769)%NHMS = 134100
      TES(2770)%NHMS = 142800
      TES(2771)%NHMS = 142900
      TES(2772)%NHMS = 143000
      TES(2773)%NHMS = 152200
      TES(2774)%NHMS = 152300
      TES(2775)%NHMS = 152300
      TES(2776)%NHMS = 152400
      TES(2777)%NHMS = 152400
      TES(2778)%NHMS = 152500
      TES(2779)%NHMS = 152500
      TES(2780)%NHMS = 152600
      TES(2781)%NHMS = 165700
      TES(2782)%NHMS = 165800
      TES(2783)%NHMS = 165800
      TES(2784)%NHMS = 165900
      TES(2785)%NHMS = 170200
      TES(2786)%NHMS = 170300
      TES(2787)%NHMS = 173500
      TES(2788)%NHMS = 173600
      TES(2789)%NHMS = 173600
      TES(2790)%NHMS = 173700
      TES(2791)%NHMS = 174000
      TES(2792)%NHMS = 174100
      TES(2793)%NHMS = 174200
      TES(2794)%NHMS = 174300
      TES(2795)%NHMS = 174300
      TES(2796)%NHMS = 174400
      TES(2797)%NHMS = 174500
      TES(2798)%NHMS = 174500
      TES(2799)%NHMS = 174600
      TES(2800)%NHMS = 174600
      TES(2801)%NHMS = 174700
      TES(2802)%NHMS = 175600
      TES(2803)%NHMS = 175700
      TES(2804)%NHMS = 175700
      TES(2805)%NHMS = 175800
      TES(2806)%NHMS = 175800
      TES(2807)%NHMS = 175900
      TES(2808)%NHMS = 175900
      TES(2809)%NHMS = 180000
      TES(2810)%NHMS = 181900
      TES(2811)%NHMS = 181900
      TES(2812)%NHMS = 182000
      TES(2813)%NHMS = 182500
      TES(2814)%NHMS = 182600
      TES(2815)%NHMS = 190900
      TES(2816)%NHMS = 190900
      TES(2817)%NHMS = 191000
      TES(2818)%NHMS = 192500
      TES(2819)%NHMS = 192600
      TES(2820)%NHMS = 192800
      TES(2821)%NHMS = 193100
      TES(2822)%NHMS = 193100
      TES(2823)%NHMS = 193200
      TES(2824)%NHMS = 193300
      TES(2825)%NHMS = 193300
      TES(2826)%NHMS = 193400
      TES(2827)%NHMS = 193500
      TES(2828)%NHMS = 193500
      TES(2829)%NHMS = 193700
      TES(2830)%NHMS = 193900
      TES(2831)%NHMS = 195800
      TES(2832)%NHMS = 200000
      TES(2833)%NHMS = 200100
      TES(2834)%NHMS = 200100
      TES(2835)%NHMS = 200200
      TES(2836)%NHMS = 200300
      TES(2837)%NHMS = 200400
      TES(2838)%NHMS = 200500
      TES(2839)%NHMS = 211000
      TES(2840)%NHMS = 211100
      TES(2841)%NHMS = 211500
      TES(2842)%NHMS = 211600
      TES(2843)%NHMS = 213500
      TES(2844)%NHMS = 213600
      TES(2845)%NHMS = 213600
      TES(2846)%NHMS = 213600
      TES(2847)%NHMS = 213700
      TES(2848)%NHMS = 213800
      TES(2849)%NHMS = 213800
      TES(2850)%NHMS = 214000
      TES(2851)%NHMS = 214100
      TES(2852)%NHMS = 231300
      TES(2853)%NHMS = 231400
      TES(2854)%NHMS = 231500
      TES(2855)%NHMS = 231600
      TES(2856)%NHMS = 231700
      TES(2857)%NHMS = 231800
      TES(2858)%NHMS = 231900
      TES(2859)%NHMS = 232300
      TES(2860)%NHMS = 232500
      TES(2861)%NHMS = 232800
      TES(2862)%NHMS = 232800
      TES(2863)%NHMS = 232800
      TES(2864)%NHMS = 232900
      TES(2865)%NHMS = 233000
      TES(2866)%NHMS = 233100
      TES(2867)%NHMS = 233100
      TES(2868)%NHMS = 233200
      TES(2869)%NHMS = 233200
      TES(2870)%NHMS = 233300
      TES(2871)%NHMS = 233400
      TES(2872)%NHMS = 233500
      TES(2873)%NHMS = 233500
      TES(2874)%NHMS = 233600
      TES(2875)%NHMS = 233700
      TES(2876)%NHMS = 233800
      TES(2877)%NHMS = 233800
      TES(2878)%NHMS = 233900
      TES(2879)%NHMS = 005300
      TES(2880)%NHMS = 005300
      TES(2881)%NHMS = 005400
      TES(2882)%NHMS = 005500
      TES(2883)%NHMS = 005500
      TES(2884)%NHMS = 005600
      TES(2885)%NHMS = 005600
      TES(2886)%NHMS = 010100
      TES(2887)%NHMS = 010100
      TES(2888)%NHMS = 010200
      TES(2889)%NHMS = 010300
      TES(2890)%NHMS = 010700
      TES(2891)%NHMS = 010700
      TES(2892)%NHMS = 010800
      TES(2893)%NHMS = 010900
      TES(2894)%NHMS = 021400
      TES(2895)%NHMS = 023200
      TES(2896)%NHMS = 023300
      TES(2897)%NHMS = 023400
      TES(2898)%NHMS = 023600
      TES(2899)%NHMS = 023600
      TES(2900)%NHMS = 023600
      TES(2901)%NHMS = 023800
      TES(2902)%NHMS = 024100
      TES(2903)%NHMS = 024200
      TES(2904)%NHMS = 024200
      TES(2905)%NHMS = 024300
      TES(2906)%NHMS = 024400
      TES(2907)%NHMS = 024400
      TES(2908)%NHMS = 024500
      TES(2909)%NHMS = 035000
      TES(2910)%NHMS = 035100
      TES(2911)%NHMS = 035200
      TES(2912)%NHMS = 042700
      TES(2913)%NHMS = 042800
      TES(2914)%NHMS = 042900
      TES(2915)%NHMS = 043000
      TES(2916)%NHMS = 043100
      TES(2917)%NHMS = 043200
      TES(2918)%NHMS = 050700
      TES(2919)%NHMS = 051500
      TES(2920)%NHMS = 051900
      TES(2921)%NHMS = 052400
      TES(2922)%NHMS = 052500
      TES(2923)%NHMS = 052900
      TES(2924)%NHMS = 052900
      TES(2925)%NHMS = 052900
      TES(2926)%NHMS = 053000
      TES(2927)%NHMS = 053100
      TES(2928)%NHMS = 053200
      TES(2929)%NHMS = 060100
      TES(2930)%NHMS = 060200
      TES(2931)%NHMS = 060300
      TES(2932)%NHMS = 060300
      TES(2933)%NHMS = 060400
      TES(2934)%NHMS = 060400
      TES(2935)%NHMS = 060500
      TES(2936)%NHMS = 060600
      TES(2937)%NHMS = 060600
      TES(2938)%NHMS = 060900
      TES(2939)%NHMS = 060900
      TES(2940)%NHMS = 061000
      TES(2941)%NHMS = 061100
      TES(2942)%NHMS = 061100
      TES(2943)%NHMS = 061300
      TES(2944)%NHMS = 061300
      TES(2945)%NHMS = 061400
      TES(2946)%NHMS = 065300
      TES(2947)%NHMS = 065400
      TES(2948)%NHMS = 065500
      TES(2949)%NHMS = 065900
      TES(2950)%NHMS = 070000
      TES(2951)%NHMS = 070000
      TES(2952)%NHMS = 070200
      TES(2953)%NHMS = 070600
      TES(2954)%NHMS = 070700
      TES(2955)%NHMS = 070800
      TES(2956)%NHMS = 070900
      TES(2957)%NHMS = 071000
      TES(2958)%NHMS = 072700
      TES(2959)%NHMS = 072800
      TES(2960)%NHMS = 073100
      TES(2961)%NHMS = 073200
      TES(2962)%NHMS = 073200
      TES(2963)%NHMS = 073300
      TES(2964)%NHMS = 073400
      TES(2965)%NHMS = 073500
      TES(2966)%NHMS = 073500
      TES(2967)%NHMS = 073500
      TES(2968)%NHMS = 073800
      TES(2969)%NHMS = 074000
      TES(2970)%NHMS = 074000
      TES(2971)%NHMS = 083600
      TES(2972)%NHMS = 083700
      TES(2973)%NHMS = 083700
      TES(2974)%NHMS = 083800
      TES(2975)%NHMS = 084100
      TES(2976)%NHMS = 084100
      TES(2977)%NHMS = 084500
      TES(2978)%NHMS = 085000
      TES(2979)%NHMS = 090700
      TES(2980)%NHMS = 090800
      TES(2981)%NHMS = 090900
      TES(2982)%NHMS = 090900
      TES(2983)%NHMS = 091000
      TES(2984)%NHMS = 091000
      TES(2985)%NHMS = 091100
      TES(2986)%NHMS = 091100
      TES(2987)%NHMS = 091200
      TES(2988)%NHMS = 091200
      TES(2989)%NHMS = 091300
      TES(2990)%NHMS = 091400
      TES(2991)%NHMS = 091500
      TES(2992)%NHMS = 101700
      TES(2993)%NHMS = 102100
      TES(2994)%NHMS = 102200
      TES(2995)%NHMS = 102400
      TES(2996)%NHMS = 102600
      TES(2997)%NHMS = 102700
      TES(2998)%NHMS = 102700
      TES(2999)%NHMS = 102700
      TES(3000)%NHMS = 102800
      TES(3001)%NHMS = 102900
      TES(3002)%NHMS = 104600
      TES(3003)%NHMS = 104600
      TES(3004)%NHMS = 104700
      TES(3005)%NHMS = 104700
      TES(3006)%NHMS = 104800
      TES(3007)%NHMS = 093000
      TES(3008)%NHMS = 093100
      TES(3009)%NHMS = 093200
      TES(3010)%NHMS = 093200
      TES(3011)%NHMS = 093300
      TES(3012)%NHMS = 093300
      TES(3013)%NHMS = 095000
      TES(3014)%NHMS = 095000
      TES(3015)%NHMS = 095000
      TES(3016)%NHMS = 095100
      TES(3017)%NHMS = 095100
      TES(3018)%NHMS = 095200
      TES(3019)%NHMS = 095200
      TES(3020)%NHMS = 095300
      TES(3021)%NHMS = 095300
      TES(3022)%NHMS = 095400
      TES(3023)%NHMS = 095400
      TES(3024)%NHMS = 095500
      TES(3025)%NHMS = 095500
      TES(3026)%NHMS = 095600
      TES(3027)%NHMS = 095600
      TES(3028)%NHMS = 095700
      TES(3029)%NHMS = 105600
      TES(3030)%NHMS = 105700
      TES(3031)%NHMS = 105700
      TES(3032)%NHMS = 110700
      TES(3033)%NHMS = 111100
      TES(3034)%NHMS = 111200
      TES(3035)%NHMS = 112800
      TES(3036)%NHMS = 112900
      TES(3037)%NHMS = 122500
      TES(3038)%NHMS = 122600
      TES(3039)%NHMS = 122600
      TES(3040)%NHMS = 122700
      TES(3041)%NHMS = 122700
      TES(3042)%NHMS = 122700
      TES(3043)%NHMS = 122800
      TES(3044)%NHMS = 122900
      TES(3045)%NHMS = 123000
      TES(3046)%NHMS = 123100
      TES(3047)%NHMS = 123100
      TES(3048)%NHMS = 123200
      TES(3049)%NHMS = 123200
      TES(3050)%NHMS = 123300
      TES(3051)%NHMS = 123300
      TES(3052)%NHMS = 123500
      TES(3053)%NHMS = 123600
      TES(3054)%NHMS = 123600
      TES(3055)%NHMS = 123700
      TES(3056)%NHMS = 123700
      TES(3057)%NHMS = 123700
      TES(3058)%NHMS = 123900
      TES(3059)%NHMS = 132900
      TES(3060)%NHMS = 133400
      TES(3061)%NHMS = 133500
      TES(3062)%NHMS = 141400
      TES(3063)%NHMS = 141900
      TES(3064)%NHMS = 164200
      TES(3065)%NHMS = 164500
      TES(3066)%NHMS = 164500
      TES(3067)%NHMS = 164600
      TES(3068)%NHMS = 165000
      TES(3069)%NHMS = 172400
      TES(3070)%NHMS = 172600
      TES(3071)%NHMS = 172600
      TES(3072)%NHMS = 172700
      TES(3073)%NHMS = 172700
      TES(3074)%NHMS = 172700
      TES(3075)%NHMS = 172800
      TES(3076)%NHMS = 172800
      TES(3077)%NHMS = 172900
      TES(3078)%NHMS = 173100
      TES(3079)%NHMS = 173100
      TES(3080)%NHMS = 173200
      TES(3081)%NHMS = 173200
      TES(3082)%NHMS = 174500
      TES(3083)%NHMS = 174700
      TES(3084)%NHMS = 174700
      TES(3085)%NHMS = 174700
      TES(3086)%NHMS = 180400
      TES(3087)%NHMS = 180400
      TES(3088)%NHMS = 180500
      TES(3089)%NHMS = 180500
      TES(3090)%NHMS = 180600
      TES(3091)%NHMS = 180700
      TES(3092)%NHMS = 181000
      TES(3093)%NHMS = 181200
      TES(3094)%NHMS = 181200
      TES(3095)%NHMS = 182100
      TES(3096)%NHMS = 182100
      TES(3097)%NHMS = 185700
      TES(3098)%NHMS = 185700
      TES(3099)%NHMS = 185800
      TES(3100)%NHMS = 185800
      TES(3101)%NHMS = 185800
      TES(3102)%NHMS = 185900
      TES(3103)%NHMS = 185900
      TES(3104)%NHMS = 190000
      TES(3105)%NHMS = 190800
      TES(3106)%NHMS = 190800
      TES(3107)%NHMS = 190900
      TES(3108)%NHMS = 191200
      TES(3109)%NHMS = 191400
      TES(3110)%NHMS = 191800
      TES(3111)%NHMS = 191900
      TES(3112)%NHMS = 192100
      TES(3113)%NHMS = 192300
      TES(3114)%NHMS = 192400
      TES(3115)%NHMS = 192500
      TES(3116)%NHMS = 192500
      TES(3117)%NHMS = 194300
      TES(3118)%NHMS = 194500
      TES(3119)%NHMS = 194500
      TES(3120)%NHMS = 194600
      TES(3121)%NHMS = 194900
      TES(3122)%NHMS = 194900
      TES(3123)%NHMS = 195000
      TES(3124)%NHMS = 195100
      TES(3125)%NHMS = 205900
      TES(3126)%NHMS = 205900
      TES(3127)%NHMS = 210000
      TES(3128)%NHMS = 210200
      TES(3129)%NHMS = 210200
      TES(3130)%NHMS = 210300
      TES(3131)%NHMS = 210300
      TES(3132)%NHMS = 210400
      TES(3133)%NHMS = 210400
      TES(3134)%NHMS = 212200
      TES(3135)%NHMS = 212300
      TES(3136)%NHMS = 212400
      TES(3137)%NHMS = 212500
      TES(3138)%NHMS = 212500
      TES(3139)%NHMS = 212600
      TES(3140)%NHMS = 212800
      TES(3141)%NHMS = 212900
      TES(3142)%NHMS = 213000
      TES(3143)%NHMS = 213000
      TES(3144)%NHMS = 213100
      TES(3145)%NHMS = 230100
      TES(3146)%NHMS = 230200
      TES(3147)%NHMS = 230300
      TES(3148)%NHMS = 230300
      TES(3149)%NHMS = 230700
      TES(3150)%NHMS = 231400
      TES(3151)%NHMS = 231500
      TES(3152)%NHMS = 231600
      TES(3153)%NHMS = 231700
      TES(3154)%NHMS = 231700
      TES(3155)%NHMS = 231700
      TES(3156)%NHMS = 231800
      TES(3157)%NHMS = 231900
      TES(3158)%NHMS = 232000
      TES(3159)%NHMS = 232100
      TES(3160)%NHMS = 232100
      TES(3161)%NHMS = 232200
      TES(3162)%NHMS = 232200
      TES(3163)%NHMS = 001300
      TES(3164)%NHMS = 004000
      TES(3165)%NHMS = 004000
      TES(3166)%NHMS = 004200
      TES(3167)%NHMS = 004300
      TES(3168)%NHMS = 004400
      TES(3169)%NHMS = 004400
      TES(3170)%NHMS = 004500
      TES(3171)%NHMS = 004500
      TES(3172)%NHMS = 004800
      TES(3173)%NHMS = 004900
      TES(3174)%NHMS = 004900
      TES(3175)%NHMS = 005300
      TES(3176)%NHMS = 005300
      TES(3177)%NHMS = 005500
      TES(3178)%NHMS = 005600
      TES(3179)%NHMS = 005600
      TES(3180)%NHMS = 005700
      TES(3181)%NHMS = 020200
      TES(3182)%NHMS = 022100
      TES(3183)%NHMS = 022100
      TES(3184)%NHMS = 022600
      TES(3185)%NHMS = 022900
      TES(3186)%NHMS = 023100
      TES(3187)%NHMS = 023200
      TES(3188)%NHMS = 023200
      TES(3189)%NHMS = 033800
      TES(3190)%NHMS = 041600
      TES(3191)%NHMS = 041600
      TES(3192)%NHMS = 041600
      TES(3193)%NHMS = 041700
      TES(3194)%NHMS = 041800
      TES(3195)%NHMS = 041800
      TES(3196)%NHMS = 041900
      TES(3197)%NHMS = 041900
      TES(3198)%NHMS = 042000
      TES(3199)%NHMS = 045300
      TES(3200)%NHMS = 045400
      TES(3201)%NHMS = 045500
      TES(3202)%NHMS = 045900
      TES(3203)%NHMS = 051300
      TES(3204)%NHMS = 051400
      TES(3205)%NHMS = 051500
      TES(3206)%NHMS = 051700
      TES(3207)%NHMS = 051700
      TES(3208)%NHMS = 051800
      TES(3209)%NHMS = 051900
      TES(3210)%NHMS = 055000
      TES(3211)%NHMS = 055100
      TES(3212)%NHMS = 055100
      TES(3213)%NHMS = 055300
      TES(3214)%NHMS = 055500
      TES(3215)%NHMS = 055600
      TES(3216)%NHMS = 055700
      TES(3217)%NHMS = 055800
      TES(3218)%NHMS = 060000
      TES(3219)%NHMS = 060200
      TES(3220)%NHMS = 060200
      TES(3221)%NHMS = 060300
      TES(3222)%NHMS = 060300
      TES(3223)%NHMS = 060300
      TES(3224)%NHMS = 060400
      TES(3225)%NHMS = 060400
      TES(3226)%NHMS = 060500
      TES(3227)%NHMS = 060500
      TES(3228)%NHMS = 060600
      TES(3229)%NHMS = 064000
      TES(3230)%NHMS = 064600
      TES(3231)%NHMS = 064600
      TES(3232)%NHMS = 064900
      TES(3233)%NHMS = 065200
      TES(3234)%NHMS = 065300
      TES(3235)%NHMS = 065400
      TES(3236)%NHMS = 065400
      TES(3237)%NHMS = 065400
      TES(3238)%NHMS = 065500
      TES(3239)%NHMS = 065600
      TES(3240)%NHMS = 065700
      TES(3241)%NHMS = 065700
      TES(3242)%NHMS = 065800
      TES(3243)%NHMS = 071500
      TES(3244)%NHMS = 071600
      TES(3245)%NHMS = 071700
      TES(3246)%NHMS = 071800
      TES(3247)%NHMS = 071900
      TES(3248)%NHMS = 072000
      TES(3249)%NHMS = 072200
      TES(3250)%NHMS = 072400
      TES(3251)%NHMS = 072900
      TES(3252)%NHMS = 082300
      TES(3253)%NHMS = 082300
      TES(3254)%NHMS = 082500
      TES(3255)%NHMS = 082500
      TES(3256)%NHMS = 082600
      TES(3257)%NHMS = 083200
      TES(3258)%NHMS = 083300
      TES(3259)%NHMS = 083300
      TES(3260)%NHMS = 083500
      TES(3261)%NHMS = 083500
      TES(3262)%NHMS = 083600
      TES(3263)%NHMS = 083700
      TES(3264)%NHMS = 083700
      TES(3265)%NHMS = 085400
      TES(3266)%NHMS = 085500
      TES(3267)%NHMS = 085500
      TES(3268)%NHMS = 085600
      TES(3269)%NHMS = 085700
      TES(3270)%NHMS = 085900
      TES(3271)%NHMS = 085900
      TES(3272)%NHMS = 090200
      TES(3273)%NHMS = 090200
      TES(3274)%NHMS = 090300
      TES(3275)%NHMS = 090300
      TES(3276)%NHMS = 090400
      TES(3277)%NHMS = 100600
      TES(3278)%NHMS = 100900
      TES(3279)%NHMS = 101000
      TES(3280)%NHMS = 101400
      TES(3281)%NHMS = 101500
      TES(3282)%NHMS = 101500
      TES(3283)%NHMS = 101600
      TES(3284)%NHMS = 103300
      TES(3285)%NHMS = 103300
      TES(3286)%NHMS = 103400
      TES(3287)%NHMS = 103500
      TES(3288)%NHMS = 103500
      TES(3289)%NHMS = 103600
      TES(3290)%NHMS = 091900
      TES(3291)%NHMS = 092100
      TES(3292)%NHMS = 093800
      TES(3293)%NHMS = 093800
      TES(3294)%NHMS = 094100
      TES(3295)%NHMS = 094200
      TES(3296)%NHMS = 094300
      TES(3297)%NHMS = 094400
      TES(3298)%NHMS = 104400
      TES(3299)%NHMS = 104500
      TES(3300)%NHMS = 104500
      TES(3301)%NHMS = 104800
      TES(3302)%NHMS = 105400
      TES(3303)%NHMS = 105500
      TES(3304)%NHMS = 105600
      TES(3305)%NHMS = 110000
      TES(3306)%NHMS = 111600
      TES(3307)%NHMS = 111700
      TES(3308)%NHMS = 121300
      TES(3309)%NHMS = 121400
      TES(3310)%NHMS = 121400
      TES(3311)%NHMS = 121500
      TES(3312)%NHMS = 121600
      TES(3313)%NHMS = 121600
      TES(3314)%NHMS = 121700
      TES(3315)%NHMS = 121800
      TES(3316)%NHMS = 121800
      TES(3317)%NHMS = 121900
      TES(3318)%NHMS = 122000
      TES(3319)%NHMS = 122100
      TES(3320)%NHMS = 122100
      TES(3321)%NHMS = 122200
      TES(3322)%NHMS = 122300
      TES(3323)%NHMS = 122700
      TES(3324)%NHMS = 122900
      TES(3325)%NHMS = 123200
      TES(3326)%NHMS = 123400
      TES(3327)%NHMS = 123500
      TES(3328)%NHMS = 123500
      TES(3329)%NHMS = 123600
      TES(3330)%NHMS = 123600
      TES(3331)%NHMS = 123600
      TES(3332)%NHMS = 140300
      TES(3333)%NHMS = 140300
      TES(3334)%NHMS = 140400
      TES(3335)%NHMS = 140500
      TES(3336)%NHMS = 144700
      TES(3337)%NHMS = 163400
      TES(3338)%NHMS = 163400
      TES(3339)%NHMS = 163800
      TES(3340)%NHMS = 171300
      TES(3341)%NHMS = 171400
      TES(3342)%NHMS = 171500
      TES(3343)%NHMS = 171600
      TES(3344)%NHMS = 171700
      TES(3345)%NHMS = 171800
      TES(3346)%NHMS = 171900
      TES(3347)%NHMS = 173100
      TES(3348)%NHMS = 173500
      TES(3349)%NHMS = 173500
      TES(3350)%NHMS = 175200
      TES(3351)%NHMS = 175300
      TES(3352)%NHMS = 175300
      TES(3353)%NHMS = 175400
      TES(3354)%NHMS = 175500
      TES(3355)%NHMS = 175500
      TES(3356)%NHMS = 175600
      TES(3357)%NHMS = 175600
      TES(3358)%NHMS = 175700
      TES(3359)%NHMS = 180900
      TES(3360)%NHMS = 181000
      TES(3361)%NHMS = 184500
      TES(3362)%NHMS = 184600
      TES(3363)%NHMS = 184600
      TES(3364)%NHMS = 184700
      TES(3365)%NHMS = 184700
      TES(3366)%NHMS = 184800
      TES(3367)%NHMS = 184900
      TES(3368)%NHMS = 185400
      TES(3369)%NHMS = 185400
      TES(3370)%NHMS = 185600
      TES(3371)%NHMS = 185600
      TES(3372)%NHMS = 185700
      TES(3373)%NHMS = 185700
      TES(3374)%NHMS = 190700
      TES(3375)%NHMS = 190700
      TES(3376)%NHMS = 190800
      TES(3377)%NHMS = 190800
      TES(3378)%NHMS = 190900
      TES(3379)%NHMS = 191200
      TES(3380)%NHMS = 191300
      TES(3381)%NHMS = 193100
      TES(3382)%NHMS = 193100
      TES(3383)%NHMS = 193200
      TES(3384)%NHMS = 193200
      TES(3385)%NHMS = 193300
      TES(3386)%NHMS = 193400
      TES(3387)%NHMS = 193400
      TES(3388)%NHMS = 193600
      TES(3389)%NHMS = 193700
      TES(3390)%NHMS = 193700
      TES(3391)%NHMS = 193700
      TES(3392)%NHMS = 194400
      TES(3393)%NHMS = 204500
      TES(3394)%NHMS = 204600
      TES(3395)%NHMS = 204900
      TES(3396)%NHMS = 204900
      TES(3397)%NHMS = 205000
      TES(3398)%NHMS = 205000
      TES(3399)%NHMS = 205100
      TES(3400)%NHMS = 205100
      TES(3401)%NHMS = 205200
      TES(3402)%NHMS = 205300
      TES(3403)%NHMS = 210900
      TES(3404)%NHMS = 211000
      TES(3405)%NHMS = 211000
      TES(3406)%NHMS = 211100
      TES(3407)%NHMS = 211200
      TES(3408)%NHMS = 211200
      TES(3409)%NHMS = 211300
      TES(3410)%NHMS = 211300
      TES(3411)%NHMS = 211400
      TES(3412)%NHMS = 211400
      TES(3413)%NHMS = 211500
      TES(3414)%NHMS = 211700
      TES(3415)%NHMS = 211700
      TES(3416)%NHMS = 211800
      TES(3417)%NHMS = 224800
      TES(3418)%NHMS = 224900
      TES(3419)%NHMS = 224900
      TES(3420)%NHMS = 225000
      TES(3421)%NHMS = 225000
      TES(3422)%NHMS = 225500
      TES(3423)%NHMS = 225900
      TES(3424)%NHMS = 230000
      TES(3425)%NHMS = 230000
      TES(3426)%NHMS = 230100
      TES(3427)%NHMS = 230200
      TES(3428)%NHMS = 230200
      TES(3429)%NHMS = 230300
      TES(3430)%NHMS = 230300
      TES(3431)%NHMS = 230400
      TES(3432)%NHMS = 230500
      TES(3433)%NHMS = 230500
      TES(3434)%NHMS = 230600
      TES(3435)%NHMS = 230600
      TES(3436)%NHMS = 230700
      TES(3437)%NHMS = 230700
      TES(3438)%NHMS = 230800
      TES(3439)%NHMS = 230800
      TES(3440)%NHMS = 230900
      TES(3441)%NHMS = 231000
      TES(3442)%NHMS = 231100
      TES(3443)%NHMS = 002900
      TES(3444)%NHMS = 002900
      TES(3445)%NHMS = 003000
      TES(3446)%NHMS = 003000
      TES(3447)%NHMS = 003200
      TES(3448)%NHMS = 003200
      TES(3449)%NHMS = 003400
      TES(3450)%NHMS = 003600
      TES(3451)%NHMS = 003700
      TES(3452)%NHMS = 003700
      TES(3453)%NHMS = 003800
      TES(3454)%NHMS = 004100
      TES(3455)%NHMS = 004100
      TES(3456)%NHMS = 004300
      TES(3457)%NHMS = 004300
      TES(3458)%NHMS = 004400
      TES(3459)%NHMS = 004500
      TES(3460)%NHMS = 004500
      TES(3461)%NHMS = 004600
      TES(3462)%NHMS = 004700
      TES(3463)%NHMS = 004800
      TES(3464)%NHMS = 020600
      TES(3465)%NHMS = 020800
      TES(3466)%NHMS = 020800
      TES(3467)%NHMS = 020900
      TES(3468)%NHMS = 021000
      TES(3469)%NHMS = 021000
      TES(3470)%NHMS = 021100
      TES(3471)%NHMS = 021500
      TES(3472)%NHMS = 021600
      TES(3473)%NHMS = 021600
      TES(3474)%NHMS = 021800
      TES(3475)%NHMS = 022000
      TES(3476)%NHMS = 025900
      TES(3477)%NHMS = 030700
      TES(3478)%NHMS = 040300
      TES(3479)%NHMS = 040400
      TES(3480)%NHMS = 044100
      TES(3481)%NHMS = 044200
      TES(3482)%NHMS = 044500
      TES(3483)%NHMS = 044600
      TES(3484)%NHMS = 044600
      TES(3485)%NHMS = 044700
      TES(3486)%NHMS = 045000
      TES(3487)%NHMS = 050200
      TES(3488)%NHMS = 050200
      TES(3489)%NHMS = 050300
      TES(3490)%NHMS = 050300
      TES(3491)%NHMS = 050400
      TES(3492)%NHMS = 050500
      TES(3493)%NHMS = 050500
      TES(3494)%NHMS = 050600
      TES(3495)%NHMS = 054000
      TES(3496)%NHMS = 054000
      TES(3497)%NHMS = 054100
      TES(3498)%NHMS = 054100
      TES(3499)%NHMS = 054200
      TES(3500)%NHMS = 054200
      TES(3501)%NHMS = 054300
      TES(3502)%NHMS = 054300
      TES(3503)%NHMS = 054400
      TES(3504)%NHMS = 054400
      TES(3505)%NHMS = 054500
      TES(3506)%NHMS = 054800
      TES(3507)%NHMS = 054900
      TES(3508)%NHMS = 055000
      TES(3509)%NHMS = 055000
      TES(3510)%NHMS = 055100
      TES(3511)%NHMS = 055100
      TES(3512)%NHMS = 055200
      TES(3513)%NHMS = 055200
      TES(3514)%NHMS = 055200
      TES(3515)%NHMS = 055300
      TES(3516)%NHMS = 062000
      TES(3517)%NHMS = 062100
      TES(3518)%NHMS = 062200
      TES(3519)%NHMS = 062700
      TES(3520)%NHMS = 062900
      TES(3521)%NHMS = 063000
      TES(3522)%NHMS = 063700
      TES(3523)%NHMS = 063700
      TES(3524)%NHMS = 063700
      TES(3525)%NHMS = 063900
      TES(3526)%NHMS = 064000
      TES(3527)%NHMS = 064200
      TES(3528)%NHMS = 064300
      TES(3529)%NHMS = 064400
      TES(3530)%NHMS = 064500
      TES(3531)%NHMS = 064600
      TES(3532)%NHMS = 064600
      TES(3533)%NHMS = 070300
      TES(3534)%NHMS = 070400
      TES(3535)%NHMS = 070500
      TES(3536)%NHMS = 070500
      TES(3537)%NHMS = 070600
      TES(3538)%NHMS = 070700
      TES(3539)%NHMS = 070700
      TES(3540)%NHMS = 070800
      TES(3541)%NHMS = 071300
      TES(3542)%NHMS = 071500
      TES(3543)%NHMS = 081300
      TES(3544)%NHMS = 081400
      TES(3545)%NHMS = 082000
      TES(3546)%NHMS = 082100
      TES(3547)%NHMS = 082100
      TES(3548)%NHMS = 082200
      TES(3549)%NHMS = 082300
      TES(3550)%NHMS = 082400
      TES(3551)%NHMS = 082400
      TES(3552)%NHMS = 082400
      TES(3553)%NHMS = 084200
      TES(3554)%NHMS = 084300
      TES(3555)%NHMS = 084400
      TES(3556)%NHMS = 084500
      TES(3557)%NHMS = 084500
      TES(3558)%NHMS = 084600
      TES(3559)%NHMS = 084700
      TES(3560)%NHMS = 084700
      TES(3561)%NHMS = 084800
      TES(3562)%NHMS = 084900
      TES(3563)%NHMS = 084900
      TES(3564)%NHMS = 085000
      TES(3565)%NHMS = 085100
      TES(3566)%NHMS = 085200
      TES(3567)%NHMS = 085200
      TES(3568)%NHMS = 100400
      TES(3569)%NHMS = 100400
      TES(3570)%NHMS = 102000
      TES(3571)%NHMS = 102100
      TES(3572)%NHMS = 102100
      TES(3573)%NHMS = 102200
      TES(3574)%NHMS = 102200
      TES(3575)%NHMS = 103300
      TES(3576)%NHMS = 103400
      TES(3577)%NHMS = 103500
      TES(3578)%NHMS = 104200
      TES(3579)%NHMS = 104300
      TES(3580)%NHMS = 104500
      TES(3581)%NHMS = 104700
      TES(3582)%NHMS = 110500
      TES(3583)%NHMS = 110500
      TES(3584)%NHMS = 120100
      TES(3585)%NHMS = 120200
      TES(3586)%NHMS = 120300
      TES(3587)%NHMS = 120300
      TES(3588)%NHMS = 120400
      TES(3589)%NHMS = 120400
      TES(3590)%NHMS = 120500
      TES(3591)%NHMS = 120600
      TES(3592)%NHMS = 120700
      TES(3593)%NHMS = 120800
      TES(3594)%NHMS = 120900
      TES(3595)%NHMS = 121000
      TES(3596)%NHMS = 121100
      TES(3597)%NHMS = 121400
      TES(3598)%NHMS = 122100
      TES(3599)%NHMS = 122100
      TES(3600)%NHMS = 122200
      TES(3601)%NHMS = 122300
      TES(3602)%NHMS = 122400
      TES(3603)%NHMS = 122500
      TES(3604)%NHMS = 135000
      TES(3605)%NHMS = 135200
      TES(3606)%NHMS = 135300
      TES(3607)%NHMS = 143600
      TES(3608)%NHMS = 160100
      TES(3609)%NHMS = 160200
      TES(3610)%NHMS = 160200
      TES(3611)%NHMS = 160300
      TES(3612)%NHMS = 162100
      TES(3613)%NHMS = 162200
      TES(3614)%NHMS = 162200
      TES(3615)%NHMS = 162300
      TES(3616)%NHMS = 162600
      TES(3617)%NHMS = 170000
      TES(3618)%NHMS = 170100
      TES(3619)%NHMS = 170100
      TES(3620)%NHMS = 170100
      TES(3621)%NHMS = 170300
      TES(3622)%NHMS = 170500
      TES(3623)%NHMS = 170500
      TES(3624)%NHMS = 170600
      TES(3625)%NHMS = 170600
      TES(3626)%NHMS = 171900
      TES(3627)%NHMS = 172200
      TES(3628)%NHMS = 174100
      TES(3629)%NHMS = 174100
      TES(3630)%NHMS = 174100
      TES(3631)%NHMS = 174200
      TES(3632)%NHMS = 174300
      TES(3633)%NHMS = 174400
      TES(3634)%NHMS = 174400
      TES(3635)%NHMS = 174500
      TES(3636)%NHMS = 174500
      TES(3637)%NHMS = 175100
      TES(3638)%NHMS = 175100
      TES(3639)%NHMS = 175800
      TES(3640)%NHMS = 183300
      TES(3641)%NHMS = 183400
      TES(3642)%NHMS = 183400
      TES(3643)%NHMS = 183500
      TES(3644)%NHMS = 183600
      TES(3645)%NHMS = 183600
      TES(3646)%NHMS = 183700
      TES(3647)%NHMS = 183800
      TES(3648)%NHMS = 183800
      TES(3649)%NHMS = 184100
      TES(3650)%NHMS = 184200
      TES(3651)%NHMS = 184300
      TES(3652)%NHMS = 184500
      TES(3653)%NHMS = 184600
      TES(3654)%NHMS = 185100
      TES(3655)%NHMS = 185300
      TES(3656)%NHMS = 185300
      TES(3657)%NHMS = 185400
      TES(3658)%NHMS = 185600
      TES(3659)%NHMS = 185600
      TES(3660)%NHMS = 185700
      TES(3661)%NHMS = 185700
      TES(3662)%NHMS = 185800
      TES(3663)%NHMS = 185900
      TES(3664)%NHMS = 185900
      TES(3665)%NHMS = 190000
      TES(3666)%NHMS = 190100
      TES(3667)%NHMS = 191800
      TES(3668)%NHMS = 192100
      TES(3669)%NHMS = 192700
      TES(3670)%NHMS = 192800
      TES(3671)%NHMS = 193000
      TES(3672)%NHMS = 202900
      TES(3673)%NHMS = 203000
      TES(3674)%NHMS = 203000
      TES(3675)%NHMS = 203500
      TES(3676)%NHMS = 203600
      TES(3677)%NHMS = 203600
      TES(3678)%NHMS = 203800
      TES(3679)%NHMS = 203800
      TES(3680)%NHMS = 203900
      TES(3681)%NHMS = 203900
      TES(3682)%NHMS = 204000
      TES(3683)%NHMS = 205700
      TES(3684)%NHMS = 205800
      TES(3685)%NHMS = 205900
      TES(3686)%NHMS = 205900
      TES(3687)%NHMS = 210000
      TES(3688)%NHMS = 210100
      TES(3689)%NHMS = 210200
      TES(3690)%NHMS = 210300
      TES(3691)%NHMS = 210400
      TES(3692)%NHMS = 210700
      TES(3693)%NHMS = 223600
      TES(3694)%NHMS = 223600
      TES(3695)%NHMS = 223800
      TES(3696)%NHMS = 223800
      TES(3697)%NHMS = 223900
      TES(3698)%NHMS = 224400
      TES(3699)%NHMS = 224400
      TES(3700)%NHMS = 224700
      TES(3701)%NHMS = 224800
      TES(3702)%NHMS = 225100
      TES(3703)%NHMS = 225600
      TES(3704)%NHMS = 235700
      TES(3705)%NHMS = 001500
      TES(3706)%NHMS = 001500
      TES(3707)%NHMS = 001600
      TES(3708)%NHMS = 002000
      TES(3709)%NHMS = 002000
      TES(3710)%NHMS = 002100
      TES(3711)%NHMS = 002300
      TES(3712)%NHMS = 002400
      TES(3713)%NHMS = 002500
      TES(3714)%NHMS = 003000
      TES(3715)%NHMS = 003200
      TES(3716)%NHMS = 003200
      TES(3717)%NHMS = 003200
      TES(3718)%NHMS = 003300
      TES(3719)%NHMS = 003400
      TES(3720)%NHMS = 003500
      TES(3721)%NHMS = 003700
      TES(3722)%NHMS = 003700
      TES(3723)%NHMS = 015500
      TES(3724)%NHMS = 015600
      TES(3725)%NHMS = 015700
      TES(3726)%NHMS = 015800
      TES(3727)%NHMS = 020100
      TES(3728)%NHMS = 020300
      TES(3729)%NHMS = 020400
      TES(3730)%NHMS = 020500
      TES(3731)%NHMS = 020600
      TES(3732)%NHMS = 020700
      TES(3733)%NHMS = 020800
      TES(3734)%NHMS = 020800
      TES(3735)%NHMS = 020800
      TES(3736)%NHMS = 024800
      TES(3737)%NHMS = 031200
      TES(3738)%NHMS = 042700
      TES(3739)%NHMS = 043200
      TES(3740)%NHMS = 043200
      TES(3741)%NHMS = 043300
      TES(3742)%NHMS = 043300
      TES(3743)%NHMS = 043800
      TES(3744)%NHMS = 044800
      TES(3745)%NHMS = 044900
      TES(3746)%NHMS = 045000
      TES(3747)%NHMS = 045000
      TES(3748)%NHMS = 045100
      TES(3749)%NHMS = 045100
      TES(3750)%NHMS = 045200
      TES(3751)%NHMS = 045400
      TES(3752)%NHMS = 045400
      TES(3753)%NHMS = 052800
      TES(3754)%NHMS = 052900
      TES(3755)%NHMS = 053100
      TES(3756)%NHMS = 053100
      TES(3757)%NHMS = 053200
      TES(3758)%NHMS = 053300
      TES(3759)%NHMS = 053400
      TES(3760)%NHMS = 053400
      TES(3761)%NHMS = 053500
      TES(3762)%NHMS = 053500
      TES(3763)%NHMS = 053600
      TES(3764)%NHMS = 053700
      TES(3765)%NHMS = 053700
      TES(3766)%NHMS = 053800
      TES(3767)%NHMS = 053800
      TES(3768)%NHMS = 053900
      TES(3769)%NHMS = 053900
      TES(3770)%NHMS = 054000
      TES(3771)%NHMS = 054000
      TES(3772)%NHMS = 060800
      TES(3773)%NHMS = 060900
      TES(3774)%NHMS = 061200
      TES(3775)%NHMS = 061500
      TES(3776)%NHMS = 061700
      TES(3777)%NHMS = 061700
      TES(3778)%NHMS = 061800
      TES(3779)%NHMS = 061800
      TES(3780)%NHMS = 062300
      TES(3781)%NHMS = 062400
      TES(3782)%NHMS = 062500
      TES(3783)%NHMS = 062800
      TES(3784)%NHMS = 063300
      TES(3785)%NHMS = 063300
      TES(3786)%NHMS = 063400
      TES(3787)%NHMS = 065100
      TES(3788)%NHMS = 065100
      TES(3789)%NHMS = 070000
      TES(3790)%NHMS = 070400
      TES(3791)%NHMS = 070700
      TES(3792)%NHMS = 070700
      TES(3793)%NHMS = 080400
      TES(3794)%NHMS = 080800
      TES(3795)%NHMS = 080800
      TES(3796)%NHMS = 080900
      TES(3797)%NHMS = 081000
      TES(3798)%NHMS = 081100
      TES(3799)%NHMS = 081200
      TES(3800)%NHMS = 081300
      TES(3801)%NHMS = 083000
      TES(3802)%NHMS = 083100
      TES(3803)%NHMS = 083100
      TES(3804)%NHMS = 083200
      TES(3805)%NHMS = 083400
      TES(3806)%NHMS = 083700
      TES(3807)%NHMS = 083700
      TES(3808)%NHMS = 083800
      TES(3809)%NHMS = 083900
      TES(3810)%NHMS = 083900
      TES(3811)%NHMS = 083900
      TES(3812)%NHMS = 084000
      TES(3813)%NHMS = 094700
      TES(3814)%NHMS = 094900
      TES(3815)%NHMS = 095000
      TES(3816)%NHMS = 095100
      TES(3817)%NHMS = 100800
      TES(3818)%NHMS = 101000
      TES(3819)%NHMS = 101000
      TES(3820)%NHMS = 101000
      TES(3821)%NHMS = 101100
      TES(3822)%NHMS = 101200
      TES(3823)%NHMS = 101200
      TES(3824)%NHMS = 101300
      TES(3825)%NHMS = 111300
      TES(3826)%NHMS = 111400
      TES(3827)%NHMS = 111400
      TES(3828)%NHMS = 111500
      TES(3829)%NHMS = 111600
      TES(3830)%NHMS = 111600
      TES(3831)%NHMS = 111600
      TES(3832)%NHMS = 111700
      TES(3833)%NHMS = 112600
      TES(3834)%NHMS = 112600
      TES(3835)%NHMS = 112800
      TES(3836)%NHMS = 112800
      TES(3837)%NHMS = 112900
      TES(3838)%NHMS = 112900
      TES(3839)%NHMS = 113000
      TES(3840)%NHMS = 114700
      TES(3841)%NHMS = 102100
      TES(3842)%NHMS = 102500
      TES(3843)%NHMS = 102700
      TES(3844)%NHMS = 102800
      TES(3845)%NHMS = 102900
      TES(3846)%NHMS = 103100
      TES(3847)%NHMS = 103100
      TES(3848)%NHMS = 103200
      TES(3849)%NHMS = 103200
      TES(3850)%NHMS = 103300
      TES(3851)%NHMS = 103300
      TES(3852)%NHMS = 103400
      TES(3853)%NHMS = 103400
      TES(3854)%NHMS = 103500
      TES(3855)%NHMS = 105100
      TES(3856)%NHMS = 105200
      TES(3857)%NHMS = 105200
      TES(3858)%NHMS = 105300
      TES(3859)%NHMS = 105300
      TES(3860)%NHMS = 105400
      TES(3861)%NHMS = 115100
      TES(3862)%NHMS = 115100
      TES(3863)%NHMS = 115200
      TES(3864)%NHMS = 115300
      TES(3865)%NHMS = 115400
      TES(3866)%NHMS = 115500
      TES(3867)%NHMS = 115500
      TES(3868)%NHMS = 115500
      TES(3869)%NHMS = 115600
      TES(3870)%NHMS = 115600
      TES(3871)%NHMS = 115900
      TES(3872)%NHMS = 115900
      TES(3873)%NHMS = 120400
      TES(3874)%NHMS = 120500
      TES(3875)%NHMS = 120500
      TES(3876)%NHMS = 120900
      TES(3877)%NHMS = 121000
      TES(3878)%NHMS = 121000
      TES(3879)%NHMS = 121100
      TES(3880)%NHMS = 121200
      TES(3881)%NHMS = 121400
      TES(3882)%NHMS = 133800
      TES(3883)%NHMS = 134100
      TES(3884)%NHMS = 134100
      TES(3885)%NHMS = 154900
      TES(3886)%NHMS = 160800
      TES(3887)%NHMS = 160900
      TES(3888)%NHMS = 160900
      TES(3889)%NHMS = 161000
      TES(3890)%NHMS = 161000
      TES(3891)%NHMS = 161100
      TES(3892)%NHMS = 161300
      TES(3893)%NHMS = 161400
      TES(3894)%NHMS = 164800
      TES(3895)%NHMS = 164900
      TES(3896)%NHMS = 164900
      TES(3897)%NHMS = 164900
      TES(3898)%NHMS = 165000
      TES(3899)%NHMS = 165000
      TES(3900)%NHMS = 165100
      TES(3901)%NHMS = 165100
      TES(3902)%NHMS = 165200
      TES(3903)%NHMS = 165200
      TES(3904)%NHMS = 165300
      TES(3905)%NHMS = 170800
      TES(3906)%NHMS = 170800
      TES(3907)%NHMS = 170900
      TES(3908)%NHMS = 171000
      TES(3909)%NHMS = 172700
      TES(3910)%NHMS = 172800
      TES(3911)%NHMS = 172900
      TES(3912)%NHMS = 172900
      TES(3913)%NHMS = 173000
      TES(3914)%NHMS = 173100
      TES(3915)%NHMS = 173300
      TES(3916)%NHMS = 174400
      TES(3917)%NHMS = 174400
      TES(3918)%NHMS = 174900
      TES(3919)%NHMS = 175000
      TES(3920)%NHMS = 182300
      TES(3921)%NHMS = 182400
      TES(3922)%NHMS = 182400
      TES(3923)%NHMS = 182500
      TES(3924)%NHMS = 182500
      TES(3925)%NHMS = 182600
      TES(3926)%NHMS = 182600
      TES(3927)%NHMS = 182800
      TES(3928)%NHMS = 182900
      TES(3929)%NHMS = 183000
      TES(3930)%NHMS = 183000
      TES(3931)%NHMS = 183100
      TES(3932)%NHMS = 183100
      TES(3933)%NHMS = 183200
      TES(3934)%NHMS = 183200
      TES(3935)%NHMS = 183300
      TES(3936)%NHMS = 183500
      TES(3937)%NHMS = 184000
      TES(3938)%NHMS = 184000
      TES(3939)%NHMS = 184400
      TES(3940)%NHMS = 184400
      TES(3941)%NHMS = 184500
      TES(3942)%NHMS = 184500
      TES(3943)%NHMS = 184600
      TES(3944)%NHMS = 184600
      TES(3945)%NHMS = 184600
      TES(3946)%NHMS = 184800
      TES(3947)%NHMS = 190600
      TES(3948)%NHMS = 190600
      TES(3949)%NHMS = 190700
      TES(3950)%NHMS = 190700
      TES(3951)%NHMS = 190800
      TES(3952)%NHMS = 190900
      TES(3953)%NHMS = 191200
      TES(3954)%NHMS = 192100
      TES(3955)%NHMS = 201700
      TES(3956)%NHMS = 201700
      TES(3957)%NHMS = 201800
      TES(3958)%NHMS = 201900
      TES(3959)%NHMS = 201900
      TES(3960)%NHMS = 202000
      TES(3961)%NHMS = 202200
      TES(3962)%NHMS = 202300
      TES(3963)%NHMS = 202500
      TES(3964)%NHMS = 202700
      TES(3965)%NHMS = 202700
      TES(3966)%NHMS = 202700
      TES(3967)%NHMS = 204500
      TES(3968)%NHMS = 204500
      TES(3969)%NHMS = 204600
      TES(3970)%NHMS = 204600
      TES(3971)%NHMS = 204700
      TES(3972)%NHMS = 204700
      TES(3973)%NHMS = 204700
      TES(3974)%NHMS = 204800
      TES(3975)%NHMS = 204900
      TES(3976)%NHMS = 204900
      TES(3977)%NHMS = 205200
      TES(3978)%NHMS = 205400
      TES(3979)%NHMS = 205600
      TES(3980)%NHMS = 205600
      TES(3981)%NHMS = 205700
      TES(3982)%NHMS = 220700
      TES(3983)%NHMS = 222500
      TES(3984)%NHMS = 222600
      TES(3985)%NHMS = 222800
      TES(3986)%NHMS = 222800
      TES(3987)%NHMS = 223000
      TES(3988)%NHMS = 223200
      TES(3989)%NHMS = 223200
      TES(3990)%NHMS = 223700
      TES(3991)%NHMS = 223800
      
      
      TES(1)%FILENAME = TRIM('retv_vars.10642_0018_003.cdf')
      TES(2)%FILENAME = TRIM('retv_vars.10642_0018_004.cdf')
      TES(3)%FILENAME = TRIM('retv_vars.10642_0019_003.cdf')
      TES(4)%FILENAME = TRIM('retv_vars.10642_0020_003.cdf')
      TES(5)%FILENAME = TRIM('retv_vars.10642_0021_002.cdf')
      TES(6)%FILENAME = TRIM('retv_vars.10642_0021_003.cdf')
      TES(7)%FILENAME = TRIM('retv_vars.10642_0021_004.cdf')
      TES(8)%FILENAME = TRIM('retv_vars.10642_0022_002.cdf')
      TES(9)%FILENAME = TRIM('retv_vars.10642_0022_003.cdf')
      TES(10)%FILENAME = TRIM('retv_vars.10642_0022_004.cdf')
      TES(11)%FILENAME = TRIM('retv_vars.10642_0023_002.cdf')
      TES(12)%FILENAME = TRIM('retv_vars.10642_0028_002.cdf')
      TES(13)%FILENAME = TRIM('retv_vars.10642_0028_003.cdf')
      TES(14)%FILENAME = TRIM('retv_vars.10642_0029_002.cdf')
      TES(15)%FILENAME = TRIM('retv_vars.10642_0055_002.cdf')
      TES(16)%FILENAME = TRIM('retv_vars.10642_0057_002.cdf')
      TES(17)%FILENAME = TRIM('retv_vars.10642_0057_003.cdf')
      TES(18)%FILENAME = TRIM('retv_vars.10642_0058_002.cdf')
      TES(19)%FILENAME = TRIM('retv_vars.10642_0060_003.cdf')
      TES(20)%FILENAME = TRIM('retv_vars.10642_0061_002.cdf')
      TES(21)%FILENAME = TRIM('retv_vars.10642_0063_002.cdf')
      TES(22)%FILENAME = TRIM('retv_vars.10642_0064_004.cdf')
      TES(23)%FILENAME = TRIM('retv_vars.10642_0066_003.cdf')
      TES(24)%FILENAME = TRIM('retv_vars.10642_0066_004.cdf')
      TES(25)%FILENAME = TRIM('retv_vars.10642_0067_003.cdf')
      TES(26)%FILENAME = TRIM('retv_vars.10642_0067_004.cdf')
      TES(27)%FILENAME = TRIM('retv_vars.10642_0068_003.cdf')
      TES(28)%FILENAME = TRIM('retv_vars.10642_0069_002.cdf')
      TES(29)%FILENAME = TRIM('retv_vars.10642_0069_003.cdf')
      TES(30)%FILENAME = TRIM('retv_vars.10642_0069_004.cdf')
      TES(31)%FILENAME = TRIM('retv_vars.10642_0070_004.cdf')
      TES(32)%FILENAME = TRIM('retv_vars.10642_0108_004.cdf')
      TES(33)%FILENAME = TRIM('retv_vars.10642_0109_002.cdf')
      TES(34)%FILENAME = TRIM('retv_vars.10642_0109_003.cdf')
      TES(35)%FILENAME = TRIM('retv_vars.10642_0185_002.cdf')
      TES(36)%FILENAME = TRIM('retv_vars.10642_0187_003.cdf')
      TES(37)%FILENAME = TRIM('retv_vars.10642_0187_004.cdf')
      TES(38)%FILENAME = TRIM('retv_vars.10642_0190_003.cdf')
      TES(39)%FILENAME = TRIM('retv_vars.10642_0190_004.cdf')
      TES(40)%FILENAME = TRIM('retv_vars.10642_0199_004.cdf')
      TES(41)%FILENAME = TRIM('retv_vars.10642_0200_002.cdf')
      TES(42)%FILENAME = TRIM('retv_vars.10642_0200_003.cdf')
      TES(43)%FILENAME = TRIM('retv_vars.10642_0200_004.cdf')
      TES(44)%FILENAME = TRIM('retv_vars.10642_0201_002.cdf')
      TES(45)%FILENAME = TRIM('retv_vars.10642_0201_003.cdf')
      TES(46)%FILENAME = TRIM('retv_vars.10642_0201_004.cdf')
      TES(47)%FILENAME = TRIM('retv_vars.10642_0212_003.cdf')
      TES(48)%FILENAME = TRIM('retv_vars.10642_0212_004.cdf')
      TES(49)%FILENAME = TRIM('retv_vars.10642_0213_002.cdf')
      TES(50)%FILENAME = TRIM('retv_vars.10642_0213_003.cdf')
      TES(51)%FILENAME = TRIM('retv_vars.10642_0213_004.cdf')
      TES(52)%FILENAME = TRIM('retv_vars.10642_0220_002.cdf')
      TES(53)%FILENAME = TRIM('retv_vars.10642_0220_003.cdf')
      TES(54)%FILENAME = TRIM('retv_vars.10642_0220_004.cdf')
      TES(55)%FILENAME = TRIM('retv_vars.10642_0224_002.cdf')
      TES(56)%FILENAME = TRIM('retv_vars.10642_0229_002.cdf')
      TES(57)%FILENAME = TRIM('retv_vars.10642_0229_003.cdf')
      TES(58)%FILENAME = TRIM('retv_vars.10642_0229_004.cdf')
      TES(59)%FILENAME = TRIM('retv_vars.10642_0235_003.cdf')
      TES(60)%FILENAME = TRIM('retv_vars.10642_0238_002.cdf')
      TES(61)%FILENAME = TRIM('retv_vars.10642_0243_003.cdf')
      TES(62)%FILENAME = TRIM('retv_vars.10642_0244_002.cdf')
      TES(63)%FILENAME = TRIM('retv_vars.10642_0244_003.cdf')
      TES(64)%FILENAME = TRIM('retv_vars.10642_0244_004.cdf')
      TES(65)%FILENAME = TRIM('retv_vars.10642_0245_004.cdf')
      TES(66)%FILENAME = TRIM('retv_vars.10642_0246_002.cdf')
      TES(67)%FILENAME = TRIM('retv_vars.10642_0246_003.cdf')
      TES(68)%FILENAME = TRIM('retv_vars.10642_0246_004.cdf')
      TES(69)%FILENAME = TRIM('retv_vars.10642_0247_004.cdf')
      TES(70)%FILENAME = TRIM('retv_vars.10642_0248_002.cdf')
      TES(71)%FILENAME = TRIM('retv_vars.10642_0248_003.cdf')
      TES(72)%FILENAME = TRIM('retv_vars.10642_0248_004.cdf')
      TES(73)%FILENAME = TRIM('retv_vars.10642_0249_002.cdf')
      TES(74)%FILENAME = TRIM('retv_vars.10642_0249_003.cdf')
      TES(75)%FILENAME = TRIM('retv_vars.10642_0249_004.cdf')
      TES(76)%FILENAME = TRIM('retv_vars.10642_0250_002.cdf')
      TES(77)%FILENAME = TRIM('retv_vars.10642_0250_003.cdf')
      TES(78)%FILENAME = TRIM('retv_vars.10642_0250_004.cdf')
      TES(79)%FILENAME = TRIM('retv_vars.10642_0251_003.cdf')
      TES(80)%FILENAME = TRIM('retv_vars.10642_0252_004.cdf')
      TES(81)%FILENAME = TRIM('retv_vars.10642_0253_002.cdf')
      TES(82)%FILENAME = TRIM('retv_vars.10642_0254_004.cdf')
      TES(83)%FILENAME = TRIM('retv_vars.10642_0258_002.cdf')
      TES(84)%FILENAME = TRIM('retv_vars.10642_0258_003.cdf')
      TES(85)%FILENAME = TRIM('retv_vars.10642_0258_004.cdf')
      TES(86)%FILENAME = TRIM('retv_vars.10642_0259_003.cdf')
      TES(87)%FILENAME = TRIM('retv_vars.10642_0259_004.cdf')
      TES(88)%FILENAME = TRIM('retv_vars.10642_0261_003.cdf')
      TES(89)%FILENAME = TRIM('retv_vars.10642_0261_004.cdf')
      TES(90)%FILENAME = TRIM('retv_vars.10642_0262_002.cdf')
      TES(91)%FILENAME = TRIM('retv_vars.10642_0267_002.cdf')
      TES(92)%FILENAME = TRIM('retv_vars.10642_0267_003.cdf')
      TES(93)%FILENAME = TRIM('retv_vars.10642_0267_004.cdf')
      TES(94)%FILENAME = TRIM('retv_vars.10642_0268_002.cdf')
      TES(95)%FILENAME = TRIM('retv_vars.10642_0268_003.cdf')
      TES(96)%FILENAME = TRIM('retv_vars.10642_0268_004.cdf')
      TES(97)%FILENAME = TRIM('retv_vars.10642_0269_002.cdf')
      TES(98)%FILENAME = TRIM('retv_vars.10642_0269_004.cdf')
      TES(99)%FILENAME = TRIM('retv_vars.10642_0274_003.cdf')
      TES(100)%FILENAME = TRIM('retv_vars.10642_0275_002.cdf')
      TES(101)%FILENAME = TRIM('retv_vars.10642_0275_003.cdf')
      TES(102)%FILENAME = TRIM('retv_vars.10642_0275_004.cdf')
      TES(103)%FILENAME = TRIM('retv_vars.10642_0277_004.cdf')
      TES(104)%FILENAME = TRIM('retv_vars.10642_0278_002.cdf')
      TES(105)%FILENAME = TRIM('retv_vars.10642_0302_003.cdf')
      TES(106)%FILENAME = TRIM('retv_vars.10642_0302_004.cdf')
      TES(107)%FILENAME = TRIM('retv_vars.10642_0303_004.cdf')
      TES(108)%FILENAME = TRIM('retv_vars.10642_0304_004.cdf')
      TES(109)%FILENAME = TRIM('retv_vars.10642_0305_002.cdf')
      TES(110)%FILENAME = TRIM('retv_vars.10642_0305_004.cdf')
      TES(111)%FILENAME = TRIM('retv_vars.10642_0306_002.cdf')
      TES(112)%FILENAME = TRIM('retv_vars.10642_0306_004.cdf')
      TES(113)%FILENAME = TRIM('retv_vars.10642_0308_002.cdf')
      TES(114)%FILENAME = TRIM('retv_vars.10642_0309_002.cdf')
      TES(115)%FILENAME = TRIM('retv_vars.10642_0309_004.cdf')
      TES(116)%FILENAME = TRIM('retv_vars.10642_0310_002.cdf')
      TES(117)%FILENAME = TRIM('retv_vars.10642_0310_004.cdf')
      TES(118)%FILENAME = TRIM('retv_vars.10642_0315_002.cdf')
      TES(119)%FILENAME = TRIM('retv_vars.10642_0315_003.cdf')
      TES(120)%FILENAME = TRIM('retv_vars.10642_0315_004.cdf')
      TES(121)%FILENAME = TRIM('retv_vars.10642_0317_002.cdf')
      TES(122)%FILENAME = TRIM('retv_vars.10642_0317_003.cdf')
      TES(123)%FILENAME = TRIM('retv_vars.10642_0317_004.cdf')
      TES(124)%FILENAME = TRIM('retv_vars.10642_0318_002.cdf')
      TES(125)%FILENAME = TRIM('retv_vars.10642_0318_003.cdf')
      TES(126)%FILENAME = TRIM('retv_vars.10642_0318_004.cdf')
      TES(127)%FILENAME = TRIM('retv_vars.10642_0321_004.cdf')
      TES(128)%FILENAME = TRIM('retv_vars.10642_0322_003.cdf')
      TES(129)%FILENAME = TRIM('retv_vars.10642_0324_003.cdf')
      TES(130)%FILENAME = TRIM('retv_vars.10642_0358_002.cdf')
      TES(131)%FILENAME = TRIM('retv_vars.10642_0364_002.cdf')
      TES(132)%FILENAME = TRIM('retv_vars.10642_0364_003.cdf')
      TES(133)%FILENAME = TRIM('retv_vars.10642_0367_004.cdf')
      TES(134)%FILENAME = TRIM('retv_vars.10642_0378_004.cdf')
      TES(135)%FILENAME = TRIM('retv_vars.10642_0379_002.cdf')
      TES(136)%FILENAME = TRIM('retv_vars.10642_0379_003.cdf')
      TES(137)%FILENAME = TRIM('retv_vars.10642_0379_004.cdf')
      TES(138)%FILENAME = TRIM('retv_vars.10642_0406_003.cdf')
      TES(139)%FILENAME = TRIM('retv_vars.10642_0407_002.cdf')
      TES(140)%FILENAME = TRIM('retv_vars.10642_0411_002.cdf')
      TES(141)%FILENAME = TRIM('retv_vars.10642_0411_003.cdf')
      TES(142)%FILENAME = TRIM('retv_vars.10642_0411_004.cdf')
      TES(143)%FILENAME = TRIM('retv_vars.10642_0412_002.cdf')
      TES(144)%FILENAME = TRIM('retv_vars.10642_0412_004.cdf')
      TES(145)%FILENAME = TRIM('retv_vars.10642_0413_003.cdf')
      TES(146)%FILENAME = TRIM('retv_vars.10642_0414_004.cdf')
      TES(147)%FILENAME = TRIM('retv_vars.10642_0415_002.cdf')
      TES(148)%FILENAME = TRIM('retv_vars.10642_0415_003.cdf')
      TES(149)%FILENAME = TRIM('retv_vars.10642_0418_004.cdf')
      TES(150)%FILENAME = TRIM('retv_vars.10642_0419_002.cdf')
      TES(151)%FILENAME = TRIM('retv_vars.10642_0419_003.cdf')
      TES(152)%FILENAME = TRIM('retv_vars.10642_0420_002.cdf')
      TES(153)%FILENAME = TRIM('retv_vars.10642_0421_004.cdf')
      TES(154)%FILENAME = TRIM('retv_vars.10642_0422_002.cdf')
      TES(155)%FILENAME = TRIM('retv_vars.10642_0422_004.cdf')
      TES(156)%FILENAME = TRIM('retv_vars.10642_0423_003.cdf')
      TES(157)%FILENAME = TRIM('retv_vars.10642_0423_004.cdf')
      TES(158)%FILENAME = TRIM('retv_vars.10642_0424_002.cdf')
      TES(159)%FILENAME = TRIM('retv_vars.10642_0424_003.cdf')
      TES(160)%FILENAME = TRIM('retv_vars.10642_0424_004.cdf')
      TES(161)%FILENAME = TRIM('retv_vars.10642_0425_002.cdf')
      TES(162)%FILENAME = TRIM('retv_vars.10642_0425_004.cdf')
      TES(163)%FILENAME = TRIM('retv_vars.10642_0426_003.cdf')
      TES(164)%FILENAME = TRIM('retv_vars.10642_0459_002.cdf')
      TES(165)%FILENAME = TRIM('retv_vars.10642_0461_003.cdf')
      TES(166)%FILENAME = TRIM('retv_vars.10642_0462_002.cdf')
      TES(167)%FILENAME = TRIM('retv_vars.10642_0462_003.cdf')
      TES(168)%FILENAME = TRIM('retv_vars.10642_0463_004.cdf')
      TES(169)%FILENAME = TRIM('retv_vars.10642_0468_004.cdf')
      TES(170)%FILENAME = TRIM('retv_vars.10642_0469_002.cdf')
      TES(171)%FILENAME = TRIM('retv_vars.10642_0469_003.cdf')
      TES(172)%FILENAME = TRIM('retv_vars.10642_0469_004.cdf')
      TES(173)%FILENAME = TRIM('retv_vars.10642_0507_003.cdf')
      TES(174)%FILENAME = TRIM('retv_vars.10642_0508_003.cdf')
      TES(175)%FILENAME = TRIM('retv_vars.10642_0513_003.cdf')
      TES(176)%FILENAME = TRIM('retv_vars.10642_0532_003.cdf')
      TES(177)%FILENAME = TRIM('retv_vars.10642_0532_004.cdf')
      TES(178)%FILENAME = TRIM('retv_vars.10642_0533_003.cdf')
      TES(179)%FILENAME = TRIM('retv_vars.10642_0533_004.cdf')
      TES(180)%FILENAME = TRIM('retv_vars.10642_0534_002.cdf')
      TES(181)%FILENAME = TRIM('retv_vars.10642_0534_003.cdf')
      TES(182)%FILENAME = TRIM('retv_vars.10642_0548_002.cdf')
      TES(183)%FILENAME = TRIM('retv_vars.10642_0548_004.cdf')
      TES(184)%FILENAME = TRIM('retv_vars.10642_0549_002.cdf')
      TES(185)%FILENAME = TRIM('retv_vars.10642_0549_004.cdf')
      TES(186)%FILENAME = TRIM('retv_vars.10642_0550_002.cdf')
      TES(187)%FILENAME = TRIM('retv_vars.10642_0567_003.cdf')
      TES(188)%FILENAME = TRIM('retv_vars.10642_0567_004.cdf')
      TES(189)%FILENAME = TRIM('retv_vars.10642_0568_002.cdf')
      TES(190)%FILENAME = TRIM('retv_vars.10642_0568_003.cdf')
      TES(191)%FILENAME = TRIM('retv_vars.10642_0568_004.cdf')
      TES(192)%FILENAME = TRIM('retv_vars.10642_0569_003.cdf')
      TES(193)%FILENAME = TRIM('retv_vars.10642_0569_004.cdf')
      TES(194)%FILENAME = TRIM('retv_vars.10642_0570_002.cdf')
      TES(195)%FILENAME = TRIM('retv_vars.10642_0570_003.cdf')
      TES(196)%FILENAME = TRIM('retv_vars.10642_0572_003.cdf')
      TES(197)%FILENAME = TRIM('retv_vars.10642_0573_002.cdf')
      TES(198)%FILENAME = TRIM('retv_vars.10642_0573_003.cdf')
      TES(199)%FILENAME = TRIM('retv_vars.10642_0573_004.cdf')
      TES(200)%FILENAME = TRIM('retv_vars.10642_0583_003.cdf')
      TES(201)%FILENAME = TRIM('retv_vars.10642_0586_003.cdf')
      TES(202)%FILENAME = TRIM('retv_vars.10642_0587_002.cdf')
      TES(203)%FILENAME = TRIM('retv_vars.10642_0587_004.cdf')
      TES(204)%FILENAME = TRIM('retv_vars.10642_0588_002.cdf')
      TES(205)%FILENAME = TRIM('retv_vars.10642_0596_003.cdf')
      TES(206)%FILENAME = TRIM('retv_vars.10642_0596_004.cdf')
      TES(207)%FILENAME = TRIM('retv_vars.10642_0598_003.cdf')
      TES(208)%FILENAME = TRIM('retv_vars.10642_0598_004.cdf')
      TES(209)%FILENAME = TRIM('retv_vars.10642_0599_002.cdf')
      TES(210)%FILENAME = TRIM('retv_vars.10642_0604_002.cdf')
      TES(211)%FILENAME = TRIM('retv_vars.10642_0604_003.cdf')
      TES(212)%FILENAME = TRIM('retv_vars.10642_0604_004.cdf')
      TES(213)%FILENAME = TRIM('retv_vars.10642_0605_004.cdf')
      TES(214)%FILENAME = TRIM('retv_vars.10642_0614_002.cdf')
      TES(215)%FILENAME = TRIM('retv_vars.10642_0615_003.cdf')
      TES(216)%FILENAME = TRIM('retv_vars.10642_0616_002.cdf')
      TES(217)%FILENAME = TRIM('retv_vars.10642_0616_003.cdf')
      TES(218)%FILENAME = TRIM('retv_vars.10642_0616_004.cdf')
      TES(219)%FILENAME = TRIM('retv_vars.10642_0617_002.cdf')
      TES(220)%FILENAME = TRIM('retv_vars.10642_0640_002.cdf')
      TES(221)%FILENAME = TRIM('retv_vars.10642_0640_004.cdf')
      TES(222)%FILENAME = TRIM('retv_vars.10642_0641_004.cdf')
      TES(223)%FILENAME = TRIM('retv_vars.10642_0643_003.cdf')
      TES(224)%FILENAME = TRIM('retv_vars.10642_0643_004.cdf')
      TES(225)%FILENAME = TRIM('retv_vars.10642_0644_004.cdf')
      TES(226)%FILENAME = TRIM('retv_vars.10642_0645_002.cdf')
      TES(227)%FILENAME = TRIM('retv_vars.10642_0645_004.cdf')
      TES(228)%FILENAME = TRIM('retv_vars.10642_0647_002.cdf')
      TES(229)%FILENAME = TRIM('retv_vars.10642_0652_003.cdf')
      TES(230)%FILENAME = TRIM('retv_vars.10642_0653_002.cdf')
      TES(231)%FILENAME = TRIM('retv_vars.10642_0654_004.cdf')
      TES(232)%FILENAME = TRIM('retv_vars.10642_0655_002.cdf')
      TES(233)%FILENAME = TRIM('retv_vars.10642_0655_003.cdf')
      TES(234)%FILENAME = TRIM('retv_vars.10642_0655_004.cdf')
      TES(235)%FILENAME = TRIM('retv_vars.10642_0656_002.cdf')
      TES(236)%FILENAME = TRIM('retv_vars.10642_0656_003.cdf')
      TES(237)%FILENAME = TRIM('retv_vars.10642_0659_004.cdf')
      TES(238)%FILENAME = TRIM('retv_vars.10642_0690_002.cdf')
      TES(239)%FILENAME = TRIM('retv_vars.10642_0690_003.cdf')
      TES(240)%FILENAME = TRIM('retv_vars.10642_0691_004.cdf')
      TES(241)%FILENAME = TRIM('retv_vars.10642_0692_003.cdf')
      TES(242)%FILENAME = TRIM('retv_vars.10642_0692_004.cdf')
      TES(243)%FILENAME = TRIM('retv_vars.10642_0693_002.cdf')
      TES(244)%FILENAME = TRIM('retv_vars.10642_0693_003.cdf')
      TES(245)%FILENAME = TRIM('retv_vars.10642_0693_004.cdf')
      TES(246)%FILENAME = TRIM('retv_vars.10642_0694_003.cdf')
      TES(247)%FILENAME = TRIM('retv_vars.10642_0694_004.cdf')
      TES(248)%FILENAME = TRIM('retv_vars.10642_0695_002.cdf')
      TES(249)%FILENAME = TRIM('retv_vars.10642_0699_004.cdf')
      TES(250)%FILENAME = TRIM('retv_vars.10642_0700_003.cdf')
      TES(251)%FILENAME = TRIM('retv_vars.10642_0700_004.cdf')
      TES(252)%FILENAME = TRIM('retv_vars.10642_0701_002.cdf')
      TES(253)%FILENAME = TRIM('retv_vars.10642_0701_003.cdf')
      TES(254)%FILENAME = TRIM('retv_vars.10642_0701_004.cdf')
      TES(255)%FILENAME = TRIM('retv_vars.10642_0702_002.cdf')
      TES(256)%FILENAME = TRIM('retv_vars.10642_0702_003.cdf')
      TES(257)%FILENAME = TRIM('retv_vars.10642_0702_004.cdf')
      TES(258)%FILENAME = TRIM('retv_vars.10642_0703_002.cdf')
      TES(259)%FILENAME = TRIM('retv_vars.10642_0703_004.cdf')
      TES(260)%FILENAME = TRIM('retv_vars.10642_0704_003.cdf')
      TES(261)%FILENAME = TRIM('retv_vars.10642_0704_004.cdf')
      TES(262)%FILENAME = TRIM('retv_vars.10642_0728_002.cdf')
      TES(263)%FILENAME = TRIM('retv_vars.10642_0731_004.cdf')
      TES(264)%FILENAME = TRIM('retv_vars.10642_0732_002.cdf')
      TES(265)%FILENAME = TRIM('retv_vars.10642_0733_002.cdf')
      TES(266)%FILENAME = TRIM('retv_vars.10642_0738_003.cdf')
      TES(267)%FILENAME = TRIM('retv_vars.10642_0739_002.cdf')
      TES(268)%FILENAME = TRIM('retv_vars.10642_0739_003.cdf')
      TES(269)%FILENAME = TRIM('retv_vars.10642_0740_003.cdf')
      TES(270)%FILENAME = TRIM('retv_vars.10642_0743_002.cdf')
      TES(271)%FILENAME = TRIM('retv_vars.10642_0747_003.cdf')
      TES(272)%FILENAME = TRIM('retv_vars.10642_0747_004.cdf')
      TES(273)%FILENAME = TRIM('retv_vars.10642_0762_004.cdf')
      TES(274)%FILENAME = TRIM('retv_vars.10647_0018_002.cdf')
      TES(275)%FILENAME = TRIM('retv_vars.10647_0019_002.cdf')
      TES(276)%FILENAME = TRIM('retv_vars.10647_0019_003.cdf')
      TES(277)%FILENAME = TRIM('retv_vars.10647_0021_003.cdf')
      TES(278)%FILENAME = TRIM('retv_vars.10647_0021_004.cdf')
      TES(279)%FILENAME = TRIM('retv_vars.10647_0022_002.cdf')
      TES(280)%FILENAME = TRIM('retv_vars.10647_0022_004.cdf')
      TES(281)%FILENAME = TRIM('retv_vars.10647_0023_002.cdf')
      TES(282)%FILENAME = TRIM('retv_vars.10647_0027_003.cdf')
      TES(283)%FILENAME = TRIM('retv_vars.10647_0028_002.cdf')
      TES(284)%FILENAME = TRIM('retv_vars.10647_0028_003.cdf')
      TES(285)%FILENAME = TRIM('retv_vars.10647_0029_003.cdf')
      TES(286)%FILENAME = TRIM('retv_vars.10647_0058_004.cdf')
      TES(287)%FILENAME = TRIM('retv_vars.10647_0059_002.cdf')
      TES(288)%FILENAME = TRIM('retv_vars.10647_0059_003.cdf')
      TES(289)%FILENAME = TRIM('retv_vars.10647_0063_002.cdf')
      TES(290)%FILENAME = TRIM('retv_vars.10647_0063_004.cdf')
      TES(291)%FILENAME = TRIM('retv_vars.10647_0064_002.cdf')
      TES(292)%FILENAME = TRIM('retv_vars.10647_0067_002.cdf')
      TES(293)%FILENAME = TRIM('retv_vars.10647_0067_003.cdf')
      TES(294)%FILENAME = TRIM('retv_vars.10647_0067_004.cdf')
      TES(295)%FILENAME = TRIM('retv_vars.10647_0068_003.cdf')
      TES(296)%FILENAME = TRIM('retv_vars.10647_0068_004.cdf')
      TES(297)%FILENAME = TRIM('retv_vars.10647_0069_002.cdf')
      TES(298)%FILENAME = TRIM('retv_vars.10647_0070_002.cdf')
      TES(299)%FILENAME = TRIM('retv_vars.10647_0109_003.cdf')
      TES(300)%FILENAME = TRIM('retv_vars.10647_0115_004.cdf')
      TES(301)%FILENAME = TRIM('retv_vars.10647_0118_002.cdf')
      TES(302)%FILENAME = TRIM('retv_vars.10647_0124_003.cdf')
      TES(303)%FILENAME = TRIM('retv_vars.10647_0185_003.cdf')
      TES(304)%FILENAME = TRIM('retv_vars.10647_0187_003.cdf')
      TES(305)%FILENAME = TRIM('retv_vars.10647_0187_004.cdf')
      TES(306)%FILENAME = TRIM('retv_vars.10647_0188_003.cdf')
      TES(307)%FILENAME = TRIM('retv_vars.10647_0189_004.cdf')
      TES(308)%FILENAME = TRIM('retv_vars.10647_0190_004.cdf')
      TES(309)%FILENAME = TRIM('retv_vars.10647_0200_003.cdf')
      TES(310)%FILENAME = TRIM('retv_vars.10647_0200_004.cdf')
      TES(311)%FILENAME = TRIM('retv_vars.10647_0201_003.cdf')
      TES(312)%FILENAME = TRIM('retv_vars.10647_0201_004.cdf')
      TES(313)%FILENAME = TRIM('retv_vars.10647_0202_002.cdf')
      TES(314)%FILENAME = TRIM('retv_vars.10647_0212_004.cdf')
      TES(315)%FILENAME = TRIM('retv_vars.10647_0221_002.cdf')
      TES(316)%FILENAME = TRIM('retv_vars.10647_0231_002.cdf')
      TES(317)%FILENAME = TRIM('retv_vars.10647_0231_004.cdf')
      TES(318)%FILENAME = TRIM('retv_vars.10647_0234_003.cdf')
      TES(319)%FILENAME = TRIM('retv_vars.10647_0237_003.cdf')
      TES(320)%FILENAME = TRIM('retv_vars.10647_0237_004.cdf')
      TES(321)%FILENAME = TRIM('retv_vars.10647_0244_002.cdf')
      TES(322)%FILENAME = TRIM('retv_vars.10647_0244_003.cdf')
      TES(323)%FILENAME = TRIM('retv_vars.10647_0244_004.cdf')
      TES(324)%FILENAME = TRIM('retv_vars.10647_0245_002.cdf')
      TES(325)%FILENAME = TRIM('retv_vars.10647_0245_004.cdf')
      TES(326)%FILENAME = TRIM('retv_vars.10647_0246_002.cdf')
      TES(327)%FILENAME = TRIM('retv_vars.10647_0246_003.cdf')
      TES(328)%FILENAME = TRIM('retv_vars.10647_0248_002.cdf')
      TES(329)%FILENAME = TRIM('retv_vars.10647_0250_002.cdf')
      TES(330)%FILENAME = TRIM('retv_vars.10647_0250_004.cdf')
      TES(331)%FILENAME = TRIM('retv_vars.10647_0251_002.cdf')
      TES(332)%FILENAME = TRIM('retv_vars.10647_0251_003.cdf')
      TES(333)%FILENAME = TRIM('retv_vars.10647_0252_002.cdf')
      TES(334)%FILENAME = TRIM('retv_vars.10647_0252_004.cdf')
      TES(335)%FILENAME = TRIM('retv_vars.10647_0259_002.cdf')
      TES(336)%FILENAME = TRIM('retv_vars.10647_0259_003.cdf')
      TES(337)%FILENAME = TRIM('retv_vars.10647_0260_003.cdf')
      TES(338)%FILENAME = TRIM('retv_vars.10647_0261_003.cdf')
      TES(339)%FILENAME = TRIM('retv_vars.10647_0268_002.cdf')
      TES(340)%FILENAME = TRIM('retv_vars.10647_0268_003.cdf')
      TES(341)%FILENAME = TRIM('retv_vars.10647_0271_004.cdf')
      TES(342)%FILENAME = TRIM('retv_vars.10647_0272_003.cdf')
      TES(343)%FILENAME = TRIM('retv_vars.10647_0273_004.cdf')
      TES(344)%FILENAME = TRIM('retv_vars.10647_0276_004.cdf')
      TES(345)%FILENAME = TRIM('retv_vars.10647_0277_002.cdf')
      TES(346)%FILENAME = TRIM('retv_vars.10647_0279_003.cdf')
      TES(347)%FILENAME = TRIM('retv_vars.10647_0279_004.cdf')
      TES(348)%FILENAME = TRIM('retv_vars.10647_0302_003.cdf')
      TES(349)%FILENAME = TRIM('retv_vars.10647_0305_002.cdf')
      TES(350)%FILENAME = TRIM('retv_vars.10647_0305_003.cdf')
      TES(351)%FILENAME = TRIM('retv_vars.10647_0305_004.cdf')
      TES(352)%FILENAME = TRIM('retv_vars.10647_0306_002.cdf')
      TES(353)%FILENAME = TRIM('retv_vars.10647_0307_002.cdf')
      TES(354)%FILENAME = TRIM('retv_vars.10647_0309_002.cdf')
      TES(355)%FILENAME = TRIM('retv_vars.10647_0310_003.cdf')
      TES(356)%FILENAME = TRIM('retv_vars.10647_0315_002.cdf')
      TES(357)%FILENAME = TRIM('retv_vars.10647_0315_003.cdf')
      TES(358)%FILENAME = TRIM('retv_vars.10647_0315_004.cdf')
      TES(359)%FILENAME = TRIM('retv_vars.10647_0316_002.cdf')
      TES(360)%FILENAME = TRIM('retv_vars.10647_0316_003.cdf')
      TES(361)%FILENAME = TRIM('retv_vars.10647_0316_004.cdf')
      TES(362)%FILENAME = TRIM('retv_vars.10647_0317_004.cdf')
      TES(363)%FILENAME = TRIM('retv_vars.10647_0318_003.cdf')
      TES(364)%FILENAME = TRIM('retv_vars.10647_0319_004.cdf')
      TES(365)%FILENAME = TRIM('retv_vars.10647_0322_002.cdf')
      TES(366)%FILENAME = TRIM('retv_vars.10647_0322_003.cdf')
      TES(367)%FILENAME = TRIM('retv_vars.10647_0323_003.cdf')
      TES(368)%FILENAME = TRIM('retv_vars.10647_0358_002.cdf')
      TES(369)%FILENAME = TRIM('retv_vars.10647_0358_003.cdf')
      TES(370)%FILENAME = TRIM('retv_vars.10647_0358_004.cdf')
      TES(371)%FILENAME = TRIM('retv_vars.10647_0359_002.cdf')
      TES(372)%FILENAME = TRIM('retv_vars.10647_0359_003.cdf')
      TES(373)%FILENAME = TRIM('retv_vars.10647_0364_002.cdf')
      TES(374)%FILENAME = TRIM('retv_vars.10647_0364_003.cdf')
      TES(375)%FILENAME = TRIM('retv_vars.10647_0366_002.cdf')
      TES(376)%FILENAME = TRIM('retv_vars.10647_0368_002.cdf')
      TES(377)%FILENAME = TRIM('retv_vars.10647_0368_004.cdf')
      TES(378)%FILENAME = TRIM('retv_vars.10647_0369_002.cdf')
      TES(379)%FILENAME = TRIM('retv_vars.10647_0378_003.cdf')
      TES(380)%FILENAME = TRIM('retv_vars.10647_0412_002.cdf')
      TES(381)%FILENAME = TRIM('retv_vars.10647_0412_003.cdf')
      TES(382)%FILENAME = TRIM('retv_vars.10647_0412_004.cdf')
      TES(383)%FILENAME = TRIM('retv_vars.10647_0413_002.cdf')
      TES(384)%FILENAME = TRIM('retv_vars.10647_0413_003.cdf')
      TES(385)%FILENAME = TRIM('retv_vars.10647_0414_002.cdf')
      TES(386)%FILENAME = TRIM('retv_vars.10647_0414_004.cdf')
      TES(387)%FILENAME = TRIM('retv_vars.10647_0415_002.cdf')
      TES(388)%FILENAME = TRIM('retv_vars.10647_0415_003.cdf')
      TES(389)%FILENAME = TRIM('retv_vars.10647_0415_004.cdf')
      TES(390)%FILENAME = TRIM('retv_vars.10647_0416_002.cdf')
      TES(391)%FILENAME = TRIM('retv_vars.10647_0416_003.cdf')
      TES(392)%FILENAME = TRIM('retv_vars.10647_0417_003.cdf')
      TES(393)%FILENAME = TRIM('retv_vars.10647_0419_002.cdf')
      TES(394)%FILENAME = TRIM('retv_vars.10647_0419_003.cdf')
      TES(395)%FILENAME = TRIM('retv_vars.10647_0420_004.cdf')
      TES(396)%FILENAME = TRIM('retv_vars.10647_0421_002.cdf')
      TES(397)%FILENAME = TRIM('retv_vars.10647_0421_003.cdf')
      TES(398)%FILENAME = TRIM('retv_vars.10647_0422_002.cdf')
      TES(399)%FILENAME = TRIM('retv_vars.10647_0422_003.cdf')
      TES(400)%FILENAME = TRIM('retv_vars.10647_0422_004.cdf')
      TES(401)%FILENAME = TRIM('retv_vars.10647_0423_002.cdf')
      TES(402)%FILENAME = TRIM('retv_vars.10647_0423_003.cdf')
      TES(403)%FILENAME = TRIM('retv_vars.10647_0423_004.cdf')
      TES(404)%FILENAME = TRIM('retv_vars.10647_0424_002.cdf')
      TES(405)%FILENAME = TRIM('retv_vars.10647_0424_004.cdf')
      TES(406)%FILENAME = TRIM('retv_vars.10647_0425_002.cdf')
      TES(407)%FILENAME = TRIM('retv_vars.10647_0425_003.cdf')
      TES(408)%FILENAME = TRIM('retv_vars.10647_0425_004.cdf')
      TES(409)%FILENAME = TRIM('retv_vars.10647_0429_004.cdf')
      TES(410)%FILENAME = TRIM('retv_vars.10647_0430_002.cdf')
      TES(411)%FILENAME = TRIM('retv_vars.10647_0461_003.cdf')
      TES(412)%FILENAME = TRIM('retv_vars.10647_0461_004.cdf')
      TES(413)%FILENAME = TRIM('retv_vars.10647_0462_004.cdf')
      TES(414)%FILENAME = TRIM('retv_vars.10647_0464_003.cdf')
      TES(415)%FILENAME = TRIM('retv_vars.10647_0465_003.cdf')
      TES(416)%FILENAME = TRIM('retv_vars.10647_0467_003.cdf')
      TES(417)%FILENAME = TRIM('retv_vars.10647_0468_004.cdf')
      TES(418)%FILENAME = TRIM('retv_vars.10647_0493_002.cdf')
      TES(419)%FILENAME = TRIM('retv_vars.10647_0501_004.cdf')
      TES(420)%FILENAME = TRIM('retv_vars.10647_0502_002.cdf')
      TES(421)%FILENAME = TRIM('retv_vars.10647_0507_003.cdf')
      TES(422)%FILENAME = TRIM('retv_vars.10647_0507_004.cdf')
      TES(423)%FILENAME = TRIM('retv_vars.10647_0533_002.cdf')
      TES(424)%FILENAME = TRIM('retv_vars.10647_0533_003.cdf')
      TES(425)%FILENAME = TRIM('retv_vars.10647_0533_004.cdf')
      TES(426)%FILENAME = TRIM('retv_vars.10647_0534_002.cdf')
      TES(427)%FILENAME = TRIM('retv_vars.10647_0537_003.cdf')
      TES(428)%FILENAME = TRIM('retv_vars.10647_0548_004.cdf')
      TES(429)%FILENAME = TRIM('retv_vars.10647_0549_002.cdf')
      TES(430)%FILENAME = TRIM('retv_vars.10647_0549_003.cdf')
      TES(431)%FILENAME = TRIM('retv_vars.10647_0549_004.cdf')
      TES(432)%FILENAME = TRIM('retv_vars.10647_0550_002.cdf')
      TES(433)%FILENAME = TRIM('retv_vars.10647_0550_003.cdf')
      TES(434)%FILENAME = TRIM('retv_vars.10647_0550_004.cdf')
      TES(435)%FILENAME = TRIM('retv_vars.10647_0551_003.cdf')
      TES(436)%FILENAME = TRIM('retv_vars.10647_0567_004.cdf')
      TES(437)%FILENAME = TRIM('retv_vars.10647_0568_003.cdf')
      TES(438)%FILENAME = TRIM('retv_vars.10647_0569_002.cdf')
      TES(439)%FILENAME = TRIM('retv_vars.10647_0569_003.cdf')
      TES(440)%FILENAME = TRIM('retv_vars.10647_0571_002.cdf')
      TES(441)%FILENAME = TRIM('retv_vars.10647_0572_002.cdf')
      TES(442)%FILENAME = TRIM('retv_vars.10647_0572_003.cdf')
      TES(443)%FILENAME = TRIM('retv_vars.10647_0572_004.cdf')
      TES(444)%FILENAME = TRIM('retv_vars.10647_0573_002.cdf')
      TES(445)%FILENAME = TRIM('retv_vars.10647_0573_003.cdf')
      TES(446)%FILENAME = TRIM('retv_vars.10647_0573_004.cdf')
      TES(447)%FILENAME = TRIM('retv_vars.10647_0585_004.cdf')
      TES(448)%FILENAME = TRIM('retv_vars.10647_0589_002.cdf')
      TES(449)%FILENAME = TRIM('retv_vars.10647_0592_003.cdf')
      TES(450)%FILENAME = TRIM('retv_vars.10647_0593_002.cdf')
      TES(451)%FILENAME = TRIM('retv_vars.10647_0593_004.cdf')
      TES(452)%FILENAME = TRIM('retv_vars.10647_0595_002.cdf')
      TES(453)%FILENAME = TRIM('retv_vars.10647_0595_004.cdf')
      TES(454)%FILENAME = TRIM('retv_vars.10647_0596_003.cdf')
      TES(455)%FILENAME = TRIM('retv_vars.10647_0597_003.cdf')
      TES(456)%FILENAME = TRIM('retv_vars.10647_0597_004.cdf')
      TES(457)%FILENAME = TRIM('retv_vars.10647_0598_002.cdf')
      TES(458)%FILENAME = TRIM('retv_vars.10647_0598_003.cdf')
      TES(459)%FILENAME = TRIM('retv_vars.10647_0598_004.cdf')
      TES(460)%FILENAME = TRIM('retv_vars.10647_0599_002.cdf')
      TES(461)%FILENAME = TRIM('retv_vars.10647_0611_003.cdf')
      TES(462)%FILENAME = TRIM('retv_vars.10647_0613_002.cdf')
      TES(463)%FILENAME = TRIM('retv_vars.10647_0615_004.cdf')
      TES(464)%FILENAME = TRIM('retv_vars.10647_0616_003.cdf')
      TES(465)%FILENAME = TRIM('retv_vars.10647_0616_004.cdf')
      TES(466)%FILENAME = TRIM('retv_vars.10647_0617_004.cdf')
      TES(467)%FILENAME = TRIM('retv_vars.10647_0618_002.cdf')
      TES(468)%FILENAME = TRIM('retv_vars.10647_0642_002.cdf')
      TES(469)%FILENAME = TRIM('retv_vars.10647_0642_004.cdf')
      TES(470)%FILENAME = TRIM('retv_vars.10647_0645_003.cdf')
      TES(471)%FILENAME = TRIM('retv_vars.10647_0645_004.cdf')
      TES(472)%FILENAME = TRIM('retv_vars.10647_0646_002.cdf')
      TES(473)%FILENAME = TRIM('retv_vars.10647_0646_003.cdf')
      TES(474)%FILENAME = TRIM('retv_vars.10647_0646_004.cdf')
      TES(475)%FILENAME = TRIM('retv_vars.10647_0647_002.cdf')
      TES(476)%FILENAME = TRIM('retv_vars.10647_0651_002.cdf')
      TES(477)%FILENAME = TRIM('retv_vars.10647_0651_004.cdf')
      TES(478)%FILENAME = TRIM('retv_vars.10647_0653_002.cdf')
      TES(479)%FILENAME = TRIM('retv_vars.10647_0653_003.cdf')
      TES(480)%FILENAME = TRIM('retv_vars.10647_0654_002.cdf')
      TES(481)%FILENAME = TRIM('retv_vars.10647_0654_003.cdf')
      TES(482)%FILENAME = TRIM('retv_vars.10647_0654_004.cdf')
      TES(483)%FILENAME = TRIM('retv_vars.10647_0655_004.cdf')
      TES(484)%FILENAME = TRIM('retv_vars.10647_0656_003.cdf')
      TES(485)%FILENAME = TRIM('retv_vars.10647_0659_004.cdf')
      TES(486)%FILENAME = TRIM('retv_vars.10647_0660_002.cdf')
      TES(487)%FILENAME = TRIM('retv_vars.10647_0689_004.cdf')
      TES(488)%FILENAME = TRIM('retv_vars.10647_0693_002.cdf')
      TES(489)%FILENAME = TRIM('retv_vars.10647_0694_003.cdf')
      TES(490)%FILENAME = TRIM('retv_vars.10647_0694_004.cdf')
      TES(491)%FILENAME = TRIM('retv_vars.10647_0695_002.cdf')
      TES(492)%FILENAME = TRIM('retv_vars.10647_0700_003.cdf')
      TES(493)%FILENAME = TRIM('retv_vars.10647_0700_004.cdf')
      TES(494)%FILENAME = TRIM('retv_vars.10647_0701_003.cdf')
      TES(495)%FILENAME = TRIM('retv_vars.10647_0701_004.cdf')
      TES(496)%FILENAME = TRIM('retv_vars.10647_0702_003.cdf')
      TES(497)%FILENAME = TRIM('retv_vars.10647_0702_004.cdf')
      TES(498)%FILENAME = TRIM('retv_vars.10647_0703_002.cdf')
      TES(499)%FILENAME = TRIM('retv_vars.10647_0704_002.cdf')
      TES(500)%FILENAME = TRIM('retv_vars.10647_0704_003.cdf')
      TES(501)%FILENAME = TRIM('retv_vars.10647_0705_003.cdf')
      TES(502)%FILENAME = TRIM('retv_vars.10647_0733_004.cdf')
      TES(503)%FILENAME = TRIM('retv_vars.10647_0734_003.cdf')
      TES(504)%FILENAME = TRIM('retv_vars.10647_0738_004.cdf')
      TES(505)%FILENAME = TRIM('retv_vars.10647_0739_003.cdf')
      TES(506)%FILENAME = TRIM('retv_vars.10647_0740_004.cdf')
      TES(507)%FILENAME = TRIM('retv_vars.10647_0741_002.cdf')
      TES(508)%FILENAME = TRIM('retv_vars.10647_0741_003.cdf')
      TES(509)%FILENAME = TRIM('retv_vars.10647_0741_004.cdf')
      TES(510)%FILENAME = TRIM('retv_vars.10647_0742_002.cdf')
      TES(511)%FILENAME = TRIM('retv_vars.10647_0742_003.cdf')
      TES(512)%FILENAME = TRIM('retv_vars.10647_0742_004.cdf')
      TES(513)%FILENAME = TRIM('retv_vars.10647_0743_002.cdf')
      TES(514)%FILENAME = TRIM('retv_vars.10647_0743_003.cdf')
      TES(515)%FILENAME = TRIM('retv_vars.10647_0747_004.cdf')
      TES(516)%FILENAME = TRIM('retv_vars.10647_0748_002.cdf')
      TES(517)%FILENAME = TRIM('retv_vars.10647_0748_003.cdf')
      TES(518)%FILENAME = TRIM('retv_vars.10649_0021_003.cdf')
      TES(519)%FILENAME = TRIM('retv_vars.10649_0022_003.cdf')
      TES(520)%FILENAME = TRIM('retv_vars.10649_0027_002.cdf')
      TES(521)%FILENAME = TRIM('retv_vars.10649_0027_004.cdf')
      TES(522)%FILENAME = TRIM('retv_vars.10649_0028_002.cdf')
      TES(523)%FILENAME = TRIM('retv_vars.10649_0028_003.cdf')
      TES(524)%FILENAME = TRIM('retv_vars.10649_0028_004.cdf')
      TES(525)%FILENAME = TRIM('retv_vars.10649_0029_002.cdf')
      TES(526)%FILENAME = TRIM('retv_vars.10649_0029_003.cdf')
      TES(527)%FILENAME = TRIM('retv_vars.10649_0030_002.cdf')
      TES(528)%FILENAME = TRIM('retv_vars.10649_0030_004.cdf')
      TES(529)%FILENAME = TRIM('retv_vars.10649_0043_003.cdf')
      TES(530)%FILENAME = TRIM('retv_vars.10649_0058_002.cdf')
      TES(531)%FILENAME = TRIM('retv_vars.10649_0059_003.cdf')
      TES(532)%FILENAME = TRIM('retv_vars.10649_0059_004.cdf')
      TES(533)%FILENAME = TRIM('retv_vars.10649_0060_003.cdf')
      TES(534)%FILENAME = TRIM('retv_vars.10649_0060_004.cdf')
      TES(535)%FILENAME = TRIM('retv_vars.10649_0065_002.cdf')
      TES(536)%FILENAME = TRIM('retv_vars.10649_0066_004.cdf')
      TES(537)%FILENAME = TRIM('retv_vars.10649_0067_002.cdf')
      TES(538)%FILENAME = TRIM('retv_vars.10649_0067_003.cdf')
      TES(539)%FILENAME = TRIM('retv_vars.10649_0067_004.cdf')
      TES(540)%FILENAME = TRIM('retv_vars.10649_0068_003.cdf')
      TES(541)%FILENAME = TRIM('retv_vars.10649_0069_003.cdf')
      TES(542)%FILENAME = TRIM('retv_vars.10649_0069_004.cdf')
      TES(543)%FILENAME = TRIM('retv_vars.10649_0070_002.cdf')
      TES(544)%FILENAME = TRIM('retv_vars.10649_0075_002.cdf')
      TES(545)%FILENAME = TRIM('retv_vars.10649_0100_002.cdf')
      TES(546)%FILENAME = TRIM('retv_vars.10649_0102_003.cdf')
      TES(547)%FILENAME = TRIM('retv_vars.10649_0102_004.cdf')
      TES(548)%FILENAME = TRIM('retv_vars.10649_0103_003.cdf')
      TES(549)%FILENAME = TRIM('retv_vars.10649_0104_004.cdf')
      TES(550)%FILENAME = TRIM('retv_vars.10649_0105_002.cdf')
      TES(551)%FILENAME = TRIM('retv_vars.10649_0106_004.cdf')
      TES(552)%FILENAME = TRIM('retv_vars.10649_0107_002.cdf')
      TES(553)%FILENAME = TRIM('retv_vars.10649_0107_003.cdf')
      TES(554)%FILENAME = TRIM('retv_vars.10649_0107_004.cdf')
      TES(555)%FILENAME = TRIM('retv_vars.10649_0108_003.cdf')
      TES(556)%FILENAME = TRIM('retv_vars.10649_0108_004.cdf')
      TES(557)%FILENAME = TRIM('retv_vars.10649_0114_002.cdf')
      TES(558)%FILENAME = TRIM('retv_vars.10649_0114_003.cdf')
      TES(559)%FILENAME = TRIM('retv_vars.10649_0115_002.cdf')
      TES(560)%FILENAME = TRIM('retv_vars.10649_0116_002.cdf')
      TES(561)%FILENAME = TRIM('retv_vars.10649_0117_002.cdf')
      TES(562)%FILENAME = TRIM('retv_vars.10649_0117_003.cdf')
      TES(563)%FILENAME = TRIM('retv_vars.10649_0117_004.cdf')
      TES(564)%FILENAME = TRIM('retv_vars.10649_0118_003.cdf')
      TES(565)%FILENAME = TRIM('retv_vars.10649_0118_004.cdf')
      TES(566)%FILENAME = TRIM('retv_vars.10649_0137_004.cdf')
      TES(567)%FILENAME = TRIM('retv_vars.10649_0144_003.cdf')
      TES(568)%FILENAME = TRIM('retv_vars.10649_0144_004.cdf')
      TES(569)%FILENAME = TRIM('retv_vars.10649_0187_004.cdf')
      TES(570)%FILENAME = TRIM('retv_vars.10649_0188_002.cdf')
      TES(571)%FILENAME = TRIM('retv_vars.10649_0189_002.cdf')
      TES(572)%FILENAME = TRIM('retv_vars.10649_0189_003.cdf')
      TES(573)%FILENAME = TRIM('retv_vars.10649_0189_004.cdf')
      TES(574)%FILENAME = TRIM('retv_vars.10649_0191_004.cdf')
      TES(575)%FILENAME = TRIM('retv_vars.10649_0201_002.cdf')
      TES(576)%FILENAME = TRIM('retv_vars.10649_0221_003.cdf')
      TES(577)%FILENAME = TRIM('retv_vars.10649_0222_003.cdf')
      TES(578)%FILENAME = TRIM('retv_vars.10649_0223_003.cdf')
      TES(579)%FILENAME = TRIM('retv_vars.10649_0223_004.cdf')
      TES(580)%FILENAME = TRIM('retv_vars.10649_0231_002.cdf')
      TES(581)%FILENAME = TRIM('retv_vars.10649_0231_003.cdf')
      TES(582)%FILENAME = TRIM('retv_vars.10649_0234_002.cdf')
      TES(583)%FILENAME = TRIM('retv_vars.10649_0236_002.cdf')
      TES(584)%FILENAME = TRIM('retv_vars.10649_0237_004.cdf')
      TES(585)%FILENAME = TRIM('retv_vars.10649_0244_002.cdf')
      TES(586)%FILENAME = TRIM('retv_vars.10649_0244_003.cdf')
      TES(587)%FILENAME = TRIM('retv_vars.10649_0244_004.cdf')
      TES(588)%FILENAME = TRIM('retv_vars.10649_0245_004.cdf')
      TES(589)%FILENAME = TRIM('retv_vars.10649_0246_003.cdf')
      TES(590)%FILENAME = TRIM('retv_vars.10649_0246_004.cdf')
      TES(591)%FILENAME = TRIM('retv_vars.10649_0247_002.cdf')
      TES(592)%FILENAME = TRIM('retv_vars.10649_0247_003.cdf')
      TES(593)%FILENAME = TRIM('retv_vars.10649_0248_002.cdf')
      TES(594)%FILENAME = TRIM('retv_vars.10649_0248_003.cdf')
      TES(595)%FILENAME = TRIM('retv_vars.10649_0248_004.cdf')
      TES(596)%FILENAME = TRIM('retv_vars.10649_0249_003.cdf')
      TES(597)%FILENAME = TRIM('retv_vars.10649_0249_004.cdf')
      TES(598)%FILENAME = TRIM('retv_vars.10649_0250_002.cdf')
      TES(599)%FILENAME = TRIM('retv_vars.10649_0252_002.cdf')
      TES(600)%FILENAME = TRIM('retv_vars.10649_0252_003.cdf')
      TES(601)%FILENAME = TRIM('retv_vars.10649_0252_004.cdf')
      TES(602)%FILENAME = TRIM('retv_vars.10649_0260_002.cdf')
      TES(603)%FILENAME = TRIM('retv_vars.10649_0260_003.cdf')
      TES(604)%FILENAME = TRIM('retv_vars.10649_0261_002.cdf')
      TES(605)%FILENAME = TRIM('retv_vars.10649_0261_004.cdf')
      TES(606)%FILENAME = TRIM('retv_vars.10649_0267_002.cdf')
      TES(607)%FILENAME = TRIM('retv_vars.10649_0269_003.cdf')
      TES(608)%FILENAME = TRIM('retv_vars.10649_0270_004.cdf')
      TES(609)%FILENAME = TRIM('retv_vars.10649_0271_004.cdf')
      TES(610)%FILENAME = TRIM('retv_vars.10649_0273_002.cdf')
      TES(611)%FILENAME = TRIM('retv_vars.10649_0273_003.cdf')
      TES(612)%FILENAME = TRIM('retv_vars.10649_0274_002.cdf')
      TES(613)%FILENAME = TRIM('retv_vars.10649_0274_004.cdf')
      TES(614)%FILENAME = TRIM('retv_vars.10649_0276_003.cdf')
      TES(615)%FILENAME = TRIM('retv_vars.10649_0280_002.cdf')
      TES(616)%FILENAME = TRIM('retv_vars.10649_0302_002.cdf')
      TES(617)%FILENAME = TRIM('retv_vars.10649_0302_003.cdf')
      TES(618)%FILENAME = TRIM('retv_vars.10649_0305_003.cdf')
      TES(619)%FILENAME = TRIM('retv_vars.10649_0307_002.cdf')
      TES(620)%FILENAME = TRIM('retv_vars.10649_0307_003.cdf')
      TES(621)%FILENAME = TRIM('retv_vars.10649_0307_004.cdf')
      TES(622)%FILENAME = TRIM('retv_vars.10649_0308_002.cdf')
      TES(623)%FILENAME = TRIM('retv_vars.10649_0308_003.cdf')
      TES(624)%FILENAME = TRIM('retv_vars.10649_0309_003.cdf')
      TES(625)%FILENAME = TRIM('retv_vars.10649_0310_003.cdf')
      TES(626)%FILENAME = TRIM('retv_vars.10649_0315_002.cdf')
      TES(627)%FILENAME = TRIM('retv_vars.10649_0315_003.cdf')
      TES(628)%FILENAME = TRIM('retv_vars.10649_0315_004.cdf')
      TES(629)%FILENAME = TRIM('retv_vars.10649_0316_003.cdf')
      TES(630)%FILENAME = TRIM('retv_vars.10649_0316_004.cdf')
      TES(631)%FILENAME = TRIM('retv_vars.10649_0318_003.cdf')
      TES(632)%FILENAME = TRIM('retv_vars.10649_0319_003.cdf')
      TES(633)%FILENAME = TRIM('retv_vars.10649_0319_004.cdf')
      TES(634)%FILENAME = TRIM('retv_vars.10649_0320_003.cdf')
      TES(635)%FILENAME = TRIM('retv_vars.10649_0321_004.cdf')
      TES(636)%FILENAME = TRIM('retv_vars.10649_0322_002.cdf')
      TES(637)%FILENAME = TRIM('retv_vars.10649_0322_003.cdf')
      TES(638)%FILENAME = TRIM('retv_vars.10649_0322_004.cdf')
      TES(639)%FILENAME = TRIM('retv_vars.10649_0355_002.cdf')
      TES(640)%FILENAME = TRIM('retv_vars.10649_0355_004.cdf')
      TES(641)%FILENAME = TRIM('retv_vars.10649_0357_002.cdf')
      TES(642)%FILENAME = TRIM('retv_vars.10649_0357_003.cdf')
      TES(643)%FILENAME = TRIM('retv_vars.10649_0358_002.cdf')
      TES(644)%FILENAME = TRIM('retv_vars.10649_0358_003.cdf')
      TES(645)%FILENAME = TRIM('retv_vars.10649_0363_004.cdf')
      TES(646)%FILENAME = TRIM('retv_vars.10649_0364_003.cdf')
      TES(647)%FILENAME = TRIM('retv_vars.10649_0364_004.cdf')
      TES(648)%FILENAME = TRIM('retv_vars.10649_0365_002.cdf')
      TES(649)%FILENAME = TRIM('retv_vars.10649_0365_004.cdf')
      TES(650)%FILENAME = TRIM('retv_vars.10649_0368_004.cdf')
      TES(651)%FILENAME = TRIM('retv_vars.10649_0369_004.cdf')
      TES(652)%FILENAME = TRIM('retv_vars.10649_0407_002.cdf')
      TES(653)%FILENAME = TRIM('retv_vars.10649_0411_003.cdf')
      TES(654)%FILENAME = TRIM('retv_vars.10649_0413_003.cdf')
      TES(655)%FILENAME = TRIM('retv_vars.10649_0413_004.cdf')
      TES(656)%FILENAME = TRIM('retv_vars.10649_0414_002.cdf')
      TES(657)%FILENAME = TRIM('retv_vars.10649_0414_003.cdf')
      TES(658)%FILENAME = TRIM('retv_vars.10649_0414_004.cdf')
      TES(659)%FILENAME = TRIM('retv_vars.10649_0415_002.cdf')
      TES(660)%FILENAME = TRIM('retv_vars.10649_0415_003.cdf')
      TES(661)%FILENAME = TRIM('retv_vars.10649_0421_002.cdf')
      TES(662)%FILENAME = TRIM('retv_vars.10649_0422_002.cdf')
      TES(663)%FILENAME = TRIM('retv_vars.10649_0422_003.cdf')
      TES(664)%FILENAME = TRIM('retv_vars.10649_0422_004.cdf')
      TES(665)%FILENAME = TRIM('retv_vars.10649_0423_003.cdf')
      TES(666)%FILENAME = TRIM('retv_vars.10649_0423_004.cdf')
      TES(667)%FILENAME = TRIM('retv_vars.10649_0424_002.cdf')
      TES(668)%FILENAME = TRIM('retv_vars.10649_0424_003.cdf')
      TES(669)%FILENAME = TRIM('retv_vars.10649_0424_004.cdf')
      TES(670)%FILENAME = TRIM('retv_vars.10649_0425_002.cdf')
      TES(671)%FILENAME = TRIM('retv_vars.10649_0425_003.cdf')
      TES(672)%FILENAME = TRIM('retv_vars.10649_0425_004.cdf')
      TES(673)%FILENAME = TRIM('retv_vars.10649_0426_002.cdf')
      TES(674)%FILENAME = TRIM('retv_vars.10649_0426_004.cdf')
      TES(675)%FILENAME = TRIM('retv_vars.10649_0428_002.cdf')
      TES(676)%FILENAME = TRIM('retv_vars.10649_0429_003.cdf')
      TES(677)%FILENAME = TRIM('retv_vars.10649_0459_002.cdf')
      TES(678)%FILENAME = TRIM('retv_vars.10649_0459_004.cdf')
      TES(679)%FILENAME = TRIM('retv_vars.10649_0460_002.cdf')
      TES(680)%FILENAME = TRIM('retv_vars.10649_0460_003.cdf')
      TES(681)%FILENAME = TRIM('retv_vars.10649_0460_004.cdf')
      TES(682)%FILENAME = TRIM('retv_vars.10649_0461_003.cdf')
      TES(683)%FILENAME = TRIM('retv_vars.10649_0461_004.cdf')
      TES(684)%FILENAME = TRIM('retv_vars.10649_0463_002.cdf')
      TES(685)%FILENAME = TRIM('retv_vars.10649_0463_003.cdf')
      TES(686)%FILENAME = TRIM('retv_vars.10649_0466_003.cdf')
      TES(687)%FILENAME = TRIM('retv_vars.10649_0467_002.cdf')
      TES(688)%FILENAME = TRIM('retv_vars.10649_0467_003.cdf')
      TES(689)%FILENAME = TRIM('retv_vars.10649_0469_002.cdf')
      TES(690)%FILENAME = TRIM('retv_vars.10649_0469_003.cdf')
      TES(691)%FILENAME = TRIM('retv_vars.10649_0469_004.cdf')
      TES(692)%FILENAME = TRIM('retv_vars.10649_0470_002.cdf')
      TES(693)%FILENAME = TRIM('retv_vars.10649_0502_002.cdf')
      TES(694)%FILENAME = TRIM('retv_vars.10649_0502_003.cdf')
      TES(695)%FILENAME = TRIM('retv_vars.10649_0510_003.cdf')
      TES(696)%FILENAME = TRIM('retv_vars.10649_0510_004.cdf')
      TES(697)%FILENAME = TRIM('retv_vars.10649_0511_002.cdf')
      TES(698)%FILENAME = TRIM('retv_vars.10649_0511_003.cdf')
      TES(699)%FILENAME = TRIM('retv_vars.10649_0514_003.cdf')
      TES(700)%FILENAME = TRIM('retv_vars.10649_0515_002.cdf')
      TES(701)%FILENAME = TRIM('retv_vars.10649_0516_004.cdf')
      TES(702)%FILENAME = TRIM('retv_vars.10649_0538_002.cdf')
      TES(703)%FILENAME = TRIM('retv_vars.10649_0546_004.cdf')
      TES(704)%FILENAME = TRIM('retv_vars.10649_0547_002.cdf')
      TES(705)%FILENAME = TRIM('retv_vars.10649_0548_004.cdf')
      TES(706)%FILENAME = TRIM('retv_vars.10649_0549_002.cdf')
      TES(707)%FILENAME = TRIM('retv_vars.10649_0549_003.cdf')
      TES(708)%FILENAME = TRIM('retv_vars.10649_0549_004.cdf')
      TES(709)%FILENAME = TRIM('retv_vars.10649_0550_002.cdf')
      TES(710)%FILENAME = TRIM('retv_vars.10649_0550_003.cdf')
      TES(711)%FILENAME = TRIM('retv_vars.10649_0551_002.cdf')
      TES(712)%FILENAME = TRIM('retv_vars.10649_0568_002.cdf')
      TES(713)%FILENAME = TRIM('retv_vars.10649_0568_004.cdf')
      TES(714)%FILENAME = TRIM('retv_vars.10649_0569_002.cdf')
      TES(715)%FILENAME = TRIM('retv_vars.10649_0569_003.cdf')
      TES(716)%FILENAME = TRIM('retv_vars.10649_0570_003.cdf')
      TES(717)%FILENAME = TRIM('retv_vars.10649_0571_002.cdf')
      TES(718)%FILENAME = TRIM('retv_vars.10649_0571_003.cdf')
      TES(719)%FILENAME = TRIM('retv_vars.10649_0572_002.cdf')
      TES(720)%FILENAME = TRIM('retv_vars.10649_0573_002.cdf')
      TES(721)%FILENAME = TRIM('retv_vars.10649_0582_002.cdf')
      TES(722)%FILENAME = TRIM('retv_vars.10649_0583_003.cdf')
      TES(723)%FILENAME = TRIM('retv_vars.10649_0583_004.cdf')
      TES(724)%FILENAME = TRIM('retv_vars.10649_0588_003.cdf')
      TES(725)%FILENAME = TRIM('retv_vars.10649_0588_004.cdf')
      TES(726)%FILENAME = TRIM('retv_vars.10649_0589_002.cdf')
      TES(727)%FILENAME = TRIM('retv_vars.10649_0589_003.cdf')
      TES(728)%FILENAME = TRIM('retv_vars.10649_0590_002.cdf')
      TES(729)%FILENAME = TRIM('retv_vars.10649_0590_003.cdf')
      TES(730)%FILENAME = TRIM('retv_vars.10649_0592_002.cdf')
      TES(731)%FILENAME = TRIM('retv_vars.10649_0592_004.cdf')
      TES(732)%FILENAME = TRIM('retv_vars.10649_0593_002.cdf')
      TES(733)%FILENAME = TRIM('retv_vars.10649_0593_004.cdf')
      TES(734)%FILENAME = TRIM('retv_vars.10649_0594_002.cdf')
      TES(735)%FILENAME = TRIM('retv_vars.10649_0594_003.cdf')
      TES(736)%FILENAME = TRIM('retv_vars.10649_0594_004.cdf')
      TES(737)%FILENAME = TRIM('retv_vars.10649_0596_002.cdf')
      TES(738)%FILENAME = TRIM('retv_vars.10649_0596_003.cdf')
      TES(739)%FILENAME = TRIM('retv_vars.10649_0596_004.cdf')
      TES(740)%FILENAME = TRIM('retv_vars.10649_0597_002.cdf')
      TES(741)%FILENAME = TRIM('retv_vars.10649_0597_003.cdf')
      TES(742)%FILENAME = TRIM('retv_vars.10649_0597_004.cdf')
      TES(743)%FILENAME = TRIM('retv_vars.10649_0598_003.cdf')
      TES(744)%FILENAME = TRIM('retv_vars.10649_0599_002.cdf')
      TES(745)%FILENAME = TRIM('retv_vars.10649_0605_003.cdf')
      TES(746)%FILENAME = TRIM('retv_vars.10649_0613_002.cdf')
      TES(747)%FILENAME = TRIM('retv_vars.10649_0615_004.cdf')
      TES(748)%FILENAME = TRIM('retv_vars.10649_0616_003.cdf')
      TES(749)%FILENAME = TRIM('retv_vars.10649_0616_004.cdf')
      TES(750)%FILENAME = TRIM('retv_vars.10649_0617_002.cdf')
      TES(751)%FILENAME = TRIM('retv_vars.10649_0617_003.cdf')
      TES(752)%FILENAME = TRIM('retv_vars.10649_0618_002.cdf')
      TES(753)%FILENAME = TRIM('retv_vars.10649_0618_003.cdf')
      TES(754)%FILENAME = TRIM('retv_vars.10649_0634_004.cdf')
      TES(755)%FILENAME = TRIM('retv_vars.10649_0635_003.cdf')
      TES(756)%FILENAME = TRIM('retv_vars.10649_0639_004.cdf')
      TES(757)%FILENAME = TRIM('retv_vars.10649_0644_004.cdf')
      TES(758)%FILENAME = TRIM('retv_vars.10649_0645_002.cdf')
      TES(759)%FILENAME = TRIM('retv_vars.10649_0645_003.cdf')
      TES(760)%FILENAME = TRIM('retv_vars.10649_0645_004.cdf')
      TES(761)%FILENAME = TRIM('retv_vars.10649_0646_002.cdf')
      TES(762)%FILENAME = TRIM('retv_vars.10649_0646_003.cdf')
      TES(763)%FILENAME = TRIM('retv_vars.10649_0646_004.cdf')
      TES(764)%FILENAME = TRIM('retv_vars.10649_0647_002.cdf')
      TES(765)%FILENAME = TRIM('retv_vars.10649_0651_002.cdf')
      TES(766)%FILENAME = TRIM('retv_vars.10649_0651_003.cdf')
      TES(767)%FILENAME = TRIM('retv_vars.10649_0651_004.cdf')
      TES(768)%FILENAME = TRIM('retv_vars.10649_0652_002.cdf')
      TES(769)%FILENAME = TRIM('retv_vars.10649_0652_003.cdf')
      TES(770)%FILENAME = TRIM('retv_vars.10649_0652_004.cdf')
      TES(771)%FILENAME = TRIM('retv_vars.10649_0653_002.cdf')
      TES(772)%FILENAME = TRIM('retv_vars.10649_0653_004.cdf')
      TES(773)%FILENAME = TRIM('retv_vars.10649_0654_002.cdf')
      TES(774)%FILENAME = TRIM('retv_vars.10649_0654_003.cdf')
      TES(775)%FILENAME = TRIM('retv_vars.10649_0655_002.cdf')
      TES(776)%FILENAME = TRIM('retv_vars.10649_0655_003.cdf')
      TES(777)%FILENAME = TRIM('retv_vars.10649_0656_003.cdf')
      TES(778)%FILENAME = TRIM('retv_vars.10649_0659_002.cdf')
      TES(779)%FILENAME = TRIM('retv_vars.10649_0659_004.cdf')
      TES(780)%FILENAME = TRIM('retv_vars.10649_0690_002.cdf')
      TES(781)%FILENAME = TRIM('retv_vars.10649_0693_004.cdf')
      TES(782)%FILENAME = TRIM('retv_vars.10649_0694_003.cdf')
      TES(783)%FILENAME = TRIM('retv_vars.10649_0699_002.cdf')
      TES(784)%FILENAME = TRIM('retv_vars.10649_0699_004.cdf')
      TES(785)%FILENAME = TRIM('retv_vars.10649_0700_002.cdf')
      TES(786)%FILENAME = TRIM('retv_vars.10649_0700_003.cdf')
      TES(787)%FILENAME = TRIM('retv_vars.10649_0700_004.cdf')
      TES(788)%FILENAME = TRIM('retv_vars.10649_0701_002.cdf')
      TES(789)%FILENAME = TRIM('retv_vars.10649_0701_004.cdf')
      TES(790)%FILENAME = TRIM('retv_vars.10649_0702_004.cdf')
      TES(791)%FILENAME = TRIM('retv_vars.10649_0703_002.cdf')
      TES(792)%FILENAME = TRIM('retv_vars.10649_0703_004.cdf')
      TES(793)%FILENAME = TRIM('retv_vars.10649_0704_003.cdf')
      TES(794)%FILENAME = TRIM('retv_vars.10649_0704_004.cdf')
      TES(795)%FILENAME = TRIM('retv_vars.10649_0732_004.cdf')
      TES(796)%FILENAME = TRIM('retv_vars.10649_0734_004.cdf')
      TES(797)%FILENAME = TRIM('retv_vars.10649_0739_002.cdf')
      TES(798)%FILENAME = TRIM('retv_vars.10649_0740_002.cdf')
      TES(799)%FILENAME = TRIM('retv_vars.10649_0741_002.cdf')
      TES(800)%FILENAME = TRIM('retv_vars.10649_0741_003.cdf')
      TES(801)%FILENAME = TRIM('retv_vars.10649_0741_004.cdf')
      TES(802)%FILENAME = TRIM('retv_vars.10649_0743_003.cdf')
      TES(803)%FILENAME = TRIM('retv_vars.10649_0748_002.cdf')
      TES(804)%FILENAME = TRIM('retv_vars.10656_0018_003.cdf')
      TES(805)%FILENAME = TRIM('retv_vars.10656_0021_003.cdf')
      TES(806)%FILENAME = TRIM('retv_vars.10656_0021_004.cdf')
      TES(807)%FILENAME = TRIM('retv_vars.10656_0022_002.cdf')
      TES(808)%FILENAME = TRIM('retv_vars.10656_0022_003.cdf')
      TES(809)%FILENAME = TRIM('retv_vars.10656_0022_004.cdf')
      TES(810)%FILENAME = TRIM('retv_vars.10656_0023_002.cdf')
      TES(811)%FILENAME = TRIM('retv_vars.10656_0027_003.cdf')
      TES(812)%FILENAME = TRIM('retv_vars.10656_0028_003.cdf')
      TES(813)%FILENAME = TRIM('retv_vars.10656_0029_002.cdf')
      TES(814)%FILENAME = TRIM('retv_vars.10656_0029_003.cdf')
      TES(815)%FILENAME = TRIM('retv_vars.10656_0029_004.cdf')
      TES(816)%FILENAME = TRIM('retv_vars.10656_0030_002.cdf')
      TES(817)%FILENAME = TRIM('retv_vars.10656_0030_003.cdf')
      TES(818)%FILENAME = TRIM('retv_vars.10656_0030_004.cdf')
      TES(819)%FILENAME = TRIM('retv_vars.10656_0031_002.cdf')
      TES(820)%FILENAME = TRIM('retv_vars.10656_0054_003.cdf')
      TES(821)%FILENAME = TRIM('retv_vars.10656_0055_002.cdf')
      TES(822)%FILENAME = TRIM('retv_vars.10656_0055_003.cdf')
      TES(823)%FILENAME = TRIM('retv_vars.10656_0055_004.cdf')
      TES(824)%FILENAME = TRIM('retv_vars.10656_0059_004.cdf')
      TES(825)%FILENAME = TRIM('retv_vars.10656_0060_002.cdf')
      TES(826)%FILENAME = TRIM('retv_vars.10656_0060_003.cdf')
      TES(827)%FILENAME = TRIM('retv_vars.10656_0066_004.cdf')
      TES(828)%FILENAME = TRIM('retv_vars.10656_0067_002.cdf')
      TES(829)%FILENAME = TRIM('retv_vars.10656_0069_003.cdf')
      TES(830)%FILENAME = TRIM('retv_vars.10656_0069_004.cdf')
      TES(831)%FILENAME = TRIM('retv_vars.10656_0070_002.cdf')
      TES(832)%FILENAME = TRIM('retv_vars.10656_0070_003.cdf')
      TES(833)%FILENAME = TRIM('retv_vars.10656_0070_004.cdf')
      TES(834)%FILENAME = TRIM('retv_vars.10656_0075_002.cdf')
      TES(835)%FILENAME = TRIM('retv_vars.10656_0075_003.cdf')
      TES(836)%FILENAME = TRIM('retv_vars.10656_0075_004.cdf')
      TES(837)%FILENAME = TRIM('retv_vars.10656_0100_004.cdf')
      TES(838)%FILENAME = TRIM('retv_vars.10656_0101_002.cdf')
      TES(839)%FILENAME = TRIM('retv_vars.10656_0101_003.cdf')
      TES(840)%FILENAME = TRIM('retv_vars.10656_0102_004.cdf')
      TES(841)%FILENAME = TRIM('retv_vars.10656_0103_004.cdf')
      TES(842)%FILENAME = TRIM('retv_vars.10656_0104_002.cdf')
      TES(843)%FILENAME = TRIM('retv_vars.10656_0104_003.cdf')
      TES(844)%FILENAME = TRIM('retv_vars.10656_0105_002.cdf')
      TES(845)%FILENAME = TRIM('retv_vars.10656_0105_004.cdf')
      TES(846)%FILENAME = TRIM('retv_vars.10656_0106_002.cdf')
      TES(847)%FILENAME = TRIM('retv_vars.10656_0106_003.cdf')
      TES(848)%FILENAME = TRIM('retv_vars.10656_0106_004.cdf')
      TES(849)%FILENAME = TRIM('retv_vars.10656_0110_003.cdf')
      TES(850)%FILENAME = TRIM('retv_vars.10656_0110_004.cdf')
      TES(851)%FILENAME = TRIM('retv_vars.10656_0114_003.cdf')
      TES(852)%FILENAME = TRIM('retv_vars.10656_0116_003.cdf')
      TES(853)%FILENAME = TRIM('retv_vars.10656_0123_003.cdf')
      TES(854)%FILENAME = TRIM('retv_vars.10656_0143_002.cdf')
      TES(855)%FILENAME = TRIM('retv_vars.10656_0143_003.cdf')
      TES(856)%FILENAME = TRIM('retv_vars.10656_0157_004.cdf')
      TES(857)%FILENAME = TRIM('retv_vars.10656_0158_002.cdf')
      TES(858)%FILENAME = TRIM('retv_vars.10656_0233_003.cdf')
      TES(859)%FILENAME = TRIM('retv_vars.10656_0247_003.cdf')
      TES(860)%FILENAME = TRIM('retv_vars.10656_0247_004.cdf')
      TES(861)%FILENAME = TRIM('retv_vars.10656_0248_003.cdf')
      TES(862)%FILENAME = TRIM('retv_vars.10656_0248_004.cdf')
      TES(863)%FILENAME = TRIM('retv_vars.10656_0250_003.cdf')
      TES(864)%FILENAME = TRIM('retv_vars.10656_0251_003.cdf')
      TES(865)%FILENAME = TRIM('retv_vars.10656_0251_004.cdf')
      TES(866)%FILENAME = TRIM('retv_vars.10656_0252_002.cdf')
      TES(867)%FILENAME = TRIM('retv_vars.10656_0260_002.cdf')
      TES(868)%FILENAME = TRIM('retv_vars.10656_0260_003.cdf')
      TES(869)%FILENAME = TRIM('retv_vars.10656_0260_004.cdf')
      TES(870)%FILENAME = TRIM('retv_vars.10656_0261_002.cdf')
      TES(871)%FILENAME = TRIM('retv_vars.10656_0261_003.cdf')
      TES(872)%FILENAME = TRIM('retv_vars.10656_0261_004.cdf')
      TES(873)%FILENAME = TRIM('retv_vars.10656_0262_002.cdf')
      TES(874)%FILENAME = TRIM('retv_vars.10656_0262_003.cdf')
      TES(875)%FILENAME = TRIM('retv_vars.10656_0263_002.cdf')
      TES(876)%FILENAME = TRIM('retv_vars.10656_0267_002.cdf')
      TES(877)%FILENAME = TRIM('retv_vars.10656_0267_003.cdf')
      TES(878)%FILENAME = TRIM('retv_vars.10656_0268_004.cdf')
      TES(879)%FILENAME = TRIM('retv_vars.10656_0269_002.cdf')
      TES(880)%FILENAME = TRIM('retv_vars.10656_0269_003.cdf')
      TES(881)%FILENAME = TRIM('retv_vars.10656_0269_004.cdf')
      TES(882)%FILENAME = TRIM('retv_vars.10656_0271_004.cdf')
      TES(883)%FILENAME = TRIM('retv_vars.10656_0273_003.cdf')
      TES(884)%FILENAME = TRIM('retv_vars.10656_0273_004.cdf')
      TES(885)%FILENAME = TRIM('retv_vars.10656_0274_003.cdf')
      TES(886)%FILENAME = TRIM('retv_vars.10656_0278_004.cdf')
      TES(887)%FILENAME = TRIM('retv_vars.10656_0279_002.cdf')
      TES(888)%FILENAME = TRIM('retv_vars.10656_0279_003.cdf')
      TES(889)%FILENAME = TRIM('retv_vars.10656_0279_004.cdf')
      TES(890)%FILENAME = TRIM('retv_vars.10656_0289_002.cdf')
      TES(891)%FILENAME = TRIM('retv_vars.10656_0289_003.cdf')
      TES(892)%FILENAME = TRIM('retv_vars.10656_0289_004.cdf')
      TES(893)%FILENAME = TRIM('retv_vars.10656_0290_004.cdf')
      TES(894)%FILENAME = TRIM('retv_vars.10656_0291_002.cdf')
      TES(895)%FILENAME = TRIM('retv_vars.10656_0300_004.cdf')
      TES(896)%FILENAME = TRIM('retv_vars.10656_0302_002.cdf')
      TES(897)%FILENAME = TRIM('retv_vars.10656_0305_003.cdf')
      TES(898)%FILENAME = TRIM('retv_vars.10656_0305_004.cdf')
      TES(899)%FILENAME = TRIM('retv_vars.10656_0306_004.cdf')
      TES(900)%FILENAME = TRIM('retv_vars.10656_0307_002.cdf')
      TES(901)%FILENAME = TRIM('retv_vars.10656_0307_003.cdf')
      TES(902)%FILENAME = TRIM('retv_vars.10656_0307_004.cdf')
      TES(903)%FILENAME = TRIM('retv_vars.10656_0308_003.cdf')
      TES(904)%FILENAME = TRIM('retv_vars.10656_0308_004.cdf')
      TES(905)%FILENAME = TRIM('retv_vars.10656_0315_003.cdf')
      TES(906)%FILENAME = TRIM('retv_vars.10656_0315_004.cdf')
      TES(907)%FILENAME = TRIM('retv_vars.10656_0316_002.cdf')
      TES(908)%FILENAME = TRIM('retv_vars.10656_0316_003.cdf')
      TES(909)%FILENAME = TRIM('retv_vars.10656_0318_003.cdf')
      TES(910)%FILENAME = TRIM('retv_vars.10656_0319_004.cdf')
      TES(911)%FILENAME = TRIM('retv_vars.10656_0320_004.cdf')
      TES(912)%FILENAME = TRIM('retv_vars.10656_0321_002.cdf')
      TES(913)%FILENAME = TRIM('retv_vars.10656_0322_004.cdf')
      TES(914)%FILENAME = TRIM('retv_vars.10656_0355_004.cdf')
      TES(915)%FILENAME = TRIM('retv_vars.10656_0356_003.cdf')
      TES(916)%FILENAME = TRIM('retv_vars.10656_0356_004.cdf')
      TES(917)%FILENAME = TRIM('retv_vars.10656_0357_002.cdf')
      TES(918)%FILENAME = TRIM('retv_vars.10656_0357_003.cdf')
      TES(919)%FILENAME = TRIM('retv_vars.10656_0363_003.cdf')
      TES(920)%FILENAME = TRIM('retv_vars.10656_0363_004.cdf')
      TES(921)%FILENAME = TRIM('retv_vars.10656_0364_003.cdf')
      TES(922)%FILENAME = TRIM('retv_vars.10656_0365_002.cdf')
      TES(923)%FILENAME = TRIM('retv_vars.10656_0365_003.cdf')
      TES(924)%FILENAME = TRIM('retv_vars.10656_0365_004.cdf')
      TES(925)%FILENAME = TRIM('retv_vars.10656_0366_002.cdf')
      TES(926)%FILENAME = TRIM('retv_vars.10656_0366_003.cdf')
      TES(927)%FILENAME = TRIM('retv_vars.10656_0366_004.cdf')
      TES(928)%FILENAME = TRIM('retv_vars.10656_0367_002.cdf')
      TES(929)%FILENAME = TRIM('retv_vars.10656_0367_004.cdf')
      TES(930)%FILENAME = TRIM('retv_vars.10656_0368_003.cdf')
      TES(931)%FILENAME = TRIM('retv_vars.10656_0411_002.cdf')
      TES(932)%FILENAME = TRIM('retv_vars.10656_0411_003.cdf')
      TES(933)%FILENAME = TRIM('retv_vars.10656_0411_004.cdf')
      TES(934)%FILENAME = TRIM('retv_vars.10656_0412_004.cdf')
      TES(935)%FILENAME = TRIM('retv_vars.10656_0413_002.cdf')
      TES(936)%FILENAME = TRIM('retv_vars.10656_0415_002.cdf')
      TES(937)%FILENAME = TRIM('retv_vars.10656_0415_004.cdf')
      TES(938)%FILENAME = TRIM('retv_vars.10656_0416_002.cdf')
      TES(939)%FILENAME = TRIM('retv_vars.10656_0421_002.cdf')
      TES(940)%FILENAME = TRIM('retv_vars.10656_0421_003.cdf')
      TES(941)%FILENAME = TRIM('retv_vars.10656_0421_004.cdf')
      TES(942)%FILENAME = TRIM('retv_vars.10656_0422_002.cdf')
      TES(943)%FILENAME = TRIM('retv_vars.10656_0422_003.cdf')
      TES(944)%FILENAME = TRIM('retv_vars.10656_0423_002.cdf')
      TES(945)%FILENAME = TRIM('retv_vars.10656_0424_002.cdf')
      TES(946)%FILENAME = TRIM('retv_vars.10656_0424_004.cdf')
      TES(947)%FILENAME = TRIM('retv_vars.10656_0425_002.cdf')
      TES(948)%FILENAME = TRIM('retv_vars.10656_0425_004.cdf')
      TES(949)%FILENAME = TRIM('retv_vars.10656_0426_004.cdf')
      TES(950)%FILENAME = TRIM('retv_vars.10656_0427_002.cdf')
      TES(951)%FILENAME = TRIM('retv_vars.10656_0429_002.cdf')
      TES(952)%FILENAME = TRIM('retv_vars.10656_0447_003.cdf')
      TES(953)%FILENAME = TRIM('retv_vars.10656_0459_002.cdf')
      TES(954)%FILENAME = TRIM('retv_vars.10656_0459_003.cdf')
      TES(955)%FILENAME = TRIM('retv_vars.10656_0459_004.cdf')
      TES(956)%FILENAME = TRIM('retv_vars.10656_0460_002.cdf')
      TES(957)%FILENAME = TRIM('retv_vars.10656_0460_003.cdf')
      TES(958)%FILENAME = TRIM('retv_vars.10656_0460_004.cdf')
      TES(959)%FILENAME = TRIM('retv_vars.10656_0461_002.cdf')
      TES(960)%FILENAME = TRIM('retv_vars.10656_0461_003.cdf')
      TES(961)%FILENAME = TRIM('retv_vars.10656_0461_004.cdf')
      TES(962)%FILENAME = TRIM('retv_vars.10656_0462_002.cdf')
      TES(963)%FILENAME = TRIM('retv_vars.10656_0462_003.cdf')
      TES(964)%FILENAME = TRIM('retv_vars.10656_0462_004.cdf')
      TES(965)%FILENAME = TRIM('retv_vars.10656_0463_002.cdf')
      TES(966)%FILENAME = TRIM('retv_vars.10656_0463_003.cdf')
      TES(967)%FILENAME = TRIM('retv_vars.10656_0465_003.cdf')
      TES(968)%FILENAME = TRIM('retv_vars.10656_0465_004.cdf')
      TES(969)%FILENAME = TRIM('retv_vars.10656_0466_002.cdf')
      TES(970)%FILENAME = TRIM('retv_vars.10656_0466_003.cdf')
      TES(971)%FILENAME = TRIM('retv_vars.10656_0467_004.cdf')
      TES(972)%FILENAME = TRIM('retv_vars.10656_0468_002.cdf')
      TES(973)%FILENAME = TRIM('retv_vars.10656_0470_004.cdf')
      TES(974)%FILENAME = TRIM('retv_vars.10656_0471_002.cdf')
      TES(975)%FILENAME = TRIM('retv_vars.10656_0471_003.cdf')
      TES(976)%FILENAME = TRIM('retv_vars.10656_0471_004.cdf')
      TES(977)%FILENAME = TRIM('retv_vars.10656_0472_002.cdf')
      TES(978)%FILENAME = TRIM('retv_vars.10656_0486_004.cdf')
      TES(979)%FILENAME = TRIM('retv_vars.10656_0507_002.cdf')
      TES(980)%FILENAME = TRIM('retv_vars.10656_0509_003.cdf')
      TES(981)%FILENAME = TRIM('retv_vars.10656_0509_004.cdf')
      TES(982)%FILENAME = TRIM('retv_vars.10656_0510_003.cdf')
      TES(983)%FILENAME = TRIM('retv_vars.10656_0511_002.cdf')
      TES(984)%FILENAME = TRIM('retv_vars.10656_0511_003.cdf')
      TES(985)%FILENAME = TRIM('retv_vars.10656_0511_004.cdf')
      TES(986)%FILENAME = TRIM('retv_vars.10656_0515_004.cdf')
      TES(987)%FILENAME = TRIM('retv_vars.10656_0516_003.cdf')
      TES(988)%FILENAME = TRIM('retv_vars.10656_0516_004.cdf')
      TES(989)%FILENAME = TRIM('retv_vars.10656_0517_003.cdf')
      TES(990)%FILENAME = TRIM('retv_vars.10656_0549_002.cdf')
      TES(991)%FILENAME = TRIM('retv_vars.10656_0549_003.cdf')
      TES(992)%FILENAME = TRIM('retv_vars.10656_0550_002.cdf')
      TES(993)%FILENAME = TRIM('retv_vars.10656_0551_002.cdf')
      TES(994)%FILENAME = TRIM('retv_vars.10656_0568_003.cdf')
      TES(995)%FILENAME = TRIM('retv_vars.10656_0569_002.cdf')
      TES(996)%FILENAME = TRIM('retv_vars.10656_0569_003.cdf')
      TES(997)%FILENAME = TRIM('retv_vars.10656_0569_004.cdf')
      TES(998)%FILENAME = TRIM('retv_vars.10656_0570_002.cdf')
      TES(999)%FILENAME = TRIM('retv_vars.10656_0570_003.cdf')
      TES(1000)%FILENAME = TRIM('retv_vars.10656_0570_004.cdf')
      TES(1001)%FILENAME = TRIM('retv_vars.10656_0571_002.cdf')
      TES(1002)%FILENAME = TRIM('retv_vars.10656_0581_002.cdf')
      TES(1003)%FILENAME = TRIM('retv_vars.10656_0583_002.cdf')
      TES(1004)%FILENAME = TRIM('retv_vars.10656_0583_004.cdf')
      TES(1005)%FILENAME = TRIM('retv_vars.10656_0594_003.cdf')
      TES(1006)%FILENAME = TRIM('retv_vars.10656_0597_002.cdf')
      TES(1007)%FILENAME = TRIM('retv_vars.10656_0597_003.cdf')
      TES(1008)%FILENAME = TRIM('retv_vars.10656_0598_002.cdf')
      TES(1009)%FILENAME = TRIM('retv_vars.10656_0599_002.cdf')
      TES(1010)%FILENAME = TRIM('retv_vars.10656_0613_003.cdf')
      TES(1011)%FILENAME = TRIM('retv_vars.10656_0614_002.cdf')
      TES(1012)%FILENAME = TRIM('retv_vars.10656_0615_004.cdf')
      TES(1013)%FILENAME = TRIM('retv_vars.10656_0616_002.cdf')
      TES(1014)%FILENAME = TRIM('retv_vars.10656_0616_003.cdf')
      TES(1015)%FILENAME = TRIM('retv_vars.10656_0617_004.cdf')
      TES(1016)%FILENAME = TRIM('retv_vars.10656_0618_002.cdf')
      TES(1017)%FILENAME = TRIM('retv_vars.10656_0619_002.cdf')
      TES(1018)%FILENAME = TRIM('retv_vars.10656_0619_003.cdf')
      TES(1019)%FILENAME = TRIM('retv_vars.10656_0620_004.cdf')
      TES(1020)%FILENAME = TRIM('retv_vars.10656_0621_004.cdf')
      TES(1021)%FILENAME = TRIM('retv_vars.10656_0622_002.cdf')
      TES(1022)%FILENAME = TRIM('retv_vars.10656_0622_003.cdf')
      TES(1023)%FILENAME = TRIM('retv_vars.10656_0622_004.cdf')
      TES(1024)%FILENAME = TRIM('retv_vars.10656_0623_002.cdf')
      TES(1025)%FILENAME = TRIM('retv_vars.10656_0623_003.cdf')
      TES(1026)%FILENAME = TRIM('retv_vars.10656_0623_004.cdf')
      TES(1027)%FILENAME = TRIM('retv_vars.10656_0624_002.cdf')
      TES(1028)%FILENAME = TRIM('retv_vars.10656_0633_004.cdf')
      TES(1029)%FILENAME = TRIM('retv_vars.10656_0634_003.cdf')
      TES(1030)%FILENAME = TRIM('retv_vars.10656_0637_002.cdf')
      TES(1031)%FILENAME = TRIM('retv_vars.10656_0638_002.cdf')
      TES(1032)%FILENAME = TRIM('retv_vars.10656_0639_002.cdf')
      TES(1033)%FILENAME = TRIM('retv_vars.10656_0640_002.cdf')
      TES(1034)%FILENAME = TRIM('retv_vars.10656_0641_002.cdf')
      TES(1035)%FILENAME = TRIM('retv_vars.10656_0642_004.cdf')
      TES(1036)%FILENAME = TRIM('retv_vars.10656_0643_002.cdf')
      TES(1037)%FILENAME = TRIM('retv_vars.10656_0644_003.cdf')
      TES(1038)%FILENAME = TRIM('retv_vars.10656_0645_002.cdf')
      TES(1039)%FILENAME = TRIM('retv_vars.10656_0645_004.cdf')
      TES(1040)%FILENAME = TRIM('retv_vars.10656_0646_002.cdf')
      TES(1041)%FILENAME = TRIM('retv_vars.10656_0646_003.cdf')
      TES(1042)%FILENAME = TRIM('retv_vars.10656_0646_004.cdf')
      TES(1043)%FILENAME = TRIM('retv_vars.10656_0647_002.cdf')
      TES(1044)%FILENAME = TRIM('retv_vars.10656_0652_002.cdf')
      TES(1045)%FILENAME = TRIM('retv_vars.10656_0652_003.cdf')
      TES(1046)%FILENAME = TRIM('retv_vars.10656_0652_004.cdf')
      TES(1047)%FILENAME = TRIM('retv_vars.10656_0653_002.cdf')
      TES(1048)%FILENAME = TRIM('retv_vars.10656_0653_004.cdf')
      TES(1049)%FILENAME = TRIM('retv_vars.10656_0654_002.cdf')
      TES(1050)%FILENAME = TRIM('retv_vars.10656_0654_003.cdf')
      TES(1051)%FILENAME = TRIM('retv_vars.10656_0654_004.cdf')
      TES(1052)%FILENAME = TRIM('retv_vars.10656_0655_002.cdf')
      TES(1053)%FILENAME = TRIM('retv_vars.10656_0655_004.cdf')
      TES(1054)%FILENAME = TRIM('retv_vars.10656_0656_002.cdf')
      TES(1055)%FILENAME = TRIM('retv_vars.10656_0658_002.cdf')
      TES(1056)%FILENAME = TRIM('retv_vars.10656_0684_004.cdf')
      TES(1057)%FILENAME = TRIM('retv_vars.10656_0687_004.cdf')
      TES(1058)%FILENAME = TRIM('retv_vars.10656_0688_002.cdf')
      TES(1059)%FILENAME = TRIM('retv_vars.10656_0688_003.cdf')
      TES(1060)%FILENAME = TRIM('retv_vars.10656_0689_002.cdf')
      TES(1061)%FILENAME = TRIM('retv_vars.10656_0691_002.cdf')
      TES(1062)%FILENAME = TRIM('retv_vars.10656_0693_002.cdf')
      TES(1063)%FILENAME = TRIM('retv_vars.10656_0693_003.cdf')
      TES(1064)%FILENAME = TRIM('retv_vars.10656_0694_003.cdf')
      TES(1065)%FILENAME = TRIM('retv_vars.10656_0701_002.cdf')
      TES(1066)%FILENAME = TRIM('retv_vars.10656_0701_004.cdf')
      TES(1067)%FILENAME = TRIM('retv_vars.10656_0703_003.cdf')
      TES(1068)%FILENAME = TRIM('retv_vars.10656_0704_002.cdf')
      TES(1069)%FILENAME = TRIM('retv_vars.10656_0705_002.cdf')
      TES(1070)%FILENAME = TRIM('retv_vars.10656_0705_003.cdf')
      TES(1071)%FILENAME = TRIM('retv_vars.10656_0705_004.cdf')
      TES(1072)%FILENAME = TRIM('retv_vars.10656_0706_002.cdf')
      TES(1073)%FILENAME = TRIM('retv_vars.10656_0735_004.cdf')
      TES(1074)%FILENAME = TRIM('retv_vars.10656_0736_003.cdf')
      TES(1075)%FILENAME = TRIM('retv_vars.10656_0737_004.cdf')
      TES(1076)%FILENAME = TRIM('retv_vars.10656_0738_004.cdf')
      TES(1077)%FILENAME = TRIM('retv_vars.10656_0739_002.cdf')
      TES(1078)%FILENAME = TRIM('retv_vars.10656_0739_003.cdf')
      TES(1079)%FILENAME = TRIM('retv_vars.10656_0740_002.cdf')
      TES(1080)%FILENAME = TRIM('retv_vars.10656_0741_003.cdf')
      TES(1081)%FILENAME = TRIM('retv_vars.10656_0741_004.cdf')
      TES(1082)%FILENAME = TRIM('retv_vars.10656_0742_002.cdf')
      TES(1083)%FILENAME = TRIM('retv_vars.10656_0742_003.cdf')
      TES(1084)%FILENAME = TRIM('retv_vars.10656_0747_003.cdf')
      TES(1085)%FILENAME = TRIM('retv_vars.10656_0747_004.cdf')
      TES(1086)%FILENAME = TRIM('retv_vars.10656_0748_002.cdf')
      TES(1087)%FILENAME = TRIM('retv_vars.10656_0748_004.cdf')
      TES(1088)%FILENAME = TRIM('retv_vars.10656_0749_003.cdf')
      TES(1089)%FILENAME = TRIM('retv_vars.10656_0749_004.cdf')
      TES(1090)%FILENAME = TRIM('retv_vars.10656_0750_002.cdf')
      TES(1091)%FILENAME = TRIM('retv_vars.10658_0020_003.cdf')
      TES(1092)%FILENAME = TRIM('retv_vars.10658_0020_004.cdf')
      TES(1093)%FILENAME = TRIM('retv_vars.10658_0021_004.cdf')
      TES(1094)%FILENAME = TRIM('retv_vars.10658_0022_004.cdf')
      TES(1095)%FILENAME = TRIM('retv_vars.10658_0023_002.cdf')
      TES(1096)%FILENAME = TRIM('retv_vars.10658_0027_002.cdf')
      TES(1097)%FILENAME = TRIM('retv_vars.10658_0027_003.cdf')
      TES(1098)%FILENAME = TRIM('retv_vars.10658_0028_003.cdf')
      TES(1099)%FILENAME = TRIM('retv_vars.10658_0028_004.cdf')
      TES(1100)%FILENAME = TRIM('retv_vars.10658_0030_002.cdf')
      TES(1101)%FILENAME = TRIM('retv_vars.10658_0030_003.cdf')
      TES(1102)%FILENAME = TRIM('retv_vars.10658_0031_002.cdf')
      TES(1103)%FILENAME = TRIM('retv_vars.10658_0031_003.cdf')
      TES(1104)%FILENAME = TRIM('retv_vars.10658_0031_004.cdf')
      TES(1105)%FILENAME = TRIM('retv_vars.10658_0032_002.cdf')
      TES(1106)%FILENAME = TRIM('retv_vars.10658_0033_002.cdf')
      TES(1107)%FILENAME = TRIM('retv_vars.10658_0059_004.cdf')
      TES(1108)%FILENAME = TRIM('retv_vars.10658_0060_003.cdf')
      TES(1109)%FILENAME = TRIM('retv_vars.10658_0061_002.cdf')
      TES(1110)%FILENAME = TRIM('retv_vars.10658_0064_004.cdf')
      TES(1111)%FILENAME = TRIM('retv_vars.10658_0066_002.cdf')
      TES(1112)%FILENAME = TRIM('retv_vars.10658_0066_004.cdf')
      TES(1113)%FILENAME = TRIM('retv_vars.10658_0068_002.cdf')
      TES(1114)%FILENAME = TRIM('retv_vars.10658_0068_004.cdf')
      TES(1115)%FILENAME = TRIM('retv_vars.10658_0069_002.cdf')
      TES(1116)%FILENAME = TRIM('retv_vars.10658_0069_003.cdf')
      TES(1117)%FILENAME = TRIM('retv_vars.10658_0070_002.cdf')
      TES(1118)%FILENAME = TRIM('retv_vars.10658_0070_003.cdf')
      TES(1119)%FILENAME = TRIM('retv_vars.10658_0075_002.cdf')
      TES(1120)%FILENAME = TRIM('retv_vars.10658_0075_003.cdf')
      TES(1121)%FILENAME = TRIM('retv_vars.10658_0076_002.cdf')
      TES(1122)%FILENAME = TRIM('retv_vars.10658_0076_003.cdf')
      TES(1123)%FILENAME = TRIM('retv_vars.10658_0101_002.cdf')
      TES(1124)%FILENAME = TRIM('retv_vars.10658_0101_003.cdf')
      TES(1125)%FILENAME = TRIM('retv_vars.10658_0102_004.cdf')
      TES(1126)%FILENAME = TRIM('retv_vars.10658_0103_004.cdf')
      TES(1127)%FILENAME = TRIM('retv_vars.10658_0104_002.cdf')
      TES(1128)%FILENAME = TRIM('retv_vars.10658_0104_003.cdf')
      TES(1129)%FILENAME = TRIM('retv_vars.10658_0104_004.cdf')
      TES(1130)%FILENAME = TRIM('retv_vars.10658_0105_004.cdf')
      TES(1131)%FILENAME = TRIM('retv_vars.10658_0106_002.cdf')
      TES(1132)%FILENAME = TRIM('retv_vars.10658_0106_003.cdf')
      TES(1133)%FILENAME = TRIM('retv_vars.10658_0108_003.cdf')
      TES(1134)%FILENAME = TRIM('retv_vars.10658_0108_004.cdf')
      TES(1135)%FILENAME = TRIM('retv_vars.10658_0109_002.cdf')
      TES(1136)%FILENAME = TRIM('retv_vars.10658_0109_003.cdf')
      TES(1137)%FILENAME = TRIM('retv_vars.10658_0112_004.cdf')
      TES(1138)%FILENAME = TRIM('retv_vars.10658_0114_004.cdf')
      TES(1139)%FILENAME = TRIM('retv_vars.10658_0115_004.cdf')
      TES(1140)%FILENAME = TRIM('retv_vars.10658_0116_003.cdf')
      TES(1141)%FILENAME = TRIM('retv_vars.10658_0117_003.cdf')
      TES(1142)%FILENAME = TRIM('retv_vars.10658_0232_004.cdf')
      TES(1143)%FILENAME = TRIM('retv_vars.10658_0234_003.cdf')
      TES(1144)%FILENAME = TRIM('retv_vars.10658_0235_004.cdf')
      TES(1145)%FILENAME = TRIM('retv_vars.10658_0237_003.cdf')
      TES(1146)%FILENAME = TRIM('retv_vars.10658_0247_002.cdf')
      TES(1147)%FILENAME = TRIM('retv_vars.10658_0247_003.cdf')
      TES(1148)%FILENAME = TRIM('retv_vars.10658_0247_004.cdf')
      TES(1149)%FILENAME = TRIM('retv_vars.10658_0248_002.cdf')
      TES(1150)%FILENAME = TRIM('retv_vars.10658_0248_003.cdf')
      TES(1151)%FILENAME = TRIM('retv_vars.10658_0248_004.cdf')
      TES(1152)%FILENAME = TRIM('retv_vars.10658_0249_002.cdf')
      TES(1153)%FILENAME = TRIM('retv_vars.10658_0249_003.cdf')
      TES(1154)%FILENAME = TRIM('retv_vars.10658_0249_004.cdf')
      TES(1155)%FILENAME = TRIM('retv_vars.10658_0250_002.cdf')
      TES(1156)%FILENAME = TRIM('retv_vars.10658_0250_003.cdf')
      TES(1157)%FILENAME = TRIM('retv_vars.10658_0250_004.cdf')
      TES(1158)%FILENAME = TRIM('retv_vars.10658_0251_002.cdf')
      TES(1159)%FILENAME = TRIM('retv_vars.10658_0251_003.cdf')
      TES(1160)%FILENAME = TRIM('retv_vars.10658_0251_004.cdf')
      TES(1161)%FILENAME = TRIM('retv_vars.10658_0260_003.cdf')
      TES(1162)%FILENAME = TRIM('retv_vars.10658_0260_004.cdf')
      TES(1163)%FILENAME = TRIM('retv_vars.10658_0262_004.cdf')
      TES(1164)%FILENAME = TRIM('retv_vars.10658_0267_004.cdf')
      TES(1165)%FILENAME = TRIM('retv_vars.10658_0268_002.cdf')
      TES(1166)%FILENAME = TRIM('retv_vars.10658_0268_004.cdf')
      TES(1167)%FILENAME = TRIM('retv_vars.10658_0269_003.cdf')
      TES(1168)%FILENAME = TRIM('retv_vars.10658_0269_004.cdf')
      TES(1169)%FILENAME = TRIM('retv_vars.10658_0270_002.cdf')
      TES(1170)%FILENAME = TRIM('retv_vars.10658_0271_002.cdf')
      TES(1171)%FILENAME = TRIM('retv_vars.10658_0273_003.cdf')
      TES(1172)%FILENAME = TRIM('retv_vars.10658_0278_003.cdf')
      TES(1173)%FILENAME = TRIM('retv_vars.10658_0278_004.cdf')
      TES(1174)%FILENAME = TRIM('retv_vars.10658_0279_002.cdf')
      TES(1175)%FILENAME = TRIM('retv_vars.10658_0279_003.cdf')
      TES(1176)%FILENAME = TRIM('retv_vars.10658_0280_004.cdf')
      TES(1177)%FILENAME = TRIM('retv_vars.10658_0291_002.cdf')
      TES(1178)%FILENAME = TRIM('retv_vars.10658_0291_003.cdf')
      TES(1179)%FILENAME = TRIM('retv_vars.10658_0292_002.cdf')
      TES(1180)%FILENAME = TRIM('retv_vars.10658_0293_002.cdf')
      TES(1181)%FILENAME = TRIM('retv_vars.10658_0296_002.cdf')
      TES(1182)%FILENAME = TRIM('retv_vars.10658_0297_003.cdf')
      TES(1183)%FILENAME = TRIM('retv_vars.10658_0297_004.cdf')
      TES(1184)%FILENAME = TRIM('retv_vars.10658_0298_002.cdf')
      TES(1185)%FILENAME = TRIM('retv_vars.10658_0298_003.cdf')
      TES(1186)%FILENAME = TRIM('retv_vars.10658_0298_004.cdf')
      TES(1187)%FILENAME = TRIM('retv_vars.10658_0299_002.cdf')
      TES(1188)%FILENAME = TRIM('retv_vars.10658_0299_003.cdf')
      TES(1189)%FILENAME = TRIM('retv_vars.10658_0303_004.cdf')
      TES(1190)%FILENAME = TRIM('retv_vars.10658_0305_004.cdf')
      TES(1191)%FILENAME = TRIM('retv_vars.10658_0307_003.cdf')
      TES(1192)%FILENAME = TRIM('retv_vars.10658_0308_003.cdf')
      TES(1193)%FILENAME = TRIM('retv_vars.10658_0309_002.cdf')
      TES(1194)%FILENAME = TRIM('retv_vars.10658_0309_003.cdf')
      TES(1195)%FILENAME = TRIM('retv_vars.10658_0315_002.cdf')
      TES(1196)%FILENAME = TRIM('retv_vars.10658_0315_004.cdf')
      TES(1197)%FILENAME = TRIM('retv_vars.10658_0316_003.cdf')
      TES(1198)%FILENAME = TRIM('retv_vars.10658_0316_004.cdf')
      TES(1199)%FILENAME = TRIM('retv_vars.10658_0317_004.cdf')
      TES(1200)%FILENAME = TRIM('retv_vars.10658_0322_002.cdf')
      TES(1201)%FILENAME = TRIM('retv_vars.10658_0322_003.cdf')
      TES(1202)%FILENAME = TRIM('retv_vars.10658_0323_002.cdf')
      TES(1203)%FILENAME = TRIM('retv_vars.10658_0324_004.cdf')
      TES(1204)%FILENAME = TRIM('retv_vars.10658_0325_002.cdf')
      TES(1205)%FILENAME = TRIM('retv_vars.10658_0353_002.cdf')
      TES(1206)%FILENAME = TRIM('retv_vars.10658_0353_003.cdf')
      TES(1207)%FILENAME = TRIM('retv_vars.10658_0354_002.cdf')
      TES(1208)%FILENAME = TRIM('retv_vars.10658_0354_004.cdf')
      TES(1209)%FILENAME = TRIM('retv_vars.10658_0355_004.cdf')
      TES(1210)%FILENAME = TRIM('retv_vars.10658_0357_002.cdf')
      TES(1211)%FILENAME = TRIM('retv_vars.10658_0357_003.cdf')
      TES(1212)%FILENAME = TRIM('retv_vars.10658_0358_002.cdf')
      TES(1213)%FILENAME = TRIM('retv_vars.10658_0359_002.cdf')
      TES(1214)%FILENAME = TRIM('retv_vars.10658_0359_003.cdf')
      TES(1215)%FILENAME = TRIM('retv_vars.10658_0363_003.cdf')
      TES(1216)%FILENAME = TRIM('retv_vars.10658_0364_002.cdf')
      TES(1217)%FILENAME = TRIM('retv_vars.10658_0365_002.cdf')
      TES(1218)%FILENAME = TRIM('retv_vars.10658_0366_002.cdf')
      TES(1219)%FILENAME = TRIM('retv_vars.10658_0366_004.cdf')
      TES(1220)%FILENAME = TRIM('retv_vars.10658_0368_002.cdf')
      TES(1221)%FILENAME = TRIM('retv_vars.10658_0368_004.cdf')
      TES(1222)%FILENAME = TRIM('retv_vars.10658_0369_002.cdf')
      TES(1223)%FILENAME = TRIM('retv_vars.10658_0371_002.cdf')
      TES(1224)%FILENAME = TRIM('retv_vars.10658_0411_002.cdf')
      TES(1225)%FILENAME = TRIM('retv_vars.10658_0411_003.cdf')
      TES(1226)%FILENAME = TRIM('retv_vars.10658_0411_004.cdf')
      TES(1227)%FILENAME = TRIM('retv_vars.10658_0412_002.cdf')
      TES(1228)%FILENAME = TRIM('retv_vars.10658_0412_003.cdf')
      TES(1229)%FILENAME = TRIM('retv_vars.10658_0412_004.cdf')
      TES(1230)%FILENAME = TRIM('retv_vars.10658_0418_002.cdf')
      TES(1231)%FILENAME = TRIM('retv_vars.10658_0418_004.cdf')
      TES(1232)%FILENAME = TRIM('retv_vars.10658_0419_002.cdf')
      TES(1233)%FILENAME = TRIM('retv_vars.10658_0419_004.cdf')
      TES(1234)%FILENAME = TRIM('retv_vars.10658_0421_004.cdf')
      TES(1235)%FILENAME = TRIM('retv_vars.10658_0422_002.cdf')
      TES(1236)%FILENAME = TRIM('retv_vars.10658_0422_003.cdf')
      TES(1237)%FILENAME = TRIM('retv_vars.10658_0423_003.cdf')
      TES(1238)%FILENAME = TRIM('retv_vars.10658_0423_004.cdf')
      TES(1239)%FILENAME = TRIM('retv_vars.10658_0424_003.cdf')
      TES(1240)%FILENAME = TRIM('retv_vars.10658_0424_004.cdf')
      TES(1241)%FILENAME = TRIM('retv_vars.10658_0425_002.cdf')
      TES(1242)%FILENAME = TRIM('retv_vars.10658_0425_003.cdf')
      TES(1243)%FILENAME = TRIM('retv_vars.10658_0426_002.cdf')
      TES(1244)%FILENAME = TRIM('retv_vars.10658_0426_003.cdf')
      TES(1245)%FILENAME = TRIM('retv_vars.10658_0426_004.cdf')
      TES(1246)%FILENAME = TRIM('retv_vars.10658_0427_004.cdf')
      TES(1247)%FILENAME = TRIM('retv_vars.10658_0459_002.cdf')
      TES(1248)%FILENAME = TRIM('retv_vars.10658_0459_003.cdf')
      TES(1249)%FILENAME = TRIM('retv_vars.10658_0459_004.cdf')
      TES(1250)%FILENAME = TRIM('retv_vars.10658_0460_002.cdf')
      TES(1251)%FILENAME = TRIM('retv_vars.10658_0460_004.cdf')
      TES(1252)%FILENAME = TRIM('retv_vars.10658_0461_002.cdf')
      TES(1253)%FILENAME = TRIM('retv_vars.10658_0462_003.cdf')
      TES(1254)%FILENAME = TRIM('retv_vars.10658_0462_004.cdf')
      TES(1255)%FILENAME = TRIM('retv_vars.10658_0465_003.cdf')
      TES(1256)%FILENAME = TRIM('retv_vars.10658_0465_004.cdf')
      TES(1257)%FILENAME = TRIM('retv_vars.10658_0466_002.cdf')
      TES(1258)%FILENAME = TRIM('retv_vars.10658_0466_003.cdf')
      TES(1259)%FILENAME = TRIM('retv_vars.10658_0467_002.cdf')
      TES(1260)%FILENAME = TRIM('retv_vars.10658_0469_003.cdf')
      TES(1261)%FILENAME = TRIM('retv_vars.10658_0469_004.cdf')
      TES(1262)%FILENAME = TRIM('retv_vars.10658_0470_002.cdf')
      TES(1263)%FILENAME = TRIM('retv_vars.10658_0472_002.cdf')
      TES(1264)%FILENAME = TRIM('retv_vars.10658_0473_002.cdf')
      TES(1265)%FILENAME = TRIM('retv_vars.10658_0473_004.cdf')
      TES(1266)%FILENAME = TRIM('retv_vars.10658_0474_002.cdf')
      TES(1267)%FILENAME = TRIM('retv_vars.10658_0474_003.cdf')
      TES(1268)%FILENAME = TRIM('retv_vars.10658_0507_002.cdf')
      TES(1269)%FILENAME = TRIM('retv_vars.10658_0509_004.cdf')
      TES(1270)%FILENAME = TRIM('retv_vars.10658_0510_002.cdf')
      TES(1271)%FILENAME = TRIM('retv_vars.10658_0512_002.cdf')
      TES(1272)%FILENAME = TRIM('retv_vars.10658_0513_002.cdf')
      TES(1273)%FILENAME = TRIM('retv_vars.10658_0518_002.cdf')
      TES(1274)%FILENAME = TRIM('retv_vars.10658_0535_003.cdf')
      TES(1275)%FILENAME = TRIM('retv_vars.10658_0569_002.cdf')
      TES(1276)%FILENAME = TRIM('retv_vars.10658_0579_002.cdf')
      TES(1277)%FILENAME = TRIM('retv_vars.10658_0582_003.cdf')
      TES(1278)%FILENAME = TRIM('retv_vars.10658_0582_004.cdf')
      TES(1279)%FILENAME = TRIM('retv_vars.10658_0583_002.cdf')
      TES(1280)%FILENAME = TRIM('retv_vars.10658_0583_003.cdf')
      TES(1281)%FILENAME = TRIM('retv_vars.10658_0584_003.cdf')
      TES(1282)%FILENAME = TRIM('retv_vars.10658_0596_003.cdf')
      TES(1283)%FILENAME = TRIM('retv_vars.10658_0599_002.cdf')
      TES(1284)%FILENAME = TRIM('retv_vars.10658_0614_002.cdf')
      TES(1285)%FILENAME = TRIM('retv_vars.10658_0614_003.cdf')
      TES(1286)%FILENAME = TRIM('retv_vars.10658_0614_004.cdf')
      TES(1287)%FILENAME = TRIM('retv_vars.10658_0615_002.cdf')
      TES(1288)%FILENAME = TRIM('retv_vars.10658_0615_003.cdf')
      TES(1289)%FILENAME = TRIM('retv_vars.10658_0615_004.cdf')
      TES(1290)%FILENAME = TRIM('retv_vars.10658_0616_002.cdf')
      TES(1291)%FILENAME = TRIM('retv_vars.10658_0616_003.cdf')
      TES(1292)%FILENAME = TRIM('retv_vars.10658_0616_004.cdf')
      TES(1293)%FILENAME = TRIM('retv_vars.10658_0617_002.cdf')
      TES(1294)%FILENAME = TRIM('retv_vars.10658_0617_004.cdf')
      TES(1295)%FILENAME = TRIM('retv_vars.10658_0618_004.cdf')
      TES(1296)%FILENAME = TRIM('retv_vars.10658_0620_004.cdf')
      TES(1297)%FILENAME = TRIM('retv_vars.10658_0621_002.cdf')
      TES(1298)%FILENAME = TRIM('retv_vars.10658_0621_003.cdf')
      TES(1299)%FILENAME = TRIM('retv_vars.10658_0621_004.cdf')
      TES(1300)%FILENAME = TRIM('retv_vars.10658_0622_002.cdf')
      TES(1301)%FILENAME = TRIM('retv_vars.10658_0622_003.cdf')
      TES(1302)%FILENAME = TRIM('retv_vars.10658_0622_004.cdf')
      TES(1303)%FILENAME = TRIM('retv_vars.10658_0623_003.cdf')
      TES(1304)%FILENAME = TRIM('retv_vars.10658_0623_004.cdf')
      TES(1305)%FILENAME = TRIM('retv_vars.10658_0624_002.cdf')
      TES(1306)%FILENAME = TRIM('retv_vars.10658_0624_003.cdf')
      TES(1307)%FILENAME = TRIM('retv_vars.10658_0628_002.cdf')
      TES(1308)%FILENAME = TRIM('retv_vars.10658_0628_003.cdf')
      TES(1309)%FILENAME = TRIM('retv_vars.10658_0628_004.cdf')
      TES(1310)%FILENAME = TRIM('retv_vars.10658_0629_002.cdf')
      TES(1311)%FILENAME = TRIM('retv_vars.10658_0629_003.cdf')
      TES(1312)%FILENAME = TRIM('retv_vars.10658_0630_002.cdf')
      TES(1313)%FILENAME = TRIM('retv_vars.10658_0630_003.cdf')
      TES(1314)%FILENAME = TRIM('retv_vars.10658_0634_004.cdf')
      TES(1315)%FILENAME = TRIM('retv_vars.10658_0635_002.cdf')
      TES(1316)%FILENAME = TRIM('retv_vars.10658_0637_004.cdf')
      TES(1317)%FILENAME = TRIM('retv_vars.10658_0639_004.cdf')
      TES(1318)%FILENAME = TRIM('retv_vars.10658_0640_002.cdf')
      TES(1319)%FILENAME = TRIM('retv_vars.10658_0640_003.cdf')
      TES(1320)%FILENAME = TRIM('retv_vars.10658_0641_002.cdf')
      TES(1321)%FILENAME = TRIM('retv_vars.10658_0644_004.cdf')
      TES(1322)%FILENAME = TRIM('retv_vars.10658_0645_002.cdf')
      TES(1323)%FILENAME = TRIM('retv_vars.10658_0645_003.cdf')
      TES(1324)%FILENAME = TRIM('retv_vars.10658_0645_004.cdf')
      TES(1325)%FILENAME = TRIM('retv_vars.10658_0646_002.cdf')
      TES(1326)%FILENAME = TRIM('retv_vars.10658_0646_003.cdf')
      TES(1327)%FILENAME = TRIM('retv_vars.10658_0646_004.cdf')
      TES(1328)%FILENAME = TRIM('retv_vars.10658_0647_002.cdf')
      TES(1329)%FILENAME = TRIM('retv_vars.10658_0654_004.cdf')
      TES(1330)%FILENAME = TRIM('retv_vars.10658_0658_002.cdf')
      TES(1331)%FILENAME = TRIM('retv_vars.10658_0658_004.cdf')
      TES(1332)%FILENAME = TRIM('retv_vars.10658_0687_004.cdf')
      TES(1333)%FILENAME = TRIM('retv_vars.10658_0688_002.cdf')
      TES(1334)%FILENAME = TRIM('retv_vars.10658_0688_004.cdf')
      TES(1335)%FILENAME = TRIM('retv_vars.10658_0689_002.cdf')
      TES(1336)%FILENAME = TRIM('retv_vars.10658_0691_003.cdf')
      TES(1337)%FILENAME = TRIM('retv_vars.10658_0691_004.cdf')
      TES(1338)%FILENAME = TRIM('retv_vars.10658_0692_003.cdf')
      TES(1339)%FILENAME = TRIM('retv_vars.10658_0693_002.cdf')
      TES(1340)%FILENAME = TRIM('retv_vars.10658_0693_003.cdf')
      TES(1341)%FILENAME = TRIM('retv_vars.10658_0694_004.cdf')
      TES(1342)%FILENAME = TRIM('retv_vars.10658_0695_002.cdf')
      TES(1343)%FILENAME = TRIM('retv_vars.10658_0700_003.cdf')
      TES(1344)%FILENAME = TRIM('retv_vars.10658_0700_004.cdf')
      TES(1345)%FILENAME = TRIM('retv_vars.10658_0701_002.cdf')
      TES(1346)%FILENAME = TRIM('retv_vars.10658_0702_002.cdf')
      TES(1347)%FILENAME = TRIM('retv_vars.10658_0702_004.cdf')
      TES(1348)%FILENAME = TRIM('retv_vars.10658_0704_002.cdf')
      TES(1349)%FILENAME = TRIM('retv_vars.10658_0704_003.cdf')
      TES(1350)%FILENAME = TRIM('retv_vars.10658_0704_004.cdf')
      TES(1351)%FILENAME = TRIM('retv_vars.10658_0705_003.cdf')
      TES(1352)%FILENAME = TRIM('retv_vars.10658_0705_004.cdf')
      TES(1353)%FILENAME = TRIM('retv_vars.10658_0706_002.cdf')
      TES(1354)%FILENAME = TRIM('retv_vars.10658_0706_003.cdf')
      TES(1355)%FILENAME = TRIM('retv_vars.10658_0706_004.cdf')
      TES(1356)%FILENAME = TRIM('retv_vars.10658_0741_003.cdf')
      TES(1357)%FILENAME = TRIM('retv_vars.10658_0741_004.cdf')
      TES(1358)%FILENAME = TRIM('retv_vars.10658_0742_002.cdf')
      TES(1359)%FILENAME = TRIM('retv_vars.10658_0742_004.cdf')
      TES(1360)%FILENAME = TRIM('retv_vars.10658_0743_002.cdf')
      TES(1361)%FILENAME = TRIM('retv_vars.10658_0747_003.cdf')
      TES(1362)%FILENAME = TRIM('retv_vars.10658_0747_004.cdf')
      TES(1363)%FILENAME = TRIM('retv_vars.10658_0748_002.cdf')
      TES(1364)%FILENAME = TRIM('retv_vars.10658_0748_003.cdf')
      TES(1365)%FILENAME = TRIM('retv_vars.10658_0748_004.cdf')
      TES(1366)%FILENAME = TRIM('retv_vars.10658_0749_002.cdf')
      TES(1367)%FILENAME = TRIM('retv_vars.10658_0749_004.cdf')
      TES(1368)%FILENAME = TRIM('retv_vars.10658_0750_002.cdf')
      TES(1369)%FILENAME = TRIM('retv_vars.10658_0750_003.cdf')
      TES(1370)%FILENAME = TRIM('retv_vars.10658_0750_004.cdf')
      TES(1371)%FILENAME = TRIM('retv_vars.10658_0751_002.cdf')
      TES(1372)%FILENAME = TRIM('retv_vars.10658_0751_003.cdf')
      TES(1373)%FILENAME = TRIM('retv_vars.10666_0012_003.cdf')
      TES(1374)%FILENAME = TRIM('retv_vars.10666_0013_002.cdf')
      TES(1375)%FILENAME = TRIM('retv_vars.10666_0020_004.cdf')
      TES(1376)%FILENAME = TRIM('retv_vars.10666_0021_003.cdf')
      TES(1377)%FILENAME = TRIM('retv_vars.10666_0021_004.cdf')
      TES(1378)%FILENAME = TRIM('retv_vars.10666_0022_003.cdf')
      TES(1379)%FILENAME = TRIM('retv_vars.10666_0027_004.cdf')
      TES(1380)%FILENAME = TRIM('retv_vars.10666_0028_002.cdf')
      TES(1381)%FILENAME = TRIM('retv_vars.10666_0053_002.cdf')
      TES(1382)%FILENAME = TRIM('retv_vars.10666_0053_003.cdf')
      TES(1383)%FILENAME = TRIM('retv_vars.10666_0053_004.cdf')
      TES(1384)%FILENAME = TRIM('retv_vars.10666_0054_002.cdf')
      TES(1385)%FILENAME = TRIM('retv_vars.10666_0054_003.cdf')
      TES(1386)%FILENAME = TRIM('retv_vars.10666_0054_004.cdf')
      TES(1387)%FILENAME = TRIM('retv_vars.10666_0055_003.cdf')
      TES(1388)%FILENAME = TRIM('retv_vars.10666_0055_004.cdf')
      TES(1389)%FILENAME = TRIM('retv_vars.10666_0056_002.cdf')
      TES(1390)%FILENAME = TRIM('retv_vars.10666_0056_003.cdf')
      TES(1391)%FILENAME = TRIM('retv_vars.10666_0057_004.cdf')
      TES(1392)%FILENAME = TRIM('retv_vars.10666_0059_002.cdf')
      TES(1393)%FILENAME = TRIM('retv_vars.10666_0059_004.cdf')
      TES(1394)%FILENAME = TRIM('retv_vars.10666_0060_003.cdf')
      TES(1395)%FILENAME = TRIM('retv_vars.10666_0067_003.cdf')
      TES(1396)%FILENAME = TRIM('retv_vars.10666_0067_004.cdf')
      TES(1397)%FILENAME = TRIM('retv_vars.10666_0068_004.cdf')
      TES(1398)%FILENAME = TRIM('retv_vars.10666_0069_003.cdf')
      TES(1399)%FILENAME = TRIM('retv_vars.10666_0069_004.cdf')
      TES(1400)%FILENAME = TRIM('retv_vars.10666_0070_003.cdf')
      TES(1401)%FILENAME = TRIM('retv_vars.10666_0109_003.cdf')
      TES(1402)%FILENAME = TRIM('retv_vars.10666_0110_003.cdf')
      TES(1403)%FILENAME = TRIM('retv_vars.10666_0113_004.cdf')
      TES(1404)%FILENAME = TRIM('retv_vars.10666_0133_004.cdf')
      TES(1405)%FILENAME = TRIM('retv_vars.10666_0171_003.cdf')
      TES(1406)%FILENAME = TRIM('retv_vars.10666_0172_002.cdf')
      TES(1407)%FILENAME = TRIM('retv_vars.10666_0172_003.cdf')
      TES(1408)%FILENAME = TRIM('retv_vars.10666_0172_004.cdf')
      TES(1409)%FILENAME = TRIM('retv_vars.10666_0184_004.cdf')
      TES(1410)%FILENAME = TRIM('retv_vars.10666_0185_002.cdf')
      TES(1411)%FILENAME = TRIM('retv_vars.10666_0186_004.cdf')
      TES(1412)%FILENAME = TRIM('retv_vars.10666_0187_002.cdf')
      TES(1413)%FILENAME = TRIM('retv_vars.10666_0187_003.cdf')
      TES(1414)%FILENAME = TRIM('retv_vars.10666_0187_004.cdf')
      TES(1415)%FILENAME = TRIM('retv_vars.10666_0190_002.cdf')
      TES(1416)%FILENAME = TRIM('retv_vars.10666_0198_003.cdf')
      TES(1417)%FILENAME = TRIM('retv_vars.10666_0198_004.cdf')
      TES(1418)%FILENAME = TRIM('retv_vars.10666_0199_002.cdf')
      TES(1419)%FILENAME = TRIM('retv_vars.10666_0199_003.cdf')
      TES(1420)%FILENAME = TRIM('retv_vars.10666_0199_004.cdf')
      TES(1421)%FILENAME = TRIM('retv_vars.10666_0201_002.cdf')
      TES(1422)%FILENAME = TRIM('retv_vars.10666_0201_004.cdf')
      TES(1423)%FILENAME = TRIM('retv_vars.10666_0202_002.cdf')
      TES(1424)%FILENAME = TRIM('retv_vars.10666_0202_003.cdf')
      TES(1425)%FILENAME = TRIM('retv_vars.10666_0202_004.cdf')
      TES(1426)%FILENAME = TRIM('retv_vars.10666_0212_002.cdf')
      TES(1427)%FILENAME = TRIM('retv_vars.10666_0213_004.cdf')
      TES(1428)%FILENAME = TRIM('retv_vars.10666_0219_002.cdf')
      TES(1429)%FILENAME = TRIM('retv_vars.10666_0219_003.cdf')
      TES(1430)%FILENAME = TRIM('retv_vars.10666_0221_003.cdf')
      TES(1431)%FILENAME = TRIM('retv_vars.10666_0222_003.cdf')
      TES(1432)%FILENAME = TRIM('retv_vars.10666_0230_004.cdf')
      TES(1433)%FILENAME = TRIM('retv_vars.10666_0232_004.cdf')
      TES(1434)%FILENAME = TRIM('retv_vars.10666_0242_003.cdf')
      TES(1435)%FILENAME = TRIM('retv_vars.10666_0243_002.cdf')
      TES(1436)%FILENAME = TRIM('retv_vars.10666_0243_003.cdf')
      TES(1437)%FILENAME = TRIM('retv_vars.10666_0243_004.cdf')
      TES(1438)%FILENAME = TRIM('retv_vars.10666_0244_002.cdf')
      TES(1439)%FILENAME = TRIM('retv_vars.10666_0244_003.cdf')
      TES(1440)%FILENAME = TRIM('retv_vars.10666_0244_004.cdf')
      TES(1441)%FILENAME = TRIM('retv_vars.10666_0247_004.cdf')
      TES(1442)%FILENAME = TRIM('retv_vars.10666_0248_002.cdf')
      TES(1443)%FILENAME = TRIM('retv_vars.10666_0249_003.cdf')
      TES(1444)%FILENAME = TRIM('retv_vars.10666_0249_004.cdf')
      TES(1445)%FILENAME = TRIM('retv_vars.10666_0250_004.cdf')
      TES(1446)%FILENAME = TRIM('retv_vars.10666_0251_003.cdf')
      TES(1447)%FILENAME = TRIM('retv_vars.10666_0251_004.cdf')
      TES(1448)%FILENAME = TRIM('retv_vars.10666_0252_004.cdf')
      TES(1449)%FILENAME = TRIM('retv_vars.10666_0256_003.cdf')
      TES(1450)%FILENAME = TRIM('retv_vars.10666_0256_004.cdf')
      TES(1451)%FILENAME = TRIM('retv_vars.10666_0257_002.cdf')
      TES(1452)%FILENAME = TRIM('retv_vars.10666_0258_004.cdf')
      TES(1453)%FILENAME = TRIM('retv_vars.10666_0259_002.cdf')
      TES(1454)%FILENAME = TRIM('retv_vars.10666_0259_004.cdf')
      TES(1455)%FILENAME = TRIM('retv_vars.10666_0260_002.cdf')
      TES(1456)%FILENAME = TRIM('retv_vars.10666_0260_003.cdf')
      TES(1457)%FILENAME = TRIM('retv_vars.10666_0261_002.cdf')
      TES(1458)%FILENAME = TRIM('retv_vars.10666_0261_003.cdf')
      TES(1459)%FILENAME = TRIM('retv_vars.10666_0261_004.cdf')
      TES(1460)%FILENAME = TRIM('retv_vars.10666_0262_002.cdf')
      TES(1461)%FILENAME = TRIM('retv_vars.10666_0262_003.cdf')
      TES(1462)%FILENAME = TRIM('retv_vars.10666_0267_004.cdf')
      TES(1463)%FILENAME = TRIM('retv_vars.10666_0268_003.cdf')
      TES(1464)%FILENAME = TRIM('retv_vars.10666_0268_004.cdf')
      TES(1465)%FILENAME = TRIM('retv_vars.10666_0269_002.cdf')
      TES(1466)%FILENAME = TRIM('retv_vars.10666_0269_003.cdf')
      TES(1467)%FILENAME = TRIM('retv_vars.10666_0272_002.cdf')
      TES(1468)%FILENAME = TRIM('retv_vars.10666_0272_004.cdf')
      TES(1469)%FILENAME = TRIM('retv_vars.10666_0273_002.cdf')
      TES(1470)%FILENAME = TRIM('retv_vars.10666_0273_003.cdf')
      TES(1471)%FILENAME = TRIM('retv_vars.10666_0274_002.cdf')
      TES(1472)%FILENAME = TRIM('retv_vars.10666_0274_004.cdf')
      TES(1473)%FILENAME = TRIM('retv_vars.10666_0275_002.cdf')
      TES(1474)%FILENAME = TRIM('retv_vars.10666_0303_002.cdf')
      TES(1475)%FILENAME = TRIM('retv_vars.10666_0303_003.cdf')
      TES(1476)%FILENAME = TRIM('retv_vars.10666_0303_004.cdf')
      TES(1477)%FILENAME = TRIM('retv_vars.10666_0305_002.cdf')
      TES(1478)%FILENAME = TRIM('retv_vars.10666_0306_002.cdf')
      TES(1479)%FILENAME = TRIM('retv_vars.10666_0306_004.cdf')
      TES(1480)%FILENAME = TRIM('retv_vars.10666_0307_003.cdf')
      TES(1481)%FILENAME = TRIM('retv_vars.10666_0309_004.cdf')
      TES(1482)%FILENAME = TRIM('retv_vars.10666_0310_002.cdf')
      TES(1483)%FILENAME = TRIM('retv_vars.10666_0315_004.cdf')
      TES(1484)%FILENAME = TRIM('retv_vars.10666_0316_004.cdf')
      TES(1485)%FILENAME = TRIM('retv_vars.10666_0317_003.cdf')
      TES(1486)%FILENAME = TRIM('retv_vars.10666_0319_002.cdf')
      TES(1487)%FILENAME = TRIM('retv_vars.10666_0320_004.cdf')
      TES(1488)%FILENAME = TRIM('retv_vars.10666_0322_003.cdf')
      TES(1489)%FILENAME = TRIM('retv_vars.10666_0363_003.cdf')
      TES(1490)%FILENAME = TRIM('retv_vars.10666_0363_004.cdf')
      TES(1491)%FILENAME = TRIM('retv_vars.10666_0364_002.cdf')
      TES(1492)%FILENAME = TRIM('retv_vars.10666_0364_003.cdf')
      TES(1493)%FILENAME = TRIM('retv_vars.10666_0364_004.cdf')
      TES(1494)%FILENAME = TRIM('retv_vars.10666_0365_003.cdf')
      TES(1495)%FILENAME = TRIM('retv_vars.10666_0366_004.cdf')
      TES(1496)%FILENAME = TRIM('retv_vars.10666_0369_002.cdf')
      TES(1497)%FILENAME = TRIM('retv_vars.10666_0370_003.cdf')
      TES(1498)%FILENAME = TRIM('retv_vars.10666_0370_004.cdf')
      TES(1499)%FILENAME = TRIM('retv_vars.10666_0372_003.cdf')
      TES(1500)%FILENAME = TRIM('retv_vars.10666_0374_004.cdf')
      TES(1501)%FILENAME = TRIM('retv_vars.10666_0375_002.cdf')
      TES(1502)%FILENAME = TRIM('retv_vars.10666_0378_003.cdf')
      TES(1503)%FILENAME = TRIM('retv_vars.10666_0406_002.cdf')
      TES(1504)%FILENAME = TRIM('retv_vars.10666_0412_003.cdf')
      TES(1505)%FILENAME = TRIM('retv_vars.10666_0415_002.cdf')
      TES(1506)%FILENAME = TRIM('retv_vars.10666_0415_003.cdf')
      TES(1507)%FILENAME = TRIM('retv_vars.10666_0417_004.cdf')
      TES(1508)%FILENAME = TRIM('retv_vars.10666_0418_002.cdf')
      TES(1509)%FILENAME = TRIM('retv_vars.10666_0418_003.cdf')
      TES(1510)%FILENAME = TRIM('retv_vars.10666_0419_002.cdf')
      TES(1511)%FILENAME = TRIM('retv_vars.10666_0420_004.cdf')
      TES(1512)%FILENAME = TRIM('retv_vars.10666_0421_002.cdf')
      TES(1513)%FILENAME = TRIM('retv_vars.10666_0421_003.cdf')
      TES(1514)%FILENAME = TRIM('retv_vars.10666_0421_004.cdf')
      TES(1515)%FILENAME = TRIM('retv_vars.10666_0423_002.cdf')
      TES(1516)%FILENAME = TRIM('retv_vars.10666_0423_003.cdf')
      TES(1517)%FILENAME = TRIM('retv_vars.10666_0423_004.cdf')
      TES(1518)%FILENAME = TRIM('retv_vars.10666_0424_003.cdf')
      TES(1519)%FILENAME = TRIM('retv_vars.10666_0424_004.cdf')
      TES(1520)%FILENAME = TRIM('retv_vars.10666_0425_002.cdf')
      TES(1521)%FILENAME = TRIM('retv_vars.10666_0426_004.cdf')
      TES(1522)%FILENAME = TRIM('retv_vars.10666_0427_002.cdf')
      TES(1523)%FILENAME = TRIM('retv_vars.10666_0459_002.cdf')
      TES(1524)%FILENAME = TRIM('retv_vars.10666_0460_002.cdf')
      TES(1525)%FILENAME = TRIM('retv_vars.10666_0460_004.cdf')
      TES(1526)%FILENAME = TRIM('retv_vars.10666_0462_002.cdf')
      TES(1527)%FILENAME = TRIM('retv_vars.10666_0466_003.cdf')
      TES(1528)%FILENAME = TRIM('retv_vars.10666_0466_004.cdf')
      TES(1529)%FILENAME = TRIM('retv_vars.10666_0467_003.cdf')
      TES(1530)%FILENAME = TRIM('retv_vars.10666_0469_002.cdf')
      TES(1531)%FILENAME = TRIM('retv_vars.10666_0469_003.cdf')
      TES(1532)%FILENAME = TRIM('retv_vars.10666_0469_004.cdf')
      TES(1533)%FILENAME = TRIM('retv_vars.10666_0482_003.cdf')
      TES(1534)%FILENAME = TRIM('retv_vars.10666_0500_003.cdf')
      TES(1535)%FILENAME = TRIM('retv_vars.10666_0530_003.cdf')
      TES(1536)%FILENAME = TRIM('retv_vars.10666_0530_004.cdf')
      TES(1537)%FILENAME = TRIM('retv_vars.10666_0532_002.cdf')
      TES(1538)%FILENAME = TRIM('retv_vars.10666_0532_003.cdf')
      TES(1539)%FILENAME = TRIM('retv_vars.10666_0532_004.cdf')
      TES(1540)%FILENAME = TRIM('retv_vars.10666_0533_002.cdf')
      TES(1541)%FILENAME = TRIM('retv_vars.10666_0534_003.cdf')
      TES(1542)%FILENAME = TRIM('retv_vars.10666_0535_002.cdf')
      TES(1543)%FILENAME = TRIM('retv_vars.10666_0535_003.cdf')
      TES(1544)%FILENAME = TRIM('retv_vars.10666_0537_004.cdf')
      TES(1545)%FILENAME = TRIM('retv_vars.10666_0546_002.cdf')
      TES(1546)%FILENAME = TRIM('retv_vars.10666_0546_003.cdf')
      TES(1547)%FILENAME = TRIM('retv_vars.10666_0550_004.cdf')
      TES(1548)%FILENAME = TRIM('retv_vars.10666_0566_003.cdf')
      TES(1549)%FILENAME = TRIM('retv_vars.10666_0567_002.cdf')
      TES(1550)%FILENAME = TRIM('retv_vars.10666_0567_003.cdf')
      TES(1551)%FILENAME = TRIM('retv_vars.10666_0567_004.cdf')
      TES(1552)%FILENAME = TRIM('retv_vars.10666_0568_004.cdf')
      TES(1553)%FILENAME = TRIM('retv_vars.10666_0569_002.cdf')
      TES(1554)%FILENAME = TRIM('retv_vars.10666_0569_004.cdf')
      TES(1555)%FILENAME = TRIM('retv_vars.10666_0571_002.cdf')
      TES(1556)%FILENAME = TRIM('retv_vars.10666_0572_002.cdf')
      TES(1557)%FILENAME = TRIM('retv_vars.10666_0572_003.cdf')
      TES(1558)%FILENAME = TRIM('retv_vars.10666_0572_004.cdf')
      TES(1559)%FILENAME = TRIM('retv_vars.10666_0573_002.cdf')
      TES(1560)%FILENAME = TRIM('retv_vars.10666_0573_004.cdf')
      TES(1561)%FILENAME = TRIM('retv_vars.10666_0574_002.cdf')
      TES(1562)%FILENAME = TRIM('retv_vars.10666_0574_003.cdf')
      TES(1563)%FILENAME = TRIM('retv_vars.10666_0574_004.cdf')
      TES(1564)%FILENAME = TRIM('retv_vars.10666_0575_002.cdf')
      TES(1565)%FILENAME = TRIM('retv_vars.10666_0575_003.cdf')
      TES(1566)%FILENAME = TRIM('retv_vars.10666_0575_004.cdf')
      TES(1567)%FILENAME = TRIM('retv_vars.10666_0576_003.cdf')
      TES(1568)%FILENAME = TRIM('retv_vars.10666_0580_002.cdf')
      TES(1569)%FILENAME = TRIM('retv_vars.10666_0581_003.cdf')
      TES(1570)%FILENAME = TRIM('retv_vars.10666_0582_003.cdf')
      TES(1571)%FILENAME = TRIM('retv_vars.10666_0582_004.cdf')
      TES(1572)%FILENAME = TRIM('retv_vars.10666_0586_003.cdf')
      TES(1573)%FILENAME = TRIM('retv_vars.10666_0586_004.cdf')
      TES(1574)%FILENAME = TRIM('retv_vars.10666_0587_003.cdf')
      TES(1575)%FILENAME = TRIM('retv_vars.10666_0591_004.cdf')
      TES(1576)%FILENAME = TRIM('retv_vars.10666_0592_002.cdf')
      TES(1577)%FILENAME = TRIM('retv_vars.10666_0592_003.cdf')
      TES(1578)%FILENAME = TRIM('retv_vars.10666_0592_004.cdf')
      TES(1579)%FILENAME = TRIM('retv_vars.10666_0594_003.cdf')
      TES(1580)%FILENAME = TRIM('retv_vars.10666_0596_002.cdf')
      TES(1581)%FILENAME = TRIM('retv_vars.10666_0597_003.cdf')
      TES(1582)%FILENAME = TRIM('retv_vars.10666_0598_002.cdf')
      TES(1583)%FILENAME = TRIM('retv_vars.10666_0598_003.cdf')
      TES(1584)%FILENAME = TRIM('retv_vars.10666_0598_004.cdf')
      TES(1585)%FILENAME = TRIM('retv_vars.10666_0599_002.cdf')
      TES(1586)%FILENAME = TRIM('retv_vars.10666_0605_004.cdf')
      TES(1587)%FILENAME = TRIM('retv_vars.10666_0606_002.cdf')
      TES(1588)%FILENAME = TRIM('retv_vars.10666_0610_003.cdf')
      TES(1589)%FILENAME = TRIM('retv_vars.10666_0613_003.cdf')
      TES(1590)%FILENAME = TRIM('retv_vars.10666_0615_002.cdf')
      TES(1591)%FILENAME = TRIM('retv_vars.10666_0615_003.cdf')
      TES(1592)%FILENAME = TRIM('retv_vars.10666_0616_002.cdf')
      TES(1593)%FILENAME = TRIM('retv_vars.10666_0639_004.cdf')
      TES(1594)%FILENAME = TRIM('retv_vars.10666_0640_002.cdf')
      TES(1595)%FILENAME = TRIM('retv_vars.10666_0640_004.cdf')
      TES(1596)%FILENAME = TRIM('retv_vars.10666_0642_002.cdf')
      TES(1597)%FILENAME = TRIM('retv_vars.10666_0644_002.cdf')
      TES(1598)%FILENAME = TRIM('retv_vars.10666_0644_003.cdf')
      TES(1599)%FILENAME = TRIM('retv_vars.10666_0644_004.cdf')
      TES(1600)%FILENAME = TRIM('retv_vars.10666_0645_002.cdf')
      TES(1601)%FILENAME = TRIM('retv_vars.10666_0645_003.cdf')
      TES(1602)%FILENAME = TRIM('retv_vars.10666_0646_004.cdf')
      TES(1603)%FILENAME = TRIM('retv_vars.10666_0647_002.cdf')
      TES(1604)%FILENAME = TRIM('retv_vars.10666_0651_004.cdf')
      TES(1605)%FILENAME = TRIM('retv_vars.10666_0652_002.cdf')
      TES(1606)%FILENAME = TRIM('retv_vars.10666_0652_003.cdf')
      TES(1607)%FILENAME = TRIM('retv_vars.10666_0653_002.cdf')
      TES(1608)%FILENAME = TRIM('retv_vars.10666_0653_003.cdf')
      TES(1609)%FILENAME = TRIM('retv_vars.10666_0655_004.cdf')
      TES(1610)%FILENAME = TRIM('retv_vars.10666_0656_003.cdf')
      TES(1611)%FILENAME = TRIM('retv_vars.10666_0657_002.cdf')
      TES(1612)%FILENAME = TRIM('retv_vars.10666_0657_003.cdf')
      TES(1613)%FILENAME = TRIM('retv_vars.10666_0657_004.cdf')
      TES(1614)%FILENAME = TRIM('retv_vars.10666_0658_003.cdf')
      TES(1615)%FILENAME = TRIM('retv_vars.10666_0658_004.cdf')
      TES(1616)%FILENAME = TRIM('retv_vars.10666_0690_003.cdf')
      TES(1617)%FILENAME = TRIM('retv_vars.10666_0694_004.cdf')
      TES(1618)%FILENAME = TRIM('retv_vars.10666_0695_002.cdf')
      TES(1619)%FILENAME = TRIM('retv_vars.10666_0699_002.cdf')
      TES(1620)%FILENAME = TRIM('retv_vars.10666_0699_003.cdf')
      TES(1621)%FILENAME = TRIM('retv_vars.10666_0701_002.cdf')
      TES(1622)%FILENAME = TRIM('retv_vars.10666_0701_004.cdf')
      TES(1623)%FILENAME = TRIM('retv_vars.10666_0702_003.cdf')
      TES(1624)%FILENAME = TRIM('retv_vars.10666_0702_004.cdf')
      TES(1625)%FILENAME = TRIM('retv_vars.10666_0703_004.cdf')
      TES(1626)%FILENAME = TRIM('retv_vars.10666_0704_002.cdf')
      TES(1627)%FILENAME = TRIM('retv_vars.10666_0731_003.cdf')
      TES(1628)%FILENAME = TRIM('retv_vars.10666_0731_004.cdf')
      TES(1629)%FILENAME = TRIM('retv_vars.10666_0732_002.cdf')
      TES(1630)%FILENAME = TRIM('retv_vars.10666_0733_002.cdf')
      TES(1631)%FILENAME = TRIM('retv_vars.10666_0733_004.cdf')
      TES(1632)%FILENAME = TRIM('retv_vars.10666_0734_002.cdf')
      TES(1633)%FILENAME = TRIM('retv_vars.10666_0737_003.cdf')
      TES(1634)%FILENAME = TRIM('retv_vars.10666_0738_004.cdf')
      TES(1635)%FILENAME = TRIM('retv_vars.10666_0740_003.cdf')
      TES(1636)%FILENAME = TRIM('retv_vars.10666_0740_004.cdf')
      TES(1637)%FILENAME = TRIM('retv_vars.10666_0741_002.cdf')
      TES(1638)%FILENAME = TRIM('retv_vars.10666_0741_003.cdf')
      TES(1639)%FILENAME = TRIM('retv_vars.10666_0741_004.cdf')
      TES(1640)%FILENAME = TRIM('retv_vars.10666_0742_003.cdf')
      TES(1641)%FILENAME = TRIM('retv_vars.10666_0742_004.cdf')
      TES(1642)%FILENAME = TRIM('retv_vars.10666_0747_003.cdf')
      TES(1643)%FILENAME = TRIM('retv_vars.10674_0020_002.cdf')
      TES(1644)%FILENAME = TRIM('retv_vars.10674_0020_003.cdf')
      TES(1645)%FILENAME = TRIM('retv_vars.10674_0020_004.cdf')
      TES(1646)%FILENAME = TRIM('retv_vars.10674_0021_002.cdf')
      TES(1647)%FILENAME = TRIM('retv_vars.10674_0021_004.cdf')
      TES(1648)%FILENAME = TRIM('retv_vars.10674_0022_004.cdf')
      TES(1649)%FILENAME = TRIM('retv_vars.10674_0023_002.cdf')
      TES(1650)%FILENAME = TRIM('retv_vars.10674_0027_002.cdf')
      TES(1651)%FILENAME = TRIM('retv_vars.10674_0027_003.cdf')
      TES(1652)%FILENAME = TRIM('retv_vars.10674_0027_004.cdf')
      TES(1653)%FILENAME = TRIM('retv_vars.10674_0028_002.cdf')
      TES(1654)%FILENAME = TRIM('retv_vars.10674_0028_003.cdf')
      TES(1655)%FILENAME = TRIM('retv_vars.10674_0028_004.cdf')
      TES(1656)%FILENAME = TRIM('retv_vars.10674_0055_002.cdf')
      TES(1657)%FILENAME = TRIM('retv_vars.10674_0055_003.cdf')
      TES(1658)%FILENAME = TRIM('retv_vars.10674_0056_002.cdf')
      TES(1659)%FILENAME = TRIM('retv_vars.10674_0056_003.cdf')
      TES(1660)%FILENAME = TRIM('retv_vars.10674_0057_002.cdf')
      TES(1661)%FILENAME = TRIM('retv_vars.10674_0057_003.cdf')
      TES(1662)%FILENAME = TRIM('retv_vars.10674_0057_004.cdf')
      TES(1663)%FILENAME = TRIM('retv_vars.10674_0059_004.cdf')
      TES(1664)%FILENAME = TRIM('retv_vars.10674_0060_003.cdf')
      TES(1665)%FILENAME = TRIM('retv_vars.10674_0061_002.cdf')
      TES(1666)%FILENAME = TRIM('retv_vars.10674_0061_004.cdf')
      TES(1667)%FILENAME = TRIM('retv_vars.10674_0062_004.cdf')
      TES(1668)%FILENAME = TRIM('retv_vars.10674_0063_004.cdf')
      TES(1669)%FILENAME = TRIM('retv_vars.10674_0064_003.cdf')
      TES(1670)%FILENAME = TRIM('retv_vars.10674_0067_004.cdf')
      TES(1671)%FILENAME = TRIM('retv_vars.10674_0068_002.cdf')
      TES(1672)%FILENAME = TRIM('retv_vars.10674_0068_003.cdf')
      TES(1673)%FILENAME = TRIM('retv_vars.10674_0068_004.cdf')
      TES(1674)%FILENAME = TRIM('retv_vars.10674_0069_002.cdf')
      TES(1675)%FILENAME = TRIM('retv_vars.10674_0069_003.cdf')
      TES(1676)%FILENAME = TRIM('retv_vars.10674_0108_003.cdf')
      TES(1677)%FILENAME = TRIM('retv_vars.10674_0110_002.cdf')
      TES(1678)%FILENAME = TRIM('retv_vars.10674_0110_004.cdf')
      TES(1679)%FILENAME = TRIM('retv_vars.10674_0114_003.cdf')
      TES(1680)%FILENAME = TRIM('retv_vars.10674_0115_004.cdf')
      TES(1681)%FILENAME = TRIM('retv_vars.10674_0171_004.cdf')
      TES(1682)%FILENAME = TRIM('retv_vars.10674_0172_002.cdf')
      TES(1683)%FILENAME = TRIM('retv_vars.10674_0185_004.cdf')
      TES(1684)%FILENAME = TRIM('retv_vars.10674_0186_004.cdf')
      TES(1685)%FILENAME = TRIM('retv_vars.10674_0189_004.cdf')
      TES(1686)%FILENAME = TRIM('retv_vars.10674_0190_003.cdf')
      TES(1687)%FILENAME = TRIM('retv_vars.10674_0199_002.cdf')
      TES(1688)%FILENAME = TRIM('retv_vars.10674_0199_003.cdf')
      TES(1689)%FILENAME = TRIM('retv_vars.10674_0199_004.cdf')
      TES(1690)%FILENAME = TRIM('retv_vars.10674_0200_002.cdf')
      TES(1691)%FILENAME = TRIM('retv_vars.10674_0200_003.cdf')
      TES(1692)%FILENAME = TRIM('retv_vars.10674_0200_004.cdf')
      TES(1693)%FILENAME = TRIM('retv_vars.10674_0201_003.cdf')
      TES(1694)%FILENAME = TRIM('retv_vars.10674_0202_003.cdf')
      TES(1695)%FILENAME = TRIM('retv_vars.10674_0213_002.cdf')
      TES(1696)%FILENAME = TRIM('retv_vars.10674_0214_002.cdf')
      TES(1697)%FILENAME = TRIM('retv_vars.10674_0219_002.cdf')
      TES(1698)%FILENAME = TRIM('retv_vars.10674_0219_003.cdf')
      TES(1699)%FILENAME = TRIM('retv_vars.10674_0221_003.cdf')
      TES(1700)%FILENAME = TRIM('retv_vars.10674_0221_004.cdf')
      TES(1701)%FILENAME = TRIM('retv_vars.10674_0222_002.cdf')
      TES(1702)%FILENAME = TRIM('retv_vars.10674_0222_003.cdf')
      TES(1703)%FILENAME = TRIM('retv_vars.10674_0225_004.cdf')
      TES(1704)%FILENAME = TRIM('retv_vars.10674_0231_003.cdf')
      TES(1705)%FILENAME = TRIM('retv_vars.10674_0231_004.cdf')
      TES(1706)%FILENAME = TRIM('retv_vars.10674_0232_002.cdf')
      TES(1707)%FILENAME = TRIM('retv_vars.10674_0235_003.cdf')
      TES(1708)%FILENAME = TRIM('retv_vars.10674_0235_004.cdf')
      TES(1709)%FILENAME = TRIM('retv_vars.10674_0243_002.cdf')
      TES(1710)%FILENAME = TRIM('retv_vars.10674_0243_003.cdf')
      TES(1711)%FILENAME = TRIM('retv_vars.10674_0243_004.cdf')
      TES(1712)%FILENAME = TRIM('retv_vars.10674_0244_004.cdf')
      TES(1713)%FILENAME = TRIM('retv_vars.10674_0245_002.cdf')
      TES(1714)%FILENAME = TRIM('retv_vars.10674_0245_004.cdf')
      TES(1715)%FILENAME = TRIM('retv_vars.10674_0248_003.cdf')
      TES(1716)%FILENAME = TRIM('retv_vars.10674_0248_004.cdf')
      TES(1717)%FILENAME = TRIM('retv_vars.10674_0249_002.cdf')
      TES(1718)%FILENAME = TRIM('retv_vars.10674_0249_003.cdf')
      TES(1719)%FILENAME = TRIM('retv_vars.10674_0249_004.cdf')
      TES(1720)%FILENAME = TRIM('retv_vars.10674_0250_002.cdf')
      TES(1721)%FILENAME = TRIM('retv_vars.10674_0250_003.cdf')
      TES(1722)%FILENAME = TRIM('retv_vars.10674_0250_004.cdf')
      TES(1723)%FILENAME = TRIM('retv_vars.10674_0252_004.cdf')
      TES(1724)%FILENAME = TRIM('retv_vars.10674_0253_002.cdf')
      TES(1725)%FILENAME = TRIM('retv_vars.10674_0255_002.cdf')
      TES(1726)%FILENAME = TRIM('retv_vars.10674_0256_003.cdf')
      TES(1727)%FILENAME = TRIM('retv_vars.10674_0257_004.cdf')
      TES(1728)%FILENAME = TRIM('retv_vars.10674_0258_003.cdf')
      TES(1729)%FILENAME = TRIM('retv_vars.10674_0259_003.cdf')
      TES(1730)%FILENAME = TRIM('retv_vars.10674_0260_004.cdf')
      TES(1731)%FILENAME = TRIM('retv_vars.10674_0261_002.cdf')
      TES(1732)%FILENAME = TRIM('retv_vars.10674_0262_003.cdf')
      TES(1733)%FILENAME = TRIM('retv_vars.10674_0267_002.cdf')
      TES(1734)%FILENAME = TRIM('retv_vars.10674_0267_003.cdf')
      TES(1735)%FILENAME = TRIM('retv_vars.10674_0267_004.cdf')
      TES(1736)%FILENAME = TRIM('retv_vars.10674_0268_002.cdf')
      TES(1737)%FILENAME = TRIM('retv_vars.10674_0268_003.cdf')
      TES(1738)%FILENAME = TRIM('retv_vars.10674_0268_004.cdf')
      TES(1739)%FILENAME = TRIM('retv_vars.10674_0269_002.cdf')
      TES(1740)%FILENAME = TRIM('retv_vars.10674_0269_003.cdf')
      TES(1741)%FILENAME = TRIM('retv_vars.10674_0271_004.cdf')
      TES(1742)%FILENAME = TRIM('retv_vars.10674_0272_004.cdf')
      TES(1743)%FILENAME = TRIM('retv_vars.10674_0302_004.cdf')
      TES(1744)%FILENAME = TRIM('retv_vars.10674_0303_002.cdf')
      TES(1745)%FILENAME = TRIM('retv_vars.10674_0303_004.cdf')
      TES(1746)%FILENAME = TRIM('retv_vars.10674_0304_002.cdf')
      TES(1747)%FILENAME = TRIM('retv_vars.10674_0306_004.cdf')
      TES(1748)%FILENAME = TRIM('retv_vars.10674_0307_003.cdf')
      TES(1749)%FILENAME = TRIM('retv_vars.10674_0308_002.cdf')
      TES(1750)%FILENAME = TRIM('retv_vars.10674_0308_003.cdf')
      TES(1751)%FILENAME = TRIM('retv_vars.10674_0308_004.cdf')
      TES(1752)%FILENAME = TRIM('retv_vars.10674_0309_002.cdf')
      TES(1753)%FILENAME = TRIM('retv_vars.10674_0309_004.cdf')
      TES(1754)%FILENAME = TRIM('retv_vars.10674_0310_002.cdf')
      TES(1755)%FILENAME = TRIM('retv_vars.10674_0310_003.cdf')
      TES(1756)%FILENAME = TRIM('retv_vars.10674_0310_004.cdf')
      TES(1757)%FILENAME = TRIM('retv_vars.10674_0311_002.cdf')
      TES(1758)%FILENAME = TRIM('retv_vars.10674_0315_002.cdf')
      TES(1759)%FILENAME = TRIM('retv_vars.10674_0315_003.cdf')
      TES(1760)%FILENAME = TRIM('retv_vars.10674_0315_004.cdf')
      TES(1761)%FILENAME = TRIM('retv_vars.10674_0316_002.cdf')
      TES(1762)%FILENAME = TRIM('retv_vars.10674_0316_004.cdf')
      TES(1763)%FILENAME = TRIM('retv_vars.10674_0317_002.cdf')
      TES(1764)%FILENAME = TRIM('retv_vars.10674_0317_003.cdf')
      TES(1765)%FILENAME = TRIM('retv_vars.10674_0318_002.cdf')
      TES(1766)%FILENAME = TRIM('retv_vars.10674_0318_003.cdf')
      TES(1767)%FILENAME = TRIM('retv_vars.10674_0318_004.cdf')
      TES(1768)%FILENAME = TRIM('retv_vars.10674_0320_004.cdf')
      TES(1769)%FILENAME = TRIM('retv_vars.10674_0321_004.cdf')
      TES(1770)%FILENAME = TRIM('retv_vars.10674_0322_002.cdf')
      TES(1771)%FILENAME = TRIM('retv_vars.10674_0363_003.cdf')
      TES(1772)%FILENAME = TRIM('retv_vars.10674_0363_004.cdf')
      TES(1773)%FILENAME = TRIM('retv_vars.10674_0364_003.cdf')
      TES(1774)%FILENAME = TRIM('retv_vars.10674_0364_004.cdf')
      TES(1775)%FILENAME = TRIM('retv_vars.10674_0365_002.cdf')
      TES(1776)%FILENAME = TRIM('retv_vars.10674_0365_003.cdf')
      TES(1777)%FILENAME = TRIM('retv_vars.10674_0365_004.cdf')
      TES(1778)%FILENAME = TRIM('retv_vars.10674_0366_002.cdf')
      TES(1779)%FILENAME = TRIM('retv_vars.10674_0366_003.cdf')
      TES(1780)%FILENAME = TRIM('retv_vars.10674_0366_004.cdf')
      TES(1781)%FILENAME = TRIM('retv_vars.10674_0367_003.cdf')
      TES(1782)%FILENAME = TRIM('retv_vars.10674_0369_003.cdf')
      TES(1783)%FILENAME = TRIM('retv_vars.10674_0373_004.cdf')
      TES(1784)%FILENAME = TRIM('retv_vars.10674_0374_003.cdf')
      TES(1785)%FILENAME = TRIM('retv_vars.10674_0411_002.cdf')
      TES(1786)%FILENAME = TRIM('retv_vars.10674_0411_003.cdf')
      TES(1787)%FILENAME = TRIM('retv_vars.10674_0412_003.cdf')
      TES(1788)%FILENAME = TRIM('retv_vars.10674_0412_004.cdf')
      TES(1789)%FILENAME = TRIM('retv_vars.10674_0413_002.cdf')
      TES(1790)%FILENAME = TRIM('retv_vars.10674_0413_003.cdf')
      TES(1791)%FILENAME = TRIM('retv_vars.10674_0413_004.cdf')
      TES(1792)%FILENAME = TRIM('retv_vars.10674_0415_002.cdf')
      TES(1793)%FILENAME = TRIM('retv_vars.10674_0417_002.cdf')
      TES(1794)%FILENAME = TRIM('retv_vars.10674_0417_003.cdf')
      TES(1795)%FILENAME = TRIM('retv_vars.10674_0418_003.cdf')
      TES(1796)%FILENAME = TRIM('retv_vars.10674_0420_002.cdf')
      TES(1797)%FILENAME = TRIM('retv_vars.10674_0421_002.cdf')
      TES(1798)%FILENAME = TRIM('retv_vars.10674_0421_004.cdf')
      TES(1799)%FILENAME = TRIM('retv_vars.10674_0422_002.cdf')
      TES(1800)%FILENAME = TRIM('retv_vars.10674_0422_004.cdf')
      TES(1801)%FILENAME = TRIM('retv_vars.10674_0423_002.cdf')
      TES(1802)%FILENAME = TRIM('retv_vars.10674_0423_003.cdf')
      TES(1803)%FILENAME = TRIM('retv_vars.10674_0423_004.cdf')
      TES(1804)%FILENAME = TRIM('retv_vars.10674_0424_003.cdf')
      TES(1805)%FILENAME = TRIM('retv_vars.10674_0425_002.cdf')
      TES(1806)%FILENAME = TRIM('retv_vars.10674_0425_003.cdf')
      TES(1807)%FILENAME = TRIM('retv_vars.10674_0427_003.cdf')
      TES(1808)%FILENAME = TRIM('retv_vars.10674_0427_004.cdf')
      TES(1809)%FILENAME = TRIM('retv_vars.10674_0428_002.cdf')
      TES(1810)%FILENAME = TRIM('retv_vars.10674_0428_004.cdf')
      TES(1811)%FILENAME = TRIM('retv_vars.10674_0459_002.cdf')
      TES(1812)%FILENAME = TRIM('retv_vars.10674_0459_003.cdf')
      TES(1813)%FILENAME = TRIM('retv_vars.10674_0460_002.cdf')
      TES(1814)%FILENAME = TRIM('retv_vars.10674_0460_003.cdf')
      TES(1815)%FILENAME = TRIM('retv_vars.10674_0460_004.cdf')
      TES(1816)%FILENAME = TRIM('retv_vars.10674_0461_002.cdf')
      TES(1817)%FILENAME = TRIM('retv_vars.10674_0463_004.cdf')
      TES(1818)%FILENAME = TRIM('retv_vars.10674_0467_003.cdf')
      TES(1819)%FILENAME = TRIM('retv_vars.10674_0469_002.cdf')
      TES(1820)%FILENAME = TRIM('retv_vars.10674_0469_004.cdf')
      TES(1821)%FILENAME = TRIM('retv_vars.10674_0532_002.cdf')
      TES(1822)%FILENAME = TRIM('retv_vars.10674_0532_003.cdf')
      TES(1823)%FILENAME = TRIM('retv_vars.10674_0532_004.cdf')
      TES(1824)%FILENAME = TRIM('retv_vars.10674_0533_003.cdf')
      TES(1825)%FILENAME = TRIM('retv_vars.10674_0534_002.cdf')
      TES(1826)%FILENAME = TRIM('retv_vars.10674_0534_004.cdf')
      TES(1827)%FILENAME = TRIM('retv_vars.10674_0535_002.cdf')
      TES(1828)%FILENAME = TRIM('retv_vars.10674_0535_003.cdf')
      TES(1829)%FILENAME = TRIM('retv_vars.10674_0535_004.cdf')
      TES(1830)%FILENAME = TRIM('retv_vars.10674_0536_002.cdf')
      TES(1831)%FILENAME = TRIM('retv_vars.10674_0536_003.cdf')
      TES(1832)%FILENAME = TRIM('retv_vars.10674_0536_004.cdf')
      TES(1833)%FILENAME = TRIM('retv_vars.10674_0538_002.cdf')
      TES(1834)%FILENAME = TRIM('retv_vars.10674_0538_003.cdf')
      TES(1835)%FILENAME = TRIM('retv_vars.10674_0547_004.cdf')
      TES(1836)%FILENAME = TRIM('retv_vars.10674_0548_002.cdf')
      TES(1837)%FILENAME = TRIM('retv_vars.10674_0550_003.cdf')
      TES(1838)%FILENAME = TRIM('retv_vars.10674_0551_002.cdf')
      TES(1839)%FILENAME = TRIM('retv_vars.10674_0567_003.cdf')
      TES(1840)%FILENAME = TRIM('retv_vars.10674_0568_003.cdf')
      TES(1841)%FILENAME = TRIM('retv_vars.10674_0569_002.cdf')
      TES(1842)%FILENAME = TRIM('retv_vars.10674_0569_003.cdf')
      TES(1843)%FILENAME = TRIM('retv_vars.10674_0569_004.cdf')
      TES(1844)%FILENAME = TRIM('retv_vars.10674_0570_002.cdf')
      TES(1845)%FILENAME = TRIM('retv_vars.10674_0570_003.cdf')
      TES(1846)%FILENAME = TRIM('retv_vars.10674_0570_004.cdf')
      TES(1847)%FILENAME = TRIM('retv_vars.10674_0571_002.cdf')
      TES(1848)%FILENAME = TRIM('retv_vars.10674_0571_003.cdf')
      TES(1849)%FILENAME = TRIM('retv_vars.10674_0571_004.cdf')
      TES(1850)%FILENAME = TRIM('retv_vars.10674_0573_002.cdf')
      TES(1851)%FILENAME = TRIM('retv_vars.10674_0573_003.cdf')
      TES(1852)%FILENAME = TRIM('retv_vars.10674_0581_002.cdf')
      TES(1853)%FILENAME = TRIM('retv_vars.10674_0585_002.cdf')
      TES(1854)%FILENAME = TRIM('retv_vars.10674_0586_003.cdf')
      TES(1855)%FILENAME = TRIM('retv_vars.10674_0586_004.cdf')
      TES(1856)%FILENAME = TRIM('retv_vars.10674_0587_003.cdf')
      TES(1857)%FILENAME = TRIM('retv_vars.10674_0587_004.cdf')
      TES(1858)%FILENAME = TRIM('retv_vars.10674_0592_002.cdf')
      TES(1859)%FILENAME = TRIM('retv_vars.10674_0592_004.cdf')
      TES(1860)%FILENAME = TRIM('retv_vars.10674_0593_003.cdf')
      TES(1861)%FILENAME = TRIM('retv_vars.10674_0597_002.cdf')
      TES(1862)%FILENAME = TRIM('retv_vars.10674_0597_003.cdf')
      TES(1863)%FILENAME = TRIM('retv_vars.10674_0597_004.cdf')
      TES(1864)%FILENAME = TRIM('retv_vars.10674_0598_002.cdf')
      TES(1865)%FILENAME = TRIM('retv_vars.10674_0598_003.cdf')
      TES(1866)%FILENAME = TRIM('retv_vars.10674_0604_003.cdf')
      TES(1867)%FILENAME = TRIM('retv_vars.10674_0604_004.cdf')
      TES(1868)%FILENAME = TRIM('retv_vars.10674_0605_004.cdf')
      TES(1869)%FILENAME = TRIM('retv_vars.10674_0611_002.cdf')
      TES(1870)%FILENAME = TRIM('retv_vars.10674_0613_002.cdf')
      TES(1871)%FILENAME = TRIM('retv_vars.10674_0613_003.cdf')
      TES(1872)%FILENAME = TRIM('retv_vars.10674_0613_004.cdf')
      TES(1873)%FILENAME = TRIM('retv_vars.10674_0614_002.cdf')
      TES(1874)%FILENAME = TRIM('retv_vars.10674_0614_003.cdf')
      TES(1875)%FILENAME = TRIM('retv_vars.10674_0639_002.cdf')
      TES(1876)%FILENAME = TRIM('retv_vars.10674_0640_003.cdf')
      TES(1877)%FILENAME = TRIM('retv_vars.10674_0643_003.cdf')
      TES(1878)%FILENAME = TRIM('retv_vars.10674_0644_003.cdf')
      TES(1879)%FILENAME = TRIM('retv_vars.10674_0645_003.cdf')
      TES(1880)%FILENAME = TRIM('retv_vars.10674_0645_004.cdf')
      TES(1881)%FILENAME = TRIM('retv_vars.10674_0646_003.cdf')
      TES(1882)%FILENAME = TRIM('retv_vars.10674_0646_004.cdf')
      TES(1883)%FILENAME = TRIM('retv_vars.10674_0647_002.cdf')
      TES(1884)%FILENAME = TRIM('retv_vars.10674_0654_002.cdf')
      TES(1885)%FILENAME = TRIM('retv_vars.10674_0654_003.cdf')
      TES(1886)%FILENAME = TRIM('retv_vars.10674_0654_004.cdf')
      TES(1887)%FILENAME = TRIM('retv_vars.10674_0655_002.cdf')
      TES(1888)%FILENAME = TRIM('retv_vars.10674_0656_002.cdf')
      TES(1889)%FILENAME = TRIM('retv_vars.10674_0656_003.cdf')
      TES(1890)%FILENAME = TRIM('retv_vars.10674_0656_004.cdf')
      TES(1891)%FILENAME = TRIM('retv_vars.10674_0657_002.cdf')
      TES(1892)%FILENAME = TRIM('retv_vars.10674_0659_002.cdf')
      TES(1893)%FILENAME = TRIM('retv_vars.10674_0659_003.cdf')
      TES(1894)%FILENAME = TRIM('retv_vars.10674_0659_004.cdf')
      TES(1895)%FILENAME = TRIM('retv_vars.10674_0690_003.cdf')
      TES(1896)%FILENAME = TRIM('retv_vars.10674_0691_004.cdf')
      TES(1897)%FILENAME = TRIM('retv_vars.10674_0693_004.cdf')
      TES(1898)%FILENAME = TRIM('retv_vars.10674_0694_004.cdf')
      TES(1899)%FILENAME = TRIM('retv_vars.10674_0695_002.cdf')
      TES(1900)%FILENAME = TRIM('retv_vars.10674_0699_002.cdf')
      TES(1901)%FILENAME = TRIM('retv_vars.10674_0699_003.cdf')
      TES(1902)%FILENAME = TRIM('retv_vars.10674_0699_004.cdf')
      TES(1903)%FILENAME = TRIM('retv_vars.10674_0700_002.cdf')
      TES(1904)%FILENAME = TRIM('retv_vars.10674_0700_003.cdf')
      TES(1905)%FILENAME = TRIM('retv_vars.10674_0700_004.cdf')
      TES(1906)%FILENAME = TRIM('retv_vars.10674_0701_002.cdf')
      TES(1907)%FILENAME = TRIM('retv_vars.10674_0701_003.cdf')
      TES(1908)%FILENAME = TRIM('retv_vars.10674_0701_004.cdf')
      TES(1909)%FILENAME = TRIM('retv_vars.10674_0702_002.cdf')
      TES(1910)%FILENAME = TRIM('retv_vars.10674_0702_003.cdf')
      TES(1911)%FILENAME = TRIM('retv_vars.10674_0702_004.cdf')
      TES(1912)%FILENAME = TRIM('retv_vars.10674_0703_002.cdf')
      TES(1913)%FILENAME = TRIM('retv_vars.10674_0703_003.cdf')
      TES(1914)%FILENAME = TRIM('retv_vars.10674_0704_003.cdf')
      TES(1915)%FILENAME = TRIM('retv_vars.10674_0727_002.cdf')
      TES(1916)%FILENAME = TRIM('retv_vars.10674_0731_003.cdf')
      TES(1917)%FILENAME = TRIM('retv_vars.10674_0731_004.cdf')
      TES(1918)%FILENAME = TRIM('retv_vars.10674_0732_002.cdf')
      TES(1919)%FILENAME = TRIM('retv_vars.10674_0734_002.cdf')
      TES(1920)%FILENAME = TRIM('retv_vars.10674_0734_003.cdf')
      TES(1921)%FILENAME = TRIM('retv_vars.10674_0740_004.cdf')
      TES(1922)%FILENAME = TRIM('retv_vars.10674_0741_003.cdf')
      TES(1923)%FILENAME = TRIM('retv_vars.10674_0742_002.cdf')
      TES(1924)%FILENAME = TRIM('retv_vars.10674_0743_002.cdf')
      TES(1925)%FILENAME = TRIM('retv_vars.10674_0743_003.cdf')
      TES(1926)%FILENAME = TRIM('retv_vars.10679_0019_002.cdf')
      TES(1927)%FILENAME = TRIM('retv_vars.10679_0019_004.cdf')
      TES(1928)%FILENAME = TRIM('retv_vars.10679_0021_003.cdf')
      TES(1929)%FILENAME = TRIM('retv_vars.10679_0021_004.cdf')
      TES(1930)%FILENAME = TRIM('retv_vars.10679_0022_002.cdf')
      TES(1931)%FILENAME = TRIM('retv_vars.10679_0022_003.cdf')
      TES(1932)%FILENAME = TRIM('retv_vars.10679_0022_004.cdf')
      TES(1933)%FILENAME = TRIM('retv_vars.10679_0023_002.cdf')
      TES(1934)%FILENAME = TRIM('retv_vars.10679_0027_003.cdf')
      TES(1935)%FILENAME = TRIM('retv_vars.10679_0027_004.cdf')
      TES(1936)%FILENAME = TRIM('retv_vars.10679_0028_002.cdf')
      TES(1937)%FILENAME = TRIM('retv_vars.10679_0028_003.cdf')
      TES(1938)%FILENAME = TRIM('retv_vars.10679_0028_004.cdf')
      TES(1939)%FILENAME = TRIM('retv_vars.10679_0029_002.cdf')
      TES(1940)%FILENAME = TRIM('retv_vars.10679_0054_003.cdf')
      TES(1941)%FILENAME = TRIM('retv_vars.10679_0054_004.cdf')
      TES(1942)%FILENAME = TRIM('retv_vars.10679_0055_002.cdf')
      TES(1943)%FILENAME = TRIM('retv_vars.10679_0056_002.cdf')
      TES(1944)%FILENAME = TRIM('retv_vars.10679_0056_003.cdf')
      TES(1945)%FILENAME = TRIM('retv_vars.10679_0057_002.cdf')
      TES(1946)%FILENAME = TRIM('retv_vars.10679_0058_003.cdf')
      TES(1947)%FILENAME = TRIM('retv_vars.10679_0058_004.cdf')
      TES(1948)%FILENAME = TRIM('retv_vars.10679_0059_002.cdf')
      TES(1949)%FILENAME = TRIM('retv_vars.10679_0060_003.cdf')
      TES(1950)%FILENAME = TRIM('retv_vars.10679_0061_002.cdf')
      TES(1951)%FILENAME = TRIM('retv_vars.10679_0062_003.cdf')
      TES(1952)%FILENAME = TRIM('retv_vars.10679_0066_003.cdf')
      TES(1953)%FILENAME = TRIM('retv_vars.10679_0066_004.cdf')
      TES(1954)%FILENAME = TRIM('retv_vars.10679_0067_002.cdf')
      TES(1955)%FILENAME = TRIM('retv_vars.10679_0067_003.cdf')
      TES(1956)%FILENAME = TRIM('retv_vars.10679_0067_004.cdf')
      TES(1957)%FILENAME = TRIM('retv_vars.10679_0068_003.cdf')
      TES(1958)%FILENAME = TRIM('retv_vars.10679_0069_002.cdf')
      TES(1959)%FILENAME = TRIM('retv_vars.10679_0069_004.cdf')
      TES(1960)%FILENAME = TRIM('retv_vars.10679_0070_003.cdf')
      TES(1961)%FILENAME = TRIM('retv_vars.10679_0070_004.cdf')
      TES(1962)%FILENAME = TRIM('retv_vars.10679_0071_002.cdf')
      TES(1963)%FILENAME = TRIM('retv_vars.10679_0108_004.cdf')
      TES(1964)%FILENAME = TRIM('retv_vars.10679_0109_002.cdf')
      TES(1965)%FILENAME = TRIM('retv_vars.10679_0109_004.cdf')
      TES(1966)%FILENAME = TRIM('retv_vars.10679_0112_003.cdf')
      TES(1967)%FILENAME = TRIM('retv_vars.10679_0115_003.cdf')
      TES(1968)%FILENAME = TRIM('retv_vars.10679_0115_004.cdf')
      TES(1969)%FILENAME = TRIM('retv_vars.10679_0117_004.cdf')
      TES(1970)%FILENAME = TRIM('retv_vars.10679_0118_002.cdf')
      TES(1971)%FILENAME = TRIM('retv_vars.10679_0187_004.cdf')
      TES(1972)%FILENAME = TRIM('retv_vars.10679_0190_002.cdf')
      TES(1973)%FILENAME = TRIM('retv_vars.10679_0190_003.cdf')
      TES(1974)%FILENAME = TRIM('retv_vars.10679_0213_003.cdf')
      TES(1975)%FILENAME = TRIM('retv_vars.10679_0221_002.cdf')
      TES(1976)%FILENAME = TRIM('retv_vars.10679_0221_004.cdf')
      TES(1977)%FILENAME = TRIM('retv_vars.10679_0224_003.cdf')
      TES(1978)%FILENAME = TRIM('retv_vars.10679_0235_003.cdf')
      TES(1979)%FILENAME = TRIM('retv_vars.10679_0236_004.cdf')
      TES(1980)%FILENAME = TRIM('retv_vars.10679_0237_002.cdf')
      TES(1981)%FILENAME = TRIM('retv_vars.10679_0237_004.cdf')
      TES(1982)%FILENAME = TRIM('retv_vars.10679_0238_002.cdf')
      TES(1983)%FILENAME = TRIM('retv_vars.10679_0244_004.cdf')
      TES(1984)%FILENAME = TRIM('retv_vars.10679_0245_002.cdf')
      TES(1985)%FILENAME = TRIM('retv_vars.10679_0246_003.cdf')
      TES(1986)%FILENAME = TRIM('retv_vars.10679_0247_002.cdf')
      TES(1987)%FILENAME = TRIM('retv_vars.10679_0247_003.cdf')
      TES(1988)%FILENAME = TRIM('retv_vars.10679_0247_004.cdf')
      TES(1989)%FILENAME = TRIM('retv_vars.10679_0248_003.cdf')
      TES(1990)%FILENAME = TRIM('retv_vars.10679_0249_002.cdf')
      TES(1991)%FILENAME = TRIM('retv_vars.10679_0249_003.cdf')
      TES(1992)%FILENAME = TRIM('retv_vars.10679_0249_004.cdf')
      TES(1993)%FILENAME = TRIM('retv_vars.10679_0250_002.cdf')
      TES(1994)%FILENAME = TRIM('retv_vars.10679_0250_003.cdf')
      TES(1995)%FILENAME = TRIM('retv_vars.10679_0250_004.cdf')
      TES(1996)%FILENAME = TRIM('retv_vars.10679_0251_003.cdf')
      TES(1997)%FILENAME = TRIM('retv_vars.10679_0252_003.cdf')
      TES(1998)%FILENAME = TRIM('retv_vars.10679_0254_004.cdf')
      TES(1999)%FILENAME = TRIM('retv_vars.10679_0261_003.cdf')
      TES(2000)%FILENAME = TRIM('retv_vars.10679_0261_004.cdf')
      TES(2001)%FILENAME = TRIM('retv_vars.10679_0262_002.cdf')
      TES(2002)%FILENAME = TRIM('retv_vars.10679_0267_002.cdf')
      TES(2003)%FILENAME = TRIM('retv_vars.10679_0267_003.cdf')
      TES(2004)%FILENAME = TRIM('retv_vars.10679_0267_004.cdf')
      TES(2005)%FILENAME = TRIM('retv_vars.10679_0268_002.cdf')
      TES(2006)%FILENAME = TRIM('retv_vars.10679_0268_003.cdf')
      TES(2007)%FILENAME = TRIM('retv_vars.10679_0268_004.cdf')
      TES(2008)%FILENAME = TRIM('retv_vars.10679_0270_004.cdf')
      TES(2009)%FILENAME = TRIM('retv_vars.10679_0273_003.cdf')
      TES(2010)%FILENAME = TRIM('retv_vars.10679_0273_004.cdf')
      TES(2011)%FILENAME = TRIM('retv_vars.10679_0274_002.cdf')
      TES(2012)%FILENAME = TRIM('retv_vars.10679_0274_003.cdf')
      TES(2013)%FILENAME = TRIM('retv_vars.10679_0278_004.cdf')
      TES(2014)%FILENAME = TRIM('retv_vars.10679_0279_002.cdf')
      TES(2015)%FILENAME = TRIM('retv_vars.10679_0279_003.cdf')
      TES(2016)%FILENAME = TRIM('retv_vars.10679_0302_003.cdf')
      TES(2017)%FILENAME = TRIM('retv_vars.10679_0302_004.cdf')
      TES(2018)%FILENAME = TRIM('retv_vars.10679_0303_003.cdf')
      TES(2019)%FILENAME = TRIM('retv_vars.10679_0304_004.cdf')
      TES(2020)%FILENAME = TRIM('retv_vars.10679_0305_004.cdf')
      TES(2021)%FILENAME = TRIM('retv_vars.10679_0306_003.cdf')
      TES(2022)%FILENAME = TRIM('retv_vars.10679_0306_004.cdf')
      TES(2023)%FILENAME = TRIM('retv_vars.10679_0307_002.cdf')
      TES(2024)%FILENAME = TRIM('retv_vars.10679_0308_002.cdf')
      TES(2025)%FILENAME = TRIM('retv_vars.10679_0309_002.cdf')
      TES(2026)%FILENAME = TRIM('retv_vars.10679_0309_004.cdf')
      TES(2027)%FILENAME = TRIM('retv_vars.10679_0310_002.cdf')
      TES(2028)%FILENAME = TRIM('retv_vars.10679_0310_003.cdf')
      TES(2029)%FILENAME = TRIM('retv_vars.10679_0310_004.cdf')
      TES(2030)%FILENAME = TRIM('retv_vars.10679_0311_002.cdf')
      TES(2031)%FILENAME = TRIM('retv_vars.10679_0315_002.cdf')
      TES(2032)%FILENAME = TRIM('retv_vars.10679_0315_004.cdf')
      TES(2033)%FILENAME = TRIM('retv_vars.10679_0316_003.cdf')
      TES(2034)%FILENAME = TRIM('retv_vars.10679_0316_004.cdf')
      TES(2035)%FILENAME = TRIM('retv_vars.10679_0317_004.cdf')
      TES(2036)%FILENAME = TRIM('retv_vars.10679_0318_002.cdf')
      TES(2037)%FILENAME = TRIM('retv_vars.10679_0318_003.cdf')
      TES(2038)%FILENAME = TRIM('retv_vars.10679_0318_004.cdf')
      TES(2039)%FILENAME = TRIM('retv_vars.10679_0319_004.cdf')
      TES(2040)%FILENAME = TRIM('retv_vars.10679_0321_003.cdf')
      TES(2041)%FILENAME = TRIM('retv_vars.10679_0321_004.cdf')
      TES(2042)%FILENAME = TRIM('retv_vars.10679_0322_002.cdf')
      TES(2043)%FILENAME = TRIM('retv_vars.10679_0323_002.cdf')
      TES(2044)%FILENAME = TRIM('retv_vars.10679_0324_002.cdf')
      TES(2045)%FILENAME = TRIM('retv_vars.10679_0358_003.cdf')
      TES(2046)%FILENAME = TRIM('retv_vars.10679_0363_003.cdf')
      TES(2047)%FILENAME = TRIM('retv_vars.10679_0363_004.cdf')
      TES(2048)%FILENAME = TRIM('retv_vars.10679_0364_002.cdf')
      TES(2049)%FILENAME = TRIM('retv_vars.10679_0364_003.cdf')
      TES(2050)%FILENAME = TRIM('retv_vars.10679_0364_004.cdf')
      TES(2051)%FILENAME = TRIM('retv_vars.10679_0365_002.cdf')
      TES(2052)%FILENAME = TRIM('retv_vars.10679_0365_003.cdf')
      TES(2053)%FILENAME = TRIM('retv_vars.10679_0366_003.cdf')
      TES(2054)%FILENAME = TRIM('retv_vars.10679_0367_004.cdf')
      TES(2055)%FILENAME = TRIM('retv_vars.10679_0369_004.cdf')
      TES(2056)%FILENAME = TRIM('retv_vars.10679_0370_003.cdf')
      TES(2057)%FILENAME = TRIM('retv_vars.10679_0370_004.cdf')
      TES(2058)%FILENAME = TRIM('retv_vars.10679_0378_004.cdf')
      TES(2059)%FILENAME = TRIM('retv_vars.10679_0379_002.cdf')
      TES(2060)%FILENAME = TRIM('retv_vars.10679_0379_003.cdf')
      TES(2061)%FILENAME = TRIM('retv_vars.10679_0379_004.cdf')
      TES(2062)%FILENAME = TRIM('retv_vars.10679_0380_002.cdf')
      TES(2063)%FILENAME = TRIM('retv_vars.10679_0380_003.cdf')
      TES(2064)%FILENAME = TRIM('retv_vars.10679_0406_003.cdf')
      TES(2065)%FILENAME = TRIM('retv_vars.10679_0411_002.cdf')
      TES(2066)%FILENAME = TRIM('retv_vars.10679_0411_003.cdf')
      TES(2067)%FILENAME = TRIM('retv_vars.10679_0412_002.cdf')
      TES(2068)%FILENAME = TRIM('retv_vars.10679_0412_003.cdf')
      TES(2069)%FILENAME = TRIM('retv_vars.10679_0415_002.cdf')
      TES(2070)%FILENAME = TRIM('retv_vars.10679_0415_003.cdf')
      TES(2071)%FILENAME = TRIM('retv_vars.10679_0416_004.cdf')
      TES(2072)%FILENAME = TRIM('retv_vars.10679_0417_004.cdf')
      TES(2073)%FILENAME = TRIM('retv_vars.10679_0418_003.cdf')
      TES(2074)%FILENAME = TRIM('retv_vars.10679_0418_004.cdf')
      TES(2075)%FILENAME = TRIM('retv_vars.10679_0419_002.cdf')
      TES(2076)%FILENAME = TRIM('retv_vars.10679_0419_003.cdf')
      TES(2077)%FILENAME = TRIM('retv_vars.10679_0419_004.cdf')
      TES(2078)%FILENAME = TRIM('retv_vars.10679_0422_002.cdf')
      TES(2079)%FILENAME = TRIM('retv_vars.10679_0422_003.cdf')
      TES(2080)%FILENAME = TRIM('retv_vars.10679_0423_002.cdf')
      TES(2081)%FILENAME = TRIM('retv_vars.10679_0423_003.cdf')
      TES(2082)%FILENAME = TRIM('retv_vars.10679_0423_004.cdf')
      TES(2083)%FILENAME = TRIM('retv_vars.10679_0424_002.cdf')
      TES(2084)%FILENAME = TRIM('retv_vars.10679_0424_004.cdf')
      TES(2085)%FILENAME = TRIM('retv_vars.10679_0425_004.cdf')
      TES(2086)%FILENAME = TRIM('retv_vars.10679_0427_002.cdf')
      TES(2087)%FILENAME = TRIM('retv_vars.10679_0428_004.cdf')
      TES(2088)%FILENAME = TRIM('retv_vars.10679_0459_002.cdf')
      TES(2089)%FILENAME = TRIM('retv_vars.10679_0461_002.cdf')
      TES(2090)%FILENAME = TRIM('retv_vars.10679_0462_003.cdf')
      TES(2091)%FILENAME = TRIM('retv_vars.10679_0463_004.cdf')
      TES(2092)%FILENAME = TRIM('retv_vars.10679_0464_003.cdf')
      TES(2093)%FILENAME = TRIM('retv_vars.10679_0467_003.cdf')
      TES(2094)%FILENAME = TRIM('retv_vars.10679_0468_004.cdf')
      TES(2095)%FILENAME = TRIM('retv_vars.10679_0493_003.cdf')
      TES(2096)%FILENAME = TRIM('retv_vars.10679_0532_003.cdf')
      TES(2097)%FILENAME = TRIM('retv_vars.10679_0532_004.cdf')
      TES(2098)%FILENAME = TRIM('retv_vars.10679_0533_002.cdf')
      TES(2099)%FILENAME = TRIM('retv_vars.10679_0533_003.cdf')
      TES(2100)%FILENAME = TRIM('retv_vars.10679_0534_004.cdf')
      TES(2101)%FILENAME = TRIM('retv_vars.10679_0548_004.cdf')
      TES(2102)%FILENAME = TRIM('retv_vars.10679_0549_002.cdf')
      TES(2103)%FILENAME = TRIM('retv_vars.10679_0567_003.cdf')
      TES(2104)%FILENAME = TRIM('retv_vars.10679_0568_004.cdf')
      TES(2105)%FILENAME = TRIM('retv_vars.10679_0569_002.cdf')
      TES(2106)%FILENAME = TRIM('retv_vars.10679_0569_003.cdf')
      TES(2107)%FILENAME = TRIM('retv_vars.10679_0569_004.cdf')
      TES(2108)%FILENAME = TRIM('retv_vars.10679_0570_002.cdf')
      TES(2109)%FILENAME = TRIM('retv_vars.10679_0571_002.cdf')
      TES(2110)%FILENAME = TRIM('retv_vars.10679_0571_004.cdf')
      TES(2111)%FILENAME = TRIM('retv_vars.10679_0572_002.cdf')
      TES(2112)%FILENAME = TRIM('retv_vars.10679_0572_004.cdf')
      TES(2113)%FILENAME = TRIM('retv_vars.10679_0573_002.cdf')
      TES(2114)%FILENAME = TRIM('retv_vars.10679_0573_004.cdf')
      TES(2115)%FILENAME = TRIM('retv_vars.10679_0580_003.cdf')
      TES(2116)%FILENAME = TRIM('retv_vars.10679_0580_004.cdf')
      TES(2117)%FILENAME = TRIM('retv_vars.10679_0583_003.cdf')
      TES(2118)%FILENAME = TRIM('retv_vars.10679_0585_004.cdf')
      TES(2119)%FILENAME = TRIM('retv_vars.10679_0586_003.cdf')
      TES(2120)%FILENAME = TRIM('retv_vars.10679_0587_002.cdf')
      TES(2121)%FILENAME = TRIM('retv_vars.10679_0587_004.cdf')
      TES(2122)%FILENAME = TRIM('retv_vars.10679_0588_002.cdf')
      TES(2123)%FILENAME = TRIM('retv_vars.10679_0591_004.cdf')
      TES(2124)%FILENAME = TRIM('retv_vars.10679_0592_003.cdf')
      TES(2125)%FILENAME = TRIM('retv_vars.10679_0594_002.cdf')
      TES(2126)%FILENAME = TRIM('retv_vars.10679_0594_003.cdf')
      TES(2127)%FILENAME = TRIM('retv_vars.10679_0596_003.cdf')
      TES(2128)%FILENAME = TRIM('retv_vars.10679_0596_004.cdf')
      TES(2129)%FILENAME = TRIM('retv_vars.10679_0614_002.cdf')
      TES(2130)%FILENAME = TRIM('retv_vars.10679_0615_002.cdf')
      TES(2131)%FILENAME = TRIM('retv_vars.10679_0615_003.cdf')
      TES(2132)%FILENAME = TRIM('retv_vars.10679_0616_002.cdf')
      TES(2133)%FILENAME = TRIM('retv_vars.10679_0616_003.cdf')
      TES(2134)%FILENAME = TRIM('retv_vars.10679_0617_002.cdf')
      TES(2135)%FILENAME = TRIM('retv_vars.10679_0617_003.cdf')
      TES(2136)%FILENAME = TRIM('retv_vars.10679_0639_003.cdf')
      TES(2137)%FILENAME = TRIM('retv_vars.10679_0639_004.cdf')
      TES(2138)%FILENAME = TRIM('retv_vars.10679_0640_002.cdf')
      TES(2139)%FILENAME = TRIM('retv_vars.10679_0643_004.cdf')
      TES(2140)%FILENAME = TRIM('retv_vars.10679_0644_003.cdf')
      TES(2141)%FILENAME = TRIM('retv_vars.10679_0644_004.cdf')
      TES(2142)%FILENAME = TRIM('retv_vars.10679_0645_002.cdf')
      TES(2143)%FILENAME = TRIM('retv_vars.10679_0645_003.cdf')
      TES(2144)%FILENAME = TRIM('retv_vars.10679_0646_002.cdf')
      TES(2145)%FILENAME = TRIM('retv_vars.10679_0646_004.cdf')
      TES(2146)%FILENAME = TRIM('retv_vars.10679_0647_002.cdf')
      TES(2147)%FILENAME = TRIM('retv_vars.10679_0651_002.cdf')
      TES(2148)%FILENAME = TRIM('retv_vars.10679_0652_002.cdf')
      TES(2149)%FILENAME = TRIM('retv_vars.10679_0652_003.cdf')
      TES(2150)%FILENAME = TRIM('retv_vars.10679_0654_003.cdf')
      TES(2151)%FILENAME = TRIM('retv_vars.10679_0655_004.cdf')
      TES(2152)%FILENAME = TRIM('retv_vars.10679_0656_003.cdf')
      TES(2153)%FILENAME = TRIM('retv_vars.10679_0689_003.cdf')
      TES(2154)%FILENAME = TRIM('retv_vars.10679_0691_004.cdf')
      TES(2155)%FILENAME = TRIM('retv_vars.10679_0692_002.cdf')
      TES(2156)%FILENAME = TRIM('retv_vars.10679_0692_004.cdf')
      TES(2157)%FILENAME = TRIM('retv_vars.10679_0693_004.cdf')
      TES(2158)%FILENAME = TRIM('retv_vars.10679_0694_002.cdf')
      TES(2159)%FILENAME = TRIM('retv_vars.10679_0694_003.cdf')
      TES(2160)%FILENAME = TRIM('retv_vars.10679_0694_004.cdf')
      TES(2161)%FILENAME = TRIM('retv_vars.10679_0695_002.cdf')
      TES(2162)%FILENAME = TRIM('retv_vars.10679_0699_002.cdf')
      TES(2163)%FILENAME = TRIM('retv_vars.10679_0699_003.cdf')
      TES(2164)%FILENAME = TRIM('retv_vars.10679_0699_004.cdf')
      TES(2165)%FILENAME = TRIM('retv_vars.10679_0700_002.cdf')
      TES(2166)%FILENAME = TRIM('retv_vars.10679_0700_003.cdf')
      TES(2167)%FILENAME = TRIM('retv_vars.10679_0701_003.cdf')
      TES(2168)%FILENAME = TRIM('retv_vars.10679_0701_004.cdf')
      TES(2169)%FILENAME = TRIM('retv_vars.10679_0702_002.cdf')
      TES(2170)%FILENAME = TRIM('retv_vars.10679_0702_003.cdf')
      TES(2171)%FILENAME = TRIM('retv_vars.10679_0703_003.cdf')
      TES(2172)%FILENAME = TRIM('retv_vars.10679_0703_004.cdf')
      TES(2173)%FILENAME = TRIM('retv_vars.10679_0727_003.cdf')
      TES(2174)%FILENAME = TRIM('retv_vars.10679_0728_002.cdf')
      TES(2175)%FILENAME = TRIM('retv_vars.10679_0731_004.cdf')
      TES(2176)%FILENAME = TRIM('retv_vars.10679_0732_004.cdf')
      TES(2177)%FILENAME = TRIM('retv_vars.10679_0738_003.cdf')
      TES(2178)%FILENAME = TRIM('retv_vars.10679_0738_004.cdf')
      TES(2179)%FILENAME = TRIM('retv_vars.10679_0739_002.cdf')
      TES(2180)%FILENAME = TRIM('retv_vars.10679_0740_004.cdf')
      TES(2181)%FILENAME = TRIM('retv_vars.10679_0741_002.cdf')
      TES(2182)%FILENAME = TRIM('retv_vars.10679_0741_003.cdf')
      TES(2183)%FILENAME = TRIM('retv_vars.10679_0742_002.cdf')
      TES(2184)%FILENAME = TRIM('retv_vars.10679_0742_003.cdf')
      TES(2185)%FILENAME = TRIM('retv_vars.10679_0742_004.cdf')
      TES(2186)%FILENAME = TRIM('retv_vars.10679_0743_002.cdf')
      TES(2187)%FILENAME = TRIM('retv_vars.10679_0747_004.cdf')
      TES(2188)%FILENAME = TRIM('retv_vars.10679_0762_004.cdf')
      TES(2189)%FILENAME = TRIM('retv_vars.10684_0018_003.cdf')
      TES(2190)%FILENAME = TRIM('retv_vars.10684_0019_003.cdf')
      TES(2191)%FILENAME = TRIM('retv_vars.10684_0021_002.cdf')
      TES(2192)%FILENAME = TRIM('retv_vars.10684_0021_003.cdf')
      TES(2193)%FILENAME = TRIM('retv_vars.10684_0021_004.cdf')
      TES(2194)%FILENAME = TRIM('retv_vars.10684_0022_002.cdf')
      TES(2195)%FILENAME = TRIM('retv_vars.10684_0022_003.cdf')
      TES(2196)%FILENAME = TRIM('retv_vars.10684_0022_004.cdf')
      TES(2197)%FILENAME = TRIM('retv_vars.10684_0023_002.cdf')
      TES(2198)%FILENAME = TRIM('retv_vars.10684_0027_002.cdf')
      TES(2199)%FILENAME = TRIM('retv_vars.10684_0027_003.cdf')
      TES(2200)%FILENAME = TRIM('retv_vars.10684_0028_003.cdf')
      TES(2201)%FILENAME = TRIM('retv_vars.10684_0028_004.cdf')
      TES(2202)%FILENAME = TRIM('retv_vars.10684_0029_002.cdf')
      TES(2203)%FILENAME = TRIM('retv_vars.10684_0029_004.cdf')
      TES(2204)%FILENAME = TRIM('retv_vars.10684_0030_004.cdf')
      TES(2205)%FILENAME = TRIM('retv_vars.10684_0055_004.cdf')
      TES(2206)%FILENAME = TRIM('retv_vars.10684_0056_002.cdf')
      TES(2207)%FILENAME = TRIM('retv_vars.10684_0057_003.cdf')
      TES(2208)%FILENAME = TRIM('retv_vars.10684_0057_004.cdf')
      TES(2209)%FILENAME = TRIM('retv_vars.10684_0058_002.cdf')
      TES(2210)%FILENAME = TRIM('retv_vars.10684_0058_003.cdf')
      TES(2211)%FILENAME = TRIM('retv_vars.10684_0058_004.cdf')
      TES(2212)%FILENAME = TRIM('retv_vars.10684_0059_002.cdf')
      TES(2213)%FILENAME = TRIM('retv_vars.10684_0059_004.cdf')
      TES(2214)%FILENAME = TRIM('retv_vars.10684_0060_002.cdf')
      TES(2215)%FILENAME = TRIM('retv_vars.10684_0061_003.cdf')
      TES(2216)%FILENAME = TRIM('retv_vars.10684_0062_002.cdf')
      TES(2217)%FILENAME = TRIM('retv_vars.10684_0063_004.cdf')
      TES(2218)%FILENAME = TRIM('retv_vars.10684_0064_002.cdf')
      TES(2219)%FILENAME = TRIM('retv_vars.10684_0067_002.cdf')
      TES(2220)%FILENAME = TRIM('retv_vars.10684_0067_003.cdf')
      TES(2221)%FILENAME = TRIM('retv_vars.10684_0068_002.cdf')
      TES(2222)%FILENAME = TRIM('retv_vars.10684_0068_003.cdf')
      TES(2223)%FILENAME = TRIM('retv_vars.10684_0106_004.cdf')
      TES(2224)%FILENAME = TRIM('retv_vars.10684_0108_004.cdf')
      TES(2225)%FILENAME = TRIM('retv_vars.10684_0115_004.cdf')
      TES(2226)%FILENAME = TRIM('retv_vars.10684_0188_003.cdf')
      TES(2227)%FILENAME = TRIM('retv_vars.10684_0188_004.cdf')
      TES(2228)%FILENAME = TRIM('retv_vars.10684_0200_004.cdf')
      TES(2229)%FILENAME = TRIM('retv_vars.10684_0202_002.cdf')
      TES(2230)%FILENAME = TRIM('retv_vars.10684_0221_002.cdf')
      TES(2231)%FILENAME = TRIM('retv_vars.10684_0231_002.cdf')
      TES(2232)%FILENAME = TRIM('retv_vars.10684_0234_003.cdf')
      TES(2233)%FILENAME = TRIM('retv_vars.10684_0236_004.cdf')
      TES(2234)%FILENAME = TRIM('retv_vars.10684_0237_002.cdf')
      TES(2235)%FILENAME = TRIM('retv_vars.10684_0237_003.cdf')
      TES(2236)%FILENAME = TRIM('retv_vars.10684_0237_004.cdf')
      TES(2237)%FILENAME = TRIM('retv_vars.10684_0244_004.cdf')
      TES(2238)%FILENAME = TRIM('retv_vars.10684_0245_002.cdf')
      TES(2239)%FILENAME = TRIM('retv_vars.10684_0246_003.cdf')
      TES(2240)%FILENAME = TRIM('retv_vars.10684_0247_002.cdf')
      TES(2241)%FILENAME = TRIM('retv_vars.10684_0247_003.cdf')
      TES(2242)%FILENAME = TRIM('retv_vars.10684_0247_004.cdf')
      TES(2243)%FILENAME = TRIM('retv_vars.10684_0248_002.cdf')
      TES(2244)%FILENAME = TRIM('retv_vars.10684_0248_003.cdf')
      TES(2245)%FILENAME = TRIM('retv_vars.10684_0248_004.cdf')
      TES(2246)%FILENAME = TRIM('retv_vars.10684_0249_002.cdf')
      TES(2247)%FILENAME = TRIM('retv_vars.10684_0249_003.cdf')
      TES(2248)%FILENAME = TRIM('retv_vars.10684_0249_004.cdf')
      TES(2249)%FILENAME = TRIM('retv_vars.10684_0250_002.cdf')
      TES(2250)%FILENAME = TRIM('retv_vars.10684_0250_003.cdf')
      TES(2251)%FILENAME = TRIM('retv_vars.10684_0250_004.cdf')
      TES(2252)%FILENAME = TRIM('retv_vars.10684_0251_002.cdf')
      TES(2253)%FILENAME = TRIM('retv_vars.10684_0252_003.cdf')
      TES(2254)%FILENAME = TRIM('retv_vars.10684_0259_002.cdf')
      TES(2255)%FILENAME = TRIM('retv_vars.10684_0259_003.cdf')
      TES(2256)%FILENAME = TRIM('retv_vars.10684_0259_004.cdf')
      TES(2257)%FILENAME = TRIM('retv_vars.10684_0260_002.cdf')
      TES(2258)%FILENAME = TRIM('retv_vars.10684_0260_003.cdf')
      TES(2259)%FILENAME = TRIM('retv_vars.10684_0260_004.cdf')
      TES(2260)%FILENAME = TRIM('retv_vars.10684_0261_003.cdf')
      TES(2261)%FILENAME = TRIM('retv_vars.10684_0267_002.cdf')
      TES(2262)%FILENAME = TRIM('retv_vars.10684_0267_004.cdf')
      TES(2263)%FILENAME = TRIM('retv_vars.10684_0268_002.cdf')
      TES(2264)%FILENAME = TRIM('retv_vars.10684_0269_002.cdf')
      TES(2265)%FILENAME = TRIM('retv_vars.10684_0269_003.cdf')
      TES(2266)%FILENAME = TRIM('retv_vars.10684_0271_002.cdf')
      TES(2267)%FILENAME = TRIM('retv_vars.10684_0272_003.cdf')
      TES(2268)%FILENAME = TRIM('retv_vars.10684_0272_004.cdf')
      TES(2269)%FILENAME = TRIM('retv_vars.10684_0273_002.cdf')
      TES(2270)%FILENAME = TRIM('retv_vars.10684_0273_003.cdf')
      TES(2271)%FILENAME = TRIM('retv_vars.10684_0276_003.cdf')
      TES(2272)%FILENAME = TRIM('retv_vars.10684_0276_004.cdf')
      TES(2273)%FILENAME = TRIM('retv_vars.10684_0277_002.cdf')
      TES(2274)%FILENAME = TRIM('retv_vars.10684_0278_004.cdf')
      TES(2275)%FILENAME = TRIM('retv_vars.10684_0279_003.cdf')
      TES(2276)%FILENAME = TRIM('retv_vars.10684_0279_004.cdf')
      TES(2277)%FILENAME = TRIM('retv_vars.10684_0305_004.cdf')
      TES(2278)%FILENAME = TRIM('retv_vars.10684_0306_004.cdf')
      TES(2279)%FILENAME = TRIM('retv_vars.10684_0308_003.cdf')
      TES(2280)%FILENAME = TRIM('retv_vars.10684_0309_004.cdf')
      TES(2281)%FILENAME = TRIM('retv_vars.10684_0310_002.cdf')
      TES(2282)%FILENAME = TRIM('retv_vars.10684_0310_003.cdf')
      TES(2283)%FILENAME = TRIM('retv_vars.10684_0315_003.cdf')
      TES(2284)%FILENAME = TRIM('retv_vars.10684_0315_004.cdf')
      TES(2285)%FILENAME = TRIM('retv_vars.10684_0316_002.cdf')
      TES(2286)%FILENAME = TRIM('retv_vars.10684_0316_003.cdf')
      TES(2287)%FILENAME = TRIM('retv_vars.10684_0316_004.cdf')
      TES(2288)%FILENAME = TRIM('retv_vars.10684_0317_004.cdf')
      TES(2289)%FILENAME = TRIM('retv_vars.10684_0320_002.cdf')
      TES(2290)%FILENAME = TRIM('retv_vars.10684_0358_002.cdf')
      TES(2291)%FILENAME = TRIM('retv_vars.10684_0363_004.cdf')
      TES(2292)%FILENAME = TRIM('retv_vars.10684_0364_002.cdf')
      TES(2293)%FILENAME = TRIM('retv_vars.10684_0364_003.cdf')
      TES(2294)%FILENAME = TRIM('retv_vars.10684_0365_002.cdf')
      TES(2295)%FILENAME = TRIM('retv_vars.10684_0365_004.cdf')
      TES(2296)%FILENAME = TRIM('retv_vars.10684_0366_002.cdf')
      TES(2297)%FILENAME = TRIM('retv_vars.10684_0367_004.cdf')
      TES(2298)%FILENAME = TRIM('retv_vars.10684_0371_003.cdf')
      TES(2299)%FILENAME = TRIM('retv_vars.10684_0378_004.cdf')
      TES(2300)%FILENAME = TRIM('retv_vars.10684_0411_002.cdf')
      TES(2301)%FILENAME = TRIM('retv_vars.10684_0411_003.cdf')
      TES(2302)%FILENAME = TRIM('retv_vars.10684_0412_003.cdf')
      TES(2303)%FILENAME = TRIM('retv_vars.10684_0413_002.cdf')
      TES(2304)%FILENAME = TRIM('retv_vars.10684_0413_004.cdf')
      TES(2305)%FILENAME = TRIM('retv_vars.10684_0414_002.cdf')
      TES(2306)%FILENAME = TRIM('retv_vars.10684_0414_004.cdf')
      TES(2307)%FILENAME = TRIM('retv_vars.10684_0415_002.cdf')
      TES(2308)%FILENAME = TRIM('retv_vars.10684_0415_003.cdf')
      TES(2309)%FILENAME = TRIM('retv_vars.10684_0415_004.cdf')
      TES(2310)%FILENAME = TRIM('retv_vars.10684_0416_003.cdf')
      TES(2311)%FILENAME = TRIM('retv_vars.10684_0419_004.cdf')
      TES(2312)%FILENAME = TRIM('retv_vars.10684_0421_002.cdf')
      TES(2313)%FILENAME = TRIM('retv_vars.10684_0421_003.cdf')
      TES(2314)%FILENAME = TRIM('retv_vars.10684_0422_002.cdf')
      TES(2315)%FILENAME = TRIM('retv_vars.10684_0422_003.cdf')
      TES(2316)%FILENAME = TRIM('retv_vars.10684_0423_003.cdf')
      TES(2317)%FILENAME = TRIM('retv_vars.10684_0423_004.cdf')
      TES(2318)%FILENAME = TRIM('retv_vars.10684_0424_003.cdf')
      TES(2319)%FILENAME = TRIM('retv_vars.10684_0424_004.cdf')
      TES(2320)%FILENAME = TRIM('retv_vars.10684_0425_002.cdf')
      TES(2321)%FILENAME = TRIM('retv_vars.10684_0425_003.cdf')
      TES(2322)%FILENAME = TRIM('retv_vars.10684_0425_004.cdf')
      TES(2323)%FILENAME = TRIM('retv_vars.10684_0426_002.cdf')
      TES(2324)%FILENAME = TRIM('retv_vars.10684_0429_004.cdf')
      TES(2325)%FILENAME = TRIM('retv_vars.10684_0460_004.cdf')
      TES(2326)%FILENAME = TRIM('retv_vars.10684_0461_002.cdf')
      TES(2327)%FILENAME = TRIM('retv_vars.10684_0461_003.cdf')
      TES(2328)%FILENAME = TRIM('retv_vars.10684_0461_004.cdf')
      TES(2329)%FILENAME = TRIM('retv_vars.10684_0462_002.cdf')
      TES(2330)%FILENAME = TRIM('retv_vars.10684_0464_003.cdf')
      TES(2331)%FILENAME = TRIM('retv_vars.10684_0465_003.cdf')
      TES(2332)%FILENAME = TRIM('retv_vars.10684_0469_002.cdf')
      TES(2333)%FILENAME = TRIM('retv_vars.10684_0501_004.cdf')
      TES(2334)%FILENAME = TRIM('retv_vars.10684_0502_002.cdf')
      TES(2335)%FILENAME = TRIM('retv_vars.10684_0507_003.cdf')
      TES(2336)%FILENAME = TRIM('retv_vars.10684_0515_002.cdf')
      TES(2337)%FILENAME = TRIM('retv_vars.10684_0533_003.cdf')
      TES(2338)%FILENAME = TRIM('retv_vars.10684_0533_004.cdf')
      TES(2339)%FILENAME = TRIM('retv_vars.10684_0534_002.cdf')
      TES(2340)%FILENAME = TRIM('retv_vars.10684_0546_003.cdf')
      TES(2341)%FILENAME = TRIM('retv_vars.10684_0548_002.cdf')
      TES(2342)%FILENAME = TRIM('retv_vars.10684_0548_003.cdf')
      TES(2343)%FILENAME = TRIM('retv_vars.10684_0548_004.cdf')
      TES(2344)%FILENAME = TRIM('retv_vars.10684_0549_002.cdf')
      TES(2345)%FILENAME = TRIM('retv_vars.10684_0551_002.cdf')
      TES(2346)%FILENAME = TRIM('retv_vars.10684_0567_004.cdf')
      TES(2347)%FILENAME = TRIM('retv_vars.10684_0568_004.cdf')
      TES(2348)%FILENAME = TRIM('retv_vars.10684_0569_003.cdf')
      TES(2349)%FILENAME = TRIM('retv_vars.10684_0569_004.cdf')
      TES(2350)%FILENAME = TRIM('retv_vars.10684_0570_002.cdf')
      TES(2351)%FILENAME = TRIM('retv_vars.10684_0570_003.cdf')
      TES(2352)%FILENAME = TRIM('retv_vars.10684_0570_004.cdf')
      TES(2353)%FILENAME = TRIM('retv_vars.10684_0571_002.cdf')
      TES(2354)%FILENAME = TRIM('retv_vars.10684_0571_003.cdf')
      TES(2355)%FILENAME = TRIM('retv_vars.10684_0571_004.cdf')
      TES(2356)%FILENAME = TRIM('retv_vars.10684_0572_003.cdf')
      TES(2357)%FILENAME = TRIM('retv_vars.10684_0573_002.cdf')
      TES(2358)%FILENAME = TRIM('retv_vars.10684_0573_004.cdf')
      TES(2359)%FILENAME = TRIM('retv_vars.10684_0580_004.cdf')
      TES(2360)%FILENAME = TRIM('retv_vars.10684_0586_002.cdf')
      TES(2361)%FILENAME = TRIM('retv_vars.10684_0591_004.cdf')
      TES(2362)%FILENAME = TRIM('retv_vars.10684_0592_004.cdf')
      TES(2363)%FILENAME = TRIM('retv_vars.10684_0593_003.cdf')
      TES(2364)%FILENAME = TRIM('retv_vars.10684_0593_004.cdf')
      TES(2365)%FILENAME = TRIM('retv_vars.10684_0594_003.cdf')
      TES(2366)%FILENAME = TRIM('retv_vars.10684_0594_004.cdf')
      TES(2367)%FILENAME = TRIM('retv_vars.10684_0595_002.cdf')
      TES(2368)%FILENAME = TRIM('retv_vars.10684_0595_004.cdf')
      TES(2369)%FILENAME = TRIM('retv_vars.10684_0596_003.cdf')
      TES(2370)%FILENAME = TRIM('retv_vars.10684_0597_004.cdf')
      TES(2371)%FILENAME = TRIM('retv_vars.10684_0598_002.cdf')
      TES(2372)%FILENAME = TRIM('retv_vars.10684_0599_002.cdf')
      TES(2373)%FILENAME = TRIM('retv_vars.10684_0604_003.cdf')
      TES(2374)%FILENAME = TRIM('retv_vars.10684_0604_004.cdf')
      TES(2375)%FILENAME = TRIM('retv_vars.10684_0611_003.cdf')
      TES(2376)%FILENAME = TRIM('retv_vars.10684_0613_002.cdf')
      TES(2377)%FILENAME = TRIM('retv_vars.10684_0613_003.cdf')
      TES(2378)%FILENAME = TRIM('retv_vars.10684_0613_004.cdf')
      TES(2379)%FILENAME = TRIM('retv_vars.10684_0614_003.cdf')
      TES(2380)%FILENAME = TRIM('retv_vars.10684_0614_004.cdf')
      TES(2381)%FILENAME = TRIM('retv_vars.10684_0615_003.cdf')
      TES(2382)%FILENAME = TRIM('retv_vars.10684_0616_002.cdf')
      TES(2383)%FILENAME = TRIM('retv_vars.10684_0616_003.cdf')
      TES(2384)%FILENAME = TRIM('retv_vars.10684_0616_004.cdf')
      TES(2385)%FILENAME = TRIM('retv_vars.10684_0617_002.cdf')
      TES(2386)%FILENAME = TRIM('retv_vars.10684_0617_003.cdf')
      TES(2387)%FILENAME = TRIM('retv_vars.10684_0639_003.cdf')
      TES(2388)%FILENAME = TRIM('retv_vars.10684_0639_004.cdf')
      TES(2389)%FILENAME = TRIM('retv_vars.10684_0640_002.cdf')
      TES(2390)%FILENAME = TRIM('retv_vars.10684_0642_004.cdf')
      TES(2391)%FILENAME = TRIM('retv_vars.10684_0643_004.cdf')
      TES(2392)%FILENAME = TRIM('retv_vars.10684_0644_004.cdf')
      TES(2393)%FILENAME = TRIM('retv_vars.10684_0645_002.cdf')
      TES(2394)%FILENAME = TRIM('retv_vars.10684_0645_003.cdf')
      TES(2395)%FILENAME = TRIM('retv_vars.10684_0645_004.cdf')
      TES(2396)%FILENAME = TRIM('retv_vars.10684_0646_003.cdf')
      TES(2397)%FILENAME = TRIM('retv_vars.10684_0646_004.cdf')
      TES(2398)%FILENAME = TRIM('retv_vars.10684_0647_002.cdf')
      TES(2399)%FILENAME = TRIM('retv_vars.10684_0651_004.cdf')
      TES(2400)%FILENAME = TRIM('retv_vars.10684_0652_003.cdf')
      TES(2401)%FILENAME = TRIM('retv_vars.10684_0652_004.cdf')
      TES(2402)%FILENAME = TRIM('retv_vars.10684_0653_003.cdf')
      TES(2403)%FILENAME = TRIM('retv_vars.10684_0654_002.cdf')
      TES(2404)%FILENAME = TRIM('retv_vars.10684_0654_004.cdf')
      TES(2405)%FILENAME = TRIM('retv_vars.10684_0655_002.cdf')
      TES(2406)%FILENAME = TRIM('retv_vars.10684_0655_004.cdf')
      TES(2407)%FILENAME = TRIM('retv_vars.10684_0656_002.cdf')
      TES(2408)%FILENAME = TRIM('retv_vars.10684_0656_003.cdf')
      TES(2409)%FILENAME = TRIM('retv_vars.10684_0656_004.cdf')
      TES(2410)%FILENAME = TRIM('retv_vars.10684_0659_003.cdf')
      TES(2411)%FILENAME = TRIM('retv_vars.10684_0659_004.cdf')
      TES(2412)%FILENAME = TRIM('retv_vars.10684_0660_002.cdf')
      TES(2413)%FILENAME = TRIM('retv_vars.10684_0693_004.cdf')
      TES(2414)%FILENAME = TRIM('retv_vars.10684_0694_002.cdf')
      TES(2415)%FILENAME = TRIM('retv_vars.10684_0699_003.cdf')
      TES(2416)%FILENAME = TRIM('retv_vars.10684_0700_002.cdf')
      TES(2417)%FILENAME = TRIM('retv_vars.10684_0700_003.cdf')
      TES(2418)%FILENAME = TRIM('retv_vars.10684_0700_004.cdf')
      TES(2419)%FILENAME = TRIM('retv_vars.10684_0701_002.cdf')
      TES(2420)%FILENAME = TRIM('retv_vars.10684_0701_004.cdf')
      TES(2421)%FILENAME = TRIM('retv_vars.10684_0702_003.cdf')
      TES(2422)%FILENAME = TRIM('retv_vars.10684_0703_002.cdf')
      TES(2423)%FILENAME = TRIM('retv_vars.10684_0703_004.cdf')
      TES(2424)%FILENAME = TRIM('retv_vars.10684_0704_002.cdf')
      TES(2425)%FILENAME = TRIM('retv_vars.10684_0704_003.cdf')
      TES(2426)%FILENAME = TRIM('retv_vars.10684_0705_003.cdf')
      TES(2427)%FILENAME = TRIM('retv_vars.10684_0732_002.cdf')
      TES(2428)%FILENAME = TRIM('retv_vars.10684_0739_002.cdf')
      TES(2429)%FILENAME = TRIM('retv_vars.10684_0742_003.cdf')
      TES(2430)%FILENAME = TRIM('retv_vars.10684_0747_003.cdf')
      TES(2431)%FILENAME = TRIM('retv_vars.10684_0747_004.cdf')
      TES(2432)%FILENAME = TRIM('retv_vars.10686_0021_003.cdf')
      TES(2433)%FILENAME = TRIM('retv_vars.10686_0027_002.cdf')
      TES(2434)%FILENAME = TRIM('retv_vars.10686_0027_003.cdf')
      TES(2435)%FILENAME = TRIM('retv_vars.10686_0027_004.cdf')
      TES(2436)%FILENAME = TRIM('retv_vars.10686_0028_004.cdf')
      TES(2437)%FILENAME = TRIM('retv_vars.10686_0029_002.cdf')
      TES(2438)%FILENAME = TRIM('retv_vars.10686_0029_003.cdf')
      TES(2439)%FILENAME = TRIM('retv_vars.10686_0030_004.cdf')
      TES(2440)%FILENAME = TRIM('retv_vars.10686_0055_004.cdf')
      TES(2441)%FILENAME = TRIM('retv_vars.10686_0056_003.cdf')
      TES(2442)%FILENAME = TRIM('retv_vars.10686_0057_003.cdf')
      TES(2443)%FILENAME = TRIM('retv_vars.10686_0057_004.cdf')
      TES(2444)%FILENAME = TRIM('retv_vars.10686_0058_002.cdf')
      TES(2445)%FILENAME = TRIM('retv_vars.10686_0058_003.cdf')
      TES(2446)%FILENAME = TRIM('retv_vars.10686_0061_002.cdf')
      TES(2447)%FILENAME = TRIM('retv_vars.10686_0061_003.cdf')
      TES(2448)%FILENAME = TRIM('retv_vars.10686_0061_004.cdf')
      TES(2449)%FILENAME = TRIM('retv_vars.10686_0063_003.cdf')
      TES(2450)%FILENAME = TRIM('retv_vars.10686_0065_002.cdf')
      TES(2451)%FILENAME = TRIM('retv_vars.10686_0066_003.cdf')
      TES(2452)%FILENAME = TRIM('retv_vars.10686_0066_004.cdf')
      TES(2453)%FILENAME = TRIM('retv_vars.10686_0067_002.cdf')
      TES(2454)%FILENAME = TRIM('retv_vars.10686_0067_003.cdf')
      TES(2455)%FILENAME = TRIM('retv_vars.10686_0067_004.cdf')
      TES(2456)%FILENAME = TRIM('retv_vars.10686_0068_002.cdf')
      TES(2457)%FILENAME = TRIM('retv_vars.10686_0068_003.cdf')
      TES(2458)%FILENAME = TRIM('retv_vars.10686_0068_004.cdf')
      TES(2459)%FILENAME = TRIM('retv_vars.10686_0069_002.cdf')
      TES(2460)%FILENAME = TRIM('retv_vars.10686_0069_003.cdf')
      TES(2461)%FILENAME = TRIM('retv_vars.10686_0071_002.cdf')
      TES(2462)%FILENAME = TRIM('retv_vars.10686_0100_002.cdf')
      TES(2463)%FILENAME = TRIM('retv_vars.10686_0100_004.cdf')
      TES(2464)%FILENAME = TRIM('retv_vars.10686_0101_004.cdf')
      TES(2465)%FILENAME = TRIM('retv_vars.10686_0102_003.cdf')
      TES(2466)%FILENAME = TRIM('retv_vars.10686_0102_004.cdf')
      TES(2467)%FILENAME = TRIM('retv_vars.10686_0103_003.cdf')
      TES(2468)%FILENAME = TRIM('retv_vars.10686_0103_004.cdf')
      TES(2469)%FILENAME = TRIM('retv_vars.10686_0104_003.cdf')
      TES(2470)%FILENAME = TRIM('retv_vars.10686_0104_004.cdf')
      TES(2471)%FILENAME = TRIM('retv_vars.10686_0105_003.cdf')
      TES(2472)%FILENAME = TRIM('retv_vars.10686_0106_004.cdf')
      TES(2473)%FILENAME = TRIM('retv_vars.10686_0107_003.cdf')
      TES(2474)%FILENAME = TRIM('retv_vars.10686_0107_004.cdf')
      TES(2475)%FILENAME = TRIM('retv_vars.10686_0108_002.cdf')
      TES(2476)%FILENAME = TRIM('retv_vars.10686_0108_004.cdf')
      TES(2477)%FILENAME = TRIM('retv_vars.10686_0109_002.cdf')
      TES(2478)%FILENAME = TRIM('retv_vars.10686_0115_002.cdf')
      TES(2479)%FILENAME = TRIM('retv_vars.10686_0115_004.cdf')
      TES(2480)%FILENAME = TRIM('retv_vars.10686_0116_002.cdf')
      TES(2481)%FILENAME = TRIM('retv_vars.10686_0116_003.cdf')
      TES(2482)%FILENAME = TRIM('retv_vars.10686_0137_004.cdf')
      TES(2483)%FILENAME = TRIM('retv_vars.10686_0144_004.cdf')
      TES(2484)%FILENAME = TRIM('retv_vars.10686_0188_002.cdf')
      TES(2485)%FILENAME = TRIM('retv_vars.10686_0188_003.cdf')
      TES(2486)%FILENAME = TRIM('retv_vars.10686_0188_004.cdf')
      TES(2487)%FILENAME = TRIM('retv_vars.10686_0190_004.cdf')
      TES(2488)%FILENAME = TRIM('retv_vars.10686_0191_004.cdf')
      TES(2489)%FILENAME = TRIM('retv_vars.10686_0231_002.cdf')
      TES(2490)%FILENAME = TRIM('retv_vars.10686_0234_002.cdf')
      TES(2491)%FILENAME = TRIM('retv_vars.10686_0235_002.cdf')
      TES(2492)%FILENAME = TRIM('retv_vars.10686_0237_003.cdf')
      TES(2493)%FILENAME = TRIM('retv_vars.10686_0237_004.cdf')
      TES(2494)%FILENAME = TRIM('retv_vars.10686_0246_002.cdf')
      TES(2495)%FILENAME = TRIM('retv_vars.10686_0246_003.cdf')
      TES(2496)%FILENAME = TRIM('retv_vars.10686_0247_003.cdf')
      TES(2497)%FILENAME = TRIM('retv_vars.10686_0247_004.cdf')
      TES(2498)%FILENAME = TRIM('retv_vars.10686_0248_003.cdf')
      TES(2499)%FILENAME = TRIM('retv_vars.10686_0248_004.cdf')
      TES(2500)%FILENAME = TRIM('retv_vars.10686_0249_004.cdf')
      TES(2501)%FILENAME = TRIM('retv_vars.10686_0250_003.cdf')
      TES(2502)%FILENAME = TRIM('retv_vars.10686_0250_004.cdf')
      TES(2503)%FILENAME = TRIM('retv_vars.10686_0251_002.cdf')
      TES(2504)%FILENAME = TRIM('retv_vars.10686_0251_003.cdf')
      TES(2505)%FILENAME = TRIM('retv_vars.10686_0251_004.cdf')
      TES(2506)%FILENAME = TRIM('retv_vars.10686_0252_002.cdf')
      TES(2507)%FILENAME = TRIM('retv_vars.10686_0252_003.cdf')
      TES(2508)%FILENAME = TRIM('retv_vars.10686_0252_004.cdf')
      TES(2509)%FILENAME = TRIM('retv_vars.10686_0253_002.cdf')
      TES(2510)%FILENAME = TRIM('retv_vars.10686_0260_002.cdf')
      TES(2511)%FILENAME = TRIM('retv_vars.10686_0260_003.cdf')
      TES(2512)%FILENAME = TRIM('retv_vars.10686_0260_004.cdf')
      TES(2513)%FILENAME = TRIM('retv_vars.10686_0261_002.cdf')
      TES(2514)%FILENAME = TRIM('retv_vars.10686_0261_003.cdf')
      TES(2515)%FILENAME = TRIM('retv_vars.10686_0261_004.cdf')
      TES(2516)%FILENAME = TRIM('retv_vars.10686_0262_002.cdf')
      TES(2517)%FILENAME = TRIM('retv_vars.10686_0267_003.cdf')
      TES(2518)%FILENAME = TRIM('retv_vars.10686_0267_004.cdf')
      TES(2519)%FILENAME = TRIM('retv_vars.10686_0268_002.cdf')
      TES(2520)%FILENAME = TRIM('retv_vars.10686_0268_003.cdf')
      TES(2521)%FILENAME = TRIM('retv_vars.10686_0268_004.cdf')
      TES(2522)%FILENAME = TRIM('retv_vars.10686_0269_002.cdf')
      TES(2523)%FILENAME = TRIM('retv_vars.10686_0269_004.cdf')
      TES(2524)%FILENAME = TRIM('retv_vars.10686_0270_004.cdf')
      TES(2525)%FILENAME = TRIM('retv_vars.10686_0271_002.cdf')
      TES(2526)%FILENAME = TRIM('retv_vars.10686_0271_003.cdf')
      TES(2527)%FILENAME = TRIM('retv_vars.10686_0271_004.cdf')
      TES(2528)%FILENAME = TRIM('retv_vars.10686_0273_003.cdf')
      TES(2529)%FILENAME = TRIM('retv_vars.10686_0273_004.cdf')
      TES(2530)%FILENAME = TRIM('retv_vars.10686_0274_002.cdf')
      TES(2531)%FILENAME = TRIM('retv_vars.10686_0274_003.cdf')
      TES(2532)%FILENAME = TRIM('retv_vars.10686_0274_004.cdf')
      TES(2533)%FILENAME = TRIM('retv_vars.10686_0275_002.cdf')
      TES(2534)%FILENAME = TRIM('retv_vars.10686_0275_003.cdf')
      TES(2535)%FILENAME = TRIM('retv_vars.10686_0276_002.cdf')
      TES(2536)%FILENAME = TRIM('retv_vars.10686_0289_002.cdf')
      TES(2537)%FILENAME = TRIM('retv_vars.10686_0289_003.cdf')
      TES(2538)%FILENAME = TRIM('retv_vars.10686_0302_002.cdf')
      TES(2539)%FILENAME = TRIM('retv_vars.10686_0302_003.cdf')
      TES(2540)%FILENAME = TRIM('retv_vars.10686_0302_004.cdf')
      TES(2541)%FILENAME = TRIM('retv_vars.10686_0306_002.cdf')
      TES(2542)%FILENAME = TRIM('retv_vars.10686_0306_003.cdf')
      TES(2543)%FILENAME = TRIM('retv_vars.10686_0308_002.cdf')
      TES(2544)%FILENAME = TRIM('retv_vars.10686_0308_003.cdf')
      TES(2545)%FILENAME = TRIM('retv_vars.10686_0308_004.cdf')
      TES(2546)%FILENAME = TRIM('retv_vars.10686_0309_003.cdf')
      TES(2547)%FILENAME = TRIM('retv_vars.10686_0315_003.cdf')
      TES(2548)%FILENAME = TRIM('retv_vars.10686_0315_004.cdf')
      TES(2549)%FILENAME = TRIM('retv_vars.10686_0316_002.cdf')
      TES(2550)%FILENAME = TRIM('retv_vars.10686_0316_003.cdf')
      TES(2551)%FILENAME = TRIM('retv_vars.10686_0316_004.cdf')
      TES(2552)%FILENAME = TRIM('retv_vars.10686_0318_003.cdf')
      TES(2553)%FILENAME = TRIM('retv_vars.10686_0319_002.cdf')
      TES(2554)%FILENAME = TRIM('retv_vars.10686_0319_003.cdf')
      TES(2555)%FILENAME = TRIM('retv_vars.10686_0319_004.cdf')
      TES(2556)%FILENAME = TRIM('retv_vars.10686_0320_003.cdf')
      TES(2557)%FILENAME = TRIM('retv_vars.10686_0320_004.cdf')
      TES(2558)%FILENAME = TRIM('retv_vars.10686_0321_004.cdf')
      TES(2559)%FILENAME = TRIM('retv_vars.10686_0354_002.cdf')
      TES(2560)%FILENAME = TRIM('retv_vars.10686_0354_003.cdf')
      TES(2561)%FILENAME = TRIM('retv_vars.10686_0356_004.cdf')
      TES(2562)%FILENAME = TRIM('retv_vars.10686_0357_003.cdf')
      TES(2563)%FILENAME = TRIM('retv_vars.10686_0358_002.cdf')
      TES(2564)%FILENAME = TRIM('retv_vars.10686_0358_004.cdf')
      TES(2565)%FILENAME = TRIM('retv_vars.10686_0359_002.cdf')
      TES(2566)%FILENAME = TRIM('retv_vars.10686_0363_004.cdf')
      TES(2567)%FILENAME = TRIM('retv_vars.10686_0364_002.cdf')
      TES(2568)%FILENAME = TRIM('retv_vars.10686_0364_003.cdf')
      TES(2569)%FILENAME = TRIM('retv_vars.10686_0364_004.cdf')
      TES(2570)%FILENAME = TRIM('retv_vars.10686_0365_002.cdf')
      TES(2571)%FILENAME = TRIM('retv_vars.10686_0365_003.cdf')
      TES(2572)%FILENAME = TRIM('retv_vars.10686_0366_002.cdf')
      TES(2573)%FILENAME = TRIM('retv_vars.10686_0366_004.cdf')
      TES(2574)%FILENAME = TRIM('retv_vars.10686_0367_002.cdf')
      TES(2575)%FILENAME = TRIM('retv_vars.10686_0367_003.cdf')
      TES(2576)%FILENAME = TRIM('retv_vars.10686_0367_004.cdf')
      TES(2577)%FILENAME = TRIM('retv_vars.10686_0369_004.cdf')
      TES(2578)%FILENAME = TRIM('retv_vars.10686_0370_002.cdf')
      TES(2579)%FILENAME = TRIM('retv_vars.10686_0413_002.cdf')
      TES(2580)%FILENAME = TRIM('retv_vars.10686_0413_003.cdf')
      TES(2581)%FILENAME = TRIM('retv_vars.10686_0413_004.cdf')
      TES(2582)%FILENAME = TRIM('retv_vars.10686_0414_002.cdf')
      TES(2583)%FILENAME = TRIM('retv_vars.10686_0414_003.cdf')
      TES(2584)%FILENAME = TRIM('retv_vars.10686_0415_002.cdf')
      TES(2585)%FILENAME = TRIM('retv_vars.10686_0415_003.cdf')
      TES(2586)%FILENAME = TRIM('retv_vars.10686_0420_002.cdf')
      TES(2587)%FILENAME = TRIM('retv_vars.10686_0420_004.cdf')
      TES(2588)%FILENAME = TRIM('retv_vars.10686_0421_002.cdf')
      TES(2589)%FILENAME = TRIM('retv_vars.10686_0421_003.cdf')
      TES(2590)%FILENAME = TRIM('retv_vars.10686_0421_004.cdf')
      TES(2591)%FILENAME = TRIM('retv_vars.10686_0422_002.cdf')
      TES(2592)%FILENAME = TRIM('retv_vars.10686_0422_004.cdf')
      TES(2593)%FILENAME = TRIM('retv_vars.10686_0423_003.cdf')
      TES(2594)%FILENAME = TRIM('retv_vars.10686_0424_002.cdf')
      TES(2595)%FILENAME = TRIM('retv_vars.10686_0424_004.cdf')
      TES(2596)%FILENAME = TRIM('retv_vars.10686_0425_002.cdf')
      TES(2597)%FILENAME = TRIM('retv_vars.10686_0426_003.cdf')
      TES(2598)%FILENAME = TRIM('retv_vars.10686_0427_002.cdf')
      TES(2599)%FILENAME = TRIM('retv_vars.10686_0427_003.cdf')
      TES(2600)%FILENAME = TRIM('retv_vars.10686_0428_002.cdf')
      TES(2601)%FILENAME = TRIM('retv_vars.10686_0429_002.cdf')
      TES(2602)%FILENAME = TRIM('retv_vars.10686_0429_003.cdf')
      TES(2603)%FILENAME = TRIM('retv_vars.10686_0429_004.cdf')
      TES(2604)%FILENAME = TRIM('retv_vars.10686_0459_004.cdf')
      TES(2605)%FILENAME = TRIM('retv_vars.10686_0460_004.cdf')
      TES(2606)%FILENAME = TRIM('retv_vars.10686_0461_002.cdf')
      TES(2607)%FILENAME = TRIM('retv_vars.10686_0461_003.cdf')
      TES(2608)%FILENAME = TRIM('retv_vars.10686_0462_002.cdf')
      TES(2609)%FILENAME = TRIM('retv_vars.10686_0463_002.cdf')
      TES(2610)%FILENAME = TRIM('retv_vars.10686_0463_003.cdf')
      TES(2611)%FILENAME = TRIM('retv_vars.10686_0466_003.cdf')
      TES(2612)%FILENAME = TRIM('retv_vars.10686_0466_004.cdf')
      TES(2613)%FILENAME = TRIM('retv_vars.10686_0467_002.cdf')
      TES(2614)%FILENAME = TRIM('retv_vars.10686_0467_003.cdf')
      TES(2615)%FILENAME = TRIM('retv_vars.10686_0467_004.cdf')
      TES(2616)%FILENAME = TRIM('retv_vars.10686_0469_002.cdf')
      TES(2617)%FILENAME = TRIM('retv_vars.10686_0469_003.cdf')
      TES(2618)%FILENAME = TRIM('retv_vars.10686_0469_004.cdf')
      TES(2619)%FILENAME = TRIM('retv_vars.10686_0471_002.cdf')
      TES(2620)%FILENAME = TRIM('retv_vars.10686_0502_002.cdf')
      TES(2621)%FILENAME = TRIM('retv_vars.10686_0502_003.cdf')
      TES(2622)%FILENAME = TRIM('retv_vars.10686_0502_004.cdf')
      TES(2623)%FILENAME = TRIM('retv_vars.10686_0508_003.cdf')
      TES(2624)%FILENAME = TRIM('retv_vars.10686_0510_004.cdf')
      TES(2625)%FILENAME = TRIM('retv_vars.10686_0511_002.cdf')
      TES(2626)%FILENAME = TRIM('retv_vars.10686_0511_003.cdf')
      TES(2627)%FILENAME = TRIM('retv_vars.10686_0515_003.cdf')
      TES(2628)%FILENAME = TRIM('retv_vars.10686_0515_004.cdf')
      TES(2629)%FILENAME = TRIM('retv_vars.10686_0516_002.cdf')
      TES(2630)%FILENAME = TRIM('retv_vars.10686_0516_004.cdf')
      TES(2631)%FILENAME = TRIM('retv_vars.10686_0542_003.cdf')
      TES(2632)%FILENAME = TRIM('retv_vars.10686_0549_002.cdf')
      TES(2633)%FILENAME = TRIM('retv_vars.10686_0549_003.cdf')
      TES(2634)%FILENAME = TRIM('retv_vars.10686_0549_004.cdf')
      TES(2635)%FILENAME = TRIM('retv_vars.10686_0550_002.cdf')
      TES(2636)%FILENAME = TRIM('retv_vars.10686_0550_003.cdf')
      TES(2637)%FILENAME = TRIM('retv_vars.10686_0551_002.cdf')
      TES(2638)%FILENAME = TRIM('retv_vars.10686_0568_002.cdf')
      TES(2639)%FILENAME = TRIM('retv_vars.10686_0568_004.cdf')
      TES(2640)%FILENAME = TRIM('retv_vars.10686_0569_002.cdf')
      TES(2641)%FILENAME = TRIM('retv_vars.10686_0569_004.cdf')
      TES(2642)%FILENAME = TRIM('retv_vars.10686_0570_003.cdf')
      TES(2643)%FILENAME = TRIM('retv_vars.10686_0570_004.cdf')
      TES(2644)%FILENAME = TRIM('retv_vars.10686_0571_002.cdf')
      TES(2645)%FILENAME = TRIM('retv_vars.10686_0571_003.cdf')
      TES(2646)%FILENAME = TRIM('retv_vars.10686_0572_002.cdf')
      TES(2647)%FILENAME = TRIM('retv_vars.10686_0572_003.cdf')
      TES(2648)%FILENAME = TRIM('retv_vars.10686_0582_004.cdf')
      TES(2649)%FILENAME = TRIM('retv_vars.10686_0583_002.cdf')
      TES(2650)%FILENAME = TRIM('retv_vars.10686_0583_004.cdf')
      TES(2651)%FILENAME = TRIM('retv_vars.10686_0587_002.cdf')
      TES(2652)%FILENAME = TRIM('retv_vars.10686_0588_003.cdf')
      TES(2653)%FILENAME = TRIM('retv_vars.10686_0589_002.cdf')
      TES(2654)%FILENAME = TRIM('retv_vars.10686_0589_003.cdf')
      TES(2655)%FILENAME = TRIM('retv_vars.10686_0590_002.cdf')
      TES(2656)%FILENAME = TRIM('retv_vars.10686_0592_004.cdf')
      TES(2657)%FILENAME = TRIM('retv_vars.10686_0594_002.cdf')
      TES(2658)%FILENAME = TRIM('retv_vars.10686_0594_004.cdf')
      TES(2659)%FILENAME = TRIM('retv_vars.10686_0595_004.cdf')
      TES(2660)%FILENAME = TRIM('retv_vars.10686_0597_002.cdf')
      TES(2661)%FILENAME = TRIM('retv_vars.10686_0598_003.cdf')
      TES(2662)%FILENAME = TRIM('retv_vars.10686_0598_004.cdf')
      TES(2663)%FILENAME = TRIM('retv_vars.10686_0599_002.cdf')
      TES(2664)%FILENAME = TRIM('retv_vars.10686_0605_002.cdf')
      TES(2665)%FILENAME = TRIM('retv_vars.10686_0613_004.cdf')
      TES(2666)%FILENAME = TRIM('retv_vars.10686_0614_002.cdf')
      TES(2667)%FILENAME = TRIM('retv_vars.10686_0615_002.cdf')
      TES(2668)%FILENAME = TRIM('retv_vars.10686_0616_003.cdf')
      TES(2669)%FILENAME = TRIM('retv_vars.10686_0617_002.cdf')
      TES(2670)%FILENAME = TRIM('retv_vars.10686_0617_003.cdf')
      TES(2671)%FILENAME = TRIM('retv_vars.10686_0618_002.cdf')
      TES(2672)%FILENAME = TRIM('retv_vars.10686_0618_003.cdf')
      TES(2673)%FILENAME = TRIM('retv_vars.10686_0635_003.cdf')
      TES(2674)%FILENAME = TRIM('retv_vars.10686_0639_004.cdf')
      TES(2675)%FILENAME = TRIM('retv_vars.10686_0640_003.cdf')
      TES(2676)%FILENAME = TRIM('retv_vars.10686_0640_004.cdf')
      TES(2677)%FILENAME = TRIM('retv_vars.10686_0645_002.cdf')
      TES(2678)%FILENAME = TRIM('retv_vars.10686_0645_003.cdf')
      TES(2679)%FILENAME = TRIM('retv_vars.10686_0645_004.cdf')
      TES(2680)%FILENAME = TRIM('retv_vars.10686_0651_002.cdf')
      TES(2681)%FILENAME = TRIM('retv_vars.10686_0651_003.cdf')
      TES(2682)%FILENAME = TRIM('retv_vars.10686_0651_004.cdf')
      TES(2683)%FILENAME = TRIM('retv_vars.10686_0652_002.cdf')
      TES(2684)%FILENAME = TRIM('retv_vars.10686_0652_004.cdf')
      TES(2685)%FILENAME = TRIM('retv_vars.10686_0653_002.cdf')
      TES(2686)%FILENAME = TRIM('retv_vars.10686_0653_003.cdf')
      TES(2687)%FILENAME = TRIM('retv_vars.10686_0653_004.cdf')
      TES(2688)%FILENAME = TRIM('retv_vars.10686_0654_002.cdf')
      TES(2689)%FILENAME = TRIM('retv_vars.10686_0654_003.cdf')
      TES(2690)%FILENAME = TRIM('retv_vars.10686_0655_004.cdf')
      TES(2691)%FILENAME = TRIM('retv_vars.10686_0656_003.cdf')
      TES(2692)%FILENAME = TRIM('retv_vars.10686_0656_004.cdf')
      TES(2693)%FILENAME = TRIM('retv_vars.10686_0657_002.cdf')
      TES(2694)%FILENAME = TRIM('retv_vars.10686_0659_003.cdf')
      TES(2695)%FILENAME = TRIM('retv_vars.10686_0659_004.cdf')
      TES(2696)%FILENAME = TRIM('retv_vars.10686_0660_002.cdf')
      TES(2697)%FILENAME = TRIM('retv_vars.10686_0694_003.cdf')
      TES(2698)%FILENAME = TRIM('retv_vars.10686_0694_004.cdf')
      TES(2699)%FILENAME = TRIM('retv_vars.10686_0695_002.cdf')
      TES(2700)%FILENAME = TRIM('retv_vars.10686_0700_004.cdf')
      TES(2701)%FILENAME = TRIM('retv_vars.10686_0701_002.cdf')
      TES(2702)%FILENAME = TRIM('retv_vars.10686_0701_003.cdf')
      TES(2703)%FILENAME = TRIM('retv_vars.10686_0701_004.cdf')
      TES(2704)%FILENAME = TRIM('retv_vars.10686_0702_002.cdf')
      TES(2705)%FILENAME = TRIM('retv_vars.10686_0702_003.cdf')
      TES(2706)%FILENAME = TRIM('retv_vars.10686_0703_002.cdf')
      TES(2707)%FILENAME = TRIM('retv_vars.10686_0703_003.cdf')
      TES(2708)%FILENAME = TRIM('retv_vars.10686_0703_004.cdf')
      TES(2709)%FILENAME = TRIM('retv_vars.10686_0704_003.cdf')
      TES(2710)%FILENAME = TRIM('retv_vars.10686_0704_004.cdf')
      TES(2711)%FILENAME = TRIM('retv_vars.10686_0705_003.cdf')
      TES(2712)%FILENAME = TRIM('retv_vars.10686_0706_002.cdf')
      TES(2713)%FILENAME = TRIM('retv_vars.10686_0733_002.cdf')
      TES(2714)%FILENAME = TRIM('retv_vars.10686_0739_003.cdf')
      TES(2715)%FILENAME = TRIM('retv_vars.10686_0740_002.cdf')
      TES(2716)%FILENAME = TRIM('retv_vars.10686_0740_004.cdf')
      TES(2717)%FILENAME = TRIM('retv_vars.10686_0741_002.cdf')
      TES(2718)%FILENAME = TRIM('retv_vars.10686_0741_003.cdf')
      TES(2719)%FILENAME = TRIM('retv_vars.10686_0747_003.cdf')
      TES(2720)%FILENAME = TRIM('retv_vars.10686_0748_002.cdf')
      TES(2721)%FILENAME = TRIM('retv_vars.10686_0748_003.cdf')
      TES(2722)%FILENAME = TRIM('retv_vars.10686_0748_004.cdf')
      TES(2723)%FILENAME = TRIM('retv_vars.10686_0749_002.cdf')
      TES(2724)%FILENAME = TRIM('retv_vars.10688_0016_004.cdf')
      TES(2725)%FILENAME = TRIM('retv_vars.10688_0021_002.cdf')
      TES(2726)%FILENAME = TRIM('retv_vars.10688_0021_003.cdf')
      TES(2727)%FILENAME = TRIM('retv_vars.10688_0022_003.cdf')
      TES(2728)%FILENAME = TRIM('retv_vars.10688_0027_003.cdf')
      TES(2729)%FILENAME = TRIM('retv_vars.10688_0027_004.cdf')
      TES(2730)%FILENAME = TRIM('retv_vars.10688_0028_003.cdf')
      TES(2731)%FILENAME = TRIM('retv_vars.10688_0028_004.cdf')
      TES(2732)%FILENAME = TRIM('retv_vars.10688_0029_002.cdf')
      TES(2733)%FILENAME = TRIM('retv_vars.10688_0029_003.cdf')
      TES(2734)%FILENAME = TRIM('retv_vars.10688_0030_002.cdf')
      TES(2735)%FILENAME = TRIM('retv_vars.10688_0030_003.cdf')
      TES(2736)%FILENAME = TRIM('retv_vars.10688_0030_004.cdf')
      TES(2737)%FILENAME = TRIM('retv_vars.10688_0031_002.cdf')
      TES(2738)%FILENAME = TRIM('retv_vars.10688_0054_004.cdf')
      TES(2739)%FILENAME = TRIM('retv_vars.10688_0060_003.cdf')
      TES(2740)%FILENAME = TRIM('retv_vars.10688_0060_004.cdf')
      TES(2741)%FILENAME = TRIM('retv_vars.10688_0061_003.cdf')
      TES(2742)%FILENAME = TRIM('retv_vars.10688_0068_002.cdf')
      TES(2743)%FILENAME = TRIM('retv_vars.10688_0068_003.cdf')
      TES(2744)%FILENAME = TRIM('retv_vars.10688_0068_004.cdf')
      TES(2745)%FILENAME = TRIM('retv_vars.10688_0069_002.cdf')
      TES(2746)%FILENAME = TRIM('retv_vars.10688_0069_003.cdf')
      TES(2747)%FILENAME = TRIM('retv_vars.10688_0069_004.cdf')
      TES(2748)%FILENAME = TRIM('retv_vars.10688_0070_003.cdf')
      TES(2749)%FILENAME = TRIM('retv_vars.10688_0100_003.cdf')
      TES(2750)%FILENAME = TRIM('retv_vars.10688_0100_004.cdf')
      TES(2751)%FILENAME = TRIM('retv_vars.10688_0103_002.cdf')
      TES(2752)%FILENAME = TRIM('retv_vars.10688_0103_004.cdf')
      TES(2753)%FILENAME = TRIM('retv_vars.10688_0104_003.cdf')
      TES(2754)%FILENAME = TRIM('retv_vars.10688_0104_004.cdf')
      TES(2755)%FILENAME = TRIM('retv_vars.10688_0105_002.cdf')
      TES(2756)%FILENAME = TRIM('retv_vars.10688_0105_003.cdf')
      TES(2757)%FILENAME = TRIM('retv_vars.10688_0105_004.cdf')
      TES(2758)%FILENAME = TRIM('retv_vars.10688_0106_003.cdf')
      TES(2759)%FILENAME = TRIM('retv_vars.10688_0106_004.cdf')
      TES(2760)%FILENAME = TRIM('retv_vars.10688_0108_003.cdf')
      TES(2761)%FILENAME = TRIM('retv_vars.10688_0108_004.cdf')
      TES(2762)%FILENAME = TRIM('retv_vars.10688_0110_003.cdf')
      TES(2763)%FILENAME = TRIM('retv_vars.10688_0110_004.cdf')
      TES(2764)%FILENAME = TRIM('retv_vars.10688_0117_002.cdf')
      TES(2765)%FILENAME = TRIM('retv_vars.10688_0117_004.cdf')
      TES(2766)%FILENAME = TRIM('retv_vars.10688_0118_002.cdf')
      TES(2767)%FILENAME = TRIM('retv_vars.10688_0119_002.cdf')
      TES(2768)%FILENAME = TRIM('retv_vars.10688_0124_002.cdf')
      TES(2769)%FILENAME = TRIM('retv_vars.10688_0138_004.cdf')
      TES(2770)%FILENAME = TRIM('retv_vars.10688_0157_002.cdf')
      TES(2771)%FILENAME = TRIM('retv_vars.10688_0157_003.cdf')
      TES(2772)%FILENAME = TRIM('retv_vars.10688_0158_004.cdf')
      TES(2773)%FILENAME = TRIM('retv_vars.10688_0188_003.cdf')
      TES(2774)%FILENAME = TRIM('retv_vars.10688_0189_002.cdf')
      TES(2775)%FILENAME = TRIM('retv_vars.10688_0189_003.cdf')
      TES(2776)%FILENAME = TRIM('retv_vars.10688_0189_004.cdf')
      TES(2777)%FILENAME = TRIM('retv_vars.10688_0190_002.cdf')
      TES(2778)%FILENAME = TRIM('retv_vars.10688_0190_003.cdf')
      TES(2779)%FILENAME = TRIM('retv_vars.10688_0190_004.cdf')
      TES(2780)%FILENAME = TRIM('retv_vars.10688_0191_004.cdf')
      TES(2781)%FILENAME = TRIM('retv_vars.10688_0233_003.cdf')
      TES(2782)%FILENAME = TRIM('retv_vars.10688_0234_002.cdf')
      TES(2783)%FILENAME = TRIM('retv_vars.10688_0234_003.cdf')
      TES(2784)%FILENAME = TRIM('retv_vars.10688_0234_004.cdf')
      TES(2785)%FILENAME = TRIM('retv_vars.10688_0237_002.cdf')
      TES(2786)%FILENAME = TRIM('retv_vars.10688_0237_003.cdf')
      TES(2787)%FILENAME = TRIM('retv_vars.10688_0244_004.cdf')
      TES(2788)%FILENAME = TRIM('retv_vars.10688_0245_003.cdf')
      TES(2789)%FILENAME = TRIM('retv_vars.10688_0245_004.cdf')
      TES(2790)%FILENAME = TRIM('retv_vars.10688_0246_002.cdf')
      TES(2791)%FILENAME = TRIM('retv_vars.10688_0248_003.cdf')
      TES(2792)%FILENAME = TRIM('retv_vars.10688_0249_002.cdf')
      TES(2793)%FILENAME = TRIM('retv_vars.10688_0250_002.cdf')
      TES(2794)%FILENAME = TRIM('retv_vars.10688_0250_003.cdf')
      TES(2795)%FILENAME = TRIM('retv_vars.10688_0250_004.cdf')
      TES(2796)%FILENAME = TRIM('retv_vars.10688_0251_004.cdf')
      TES(2797)%FILENAME = TRIM('retv_vars.10688_0252_002.cdf')
      TES(2798)%FILENAME = TRIM('retv_vars.10688_0252_003.cdf')
      TES(2799)%FILENAME = TRIM('retv_vars.10688_0252_004.cdf')
      TES(2800)%FILENAME = TRIM('retv_vars.10688_0253_002.cdf')
      TES(2801)%FILENAME = TRIM('retv_vars.10688_0253_003.cdf')
      TES(2802)%FILENAME = TRIM('retv_vars.10688_0260_002.cdf')
      TES(2803)%FILENAME = TRIM('retv_vars.10688_0260_004.cdf')
      TES(2804)%FILENAME = TRIM('retv_vars.10688_0261_002.cdf')
      TES(2805)%FILENAME = TRIM('retv_vars.10688_0261_003.cdf')
      TES(2806)%FILENAME = TRIM('retv_vars.10688_0261_004.cdf')
      TES(2807)%FILENAME = TRIM('retv_vars.10688_0262_003.cdf')
      TES(2808)%FILENAME = TRIM('retv_vars.10688_0262_004.cdf')
      TES(2809)%FILENAME = TRIM('retv_vars.10688_0263_002.cdf')
      TES(2810)%FILENAME = TRIM('retv_vars.10688_0268_004.cdf')
      TES(2811)%FILENAME = TRIM('retv_vars.10688_0269_002.cdf')
      TES(2812)%FILENAME = TRIM('retv_vars.10688_0269_004.cdf')
      TES(2813)%FILENAME = TRIM('retv_vars.10688_0273_004.cdf')
      TES(2814)%FILENAME = TRIM('retv_vars.10688_0274_003.cdf')
      TES(2815)%FILENAME = TRIM('retv_vars.10688_0289_002.cdf')
      TES(2816)%FILENAME = TRIM('retv_vars.10688_0289_004.cdf')
      TES(2817)%FILENAME = TRIM('retv_vars.10688_0290_002.cdf')
      TES(2818)%FILENAME = TRIM('retv_vars.10688_0301_003.cdf')
      TES(2819)%FILENAME = TRIM('retv_vars.10688_0302_002.cdf')
      TES(2820)%FILENAME = TRIM('retv_vars.10688_0303_002.cdf')
      TES(2821)%FILENAME = TRIM('retv_vars.10688_0305_003.cdf')
      TES(2822)%FILENAME = TRIM('retv_vars.10688_0305_004.cdf')
      TES(2823)%FILENAME = TRIM('retv_vars.10688_0306_002.cdf')
      TES(2824)%FILENAME = TRIM('retv_vars.10688_0306_004.cdf')
      TES(2825)%FILENAME = TRIM('retv_vars.10688_0307_002.cdf')
      TES(2826)%FILENAME = TRIM('retv_vars.10688_0307_004.cdf')
      TES(2827)%FILENAME = TRIM('retv_vars.10688_0308_002.cdf')
      TES(2828)%FILENAME = TRIM('retv_vars.10688_0308_004.cdf')
      TES(2829)%FILENAME = TRIM('retv_vars.10688_0309_004.cdf')
      TES(2830)%FILENAME = TRIM('retv_vars.10688_0311_002.cdf')
      TES(2831)%FILENAME = TRIM('retv_vars.10688_0317_002.cdf')
      TES(2832)%FILENAME = TRIM('retv_vars.10688_0318_003.cdf')
      TES(2833)%FILENAME = TRIM('retv_vars.10688_0319_003.cdf')
      TES(2834)%FILENAME = TRIM('retv_vars.10688_0319_004.cdf')
      TES(2835)%FILENAME = TRIM('retv_vars.10688_0320_003.cdf')
      TES(2836)%FILENAME = TRIM('retv_vars.10688_0321_002.cdf')
      TES(2837)%FILENAME = TRIM('retv_vars.10688_0321_004.cdf')
      TES(2838)%FILENAME = TRIM('retv_vars.10688_0322_002.cdf')
      TES(2839)%FILENAME = TRIM('retv_vars.10688_0353_004.cdf')
      TES(2840)%FILENAME = TRIM('retv_vars.10688_0354_004.cdf')
      TES(2841)%FILENAME = TRIM('retv_vars.10688_0357_002.cdf')
      TES(2842)%FILENAME = TRIM('retv_vars.10688_0358_003.cdf')
      TES(2843)%FILENAME = TRIM('retv_vars.10688_0363_004.cdf')
      TES(2844)%FILENAME = TRIM('retv_vars.10688_0364_003.cdf')
      TES(2845)%FILENAME = TRIM('retv_vars.10688_0364_004.cdf')
      TES(2846)%FILENAME = TRIM('retv_vars.10688_0365_002.cdf')
      TES(2847)%FILENAME = TRIM('retv_vars.10688_0365_003.cdf')
      TES(2848)%FILENAME = TRIM('retv_vars.10688_0366_002.cdf')
      TES(2849)%FILENAME = TRIM('retv_vars.10688_0366_003.cdf')
      TES(2850)%FILENAME = TRIM('retv_vars.10688_0367_003.cdf')
      TES(2851)%FILENAME = TRIM('retv_vars.10688_0368_002.cdf')
      TES(2852)%FILENAME = TRIM('retv_vars.10688_0411_003.cdf')
      TES(2853)%FILENAME = TRIM('retv_vars.10688_0411_004.cdf')
      TES(2854)%FILENAME = TRIM('retv_vars.10688_0412_003.cdf')
      TES(2855)%FILENAME = TRIM('retv_vars.10688_0413_002.cdf')
      TES(2856)%FILENAME = TRIM('retv_vars.10688_0414_002.cdf')
      TES(2857)%FILENAME = TRIM('retv_vars.10688_0415_002.cdf')
      TES(2858)%FILENAME = TRIM('retv_vars.10688_0415_004.cdf')
      TES(2859)%FILENAME = TRIM('retv_vars.10688_0418_004.cdf')
      TES(2860)%FILENAME = TRIM('retv_vars.10688_0420_002.cdf')
      TES(2861)%FILENAME = TRIM('retv_vars.10688_0421_004.cdf')
      TES(2862)%FILENAME = TRIM('retv_vars.10688_0422_002.cdf')
      TES(2863)%FILENAME = TRIM('retv_vars.10688_0422_003.cdf')
      TES(2864)%FILENAME = TRIM('retv_vars.10688_0422_004.cdf')
      TES(2865)%FILENAME = TRIM('retv_vars.10688_0423_004.cdf')
      TES(2866)%FILENAME = TRIM('retv_vars.10688_0424_002.cdf')
      TES(2867)%FILENAME = TRIM('retv_vars.10688_0424_003.cdf')
      TES(2868)%FILENAME = TRIM('retv_vars.10688_0424_004.cdf')
      TES(2869)%FILENAME = TRIM('retv_vars.10688_0425_002.cdf')
      TES(2870)%FILENAME = TRIM('retv_vars.10688_0425_003.cdf')
      TES(2871)%FILENAME = TRIM('retv_vars.10688_0426_004.cdf')
      TES(2872)%FILENAME = TRIM('retv_vars.10688_0427_002.cdf')
      TES(2873)%FILENAME = TRIM('retv_vars.10688_0427_003.cdf')
      TES(2874)%FILENAME = TRIM('retv_vars.10688_0427_004.cdf')
      TES(2875)%FILENAME = TRIM('retv_vars.10688_0428_003.cdf')
      TES(2876)%FILENAME = TRIM('retv_vars.10688_0429_002.cdf')
      TES(2877)%FILENAME = TRIM('retv_vars.10688_0429_003.cdf')
      TES(2878)%FILENAME = TRIM('retv_vars.10688_0429_004.cdf')
      TES(2879)%FILENAME = TRIM('retv_vars.10688_0459_004.cdf')
      TES(2880)%FILENAME = TRIM('retv_vars.10688_0460_002.cdf')
      TES(2881)%FILENAME = TRIM('retv_vars.10688_0460_003.cdf')
      TES(2882)%FILENAME = TRIM('retv_vars.10688_0461_002.cdf')
      TES(2883)%FILENAME = TRIM('retv_vars.10688_0461_004.cdf')
      TES(2884)%FILENAME = TRIM('retv_vars.10688_0462_002.cdf')
      TES(2885)%FILENAME = TRIM('retv_vars.10688_0462_003.cdf')
      TES(2886)%FILENAME = TRIM('retv_vars.10688_0465_004.cdf')
      TES(2887)%FILENAME = TRIM('retv_vars.10688_0466_002.cdf')
      TES(2888)%FILENAME = TRIM('retv_vars.10688_0466_003.cdf')
      TES(2889)%FILENAME = TRIM('retv_vars.10688_0467_002.cdf')
      TES(2890)%FILENAME = TRIM('retv_vars.10688_0470_002.cdf')
      TES(2891)%FILENAME = TRIM('retv_vars.10688_0470_003.cdf')
      TES(2892)%FILENAME = TRIM('retv_vars.10688_0471_002.cdf')
      TES(2893)%FILENAME = TRIM('retv_vars.10688_0471_003.cdf')
      TES(2894)%FILENAME = TRIM('retv_vars.10688_0502_004.cdf')
      TES(2895)%FILENAME = TRIM('retv_vars.10688_0508_003.cdf')
      TES(2896)%FILENAME = TRIM('retv_vars.10688_0508_004.cdf')
      TES(2897)%FILENAME = TRIM('retv_vars.10688_0509_003.cdf')
      TES(2898)%FILENAME = TRIM('retv_vars.10688_0510_004.cdf')
      TES(2899)%FILENAME = TRIM('retv_vars.10688_0511_002.cdf')
      TES(2900)%FILENAME = TRIM('retv_vars.10688_0511_003.cdf')
      TES(2901)%FILENAME = TRIM('retv_vars.10688_0512_004.cdf')
      TES(2902)%FILENAME = TRIM('retv_vars.10688_0515_002.cdf')
      TES(2903)%FILENAME = TRIM('retv_vars.10688_0515_003.cdf')
      TES(2904)%FILENAME = TRIM('retv_vars.10688_0515_004.cdf')
      TES(2905)%FILENAME = TRIM('retv_vars.10688_0516_003.cdf')
      TES(2906)%FILENAME = TRIM('retv_vars.10688_0516_004.cdf')
      TES(2907)%FILENAME = TRIM('retv_vars.10688_0517_002.cdf')
      TES(2908)%FILENAME = TRIM('retv_vars.10688_0517_003.cdf')
      TES(2909)%FILENAME = TRIM('retv_vars.10688_0549_002.cdf')
      TES(2910)%FILENAME = TRIM('retv_vars.10688_0549_004.cdf')
      TES(2911)%FILENAME = TRIM('retv_vars.10688_0550_002.cdf')
      TES(2912)%FILENAME = TRIM('retv_vars.10688_0568_002.cdf')
      TES(2913)%FILENAME = TRIM('retv_vars.10688_0568_004.cdf')
      TES(2914)%FILENAME = TRIM('retv_vars.10688_0569_004.cdf')
      TES(2915)%FILENAME = TRIM('retv_vars.10688_0570_003.cdf')
      TES(2916)%FILENAME = TRIM('retv_vars.10688_0571_002.cdf')
      TES(2917)%FILENAME = TRIM('retv_vars.10688_0571_004.cdf')
      TES(2918)%FILENAME = TRIM('retv_vars.10688_0580_003.cdf')
      TES(2919)%FILENAME = TRIM('retv_vars.10688_0586_004.cdf')
      TES(2920)%FILENAME = TRIM('retv_vars.10688_0589_003.cdf')
      TES(2921)%FILENAME = TRIM('retv_vars.10688_0593_003.cdf')
      TES(2922)%FILENAME = TRIM('retv_vars.10688_0594_002.cdf')
      TES(2923)%FILENAME = TRIM('retv_vars.10688_0596_003.cdf')
      TES(2924)%FILENAME = TRIM('retv_vars.10688_0596_004.cdf')
      TES(2925)%FILENAME = TRIM('retv_vars.10688_0597_002.cdf')
      TES(2926)%FILENAME = TRIM('retv_vars.10688_0597_003.cdf')
      TES(2927)%FILENAME = TRIM('retv_vars.10688_0598_002.cdf')
      TES(2928)%FILENAME = TRIM('retv_vars.10688_0599_002.cdf')
      TES(2929)%FILENAME = TRIM('retv_vars.10688_0612_002.cdf')
      TES(2930)%FILENAME = TRIM('retv_vars.10688_0613_002.cdf')
      TES(2931)%FILENAME = TRIM('retv_vars.10688_0613_003.cdf')
      TES(2932)%FILENAME = TRIM('retv_vars.10688_0613_004.cdf')
      TES(2933)%FILENAME = TRIM('retv_vars.10688_0614_002.cdf')
      TES(2934)%FILENAME = TRIM('retv_vars.10688_0614_003.cdf')
      TES(2935)%FILENAME = TRIM('retv_vars.10688_0615_003.cdf')
      TES(2936)%FILENAME = TRIM('retv_vars.10688_0615_004.cdf')
      TES(2937)%FILENAME = TRIM('retv_vars.10688_0616_002.cdf')
      TES(2938)%FILENAME = TRIM('retv_vars.10688_0617_004.cdf')
      TES(2939)%FILENAME = TRIM('retv_vars.10688_0618_003.cdf')
      TES(2940)%FILENAME = TRIM('retv_vars.10688_0619_002.cdf')
      TES(2941)%FILENAME = TRIM('retv_vars.10688_0619_003.cdf')
      TES(2942)%FILENAME = TRIM('retv_vars.10688_0619_004.cdf')
      TES(2943)%FILENAME = TRIM('retv_vars.10688_0620_004.cdf')
      TES(2944)%FILENAME = TRIM('retv_vars.10688_0621_002.cdf')
      TES(2945)%FILENAME = TRIM('retv_vars.10688_0621_003.cdf')
      TES(2946)%FILENAME = TRIM('retv_vars.10688_0634_002.cdf')
      TES(2947)%FILENAME = TRIM('retv_vars.10688_0634_003.cdf')
      TES(2948)%FILENAME = TRIM('retv_vars.10688_0635_003.cdf')
      TES(2949)%FILENAME = TRIM('retv_vars.10688_0638_003.cdf')
      TES(2950)%FILENAME = TRIM('retv_vars.10688_0638_004.cdf')
      TES(2951)%FILENAME = TRIM('retv_vars.10688_0639_003.cdf')
      TES(2952)%FILENAME = TRIM('retv_vars.10688_0640_003.cdf')
      TES(2953)%FILENAME = TRIM('retv_vars.10688_0643_003.cdf')
      TES(2954)%FILENAME = TRIM('retv_vars.10688_0644_003.cdf')
      TES(2955)%FILENAME = TRIM('retv_vars.10688_0644_004.cdf')
      TES(2956)%FILENAME = TRIM('retv_vars.10688_0645_003.cdf')
      TES(2957)%FILENAME = TRIM('retv_vars.10688_0646_003.cdf')
      TES(2958)%FILENAME = TRIM('retv_vars.10688_0651_002.cdf')
      TES(2959)%FILENAME = TRIM('retv_vars.10688_0651_004.cdf')
      TES(2960)%FILENAME = TRIM('retv_vars.10688_0654_002.cdf')
      TES(2961)%FILENAME = TRIM('retv_vars.10688_0654_003.cdf')
      TES(2962)%FILENAME = TRIM('retv_vars.10688_0654_004.cdf')
      TES(2963)%FILENAME = TRIM('retv_vars.10688_0655_003.cdf')
      TES(2964)%FILENAME = TRIM('retv_vars.10688_0655_004.cdf')
      TES(2965)%FILENAME = TRIM('retv_vars.10688_0656_003.cdf')
      TES(2966)%FILENAME = TRIM('retv_vars.10688_0656_004.cdf')
      TES(2967)%FILENAME = TRIM('retv_vars.10688_0657_002.cdf')
      TES(2968)%FILENAME = TRIM('retv_vars.10688_0658_004.cdf')
      TES(2969)%FILENAME = TRIM('retv_vars.10688_0660_002.cdf')
      TES(2970)%FILENAME = TRIM('retv_vars.10688_0660_003.cdf')
      TES(2971)%FILENAME = TRIM('retv_vars.10688_0685_002.cdf')
      TES(2972)%FILENAME = TRIM('retv_vars.10688_0685_004.cdf')
      TES(2973)%FILENAME = TRIM('retv_vars.10688_0686_002.cdf')
      TES(2974)%FILENAME = TRIM('retv_vars.10688_0686_004.cdf')
      TES(2975)%FILENAME = TRIM('retv_vars.10688_0688_004.cdf')
      TES(2976)%FILENAME = TRIM('retv_vars.10688_0689_002.cdf')
      TES(2977)%FILENAME = TRIM('retv_vars.10688_0691_003.cdf')
      TES(2978)%FILENAME = TRIM('retv_vars.10688_0695_002.cdf')
      TES(2979)%FILENAME = TRIM('retv_vars.10688_0700_002.cdf')
      TES(2980)%FILENAME = TRIM('retv_vars.10688_0700_003.cdf')
      TES(2981)%FILENAME = TRIM('retv_vars.10688_0701_002.cdf')
      TES(2982)%FILENAME = TRIM('retv_vars.10688_0701_003.cdf')
      TES(2983)%FILENAME = TRIM('retv_vars.10688_0701_004.cdf')
      TES(2984)%FILENAME = TRIM('retv_vars.10688_0702_002.cdf')
      TES(2985)%FILENAME = TRIM('retv_vars.10688_0702_003.cdf')
      TES(2986)%FILENAME = TRIM('retv_vars.10688_0703_002.cdf')
      TES(2987)%FILENAME = TRIM('retv_vars.10688_0703_003.cdf')
      TES(2988)%FILENAME = TRIM('retv_vars.10688_0703_004.cdf')
      TES(2989)%FILENAME = TRIM('retv_vars.10688_0704_003.cdf')
      TES(2990)%FILENAME = TRIM('retv_vars.10688_0704_004.cdf')
      TES(2991)%FILENAME = TRIM('retv_vars.10688_0705_003.cdf')
      TES(2992)%FILENAME = TRIM('retv_vars.10688_0734_003.cdf')
      TES(2993)%FILENAME = TRIM('retv_vars.10688_0737_004.cdf')
      TES(2994)%FILENAME = TRIM('retv_vars.10688_0738_004.cdf')
      TES(2995)%FILENAME = TRIM('retv_vars.10688_0740_002.cdf')
      TES(2996)%FILENAME = TRIM('retv_vars.10688_0741_003.cdf')
      TES(2997)%FILENAME = TRIM('retv_vars.10688_0741_004.cdf')
      TES(2998)%FILENAME = TRIM('retv_vars.10688_0742_002.cdf')
      TES(2999)%FILENAME = TRIM('retv_vars.10688_0742_003.cdf')
      TES(3000)%FILENAME = TRIM('retv_vars.10688_0742_004.cdf')
      TES(3001)%FILENAME = TRIM('retv_vars.10688_0743_003.cdf')
      TES(3002)%FILENAME = TRIM('retv_vars.10688_0747_004.cdf')
      TES(3003)%FILENAME = TRIM('retv_vars.10688_0748_002.cdf')
      TES(3004)%FILENAME = TRIM('retv_vars.10688_0748_003.cdf')
      TES(3005)%FILENAME = TRIM('retv_vars.10688_0748_004.cdf')
      TES(3006)%FILENAME = TRIM('retv_vars.10688_0749_002.cdf')
      TES(3007)%FILENAME = TRIM('retv_vars.10693_0020_003.cdf')
      TES(3008)%FILENAME = TRIM('retv_vars.10693_0021_003.cdf')
      TES(3009)%FILENAME = TRIM('retv_vars.10693_0022_002.cdf')
      TES(3010)%FILENAME = TRIM('retv_vars.10693_0022_003.cdf')
      TES(3011)%FILENAME = TRIM('retv_vars.10693_0022_004.cdf')
      TES(3012)%FILENAME = TRIM('retv_vars.10693_0023_002.cdf')
      TES(3013)%FILENAME = TRIM('retv_vars.10693_0027_002.cdf')
      TES(3014)%FILENAME = TRIM('retv_vars.10693_0027_003.cdf')
      TES(3015)%FILENAME = TRIM('retv_vars.10693_0027_004.cdf')
      TES(3016)%FILENAME = TRIM('retv_vars.10693_0028_002.cdf')
      TES(3017)%FILENAME = TRIM('retv_vars.10693_0028_003.cdf')
      TES(3018)%FILENAME = TRIM('retv_vars.10693_0028_004.cdf')
      TES(3019)%FILENAME = TRIM('retv_vars.10693_0029_002.cdf')
      TES(3020)%FILENAME = TRIM('retv_vars.10693_0029_003.cdf')
      TES(3021)%FILENAME = TRIM('retv_vars.10693_0029_004.cdf')
      TES(3022)%FILENAME = TRIM('retv_vars.10693_0030_002.cdf')
      TES(3023)%FILENAME = TRIM('retv_vars.10693_0030_003.cdf')
      TES(3024)%FILENAME = TRIM('retv_vars.10693_0031_002.cdf')
      TES(3025)%FILENAME = TRIM('retv_vars.10693_0031_003.cdf')
      TES(3026)%FILENAME = TRIM('retv_vars.10693_0031_004.cdf')
      TES(3027)%FILENAME = TRIM('retv_vars.10693_0032_002.cdf')
      TES(3028)%FILENAME = TRIM('retv_vars.10693_0032_003.cdf')
      TES(3029)%FILENAME = TRIM('retv_vars.10693_0059_003.cdf')
      TES(3030)%FILENAME = TRIM('retv_vars.10693_0060_002.cdf')
      TES(3031)%FILENAME = TRIM('retv_vars.10693_0060_003.cdf')
      TES(3032)%FILENAME = TRIM('retv_vars.10693_0067_003.cdf')
      TES(3033)%FILENAME = TRIM('retv_vars.10693_0070_003.cdf')
      TES(3034)%FILENAME = TRIM('retv_vars.10693_0071_002.cdf')
      TES(3035)%FILENAME = TRIM('retv_vars.10693_0075_002.cdf')
      TES(3036)%FILENAME = TRIM('retv_vars.10693_0075_003.cdf')
      TES(3037)%FILENAME = TRIM('retv_vars.10693_0100_003.cdf')
      TES(3038)%FILENAME = TRIM('retv_vars.10693_0100_004.cdf')
      TES(3039)%FILENAME = TRIM('retv_vars.10693_0101_002.cdf')
      TES(3040)%FILENAME = TRIM('retv_vars.10693_0101_003.cdf')
      TES(3041)%FILENAME = TRIM('retv_vars.10693_0101_004.cdf')
      TES(3042)%FILENAME = TRIM('retv_vars.10693_0102_002.cdf')
      TES(3043)%FILENAME = TRIM('retv_vars.10693_0102_004.cdf')
      TES(3044)%FILENAME = TRIM('retv_vars.10693_0103_002.cdf')
      TES(3045)%FILENAME = TRIM('retv_vars.10693_0103_004.cdf')
      TES(3046)%FILENAME = TRIM('retv_vars.10693_0104_003.cdf')
      TES(3047)%FILENAME = TRIM('retv_vars.10693_0104_004.cdf')
      TES(3048)%FILENAME = TRIM('retv_vars.10693_0105_002.cdf')
      TES(3049)%FILENAME = TRIM('retv_vars.10693_0105_004.cdf')
      TES(3050)%FILENAME = TRIM('retv_vars.10693_0106_002.cdf')
      TES(3051)%FILENAME = TRIM('retv_vars.10693_0106_003.cdf')
      TES(3052)%FILENAME = TRIM('retv_vars.10693_0107_004.cdf')
      TES(3053)%FILENAME = TRIM('retv_vars.10693_0108_002.cdf')
      TES(3054)%FILENAME = TRIM('retv_vars.10693_0108_003.cdf')
      TES(3055)%FILENAME = TRIM('retv_vars.10693_0108_004.cdf')
      TES(3056)%FILENAME = TRIM('retv_vars.10693_0109_002.cdf')
      TES(3057)%FILENAME = TRIM('retv_vars.10693_0109_003.cdf')
      TES(3058)%FILENAME = TRIM('retv_vars.10693_0110_004.cdf')
      TES(3059)%FILENAME = TRIM('retv_vars.10693_0139_003.cdf')
      TES(3060)%FILENAME = TRIM('retv_vars.10693_0143_002.cdf')
      TES(3061)%FILENAME = TRIM('retv_vars.10693_0143_003.cdf')
      TES(3062)%FILENAME = TRIM('retv_vars.10693_0156_002.cdf')
      TES(3063)%FILENAME = TRIM('retv_vars.10693_0159_002.cdf')
      TES(3064)%FILENAME = TRIM('retv_vars.10693_0231_004.cdf')
      TES(3065)%FILENAME = TRIM('retv_vars.10693_0233_003.cdf')
      TES(3066)%FILENAME = TRIM('retv_vars.10693_0233_004.cdf')
      TES(3067)%FILENAME = TRIM('retv_vars.10693_0234_003.cdf')
      TES(3068)%FILENAME = TRIM('retv_vars.10693_0237_003.cdf')
      TES(3069)%FILENAME = TRIM('retv_vars.10693_0245_004.cdf')
      TES(3070)%FILENAME = TRIM('retv_vars.10693_0247_002.cdf')
      TES(3071)%FILENAME = TRIM('retv_vars.10693_0247_003.cdf')
      TES(3072)%FILENAME = TRIM('retv_vars.10693_0247_004.cdf')
      TES(3073)%FILENAME = TRIM('retv_vars.10693_0248_002.cdf')
      TES(3074)%FILENAME = TRIM('retv_vars.10693_0248_003.cdf')
      TES(3075)%FILENAME = TRIM('retv_vars.10693_0248_004.cdf')
      TES(3076)%FILENAME = TRIM('retv_vars.10693_0249_002.cdf')
      TES(3077)%FILENAME = TRIM('retv_vars.10693_0249_003.cdf')
      TES(3078)%FILENAME = TRIM('retv_vars.10693_0250_004.cdf')
      TES(3079)%FILENAME = TRIM('retv_vars.10693_0251_002.cdf')
      TES(3080)%FILENAME = TRIM('retv_vars.10693_0251_003.cdf')
      TES(3081)%FILENAME = TRIM('retv_vars.10693_0251_004.cdf')
      TES(3082)%FILENAME = TRIM('retv_vars.10693_0261_003.cdf')
      TES(3083)%FILENAME = TRIM('retv_vars.10693_0262_003.cdf')
      TES(3084)%FILENAME = TRIM('retv_vars.10693_0262_004.cdf')
      TES(3085)%FILENAME = TRIM('retv_vars.10693_0263_002.cdf')
      TES(3086)%FILENAME = TRIM('retv_vars.10693_0267_002.cdf')
      TES(3087)%FILENAME = TRIM('retv_vars.10693_0267_003.cdf')
      TES(3088)%FILENAME = TRIM('retv_vars.10693_0267_004.cdf')
      TES(3089)%FILENAME = TRIM('retv_vars.10693_0268_002.cdf')
      TES(3090)%FILENAME = TRIM('retv_vars.10693_0268_003.cdf')
      TES(3091)%FILENAME = TRIM('retv_vars.10693_0269_003.cdf')
      TES(3092)%FILENAME = TRIM('retv_vars.10693_0271_003.cdf')
      TES(3093)%FILENAME = TRIM('retv_vars.10693_0272_004.cdf')
      TES(3094)%FILENAME = TRIM('retv_vars.10693_0273_002.cdf')
      TES(3095)%FILENAME = TRIM('retv_vars.10693_0279_003.cdf')
      TES(3096)%FILENAME = TRIM('retv_vars.10693_0279_004.cdf')
      TES(3097)%FILENAME = TRIM('retv_vars.10693_0289_003.cdf')
      TES(3098)%FILENAME = TRIM('retv_vars.10693_0289_004.cdf')
      TES(3099)%FILENAME = TRIM('retv_vars.10693_0290_002.cdf')
      TES(3100)%FILENAME = TRIM('retv_vars.10693_0290_003.cdf')
      TES(3101)%FILENAME = TRIM('retv_vars.10693_0290_004.cdf')
      TES(3102)%FILENAME = TRIM('retv_vars.10693_0291_002.cdf')
      TES(3103)%FILENAME = TRIM('retv_vars.10693_0291_003.cdf')
      TES(3104)%FILENAME = TRIM('retv_vars.10693_0291_004.cdf')
      TES(3105)%FILENAME = TRIM('retv_vars.10693_0297_004.cdf')
      TES(3106)%FILENAME = TRIM('retv_vars.10693_0298_002.cdf')
      TES(3107)%FILENAME = TRIM('retv_vars.10693_0298_004.cdf')
      TES(3108)%FILENAME = TRIM('retv_vars.10693_0300_004.cdf')
      TES(3109)%FILENAME = TRIM('retv_vars.10693_0302_002.cdf')
      TES(3110)%FILENAME = TRIM('retv_vars.10693_0305_003.cdf')
      TES(3111)%FILENAME = TRIM('retv_vars.10693_0305_004.cdf')
      TES(3112)%FILENAME = TRIM('retv_vars.10693_0307_003.cdf')
      TES(3113)%FILENAME = TRIM('retv_vars.10693_0308_003.cdf')
      TES(3114)%FILENAME = TRIM('retv_vars.10693_0309_003.cdf')
      TES(3115)%FILENAME = TRIM('retv_vars.10693_0310_002.cdf')
      TES(3116)%FILENAME = TRIM('retv_vars.10693_0310_003.cdf')
      TES(3117)%FILENAME = TRIM('retv_vars.10693_0315_002.cdf')
      TES(3118)%FILENAME = TRIM('retv_vars.10693_0316_004.cdf')
      TES(3119)%FILENAME = TRIM('retv_vars.10693_0317_002.cdf')
      TES(3120)%FILENAME = TRIM('retv_vars.10693_0317_003.cdf')
      TES(3121)%FILENAME = TRIM('retv_vars.10693_0319_003.cdf')
      TES(3122)%FILENAME = TRIM('retv_vars.10693_0320_002.cdf')
      TES(3123)%FILENAME = TRIM('retv_vars.10693_0320_004.cdf')
      TES(3124)%FILENAME = TRIM('retv_vars.10693_0321_002.cdf')
      TES(3125)%FILENAME = TRIM('retv_vars.10693_0354_003.cdf')
      TES(3126)%FILENAME = TRIM('retv_vars.10693_0354_004.cdf')
      TES(3127)%FILENAME = TRIM('retv_vars.10693_0355_003.cdf')
      TES(3128)%FILENAME = TRIM('retv_vars.10693_0356_004.cdf')
      TES(3129)%FILENAME = TRIM('retv_vars.10693_0357_002.cdf')
      TES(3130)%FILENAME = TRIM('retv_vars.10693_0357_003.cdf')
      TES(3131)%FILENAME = TRIM('retv_vars.10693_0357_004.cdf')
      TES(3132)%FILENAME = TRIM('retv_vars.10693_0358_002.cdf')
      TES(3133)%FILENAME = TRIM('retv_vars.10693_0358_003.cdf')
      TES(3134)%FILENAME = TRIM('retv_vars.10693_0363_003.cdf')
      TES(3135)%FILENAME = TRIM('retv_vars.10693_0364_003.cdf')
      TES(3136)%FILENAME = TRIM('retv_vars.10693_0365_002.cdf')
      TES(3137)%FILENAME = TRIM('retv_vars.10693_0365_003.cdf')
      TES(3138)%FILENAME = TRIM('retv_vars.10693_0365_004.cdf')
      TES(3139)%FILENAME = TRIM('retv_vars.10693_0366_004.cdf')
      TES(3140)%FILENAME = TRIM('retv_vars.10693_0367_004.cdf')
      TES(3141)%FILENAME = TRIM('retv_vars.10693_0368_004.cdf')
      TES(3142)%FILENAME = TRIM('retv_vars.10693_0369_003.cdf')
      TES(3143)%FILENAME = TRIM('retv_vars.10693_0369_004.cdf')
      TES(3144)%FILENAME = TRIM('retv_vars.10693_0370_003.cdf')
      TES(3145)%FILENAME = TRIM('retv_vars.10693_0411_003.cdf')
      TES(3146)%FILENAME = TRIM('retv_vars.10693_0412_002.cdf')
      TES(3147)%FILENAME = TRIM('retv_vars.10693_0412_004.cdf')
      TES(3148)%FILENAME = TRIM('retv_vars.10693_0413_002.cdf')
      TES(3149)%FILENAME = TRIM('retv_vars.10693_0415_004.cdf')
      TES(3150)%FILENAME = TRIM('retv_vars.10693_0421_002.cdf')
      TES(3151)%FILENAME = TRIM('retv_vars.10693_0421_004.cdf')
      TES(3152)%FILENAME = TRIM('retv_vars.10693_0422_003.cdf')
      TES(3153)%FILENAME = TRIM('retv_vars.10693_0422_004.cdf')
      TES(3154)%FILENAME = TRIM('retv_vars.10693_0423_002.cdf')
      TES(3155)%FILENAME = TRIM('retv_vars.10693_0423_003.cdf')
      TES(3156)%FILENAME = TRIM('retv_vars.10693_0423_004.cdf')
      TES(3157)%FILENAME = TRIM('retv_vars.10693_0424_003.cdf')
      TES(3158)%FILENAME = TRIM('retv_vars.10693_0425_003.cdf')
      TES(3159)%FILENAME = TRIM('retv_vars.10693_0425_004.cdf')
      TES(3160)%FILENAME = TRIM('retv_vars.10693_0426_002.cdf')
      TES(3161)%FILENAME = TRIM('retv_vars.10693_0426_003.cdf')
      TES(3162)%FILENAME = TRIM('retv_vars.10693_0426_004.cdf')
      TES(3163)%FILENAME = TRIM('retv_vars.10693_0447_003.cdf')
      TES(3164)%FILENAME = TRIM('retv_vars.10693_0459_003.cdf')
      TES(3165)%FILENAME = TRIM('retv_vars.10693_0459_004.cdf')
      TES(3166)%FILENAME = TRIM('retv_vars.10693_0461_002.cdf')
      TES(3167)%FILENAME = TRIM('retv_vars.10693_0461_004.cdf')
      TES(3168)%FILENAME = TRIM('retv_vars.10693_0462_003.cdf')
      TES(3169)%FILENAME = TRIM('retv_vars.10693_0462_004.cdf')
      TES(3170)%FILENAME = TRIM('retv_vars.10693_0463_002.cdf')
      TES(3171)%FILENAME = TRIM('retv_vars.10693_0463_003.cdf')
      TES(3172)%FILENAME = TRIM('retv_vars.10693_0465_003.cdf')
      TES(3173)%FILENAME = TRIM('retv_vars.10693_0466_002.cdf')
      TES(3174)%FILENAME = TRIM('retv_vars.10693_0466_003.cdf')
      TES(3175)%FILENAME = TRIM('retv_vars.10693_0469_002.cdf')
      TES(3176)%FILENAME = TRIM('retv_vars.10693_0469_003.cdf')
      TES(3177)%FILENAME = TRIM('retv_vars.10693_0470_004.cdf')
      TES(3178)%FILENAME = TRIM('retv_vars.10693_0471_002.cdf')
      TES(3179)%FILENAME = TRIM('retv_vars.10693_0471_003.cdf')
      TES(3180)%FILENAME = TRIM('retv_vars.10693_0471_004.cdf')
      TES(3181)%FILENAME = TRIM('retv_vars.10693_0503_002.cdf')
      TES(3182)%FILENAME = TRIM('retv_vars.10693_0509_002.cdf')
      TES(3183)%FILENAME = TRIM('retv_vars.10693_0509_003.cdf')
      TES(3184)%FILENAME = TRIM('retv_vars.10693_0512_004.cdf')
      TES(3185)%FILENAME = TRIM('retv_vars.10693_0515_002.cdf')
      TES(3186)%FILENAME = TRIM('retv_vars.10693_0516_003.cdf')
      TES(3187)%FILENAME = TRIM('retv_vars.10693_0517_002.cdf')
      TES(3188)%FILENAME = TRIM('retv_vars.10693_0517_003.cdf')
      TES(3189)%FILENAME = TRIM('retv_vars.10693_0549_003.cdf')
      TES(3190)%FILENAME = TRIM('retv_vars.10693_0568_004.cdf')
      TES(3191)%FILENAME = TRIM('retv_vars.10693_0569_002.cdf')
      TES(3192)%FILENAME = TRIM('retv_vars.10693_0569_003.cdf')
      TES(3193)%FILENAME = TRIM('retv_vars.10693_0569_004.cdf')
      TES(3194)%FILENAME = TRIM('retv_vars.10693_0570_003.cdf')
      TES(3195)%FILENAME = TRIM('retv_vars.10693_0570_004.cdf')
      TES(3196)%FILENAME = TRIM('retv_vars.10693_0571_002.cdf')
      TES(3197)%FILENAME = TRIM('retv_vars.10693_0571_003.cdf')
      TES(3198)%FILENAME = TRIM('retv_vars.10693_0571_004.cdf')
      TES(3199)%FILENAME = TRIM('retv_vars.10693_0579_004.cdf')
      TES(3200)%FILENAME = TRIM('retv_vars.10693_0580_002.cdf')
      TES(3201)%FILENAME = TRIM('retv_vars.10693_0581_002.cdf')
      TES(3202)%FILENAME = TRIM('retv_vars.10693_0584_002.cdf')
      TES(3203)%FILENAME = TRIM('retv_vars.10693_0594_003.cdf')
      TES(3204)%FILENAME = TRIM('retv_vars.10693_0594_004.cdf')
      TES(3205)%FILENAME = TRIM('retv_vars.10693_0595_003.cdf')
      TES(3206)%FILENAME = TRIM('retv_vars.10693_0596_004.cdf')
      TES(3207)%FILENAME = TRIM('retv_vars.10693_0597_002.cdf')
      TES(3208)%FILENAME = TRIM('retv_vars.10693_0597_003.cdf')
      TES(3209)%FILENAME = TRIM('retv_vars.10693_0598_003.cdf')
      TES(3210)%FILENAME = TRIM('retv_vars.10693_0613_003.cdf')
      TES(3211)%FILENAME = TRIM('retv_vars.10693_0613_004.cdf')
      TES(3212)%FILENAME = TRIM('retv_vars.10693_0614_002.cdf')
      TES(3213)%FILENAME = TRIM('retv_vars.10693_0615_002.cdf')
      TES(3214)%FILENAME = TRIM('retv_vars.10693_0616_004.cdf')
      TES(3215)%FILENAME = TRIM('retv_vars.10693_0617_004.cdf')
      TES(3216)%FILENAME = TRIM('retv_vars.10693_0618_003.cdf')
      TES(3217)%FILENAME = TRIM('retv_vars.10693_0619_002.cdf')
      TES(3218)%FILENAME = TRIM('retv_vars.10693_0620_004.cdf')
      TES(3219)%FILENAME = TRIM('retv_vars.10693_0621_004.cdf')
      TES(3220)%FILENAME = TRIM('retv_vars.10693_0622_002.cdf')
      TES(3221)%FILENAME = TRIM('retv_vars.10693_0622_003.cdf')
      TES(3222)%FILENAME = TRIM('retv_vars.10693_0622_004.cdf')
      TES(3223)%FILENAME = TRIM('retv_vars.10693_0623_002.cdf')
      TES(3224)%FILENAME = TRIM('retv_vars.10693_0623_003.cdf')
      TES(3225)%FILENAME = TRIM('retv_vars.10693_0623_004.cdf')
      TES(3226)%FILENAME = TRIM('retv_vars.10693_0624_002.cdf')
      TES(3227)%FILENAME = TRIM('retv_vars.10693_0624_003.cdf')
      TES(3228)%FILENAME = TRIM('retv_vars.10693_0624_004.cdf')
      TES(3229)%FILENAME = TRIM('retv_vars.10693_0633_003.cdf')
      TES(3230)%FILENAME = TRIM('retv_vars.10693_0637_004.cdf')
      TES(3231)%FILENAME = TRIM('retv_vars.10693_0638_002.cdf')
      TES(3232)%FILENAME = TRIM('retv_vars.10693_0639_004.cdf')
      TES(3233)%FILENAME = TRIM('retv_vars.10693_0642_002.cdf')
      TES(3234)%FILENAME = TRIM('retv_vars.10693_0642_004.cdf')
      TES(3235)%FILENAME = TRIM('retv_vars.10693_0643_003.cdf')
      TES(3236)%FILENAME = TRIM('retv_vars.10693_0643_004.cdf')
      TES(3237)%FILENAME = TRIM('retv_vars.10693_0644_002.cdf')
      TES(3238)%FILENAME = TRIM('retv_vars.10693_0644_003.cdf')
      TES(3239)%FILENAME = TRIM('retv_vars.10693_0645_003.cdf')
      TES(3240)%FILENAME = TRIM('retv_vars.10693_0645_004.cdf')
      TES(3241)%FILENAME = TRIM('retv_vars.10693_0646_002.cdf')
      TES(3242)%FILENAME = TRIM('retv_vars.10693_0646_004.cdf')
      TES(3243)%FILENAME = TRIM('retv_vars.10693_0651_003.cdf')
      TES(3244)%FILENAME = TRIM('retv_vars.10693_0651_004.cdf')
      TES(3245)%FILENAME = TRIM('retv_vars.10693_0652_004.cdf')
      TES(3246)%FILENAME = TRIM('retv_vars.10693_0653_002.cdf')
      TES(3247)%FILENAME = TRIM('retv_vars.10693_0654_002.cdf')
      TES(3248)%FILENAME = TRIM('retv_vars.10693_0655_002.cdf')
      TES(3249)%FILENAME = TRIM('retv_vars.10693_0656_002.cdf')
      TES(3250)%FILENAME = TRIM('retv_vars.10693_0657_004.cdf')
      TES(3251)%FILENAME = TRIM('retv_vars.10693_0661_002.cdf')
      TES(3252)%FILENAME = TRIM('retv_vars.10693_0684_003.cdf')
      TES(3253)%FILENAME = TRIM('retv_vars.10693_0684_004.cdf')
      TES(3254)%FILENAME = TRIM('retv_vars.10693_0685_004.cdf')
      TES(3255)%FILENAME = TRIM('retv_vars.10693_0686_002.cdf')
      TES(3256)%FILENAME = TRIM('retv_vars.10693_0687_002.cdf')
      TES(3257)%FILENAME = TRIM('retv_vars.10693_0691_003.cdf')
      TES(3258)%FILENAME = TRIM('retv_vars.10693_0691_004.cdf')
      TES(3259)%FILENAME = TRIM('retv_vars.10693_0692_002.cdf')
      TES(3260)%FILENAME = TRIM('retv_vars.10693_0693_002.cdf')
      TES(3261)%FILENAME = TRIM('retv_vars.10693_0693_003.cdf')
      TES(3262)%FILENAME = TRIM('retv_vars.10693_0694_003.cdf')
      TES(3263)%FILENAME = TRIM('retv_vars.10693_0694_004.cdf')
      TES(3264)%FILENAME = TRIM('retv_vars.10693_0695_002.cdf')
      TES(3265)%FILENAME = TRIM('retv_vars.10693_0699_003.cdf')
      TES(3266)%FILENAME = TRIM('retv_vars.10693_0700_002.cdf')
      TES(3267)%FILENAME = TRIM('retv_vars.10693_0700_003.cdf')
      TES(3268)%FILENAME = TRIM('retv_vars.10693_0700_004.cdf')
      TES(3269)%FILENAME = TRIM('retv_vars.10693_0701_004.cdf')
      TES(3270)%FILENAME = TRIM('retv_vars.10693_0702_004.cdf')
      TES(3271)%FILENAME = TRIM('retv_vars.10693_0703_002.cdf')
      TES(3272)%FILENAME = TRIM('retv_vars.10693_0705_002.cdf')
      TES(3273)%FILENAME = TRIM('retv_vars.10693_0705_003.cdf')
      TES(3274)%FILENAME = TRIM('retv_vars.10693_0705_004.cdf')
      TES(3275)%FILENAME = TRIM('retv_vars.10693_0706_002.cdf')
      TES(3276)%FILENAME = TRIM('retv_vars.10693_0706_003.cdf')
      TES(3277)%FILENAME = TRIM('retv_vars.10693_0735_003.cdf')
      TES(3278)%FILENAME = TRIM('retv_vars.10693_0737_004.cdf')
      TES(3279)%FILENAME = TRIM('retv_vars.10693_0738_004.cdf')
      TES(3280)%FILENAME = TRIM('retv_vars.10693_0741_003.cdf')
      TES(3281)%FILENAME = TRIM('retv_vars.10693_0742_002.cdf')
      TES(3282)%FILENAME = TRIM('retv_vars.10693_0742_003.cdf')
      TES(3283)%FILENAME = TRIM('retv_vars.10693_0743_002.cdf')
      TES(3284)%FILENAME = TRIM('retv_vars.10693_0747_003.cdf')
      TES(3285)%FILENAME = TRIM('retv_vars.10693_0747_004.cdf')
      TES(3286)%FILENAME = TRIM('retv_vars.10693_0748_002.cdf')
      TES(3287)%FILENAME = TRIM('retv_vars.10693_0748_004.cdf')
      TES(3288)%FILENAME = TRIM('retv_vars.10693_0749_002.cdf')
      TES(3289)%FILENAME = TRIM('retv_vars.10693_0749_004.cdf')
      TES(3290)%FILENAME = TRIM('retv_vars.10695_0021_004.cdf')
      TES(3291)%FILENAME = TRIM('retv_vars.10695_0023_002.cdf')
      TES(3292)%FILENAME = TRIM('retv_vars.10695_0027_003.cdf')
      TES(3293)%FILENAME = TRIM('retv_vars.10695_0027_004.cdf')
      TES(3294)%FILENAME = TRIM('retv_vars.10695_0030_002.cdf')
      TES(3295)%FILENAME = TRIM('retv_vars.10695_0030_004.cdf')
      TES(3296)%FILENAME = TRIM('retv_vars.10695_0031_003.cdf')
      TES(3297)%FILENAME = TRIM('retv_vars.10695_0031_004.cdf')
      TES(3298)%FILENAME = TRIM('retv_vars.10695_0059_003.cdf')
      TES(3299)%FILENAME = TRIM('retv_vars.10695_0060_003.cdf')
      TES(3300)%FILENAME = TRIM('retv_vars.10695_0060_004.cdf')
      TES(3301)%FILENAME = TRIM('retv_vars.10695_0062_003.cdf')
      TES(3302)%FILENAME = TRIM('retv_vars.10695_0067_002.cdf')
      TES(3303)%FILENAME = TRIM('retv_vars.10695_0068_002.cdf')
      TES(3304)%FILENAME = TRIM('retv_vars.10695_0068_003.cdf')
      TES(3305)%FILENAME = TRIM('retv_vars.10695_0071_002.cdf')
      TES(3306)%FILENAME = TRIM('retv_vars.10695_0075_002.cdf')
      TES(3307)%FILENAME = TRIM('retv_vars.10695_0075_004.cdf')
      TES(3308)%FILENAME = TRIM('retv_vars.10695_0100_004.cdf')
      TES(3309)%FILENAME = TRIM('retv_vars.10695_0101_002.cdf')
      TES(3310)%FILENAME = TRIM('retv_vars.10695_0101_003.cdf')
      TES(3311)%FILENAME = TRIM('retv_vars.10695_0101_004.cdf')
      TES(3312)%FILENAME = TRIM('retv_vars.10695_0102_003.cdf')
      TES(3313)%FILENAME = TRIM('retv_vars.10695_0102_004.cdf')
      TES(3314)%FILENAME = TRIM('retv_vars.10695_0103_004.cdf')
      TES(3315)%FILENAME = TRIM('retv_vars.10695_0104_002.cdf')
      TES(3316)%FILENAME = TRIM('retv_vars.10695_0104_003.cdf')
      TES(3317)%FILENAME = TRIM('retv_vars.10695_0104_004.cdf')
      TES(3318)%FILENAME = TRIM('retv_vars.10695_0105_003.cdf')
      TES(3319)%FILENAME = TRIM('retv_vars.10695_0106_002.cdf')
      TES(3320)%FILENAME = TRIM('retv_vars.10695_0106_003.cdf')
      TES(3321)%FILENAME = TRIM('retv_vars.10695_0107_003.cdf')
      TES(3322)%FILENAME = TRIM('retv_vars.10695_0107_004.cdf')
      TES(3323)%FILENAME = TRIM('retv_vars.10695_0110_004.cdf')
      TES(3324)%FILENAME = TRIM('retv_vars.10695_0112_003.cdf')
      TES(3325)%FILENAME = TRIM('retv_vars.10695_0114_004.cdf')
      TES(3326)%FILENAME = TRIM('retv_vars.10695_0116_002.cdf')
      TES(3327)%FILENAME = TRIM('retv_vars.10695_0116_003.cdf')
      TES(3328)%FILENAME = TRIM('retv_vars.10695_0116_004.cdf')
      TES(3329)%FILENAME = TRIM('retv_vars.10695_0117_002.cdf')
      TES(3330)%FILENAME = TRIM('retv_vars.10695_0117_003.cdf')
      TES(3331)%FILENAME = TRIM('retv_vars.10695_0117_004.cdf')
      TES(3332)%FILENAME = TRIM('retv_vars.10695_0156_004.cdf')
      TES(3333)%FILENAME = TRIM('retv_vars.10695_0157_002.cdf')
      TES(3334)%FILENAME = TRIM('retv_vars.10695_0157_003.cdf')
      TES(3335)%FILENAME = TRIM('retv_vars.10695_0158_003.cdf')
      TES(3336)%FILENAME = TRIM('retv_vars.10695_0181_002.cdf')
      TES(3337)%FILENAME = TRIM('retv_vars.10695_0234_003.cdf')
      TES(3338)%FILENAME = TRIM('retv_vars.10695_0234_004.cdf')
      TES(3339)%FILENAME = TRIM('retv_vars.10695_0237_003.cdf')
      TES(3340)%FILENAME = TRIM('retv_vars.10695_0247_002.cdf')
      TES(3341)%FILENAME = TRIM('retv_vars.10695_0247_004.cdf')
      TES(3342)%FILENAME = TRIM('retv_vars.10695_0248_003.cdf')
      TES(3343)%FILENAME = TRIM('retv_vars.10695_0249_003.cdf')
      TES(3344)%FILENAME = TRIM('retv_vars.10695_0249_004.cdf')
      TES(3345)%FILENAME = TRIM('retv_vars.10695_0250_003.cdf')
      TES(3346)%FILENAME = TRIM('retv_vars.10695_0251_002.cdf')
      TES(3347)%FILENAME = TRIM('retv_vars.10695_0260_002.cdf')
      TES(3348)%FILENAME = TRIM('retv_vars.10695_0262_004.cdf')
      TES(3349)%FILENAME = TRIM('retv_vars.10695_0263_002.cdf')
      TES(3350)%FILENAME = TRIM('retv_vars.10695_0267_004.cdf')
      TES(3351)%FILENAME = TRIM('retv_vars.10695_0268_002.cdf')
      TES(3352)%FILENAME = TRIM('retv_vars.10695_0268_003.cdf')
      TES(3353)%FILENAME = TRIM('retv_vars.10695_0269_002.cdf')
      TES(3354)%FILENAME = TRIM('retv_vars.10695_0269_003.cdf')
      TES(3355)%FILENAME = TRIM('retv_vars.10695_0269_004.cdf')
      TES(3356)%FILENAME = TRIM('retv_vars.10695_0270_003.cdf')
      TES(3357)%FILENAME = TRIM('retv_vars.10695_0270_004.cdf')
      TES(3358)%FILENAME = TRIM('retv_vars.10695_0271_002.cdf')
      TES(3359)%FILENAME = TRIM('retv_vars.10695_0279_004.cdf')
      TES(3360)%FILENAME = TRIM('retv_vars.10695_0280_004.cdf')
      TES(3361)%FILENAME = TRIM('retv_vars.10695_0290_002.cdf')
      TES(3362)%FILENAME = TRIM('retv_vars.10695_0290_003.cdf')
      TES(3363)%FILENAME = TRIM('retv_vars.10695_0290_004.cdf')
      TES(3364)%FILENAME = TRIM('retv_vars.10695_0291_002.cdf')
      TES(3365)%FILENAME = TRIM('retv_vars.10695_0291_004.cdf')
      TES(3366)%FILENAME = TRIM('retv_vars.10695_0292_002.cdf')
      TES(3367)%FILENAME = TRIM('retv_vars.10695_0293_002.cdf')
      TES(3368)%FILENAME = TRIM('retv_vars.10695_0296_003.cdf')
      TES(3369)%FILENAME = TRIM('retv_vars.10695_0296_004.cdf')
      TES(3370)%FILENAME = TRIM('retv_vars.10695_0297_004.cdf')
      TES(3371)%FILENAME = TRIM('retv_vars.10695_0298_002.cdf')
      TES(3372)%FILENAME = TRIM('retv_vars.10695_0298_003.cdf')
      TES(3373)%FILENAME = TRIM('retv_vars.10695_0298_004.cdf')
      TES(3374)%FILENAME = TRIM('retv_vars.10695_0306_002.cdf')
      TES(3375)%FILENAME = TRIM('retv_vars.10695_0306_003.cdf')
      TES(3376)%FILENAME = TRIM('retv_vars.10695_0306_004.cdf')
      TES(3377)%FILENAME = TRIM('retv_vars.10695_0307_002.cdf')
      TES(3378)%FILENAME = TRIM('retv_vars.10695_0307_004.cdf')
      TES(3379)%FILENAME = TRIM('retv_vars.10695_0310_002.cdf')
      TES(3380)%FILENAME = TRIM('retv_vars.10695_0310_004.cdf')
      TES(3381)%FILENAME = TRIM('retv_vars.10695_0315_003.cdf')
      TES(3382)%FILENAME = TRIM('retv_vars.10695_0315_004.cdf')
      TES(3383)%FILENAME = TRIM('retv_vars.10695_0316_003.cdf')
      TES(3384)%FILENAME = TRIM('retv_vars.10695_0316_004.cdf')
      TES(3385)%FILENAME = TRIM('retv_vars.10695_0317_003.cdf')
      TES(3386)%FILENAME = TRIM('retv_vars.10695_0317_004.cdf')
      TES(3387)%FILENAME = TRIM('retv_vars.10695_0318_002.cdf')
      TES(3388)%FILENAME = TRIM('retv_vars.10695_0319_003.cdf')
      TES(3389)%FILENAME = TRIM('retv_vars.10695_0319_004.cdf')
      TES(3390)%FILENAME = TRIM('retv_vars.10695_0320_002.cdf')
      TES(3391)%FILENAME = TRIM('retv_vars.10695_0320_003.cdf')
      TES(3392)%FILENAME = TRIM('retv_vars.10695_0325_002.cdf')
      TES(3393)%FILENAME = TRIM('retv_vars.10695_0353_003.cdf')
      TES(3394)%FILENAME = TRIM('retv_vars.10695_0354_003.cdf')
      TES(3395)%FILENAME = TRIM('retv_vars.10695_0356_003.cdf')
      TES(3396)%FILENAME = TRIM('retv_vars.10695_0356_004.cdf')
      TES(3397)%FILENAME = TRIM('retv_vars.10695_0357_002.cdf')
      TES(3398)%FILENAME = TRIM('retv_vars.10695_0357_003.cdf')
      TES(3399)%FILENAME = TRIM('retv_vars.10695_0357_004.cdf')
      TES(3400)%FILENAME = TRIM('retv_vars.10695_0358_002.cdf')
      TES(3401)%FILENAME = TRIM('retv_vars.10695_0358_003.cdf')
      TES(3402)%FILENAME = TRIM('retv_vars.10695_0359_002.cdf')
      TES(3403)%FILENAME = TRIM('retv_vars.10695_0363_003.cdf')
      TES(3404)%FILENAME = TRIM('retv_vars.10695_0363_004.cdf')
      TES(3405)%FILENAME = TRIM('retv_vars.10695_0364_002.cdf')
      TES(3406)%FILENAME = TRIM('retv_vars.10695_0364_003.cdf')
      TES(3407)%FILENAME = TRIM('retv_vars.10695_0365_002.cdf')
      TES(3408)%FILENAME = TRIM('retv_vars.10695_0365_003.cdf')
      TES(3409)%FILENAME = TRIM('retv_vars.10695_0365_004.cdf')
      TES(3410)%FILENAME = TRIM('retv_vars.10695_0366_002.cdf')
      TES(3411)%FILENAME = TRIM('retv_vars.10695_0366_003.cdf')
      TES(3412)%FILENAME = TRIM('retv_vars.10695_0367_002.cdf')
      TES(3413)%FILENAME = TRIM('retv_vars.10695_0367_003.cdf')
      TES(3414)%FILENAME = TRIM('retv_vars.10695_0368_004.cdf')
      TES(3415)%FILENAME = TRIM('retv_vars.10695_0369_002.cdf')
      TES(3416)%FILENAME = TRIM('retv_vars.10695_0369_004.cdf')
      TES(3417)%FILENAME = TRIM('retv_vars.10695_0411_002.cdf')
      TES(3418)%FILENAME = TRIM('retv_vars.10695_0411_003.cdf')
      TES(3419)%FILENAME = TRIM('retv_vars.10695_0411_004.cdf')
      TES(3420)%FILENAME = TRIM('retv_vars.10695_0412_002.cdf')
      TES(3421)%FILENAME = TRIM('retv_vars.10695_0412_003.cdf')
      TES(3422)%FILENAME = TRIM('retv_vars.10695_0415_004.cdf')
      TES(3423)%FILENAME = TRIM('retv_vars.10695_0419_002.cdf')
      TES(3424)%FILENAME = TRIM('retv_vars.10695_0419_003.cdf')
      TES(3425)%FILENAME = TRIM('retv_vars.10695_0419_004.cdf')
      TES(3426)%FILENAME = TRIM('retv_vars.10695_0420_004.cdf')
      TES(3427)%FILENAME = TRIM('retv_vars.10695_0421_002.cdf')
      TES(3428)%FILENAME = TRIM('retv_vars.10695_0421_003.cdf')
      TES(3429)%FILENAME = TRIM('retv_vars.10695_0421_004.cdf')
      TES(3430)%FILENAME = TRIM('retv_vars.10695_0422_002.cdf')
      TES(3431)%FILENAME = TRIM('retv_vars.10695_0422_004.cdf')
      TES(3432)%FILENAME = TRIM('retv_vars.10695_0423_002.cdf')
      TES(3433)%FILENAME = TRIM('retv_vars.10695_0423_003.cdf')
      TES(3434)%FILENAME = TRIM('retv_vars.10695_0423_004.cdf')
      TES(3435)%FILENAME = TRIM('retv_vars.10695_0424_003.cdf')
      TES(3436)%FILENAME = TRIM('retv_vars.10695_0424_004.cdf')
      TES(3437)%FILENAME = TRIM('retv_vars.10695_0425_002.cdf')
      TES(3438)%FILENAME = TRIM('retv_vars.10695_0425_003.cdf')
      TES(3439)%FILENAME = TRIM('retv_vars.10695_0425_004.cdf')
      TES(3440)%FILENAME = TRIM('retv_vars.10695_0426_002.cdf')
      TES(3441)%FILENAME = TRIM('retv_vars.10695_0426_004.cdf')
      TES(3442)%FILENAME = TRIM('retv_vars.10695_0427_004.cdf')
      TES(3443)%FILENAME = TRIM('retv_vars.10695_0460_003.cdf')
      TES(3444)%FILENAME = TRIM('retv_vars.10695_0460_004.cdf')
      TES(3445)%FILENAME = TRIM('retv_vars.10695_0461_002.cdf')
      TES(3446)%FILENAME = TRIM('retv_vars.10695_0461_003.cdf')
      TES(3447)%FILENAME = TRIM('retv_vars.10695_0462_003.cdf')
      TES(3448)%FILENAME = TRIM('retv_vars.10695_0462_004.cdf')
      TES(3449)%FILENAME = TRIM('retv_vars.10695_0464_002.cdf')
      TES(3450)%FILENAME = TRIM('retv_vars.10695_0465_003.cdf')
      TES(3451)%FILENAME = TRIM('retv_vars.10695_0466_002.cdf')
      TES(3452)%FILENAME = TRIM('retv_vars.10695_0466_004.cdf')
      TES(3453)%FILENAME = TRIM('retv_vars.10695_0467_003.cdf')
      TES(3454)%FILENAME = TRIM('retv_vars.10695_0469_002.cdf')
      TES(3455)%FILENAME = TRIM('retv_vars.10695_0469_003.cdf')
      TES(3456)%FILENAME = TRIM('retv_vars.10695_0470_004.cdf')
      TES(3457)%FILENAME = TRIM('retv_vars.10695_0471_002.cdf')
      TES(3458)%FILENAME = TRIM('retv_vars.10695_0471_004.cdf')
      TES(3459)%FILENAME = TRIM('retv_vars.10695_0472_002.cdf')
      TES(3460)%FILENAME = TRIM('retv_vars.10695_0472_003.cdf')
      TES(3461)%FILENAME = TRIM('retv_vars.10695_0472_004.cdf')
      TES(3462)%FILENAME = TRIM('retv_vars.10695_0474_002.cdf')
      TES(3463)%FILENAME = TRIM('retv_vars.10695_0474_003.cdf')
      TES(3464)%FILENAME = TRIM('retv_vars.10695_0507_002.cdf')
      TES(3465)%FILENAME = TRIM('retv_vars.10695_0508_004.cdf')
      TES(3466)%FILENAME = TRIM('retv_vars.10695_0509_002.cdf')
      TES(3467)%FILENAME = TRIM('retv_vars.10695_0509_004.cdf')
      TES(3468)%FILENAME = TRIM('retv_vars.10695_0510_002.cdf')
      TES(3469)%FILENAME = TRIM('retv_vars.10695_0510_003.cdf')
      TES(3470)%FILENAME = TRIM('retv_vars.10695_0510_004.cdf')
      TES(3471)%FILENAME = TRIM('retv_vars.10695_0513_004.cdf')
      TES(3472)%FILENAME = TRIM('retv_vars.10695_0514_003.cdf')
      TES(3473)%FILENAME = TRIM('retv_vars.10695_0514_004.cdf')
      TES(3474)%FILENAME = TRIM('retv_vars.10695_0516_002.cdf')
      TES(3475)%FILENAME = TRIM('retv_vars.10695_0517_004.cdf')
      TES(3476)%FILENAME = TRIM('retv_vars.10695_0530_002.cdf')
      TES(3477)%FILENAME = TRIM('retv_vars.10695_0535_003.cdf')
      TES(3478)%FILENAME = TRIM('retv_vars.10695_0568_004.cdf')
      TES(3479)%FILENAME = TRIM('retv_vars.10695_0569_002.cdf')
      TES(3480)%FILENAME = TRIM('retv_vars.10695_0579_003.cdf')
      TES(3481)%FILENAME = TRIM('retv_vars.10695_0580_003.cdf')
      TES(3482)%FILENAME = TRIM('retv_vars.10695_0582_003.cdf')
      TES(3483)%FILENAME = TRIM('retv_vars.10695_0583_002.cdf')
      TES(3484)%FILENAME = TRIM('retv_vars.10695_0583_003.cdf')
      TES(3485)%FILENAME = TRIM('retv_vars.10695_0584_003.cdf')
      TES(3486)%FILENAME = TRIM('retv_vars.10695_0586_003.cdf')
      TES(3487)%FILENAME = TRIM('retv_vars.10695_0595_002.cdf')
      TES(3488)%FILENAME = TRIM('retv_vars.10695_0595_003.cdf')
      TES(3489)%FILENAME = TRIM('retv_vars.10695_0595_004.cdf')
      TES(3490)%FILENAME = TRIM('retv_vars.10695_0596_002.cdf')
      TES(3491)%FILENAME = TRIM('retv_vars.10695_0596_004.cdf')
      TES(3492)%FILENAME = TRIM('retv_vars.10695_0597_002.cdf')
      TES(3493)%FILENAME = TRIM('retv_vars.10695_0597_003.cdf')
      TES(3494)%FILENAME = TRIM('retv_vars.10695_0597_004.cdf')
      TES(3495)%FILENAME = TRIM('retv_vars.10695_0614_004.cdf')
      TES(3496)%FILENAME = TRIM('retv_vars.10695_0615_002.cdf')
      TES(3497)%FILENAME = TRIM('retv_vars.10695_0615_003.cdf')
      TES(3498)%FILENAME = TRIM('retv_vars.10695_0615_004.cdf')
      TES(3499)%FILENAME = TRIM('retv_vars.10695_0616_002.cdf')
      TES(3500)%FILENAME = TRIM('retv_vars.10695_0616_004.cdf')
      TES(3501)%FILENAME = TRIM('retv_vars.10695_0617_002.cdf')
      TES(3502)%FILENAME = TRIM('retv_vars.10695_0617_003.cdf')
      TES(3503)%FILENAME = TRIM('retv_vars.10695_0617_004.cdf')
      TES(3504)%FILENAME = TRIM('retv_vars.10695_0618_002.cdf')
      TES(3505)%FILENAME = TRIM('retv_vars.10695_0618_003.cdf')
      TES(3506)%FILENAME = TRIM('retv_vars.10695_0621_002.cdf')
      TES(3507)%FILENAME = TRIM('retv_vars.10695_0621_003.cdf')
      TES(3508)%FILENAME = TRIM('retv_vars.10695_0622_002.cdf')
      TES(3509)%FILENAME = TRIM('retv_vars.10695_0622_003.cdf')
      TES(3510)%FILENAME = TRIM('retv_vars.10695_0622_004.cdf')
      TES(3511)%FILENAME = TRIM('retv_vars.10695_0623_002.cdf')
      TES(3512)%FILENAME = TRIM('retv_vars.10695_0623_003.cdf')
      TES(3513)%FILENAME = TRIM('retv_vars.10695_0623_004.cdf')
      TES(3514)%FILENAME = TRIM('retv_vars.10695_0624_002.cdf')
      TES(3515)%FILENAME = TRIM('retv_vars.10695_0624_003.cdf')
      TES(3516)%FILENAME = TRIM('retv_vars.10695_0628_002.cdf')
      TES(3517)%FILENAME = TRIM('retv_vars.10695_0628_004.cdf')
      TES(3518)%FILENAME = TRIM('retv_vars.10695_0629_002.cdf')
      TES(3519)%FILENAME = TRIM('retv_vars.10695_0633_003.cdf')
      TES(3520)%FILENAME = TRIM('retv_vars.10695_0634_004.cdf')
      TES(3521)%FILENAME = TRIM('retv_vars.10695_0635_002.cdf')
      TES(3522)%FILENAME = TRIM('retv_vars.10695_0640_002.cdf')
      TES(3523)%FILENAME = TRIM('retv_vars.10695_0640_003.cdf')
      TES(3524)%FILENAME = TRIM('retv_vars.10695_0640_004.cdf')
      TES(3525)%FILENAME = TRIM('retv_vars.10695_0641_004.cdf')
      TES(3526)%FILENAME = TRIM('retv_vars.10695_0642_003.cdf')
      TES(3527)%FILENAME = TRIM('retv_vars.10695_0644_002.cdf')
      TES(3528)%FILENAME = TRIM('retv_vars.10695_0644_004.cdf')
      TES(3529)%FILENAME = TRIM('retv_vars.10695_0645_003.cdf')
      TES(3530)%FILENAME = TRIM('retv_vars.10695_0646_003.cdf')
      TES(3531)%FILENAME = TRIM('retv_vars.10695_0646_004.cdf')
      TES(3532)%FILENAME = TRIM('retv_vars.10695_0647_002.cdf')
      TES(3533)%FILENAME = TRIM('retv_vars.10695_0651_003.cdf')
      TES(3534)%FILENAME = TRIM('retv_vars.10695_0652_002.cdf')
      TES(3535)%FILENAME = TRIM('retv_vars.10695_0652_004.cdf')
      TES(3536)%FILENAME = TRIM('retv_vars.10695_0653_002.cdf')
      TES(3537)%FILENAME = TRIM('retv_vars.10695_0653_003.cdf')
      TES(3538)%FILENAME = TRIM('retv_vars.10695_0654_002.cdf')
      TES(3539)%FILENAME = TRIM('retv_vars.10695_0654_003.cdf')
      TES(3540)%FILENAME = TRIM('retv_vars.10695_0654_004.cdf')
      TES(3541)%FILENAME = TRIM('retv_vars.10695_0658_004.cdf')
      TES(3542)%FILENAME = TRIM('retv_vars.10695_0660_003.cdf')
      TES(3543)%FILENAME = TRIM('retv_vars.10695_0686_003.cdf')
      TES(3544)%FILENAME = TRIM('retv_vars.10695_0687_003.cdf')
      TES(3545)%FILENAME = TRIM('retv_vars.10695_0691_003.cdf')
      TES(3546)%FILENAME = TRIM('retv_vars.10695_0692_002.cdf')
      TES(3547)%FILENAME = TRIM('retv_vars.10695_0692_003.cdf')
      TES(3548)%FILENAME = TRIM('retv_vars.10695_0693_002.cdf')
      TES(3549)%FILENAME = TRIM('retv_vars.10695_0693_004.cdf')
      TES(3550)%FILENAME = TRIM('retv_vars.10695_0694_002.cdf')
      TES(3551)%FILENAME = TRIM('retv_vars.10695_0694_003.cdf')
      TES(3552)%FILENAME = TRIM('retv_vars.10695_0694_004.cdf')
      TES(3553)%FILENAME = TRIM('retv_vars.10695_0699_004.cdf')
      TES(3554)%FILENAME = TRIM('retv_vars.10695_0700_002.cdf')
      TES(3555)%FILENAME = TRIM('retv_vars.10695_0701_003.cdf')
      TES(3556)%FILENAME = TRIM('retv_vars.10695_0701_004.cdf')
      TES(3557)%FILENAME = TRIM('retv_vars.10695_0702_002.cdf')
      TES(3558)%FILENAME = TRIM('retv_vars.10695_0702_003.cdf')
      TES(3559)%FILENAME = TRIM('retv_vars.10695_0703_002.cdf')
      TES(3560)%FILENAME = TRIM('retv_vars.10695_0703_003.cdf')
      TES(3561)%FILENAME = TRIM('retv_vars.10695_0704_002.cdf')
      TES(3562)%FILENAME = TRIM('retv_vars.10695_0704_003.cdf')
      TES(3563)%FILENAME = TRIM('retv_vars.10695_0705_002.cdf')
      TES(3564)%FILENAME = TRIM('retv_vars.10695_0705_003.cdf')
      TES(3565)%FILENAME = TRIM('retv_vars.10695_0706_003.cdf')
      TES(3566)%FILENAME = TRIM('retv_vars.10695_0706_004.cdf')
      TES(3567)%FILENAME = TRIM('retv_vars.10695_0707_002.cdf')
      TES(3568)%FILENAME = TRIM('retv_vars.10695_0743_002.cdf')
      TES(3569)%FILENAME = TRIM('retv_vars.10695_0743_003.cdf')
      TES(3570)%FILENAME = TRIM('retv_vars.10695_0747_003.cdf')
      TES(3571)%FILENAME = TRIM('retv_vars.10695_0747_004.cdf')
      TES(3572)%FILENAME = TRIM('retv_vars.10695_0748_002.cdf')
      TES(3573)%FILENAME = TRIM('retv_vars.10695_0748_003.cdf')
      TES(3574)%FILENAME = TRIM('retv_vars.10695_0748_004.cdf')
      TES(3575)%FILENAME = TRIM('retv_vars.10703_0012_003.cdf')
      TES(3576)%FILENAME = TRIM('retv_vars.10703_0013_002.cdf')
      TES(3577)%FILENAME = TRIM('retv_vars.10703_0013_004.cdf')
      TES(3578)%FILENAME = TRIM('retv_vars.10703_0019_002.cdf')
      TES(3579)%FILENAME = TRIM('retv_vars.10703_0020_002.cdf')
      TES(3580)%FILENAME = TRIM('retv_vars.10703_0021_004.cdf')
      TES(3581)%FILENAME = TRIM('retv_vars.10703_0022_004.cdf')
      TES(3582)%FILENAME = TRIM('retv_vars.10703_0027_004.cdf')
      TES(3583)%FILENAME = TRIM('retv_vars.10703_0028_002.cdf')
      TES(3584)%FILENAME = TRIM('retv_vars.10703_0053_002.cdf')
      TES(3585)%FILENAME = TRIM('retv_vars.10703_0053_003.cdf')
      TES(3586)%FILENAME = TRIM('retv_vars.10703_0054_002.cdf')
      TES(3587)%FILENAME = TRIM('retv_vars.10703_0054_003.cdf')
      TES(3588)%FILENAME = TRIM('retv_vars.10703_0054_004.cdf')
      TES(3589)%FILENAME = TRIM('retv_vars.10703_0055_002.cdf')
      TES(3590)%FILENAME = TRIM('retv_vars.10703_0055_004.cdf')
      TES(3591)%FILENAME = TRIM('retv_vars.10703_0056_003.cdf')
      TES(3592)%FILENAME = TRIM('retv_vars.10703_0057_003.cdf')
      TES(3593)%FILENAME = TRIM('retv_vars.10703_0058_002.cdf')
      TES(3594)%FILENAME = TRIM('retv_vars.10703_0058_003.cdf')
      TES(3595)%FILENAME = TRIM('retv_vars.10703_0059_003.cdf')
      TES(3596)%FILENAME = TRIM('retv_vars.10703_0059_004.cdf')
      TES(3597)%FILENAME = TRIM('retv_vars.10703_0062_003.cdf')
      TES(3598)%FILENAME = TRIM('retv_vars.10703_0067_003.cdf')
      TES(3599)%FILENAME = TRIM('retv_vars.10703_0067_004.cdf')
      TES(3600)%FILENAME = TRIM('retv_vars.10703_0068_003.cdf')
      TES(3601)%FILENAME = TRIM('retv_vars.10703_0069_002.cdf')
      TES(3602)%FILENAME = TRIM('retv_vars.10703_0069_004.cdf')
      TES(3603)%FILENAME = TRIM('retv_vars.10703_0070_002.cdf')
      TES(3604)%FILENAME = TRIM('retv_vars.10703_0108_003.cdf')
      TES(3605)%FILENAME = TRIM('retv_vars.10703_0109_003.cdf')
      TES(3606)%FILENAME = TRIM('retv_vars.10703_0110_003.cdf')
      TES(3607)%FILENAME = TRIM('retv_vars.10703_0133_004.cdf')
      TES(3608)%FILENAME = TRIM('retv_vars.10703_0172_002.cdf')
      TES(3609)%FILENAME = TRIM('retv_vars.10703_0172_003.cdf')
      TES(3610)%FILENAME = TRIM('retv_vars.10703_0172_004.cdf')
      TES(3611)%FILENAME = TRIM('retv_vars.10703_0173_002.cdf')
      TES(3612)%FILENAME = TRIM('retv_vars.10703_0186_004.cdf')
      TES(3613)%FILENAME = TRIM('retv_vars.10703_0187_002.cdf')
      TES(3614)%FILENAME = TRIM('retv_vars.10703_0187_003.cdf')
      TES(3615)%FILENAME = TRIM('retv_vars.10703_0187_004.cdf')
      TES(3616)%FILENAME = TRIM('retv_vars.10703_0190_002.cdf')
      TES(3617)%FILENAME = TRIM('retv_vars.10703_0198_003.cdf')
      TES(3618)%FILENAME = TRIM('retv_vars.10703_0198_004.cdf')
      TES(3619)%FILENAME = TRIM('retv_vars.10703_0199_002.cdf')
      TES(3620)%FILENAME = TRIM('retv_vars.10703_0199_003.cdf')
      TES(3621)%FILENAME = TRIM('retv_vars.10703_0200_004.cdf')
      TES(3622)%FILENAME = TRIM('retv_vars.10703_0201_004.cdf')
      TES(3623)%FILENAME = TRIM('retv_vars.10703_0202_002.cdf')
      TES(3624)%FILENAME = TRIM('retv_vars.10703_0202_004.cdf')
      TES(3625)%FILENAME = TRIM('retv_vars.10703_0203_002.cdf')
      TES(3626)%FILENAME = TRIM('retv_vars.10703_0212_002.cdf')
      TES(3627)%FILENAME = TRIM('retv_vars.10703_0214_004.cdf')
      TES(3628)%FILENAME = TRIM('retv_vars.10703_0220_002.cdf')
      TES(3629)%FILENAME = TRIM('retv_vars.10703_0220_003.cdf')
      TES(3630)%FILENAME = TRIM('retv_vars.10703_0220_004.cdf')
      TES(3631)%FILENAME = TRIM('retv_vars.10703_0221_002.cdf')
      TES(3632)%FILENAME = TRIM('retv_vars.10703_0221_004.cdf')
      TES(3633)%FILENAME = TRIM('retv_vars.10703_0222_003.cdf')
      TES(3634)%FILENAME = TRIM('retv_vars.10703_0222_004.cdf')
      TES(3635)%FILENAME = TRIM('retv_vars.10703_0223_002.cdf')
      TES(3636)%FILENAME = TRIM('retv_vars.10703_0223_003.cdf')
      TES(3637)%FILENAME = TRIM('retv_vars.10703_0227_003.cdf')
      TES(3638)%FILENAME = TRIM('retv_vars.10703_0228_002.cdf')
      TES(3639)%FILENAME = TRIM('retv_vars.10703_0232_004.cdf')
      TES(3640)%FILENAME = TRIM('retv_vars.10703_0242_003.cdf')
      TES(3641)%FILENAME = TRIM('retv_vars.10703_0242_004.cdf')
      TES(3642)%FILENAME = TRIM('retv_vars.10703_0243_002.cdf')
      TES(3643)%FILENAME = TRIM('retv_vars.10703_0243_003.cdf')
      TES(3644)%FILENAME = TRIM('retv_vars.10703_0244_002.cdf')
      TES(3645)%FILENAME = TRIM('retv_vars.10703_0244_003.cdf')
      TES(3646)%FILENAME = TRIM('retv_vars.10703_0244_004.cdf')
      TES(3647)%FILENAME = TRIM('retv_vars.10703_0245_004.cdf')
      TES(3648)%FILENAME = TRIM('retv_vars.10703_0246_002.cdf')
      TES(3649)%FILENAME = TRIM('retv_vars.10703_0247_004.cdf')
      TES(3650)%FILENAME = TRIM('retv_vars.10703_0249_002.cdf')
      TES(3651)%FILENAME = TRIM('retv_vars.10703_0249_004.cdf')
      TES(3652)%FILENAME = TRIM('retv_vars.10703_0250_004.cdf')
      TES(3653)%FILENAME = TRIM('retv_vars.10703_0251_004.cdf')
      TES(3654)%FILENAME = TRIM('retv_vars.10703_0255_003.cdf')
      TES(3655)%FILENAME = TRIM('retv_vars.10703_0256_004.cdf')
      TES(3656)%FILENAME = TRIM('retv_vars.10703_0257_002.cdf')
      TES(3657)%FILENAME = TRIM('retv_vars.10703_0257_003.cdf')
      TES(3658)%FILENAME = TRIM('retv_vars.10703_0258_004.cdf')
      TES(3659)%FILENAME = TRIM('retv_vars.10703_0259_002.cdf')
      TES(3660)%FILENAME = TRIM('retv_vars.10703_0259_004.cdf')
      TES(3661)%FILENAME = TRIM('retv_vars.10703_0260_002.cdf')
      TES(3662)%FILENAME = TRIM('retv_vars.10703_0260_003.cdf')
      TES(3663)%FILENAME = TRIM('retv_vars.10703_0261_002.cdf')
      TES(3664)%FILENAME = TRIM('retv_vars.10703_0261_003.cdf')
      TES(3665)%FILENAME = TRIM('retv_vars.10703_0262_002.cdf')
      TES(3666)%FILENAME = TRIM('retv_vars.10703_0262_003.cdf')
      TES(3667)%FILENAME = TRIM('retv_vars.10703_0267_003.cdf')
      TES(3668)%FILENAME = TRIM('retv_vars.10703_0269_002.cdf')
      TES(3669)%FILENAME = TRIM('retv_vars.10703_0274_002.cdf')
      TES(3670)%FILENAME = TRIM('retv_vars.10703_0274_003.cdf')
      TES(3671)%FILENAME = TRIM('retv_vars.10703_0275_004.cdf')
      TES(3672)%FILENAME = TRIM('retv_vars.10703_0303_002.cdf')
      TES(3673)%FILENAME = TRIM('retv_vars.10703_0303_003.cdf')
      TES(3674)%FILENAME = TRIM('retv_vars.10703_0303_004.cdf')
      TES(3675)%FILENAME = TRIM('retv_vars.10703_0307_003.cdf')
      TES(3676)%FILENAME = TRIM('retv_vars.10703_0307_004.cdf')
      TES(3677)%FILENAME = TRIM('retv_vars.10703_0308_002.cdf')
      TES(3678)%FILENAME = TRIM('retv_vars.10703_0309_002.cdf')
      TES(3679)%FILENAME = TRIM('retv_vars.10703_0309_004.cdf')
      TES(3680)%FILENAME = TRIM('retv_vars.10703_0310_002.cdf')
      TES(3681)%FILENAME = TRIM('retv_vars.10703_0310_003.cdf')
      TES(3682)%FILENAME = TRIM('retv_vars.10703_0310_004.cdf')
      TES(3683)%FILENAME = TRIM('retv_vars.10703_0315_003.cdf')
      TES(3684)%FILENAME = TRIM('retv_vars.10703_0316_002.cdf')
      TES(3685)%FILENAME = TRIM('retv_vars.10703_0316_003.cdf')
      TES(3686)%FILENAME = TRIM('retv_vars.10703_0316_004.cdf')
      TES(3687)%FILENAME = TRIM('retv_vars.10703_0317_003.cdf')
      TES(3688)%FILENAME = TRIM('retv_vars.10703_0318_002.cdf')
      TES(3689)%FILENAME = TRIM('retv_vars.10703_0318_004.cdf')
      TES(3690)%FILENAME = TRIM('retv_vars.10703_0319_004.cdf')
      TES(3691)%FILENAME = TRIM('retv_vars.10703_0320_003.cdf')
      TES(3692)%FILENAME = TRIM('retv_vars.10703_0322_003.cdf')
      TES(3693)%FILENAME = TRIM('retv_vars.10703_0363_003.cdf')
      TES(3694)%FILENAME = TRIM('retv_vars.10703_0363_004.cdf')
      TES(3695)%FILENAME = TRIM('retv_vars.10703_0364_004.cdf')
      TES(3696)%FILENAME = TRIM('retv_vars.10703_0365_002.cdf')
      TES(3697)%FILENAME = TRIM('retv_vars.10703_0365_004.cdf')
      TES(3698)%FILENAME = TRIM('retv_vars.10703_0369_002.cdf')
      TES(3699)%FILENAME = TRIM('retv_vars.10703_0369_003.cdf')
      TES(3700)%FILENAME = TRIM('retv_vars.10703_0371_003.cdf')
      TES(3701)%FILENAME = TRIM('retv_vars.10703_0372_003.cdf')
      TES(3702)%FILENAME = TRIM('retv_vars.10703_0374_004.cdf')
      TES(3703)%FILENAME = TRIM('retv_vars.10703_0378_003.cdf')
      TES(3704)%FILENAME = TRIM('retv_vars.10703_0406_002.cdf')
      TES(3705)%FILENAME = TRIM('retv_vars.10703_0411_002.cdf')
      TES(3706)%FILENAME = TRIM('retv_vars.10703_0411_003.cdf')
      TES(3707)%FILENAME = TRIM('retv_vars.10703_0411_004.cdf')
      TES(3708)%FILENAME = TRIM('retv_vars.10703_0414_004.cdf')
      TES(3709)%FILENAME = TRIM('retv_vars.10703_0415_002.cdf')
      TES(3710)%FILENAME = TRIM('retv_vars.10703_0415_003.cdf')
      TES(3711)%FILENAME = TRIM('retv_vars.10703_0417_003.cdf')
      TES(3712)%FILENAME = TRIM('retv_vars.10703_0417_004.cdf')
      TES(3713)%FILENAME = TRIM('retv_vars.10703_0418_003.cdf')
      TES(3714)%FILENAME = TRIM('retv_vars.10703_0422_003.cdf')
      TES(3715)%FILENAME = TRIM('retv_vars.10703_0423_003.cdf')
      TES(3716)%FILENAME = TRIM('retv_vars.10703_0423_004.cdf')
      TES(3717)%FILENAME = TRIM('retv_vars.10703_0424_002.cdf')
      TES(3718)%FILENAME = TRIM('retv_vars.10703_0424_003.cdf')
      TES(3719)%FILENAME = TRIM('retv_vars.10703_0425_002.cdf')
      TES(3720)%FILENAME = TRIM('retv_vars.10703_0425_004.cdf')
      TES(3721)%FILENAME = TRIM('retv_vars.10703_0427_002.cdf')
      TES(3722)%FILENAME = TRIM('retv_vars.10703_0427_003.cdf')
      TES(3723)%FILENAME = TRIM('retv_vars.10703_0460_002.cdf')
      TES(3724)%FILENAME = TRIM('retv_vars.10703_0461_002.cdf')
      TES(3725)%FILENAME = TRIM('retv_vars.10703_0461_004.cdf')
      TES(3726)%FILENAME = TRIM('retv_vars.10703_0462_002.cdf')
      TES(3727)%FILENAME = TRIM('retv_vars.10703_0464_004.cdf')
      TES(3728)%FILENAME = TRIM('retv_vars.10703_0466_002.cdf')
      TES(3729)%FILENAME = TRIM('retv_vars.10703_0467_002.cdf')
      TES(3730)%FILENAME = TRIM('retv_vars.10703_0467_003.cdf')
      TES(3731)%FILENAME = TRIM('retv_vars.10703_0468_003.cdf')
      TES(3732)%FILENAME = TRIM('retv_vars.10703_0469_002.cdf')
      TES(3733)%FILENAME = TRIM('retv_vars.10703_0469_003.cdf')
      TES(3734)%FILENAME = TRIM('retv_vars.10703_0469_004.cdf')
      TES(3735)%FILENAME = TRIM('retv_vars.10703_0470_002.cdf')
      TES(3736)%FILENAME = TRIM('retv_vars.10703_0482_003.cdf')
      TES(3737)%FILENAME = TRIM('retv_vars.10703_0500_003.cdf')
      TES(3738)%FILENAME = TRIM('retv_vars.10703_0530_004.cdf')
      TES(3739)%FILENAME = TRIM('retv_vars.10703_0534_003.cdf')
      TES(3740)%FILENAME = TRIM('retv_vars.10703_0534_004.cdf')
      TES(3741)%FILENAME = TRIM('retv_vars.10703_0535_002.cdf')
      TES(3742)%FILENAME = TRIM('retv_vars.10703_0535_003.cdf')
      TES(3743)%FILENAME = TRIM('retv_vars.10703_0538_004.cdf')
      TES(3744)%FILENAME = TRIM('retv_vars.10703_0546_002.cdf')
      TES(3745)%FILENAME = TRIM('retv_vars.10703_0546_004.cdf')
      TES(3746)%FILENAME = TRIM('retv_vars.10703_0547_003.cdf')
      TES(3747)%FILENAME = TRIM('retv_vars.10703_0547_004.cdf')
      TES(3748)%FILENAME = TRIM('retv_vars.10703_0548_003.cdf')
      TES(3749)%FILENAME = TRIM('retv_vars.10703_0548_004.cdf')
      TES(3750)%FILENAME = TRIM('retv_vars.10703_0549_003.cdf')
      TES(3751)%FILENAME = TRIM('retv_vars.10703_0550_003.cdf')
      TES(3752)%FILENAME = TRIM('retv_vars.10703_0550_004.cdf')
      TES(3753)%FILENAME = TRIM('retv_vars.10703_0567_004.cdf')
      TES(3754)%FILENAME = TRIM('retv_vars.10703_0568_002.cdf')
      TES(3755)%FILENAME = TRIM('retv_vars.10703_0569_004.cdf')
      TES(3756)%FILENAME = TRIM('retv_vars.10703_0570_002.cdf')
      TES(3757)%FILENAME = TRIM('retv_vars.10703_0570_003.cdf')
      TES(3758)%FILENAME = TRIM('retv_vars.10703_0571_003.cdf')
      TES(3759)%FILENAME = TRIM('retv_vars.10703_0571_004.cdf')
      TES(3760)%FILENAME = TRIM('retv_vars.10703_0572_002.cdf')
      TES(3761)%FILENAME = TRIM('retv_vars.10703_0572_003.cdf')
      TES(3762)%FILENAME = TRIM('retv_vars.10703_0572_004.cdf')
      TES(3763)%FILENAME = TRIM('retv_vars.10703_0573_002.cdf')
      TES(3764)%FILENAME = TRIM('retv_vars.10703_0574_002.cdf')
      TES(3765)%FILENAME = TRIM('retv_vars.10703_0574_003.cdf')
      TES(3766)%FILENAME = TRIM('retv_vars.10703_0574_004.cdf')
      TES(3767)%FILENAME = TRIM('retv_vars.10703_0575_002.cdf')
      TES(3768)%FILENAME = TRIM('retv_vars.10703_0575_003.cdf')
      TES(3769)%FILENAME = TRIM('retv_vars.10703_0575_004.cdf')
      TES(3770)%FILENAME = TRIM('retv_vars.10703_0576_002.cdf')
      TES(3771)%FILENAME = TRIM('retv_vars.10703_0576_003.cdf')
      TES(3772)%FILENAME = TRIM('retv_vars.10703_0580_002.cdf')
      TES(3773)%FILENAME = TRIM('retv_vars.10703_0580_004.cdf')
      TES(3774)%FILENAME = TRIM('retv_vars.10703_0582_004.cdf')
      TES(3775)%FILENAME = TRIM('retv_vars.10703_0585_003.cdf')
      TES(3776)%FILENAME = TRIM('retv_vars.10703_0586_003.cdf')
      TES(3777)%FILENAME = TRIM('retv_vars.10703_0586_004.cdf')
      TES(3778)%FILENAME = TRIM('retv_vars.10703_0587_002.cdf')
      TES(3779)%FILENAME = TRIM('retv_vars.10703_0587_003.cdf')
      TES(3780)%FILENAME = TRIM('retv_vars.10703_0591_003.cdf')
      TES(3781)%FILENAME = TRIM('retv_vars.10703_0591_004.cdf')
      TES(3782)%FILENAME = TRIM('retv_vars.10703_0592_004.cdf')
      TES(3783)%FILENAME = TRIM('retv_vars.10703_0595_002.cdf')
      TES(3784)%FILENAME = TRIM('retv_vars.10703_0598_003.cdf')
      TES(3785)%FILENAME = TRIM('retv_vars.10703_0598_004.cdf')
      TES(3786)%FILENAME = TRIM('retv_vars.10703_0599_002.cdf')
      TES(3787)%FILENAME = TRIM('retv_vars.10703_0603_003.cdf')
      TES(3788)%FILENAME = TRIM('retv_vars.10703_0603_004.cdf')
      TES(3789)%FILENAME = TRIM('retv_vars.10703_0610_003.cdf')
      TES(3790)%FILENAME = TRIM('retv_vars.10703_0613_003.cdf')
      TES(3791)%FILENAME = TRIM('retv_vars.10703_0615_002.cdf')
      TES(3792)%FILENAME = TRIM('retv_vars.10703_0615_003.cdf')
      TES(3793)%FILENAME = TRIM('retv_vars.10703_0640_003.cdf')
      TES(3794)%FILENAME = TRIM('retv_vars.10703_0643_003.cdf')
      TES(3795)%FILENAME = TRIM('retv_vars.10703_0643_004.cdf')
      TES(3796)%FILENAME = TRIM('retv_vars.10703_0644_004.cdf')
      TES(3797)%FILENAME = TRIM('retv_vars.10703_0645_002.cdf')
      TES(3798)%FILENAME = TRIM('retv_vars.10703_0645_004.cdf')
      TES(3799)%FILENAME = TRIM('retv_vars.10703_0646_003.cdf')
      TES(3800)%FILENAME = TRIM('retv_vars.10703_0647_002.cdf')
      TES(3801)%FILENAME = TRIM('retv_vars.10703_0652_002.cdf')
      TES(3802)%FILENAME = TRIM('retv_vars.10703_0652_003.cdf')
      TES(3803)%FILENAME = TRIM('retv_vars.10703_0652_004.cdf')
      TES(3804)%FILENAME = TRIM('retv_vars.10703_0653_003.cdf')
      TES(3805)%FILENAME = TRIM('retv_vars.10703_0654_003.cdf')
      TES(3806)%FILENAME = TRIM('retv_vars.10703_0656_004.cdf')
      TES(3807)%FILENAME = TRIM('retv_vars.10703_0657_002.cdf')
      TES(3808)%FILENAME = TRIM('retv_vars.10703_0657_004.cdf')
      TES(3809)%FILENAME = TRIM('retv_vars.10703_0658_002.cdf')
      TES(3810)%FILENAME = TRIM('retv_vars.10703_0658_003.cdf')
      TES(3811)%FILENAME = TRIM('retv_vars.10703_0658_004.cdf')
      TES(3812)%FILENAME = TRIM('retv_vars.10703_0659_002.cdf')
      TES(3813)%FILENAME = TRIM('retv_vars.10703_0692_002.cdf')
      TES(3814)%FILENAME = TRIM('retv_vars.10703_0693_003.cdf')
      TES(3815)%FILENAME = TRIM('retv_vars.10703_0694_003.cdf')
      TES(3816)%FILENAME = TRIM('retv_vars.10703_0694_004.cdf')
      TES(3817)%FILENAME = TRIM('retv_vars.10703_0699_002.cdf')
      TES(3818)%FILENAME = TRIM('retv_vars.10703_0700_003.cdf')
      TES(3819)%FILENAME = TRIM('retv_vars.10703_0700_004.cdf')
      TES(3820)%FILENAME = TRIM('retv_vars.10703_0701_002.cdf')
      TES(3821)%FILENAME = TRIM('retv_vars.10703_0701_004.cdf')
      TES(3822)%FILENAME = TRIM('retv_vars.10703_0702_002.cdf')
      TES(3823)%FILENAME = TRIM('retv_vars.10703_0702_003.cdf')
      TES(3824)%FILENAME = TRIM('retv_vars.10703_0702_004.cdf')
      TES(3825)%FILENAME = TRIM('retv_vars.10703_0730_004.cdf')
      TES(3826)%FILENAME = TRIM('retv_vars.10703_0731_002.cdf')
      TES(3827)%FILENAME = TRIM('retv_vars.10703_0731_003.cdf')
      TES(3828)%FILENAME = TRIM('retv_vars.10703_0731_004.cdf')
      TES(3829)%FILENAME = TRIM('retv_vars.10703_0732_003.cdf')
      TES(3830)%FILENAME = TRIM('retv_vars.10703_0732_004.cdf')
      TES(3831)%FILENAME = TRIM('retv_vars.10703_0733_002.cdf')
      TES(3832)%FILENAME = TRIM('retv_vars.10703_0733_004.cdf')
      TES(3833)%FILENAME = TRIM('retv_vars.10703_0740_002.cdf')
      TES(3834)%FILENAME = TRIM('retv_vars.10703_0740_003.cdf')
      TES(3835)%FILENAME = TRIM('retv_vars.10703_0741_003.cdf')
      TES(3836)%FILENAME = TRIM('retv_vars.10703_0741_004.cdf')
      TES(3837)%FILENAME = TRIM('retv_vars.10703_0742_002.cdf')
      TES(3838)%FILENAME = TRIM('retv_vars.10703_0742_003.cdf')
      TES(3839)%FILENAME = TRIM('retv_vars.10703_0742_004.cdf')
      TES(3840)%FILENAME = TRIM('retv_vars.10703_0747_003.cdf')
      TES(3841)%FILENAME = TRIM('retv_vars.10708_0013_002.cdf')
      TES(3842)%FILENAME = TRIM('retv_vars.10708_0016_002.cdf')
      TES(3843)%FILENAME = TRIM('retv_vars.10708_0017_003.cdf')
      TES(3844)%FILENAME = TRIM('retv_vars.10708_0017_004.cdf')
      TES(3845)%FILENAME = TRIM('retv_vars.10708_0018_004.cdf')
      TES(3846)%FILENAME = TRIM('retv_vars.10708_0020_002.cdf')
      TES(3847)%FILENAME = TRIM('retv_vars.10708_0020_003.cdf')
      TES(3848)%FILENAME = TRIM('retv_vars.10708_0020_004.cdf')
      TES(3849)%FILENAME = TRIM('retv_vars.10708_0021_002.cdf')
      TES(3850)%FILENAME = TRIM('retv_vars.10708_0021_003.cdf')
      TES(3851)%FILENAME = TRIM('retv_vars.10708_0021_004.cdf')
      TES(3852)%FILENAME = TRIM('retv_vars.10708_0022_002.cdf')
      TES(3853)%FILENAME = TRIM('retv_vars.10708_0022_004.cdf')
      TES(3854)%FILENAME = TRIM('retv_vars.10708_0023_002.cdf')
      TES(3855)%FILENAME = TRIM('retv_vars.10708_0027_002.cdf')
      TES(3856)%FILENAME = TRIM('retv_vars.10708_0027_003.cdf')
      TES(3857)%FILENAME = TRIM('retv_vars.10708_0027_004.cdf')
      TES(3858)%FILENAME = TRIM('retv_vars.10708_0028_002.cdf')
      TES(3859)%FILENAME = TRIM('retv_vars.10708_0028_003.cdf')
      TES(3860)%FILENAME = TRIM('retv_vars.10708_0028_004.cdf')
      TES(3861)%FILENAME = TRIM('retv_vars.10708_0054_003.cdf')
      TES(3862)%FILENAME = TRIM('retv_vars.10708_0054_004.cdf')
      TES(3863)%FILENAME = TRIM('retv_vars.10708_0055_003.cdf')
      TES(3864)%FILENAME = TRIM('retv_vars.10708_0055_004.cdf')
      TES(3865)%FILENAME = TRIM('retv_vars.10708_0056_003.cdf')
      TES(3866)%FILENAME = TRIM('retv_vars.10708_0057_002.cdf')
      TES(3867)%FILENAME = TRIM('retv_vars.10708_0057_003.cdf')
      TES(3868)%FILENAME = TRIM('retv_vars.10708_0057_004.cdf')
      TES(3869)%FILENAME = TRIM('retv_vars.10708_0058_002.cdf')
      TES(3870)%FILENAME = TRIM('retv_vars.10708_0058_003.cdf')
      TES(3871)%FILENAME = TRIM('retv_vars.10708_0060_002.cdf')
      TES(3872)%FILENAME = TRIM('retv_vars.10708_0060_003.cdf')
      TES(3873)%FILENAME = TRIM('retv_vars.10708_0064_002.cdf')
      TES(3874)%FILENAME = TRIM('retv_vars.10708_0064_003.cdf')
      TES(3875)%FILENAME = TRIM('retv_vars.10708_0064_004.cdf')
      TES(3876)%FILENAME = TRIM('retv_vars.10708_0067_004.cdf')
      TES(3877)%FILENAME = TRIM('retv_vars.10708_0068_002.cdf')
      TES(3878)%FILENAME = TRIM('retv_vars.10708_0068_003.cdf')
      TES(3879)%FILENAME = TRIM('retv_vars.10708_0069_003.cdf')
      TES(3880)%FILENAME = TRIM('retv_vars.10708_0069_004.cdf')
      TES(3881)%FILENAME = TRIM('retv_vars.10708_0071_002.cdf')
      TES(3882)%FILENAME = TRIM('retv_vars.10708_0108_003.cdf')
      TES(3883)%FILENAME = TRIM('retv_vars.10708_0110_004.cdf')
      TES(3884)%FILENAME = TRIM('retv_vars.10708_0111_002.cdf')
      TES(3885)%FILENAME = TRIM('retv_vars.10708_0172_002.cdf')
      TES(3886)%FILENAME = TRIM('retv_vars.10708_0186_002.cdf')
      TES(3887)%FILENAME = TRIM('retv_vars.10708_0186_004.cdf')
      TES(3888)%FILENAME = TRIM('retv_vars.10708_0187_002.cdf')
      TES(3889)%FILENAME = TRIM('retv_vars.10708_0187_003.cdf')
      TES(3890)%FILENAME = TRIM('retv_vars.10708_0187_004.cdf')
      TES(3891)%FILENAME = TRIM('retv_vars.10708_0188_002.cdf')
      TES(3892)%FILENAME = TRIM('retv_vars.10708_0190_002.cdf')
      TES(3893)%FILENAME = TRIM('retv_vars.10708_0190_003.cdf')
      TES(3894)%FILENAME = TRIM('retv_vars.10708_0198_004.cdf')
      TES(3895)%FILENAME = TRIM('retv_vars.10708_0199_002.cdf')
      TES(3896)%FILENAME = TRIM('retv_vars.10708_0199_003.cdf')
      TES(3897)%FILENAME = TRIM('retv_vars.10708_0199_004.cdf')
      TES(3898)%FILENAME = TRIM('retv_vars.10708_0200_002.cdf')
      TES(3899)%FILENAME = TRIM('retv_vars.10708_0200_003.cdf')
      TES(3900)%FILENAME = TRIM('retv_vars.10708_0200_004.cdf')
      TES(3901)%FILENAME = TRIM('retv_vars.10708_0201_002.cdf')
      TES(3902)%FILENAME = TRIM('retv_vars.10708_0201_003.cdf')
      TES(3903)%FILENAME = TRIM('retv_vars.10708_0201_004.cdf')
      TES(3904)%FILENAME = TRIM('retv_vars.10708_0202_002.cdf')
      TES(3905)%FILENAME = TRIM('retv_vars.10708_0213_002.cdf')
      TES(3906)%FILENAME = TRIM('retv_vars.10708_0213_003.cdf')
      TES(3907)%FILENAME = TRIM('retv_vars.10708_0213_004.cdf')
      TES(3908)%FILENAME = TRIM('retv_vars.10708_0214_003.cdf')
      TES(3909)%FILENAME = TRIM('retv_vars.10708_0219_002.cdf')
      TES(3910)%FILENAME = TRIM('retv_vars.10708_0220_002.cdf')
      TES(3911)%FILENAME = TRIM('retv_vars.10708_0220_003.cdf')
      TES(3912)%FILENAME = TRIM('retv_vars.10708_0220_004.cdf')
      TES(3913)%FILENAME = TRIM('retv_vars.10708_0221_003.cdf')
      TES(3914)%FILENAME = TRIM('retv_vars.10708_0222_003.cdf')
      TES(3915)%FILENAME = TRIM('retv_vars.10708_0223_004.cdf')
      TES(3916)%FILENAME = TRIM('retv_vars.10708_0231_003.cdf')
      TES(3917)%FILENAME = TRIM('retv_vars.10708_0231_004.cdf')
      TES(3918)%FILENAME = TRIM('retv_vars.10708_0235_003.cdf')
      TES(3919)%FILENAME = TRIM('retv_vars.10708_0235_004.cdf')
      TES(3920)%FILENAME = TRIM('retv_vars.10708_0244_002.cdf')
      TES(3921)%FILENAME = TRIM('retv_vars.10708_0244_003.cdf')
      TES(3922)%FILENAME = TRIM('retv_vars.10708_0244_004.cdf')
      TES(3923)%FILENAME = TRIM('retv_vars.10708_0245_002.cdf')
      TES(3924)%FILENAME = TRIM('retv_vars.10708_0245_003.cdf')
      TES(3925)%FILENAME = TRIM('retv_vars.10708_0245_004.cdf')
      TES(3926)%FILENAME = TRIM('retv_vars.10708_0246_003.cdf')
      TES(3927)%FILENAME = TRIM('retv_vars.10708_0247_003.cdf')
      TES(3928)%FILENAME = TRIM('retv_vars.10708_0248_003.cdf')
      TES(3929)%FILENAME = TRIM('retv_vars.10708_0248_004.cdf')
      TES(3930)%FILENAME = TRIM('retv_vars.10708_0249_002.cdf')
      TES(3931)%FILENAME = TRIM('retv_vars.10708_0249_003.cdf')
      TES(3932)%FILENAME = TRIM('retv_vars.10708_0249_004.cdf')
      TES(3933)%FILENAME = TRIM('retv_vars.10708_0250_003.cdf')
      TES(3934)%FILENAME = TRIM('retv_vars.10708_0250_004.cdf')
      TES(3935)%FILENAME = TRIM('retv_vars.10708_0251_002.cdf')
      TES(3936)%FILENAME = TRIM('retv_vars.10708_0252_004.cdf')
      TES(3937)%FILENAME = TRIM('retv_vars.10708_0256_002.cdf')
      TES(3938)%FILENAME = TRIM('retv_vars.10708_0256_003.cdf')
      TES(3939)%FILENAME = TRIM('retv_vars.10708_0259_002.cdf')
      TES(3940)%FILENAME = TRIM('retv_vars.10708_0259_003.cdf')
      TES(3941)%FILENAME = TRIM('retv_vars.10708_0259_004.cdf')
      TES(3942)%FILENAME = TRIM('retv_vars.10708_0260_002.cdf')
      TES(3943)%FILENAME = TRIM('retv_vars.10708_0260_003.cdf')
      TES(3944)%FILENAME = TRIM('retv_vars.10708_0260_004.cdf')
      TES(3945)%FILENAME = TRIM('retv_vars.10708_0261_002.cdf')
      TES(3946)%FILENAME = TRIM('retv_vars.10708_0262_002.cdf')
      TES(3947)%FILENAME = TRIM('retv_vars.10708_0267_002.cdf')
      TES(3948)%FILENAME = TRIM('retv_vars.10708_0267_003.cdf')
      TES(3949)%FILENAME = TRIM('retv_vars.10708_0268_002.cdf')
      TES(3950)%FILENAME = TRIM('retv_vars.10708_0268_003.cdf')
      TES(3951)%FILENAME = TRIM('retv_vars.10708_0268_004.cdf')
      TES(3952)%FILENAME = TRIM('retv_vars.10708_0269_003.cdf')
      TES(3953)%FILENAME = TRIM('retv_vars.10708_0272_002.cdf')
      TES(3954)%FILENAME = TRIM('retv_vars.10708_0278_003.cdf')
      TES(3955)%FILENAME = TRIM('retv_vars.10708_0302_004.cdf')
      TES(3956)%FILENAME = TRIM('retv_vars.10708_0303_003.cdf')
      TES(3957)%FILENAME = TRIM('retv_vars.10708_0304_002.cdf')
      TES(3958)%FILENAME = TRIM('retv_vars.10708_0304_003.cdf')
      TES(3959)%FILENAME = TRIM('retv_vars.10708_0304_004.cdf')
      TES(3960)%FILENAME = TRIM('retv_vars.10708_0305_002.cdf')
      TES(3961)%FILENAME = TRIM('retv_vars.10708_0306_003.cdf')
      TES(3962)%FILENAME = TRIM('retv_vars.10708_0307_004.cdf')
      TES(3963)%FILENAME = TRIM('retv_vars.10708_0308_004.cdf')
      TES(3964)%FILENAME = TRIM('retv_vars.10708_0310_002.cdf')
      TES(3965)%FILENAME = TRIM('retv_vars.10708_0310_003.cdf')
      TES(3966)%FILENAME = TRIM('retv_vars.10708_0310_004.cdf')
      TES(3967)%FILENAME = TRIM('retv_vars.10708_0315_003.cdf')
      TES(3968)%FILENAME = TRIM('retv_vars.10708_0315_004.cdf')
      TES(3969)%FILENAME = TRIM('retv_vars.10708_0316_002.cdf')
      TES(3970)%FILENAME = TRIM('retv_vars.10708_0316_003.cdf')
      TES(3971)%FILENAME = TRIM('retv_vars.10708_0316_004.cdf')
      TES(3972)%FILENAME = TRIM('retv_vars.10708_0317_002.cdf')
      TES(3973)%FILENAME = TRIM('retv_vars.10708_0317_003.cdf')
      TES(3974)%FILENAME = TRIM('retv_vars.10708_0318_002.cdf')
      TES(3975)%FILENAME = TRIM('retv_vars.10708_0318_003.cdf')
      TES(3976)%FILENAME = TRIM('retv_vars.10708_0318_004.cdf')
      TES(3977)%FILENAME = TRIM('retv_vars.10708_0320_004.cdf')
      TES(3978)%FILENAME = TRIM('retv_vars.10708_0322_002.cdf')
      TES(3979)%FILENAME = TRIM('retv_vars.10708_0323_003.cdf')
      TES(3980)%FILENAME = TRIM('retv_vars.10708_0323_004.cdf')
      TES(3981)%FILENAME = TRIM('retv_vars.10708_0324_003.cdf')
      TES(3982)%FILENAME = TRIM('retv_vars.10708_0359_003.cdf')
      TES(3983)%FILENAME = TRIM('retv_vars.10708_0364_003.cdf')
      TES(3984)%FILENAME = TRIM('retv_vars.10708_0365_003.cdf')
      TES(3985)%FILENAME = TRIM('retv_vars.10708_0366_003.cdf')
      TES(3986)%FILENAME = TRIM('retv_vars.10708_0366_004.cdf')
      TES(3987)%FILENAME = TRIM('retv_vars.10708_0368_002.cdf')
      TES(3988)%FILENAME = TRIM('retv_vars.10708_0369_003.cdf')
      TES(3989)%FILENAME = TRIM('retv_vars.10708_0369_004.cdf')
      TES(3990)%FILENAME = TRIM('retv_vars.10708_0373_003.cdf')
      TES(3991)%FILENAME = TRIM('retv_vars.10708_0373_004.cdf')
				

      ! for BLVMR: make date same as 2008
      !print*, ' subtract 1 year from TES date '
      !TES(:)%NYMD = TES(:)%NYMD - 10000

      ! Return to calling program 
      END SUBROUTINE INIT_TES_NH3
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TES_NH3
!
!*****************************************************************************
!  Subroutine CLEANUP_TES_NH3 deallocates all module arrays. (dkh, 02/15/09) 
!        
!  NOTES:
!
!******************************************************************************
!     

      IF ( ALLOCATED( NH3_SAVE ) )      DEALLOCATE( NH3_SAVE )


      ! Return to calling program 
      END SUBROUTINE CLEANUP_TES_NH3
!------------------------------------------------------------------------------

      END MODULE TES_NH3_MOD
