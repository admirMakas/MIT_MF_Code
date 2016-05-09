! $Id: diag48_mod.f,v 1.2 2009/11/18 07:09:33 daven Exp $
      MODULE DIAG48_MOD
!
!******************************************************************************
!  Module DIAG48_MOD contains variables and routines to save out 3-D 
!  timeseries output to disk (bmy, 7/20/04, 10/7/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) DO_SAVE_DIAG48 (LOGICAL ) : Switch to turn ND49 timeseries on/off 
!  (2 ) I0             (INTEGER ) : Lon offset between global & nested grid
!
!  Module Routines:
!  ============================================================================
!  (1 ) DIAG48                 : Main driver routine
!  (2 ) ITS_TIME_FOR_DIAG48    : Returns TRUE if it's time to save to disk
!  (3 ) INIT_DIAG48            : Gets variable values from "input_mod.f"
!
!  GEOS-CHEM modules referenced by "diag48_mod.f" 
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f      : Module w/ arrays for DAO met fields
!  (3 ) file_mod.f     : Module w/ file unit numbers & error checks
!  (4 ) grid_mod.f     : Module w/ horizontal grid information   
!  (5 ) pbl_mix_mod.f  : Module w/ routines for PBL height & mixing
!  (6 ) pressure_mod.f : Module w/ routines to compute P(I,J,L)
!  (7 ) time_mod.f     : Module w/ routines for computing time & date 
!  (8 ) tracer_mod.f   : Module w/ GEOS-CHEM tracer array STT etc.  
!  (9 ) tracerid_mod.f : Module w/ pointers to tracers & emissions
!
!  ND48 tracer numbers:
!  ============================================================================
!  1 - N_TRACERS : GEOS-CHEM transported tracers            [v/v      ]
!  74            : OH concentration                         [molec/cm3]
!  75            : NO2 concentration                        [v/v      ]
!  76            : PBL heights                              [m        ]
!  77            : PBL heights                              [levels   ]
!  78            : Air density                              [molec/cm3]
!  79            : 3-D Cloud fractions                      [unitless ]
!  80            : Column optical depths                    [unitless ]
!  81            : Cloud top heights                        [hPa      ]
!  82            : Sulfate aerosol optical depth            [unitless ]
!  83            : Black carbon aerosol optical depth       [unitless ]
!  84            : Organic carbon aerosol optical depth     [unitless ]
!  85            : Accumulation mode seasalt optical depth  [unitless ]
!  86            : Coarse mode seasalt optical depth        [unitless ]
!  87            : Total dust optical depth                 [unitless ]
!  88            : Total seasalt tracer concentration       [unitless ]
!  89            : Pure O3 (not Ox) concentration           [v/v      ]
!  90            : NO concentration                         [v/v      ]
!  91            : NOy concentration                        [v/v      ]
!  92            : RESERVED FOR FUTURE USE
!  93            : Grid box heights                         [m        ]
!  94            : Relative humidity                        [%        ]
!  95            : Sea level pressure                       [hPa      ]
!  96            : Zonal wind (a.k.a. U-wind)               [m/s      ]
!  97            : Meridional wind (a.k.a. V-wind)          [m/s      ]
!  98            : P(surface) - PTOP                        [hPa      ]
!  99            : Temperature                              [K        ]
!  
!  NOTES:
!  (1 ) Now save out cld frac and grid box heights (bmy, 4/20/05)
!  (2 ) Remove TRCOFFSET because it's always zero.  Now call GET_HALFPOLAR
!        to get the value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Minor bug fixes in DIAG48 (cdh, bmy, 2/11/08)
!  (6 ) Bug fix: replace "PS-PTOP" with "PEDGE-$" (phs, bmy, 10/7/08)
!  (7 ) Modified to archive O3, NO, NOy as tracers 89, 90, 91  (tmf, 10/22/07)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag48_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE
 
      ! ... except these variables ...
      PUBLIC :: DO_SAVE_DIAG48
      PUBLIC :: ND48_MAX_STATIONS

      ! ... and these routines
      PUBLIC :: CLEANUP_DIAG48
      PUBLIC :: DIAG48
      PUBLIC :: INIT_DIAG48 
      PUBLIC :: ITS_TIME_FOR_DIAG48

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      LOGICAL              :: DO_SAVE_DIAG48
      INTEGER, PARAMETER   :: ND48_MAX_STATIONS=1000
      INTEGER              :: ND48_FREQ
      INTEGER              :: ND48_N_STATIONS
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180=1 
      TYPE (XPLEX)               :: LONRES
      TYPE (XPLEX)               :: LATRES
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = 'ND48 station timeseries'
      CHARACTER(LEN=80)    :: TITLE
      CHARACTER(LEN=255)   :: ND48_OUTPUT_FILE

      ! Arrays
      INTEGER, ALLOCATABLE :: ND48_I(:)
      INTEGER, ALLOCATABLE :: ND48_J(:)
      INTEGER, ALLOCATABLE :: ND48_L(:)
      INTEGER, ALLOCATABLE :: ND48_N(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DIAG48 
!
!******************************************************************************
!  Subroutine DIAG48 saves station time series diagnostics to disk.
!  (bmy, bey, amf, 6/1/99, 10/7/08)
!
!  NOTES:
!  (1 ) Remove reference to "CMN".  Also now get PBL heights in meters and
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (2 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D cloud 
!        fraction as tracer #79 and box height as tracer #93.  Now remove
!        reference to PBL from "dao_mod.f" (bmy, 4/20/05)
!  (3 ) Remove references to TRCOFFSET because it's always zero.  Now call 
!        GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag value for 
!        GEOS or GCAP grids.  (bmy, 6/28/05)
!  (4 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (7 ) Bug fix: unit for tracer #77 should be "layers".  Also RH should be 
!        tracer #17 under "TIME-SER" category. (cdh, bmy, 2/11/08)
!  (8 ) Bug fix: replace "PS-PTOP" with "PEDGE-$" (phs, bmy, 10/7/08)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : BPCH2
      USE DAO_MOD,      ONLY : AD,      AIRDEN, BXHEIGHT, CLDF 
      USE DAO_MOD,      ONLY : CLDTOPS, OPTD,   RH,       SLP   
      USE DAO_MOD,      ONLY : T,       UWND,   VWND 
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_ND48
      USE GRID_MOD,     ONLY : GET_XOFFSET,        GET_YOFFSET
      USE PBL_MIX_MOD,  ONLY : GET_PBL_TOP_L,      GET_PBL_TOP_m
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_LOCALTIME,      GET_NYMD
      USE TIME_MOD,     ONLY : GET_NHMS,           GET_TAU
      USE TIME_MOD,     ONLY : EXPAND_DATE,        ITS_A_NEW_DAY
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM, N_TRACERS
      USE TRACER_MOD,   ONLY : STT,                TCVV  
      USE TRACER_MOD,   ONLY : XNUMOLAIR
      USE TRACERID_MOD, ONLY : IDTHNO3, IDTN2O5, IDTHNO4, IDTNOX
      USE TRACERID_MOD, ONLY : IDTPAN,  IDTPMN,  IDTPPN,  IDTOX   
      USE TRACERID_MOD, ONLY : IDTR4N2, IDTSALA, IDTSALC 

#     include "cmn_fj.h"    ! Size parameters + FAST-J stuff
#     include "jv_cmn.h"    ! ODAER, ODMDUST
#     include "CMN_O3"      ! XNUMOLAIR
#     include "CMN_GCTM"    ! SCALE_HEIGHT

      ! Local variables
      LOGICAL, SAVE      :: IS_CLDTOPS,  IS_OPTD, IS_SEASALT, IS_SLP
      LOGICAL, SAVE      :: IS_FULLCHEM, IS_Ox,   IS_NOx,     IS_NOy
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: GMTRC, I, I0, J, J0, L, N, K, R, TMP, W
      TYPE (XPLEX)             :: LT, TAU, Q(LLPAR), SCALE400nm
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=255) :: FILENAME

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER            :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      !=================================================================
      ! DIAG48 begins here!
      !=================================================================

      ! Get grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Set logical flags on first timestep
         IS_CLDTOPS  = ALLOCATED( CLDTOPS )
         IS_OPTD     = ALLOCATED( OPTD    )
         IS_SLP      = ALLOCATED( SLP     )
         IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
         IS_SEASALT  = ( IDTSALA > 0 .and. IDTSALC > 0 )
         IS_Ox       = ( IS_FULLCHEM .and. IDTOX   > 0 )
         IS_NOx      = ( IS_FULLCHEM .and. IDTNOX  > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. 
     &                   IDTNOX  > 0 .and. IDTPAN  > 0 .and.
     &                   IDTHNO3 > 0 .and. IDTPMN  > 0 .and.
     &                   IDTPPN  > 0 .and. IDTR4N2 > 0 .and.
     &                   IDTN2O5 > 0 .and. IDTHNO4 > 0 ) 

         ! Replace YYYYMMDD tokens in filename
         FILENAME = TRIM( ND48_OUTPUT_FILE )
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

         ! Open file for writing
         CALL OPEN_ND48_FILE( FILENAME )

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Store station timeseries data in the Q array
      !=================================================================

      ! Get TAU value
      TAU = GET_TAU()

      ! Loop over each station
      DO W = 1, ND48_N_STATIONS

         ! Get lon, lat, alt, tracer indices
         I = ND48_I(W)             
         J = ND48_J(W)  
         K = ND48_L(W) 
         N = ND48_N(W)

         ! Get local time
         LT = GET_LOCALTIME( I )

         ! Initialize
         Q(:) = 0

         IF ( N <= N_TRACERS ) THEN

            !------------------------------------
            ! GEOS-CHEM tracers [v/v]
            !------------------------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick the unit
            GMTRC    = N

            DO L = 1, K
               Q(L) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
            ENDDO

         ELSE IF ( N == 89 .and. IS_Ox ) THEN
 
            !------------------------------------
            ! PURE O3 CONCENTRATION [v/v]
            !------------------------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick the unit
            GMTRC    = N_TRACERS + 1

            DO L = 1, K
               Q(L) = STT(I,J,L,IDTOX) * TCVV(IDTOX)  / 
     &                AD(I,J,L)        * FRACO3(I,J,L)
            ENDDO
               
         ELSE IF ( N == 90 .and. IS_NOx ) THEN

            !------------------------------------
            ! NO CONCENTRATION [v/v]
            !------------------------------------ 
            CATEGORY = 'TIME-SER'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMTRC    = 9

            DO L = 1, K
               Q(L) = STT(I,J,L,IDTNOX) * TCVV(IDTNOX) * 
     &                FRACNO(I,J,L)     / AD(I,J,L)
            ENDDO

         ELSE IF ( N == 91 .and. IS_NOy ) THEN

            !-------------------------------------
            ! NOy CONCENTRATION [v/v]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMTRC    = 3

            DO L = 1, K

               ! Initialize
               Q(L) = 0d0

               ! NOx
               Q(L) = Q(L) + ( TCVV(IDTNOX)        * 
     &                         STT(I,J,L,IDTNOX)   / AD(I,J,L) )

               ! PAN
               Q(L) = Q(L) + ( TCVV(IDTPAN)        * 
     &                         STT(I,J,L,IDTPAN)   / AD(I,J,L) )

               ! HNO3
               Q(L) = Q(L) + ( TCVV(IDTHNO3)       *
     &                         STT(I,J,L,IDTHNO3)  / AD(I,J,L) )
            
               ! PMN
               Q(L) = Q(L) + ( TCVV(IDTPMN)        *
     &                         STT(I,J,L,IDTPMN)   / AD(I,J,L) )

               ! PPN
               Q(L) = Q(L) + ( TCVV(IDTPPN)        * 
     &                         STT(I,J,L,IDTPPN)   / AD(I,J,L) )
 
               ! R4N2
               Q(L) = Q(L) + ( TCVV(IDTR4N2)       *
     &                         STT(I,J,L,IDTR4N2)  / AD(I,J,L) )
            
               ! N2O5
               Q(L) = Q(L) + ( 2d0 * TCVV(IDTN2O5) *
     &                         STT(I,J,L,IDTN2O5)  / AD(I,J,L) )
                        
               ! HNO4
               Q(L) = Q(L) + ( TCVV(IDTHNO4)       *
     &                         STT(I,J,L,IDTHNO4)  / AD(I,J,L) )
            ENDDO

         ELSE IF ( N == 74 .and. IS_FULLCHEM ) THEN

            !------------------------------------
            ! OH CONCENTRATION [molec/cm3]
            !------------------------------------ 
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMTRC    = 2

            DO L = 1, K
               Q(L) = SAVEOH(I,J,L)
            ENDDO


         ELSE IF ( N == 75 .and. IS_FULLCHEM ) THEN

            !------------------------------------
            ! NO2 CONCENTRATION [molec/cm3]
            !------------------------------------ 
            CATEGORY = 'TIME-SER'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMTRC    = 19

            DO L = 1, K
               Q(L) = SAVENO2(I,J,L)
            ENDDO

         ELSE IF ( N == 76 ) THEN

            !-------------------------------------
            ! PBL HEIGHTS [m] 
            !-------------------------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'  
            GMTRC    = 1

            IF ( K == 1 ) THEN
               Q(1) = GET_PBL_TOP_m( I, J )
            ENDIF
               
         ELSE IF ( N == 77 ) THEN

            !-------------------------------------
            ! PBL HEIGHTS [layers] 
            !-------------------------------------    
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'
            GMTRC    = 2

            IF ( K == 1 ) THEN
               Q(1) = GET_PBL_TOP_L( I, J )
            ENDIF

         ELSE IF ( N == 78 ) THEN 

            !-------------------------------------
            ! AIR DENSITY [molec/cm3]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMTRC    = 22

            DO L = 1, K
               Q(L) = AIRDEN(L,I,J) * XNUMOLAIR * 1d-6
            ENDDO

         ELSE IF ( N == 79 ) THEN

            !-------------------------------------
            ! CLOUD FRACTIONS [unitless]
            !------------------------------------- 
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMTRC    = 19

            DO L = 1, K
               Q(L) = CLDF(L,I,J)
            ENDDO

         ELSE IF ( N == 80 .and. IS_OPTD ) THEN
            
            !---------------------------------------
            ! COLUMN OPTICAL DEPTHS [unitless]
            !---------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMTRC    = 20

            Q(1) = SUM( OPTD(:,I,J) )

         ELSE IF ( N == 81 .and. IS_CLDTOPS ) THEN

            !---------------------------------------
            ! CLOUD TOP HEIGHTS [hPa]
            !---------------------------------------
            CATEGORY = 'TIME_SER'
            UNIT     = 'hPa'
            GMTRC    = 21

            IF ( K == 1 ) THEN
               Q(1) = GET_PEDGE( I, J, CLDTOPS(I,J) )
            ENDIF

         ELSE IF ( N == 82 ) THEN

            !---------------------------------------
            ! SULFATE AOD @ 400nm [unitless]
            !--------------------------------------- 
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 6

            DO R = 1, NRH

               ! Scaling factor for 400nm
               SCALE400nm = QAA(2,IND(1)+R-1) / QAA(4,IND(1)+R-1) 

               DO L = 1, K
                  Q(L) = ODAER(I,J,L,R) * SCALE400nm
               ENDDO
            ENDDO

         ELSE IF ( N == 83 ) THEN

            !-------------------------------------
            ! BLACK CARBON AOD @ 400nm [unitless]
            !-------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 9            

            DO R = 1, NRH

               ! Scaling factor for 400nm
               SCALE400nm = QAA(2,IND(2)+R-1) / QAA(4,IND(2)+R-1) 

               DO L = 1, K
                  Q(L) = ODAER(I,J,L,NRH+R) * SCALE400nm
               ENDDO
            ENDDO

         ELSE IF ( N == 84 ) THEN

            !-----------------------------------
            ! ORGANIC CARBON AOD [unitless]
            !-----------------------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 12

            DO R = 1, NRH

               ! Scaling factor for 400nm
               SCALE400nm = QAA(2,IND(3)+R-1) / QAA(4,IND(3)+R-1) 

               DO L = 1, K
                  Q(L) = ODAER(I,J,L,2*NRH+R) * SCALE400nm
               ENDDO
            ENDDO

         ELSE IF ( N == 85 ) THEN
            
            !-------------------------------------
            ! ACCUM SEASALT AOD @ 400nm [unitless]
            !-------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 15

            DO R = 1, NRH

               ! Scaling factor for 400nm
               SCALE400nm = QAA(2,IND(4)+R-1) / QAA(4,IND(4)+R-1) 

               DO L = 1, K
                  Q(L) = ODAER(I,J,L,3*NRH+R) * SCALE400nm
               ENDDO
            ENDDO

         ELSE IF ( N == 86 ) THEN

            !-------------------------------------
            ! COARSE SEASALT AOD @ 40nm [unitless]
            !-------------------------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 18

            DO R = 1, NRH

               ! Scaling factor for 400nm
               SCALE400nm = QAA(2,IND(5)+R-1) / QAA(4,IND(5)+R-1) 

               DO L = 1, K
                  Q(L) = ODAER(I,J,L,4*NRH+R) * SCALE400nm
               ENDDO
            ENDDO

         ELSE IF ( N == 87 ) THEN

            !-----------------------------------
            ! TOTAL DUST OPT DEPTH [unitless]
            !-----------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMTRC     = 4

            DO R = 1, NDUST

               ! Scaling factor for 400nm
               SCALE400nm = QAA(2,IND(6)+R-1) / QAA(4,IND(6)+R-1) 

               DO L = 1, K
                  !Q(L) = ODMDUST(I,J,L,R)
                  Q(L) = ODMDUST(I,J,L,R) + SCALE400nm
               ENDDO
            ENDDO

         ELSE IF ( N == 88 .and. IS_SEASALT ) THEN

            !-----------------------------------
            ! TOTAL SEASALT TRACER [v/v]
            !-----------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMTRC    = 24

            DO L = 1, K
               Q(L) = ( STT(I,J,L,IDTSALA) + STT(I,J,L,IDTSALC) ) *
     &                  TCVV(IDTSALA)      / AD(I,J,L) 
            ENDDO

         ELSE IF ( N == 93 ) THEN 

            !-----------------------------------
            ! GRID BOX HEIGHTS [m]
            !----------------------------------- 
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMTRC    = 1

            DO L = 1, K
               Q(L) = BXHEIGHT(I,J,L)
            ENDDO

         ELSE IF ( N == 94 ) THEN 

            !-----------------------------------
            ! RELATIVE HUMIDITY [%]
            !----------------------------------- 
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMTRC    = 17

            DO L = 1, K
               Q(L) = RH(I,J,L)
            ENDDO

         ELSE IF ( N == 95 .and. IS_SLP ) THEN

            !-----------------------------------
            ! SEA LEVEL PRESSURE [hPa]
            !----------------------------------- 
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMTRC    = 21            

            IF ( K == 1 ) THEN
               Q(1) = SLP(I,J)
            ENDIF

         ELSE IF ( N == 96 ) THEN

            !-----------------------------------
            ! ZONAL (U) WIND [M/S]
            !----------------------------------- 
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMTRC    = 1

            DO L = 1, K
               Q(L) = UWND(I,J,L)
            ENDDO
            
         ELSE IF ( N == 97 ) THEN 

            !-----------------------------------
            ! ZONAL (V) WIND [M/S]
            !----------------------------------- 
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMTRC    = 2

            DO L = 1, K
               Q(L) = VWND(I,J,L)
            ENDDO

         ELSE IF ( N == 98 ) THEN

            !-----------------------------------
            ! PSURFACE - PTOP [hPa]
            !----------------------------------- 
            CATEGORY = 'PEDGE-$'
            UNIT     = 'hPa'
            GMTRC    = 1            
            
            IF ( K == 1 ) THEN
               Q(1) = GET_PEDGE(I,J,1) - PTOP
            ENDIF

         ELSE IF ( N == 99 ) THEN

            !-----------------------------------
            ! TEMPERATURE [K]
            !----------------------------------- 
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'K'
            GMTRC    = 3            

            DO L = 1, K
               Q(L) = T(I,J,L)
            ENDDO

         ELSE

            ! Skip other tracers
            CYCLE

         ENDIF
         
         !==============================================================
         ! Write each station to a bpch file
         !==============================================================
         CALL BPCH2( IU_ND48,   MODELNAME, LONRES,   LATRES,       
     &               HALFPOLAR, CENTER180, CATEGORY, GMTRC,        
     &               UNIT,      TAU,       TAU,      RESERVED,  
     &               1,         1,         K,        I+I0,         
     &               J+J0,      1,      ( Q(1:K) ) )

      ENDDO

      !=================================================================
      ! Close the file at the proper time
      !=================================================================
!      IF ( ITS_TIME_TO_CLOSE_FILE() ) THEN
!
!         ! Expand date tokens in the file name
!         FILENAME = TRIM( ND49_OUTPUT_FILE )
!         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
!
!         ! Echo info
!         WRITE( 6, 120 ) TRIM( FILENAME )
! 120     FORMAT( '     - DIAG49: Closing file : ', a )
!
!         ! Close file
!         CLOSE( IU_ND49 ) 
!      ENDIF

      ! Flush the file once per day
      IF ( ITS_A_NEW_DAY() ) CALL FLUSH( IU_ND48 ) 

      ! Return to calling program
      END SUBROUTINE DIAG48

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_ND48_FILE( FILENAME )
!
!******************************************************************************
!  Subroutine OPEN_ND48_FILE opens a binary punch file (version 2.0)
!  for writing GEOS-CHEM ND48 station timeseries data. (bmy, 7/30/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of the file to be opened
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_ND48, IOERROR

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

      ! Local variables
      INTEGER                      :: IOS
      CHARACTER(LEN=40)            :: FTI
      CHARACTER(LEN=80)            :: TITLE

      !=================================================================
      ! OPEN_ND48_FOR_WRITE begins here!
      !=================================================================

      ! Initialize
      FTI   = 'CTM bin 02 -- GEOS-CHEM station ts' 
      TITLE = 'GEOS-CHEM ND48 station timeseries diagnostic'

      ! Open file for output
      OPEN( IU_ND48,    FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_ND48, 'open_nd48_file:1' )

      ! Write file type identifier
      WRITE ( IU_ND48, IOSTAT=IOS ) FTI
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_ND48, 'open_nd48_file:2' )

      ! Write top title
      WRITE ( IU_ND48, IOSTAT=IOS ) TITLE
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_ND48, 'open_nd48_file:3' )

      ! Return to calling program
      END SUBROUTINE OPEN_ND48_FILE

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_DIAG48() RESULT( ITS_TIME )
!
!******************************************************************************
!  Function ITS_TIME_FOR_DIAG48 returns TRUE if it's time for the next
!  timeseries data write. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_ELAPSED_MIN

      ! Local variables
      LOGICAL :: ITS_TIME
      INTEGER :: XMIN

      !=================================================================
      ! ITS_TIME_FOR_DIAG48 begins here!
      !=================================================================
      
      ! Get elapsed minutes
      XMIN     = GET_ELAPSED_MIN()

      ! Is it time to save the next timeseries station?
      ITS_TIME = ( DO_SAVE_DIAG48 .and. MOD( XMIN, ND48_FREQ ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_DIAG48

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG48( DO_ND48, FREQ, N_STA, IARR, 
     &                        JARR,    LARR, NARR,  FILE )
!
!******************************************************************************
!  Subroutine INIT_DIAG48 allocates and zeroes all module arrays (bmy, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DO_ND48 (LOGICAL ) : Switch to turn on ND49 timeseries diagnostic
!  (2 ) FREQ    (INTEGER ) : Frequency for saving to disk [min]
!  (3 ) N_STA   (INTEGER ) : Number of ND48 stations from "input_mod.f"
!  (4 ) IARR    (INTEGER ) : Array w/ ND48 lon    indices from "input_mod.f"
!  (5 ) JARR    (INTEGER ) : Array w/ ND48 lat    indices from "input_mod.f"
!  (6 ) LARR    (INTEGER ) : Array w/ ND48 alt    indices from "input_mod.f"
!  (7 ) NARR    (INTEGER ) : Array w/ ND48 tracer indices from "input_mod.f"
!  (8 ) FILE    (CHAR*255) : ND48 output file name read by "input_mod.f"
! 
!  NOTES:
!  (1 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_MODELNAME, GET_HALFPOLAR
      USE ERROR_MOD, ONLY : ALLOC_ERR, ERROR_STOP

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      LOGICAL,            INTENT(IN) :: DO_ND48
      INTEGER,            INTENT(IN) :: FREQ
      INTEGER,            INTENT(IN) :: N_STA
      INTEGER,            INTENT(IN) :: IARR(ND48_MAX_STATIONS)
      INTEGER,            INTENT(IN) :: JARR(ND48_MAX_STATIONS)
      INTEGER,            INTENT(IN) :: LARR(ND48_MAX_STATIONS)
      INTEGER,            INTENT(IN) :: NARR(ND48_MAX_STATIONS)
      CHARACTER(LEN=255), INTENT(IN) :: FILE

      ! Local variables
      INTEGER                        :: AS, N
      CHARACTER(LEN=255)             :: LOCATION

      !=================================================================
      ! INIT_DIAG48 begins here!
      !=================================================================

      ! Location string
      LOCATION         = 'INIT_DIAG48 ("diag48_mod.f")'
     
      !-----------------------------
      ! Get values from input_mod.f
      !-----------------------------
      DO_SAVE_DIAG48   = DO_ND48
      ND48_FREQ        = FREQ
      ND48_N_STATIONS  = N_STA
      ND48_OUTPUT_FILE = TRIM( FILE )

      ! Error check
      IF ( ND48_N_STATIONS > ND48_MAX_STATIONS ) THEN
         CALL ERROR_STOP( 'Too many ND48 stations!', LOCATION )
      ENDIF

      !------------------------------
      ! Allocate module arrays
      !------------------------------
      ALLOCATE( ND48_I( ND48_N_STATIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_I' )

      ALLOCATE( ND48_J( ND48_N_STATIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_J' )

      ALLOCATE( ND48_L( ND48_N_STATIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_L' )

      ALLOCATE( ND48_N( ND48_N_STATIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_N' )

      !-------------------------------
      ! Copy values & error check
      !-------------------------------
      DO N = 1, ND48_N_STATIONS

         ! Error check longitude
         IF ( IARR(N) < 1 .or. IARR(N) > IIPAR ) THEN
            CALL ERROR_STOP( 'Bad longitude index!', LOCATION )
         ELSE
            ND48_I(N) = IARR(N)
         ENDIF

         ! Error check latitude
         IF ( JARR(N) < 1 .or. JARR(N) > JJPAR ) THEN
            CALL ERROR_STOP( 'Bad latitude index!', LOCATION )
         ELSE
            ND48_J(N) = JARR(N)
         ENDIF

         ! Error check longitude
         IF ( LARR(N) < 1 .or. LARR(N) > LLPAR ) THEN
            CALL ERROR_STOP( 'Bad altitude index!', LOCATION )
         ELSE
            ND48_L(N) = LARR(N)
         ENDIF

         ! Tracer array
         ND48_N(N) = NARR(N)
      ENDDO

      !-------------------------------
      ! For bpch file output
      !-------------------------------
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()

      ! Return to calling program
      END SUBROUTINE INIT_DIAG48

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG48
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG48 deallocates all module arrays (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG48 begins here!
      !=================================================================
      IF ( ALLOCATED( ND48_I ) ) DEALLOCATE( ND48_I )
      IF ( ALLOCATED( ND48_J ) ) DEALLOCATE( ND48_J )
      IF ( ALLOCATED( ND48_L ) ) DEALLOCATE( ND48_L )
      IF ( ALLOCATED( ND48_N ) ) DEALLOCATE( ND48_N )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG48

!------------------------------------------------------------------------------

      ! Return to calling program
      END MODULE DIAG48_MOD
