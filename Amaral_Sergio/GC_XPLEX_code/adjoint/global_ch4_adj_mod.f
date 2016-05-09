!$Id: global_ch4_adj_mod.f,v 1.3 2012/03/02 06:10:57 daven Exp $
      MODULE GLOBAL_CH4_ADJ_MOD
!
!******************************************************************************
!  Module GLOBAL_CH4_ADJ_MOD contains variables and routines used for the 
!  adjoint CH4 simulation (adj_group, kjw, 2/22/10, adj32_023)
!
!  To perform identical twin tests using TES pseudo-observations, I made some
!  rather ugly work-arounds in the standard adjoint code distribution. These
!  changes are not incorporated into the standard code because they would
!  make the code messy and slower. If you want a copy of the v8 adjoint code 
!  used to perform identical twin tests, contact Kevin Wecht 
!  (wecht-at-fas.harvard.edu) (kjw, 7/06/11)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BOH          : Array to hold monthly mean OH concentrations
!  (2 ) CH4_EMIS_ADJ : Array to hold methane emissions
!  (3 ) COPROD       : Array to hold CH4 loss from stratosphere
!  (4 ) TAVG_ADJ     : Array to hold average daily temperature
!  (5 ) BAIRDENS     : Array to hold density of air
!  (6 ) FMOL_CH4     : Molecular weight of CH4 [kg / mol]
!  (7 ) XNUMOL_CH4   : molec CH4 / kg CH4
!  (8 ) IU_A6_CH4_ADJ: file unit number
!
!  Module Routines:
!  ============================================================================
!  (1 ) EMISSCH4_ADJ         : Adjoint of CH4 emissions
!  (2 ) CHEMCH4_ADJ          : Adjoint of CH4 chemistry
!  (3 ) CH4_DECAY_ADJ        : Adjoint of decay rate of CH4 by OH.
!  (4 ) CH4_STRAT_ADJ        : Adjoint of loss of CH4 in the stratosphere
!  (5 ) READ_COPROD          : Reads prescribed zonal CH4 loss from stratosphere
!  (6 ) CH4_AVGTP_AVG        : Gets 24h avg temp and pressure for CHEMCH4_ADJ
!  (7 ) OPEN_A6_CH4_ADJ      : Opens A6 met files for use by CH4_AVGTP_ADJ
!  (8 ) READ_A6_CH4_ADJ      : Reads A6 met files for use by CH4_AVGTP_ADJ
!  (9 ) FIND_CLOSEST_A6      : Finds date and time of nearest A6 met field
!  (10) GET_SCALE_GROUP      : Determines which temporal/spatial scaling index to use.
!  (11) INIT_CH4_ADJ         : Allocates and initializes module arrays
!  (12) CLEANUP_CH4_ADJ      : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by global_ch4_mod.f
!  ============================================================================
!  (1 ) biofuel_mod.f    : Module w/ routines to read biofuel emissions
!  (2 ) biomass_mod.f    : Module w/ routines to read biomass emissions
!  (3 ) bpch2_mod.f      : Module w/ routines for binary punch file I/O
!  (4 ) dao_mod.f        : Module w/ arrays for DAO met fields
!  (5 ) diag_mod.f       : Module w/ GEOS-CHEM diagnostic arrays
!  (6 ) diag_pl_mod.f    : Module w/ routines for prod & loss diag's
!  (7 ) directory_mod.f  : Module w/ GEOS-CHEM data & met field dirs
!  (8 ) error_mod.f      : Module w/ I/O error and NaN check routines
!  (9 ) geia_mod         : Module w/ routines to read anthro emissions
!  (10) global_oh_mod.f  : Module w/ routines to read 3-D OH field
!  (11) global_nox_mod.f : Module w/ routines to read 3-D NOx field
!  (12) grid_mod.f       : Module w/ horizontal grid information
!  (13) logical_mod.f    : Module w/ GEOS-CHEM logical switches
!  (14) pressure_mod.f   : Module w/ routines to compute P(I,J,L)
!  (15) time_mod.f       : Module w/ routines for computing time & date
!  (16) tracer_mod.f     : Module w/ GEOS-CHEM tracer array STT etc.
!  (17) tropopause_mod.f : Module w/ routines to read ann mean tropopause
!  (18) logical_adj_mod.f: Module w/ adj logical flags
!
!  NOTES:
!******************************************************************************
!
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep all variables and routines 
      ! from being seen outside "global_ch4_adj_mod.f"
      ! Exceptions for those listed under "PUBLIC ROUTINES"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE

      ! PUBLIC ROUTINES
      PUBLIC :: CHEMCH4_ADJ
      PUBLIC :: EMISSCH4_ADJ
      PUBLIC :: CLEANUP_GLOBAL_CH4_ADJ
      
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      TYPE (XPLEX),  ALLOCATABLE   :: CH4_EMIS_ADJ(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE   :: BOH(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE   :: CH4LOSS(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE   :: COPROD(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE   :: TAVG_ADJ(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE   :: BAIRDENS(:,:,:)

      ! FMOL_CH4     - kg CH4 / mole CH4
      ! XNUMOL_CH4   - molecules CH4 / kg CH4
      TYPE (XPLEX),  PARAMETER     :: FMOL_CH4     = xplex(16d-3,0d0)
      TYPE (XPLEX),PARAMETER::XNUMOL_CH4=xplex(6.022d+23/FMOL_CH4%r,0d0)
      INTEGER                :: IU_A6_CH4_ADJ

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS


!------------------------------------------------------------------------------

      SUBROUTINE EMISSCH4_ADJ
!
!******************************************************************************
!  Subroutine EMISSCH4_ADJ does adjoint of CH4 emissions
!  (adj_group, kjw, 2/22/10)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,   ONLY : ADCH4EMS, EMS_SF_ADJ, STT_ADJ
      USE GRID_MOD,         ONLY : GET_XOFFSET, GET_YOFFSET
      USE GRID_MOD,         ONLY : GET_AREA_CM2
      USE TIME_MOD,         ONLY : GET_TS_EMIS
      USE TIME_MOD,         ONLY : GET_MONTH, GET_YEAR
      USE TIME_MOD,         ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE GLOBAL_CH4_MOD,   ONLY : WETLAND_EMIS, BIOBURN_EMIS
      USE GLOBAL_CH4_MOD,   ONLY : RICE_EMIS!,    BIOFUEL_EMIS
      USE GLOBAL_CH4_MOD,   ONLY : ASEASONAL_ANTHRO_EMIS
      USE GLOBAL_CH4_MOD,   ONLY : ASEASONAL_NATURAL_EMIS
      USE TIME_MOD,         ONLY : GET_NYMD,     GET_NHMS
      
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      !LOGICAL, SAVE          :: FIRSTEMIS = .TRUE.
      LOGICAL, SAVE          :: LASTOFMONTH = .TRUE.
      LOGICAL, SAVE          :: LASTOFYEAR  = .TRUE.
      INTEGER                :: I, J, I0, J0, M, IREF, JREF

      ! Local variables
      TYPE (XPLEX)                 :: E_CH4, DTSRCE, AREA_CM2
      TYPE (XPLEX)                 :: CH4_EMIS_for_SF(IIPAR,JJPAR)


      !=================================================================
      ! EMISSCH4_ADJ begins here!
      !=================================================================

      WRITE(6,*) '% --- ENTERING EMISSCH4_ADJ! ---'

      ! Initialize GLOBAL_CH4_ADJ_MOD variables
      ! Do here because CHEMCH4 isn't called at first midnight
      ! NO. Initialize global_ch4_adj_mod variables in chemch4_adj 
      ! because chemistry and emissions now have the same time step 
      ! and chemistry is called first (kjw, 12/2/2011).
      !IF ( FIRSTEMIS ) THEN
         !CALL INIT_GLOBAL_CH4_ADJ
         !FIRSTEMIS=.FALSE.
      !ENDIF

      ! Determine group (temporal)
      M = GET_SCALE_GROUP()
      ! Print out scaling info
      WRITE(6,*) '    - READ / RESCALE CHEMISTRY: 
     &     use SCALE_GROUP ',  M 

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()


      !===================================================================
      ! Emissions are read or calculated at the first of every:
      !   1) Emission time step - Natural Wetlands (from J Kaplan)
      !   2) Month              - Biomass Burning and Rice
      !   3) Year               - All other sources
      !
      ! Emissions are stored in CH4_EMIS_ADJ(IIPAR,JJPAR,N).
      !     Where N = 1:12
      !       1. Total Emissions (including soil absorption, counted neg.)
      !       2. Oil and Gas Processing
      !       3. Coal Mining
      !       4. Livestock
      !       5. Waste
      !       6. Biofuel
      !       7. Rice
      !       8. Other Anthropogenic
      !       9. Biomass Burning
      !       10. Wetlands
      !       11. Soil Absorption
      !       12. Other Natural
      !
      ! Emissions are then summed 
      !                                                      (kjw, 6/4/09)
      !===================================================================

      ! Do Adjoint CH4 emissions!

      !4.1 Wetland Emissions (CH4_WTL, #10)
      CALL WETLAND_EMIS( CH4_EMIS_ADJ )


      IF ( LASTOFMONTH ) THEN

         !4.2 Biomass Burning emissions (CH4_BBN, #9)
         CALL BIOBURN_EMIS( CH4_EMIS_ADJ )

         !4.3 Rice emissions (CH4_RIC, #7)
         CALL RICE_EMIS( CH4_EMIS_ADJ )

      ENDIF


      IF ( LASTOFYEAR ) THEN

        !4.4 Biofuel emissions (CH4_BFL, #6)
        !kjw replace with EDGARv4 biofuels in ASEASONAL_ANTHRO_EMIS
        !    (kjw, 11/17/11)
        !CALL BIOFUEL_EMIS( CH4_EMIS_ADJ )

        !4.5 Aseasonal Anthropogenic emissions
        ! (CH4_OAG, #2; CH4_COL, #3; CH4_LIV, #4; CH4_WST, #5; CH4_OTA, #8)
        CALL ASEASONAL_ANTHRO_EMIS( CH4_EMIS_ADJ )

        !4.6 Aseasonal Natural emissions (CH4_SAB, #11; CH4_OTN, #12)
        CALL ASEASONAL_NATURAL_EMIS( CH4_EMIS_ADJ )

      ENDIF


      ! Total emission: sum of all emissions - (2*soil absorption)
      ! We have to substract soil absorption twice because it is added 
      ! to other emissions in the SUM function. (ccc, 7/23/09)
      CH4_EMIS_ADJ(:,:,1) = 0d0
      CH4_EMIS_ADJ(:,:,1)=xplx(SUM(CH4_EMIS_ADJ%r,3),
     & SUM(CH4_EMIS_ADJ%i,3))-(2*CH4_EMIS_ADJ(:,:,11))

      ! Select emissions to be optimized (all but soil absorption).
      ! Exclude soil abs to prevent having negative scaling factors.
      ! To do this, add the magnitude of soil absorption back to the total
      CH4_EMIS_for_SF(:,:) = CH4_EMIS_ADJ(:,:,1) 
     &                             + CH4_EMIS_ADJ(:,:,11)

      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Accumulate gradients
      DO J = 1, JJPAR
         JREF = J + J0

         ! Get area [cm2] of each box for unit conversion
         AREA_CM2 = GET_AREA_CM2( J )

         DO I = 1, IIPAR
            IREF = I + I0

            ! Convert from [molec cm-2 s-1] --> [kg CH4]
            E_CH4 = CH4_EMIS_for_SF(I,J) * DTSRCE
     &                              * AREA_CM2 / XNUMOL_CH4

            ! Calculate Gradients
            EMS_SF_ADJ(I,J,M,ADCH4EMS) = EMS_SF_ADJ(I,J,M,ADCH4EMS) + 
     &           STT_ADJ(I,J,1,1) * E_CH4

         ENDDO
      ENDDO


      ! RESET LASTOF logicals
      ! If 12am on Jan1, next time step in adjoint will be last of preceeding year
      IF ( ITS_A_NEW_YEAR()  ) THEN
         LASTOFYEAR  = .TRUE.
      ELSE
         LASTOFYEAR  = .FALSE.
      ENDIF

      ! If 12am on 1st of month, next time step in adjoint will be last of preceeding month
      IF ( ITS_A_NEW_MONTH() ) THEN
         LASTOFMONTH = .TRUE.
      ELSE
         LASTOFMONTH = .FALSE.
      ENDIF


      WRITE(6,'(a)') '% --- EXITING EMISSCH4_ADJ! ---'

      ! Return to calling program
      END SUBROUTINE EMISSCH4_ADJ



!------------------------------------------------------------------------------

      SUBROUTINE CHEMCH4_ADJ
!
!******************************************************************************
!  Subroutine CHEMCH4_ADJ does adjoint of CH4 chemistry
!  (adj_group, kjw, 2/22/10)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,    ONLY : OH_DIR
      USE BPCH2_MOD,        ONLY : GET_TAU0,     READ_BPCH2
      USE BPCH2_MOD,        ONLY : GET_NAME_EXT, GET_RES_EXT
      USE TIME_MOD,         ONLY : GET_MONTH
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D, TRANSFER_3D
      USE ERROR_MOD,        ONLY : GEOS_CHEM_STOP
      USE GLOBAL_OH_MOD,    ONLY : GET_GLOBAL_OH, OH
      USE GLOBAL_CH4_MOD,   ONLY : CH4LOSS

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! LPAUSE

      ! Local variables
      CHARACTER(LEN=255)     :: FILENAME
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER                :: L, NOHDO, LMN, M, I, J
      TYPE (XPLEX)                 :: ARRAY(IIPAR,JJPAR,LGLOB)
      TYPE (XPLEX)                 :: XTAU


      !=================================================================
      ! CHEMCH4_ADJ begins here!
      !=================================================================
      WRITE( 6, '(a)' ) '% --- ENTERING CHEMCH4_ADJ! ---'


      ! Initialize GLOBAL_CH4_ADJ_MOD variables
      ! Do here because chemistry called before emissions and both
      ! have the same time step (kjw, 12/2/2011)
      IF ( FIRSTCHEM ) THEN
         CALL INIT_GLOBAL_CH4_ADJ

         ! Read Stratospheric loss rates
         CH4LOSS(:,:,:,:) = 0d0
         CALL READ_CH4LOSS

      ENDIF


      ! Get Average temp for the preceeding day
      !CALL CH4_AVGTP_ADJ

      !================================================================
      ! (1) get parameterized OH fields or monthly mean fields.
      !
      ! Variables of note:
      ! ---------------------------------------------------------------
      ! (1) BOH = storage array for OH fields.
      !
      ! (2) NOHDO = switch
      !       ONLY USE CASE 1 as of 5/28/08 (kjw)
      !       = 1 : Get GEOS-Chem OH (v5-07-08) (kjw, 5/28/08)
      !
      ! (3) TROPP =  the vertical level of the tropopause.  Above this
      !     level, no [OH] is calculated.  The user can feed this
      !     SR a high value for LPAUSE which effectively turns this 
      !     option off (i.e., LPAUSE > MVRTBX). If the [OH] = -999 
      !     then the [OH] was not calculated.
      !================================================================

      ! 3D OH Field
      BOH(:,:,:,:) = 0d0


      ! Change value of NOHDO as listed above
      NOHDO = 1

      SELECT CASE ( NOHDO )

         ! NOHDO = 1: GEOS-Chem OH v5-07-08
         CASE ( 1 )

            ! If first of month, read monthly mean OH
            IF ( FIRSTCHEM ) THEN

               ! Clear 3D OH field
               BOH(:,:,:,:) = 0d0
               LMN = GET_MONTH()
               
               ! Loop over each month, reading OH
               DO M=1,12
               
                  ! Global OH 
                  CALL GET_GLOBAL_OH( M )

                  ! Assign to module variable BOH
                  BOH(:,:,:,M) = OH(:,:,:)
               
               ENDDO

            ENDIF
        
         CASE DEFAULT
            WRITE( 6, '(a)' ) 'Invalid selection for NOHDO!'
            WRITE( 6, '(a)' ) 'Halting execution in CHEMCH4!'
            CALL GEOS_CHEM_STOP
            
      END SELECT


      !=================================================================
      ! (3) adjoint of CH4 chemistry in layers above tropopause.
      !=================================================================
      CALL CH4_STRAT_ADJ

      !=================================================================
      ! (3) adjoint of rate of decay of CH4 by OH oxidation.
      !=================================================================
      CALL CH4_DECAY_ADJ



      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.


      ! Return to calling program
      END SUBROUTINE CHEMCH4_ADJ


!------------------------------------------------------------------------------

      SUBROUTINE CH4_DECAY_ADJ
!
!******************************************************************************
!  Subroutine CH4_DECAY_ADJ is the adjoint of decay rate of CH4 by OH.  OH is the 
!  only sink for CH4 considered here. (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
!  We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
!  (bmy, 4/18/00)  
!
!  Monthly loss of CH4 is summed in TCH4(3)
!     TCH4(3)  = CH4 sink by OH
!
!  Module Variables:
!  ============================================================================
!  (1) BOH        (TYPE (XPLEX)) : Array holding global OH concentrations
!  (2) XNUMOL_CH4 (TYPE (XPLEX)) : Molec CH4 / kg CH4
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_DECAY is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,   ONLY : ADCH4EMS, EMS_SF_ADJ, STT_ADJ
      USE DAO_MOD,          ONLY : AIRVOL, AD, CONVERT_UNITS, T
      USE TIME_MOD,         ONLY : GET_TS_CHEM, GET_NYMD, GET_NHMS
      USE TRACER_MOD,       ONLY : TCVV, N_TRACERS
      USE TIME_MOD,         ONLY : GET_MONTH
      

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! LPAUSE

      ! Local variables
      INTEGER          :: I, J, L, LMN
      TYPE (XPLEX)           :: DT, GCH4_ADJ, STT2GCH4, KRATE, TROPCH4


      !=================================================================
      ! CH4_DECAY_ADJ begins here!
      !=================================================================

      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM() * 60d0

      ! Current month
      LMN = GET_MONTH()

      !=================================================================
      ! Compute decay of CH4 by OH in the troposphere
      !
      ! The decay for CH4 is calculated by:
      !    OH + CH4 -> CH3 + H2O 
      !    k = 2.45E-12 exp(-1775/T)
      !
      ! This is from JPL '97.
      ! JPL '00 & '06 do not revise '97 value. (jsw, kjw)
      !=================================================================

      ! Convert STT_ADJ from [v/v] --> [kg]
      CALL CONVERT_UNITS( 2,  N_TRACERS, TCVV, AD, STT_ADJ )


      TROPCH4 = 0d0

      DO L = 1, MAXVAL( LPAUSE )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( L < LPAUSE(I,J) ) THEN 

            ! Use 24-hr avg temperature to calc. rate coeff.
            ! citation needed
            KRATE = 2.45d-12 * EXP( -1775d0 / T(I,J,L) )  

            ! Conversion from [kg/box] --> [molec/cm3]
            ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4 

            ! CH4 in [molec/cm3]
            GCH4_ADJ = STT_ADJ(I,J,L,1) * STT2GCH4

            ! Calculate new CH4 value: [CH4]=[CH4](1-k*[OH]*delta) 
            GCH4_ADJ = GCH4_ADJ * ( 1d0 - KRATE * BOH(I,J,L,LMN) * DT )
		
            ! Convert back from [molec/cm3] --> [kg/box]
            STT_ADJ(I,J,L,1) = GCH4_ADJ / STT2GCH4

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Convert STT_ADJ back from [kg] --> [v/v]
      CALL CONVERT_UNITS( 1,  N_TRACERS, TCVV, AD, STT_ADJ )


      ! Return to calling program
      END SUBROUTINE CH4_DECAY_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE READ_CH4LOSS
!
!*****************************************************************************
!  Subroutine READ_CH4LOSS reads CH4 loss frequencies in the stratosphere. 
!  These values constitute a linearized stratospheric CH4 chemistry scheme.
!  Loss frequencies from 4x5 degree output from the GMI model. Thanks to Lee 
!  Murray for the ch4 loss frequencies. (kjw, 11/19/2011)
!
!  Module Variables:
!  ===========================================================================
!  (1) CH4LOSS (TYPE (XPLEX)) : Array containing ch4 loss frequencies for all 12 months [1/s]
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_CH4LOSS is independent of "F77_CMN_OH", "F77_CMN_CO", and "F77_CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JJPAR,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Treat MERRA in the same way as for GEOS-5 (bmy, 8/13/10)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT,    GET_MODELNAME
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
#     include "define.h"
#     include "CMN_SIZE"

      ! Local variables
      INTEGER            :: I, J, L, M
      TYPE (XPLEX)             :: ARRAY(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)             :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_CH4LOSS begins here!
      ! 
      ! Read P(CO) for all 12 months
      !=================================================================
      ! Construct filename
!      FILENAME = TRIM( DATA_DIR ) // ...
      FILENAME = '/data/gc/CH4/gmi.ch4loss.' //
     &           'geos5_47L.' // get_res_ext() // '.bpch'
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = '/data/gc/CH4/gmi.ch4loss.' //
     &           'geos5_47L.05x0666_NA.bpch'
#endif

      WRITE( 6, 93 ) TRIM ( FILENAME )
 93   FORMAT( '     - READ_CH4LOSS: Reading Ch4loss: ', a )
      CALL FLUSH( 6 )

      ! Read data for each month
      DO M = 1, 12 

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Read Loss frequencies in units of [1/s].  drevet.
         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,     
     &                    XTAU,      IIPAR,     JJPAR,     
     &                    LLPAR,     ARRAY,     QUIET=.TRUE. )

         ! Place array into CH4LOSS module variable
         CH4LOSS(:,:,:,M) = ARRAY(:,:,:)

      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_CH4LOSS

!------------------------------------------------------------------------------

      SUBROUTINE CH4_STRAT_ADJ
!
!*****************************************************************************
!  Subroutine CH4_STRAT_ADJ is adjonit of loss of CH4 above tropopause. 
!
!  Production (mixing ratio/sec) rate provided by Dylan Jones.  
!  Only production by CH4 + OH is considered.
!  
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
! 
!  Levels           1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/18/00)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_STRAT is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed LMN from the arg list and made it a local variable.  Now use 
!        functions GET_MONTH and GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : AIRVOL
      USE TIME_MOD,   ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD, ONLY : STT
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE

      ! Local variables
      INTEGER             :: I, J, L, LMN
      TYPE (XPLEX)              :: DT, GCH4, STT2GCH4, LRATE
      CHARACTER*20        :: STT_TEST
      CHARACTER*20        :: STT2GCH4_CHAR

      ! External functions
      TYPE (XPLEX), EXTERNAL    :: BOXVL

      !=================================================================
      ! CH4_STRAT_ADJ begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DT  = GET_TS_CHEM() * 60d0

      ! Current month
      LMN = GET_MONTH()

      !=================================================================
      ! Loop over stratospheric boxes only
      !=================================================================
      DO L = MINVAL( LPAUSE ), LLPAR 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( L >= LPAUSE(I,J) ) THEN

            ! Conversion factor [kg/box] --> [molec/cm3]
            ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4

            ! CH4 in [molec/cm3]
            GCH4 = STT_ADJ(I,J,L,1) * STT2GCH4

            ! Loss rate [molec/cm3/s]
            LRATE = GCH4 * CH4LOSS( I,J,L,LMN )

            ! CH4 in [molec/cm3]
            GCH4 = GCH4 - ( LRATE * DT )

!kjw. Update stratospheric chem to use linearized CH4 loss frequencies
!  (kjw, 11/19/11)
!            ! Sum loss in TCH4(3) [molec CH4/box] in the stratosphere
!            ! [molec/cm3] * [v/v/s] * [s] * [cm3/box] = [molec CH4/box]
!            TCH4(I,J,L,3) = TCH4(I,J,L,3) + 
!     &                      ( BAIRDENS(I,J,L) * COPROD(J,L,LMN) *
!     &                        DT              * BOXVL(I,J,L)    )
!
!            ! Calculate new CH4 value [molec CH4/cm3] in the stratosphere
!            ! [v/v/s] * [s] * [molec/cm3] = [molec CH4/cm3] 
!            GCH4 = GCH4 - ( COPROD(J,L,LMN) * DT * BAIRDENS(I,J,L) )
!kjw

            ! Convert back from [molec CH4/cm3] --> [kg/box] 
            STT_ADJ(I,J,L,1) = GCH4 / STT2GCH4


!kjw. With new linearized chemistry, STT should never be negative
!  (kjw, 11/19/11)		
!	    IF ( STT(I,J,L,1) < 0 ) THEN
!		STT(I,J,L,1)=0
!	    ENDIF
!kjw

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CH4_STRAT_ADJ

!------------------------------------------------------------------------------

!      SUBROUTINE CH4_STRAT
!!
!!*****************************************************************************
!!  Subroutine CH4_STRAT calculates uses production rates for CH4 to 
!!  calculate loss of CH4 in above the tropopause. 
!!  (jsw, bnd, bmy, 1/16/01, 7/20/04)
!!  DO NOT CALL IN ADJOINT SIMULATION (kjw, 7/6/11)
!!
!!  Production (mixing ratio/sec) rate provided by Dylan Jones.  
!!  Only production by CH4 + OH is considered.
!!  
!!  The annual mean tropopause is stored in the LPAUSE array 
!!  (from header file "CMN").  LPAUSE is defined such that: 
!! 
!!  Levels           1 <= L <  LPAUSE(I,J) - 1 are tropospheric
!!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/18/00)
!!
!!  NOTES:
!!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!!        by Bob Yantosca. (bmy, 1/16/01)
!!  (2 ) CH4_STRAT is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!!        (bmy, 1/16/01)
!!  (3 ) Removed LMN from the arg list and made it a local variable.  Now use 
!!        functions GET_MONTH and GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!!*****************************************************************************
!!
!      ! References to F90 modules
!      USE ADJ_ARRAYS_MOD,   ONLY : ADCH4EMS,  EMS_SF_ADJ,  STT_ADJ
!      USE DAO_MOD,          ONLY : AIRVOL,    AD
!      USE TIME_MOD,         ONLY : GET_MONTH, GET_TS_CHEM
!
!#     include "CMN_SIZE"       ! Size parameters
!#     include "CMN"            ! LPAUSE
!
!      ! Local variables
!      LOGICAL, SAVE          :: FIRSTCHEM
!      INTEGER                :: I, J, L, LMN
!      TYPE (XPLEX)                 :: DT, GCH4_ADJ, STT2GCH4
!      TYPE (XPLEX), PARAMETER      :: WTAIR = 28.966d0
!
!      ! External functions
!      TYPE (XPLEX), EXTERNAL       :: BOXVL
!
!      !=================================================================
!      ! CH4_STRAT begins here!
!      !=================================================================
!
!      !=================================================================
!      ! (1) If first time step, read LCO data
!      !=================================================================
!      IF ( FIRSTCHEM ) THEN
!
!         ! Zero CO Production array
!         COPROD(:,:,:) = 0d0
!
!         ! Read zonally-averaged CO production [v/v/s]
!         CALL READ_COPROD
! 
!      ENDIF
!
!      ! Chemistry timestep [s]
!      DT  = GET_TS_CHEM() * 60d0
!
!      ! Current month
!      LMN = GET_MONTH()
!
!
!      !=================================================================
!      ! (2) Calculate each box's air density [molec/cm3]
!      !=================================================================
!
!      DO L = 1, LLPAR
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!         BAIRDENS(I,J,L) = AD(I,J,L) * 1000d0   / BOXVL(I,J,L) * 
!     &                                 6.023D23 / WTAIR
!      ENDDO
!      ENDDO
!      ENDDO
!
!
!      !=================================================================
!      ! (3) Calculate stratospheric CH4 loss from COPRODuction
!      !        Loop over stratospheric boxes only
!      !=================================================================
!      DO L = MINVAL( LPAUSE ), LLPAR 
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!         IF ( L >= LPAUSE(I,J) ) THEN
!
!            ! Conversion factor [kg/box] --> [molec/cm3]
!            ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
!            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4
!
!            ! CH4 in [molec/cm3]
!            GCH4_ADJ = STT_ADJ(I,J,L,1) * STT2GCH4
!
!            ! Calculate new CH4 value [molec CH4/cm3] in the stratosphere
!            ! [v/v/s] * [s] * [molec/cm3] = [molec CH4/cm3] 
!            GCH4_ADJ = GCH4_ADJ - 
!     &           ( COPROD(J,L,LMN) * DT * BAIRDENS(I,J,L) )
!
!            ! Convert back from [molec CH4/cm3] --> [kg/box] 
!            STT_ADJ(I,J,L,1) = GCH4_ADJ / STT2GCH4
!
!         ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!
!      ! Set FIRSTCHEM to FALSE
!      FIRSTCHEM = .FALSE.
!
!      ! Return to calling program
!      END SUBROUTINE CH4_STRAT
!
!
!
!------------------------------------------------------------------------------

      SUBROUTINE READ_COPROD
!
!*****************************************************************************
!  Subroutine READ_COPROD reads production and destruction rates for CO in 
!  the stratosphere. (bnd, bmy, 1/17/01, 10/3/05)
!
!  Module Variables:
!  ===========================================================================
!  (1) COPROD (TYPE (XPLEX)) : Array containing P(CO) for all 12 months [v/v/s]
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_COPROD is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JGLOB,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT,    GET_MODELNAME
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_ZONAL
        
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters


      ! Local variables
      CHARACTER(LEN=255) :: FILENAME
      INTEGER            :: M
      TYPE (XPLEX)             :: ARRAY(1,JGLOB,LGLOB)
      TYPE (XPLEX)             :: DUMMY_IN(JGLOB,LGLOB)
      TYPE (XPLEX)             :: XTAU
      TYPE (XPLEX)             :: DUMMY_OUT(JGLOB,LGLOB)


      !=================================================================
      ! READ_COPROD begins here!
      ! 
      ! Read P(CO) for all 12 months
      !=================================================================
      DO M = 1, 12 

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Construct filename
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/' //
     &   'COprod.' // GET_NAME_EXT() // '.' // GET_RES_EXT()


         WRITE( 6, 93 ) TRIM ( FILENAME )
 93      FORMAT( '     - READ_COPROD: Reading COprod: ', a )
	 CALL FLUSH( 6 )


cdrevet
         ! Read P(CO) in units of [v/v/s]
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,     
     &                    XTAU,      1,         JGLOB,     
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )
cdrevet

         ! use 2D arrays for TRANSFER ZONAL
         DUMMY_IN(:,:) = ARRAY(1,:,:)

         ! Copy TYPE (XPLEX) to COMPLEX*16 data, and resize from (JGLOB,LGLOB) 
         ! to (JJPAR,LLPAR) -- vertically regrid if necessary
         CALL TRANSFER_ZONAL( DUMMY_IN, DUMMY_OUT )

         COPROD(:,:,M) = DUMMY_OUT(:,:)

      ENDDO


      ! Return to calling program
      END SUBROUTINE READ_COPROD


!------------------------------------------------------------------------------

      SUBROUTINE CH4_AVGTP_ADJ
!
!******************************************************************************
!  Subroutine CH4_AVGTP gets the 24-h average surface pressure and temperature
!  needed for the CH4 simulation. (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry and
!        placed into module "global_ch4_mod.f" by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_AVGTP is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed duplicate definition for NTDT, NMIN (bmy, 11/15/01)
!  (4 ) Removed PS from argument list.  Now use P(I,J)+PTOP instead of
!        PS, this ensures that we have consistency between P and AD.
!        (bmy, 4/11/02)
!  (5 ) Removed obsolete code (bmy, 6/27/02)
!  (6 ) Now uses GET_PCENTER from "pressure_mod.f" to return the pressure
!        at the midpoint of the box (I,J,L).  Also added parallel DO-loops.
!        Updated comments. (dsa, bdf, bmy, 8/21/02)
!  (7 ) Now reference T from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f" (bmy, 10/15/02)
!  (8 ) Removed NTDT, NMIN from the arg list.  Now uses functions GET_TS_DYN,
!        GET_TS_CHEM, and GET_ELAPSED_MIN from "time_mod.f" (bmy, 3/27/03)
!  (9 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!  (10) Modify from CH4_AVGTP_ADJ for compatability with the adjoint (kjw, 2/22/10)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,     ONLY : GET_TS_DYN, GET_TS_CHEM
      USE TIME_MOD,     ONLY : GET_NYMD,   GET_NYMDe,   GET_NYMDb
      USE TIME_MOD,     ONLY : GET_NHMS,   GET_TIME_BEHIND_ADJ
      USE DAO_MOD,      ONLY : T

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      CHARACTER(LEN=255)  :: PATH
      INTEGER             :: I, NYMD, NHMS, NDYN, TAG, INFO(3)
      INTEGER             :: NYMDp, NHMSp, countit
      TYPE (XPLEX)              :: Temp(IIPAR, JJPAR, LLPAR)
      TYPE (XPLEX)              :: result(2)
      
      !=================================================================
      ! CH4_AVGTP_ADJ begins here!
      !=================================================================

      WRITE(6,'(a)') '   % CH4_AVGTP_ADJ begins'

      ! Initialize Tavg_adj for the current time step
      Tavg_adj(:,:,:) = 0d0

      ! NDYN = # of dynamics time steps in each chemical time step
      NDYN = GET_TS_CHEM() / GET_TS_DYN()

      ! If 1 day from beginning of the simulation
      IF ( GET_NYMD() .EQ. GET_NYMDb()+1 ) NDYN = NDYN + 1

      countit=0

      DO I = 1,NDYN

         ! Initialize T
         !T(:,:,:) = 0d0

         ! Get time stamps at current and every dyn time step during the day
         result = GET_TIME_BEHIND_ADJ( GET_TS_DYN()*(I-1) )

         ! NYMD, NHMS, TAG
         NYMD = result(1)
         NHMS = result(2)

         ! Get file unit number for A6 file to be opened
         ! Such that 20 <= IU_NUM <= 64 .OR. IU_NUM => 100 (kjw, 2/23/10)
         ! Check available unit numbers in file_mod.f
         ! IU_NUM = 52

         ! Find date of closest A-6 file to open and which occurence of Temp to use
         INFO  = FIND_CLOSEST_A6( NYMD, NHMS )
         NYMDp = INFO(1)
         NHMSp = INFO(2)
         TAG   = INFO(3)


         ! Open A6 file to read Temperature data
         CALL OPEN_A6_CH4_ADJ( NYMDp, NHMSp )


         ! If the desired A-6 file is already in use, use DAO_MOD, ONLY : T
         IF ( IU_A6_CH4_ADJ == 72 ) THEN

            Tavg_adj(:,:,:) = Tavg_adj(:,:,:) + T(:,:,:)
            countit=countit+1

         ELSE IF ( IU_A6_CH4_ADJ == 52 ) THEN

            ! READ A6 fields with temp data
            CALL READ_A6_CH4_ADJ( NYMDp, NHMSp, TAG, Temp )

            ! Collect temperature data
            Tavg_adj(:,:,:) = Tavg_adj(:,:,:) + Temp(:,:,:)
            countit=countit+1

            ! Close file we just opened
            CLOSE( IU_A6_CH4_ADJ )

         ENDIF

      ENDDO

      ! Average Temperature information
      Tavg_adj(:,:,:) = Tavg_adj(:,:,:) / NDYN


      WRITE(6,'(a)') '   % CH4_AVGTP_ADJ ends'

      ! Return to calling program
      END SUBROUTINE CH4_AVGTP_ADJ



!------------------------------------------------------------------------------

!      SUBROUTINE UPDATE_LASTOF
!
!********************************************************************************

! Subroutine UPDATE_LASTOF determines whether the next time step will be the last
! of a month (kjw, 2/22/10)
!
! NOTES
!
!********************************************************************************


      ! If 12am on Jan1, next time step in adjoint will be last of preceeding year
!      IF ( ITS_A_NEW_YEAR()  ) THEN
!         LASTOFYEAR  = .TRUE.
!      ELSE
!         LASTOFYEAR  = .FALSE.
!      ENDIF

      ! If 12am on 1st of month, next time step in adjoint will be last of preceeding month
!      IF ( TIS_A_NEW_MONTH() ) THEN
!         LASTOFMONTH = .TRUE.
!      ELSE
!         LASTOFMONTH = .TRUE.
!      ENDIF



      ! Return to calling program
!      END SUBROUTINE UPDATE_LASTOF




!------------------------------------------------------------------------------

      SUBROUTINE OPEN_A6_CH4_ADJ( NYMDp, NHMSp )
!
!******************************************************************************
!  Subroutine OPEN_A6_CH4_ADJ opens the A-6 met fields file for date NYMD and 
!  time NHMS for use by CH4 adjoint simulation. Based on GET_A6_FIELDS.
!  (bmy, bdf, 6/15/98, 2/12/09), (kjw, 2/23/10)
!
!  Difference with OPEN_A6_FIELDS is that this uses a different file unit
!  number than IU_A6. File unit # is a parameter, IU_A6_CH4_ADJ
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER)   : Current value of YYYYMMDD
!  (2 ) NHMS (INTEGER)   : Current value of HHMMSS
!
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (4 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (5 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit IU_A6 
!        refers to a valid file on disk (bmy, 3/23/05)
!  (6 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (9 ) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!  (10) Set N_A6_FIELDS=21 for GEOS-5 and IN_CLOUD_OD (jmao, bmy, 2/12/09)
!******************************************************************************


      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LUNZIP
      USE FILE_MOD,      ONLY : IOERROR, FILE_EXISTS
      USE TIME_MOD,      ONLY : EXPAND_DATE, GET_NYMD, GET_NHMS

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMDp, NHMSp

      ! Local variables
      INTEGER                :: IOS
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: A6_FILE
      CHARACTER(LEN=255)     :: A6_NOW
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_A6_FIELDS begins here!
      !=================================================================


      ! Get Filename
      ! ----------------------------------------------------------------
#if defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         A6_FILE  = 'YYYYMMDD.a6.' // GET_RES_EXT()
         A6_NOW   = 'YYYYMMDD.a6.' // GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         A6_FILE  = 'YYYYMMDD.a6.' // GET_RES_EXT()
         A6_NOW   = 'YYYYMMDD.a6.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMDp, NHMSp )
         CALL EXPAND_DATE( A6_FILE,  NYMDp, NHMSp )
         CALL EXPAND_DATE( A6_NOW,   GET_NYMD(), GET_NHMS() )

         print*,'NYMD and NHMS right now ',GET_NYMD(),GET_NHMS()
         print*,'NYMDp and NHMSp of file to open ', NYMDp,NHMSp

         ! If the A-6 file is already open, return to calling program
         IF ( TRIM( A6_FILE ) == TRIM( A6_NOW ) ) THEN
            print*,'This file is already open',TRIM(A6_NOW)
            print*,'Using previously opened A-6 file...'
            IU_A6_CH4_ADJ = 72
            RETURN
         ELSE
            print*,'We have to open a new A-6 file: ',TRIM(A6_NOW)
            IU_A6_CH4_ADJ = 52
         ENDIF

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( A6_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( A6_FILE )
         ENDIF

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_A6_CH4_ADJ ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_A6_FIELDS (a6_read_mod.f)' )
         ENDIF

      ! Open the file
      ! ----------------------------------------------------------------
         ! Hardwire unit number to not conflict with current IU_A6 (kjw, 2/23/10)
         OPEN( UNIT   = IU_A6_CH4_ADJ, FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A6_CH4_ADJ, 'open_a6_fields:1' )
         ENDIF

         ! Skip past the ident string
         READ( IU_A6_CH4_ADJ, IOSTAT=IOS ) IDENT

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a ) 


      ! Return to calling program
      END SUBROUTINE OPEN_A6_CH4_ADJ


!-----------------------------------------------------------------------------

      SUBROUTINE READ_A6_CH4_ADJ( NYMDp, NHMSp, TAG, Temp )
!
!******************************************************************************
!  Subroutine READ_A6_CH4_ADJ reads A-6 (avg 6-hr) met fields from disk. 
!  (bmy, 6/5/98, 3/28/08)
!
!  For CH4 adjoint simulation, hardwire file unit # as a parameter
! 
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD     : YYYYMMDD
!  (2 ) NHMS     :  and HHMMSS of A-6 met fields to be accessed
!
!  A-6 Met Fields as Output (Optional Arguments):
!  ============================================================================
!  (1) T         : (3-D) Temperature                           [K]
!
!  NOTES:
!  (1 ) Adapted from READ_A6 of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!  (3 ) Now compute CLDTOPS using ZMMU for GEOS-4 (bmy, 3/4/04)
!  (4 ) Now modified for GEOS-5 and GCAP fields.  Added DETRAINE, 
!        DETRAINN, DNDE, DNDN, ENTRAIN, UPDE, UPDN as optional arguments.
!        Now references "CMN_DIAG". (swu, bmy, 5/25/05)
!  (5 ) Bug fix in ND66 diagnostic for GEOS-4 (bmy, 2/1/06)
!  (6 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (7 ) Now set negative SPHU to a small positive # (1d-32) instead of zero, 
!        so as not to blow up logarithms (bmy, 9/8/06)
!  (8 ) Add CMFMC, DQIDTMST, DQLDTMST, DQRCON, DQRLSC, DQVDTMST, MFXC, MFYC, 
!        MFZ, PLE, PV, RH, TAUCLI, and TAUCLW as optional arguments.  Also 
!        update the CASE statement accordingly for GEOS-5 met fields. 
!        Now reference TRANSFER_3D_Lp1 from "transfer_mod.f".  Now convert
!        GEOS-5 specific humidity from [kg/kg] to [g/kg] for compatibility
!        with existing routines.  Also recognize EPV, which is an alternate 
!        name for PV.  Bug fix: convert GEOS-5 RH from unitless to %.
!        (phs, bmy, 3/28/08)
!  (8 ) Now get the # of A-6 fields from the file ident string (bmy, 10/7/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD66,        AD67
      USE FILE_MOD,     ONLY : IOERROR
      USE TIME_MOD,     ONLY : SET_CT_A6,   TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_A6, TRANSFER_3D_Lp1
      USE TRANSFER_MOD, ONLY : TRANSFER_3D, TRANSFER_G5_PLE

#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! ND66, ND67
#     include "CMN_GCTM"             ! g0

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMDp, NHMSp, TAG
      TYPE (XPLEX),  INTENT(OUT)           :: Temp(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                        :: IOS
      TYPE (XPLEX)                         :: D(IGLOB,JGLOB,LGLOB)
      CHARACTER(LEN=8)               :: NAME
      INTEGER                        :: XYMD, XHMS, NFOUND

      !=================================================================
      ! READ_A6 begins here!      
      !=================================================================

      ! Number of A-6 fields.
      ! We only want 1: temperature
      !N_A6 = 1  DON'T NEED

      ! Zero number of fields that we have found
      !NFOUND = 0   DON'T NEED

      !=================================================================
      ! Read the A-6 fields from disk
      !=================================================================

      ! Count # of times we find temperature
      NFOUND = 0

      ! Read each available data set in the file, but only save Temperature
      DO

         ! Read A-6 field name
         READ( IU_A6_CH4_ADJ, IOSTAT=IOS ) NAME
         !print*, 'A6 NAME : ', NAME

         ! IOS < 0: End-of-file; make sure we've found 
         ! all the A-6 fields before exiting this loop
         IF ( IOS < 0 ) EXIT

         ! IOS > 0: True I/O Error, stop w/ error msg 
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_A6_CH4_ADJ, 
     &                               'read_a6_ch4_adj:1' )

         ! Examine the name of field
         SELECT CASE ( TRIM( NAME ) )

            ! If we've found temperature
            CASE( 'T' )

               ! Increase count
               NFOUND = NFOUND + 1

               IF ( NFOUND == TAG ) THEN

                  print*,'% --- READ_A6_CH4_ADJ : Found T field desired'

                  ! Read into array
                  READ( IU_A6_CH4_ADJ, IOSTAT=IOS ) XYMD, XHMS, D
                  IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6_CH4_ADJ,  
     &                                      'read_a6_ch4_adj:29' )

                  !IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  !   IF ( PRESENT( T ) ) CALL TRANSFER_3D( D, T )
                  CALL TRANSFER_3D( D, Temp )
                  !ENDIF

                  ! Return to Calling Program
                  RETURN

               ELSE

                  print*,' % --- READ_A6_CH4_ADJ : Not Correct T field'
                  
                  ! Read into array
                  READ( IU_A6_CH4_ADJ, IOSTAT=IOS ) XYMD, XHMS, D
                  IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6_CH4_ADJ,  
     &                                      'read_a6_ch4_adj:29' )

               ENDIF

            ! If field is not temperature -- skip over
            CASE DEFAULT
               !WRITE( 6, '(a)' ) 'Searching for next A-6 field!'
               !WRITE( 6, '(2a)' ) 'THIS name = ',TRIM(NAME)
               !print*,'LLPAR = ',LLPAR
               !print*,'LGLOB = ',LGLOB
               READ( IU_A6_CH4_ADJ, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6_CH4_ADJ,  
     &                                      'read_a6_ch4_adj:41' )

         END SELECT

      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_A6_CH4_ADJ



!------------------------------------------------------------------------------

      FUNCTION FIND_CLOSEST_A6( NYMD, NHMS ) RESULT( INFO )
!
!********************************************************************************
! Subroutine FIND_CLOSEST_A6 finds the date and time of the nearest A-6 met field
!  (kjw, 2/24/10)
!
! NOTES
!    time tag, A6_TIME(3), tells us which occurence of A-6 Temperature we want (1st-4th)
!
!********************************************************************************

      ! Reference to f90 modules
      USE TIME_MOD,         ONLY : GET_NYMD, GET_NHMS


#     include "CMN_SIZE" ! Size

      ! Arguments
      INTEGER, INTENT(IN)   :: NYMD, NHMS

      ! Output
      INTEGER               :: INFO(3)

      ! Local Variables
      CHARACTER(LEN=6)      :: CNHMS
!      CHARACTER(LEN=2)      :: CHH
!      CHARACTER(LEN=2)      :: CMM
!      INTEGER               :: HH,MM,NMIN


      !============================================================
      ! FIND_CLOSEST_A6 begins here!
      !============================================================

      ! If 12-2:59am
      IF ( ( NHMS >= 000000 ) .AND. ( NHMS <= 025959 ) ) THEN

         ! NYMD remains the same
         INFO(1) = NYMD

         ! NHMS = 6am
         INFO(2) = 000000

         ! Set time tag
         INFO(3) = 1

      ! If 3-8:59am
      ELSE IF ( ( NHMS >= 030000 ) .AND. ( NHMS <= 085959 ) ) THEN

         ! NYMD remains the same
         INFO(1) = NYMD

         ! NHMS = 6am
         INFO(2) = 060000

         ! Set time tag
         INFO(3) = 2

      ! If 9am-2:59pm
      ELSE IF ( ( NHMS >= 090000 ) .AND. ( NHMS <= 145959 ) ) THEN

         ! NYMD remains the same
         INFO(1) = NYMD

         ! NHMS = 12pm
         INFO(2) = 120000

         ! Set time tag
         INFO(3) = 3

      ! If 3pm-8:59pm
      ELSE IF ( ( NHMS >= 150000 ) .AND. ( NHMS <= 205959 ) ) THEN

         ! NYMD remains the same
         INFO(1) = NYMD

         ! NHMS = 12pm
         INFO(2) = 180000

         ! Set time tag
         INFO(3) = 4

      ! If 9pm-11:59pm
      ELSE IF ( ( NHMS >= 210000 ) .AND. ( NHMS <= 235959 ) ) THEN

         ! Since calling at midnight, these values should be current time
         INFO(1) = GET_NYMD()
         INFO(2) = GET_NHMS()

         ! Set time tag
         INFO(3) = 1

         ! We should find how many minutes behind current time is NYMD and NHMS 
         ! Turn NHMS into a string
         !WRITE( CNHMS, '(i6)' ) NHMS

         ! Get Hour and Minute values from this
         !CHH = CNHMS(1:2)
         !CMM = CNHMS(3:4)
         !READ( CHH, * ) HH
         !READ( CMM, * ) MM

         ! Get number of minutes from midnight
         !NMIN = 60 * (24 - HH) + (60 - MM) - 60

         ! Get proper date stamp for midnight on the next day
         !INFO = GET_TIME_AHEAD( NMIN )

      ENDIF
      RETURN


      ! Return to calling program
      END FUNCTION FIND_CLOSEST_A6


!------------------------------------------------------------------------------

      FUNCTION GET_SCALE_GROUP( ) RESULT( CURRENT_GROUP )
!
!********************************************************************************
! Subroutine GET_SCALE_GROUP determines which predifined scaling index corresponds
! to the current time and location  (dkh, 12/02/04)
!
! NOTES
! (1 ) CURRENT_GROUP is currently only a function of TAU
! (2 ) Get rid of I,J as argument. (dkh, 03/28/05)
!
!********************************************************************************

      ! Reference to f90 modules
      USE TIME_MOD, ONLY : GET_TAU, GET_TAUe, GET_TAUb, GET_MONTH
      USE ADJ_ARRAYS_MOD, ONLY: MMSCL

#     include "CMN_SIZE" ! Size stuff

      ! Arguments
      INTEGER      :: I, J

      ! Local Variables
      TYPE (XPLEX)       :: TOTAL_HR, CURRENT_HR, GROUP_LENGTH
      TYPE (XPLEX)       :: TAU, TAUe, TAUb

      ! Function variable
      INTEGER      :: CURRENT_GROUP
      LOGICAL, SAVE :: MONTHLY = .TRUE.
      INTEGER, SAVE :: MONTH_SAVE
      INTEGER, SAVE :: GROUP_SAVE
      LOGICAL, SAVE :: FIRST = .TRUE.

      !============================================================
      ! GET_SCALE_GROUP begins here!
      !============================================================

      ! Currently there is no spatial grouping

      ! Determine temporal grouping
      IF ( MMSCL == 1 ) THEN
         CURRENT_GROUP = 1
         RETURN
      ENDIF

      IF ( MONTHLY ) THEN
         IF (FIRST) THEN 
            MONTH_SAVE = GET_MONTH()
            CURRENT_GROUP = MMSCL
            GROUP_SAVE = MMSCL
            FIRST = .FALSE.
         ENDIF
         IF ( MONTH_SAVE /= GET_MONTH() ) THEN 
            MONTH_SAVE = GET_MONTH()
            GROUP_SAVE = GROUP_SAVE - 1
            CURRENT_GROUP = GROUP_SAVE
         ELSE
               CURRENT_GROUP = GROUP_SAVE
         ENDIF

      ELSE
         ! Retrieve time parameters
         TAUe       = GET_TAUe()
         TAUb       = GET_TAUb()
         TAU        = GET_TAU()
         TOTAL_HR   = TAUe - TAUb
         CURRENT_HR = TAU  - TAUb


         ! The last time step always belongs to the last group
         IF ( TAU == TAUe ) THEN
            CURRENT_GROUP = MMSCL
            RETURN
         ELSE    

            ! Determine the length of each group
            GROUP_LENGTH = ( TOTAL_HR / MMSCL )

            ! Index is the current time divided by the group length, plus one
            CURRENT_GROUP = ( CURRENT_HR / GROUP_LENGTH ) + 1

         ENDIF

      ENDIF

      END FUNCTION GET_SCALE_GROUP


!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_CH4_ADJ
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_CH4 allocates and zeroes module arrays. 
!  (bmy, 1/16/01, 10/15/02)
!
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!      
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"
#     include "CMN_DIAG"
      
      ! Local variables
      INTEGER         :: AS
      LOGICAL, SAVE   :: FIRST = .TRUE.


      ! If NOT first, return
      IF ( FIRST==.FALSE. ) RETURN


      ALLOCATE( BAIRDENS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BAIRDENS' )
      BAIRDENS = 0d0

      ALLOCATE( BOH( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BOH' )
      BOH = 0d0

      ALLOCATE( CH4LOSS( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4LOSS' )
      CH4LOSS = 0d0

      ALLOCATE( COPROD( JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
      COPROD = 0d0

      ALLOCATE( TAVG_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAVG_ADJ' )
      TAVG_ADJ = 0d0     

      ALLOCATE( CH4_EMIS_ADJ( IIPAR, JJPAR, PD58), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_EMIS_ADJ' )
      CH4_EMIS_ADJ = 0d0


      ! We've now initialized, do not attempt again!
      FIRST = .FALSE.

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_CH4_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_CH4_ADJ
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_CH4 deallocates module arrays. (bmy, 1/16/01)
!******************************************************************************


      IF ( ALLOCATED( BAIRDENS  ) ) DEALLOCATE( BAIRDENS  )
      IF ( ALLOCATED( BOH       ) ) DEALLOCATE( BOH       )
      IF ( ALLOCATED( CH4LOSS   ) ) DEALLOCATE( CH4LOSS   )
      IF ( ALLOCATED( COPROD    ) ) DEALLOCATE( COPROD    )
      IF ( ALLOCATED( Tavg_adj  ) ) DEALLOCATE( Tavg_adj  )
      IF ( ALLOCATED( CH4_EMIS_ADJ  ) ) DEALLOCATE( CH4_EMIS_ADJ  )


      END SUBROUTINE CLEANUP_GLOBAL_CH4_ADJ

!------------------------------------------------------------------------------
     

      ! End of module
      END MODULE GLOBAL_CH4_ADJ_MOD
