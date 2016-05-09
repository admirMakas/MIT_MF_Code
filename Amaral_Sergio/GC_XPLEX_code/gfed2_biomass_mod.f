! $Id: gfed2_biomass_mod.f,v 1.5 2010/05/07 20:39:47 daven Exp $
      MODULE GFED2_BIOMASS_MOD
!
!******************************************************************************
!  Module GFED2_BIOMASS_MOD contains variables and routines to compute the
!  GFED2 biomass burning emissions. (psk, bmy, 4/20/06, 7/8/09)
!
!     Monthly/8-day/3-hr emissions of C are read from disk and then
!     multiplied by the appropriate emission factors to produce biomass
!     burning emissions on a "generic" 1x1 grid.  The emissions are then
!     regridded to the current GEOS-Chem or GCAP grid (1x1, 2x25, or
!     4x5).
!     If several GFED2 options are switched on, the smaller period
!     product is used: 3-hr before 8-day before monthly.
!
!  GFED2 biomass burning emissions are computed for the following gas-phase 
!  and aerosol-phase species:
!
!     (1 ) NOx  [  molec/cm2/s]     (9 ) CH2O [  molec/cm2/s]
!     (2 ) CO   [  molec/cm2/s]     (10) C2H6 [atoms C/cm2/s]                  
!     (3 ) ALK4 [atoms C/cm2/s]     (11) SO2  [  molec/cm2/s]    
!     (4 ) ACET [atoms C/cm2/s]     (12) NH3  [  molec/cm2/s]    
!     (5 ) MEK  [atoms C/cm2/s]     (13) BC   [atoms C/cm2/s]  
!     (6 ) ALD2 [atoms C/cm2/s]     (14) OC   [atoms C/cm2/s]     
!     (7 ) PRPE [atoms C/cm2/s]     (15) CO2  [  molec/cm2/s]
!     (8 ) C3H8 [atoms C/cm2/s]   
!
!  Module Variables:
!  ============================================================================
!  (1 ) IDBNOx          (INTEGER) : Local index for NOx  in BIOM_OUT array
!  (2 ) IDBCO           (INTEGER) : Local index for CO   in BIOM_OUT array
!  (3 ) IDBALK4         (INTEGER) : Local index for ALK4 in BIOM_OUT array
!  (4 ) IDBACET         (INTEGER) : Local index for ACET in BIOM_OUT array
!  (5 ) IDBMEK          (INTEGER) : Local index for MEK  in BIOM_OUT array
!  (6 ) IDBALD2         (INTEGER) : Local index for ALD2 in BIOM_OUT array
!  (7 ) IDBPRPE         (INTEGER) : Local index for PRPE in BIOM_OUT array
!  (8 ) IDBC3H8         (INTEGER) : Local index for C3H8 in BIOM_OUT array
!  (9 ) IDBCH2O         (INTEGER) : Local index for CH2O in BIOM_OUT array
!  (10) IDBC2H6         (INTEGER) : Local index for C2H6 in BIOM_OUT array
!  (11) IDBSO2          (INTEGER) : Local index for SO2  in BIOM_OUT array
!  (12) IDBNH3          (INTEGER) : Local index for NH3  in BIOM_OUT array
!  (13) IDBBC           (INTEGER) : Local index for BC   in BIOM_OUT array
!  (14) IDBOC           (INTEGER) : Local index for OC   in BIOM_OUT array
!  (15) IDBCO2          (INTEGER) : Local index for CO2  in BIOM_OUT array
!  (11) SECONDS         (TYPE (XPLEX) ) : Number of seconds in the current month
!  (12) N_EMFAC         (INTEGER) : Number of emission factors per species
!  (13) N_SPEC          (INTEGER) : Number of species
!  (14) VEG_GEN_1x1     (TYPE (XPLEX) ) : Array for GFED2 1x1 vegetation ma
!  (15) GFED2_SPEC_NAME (CHAR*4 ) : Array for GFED2 biomass species names
!  (16) GFED2_EMFAC     (TYPE (XPLEX) ) : Array for user-defined emission factors
!  (17) BIOM_OUT        (TYPE (XPLEX) ) : Array for biomass emissions on model grid
!  (18) DOY8DAY         (INTEGER) : Day Of the Year at start of the current
!                                   8-day period. 
!  (19) T3HR            (INTEGER) : HH at start of the current 3-hr period.
!  (20) UPDATED         (LOGICAL) : flag to indicate if new data are read at
!                                   the current emission time step.
!
!  Module Routines:
!  ============================================================================
!  (1 ) GFED2_COMPUTE_BIOMASS     : Computes biomass emissions once per month
!  (2 ) GFED2_SCALE_FUTURE        : Applies IPCC future scale factors to GFED2
!  (3 ) GFED2_TOTAL_Tg            : Totals GFED2 biomass emissions [Tg/month]
!  (4 ) INIT_GFED2_BIOMASS        : Initializes arrays and reads startup data
!  (5 ) CLEANUP_GFED2_BIOMASS     : Deallocates all module arrays
!
!  GEOS-Chem modules referenced by "gfed2_biomass_mod.f":
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f        : Module w/ GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f            : Module w/ error and NaN check routines
!  (4 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (5 ) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (6 ) grid_mod.f             : Module w/ horizontal grid information
!  (7 ) time_mod.f             : Module w/ routines for computing time & date
!  (8 ) regrid_1x1_mod.f       : Module w/ routines for regrid 1x1 data
!
!  References:
!  ============================================================================
!  (1 ) Original GFED2 database from Jim Randerson:
!        http://ess1.ess.uci.edu/~jranders/data/GFED2/
!  (2 ) Giglio, L., G.R. van der Werf, J.T. Randerson, G.J. Collatz, and
!        P. Kasibhatla, "Global estimation of burned area using MODIS active
!        fire observations", Atm. Chem. Phys. Discuss, Vol 5, 11091, 2005.
!        http://www.copernicus.org/EGU/acp/acpd/5/11091/acpd-5-11091.pdf
!  (3 ) G.R. van der Werf, J.T. Randerson, L. Giglio, G.J. Collatz, 
!        P.S. Kasibhatla, and A.F. Arellano, Jr., "Interannual variability
!        in global biomass burning emissions from 1997 to 2004", Atm. Chem.
!        Phys. Discuss., submitted, 2005, 
!        http://sheba.geo.vu.nl/~gwerf/pubs/VanderWerfEA2005ACPD.pdf
!
!  NOTES:
!  (1 ) Added private routine GFED2_SCALE_FUTURE (swu, bmy, 5/30/06)
!  (2 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  (3 ) Added BC, OC, SO2, NH3, CO2 species.  Also now can read 2005 GFED2
!        C emissions from disk. (rjp, yxw, bmy, 9/25/06)
!  (4 ) 2006 is now the last year of GFED2 emissions (bmy, 1/2/08)
!  (5 ) Add routines to check if GFED2 has been or must be updated. Added
!        module variable UPDATED, DOY8DAY, and T3HR. Now choice between 3-hr,
!        8-day or monthly data (phs, psk, yc, 12/18/08) 
!  (6 ) Added 9 gaseous biomass burning emissions (tmf, 1/7/09)
!  (7 ) The value of N_SPEC is now determined in INIT_GFED2_BIOMASS and depend
!        on the tracers used in the simulation. (ccc, 4/23/09)
!  (6 ) Adjust call to GFED2_AVAILABLE to reflect that monthly data for
!        2008 is now available on disk (bmy, 7/8/09)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gfed2_biomass_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: GFED2_COMPUTE_BIOMASS
      PUBLIC :: CLEANUP_GFED2_BIOMASS
      PUBLIC :: GFED2_IS_NEW

      !==================================================================
      ! MODULE VARIABLES
      !==================================================================

      ! Scalars
      INTEGER                       :: IDBNOx,  IDBCO,   IDBALK4
      INTEGER                       :: IDBACET, IDBMEK,  IDBALD2
      INTEGER                       :: IDBPRPE, IDBC3H8, IDBCH2O
      INTEGER                       :: IDBC2H6, IDBBC,   IDBOC
      INTEGER                       :: IDBSO2,  IDBNH3,  IDBCO2
      INTEGER                       :: IDBBENZ, IDBTOLU, IDBXYLE
      INTEGER                       :: IDBC2H2, IDBC2H4, IDBGLYX
      INTEGER                       :: IDBMGLY, IDBGLYC, IDBHAC
      INTEGER                       :: DOY8DAY, T3HR
      LOGICAL                       :: UPDATED
      TYPE (XPLEX)                        :: SECONDS

      ! Parameters
      INTEGER,          PARAMETER   :: N_EMFAC = 3
!------------------------------------------------------------------------
      ! Why is this hardwired?  (dkh, 09/20/09) 
      INTEGER,          PARAMETER   :: N_SPEC  = 24
      !INTEGER,          PARAMETER   :: N_SPEC  = 1

      ! Arrays
      INTEGER,          ALLOCATABLE :: VEG_GEN_1x1(:,:)
      TYPE (XPLEX),           ALLOCATABLE :: GFED2_EMFAC(:,:)
      REAL*8, ALLOCATABLE :: D_GFED2_EMFAC(:,:)
      TYPE (XPLEX),           ALLOCATABLE :: GFED2_SPEC_MOLWT(:)
      CHARACTER(LEN=4), ALLOCATABLE :: GFED2_SPEC_NAME(:)
      CHARACTER(LEN=6), ALLOCATABLE :: GFED2_SPEC_UNIT(:)
      TYPE (XPLEX),           ALLOCATABLE :: GFED2_BIOMASS(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GFED2_IS_NEW( ) RESULT( IS_UPDATED )
!
!******************************************************************************
!  Function GFED2_IS_NEW returns TRUE if GFED2 emissions have been updated.
!  (phs, 12/18/08)
!
!  NOTES:
!     (1 ) Used in carbon_mod.f and sulfate_mod.f
!******************************************************************************
!
      ! Function value
      LOGICAL    :: IS_UPDATED

      IS_UPDATED = UPDATED      

      ! Return to calling program
      END FUNCTION GFED2_IS_NEW

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_GFED2( DOY, HH )
!
!******************************************************************************
!     Subroutine GFED2_UPDATE checks if we entered a new GFED period
!     since last emission timestep (ie, last call). The result depends
!     on the emissions time step, and the GFED time period used, as well
!     as MMDDHH at beginning of the GEOS-Chem run (phs, 12/18/08)
!
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DOY  (INTEGER) : Day of the Year (0-366) 
!  (2 ) HH   (INTEGER) : Current hour of the day (0-23)
!
!  NOTES:
!  (1 ) the routine computes the DOY (resp. HOUR) at start of the 8-day (resp.
!       3-hour) period we are in, if the 8-day (resp. 3-hr or synoptic) GFED2
!       option is on. Result is compared to previous value to indicate if new
!       data should be read.
!******************************************************************************
!
      USE LOGICAL_MOD, ONLY : LGFED2BB, L8DAYBB, L3HRBB, LSYNOPBB
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH

      ! Arguments
      INTEGER, INTENT(IN) :: DOY, HH

      ! Local
      INTEGER             :: NEW_T3HR, NEW_DOY8DAY
      
      ! Reset to default
      UPDATED = .FALSE.

      ! Check if we enter a new 3hr GFED period (we assume that
      ! emissions time step is less than a day)
      IF ( L3HRBB .OR. LSYNOPBB ) THEN

         NEW_T3HR = INT( HH / 3 ) * 3

         IF ( NEW_T3HR .NE. T3HR ) THEN
            UPDATED = .TRUE.
            T3HR    = NEW_T3HR
         ENDIF         

      ! or a new 8-day GFED period
      ELSE IF ( L8DAYBB ) THEN

         NEW_DOY8DAY = DOY - MOD( DOY - 1, 8 )

         IF ( NEW_DOY8DAY .NE. DOY8DAY ) THEN
            UPDATED = .TRUE.
            DOY8DAY = NEW_DOY8DAY
         ENDIF

      ! or a new month (we assume that we always do emissions on
      ! 1st day 00 GMT of each month - except for the month the
      ! run starts, for which it is not required)
      ELSE IF ( LGFED2BB ) THEN 

         IF ( ITS_A_NEW_MONTH() ) UPDATED = .TRUE.
      
      ENDIF
      
      
      END SUBROUTINE CHECK_GFED2
      
!------------------------------------------------------------------------------

      SUBROUTINE GFED2_AVAILABLE( YYYY, YMIN, YMAX, MM, MMIN, MMAX )
!     
!******************************************************************************
!     Function GFED2_AVAILABLE checks if data are available for input YYYY/MM
!     date, and constrains the later if needed (phs, 1/5/08)
!     
!     NOTES:
!     (1 ) 
!******************************************************************************
!     
      ! Arguments 
      INTEGER, INTENT(INOUT)           :: YYYY
      INTEGER, INTENT(IN)              :: YMIN, YMAX
      INTEGER, INTENT(INOUT), OPTIONAL :: MM
      INTEGER, INTENT(IN),    OPTIONAL :: MMIN, MMAX


      ! Check year
      IF ( YYYY > YMAX .OR. YYYY < YMIN ) THEN
         
         YYYY = MAX( YMIN, MIN( YYYY, YMAX) )
         
         WRITE( 6, 100 ) YMAX, YMIN, YYYY
 100     FORMAT( 'YEAR > ', i4, ' or YEAR < ', i4, 
     $        '. Using GFED2 biomass for ', i4)
      ENDIF
      

      ! Check month
      IF ( PRESENT( MM ) ) THEN 
         IF ( MM > MMAX .OR. MM < MMIN ) THEN

            MM = MAX( MMIN, MIN( MM, MMAX) )
            
            WRITE( 6, 200 ) MMIN, MMAX, MM
 200        FORMAT( ' ** WARNING ** : MONTH is not within ', i2,'-',
     $              i2, '. Using GFED2 biomass for month #', i2)
         ENDIF
      ENDIF

      ! Return to calling program
      END subroutine GFED2_AVAILABLE

!------------------------------------------------------------------------------

      SUBROUTINE GFED2_COMPUTE_BIOMASS( THIS_YYYY, THIS_MM, BIOM_OUT )
!
!******************************************************************************
!  Subroutine GFED2_COMPUTE_BIOMASS computes the monthly GFED2 biomass burning
!  emissions for a given year and month. (psk, bmy, 4/20/06, 1/2/08)
!
!  This routine has to be called on EVERY emissions-timestep if you use one
!  of the GFED2 options.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THIS_YYYY (INTEGER) : Current year 
!  (2 ) THIS_MM   (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Now references LFUTURE from "logical_mod.f".  Now call private routine
!        GFED2_SCALE_FUTURE to compute future biomass emissions, if necessary. 
!        (swu, bmy, 5/30/06)
!  (2 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  (3 ) 2005 is now the last year of available data (bmy, 10/16/06)
!  (4 ) 2006 is now the last year of available data (bmy, 1/2/08)
!  (5 ) Biomass emissions array BIOM_OUT is now INOUT; automatically update
!        BIOMASS array if needed, account for different GFED2 products:
!        monthly, 8-day, 3hr, synoptic (phs, yc, psk, 12/18/08)
!  (6 ) Adjust call to GFED2_AVAILABLE to reflect that monthly data for
!        2008 is now available on disk (bmy, 7/8/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2,    GET_TAU0
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE JULDAY_MOD,     ONLY : JULDAY, CALDATE
      USE LOGICAL_MOD,    ONLY : LFUTURE
      USE LOGICAL_MOD,    ONLY : L8DAYBB, L3HRBB, LSYNOPBB, LGFED2BB
      USE TIME_MOD,       ONLY : EXPAND_DATE,   TIMESTAMP_STRING
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE TIME_MOD,       ONLY : GET_DAY, GET_HOUR, GET_DAY_OF_YEAR
      USE TIME_MOD,       ONLY : ITS_A_LEAPYEAR
      ! adj_group
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE TRACER_MOD,     ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM


#     include "CMN_SIZE"       ! Size parameters

      ! Arguments 
      INTEGER,             INTENT(IN)    :: THIS_YYYY
      INTEGER,             INTENT(IN)    :: THIS_MM
      TYPE (XPLEX),        INTENT(INOUT) :: BIOM_OUT(IIPAR,JJPAR,N_SPEC)

      ! Local variables
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: I,    J,  N,   N_VEG 
      INTEGER                 :: YYYY, MM, MM1, YYYY1
      INTEGER                 :: YYYYMMDD, HHMMSS
      TYPE (XPLEX)                  :: DM_GEN_1x1(I1x1,J1x1-1)
      TYPE (XPLEX)                  :: BIOM_GEN_1x1(I1x1,J1x1-1,N_SPEC)
      TYPE (XPLEX)                  :: BIOM_GEOS_1x1(I1x1,J1x1,N_SPEC)
      TYPE (XPLEX)                  :: TAU0, TAU1, JD8DAY
      TYPE (XPLEX)                  :: TMP
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=16 )      :: TIME_STR
      INTEGER                 :: DD, HH, DOY

      
      !=================================================================
      ! GFED2_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_GFED2_BIOMASS
         FIRST = .FALSE.

         ! adj_group: kludgy check value of N_SPEC (dkh, 09/20/09).
         ! Forward code should adjust N_SPEC automatically, see note (7), 
         ! but does not yet.  N_SPEC is hardwired at the top of this 
         ! module. 
         ! Since applying TPCORE patch, maybe don't need this. 
!         IF ( ITS_A_TAGCO_SIM() .and. N_SPEC /= 1 ) THEN 
!            CALL ERROR_STOP('N_SPEC needs to be 1 for adj tag co',
!     &                      'GFED2_COMPUTE_BIOMASS')
!         ENDIF 
!         IF ( ITS_A_FULLCHEM_SIM() .and. N_SPEC /= 24 ) THEN 
!            CALL ERROR_STOP('N_SPEC needs to be 24 for adj fullchem',
!     &                      'GFED2_COMPUTE_BIOMASS')
!         ENDIF 

      ENDIF

      ! Save in local variables
      YYYY = THIS_YYYY
      MM   = THIS_MM
      DD   = GET_DAY()
      HH   = GET_HOUR()
      DOY  = GET_DAY_OF_YEAR()
      
      ! Check if we need to update GFED2 (phs, 18/12/08)
      CALL CHECK_GFED2( DOY, HH )
      
      IF ( UPDATED ) THEN
         GFED2_BIOMASS  = 0D0
      ELSE
         BIOM_OUT = GFED2_BIOMASS
         RETURN
      ENDIF
      
      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 
     &  'G F E D 2   B I O M A S S   B U R N I N G   E M I S S I O N S'

      
      !=================================================================
      ! Check GFED2 availability & get YYYYMMDD of data to read.
      !=================================================================
         
      ! Availability of 3-HR data
      !-------------------------------
      IF ( L3HRBB .OR. LSYNOPBB ) THEN

         CALL GFED2_AVAILABLE( YYYY, 2004, 2004, MM, 6, 9 )

         ! Kludge until we have a full year of data: make sure that
         ! DD .ne. 31 if data for MM are not available
         IF ( MM /= THIS_MM ) DD = MIN( 30, DD )
         
         YYYYMMDD = YYYY * 10000 + MM * 100 + DD
         HHMMSS   = T3HR * 10000

         TIME_STR = TIMESTAMP_STRING( YYYYMMDD, HHMMSS )
         
         WRITE( 6, 210 ) TIME_STR
 210     FORMAT( 'for 3-hr period starting: ', a16, / )

         
      ! Availability of 8-DAY data
      !-------------------------------
      ELSE IF ( L8DAYBB ) THEN
         
         CALL GFED2_AVAILABLE( YYYY, 2001, 2007 )
         
         ! Get Julian day at start of 8-day period & its YYYYMMDD
         JD8DAY = JULDAY( YYYY, 1, 0d0) + dble( DOY8DAY )
         CALL CALDATE( JD8DAY, YYYYMMDD, HHMMSS )

         TIME_STR = TIMESTAMP_STRING( YYYYMMDD, 0 )
         
         WRITE( 6, 310 ) TIME_STR
 310     FORMAT( 'for 8-day period starting: ', a16, / )


      ! Availability of MONTHLY data
      !-------------------------------
      ELSE IF ( LGFED2BB ) THEN
         
         !-----------------------------------------------------------
         ! Prior to 7/8/09:
         ! GFED2 2008 monthly data is now available (bmy, 7/8/09)
         !CALL GFED2_AVAILABLE( YYYY, 1997, 2007 )
         !-----------------------------------------------------------
         CALL GFED2_AVAILABLE( YYYY, 1997, 2008 )

         WRITE( 6, 410 ) YYYY, MM
 410     FORMAT( 'for year and month: ', i4, '/', i2.2, / )

         ! Create YYYYMMDD integer value
         YYYYMMDD = YYYY*10000 + MM*100 + 01

      ENDIF

     
      !=================================================================
      ! Filename, TAU0 and number of seconds
      !=================================================================
      
      ! for 3-HR data
      !-------------------------------
      IF ( L3HRBB .OR. LSYNOPBB ) THEN
         
         TAU0 = GET_TAU0( MM, DD, YYYY, T3HR )
        
         IF ( LSYNOPBB ) THEN
            FILENAME = TRIM( DATA_DIR_1x1 )                  //
     &                 'GFED2_3hr_200901/YYYY/'              //
     $                 'GFED2.synoptic.C_YYYYMM.generic.1x1'
         ELSE 
            FILENAME = TRIM( DATA_DIR_1x1 )                  //
     &                 'GFED2_3hr_200901/YYYY/'              //
     $                 'GFED2.3hr.C_YYYYMM.generic.1x1'
         ENDIF

         SECONDS = 3 * 3600d0
        
      
      ! for 8-day data
      !-------------------------------
      ELSE IF (L8DAYBB) THEN

         ! get TAU0 from two JDs: start of 8-day period and 1985/1/1 
         TMP      = ( JD8DAY - 2446066.5d0 ) * 24e0
         TAU0     = ( TMP )
         
         FILENAME = TRIM( DATA_DIR_1x1 )                 //
     &              'GFED2_8day_200712/YYYY/'            //
     $              'GFED2_8day_C_YYYYMMDD.generic.1x1'

         SECONDS  = 8 * 24 * 3600d0

         ! check for last period of the year
         IF ( DOY > 360) THEN
            IF ( ITS_A_LEAPYEAR( YYYY ) ) THEN
               SECONDS = 6*24* 3600d0
            ELSE
               SECONDS = 5*24* 3600d0
            ENDIF
         ENDIF 

         
      ! for monthly data
      !-------------------------------
      ELSE IF ( LGFED2BB ) THEN 
      
         ! TAU value at start of YYYY/MM
         TAU0     = GET_TAU0( MM, 1, YYYY )

         ! Get YYYY/MM value for next month
         MM1      = MM + 1
         YYYY1    = YYYY

         ! Increment year if necessary
         IF ( MM1 == 13 ) THEN
            MM1   = 1
            YYYY1 = YYYY + 1
         ENDIF

         ! TAU value at start of next month
         TAU1     = GET_TAU0( MM1, 1, YYYY1 )

         ! Number of seconds in this month 
         ! (NOTE: its value will be saved until the next month)
         SECONDS  = ( TAU1 - TAU0 ) * 3600d0

         ! File name with GFED2 C emissions
         FILENAME = TRIM( DATA_DIR_1x1 ) //
     &              'GFED2_200601/YYYY/GFED2_C_YYYYMM.generic.1x1'

      ENDIF
      
      !=================================================================
      ! Read GFED2 C emissions [g/m2/month, g/m2/8-day, or g/m2/3-hr]
      !=================================================================
      
      ! Replace YYYY/MM in the file name
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read GFED2 C emissions [g C/m2/month]
      CALL READ_BPCH2( FILENAME, 'GFED2-BB',   99, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         DM_GEN_1x1,  QUIET=.TRUE. ) 

      !=================================================================
      ! Convert C [g/m2/month] to dry matter burned [kg/cm2/month]
      !
      ! Unit Conversions:
      ! (1) C    to DM    --> Divide by 0.45  
      ! (2) g    to kg    --> Divide by 1000  
      ! (3) 1/m2 to 1/cm2 --> Divide by 10000 
      !=================================================================

      ! Loop over GENERIC 1x1 GRID
      DO J = 1, J1x1-1
      DO I = 1, I1x1

         ! Set negatives to zero
         DM_GEN_1x1(I,J) = MAX( DM_GEN_1x1(I,J), 0d0 )

         ! Convert [g C/m2/month] to [kg DM/cm2/month]
         DM_GEN_1x1(I,J) = DM_GEN_1x1(I,J) / ( 0.45d0 * 1d3 * 1d4 )

      ENDDO
      ENDDO

      !=================================================================
      ! Calculate biomass species emissions on 1x1 emissions grid
      !
      ! Emission factors convert from [kg/cm2/month] to either
      ! [molec/cm2/month] or [atoms C/cm2/month]
      !
      ! Units:
      !  [  molec/cm2/month] : NOx,  CO,   CH2O, SO2,  NH3,  CO2
      !  [atoms C/cm2/month] : ALK4, ACET, MEK,  ALD2, PRPE, C3H8, 
      !                        C2H6, BC,   OC
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N_VEG )
         DO J = 1, J1x1-1
         DO I = 1, I1x1
 
            ! Vegetation type index
            N_VEG = VEG_GEN_1x1(I,J)
            
            ! Multiply DM * EMISSION FACTOR to get biomass emissions
            ! for each species on the GENERIC 1x1 GRID 
            SELECT CASE( N_VEG )

               ! Ocean 
               CASE( 0 ) 
                  BIOM_GEN_1x1(I,J,N) = 0d0

               ! Land
               CASE( 1:3 )
                  BIOM_GEN_1x1(I,J,N) = DM_GEN_1x1(I,J) * 
     &                                  GFED2_EMFAC(N,N_VEG)

               ! Otherwise
               CASE DEFAULT
                  ! Nothing

            END SELECT
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Regrid each species from GENERIC 1x1 GRID to GEOS-Chem 1x1 GRID
         CALL DO_REGRID_G2G_1x1( 'molec/cm2',
     &                            BIOM_GEN_1x1(:,:,N), 
     &                            BIOM_GEOS_1x1(:,:,N) )
      ENDDO

      ! Regrid from GEOS 1x1 grid to current grid.  (The unit 'molec/cm2' 
      ! is just used to denote that the quantity is per unit area.)
      CALL DO_REGRID_1x1( N_SPEC,       'molec/cm2', 
     &                    BIOM_GEOS_1x1, GFED2_BIOMASS ) 

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL GFED2_SCALE_FUTURE( GFED2_BIOMASS )
      ENDIF

      ! Print totals in Tg/month
      CALL GFED2_TOTAL_Tg( THIS_YYYY, THIS_MM )

      ! Convert from [molec/cm2/month], [molec/cm2/8day] or
      ! [molec/cm2/3hr] to [molec/cm2/s]
      GFED2_BIOMASS = GFED2_BIOMASS / SECONDS

      ! set output
      BIOM_OUT = GFED2_BIOMASS
      
      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE GFED2_COMPUTE_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE GFED2_SCALE_FUTURE( BB )
!
!******************************************************************************
!  Subroutine GFED2_SCALE_FUTURE applies the IPCC future emissions scale 
!  factors to the GFED2 biomass burning emisisons in order to compute the 
!  future emissions of biomass burning for NOx, CO, and VOC's.  
!  (swu, bmy, 5/30/06, 9/25/06)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) BB (TYPE (XPLEX)) : Array w/ biomass burning emisisons [molec/cm2]
!
!  NOTES:
!  (1 ) Now scale to IPCC future scenario for BC, OC, SO2, NH3 (bmy, 9/25/03)
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_BCbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_CObb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NH3bb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NOxbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_OCbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_SO2bb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_VOCbb
      USE TRACER_MOD,             ONLY : ITS_A_CO2_SIM       

#     include "CMN_SIZE"               ! Size parameters

      ! Arguments
      TYPE (XPLEX),           INTENT(INOUT) :: BB(IIPAR,JJPAR,N_SPEC)

      ! Local variables
      LOGICAL                         :: ITS_CO2
      INTEGER                         :: I, J, N
      
      !=================================================================
      ! GFED2_SCALE_FUTURE begins here!
      !=================================================================

      ! Test if it's a CO2 simulation outside of the loop
      ITS_CO2 = ITS_A_CO2_SIM()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )

      ! Loop over species and grid boxes
      DO N = 1, N_SPEC
      DO J = 1, JJPAR
      DO I = 1, IIPAR 

         ! Scale each species to IPCC future scenario
         IF ( N == IDBNOx ) THEN

            ! Future biomass NOx [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NOxbb( I, J )

         ELSE IF ( N == IDBCO ) THEN

            ! Future biomass CO [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_CObb( I, J )

         ELSE IF ( N == IDBSO2 ) THEN

            ! Future biomass SO2 [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_SO2bb( I, J )

         ELSE IF ( N == IDBNH3 ) THEN

            ! Future biomass NH3 [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NH3bb( I, J )

         ELSE IF ( N == IDBBC ) THEN

            ! Future biomass BC [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_BCbb( I, J )

         ELSE IF ( N == IDBOC ) THEN

            ! Future biomass OC [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_OCbb( I, J )

         ELSE IF ( ITS_CO2 ) THEN

            ! Nothing

         ELSE

            ! Future biomass Hydrocarbons [atoms C/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_VOCbb( I, J )

         ENDIF
         
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GFED2_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE GFED2_TOTAL_Tg( YYYY, MM )
!
!******************************************************************************
!  Subroutine TOTAL_BIOMASS_TG prints the amount of biomass burning emissions 
!  that are emitted each month/8-day/3-hr in Tg or Tg C. (bmy, 3/20/01,
!  12/23/08)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYY    (INTEGER) : Current year
!  (2 ) MM      (INTEGER) : Currrent month
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,   ONLY : GET_AREA_CM2

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYY, MM

      ! Local variables
      INTEGER             :: I,    J,     N
      TYPE (XPLEX)              :: CONV, MOLWT, TOTAL
      CHARACTER(LEN=4)    :: NAME
      CHARACTER(LEN=6)    :: UNIT

      !=================================================================
      ! GFED2_TOTAL_Tg begins here!
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

         ! Initialize
         NAME  = GFED2_SPEC_NAME(N)
         MOLWT = GFED2_SPEC_MOLWT(N)
         UNIT  = GFED2_SPEC_UNIT(N)
         TOTAL = 0d0

         ! Loop over latitudes
         DO J = 1, JJPAR
         
            ! Convert to [Tg/gfed-period] (or [Tg C/gfed-period] for HC's)
            CONV = GET_AREA_CM2( J ) * ( MOLWT / 6.023d23 ) * 1d-9

            ! Loop over longitudes
            DO I = 1, IIPAR
               TOTAL = TOTAL + ( GFED2_BIOMASS(I,J,N) * CONV )
            ENDDO
         ENDDO
     
         ! Write totals
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', 2f9.4, 1x, a6 )
      ENDDO

      ! Return to calling program
      END SUBROUTINE GFED2_TOTAL_Tg

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GFED2_BIOMASS
!
!******************************************************************************
!  Subroutine INIT_GFED2_BIOMASS allocates all module arrays.  It also reads
!  the emission factors and vegetation map files at the start of a GEOS-Chem
!  simulation. (psk, bmy, 4/20/06, 9/25/06)
!
!  NOTES:
!  (1 ) Now initialize for BC, OC, SO2, NH3, CO2 (bmy, 9/25/06)
!  (2 ) Bug fix: IDBSO2, IDBNH3, IDBOC, and IDBCO2 are correctly 
!        set (phs, 3/18/08)
!  (3 ) Add 9 gaseous biomass burning emissions (tmf, 1/7/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE FILE_MOD,      ONLY : IOERROR, IU_FILE
      USE LOGICAL_MOD,   ONLY : LDICARB


#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: AS, IOS, M, N, NDUM
      TYPE (XPLEX)                 :: ARRAY(I1x1,J1x1-1,1)
      CHARACTER(LEN=255)     :: FILENAME
      
      !=================================================================
      ! INIT_GFED2_BIOMASS begins here!
      !=================================================================

      ! Allocate array to hold emissions
      ALLOCATE( GFED2_BIOMASS( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_BIOMASS' )
      GFED2_BIOMASS = 0d0

      ! Allocate array for emission factors
      ALLOCATE( GFED2_EMFAC( N_SPEC, N_EMFAC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_EMFAC' )
      GFED2_EMFAC = 0d0

      ALLOCATE( D_GFED2_EMFAC( N_SPEC, N_EMFAC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'D_GFED2_EMFAC' )
      D_GFED2_EMFAC = 0d0      

      ! Allocate array for species molecular weight
      ALLOCATE( GFED2_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_SPEC_MOLWT' )
      GFED2_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( GFED2_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_SPEC_NAME' )
      GFED2_SPEC_NAME = ''

      ! Allocate array for species molecular weight
      ALLOCATE( GFED2_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_SPEC_UNIT' )
      GFED2_SPEC_UNIT = ''

      ! Allocate array for vegetation map
      ALLOCATE( VEG_GEN_1x1( I1x1, J1x1-1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VEG_GEN_1x1' )

      ! Set default values for module variables
      T3HR    = -1
      DOY8DAY = -1

      !=================================================================
      ! Read emission factors (which convert from kg DM to 
      ! either [molec species] or [atoms C]) from bpch file
      !=================================================================
     
      ! File name
      FILENAME = TRIM( DATA_DIR_1x1) // 
     &           'GFED2_200601/GFED2_emission_factors_73t.txt'

      ! Open emission factor file (ASCII format)
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed2:1' )

      ! Skip header lines
      DO N = 1, 6 
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed2:2' )
      ENDDO

      ! Read emission factors for each species
      DO N = 1, N_SPEC
         READ( IU_FILE, 100, IOSTAT=IOS ) 
     &       NDUM, GFED2_SPEC_NAME(N),(D_GFED2_EMFAC(N,M), M=1,N_EMFAC )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed2:3' )
      ENDDO
      GFED2_EMFAC(:,:) = (D_GFED2_EMFAC(:,:))
      ! FORMAT string
 100  FORMAT( 1x, i2, 1x, a4, 3(3x,es14.6) )

      ! Close file
      CLOSE( IU_FILE )
      
      !=================================================================
      ! Read GFED2 vegetation map from bpch file
      ! 
      ! Values:  3 = boreal forest 
      !          2 = tropical forest; 
      !          1 = savanna / herb / other land
      !          0 = water
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'GFED2_200601/GFED2_vegmap.generic.1x1'

      ! Read GFED2 veg map 
      CALL READ_BPCH2( FILENAME, 'LANDMAP',  1, 
     &                xplex(0d0,0d0),       I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from TYPE (XPLEX) to INTEGER
      VEG_GEN_1x1(:,:) = ARRAY(:,:,1)

      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the GFED2 biomass species
      !=================================================================
      
      ! Initialize 
      IDBNOx  = 0  
      IDBCO   = 0
      IDBALK4 = 0
      IDBACET = 0 
      IDBMEK  = 0 
      IDBALD2 = 0
      IDBPRPE = 0
      IDBC3H8 = 0
      IDBCH2O = 0
      IDBC2H6 = 0
      IDBBC   = 0
      IDBOC   = 0
      IDBSO2  = 0
      IDBNH3  = 0
      IDBCO2  = 0
      IDBGLYX = 0
      IDBMGLY = 0
      IDBBENZ = 0
      IDBTOLU = 0   
      IDBXYLE = 0
      IDBC2H4 = 0
      IDBC2H2 = 0 
      IDBGLYC = 0
      IDBHAC  = 0
 
      ! Save species # in IDBxxxx flags for future reference
      ! and also initialize arrays for mol wts and units
      DO N = 1, N_SPEC
         SELECT CASE ( TRIM( GFED2_SPEC_NAME(N) ) ) 
            CASE( 'NOx'  )
               IDBNOx              = N
               GFED2_SPEC_MOLWT(N) = 14d-3
               GFED2_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'CO'   )
               IDBCO               = N
               GFED2_SPEC_MOLWT(N) = 28d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'ALK4' )
               IDBALK4             = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ACET' )
               IDBACET = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'MEK'  )
               IDBMEK  = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ALD2' )
               IDBALD2 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'PRPE' )
               IDBPRPE = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C3H8' )
               IDBC3H8 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'CH2O' )
               IDBCH2O = N
               GFED2_SPEC_MOLWT(N) = 30d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'C2H6' )
               IDBC2H6 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'SO2'  )
               IDBSO2 = N
               GFED2_SPEC_MOLWT(N) = 64d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'NH3'  )
               IDBNH3 = N
               GFED2_SPEC_MOLWT(N) = 17d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BC'   ) 
               IDBBC = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'OC'   )
               IDBOC = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYX' )
               IDBGLYX = N
               GFED2_SPEC_MOLWT(N) = 58d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'MGLY' )
               IDBMGLY = N
               GFED2_SPEC_MOLWT(N) = 72d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BENZ' )
               IDBBENZ = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'TOLU' )
               IDBTOLU = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'XYLE' )
               IDBXYLE = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H4' )
               IDBC2H4 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H2' )
               IDBC2H2 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYC' )
               IDBGLYC = N
               GFED2_SPEC_MOLWT(N) = 60d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'HAC' )
               IDBHAC  = N
               GFED2_SPEC_MOLWT(N) = 74d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CO2'  )
               IDBCO2 = N
               GFED2_SPEC_MOLWT(N) = 44d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE DEFAULT
               ! Nothing
         END SELECT
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_GFED2_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GFED2_BIOMASS
!
!******************************************************************************
!  Subroutine CLEANUP_GFED2_BIOMASS deallocates all module arrays.
!  (psk, bmy, 4/20/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_GFED2_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( GFED2_EMFAC      ) ) DEALLOCATE( GFED2_EMFAC     )
      IF ( ALLOCATED( D_GFED2_EMFAC    ) ) DEALLOCATE( D_GFED2_EMFAC   )
      IF ( ALLOCATED( GFED2_SPEC_MOLWT ) ) DEALLOCATE( GFED2_SPEC_MOLWT)
      IF ( ALLOCATED( GFED2_SPEC_NAME  ) ) DEALLOCATE( GFED2_SPEC_NAME )
      IF ( ALLOCATED( VEG_GEN_1x1      ) ) DEALLOCATE( VEG_GEN_1x1     )
      IF ( ALLOCATED( GFED2_BIOMASS    ) ) DEALLOCATE( GFED2_BIOMASS   )
      
      ! Return to calling program
      END SUBROUTINE CLEANUP_GFED2_BIOMASS

!------------------------------------------------------------------------------

      ! End of module 
      END MODULE GFED2_BIOMASS_MOD
