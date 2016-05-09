!$Id: chemdr_adj.f,v 1.9 2012/03/01 22:00:26 daven Exp $
      SUBROUTINE CHEMDR_ADJ
!
!******************************************************************************
!  Subroutine CHEMDR_ADJ is the driver subroutine for full chemistry adjoint. 
!  Adapted from original code by lwh, jyl, gmg, djj. (bmy, 11/15/01, 6/3/08)
!  Adjoint developed by dkh, ks, 07/30/09
!
!  Important input variables from "dao_mod.f" and "uvalbedo_mod.f":
!  ============================================================================
!  ALBD        : DAO visible albedo                         [unitless]
!  AVGW        : Mixing ratio of water vapor                [v/v] 
!  BXHEIGHT    : Grid box heights                           [m]
!  OPTD        : DAO grid-box optical depths (for FAST-J)   [unitless]
!  SUNCOS      : Cosine of solar zenith angle               [unitless]
!  SUNCOSB     : Cosine of solar zenith angle 1 hr from now [unitless]
!  UVALBEDO    : TOMS UV albedo 340-380 nm (for FAST-J)     [unitless]
!
!  Important input variables from "comode.h" or "comode_mod.f":
!  ============================================================================
!  NPTS        : Number of points (grid-boxes) to calculate
!  REMIS       : Emission rates                             [molec/cm3/s-1]
!  RAERSOL     : Frequency of gas-aerosol collisions        [s-1]
!  PRESS       : Pressure                                   [Pa]
!  TMPK        : Temperature                                [K]
!  ABSHUM      : Absolute humidity                          [molec/cm3]
!  CSPEC       : Initial species concentrations             [molec/cm3]
!
!  Important output variables in "comode.h" etc.
!  ============================================================================
!  NAMESPEC    : Character array of species names
!  NNSPEC      : # of ACTIVE + INACTIVE (not DEAD) species
!  CSPEC       : Final species concentrations               [molec/cm3]
!
!  Other Important Variables
!  ============================================================================
!  MAXPTS      : Maximum number of points or grid-boxes (in "comsol.h")
!                (NPTS must be <= MAXPTS, for SLOW-J only)
!  MAXDEP      : Maximum number of deposition species (note # of
!                depositing species listed in tracer.dat must be <= MAXDEP)
!  IGAS        : Maximum number of gases, ACTIVE + INACTIVE
!  IO93        : I/O unit for output for "ctm.chem" file
!
!  Input files for SMVGEAR II:
!  ============================================================================
!   mglob.dat  : control switches                       (read in "reader.f")
!  tracer.dat  : list of tracers, emitting species      (read in "reader.f")
!                and depositing species
! globchem.dat : species list, reaction list,           (read in "chemset.f")
!                photolysis reaction list
!
!  Input files for FAST-J photolysis:
!  ============================================================================
!     ratj.d   : Lists photo species, branching ratios  (read in "rd_js.f")
! jv_atms.dat  : Climatology of T and O3                (read in "rd_prof.f")
! jv_spec.dat  : Cross-sections for each species        (read in "RD_TJPL.f")
!
!  Input files for SLOW-J photolysis:
!  ============================================================================
!  jvalue.dat  : Solar flux data, standard T and O3     (read in "jvaluein.f")
!                profiles, aerosol optical depths 
!    8col.dat  : SLOW-J cross-section data              (read in "jvaluein.f")
!  chemga.dat  : Aerosol data
!    o3du.dat  : O3 in Dobson units, cloud data         (read in "jvaluein.f")
!
!  NOTES:
!  (1 ) Cleaned up a lot of stuff.  SUNCOS, OPTD, ALBD, and AVGW are now 
!        referenced from dao_mod.f.  IREF and JREF are obsolete.  Also 
!        updated comments. (bmy, 9/27/01)
!  (2 ) Do not declare LPRT or set LPRT = .FALSE. in "chemdr.f".  LPRT is 
!        included via "CMN" and is defined in "main.f". (bmy, 10/9/01)
!  (3 ) Removed obsolete data from 9/01 (bmy, 10/23/01)
!  (4 ) ERADIUS(JLOOP) is now ERADIUS(JLOOP,1) and TAREA(JLOOP) is now
!        TAREA(JLOOP,1) for sulfate aerosol.  Updated comments. (bmy, 11/15/01)
!  (5 ) Renamed routine PAFTOP to DEBUG_MSG.  Also deleted obsolete code
!        from 11/01.  Enhanced debug output via DEBUG_MSG.  Also reference
!        the UVALBEDO array directly from "uvalbedo_mod.f".  Remove UVALBEDO
!        from the argument list.  Removed obsolete comments. (bmy, 1/15/02)
!  (6 ) Now pass LPAUSE to "initgas.f" via the arg list (bmy, 2/14/02)
!  (7 ) Now call "rdaer.f" instead of RDAEROSOL to read the aerosol and dust 
!        fields from disk.  Also, ignore hygroscopic growth for dust.  Now
!        pass SAVEHO2 and FRACNO2 arrays in the arg list to "ohsave.f"; these 
!        return HO2 conc.'s and NO2 fraction.  Delete NTRACE from call
!        to "ohsave.f", it's obsolete.  Delete reference to DARSFCA from
!        "comode_mod.f", it's obsolete. (rvm, bmy, 2/27/02)
!  (8 ) Removed obsolete code from 2/02. (bmy, 4/15/02)
!  (9 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (10) Now reference IU_CTMCHEM from "file_mod.f".  Assign the value of
!        IU_CTMCHEM (which =93) to IO93 for SMVGEAR routines.  Also open 
!        "ctm.chem" file on the first call as file unit # IO93.  Add
!        references to "CMN_DIAG" and "planeflight_mod.f".  Call routine
!        SETUP_PLANEFLIGHT to initialize the plane track diagnostic
!        after reading the "chem.dat" file.  (bmy, 7/2/02)
!  (11) Now reference AD, T and BXHEIGHT from "dao_mod.f".  Also removed 
!        obsolete commented out code in various sections below.  Now also
!        references "tracerid_mod.f".  Also remove reference to BIOTRCE, since
!        this is now obsolete.  Now make FIRSTCHEM a local SAVED variable
!        instead of an argument.  Now calls MAKE_AVGW, which was formerly
!        called in "main.f". (bmy, 11/15/02)
!  (12) Now reference "AIRVOL" from "dao_mod.f".  Now declare local array
!        SO4_NH4_NIT, which will contain lumped SO4, NH3, NIT aerosol.  Now
!        pass SO4_NH4_NIT to "rdaer.f" via the arg list if sulfate chemistry
!        is turned on.  Now also references CMN_SETUP. (rjp, bmy, 3/23/03)
!  (13) Removed ITAU from the arg list.  Removed reference to IHOUR.  Now use
!        functions GET_MONTH, GET_YEAR from "time_mod.f" (bmy, 3/27/03)
!  (14) Remove KYEAR and TWO_PI, these are now obsolete for SMVGEAR II.  Now 
!        open unit #93 and call READER in the same FIRSTCHEM if-block.  Now
!        Replace call to CHEMSET with call to READCHEM.  JPARSE is now called 
!        from w/in READCHEM.  Replace call to INITGAS w/ call to GASCONC.
!        Removed reference to "file_mod.f".  Remove call to SETPL, we now must
!        call this in "readchem.f" before the call to JSPARSE. 
!        (bdf, ljm, bmy, 5/8/03)
!  (15) Now reference routine GET_GLOBAL_CH4 from "global_ch4_mod.f".  Also
!        added CH4_YEAR as a SAVEd variable. (bnd, bmy, 7/1/03)
!  (16) Remove references to MONTHP, IMIN, ISEC; they are obsolete and not 
!        defined anywhere. (bmy, 7/16/03)
!  (17) Now reference SUNCOSB from "dao_mod.f".  Now pass SUNCOSB to "chem.f". 
!        Also remove LSAMERAD from call to CHEM, since it's obsolete. 
!        (gcc, bmy, 7/30/03)
!  (18) Added BCPO, BCPI, OCPO, OCPI, and SOILDUST arrays.  Now pass SOILDUST
!       to RDUST_ONLINE (in "dust_mod.f").  Now pass PIEC, POEC, PIOC, POOC to
!       "rdaer.f".  Now references "dust_mod.f". (rjp, tdf, bmy, 4/1/04)
!  (19) Added SALA and SALC arrays for passing seasalt to rdaer.f.  Now
!        rearranged the DO loop for computational efficiency. (bmy, 4/20/04)
!  (20) Added OCF parameter to account for the other chemical components that
!        are attached to OC.  Also now handle hydrophilic OC differently for
!        online & offline SOA. (rjp, bmy, 7/15/04)
!  (21) Now reference "logical_mod.f".  Now reference STT and N_TRACERS from
!        "tracer_mod.f".  Now references DO_DIAG_PL from "diag_pl_mod.f".
!        Now references DO_DIAG_OH from "diag_oh_mod.f".  Now references
!        AEROSOL_CONC, RDAER, & SOILDUST from "aerosol_mod.f" (bmy, 7/20/04)
!  (22) Now references ITS_A_NEW_DAY from "time_mod.f".  Now calls routine
!        SETUP_PLANEFLIGHT at the start of each new day. (bmy, 3/24/05)
!  (23) FAST-J is now the default, so we don't need the LFASTJ C-preprocessor 
!        switch any more (bmy, 6/23/05)
!  (24) Now remove LPAUSE from the arg list to "ruralbox.f" and "gasconc.f".
!        (bmy, 8/22/05)
!  (25) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (26) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (27) Remove more obsolete SLOW-J code references.  Also now move function
!        calls from subroutine "chem.f" into "chemdr.f".  Remove obsolete
!        arguments from call to RURALBOX. (bmy, 4/10/06) 
!  (28) Remove reference to "global_ch4_mod.f".  Add error check for LISOPOH
!        when using the online SOA tracers. (dkh, bmy, 6/1/06)
!  (29) Now support variable tropopause (bdf, phs, bmy, 10/3/06)
!  (30) Now get CH4 concentrations for FUTURE_YEAR when using the future
!        emissions scale factors (swu, havala, bmy, 1/28/04)
!  (31) Now call "save_full_trop" at the end to account for "do_diag_pl" 
!        resetting some of CSPEC elements (phs, 6/3/08)
!  (32) Reading the CSPEC_FULL restart file if asked.(dkh, hotp, ccc 2/26/09)
!  (33) Now use GET_DIRECTION
!  (34) LVARTROP treated correctly (dkh, 01/26/11)
!  (35) Add support for strat chem adj LADJ_STRAT and check to make sure that 
!        FD location is in the trop prior to printing debug info for LPRINTFD.
!        (hml, dkh, 02/14/12, adj32_025) 
!******************************************************************************
!
      ! References to F90 modules
      USE AEROSOL_MOD,          ONLY : AEROSOL_CONC, RDAER, SOILDUST
      USE COMODE_MOD,           ONLY : ABSHUM, CSPEC, ERADIUS, TAREA,
     &                                 CSPEC_FOR_KPP, JLOP,    R_KPP
      USE DAO_MOD,              ONLY : AD,       AIRVOL,    ALBD, AVGW   
      USE DAO_MOD,              ONLY : BXHEIGHT, MAKE_AVGW, OPTD, SUNCOS  
      USE DAO_MOD,              ONLY : SUNCOSB,  T
      USE DIAG_OH_MOD,          ONLY : DO_DIAG_OH
      USE DIAG_PL_MOD,          ONLY : DO_DIAG_PL
      USE DUST_MOD,             ONLY : RDUST_ONLINE, RDUST_OFFLINE
      USE ERROR_MOD,            ONLY : DEBUG_MSG,    ERROR_STOP
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_YEAR
      USE LOGICAL_MOD,          ONLY : LCARB,        LDUST,     LEMBED
      USE LOGICAL_MOD,          ONLY : LPRT,         LSSALT,    LSULF  
      USE LOGICAL_MOD,          ONLY : LSOA,         LVARTROP,  LFUTURE
      USE LOGICAL_MOD,          ONLY : LEMIS
      USE PLANEFLIGHT_MOD,      ONLY : SETUP_PLANEFLIGHT
      USE TIME_MOD,             ONLY : GET_MONTH,    GET_YEAR
      USE TIME_MOD,             ONLY : ITS_A_NEW_DAY
      USE TRACER_MOD,           ONLY : STT,          N_TRACERS, XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNOX,       IDTOX,     SETTRACE
      USE TROPOPAUSE_MOD,       ONLY : SAVE_FULL_TROP
      USE UVALBEDO_MOD,         ONLY : UVALBEDO
      ! To use CSPEC_FULL restart (dkh, 02/12/09
      USE RESTART_MOD,          ONLY : READ_CSPEC_FILE 
      USE TIME_MOD,             ONLY : GET_NYMD,     GET_NHMS,  GET_TAU
      USE LOGICAL_MOD,          ONLY : LSVCSPEC
      ! add for adjoint (dkh, 07/31/09) 
      USE ADJ_ARRAYS_MOD,       ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : IFD,          JFD,       LFD
      USE ADJ_ARRAYS_MOD,       ONLY : NFD 
      USE ADJ_ARRAYS_MOD,       ONLY : CHECK_STT_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : ID2C
      USE ADJ_ARRAYS_MOD,       ONLY : IDCSPEC_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : NOBS_CSPEC 
      USE TIME_MOD,             ONLY : GET_DIRECTION
      USE CHECKPT_MOD,          ONLY : CHK_STT_BEFCHEM 
      USE CHEMISTRY_MOD,        ONLY : GCKPP_ADJ_DRIVER
      USE COMODE_MOD,           ONLY : CSPEC_FOR_KPP, CSPEC_ADJ
      !USE COMODE_MOD,           ONLY : NO2_AFTER_CHEM_ADJ
      !USE COMODE_MOD,           ONLY : CSPEC_ADJ_FORCE 
      USE COMODE_MOD,           ONLY : CSPEC_AFTER_CHEM_ADJ
      USE COMODE_MOD,           ONLY : CHK_CSPEC
      USE COMODE_MOD,           ONLY : CSPEC_ORIG
      USE DRYDEP_MOD,           ONLY : NUMDEP
      USE GCKPP_ADJ_Global,     ONLY : NTT
      USE GRID_MOD,             ONLY : GET_XMID, GET_YMID
      USE LOGICAL_MOD,          ONLY : LSCHEM
      USE LOGICAL_ADJ_MOD,      ONLY : LPRINTFD,     LFDTEST
      USE LOGICAL_ADJ_MOD,      ONLY : LCSPEC_OBS
      USE STRAT_CHEM_ADJ_MOD,   ONLY : DO_STRAT_CHEM_ADJ

      USE TIME_MOD,             ONLY : GET_NYMDb,    GET_NHMSb
      USE TIME_MOD,             ONLY : GET_LOCALTIME
      USE TRACERID_MOD,         ONLY : IDTSO2,       IDTSO4

      ! dkh debug 
      USE TRACERID_MOD,         ONLY : IDNO
      USE PARTNER_SAVE_MOD,    ONLY :  SAVE_CSPEC_NETCDF_FILE
      USE PARTNER_SAVE_MOD,    ONLY :  SAVE_NETCDF_FILE
      USE TIME_MOD,              ONLY : EXPAND_DATE


      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"        ! Size parameters
#     include "CMN"             ! IEBD1, IEBD2, etc.
#     include "CMN_O3"          ! EMISRRN, EMISRR
#     include "CMN_NOX"         ! SLBASE
#     include "comode.h"        ! SMVGEAR variables
#     include "CMN_DEP"         ! FRCLND
#     include "CMN_DIAG"        ! ND40
#     include "define_adj.h"    ! OBS operators 

      ! Local variables
      LOGICAL, SAVE            :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE            :: CH4_YEAR  = -1
      INTEGER                  :: I, J, JLOOP, L, NPTS, N, MONTH, YEAR
      ! Now use GET_DIRECTION (dkh, 01/26/10)
      !INTEGER                  :: DIRECTION

      
      ! To use CSPEC_FULL restart (dkh, 02/12/09) 
      LOGICAL                  :: IT_EXISTS 
      
      ! ADJ_GROUP
      INTEGER                  :: NYMD, NHMS
      INTEGER                  :: JLOOPTMP
      INTEGER                  :: IDCSPEC 
      TYPE (XPLEX)                   :: TAU
      TYPE (XPLEX)                :: ADJ_SO4_NH4_NIT(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                   :: ADJ_BCPI(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                   :: ADJ_BCPO(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                   :: ADJ_OCPI(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                   :: ADJ_OCPO(IIPAR,JJPAR,LLPAR)
      ! (dkh, 01/06/12, adj32_006) 
      INTEGER                  :: JJ, NK
      CHARACTER(LEN=255)       :: FILENAME_NC

      !=================================================================
      ! CHEMDR_ADJ begins here!
      !=================================================================

      ! Set some size variables
      NLAT   = JJPAR
      NLONG  = IIPAR
      NVERT  = IVERT 
      NPVERT = NVERT
      NPVERT = NVERT + IPLUME

      ! Get month and year
      MONTH  = GET_MONTH()
      YEAR   = GET_YEAR()

      !=================================================================
      ! Compute AVGW, the mixing ratio of water vapor
      !=================================================================
      CALL MAKE_AVGW

      ! All the FIRSTCHEM stuff will have been done during the forward run, 
      ! Only redo this on the final adjoint step. 
      IF ( GET_NYMD() == GET_NYMDb() .AND.
     &     GET_NHMS() == GET_NHMSb()       ) THEN
         ! dkh debug
         print*, ' FIRSTCHEM = TRUE '  
         FIRSTCHEM = .TRUE.
      ELSE
         FIRSTCHEM = .FALSE.
      ENDIF

      !=================================================================
      ! Open "smv2.log" output file and read chem mechanism switches
      !=================================================================
      IF ( FIRSTCHEM ) THEN
         
         ! Read from data file mglob.dat
         CALL READER( FIRSTCHEM )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after READER' )

         ! Set NCS for urban chemistry only (since that is where we
         ! have defined the GEOS-CHEM mechanism) (bdf, bmy, 4/21/03)
         NCS = NCSURBAN
      ENDIF

      !=================================================================      
      ! Call RURALBOX, which defines tropospheric boxes to be sent to
      ! the SMVGEAR solver, as well as setting up some SMVGEAR arrays.
      !=================================================================      

      ! Redefine NTLOOP since READER defines it initially (bmy, 9/28/04)
      NLOOP  = NLAT  * NLONG
      NTLOOP = NLOOP * NVERT

      CALL RURALBOX( AD,     T,     AVGW,  ALBD,  SUNCOS, 
     &               LEMBED, IEBD1, IEBD2, JEBD1, JEBD2 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after RURALBOX' ) 

      ! Reset NTTLOOP, the # of tropospheric grid boxes
      NTTLOOP = NTLOOP

      !=================================================================
      ! Call SETMODEL which defines number of grid-blocks in calculation,
      ! and copies meteorological parameters into local variables 
      !=================================================================
      CALL SETMODEL

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after SETMODEL' )

      !=================================================================
      ! Do the following only on the first call ...
      !=================================================================
      IF ( FIRSTCHEM ) THEN

         !---------------------------------
         ! Initialize chemistry mechanism
         !---------------------------------

         
         NEMIS(NCSURBAN) = 0 
         NNADDV(NCSURBAN) = 0 
         NNADDA(NCSURBAN) = 0 
         NNADDB(NCSURBAN) = 0 
         NNADDC(NCSURBAN) = 0 
         NNADDD(NCSURBAN) = 0 
         NNADDF(NCSURBAN) = 0 
         NNADDH(NCSURBAN) = 0 
         NNADDG(NCSURBAN) = 0 
         ! dkh debug: try not reading this during adj integration 
         ! to avoid over incrementing NEMIS (dkh, 07/31/09) 
         ! Read "globchem.dat" chemistry mechanism
         CALL READCHEM

         ! Set NCS=NCSURBAN here since we have defined our tropospheric
         ! chemistry mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
         NCS = NCSURBAN

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after READCHEM' )

         !---------------------------------
         ! Check for LISOPOH for SOA
         !---------------------------------
         IF ( LSOA .and. ILISOPOH == 0 ) THEN
            CALL ERROR_STOP( 'LISOPOH needs to be defined for SOA!',
     &                       'chemdr.f' )
         ENDIF

         !---------------------------------
         ! Set global concentration of CH4
         !---------------------------------
         IF ( ICH4 > 0 .and. ( CH4_YEAR /= GET_YEAR() ) ) THEN

            ! If CH4 is a SMVGEAR II species, then call GET_GLOBAL_CH4
            ! to return the globally-varying CH4 conc. as a function of
            ! year and latitude bin.  (ICH4 is defined in READCHEM.)
            ! (bnd, bmy, 7/1/03)
            !
            ! If we are using the future emissions, then get the CH4
            ! concentrations for FUTURE_YEAR.  Otherwise get the CH4
            ! concentrations for the current met field year. 
            ! (swu, havala, bmy, 1/24/08)
            IF ( LFUTURE ) THEN
               CH4_YEAR = GET_FUTURE_YEAR()
            ELSE
               CH4_YEAR = GET_YEAR()
            ENDIF

            ! Get CH4 [ppbv] in 4 latitude bins for each year
            CALL GET_GLOBAL_CH4( CH4_YEAR, .TRUE., C3090S,
     &                           C0030S,   C0030N, C3090N )
         ENDIF

         !-------------------------------
         ! Initialize FAST-J photolysis
         !-------------------------------
         CALL INPHOT( LLTROP, NPHOT ) 
         
         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after INPHOT' )        

         !-------------------------------
         ! Flag certain chemical species
         !-------------------------------
         CALL SETTRACE

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after SETTRACE' )

         !-------------------------------
         ! Flag emission & drydep rxns
         !-------------------------------
         CALL SETEMDEP( N_TRACERS )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after SETEMDEP' )

      ENDIF

      !=================================================================
      ! At the beginning of each new day, call SETUP_PLANEFLIGHT
      ! to see if there are any plane flight points to be processed
      !=================================================================
      IF ( ND40 > 0 .and. ITS_A_NEW_DAY() ) THEN
         CALL SETUP_PLANEFLIGHT
      ENDIF

      !================================================================
      ! Get concentrations of aerosols in [kg/m3] 
      ! for FAST-J and optical depth diagnostics
      !=================================================================
      IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

         ! Skip this section if all these are turned off
         CALL AEROSOL_CONC

      ENDIF

! Now this is done at the end of DO_WETDEP_ADJ
!      ! SO2 and SO4 may have changed during DO_ADJ_WETDEP, so 
!      ! reload their values here. (dkh, 2006)
!      IF ( GET_DIRECTION() < 0 ) THEN
!         STT(:,:,:,IDTSO2) = CHK_STT_BEFCHEM(:,:,:,IDTSO2)
!         STT(:,:,:,IDTSO4) = CHK_STT_BEFCHEM(:,:,:,IDTSO4)
!      ENDIF


      !=================================================================
      ! Call GASCONC which initializes gas concentrations and sets 
      ! miscellaneous parameters.  GASCONC also calls PARTITION, which
      ! splits up family tracers like NOx and Ox into individual
      ! chemical species for SMVGEAR.
      ! NOTE:
      ! (1) The call to GASCONC is modified to use CSPEC_FULL restart 
      !     file (dkh, hotp, ccc,2/26/09)
      !=================================================================
      IT_EXISTS = .FALSE.
      IF ( FIRSTCHEM .AND. LSVCSPEC ) THEN 

         CALL READ_CSPEC_FILE( GET_NYMD(), GET_NHMS(), IT_EXISTS ) 
   
         IF ( .not. IT_EXISTS ) THEN 
             
            ! Use default background values 
            WRITE(6,*) 
     &   '    - CHEMDR: CSPEC restart not found, use background values'
 
            CALL GASCONC( FIRSTCHEM, N_TRACERS, STT, XNUMOL, FRCLND,
     &                    IT_EXISTS )
      
         ELSE 

            ! Use default background values 
            WRITE(6,*) 
     &   '    - CHEMDR: using CSPEC values from restart file'                  

            ! Call GASCONC but don't reset CSPEC values
            CALL GASCONC( .FALSE., N_TRACERS, STT, XNUMOL, FRCLND,
     &                    IT_EXISTS )

         ENDIF 

      ELSE 
         
         ! dkh debug
         !IF ( LPRINTFD ) THEN 
         IF ( LPRINTFD .and. JLOP(IFD,JFD,LFD) > 0 ) THEN 
            print*, ' CSPEC before partition adj = ', 
     &         CSPEC(JLOP(IFD,JFD,LFD),:) 
         ENDIF 

         ! Copy CSPEC top CSPEC_FULL so that CSPEC doesn't get 
         ! overwritten with an old CSPEC_FULL in GASCONC (dkh, 08/04/09) 
          ! LVARTROP support for adj (dkh, 01/26/11)
          ! don't need to do this now becuase we actually checkpt CSPEC_FULL
!         IF ( LVARTROP ) CALL SAVE_FULL_TROP 
         CALL GASCONC( FIRSTCHEM, N_TRACERS, STT, XNUMOL, FRCLND,
     &                 IT_EXISTS )

      ENDIF 
      IT_EXISTS = .FALSE.

      !ADJ_GROUP: Saving CSPEC for KPP calculation
      CSPEC_FOR_KPP(:,:) = CSPEC(:,:)
      CSPEC_ORIG(:,:)    = CSPEC(:,:) 

      !IF ( LPRINTFD ) THEN
      IF ( LPRINTFD .and. JLOP(IFD,JFD,LFD) > 0 ) THEN

         WRITE(6,*) 'CSPEC(FD) after GASCONC = ',
     &               CSPEC(JLOP(IFD,JFD,LFD),:)
         print*, 'STT_ADJ in chemdr_adj', STT_ADJ(IFD,JFD,LFD,:)
      ENDIF

      FILENAME_NC='B_LUMP_CSPEC_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, CSPEC_ADJ)
      FILENAME_NC='B_LUMP_STT_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_NETCDF_FILE(FILENAME_NC, STT_ADJ)

      ! Use dkh's adjoint routine for lumping and partioning 
      CALL LUMP_ADJ( N_TRACERS, XNUMOL, STT_ADJ )

      FILENAME_NC='A_LUMP_CSPEC_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, CSPEC_ADJ)
      FILENAME_NC='A_LUMP_STT_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_NETCDF_FILE(FILENAME_NC, STT_ADJ)
      !IF ( LPRINTFD ) THEN
      IF ( LPRINTFD .and. JLOP(IFD,JFD,LFD) > 0 ) THEN
         print*, 'STT_ADJ after lump_adj', STT_ADJ(IFD,JFD,LFD,:)
      ENDIF

      
      ! Update for new strat chem (hml dkh, 02/14/12, adj32_025) 
      !! SCHEM applies a simplified strat chemistry in order
      !! to prevent stuff from building up in the stratosphere
      !!CALL SCHEM_ADJ
      ! Do stratospheric chemistry adjoint
      IF ( LSCHEM ) CALL DO_STRAT_CHEM_ADJ

      !### Debug
      IF ( LPRT ) 
     &   CALL DEBUG_MSG( '### CHEMDR_ADJ: after STRAT_CHEM_ADJ' )

      ! Reset dry dep adjoints (fp, dkh, 01/06/12, adj32_006) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( N, NK, JJ )
      DO N = 1,NUMDEP
         NK = NTDEP(N)
         IF (NK.NE.0) THEN
            JJ = IRM(NPRODLO+1,NK,NCS)
            IF (JJ.GT.0) THEN
               CSPEC_ADJ(:,JJ) = 0.0D0
            ENDIF
         ENDIF
      ENDDO
!$OMP END PARALLEL DO


      ! Apply forcing from observation (or sensitivyt w.r.t) 
      ! of CSPEC species. (dkh, 10/25/07)
      ! Now use CSPEC_AFTER_CHEM_ADJ (dkh, 02/09/11) 
      IF ( LCSPEC_OBS ) THEN 
         DO N = 1, NOBS_CSPEC

            IDCSPEC = IDCSPEC_ADJ(N) 

            CSPEC_ADJ(:,IDCSPEC) = CSPEC_ADJ(:,IDCSPEC) 
     &                           + CSPEC_AFTER_CHEM_ADJ(:,N)

            CSPEC_AFTER_CHEM_ADJ(:,N) = 0d0

         ENDDO 
      ENDIF 

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after GASCONC' )

      !=================================================================  
      ! Call SETEMIS which sets emission rates REMIS 
      !=================================================================
      CALL SETEMIS( EMISRR, EMISRRN )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after SETEMIS' )

      !=================================================================
      ! Call RDAER -- computes aerosol optical depths
      !=================================================================
      CALL RDAER( MONTH, YEAR )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after RDAER' )

      !=================================================================
      ! If LDUST is turned on, then we have online dust aerosol in
      ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
      ! compute aerosol optical depth for FAST-J, etc.
      !
      ! If LDUST is turned off, then we do not have online dust aerosol
      ! in GEOS-CHEM...so read monthly-mean dust files from disk.
      ! (rjp, tdf, bmy, 4/1/04)
      !=================================================================
      IF ( LDUST ) THEN
         CALL RDUST_ONLINE( SOILDUST )
      ELSE
         CALL RDUST_OFFLINE( MONTH, YEAR )
      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after RDUST' )

      NPTS = NTTLOOP

      ! At present, we are only doing tropospheric chemistry, which 
      ! for the moment we are storing in SMVGEAR II's "urban" slot
      NCS = NCSURBAN

      !=================================================================
      ! Call photolysis routine to compute J-Values
      !=================================================================
      CALL FAST_J( SUNCOS, OPTD, UVALBEDO )              

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after FAST-J' )

      !================================================================
      ! Call chemistry routines
      !================================================================

      ! PHYSPROC calls both CALCRATE, which computes rxn rates 
      ! and SMVGEAR, which is the chemistry solver
      CALL PHYSPROC( SUNCOS, SUNCOSB )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR_ADJ: after PHYSPROC' )

      NHMS        = GET_NHMS()
      NYMD        = GET_NYMD()
      TAU         = GET_TAU()

      NTT = NTTLOOP
      
      !================================================================
      !  Call KPP generated chemical solver. DIRECTION = -1 is adjoint
      !================================================================
      FILENAME_NC='B_KPP_CSPEC_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, CSPEC_ADJ)
      FILENAME_NC='B_KPP_STT_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_NETCDF_FILE(FILENAME_NC, STT_ADJ)
         print*, 'ddd jk STT_ADJ(FD) before gckpp_adj_driver',
     &      STT_ADJ(IFD,JFD,LFD,NFD)
         print*, 'ddd jk CSPEC_ADJ before gckpp_adj_driver',
     &              CSPEC_ADJ(JLOP(IFD,JFD,LFD),85)

      CALL GCKPP_ADJ_DRIVER( GET_DIRECTION() )

         print*, 'ddd jk CSPEC_ADJ after gckpp_adj_driver',
     &              CSPEC_ADJ(JLOP(IFD,JFD,LFD),85)
         print*, 'ddd jk STT_ADJ(FD) after gckpp_adj_driver',
     &      STT_ADJ(IFD,JFD,LFD,NFD)
      FILENAME_NC='A_KPP_CSPEC_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, CSPEC_ADJ)
      FILENAME_NC='A_KPP_STT_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_NETCDF_FILE(FILENAME_NC, STT_ADJ)

      !### Debug
      IF ( LPRT ) 
     &   CALL DEBUG_MSG( '### CHEMDR_ADJ: after GCKPP_ADJ_DRIVER')

      ! Can compare KPP to SMVGEAR side-by-side (dkh, 06/20/05)
      !IF ( LPRINTFD ) THEN
      IF ( LPRINTFD .and. JLOP(IFD,JFD,LFD) > 0 ) THEN

      
         !! Display times (SMVGEAR time will be padded by time to do CALCRATE )
         !WRITE(6,*) ' SMVGEAR TIME : ', TIME2 - TIME1
         !WRITE(6,*) ' ROSENBK TIME : ', TIME3 - TIME2
         
         ! Display comparison for a particular cell 
         JLOOPTMP =  JLOP(IFD,JFD,LFD)
         WRITE(6,*)  ' Spot test in cell:', JLOOPTMP
         WRITE(6,*)  ' LON = ', GET_XMID(IFD)
         WRITE(6,*)  ' LAT = ', GET_YMID(JFD)
         WRITE(6,*)  ' LOCAL TIME = ', GET_LOCALTIME(IFD)
         
         WRITE(6,*)  ' Species        SMVGEAR         ROS          R/S      
     &         ORIG '
         WRITE(6,69) (NAMEGAS(I),CSPEC(JLOOPTMP,I),
     &         CSPEC_FOR_KPP(JLOOPTMP,I),
     &         CSPEC_FOR_KPP(JLOOPTMP,I) /  CSPEC(JLOOPTMP,I),
     &         CSPEC_ORIG(JLOOPTMP,I),
     &         I=1,87)
 69      FORMAT(A10,1X,F20.2,1X,F20.2,1X,1PE10.2,1X,F20.2)

      ENDIF



      ! dkh debug
      WRITE(6,*) '     - CHECK_STT_ADJ after GCKPP_ADJ_DRIVER'
      CALL CHECK_STT_ADJ('AFTER GCKPP_ADJ_DRIVER')

      ! We don't call any adjoint of PHYSPROC.  We just directly call CALCRATE_ADJ
      ! from within GCKPP_ADJ_DRIVER.  That saves us some RAM. 
      

      !=================================================================
      ! Do adjoint of rdaer
      !=================================================================
      print*, ' NEED to update ADJ_RDAER ' 
      print*, ' NEED to update ADJ_RDAER ' 
      print*, ' NEED to update ADJ_RDAER ' 
      print*, ' NEED to update RDAER_ADJ ' 
!      CALL RDAER_ADJ( SO4_NH4_NIT, BCPI, BCPO, OCPI, OCPO,
!     &                ADJ_SO4_NH4_NIT, ADJ_BCPI, ADJ_BCPO,
!     &                ADJ_OCPI, ADJ_OCPO )

      ! For now, don't include these in FDTESTs 
      IF ( LFDTEST ) THEN
         ADJ_SO4_NH4_NIT = 0d0
         ADJ_BCPI        = 0D0
         ADJ_BCPO        = 0D0
         ADJ_OCPI        = 0D0
         ADJ_OCPO        = 0D0
      ENDIF

      !=================================================================
      ! Do adjoint of setemis
      !=================================================================
      CALL SETEMIS_ADJ

      !=================================================================
      ! Do adjoint of partitioning 
      !=================================================================
      ! To do this we need STT = STT_BEFCHEM and
      ! CSPEC = CSPEC_CHK = CSPEC from after chem of step n  - 1. 
      ! CSPEC neads to be reloaded, it 
      ! was overwritten in the above call to PARTITION
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N  )    
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         STT(I,J,L,N) = CHK_STT_BEFCHEM(I,J,L,N)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP, N )    
      DO N     = 1, IGAS
      DO JLOOP = 1, ITLOOP

         CSPEC(JLOOP,N) = CHK_CSPEC(JLOOP,N)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !IF ( LPRINTFD ) THEN
      IF ( LPRINTFD  .and. JLOP(IFD,JFD,LFD) > 0 )  THEN
         print*, 'CSPEC_ADJ before partition_adj',
     &      CSPEC_ADJ(JLOP(IFD,JFD,LFD),:)
      ENDIF
         print*, 'ddd jk CSPEC_ADJ before partition',
     &              CSPEC_ADJ(JLOP(IFD,JFD,LFD),85)
         print*, 'ddd jk STT_ADJ before partition',
     &                   STT_ADJ(IFD,JFD,LFD,NFD)
      
      NHMS        = GET_NHMS()
      NYMD        = GET_NYMD()
      FILENAME_NC='B_PART_CSPEC_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, CSPEC_ADJ)
      FILENAME_NC='B_PART_STT_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_NETCDF_FILE(FILENAME_NC, STT_ADJ)

      ! Use partition adjoint from dkh 
      CALL PARTITION_ADJ( STT_ADJ, STT, N_TRACERS, XNUMOL )
      FILENAME_NC='A_PART_CSPEC_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_CSPEC_NETCDF_FILE(FILENAME_NC, CSPEC_ADJ)
      FILENAME_NC='A_PART_STT_ADJ.YYYYMMDD.hhmm.nc'
      CALL EXPAND_DATE( FILENAME_NC, GET_NYMD(), GET_NHMS() )
!      CALL SAVE_NETCDF_FILE(FILENAME_NC, STT_ADJ)

      ! Now use insead of CSPEC_AFTER_CHEM_ADJ (nb, dkh, 01/06/12, adj32_003)
!#if      defined( SCIA_KNMI_NO2_OBS ) || defined( SCIA_DAL_NO2_OBS )
!      ! Apply forcing from satellite observations
!      CSPEC_ADJ(:,IDNO2)   = CSPEC_ADJ(:,IDNO2) + CSPEC_NO2_ADJ(:)
!      CSPEC_NO2_ADJ(:)     = 0d0
!#endif

         print*, 'ddd jk CSPEC_ADJ afte partition',
     &              CSPEC_ADJ(JLOP(IFD,JFD,LFD),85)
         print*, 'ddd jk STT_ADJ after partition',
     &                   STT_ADJ(IFD,JFD,LFD,NFD)

      !### Debug
      IF ( LPRT )
     &   CALL DEBUG_MSG( '### CHEMDR_ADJ: after PARTITION_ADJ' )

      !IF ( LPRINTFD ) THEN
      IF ( LPRINTFD .and. JLOP(IFD,JFD,LFD) > 0 ) THEN
         print*, 'CSPEC_ADJ after partition_adj',
     &      CSPEC_ADJ(JLOP(IFD,JFD,LFD),:)
         print*, 'STT_ADJ after partition_adj',
     &      STT_ADJ(IFD,JFD,LFD,:)
      ENDIF

      ! dkh debug
      WRITE(6,*) 'CHECK_STT_ADJ after PARTITION_ADJ'
      CALL CHECK_STT_ADJ('after partition_adj')


      ! Adjoint of AEROSOL_CONT
      !=================================================================
      IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

         ! Skip this section if all these are turned off
         print*, ' NEED to updae AEROSOL_CONC_ADJ' 
         print*, ' NEED to updae AEROSOL_CONC_ADJ' 
         print*, ' NEED to updae AEROSOL_CONC_ADJ' 
         print*, ' NEED to updae AEROSOL_CONC_ADJ' 
         !CALL AEROSOL_CONC_ADJ

      ENDIF


      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### Now exiting CHEMDR_ADJ!' )

      ! Return to calling program
      END SUBROUTINE CHEMDR_ADJ




