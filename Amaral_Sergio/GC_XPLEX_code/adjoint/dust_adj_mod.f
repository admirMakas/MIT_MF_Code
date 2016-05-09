! $Id: dust_adj_mod.f,v 1.1 2012/03/01 22:00:26 daven Exp $
      MODULE DUST_ADJ_MOD
!
!******************************************************************************
!  Module DUST_ADJ_MOD contains arrays and routines for performing mineral 
!  dust adjoint simulation. Original code taken from forward model routines
!  in DUST_MOD and modified accordingly. (xxu, 5/20/11)
!  Added to adjoint standard code (xxu, dkh, 01/09/12, adj32_011)
!
!  Module Variables:
!  ============================================================================
!
!  Module Routines:
!  ============================================================================
!  (1 ) EMISSDUST_ADJ      : Driver routine for adjoint dust 
!  (2 ) SRC_DUST_DEAD_ADJ  : Adjoint of DEAD dust emits
!  (3 ) CHEMDUST_ADJ       : Adjoint of DEAD dust emits
!  (4 ) DRY_SETTLING_ADJ   : Adjoint of dust settling 
!  (5 ) DRY_DEPOSITION_ADJ : Adjoint of dust dry deposition
!
!  GEOS-CHEM modules referenced by "dust_mod.f"
!  ============================================================================
!  (1 ) dao_mod.f          : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f         : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) directory_mod.f    : Module containing GEOS-CHEM data & met field dirs
!  (4 ) drydep_mod.f       : Module containing GEOS-CHEM drydep routines
!  (5 ) dust_dead_mod.f    : Module containing Zender's DEAD dust routines
!  (6 ) error_mod.f        : Module containing I/O error and NaN check routines
!  (7 ) file_mod.f         : Contains file unit numbers and error checks
!  (8 ) grid_mod.f         : Module containing horizontal grid information
!  (9 ) logical_mod.f      : Module containing GEOS-CHEM logical switches
!  (10) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!  (11) time_mod.f         : Module containing routines for computing time/date
!  (12) tracer_mod.f       : Module containing GEOS-CHEM tracer array STT etc.
!  (13) tracerid_mod.f     : Module containing pointers to tracers & emissions
!  
!  NOTES:
!  (1 ) See forward model module for complete documentation. (xxu, 5/20/11)
!  (2 ) Implemented the CHEMDUST_ADJ for adjoint change of dry deposition & 
!       settling. (dkh, 1/10/12) 
!  
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "dust_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC                 :: EMISSDUST_ADJ
      PUBLIC                 :: CHEMDUST_ADJ

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE EMISSDUST_ADJ
!
!******************************************************************************
!  Subroutine EMISSDUST_ADJ is the driver routine for the adjoint of 
!  the mineral dust emission. (xxu, 5/20/11)
! 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_MOD,       ONLY : LDEAD
      USE TRACERID_MOD,      ONLY : IDTDST1, IDTDST4
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ

      !================================================================
      ! EMISSDUST_ADJ begins here!
      !================================================================

      ! Check the selected dust scheme
      IF ( LDEAD ) THEN

         ! Adjoint of Zender's DEAD dust source
         CALL SRC_DUST_DEAD_ADJ( STT_ADJ(:,:,:,IDTDST1:IDTDST4) )

      ELSE

         ! Adjoint of Ginoux dust source not yet supported
         CALL ERROR_STOP(' Adjoint of Ginoux dust not yet supported', 
     &                   'dust_adj_mod.f')

      ENDIF

      ! Return to calling program
      END SUBROUTINE EMISSDUST_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE SRC_DUST_DEAD_ADJ( TC )
!
!******************************************************************************
!  Subroutine SRC_DUST_DEAD_ADJ is the adjoint routine for DEAD dust
!  emissions. (xxu, 5/20/11)
!  Based on forward model code. (tdf, bmy, 4/8/04, 1/23/07)
!  
!  NOTES
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : BXHEIGHT,     GWETTOP,   LWI
      USE DAO_MOD,       ONLY : SNOW,         SPHU,      T
      USE DAO_MOD,       ONLY : TS,           UWND,      VWND
      USE DAO_MOD,       ONLY : SNOMAS
      USE DUST_DEAD_MOD, ONLY : GET_TIME_INVARIANT_DATA, GET_ORO
      USE DUST_DEAD_MOD, ONLY : GET_MONTHLY_DATA,        DST_MBL
      USE DIAG_MOD,      ONLY : AD06
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE GRID_MOD,      ONLY : GET_YMID_R
      USE PRESSURE_MOD,  ONLY : GET_PEDGE,       GET_PCENTER
      USE TIME_MOD,      ONLY : GET_TS_EMIS,     GET_MONTH
      USE TIME_MOD,      ONLY : GET_DAY_OF_YEAR, ITS_A_NEW_MONTH
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      ! adj_group: add for emissions scale factors (xxu, 11/02/10)
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EDST1, IDADJ_EDST2
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_EDST3, IDADJ_EDST4
      USE ADJ_ARRAYS_MOD,    ONLY : IS_DUST_EMS_ADJ
      USE DUST_MOD,          ONLY : GET_SCALE_GROUP 
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND06
#     include "CMN_GCTM"      ! g0

      !----------------
      ! Arguments
      !----------------
      TYPE (XPLEX),  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      !-----------------
      ! Local variables
      !-----------------

      ! Scalars
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,      J,      L,       N
      INTEGER                :: M,      IOS,    INC,     LAT_IDX
      INTEGER                :: NDB,    NSTEP,  MM
      TYPE (XPLEX)                 :: W10M,   DEN,    DIAM,    U_TS0
      TYPE (XPLEX)                 :: U_TS,   SRCE_P, Reynol,  YMID_R
      TYPE (XPLEX)                 :: ALPHA,  BETA,   GAMMA,   CW
      TYPE (XPLEX)                 :: DTSRCE, XTAU,   P1,      P2
      TYPE (XPLEX)                 :: DOY
      CHARACTER(LEN=255)     :: FILENAME

      ! Arrays
      INTEGER                :: OROGRAPHY(IIPAR,JJPAR)
      TYPE (XPLEX)                 :: PSLON(IIPAR)         ! surface pressure
      TYPE (XPLEX)                 :: PTHICK(IIPAR)        ! delta P (L=1)
      TYPE (XPLEX)                 :: PMID(IIPAR)          ! mid layer P (L=1)
      TYPE (XPLEX)                 :: TLON(IIPAR)          ! temperature (L=1)
      TYPE (XPLEX)                 :: THLON(IIPAR)         ! pot. temp. (L=1)
      TYPE (XPLEX)                 :: ULON(IIPAR)          ! U-wind (L=1)
      TYPE (XPLEX)                 :: VLON(IIPAR)          ! V-wind (L=1)
      TYPE (XPLEX)                 :: BHT2(IIPAR)          ! half box height (L=1)
      TYPE (XPLEX)                 :: Q_H2O(IIPAR)         ! specific humidity (L=1)
      TYPE (XPLEX)                 :: ORO(IIPAR)           ! "orography" 
      TYPE (XPLEX)                 :: SNW_HGT_LQD(IIPAR)   ! equivalent snow ht.
      TYPE (XPLEX)                 :: DSRC(IIPAR,NDSTBIN)  ! dust mixing ratio incr.

      !----------------
      ! Parameters
      !----------------      
      TYPE (XPLEX), PARAMETER      :: Ch_dust = xplex(9.375d-10,0d0)
      TYPE (XPLEX), PARAMETER      :: G       = xplex(g0%r * 1.D2,0d0)
      TYPE (XPLEX), PARAMETER      :: RHOA    = xplex(1.25D-3,0d0)
      TYPE (XPLEX), PARAMETER      :: CP      = xplex(1004.16d0,0d0)
      TYPE (XPLEX), PARAMETER      :: RGAS=xplex(8314.3d0 / 28.97d0,0d0)
      TYPE (XPLEX), PARAMETER      :: AKAP    = xplex(RGAS%r/CP%r,0d0)
      TYPE (XPLEX), PARAMETER      :: P1000   = xplex(1000d0,0d0)

      ! External functions
      TYPE (XPLEX),  EXTERNAL       :: SFCWINDSQR

      !=================================================================
      ! SRC_DUST_DEAD begins here!
      !=================================================================      

      ! DTSRCE is the emission timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

      ! DOY is the day of year (0-365 or 0-366)
      DOY    = XPLX( GET_DAY_OF_YEAR() )

! fwd cide: 
! we don't need to read in the time invariant data again 
!      !=================================================================
!      ! Read data fields for the DEAD model from disk
!      !=================================================================      
!      IF ( FIRST ) THEN
!
!         ! Echo info
!         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!         WRITE( 6, 100   )
!         WRITE( 6, 110   )
!         WRITE( 6, 120   )
!         WRITE( 6, 130   ) 
!         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!
!         ! FORMAT strings
! 100     FORMAT( 'D E A D   D U S T   M O B I L I Z A T I O N'         )
! 110     FORMAT( 'Routines from DEAD model by Charlie Zender et al'    )
! 120     FORMAT( 'Modified for GEOS-CHEM by D. Fairlie and R. Yantosca')
! 130     FORMAT( 'Last Modification Date: 1/23/07'                     )
!
!         ! Read fields for DEAD that are time-invariant
!         CALL GET_TIME_INVARIANT_DATA
!
!         ! Reset first-time flag
!         FIRST = .FALSE.
!      ENDIF

      ! Read monthly data for DEAD
      IF ( ITS_A_NEW_MONTH() ) THEN
         CALL GET_MONTHLY_DATA
      ENDIF

      ! Determine group (temporal)
      MM = GET_SCALE_GROUP()
      ! Print out scaling info
      WRITE(6,*) ' - READ / RESCALE DUST EMISSION: use SCALE_GROUP ', MM

      !=================================================================
      ! Call dust mobilization scheme
      !=================================================================

      ! Make OROGRAPHY array from GEOS-CHEM LWI
      CALL GET_ORO( OROGRAPHY )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,      P1,   P2,   PTHICK, PMID, TLON        )
!$OMP+PRIVATE( THLON, ULON,   VLON, BHT2, Q_H2O,  ORO,  SNW_HGT_LQD )
!$OMP+PRIVATE( N,     YMID_R, DSRC                                  )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Pressure [hPa] at bottom and top edge of level 1
            P1             = GET_PEDGE(I,J,1)
            P2             = GET_PEDGE(I,J,2)

            ! Pressure thickness of 1st layer [Pa]
            PTHICK(I)      = ( P1 - P2 ) * 100d0

            ! Pressure at midpt of surface layer [Pa]
            PMID(I)        = GET_PCENTER(I,J,1) * 100d0

            ! Temperature [K] at surface layer
            TLON(I)        = T(I,J,1)

            ! Potential temperature [K] at surface layer
            THLON(I)       = TLON(I) * ( P1000 / PMID(I) )**AKAP

            ! U and V winds at surface [m/s]
            ULON(I)        = UWND(I,J,1)
            VLON(I)        = VWND(I,J,1)

            ! Half box height at surface [m]
            BHT2(I)        = BXHEIGHT(I,J,1) / 2.d0

            ! Specific humidity at surface [kg H2O/kg air]
            Q_H2O(I)       = SPHU(I,J,1) / 1000.d0

            ! Orography at surface
            ! Ocean is 0; land is 1; ice is 2
            ORO(I)         = OROGRAPHY(I,J)

            ! Snow height [m H2O]
#if   defined( GEOS_5  )
            SNW_HGT_LQD(I) = SNOMAS(I,J)
#else
            SNW_HGT_LQD(I) = SNOW(I,J) / 1000d0
#endif
            ! Dust tracer and increments
            DO N = 1, NDSTBIN
               DSRC(I,N)   = 0.0d0
            ENDDO
         ENDDO 

         !==============================================================
         ! Call dust mobilization driver (DST_MBL) for latitude J
         !==============================================================

         ! Latitude in RADIANS
         YMID_R = GET_YMID_R(J)

         ! Call DEAD dust mobilization
         CALL DST_MBL( DOY,    BHT2,  J,     YMID_R, ORO,
     &                 PTHICK, PMID,  Q_H2O, DSRC,   SNW_HGT_LQD,
     &                 DTSRCE, TLON,  THLON, VLON,   ULON,
     &                 FIRST,  J )

         ! Update
         DO N = 1, NDSTBIN
         DO I = 1, IIPAR

            ! fwd code:
            !! Include dust adjoint scale factor (xxu, 11/02/10)
            !IF ( LADJ_EMS .and. IS_DUST_EMS_ADJ) THEN
            !
            !   IF(N==1) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST1)
            !   IF(N==2) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST2)
            !   IF(N==3) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST3)
            !   IF(N==4) DSRC(I,N) = DSRC(I,N)*EMS_SF(I,J,MM,IDADJ_EDST4)
            !
            !ENDIF
            !
            !! Add dust emissions into tracer array [kg]
            !TC(I,J,1,N) = TC(I,J,1,N) + DSRC(I,N)

            ! adj code:
            ! Include dust adjoint scale factor (xxu, 5/20/11) 
            IF ( LADJ_EMS .and. IS_DUST_EMS_ADJ) THEN

               IF (N==1) EMS_SF_ADJ(I,J,MM,IDADJ_EDST1)
     &                    = EMS_SF_ADJ(I,J,MM,IDADJ_EDST1) 
     &                    + DSRC(I,N) * TC(I,J,1,N)
               IF (N==2) EMS_SF_ADJ(I,J,MM,IDADJ_EDST2)
     &                    = EMS_SF_ADJ(I,J,MM,IDADJ_EDST2)
     &                    + DSRC(I,N) * TC(I,J,1,N)
               IF (N==3) EMS_SF_ADJ(I,J,MM,IDADJ_EDST3)
     &                    = EMS_SF_ADJ(I,J,MM,IDADJ_EDST3)
     &                    + DSRC(I,N) * TC(I,J,1,N)
               IF (N==4) EMS_SF_ADJ(I,J,MM,IDADJ_EDST4)
     &                    = EMS_SF_ADJ(I,J,MM,IDADJ_EDST4)
     &                    + DSRC(I,N) * TC(I,J,1,N)

            ENDIF         

            ! ND19 diagnostics [kg]
            IF ( ND06 > 0 ) THEN
               AD06(I,J,N) = AD06(I,J,N) + DSRC(I,N)
            ENDIF
         ENDDO
         ENDDO

      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRC_DUST_DEAD_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE CHEMDUST_ADJ
!        
!******************************************************************************
!  Subroutine CHEMDUST_ADJ is the adjoint of CHEMDUST.  Based on the forward
!  model routine (tdf, bmy, 3/30/04, 5/23/06).   (dkh, 01/10/12, adj32_011) 
!           
!  NOTES:   
!           
!******************************************************************************
!           
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE LOGICAL_MOD,    ONLY : LDRYD,   LDUST
      USE DRYDEP_MOD,     ONLY : DEPNAME, NUMDEP
      USE TRACERID_MOD,   ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4

#     include "CMN_SIZE"       ! Size parameters
            
      ! Local variables
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: N
     
      !=================================================================
      ! CHEMDUST_ADJ begins here!
      !=================================================================
     
      !=================================================================
      ! Do dust settling & deposition adjoints
      !=================================================================
               
      ! Dust deposition.  
      IF ( LDRYD ) THEN   
         ! fwd:
         !CALL DRY_DEPOSITION( STT(:,:,:,IDTDST1:IDTDST4) )
         ! adj:
         CALL DRY_DEPOSITION_ADJ( STT_ADJ(:,:,:,IDTDST1:IDTDST4) )
      ENDIF 
               
      ! Dust settling
      ! fwd:
      !CALL DRY_SETTLING(   STT(:,:,:,IDTDST1:IDTDST4) )
      ! adj:
      CALL DRY_SETTLING_ADJ( STT_ADJ(:,:,:,IDTDST1:IDTDST4) )

      ! Return to calling program
      END SUBROUTINE CHEMDUST_ADJ
      
!------------------------------------------------------------------------------
      SUBROUTINE DRY_SETTLING_ADJ( TC )
!
!******************************************************************************
!  Subroutine DRY_SETTLING_ADJ computes the adjoint dry settling of dust tracers.
!
!  Based on the forward model routine (tdf, bmy, 3/30/04, 10/25/05), modified
!  here to calculate the adjoint (dkh, 01/10/12, adj32_011), and with taking
!  out superfluous diagnostic code. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Dust tracer adjoint array 
!
!  NOTES
!
!******************************************************************************
! 
      USE DAO_MOD,      ONLY : T, BXHEIGHT
      USE DUST_MOD,     ONLY : DUSTREFF, DUSTDEN
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTDST1

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! g0

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! Local variables
      INTEGER               :: I, J, L, N
      TYPE (XPLEX)                :: DT_SETTL, DELZ,  DELZ1
      TYPE (XPLEX)                :: REFF,     DEN,   CONST
      TYPE (XPLEX)                :: NUM,      LAMDA
      TYPE (XPLEX)                :: AREA_CM2, TC0(LLPAR)

      ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa      
      TYPE (XPLEX)                :: P

      ! Diameter of aerosol [um]
      TYPE (XPLEX)                :: Dp

      ! Pressure * DP
      TYPE (XPLEX)                :: PDp

      ! Temperature (K)    
      TYPE (XPLEX)                :: TEMP

      ! Slip correction factor
      TYPE (XPLEX)                :: Slip

      ! Viscosity of air (Pa s)
      TYPE (XPLEX)                :: Visc

      ! Settling velocity of particle (m/s)
      TYPE (XPLEX)                :: VTS(LLPAR)

      ! Parameters
      TYPE (XPLEX),  PARAMETER    :: C1 = xplex( 0.7674D0,0d0)
      TYPE (XPLEX),  PARAMETER    :: C2 = xplex( 3.079d0,0d0)
      TYPE (XPLEX),  PARAMETER    :: C3 = xplex( 2.573D-11,0d0)
      TYPE (XPLEX),  PARAMETER    :: C4 = xplex(-1.424d0,0d0)

      !=================================================================
      ! DRY_SETTLING_ADJ begins here!
      !=================================================================

      ! Dust settling timestep [s]
      DT_SETTL = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,        L,    N,     DEN,  REFF, DP    )
!$OMP+PRIVATE( CONST, VTS,      TEMP,  P,    PDP,  SLIP  )
!$OMP+PRIVATE( VISC,  TC0,      DELZ, DELZ1 )

      ! Loop over dust bins
      DO N = 1, NDSTBIN

         ! Initialize
         DEN   = DUSTDEN(N)
         REFF  = DUSTREFF(N)
         DP    = 2D0 * REFF * 1.D6              ! Dp [um] = particle diameter
         CONST = 2D0 * DEN * REFF**2 * G0 / 9D0

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Initialize settling velocity
               DO L = 1, LLPAR
                  VTS(L) = 0d0
               ENDDO

               ! Loop over levels
               DO L = 1, LLPAR

                  ! Get P [kPa], T [K], and P*DP
                  P    = GET_PCENTER(I,J,L) * 0.1d0
                  TEMP = T(I,J,L)
                  PDP  = P * DP

                  !=====================================================
                  ! # air molecule number density
                  ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
                  !
                  ! # gas mean free path
                  ! lamda = 1.d6 / 
                  !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
                  !
                  ! # Slip correction
                  ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * 
                  !      &  exp( -1.1 * Dp / (2. * lamda))) / Dp
                  !=====================================================
                  ! NOTE, Slip correction factor calculations following 
                  !       Seinfeld, pp464 which is thought to be more 
                  !       accurate but more computation required.
                  !=====================================================

                  ! Slip correction factor as function of (P*dp)
                  SLIP = 1d0 +
     &                   ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP

                  !=====================================================
                  ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                  ! which produce slip correction factor with small 
                  ! error compared to the above with less computation.
                  !=====================================================

                  ! Viscosity [Pa s] of air as a function of temp (K)
                  VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

                  ! Settling velocity [m/s]
                  VTS(L) = CONST * SLIP / VISC

               ENDDO

               ! fwd:
               !L           = LLTROP
               !DELZ        = BXHEIGHT(I,J,L)
               !TC(I,J,L,N) = TC(I,J,L,N) / 
               !               ( 1.d0 + DT_SETTL * VTS(L) / DELZ )
               !
               !DO L = LLTROP-1, 1, -1
               !   DELZ        = BXHEIGHT(I,J,L)
               !   DELZ1       = BXHEIGHT(I,J,L+1)
               !   TC(I,J,L,N) = 1.d0 / 
               !                 ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
               !        * (TC(I,J,L,N)  + DT_SETTL * VTS(L+1) / DELZ1
               !        *  TC(I,J,L+1,N) )
               !ENDDO
               ! adj:

               ! Save a copy of the input TC as TC0 and 
               ! intialize the output TC to 0
               TC0(:)             = TC(I,J,:,N)
               TC(I,J,1:LLTROP,N) = 0d0

               DO L = 1, LLTROP -1

                  DELZ          = BXHEIGHT(I,J,L)
                  DELZ1         = BXHEIGHT(I,J,L+1)

                  TC(I,J,L+1,N) = 1d0
     &                          / ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
     &                          * DT_SETTL * VTS(L+1) / DELZ1
     &                          * TC0(L)

                  TC(I,J,L,N)   = 1d0
     &                          / ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
     &                          * TC0(L)
     &                          + TC(I,J,L,N)

               ENDDO

               L             = LLTROP
               DELZ          = BXHEIGHT(I,J,L)
               TC(I,J,L,N)   = 1d0
     &                       / ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
     &                       * TC0(L)
     &                       + TC(I,J,L,N)


            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRY_SETTLING_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE DRY_DEPOSITION_ADJ( TC )
!
!******************************************************************************
!  Subroutine DRY_DEPOSITION_ADJ computes the adjoint of dust deposition. 
!  Deposition is linear and thus self adjoint, so we simply use the forward
!  model routine (tdf, bmy, 3/30/04, 10/25/05) modified slightly to skip 
!  the diagnostics for efficiency (dkh, 01/10/12, adj32_011) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (TYPE (XPLEX)) : Dust adjoint tracer array   
!
!  NOTES: 
!  
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE DUST_MOD,     ONLY : IDDEP
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTDST1

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      TYPE (XPLEX), INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! local variables
      INTEGER               :: I,   J,   L,      N
      TYPE (XPLEX)                :: OLD, NEW, DTCHEM

      !=================================================================
      ! DRY_DEPOSITION_ADJ begins here!
      !=================================================================

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over dust bins
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, OLD, NEW )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over dust bins
      DO N = 1, NDSTBIN

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Original dust concentration at surface
               OLD = TC(I,J,1,N)

               ! Dust left after dry deposition
               NEW = OLD * EXP( -DEPSAV(I,J,IDDEP(N)) * DTCHEM  )

               ! Save back into STT_ADJ
               TC(I,J,1,N) = NEW
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRY_DEPOSITION_ADJ

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DUST_ADJ_MOD
