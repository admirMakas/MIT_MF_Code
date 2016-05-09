!$Id: sciabr_co_obs_mod.f,v 1.3 2012/03/01 22:00:27 daven Exp $
      MODULE SCIAbr_CO_OBS_MOD

!******************************************************************************
! Module SCIAbr_CO_OBS_MOD contains subroutines necessary to 
! 1. Read SCIA Bremen (ASCII) file with CO observations (monthly data files)
! 2. Determine when SCIA CO obs are available
! 3. Transform CHK_STT into SCIA levels and then columns
! 4. Compute adjoint forcing, including transforming the difference between
!    model and MOPITT back to model space using the adjoint of averaging 
!    kernel and interpolation code.
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_SCIA_CO_FILE    : Read SCIA ASCII file
!  (2 ) ITS_TIME_FOR_SCIA_CO_OBS: function checks model time vs. OBS_HOUR array
!  (3 ) CALC_OBS_HOUR        : calculates OBS_HOUR_SCIA_CO
!  (4 ) COMPUTE_COLUMN       : vertical gridding of CHK_STT, bin data
!  (5 ) READ_ERROR_VARIANCE  : Reads error variance file
!  (6 ) CALC_SCIA_CO_FORCE   : Calculates cost function and ADJ_STT increments
!  (7 ) CLEANUP_SCIA         : Deallocates memory of arrays
!
!  ============================================================================
!  Module Variables:
!  ============================================================================
! (1 ) SCIACOcol(nobss)        : columns read from AIRS file
! (2 ) SCIACOcol_err(nobss)    :
! (3 ) COpressure(nlevs, nobss) :
! (4 ) Longitude(nobss)        : vector of longitudes read from AIRS file
! (5 ) Latitude(nobss)         : vector of latitudes read from AIRS file
! (6 ) COUNT_GRID(:,:,:)       : array of # obs/gridsquare (in 1 day)
! (7 ) SCIA_COL_GRID(:,:,:)    : gridded AIRS CO columns (computed and from file)
! (8 ) iday(:)                 : is a fraction of the day since beginning 
!                                 of the year
! (9 ) NObss                   : number of observation in each SCIA file
! (10 ) GRID_SCIA
! (11) SZA(nobs)                : vector of solar zenith angle read from file
! (12) mday(nobs)               : vector with day of month of obs
! (13) time_h(nobs)             : vector with hour of day of obs
! (13) local_t(nobs)            : vector with local time of obs
! (14) Cloud(nobs)              : vector with cloud-free(0) or contam (1)

!
!  NOTES: 
!  (1 ) Filter on MOPITT data: morning over pass only (in mop_mod.f) and
!       only use obs that are greater than 5e17 (mak, 11/18/05)    
!  (2 ) Now READ(read data, calculate OBS_HOUR_SCIA_CO), then 
!       CALC_SCIA_CO_FORCE(GRID_SCIA, compute adj forcing): (mak, 8/1/07)
!  (3 ) Fixed conflict with NOBS in CMN_ADJ, now we have NObss for SCIA.
!  (4 ) Move NLev to obs modules from CMN_ADJ. this variable controls the
!       vertical levels on which cost function is computed; NLev=1 indicates
!       we're comparing model-satellite columns (mak, 8/14/07)
!
!******************************************************************************

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=============================================================
      ! MODULE VARIABLES
      !=============================================================
      PRIVATE

      PUBLIC :: READ_SCIAbr_CO_FILE
      PUBLIC :: ITS_TIME_FOR_SCIABR_CO_OBS
      PUBLIC :: CALC_SCIABr_CO_FORCE

      TYPE (XPLEX), ALLOCATABLE  :: SCIACOcol_err(:)
      TYPE (XPLEX), ALLOCATABLE  :: SCIACOcol(:)
      TYPE (XPLEX), ALLOCATABLE  :: Longitude(:)
      TYPE (XPLEX), ALLOCATABLE  :: Latitude(:)
      TYPE (XPLEX), ALLOCATABLE  :: SZA(:)
      INTEGER, ALLOCATABLE :: COUNT_GRID(:,:,:)
      INTEGER, ALLOCATABLE :: iday(:)
      INTEGER, ALLOCATABLE :: mday(:)
      TYPE (XPLEX), ALLOCATABLE  :: time_h(:)
      INTEGER, ALLOCATABLE :: Cloud(:)
      TYPE (XPLEX), ALLOCATABLE  :: SCIA_COL_GRID(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: ERR_COL_GRID(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: CHK_STT_SCIA(:,:)
      INTEGER              :: NObss
      !INTEGER, PARAMETER   :: NLev  = 1
      INTEGER, PARAMETER   :: NLevs = 60
      INTEGER              :: NDays
      INTEGER              :: fId
      INTEGER              :: TRCNUM
      INTEGER, ALLOCATABLE :: OBS_HOUR_SCIA_CO(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: FRACTION(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: AirMass(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: ADJ_SCIA_ALL(:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: ERR_PERCENT(:,:)
      INTEGER, ALLOCATABLE :: DOMAIN_OBS(:,:)
 
      CONTAINS

      SUBROUTINE READ_SCIAbr_CO_FILE( YYYYMMDD, HHMMSS )

!******************************************************************************
!  Subroutine READ_SCIA_FILE reads the SCIA ASCII file and assigns OBS_HOUR
!  array based on available data. SCIA CO data are stored in a 1 month/file
!  frequency. (mak, 7/12/07)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!
      USE ErrorModule, ONLY : ReplaceNanAndInf
      USE ERROR_MOD,   ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,    ONLY : IOERROR
      USE TIME_MOD,    ONLY : EXPAND_DATE, YMD_EXTRACT, GET_MONTH

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, L, N
      CHARACTER(LEN=255) DIR_SCIA
      CHARACTER(LEN=255) FILENAME_SCIA
      CHARACTER(LEN=255) FILENAME
      
 
      integer ipx,ist,iread,it
      TYPE (XPLEX) dsr_time,t_int,lat_c,lon_c,lat_1
      TYPE (XPLEX) lon_1, lat_2,lon_2,lat_3,lon_3,lat_4,lon_4
      TYPE (XPLEX) sza_in,los,azi
      integer cld,lnd
      TYPE (XPLEX) rms,snrad,alt,h20,h20_err,ch4,ch4_err,co,co_err
      TYPE (XPLEX) co_corr, co_corr_err
      integer coq
      integer counter, count_day
      INTEGER    :: IU_FILE, IOS, IOS1
      CHARACTER(LEN=255)   :: header
      INTEGER    :: LINECOUNT
      INTEGER    :: YEAR, MONTH, DAY
      INTEGER    :: DAYOM_year(12), days_so_far

      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: LASTMONTH = -999

      !=============================
      ! FIRST CLEANUP IF NECESSARY:
      !=============================
      CALL CLEANUP_SCIA

      !========================
      ! FILENAME
      !=========================

      DIR_SCIA = '/as/data/scia/monthly_good/'
      !DIR_SCIA = '/as/home/ctm/mak/sciamachy/data_scia/monthly_good/'
      !DIR_SCIA = '/as/home/ctm/jaf/SCIA/data/'
      FILENAME_SCIA = 'SCI_WFMD_L2_w8003_YYYYMM_v0.6.was'
      IU_FILE = 15
      days_so_far = -1
      DAYOM_year(:) = -1


      MONTH = YYYYMMDD/100 - (YYYYMMDD/10000)*100
      YEAR  = YYYYMMDD/10000

      ! Determine Ndays, number of days in the given month
      SELECT CASE (MONTH)
      CASE(4,6,9,11)
         Ndays = 30
      CASE(2)
         IF (YEAR .eq. 2004) Ndays = 29
         IF (YEAR .ne. 2004) Ndays = 28
      CASE DEFAULT
         Ndays = 31
      END SELECT
      
      IF( YEAR ==2004) THEN
         DAYOM_year = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      ELSEIF (YEAR ==2005) THEN
         DAYOM_year = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      ENDIF
      days_so_far = sum(DAYOM_year(1:MONTH-1))

      CALL EXPAND_DATE( FILENAME_SCIA, YYYYMMDD, 0 )
      FILENAME = trim(DIR_SCIA)//FILENAME_SCIA

      !print*, 'filename:', FILENAME

      ! zero counters
      counter = 0
      count_day = 0
      LINECOUNT = 0

      ! Figure out how many observations to read (#lines in the file):
      CALL SYSTEM('wc -l '//trim(fileName)//' > tmp.txt')
   
      OPEN( 5, FILE='tmp.txt', IOSTAT=IOS1 )
      IF ( IOS1 /= 0 ) CALL IOERROR( IOS1, 5, 'tmp:1' )
      
      ! Read #lines
      READ( 5, *, IOSTAT=IOS1  ) linecount
      IF ( IOS1 /= 0 ) CALL IOERROR( IOS1, 5, 'tmp:2' )
   
      ! Close file
      CLOSE( 5 )

      !PRINT*, 'There are:', LINECOUNT-39, 'good observations'
      IF ( LINECOUNT == 0 ) THEN
         PRINT*, 'There are no obs available, exit'
         CALL GEOS_CHEM_STOP
      ENDIF

      OPEN( IU_FILE, FILE=fileName, IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'scia:1' )
      
      NObss = linecount-39

      ! Read SCIA header (first 41 lines, 39 lines for monthly files)
      READ(IU_FILE,'(a)') HEADER
      WRITE(*,*) TRIM( HEADER )

      ! rest of header
      do i=2,39 
         read(IU_FILE, *)
      enddo
   
      !=========================================================
      ! ALLOCATE ARRAYS: (now that we read in their dimensions!)
      !=========================================================

      ALLOCATE( iday( NObss ) )
      iday(:) = 0d0
      ALLOCATE( mday( NObss ) )
      mday(:) = 0d0
      ALLOCATE( SCIACOcol( NObss ) )
      SCIACOcol(:) = 0d0
      ALLOCATE( SCIACOcol_err( NObss ) )
      SCIACOcol_err(:) = 0d0
      ALLOCATE( Longitude( NObss ) )
      Longitude(:) = 0d0
      ALLOCATE( Latitude( NObss ) )
      Latitude(:) = 0d0
      ALLOCATE( SZA( NObss ) )
      SZA(:) = 0d0
      ALLOCATE( time_h(NObss ) ) 
      time_h(:) = 0d0
      ALLOCATE( Cloud( NObss ) )
      Cloud(:) = 0d0

      !compute iday based on the input file, meaning August 1, for july file
      ! then subtract iday_input - iday = day fraction within the given month
      ! multiply the fraction times the NDAYS (from above), get, e.g. 15.5
      ! meaning it's July 15 at noon. use the modulus command to get day and
      ! hour

      DO i = 40, NObss+39
     
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'scia:2' )

         ! Read SCIA data
         read( IU_FILE,
     &        '(i6,i4,i6,i2,f14.9,f9.5,10f10.5,3f8.3,2i4,11e13.9,i3)', 
     &   IOSTAT=IOS)  ipx,ist,iread,it,dsr_time,t_int,lat_c,lon_c,lat_1,
     &   lon_1, lat_2,lon_2,lat_3,lon_3,lat_4,lon_4,sza_in,los,azi,cld,
     &   lnd, rms,snrad,alt,h20,h20_err,ch4,ch4_err,co,co_err,co_corr, 
     &       co_corr_err,coq

         ! dsr_time is Starttime in frac.days since 1.1.2000
         ! time_h is obs hour for each obs
         if(year == 2004) then
            iday(i-39)=dsr_time-366-365*3+1
            time_h(i-39)=(dsr_time-366-365*3-iday(i-39)+1)*24
            Cloud(i-39) = cld
         ELSEIF(year == 2005 ) then
            iday(i-39)=dsr_time-366*2-365*3+1
            time_h(i-39)=(dsr_time-366*2-365*3-iday(i-39)+1)*24
            Cloud(i-39) = cld
         else
            print*, 'year is:', year
            print*, 'only looking at 2004 and 2005 data'
            CALL GEOS_CHEM_STOP
         endif
         mday(i-39)=iday(i-39)-days_so_far !gives day of the month
 
         ! if it's a good data (quality flag coq != 1) 

         SCIACOcol(i-39) = co
         SCIACOcol_err(i-39) = co_err
         Latitude(i-39) = lat_c
         Longitude(i-39) = lon_c !-180 ! for -180..180 range
         SZA(i-39) = sza_in
         
         counter = counter+1
         
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Echo info
      !PRINT*, '### number of observations (total): ', counter
      !PRINT*, 'NObss: ', NObss
!      PRINT*, '### Longitude min, max: ', MINVAL( Longitude ), 
!     & MAXVAL( Longitude )
!      PRINT*, '### Latitude min, max: ', MINVAL( Latitude ), 
!     & MAXVAL( Latitude )
      !PRINT*, '### SZA min, max: ', MINVAL( SZA ), MAXVAL( SZA )
      !PRINT*, '### COcol min, max:', MINVAL(SCIACOcol),MAXVAL(SCIACOcol)
!      PRINT*, '### COcol_err min, max:', MINVAL(SCIACOcol_err),
!     & MAXVAL(SCIACOcol_err)
      !print*, 'min/max of day of the month:', minval(mday), maxval(mday)
      !print*, 'min/max hour of day:', minval(time_h), maxval(time_h)
 
      ! Decide how many arrays to read/store, for now 2 CO col computed
      ! and CO col read from file

      TRCNUM = 1

      !=========================================================
      ! ALLOCATE ARRAYS: (now that we read in their dimensions!)
      !=========================================================

      ! Allocate 2 CO column matrices for CO columns computed and read in 
      !from file
      ALLOCATE( SCIA_COL_GRID(IIPAR, JJPAR, Ndays,TRCNUM))
      SCIA_COL_GRID(:,:,:,:) = 0d0
      ALLOCATE( ERR_COL_GRID(IIPAR, JJPAR, Ndays))
      ERR_COL_GRID(:,:,:) = 0d0
      ALLOCATE ( COUNT_GRID(IIPAR, JJPAR, Ndays))
      COUNT_GRID(:,:,:) = 0
      ALLOCATE( OBS_HOUR_SCIA_CO(IIPAR, JJPAR, NDAYS ))
      OBS_HOUR_SCIA_CO(:,:,:) = 0d0
      ALLOCATE( ERR_PERCENT(IIPAR,JJPAR) ) 
      ERR_PERCENT(:,:) = 0d0

      ! only compute SCIA_OBS_HOUR; grid when computing forcing
      CALL CALC_OBS_HOUR

      CALL INIT_DOMAIN

      ! READ ERROR FILE
      IF(GET_MONTH() .NE. LASTMONTH) THEN
         CALL READ_ERROR_VARIANCE
         LASTMONTH = GET_MONTH()
      ENDIF

      END SUBROUTINE READ_SCIAbr_CO_FILE

!----------------------------------------------------------------------

      SUBROUTINE COMPUTE_COLUMN
!*************************************************************************
! This subroutine computes SCIA column, based on the averaging kernels
! COverticality from SCIA file. (jaf, 7/07) It now includes vertical
! regridding, previously done in IDL ( mak, 8/2/07)
!
! Notes:
! (1 ) The subroutine does vertical regridding of GC column. 
!      Then reads averaging kernels, a priori guess and pressure levels
!      Then it does retrieval to get a GC column * SCIA AK.
! (2 ) CHK_STT(IIPAR,JJPAR,LLPAR,1) -> Model_CO_MR(IIPAR,JJPAR,LLPAR)->
!      ->CHK_STT_SCIA_VGRID(IIPAR,JJPAR,NLevs)->CHK_STT_STT(IIPAR,JJPAR)
!*************************************************************************

      USE FILE_MOD, ONLY    : IOERROR
      USE TIME_MOD, ONLY    : GET_HOUR, GET_DAY
      USE DAO_MOD,  ONLY    : AD
      USE TRACER_MOD, ONLY  : TCVV

#     include "CMN_SIZE"   ! Size parameters

      TYPE (XPLEX)     :: T(NLevs)
      TYPE (XPLEX)     :: delP(NLevs)
      TYPE (XPLEX)     :: Pedge(NLevs)
      TYPE (XPLEX)     :: Pa(NLevs+1)
      TYPE (XPLEX)     :: sciapress(NLevs)
      TYPE (XPLEX)     :: xa(NLevs)
      TYPE (XPLEX)     :: ap(NLevs)
      TYPE (XPLEX)     :: AK(NLevs)
      TYPE (XPLEX)     :: A(NLevs)
      
      TYPE (XPLEX)     :: alt(NLevs)
      TYPE (XPLEX)     :: temp(NLevs)
      TYPE (XPLEX)     :: sza20(NLevs), sza30(NLevs), sza40(NLevs)
      TYPE (XPLEX)     :: sza50(NLevs), sza60(NLevs), sza65(NLevs)
      TYPE (XPLEX)::sza70(NLevs),sza75(NLevs),sza80(NLevs),sza85(NLevs)
      
      INTEGER    :: L, N, I, J, K, D, w, LL
      INTEGER    :: IU_FILE, IOS
      INTEGER    :: ILON, ILAT
      CHARACTER(LEN=255) ::  HEADER
      LOGICAL, SAVE :: FIRST = .TRUE.

      TYPE (XPLEX)  :: Model_CO_MR(IIPAR, JJPAR, NLevs)

      xa(:) = 0d0
      ap(:) = 0d0
      Pa(:) = 0d0
      sciapress(:) = 0d0
      T(:) = 0d0
      delP(:) = 0d0
      PEdge(:) = 0d0

     
      IF ( FIRST ) THEN
         ALLOCATE( CHK_STT_SCIA( IIPAR,JJPAR ) )
         CHK_STT_SCIA(:,:) = 0d0
         ALLOCATE( FRACTION(IIPAR,JJPAR,LLPAR,NLevs))
         FRACTION(:,:,:,:) = 0e0
         ALLOCATE( Airmass(IIPAR,JJPAR,NLevs))
         Airmass(:,:,:) = 0d0
         ALLOCATE( ADJ_SCIA_ALL(IIPAR,JJPAR,LLPAR))
         ADJ_SCIA_ALL= 0d0
         FIRST = .FALSE.
      ELSE
         CHK_STT_SCIA(:,:) = 0d0    
         FRACTION(:,:,:,:) = 0e0   
         Airmass(:,:,:) = 0d0  
         ADJ_SCIA_ALL(:,:,:)= 0d0
      ENDIF

      Model_CO_MR(:,:,:) = 0d0

      !       xcol =    AK*geos_raw   + (I - AK)          *xapriori
      !  OR   xcol =    AK*COppbv     + (I-COverticality) *xapriori
      !  OR   model_col(1) = T.xa + A.(geos_raw - xa)
      ! where T is COverticality and xa is read from file
      ! READ AIRS FIRST GUESS (A PRIORI) FROM A TEXT FILE:

      !--------------------------------------------------------
      ! Open file with SCIA a priori (xa) and averaging kernels
      !--------------------------------------------------------

      IU_FILE = 15

      OPEN( IU_FILE, FILE='ak_co_wfmdscia_V2.dat', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'aveker:1' )
      
      ! Read SCIA AK file header (first 4 lines)
      READ(IU_FILE,'(a)') HEADER
      !WRITE(*,*) TRIM( HEADER )

      ! rest of header
      DO I = 2, 4 
         READ(IU_FILE, *)
      ENDDO
      
      DO I = 5, NLevs+4
      
         L = I - 4
         ! Read SCIA info
         ! Note: Pressure (Pa) is in hPa, a priori (xa) is ppmv!!!
         READ( IU_FILE, 100, IOSTAT=IOS ) alt(L), sciapress(L),  
     &    temp(L), ap(L), sza20(L), sza30(L), sza40(L),       
     &    sza50(L), sza60(L), sza65(L), sza70(L), sza75(L),   
     &    sza80(L), sza85(L)
         
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'aveker:2' )
 100     FORMAT(f7.2,2f8.2,f14.5,2x,10f8.3)
      
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! ak file lists first row starting w/ highest altitude
      ! Need to reverse order for Pa, xa; AK reversal done below
      Pa(1:NLevs) = sciapress(NLevs:1:-1)
      xa(1:NLevs) = ap(NLevs:1:-1)

      ! Assign top pressure (not given in AK file).
      ! For this, I assume the first pressure given is halfway
      ! between the pressure below (given) and the pressure
      ! above (not given)
      Pa(NLevs+1) = 2*Pa(NLevs) - Pa(NLevs-1)

      !---------------------------------------------------------
      ! GEOS-Chem profile on LLPAR levels, convert to NLevs =60
      !---------------------------------------------------------
      ! have: GC profile for particular day of the month and for particular
      ! hour in the day CHK_STT(I,J,L)
      CALL REGRIDV_SCIA( Model_CO_MR ) 
      print*, 'max/min of Model_CO_MR:', maxval(Model_CO_MR), 
     &     minval(Model_CO_MR)
!       print*, 'Model_CO_MR location of min',minloc(Model_CO_MR),'max ',
!     & MAXLOC(Model_CO_MR)


      ! CHK_STT is kg, since stored that way in geos_chem_mod
      ! convert to v/v from kg
c$$$!$OMP PARALLEL DO
c$$$!$OMP+DEFAULT( SHARED )
c$$$!$OMP+PRIVATE( I,   J,   L)
c$$$!$OMP+SCHEDULE( DYNAMIC )
c$$$      DO L= 1,NLevs
c$$$      DO J= 1,JJPAR
c$$$      DO I= 1,IIPAR
c$$$         Model_CO_MR(I,J,L) = CHK_STT_SCIA_VGRID(I,J,L)
c$$$     &                                   * ADJ_TCVV(1) / AD(I,J,L)
c$$$      ENDDO
c$$$      ENDDO
c$$$      ENDDO
c$$$!$OMP END PARALLEL DO

      !--------------------------------------------------------
      ! transfer function [molec/cm2/ppmv]
      !DO L = 1, NLevs
      !   t(L) = 2.12E+16*delP(L)  ! FOR SCIA (since P is in hPa and have to 
      !ENDDO                       !           convert to 1/cm2 from 1/m2
                                   !           and a priori in ppmv)      
      ! column version of the averaging kernel
      ! A in molec/cm2/ppmv (avgker: unitless)

      ! For each observation 1 to NObss:   
      DO w = 1, NObss            ! loop over # of obs

         IF(floor((time_h(w)%r)) == GET_HOUR() .AND. 
     &      mday(w) == GET_DAY() .and.
     &      Cloud(w) == 0 ) THEN

            ! Look for model grid box corresponding to the SCIA observation:
            ! Get I corresponding to PLON(IND)
            ILON = INT( ( Longitude(w) + 180d0 ) / DISIZE + 1.5d0 )
            ! Handle date line correctly (bmy, 4/23/04)
            IF (ILON > IIPAR ) ILON = ILON - IIPAR 
            ! Get J corresponding to PLAT(IND)
            ILAT = INT( ( Latitude(w) +  90d0 ) / DJSIZE + 1.5d0 )

            ! Initialize AK
            AK(:) = 0d0

            ! Get correct averaging kernel for given SZA of observation
            IF ( (SZA(w) .ge. 0   ) .AND. (SZA(w) .lt. 25  ) ) 
     &           AK = sza20(NLevs:1:-1)
            IF ( (SZA(w) .ge. 25  ) .AND. (SZA(w) .lt. 35  ) ) 
     &           AK = sza30(NLevs:1:-1)
            IF ( (SZA(w) .ge. 35  ) .AND. (SZA(w) .lt. 45  ) ) 
     &           AK = sza40(NLevs:1:-1)
            IF ( (SZA(w) .ge. 45  ) .AND. (SZA(w) .lt. 55  ) ) 
     &           AK = sza50(NLevs:1:-1)
            IF ( (SZA(w) .ge. 55  ) .AND. (SZA(w) .lt. 62.5) ) 
     &           AK = sza60(NLevs:1:-1)
            IF ( (SZA(w) .ge. 62.5) .AND. (SZA(w) .lt. 67.5) ) 
     &           AK = sza65(NLevs:1:-1)
            IF ( (SZA(w) .ge. 67.5) .AND. (SZA(w) .lt. 72.5) ) 
     &           AK = sza70(NLevs:1:-1)
            IF ( (SZA(w) .ge. 72.5) .AND. (SZA(w) .lt. 77.5) ) 
     &           AK = sza75(NLevs:1:-1)
            IF ( (SZA(w) .ge. 77.5) .AND. (SZA(w) .lt. 82.5) ) 
     &           AK = sza80(NLevs:1:-1)
            IF ( (SZA(w) .ge. 82.5) .AND. (SZA(w) .le. 90  ) ) 
     &           AK = sza85(NLevs:1:-1)

            ! Check that averaging kernel has been assigned
            IF(MAXVAL(AK) .eq. 0) PRINT*, 'No averaging kernel assigned'

            !print*, 'scia col:', SCIACOcol(w)

            DO L = 1, NLevs     ! loop over # of levels

               delP(L) = Pa(L) - Pa(L+1)

               ! Compute transfer function for particular level
               T(L) = 2.12E+16*delP(L)
 
               ! If we have all info, compute the CO column
               IF ( (T(L) .gt. 0) .AND. (AK(L) .gt. 0) ) THEN

                  A(L) = T(L) * AK(L)
                  ! SCIA a priori info (xa) is in ppmv
                  ! Model profile is in v/v
                  CHK_STT_SCIA(ILON,ILAT) =  
     &                 CHK_STT_SCIA(ILON,ILAT) 
     &                 + T(L) * xa(L) + A(L) *  
     &                 (Model_CO_MR(ILON,ILAT,L)*1e6 - xa(L) )

               ENDIF ! If we have all info
               
            ENDDO   ! Loop over levels

!            print*, 'min max airmass:', minval(airmass), maxval(airmass)
!            print*,'min max fraction:',minval(fraction),maxval(fraction)
!            print*, 'min  max ADJ_SCIA_ALL:', 
!     &           minval(ADJ_SCIA_ALL(ILON,ILAT,:)),
!     &           MAXVAL(ADJ_SCIA_ALL(ILON,ILAT,:))
!            PRINT*, 'ADJ_TCVV:', ADJ_TCVV(1)
!            print*, 'min max a:', minval(a), maxval(a)
            DO L = 1,LLPAR
               DO LL = 1,NLevs
                  ! d(CHK_STT_SCIA)/d(Model_CO_MR) = A(LL)*1e6
                  IF (AirMass(ILON,ILAT,LL) .GT. 0) THEN
                  ADJ_SCIA_ALL(ILON,ILAT,L) = ADJ_SCIA_ALL(ILON,ILAT,L)+
     &                 (A(LL)*1e6) *(TCVV(1)/AirMass(ILON,ILAT,LL))
     &                 *FRACTION(ILON,ILAT,L,LL)
                  ENDIF
               ENDDO
            ENDDO

!            PRINT*, 'ADJ_SCIA_ALL:', ADJ_SCIA_ALL(ILON,ILAT,:)

            ! SCIA column 
            SCIA_COL_GRID(ILON,ILAT,mday(w),1) = 
     &           SCIA_COL_GRID(ILON,ILAT,mday(w),1) + SCIACOcol(W)/0.9

            ! SCIA column error
            ERR_COL_GRID(ILON,ILAT,mday(w)) = 
     &           ERR_COL_GRID(ILON,ILAT,mday(w))+ SCIACOcol_err(w)
            !print*, 'model col:', CHK_STT_SCIA(ILON,ILAT)

         ENDIF !time and day of the model is the same as this obs
      ENDDO   ! Loop over observations

      !=======================================
      ! BIN OUTPUT INFO INTO MODEL GRID BOXES
      !=======================================

      ! Getting a day of the month ok if one month at a time
      D = GET_DAY()

!      print*, 'before averaging:'
!      PRINT*, '### SCIA_COL_GRID  min, max: ', 
!     & MINVAL( SCIA_COL_GRID(:,:,D,1) ), MAXVAL( SCIA_COL_GRID(:,:,D,1))
!      PRINT*, '### CHK_STT_SCIA  min, max: ', 
!     & MINVAL( CHK_STT_SCIA ), MAXVAL( CHK_STT_SCIA )
!      PRINT*, '### ERR_COL_GRID  min, max: ', 
!     & MINVAL( ERR_COL_GRID(:,:,D)  ), MAXVAL( ERR_COL_GRID(:,:,D))

      DO I = 1, IIPAR
      DO J = 1, JJPAR
         IF ( COUNT_GRID(I,J,D) .gt. 0. ) then   
            ! average SCIA
            SCIA_COL_GRID(I,J,D,1) = SCIA_COL_GRID(I,J,D,1)/
     &           COUNT_GRID(I,J,D)

            ! average SCIA error
            ERR_COL_GRID(I,J,D) = ERR_COL_GRID(I,J,D)/COUNT_GRID(I,J,D)

            ! average model
            CHK_STT_SCIA(I,J) = CHK_STT_SCIA(I,J)/COUNT_GRID(I,J,D) 

            DO L = 1,LLPAR
             ! d(CHK_STT_SCIA)/d(Model_CO_MR) = A(LL)*1e6
             ADJ_SCIA_ALL(I,J,L) = ADJ_SCIA_ALL(I,J,L)/COUNT_GRID(I,J,D)
            ENDDO
        ELSE
           SCIA_COL_GRID(I,J,D,:) = -999.
           ERR_COL_GRID(I,J,D) = -999.
           CHK_STT_SCIA(I,J) = -999.
           ADJ_SCIA_ALL(I,J,:) = 0d0
        ENDIF
      ENDDO
      ENDDO

      PRINT*, '### SCIA_COL_GRID  min, max: ', 
     & MINVAL( SCIA_COL_GRID(:,:,D,1) ), MAXVAL( SCIA_COL_GRID(:,:,D,1))
      PRINT*, '### CHK_STT_SCIA  min, max: ', 
     & MINVAL( CHK_STT_SCIA ), MAXVAL( CHK_STT_SCIA )
!      PRINT*, '### ERR_COL_GRID  min, max: ', 
!     & MINVAL( ERR_COL_GRID(:,:,D)  ), MAXVAL( ERR_COL_GRID(:,:,D))
!      PRINT*, '### ADJ_SCIA_ALL min, max:', 
!     & MINVAL( ADJ_SCIA_ALL ) , MAXVAL(ADJ_SCIA_ALL)
!      call flush(6)
!      PRINT*, '### COUNT_GRID min, max: ', 
!     & MINVAL( COUNT_GRID(:,:,:)  ), MAXVAL( COUNT_GRID(:,:,:))
!      PRINT*, 'TODAY:'
!      PRINT*, '### COUNT_GRID min, max: ', 
!     & MINVAL( COUNT_GRID(:,:,D) ), MAXVAL( COUNT_GRID(:,:,D) )

      END SUBROUTINE COMPUTE_COLUMN

!----------------------------------------------------------------------

      SUBROUTINE READ_ERROR_VARIANCE 

      USE ERROR_MOD, 	ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE BPCH2_MOD
      USE FILE_MOD,    ONLY : IOERROR
      USE TIME_MOD,    ONLY : GET_TAU

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      CHARACTER(LEN=255)     :: INPUT_FILE
      INTEGER                :: I, IOS, J,  L
      CHARACTER(LEN=255)     :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR)
      INTEGER              :: HALFPOLAR 
      INTEGER              :: CENTER180 
      INTEGER              :: NI,     NJ,     NL, k
      INTEGER              :: IFIRST, JFIRST, LFIRST
      INTEGER              :: NTRACER,   NSKIP
      TYPE (XPLEX)               :: ZTAU0,     ZTAU1, TAUTMP

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      INTEGER              :: IU_FILE
      LOGICAL :: IT_EXISTS

      INPUT_FILE = 'RRE_seasonMay1sciabrGlobal.bpch'
      IU_FILE = 66

      PRINT*, 'SET ERROR TO 30% AS WE SAVE SCIA FOR RRE CALCULATION'
      ERR_PERCENT(:,:) = 0.3
      GOTO 121

      !================================================================
      ! READ OLD RESTART FILE
      !================================================================
      FILENAME = TRIM( INPUT_FILE )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'O B S   E R R O R   F I L E'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_FILE: Reading ', a )

      ERR_PERCENT(:,:) = -999d0

      !READ SEASONAL ERRORS:
      IF (((GET_TAU() .ge. 169440.) .and.(GET_TAU().lt. 170184.)) 
     &.or.((GET_TAU() .ge. 176736.) .and.(GET_TAU().lt. 178200.)))THEN
         ! if May 2004, March 2005 or April 2005: read spring error
         TAUTMP = 169440.00
      ELSEIF((GET_TAU() .ge. 170184.).and.(GET_TAU().lt. 172392.))THEN
         ! if June 2004 through August 2004, read summer error
         TAUTMP = 170184.00
      ELSEIF((GET_TAU() .ge. 172392.).and.(GET_TAU().lt. 174576.))THEN
         ! if September 2004 through November 2004, read fall error
         TAUTMP = 172392.00
      ELSEIF((GET_TAU().ge. 174576.).and.(GET_TAU() .lt. 176736.))THEN
         TAUTMP = 174576.00
      ELSE
         PRINT*, 'missing obs error for tau:', GET_TAU()
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )
      ! Echo more output

      !=================================================================
      ! Read tracers -- store in the TRACER array
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
      
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a TYPE (XPLEX) I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_FILE,'read_file:4' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_file:5')

         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( TRACER(I,J), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_file:6')

         !==============================================================
         ! Assign data from the TRACER array to the ERR_PERCENT array. 
         !==============================================================
         IF ( ZTAU0 == TAUTMP ) THEN 
            PRINT*, 'Reading error for tau:', ztau0
         ! Make sure array dimensions are of global size
         ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            !print*, 'inside the loop'
            !print*, 'max value is:', maxval(TRACER(:,:,:))
            !print*, 'min value is:', minval(TRACER)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J , L)
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF(TRACER(I,J) .ge. 0.05)THEN
               ERR_PERCENT(I,J) = TRACER(I,J)
               
            ELSEIF((TRACER(I,J) .lt. 0.05).and.(TRACER(I,J).gt.0)) THEN
               ERR_PERCENT(I,J) = 0.05
               
            ELSE
               ERR_PERCENT(I,J) = TRACER(I,J)
   
            ENDIF
               
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ENDIF     
      ENDDO

      ! Close file
      CLOSE( IU_FILE )      

 121  CONTINUE

      print*, 'max value is:', maxval(ERR_PERCENT(:,:))
      print*, 'min value is:', minval(ERR_PERCENT)

      CALL FLUSH(6)
      END SUBROUTINE READ_ERROR_VARIANCE

!---------------------------------------------------------------------------

      SUBROUTINE GRID_SCIA
!*********************************************************************
!   GRIDS SCIA ARRAYS ONTO GEOS-CHEM GRID OF CHOICE (e.g GEOS3 2x2.5)
!*********************************************************************    

# include "CMN_SIZE"

      INTEGER             :: W, ILON, ILAT, I, J, M, D, count_day


      !=====================================================
      ! GRID_SCIA begins here
      !=====================================================

      ! at time point while looping over all obs we need to figure out
      ! what day is the obs from to assign it to the correct 3rd dimension.
      count_day = 0

      DO W = 1, NOBSS
         ! COMPUTE LONGITUDE AND LATITUDE
         !=================================================================
         ! Look for model grid box corresponding to the MOPITT observation:
         ! Get I corresponding to PLON(IND)
         ILON = INT( ( Longitude(W) + 180d0 ) / DISIZE + 1.5d0 )
         ! Handle date line correctly (bmy, 4/23/04)
         IF (ILON > IIPAR ) ILON = ILON - IIPAR 
         ! Get J corresponding to PLAT(IND)
         ILAT = INT( ( Latitude(W) +  90d0 ) / DJSIZE + 1.5d0 )

         ! IF valid column obs, increment the counter for that box
         IF(SCIACOcol(W) .gt. 0) THEN
            D = mday(w)

            DO I = 1, TRCNUM
               IF (I .eq. 1 ) THEN
                  SCIA_COL_GRID(ILON,ILAT,D,I) = 
     &                 SCIA_COL_GRID(ILON,ILAT,D,I) + SCIACOcol(W)
                  ERR_COL_GRID(ILON,ILAT,D) = 
     &                 ERR_COL_GRID(ILON,ILAT,D) + SCIACOcol_err(w)
               ELSEIF ( I .eq. 2 ) THEN
!                  CO_COL_GRID(ILON,ILAT,D,I) = 
!     &                 CO_COL_GRID(ILON,ILAT,D,I) + 
!     &                 (COcol(W) - COcol_pro(W))
               ENDIF
            ENDDO
            
         ENDIF
         
         count_day = count_day+1
      ENDDO  ! Loop over obs
 
      print*,'number of obs this month',count_day

! print*,'Bin output model and scia data into array....'

      !=======================================
      ! BIN OUTPUT INFO INTO MODEL GRID BOXES
      !=======================================

      DO D = 1, NDays
      DO I = 1, IIPAR
      DO J = 1, JJPAR
         IF ( COUNT_GRID(I,J,D) .gt. 0. ) then
            DO M = 1, TRCNUM
               SCIA_COL_GRID(I,J,D,M) = SCIA_COL_GRID(I,J,D,M)/
     &              COUNT_GRID(I,J,D)
            IF(M == 1) THEN
            ERR_COL_GRID(I,J,D) = ERR_COL_GRID(I,J,D)/COUNT_GRID(I,J,D)
            ENDIF
           ENDDO
        ELSE
           SCIA_COL_GRID(I,J,D,:) = -999.
           ERR_COL_GRID(I,J,D) = -999.
        ENDIF
      ENDDO
      ENDDO
      ENDDO ! Loop over days

      PRINT*, '### SCIA_COL_GRID  min, max: ', 
     & MINVAL( SCIA_COL_GRID(:,:,:,1) ), MAXVAL( SCIA_COL_GRID(:,:,:,1))
      PRINT*, '### ERR_COL_GRID  min, max: ', 
     & MINVAL( ERR_COL_GRID(:,:,:)  ), MAXVAL( ERR_COL_GRID(:,:,:))

      END SUBROUTINE GRID_SCIA

!-------------------------------------------------------------------------
      
      FUNCTION ITS_TIME_FOR_SCIAbr_CO_OBS( ) RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_MOPITT_OBS returns TRUE if there are observations
!  available for particular time (hour of a particular day) based on
!  the OBS_HOUR array which holds the hour of obs in each gridbox (computed
!  when file read in mop02_mod.f) (mak, 7/12/07)
!
!  NOTES:
!
!******************************************************************************
!

      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE, GET_DAY

#     include "CMN_SIZE"  ! Size params

      ! Function value
      LOGICAL :: FLAG

      INTEGER :: I,J, D

      !=================================================================
      ! ITS_TIME_FOR_SCIAbr_CO_OBS begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.

      DO J = 1,JJPAR
      DO I = 1,IIPAR
         IF(GET_HOUR() == OBS_HOUR_SCIA_CO(I,J,GET_DAY()) 
     &        .AND. GET_MINUTE() == 0  ) THEN
                print*, 'obs_hour was', GET_HOUR(), 'in box', i, j
            FLAG = .TRUE.
            GOTO 11
         ENDIF
      ENDDO 
      ENDDO
      
 11   CONTINUE
      END FUNCTION ITS_TIME_FOR_SCIAbr_CO_OBS

!---------------------------------------------------------------------------

      SUBROUTINE CALC_SCIAbr_CO_FORCE

      ! References to F90 modules
      USE ERROR_MOD, 		ONLY : IT_IS_NAN, ERROR_STOP
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TIME_MOD,             ONLY : GET_HOUR, GET_NYMDe, GET_NHMSe,
     &                                 GET_DAY
      USE ADJ_ARRAYS_MOD,       ONLY : SET_OBS,SET_MODEL,SET_MODEL_BIAS,
     &                                 SET_FORCING,COST_ARRAY,OBS_COUNT,
     &                                 STT_ADJ, COST_FUNC, IFD, JFD, 
     &                                 LFD, NFD, ADJ_FORCE,
     &                                 DAY_OF_SIM
       ! if no AK, need access to STT
      USE TRACER_MOD,           ONLY : STT
      USE LOGICAL_ADJ_MOD,      ONLY : LPRINTFD, LDCOSAT

#     include "CMN_SIZE" 	! Size parameters

      ! Internal variables 
      TYPE (XPLEX)  :: DIFF_COST
      TYPE (XPLEX)  :: NEW_COST(IIPAR,JJPAR) 
      INTEGER :: I, J, L, N, LL
      INTEGER :: ADJ_EXPLD_COUNT
      INTEGER, PARAMETER  ::  MAX_ALLOWED_EXPLD    = 10
      TYPE (XPLEX),  PARAMETER::MAX_ALLOWED_INCREASE=xplex(10D15,0d0)
      TYPE (XPLEX)  :: MAX_ADJ_TMP
      TYPE (XPLEX)  :: invSy(IIPAR,JJPAR) !error variance for column
      INTEGER :: DAYOM 
      TYPE (XPLEX)  :: Sy
      LOGICAL :: USING_AK = .TRUE. 

      !================================================================ 
      ! CALC_SCIAbr_CO_FORCE begins here!
      !================================================================
      
      !print*, 'in CALC_SCIA_CO_FORCE'

      ! Some error checking stuff
      MAX_ADJ_TMP     = MAXVAL( STT_ADJ )
      ADJ_EXPLD_COUNT = 0
      Sy = 0d0

      !initialize:
      NEW_COST(:,:) = 0d0

      ! reinitialize domain
      CALL INIT_DOMAIN

      IF ( USING_AK ) THEN
         ! grid scia and compute GC*AK (column value as CHK_STT_SCIA)
         ! COL obs from GC to compare to scia
         print*, 'using averaging kernels'
         CALL COMPUTE_COLUMN
      ELSE
         ! NO AVERAGING KERNELS
         print*, 'not using averaging kernels'
         IF(.not. (ALLOCATED(  CHK_STT_SCIA) )) THEN 
            ALLOCATE( CHK_STT_SCIA( IIPAR,JJPAR ) )
            CHK_STT_SCIA(:,:) = 0d0
         ENDIF
         CALL GRID_SCIA

         ! compute straight column, no averaging kernels for now
         !CHK_STT_SCIA(:,:) = SUM(CHK_STT,3)
         ! this gives us a GC column in kg
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N)
         DO N = 1, 1            !FOR NOW, UNTIL WE have more than just CO obs
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            CHK_STT_SCIA(I,J) = CHK_STT_SCIA(I,J) + STT(I,J,L,N)
         ENDDO
         ! kg -> molec/cm2 conversion
         ! molec/cm2 = kg * 1000g/kg / (28g/mole) * 6.02 *10^23 molec/mole
         ! * (1/GET_AREA_CM2)
         CHK_STT_SCIA(I,J) = CHK_STT_SCIA(I,J) * 1000 / 28 * 6.02 * 1e23 
     &        * (1/GET_AREA_CM2(J))
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         print*, 'min/max GC/scia column:', minval(chk_stt_scia), 
     &        maxval(chk_stt_scia)

      ENDIF !using AK

      ! DAY of the month:
      DAYOM = GET_DAY()

      ! Compute error for each day of the obs and store its inverse in invSy
      invSy(:,:) = 0d0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, Sy)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF( (SCIA_COL_GRID(I,J,DAYOM,1) .GT. 1e15) .and.
     &       (GET_HOUR() .EQ. OBS_HOUR_SCIA_CO(I,J,DAYOM)) .and.
     &        (ERR_PERCENT(I,J) .gt. 0) .and.
     &        (DOMAIN_OBS(I,J) .eq. 1) )THEN
!            invSy(I,J) = 1/((ERR_COL_GRID(I,J,DAYOM)/100.0)**2 * 
            ! 50% error
!            invSy(I,J) = 1/((50.0/100.0)**2 * 
!     &           SCIA_COL_GRID(I,J,DAYOM,1)**2)
            Sy= ERR_PERCENT(I,J)**2 * 
     &           SCIA_COL_GRID(I,J,DAYOM,1)**2
            invSy(I,J) = 1.0/Sy

            IF ( invSy(i,j) .gt. 1 ) THEN
               CALL ERROR_STOP('invSy is too big', 'scia_co_obs_mod.f')
            ENDIF                      
         ELSE
            !DOMAIN_OBS(I,J)=0

         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!      print*, 'min/max of % error:', 
!     &  minval(err_col_grid(:,:,DAYOM)), maxval(err_col_grid(:,:,DAYOM))
!      print*, 'min/max of error:', 
!     &  minval(invSy), maxval(invSy)

!       PRINT*, 'min/max of STT_ADJ, before obs:'
!       PRINT*, minval(STT_ADJ), maxval(STT_ADJ)
!       print*, 'STT_ADJ location of min', minloc(STT_ADJ),'max ',
!     & MAXLOC(STT_ADJ)
!       print*, 'adj stt(6,39,:)',adj_stt(6,39,:,1)

      !DO N = 1, NOBS
      !DO L = 1, NLEV   !for SCIA NLEV=1, since using column data
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J, LL)
!$OMP+PRIVATE( DIFF_COST )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Determine the contribution to the cost function in each grid cell
         ! from each species
         ! CO_COL_GRID is SCIA observations
         IF ( (SCIA_COL_GRID(I,J,DAYOM,1) .GT. 1e15) .and.      
     &       (GET_HOUR() .EQ. OBS_HOUR_SCIA_CO(I,J,DAYOM)) .and. 
     &       (DOMAIN_OBS(I,J) .eq. 1) )THEN

            DIFF_COST = (CHK_STT_SCIA(I,J) - SCIA_COL_GRID(I,J,DAYOM,1))
        
            ! Calculate new additions to cost function
            ! include all regions for which there are obs
            ! NOTE: a bit of a mismatch in weight_obs in vertical
            NEW_COST(I,J)  = DOMAIN_OBS(I,J) * 
            ! Updated for consistency with merged CALC_APRIOR (dkh, 01/18/12, adj32_017) 
            !  (DIFF_COST ** 2) * invSy(I,J) 
     &           0.5d0 * (DIFF_COST ** 2) * invSy(I,J) 

            ! Check for errors
!$OMP CRITICAL
            IF ( IT_IS_NAN( NEW_COST(I,J) ) ) THEN
               WRITE(6,*) ' Bad NEW_COST in ', I, J,
     &              ' from OBS, CHK, DOMAIN_OBS = ', 
!     &              OBS_STT(I,J,L,N), CHK_STT_MOP(I,J,L,N), 
     &              DOMAIN_OBS(I,J), DIFF_COST, invSy(i,j)   
               
               CALL ERROR_STOP('NEW_COST is NaN', 'adjoint_mod.f')
            ENDIF
!$OMP END CRITICAL

            ! update diagnostic arrays if we're saving these diagnostics
            IF ( LDCOSAT ) THEN
               CALL SET_MODEL(I,J,DAY_OF_SIM,2,CHK_STT_SCIA(I,J))
               CALL SET_OBS(I,J,DAY_OF_SIM,2,SCIA_COL_GRID(I,J,DAYOM,1))
               CALL SET_MODEL_BIAS(I,J,DAY_OF_SIM,2, 
     &              DIFF_COST/SCIA_COL_GRID(I,J,DAYOM,1))
               CALL SET_FORCING(I,J,DAY_OF_SIM,NEW_COST(I,J))
               PRINT*, 'MODEL:', I,J,CHK_STT_SCIA(I,J)
               PRINT*, 'SCIA:', SCIA_COL_GRID(I,J,DAYOM,1)

            ! Update cost array
               COST_ARRAY(I,J,DAY_OF_SIM) = COST_ARRAY(I,J,DAY_OF_SIM) + 
     &              NEW_COST(I,J)

            ENDIF

            OBS_COUNT(I,J) = OBS_COUNT(I,J) + 1

            !Adjoint of obs operator, LOOP over all 30 levels
            DO LL=1,LLPAR

               ! Force the adjoint variables x with dJ/dx
               ! NO AVERAGING KERNERLS
!               ADJ_FORCE(I,J,LL,1) = 2.0D0 * DOMAIN_OBS(I,J) 
!     &              * DIFF_COST * invSy(I,J) * 1000.0/28.0*6.02 * 1e23 
!     &              * (1/GET_AREA_CM2(J)) 

               ! WITH AVERAGING KERNERLS
               ! Updated for consistency with merged CALC_APRIOR (dkh, 01/18/12, adj32_017) 
               !ADJ_FORCE(I,J,LL,1) = 2.0D0 * DOMAIN_OBS(I,J) 
               ADJ_FORCE(I,J,LL,1) = DOMAIN_OBS(I,J) 
     &              * DIFF_COST * invSy(I,J)* ADJ_SCIA_ALL(I,J,LL)

               ! Update STT_ADJ, first tracer (CO)
             STT_ADJ(I,J,LL,1) = STT_ADJ(I,J,LL,1) + ADJ_FORCE(I,J,LL,1)

            ENDDO

         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      !ENDDO
      !ENDDO      

      PRINT*, 'OBS this hour:', sum(domain_obs(:,:))
      print*, 'OBS so far:', sum(obs_count)

      ! Error checking: warn of exploding adjoit values, except
      ! the first jump up from zero (MAX_ADJ_TMP = 0 first few times)
      IF ( MAXVAL(ABS(STT_ADJ)) > (MAX_ADJ_TMP * MAX_ALLOWED_INCREASE)
     &   .AND. ( MAX_ADJ_TMP > 0d0 )  ) THEN

         WRITE(6,*)' *** - WARNING: EXPLODING adjoints in ADJ'
         WRITE(6,*)' *** - MAX(STT_ADJ) before = ',MAX_ADJ_TMP
         WRITE(6,*)' *** - MAX(STT_ADJ) after  = ',MAXVAL(ABS(STT_ADJ))

         ADJ_EXPLD_COUNT = ADJ_EXPLD_COUNT + 1

         IF (ADJ_EXPLD_COUNT > MAX_ALLOWED_EXPLD )
     &      CALL ERROR_STOP('Too many exploding adjoints',
     &                       'ADJ_AEROSOL, adjoint_mod.f')

       ENDIF
       
       ! Update cost function
       !PRINT*, 'min/max of ADJ_FORCE:'
       !PRINT*, minval(ADJ_FORCE), maxval(ADJ_FORCE)
       !print*, 'location of min/max of ADJ_FORCE'
       !PRINT*, minloc(ADJ_FORCE), MAXLOC(ADJ_FORCE)
       !print*, 'adj force(6,41,:)',adj_force(6,41,:,1)

       !PRINT*, 'min/max of STT_ADJ:'
       !PRINT*, minval(STT_ADJ), maxval(STT_ADJ)
       !PRINT*, 'min/max of ADJ_EMS:'
       !PRINT*, minval(ADJ_EMS), maxval(ADJ_EMS)
       !print*, 'location of min', minloc(ADJ_EMS),'max of ADJ_EMS',
!     & MAXLOC(ADJ_EMS)
!       print*, 'adj stt(6,41,:)',adj_stt(6,41,:,1)
!       print*, 'adj_stt(9,39,:)',adj_stt(9,39,:,1)

       PRINT*, 'min/max of NEW_COST'
       PRINT*, minval(NEW_COST), maxval(NEW_COST)
       !PRINT*, 'NEW_COST(FD)=', NEW_COST(IFD,JFD,NFD)
       PRINT*, 'TOTAL NEW_COST = ', SUM(NEW_COST)
       PRINT*, 'COST_FUNC BEFORE ADDING NEW_COST=', COST_FUNC
       COST_FUNC = COST_FUNC + SUM ( NEW_COST ) 
       !COST_ARRAY(1,1,1) = COST_ARRAY(1,1,1) + SUM ( NEW_COST ) 

      ! Echo output to screen
       IF ( LPRINTFD ) THEN
          !WRITE(6,*) ' ADJ_FORCE(:) = ', ADJ_FORCE(IFD,JFD,:,NFD)
          WRITE(6,*) ' Using predicted value (CHK_STT_SCIA) = ',
     &         CHK_STT_SCIA(IFD,JFD), '[molec/cm2]' 
          WRITE(6,*) ' Using observed value  (SCIA_STT) = ', 
     &                SCIA_COL_GRID(IFD,JFD,DAYOM,1), '[molec/cm2]' 
          WRITE(6,*) ' Using WEIGHT  = ', DOMAIN_OBS (IFD,JFD) 
          WRITE(6,*) ' ADJ_FORCE = ',
     &         ADJ_FORCE(IFD,JFD,LFD,NFD) 
          WRITE(6,*) ' STT_ADJ = ',
     &         STT_ADJ(IFD,JFD,LFD,NFD)
          WRITE(6,*) ' NEW_COST = ',
     &         NEW_COST(IFD,JFD)
       ENDIF
       
      END SUBROUTINE CALC_SCIAbr_CO_FORCE

!------------------------------------------------------------------------------

      SUBROUTINE CALC_OBS_HOUR

!***************************************************************************
! Subroutine CALC_OBS_HOUR computes an array of hours for each day of obs. 
! If there is an obs in a particular gridbox on that day, it assigns the 
! hour (0..23). If there isn't, OBS_HOUR stays initialized to -1. Also,
! this subroutine computes COUNT_GRID array.
! (mak, 12/14/05)
!***************************************************************************

      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY,
     &                         GET_YEAR, GET_HOUR

#     include "CMN_SIZE"

      INTEGER  :: W, I, J, D
      INTEGER  :: ilon, ilat
      TYPE (XPLEX)   :: OBS_HOURr(IIPAR,JJPAR, NDAYS)
      integer  :: count, as

      count_grid(:,:,:) = 0d0         
      OBS_HOUR_SCIA_CO(:,:,:) = -99
      OBS_HOURr(:,:,:) = 0
      count = 0

      !print*, 'in calc_obs_hour'
      DO W = 1, NOBSS

         !=================================================================
         ! COMPUTE LONGITUDE AND LATITUDE
         !=================================================================
         ! Look for model grid box corresponding to the MOPITT observation:
         ! Get I corresponding to PLON(IND)
         ILON = INT( ( Longitude(W) + 180d0 ) / DISIZE + 1.5d0 )
         ! Handle date line correctly (bmy, 4/23/04)
         IF ( ILON > IIPAR ) ILON = ILON - IIPAR 
         ! Get J corresponding to PLAT(IND)
         ILAT = INT( ( Latitude(W) +  90d0 ) / DJSIZE + 1.5d0 )
         if ( (ilon .eq. -999) .or. (ilat .eq. -999) ) then
            print*,'ilon,ilat=',ilon,ilat
            print*,'STOP'
            stop
         endif

         ! If there's an obs, calculate the time
         IF(SCIACOcol(W) .gt. 0) THEN
 
            D = mday(w)
            
            OBS_HOURr(ILON,ILAT,D) = OBS_HOURr(ILON,ILAT,D)
     &        + time_h(w) 
            count_grid(ILON,ILAT,D) = count_grid(ILON, ILAT,D) + 1

         ENDIF

!          print*, 'obs hour in:', ilon, ilat, 'is:', obs_hour(ilon,ilat)

      ENDDO

      ! average obs_hour on the grid
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, D, count)
      DO D = 1, NDAYS
      DO J = 1, jjPAR
      DO I = 1, IIPAR

         IF ( COUNT_GRID(I,J,D) .gt. 0. ) then
           
            OBS_HOUR_SCIA_CO(I,J,D) = FLOOR((OBS_HOURr(I,J,D)/
     &           COUNT_GRID(I,J,D)))
            count = count + 1
            !IF( D == 2 ) THEN
            !   PRINT*, I,J, OBS_HOUR_SCIA_CO(I,J,D)
            !ENDIF

         ENDIF
            
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !print*, 'obs hour in 144,45 is:', obs_hour(144,45)
      !print*, 'today we have',sum(COUNT_GRID(:,:,GET_DAY()),'obs.'

      END SUBROUTINE CALC_OBS_HOUR
!---------------------------------------------------------------------------

      SUBROUTINE REGRIDV_SCIA( Model_CO_MR )
!***************************************************************************
! Subroutine REGRIDV_SCIA regrids CHK_STT from LLPAR levels of the GC to
! NLevs of SCIA retrieval. This code is a direct Fortran translation 
! of Jenny Fisher's IDL code, which was in turn constructed from
! gamap's regridv.pro It calls a subroutine REGRID_COLUMN, which is a 
! Fortran translation of IDL code, which apparently was a translation of 
! Fortran code that we can no longer locate, which is a shame. (mak, 8/8/07)
!
! NOTES:
! (1 ) Missing from the idl version (Not needed here):
!         ! Airmass on input grid is AD(I,J,L)
!         ! Convert data from [v/v] to mass 
!         ! CHK_STT is already in kg
!         ! Regrid vertically -- preserve column mass (now in kg)
!         !OutCol   = Regrid_Column( InCol, InPEdge, OutPEdge, $!         !                          No_Check=No_Check, _EXTRA=e )
!         ! OutCol is now CHK_STT_SCIA_VGRID

!***************************************************************************

      USE ERROR_MOD, ONLY   : ERROR_STOP
      USE BPCH2_MOD, ONLY   : GET_NAME_EXT,  GET_RES_EXT 
      USE FILE_MOD, ONLY    : IOERROR
      USE TIME_MOD,    ONLY : GET_DAY
      USE PRESSURE_MOD, ONLY: GET_PEDGE
      USE CHECKPT_MOD, ONLY : CHK_PSC
      USE TRACER_MOD,  ONLY : STT, TCVV

#     include "CMN_SIZE" ! PTOP, LLPAR, JJPAR, IIPAR

      ! NLevs = 60 levels of SCIA AKs and pressure levels
      TYPE (XPLEX), INTENT(INOUT):: Model_CO_MR(IIPAR, JJPAR, NLevs)
      TYPE (XPLEX)            :: CHK_STT_SCIA_VGRID(IIPAR,JJPAR,NLevs)
      TYPE (XPLEX)               :: SCIAPress(NLevs)
      TYPE (XPLEX)               :: InPEdge(LLPAR+1)
      TYPE (XPLEX)               :: OutPEdge(NLevs+1)
      TYPE (XPLEX)               :: SCIAEdge(NLevs+1) 
      TYPE (XPLEX)               :: SCIAEdgePressure(NLevs+1)
      TYPE (XPLEX)               :: surfP
      !TYPE (XPLEX)               :: FRACTION(LLPAR,NLevs)
      !TYPE (XPLEX)               :: AirMass(NLevs)

      INTEGER I, J, D, L, LL, IU_FILE, IOS, k
      TYPE (XPLEX)               :: HI, LOW, DIFF
      LOGICAL, SAVE :: FIRST = .TRUE.
      LOGICAL, SAVE :: valid = .FALSE. 

      !================================================================
      ! Read SCIA Pressure info
      !================================================================

      IU_FILE=16
      FRACTION(:,:,:,:) = 0d0

      ! Read SCIA pressures
      OPEN( IU_FILE, FILE='SCIA_pressure.dat', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'sciaP:1' )

      READ( IU_FILE, 100, IOSTAT=IOS) SCIAPress

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'sciaP:2' )
 100  FORMAT(f7.2)


      ! Close file
      CLOSE( IU_FILE )

      ! SURFACE PRESSURE: use checkpointed pressure
      !-------------------
      SCIAEdge(1:60) = SCIAPress(60:1:-1)
      !print*, 'scia pressure, max and min:'
      !print*, maxval(sciaedge(1:60)), minval(sciaedge(1:60))

      !Assume first given edge is 0.01hPa
      SCIAEdge(61) = 0.01

      !Store pressure edges
      SCIAEdgePressure = SciaEDGE

      ! Convert to sigma scale
      surfP=SCIAEdge(1)
      DO k = 1,61 
         SCIAEdge(k)=SCIAEdge(k)/surfP
      ENDDO

      !-------------------
      ! REGRID DATA
      !-------------------

      CHK_STT_SCIA_VGRID(:,:,:) = 0d0

      D = GET_DAY()

      ! Loop over surface grid boxes
      ! WARNING: parallelization screws it up. (mak, 8/15/07)
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, LL, valid, first)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( COUNT_GRID(I,J,D) .gt. 0. ) then

            !to be safe, remove junk values:
            fraction(i,j,:,:) = 0d0                  

            DO L = 1, LLPAR
               ! OutVertEdge = AIRSEdgePressure / PSurf[I,J]
               ! Pressure edges on INPUT and OUTPUT grids
               ! both in and out pressures in hPa
               InPEdge(L)  = ( GET_PEDGE(I,J,L) ) !* Psurf(I,J)) + PTOP)/100
            ENDDO
            InPEdge(LLPAR+1) = PTOP

            ! OutPEdge = ( OutVertEdge * PSurf[I,J] ) ;+ OutType.PTOP
            OutPEdge(:) = SCIAEdgePressure(:) 
            
            !=====================================================
            ! Determine fraction of each INPUT box 
            ! which contributes to each OUTPUT box
            !=====================================================
            ! LM1 = LLPAR, L = L, LM2 = NLevs, K = LL
            ! Loop over INPUT layers
            FIRST = .TRUE.
            valid = .false.
            DO L = 1, LLPAR
       
               ! Reset VALID flag
               Valid = .false.
 
               ! If the thickness of this pressure level is zero, then this 
               ! means that this pressure level lies below the surface 
               ! pressure (due to topography), as set up in the calling
               ! program.  Therefore, skip to the next INPUT level.
               ! This also helps avoid divide by zero errors. (bmy, 8/6/01)
               IF ( ( InPEdge(L) - InPedge(L+1) ) .lt. 1e-5 ) then 
                  !print*, '1. going to 12'
                  goto 12 !NextL
               ENDIF


               ! Loop over OUTPUT layers
               DO LL = 1, NLevs

                 if( OutPEdge(LL) .lt. InPEdge(L) .and.
     &                OutPEdge(LL) .lt. InPEdge(L+1) .and.
     &                 (ll .eq. 1) .and. (l.eq.1) ) THEN
                     Fraction(i,j,L,LL) = 1d0
                     !print*, 'first GC layer lower than 1st SCIA layer'
                     ! Go to next iteration
                     goto 12 !NextL
                  endif

                  !=================================================
                  ! No contribution if:
                  ! -------------------
                  ! Bottom of OUTPUT layer above Top of INPUT layer  OR
                  ! Top    of OUTPUT layer below Bottom of INPUT layer
                  ! ..unless it's the first layer in GC (mak, 8/15/07)
                  !===================================================
                  if ( OutPEdge(LL) .lt. InPEdge(L+1) .OR.
     &                 OutPEdge(LL+1) .gt. InPEdge(L) ) THEN
                     goto 13    !NextLL
                  ENDIF
                      
                  !==================================================
                  ! Contribution if: 
                  ! ----------------
                  ! Entire INPUT layer in OUTPUT layer
                  !===================================================
                  if ( OutPEdge(LL) .ge. InPEdge(L) .AND.
     &                 OutPEdge(LL+1) .le. InPEdge(L+1) ) then 
               
                     Fraction(i,j,L,LL) = 1d0

                     !Indicate a valid contribution from L to LL
                     Valid = .true.
 
                     ! Go to next iteration
                     goto 13 !NextLL
                  endif

                  !==================================================
                  ! Contribution if: 
                  ! ----------------
                  ! Top of OUTPUT layer in INPUT layer
                  !==================================================
                  if ( OutPEdge(LL+1) .le. InPEdge(L)  .AND.
     &                 OutPEdge(LL)  .ge. InPEdge(L) ) THEN 
 
                     Fraction(i,j,L,LL) =(InPEdge(L) - OutPEdge(LL+1)) / 
     &                                ( InPEdge(L) -  InPEdge(L+1) ) 
 
                     ! Indicate a valid contribution from L to LL
                     Valid = .true.
 
                     ! Go to next iteration
                     goto 13 !NextLL
                  endif
            
                  !==================================================
                  ! Contribution if: 
                  ! ----------------
                  ! Entire OUTPUT layer in INPUT layer
                  !==================================================
                  if ( OutPEdge(LL)   .le. InPEdge(L) .AND.
     &                 OutPEdge(LL+1) .ge. InPEdge(L+1) ) then 
 
                     Fraction(i,j,L,LL)=(OutPEdge(LL) - OutPEdge(LL+1))/ 
     &                                 ( InPEdge(L)  -  InPEdge(L+1) )
 
                     ! Also add the to the first OUTPUT layer the fraction
                     ! of the first INPUT layer that is below sigma = 1.0
                     ! This is a condition that can be found in GEOS-3 data.
                     if ( ( First                  )   .AND.  
     &                    ( LL .eq. 1              )   .AND.
     &                    ( InPEdge(L) .gt. OutPEdge(1) ) ) then 
 
                        Fraction(i,j,L,LL) = Fraction(i,j,L,LL) +             
     &                        ( InPEdge(L) - OutPEdge(1)  ) / 
     &                        ( InPEdge(L) - InPEdge(L+1) )                
 
                        ! We only need to do this once...
                        First = .false.
                     endif
 
                     ! Indicate a valid contribution from L to LL
                     Valid = .true.
 
                     ! Go to next iteration
                     goto 13 !NextLL
                  endif
            
                  !===================================================
                  ! Contribution if: 
                  ! ----------------
                  ! Bottom of OUTPUT layer in INPUT layer
                  !===================================================
                  if ( OutPEdge(LL)   .ge. InPEdge(L+1) .AND. 
     &                 OutPEdge(LL+1) .le. InPEdge(L+1) ) then 
                  
                  Fraction(i,j,L,LL) = ( OutPEdge(LL) - InPEdge(L+1) ) /
     &                              ( InPEdge(L)  - InPEdge(L+1) )
                 
                  ! Also add the to the first OUTPUT layer the fraction
                  ! of the first INPUT layer that is below sigma = 1.0
                  ! This is a condition that can be found in GEOS-3 data.
                  if ( ( First         )   .AND.
     &                ( LL .eq. 1      )   .AND. 
     &            ( InPEdge(L) .gt. OutPEdge(1) ) ) then 

                     Fraction(i,j,L,LL) = Fraction(i,j,L,LL) +    
     &                 ( InPEdge(L) - OutPEdge(1)   ) / 
     &                ( InPEdge(L) - InPEdge(L+1) )                
                    
                     ! We only need to do this once...
                     First = .false.
                  endif
              
 
                  ! Indicate a valid contribution from L to LL
                  Valid = .true.
 
                  ! Go to next iteration
                  goto 13       !NextLL
               endif
 
 13            CONTINUE !NextLL

            ENDDO ! LL
 
            !======================================================
            ! Consistency Check:
            ! ------------------
            ! If SUM( FRACTION(L,:) ) does not = 1, there is a problem.
            ! Test those INPUT layers (L) which make a contribution to 
            ! OUTPUT layers (LL) for this criterion.
            !
            !======================================================
            if ( Valid ) then 
               if ( Abs( 1e0 - sum( Fraction(i,j,L,:))) .ge. 1e-4 ) THEN 
                  print*, 'Fraction does not add to 1'
                  print*, L, LL,sum( Fraction(i,j,L,:) )
                  print*, 'frac(5,:):', fraction(i,j,L,:)
                  PRINT*, 'InPEdge:', InPEdge
                  print*, 'OutPEdge:', OutPEdge
    
                  CALL ERROR_STOP ('REGRIDV still sucks', 
     &                 'scia_co_obs_mod.f' )
               endif
            endif
      
 12         CONTINUE  !NextL
         ENDDO !L
      
         !==========================================================
         ! Compute "new" data -- multiply "old" data by fraction of
         ! "old" data residing in the "new" layer
         !==========================================================
         ! Map CO from GC to SCIA grid
         DO LL = 1 , NLevs
            DO L = 1 , LLPAR
               CHK_STT_SCIA_VGRID(I,J,LL) = 
     &              CHK_STT_SCIA_VGRID(I,J,LL)
     &              + STT(I,J,L,1)*FRACTION(i,j,L,LL)
            ENDDO
         ENDDO 

         !print*, 'columns before and after regridding:'
         !PRINT*,I,J,SUM(CHK_STT(I,J,:,1)),SUM(CHK_STT_SCIA_VGRID(I,J,:)) 
         IF(Abs( SUM(CHK_STT_SCIA_VGRID(I,J,:)) - SUM(STT(I,J,:,1)))
     &        /SUM(STT(I,J,:,1)) .gt. 1e-5 ) THEN
            PRINT*, 'columns before and after regrid dont add up:'
            print*, 'columns before and after regridding:'
         PRINT*,I,J,SUM(STT(I,J,:,1)),SUM(CHK_STT_SCIA_VGRID(I,J,:)) 
            PRINT*, 'InPEdge:', InPEdge
            print*, 'OutPEdge:', OutPEdge
            print*, 'chk_stt'
            print*, 'chk_stt_scia_vgrid:'
            print*, chk_stt_scia_vgrid(i,j,:)
            CALL ERROR_STOP ('REGRIDV sucks', 
     &           'scia_co_obs_mod.f' )
         ENDIF

         ! Airmass on output grid (in kg/box in each level)
         AirMass(I,J,:)  = RVR_GetAirMass( SCIAEdge, j, surfP )
         !AirMass  = RVR_GetAirMass( OutVertEdge, OutArea[I,J], surfP )

         ! Convert data from kg to [v/v]
         ! Model_CO_MR = kgCO * gair/gCO / kgair = [v/v]
         DO LL = 1, NLevs
            Model_CO_MR(I,J,LL) = CHK_STT_SCIA_VGRID(I,J,LL) *
     &           TCVV(1)/AirMass(I,J,LL)
         ENDDO
!            DO L = 1, LLPAR
!               DO LL = 1, NLev
!                  ADJ_SCIA_ALL(I,J,L) = ADJ_SCIA_ALL(I,J,L) +
!     &             A(LL)*1e6*ADJ_TCVV(1)/AirMass(I,J,LL)*FRACTION(L,LL)
!                  !ADJ_SCIA_REGRID(I,J,L,LL) = FRACTION(L,LL)
!               ENDDO
!            ENDDO


        ! Model_CO_MR = [CHK_STT(L)*FRACTION(L,LL)]*ADJ_TCVV/Airmass
        ! d(Model_CO_MR)/d(CHK_STT) = (ADJ_TCVV/Airmass)*FRACTION(L,LL)
!         DO LL = 1,NLevs
!            ADJ_SCIA_CONVERT(I,J,LL) =  ADJ_TCVV(1)/AirMass(LL)
!         ENDDO

      ENDIF
      ENDDO
      ENDDO
!!$OMP END PARALLEL DO  


      END SUBROUTINE REGRIDV_SCIA

!---------------------------------------------------------------------------
 
      FUNCTION RVR_GetAirMass( SCIAEdge, J, SurfP) RESULT ( AirMassloc )

      !====================================================================
      ! Internal function RVR_GETAIRMASS returns a column vector of air 
      ! mass given the vertical coordinates, the surface area,
      ! and surface pressure. (bmy, 12/19/03)
      !====================================================================

      USE GRID_MOD,             ONLY : GET_AREA_M2

#     include "CMN_SIZE"

      INTEGER, INTENT(IN)   :: J
      TYPE (XPLEX),  INTENT(IN)   :: SurfP
      TYPE (XPLEX)                :: AirMassloc(NLevs)
      TYPE (XPLEX),  INTENT(IN)   :: SCIAEdge(NLevs+1)
      INTEGER               :: L
      TYPE (XPLEX)                :: g100

      AirMassloc(:) = 0d0

      ! Constant 100/g 
      g100    = 100d0 / 9.8d0 
      
      ! Loop over levels
      ! airmass(L) = hPa * m2 * 1 * 100Pa/hPa * 1/(m/s2) = 
      !            = N * 1/(m/s2) = kg
      DO L = 1, NLevs
         AirMassloc(L) = SurfP * GET_AREA_M2(J) * 
     &        ( SCIAEdge(L) - SCIAEdge(L+1) ) * g100
      ENDDO
      
      END FUNCTION RVR_GetAirMass

!---------------------------------------------------------------------------
      SUBROUTINE INIT_DOMAIN

      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      INTEGER  I, J, as
      LOGICAL, SAVE ::  FIRST = .TRUE.

      IF ( FIRST ) THEN
         ALLOCATE( DOMAIN_OBS( IIPAR,JJPAR ) ,stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'DOMAIN_OBS' )
         FIRST = .FALSE.
      ENDIF

      DOMAIN_OBS(:,:) = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( 
!     &        .and. (I .ge. 93) .and. (I .le. 144) ! MOPITT TRACE-P..
!     &        .and. (J .ge. 41) .and. (J .le.  73) ! ..region 2x2.5
!     &        .and. (J .ge. 40) .and. (J .le.  72) ! ..region 2x2.5
!!     &        .and. (EMS_orig(I,J,NEMS) .NE. 0 )
!!     &        .and. L < LPAUSE(I,J)           ! Only in the troposphere
!!     &        .and.  IS_LAND(I,J)            ! Only the land species
!     &         .and.  ( MOD( I, 2 )  == 0 )   ! Only in every other cell
!!     &        .and.  J >= 10                 ! Not in antarctica
!     &         .and.  L == 8                   ! Only at ~500mb 
!     &         .and. (J .ge. 24)               ! only N.Hemisphere
     &          (J .le. 38)               ! not poleward of 60N
     &         .and. (J .ge. 9)                ! not poleward of 60S
     &                                  ) THEN
        
             DOMAIN_OBS(I,J) = 1
          ELSE
             DOMAIN_OBS(I,J) = 0
          ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      PRINT*, sum(DOMAIN_obs), 'MAX observations today'
     
      END SUBROUTINE INIT_DOMAIN
!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SCIA
      ! DEALLOCATE ALL MEMORY (DONE BEFORE READING REACH MONTHLY FILE)
   
      ! Deallocate
      IF ( ALLOCATED( SCIACOcol     ) ) DEALLOCATE( SCIACOcol     )
      IF ( ALLOCATED( SCIACOcol_err ) ) DEALLOCATE( SCIACOcol_err )
      IF ( ALLOCATED( Longitude     ) ) DEALLOCATE( Longitude     )
      IF ( ALLOCATED( Latitude      ) ) DEALLOCATE( Latitude      )
      IF ( ALLOCATED( SZA           ) ) DEALLOCATE( SZA           )
      IF ( ALLOCATED( SCIA_COL_GRID ) ) DEALLOCATE( SCIA_COL_GRID )
      IF ( ALLOCATED( ERR_COL_GRID  ) ) DEALLOCATE( ERR_COL_GRID  )
      IF ( ALLOCATED( COUNT_GRID    ) ) DEALLOCATE( COUNT_GRID    )
      IF ( ALLOCATED( iday          ) ) DEALLOCATE( iday          )
      IF ( ALLOCATED( mday          ) ) DEALLOCATE( mday          )
      IF ( ALLOCATED( time_h        ) ) DEALLOCATE( time_h        )
      IF ( ALLOCATED( OBS_HOUR_SCIA_CO ) ) DEALLOCATE( OBS_HOUR_SCIA_CO)
      IF ( ALLOCATED( ERR_PERCENT   ) ) DEALLOCATE( ERR_PERCENT   )
      IF ( ALLOCATED( Cloud         ) ) DEALLOCATE( Cloud         )

      END SUBROUTINE CLEANUP_SCIA

!---------------------------------------------------------------------------

      END MODULE SCIAbr_CO_OBS_MOD
