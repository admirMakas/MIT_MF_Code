! $Id: airs_co_obs_mod.f,v 1.3 2012/03/01 22:00:27 daven Exp $
      MODULE AIRS_CO_OBS_MOD

!******************************************************************************
! Module AIRS_CO_OBS_MOD contains subroutines necessary to 
! 1. Transform CHK_STT into AIRS space
! 2. Compute AIRS-GEOS-Chem difference, cost function and adj forcing
! 3. Transform the difference between model and AIRS back to model space
!    using the adjoint of averaging kernel and interpolation code.
!
!  Module Variables:
!  ============================================================================
! (1 ) COUNT_GRID              : number of observations in one GC gridbox
! (2 ) invtest                 : array showing if observation was invertible
! (3 ) ModelPS                 : model surface pressure array
! (2 )
!  Module Routines:
!  ============================================================================
!  (1 ) READ_AIRS_CO_FILES     : Reas AIRS hdf file
!  (2 ) ITS_TIME_FOR_AIRS_CO_OBS: FUNCTION that checks time vs. OBS_HOUR array
!  (3 ) AIRS_FWD         : Driver for fwd obs operator 
!  (4 ) ADJ_AIRS         : Computes the adjoint of observation operator
!  (5 ) READ_ERROR_VARIANCE: Reads error variance file
!  (6 ) CALC_AIRS_CO_FORCE  : Calculates cost function and STT_ADJ increments
!  (7 ) CALC_OBS_HOUR   : Calculated hour of morning obs
!
!  ============================================================================
!  NOTES: 
!  (1 ) Filter on AIRS data: morning over pass only (in airs_mod.f) and
!       only use obs that are greater than 5e17 (mak, 6/07/08)    
!  (2 ) Remove invtest array from the module variables, since it was only 
!       needed when gridding was done separately from computing the column 
!       (mak, 6/12/08)
!  (3 ) Add adjoint code (mak, 6/17/08)
!  (4 ) Update to v8 adjoint (6/20/09)
!
!******************************************************************************

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      ! Everything PRIVATE unless specified otherwise
      ! PRIVATE module variables
      ! PRIVATE module routines
      PRIVATE 

      PUBLIC :: READ_AIRS_CO_FILES, ITS_TIME_FOR_AIRS_CO_OBS
      PUBLIC :: CALC_AIRS_CO_FORCE, OBS_HOUR_AIRS_CO

      TYPE (XPLEX),  ALLOCATABLE :: ERR_PERCENT(:,:)
      INTEGER, ALLOCATABLE :: OBS_HOUR_AIRS_CO(:,:)       
      TYPE (XPLEX),  ALLOCATABLE :: AIRS_COL_GRID(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AIRSDOF_COL_GRID(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: CHK_STT_AIRS(:,:)
      TYPE (XPLEX), ALLOCATABLE :: COUNT_GRID(:,:)
      !INTEGER, ALLOCATABLE :: invtest(:)
      TYPE (XPLEX), ALLOCATABLE  :: ModelPS(:,:)
      TYPE (XPLEX), ALLOCATABLE  :: FRACTION(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE  :: ADJ_AIRS_ALL(:,:,:)

      CONTAINS

      SUBROUTINE READ_AIRS_CO_FILES( YYYYMMDD, HHMMSS )

!******************************************************************************
!  Subroutine READ_AIRS_CO_FILES reads the AIRS hdf file and assigns OBS_HOUR
!  array based on available data. AIRS data are stored in a 1 day/file
!  frequency. (mak, 7/12/07, 6/08/08)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!

      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE AIRSv5_MOD
      USE TIME_MOD,  ONLY : EXPAND_DATE
      USE FILE_MOD,  ONLY : IOERROR


#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      CHARACTER(LEN=255)  :: DIR_AIRS
      CHARACTER(LEN=255)  :: FILENAME_IN      
      CHARACTER(LEN=255)  :: file
      CHARACTER(LEN=  8)  :: YYYYMMDDs
      INTEGER             :: IU_FILE, IOS, IOS1, I, as

      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
      LOGICAL, SAVE       :: FIRST = .TRUE. 
      
      !=================================================================
      ! READ_AIRS_CO_FILES begins here!
      !=================================================================
      
      CALL CLEANUP_AIRS

      ! Set date and corresponding input AIRS filename
      DIR_AIRS = '/lustre/data/obs/airs/YYYY/MM/YYYYMMDD/'
      !DIR_AIRS = '/san/as04/home/ctm/mak/AIRS/data_airs/'
      !DIR_AIRS = '/as/data-rw/corrections/as/data/airs/'
      FILENAME_IN='input.txt'

      CALL EXPAND_DATE( DIR_AIRS, YYYYMMDD, 0 )

      print*, 'dir_airs is:', trim(dir_airs)

      ! mak debug: use just 20 files
      CALL SYSTEM('ls '//trim(DIR_AIRS)//' > input.txt')

      IU_FILE=15
      OPEN( IU_FILE, FILE=FILENAME_IN, IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'airs:1' )

      ! zero counters
      Nfiles = 0

      ! Figure out how many files to read (#lines in the file):
      CALL SYSTEM('wc -l '//trim(FILENAME_IN)//' > tmp.txt')
      
      OPEN( 5, FILE='tmp.txt', IOSTAT=IOS1 )
      IF ( IOS1 /= 0 ) CALL IOERROR( IOS1, 5, 'tmp:1' )
      
      ! Read #lines
      READ( 5, *, IOSTAT=IOS1  ) NFiles
      IF ( IOS1 /= 0 ) CALL IOERROR( IOS1, 5, 'tmp:2' )
   
      ! Close file
      CLOSE( 5 )

      ALLOCATE( iNObs (NFiles) )

      ALLOCATE( FILENAME(NFiles) )
      DO i = 1, NFiles
     
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'airs:2' )
         
         READ( IU_FILE,'(a)',IOSTAT=IOS) file

         IF (i .eq. 1) 
     &        YYYYMMDDs = file(6:9)//file(11:12)//file(14:15)

         ! on ceres: 
         FILENAME(i)=trim(DIR_AIRS)//trim(file)
         ! on prometheus and tethys:
         !FILENAME(i)=trim(file)
         print*, 'filename:', trim(filename(i))

      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      !READ (YYYYMMDDs,*) YYYYMMDD
      
      PRINT*,'Date: ',YYYYMMDD

      IF(FIRST) THEN
         ALLOCATE( OBS_HOUR_AIRS_CO( IIPAR, JJPAR ), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'OBS_HOUR_AIRS_CO' )
         ALLOCATE( ERR_PERCENT( IIPAR, JJPAR ), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'ERR_PERCENT' )
         !ALLOCATE( ADJ_FACTOR( IIPAR, JJPAR, LLPAR), stat=as )
         !IF ( as /= 0 ) CALL ALLOC_ERR( 'ADJ_FACTOR' )    
         ALLOCATE( AIRS_COL_GRID( IIPAR, JJPAR ), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'AIRS_COL_GRID' )
         ALLOCATE( AIRSDOF_COL_GRID( IIPAR, JJPAR ), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'AIRSDOF_COL_GRID' )
         ALLOCATE( CHK_STT_AIRS( IIPAR, JJPAR ), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'CHK_STT_AIRS' )
         ALLOCATE ( COUNT_GRID(IIPAR, JJPAR), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'COUNT_GRID' ) 
         ALLOCATE( ModelPS(IIPAR,JJPAR) ,stat = as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'ModelPs' )
         ALLOCATE( FRACTION(IIPAR,JJPAR,LLPAR,NLevs), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'FRACTION' ) 
         ALLOCATE( ADJ_AIRS_ALL(IIPAR,JJPAR,LLPAR), stat=as )
         IF ( as /= 0 ) CALL ALLOC_ERR( 'ADJ_AIRS_ALL' )
      
         FIRST = .FALSE.
      ENDIF

      ! initialize some arrays, the rest is initialized before 
      ! relevant calculations, every hour when we have obs

      ! Get dimensions of the arrays, allocate arrays, and read AIRS file
      CALL INIT_READ_AIRS

      ! Calculate hour of day when obs should be compared to model
      CALL CALC_OBS_HOUR


      END SUBROUTINE READ_AIRS_CO_FILES

!--------------------------------------------------------------------------

      SUBROUTINE CALC_OBS_HOUR

!***************************************************************************
! Subroutine CALC_OBS_HOUR computes an array of hours for each day of obs. 
! If there is an obs in a particular gridbox on that day, it assigns the 
! hour (0..23). If there isn't, OBS_HOUR stays initialized to -1. 
! (mak, 12/14/05, 6/10/08)
!***************************************************************************

!      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY,
     &                         GET_YEAR, GET_HOUR
      USE AIRSV5_MOD 

#     include "CMN_SIZE"

      TYPE (XPLEX)   :: tau0
      INTEGER  :: W, I, J
      INTEGER  :: ilon, ilat
      TYPE (XPLEX)   :: OBS_HOURr(IIPAR,JJPAR)
      integer  :: count
      TYPE (XPLEX)   :: AirsGMT, lon15
      
      ! Get TAU0 from the date (at 0GMT)
      tau0 = GET_TAU0(GET_MONTH(), GET_DAY(), GET_YEAR())

      OBS_HOUR_AIRS_CO(:,:) = -1
      OBS_HOURr(:,:) = 0
      COUNT_GRID(:,:) = 0d0
      count = 0

      DO W = 1, Nobss

         !============================================================
         !Only consider day time AIRS measurements
         ! am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)
         !==============================================================
         IF ( (qual(W) .eq. 0) .AND. (NTraps(W) .gt. 1) .AND. 
     &         (DNFlag(W) .eq. 'Day') .AND. 
     &         (Tsurf(W) .ge. 250) ) THEN
 

         ! Compute local time: 
         ! Local Time = GMT + ( longitude / 15 ) since each hour of time
         ! corresponds to 15 degrees of longitude on the globe
         ! so 
         ! GMT = local time - (longitude/15)
         !============================================================
         lon15 = longitude(w)/15.
         AirsGMT = hour(w)- lon15
         if (AirsGMT .lt. 0.)  AirsGMT = AirsGMT + 24
         if (AirsGMT .gt. 24.)  AirsGMT = AirsGMT - 24

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
         IF ( (COcol(W) .gt. 0.) .and. 
     &        (qual(W) .eq. 0) .AND. (NTraps(W) .gt. 1) .AND. 
     &        (DNFlag(W) .eq. 'Day') .AND. 
     &        (Tsurf(W) .ge. 250) )THEN

            COUNT_GRID(ILON,ILAT) = COUNT_GRID(ILON,ILAT) + 1.
            !Add the time of obs, to be averaged and floored later
            OBS_HOURr(ILON,ILAT) = OBS_HOURr(ILON,ILAT)
     &           + AirsGMT
!          print*, 'obs hour in:', ilon, ilat, 'is:', obs_hour(ilon,ilat)
         ENDIF
      ENDIF !morning overpass
      ENDDO

      ! average obs_hour on the grid
      DO J = 1, jjPAR
      DO I = 1, IIPAR
         IF ( COUNT_GRID(I,J) .gt. 0. ) then
           OBS_HOUR_AIRS_CO(I,J) = FLOOR(OBS_HOURr(I,J)/COUNT_GRID(I,J))

               count = count + 1
         ENDIF
      ENDDO
      ENDDO

      print*, 'today we have (globally)',count,'AIRS observations.'

      END SUBROUTINE CALC_OBS_HOUR

!-------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_AIRS_CO_OBS( ) RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_AIRS_CO_OBS returns TRUE if there are observations
!  available for particular time (hour of a particular day) based on
!  the OBS_HOUR_AIRS_CO array which holds the hour of obs in each gridbox 
!  (computed when file read in airsv5_mod.f90) (mak, 6/09/08)
!
!  NOTES: 
!  ( 1) Also like corresponding MOPITT code
!
!******************************************************************************
!

      USE TIME_MOD,   ONLY : GET_HOUR, GET_MINUTE
      
#     include "CMN_SIZE"  ! Size params

      ! Function value
      LOGICAL :: FLAG

      INTEGER :: I,J

      !=================================================================
      ! ITS_TIME_FOR_AIRS_CO_OBS begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.

      DO J = 1,JJPAR
         DO I = 1,IIPAR
            IF(GET_HOUR() == OBS_HOUR_AIRS_CO(I,J) 
     &   .AND. GET_MINUTE() == 0) THEN
                !print*, 'obs_hour was', get_hour(), 'in box', i, j
               FLAG = .TRUE.
               GOTO 11
            ENDIF
         ENDDO
      ENDDO 

 11   CONTINUE
      END FUNCTION ITS_TIME_FOR_AIRS_CO_OBS

!---------------------------------------------------------------------------

      SUBROUTINE CALC_AIRS_CO_FORCE

      ! References to F90 modules
      USE ERROR_MOD, 		ONLY : IT_IS_NAN, ERROR_STOP
      USE AIRSV5_MOD,           ONLY : DOMAIN_OBS
      USE TIME_MOD,             ONLY : GET_HOUR, GET_NYMDe, GET_NHMSe,
     &                                 GET_MONTH
      USE ADJ_ARRAYS_MOD,       ONLY : SET_FORCING, SET_MOP_MOD_DIFF, 
     &                          SET_MODEL_BIAS, SET_MODEL, SET_OBS,
     &                          GET_FORCING, COST_ARRAY, OBS_COUNT,
     &                          SET_DOFS, IFD, JFD, LFD, NFD, COST_FUNC,
     &                          NOBS, DAY_OF_SIM, ADJ_FORCE, STT_ADJ
      USE TRACER_MOD,           ONLY : N_TRACERS
      USE LOGICAL_ADJ_MOD,      ONLY : LPRINTFD


#     include "CMN_SIZE" 	! Size parameters

      ! Internal variables 
      TYPE (XPLEX)  :: DIFF
      TYPE (XPLEX)  :: DIFF_COST
      TYPE (XPLEX)  :: DIFF_ADJ(LLPAR)
      TYPE (XPLEX)  :: NEW_COST(IIPAR,JJPAR,NOBS) !column cost 
      INTEGER :: I, J, L, N, LL
      INTEGER :: ADJ_EXPLD_COUNT
      INTEGER, PARAMETER  ::  MAX_ALLOWED_EXPLD    = 10
      TYPE (XPLEX),  PARAMETER  ::  MAX_ALLOWED_INCREASE = 10D15
      TYPE (XPLEX)  :: MAX_ADJ_TMP
      TYPE (XPLEX)  :: invSy(IIPAR,JJPAR) !error variance for column
      LOGICAL, SAVE ::  FIRST= .TRUE.
      TYPE (XPLEX)  :: Sy
      INTEGER, SAVE :: LASTMONTH = -999

      !================================================================ 
      ! CALC_MOPITT_FORCE begins here!
      !================================================================

      !initialize:
      CHK_STT_AIRS(:,:) = 0d0
      AIRS_COL_GRID(:,:) = 0d0
      AIRSDOF_COL_GRID(:,:) = 0d0
      invSy(:,:) = 0d0
      NEW_COST(:,:,:) = 0d0
      
      ! column AIRS data is in COTotalColumn
      CALL AIRS_COMPUTE_COLUMN

      ! Read in the matrix with mean % variance to compute
      ! 1/(ymod*%)^2=invSy
      IF(GET_MONTH() .ne. LASTMONTH)THEN
         ! using seasonal errors, but read file once month for simplicity
         ! better than reading the file every time step (mak, 1/27/08)
         print*, 'read error variance matrix'
         CALL READ_ERROR_VARIANCE
         LASTMONTH = GET_MONTH()
      ENDIF

      print*, 'max AIRS value is:', maxval(AIRS_COL_GRID)
      print*, 'min AIRS value is:', minval(AIRS_COL_GRID)
      print*, 'max model value is:',maxval(CHK_STT_AIRS)
      print*, 'min model value is:',minval(CHK_STT_AIRS)
      !print*, 'max err value is:', maxval(ERR_PERCENT(:,:))
      !print*, 'min err value is:', minval(ERR_PERCENT)
      
      ! CHK_STT_AIRS in molec/cm2, OBS_STT in molec/cm2
      !print*, 'before loop: domain_obs is', sum(domain_obs)/30
      !print*, 'before loop, count is:', count

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, Sy)
      DO J= 1,JJPAR
      DO I= 1,IIPAR

         IF ((AIRS_COL_GRID(I,J) .gt. 0) .AND. 
     &        (OBS_HOUR_AIRS_CO(I,J)  .eq. GET_HOUR()).AND.
     &        (DOMAIN_OBS(I,J) .eq. 1) ) THEN

            Sy = ERR_PERCENT(I,J)**2 * 
     &           AIRS_COL_GRID(I,J)**2
            invSy(I,J) = 1/Sy

            OBS_COUNT(I,J) = OBS_COUNT(I,J) + 1

            IF ( invSy(i,j) .ge. 1 ) THEN
               CALL ERROR_STOP('invSy is too big', 'airsitt_obs_mod.f')
            ENDIF                      
         ELSE

            DOMAIN_OBS(I,J)=0
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      print*,'AIRS obs used this hour', sum(domain_obs)/30
      PRINT*, 'OBS_COUNT TOTAL:', SUM(OBS_COUNT)
      print*, 'min/max of invSy:', minval(invSy), maxval(invSy)

      DO N = 1, NOBS 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,    J)
!!$OMP+PRIVATE( DIFF_COST )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( (AIRS_COL_GRID(I,J) .GT. 1e15) .and.      
     &       (GET_HOUR() .EQ. OBS_HOUR_AIRS_CO(I,J)) .and. 
     &       (DOMAIN_OBS(I,J) .eq. 1))then! .and.
!     &       (CHK_STT_AIRS(I,J) .GT. 0) )THEN


         ! Determine the contribution to the cost function in each grid cell
         ! from each species
         DIFF_COST  = ( CHK_STT_AIRS(I,J) - AIRS_COL_GRID(I,J) )
        
         ! Calculate new additions to cost function
         ! include all regions for which there are obs
         ! NOTE: a bit of a mismatch in domain_obs in vertical
         NEW_COST(I,J,N)  = DOMAIN_OBS(I,J) * 
         ! Update to be consistent with merged APCOST routine (dkh, 01/18/12, adj32_017) 
         !     (DIFF_COST ** 2) * invSy(I,J) 
     &        0.5d0 * (DIFF_COST ** 2) * invSy(I,J) 

         ! Diagnostic stuff: FORCING, MOP_MOD_DIFF, MODEL_BIAS
!         CALL SET_FORCING(I,J,DAY_OF_SIM,NEW_COST(I,J,N))
!         FORCING(I,J,DAY_OQF_SIM) = FORCING(I,J,DAY_OF_SIM) 
!     &        + NEW_COST(I,J,L,N)

         if((DOMAIN_OBS(I,J) .eq. 1) )then
      CALL SET_MODEL_BIAS(I,J,DAY_OF_SIM,3,DIFF_COST/AIRS_COL_GRID(I,J))
         CALL SET_MODEL(I,J,DAY_OF_SIM,3,CHK_STT_AIRS(I,J))
         CALL SET_OBS(I,J,DAY_OF_SIM,3, AIRS_COL_GRID(I,J))
         CALL SET_DOFS(I,J,DAY_OF_SIM,3, AIRSDOF_COL_GRID(I,J))
         endif

         ! update cost array
         COST_ARRAY(I,J,DAY_OF_SIM) = COST_ARRAY(I,J,DAY_OF_SIM) + 
     &         NEW_COST(I,J,1)

         ! Check for errors
!!$OMP CRITICAL
         IF ( IT_IS_NAN( NEW_COST(I,J,N) ) ) THEN
            WRITE(6,*) ' Bad NEW_COST in ', I, J, L, N,
     &                 ' from OBS, CHK, DOMAIN_OBS = ', 
     &                 AIRS_COL_GRID(I,J), CHK_STT_AIRS(I,J), 
     &                 DOMAIN_OBS(I,J), DIFF_COST, invSy(i,j)   

            CALL ERROR_STOP('NEW_COST is NaN', 'adjoint_mod.f')
         ENDIF
!!$OMP END CRITICAL

         !LOOP over all 30 levels
         DO LL=1,LLPAR

            ! Force the adjoint variables x with dJ/dx
            ! Update to be consistent with merged APCOST routine (dkh, 01/18/12, adj32_017) 
            !ADJ_FORCE(I,J,LL,N) = 2.0D0 * DOMAIN_OBS(I,J) 
            ADJ_FORCE(I,J,LL,N) = DOMAIN_OBS(I,J) 
     &           * DIFF_COST * invSy(I,J) * ADJ_AIRS_ALL(I,J,LL)
!     &           * 1.0 * 1.0

            ! Update STT_ADJ 
            IF ( N <= N_TRACERS ) THEN 
             STT_ADJ(I,J,LL,N) = STT_ADJ(I,J,LL,N) + ADJ_FORCE(I,J,LL,N)
             !PRINT*, 'ADJ_FORCE,I,J,L:', I,J,LL,ADJ_FORCE(I,J,LL,N) 
             !print*, 'ADJ_AIRS_ALL:', ADJ_AIRS_ALL(I,J,LL)
            ENDIF
         ENDDO

         IF(I == IFD .AND. J == JFD) THEN
            !PRINT*, 'CHK_STT:', CHK_STT(I,J,:,N)
            PRINT*, 'N = ', N
            PRINT*, 'CHK_STT_AIRS:', CHK_STT_AIRS(I,J)
            PRINT*, 'OBS_STT:', AIRS_COL_GRID(I,J)
            PRINT*, 'NEW_COST', NEW_COST(I,J,N)
            PRINT*, 'DIFF_COST', DIFF_COST
            PRINT*, 'DOMAIN_OBS', DOMAIN_OBS(I,J) 
            PRINT*, 'ADJ_FORCE:', ADJ_FORCE(I,J,:,N)
            PRINT*, 'STT_ADJ:', STT_ADJ(I,J,:,N)
         ENDIF


         ENDIF
      ENDDO
      ENDDO
!!$OMP END PARALLEL DO
      ENDDO

      !have to zero the NEW_COST that is above 7th layer
      !NEW_COST(:,:,NLEV+1:LLPAR,:)=0d0

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
       
       ! Update cost array, uncomment if L=1,LFDSIZE
c$$$       DO L = 1, LFDSIZE
c$$$       DO J = 1, JJPAR
c$$$       DO I = 1, IIPAR
c$$$          ! COST_ARRAY
c$$$!          COST_ARRAY(I,J,DAY_OF_SIM) = COST_ARRAY(I,J,DAY_OF_SIM) + 
c$$$!     &         NEW_COST(I,J,1,1)
c$$$          COST_ARRAY(I,J,L) = COST_ARRAY(I,J,L) + 
c$$$     &         NEW_COST(I,J,L,1)
c$$$       ENDDO
c$$$       ENDDO
c$$$       ENDDO


       ! Update cost function
       !PRINT*, 'NEW_COST(FD)=', NEW_COST(IFD,JFD,LFD,NFD)
       PRINT*, 'TOTAL NEW_COST = ', SUM(NEW_COST)
       PRINT*, 'COST_FUNC BEFORE ADDING NEW_COST=', COST_FUNC
       COST_FUNC = COST_FUNC + SUM ( NEW_COST ) 

      ! Echo output to screen
      IF ( LPRINTFD ) THEN
         WRITE(6,*) ' ADJ_FORCE(:) = ', ADJ_FORCE(IFD,JFD,:,NFD)
         WRITE(6,*) ' Using predicted value (CHK_STT_AIRS) = '
     &              , CHK_STT_AIRS(IFD,JFD), '[molec/cm2]' 
         WRITE(6,*) ' Using observed value  (OBS_STT) = '
     &              , AIRS_COL_GRID(IFD,JFD), '[molec/cm2]' 
         WRITE(6,*) ' Using WEIGHT  = ', DOMAIN_OBS (IFD,JFD) 
         WRITE(6,*) ' ADJ_FORCE = '
     &              , ADJ_FORCE(IFD,JFD,LFD,NFD), '[1/molec/cm2]' 
         WRITE(6,*) ' STT_ADJ = '
     &              , STT_ADJ(IFD,JFD,LFD,NFD), '[1/molec/cm2]' 
         WRITE(6,*) ' NEW_COST = '
     &              , NEW_COST(IFD,JFD,NFD)
      ENDIF

      !PRINT*, 'END CALC_AIRS_CO_FORCE'
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )     

      END SUBROUTINE CALC_AIRS_CO_FORCE

!---------------------------------------------------------------------------

      FUNCTION APRIORI_COMP( H2Ocd_loc, BotLev_loc) RESULT ( Xaloc )

      !====================================================================
      ! Function APRIORI_COMP computes the a priori layer column density,
      ! based on the a priori mixing ratio profile, for a given observation,
      ! given the layer column density of water for the observation.
      !
      ! This function is based on the subroutine coch4fg.f written by Evan
      ! Manning at NASA-JPL/Caltech for AIRS Level 2 Data Processing. It
      ! has been adapted for our use at Harvard. (jaf, 11/04/07)
      !====================================================================

      USE AIRSv5_MOD

      TYPE (XPLEX), INTENT(IN)    :: H2Ocd_loc(NLevs)
      INTEGER, INTENT(IN)   :: BotLev_loc
      TYPE (XPLEX)                :: Xaloc(NLevs)
      TYPE (XPLEX)                :: COmr(NLevs)
      REAL                  :: LCDdry, c0, eps, delP
      INTEGER               :: L

      TYPE (XPLEX), PARAMETER     :: AVOGAD = 6.02214199D23 !molec/mol
      TYPE (XPLEX), PARAMETER     :: WATMOL_SI=18.0152D-3   !kg/mol
      TYPE (XPLEX), PARAMETER     :: GRAV_SI=9.80665        !m/s^2
      TYPE (XPLEX), PARAMETER     :: MDRYAIR_SI = 28.964D-3 !kg/mol

      ! This a priori is the MOPITT a priori from Eric Maddy
      ! with AFGL profile replacing the top 21 layers
      COmr = (/1751.373, 652.180, 313.65, 198.207, 132.753,  
     &      93.20174,  69.14629,  55.71019,  45.23621,  37.38374, 
     &      32.96809,  29.58296,  26.98655,  24.92291,  23.32885, 
     &      21.93244,  20.72882,  19.77256, 19.2927,  18.8556,    
     &      18.4784,  18.1747,  17.9416,  17.7607,  17.5914,      
     &      17.3647,  17.0681,  16.8321,  16.8665,  17.4449,      
     &      18.6399,  20.3268,  22.3094,  24.4428,  26.5907,      
     &      28.6583,  30.6771,  32.7125,  34.8565,  37.1747,      
     &      39.6367,  42.1871,  44.7521,  47.2370,  49.5845,      
     &      51.8020,  53.9116,  55.9452,  57.9489,  59.9825,      
     &      62.0931,  64.2530,  66.4136,  68.5153,  70.4870,      
     &      72.2682,  73.8662,  75.3098,  76.6364,  77.8896,      
     &      79.0853,  80.2262,  81.3165,  82.3502,  83.2717,      
     &      84.0038,  84.4653,  84.6851,  84.8021,  84.9779,      
     &      85.2669,  85.5724,  85.7798,  85.8779,  85.9660,      
     &      86.1554,  86.5235,  87.1354,  87.9981,  88.9903,      
     &      89.9855,  91.0270,  92.2368,  93.5934,  94.9150,      
     &      96.2908,  98.3620, 101.2763, 103.5686, 104.2298,      
     &      104.8466, 106.3291, 107.6906, 110.6079, 114.1500,      
     &      118.1522, 120.4714, 120.6016, 120.6461, 120.6100/)

      Xaloc(:)=0d0

      eps = WATMOL_SI/MDRYAIR_SI
      c0 = 1.e-2*AVOGAD/(MDRYAIR_SI*GRAV_SI)
      DO L = 1, NLevs
         IF (L .eq. 1) THEN
            delP = Plev(L) - 0.005
         ELSE
            delP = Plev(L) - Plev(L-1)
         ENDIF
         LCDdry = c0*delP - eps*H2Ocd_loc(L)
         IF (L .le. BotLev_loc) THEN
            Xaloc(L) = LCDdry*comr(L)*1.E-9
         ELSE
            Xaloc(L) = Xaloc(BotLev_loc)
         ENDIF
      ENDDO

      END FUNCTION APRIORI_COMP

!---------------------------------------------------------------------------

      FUNCTION FMATRIX_COMP( PSurf, NTloc) RESULT( Fmatloc)

!*********************************************************************
! FMATRIX_COMP CALCULATES THE F MATRIX NEEDED TO COMPARE GC TO AIRS
!*********************************************************************    
! Note: This function is translated from the MATLAB script
!       Fmatrix_comp written by Wallace McMillan (8/28/07)
!*********************************************************************    

      USE AIRSv5_MOD

      TYPE (XPLEX), INTENT(IN)   :: Psurf
      INTEGER, INTENT(IN)  :: NTloc
      TYPE (XPLEX)               :: Fmatloc(NTloc,NLevs)
      INTEGER              :: ind0(3),ind_end(3),ind(NTloc-2,4)
      INTEGER              :: I, J, P, DD
      TYPE (XPLEX)               :: PiBot

      !=====================================================
      ! FMATRIX_COMP begins here
      !=====================================================

      ! First, define the sets of indice vectors used to setup a
      ! vector for each trapezoid. The first and last indice vectors
      ! contain only three entries because the trapezoids go to 0.5
      ! at the top and bottom of the atmosphere.

      ind0 = COlev(1:3)
      DO I = 1, NTloc-2
         ind(i,:)=COlev(i:i+3)
      ENDDO
      ind_end=COlev(NTloc-1:NTloc+1)

      ! Construct f vectors, one for each trapezoid, and fill with zeros.
      ! f vectors for ALL possible trapezoids are initially setup this way.
      ! Only the f vectors corresponding to the actual trapezoids used in
      ! the retrieval will be assigned non-zero values, below.
      Fmatloc(:,:)=0d0

      ! Now, compute the value of each trapezoid used in the retrieval at
      ! each of the 100 AIRS levels.

      ! First, the top trapezoid must be defined from only the first three
      ! entries in COlev because its value at the top of the atmosphere=0.5
      Fmatloc( 1, ind0(1):ind0(2) ) = 0.5
      Fmatloc( 1, ind0(2):ind0(3) ) = 0.5 * 
     &     (1-(log(Plev(ind0(2):ind0(3)))-log(Plev(ind0(2)))) / 
     &     (log(Plev(ind0(3)))-log(Plev(ind0(2)))))

      ! Next, compute the f vectors for the middle trapezoids.  The first
      ! trapezoid that encounters the surface must be scaled to terminate 
      ! at the pressure level nearest to and above the surface.

      DO I = 1, NTloc-2
         PiBot=Plev(ind(i,4))
         IF ( PiBot .le. Psurf) THEN
            Fmatloc(i+1,ind(i,1):ind(i,2)) = 0.5 * 
     &             (1-(log(Plev(ind(i,1):ind(i,2))) - 
     &           log(Plev(ind(i,2)))) / 
     &            (log(Plev(ind(i,1))) - log(Plev(ind(i,2)))))
            Fmatloc(i+1,ind(i,2):ind(i,3)) = 0.5
            Fmatloc(i+1,ind(i,3):ind(i,4)) = 0.5 * 
     &           (1-(log(Plev(ind(i,3):ind(i,4))) - 
     &           log(Plev(ind(i,3)))) / 
     &           (log(Plev(ind(i,4))) - log(Plev(ind(i,3)))))
         ELSE
            ! J is the first index for which Plev is gt Psurf
            J = 0
            DO P = 1, NLevs
               IF ( (J .eq. 0) .AND. (PLev(P) .gt. Psurf) ) J=P
            ENDDO
            Fmatloc(i+1,ind(i,1):ind(i,2)) = 0.5 * 
     &           (1-(log(Plev(ind(i,1):ind(i,2))) - 
     &           log(Plev(ind(i,2)))) / 
     &           (log(Plev(ind(i,1))) - log(Plev(ind(i,2)))))
            Fmatloc(i+1,ind(i,2):ind(i,3)) = 0.5

            ! If the bottom trapezoid has a face that is less than one AIRS
            ! layer wide, then we must make sure the next to the bottom
            ! trapezoid goes to zero at the lowest level above the surface.
            DD = (J-1 - ind(i,3))
            IF (DD .ne. 0)  THEN
               Fmatloc(i+1,ind(i,3):J-1) = 0.5 * 
     &              (1-(log(Plev(ind(i,3):J-1)) - 
     &              log(Plev(ind(i,3)))) / 
     &              (log(Plev(J-1)) - log(Plev(ind(i,3)))))
            ELSE
               Fmatloc(i+1,ind(i,3):J-1) = 0.0
            ENDIF
         ENDIF
      ENDDO

      ! Then, compute the bottom trapezoid.  If it encounters the surface, 
      ! only set it equal to 0.5 down to the lowest level above the surface.
      ! Below the surface, it will remain at zero.  All trapezoids below the
      ! surface (i.e. not used in the retrieval) will keep zero values.
      PiBot=Plev(ind_end(3))
      IF ( PiBot .gt. Psurf) THEN
         ! J is the first index for which Plev is gt Psurf
         J = 0
         DO P = 1, NLevs
            IF ( (J .eq. 0) .AND. (PLev(P) .gt. Psurf) ) J=P
         ENDDO
         Fmatloc(NTloc,ind_end(1):ind_end(2)) = 0.5 * 
     &        (1-(log(Plev(ind_end(1):ind_end(2))) - 
     &        log(Plev(ind_end(2)))) / 
     &        (log(Plev(ind_end(1))) - log(Plev(ind_end(2)))))
         Fmatloc(NTloc,ind_end(2):J-1) = 0.5
      ELSEIF (PiBot .eq. Psurf) THEN
         J = 0
         DO P = 1, NLevs
            IF ( (J .eq. 0) .AND. (PLev(P) .eq. Psurf) ) J=P
         ENDDO
         Fmatloc(NTloc,ind_end(1):ind_end(2)) = 0.5 * 
     &        (1-(log(Plev(ind_end(1):ind_end(2))) - 
     &        log(Plev(ind_end(2)))) / 
     &        (log(Plev(ind_end(1))) - log(Plev(ind_end(2)))))
         Fmatloc(NTloc,ind_end(2):J) = 0.5
      ENDIF

      END FUNCTION FMATRIX_COMP

!------------------------------------------------------------------------------

      FUNCTION CALCtotcolden( laycolden,psurf) RESULT( totcolden)

!*********************************************************************
! CALCtotcolden CALCULATES THE TOTAL COLUMN DENSITY
!*********************************************************************    
! Note: This function is translated from the MATLAB script
!       CALCtotcolden written by Wallace McMillan (11/20/07)
!*********************************************************************    

      !This routine calculates total column density for a given input
      !layer column density profile with specified pressures for the
      !layer boundaries and an input surface pressure.  For 100 layer
      !column densities, there must be 101 pressure boundaries (top and
      !bottom of each layer).  In the case of the surface pressure lying
      !in the middle of a layer (between two pressure boundaries), only
      !a fraction of this lowest layer is used in computing the total column
      !density.  This fraction is computed in a dlogP since to account for
      !the general exponential variation of pressure with altitude.
      !W. McMillan, 11/12/03

      !Correction to if statement pbound(il-1) changed to pbound(il)
      ! W. McMillan, 8/22/06

      USE AIRSv5_MOD

      TYPE (XPLEX), INTENT(IN)   :: psurf,laycolden(NLevs)
      TYPE (XPLEX)               :: pbound(NLevs),totcolden
      TYPE (XPLEX)               :: frac
      INTEGER              :: I,IL

      pbound=Plev
      totcolden=0d0
      frac=0d0
      IL=0

      findil:      DO I=1,NLevs
         IF (pbound(I) .gt. psurf) THEN
            EXIT findil
         ENDIF
         IL=I
      ENDDO findil
      
      IF (pbound(IL) .eq. psurf) THEN
         totcolden = SUM(laycolden(1:IL+1-1))
      ELSE
         frac=(LOG(psurf)-LOG(pbound(IL)))/(LOG(pbound(IL+1))-
     &        LOG(pbound(IL)))
         totcolden = SUM(laycolden(1:IL+1-1)) + laycolden(IL+1)*frac
      ENDIF

      END FUNCTION CALCtotcolden

!------------------------------------------------------------------------------

      SUBROUTINE AIRS_COMPUTE_COLUMN

!*********************************************************************
! COMPUTES COLUMN FROM GEOS-CHEM OUTPUT WITH AIRS AVERAGING KERNELS 
!*********************************************************************    

      USE AIRSv5_MOD
      USE TIME_MOD,   ONLY : GET_HOUR
      USE GRID_MOD,             ONLY : GET_AREA_M2

# include "CMN_SIZE"

      INTEGER             :: L, W, J, I, LL
      INTEGER             :: ILON, ILAT
      INTEGER             :: NTloc
      TYPE (XPLEX)              :: Psurf
      TYPE (XPLEX), ALLOCATABLE :: Ftrans(:,:), Fpi(:,:)
      TYPE (XPLEX), ALLOCATABLE :: Fmat(:,:)
      TYPE (XPLEX), ALLOCATABLE :: temp(:,:), invtemp(:,:)
      TYPE (XPLEX), ALLOCATABLE :: Amat(:,:)
      TYPE (XPLEX)              :: FAFp(NLevs,NLevs)
      TYPE (XPLEX)              :: temp3(NLevs)
      TYPE (XPLEX)              :: temp4
      TYPE (XPLEX)              :: temp5(NLevs), temp6(NLevs)
      TYPE (XPLEX)              :: Xa(NLevs)
      TYPE (XPLEX)              :: Model_CO_layer(IIPAR, JJPAR, NLevs)
      TYPE (XPLEX)              :: Model_lcd(IIPAR, JJPAR, NLevs)
      TYPE (XPLEX)              :: COUNT_GRID_LOCAL(IIPAR,JJPAR)
      TYPE (XPLEX)              :: adj_factor(IIPAR,JJPAR,NLevs)
      INTEGER             :: ErrorFlag 

      !TYPE (XPLEX), intent(in) :: Model_CO_MR(iipar,jjpar,llpar)

      !=====================================================
      ! COMPUTE_COLUMN begins here
      !=====================================================

      Model_CO_layer(:,:,:)=0d0
      Model_lcd(:,:,:)=0d0
      Xa(:) = 0d0
      adj_factor(:,:,:) = 0d0
      COUNT_GRID_LOCAL(:,:) = 0d0

      ! Read and regrid input GEOS-Chem file on AIRS levels
      CALL REGRIDV_AIRS(Model_CO_layer)
      PRINT*,'### Model_CO_layer min, max: ',MINVAL(Model_CO_layer), 
     &      MAXVAL(Model_CO_layer)
      Model_lcd = Model_CO_layer


      print*, 'data corrected for 10% low bias in 200405'
      DO W = 1, NObss

         !print*, 'w is:', w
         !print*, 'if stmt stuf:',qual(W),NTraps(W),DNFlag(W),Tsurf(W)
         !call flush(6)
         ! Quality control flag
         ! NTraps > 1 means the averaging kernel exists
         ! Select daytime only measurements (12 hours centered
         ! around 1:30pm overpass time based on Colette's MOPITT code)
         ! LocalT is in minutes, so this is 7:30am - 7:30pm
         IF ( (qual(W) .eq. 0) .AND. (NTraps(W) .gt. 1) .AND. 
     &         (DNFlag(W) .eq. 'Day') .AND. 
     &         (Tsurf(W) .ge. 250) .and. (COcol(w) .gt. 0) )THEN

            !print*, 'got inside the if statement?'
            !call flush(6)

         !print*, 'w is:', w
         !print*, 'if stmt stuf:',qual(W),NTraps(W),DNFlag(W),Tsurf(W)
         !call flush(6)

            NTloc = NTraps(W)

             ! Look for model grid box corresponding to the observation:
            ! Get I corresponding to PLON(IND)
            ILON = INT( ( Longitude(W) + 180d0 ) / DISIZE + 1.5d0 )
            ! Handle date line correctly (bmy, 4/23/04)
            IF (ILON > IIPAR ) ILON = ILON - IIPAR 
            ! Get J corresponding to PLAT(IND)
            ILAT = INT( ( Latitude(W) +  90d0 ) / DJSIZE + 1.5d0 )

            IF(GET_HOUR() .EQ. OBS_HOUR_AIRS_CO(ILON,ILAT)) THEN
            ALLOCATE( Ftrans  (NTloc, NLevs) )
            ALLOCATE( Fpi     (NTloc, NLevs) )
            ALLOCATE( Fmat    (NLevs, NTloc) )
            ALLOCATE( temp    (NTloc, NTloc) )
            ALLOCATE( invtemp (NTloc, NTloc) )
            ALLOCATE( Amat    (NTloc, NTloc) )

            ! Initialize variables
            Ftrans(:,:) = 0d0
            Fpi(:,:)=0d0
            Fmat(:,:) = 0d0
            temp(:,:)=0d0
            invtemp(:,:)=0d0
            Amat(:,:)=0d0
            FAFp(:,:)=0d0
            Xa(:)=0d0
            temp3(:)=0d0
            temp4=0e0
            temp5(:) = 0e0
            temp6(:) = 0e0

            Psurf = SurfPressure( W )
            Ftrans = FMATRIX_COMP( Psurf, NTloc )
            Fmat = TRANSPOSE(Ftrans)

            temp = MATMUL( Ftrans, Fmat )

            CALL FINDINV( temp, invtemp, NTloc, ErrorFlag )
            IF (ErrorFlag .ne. 0) THEN
               !invtest(W)=1
               GOTO 291
            ENDIF

            Fpi = MATMUL( invtemp, Ftrans )
            Amat(:,:) = AvgKer(1:NTloc,1:NTloc,W)
            FAFp = MATMUL( MATMUL(Fmat,Amat), Fpi )

            ! Need to use LOG(X) where X is layer column density
            ! log(x')=log(x0)+FAFpi*log(x/x0)
            ! x'=exp[log(x0)+FAFpi*log(x/x0)]
            Xa = APRIORI_COMP( H2Ocd(:,W), BotLev(W) )

            DO L = 1, NLevs
               IF ( (Model_lcd(ILON,ILAT,L) .gt. 0) .AND. 
     &               (Xa(L) .gt. 0) ) THEN
                  temp3(L)=LOG(Model_lcd(ILON,ILAT,L)/Xa(L))
               ENDIF
            ENDDO

            temp3=MATMUL(FAFp,temp3)

            DO L = 1, NLevs
               IF (Xa(L) .gt. 0) THEN
                  temp3(L)=temp3(L)+LOG(Xa(L))
               ELSE
                  temp3(L)=0
               ENDIF
            ENDDO

            DO L = 1, NLevs
               IF (temp3(L) .gt. 0) THEN
                  temp3(L)=EXP(temp3(L))
               ELSE
                  temp3(L)=0
               ENDIF
            ENDDO

            ! in stand-alone code, we computed Model_col(w),
            ! in the adjoint/gc code, we compute CHK_STT_AIRS(I,J)
            !Model_col(W) = CALCtotcolden( temp3, PSurf)
            CHK_STT_AIRS(ILON,ILAT) = CHK_STT_AIRS(ILON,ILAT)
     &           + CALCtotcolden( temp3, PSurf)
       
            AIRS_COL_GRID(ILON,ILAT) = AIRS_COL_GRID(ILON,ILAT)
     &           + COcol(w)
            AIRSDOF_COL_GRID(ILON,ILAT) = AIRSDOF_COL_GRID(ILON,ILAT)
     &           + dofs(w)

            temp4 = CALCtotcolden( temp3, PSurf)
            COUNT_GRID_LOCAL(ILON,ILAT) = COUNT_GRID_LOCAL(ILON,ILAT)+1

            ! ADJOINT of AIRS retrieval (part 1 of 3) mak, 6/17/08
            DO L = 1, NLevs
               if (model_lcd(ilon,ilat,l) .gt. 0 ) then
                  temp5(L) = temp4/model_lcd(ILON,ILAT,L)
               endif
            ENDDO
            temp6 = MATMUL(FAFp,temp5) 
            adj_factor(ILON,ILAT,:) = temp6
            !print*, 'adj_factor',ilon,ilat,adj_factor(ilon,ilat,:)
            !print*, 'model col:', temp4
            !print*, 'x"/x:', temp5
            !print*, 'temp6:', temp6

291         CONTINUE

            IF ( ALLOCATED( Ftrans) ) DEALLOCATE( Ftrans  )
            IF ( ALLOCATED( Fpi   ) ) DEALLOCATE( Fpi     )
            IF ( ALLOCATED( Fmat  ) ) DEALLOCATE( Fmat    )
            IF ( ALLOCATED( temp  ) ) DEALLOCATE( temp    )
            IF ( ALLOCATED( invtemp)) DEALLOCATE( invtemp )
            IF ( ALLOCATED( Amat  ) )  DEALLOCATE( Amat    )

            ENDIF !OBS_HOUR
         ENDIF !quality flags etc

      ENDDO

      !PRINT*,'# of noninvertible matrix obs',SUM(invtest)
!      print*,'min/max of CHK_STT_AIRS:',minval(chk_stt_airs),

      !=======================================
      ! BIN OUTPUT INFO INTO MODEL GRID BOXES
      !=======================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( COUNT_GRID_LOCAL(I,J) .gt. 0. ) then   
            ! average AIRS
            AIRS_COL_GRID(I,J) = AIRS_COL_GRID(I,J)/
     &           COUNT_GRID_LOCAL(I,J)
            AIRSDOF_COL_GRID(I,J) = AIRSDOF_COL_GRID(I,J)/
     &           COUNT_GRID_LOCAL(I,J)

            ! average model
            CHK_STT_AIRS(I,J) = CHK_STT_AIRS(I,J)/COUNT_GRID_LOCAL(I,J) 
            
            ! average adjoint of AIRS retrieval (part 2 of 3)
            adj_factor(I,J,:) = adj_factor(I,J,:)/COUNT_GRID_LOCAL(I,J)
        ELSE
           AIRS_COL_GRID(I,J) = -999.
           AIRSDOF_COL_GRID(I,J) = -999.
           CHK_STT_AIRS(I,J) = -999.
           adj_factor(I,J,:) = -999.
        ENDIF
      ENDDO
      ENDDO
  
      print*,'min/max of CHK_STT_AIRS:',minval(chk_stt_airs),
     & maxval(chk_stt_airs)
      print*,'min/max of AIRS_COL_GRID:',minval(AIRS_COL_GRID),
     & maxval(airs_COL_GRID)

      ! ADJOINT of AIRS retrieval (part 3 of 3)
      ADJ_AIRS_ALL(:,:,:) = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, LL )
      DO LL = 1,NLevs
      DO L = 1,LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! d(CHK_STT_AIRS)/d(CHK_STT) = FAFp*(x'/x), then average
         ! then multiply by unit conversion factor (kg->molec/cm2)
         ! then flip the array vertically
         IF (adj_factor(I,J,1) .GE. 0) THEN
      
            ADJ_AIRS_ALL(I,J,L) = ADJ_AIRS_ALL(I,J,L) 
     &           + adj_factor(i,j,NLevs-LL+1)
     &           * ( 6.022d22 / (28.0d0 * GET_AREA_M2(J) ))
     &           * FRACTION(I,J,L,LL)
            !print*, 'adj_airs_all,i,j,l', i,j,l,adj_airs_all

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

c$$$      DO L = 1,LLPAR
c$$$      DO J = 1, JJPAR
c$$$      DO I = 1, IIPAR
c$$$         IF(ADJ_AIRS_ALL(I,J,L) .GT. 0) THEN 
c$$$            PRINT*, 'ADJ_AIRS_ALL:', I,J,L,ADJ_AIRS_ALL(I,J,L)
c$$$         ENDIF
c$$$      ENDDO
c$$$      ENDDO
c$$$      ENDDO

      END SUBROUTINE AIRS_COMPUTE_COLUMN

!--------------------------------------------------------------------------

      SUBROUTINE REGRIDV_AIRS( Model_CO_layer )

!*********************************************************************
! REGRIDS MODEL ARRAY FROM GEOS-CHEM LEVELS TO 100 AIRS LEVELS
!*********************************************************************    
! Note: This subroutine is copied and adapted slightly from Monika
!       Kopacz's REGRIDV_AIRS, which is a direct Fortran translation 
!       of my IDL code, which was in turn constructed from gamap's 
!       regridv.pro. It calls a subroutine REGRID_COLUMN, which is a 
!       Fortran translation of IDL code, which apparently was a
!       translation of Fortran code that we can no longer locate.
!       (jaf, 10/4/07)
!*********************************************************************    

      USE AIRSv5_MOD
      USE CHECKPT_MOD,    ONLY : CHK_PSC, CHK_STT
      USE GRID_MOD,             ONLY : GET_AREA_M2

#     include "CMN_SIZE" ! PTOP, LLPAR, JJPAR, IIPAR


      !TYPE (XPLEX), intent(in) :: Model_CO_MR(iipar,jjpar,llpar)
      TYPE (XPLEX), INTENT(INOUT):: Model_CO_layer(IIPAR, JJPAR, NLevs)
      TYPE (XPLEX)               :: STT_AIRS_VGRID(IIPAR,JJPAR,NLevs)
      TYPE (XPLEX)               :: InPEdge(LLPAR+1)
      TYPE (XPLEX)               :: ModelEdge(LLPAR+1)
      TYPE (XPLEX)               :: OutPEdge(NLevs+1)
      TYPE (XPLEX)               :: AIRSEdge(NLevs+1) 
      TYPE (XPLEX)               :: AIRSEdgePressure(NLevs+1)
      !TYPE (XPLEX)               :: FRACTION(IIPAR,JJPAR,LLPAR,NLevs)
      TYPE (XPLEX)               :: SurfP
      !TYPE (XPLEX)               :: STT_KG(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)               :: AirMass_AIRS(NLevs)
      TYPE (XPLEX)               :: AirMass_GC(LLPAR)
      TYPE (XPLEX)               :: TCVV

      INTEGER I, J, Y, M, L, LL, K
      LOGICAL, SAVE :: FIRST = .TRUE.
      LOGICAL, SAVE :: valid = .FALSE. 

      !=====================================================
      ! REGRIDV_AIRS begins here
      !=====================================================

      FRACTION(:,:,:,:) = 0d0
      AirMass_AIRS(:) = 0d0
      AirMass_GC(:) = 0d0
      Model_CO_layer(:,:,:)=0d0

      ModelPS(:,:) = CHK_PSC(:,:,2)

      ! TCVV is the ratio MW air / MW tracer
      TCVV   = 28.97d0 / 28.0d0

      ! Reorder pressure levels, so the first edge pressure is
      ! the surface, rather than the top of the atmosphere.
      AIRSEdge(1:NLevs) = Plev(NLevs:1:-1)

      !Assume first given edge is 0.01hPa
      AIRSEdge(NLevs+1) = 0.01

      !Store pressure edges
      AIRSEdgePressure = AIRSEDGE

      ! Convert to sigma scale
      SurfP=AIRSEdge(1)
      DO k = 1,NLevs+1
         AIRSEdge(k)=AIRSEdge(k)/SurfP
      ENDDO

      !-------------------
      ! REGRID DATA
      !-------------------
      
      !STT_KG(:,:,:)=0d0
      ! First need to convert input model v/v into kg, so need airmass
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         DO L = 1, LLPAR
            InPEdge(L)  = ( GET_PEDGE_JAF(I,J,L) )
            ! Need sigma grid value for kg conversion
            ModelEdge(L) = InPEdge(L)/ModelPS(I,J)
         ENDDO
         InPEdge(LLPAR+1) = PTOP
         ModelEdge(LLPAR+1) = InPEdge(LLPAR+1)/ModelPS(I,J)
         AirMass_GC = RVR_GetAirMass( ModelEdge, J, ModelPS(I,J), LLPAR)
         !DO L = 1, LLPAR
         !      ! STT_KG = v/v * kgair / (gair/gCO) = kg CO
         !   STT_KG(I,J,L) = CHK_STT(I,J,L,1)! * AirMass_GC(L) / TCVV
         !ENDDO
      ENDDO
      ENDDO

      STT_AIRS_VGRID(:,:,:) = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         !to be safe, remove junk values:
         fraction(i,j,:,:) = 0d0                  
            
         DO L = 1, LLPAR
            ! Pressure edges on INPUT and OUTPUT grids
            ! both in and out pressures in hPa
            InPEdge(L)  = ( GET_PEDGE_JAF(I,J,L) )
         ENDDO
         InPEdge(LLPAR+1) = PTOP
         OutPEdge(:) = AIRSEdgePressure(:) 

         !=====================================================
         ! Determine fraction of each INPUT box 
         ! which contributes to each OUTPUT box
         !=====================================================
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
            IF ( ( InPEdge(L) - InPedge(L+1) ) .lt. 1e-5 ) THEN
               goto 12   !NextL
            ENDIF

            ! Loop over OUTPUT layers
            DO LL = 1, NLevs

               IF( OutPEdge(LL) .lt. InPEdge(L) .and.    
     &              OutPEdge(LL) .lt. InPEdge(L+1) .and. 
     &              (LL .eq. 1) .and. (L.eq.1) ) THEN
                  Fraction(i,j,L,LL) = 1d0
                  ! Go to next iteration
                  goto 12  !NextL
               ENDIF
               !===================================================
               ! No contribution if:
               ! -------------------
               ! Bottom of OUTPUT layer above Top of INPUT layer OR
               ! Top    of OUTPUT layer below Bottom of INPUT layer
               ! ..unless it's the first layer in GC (mak, 8/15/07)
               !===================================================
               IF ( OutPEdge(LL) .lt. InPEdge(L+1) .OR.           
     &              OutPEdge(LL+1) .gt. InPEdge(L) ) THEN
                  goto 13   !NextLL
               ENDIF

               !===================================================
               ! Contribution if: 
               ! ----------------
               ! Entire INPUT layer in OUTPUT layer
               !===================================================
               IF  ( OutPEdge(LL) .ge. InPEdge(L) .AND.           
     &              OutPEdge(LL+1) .le. InPEdge(L+1) ) THEN
                     
                  Fraction(i,j,L,LL) = 1d0
                  !Indicate a valid contribution from L to LL
                  Valid = .true.
                     
                  ! Go to next iteration
                  goto 13       !NextLL
               ENDIF

               !==================================================
               ! Contribution if: 
               ! ----------------
               ! Top of OUTPUT layer in INPUT layer
               !==================================================
               IF ( OutPEdge(LL+1) .le. InPEdge(L)  .AND.        
     &              OutPEdge(LL)  .ge. InPEdge(L) ) THEN 
               
                  Fraction(i,j,L,LL) =(InPEdge(L) - OutPEdge(LL+1)) / 
     &                 ( InPEdge(L) -  InPEdge(L+1) ) 
                  ! Indicate a valid contribution from L to LL
                  Valid = .true.
                     
                  ! Go to next iteration
                  goto 13       !NextLL
               ENDIF
                  
               !==================================================
               ! Contribution if: 
               ! ----------------
               ! Entire OUTPUT layer in INPUT layer
               !==================================================
               IF ( OutPEdge(LL)   .le. InPEdge(L) .AND.         
     &              OutPEdge(LL+1) .ge. InPEdge(L+1) ) then 
                  
                  Fraction(i,j,L,LL)=(OutPEdge(LL) - OutPEdge(LL+1))/ 
     &                 ( InPEdge(L)  -  InPEdge(L+1) )
                     
                  ! Also add the to the first OUTPUT layer the fraction
                  ! of the first INPUT layer that is below sigma = 1.0
                  ! This is a condition that can be found in GEOS-3 data.
                  IF ( ( First                  )   .AND.        
     &                 ( LL .eq. 1              )   .AND.        
     &                 ( InPEdge(L) .gt. OutPEdge(1) ) ) THEN
                        
                     Fraction(i,j,L,LL) = Fraction(i,j,L,LL) +   
     &                    ( InPEdge(L) - OutPEdge(1)  ) /        
     &                    ( InPEdge(L) - InPEdge(L+1) )                
                        
                     ! We only need to do this once...
                     First = .false.
                  ENDIF
 
                  ! Indicate a valid contribution from L to LL
                  Valid = .true.
                     
                  ! Go to next iteration
                  goto 13       !NextLL
               ENDIF

               !===================================================
               ! Contribution if: 
               ! ----------------
               ! Bottom of OUTPUT layer in INPUT layer
               !===================================================
               IF ( OutPEdge(LL)   .ge. InPEdge(L+1) .AND.        
     &              OutPEdge(LL+1) .le. InPEdge(L+1) ) THEN
                     
                  Fraction(i,j,L,LL) = ( OutPEdge(LL) - InPEdge(L+1) ) / 
     &                 ( InPEdge(L)  - InPEdge(L+1) )
                     
                  ! Also add the to the first OUTPUT layer the fraction
                  ! of the first INPUT layer that is below sigma = 1.0
                  ! This is a condition that can be found in GEOS-3 data.
                  IF ( ( First         )   .AND.                  
     &                 ( LL .eq. 1      )   .AND.                 
     &                 ( InPEdge(L) .gt. OutPEdge(1) ) ) then 
                     
                     Fraction(i,j,L,LL) = Fraction(i,j,L,LL) +    
     &                    ( InPEdge(L) - OutPEdge(1)   ) /        
     &                    ( InPEdge(L) - InPEdge(L+1) )                
                        
                     ! We only need to do this once...
                     First = .false.
                  ENDIF
 
                  ! Indicate a valid contribution from L to LL
                  Valid = .true.
                     
                  ! Go to next iteration
                  goto 13       !NextLL
               ENDIF
 
 13            CONTINUE         !NextLL

            ENDDO               !LL

            !======================================================
            ! Consistency Check:
            ! ------------------
            ! If SUM( FRACTION(L,:) ) does not = 1, there is a problem.
            ! Test those INPUT layers (L) which make a contribution to 
            ! OUTPUT layers (LL) for this criterion.
            !
            !======================================================
            IF ( Valid ) THEN
               IF ( Abs( 1e0 - sum( Fraction(i,j,L,:))) .ge. 1e-4 ) THEN 
                  print*, 'Fraction does not add to 1'
                  print*, L, LL,sum( Fraction(i,j,L,:) )
                  print*, 'frac(5,:):', fraction(i,j,L,:)
                  print*, 'InPEdge:', InPEdge
                  print*, 'OutPEdge:', OutPEdge
               ENDIF
            ENDIF
      
 12         CONTINUE            !NextL
            
         ENDDO                  !L

         !==========================================================
         ! Compute "new" data -- multiply "old" data by fraction of
         ! "old" data residing in the "new" layer
         !==========================================================
         ! Map CO from GC to AIRS grid
         DO LL = 1 , NLevs
            DO L = 1 , LLPAR
               STT_AIRS_VGRID(I,J,LL) =                        
     &              STT_AIRS_VGRID(I,J,LL)                     
     &              + CHK_STT(I,J,L,1)*FRACTION(I,J,L,LL) !CHK_STT in kg
            ENDDO
         ENDDO

         IF(Abs( SUM(STT_AIRS_VGRID(I,J,:)) - SUM(CHK_STT(I,J,:,1))) 
     &        /SUM(CHK_STT(I,J,:,1)) .gt. 1e-5 ) THEN
         PRINT*, 'columns before and after regrid dont add up:'
         PRINT*, 'columns before and after regridding:'
         PRINT*, I,J,SUM(CHK_STT(I,J,:,1)),SUM(STT_AIRS_VGRID(I,J,:)) 
         PRINT*, 'InPEdge:', InPEdge
         PRINT*, 'OutPEdge:', OutPEdge
         PRINT*, 'CHK_STT'
         PRINT*, CHK_STT(I,J,:,1)
         PRINT*, 'STT_AIRS_VGRID:'
         PRINT*, STT_AIRS_VGRID(I,J,:)
      ENDIF
!!$
!!$            ! Airmass on output grid (in kg/box in each level)
!!$            AirMass_AIRS = RVR_GetAirMass( AIRSEdge, J, SurfP, NLevs )
!!$
!!$            ! Convert data from kg to [v/v]
!!$            ! Model_CO_MR = kgCO * gair/gCO / kgair = [v/v]
!!$            DO LL = 1, NLevs
!!$               Model_CO_MR(I,J,LL) = STT_AIRS_VGRID(I,J,LL) *  &
!!$                    TCVV/AirMass_AIRS(LL)
!!$            ENDDO

      ! Need to reverse the vertical order of the array...
      ! The surface should be L100 and TOA L1
      DO LL = 1, NLevs
         Model_CO_layer(I,J,LL) = STT_AIRS_VGRID(I,J,NLevs-LL+1)
      ENDDO

      ! Convert kg CO to layer column in molecules/cm^2
      ! Model_CO_layer = kgCO*Avogadro*(g/kg)*(m2/cm2)/[(g/mol)CO*Area]
      Model_CO_layer(I,J,:) = Model_CO_layer(I,J,:) *       
     &     6.022d22 / (28.0d0 * GET_AREA_M2(J) )
!            Model_CO_layer(I,J,:) = STT_AIRS_VGRID(I,J,:) *       &
!                 6.022d22 / (28.0d0 * Area (J) )


      ENDDO  !IIPAR
      ENDDO  !JJPAR

      END SUBROUTINE REGRIDV_AIRS

!-------------------------------------------------------------------------

      FUNCTION GET_PEDGE_JAF( I, J, L ) RESULT( PEDGE )
!
!******************************************************************************
!  Function GET_PEDGE returns the pressure at the bottom edge of level L.
!  (dsa, bmy, 8/20/02, 10/24/03)
!  This version has been slightly modified for use with AIRS regridding.
!  (jaf, 10/04/07)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) P (TYPE (XPLEX) ) : P_surface - P_top (PS-PTOP)
!  (2 ) L (INTEGER) : Pressure will be returned at the bottom edge of level L
!
!  NOTES:
!  (1 ) Bug fix: use PFLT instead of PFLT-PTOP for GEOS-4 (bmy, 10/24/03)
!******************************************************************************
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! LLPAR, PTOP

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      
      ! Local Variables
      INTEGER :: AS
      TYPE (XPLEX), ALLOCATABLE :: AP(:)
      TYPE (XPLEX), ALLOCATABLE :: BP(:)

      ! Return value
      TYPE (XPLEX)              :: PEDGE 

      !=================================================================
      ! GET_PEDGE_JAF begins here!
      !=================================================================
      ALLOCATE( AP( LLPAR + 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AP' )
      AP = 1d0

      ALLOCATE( BP( LLPAR + 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BP' )
      BP = 0d0

#if   defined( GRID30LEV )
      !----------------------
      ! GEOS-4 30 level grid
      !----------------------

      ! Ap [hPa] for 30 levels (31 edges)
      AP = (/  0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &     66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0, 
     &     198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0, 
     &     216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0, 
     &     127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &     56.387939d0,  40.175419d0,  28.367815d0,  19.791553d0, 
     &     9.292943d0,   4.076567d0,   1.650792d0,   0.616779d0, 
     &     0.211349d0,   0.066000d0,   0.010000d0 /)

      ! Bp [unitless] for 30 levels (31 edges)
      BP = (/  1.000000d0,   0.985110d0,   0.943290d0,   0.867830d0, 
     &     0.764920d0,   0.642710d0,   0.510460d0,   0.378440d0, 
     &     0.270330d0,   0.183300d0,   0.115030d0,   0.063720d0, 
     &     0.028010d0,   0.006960d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,   0.000000d0,   0.000000d0 /)

#else
      !----------------------
      ! GEOS-4 55 level grid
      !----------------------

      ! AP [hPa] for 55 levels (56 edges)
      AP = (/ 0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &     66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0, 
     &     198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0, 
     &     216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0, 
     &     127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &     66.603378d0,  56.387939d0,  47.643932d0,  40.175419d0, 
     &     33.809956d0,  28.367815d0,  23.730362d0,  19.791553d0, 
     &     16.457071d0,  13.643393d0,  11.276889d0,   9.292943d0, 
     &     7.619839d0,   6.216800d0,   5.046805d0,   4.076567d0, 
     &     3.276433d0,   2.620212d0,   2.084972d0,   1.650792d0, 
     &     1.300508d0,   1.019442d0,   0.795134d0,   0.616779d0, 
     &     0.475806d0,   0.365041d0,   0.278526d0,   0.211349d0, 
     &     0.159495d0,   0.119703d0,   0.089345d0,   0.066000d0, 
     &     0.047585d0,   0.032700d0,   0.020000d0,   0.010000d0 /)

      ! BP [unitless] for 55 levels (56 edges)
      BP = (/  1.000000d0,  0.985110d0,   0.943290d0,   0.867830d0, 
     &     0.764920d0,  0.642710d0,   0.510460d0,   0.378440d0, 
     &     0.270330d0,  0.183300d0,   0.115030d0,   0.063720d0, 
     &     0.028010d0,  0.006960d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0, 
     &     0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0 /)
#endif

      ! Here Ap is in [hPa] and Bp is [unitless].  
      ! For GEOS-4, we need to have PFLT as true surface pressure, 
      ! since Ap(1)=0 and Bp(1)=1.0.  This ensures that the true
      ! surface pressure will be returned for L=1. (bmy, 10/24/03)
      PEDGE = AP(L) + ( BP(L) * ModelPS(I,J) )

      IF ( ALLOCATED( AP   ) ) DEALLOCATE( AP   )
      IF ( ALLOCATED( BP   ) ) DEALLOCATE( BP   )

      ! Return to calling program
      END FUNCTION GET_PEDGE_JAF 

!------------------------------------------------------------------------------
 
      FUNCTION RVR_GetAirMass( Edge, J, SurfP,Levels) RESULT(AirMassloc)

      !====================================================================
      ! Internal function RVR_GETAIRMASS returns a column vector of air 
      ! mass given the vertical coordinates, the surface area,
      ! and surface pressure. (bmy, 12/19/03)
      !====================================================================

      USE GRID_MOD,             ONLY : GET_AREA_M2

#     include "CMN_SIZE"

      INTEGER, INTENT(IN)   :: J, Levels
      TYPE (XPLEX),  INTENT(IN)   :: SurfP
      TYPE (XPLEX)                :: AirMassloc(Levels)
      TYPE (XPLEX),  INTENT(IN)   :: Edge(Levels+1)
      INTEGER               :: L
      TYPE (XPLEX)                :: g100

      AirMassloc(:) = 0d0

      ! Constant 100/g 
      g100    = 100d0 / 9.8d0 

      ! Loop over levels
      ! airmass(L) = hPa * m2 * 1 * 100Pa/hPa * 1/(m/s2) = 
      !            = N * 1/(m/s2) = kg
      DO L = 1, Levels
         AirMassloc(L) = SurfP * GET_AREA_M2(J) *
     &          ( Edge(L) - Edge(L+1) ) * g100
      ENDDO
      
      END FUNCTION RVR_GetAirMass

!------------------------------------------------------------------------------

      SUBROUTINE READ_ERROR_VARIANCE 

      USE ERROR_MOD, 	ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE BPCH2_MOD
      USE FILE_MOD,    ONLY : IOERROR
      USE TIME_MOD,    ONLY : GET_TAU, EXPAND_DATE, GET_NYMD


      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

!#     include "define.h"
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      CHARACTER(LEN=255)     :: INPUT_FILE
      INTEGER                :: I, IOS, J,  L
      TYPE (XPLEX)                 :: TRACER(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: FILENAME_err

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
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
      INTEGER              :: YYYYMMDD

      INPUT_FILE = 'RRE_YYYYMMairsGlobal.bpch'

      IU_FILE = 66

      PRINT*, 'SET ERROR TO 50% AS WE SAVE AIRS FOR RRE CALCULATION'
      ERR_PERCENT(:,:) = 0.2
      GOTO 121
      !================================================================
      ! READ OLD RESTART FILE
      !================================================================
      YYYYMMDD = GET_NYMD()

      FILENAME_err = TRIM( INPUT_FILE )

      CALL EXPAND_DATE( FILENAME_err, YYYYMMDD, 0 )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'OBS ERROR FILE'
      WRITE( 6, 100 ) TRIM( FILENAME_err )
 100  FORMAT( 'READ_FILE: Reading ', a )

      ERR_PERCENT(:,:) = -999e0

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME_err )

      !=================================================================
      ! Read tracers -- store in the TRACER array
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS ) 
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
      
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
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
         PRINT*, 'Reading error for tau:', ztau0
              ! Make sure array dimensions are of global size
              ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
              !print*, 'inside the loop'
              !print*, 'max value is:', maxval(TRACER(:,:,:))
              !print*, 'min value is:', minval(TRACER)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF(TRACER(I,J) .ge. 0.05)THEN
               ERR_PERCENT(I,J) = TRACER(I,J)
               
            ELSEIF((TRACER(I,J).lt. 0.05).and.(TRACER(I,J).gt.0))THEN
               ERR_PERCENT(I,J) = 0.05
                    
            ELSE
               ERR_PERCENT(I,J) = TRACER(I,J)
                    
            ENDIF
                 
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDDO ! infinite reading loop

      ! Close file
      CLOSE( IU_FILE )      

 121  CONTINUE

      print*, 'max err value is:', maxval(ERR_PERCENT(:,:))
      print*, 'min err value is:', minval(ERR_PERCENT)

      END SUBROUTINE READ_ERROR_VARIANCE

!-----------------------------------------------------------------------------

      SUBROUTINE INIT_DOMAIN

      USE AIRSV5_MOD,           ONLY : DOMAIN_OBS

#include "CMN_SIZE"

      INTEGER  I, J, L, N

      ALLOCATE( DOMAIN_OBS( IIPAR,JJPAR ) )
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
     &          L == 1                    ! Only at the surface or col 
!!     &        .and.  J >= 10                 ! Not in antarctica
!     &         .and.  L == 8                   ! Only at ~500mb 
!     &         .and. (J .ge. 24)               ! only N.Hemisphere
     &         .and. (J .le. 38)               ! not poleward of 60N
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

      END MODULE AIRS_CO_OBS_MOD


      





