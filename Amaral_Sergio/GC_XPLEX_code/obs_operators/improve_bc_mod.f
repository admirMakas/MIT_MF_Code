!$Id: improve_bc_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE IMPROVE_BC_MOD
!
!******************************************************************************
! Mdoule IMPROVE_BC_MOD contains subroutines necessary to assimilate 
! observations of BC and OC from the IMPROVE network. 
! (yhmao, dkh, 01/13/12, adj32_013) 
!
! Notes
! (1 ) Based on the v6 adjoint improve obs operator
!
!******************************************************************************
  
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
      INTEGER, PARAMETER :: IU_IMPRV_ASCI = 901
      INTEGER, PARAMETER :: IU_IMPRV_BPCH = 902
      INTEGER, PARAMETER :: IU_AEROAVE    = 903
      INTEGER, PARAMETER :: IU_JSAVE      = 923

      INTEGER, PARAMETER :: LLAVE         = 1
 
      LOGICAL            :: DURING_IMPRV_OBS 


      TYPE (XPLEX), ALLOCATABLE :: IMPRV_BC(:,:,:)
      !TYPE (XPLEX), ALLOCATABLE :: IMPRV_OC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: AVE_BCPI(:,:)
      !TYPE (XPLEX), ALLOCATABLE :: AVE_OCPI(:,:)
      TYPE (XPLEX), ALLOCATABLE :: AVE_BCPO(:,:)
      !TYPE (XPLEX), ALLOCATABLE :: AVE_OCPO(:,:)
      TYPE (XPLEX), ALLOCATABLE :: ADJ_AVE_BCPI(:,:)
      !TYPE (XPLEX), ALLOCATABLE :: ADJ_AVE_OCPI(:,:)
      TYPE (XPLEX), ALLOCATABLE :: ADJ_AVE_BCPO(:,:)
      !TYPE (XPLEX), ALLOCATABLE :: ADJ_AVE_OCPO(:,:)
      
      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE IMPROVE_DATAPROC(YYYYMMDD)
!
!******************************************************************************
!  Subroutine IMPROVE_DATAPROC 
!  - this routine reads in raw IMPROVE data from a text file and puts in  
!    on a GEOS-Chem grid. It assumes that the text file has a row for  
!    each station, and is all of the data for a singe day.  Currently,
!    it is set for "imprv.20050730", which would by july 30th 2005. I
!    would  change this string to correspond to the name of your data text
!    file, and rerun just this subroutine.  You'll have to call it once for
!    every day of text data you need to process, so maybe put it in a loop.
!  - The columns of the text file are: YYYYMMDD, Lat, Lon, BCrate concentration,
!    BCrate uncertainty, BCrate detection limit
!    - the result is a bpch file with IMPROVE data.
!
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : Current date
!     
!  NOTES:
!  (1 ) Now add model resultion extension to binary file name. (dkh, 11/28/06) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_IJ, GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : GET_TAU
      USE TIME_MOD,    ONLY : EXPAND_DATE

#     include "CMN_SIZE"
!#     include "CMN_SETUP"   ! DATA_DIR
      
      INTEGER, INTENT(IN)   :: YYYYMMDD
      ! Arguments
    
      ! Local variables 
      INTEGER :: I, J, I0, J0, IIJJ(2), HHMMSS_dum
      CHARACTER(LEN=120) :: FILENAME, FILENAME1
      CHARACTER(LEN=120) :: READ_FILENAME, WRITE_FILENAME
      TYPE (XPLEX)  :: BC_BAR(IIPAR,JJPAR)!,  OC_BAR(IIPAR,JJPAR)
      TYPE (XPLEX)  :: SBC_MES(IIPAR,JJPAR)!, SOC_MES(IIPAR,JJPAR)
      !TYPE (XPLEX)  :: SBC_REP(IIPAR,JJPAR), SOC_REP(IIPAR,JJPAR)
      TYPE (XPLEX)  :: SBC_RES(IIPAR,JJPAR)!, SOC_RES(IIPAR,JJPAR)
      TYPE (XPLEX)  :: BC_MMDL(IIPAR,JJPAR)!, OC_MMDL(IIPAR,JJPAR)
      TYPE (XPLEX)  :: DELTA_BC!,             DELTA_OC
      INTEGER :: N(IIPAR,JJPAR), IOS
      LOGICAL :: EOF
      TYPE (XPLEX)               :: DAT(IIPAR,JJPAR,4)

      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UBC     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1 

      ! Variables read from file
      INTEGER  :: YYYYMMDD1
      TYPE (XPLEX)   :: LAT, LON
      TYPE (XPLEX)   :: BC, BC_UNC, BC_MDL
      !TYPE (XPLEX)   :: OC, OC_UNC, OC_MDL

      !=================================================================
      ! IMPROVE_DATAPROC begins here!
      !=================================================================
      BC_BAR(:,:)  = 0d0
      !OC_BAR(:,:)  = 0d0
      !SBC_REP(:,:) = 0d0 
      !SOC_REP(:,:) = 0d0
      SBC_MES(:,:) = 0d0 
      SBC_MES(:,:) = 0D0
      BC_MMDL(:,:) = 0d0 
      !OC_MMDL(:,:) = 0D0
      SBC_RES(:,:) = 0d0 
      !SOC_RES(:,:) = 0d0
      N(:,:)        = 0d0
      EOF           = .FALSE.

       print*,2 !yhmao 
      FILENAME = TRIM( 'imprv.YYYYMMDD' )
      ! #######
      ! manually enter the name of the file to process here
      ! ( also have to enter it in CALC_SRES )
      ! #######
      CALL EXPAND_DATE( FILENAME, YYYYMMDD,HHMMSS_dum  ) 
      ! Add input file suffix (confusingly named .out)
      FILENAME1='/qb6/yhmao/geos-chem/adjoint/new/gcadj_std/obsdata/'//
     &TRIM( FILENAME )
      READ_FILENAME =TRIM( FILENAME1 ) // TRIM( '.txt' )

        
      print*,3 !yhmao
      WRITE(6,*) ' IMPROVE_DATAPROC: reading file: ', READ_FILENAME

      OPEN( IU_IMPRV_ASCI,  FILE = TRIM( READ_FILENAME ),
     &      STATUS = 'OLD', IOSTAT = IOS )

      DO WHILE (.not. EOF ) 

         ! Read ascii text file
         READ( IU_IMPRV_ASCI, *, IOSTAT = IOS ) 
     &                                YYYYMMDD1, LAT,     LON, 
     &                                BC,      BC_UNC, BC_MDL!, 
      !&                                OC,      OC_UNC, OC_MDL

         ! dkh debug
         print*, ' YYYYMMDD = ', YYYYMMDD1 
         print*, ' LAT      = ', LAT      
         print*, ' LON      = ', LON      
         print*, ' BC      = ', BC      
         print*, ' BC_UN   = ', BC_UNC
         print*, ' BC_MDL  = ', BC_MDL
         !print*, ' OC      = ', OC      
         !print*, ' OC_UN   = ', OC_UNC
         !print*, ' OC_MDL  = ', OC_MDL
         print*, ' IOS      = ', IOS

         IF ( IOS < 0 ) THEN 

            EOF = .TRUE. 
 
         ELSEIF( IOS > 0 ) THEN
   
            CALL ERROR_STOP( 'Error reading improve.out file', 
     &                       'improve_mod.f' )
 
         ENDIF 

         ! Quality check
         IF (
     &         BC < 0       .or. 
      !&         OC < 0       .or. 
     &         BC < BC_MDL .or.
      !&         OC < OC_MDL .or. 
     &         BC < BC_UNC !.or. 
      !&         OC < OC_UNC 
     &                            ) CYCLE


         ! Get grid box 
         IIJJ = GET_IJ( LON, LAT)
         I    = IIJJ(1)
         J    = IIJJ(2)
   
         ! Update local count 
         N(I,J) = N(I,J) + 1
      
         ! dkh debug 
         print*, 'LON, LAT, I, J, N' , LON, LAT, I, J, N(I,J)
         print*, 'BC, BC_UNC, BC_MDL', BC, BC_UNC, BC_MDL
         !print*, 'OC, OC_UNC, OC_MDL', OC, OC_UNC, OC_MDL
     
         ! Update mean
         DELTA_BC     = BC - BC_BAR(I,J)
         !DELTA_OC     = OC - OC_BAR(I,J)
   
         BC_BAR(I,J)  = BC_BAR(I,J) + DELTA_BC / N(I,J)
         !OC_BAR(I,J)  = OC_BAR(I,J) + DELTA_OC / N(I,J)
      
         ! Update representational error
!         SBC_REP(I,J) = SBC_REP(I,J) 
!     &                 + DELTA_BC * ( BC - BC_BAR(I,J))
!         SOC_REP(I,J) = SOC_REP(I,J) 
!     &                 + DELTA_OC * ( OC - OC_BAR(I,J))

         ! Update the maximum minimum detection limit 
         BC_MMDL(I,J)  = MAX(BC_MMDL(I,J),BC_MDL)
         !OC_MMDL(I,J)  = MAX(OC_MMDL(I,J),OC_MDL)
   
         ! Update measurement error
         ! These are already variances, no need to square
         ! them. (dkh, 04/28/08) 
         !SBC_MES(I,J) = SBC_MES(I,J) + BC_UNC ** 2
         !SOC_MES(I,J) = SOC_MES(I,J) + OC_UNC ** 2 
         SBC_MES(I,J) = SBC_MES(I,J) + BC_UNC ** 2
         !SOC_MES(I,J) = SOC_MES(I,J) + OC_UNC ** 2 




         IF ( IOS < 0 ) THEN 

            EOF = .TRUE. 
 
         ELSEIF( IOS > 0 ) THEN
   
            CALL ERROR_STOP( 'Error reading improve.out file', 
     &                       'improve_mod.f' )
 
         ENDIF 

      ENDDO


      ! Close file
      CLOSE( IU_IMPRV_ASCI )

      WRITE(6,*) 'Done reading improve data file '
      WRITE(6,*) 'Number of good data points found: ', SUM(N(:,:))

      
      ! Finalize Statistical quantities
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO I = 1, IIPAR
      DO J = 1, JJPAR

         IF ( N(I,J) /= 0 ) THEN 
   
            !IF ( N(I,J) > 1 ) THEN 
            !
            !   SBC_REP(I,J) = SBC_REP(I,J) / ( N(I,J) - 1 ) 
            !   SOC_REP(I,J) = SOC_REP(I,J) / ( N(I,J) - 1 ) 
            !
            !ENDIF 
     
            !SBC_MES(I,J) = SBC_MES(I,J) / ( N(I,J) ** 2 )
            !SOC_MES(I,J) = SOC_MES(I,J) / ( N(I,J) ** 2 )
            !SBC_MES(I,J) = SQRT( SBC_MES(I,J) ) / N(I,J) 
            !SOC_MES(I,J) = SQRT( SOC_MES(I,J) ) / N(I,J)
            SBC_MES(I,J) =  SBC_MES(I,J) / N(I,J) 
            !SOC_MES(I,J) =  SOC_MES(I,J) / N(I,J)

         ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Get resolution error
      !CALL CALC_SRES(N, SBC_RES, SOC_RES)

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'Improve data file ' 
      CATEGORY = 'IJ-IMP-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UBC     = 'ug/m3'

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Add ADJ_DIR prefix to filename
      !WRITE_FILENAME = TRIM( FILENAME ) // TRIM('.bpch.') // 
      WRITE_FILENAME = TRIM( FILENAME1 ) // 
     &                 TRIM('.v2') // 
     &                 TRIM('.bpch.') // 
     &                 GET_RES_EXT()

      WRITE( 6, 100 ) TRIM( WRITE_FILENAME )
 100  FORMAT( '     - IMPROVE_DATAPROC: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_IMPRV_BPCH, WRITE_FILENAME, TITLE )

      ! BC
      ! Temporarily store data in DAT
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         DAT(I,J,1) = BC_BAR(I,J)
         !DAT(I,J,2) = SBC_REP(I,J)
         DAT(I,J,2) = BC_MMDL(I,J)
         DAT(I,J,3) = SBC_MES(I,J)
         DAT(I,J,4) = SBC_RES(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_IMPRV_BPCH,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UBC,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     4,     I0+1,
     &            J0+1,      1,         DAT )

      ! dkh debug
      print*, 'TOTAL BC : ', SUM(DAT(:,:,:))

      ! OC
      ! Temporarily store data in DAT
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
      !DO J = 1, JJPAR
      !DO I = 1, IIPAR

         !DAT(I,J,1) = OC_BAR(I,J)
         !DAT(I,J,2) = SOC_REP(I,J)
         !DAT(I,J,2) = OC_MMDL(I,J)
         !DAT(I,J,3) = SOC_MES(I,J)
         !DAT(I,J,4) = SOC_RES(I,J)

      !ENDDO
      !ENDDO
!!$OMP END PARALLEL DO

      !CALL BPCH2( IU_IMPRV_BPCH,    MODELNAME, LONRES,    LATRES,
      !&            HALFPOLAR, CENTER180, CATEGORY,  2,
      !&            UBC,      GET_TAU(), GET_TAU(), RESERVED,
      !&            IIPAR,     JJPAR,     4,     I0+1,
      !&            J0+1,      1,         DAT )

      ! dkh debug
      !print*, 'TOTAL OC : ', SUM(DAT(:,:,:))
      !print*, 'TOTAL OCa : ', SUM(OC_BAR(:,:))
      !print*, 'TOTAL OCb : ', SUM(SOC_REP(:,:))
      !print*, 'TOTAL OCm : ', SUM(OC_MMDL(:,:))
      !print*, 'TOTAL OCc : ', SUM(SOC_MES(:,:))


      ! Close file
      CLOSE( IU_IMPRV_BPCH )

      ! Return to calling program
      END SUBROUTINE IMPROVE_DATAPROC

!------------------------------------------------------------------------------



      SUBROUTINE READ_IMPRV_BPCH( YYYYMMDD )
!
!******************************************************************************
!  Subroutine READ_IMPRV_BPCH 
!  - reads in the IMPROVE bpch files that were made in SUBROUTINE  
!    IMPROVE_DATAPROC
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD, ONLY   : OPEN_BPCH2_FOR_READ, GET_RES_EXT
      USE ERROR_MOD, ONLY   : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,  ONLY   : IOERROR
      USE TIME_MOD,  ONLY   : EXPAND_DATE

#     include "CMN_SIZE"  
!#     include "CMN_SETUP"   ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN)   :: YYYYMMDD
    
      ! Local variables 
      INTEGER             :: IOS, I, J, L, HHMMSS_dum
      CHARACTER(LEN=120)  :: FILENAME, READ_FILENAME

       ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)           :: LONRES,    LATRES, TEMP4(IIPAR,JJPAR,4)
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UBC     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_IMPRV_BPCH begins here!
      !=================================================================
      FILENAME = TRIM( 'imprv.YYYYMMDD' )
 
      ! Replace token with actual year, month and day. 
      ! Also pass a dummy integer for HHMMSS, but it won't be used
      ! since the filename does not depend upon HHMMSS     
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum )
 
      READ_FILENAME = '/qb6/yhmao/geos-chem/adjoint/new/gcadj_std/
!     &                TRIM( 'improve_2006/' ) // 
!     &                TRIM( FILENAME ) // TRIM( '.bpch.' ) //
     &obsdata/' //TRIM( FILENAME ) // 
     &                TRIM( '.v2' ) //
     &                TRIM( '.bpch.' ) //
     &                GET_RES_EXT()

      print*, 'READ_IMPRV_BPCH: reading file : ', READ_FILENAME

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_IMPRV_BPCH, READ_FILENAME )

      !------------------
      ! BC
      !------------------
      READ( IU_IMPRV_BPCH, IOSTAT=IOS )
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_IMPRV_BPCH,'read_imprv_bpch:4' )

      READ( IU_IMPRV_BPCH, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UBC, ZTAU0,  ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP
  
      print*, 'ntracer = ', ntracer
      
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_IMPRV_BPCH,'read_imprv_bpch:5')

      READ( IU_IMPRV_BPCH, IOSTAT=IOS ) 
     &    ( ( (TEMP4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      IMPRV_BC(:,:,:) = TEMP4(:,:,:)

      !IF ( IOS /= 0 ) THEN 
         print*, 'ios = ', IOS
         print*, 'NI, NJ, NL = ', NI, NJ, NL 
         print*, 'SUM BC  = ', SUM(IMPRV_BC(:,:,:))
         !print*, ' BC  = ', IMPRV_BC(25,65,:)
         !print*, ' TEMP4 = ', TEMP4(25,65,:)


         !CALL IOERROR( IOS,IU_IMPRV_BPCH,'read_imprv_bpch:6')
      !ENDIF

      ! Only process improve data 
      IF ( CATEGORY(1:8) /= 'IJ-IMP-$' ) THEN
         CALL ERROR_STOP( 'Wrong data type', 'read_imprv_bpch')
      ENDIF


      !------------------
      ! OC
      !------------------
      !READ( IU_IMPRV_BPCH, IOSTAT=IOS )
      !&     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      !IF ( IOS > 0 ) 
      !&   CALL IOERROR( IOS,IU_IMPRV_BPCH,'read_checkpt_file:4' )

      !READ( IU_IMPRV_BPCH, IOSTAT=IOS )
      !&     CATEGORY, NTRACER,  UBC, ZTAU0,  ZTAU1,  RESERVED,
      !&     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
      !&     NSKIP
      
      !IF ( IOS /= 0 ) 
      !&   CALL IOERROR( IOS,IU_IMPRV_BPCH,'read_checkpt_file:5')

      !READ( IU_IMPRV_BPCH, IOSTAT=IOS ) 
      !&    ( ( (TEMP4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
 
      !IMPRV_OC(:,:,:) = TEMP4(:,:,:)

!      IF ( IOS /= 0 ) 
!     &   CALL IOERROR( IOS,IU_IMPRV_BPCH,'read_checkpt_file:6')

      ! Only process improve data 
      !IF ( CATEGORY(1:8) /= 'IJ-IMP-$' ) THEN
      !   CALL ERROR_STOP( 'Wrong data type', 'read_imprv_bpch')
      !ENDIF



!      ! Make sure array dimensions are of global size
!      ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!      CALL CHECK_DIMENSIONS( NI, NJ, NL )

      ! dkh debug
      print*, 'TOTAL BC = ', SUM(IMPRV_BC(:,:,:))
      !print*, 'TOTAL OC = ', SUM(IMPRV_OC(:,:,:))
      !print*, 'TOTAL OCa : ', SUM(IMPRV_OC(:,:,1))
      !print*, 'TOTAL OCm : ', SUM(IMPRV_OC(:,:,2))
      !print*, 'TOTAL OCc : ', SUM(IMPRV_OC(:,:,3))


      ! Return to calling program
      END SUBROUTINE READ_IMPRV_BPCH

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_AEROAVE_FILE( YYYYMMDD )
!
!******************************************************************************
!  Subroutine MAKE_AEROAVE_FILE saves daily average concentrations of 
!  OC and BC aerosol. (dkh, 11/18/06)  
!  - makes a bpch file with model values that match the time and
!    locations of the IMPROVE data.  This is called during the
!    forward part of the inverse model, from geos_chem_mod.f
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : Date of average                          [uBC]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : GET_TAU,     EXPAND_DATE

#     include "CMN_SIZE"    ! Size params
!#     include "CMN_ADJ"     ! ADJ_DIR

      ! Arguments
      INTEGER              :: YYYYMMDD
    
      ! Local variables 
      INTEGER              :: I, J, I0, J0, HHMMSS_dum
      CHARACTER(LEN=120)   :: FILENAME
      TYPE (XPLEX)               :: DAT(IIPAR,JJPAR,LLAVE)

      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UBC     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1 

      !=================================================================
      ! MAKE_AEROAVE_FILE begins here!
      !=================================================================

      FILENAME = TRIM( 'aero.ave.YYYYMMDD' )

      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum)

      ! Append data directory prefix
      FILENAME = '/qb6/yhmao/geos-chem/adjoint/new/gcadj_std/obsdata/'//
     &           TRIM( FILENAME )

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'Average Aerosol data file ' 
      CATEGORY = 'IJ-AVE-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UBC     = 'ug/m3'

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
 100  FORMAT( '     - MAKE_AEROAVE_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_AEROAVE, FILENAME, TITLE )

      ! BC
      ! Temporarily store data in DAT
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         DAT(I,J,1) = AVE_BCPI(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_AEROAVE,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UBC,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     LLAVE,     I0+1,
     &            J0+1,      1,         DAT )

      ! dkh debug
      print*, 'SUM AVE_BCPI : ', SUM(DAT(:,:,:))
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         DAT(I,J,1) = AVE_BCPO(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_AEROAVE,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  2,
     &            UBC,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     LLAVE,     I0+1,
     &            J0+1,      1,         DAT )

      ! dkh debug
      print*, 'SUM AVE_BCPO : ', SUM(DAT(:,:,:))

      ! OC
      ! Temporarily store data in DAT
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
      !DO J = 1, JJPAR
      !DO I = 1, IIPAR

        ! DAT(I,J,1) = AVE_OCPI(I,J)

      !ENDDO
      !ENDDO
!!$OMP END PARALLEL DO

      ! CALL BPCH2( IU_AEROAVE,    MODELNAME, LONRES,    LATRES,
      !&            HALFPOLAR, CENTER180, CATEGORY,  3,
      !&            UBC,      GET_TAU(), GET_TAU(), RESERVED,
      !&            IIPAR,     JJPAR,     LLAVE,     I0+1,
      !&            J0+1,      1,         DAT )

!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
      !DO J = 1, JJPAR
      !DO I = 1, IIPAR

       !  DAT(I,J,1) = AVE_OCPO(I,J)

      !ENDDO
      !ENDDO
!!$OMP END PARALLEL DO

      !CALL BPCH2( IU_AEROAVE,    MODELNAME, LONRES,    LATRES,
      !&            HALFPOLAR, CENTER180, CATEGORY,  4,
      !&            UBC,      GET_TAU(), GET_TAU(), RESERVED,
      !&            IIPAR,     JJPAR,     LLAVE,     I0+1,
      !&            J0+1,      1,         DAT )    


      ! Close file
      CLOSE( IU_AEROAVE )

      ! Return to calling program
      END SUBROUTINE MAKE_AEROAVE_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_AEROAVE_FILE( YYYYMMDD )
!
!******************************************************************************
!  Subroutine READ_AEROAVE_FILE 
!  - called during the adjoint part of the inverse model.  Reads in the  
!      24h average concentrations saved from the forward part that were  
!     written during MAKE_AEROAVE_FILE
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD, ONLY   : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD, ONLY   : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,  ONLY   : IOERROR
      USE TIME_MOD,  ONLY   : EXPAND_DATE

#     include "CMN_SIZE"    ! Size params 
!#     include "CMN_ADJ"     ! ADJ_DIR

      ! Arguments
      INTEGER, INTENT(IN)   :: YYYYMMDD
    
      ! Local variables 
      INTEGER             :: IOS, I, J, L, HHMMSS_dum
      CHARACTER(LEN=120)  :: FILENAME

       ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)         :: LONRES,    LATRES, TEMP4(IIPAR,JJPAR,3)
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UBC     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_AEROAVE_FILE begins here!
      !=================================================================
      FILENAME = TRIM( 'aero.ave.YYYYMMDD' )
 
      ! Replace token with actual year, month and day. 
      ! Also pass a dummy integer for HHMMSS, but it won't be used
      ! since the filename does not depend upon HHMMSS     
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum )
 
      FILENAME = '/qb6/yhmao/geos-chem/adjoint/new/gcadj_std/obsdata/'//
     &           TRIM( FILENAME ) 

      print*, 'READ_AEROAVE_FILE: reading file - ', FILENAME

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_AEROAVE, FILENAME )

      !------------------
      ! BC
      !------------------
      READ( IU_AEROAVE, IOSTAT=IOS )
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_AEROAVE,'read_imprv_bpch:4' )

      READ( IU_AEROAVE, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UBC, ZTAU0,  ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP
  
      print*, 'ntracer = ', ntracer
      
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_AEROAVE,'read_imprv_bpch:5')

      READ( IU_AEROAVE, IOSTAT=IOS ) 
     &    ( ( (TEMP4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      AVE_BCPI(:,:) = TEMP4(:,:,1)

      IF ( IOS /= 0 ) THEN 
         CALL IOERROR( IOS,IU_AEROAVE,'read_imprv_bpch:6')
      ENDIF

      ! Only process improve data 
      IF ( CATEGORY(1:8) /= 'IJ-AVE-$' ) THEN
         CALL ERROR_STOP( 'Wrong data type', 'read_imprv_bpch')
      ENDIF
      
      !BCPO
      READ( IU_AEROAVE, IOSTAT=IOS )
     &MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_AEROAVE,'read_imprv_bpch:4' )

      READ( IU_AEROAVE, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UBC, ZTAU0,  ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP
  
      print*, 'ntracer = ', ntracer
      
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_AEROAVE,'read_imprv_bpch:5')

      READ( IU_AEROAVE, IOSTAT=IOS ) 
     &    ( ( (TEMP4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      AVE_BCPO(:,:) = TEMP4(:,:,1)

      IF ( IOS /= 0 ) THEN 
         CALL IOERROR( IOS,IU_AEROAVE,'read_imprv_bpch:6')
      ENDIF

      ! Only process improve data 
      IF ( CATEGORY(1:8) /= 'IJ-AVE-$' ) THEN
         CALL ERROR_STOP( 'Wrong data type', 'read_imprv_bpch')
      ENDIF

      !------------------
      ! OC
      !------------------
      !READ( IU_AEROAVE, IOSTAT=IOS )
      !&     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      !IF ( IOS > 0 ) 
      !&   CALL IOERROR( IOS,IU_AEROAVE,'read_checkpt_file:4' )

      !READ( IU_AEROAVE, IOSTAT=IOS )
      !&     CATEGORY, NTRACER,  UBC, ZTAU0,  ZTAU1,  RESERVED,
      !&     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
      !&     NSKIP
      
      !IF ( IOS /= 0 ) 
      !&   CALL IOERROR( IOS,IU_AEROAVE,'read_checkpt_file:5')

      !READ( IU_AEROAVE, IOSTAT=IOS ) 
      !&    ( ( (TEMP4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
 
      !AVE_OCPI(:,:) = TEMP4(:,:,1)

      ! IF ( IOS /= 0 ) 
      !&   CALL IOERROR( IOS,IU_AEROAVE,'read_checkpt_file:6')

      ! Only process improve data 
      !IF ( CATEGORY(1:8) /= 'IJ-AVE-$' ) THEN
      !   CALL ERROR_STOP( 'Wrong data type', 'read_imprv_bpch')
      !ENDIF
      
      !OCPO
      !READ( IU_AEROAVE, IOSTAT=IOS )
      !&     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      !IF ( IOS > 0 ) 
      !&   CALL IOERROR( IOS,IU_AEROAVE,'read_checkpt_file:4' )

      !READ( IU_AEROAVE, IOSTAT=IOS )
      !&     CATEGORY, NTRACER,  UBC, ZTAU0,  ZTAU1,  RESERVED,
      !&     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
      !&     NSKIP
      
      !IF ( IOS /= 0 ) 
      !&   CALL IOERROR( IOS,IU_AEROAVE,'read_checkpt_file:5')

      !READ( IU_AEROAVE, IOSTAT=IOS ) 
      !&    ( ( (TEMP4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
 
      !AVE_OCPO(:,:) = TEMP4(:,:,1)

      !IF ( IOS /= 0 ) 
      !&   CALL IOERROR( IOS,IU_AEROAVE,'read_checkpt_file:6')

      ! Only process improve data 
      !IF ( CATEGORY(1:8) /= 'IJ-AVE-$' ) THEN
       !  CALL ERROR_STOP( 'Wrong data type', 'read_imprv_bpch')
      !ENDIF
      !------------------
 

      ! dkh debug
      print*, 'SUM AVE_BCPI = ', SUM(AVE_BCPI(:,:))
      !print*, 'SUM AVE_OCPI = ', SUM(AVE_OCPI(:,:))
      print*, 'SUM AVE_BCPO = ', SUM(AVE_BCPO(:,:))
      !print*, 'SUM AVE_OCPO = ', SUM(AVE_OCPO(:,:))


      IF ( NL /= 1 )  CALL ERROR_STOP( 'wrong dimension', 
     &                                 'read_aerosave')

      ! Close file
      CLOSE( IU_AEROAVE )

      ! Return to calling program
      END SUBROUTINE READ_AEROAVE_FILE

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_IMPRV_OBS_START( DIRECTION ) RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_IMPRV_OBS returns TRUE if it is the start of a day 
!  for which we have improve observations. (dkh, 11/18/06)  
!  - determines wether or not the simulation is at the beginning of a  
!    day where there are IMPROVE observations.  Based on the months that  
!    you choose to analyze, you will need to adjust this line:
!        IF ( MOD( DATE(1) - 20050400 + 3, 3 ) == 0 ) THEN
!                         FLAG = .TRUE.
!    so that it is looking at the correct 1-out-of-3 days for your time period.
!    - assumes the entire dataset is in the same time zone, which I think is 
!    reasonable for 24 hour average measurements over just the continental US.
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules 
      USE TIME_MOD,    ONLY : GET_TIME_AHEAD
      USE TIME_MOD,    ONLY : GET_NYMDb

      ! Function arguments
      INTEGER :: DIRECTION 

      ! Function value
      LOGICAL :: FLAG
 
      ! Local variables
      INTEGER :: DATE(2)

      !=================================================================
      ! ITS_TIME_FOR_IMPRV_OBS_START begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.
      DURING_IMPRV_OBS = .FALSE.
      
      ! Get time in midwest
      DATE = GET_TIME_AHEAD( - 7 * 60 )

      ! Check if it's midnight
      IF ( DATE(2) == 00 ) THEN 

         ! Check if its a day where we have observations
         ! for 200201:
         !IF ( MOD( DATE(1) - 20020100 + 1, 3 ) == 0 ) THEN 
         ! for 200104:
         !IF ( MOD( DATE(1) - 20010400 + 2, 3 ) == 0 ) THEN 
         ! for 200107:
         !IF ( MOD( DATE(1) - 20010400 + 3, 3 ) == 0 ) THEN
         ! for 200507:
         IF ( MOD( DATE(1) - 20060701 + 3, 3 ) == 0 ) THEN 
            FLAG = .TRUE. 

            WRITE(6,*) ' ITS_TIME_FOR_IMPRV_OBS_START at ',
     &         DATE(1), DATE(2)

            ! forward calculation
            IF ( DIRECTION > 0 ) THEN 

               ! moving into new measurement period
               DURING_IMPRV_OBS = .TRUE.
 
            ! backward calculation 
            ELSEIF ( DIRECTION < 0 ) THEN
    
               ! leaving measurement period 
               DURING_IMPRV_OBS = .FALSE.
    
            ENDIF 
       !  ELSE
       !   DURING_IMPRV_OBS = .FALSE.
         ENDIF
       !ELSE 
       !     DURING_IMPRV_OBS = .FALSE.
      ENDIF 

      print*, 'TIME_FOR_IMPRV_OBS_START: DATE = ', DATE, 
     &        DURING_IMPRV_OBS

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_IMPRV_OBS_START

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_IMPRV_OBS() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_IMPRV_OBS returns TRUE if it is during a day in which
!  we have improve measurements. (dkh, 11/18/06)  
!
!  NOTES:
!  - returns TRUE if the simulation is currently in within a day for  
!      which there are IMPROVE observations
!******************************************************************************
!
       !USE TIME_MOD,    ONLY : GET_TIME_AHEAD



      ! Local variables
      !INTEGER :: DATE(2)

      ! Function value
      LOGICAL :: FLAG
 
      !=================================================================
      ! ITS_TIME_FOR_IMPRV_OBS begins here!
      !=================================================================
      
      !FLAG = .FALSE.

      ! DATE = GET_TIME_AHEAD( - 7 * 60 )
      ! Check if it's midnight
      !IF ( DATE(2) == 00 ) THEN
      !IF ( MOD( DATE(1) - 20060102 + 3, 3 ) == 0 ) THEN
      !   FLAG = .TRUE.
      !ENDIF
      !ENDIF
       FLAG =  DURING_IMPRV_OBS 

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_IMPRV_OBS

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_IMPRV_OBS_STOP( DIRECTION ) RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_IMPRV_OBS returns TRUE if it is the end of a day for
!  which we have improve observations. (dkh, 11/18/06)  
!
!  NOTES:
!  - determines if it is the end of a day for which there are IMPROVE  
!    observations.   Like ITS_TIME_FOR_IMPRV_OBS_START, you will need to  
!    modify the exact date range here.
!  - also assumes the entire dataset is in the same time zone,which I  
!        think is reasonable for 24 hour average measurements.
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE TIME_MOD,    ONLY  : GET_TIME_AHEAD 
      USE TIME_MOD,    ONLY  : GET_NYMDb

      ! Function argument 
      INTEGER :: DIRECTION

      ! Function value
      LOGICAL :: FLAG
 
      ! Local variables
      INTEGER :: DATE(2)

      !=================================================================
      ! ITS_TIME_FOR_IMPRV_OBS_STOP begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.
      DURING_IMPRV_OBS = .FALSE.
  
      ! Get time in midwest
      DATE = GET_TIME_AHEAD( - 7 * 60 )

      ! Check if it's midnight
      IF ( DATE(2) == 00 ) THEN 

         ! Check if its a day where we have observations
         ! for 200201
         IF ( MOD( DATE(1) - 20060701+2, 3 ) == 0 ) THEN 
         ! for 200104
         !IF ( MOD( DATE(1) - 20010400 + 1, 3 ) == 0 ) THEN 
         ! for 200107
         !IF ( MOD( DATE(1) - 20010700 + 2, 3 ) == 0 .and.
!     &        DATE(1) .ne. 20010701                       ) THEN
         ! for 200507
         !IF ( MOD( DATE(1) - 20050700 + 2, 3 ) == 0 .and.
     &    !    DATE(1) .ne. 20050701                       ) THEN

            FLAG = .TRUE. 

            WRITE(6,*) ' ITS_TIME_FOR_IMPRV_OBS_STOP at ', 
     &                 DATE(1), DATE(2)

            ! Forward calculation 
            IF ( DIRECTION > 0 ) THEN 
 
               ! leaving measument period
               DURING_IMPRV_OBS = .FALSE.
 
            ! Adjoint calculation  
            ELSEIF ( DIRECTION < 0 ) THEN 
               
               ! entereing new measurement period
               DURING_IMPRV_OBS = .TRUE.

            ENDIF 

         ENDIF

      ENDIF 

      print*, 'TIME_FOR_IMPRV_OBS_STOP: DATE = ', DATE, 
     &        DURING_IMPRV_OBS

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_IMPRV_OBS_STOP

!------------------------------------------------------------------------------

      SUBROUTINE RESET_AEROAVE(  )
!
!******************************************************************************
!  Subroutine RESET_AEROAVE 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  NOTES:
!  - called to calculate the 24-hour average model concentrations,  
!      AVE_BC and AVE_OC.  You will need to modify to reference BCPI
!      and BCPO instead of BC / OC / NH4
!
!******************************************************************************
!

#     include "CMN_SIZE"    ! Size params

      ! Local variables
      INTEGER :: I, J

      !=================================================================
      ! RESET_AEROAVE begins here!
      !=================================================================
    
      WRITE(6,*) ' RESET_AEROSAVE '

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J   )    
      DO I = 1, IIPAR
      DO J = 1, JJPAR

         AVE_BCPI(I,J) = 0d0
         !AVE_OCPI(I,J) = 0d0
         AVE_BCPO(I,J) = 0d0
         !AVE_OCPO(I,J) = 0d0

 
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE RESET_AEROAVE

!------------------------------------------------------------------------------

      SUBROUTINE UPDATE_AEROAVE( BCPI, BCPO)
!
!******************************************************************************
!  Subroutine UPDATE_AEROAVE 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,     ONLY : AIRVOL
      USE TIME_MOD,    ONLY : GET_TS_CHEM

#     include  "CMN_SIZE"   ! IIPAR, JJPAR

      ! Arguments
      TYPE (XPLEX),  INTENT(IN)   :: BCPI(IIPAR,JJPAR)
      !TYPE (XPLEX),  INTENT(IN)   :: OCPI(IIPAR,JJPAR)
      TYPE (XPLEX),  INTENT(IN)   :: BCPO(IIPAR,JJPAR)
      !TYPE (XPLEX),  INTENT(IN)   :: OCPO(IIPAR,JJPAR)

    
      ! Local variables
      TYPE (XPLEX)                :: FACTOR
      INTEGER               :: I, J
    
      !=================================================================
      ! UPDATE_AEROAVE begins here!
      !=================================================================
      
      WRITE(6,*) ' UPDATE_AEROSAVE ' 

      ! Percent of a day for a given chemistry timestep multiplied by
      ! uBC conversion factor (kg --> ug)
      FACTOR = GET_TS_CHEM() / ( 1440d0 ) * 1d9 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )    
      DO I = 1, IIPAR
      DO J = 1, JJPAR

         ! Update average, and convert from ug/box to ug/m3 
         AVE_BCPI(I,J) = AVE_BCPI(I,J) + BCPI(I,J) 
     &* FACTOR / AIRVOL(I,J,1)
       !  AVE_OCPI(I,J) = AVE_OCPI(I,J) + OCPI(I,J) 
      !&* FACTOR / AIRVOL(I,J,1)
         AVE_BCPO(I,J) = AVE_BCPO(I,J) + BCPO(I,J) 
     &* FACTOR / AIRVOL(I,J,1)
      !   AVE_OCPO(I,J) = AVE_OCPO(I,J) + OCPO(I,J) 
      !&* FACTOR / AIRVOL(I,J,1)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO



      ! Return to calling program
      END SUBROUTINE UPDATE_AEROAVE

!------------------------------------------------------------------------------
      SUBROUTINE CALC_IMPRV_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_IMPRV_FORCE 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ARG (TYPE) : Description                          [uBC]
!     
!  NOTES:
!  - This is where the cost function actually gets calculated and the  
!    adjoint variables given values.  The cost function is the sum of the  
!    error weighted squared residuals.
!   - This routine should be called from within geos_chem_adj_mod
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,  ONLY  : IT_IS_NAN, ERROR_STOP
      USE TIME_MOD,   ONLY  : GET_TIME_AHEAD
      USE TIME_MOD,   ONLY  : GET_NYMD, GET_NYMDb
      USE TIME_MOD,   ONLY  : GET_NHMS, GET_NHMSb


#     include "CMN_SIZE"    ! IIPAR, JJPAR

      ! Arguments
      TYPE (XPLEX)             :: COST_FUNC
   
      ! Parameters
      TYPE (XPLEX), PARAMETER  :: OBS_REP_UNC = xplex(0.3d0,0d0) 

      ! Local variables 
      TYPE (XPLEX)             :: NEW_COST(IIPAR,JJPAR)
      TYPE (XPLEX), SAVE       :: JSAVE(IIPAR,JJPAR,5)
      INTEGER, SAVE      :: NEXCD(IIPAR,JJPAR)
      TYPE (XPLEX)             :: DIFF, OBS_ERRCOV
      INTEGER            :: I, J, DATE(2), YYYYMMDD
      LOGICAL, SAVE      :: FIRST = .TRUE. 
      TYPE (XPLEX), SAVE       :: USA_MASK(IIPAR,JJPAR)
      TYPE (XPLEX), PARAMETER  :: THRESH = xplex(0d0,0d0)
      !TYPE (XPLEX)             :: THRESH_2 

      !=================================================================
      ! CALC_IMPRV_FORCE begins here!
      !=================================================================

      ! IBC NEW_COST 
      NEW_COST(:,:) = 0d0 

      DATE = GET_TIME_AHEAD( -8 * 60 ) 

      YYYYMMDD = DATE(1) 

      ! Read BPCH file
      CALL READ_IMPRV_BPCH( YYYYMMDD )

      ! Read improve checkpt
      CALL READ_AEROAVE_FILE( YYYYMMDD )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, DIFF, OBS_ERRCOV )    
      DO I = 1, IIPAR 
      DO J = 1, JJPAR 
 
   
         !------------------------------------------
         ! BC
         !------------------------------------------
         ! Only include points above detection limit
         !IF  ( IMPRV_BC(I,J,1) > 1d-3 ) THEN 
         IF  ( IMPRV_BC(I,J,1) > IMPRV_BC(I,J,2) ) THEN 
 
            ! Difference between predicted and observed daily average 
            DIFF = AVE_BCPI(I,J) + AVE_BCPO(I,J) - IMPRV_BC(I,J,1)

            ! Calculate obs error (30% representation + reported err)
!            OBS_ERRCOV = ( OBS_REP_UNC * IMPRV_BC(I,J,1)
!     &                 + SQRT( IMPRV_BC(I,J,3) ) ) **2
            OBS_ERRCOV =  OBS_REP_UNC * IMPRV_BC(I,J,1) 
     &                 +  IMPRV_BC(I,J,3) 
!            ! Use representational error calculated with 4x5 vs 2x2.5
!            OBS_ERRCOV = IMPRV_BC(I,J,4) ** 2 + IMPRV_BC(I,J,3) ** 2

            ! Calculate new cost 
            NEW_COST(I,J)  = NEW_COST(I,J)
     &                     + 0.5d0 * DIFF ** 2 / OBS_ERRCOV
 
            ! Calculate adjoint forcing 
            ADJ_AVE_BCPI(I,J) = DIFF / OBS_ERRCOV
            ADJ_AVE_BCPO(I,J) = DIFF / OBS_ERRCOV

         ENDIF 

         !------------------------------------------
         ! OC
         !------------------------------------------
         ! Only include points above detection limit
         !IF  ( IMPRV_OC(I,J,1) > 1d-3 ) THEN 
         !IF  ( IMPRV_OC(I,J,1) >  IMPRV_OC(I,J,2) ) THEN 
 
            ! Difference between predicted and observed daily average 
          !  DIFF = AVE_OCPI(I,J) + AVE_OCPO(I,J) - IMPRV_OC(I,J,1)

            ! Calculate obs error (30% representation + reported err)
!            OBS_ERRCOV = ( OBS_REP_UNC * IMPRV_OC(I,J,1)
!     &                 + SQRT( IMPRV_OC(I,J,3) ) ) **2
           ! OBS_ERRCOV = OBS_REP_UNC * IMPRV_OC(I,J,1)
      !&                 + IMPRV_OC(I,J,3) 

            ! Calculate new cost 
      !      NEW_COST(I,J)  = NEW_COST(I,J)
      !&                     + 0.5d0 * DIFF ** 2 / OBS_ERRCOV
 
            ! Calculate adjoint forcing 
      !      ADJ_AVE_OCPI(I,J) = DIFF / OBS_ERRCOV
      !      ADJ_AVE_OCPO(I,J) = DIFF / OBS_ERRCOV
       !  ENDIF 

         !------------------------------------------
  
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! dkh debug
      print*, 'AVE_BCPI = ', AVE_BCPI(17:23,34)
      print*, 'IMPRV_BC = ', IMPRV_BC(17:23,34,1)
      print*, 'IMPRV_BC2 = ', IMPRV_BC(17:23,34,2)
      print*, 'IMPRV_BC3 = ', IMPRV_BC(17:23,34,3)
      print*, ' new cost = ', NEW_COST(17:23,34)
      print*, 'ADJ_AVE_BCPI = ', ADJ_AVE_BCPI(17:23,34)

      ! Update cost function 
      WRITE(6,*) ' CALC_IMPRV_FORCE: NEW_COST = ', SUM( NEW_COST(:,:))

      COST_FUNC = COST_FUNC + SUM(NEW_COST(:,:))

      ! Error check
      IF ( IT_IS_NAN( COST_FUNC ) ) THEN 
         CALL ERROR_STOP( 'COST_FUNC IS NaN', 'calc_imprv_force')
      ENDIF
     
      IF ( GET_NYMD() == GET_NYMDb() ) THEN 
         CALL MAKE_JSAVE_FILE( GET_NYMD(), JSAVE, COST_FUNC, NEXCD )
      ENDIF

 
      ! Return to calling program
      END SUBROUTINE CALC_IMPRV_FORCE

!------------------------------------------------------------------------------
      SUBROUTINE ADJ_UPDATE_AEROAVE( ADJ_BCPI, ADJ_BCPO )
!
!******************************************************************************
!  Subroutine ADJ_UPDATE_AEROAVE applies the adjoint of the average 
!  concentrations corresponding to IMPROVE measurements to the adjoint tracer
!  of BC and OC. (dkh, 11/19/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ADJ_BC (TYPE (XPLEX)) : Adjoint of BC aerosol
!  (2 ) ADJ_OC (TYPE (XPLEX)) : Adjoint of OC aerosol
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) ADJ_BC (TYPE (XPLEX)) : Adjoint of BC aerosol
!  (2 ) ADJ_OC (TYPE (XPLEX)) : Adjoint of OC aerosol
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,     ONLY : AIRVOL
      USE TIME_MOD,    ONLY : GET_TS_CHEM

#     include  "CMN_SIZE"   ! IIPAR, JJPAR

      ! Arguments
      TYPE (XPLEX),  INTENT(INOUT)   :: ADJ_BCPI(IIPAR,JJPAR)
      !TYPE (XPLEX),  INTENT(INOUT)   :: ADJ_OCPI(IIPAR,JJPAR)
      TYPE (XPLEX),  INTENT(INOUT)   :: ADJ_BCPO(IIPAR,JJPAR)
      !TYPE (XPLEX),  INTENT(INOUT)   :: ADJ_OCPO(IIPAR,JJPAR)
    
      ! Local variables
      TYPE (XPLEX)                   :: FACTOR
      INTEGER                  :: I, J
    
      !=================================================================
      ! ADJ_UPDATE_AEROAVE begins here!
      !=================================================================
      
      WRITE(6,*) ' ADJ_UPDATE_AEROSAVE ' 
      ! dkh debug
      print*, 'ADJ_BCPI before  = ', ADJ_BCPI(17:23,34)

      ! Percent of a day for a given chemistry timestep multiplied by
      ! uBC conversion factor (kg --> ug)
      FACTOR = GET_TS_CHEM() / ( 1440d0 ) * 1d9 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )    
      DO I = 1, IIPAR
      DO J = 1, JJPAR


         ! fwd code:
         !AVE_BC(I,J) = AVE_BC(I,J) + BC(I,J) * FACTOR / AIRVOL(I,J,1)
         !AVE_OC(I,J) = AVE_OC(I,J) + OC(I,J) * FACTOR / AIRVOL(I,J,1)
         ADJ_BCPI(I,J) = ADJ_BCPI(I,J) 
     &                + ADJ_AVE_BCPI(I,J) * FACTOR / AIRVOL(I,J,1) 
         ADJ_BCPO(I,J) = ADJ_BCPO(I,J) 
     &                + ADJ_AVE_BCPO(I,J) * FACTOR / AIRVOL(I,J,1)
      !    ADJ_OCPI(I,J) = ADJ_OCPI(I,J) 
      !&                + ADJ_AVE_OCPI(I,J) * FACTOR / AIRVOL(I,J,1)
      !    ADJ_OCPO(I,J) = ADJ_OCPO(I,J) 
      !&                + ADJ_AVE_OCPO(I,J) * FACTOR / AIRVOL(I,J,1)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ADJ_UPDATE_AEROAVE
!------------------------------------------------------------------------------

      SUBROUTINE ADJ_RESET_AEROAVE(  )
!
!******************************************************************************
!  Subroutine ADJ_RESET_AEROAVE resets the adjoint of the average aerosol 
!  concentrations that correspond to IMPROVE measurements. (dkh, 11/19/06)  
!
!  NOTES:
!
!******************************************************************************
!

#     include "CMN_SIZE"    ! Size params

      ! Local variables
      INTEGER :: I, J

      !=================================================================
      ! ADJ_RESET_AEROAVE begins here!
      !=================================================================
    
      WRITE(6,*) ' ADJ_RESET_AEROSAVE '

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J   )    
      DO I = 1, IIPAR
      DO J = 1, JJPAR

         ADJ_AVE_BCPI(I,J) = 0d0
       !  ADJ_AVE_OCPI(I,J) = 0d0
         ADJ_AVE_BCPO(I,J) = 0d0
        ! ADJ_AVE_OCPO(I,J) = 0d0
 
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ADJ_RESET_AEROAVE

!------------------------------------------------------------------------------
!      SUBROUTINE READ_USA_MASK( USA_MASK )
!
!******************************************************************************
!  Subroutine READ_USA_MASK reads the USA mask from disk.   The USA mask is
!  the fraction of the grid box (I,J) which lies w/in the continental USA.
!  (rch, bmy, 11/10/04, 10/3/05)
!   - just for diagnostic; you don't need this
!      
!  NOTES:
!  (1 ) Now can read data for GEOS and GCAP grids (bmy, 8/16/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! Reference to F90 modules
      !USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_RES_EXT
!      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      !USE DIRECTORY_MOD, ONLY : DATA_DIR
      !USE TRANSFER_MOD,  ONLY : TRANSFER_2D

!#     include "CMN_SIZE"  ! Size parameters
!#     include "CMN_SETUP" ! DATA_DIR

      ! Local variables
!     TYPE (XPLEX)             :: ARRAY(IGLOB,JGLOB,1)
!      TYPE (XPLEX)             :: XTAU
!      TYPE (XPLEX)             :: USA_MASK(IGLOB,JGLOB)
!      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_USA_MASK begins here!
      !=================================================================

      ! File name
!      FILENAME = TRIM( DATA_DIR )           //
!     &           'EPA_NEI_200411/usa_mask.' // GET_NAME_EXT_2D() //
!     &           'EPA_NEI_200411/usa_mask.geos' //
!     &           '.'                        // GET_RES_EXT()

      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_USA_MASK: Reading ', a )
      
      ! Get TAU0 for Jan 1985 
!      XTAU  = GET_TAU0( 1, 1, 1985 )
      
      ! USA mask is stored in the bpch file as #2
!      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
!     &                 XTAU,      IGLOB,    JGLOB,
!     &                 1,         ARRAY,    QUIET=.TRUE. )
      
      ! Cast to TYPE (XPLEX)
      !CALL TRANSFER_2D( ARRAY(:,:,1), USA_MASK )
!      USA_MASK(:,:) = ARRAY(:,:,1)
      
      ! Return to calling program
!      END SUBROUTINE READ_USA_MASK

!------------------------------------------------------------------------------
      SUBROUTINE MAKE_JSAVE_FILE( YYYYMMDD, JSAVE, COST_FUNC, NEXCD )
!
!******************************************************************************
!  Subroutine MAKE_JSAVE_FILE saves daily average concentrations of 
!  OC and BC aerosol. (dkh, 11/18/06)  
!  - another diagnostic.  It saves the residuals to bpch file for plotting
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : Date of average                          [uBC]
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,    ONLY : GET_TAU,     EXPAND_DATE

#     include "CMN_SIZE"    ! Size params
!#     include "CMN_ADJ"     ! ADJ_DIR
      
      ! Arguments          
      INTEGER              :: YYYYMMDD
      TYPE (XPLEX)               :: JSAVE(IIPAR,JJPAR,5)
      TYPE (XPLEX)               :: COST_FUNC
      INTEGER              :: NEXCD(IIPAR,JJPAR)
      
      ! Local variables    
      INTEGER              :: I, J, I0, J0, HHMMSS_dum, L
      CHARACTER(LEN=120)   :: FILENAME
      TYPE (XPLEX)               :: DAT(IIPAR,JJPAR,5)
      
      ! For binary punch file, version 2.0
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UBC     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      !=================================================================
      ! MAKE_JSAVE_FILE begins here!
      !=================================================================
      
      FILENAME = TRIM( 'jsave.YYYYMMDD' )
      
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum)
      
      ! Append data directory prefix
      FILENAME = '/qb6/yhmao/geos-chem/adjoint/new/gcadj_std/obsdata/'//
     &           TRIM( FILENAME )
      
      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'Average Aerosol data file '
      CATEGORY = 'IJ-AVE-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      UBC     = '%'

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
 100  FORMAT( '     - MAKE_JSAVE_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_JSAVE, FILENAME, TITLE )

      IF ( COST_FUNC > 0d0 ) THEN
         ! Temporarily store data in DAT as REAL4
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, 5
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( L < 5 ) THEN

               IF ( NEXCD(I,J) > 0 ) THEN
                  DAT(I,J,L) = JSAVE(I,J,L) / XPLX(NEXCD(I,J))
               ELSE
                  DAT(I,J,L) = 0
               ENDIF

            ELSE

               DAT(I,J,L) = JSAVE(I,J,L) / COST_FUNC * 100d0

            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE

         print*, 'COST FUNCTION IS ZERO'
         print*, 'COST FUNCTION IS ZERO'
         print*, 'COST FUNCTION IS ZERO'
         DAT = 0d0

      ENDIF

      CALL BPCH2( IU_JSAVE,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UBC,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     5,         I0+1,
     &            J0+1,      1,         DAT )

      ! Close file
      CLOSE( IU_JSAVE )

      ! Return to calling program
      END SUBROUTINE MAKE_JSAVE_FILE

!-----------------------------------------------------------------------------

      SUBROUTINE INIT_IMPROVE 
!
!*****************************************************************************
!  Subroutine INIT_IMPROVE deallocates all module arrays. (dkh, 11/16/06)  
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
      ! INIT_IMPROVE begins here
      !================================================================= 
      ALLOCATE( IMPRV_BC( IIPAR, JJPAR, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IMPRV_BC' )
      IMPRV_BC = 0d0

      !ALLOCATE( IMPRV_OC( IIPAR, JJPAR, 4 ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'IMPRV_OC' )
      !IMPRV_OC = 0d0

      ALLOCATE( AVE_BCPI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVE_BCPI' )
      AVE_BCPI = 0d0

      !ALLOCATE( AVE_OCPI( IIPAR, JJPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVE_OCPI' )
      !AVE_OCPI = 0d0
      
      ALLOCATE( AVE_BCPO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVE_BCPO' )
      AVE_BCPO = 0d0

      !ALLOCATE( AVE_OCPO( IIPAR, JJPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVE_OCPO' )
      !AVE_OCPO = 0d0

      ALLOCATE( ADJ_AVE_BCPI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ADJ_AVE_BCPI' )
      ADJ_AVE_BCPI = 0d0

      !ALLOCATE( ADJ_AVE_OCPI( IIPAR, JJPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'ADJ_AVE_OCPI' )
      !ADJ_AVE_OCPI = 0d0

      ALLOCATE( ADJ_AVE_BCPO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ADJ_AVE_BCPO' )
      ADJ_AVE_BCPO = 0d0

      !ALLOCATE( ADJ_AVE_OCPO( IIPAR, JJPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'ADJ_AVE_OCPO' )
      !ADJ_AVE_OCPO = 0d0

      ! Return to calling program 
      END SUBROUTINE INIT_IMPROVE
!------------------------------------------------------------------------------
      SUBROUTINE CLEANUP_IMPROVE 
!
!*****************************************************************************
!  Subroutine CLEANUP_IMPROVE deallocates all module arrays. (dkh, 11/16/06)  
!        
!  NOTES:
!
!******************************************************************************
!     
      IF ( ALLOCATED( IMPRV_BC  ) )      DEALLOCATE(  IMPRV_BC)
      !IF ( ALLOCATED( IMPRV_OC  ) )      DEALLOCATE(  IMPRV_OC)
      IF ( ALLOCATED( AVE_BCPI  ) )      DEALLOCATE(  AVE_BCPI)
      !IF ( ALLOCATED( AVE_OCPI  ) )      DEALLOCATE(  AVE_OCPI)
      IF ( ALLOCATED( AVE_BCPO  ) )      DEALLOCATE(  AVE_BCPO)
      !IF ( ALLOCATED( AVE_OCPO  ) )      DEALLOCATE(  AVE_OCPO)
      IF ( ALLOCATED( ADJ_AVE_BCPI  ) )      DEALLOCATE(  ADJ_AVE_BCPI)
      !IF ( ALLOCATED( ADJ_AVE_OCPI  ) )      DEALLOCATE(  ADJ_AVE_OCPI)
      IF ( ALLOCATED( ADJ_AVE_BCPO  ) )      DEALLOCATE(  ADJ_AVE_BCPO)
      !IF ( ALLOCATED( ADJ_AVE_OCPO  ) )      DEALLOCATE(  ADJ_AVE_OCPO)
      
      ! Return to calling program 
      END SUBROUTINE CLEANUP_IMPROVE
!------------------------------------------------------------------------------
      END MODULE IMPROVE_BC_MOD

